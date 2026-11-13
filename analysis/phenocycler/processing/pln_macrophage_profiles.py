"""Apply the CD68-or-CD163 macrophage anchor and Gal3 subdivision to the reviewed pLN hierarchy."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

try:
    from .pln_unmatched_review_markers import _hash, _keys, _summarize
except ImportError:
    from pln_unmatched_review_markers import _hash, _keys, _summarize

LEAVES = ("Macrophage_Gal3+", "Macrophage_Gal3-")
MAC_PARENTS = {"Macrophage", "M1_macrophage", "M2_macrophage"}
UNKNOWN = "Macrophage candidate · anchor unconfirmed"
GAL_MISSING = "Macrophage · Gal3 unassessed"
QC_STATES = {"ambiguous": "Ambiguous candidate", "unavailable": "Unavailable candidate"}
GROUP = ["tissue", "donor", "disease", "label"]


def _section(job):
    donor, disease, data, root, calls_dir, ref, threshold, raw_hash, calls_hash, sample_ids, parts = job
    formal = _keys(pq.ParquetFile(calls_dir / f"{donor}.parquet").read(
        columns=["object_id", "analysis_eligible", "candidate_state", "candidate_phenotype"]).to_pandas())
    admitted = formal[formal.analysis_eligible]
    source = _keys(pd.read_parquet(data / "pln_lineage_review_cells.parquet", filters=[("donor", "=", donor)],
                       columns=["object_id","candidate_state","pln_lineage_review","pln_lineage_review_status"]))
    if len(source)!=len(admitted) or not admitted.index.isin(source.index).all():
        raise ValueError(f"{donor}: lineage review and admitted census disagree")
    source = source.loc[admitted.index]
    state = admitted.candidate_state.astype("string")
    if not state.equals(source.candidate_state.astype("string")):
        raise ValueError(f"{donor}: candidate states disagree")
    rule = admitted.candidate_phenotype.astype("string")
    lineage = source.pln_lineage_review.astype("string")
    base = rule.where(state.eq("unique_candidate"),state.map(QC_STATES)).astype("string")
    no_match = state.eq("no_match")
    base.loc[no_match] = ("Review · "+lineage).fillna(source.pln_lineage_review_status)[no_match]
    if base.isna().any(): raise ValueError(f"{donor}: missing base leaf")
    markers=["CD68","CD163","M2Gal3","Pan-Cytokeratin"]
    entries = _keys(pq.ParquetFile(root/"06_typing"/"protein"/f"{donor}__pln_1.parquet").read(
        columns=["object_id",*["marker_call:"+m for m in markers]]).to_pandas())
    if not entries.index.equals(formal.index): raise ValueError(f"{donor}: native UUID order differs")
    native=entries.loc[admitted.index].set_axis(markers,axis=1)
    if not native.isin([-1,0,1]).all().all(): raise ValueError("Unexpected marker call")
    anchor=(native.CD68.eq(1) | native.CD163.eq(1))
    anchor_status=np.select([anchor,native[["CD68","CD163"]].lt(0).any(axis=1)], ["positive","unassessed"],default="negative")
    source_mac=base.isin(MAC_PARENTS)
    parent=source_mac | base.isin(["Immune","Myeloid","Review · Myeloid"])
    gate = parent & native["Pan-Cytokeratin"].eq(0)
    assigned=gate & anchor
    result=base.copy()
    result.loc[source_mac & ~assigned]=UNKNOWN
    result.loc[assigned]=GAL_MISSING
    result.loc[assigned & native.M2Gal3.eq(1)]=LEAVES[0]
    result.loc[assigned & native.M2Gal3.eq(0)]=LEAVES[1]
    disposition=np.select([assigned,source_mac & ~anchor, parent & anchor,anchor],
        ["Within checked lineage","Previous macrophage lacks measured anchor","Immune exclusion or missing cytokeratin",
         "Outside accepted myeloid or unresolved immune parent"],default="Outside macrophage call")
    q=pd.DataFrame({"source":base,"anchor_status":anchor_status,"cd68":native.CD68,"cd163":native.CD163,
                    "gal3":native.M2Gal3,"pan_cytokeratin":native["Pan-Cytokeratin"],"disposition":disposition,
                    "label":result}).groupby(["source","anchor_status","cd68","cd163","gal3","pan_cytokeratin","disposition","label"],
                    observed=True,dropna=False).size().rename("cells").reset_index().assign(donor=donor,disease=disease)
    assert q.cells.sum()==len(admitted)
    counts,stats = _summarize(donor,disease,root,raw_hash,calls_hash,formal,result.to_frame("label"),ref,threshold)
    counts["eligible_cells"] = len(admitted)
    sampled=result.loc[sample_ids].rename("pln_macrophage_detail").to_frame().assign(donor=donor).reset_index()
    part=pd.DataFrame({"donor":donor,"object_id":result.index,"pln_macrophage_detail":result.values})
    part.to_parquet(parts/f"{donor}.parquet",index=False,compression="zstd")
    return counts,stats,q,sampled


def build(*,data_dir,dedicated_root,calls_dir,workers=3):
    data,root,calls_dir=map(Path,(data_dir,dedicated_root,calls_dir))
    reference=pd.read_parquet(data/"donor_markers.parquet")
    reference=reference[(reference.tissue=="pLN") & (reference.level=="detail") & (reference.marker!="EpCAM")].copy()
    reference["signal"]=reference.mean_scaled*reference.n_measured
    reference=reference.groupby(["donor","disease","marker"],observed=True)[["n_positive","n_measured","signal"]].sum().reset_index()
    reference["mean_scaled"]=reference.signal/reference.n_measured
    receipt_path=root/"06_typing"/"receipt.json"; receipt=json.loads(receipt_path.read_text())
    raw_pins={r["name"]:r["sha256"] for r in receipt["inputs"]}
    pins={r["path"]:r["sha256"] for r in receipt["outputs"]}
    thr_file=root/"06_typing"/"protein_thresholds.csv"
    if _hash(thr_file)!=pins[thr_file.name]:raise ValueError("Native threshold hash differs")
    threshold=pd.read_csv(thr_file)
    obj=ad.read_h5ad(data/"pln_umap.h5ad",backed="r")
    obs=obj.obs[["donor"]].copy();obj.file.close()
    if obs.index.has_duplicates:raise ValueError("Embedding has duplicate UUIDs")
    scratch=TemporaryDirectory(prefix=".pln-macrophage-cells-",dir=data)
    parts=Path(scratch.name)
    jobs=[]
    for (d,dx),ref in reference.groupby(["donor","disease"],observed=True):
        entry=f"{d}__pln_1"
        jobs.append((d,dx,data,root,calls_dir,ref.set_index("marker"),threshold[threshold.protein_entry.eq(entry)],
            raw_pins["protein_cells:"+entry],pins["protein/"+entry+".parquet"],obs.index[obs.donor.astype(str).eq(str(d))],parts))
    with ProcessPoolExecutor(max_workers=workers) as exec:
        res=list(exec.map(_section,jobs))
    counts,stats,q,sample=[pd.concat(x,ignore_index=True) for x in zip(*res)]
    base=pd.read_parquet(data/"donor_counts.parquet")
    base=base[(base.tissue=="pLN") & (base.level=="detail")]
    unavailable=base[~base.assessable][GROUP+["cells","eligible_cells","assessable","comparable"]]
    donors=counts[GROUP[:3]+["eligible_cells"]].drop_duplicates()
    roster=set(counts.label)|set(LEAVES)|set(base.loc[base.assessable,"label"])-MAC_PARENTS-{"No rule match"}
    complete=donors.merge(pd.DataFrame({"label":sorted(roster)}),how="cross").merge(
        counts[GROUP+["cells"]],on=GROUP,how="left",validate="one_to_one")
    complete["cells"]=complete.cells.fillna(0).astype(int)
    complete=complete.assign(assessable=True,comparable=True)
    counts=pd.concat([complete,unavailable],ignore_index=True)
    counts=counts.assign(level="subtype",parent_cells=counts.eligible_cells,fraction_tissue=counts.cells/counts.eligible_cells)
    counts["fraction_parent"]=counts.fraction_tissue
    totals=counts.groupby("donor",observed=True).cells.sum()
    if not totals.eq(donors.set_index("donor").eligible_cells).all(): raise ValueError("New leaves overlap")
    stats=counts[GROUP+["cells"]].merge(reference[["marker"]].drop_duplicates(),how="cross").merge(
        stats.drop(columns="cells"),on=GROUP+["marker"],how="left",validate="one_to_one").assign(level="subtype")
    for k in ["n_measured","n_positive"]:stats[k]=stats[k].fillna(0).astype(int)
    if len(sample)!=len(obs) or not obs.index.isin(sample.object_id).all():raise ValueError("Embedding UUIDs differ")
    outputs={"pln_macrophage_counts.parquet":counts,"pln_macrophage_markers.parquet":stats,
             "pln_macrophage_transitions.parquet":q,"pln_macrophage_embedding.parquet":sample}
    with TemporaryDirectory(prefix=".pln-macrophage-",dir=data) as tmp:
        stage=Path(tmp)
        for name,frame in outputs.items():frame.to_parquet(stage/name,index=False,compression="zstd")
        cells_file="pln_macrophage_cells.parquet"
        schema=pa.schema([("donor",pa.string()),("object_id",pa.string()),("pln_macrophage_detail",pa.string())])
        with pq.ParquetWriter(stage/cells_file,schema,compression="zstd") as writer:
            for d in donors.donor:
                block=pq.ParquetFile(parts/f"{d}.parquet").read().cast(schema)
                if len(block)!=totals.loc[d]:raise ValueError(f"{d}: exported native census disagrees")
                writer.write_table(block)
        method={"rule":"In an accepted macrophage, unspecialized immune/myeloid, or corroborated review-myeloid parent, CD68 or CD163 positive and Pan-Cytokeratin negative establishes macrophage. Only acquired M2Gal3 determines its two subtypes. Other branches and excluded/uncertain cells retain the original reviewed label. Earlier macrophages without a positive anchor are separately flagged.",
                "marker_call_values":{"-1":"unassessed","0":"negative","1":"positive"},
                "cells":int(totals.sum()),"sample_cells":len(sample),"donors":len(donors),"measured_markers":stats.marker.nunique(),
                "outputs":{n:_hash(stage/n) for n in [*outputs,cells_file]},"producer_sha256":_hash(__file__),
                "source_receipt_sha256":_hash(receipt_path),"source_threshold_sha256":_hash(thr_file),
                "review_cells_sha256":_hash(data/"pln_lineage_review_cells.parquet")}
        method_name="pln_macrophage_method.json";(stage/method_name).write_text(json.dumps(method,indent=2)+"\n")
        for n in [*outputs,cells_file,method_name]:(stage/n).replace(data/n)
    scratch.cleanup()
    return method


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","dedicated-root","calls-dir"]:parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--workers",type=int,default=3)
    print(build(**vars(parser.parse_args())))
