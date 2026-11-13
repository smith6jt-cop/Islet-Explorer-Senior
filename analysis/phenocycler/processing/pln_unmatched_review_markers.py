"""Summarize full-panel pLN evidence for the no-match review or the SCIMAP candidate hierarchy."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

SUPPORTED = ("T cell", "B/plasma", "Myeloid", "Granulocyte")
GROUP_KEYS = ["tissue", "donor", "disease", "review_group", "review_status"]


def _hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _marker_name(name):
    return "PD_L1" if name == "PDL1" else name.replace("-", "_").replace(" ", "_")


def _groups(cells):
    lineage = cells["pln_lineage_review"].astype("string")
    reason = cells["pln_lineage_review_status"].astype("string")
    if reason.isna().any() or not set(lineage.dropna()) <= set(SUPPORTED):
        raise ValueError("Unexpected lineage review value or missing evidence status")
    return pd.DataFrame({"review_group": lineage.fillna(reason), "review_status": reason}, index=cells.index)


def _keys(frame):
    result = frame if frame.index.name == "object_id" else frame.set_index("object_id")
    if not result.index.is_unique or result.index.hasnans:
        raise ValueError("Each native object_id must be unique and present in its section")
    return result


def _donor(job):
    (donor, disease, root, calls_dir, data, raw_hash, calls_hash, reference, thresholds) = job
    formal = _keys(pq.ParquetFile(calls_dir/f"{donor}.parquet").read(
        columns=["object_id", "analysis_eligible", "candidate_state"]).to_pandas())
    unmatched = formal.index[formal.analysis_eligible & formal.candidate_state.eq("no_match")]
    cells = pd.read_parquet(data/"pln_lineage_review_cells.parquet", filters=[("donor", "=", donor),
                            ("candidate_state", "=", "no_match")], columns=["object_id", "pln_lineage_review",
                            "pln_lineage_review_status"])
    cells = _keys(cells)
    if len(cells) != len(unmatched) or not cells.index.isin(unmatched).all():
        raise ValueError(f"{donor}: reviewed object IDs differ from the full formal no-match census")
    cells = cells.loc[unmatched]
    grouping = _groups(cells)
    return _summarize(donor, disease, root, raw_hash, calls_hash, formal, grouping, reference, thresholds)


def _summarize(donor, disease, root, raw_hash, calls_hash, formal, grouping, reference, thresholds):
    filename = f"{donor}__pln_1.parquet"
    raw_path = root / "03_protein_measurement" / "cells" / filename
    markers_path = root / "06_typing" / "protein" / filename
    if _hash(raw_path) != raw_hash or _hash(markers_path) != calls_hash:
        raise ValueError(f"{donor}: dedicated pLN measurements do not match the typing receipt")
    row_codes, labels = pd.factorize(pd.MultiIndex.from_frame(grouping), sort=True)
    quantities = labels.set_names(grouping.columns.tolist()).to_frame(index=False)
    quantities["cells"] = np.bincount(row_codes, minlength=len(labels))
    quantities = quantities.assign(tissue="pLN", donor=donor, disease=disease)
    ranks = thresholds.set_index("marker")
    markers = [m for m in ranks.index if _marker_name(m) in reference.index]
    if len(markers) != len(reference) or len(reference) != 53:
        raise ValueError(f"{donor}: unmatched reference lacks an active protein or calibration")
    raw = _keys(pq.ParquetFile(raw_path).read(columns=["object_id"]+["marker:"+m for m in markers]).to_pandas())
    states = _keys(pq.ParquetFile(markers_path).read(columns=["object_id"]+["marker_call:"+m for m in markers]).to_pandas())
    if not raw.index.equals(formal.index) or not states.index.equals(formal.index):
        raise ValueError(f"{donor}: native measurement or marker-call identity mismatch")
    pos = formal.index.get_indexer(grouping.index)
    if (pos < 0).any():
        raise ValueError(f"{donor}: a grouped cell lacks a native measurement row")
    outputs = []
    for m in markers:
        name = _marker_name(m)
        t = ranks.loc[m]
        calibrated = raw["marker:"+m].to_numpy()[pos].astype(float) * t.shift_scale
        call = states["marker_call:"+m].to_numpy()[pos]
        valid = call>=0
        if not np.isin(call,[-1,0,1]).all() or not np.isfinite(calibrated[valid]).all():
            raise ValueError(f"{donor}/{m}: invalid measured marker evidence")
        if not (t.top>t.floor and np.isclose(t.cutoff, t.floor+.1*(t.top-t.floor))):
            raise ValueError(f"{donor}/{m}: invalid dedicated section calibration")
        if not np.array_equal(calibrated[valid]>=t.cutoff, call[valid]==1):
            raise ValueError(f"{donor}/{m}: marker positivity does not reproduce the dedicated stage")
        scaled = np.clip((calibrated[valid]-t.floor)/(t.top-t.floor),0,1).astype(np.float32)
        count = np.bincount(row_codes[valid],minlength=len(labels))
        positive = np.bincount(row_codes[call==1],minlength=len(labels))
        signal = np.bincount(row_codes[valid],weights=scaled,minlength=len(labels))
        expected = reference.loc[name]
        if (count.sum()!=expected.n_measured or positive.sum()!=expected.n_positive
                or not np.isclose(signal.sum()/count.sum(),expected.mean_scaled,rtol=1e-7,atol=1e-8)):
            raise ValueError(f"{donor}/{m}: candidate groups fail to reconstruct the original marker profile")
        mean = np.full(len(labels),np.nan)
        fraction = np.full(len(labels),np.nan)
        np.divide(signal,count,out=mean,where=count>0)
        np.divide(positive,count,out=fraction,where=count>0)
        outputs.append(quantities.assign(marker=name,n_measured=count,n_positive=positive,
                                        mean_scaled=mean,positive_fraction=fraction))
    return quantities,pd.concat(outputs,ignore_index=True)


def _scimap_donor(job):
    donor, disease, root, calls_dir, data, raw_hash, calls_hash, reference, thresholds = job
    formal = _keys(pq.ParquetFile(calls_dir/f"{donor}.parquet").read(
        columns=["object_id", "analysis_eligible", "candidate_state"]).to_pandas())
    admitted = formal.index[formal.analysis_eligible]
    keys = ["scimap_candidate_root", "scimap_candidate_detail", "scimap_evidence_status"]
    source = _keys(pd.read_parquet(data/"pln_scimap_candidate_cells.parquet",
                   filters=[("donor","=",donor)], columns=["object_id","candidate_state"]+keys))
    if len(source)!=len(admitted) or not source.index.isin(admitted).all():
        raise ValueError(f"{donor}: cells differ from formal admitted census")
    source = source.loc[admitted]
    if (not np.array_equal(source.candidate_state.astype(str),formal.loc[admitted,"candidate_state"].astype(str))
            or source[keys].isna().any().any()):
        raise ValueError(f"{donor}: identity or status mismatch")
    support, stats = _summarize(donor,disease,root,raw_hash,calls_hash,formal,source[keys],reference,thresholds)
    support["eligible_cells"] = len(admitted)
    support["fraction_tissue"] = support.cells/len(admitted)
    support["root_cells"] = support.groupby("scimap_candidate_root",observed=True).cells.transform("sum")
    support["fraction_root"] = support.cells/support.root_cells
    return support, stats


def _run(data,root,calls_dir,official,worker,workers):
    receipt_path = root/"06_typing"/"receipt.json"
    receipt = json.loads(receipt_path.read_text())
    sources = {rec["name"]: rec["sha256"] for rec in receipt["inputs"]}
    saved = {rec["path"]: rec["sha256"] for rec in receipt["outputs"]}
    thresholds_path = root/"06_typing"/"protein_thresholds.csv"
    if _hash(thresholds_path)!=saved[thresholds_path.name]:
        raise ValueError("Section marker thresholds disagree with their typing receipt")
    thresholds = pd.read_csv(thresholds_path)
    jobs = []
    for d,group in official.groupby("donor", observed=True):
        entry = f"{d}__pln_1"
        jobs.append((str(d),group.disease.unique().item(),root,calls_dir,data,
                     sources[f"protein_cells:{entry}"],saved[f"protein/{entry}.parquet"],
                     group.set_index("marker"),thresholds[thresholds.protein_entry.eq(entry)]))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(worker,jobs))
    counts = pd.concat([part[0] for part in results],ignore_index=True)
    stats = pd.concat([part[1] for part in results],ignore_index=True)
    return counts,stats,receipt_path,thresholds_path


def build_review_marker_tables(*,data_dir,dedicated_root,calls_dir,workers=3):
    """Rebuild after publishing the full review-cell table; inputs and outputs use native UUIDs."""
    data = Path(data_dir)
    root = Path(dedicated_root)
    calls_dir = Path(calls_dir)
    official = pd.read_parquet(data/"donor_markers.parquet")
    official = official[(official.tissue=="pLN") & (official.level=="detail")
                        & (official.label=="No rule match") & (official.marker!="EpCAM")]
    official_counts = pd.read_parquet(data/"donor_counts.parquet")
    official_counts = official_counts[(official_counts.tissue=="pLN") & (official_counts.level=="detail")
                                     & (official_counts.label=="No rule match")].set_index("donor")
    counts,stats,receipt_path,thresholds_path = _run(data,root,calls_dir,official,_donor,workers)
    observed = counts.groupby("donor",observed=True).cells.sum().sort_index()
    expected = official_counts.cells.sort_index()
    if not observed.index.equals(expected.index) or not np.array_equal(observed,expected):
        raise ValueError("Review-marker output differs from full donor no-match census")
    if stats.duplicated(GROUP_KEYS+["marker"]).any():
        raise ValueError("Repeated donor, review group and marker")
    metrics = ["cells","marker","n_measured","n_positive","mean_scaled","positive_fraction"]
    outputs = {"pln_unmatched_review_markers.parquet":stats[GROUP_KEYS+metrics],
               "pln_unmatched_review_counts.parquet":counts[GROUP_KEYS+["cells"]]}
    with tempfile.TemporaryDirectory(prefix=".pln-review-markers-", dir=data) as work:
        staging = Path(work)
        for filename, frame in outputs.items():
            frame.to_parquet(staging/filename,index=False,compression="zstd")
        method = {"formal_state":"no_match","measurements":"Dedicated pLN aligned protein, section-specific floor and top, clipped to [0,1]. Positivity uses the original dedicated threshold and counts measured cells only.",
                  "review":"Current corrected candidate state with prior cluster and matching own lineage evidence. A null supported lineage uses its recorded unresolved reason.",
                  "n_donors":counts.donor.nunique(),"unmatched_cells":int(counts.cells.sum()),"active_markers":int(stats.marker.nunique()),
                  "producer_sha256":_hash(__file__),"dedicated_receipt_sha256":_hash(receipt_path),
                  "thresholds_sha256":_hash(thresholds_path),"native_review_table_sha256":_hash(data/"pln_lineage_review_cells.parquet"),
                  "outputs":{name:_hash(staging/name) for name in outputs}}
        method_file = "pln_unmatched_review_marker_method.json"
        (staging/method_file).write_text(json.dumps(method,indent=2)+"\n")
        for name in [*outputs,method_file]:
            (staging/name).replace(data/name)
    return method


def build_scimap_marker_tables(*,data_dir,dedicated_root,calls_dir,workers=3):
    """Measure all admitted cells by SCIMAP candidate, with independent pooled reference checks."""
    data,root,calls_dir = map(Path,(data_dir,dedicated_root,calls_dir))
    baseline = pd.read_parquet(data/"donor_markers.parquet")
    baseline = baseline[(baseline.tissue=="pLN") & (baseline.level=="detail") & (baseline.marker!="EpCAM")].copy()
    baseline["signal"] = baseline.mean_scaled*baseline.n_measured
    original = baseline.groupby(["donor","disease","marker"],observed=True).agg(
        n_measured=("n_measured","sum"),n_positive=("n_positive","sum"),signal=("signal","sum")).reset_index()
    original["mean_scaled"] = original.signal/original.n_measured
    counts,stats,receipt_path,thresholds_path = _run(data,root,calls_dir,original,_scimap_donor,workers)
    status_rows = pd.read_parquet(data/"donor_counts.parquet")
    unavailable = status_rows[(status_rows.tissue=="pLN") & (status_rows.level=="detail") & ~status_rows.assessable]
    # These same five B-cell descendants require CD38, which is unavailable for the
    # dedicated panel. A denominator may be known even when its subtype count is unknown.
    if set(unavailable.label) != {"Naive_B_cell","Memory_B_cell","GC_B_cell","Plasma_cell","Plasmablast"}:
        raise ValueError("Re-evaluate the SCIMAP missing-marker rule coverage for this release")
    unavailable = unavailable[["tissue","donor","disease","label","eligible_cells"]].rename(columns={"label":"scimap_candidate_detail"})
    if set(unavailable.scimap_candidate_detail) & set(counts.scimap_candidate_detail):
        raise ValueError("A candidate declared unavailable has observed calls")
    # The available classifier vocabulary is the union of its observed output classes.
    # Cross-check that union against every admitted, assessable source-rule category.
    registered=set(status_rows.loc[(status_rows.tissue=="pLN") & (status_rows.level=="detail") & status_rows.assessable,"label"])
    registered-= {"No rule match","Unavailable candidate","Ambiguous candidate"}
    if registered-set(counts.scimap_candidate_detail):
        raise ValueError("A registered assessable hierarchy category lacks a verified SCIMAP root")
    class_keys=["scimap_candidate_root","scimap_candidate_detail","scimap_evidence_status"]
    roster=counts[class_keys].drop_duplicates()
    if roster.scimap_candidate_detail.duplicated().any():
        raise ValueError("A SCIMAP candidate label has conflicting roots or evidence status")
    sections=counts[["tissue","donor","disease","eligible_cells"]].drop_duplicates()
    root_sizes=counts.groupby(["donor","scimap_candidate_root"],observed=True).cells.sum().rename("root_cells")
    complete=sections.merge(roster,how="cross").merge(
        counts[["donor",*class_keys,"cells"]], on=["donor",*class_keys], how="left", validate="one_to_one")
    complete["cells"]=complete.cells.fillna(0).astype("int64")
    counts=complete.join(root_sizes,on=["donor","scimap_candidate_root"])
    counts["root_cells"]=counts.root_cells.fillna(0).astype("int64")
    counts["fraction_tissue"]=counts.cells/counts.eligible_cells
    counts["fraction_root"]=counts.cells/counts.root_cells.where(counts.root_cells.gt(0))
    counts["assessable"] = True
    unavailable = unavailable.assign(scimap_candidate_root="Immune",scimap_evidence_status="Unassessable without CD38",
                                    cells=np.nan, fraction_tissue=np.nan, fraction_root=np.nan, assessable=False)
    immune = counts[counts.scimap_candidate_root.eq("Immune")].groupby("donor",observed=True).root_cells.first()
    unavailable["root_cells"] = unavailable.donor.map(immune)
    counts=pd.concat([counts,unavailable],ignore_index=True)
    keys = ["tissue","donor","disease","scimap_candidate_root","scimap_candidate_detail","scimap_evidence_status"]
    if counts.duplicated(keys).any() or stats.duplicated(keys+["marker"]).any():
        raise ValueError("Repeated donor and SCIMAP candidate combination")
    total=counts.groupby("donor",observed=True).cells.sum()
    parent=counts.groupby("donor",observed=True).eligible_cells.first()
    if not total.eq(parent).all():
        raise ValueError("SCIMAP leaf counts fail to partition each admitted section")
    compare=pd.read_parquet(data/"pln_scimap_candidate_counts.parquet").groupby(keys[1:2]+keys[3:],observed=True).cells.sum()
    observed=counts[counts.assessable].set_index(keys[1:2]+keys[3:]).cells
    obs,expected=observed.align(compare,fill_value=0)
    if not np.array_equal(obs,expected):
        raise ValueError("SCIMAP marker strata do not match the independent classification count export")
    all_bins = counts[keys+["cells"]].merge(original[["marker"]].drop_duplicates(),how="cross")
    stats = all_bins.merge(stats.drop(columns="cells"),on=keys+["marker"],how="left",validate="one_to_one")
    for col in ("n_measured","n_positive"):
        stats[col]=stats[col].fillna(0).astype("int64")
    if (stats.n_measured.eq(0) & (stats.mean_scaled.notna() | stats.positive_fraction.notna())).any():
        raise ValueError("A missing marker measurement has been displayed as a measured result")
    outputs = {"pln_scimap_candidate_markers.parquet":stats[keys+["marker","n_measured","n_positive","mean_scaled","positive_fraction","cells"]],
               "pln_scimap_candidate_proportions.parquet":counts[keys+["cells","eligible_cells","fraction_tissue","root_cells","fraction_root","assessable"]]}
    with tempfile.TemporaryDirectory(prefix=".pln-scimap-markers-",dir=data) as scratch:
        stage=Path(scratch)
        for name,frame in outputs.items():frame.to_parquet(stage/name,index=False,compression="zstd")
        method={"cells":int(total.sum()),"donors":len(total),"markers":stats.marker.nunique(),
                "signal":"Same dedicated, calibrated [0,1] measured signal and positivity as the release dotplots. The SCIMAP phenotype scoring uses its own gate mapping.",
                "denominators":"fraction_tissue uses the full admitted donor section. fraction_root uses all outcome classes in the same SCIMAP root. All five CD38-dependent branches remain missing.",
                "source_cells_sha256":_hash(data/"pln_scimap_candidate_cells.parquet"),"formal_markers_sha256":_hash(data/"donor_markers.parquet"),
                "dedicated_receipt_sha256":_hash(receipt_path),"thresholds_sha256":_hash(thresholds_path),"producer_sha256":_hash(__file__),
                "outputs":{name:_hash(stage/name) for name in outputs}}
        note="pln_scimap_candidate_marker_method.json"
        (stage/note).write_text(json.dumps(method,indent=2)+"\n")
        for name in [*outputs,note]:(stage/name).replace(data/name)
    return method


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","dedicated-root","calls-dir"]:
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--workers",type=int,default=3)
    parser.add_argument("--mode",choices=["review","scimap"],default="review")
    args=vars(parser.parse_args())
    run=build_scimap_marker_tables if args.pop("mode")=="scimap" else build_review_marker_tables
    print(run(**args))
