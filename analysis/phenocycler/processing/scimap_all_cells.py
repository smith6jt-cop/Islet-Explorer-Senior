"""Score the assessable dedicated pLN hierarchy with unmodified SCIMAP 2.3.5."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import distribution
from importlib.util import module_from_spec, spec_from_file_location
from functools import partial
from pathlib import Path
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from multiplex5000.tissues.lymph_node import load_rule_source

DRAW, SEED = 5000, 20261113


ALT = ["scimap_candidate_root", "scimap_candidate_detail", "scimap_evidence_status"]
def scorer():
    installed = distribution("scimap")
    if installed.version != "2.3.5":
        raise ValueError("Use the original SCIMAP 2.3.5 scorer")
    path = installed.locate_file("scimap/tools/phenotype_cells.py")
    spec = spec_from_file_location("scimap_235_native_phenotype", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.phenotype_cells


def root_name(source, label):
    rules = source.by_phenotype
    while rules[label].parent in rules:
        label = rules[label].parent
    return label


METHOD = ("SCIMAP 2.3.5 phenotype_cells with the five dedicated pLN lineage branches, in source sheet order. "
          "Donor-specific calibrated floor/0.10 cutoff/99.9th percentile map linearly to 0/.5/1. "
          "EpCAM is excluded. CD38-dependent descendants are unassessable. "
          "These are whole-cell score-based candidate calls. No manual slide gates were transferred.")


def build_sheet(source, available):
    eligible_rules = [r for r in source.rules if root_name(source, r.phenotype) in source.candidate_lineage_roots]
    rules = [r for r in eligible_rules if source.inherited_markers(r.phenotype) <= available]
    sheet = pd.DataFrame([{"CellType":r.parent,"Phenotype":r.phenotype, **dict(r.constraints)} for r in rules])
    if "EpCAM" in sheet.columns:
        raise ValueError("Dedicated cell type rules unexpectedly include EpCAM")
    return sheet, eligible_rules, rules


def classify_chunk(raw, index, thresholds, sheet, score):
    valid = np.isfinite(raw).all(axis=1) & (raw >= 0).all(axis=1)
    result = pd.Series("Unassessable", index=index, dtype="string")
    if not valid.any():
        return result
    scaled = np.column_stack([np.interp(raw[valid, i]*t.shift_scale, [t.floor, t.cutoff, t.top], [0, .5, 1])
                              for i,t in enumerate(thresholds.itertuples())])
    ids = index[valid]
    obj = ad.AnnData(scaled, obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=thresholds.index))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
        score(obj, sheet.copy(), gate=.5, label="alternative", verbose=False)
    result.loc[ids] = obj.obs.alternative.astype("string")
    return result


def check_scoring(manifest):
    """Exercise the native hierarchy and continuous scorer using synthetic cell measurements."""
    source = load_rule_source(manifest, dedicated=True)
    features = set().union(*(r.markers for r in source.rules)) - {"CD38"}
    sheet, _, _ = build_sheet(source, features)
    if "Naive_B_cell" in set(sheet.Phenotype) or "B_cell" not in set(sheet.Phenotype):
        raise AssertionError("Absent CD38 must remove descendant evidence and retain the B parent")
    markers = list(sheet.columns[2:])
    gates = pd.DataFrame({"floor":0, "cutoff":.1, "top":1, "shift_scale":1},index=markers)
    sample = pd.DataFrame(.08,index=["B","T","veto","missing","weak"],columns=markers)
    sample.loc["B", ["CD20","CD79a"]] = 1
    sample.loc["T", ["CD3e","CD8","TCF-1"]] = 1
    sample.loc["veto", ["CD3e","Pan-Cytokeratin"]] = [1,.101]
    sample.loc["missing","CD3e"] = np.nan
    score = scorer()
    result = classify_chunk(sample.to_numpy(),sample.index,gates,sheet,score)
    assert result.tolist() == ["B_cell","CD8_T_cytotoxic","DN_T_cell","Unassessable","Nerve_fiber"]
    # SCIMAP can call a type from averaged group scores even if no positive marker reaches its gate.
    # The sensitivity deliberately retains that behavior; its candidate label alone is not proof of positivity.
    root_sheet = sheet.loc[sheet.CellType.eq("all")].dropna(axis=1,how="all")
    roots = list(root_sheet.columns[2:])
    result = classify_chunk(sample[roots].to_numpy(),sample.index,gates.loc[roots],root_sheet,score)
    assert result.tolist() == ["Immune","Immune","Immune","Unassessable","Neural"]


def classify(donor, *, data_dir, calls_dir, dedicated_root, manifest, parts):
    function = scorer()
    source = load_rule_source(manifest, dedicated=True)
    entry = donor + "__pln_1"
    gates = pd.read_csv(dedicated_root / "06_typing/protein_thresholds.csv")
    gates = gates.loc[gates.protein_entry.eq(entry)].set_index("marker")
    available = set(gates.index[gates.cutoff.notna()])
    sheet, eligible_rules, rules = build_sheet(source, available)
    markers = list(sheet.columns[2:])
    thresholds = gates.loc[markers]
    data = pq.ParquetFile(dedicated_root/"03_protein_measurement"/"cells"/f"{entry}.parquet").read(
        columns=["object_id", *["marker:"+m for m in markers]],use_threads=False).to_pandas().set_index("object_id")
    current = pd.read_parquet(calls_dir/f"{donor}.parquet",
                              columns=["analysis_eligible","candidate_state","candidate_phenotype"])
    current = current.loc[current.analysis_eligible, ["candidate_state","candidate_phenotype"]].copy()
    data = data.loc[current.index]
    if not data.index.is_unique or not data.index.equals(current.index):
        raise ValueError(f"{donor}: row identity mismatch")
    out = current.copy()
    out["scimap_candidate_detail"] = "Unassessable"
    out["scimap_evidence_status"] = "Incomplete measurements"
    mapping = {r.phenotype:root_name(source,r.phenotype) for r in eligible_rules}
    children = {r.parent for r in eligible_rules}
    unavailable = {r.parent for r in eligible_rules if r not in rules}
    subset_by_state = [x.sample(min(len(x), DRAW), random_state=SEED+int(donor)).index for _, x in
                       current.groupby("candidate_state", observed=True, sort=True)]
    sampled_ids = pd.Index(np.concatenate(subset_by_state))
    for begin in range(0,len(data),40000):
        chunk = data.iloc[begin:begin+40000]
        labels = classify_chunk(chunk.to_numpy(dtype=float), chunk.index, thresholds, sheet, function)
        base_labels = labels.str.removeprefix("likely-")
        if not base_labels.isin(set(mapping) | {"Unknown","Unassessable"}).all():
            raise ValueError(f"{donor}: unexpected SCIMAP output")
        status = pd.Series("Scored subtype", index=chunk.index, dtype="string")
        status.loc[base_labels.isin(children)] = "Scored parent"
        status.loc[base_labels.isin(unavailable)] = "Scored parent; child marker unavailable"
        status.loc[labels.str.startswith("likely-")] = "Likely candidate"
        status.loc[labels.eq("Unknown")] = "No broad identity scored"
        status.loc[labels.eq("Unassessable")] = "Incomplete measurements"
        out.loc[labels.index,"scimap_candidate_detail"] = labels
        out.loc[status.index,"scimap_evidence_status"] = status
    out["scimap_candidate_root"] = out.scimap_candidate_detail.str.removeprefix("likely-").map(
        {**mapping,"Unknown":"Unknown","Unassessable":"Unassessable"})
    source_root = current.candidate_phenotype.astype("string").map(mapping).fillna(current.candidate_state.astype("string"))
    replay = out.loc[sampled_ids].groupby(["candidate_state","scimap_candidate_root"],observed=True).size()
    validation = pd.read_parquet(data_dir/"pln_no_match_scimap_sensitivity.parquet",
                                 filters=[("donor","=",donor)]).replace({"scimap_root":{"Unassessable root marker":"Unassessable"}})
    baseline = validation.groupby(["candidate_state","scimap_root"],observed=True).sampled_cells.sum()
    left,right = replay.astype("int64").align(baseline.astype("int64").rename_axis(replay.index.names),fill_value=0)
    if not np.array_equal(left,right):
        raise ValueError(f"{donor}: full hierarchy failed independently scored sampled roots")
    out["source_root"] = source_root
    out["donor"] = donor
    parts.mkdir(parents=True,exist_ok=True)
    final = parts / (donor+".parquet")
    temp = final.with_suffix(".building")
    out.reset_index().to_parquet(temp,compression="zstd",index=False)
    temp.replace(final)
    print(donor,len(out),flush=True)


def assemble(donors, *, data_dir, parts):
    ad.settings.allow_write_nullable_strings = True
    outpath = data_dir/"pln_scimap_candidate_cells.parquet"
    obj = ad.read_h5ad(data_dir/"pln_umap.h5ad")
    schema = pa.schema([(c,pa.string()) for c in ["donor","object_id","candidate_state","candidate_phenotype",
                                                    "source_root",*ALT]])
    summary, shown = [], []
    with pq.ParquetWriter(outpath.with_suffix(".building"),schema,compression="zstd") as writer:
        for donor in donors:
            frame = pd.read_parquet(parts/(donor+".parquet")).set_index("object_id")
            review = pd.read_parquet(data_dir/"pln_lineage_review_cells.parquet",filters=[("donor","=",donor)],
                                     columns=["object_id","pln_lineage_review","own_compartment_root_status","candidate_state"]).set_index("object_id")
            if (not review.index.equals(frame.index) or not frame.index.is_unique or
                    not np.array_equal(frame.candidate_state.astype(str), review.candidate_state.astype(str))):
                raise ValueError(f"{donor}: strict/alternate/review membership differs")
            writer.write_table(pa.Table.from_pandas(frame.reset_index(),schema=schema,preserve_index=False))
            shown.append(frame.loc[obj.obs_names[obj.obs.donor.astype(str).eq(donor)], ALT])
            join = frame.join(review.drop(columns="candidate_state"))
            keys = ["donor","candidate_state","candidate_phenotype","source_root",*ALT,"pln_lineage_review","own_compartment_root_status"]
            result = join.groupby(keys,observed=True,dropna=False).size().rename("cells").reset_index()
            if result.cells.sum()!=len(frame):
                raise ValueError(f"{donor}: missing sample in counts")
            summary.append(result)
    annotations = pd.concat(shown).reindex(obj.obs_names)
    if annotations.isna().any().any():
        raise ValueError("Partial UMAP join")
    for c in ALT:
        obj.obs[c]=annotations[c].astype("category")
    obj.uns["pln_scimap_candidate_method"] = METHOD
    obj.write_h5ad(data_dir/"pln_umap.full.building.h5ad",compression="lzf")
    metrics = data_dir/"pln_scimap_candidate_counts.parquet"
    pd.concat(summary,ignore_index=True).to_parquet(metrics.with_suffix(".building"),index=False)
    outpath.with_suffix(".building").replace(outpath)
    metrics.with_suffix(".building").replace(metrics)
    (data_dir/"pln_umap.full.building.h5ad").replace(data_dir/"pln_umap.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--donor", action="append")
    parser.add_argument("--assemble-only",action="store_true")
    for name in ["data-dir","calls-dir","dedicated-root","manifest","parts"]:
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    check_scoring(args.manifest)
    donors = sorted(f.stem for f in args.calls_dir.glob("*.parquet"))
    run = partial(classify, data_dir=args.data_dir, calls_dir=args.calls_dir, dedicated_root=args.dedicated_root,
                  manifest=args.manifest, parts=args.parts)
    if not args.assemble_only:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(run,args.donor or donors))
    if args.donor is None:
        assemble(donors, data_dir=args.data_dir, parts=args.parts)
