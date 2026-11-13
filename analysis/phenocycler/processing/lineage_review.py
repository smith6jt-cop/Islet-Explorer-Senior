"""Attach broad, independent review evidence to strict pLN calls by native cell ID."""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq



BROAD = {"T_cell": "T cell", "B_Plasma": "B/plasma", "Myeloid": "Myeloid",
         "Granulocyte": "Granulocyte", "NK_ILC": "NK/ILC"}
STATUS = "pln_lineage_review_status"
LABEL = "pln_lineage_review"
PRIOR = ["prior_cluster_id", "prior_cluster_class", "prior_cell_top_subpanel",
         "prior_cluster_fraction_agreeing", "prior_cell_agrees"]
COMP = ["own_compartment_root_status", "own_compartment_single_immune_root"]
ROOTS = ["immune", "endothelial", "neural", "stromal", "smooth_muscle"]
METHOD = ("Only admitted strict no-match cells receive a broad immune review label. "
          "October 1 cluster-inherited immune class and the cell's highest immune subpanel must agree. "
          "The subpanel uses one representative marker per immune panel chosen in that cluster. "
          "The reported cluster fraction is agreement among evaluable historical cluster cells. "
          "Own-compartment gate replay tests the broad immune parent separately on raw measurements. "
          "No class is validated by an independent ground-truth label.")


def section(donor, compartment_root, metrics, *, calls_dir, cluster_root):
    current = pd.read_parquet(calls_dir / f"{donor}.parquet",
                              columns=["analysis_eligible", "candidate_state", "candidate_phenotype"])
    current = current.loc[current.analysis_eligible, ["candidate_state", "candidate_phenotype"]].copy()
    if not current.index.is_unique:
        raise ValueError(f"{donor}: repeated native IDs")
    current["donor"] = donor
    prior = pd.read_parquet(cluster_root / f"donor_id={donor}/roi=pln_1/data.parquet",
                            columns=["object_id", "cluster_id", "is_immune", "immune_class", "cell_argmax_subpanel",
                                     "cell_agrees_with_cluster_class"]).set_index("object_id")
    if not prior.index.is_unique:
        raise ValueError(f"{donor}: repeated cluster IDs for a cell")
    matched = current.join(prior, validate="one_to_one")
    matched["prior_cluster_id"] = matched.cluster_id.astype("string")
    matched["prior_cluster_fraction_agreeing"] = matched.prior_cluster_id.map(metrics.loc[donor])
    cluster = matched.immune_class.astype("string")
    named = cluster.isin(BROAD)
    matched["prior_cluster_class"] = cluster.mask(matched.is_immune.eq(False).fillna(False), "Non-immune cluster")
    matched["prior_cluster_class"] = matched.prior_cluster_class.mask(
        matched.is_immune.isna() & matched.prior_cluster_id.notna(), "Undecided cluster")
    matched["prior_cell_top_subpanel"] = matched.cell_argmax_subpanel.astype("string").replace("",pd.NA)
    matched["prior_cell_agrees"] = matched.cell_agrees_with_cluster_class.astype("boolean").where(named)
    agrees = matched.prior_cell_agrees.fillna(False)
    check = matched.prior_cell_top_subpanel.eq(cluster).fillna(False)
    if not agrees[named].equals(check[named].astype("boolean")):
        raise ValueError(f"{donor}: historical cell-level agreement was not reproduced")
    no_match = matched.candidate_state.eq("no_match")
    matched[LABEL] = cluster.map(BROAD).where(no_match & named & agrees)
    matched[STATUS] = matched.candidate_state.astype("string").map({
        "unique_candidate": "Strict unique candidate", "ambiguous": "Strict overlapping candidates",
        "unavailable": "Strict evidence unavailable", "no_match": "No named immune cluster"})
    matched.loc[no_match & named & ~agrees, STATUS] = "Cell and named cluster disagree"
    matched.loc[no_match & named & matched.prior_cell_agrees.isna(), STATUS] = "Cell subpanel unavailable"
    matched.loc[no_match & matched.prior_cluster_id.isna(), STATUS] = "Outside prior cluster analysis"
    matched.loc[matched[LABEL].notna(), STATUS] = "Cell and named cluster agree"
    evidence = pd.read_parquet(compartment_root / f"{donor}.parquet").set_index("object_id")
    if (not evidence.index.is_unique or not evidence.donor.astype("string").eq(donor).all()
            or len(evidence) != no_match.sum() or not matched.index[no_match].isin(evidence.index).all()):
        raise ValueError(f"{donor}: compartment audit has different no-match cells")
    roots = evidence[[f"root_{r}_own_raw" for r in ROOTS]]
    if not np.array_equal(roots.sum(axis=1), evidence.own_raw_roots):
        raise ValueError(f"{donor}: compartment root counts disagree")
    complete = evidence[[f"root_{r}_own_all_markers_measured" for r in ROOTS]].all(axis=1)
    own_status = pd.Series("No root passed", index=evidence.index, dtype="string")
    for root in ROOTS:
        only = evidence.own_raw_roots.eq(1) & evidence[f"root_{root}_own_raw"]
        own_status.loc[only] = "Single " + root.replace("_", " ") + " root"
    own_status.loc[evidence.own_raw_roots.gt(1)] = "Multiple roots passed"
    own_status.loc[~complete] = "Incomplete root measurements"
    matched["own_compartment_root_status"] = own_status.reindex(matched.index)
    matched["own_compartment_single_immune_root"] = own_status.eq("Single immune root").where(complete).reindex(matched.index)
    return matched[["donor","candidate_state","candidate_phenotype",LABEL,STATUS,*PRIOR,*COMP]]


def main(data_dir, calls_dir, cluster_root, compartment_root):
    ad.settings.allow_write_nullable_strings = True
    metrics = pd.read_csv(cluster_root / "cluster_heterogeneity.csv",
                          dtype={"donor_id":str, "cluster_id":str})
    metrics = metrics.set_index(["donor_id","cluster_id"]).cell_frac_agrees_class
    data_file = data_dir / "pln_lineage_review_cells.parquet"
    cell_temp = data_file.with_suffix(".parquet.tmp")
    h5ad_file = data_dir / "pln_umap.h5ad"
    obj = ad.read_h5ad(h5ad_file)
    qc = pd.read_parquet(data_dir / "pln_no_match_by_donor.parquet").set_index("donor")
    prior_tallies = pd.read_parquet(data_dir / "pln_no_match_prior_clusters.parquet").groupby("donor").cell_agrees_with_named_cluster.sum()
    points, summaries = [], []
    schema = pa.schema([(col, pa.string()) for col in ["object_id","donor","candidate_state","candidate_phenotype",LABEL,STATUS,*PRIOR[:3]]]
                       + [("prior_cluster_fraction_agreeing",pa.float32()),("prior_cell_agrees",pa.bool_()),
                          (COMP[0],pa.string()),(COMP[1],pa.bool_())])
    with pq.ParquetWriter(cell_temp, schema, compression="zstd") as writer:
        for donor, expected in qc.iterrows():
            frame = section(donor, compartment_root, metrics, calls_dir=calls_dir, cluster_root=cluster_root)
            if (len(frame) != expected.admitted or frame[LABEL].notna().sum() != prior_tallies[donor]
                    or frame[STATUS].isna().any()):
                raise ValueError(f"{donor}: full or review class counts differ")
            sample = obj.obs.loc[obj.obs.donor.astype("string").eq(donor)]
            selected = frame.loc[sample.index]
            if not selected.candidate_state.astype("string").equals(sample.candidate_state.astype("string")):
                raise ValueError(f"{donor}: native UMAP candidate status differs")
            points.append(selected[[LABEL,STATUS,*PRIOR,*COMP]])
            out = frame.reset_index()
            writer.write_table(pa.Table.from_pandas(out,schema=schema,preserve_index=False))
            summary = frame.groupby(["candidate_state",STATUS,LABEL,COMP[0]], dropna=False,observed=True).size().rename("cells").reset_index()
            summaries.append(summary.assign(donor=donor))
    annotations = pd.concat(points).loc[obj.obs_names]
    if not obj.obs_names.is_unique or not annotations.index.equals(obj.obs_names):
        raise ValueError("UMAP annotations do not align")
    for col in annotations:
        if col == "prior_cluster_fraction_agreeing":
            obj.obs[col] = annotations[col].astype("float32")
        elif col in ("prior_cell_agrees",COMP[1]):
            obj.obs[col] = annotations[col].astype("boolean")
        else:
            obj.obs[col] = annotations[col].astype("category")
    obj.uns["pln_lineage_review_method"] = METHOD
    summary_file = data_dir / "pln_lineage_review_counts.parquet"
    count_temp = summary_file.with_suffix(".parquet.tmp")
    pd.concat(summaries,ignore_index=True).to_parquet(count_temp,index=False)
    h5ad_temp = h5ad_file.with_suffix(".building.h5ad")
    obj.write_h5ad(h5ad_temp,compression="lzf")
    check = ad.read_h5ad(h5ad_temp,backed="r")
    if not check.obs_names.equals(obj.obs_names) or not check.obs[STATUS].equals(obj.obs[STATUS]):
        raise ValueError("Written UMAP review annotations failed verification")
    check.file.close()
    cell_temp.replace(data_file)
    count_temp.replace(summary_file)
    h5ad_temp.replace(h5ad_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","calls-dir","cluster-root","compartment-root"]:
        parser.add_argument("--"+name,type=Path,required=True)
    main(**vars(parser.parse_args()))
