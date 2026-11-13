"""Trace every dedicated pLN no-match to the measured conditions of its five lineage rules."""
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from multiplex5000.tissues.lymph_node import load_rule_source

def marker_name(name):
    return "PD_L1" if name == "PDL1" else name.replace("-", "_").replace(" ", "_")


def audit_section(donor, source, counts, evidence, *, calls_dir, dedicated_root, cluster_root, gate_table):
    entry = donor + "__pln_1"
    status = pd.read_parquet(calls_dir / f"{donor}.parquet",
                             columns=["analysis_eligible", "candidate_state"])
    thresholds = gate_table.loc[gate_table.protein_entry.eq(entry)].set_index("marker")
    features = [m for m in thresholds.index if marker_name(m) not in {"DAPI","Beta_actin"}]
    stage = dedicated_root / "06_typing" / "protein" / f"{entry}.parquet"
    calls = pq.ParquetFile(stage).read(columns=["object_id"] + ["marker_call:" + m for m in features]).to_pandas().set_index("object_id")
    if not status.index.is_unique or not status.index.equals(calls.index):
        raise ValueError(f"{donor}: unmatched native IDs")
    admitted = status.analysis_eligible.to_numpy()
    unmatched = admitted & status.candidate_state.eq("no_match").to_numpy()
    cells = int(unmatched.sum())
    census = counts.loc[counts.donor.eq(donor) & counts.level.eq("detail")]
    if census.loc[census.label.eq("No rule match"), "cells"].item() != cells or census.cells.sum() != admitted.sum():
        raise ValueError(f"{donor}: active donor census differs from native calls")
    markers = {m: calls.loc[unmatched, "marker_call:" + m].to_numpy() for m in features}
    native = dedicated_root / "03_protein_measurement" / "cells" / f"{entry}.parquet"
    data = pq.ParquetFile(native).read(columns=["object_id"] + ["marker:" + m for m in features]).to_pandas().set_index("object_id")
    if not data.index.equals(calls.index):
        raise ValueError(f"{donor}: unmatched raw and candidate native IDs")
    key = {"donor": donor, "disease": census.disease.unique().item()}
    root_rows = []
    has_complete_positive = np.zeros(cells, dtype=bool)
    has_anchor = np.zeros(cells, dtype=bool)
    root_markers = set()
    for root in source.candidate_lineage_roots:
        rule = source.by_phenotype[root]
        mandatory = [m for m, op in rule.constraints if op in {"pos", "allpos"}]
        alternatives = [m for m, op in rule.constraints if op == "anypos"]
        excludes = [m for m, op in rule.constraints if op == "allneg"]
        root_markers.update(mandatory + alternatives)
        for marker in mandatory + alternatives + excludes:
            if not np.isin(markers[marker], [0, 1]).all():
                raise ValueError(f"{donor}: no-match root includes unavailable {marker}")
        complete = np.ones(cells, dtype=bool)
        if mandatory:
            complete &= np.logical_and.reduce([markers[m] == 1 for m in mandatory])
        if alternatives:
            complete &= np.logical_or.reduce([markers[m] == 1 for m in alternatives])
        veto = np.logical_or.reduce([markers[m] == 1 for m in excludes])
        if (complete & ~veto).any():
            raise ValueError(f"{donor}: no-match satisfies full {root} root")
        anchors = np.logical_or.reduce([markers[m] == 1 for m in mandatory + alternatives])
        has_anchor |= anchors
        has_complete_positive |= complete
        root_rows.append(dict(**key, root=root, no_match=cells, required_positive="; ".join(mandatory),
                             any_positive="; ".join(alternatives), required_negative="; ".join(excludes),
                             any_root_positive=int(anchors.sum()), complete_positive_blocked=int(complete.sum()),
                             incomplete_positive=int((~complete).sum()), complete_positive_no_veto=0))
    marker_rows = []
    expected = evidence.loc[evidence.donor.eq(donor) & evidence.level.eq("detail")
                            & evidence.label.eq("No rule match")].set_index("marker")
    for marker in features:
        threshold = thresholds.loc[marker]
        raw = data.loc[unmatched, "marker:" + marker].to_numpy(dtype=float)
        aligned = raw * threshold.shift_scale
        score = (aligned - threshold.floor) / (threshold.top - threshold.floor)
        positive = markers[marker] == 1
        measured = markers[marker] >= 0
        name = marker_name(marker)
        row = expected.loc[name]
        if not np.isclose(threshold.cutoff, threshold.floor + .1 * (threshold.top - threshold.floor)):
            raise ValueError(f"{donor}: cutoff differs for {marker}")
        if not (np.array_equal(np.isfinite(raw) & (raw >= 0), measured)
                and np.array_equal(aligned[measured] >= threshold.cutoff, positive[measured])
                and positive.sum() == row.n_positive and measured.sum() == row.n_measured):
            raise ValueError(f"{donor}: threshold or summary mismatch for {marker}")
        marker_rows.append(dict(**key, marker=name, no_match=cells, measured=int(measured.sum()),
                                positive=int(positive.sum()), mean_scaled=row.mean_scaled,
                                floor=threshold.floor, top=threshold.top, cutoff=threshold.cutoff,
                                shift_scale=threshold.shift_scale, cutoff_fraction=.1,
                                at_or_above_0_20=int((score[measured] >= .2).sum()),
                                positive_root_marker=marker in root_markers))
    donor_row = dict(**key, admitted=int(admitted.sum()), no_match=cells, no_match_pct=100*cells/admitted.sum(),
                     complete_positive_blocked=int(has_complete_positive.sum()),
                     partial_positive_only=int((has_anchor & ~has_complete_positive).sum()),
                     no_root_positive=int((~has_anchor).sum()))
    cluster = pd.read_parquet(cluster_root / f"donor_id={donor}/roi=pln_1/data.parquet",
                              columns=["object_id", "is_immune", "immune_class", "cell_agrees_with_cluster_class",
                                       "call_source"]).set_index("object_id")
    if not cluster.index.is_unique or not cluster.call_source.eq("inherited from this cell's cluster").all():
        raise ValueError(f"{donor}: prior cluster identities or label basis failed verification")
    cluster_class = cluster.immune_class.astype("string")
    cluster_class = cluster_class.mask(cluster.is_immune.eq(False).fillna(False), "Non-immune cluster")
    cluster_class = cluster_class.mask(cluster.is_immune.isna(), "Undecided cluster")
    cluster_class = cluster_class.replace({"Uncalled": "Immune cluster, class uncalled"})
    if cluster_class.fillna("").eq("").any():
        raise ValueError(f"{donor}: cluster status missing")
    prior = pd.DataFrame({"cluster_class": cluster_class,
                          "cell_agrees": cluster.cell_agrees_with_cluster_class})
    comparison = prior.reindex(status.index[unmatched])
    comparison["cluster_class"] = comparison.cluster_class.fillna("Outside prior cluster analysis")
    comparison["named_immune"] = comparison.cluster_class.isin(
        ["B_Plasma", "T_cell", "Myeloid", "Granulocyte", "NK_ILC"])
    comparison["agrees_named"] = comparison.cell_agrees.eq(True) & comparison.named_immune
    comparison["measured_agreement"] = comparison.cell_agrees.notna() & comparison.named_immune
    history = comparison.groupby("cluster_class", observed=True).agg(
        cells=("named_immune", "size"), named_immune=("named_immune", "sum"),
        cell_agrees_with_named_cluster=("agrees_named", "sum"),
        evaluated_for_named_cluster_agreement=("measured_agreement", "sum")).reset_index()
    if history.cells.sum() != cells:
        raise ValueError(f"{donor}: prior cluster comparison lost native cells")
    return donor_row, root_rows, marker_rows, history.assign(**key)


def main(data_dir, calls_dir, dedicated_root, cluster_root, manifest):
    gate_table = pd.read_csv(dedicated_root/"06_typing/protein_thresholds.csv")
    source = load_rule_source(manifest, dedicated=True)
    counts = pd.read_parquet(data_dir / "donor_counts.parquet").query("tissue == 'pLN'")
    evidence = pd.read_parquet(data_dir / "donor_markers.parquet").query("tissue == 'pLN'")
    donors, roots, markers, clusters = [], [], [], []
    for donor in sorted(counts.donor.unique()):
        summary, rules, thresholds, history = audit_section(
            donor, source, counts, evidence, calls_dir=calls_dir, dedicated_root=dedicated_root,
            cluster_root=cluster_root, gate_table=gate_table)
        donors.append(summary)
        roots.extend(rules)
        markers.extend(thresholds)
        clusters.extend(history.to_dict("records"))
    for filename, rows in [("pln_no_match_by_donor", donors), ("pln_no_match_by_root", roots),
                           ("pln_no_match_marker_qc", markers), ("pln_no_match_prior_clusters", clusters)]:
        out = data_dir / (filename + ".parquet")
        temporary = out.with_suffix(".parquet.tmp")
        pd.DataFrame(rows).to_parquet(temporary, index=False)
        temporary.replace(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","calls-dir","dedicated-root","cluster-root","manifest"]:
        parser.add_argument("--"+name,type=Path,required=True)
    main(**vars(parser.parse_args()))
