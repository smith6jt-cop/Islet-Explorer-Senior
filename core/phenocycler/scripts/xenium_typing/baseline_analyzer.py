"""Arm X3 ``xenium-analyzer-baseline`` -- the instrument's own onboard clusters, annotated by the null gate.

The Xenium Analyzer (10x onboard analysis) ships ``analysis/clustering/gene_expression_graphclust/
clusters.csv`` (graph-based clusters over all genes) and ``gene_expression_kmeans_<k>_clusters/
clusters.csv`` for k = 2..10.  Those clusters are the baseline the maintainer asked for: they use
every gene (state genes included) by construction, which is recorded on this arm's row; the
ANNOTATION uses identity cores only -- ``gate.nullgate_table`` asks, per cluster, whether a lineage's
own core beats a size- and detection-matched random core, and a cluster is awarded a lineage only
when exactly one lineage clears (two or more -> Ambiguous(multi_lineage); none -> Unavailable).

Primary = graphclust; replicates = the k-means ladder (no seeds exist in the shipped output).  Cells in
the base but absent from the Analyzer table are Unavailable(not_in_analysis); if the analysis folder
has not been transferred, every cell is Unavailable(not_transferred) and the record says so.

    cd core/phenocycler && ../../.venv/bin/python scripts/xenium_typing/baseline_analyzer.py --donors 6436 6442 [--judge]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gate  # noqa: E402
import pins  # noqa: E402

ARM = pins.ARM_BASELINE
KMEANS_KS = tuple(range(2, 11))


def read_clusters(path: pathlib.Path) -> pd.Series:
    df = pd.read_csv(path)
    bc = "Barcode" if "Barcode" in df.columns else df.columns[0]
    cl = "Cluster" if "Cluster" in df.columns else df.columns[1]
    return pd.Series(df[cl].astype(str).to_numpy(), index=df[bc].astype(str))


def annotate_clusters(X, genes, clusters: np.ndarray, cores: dict, *, seed: int = 0) -> tuple:
    """Per-cluster label from the null gate: exactly one lineage clears -> that lineage."""
    long, _calls = gate.nullgate_table(X, genes, clusters, cores, seed=seed)
    label, reason = {}, {}
    for g, sub in long.groupby("group"):
        if not sub["judgeable"].any():
            label[g], reason[g] = gate.UNAVAILABLE, "cluster_below_min_cells"
            continue
        clear = sub.loc[sub["clears"], "lineage"].tolist()
        if len(clear) == 1:
            label[g], reason[g] = clear[0], ""
        elif not clear:
            label[g], reason[g] = gate.UNAVAILABLE, "no_lineage_beats_null"
        else:
            label[g], reason[g] = gate.AMBIGUOUS, "multi_lineage"
    return label, reason, long


def run(donor: str, *, seed: int = 0) -> pathlib.Path:
    t0 = time.time()
    base = gate.load_base_for_gate(donor)
    ids = base["obs_names"]
    adir = pins.analysis_dir(donor)
    gc = adir / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    cores = pins.identity_cores()
    calls = pd.DataFrame({"cell_id": ids, "lineage_native": gate.UNAVAILABLE, "celltype_fine": "",
                          "confidence": np.nan, "reason": "not_transferred", "cycling": pd.array([pd.NA] * len(ids), dtype="boolean")})
    replicates = None
    params = {"primary": "gene_expression_graphclust", "replicates": [f"kmeans_{k}" for k in KMEANS_KS],
              "annotation": "gate.nullgate_table per cluster; exactly one lineage clearing -> label",
              "clusters_use_all_genes": True, "analysis_dir": str(adir)}
    notes = "Xenium Analyzer clusters are computed over all genes (state genes included) by the instrument"
    longs = {}
    if gc.exists():
        s = read_clusters(gc).reindex(ids)
        present = s.notna().to_numpy()
        cl = np.where(present, s.fillna("__absent__").to_numpy(), "__absent__")
        lab, rsn, long = annotate_clusters(base["lognorm"][present], base["var_names"], cl[present], cores, seed=seed)
        longs["graphclust"] = long
        calls.loc[present, "lineage_native"] = [lab[c] for c in cl[present]]
        calls.loc[present, "reason"] = [rsn[c] for c in cl[present]]
        calls.loc[~present, "reason"] = "not_in_analysis"
        params["n_graphclust_clusters"] = int(len(lab))
        reps = {}
        for k in KMEANS_KS:
            pk = adir / "clustering" / f"gene_expression_kmeans_{k}_clusters" / "clusters.csv"
            if not pk.exists():
                continue
            sk = read_clusters(pk).reindex(ids)
            pk_present = sk.notna().to_numpy()
            clk = sk.fillna("__absent__").to_numpy()
            labk, _r, longk = annotate_clusters(base["lognorm"][pk_present], base["var_names"], clk[pk_present], cores, seed=seed)
            longs[f"kmeans_{k}"] = longk
            reps[f"kmeans_{k}"] = np.where(pk_present, [labk.get(c, gate.UNAVAILABLE) for c in clk], gate.UNAVAILABLE)
        if reps:
            replicates = pd.DataFrame({"cell_id": ids, **reps})
    else:
        notes += "; analysis/ not present locally -> every cell Unavailable(not_transferred)"
    called = calls["lineage_native"].isin(pins.LINEAGES)
    rec = gate.ArmRecord(
        arm=ARM, donor_id=str(donor), config_id="primary", identity_inputs={L: len(g) for L, g in cores.items()},
        n_cells=int(len(calls)), n_called=int(called.sum()), n_ambiguous=int((calls["lineage_native"] == gate.AMBIGUOUS).sum()),
        n_unavailable=int((calls["lineage_native"] == gate.UNAVAILABLE).sum()), seconds=time.time() - t0,
        env=sys.executable, stochastic=False, pinned={"base": base["path"], "panel": pins.panel_snapshot()},
        params=params, notes=notes)
    out = gate.write_calls(calls, rec, expected_ids=ids, replicates=replicates)
    d = pins.arm_dir(ARM, donor)
    for name, long in longs.items():
        long.to_parquet(d / f"nullgate__{name}.parquet", index=False)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=["6436", "6442"])
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    for d in a.donors:
        p = run(d, seed=a.seed)
        print(f"{ARM} {d}: wrote {p}")
        if a.judge:
            gate.judge(ARM, d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
