"""Step 1: marker x donor usability over the whole cohort.

Reads the run's expression stage and writes one row per (donor, marker) with the statistics
in ``usability.py`` and a fail-closed verdict.  Decides nothing about cell types; it decides
which measurements the later steps are allowed to lean on.

    python scripts/celltype/run_step1.py                 # from core/phenocycler/
    python scripts/celltype/run_step1.py --donors 6450   # one donor, for a quick check
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import pathlib

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                      # noqa: E402
import usability as U                             # noqa: E402


def donor_frame(donor: str, markers: list) -> pd.DataFrame:
    """Analysis-eligible cells for one donor, complete across the markers IT CARRIES.

    Intersecting with ``panel.donor_markers`` rather than demanding the full list is what
    lets the panel-variant markers in at all: requiring CD38 of every donor would drop all
    cells for the 7 that never acquired it, which is the opposite of "a comparison that uses
    that marker runs on fewer donors".
    """
    markers = [m for m in markers if m in set(panel.donor_markers(donor))]
    cols = ["object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid"] + markers
    # Two separate reasons this line is what it is, found by two sessions on 2026-09-17.
    #
    # panel.expression_path, not a fixed path: it resolves the stable flat root first and falls
    # back to the retired run tree per donor. A donor added after 2026-09-16 (6533, 112) exists
    # only in the flat root, so the retired run's fixed path raised FileNotFoundError for them.
    #
    # The FILE, not the dataset: a partition that also carries ``donor_id`` as a column conflicts
    # with the hive key inferred from its directory (large_string vs dictionary), and dataset
    # discovery refuses to merge the two. Reading the file skips the inference, and is identical
    # for a partition that does not carry the column.
    df = pq.ParquetFile(panel.expression_path(donor)).read(columns=cols).to_pandas()
    n_all = len(df)
    df = df[df["qc_analysis_eligible"]]
    n_elig = len(df)
    df = df[np.isfinite(df[markers].to_numpy("float32")).all(axis=1)]
    df = df[np.isfinite(df[["X_centroid", "Y_centroid"]].to_numpy()).all(axis=1)]
    df.attrs["counts"] = {"rows": n_all, "eligible": n_elig, "complete": len(df)}
    return df.reset_index(drop=True)


def stats_for_donor(donor: str, markers: list, partners: dict, controls: dict) -> list:
    df = donor_frame(donor, markers)
    markers = [m for m in markers if m in df.columns]
    c = df.attrs["counts"]
    print(f"  {donor}: {c['rows']:,} rows -> {c['eligible']:,} eligible -> "
          f"{c['complete']:,} complete", flush=True)
    X = {m: df[m].to_numpy("float64") for m in markers}
    nbr = U.neighbour_graph(df[["X_centroid", "Y_centroid"]].to_numpy(), U.BARS["morans_k"])

    # The spatial null, measured on this donor: probe markers spanning the concentration
    # range, so the bar reflects how heavy-tailed this donor's values actually are.
    conc = sorted(markers, key=lambda m: U.shape_stats(X[m])["top1_signal_share"] or 0)
    probes = [conc[0], conc[len(conc) // 2], conc[-1]][:U.BARS["morans_null_probes"]]
    null = U.morans_null([X[m] for m in probes], nbr,
                         draws=U.BARS["morans_null_draws"], q=U.BARS["morans_null_q"])
    print(f"    spatial null |I| q{U.BARS['morans_null_q']:.2f} = {null:.5f} "
          f"(probes: {', '.join(probes)})", flush=True)

    rows = []
    for m in markers:
        r = {"donor_id": donor, "marker": m,
             "n_eligible": c["eligible"], "n_complete": c["complete"], "morans_null": null}
        r |= U.shape_stats(X[m])
        r |= U.spatial_stats(X[m], nbr)

        # tails, against declared partners and against derived non-partners
        pv = {p: U.tail_overlap(X[m], X[p], U.BARS["tail_q"])
              for p, _ in partners.get(m, []) if p in X}
        cv = {q: U.tail_overlap(X[m], X[q], U.BARS["tail_q"])
              for q in controls.get(m, []) if q in X}
        pv = {k: v for k, v in pv.items() if np.isfinite(v)}
        cv = {k: v for k, v in cv.items() if np.isfinite(v)}
        r["partner_best"] = max(pv.values()) if pv else np.nan
        r["partner_best_with"] = max(pv, key=pv.get) if pv else ""
        r["control_worst"] = max(cv.values()) if cv else np.nan
        r["control_worst_with"] = max(cv, key=cv.get) if cv else ""
        r["partners"] = json.dumps({k: round(v, 2) for k, v in pv.items()})
        r["controls"] = json.dumps({k: round(v, 2) for k, v in cv.items()})
        rows.append(r)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    # The acquisition table describes the 20 cohort donors. A donor acquired outside it --
    # or a checkout without the table -- can still be measured per (donor, marker), which is
    # everything Step 3 leans on; only the acquisition-GROUP statistic needs the cohort, and
    # it is a between-group ratio that a single donor cannot carry anyway.
    try:
        coh = panel.cohort()
    except FileNotFoundError:
        coh = None
    donors = args.donors or (coh["donor_id"].tolist() if coh is not None else [])
    if not donors:
        raise SystemExit("no donor_acquisition_batch.csv on disk, so --donors is required")
    markers = panel.markers(include_variant=True)
    partners, controls = panel.partner_map(), panel.control_map()
    print(f"Step 1 — {len(donors)} donors x {len(markers)} markers, "
          f"tail at top {U.BARS['tail_q']:.1%}, Moran's I k={U.BARS['morans_k']}")

    t0 = time.time()
    rows = []
    for d in donors:
        t = time.time()
        rows += stats_for_donor(d, markers, partners, controls)
        print(f"    {time.time() - t:.0f}s", flush=True)
    df = pd.DataFrame(rows)

    # Acquisition-group structure needs every donor, so it is added after the loop.
    off_table = sorted(set(df["donor_id"]) - set(coh["donor_id"])) if coh is not None else \
        sorted(set(df["donor_id"]))
    if off_table:
        print(f"\n  acquisition-group structure NOT computed: {len(off_table)} donor(s) are "
              f"not in the acquisition table ({off_table}). The per-(donor, marker) "
              f"statistics below are unaffected; only the between-group ratio is skipped.")
    else:
        grp = coh.set_index("donor_id")["acq_group"]
        df = df.merge(
            coh[["donor_id", "batch", "epoch", "panel_variant", "acq_group", "status"]],
            on="donor_id", how="left")
        acq = []
        for m, g in df.groupby("marker"):
            acq.append({"marker": m,
                        **U.acquisition_group_ratio(g.set_index("donor_id")["p50"], grp)})
        df = df.merge(pd.DataFrame(acq), on="marker", how="left")


    out = pathlib.Path(args.out) if args.out else panel.OUT / "marker_donor_usability.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"\nwrote {out}  ({len(df):,} rows, {time.time() - t0:.0f}s)")
    print("\nmarkers whose level steps with the acquisition group "
          f"(ratio >= {U.BARS['acq_ratio_flag']}) — the one measurement Step 2 acts on:")
    r = (df.groupby("marker")["acq_ratio"].first().dropna().sort_values(ascending=False)
         if "acq_ratio" in df.columns else pd.Series(dtype=float))
    if r.empty:
        print("    not computable: the ratio needs at least two acquisition groups")
    else:
        for m, v in r[r >= U.BARS["acq_ratio_flag"]].items():
            print(f"    {m:<18} {v:8.2f}")
        rest = r[r < U.BARS["acq_ratio_flag"]]
        if not rest.empty:
            print(f"    (next highest: {rest.idxmax()} {rest.max():.2f})")

    print("\nweakest declared partner tails, pooled over donors (Step 3's input):")
    q = (df.groupby("marker")[["partner_best", "control_worst"]].median()
           .dropna().sort_values("partner_best"))
    print(q.head(8).to_string(float_format=lambda v: f"{v:9.1f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
