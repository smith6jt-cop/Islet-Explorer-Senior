#!/usr/bin/env python3
"""
Gate-6 accounting audit — prove every canonical cell is assigned exactly once, and show WHERE the
ordered residual actually decides.

Gate 6's exit criterion is "all 25,151,028 eligible canonical cells accounted for exactly once". Nothing
checks that today: ``lineage.py`` writes per-donor partitions and prints percentages, and a cell lost
between ``data/cells`` and ``data/phenotype/broad`` would show up as a slightly different percentage
that nobody could distinguish from biology.

This asserts the accounting and then reports the tree's own behaviour: which gate claimed each cell,
how many cells more than one gate could have claimed (the ordered residual's silent tie-breaks), and
what the `Other` and `Unresolved` populations are actually made of. Descriptive per donor — class
frequency is never a QC target.

    python scripts/gate5/audit_lineage.py [--config config.ini] [--outdir data/gate5]

Exit 1 if the accounting fails. Outputs ``gate6_accounting.csv``, ``gate6_tie_breaks.csv``,
``gate6_audit.png``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                # noqa: E402
from matplotlib.patches import Patch                                           # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenocycler.config import (load_config, COMPARTMENT_COLORS, COMPARTMENT_GATES,   # noqa: E402
                                COMPARTMENT_ORDER, OTHER_LABEL, UNRESOLVED_LABEL)
from build_dashboard import CLASSES, STATUS_LABEL                              # noqa: E402


def audit_donor(cfg, donor: str) -> tuple[dict, pd.DataFrame]:
    """One donor: accounting totals + the per-gate claim/tie-break table."""
    broad = ds.dataset(cfg.broad_dir / f"donor_id={donor}", format="parquet").to_table(
        columns=["object_id", "compartment", "assignment_reason", "blocked_by"]).to_pandas()
    gated = ds.dataset(cfg.restore_gated_dir / f"donor_id={donor}", format="parquet",
                       ignore_prefixes=["apply_manifest", "_", "."]).to_table(
        columns=["object_id"] + [f"{m}_state" for gate in COMPARTMENT_GATES.values() for m in gate]
    ).to_pandas()
    canonical = ds.dataset(cfg.cells_dir / f"donor_id={donor}", format="parquet").to_table(
        columns=["object_id"]).to_pandas()

    for frame in (broad, gated, canonical):
        frame["object_id"] = frame["object_id"].astype(str)

    # --- the accounting -------------------------------------------------------------------------
    n_canonical, n_broad = len(canonical), len(broad)
    dup = int(n_broad - broad["object_id"].nunique())
    canon_ids, broad_ids = set(canonical["object_id"]), set(broad["object_id"])
    unassigned = len(canon_ids - broad_ids)
    extra = len(broad_ids - canon_ids)

    # --- which gate could have claimed each cell (the tie-breaks the order resolves silently) ----
    g = gated.set_index("object_id").reindex(broad["object_id"])
    could_claim = pd.DataFrame(
        {comp: np.zeros(len(g), bool) for comp in COMPARTMENT_ORDER}, index=g.index)
    for comp, gate in COMPARTMENT_GATES.items():
        block = np.column_stack([(g[f"{m}_state"] == "positive").to_numpy() for m in gate])
        could_claim[comp] = block.any(axis=1)
    n_eligible = could_claim.sum(axis=1).to_numpy()

    rows = []
    claimed = broad["compartment"].to_numpy()
    for comp in COMPARTMENT_ORDER:
        eligible = could_claim[comp].to_numpy()
        got = claimed == comp
        rows.append({
            "donor": donor, "compartment": comp,
            "n_eligible": int(eligible.sum()),          # some gate marker positive
            "n_claimed": int(got.sum()),                # actually assigned here
            # cells this gate could have claimed but a HIGHER-priority gate took first
            "n_lost_to_higher": int((eligible & ~got & np.isin(claimed, COMPARTMENT_ORDER)).sum()),
        })
    summary = {
        "donor": donor, "n_canonical": n_canonical, "n_broad": n_broad,
        "n_duplicate_object_id": dup, "n_canonical_unassigned": unassigned,
        "n_broad_not_canonical": extra,
        "n_multi_eligible": int((n_eligible > 1).sum()),
        "n_no_gate_positive": int((n_eligible == 0).sum()),
        "n_other": int((claimed == OTHER_LABEL).sum()),
        "n_unresolved": int((claimed == UNRESOLVED_LABEL).sum()),
        "n_reason_positive": int((broad["assignment_reason"] == "positive").sum()),
        "n_reason_all_negative": int((broad["assignment_reason"] == "all_negative").sum()),
        "n_reason_blocked": int((broad["assignment_reason"] == "blocked_unavailable_gate").sum()),
    }
    print(f"[audit] {donor}: canonical {n_canonical:,} broad {n_broad:,} | "
          f"multi-eligible {summary['n_multi_eligible']:,} "
          f"({100 * summary['n_multi_eligible'] / max(n_broad, 1):.1f}%) | "
          f"Other {summary['n_other']:,} Unresolved {summary['n_unresolved']:,}", flush=True)
    return summary, pd.DataFrame(rows)


def figure(acc: pd.DataFrame, ties: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})

    ax = axes[0]
    agg = ties.groupby("compartment")[["n_eligible", "n_claimed", "n_lost_to_higher"]].sum()
    agg = agg.reindex(COMPARTMENT_ORDER)
    xs = np.arange(len(agg))
    total = acc["n_broad"].sum()
    ax.bar(xs, 100 * agg["n_eligible"] / total, width=0.72, color="#DDDDDD", zorder=2,
           label="a gate marker is positive")
    ax.bar(xs, 100 * agg["n_claimed"] / total, width=0.72, zorder=3, label="assigned to this compartment",
           color=[COMPARTMENT_COLORS[c] for c in agg.index])
    for xi, row in zip(xs, agg.itertuples()):
        lost = 100 * row.n_lost_to_higher / total
        if lost > 0.25:
            ax.text(xi, 100 * row.n_eligible / total + 0.6, f"−{lost:.1f}pp", ha="center", fontsize=8,
                    color="#CC3311")
    ax.set_xticks(xs)
    ax.set_xticklabels(agg.index, rotation=20, ha="right")
    ax.set_ylabel("% of all 25.1M cells")
    ax.set_title("A — what the ordered residual costs each compartment\n"
                 "grey = a gate marker is positive · coloured = actually assigned here\n"
                 "red = lost to a higher-priority gate", fontsize=10.5, loc="left", pad=26)
    ax.set_ylim(0, 72)
    ax.legend(fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.93))

    ax = axes[1]
    frac_multi = 100 * acc["n_multi_eligible"] / acc["n_broad"]
    frac_none = 100 * acc["n_no_gate_positive"] / acc["n_broad"]
    frac_other = 100 * acc["n_other"] / acc["n_broad"]
    xs = np.arange(len(acc))
    ax.bar(xs, frac_multi, width=0.72, color="#EE7733", zorder=3,
           label="≥2 compartments eligible (order decides)")
    ax.plot(xs, frac_none, "o-", color="#333333", ms=4, lw=1.2, zorder=4,
            label="no gate marker positive")
    ax.plot(xs, frac_other, "s--", color="#BBBBBB", ms=4, lw=1.2, zorder=4, label="assigned Other")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r.donor} {STATUS_LABEL.get(r.status, '?')}" for r in acc.itertuples()],
                       rotation=90, fontsize=8)
    ax.set_ylabel("% of the donor's cells")
    ax.set_title("B — per donor. 'No gate marker positive' and 'Other' should coincide exactly:\n"
                 "with zero Unresolved, Other IS the all-negative population.", fontsize=10.5,
                 loc="left", pad=26)
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(0.0, 0.98))

    fig.suptitle("Gate 6 — lineage accounting and where the ordered residual actually decides",
                 fontsize=13, y=1.02)
    fig.savefig(out / "gate6_audit.png")
    plt.close(fig)
    print(f"[saved] {out / 'gate6_audit.png'}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5")
    ap.add_argument("--figure-only", action="store_true",
                    help="redraw from the cached CSVs instead of re-streaming 25M cells")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    acc_csv, ties_csv = out / "gate6_accounting.csv", out / "gate6_tie_breaks.csv"

    if a.figure_only:
        if not (acc_csv.is_file() and ties_csv.is_file()):
            raise SystemExit(f"[err] --figure-only needs {acc_csv} and {ties_csv}")
        acc = pd.read_csv(acc_csv, dtype={"donor": str})
        ties = pd.read_csv(ties_csv, dtype={"donor": str})
    else:
        summaries, ties_frames = [], []
        for donor in cfg.discover_donors(cfg.broad_dir):
            s, t = audit_donor(cfg, donor)
            summaries.append(s)
            ties_frames.append(t)
        acc = pd.DataFrame(summaries)
        ties = pd.concat(ties_frames, ignore_index=True)

    # Disease status for the donor labels. `broad_lineage_composition.csv` only exists once lineage has
    # been re-run since it started writing one, so fall back to the Gate-5 composition table.
    acc["status"] = "?"
    for source in (cfg.phenotype_dir / "broad_lineage_composition.csv",
                   out / "composition_by_region.csv"):
        if source.is_file():
            table = pd.read_csv(source, dtype={"donor": str})
            status = table.drop_duplicates("donor").set_index("donor")["status"]
            mapped = acc["donor"].map(status)
            if mapped.notna().any():
                acc["status"] = mapped.fillna("?")
                break

    acc.to_csv(acc_csv, index=False)
    ties.to_csv(ties_csv, index=False)
    figure(acc, ties, out)

    # --- the assertions -------------------------------------------------------------------------
    failures = []
    if int(acc["n_duplicate_object_id"].sum()):
        failures.append(f"{int(acc['n_duplicate_object_id'].sum()):,} duplicate object_id in broad output")
    if int(acc["n_canonical_unassigned"].sum()):
        failures.append(f"{int(acc['n_canonical_unassigned'].sum()):,} canonical cells never assigned")
    if int(acc["n_broad_not_canonical"].sum()):
        failures.append(f"{int(acc['n_broad_not_canonical'].sum()):,} assigned cells are not canonical")
    mismatch = acc[acc["n_no_gate_positive"] != acc["n_other"] + acc["n_unresolved"]]
    if len(mismatch):
        failures.append(f"{len(mismatch)} donors where 'no gate positive' != Other + Unresolved "
                        f"({mismatch['donor'].tolist()})")
    reason_total = acc[["n_reason_positive", "n_reason_all_negative", "n_reason_blocked"]].sum(axis=1)
    if not (reason_total == acc["n_broad"]).all():
        failures.append("assignment_reason does not partition every cell")

    total = int(acc["n_broad"].sum())
    print(f"\n[total] {total:,} cells assigned | "
          f"{int(acc['n_multi_eligible'].sum()):,} "
          f"({100 * acc['n_multi_eligible'].sum() / total:.1f}%) had ≥2 compartments eligible, "
          "resolved silently by COMPARTMENT_ORDER")
    if failures:
        print("\n[FAIL] accounting is not exact:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"[PASS] every one of the {total:,} canonical cells is assigned exactly once; "
          "Other is exactly the all-gates-negative population.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
