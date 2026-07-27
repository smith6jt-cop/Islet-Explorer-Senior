#!/usr/bin/env python3
"""
Anatomy check — assign compartments under every divisor arm and under two Epithelial gate designs.

WHY. The Gate-4 pilot was rejected for ~15-19% Epithelial and ~65% ``Other``. The RESTORE control
definition has since been fixed, moving E-cadherin to 43-54% under ``reference_only``/``p99``, but the
compartments' own accounting implies ~75%. That accounting is monotone in the divisor -- every policy
crosses 100% at exactly one point -- so it locates an operating point without excluding compensating
errors, and cannot arbitrate on its own. The independent arbiter is anatomy: the QuPath ``Islet_N``
annotations are drawn from the image, not from any marker threshold, and they are validated as anatomy
by hormone expression (islet-core INS 407 / GCG 1004 against 10 / 71 in tissue, cohort medians).

Two hypotheses are tested in one pass, because they trade against each other:

  H1 (divisor)  the E-cadherin threshold is too high. Varied by ``divisor_statistic``.
  H2 (gate)     one E-cadherin threshold cannot cover both epithelial populations this tissue has.
                Islet-core (endocrine) cells sit BELOW acinar on E-cadherin (630 vs 762) and ABOVE it
                on EpCAM (656 vs 516), so the two markers have opposite polarity across the boundary
                that matters. Varied by gating Epithelial on E-cadherin alone vs a 4-marker union.

Nothing here is production. Every arm writes to its own scratch data root under ``data/_arm/`` whose
``cells``/``cells_redsea``/``restore_pair_validation`` are relative symlinks to the canonical
directories, so canonical ``data/`` is read-only throughout. The comparator pairs
(``PAIR_SETS['epithelial-comparator']``) are diagnostic: they are NOT in ``ACCEPTED_PAIRS`` and emit no
production call.

    python scripts/anatomy_check/run_arms.py                       # all arms, the 7 pilot donors
    python scripts/anatomy_check/run_arms.py --arm reference_only/p99 --donor 6539
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))

from phenocycler import lineage, load_config                                    # noqa: E402
from phenocycler import restore_apply as ra                                     # noqa: E402
from phenocycler.config import COMPARTMENT_GATES                                # noqa: E402
from phenocycler.restore_validation import (                                    # noqa: E402
    PAIR_SETS,
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    evaluate_locked_pair,
    load_validation_sample,
    locked_pair_metrics,
)

# The 5 policy arms. reference_and_target/max reproduces v9 (control capped at the target's own Otsu
# floor, which made the divisor BE that floor); the reference_only rows walk the divisor statistic.
ARMS = [
    "reference_and_target/max",
    "reference_only/max",
    "reference_only/p999",
    "reference_only/p99",
    "reference_only/p95",
]
# The Gate-3 pilot donors: 6539 and 6476 by standing decision, plus the weakest case per pair.
PILOT_DONORS = ["115", "6380", "6436", "6450", "6476", "6539", "6591"]

# The Epithelial gate designs under test. "single" is production today; "union" adds the three other
# epithelial markers the panel carries, mirroring how Immune already needs four gates for want of CD45.
# Each comparator target takes ONE reference here; both references are screened separately so the
# choice can be revisited from the comparator manifest.
UNION_EPITHELIAL = ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM"]
COMPARATOR_REFERENCE = "CD31"          # E-cadherin's own frozen reference -> like-for-like comparison


def arm_slug(arm: str) -> str:
    ctrl, stat = arm.split("/")
    return f"ctrl-{ctrl}_div-{stat}"


def arm_root(out_root: Path, arm: str) -> Path:
    """A scratch data root whose inputs are symlinks to canonical data/ (outputs land inside it)."""
    root = out_root / arm_slug(arm)
    root.mkdir(parents=True, exist_ok=True)
    for name in ("cells", "cells_redsea", "restore_pair_validation"):
        link = root / name
        if not link.is_symlink() and not link.exists():
            # relative, so the tree stays valid if data/ is moved (it has been, twice)
            link.symlink_to(os.path.relpath(REPO / "data" / name, root))
    return root


def apply_comparator(cfg, donors, pvc: PairValidationConfig) -> pd.DataFrame:
    """Per-cell tri-state calls for the DIAGNOSTIC epithelial-comparator pairs.

    Deliberately separate from ``restore_apply.run_apply``: that function's contract is the frozen
    ACCEPTED_PAIRS and it halts on any pair it does not recognise, which is exactly the guard that
    should not be loosened for a hypothesis test. The frozen math (``evaluate_locked_pair``) and the
    tested per-cell mapper (``restore_apply.cell_calls_from_evaluation``) are reused unchanged, so the
    comparator calls are produced by the same code path as production ones.

    Unlike production, a non-callable pair state here is RECORDED, not fatal: a comparator marker
    failing on one donor is a result about that marker, not a pipeline error.
    """
    out_dir = cfg.data_dir / "restore_gated_comparator"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for donor in donors:
        donor_df, audit = load_validation_sample(
            cfg.restore_allcells_sample_dir, cfg.cells_dir, cfg.cells_redsea_dir, donor
        )
        cols = {"object_id": donor_df["object_id"].astype(str).to_numpy(),
                "donor_id": np.full(len(donor_df), str(donor))}
        for target, reference in PAIR_SETS["epithelial-comparator"]:
            ev = evaluate_locked_pair(
                donor_df, donor, target, reference,
                config=pvc,
                expression_review=ExpressionReview.UNREVIEWED,
                pair_review=PairReview.ACCEPTED,       # diagnostic: assert nothing about sign-off
            )
            pos, norm, state = ra.cell_calls_from_evaluation(ev)
            key = f"{target}__{reference}"
            cols[f"{key}_pos"], cols[f"{key}_norm"], cols[f"{key}_state"] = pos, norm, state
            row = dict(locked_pair_metrics(ev))
            row.update(donor=str(donor), image=audit.get("image"),
                       called_positive_n=int((state == ra.CELL_POSITIVE).sum()))
            rows.append(row)
        part = out_dir / f"donor_id={donor}"
        part.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(cols).to_parquet(part / "data_0.parquet", index=False)
        print(f"[comparator {donor}] {len(donor_df):,} cells | " + " ".join(
            f"{t}<-{r}={100 * float((cols[f'{t}__{r}_state'] == ra.CELL_POSITIVE).mean()):.1f}%"
            for t, r in PAIR_SETS["epithelial-comparator"]), flush=True)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "comparator_manifest.csv", index=False)
    return manifest


def _merge_union_gated(cfg, donor: str) -> Path:
    """Gated parquet carrying the accepted targets PLUS the union Epithelial markers, for assign_donor."""
    acc = pd.read_parquet(cfg.restore_gated_dir / f"donor_id={donor}" / "data_0.parquet")
    cmp_ = pd.read_parquet(cfg.data_dir / "restore_gated_comparator" / f"donor_id={donor}" / "data_0.parquet")
    acc["object_id"] = acc["object_id"].astype(str)
    cmp_["object_id"] = cmp_["object_id"].astype(str)
    extra = [m for m in UNION_EPITHELIAL if m != "E_cadherin"]
    take = {f"{m}__{COMPARATOR_REFERENCE}_{suf}": f"{m}_{suf}"
            for m in extra for suf in ("pos", "norm", "state")}
    merged = acc.merge(
        cmp_[["object_id"] + list(take)].rename(columns=take),
        on="object_id", how="left", validate="1:1",
    )
    if merged[[f"{m}_state" for m in extra]].isna().any().any():
        raise ValueError(f"donor {donor}: comparator calls missing for some cells after the join")
    dst = cfg.data_dir / "restore_gated_union" / f"donor_id={donor}" / "data_0.parquet"
    dst.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(dst, index=False)
    return dst


def assign_both_gates(cfg, donors, arm: str) -> pd.DataFrame:
    """Run the ordered residual under the single-marker and union Epithelial gates. Returns per-cell rows."""
    frames = []
    original = dict(lineage.COMPARTMENT_GATES)      # restored after every assign_donor call
    for donor in donors:
        gated = str(cfg.restore_gated_dir / f"donor_id={donor}" / "data_0.parquet")
        union_gated = str(_merge_union_gated(cfg, donor))
        # dist_islet is not carried onto the broad parquet (assign_donor attaches only cell_region /
        # islet_num), and the radial profile needs it.
        cells = pd.read_parquet(
            sorted((cfg.cells_dir / f"donor_id={donor}").glob("*.parquet"))[0],
            columns=["object_id", "dist_islet", "X_centroid", "Y_centroid"],
        )
        cells["object_id"] = cells["object_id"].astype(str)
        for design, path, gates in (
            ("single", gated, original),
            ("union", union_gated, {**original, "Epithelial": UNION_EPITHELIAL}),
        ):
            lineage.COMPARTMENT_GATES = gates       # assign_donor reads this module global
            try:
                out = lineage.assign_donor(donor, path, cfg.cells_dir)
            finally:
                lineage.COMPARTMENT_GATES = original
            out = out.merge(cells, on="object_id", how="left", validate="1:1")
            out["arm"], out["gate_design"], out["donor"] = arm, design, str(donor)
            frames.append(out)
            comp = out["compartment"].value_counts(normalize=True).mul(100).round(1)
            print(f"[{arm} {design} {donor}] {len(out):,} cells | "
                  + " ".join(f"{k}={v}%" for k, v in comp.items()), flush=True)
    return pd.concat(frames, ignore_index=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", action="append", default=None, choices=ARMS)
    ap.add_argument("--donor", action="append", default=None)
    ap.add_argument("--out-root", type=Path, default=REPO / "data" / "_arm")
    ap.add_argument("--summary", type=Path, default=REPO / "data" / "anatomy_check")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--force", action="store_true", help="redo arms whose aggregates already exist")
    a = ap.parse_args(argv)

    arms = a.arm or ARMS
    donors = a.donor or PILOT_DONORS
    a.summary.mkdir(parents=True, exist_ok=True)

    outcomes = []
    for arm in arms:
        ctrl, stat = arm.split("/")
        done = a.summary / f"region_composition_{arm_slug(arm)}.csv"
        if done.exists() and not a.force:
            print(f"[skip] {arm} — {done.name} exists (pass --force to redo)", flush=True)
            outcomes.append(dict(arm=arm, outcome="already_done", detail=""))
            continue
        root = arm_root(a.out_root, arm)
        cfg = load_config(
            REPO / "config.ini",
            data_dir=root,
            n_jobs=a.jobs,
            restore_control_definition=ctrl,
            restore_divisor_statistic=stat,
        )
        pvc = ra.config_from_pipeline(cfg)
        print(f"\n{'=' * 78}\n=== ARM {arm}  ->  {root}\n{'=' * 78}", flush=True)

        try:
            ra.run_apply(cfg, donors=donors)
        except Exception as e:
            # restore_apply halts a production run on any technical pair state, by design. An arm that
            # cannot emit a call for every donor is not a candidate policy — that is a RESULT about the
            # arm, so record it and keep going rather than losing the remaining arms to it.
            print(f"\n[ARM {arm} NOT VIABLE] {type(e).__name__}: {e}\n", flush=True)
            outcomes.append(dict(arm=arm, outcome="halted", detail=str(e)))
            continue
        apply_comparator(cfg, donors, pvc)
        cells = assign_both_gates(cfg, donors, arm)
        outcomes.append(dict(arm=arm, outcome="ok", detail=""))

        # Per-cell tables are large and per-arm; only the aggregates the figures need are kept.
        agg = (
            cells.groupby(["arm", "gate_design", "donor", "cell_region", "compartment"], observed=True)
            .size().rename("n").reset_index()
        )
        agg.to_csv(a.summary / f"region_composition_{arm_slug(arm)}.csv", index=False)

        edges = [0, 1, 5, 10, 20, 40, 80, 160, 320, np.inf]
        radial = cells.assign(
            dist_bin=pd.cut(cells["dist_islet"], bins=edges, right=False)
        ).groupby(["arm", "gate_design", "donor", "dist_bin", "compartment"], observed=True) \
            .size().rename("n").reset_index()
        radial["dist_bin"] = radial["dist_bin"].astype(str)
        radial.to_csv(a.summary / f"radial_profile_{arm_slug(arm)}.csv", index=False)
        print(f"[arm {arm}] aggregates -> {a.summary}", flush=True)

    oc = pd.DataFrame(outcomes)
    oc.to_csv(a.summary / "arm_outcomes.csv", index=False)
    print("\n[arms]\n" + oc[["arm", "outcome"]].to_string(index=False))
    print("\n[done] anatomy-check aggregates in", a.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
