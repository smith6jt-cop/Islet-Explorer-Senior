#!/usr/bin/env python3
"""
Build the manual gold-standard annotation bundle (data/gold_standard_v1/).

WHY. There is no accuracy measurement anywhere in this project. Every QC in the pipeline is an
internal-consistency check — pair states, seed stability, divisor reproducibility, the compartment
accounting — and internal consistency is exactly what a systematic error survives. A ~55 pp epithelial
error passed 217 tests, 11 review bundles and three gates. Nothing in the repo can currently answer
"what fraction of the cells called Immune are immune cells", which is the question the study's
conclusion rests on.

WHAT THIS PRODUCES. A blinded, stratified sample of canonical cells for two independent human
annotators, plus the QuPath round-trip to annotate it in and the scorer to read it back. It performs
no annotation and asserts no accuracy — it is the instrument, not the measurement.

Sampling is stratified because the population is not: the immune compartments carry the conclusion and
are a few percent of cells, so a uniform sample of 2,500 would contain too few to estimate immune
precision. Every stratum records an inclusion weight (N_stratum / n_sampled), so the scorer reports
population-weighted rates and the stratification does not bias the answer.

    python scripts/gold_standard/build_bundle.py
    python scripts/gold_standard/build_bundle.py --donor 6539 --n-per-donor 200
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.config import COMPARTMENT_ORDER, OTHER_LABEL, UNRESOLVED_LABEL   # noqa: E402

# Five donors spanning all three disease groups AND all three acquisition-software batches AND both
# batch-2 panel variants, so a per-compartment accuracy estimate is never confounded with the
# ND/Fusion-2.4.1 aliasing (Fisher p = 0.0147; see data/donor_acquisition_batch.csv). All five are in
# the Gate-3 pilot set, so every arm's calls already exist for them.
DEFAULT_DONORS = ["6476", "6539", "6450", "6436", "6591"]
N_PER_DONOR = 500
# Fractions per stratum. Immune is over-sampled on purpose (see module docstring); the uniform stratum
# is what makes RECALL estimable — without it the sample only contains cells the pipeline already
# called, and a marker the pipeline misses entirely would score perfectly.
STRATA = {"immune_called": 0.40, "islet_proximal": 0.25, "uniform": 0.35}
ISLET_PROXIMAL_UM = 20.0

VOCAB = COMPARTMENT_ORDER + [OTHER_LABEL, "Unsure"]


# Arms the immune stratum draws from. Deliberately NOT every arm on disk: reference_only/max is not a
# viable policy (restore_apply halts on it) and reference_only/p95 over-calls every sparse marker
# 5-60x, so including either would flood the stratum with candidates no candidate policy produces.
IMMUNE_ARMS = ["ctrl-reference_and_target_div-max", "ctrl-reference_only_div-p999",
               "ctrl-reference_only_div-p99"]


def _load_arm_calls(arm_root: Path, donor: str, arms=IMMUNE_ARMS) -> pd.DataFrame:
    """Per-cell compartment for each arm/design, rebuilt from the gated parquets on disk.

    run_arms.py keeps only aggregates, so the immune stratum is defined from the gated tri-state
    calls directly: a cell is "immune-called" if ANY candidate arm called any of the three immune
    gate markers positive. That is the widest net over policies still in contention, which is the
    right frame for measuring immune PRECISION (every candidate a viable policy could produce is
    eligible for review).
    """
    frames = []
    for arm in arms:
        gated_files = sorted(glob.glob(str(arm_root / arm / "restore_gated_redsea"
                                           / f"donor_id={donor}" / "*.parquet")))
        if not gated_files:
            continue
        gated = gated_files[0]
        g = pd.read_parquet(gated, columns=["object_id", "CD3e_state", "CD68_state", "CD11b_state"])
        g["object_id"] = g["object_id"].astype(str)
        g[f"immune_{arm}"] = (g[["CD3e_state", "CD68_state", "CD11b_state"]] == "positive").any(axis=1)
        frames.append(g[["object_id", f"immune_{arm}"]].set_index("object_id"))
    if not frames:
        raise SystemExit(f"[err] no gated parquets for donor {donor} under {arm_root}; run run_arms.py first")
    wide = pd.concat(frames, axis=1)
    wide["immune_any_arm"] = wide.any(axis=1)
    return wide.reset_index()


def sample_donor(cells_dir: Path, arm_root: Path, donor: str, n: int, rng) -> pd.DataFrame:
    cells = pd.read_parquet(
        sorted((cells_dir / f"donor_id={donor}").glob("*.parquet"))[0],
        columns=["object_id", "image", "X_centroid", "Y_centroid", "cell_region", "islet_num",
                 "dist_islet", "cell_area"],
    )
    cells["object_id"] = cells["object_id"].astype(str)
    calls = _load_arm_calls(arm_root, donor)
    d = cells.merge(calls[["object_id", "immune_any_arm"]], on="object_id", how="left", validate="1:1")
    d["immune_any_arm"] = d["immune_any_arm"].fillna(False)

    # Strata are disjoint and exhaustive, assigned in priority order.
    stratum = np.full(len(d), "uniform", dtype=object)
    stratum[(d["dist_islet"].to_numpy() <= ISLET_PROXIMAL_UM)] = "islet_proximal"
    stratum[d["immune_any_arm"].to_numpy()] = "immune_called"
    d["stratum"] = stratum

    picks = []
    for name, frac in STRATA.items():
        pool = d.index[d["stratum"] == name].to_numpy()
        want = int(round(n * frac))
        take = min(want, len(pool))
        if take < want:
            print(f"[warn] {donor}/{name}: wanted {want}, only {len(pool)} available")
        if take == 0:
            continue
        chosen = rng.choice(pool, size=take, replace=False)
        sub = d.loc[chosen].copy()
        sub["stratum_N"] = len(pool)               # population size of this stratum
        sub["stratum_n"] = take                    # how many were drawn
        sub["weight"] = len(pool) / take           # each drawn cell stands for this many cells
        picks.append(sub)
    s = pd.concat(picks, ignore_index=True)
    s["donor"] = str(donor)
    s["donor_N"] = len(d)
    return s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor", action="append", default=None)
    ap.add_argument("--n-per-donor", type=int, default=N_PER_DONOR)
    ap.add_argument("--cells", type=Path, default=REPO / "data" / "cells")
    ap.add_argument("--arm-root", type=Path, default=REPO / "data" / "_arm")
    ap.add_argument("--out", type=Path, default=REPO / "data" / "gold_standard_v1")
    ap.add_argument("--annotators", nargs="*", default=["A", "B"])
    ap.add_argument("--seed", type=int, default=20260727)
    a = ap.parse_args(argv)

    donors = a.donor or DEFAULT_DONORS
    rng = np.random.default_rng(a.seed)
    a.out.mkdir(parents=True, exist_ok=True)

    frames = [sample_donor(a.cells, a.arm_root, d, a.n_per_donor, rng) for d in donors]
    s = pd.concat(frames, ignore_index=True)

    # Blind: shuffle so neither donor nor stratum can be read off the row order, then assign ids.
    s = s.sample(frac=1.0, random_state=a.seed).reset_index(drop=True)
    s.insert(0, "review_id", [f"GS{i:05d}" for i in range(len(s))])

    # 1. The key — carries stratum, weight and every pipeline call. NEVER given to an annotator.
    key_cols = ["review_id", "donor", "object_id", "image", "X_centroid", "Y_centroid",
                "cell_region", "islet_num", "dist_islet", "cell_area",
                "stratum", "stratum_N", "stratum_n", "weight", "immune_any_arm", "donor_N"]
    s[key_cols].to_csv(a.out / "sample_key.csv", index=False)

    # 2. Blinded annotation sheets — location only. No pipeline call, no stratum, no marker value.
    sheet = s[["review_id", "donor", "image", "X_centroid", "Y_centroid"]].copy()
    sheet["compartment"] = ""          # <- the annotator fills this from VOCAB
    sheet["confidence"] = ""           # high | medium | low
    sheet["notes"] = ""
    for who in a.annotators:
        sheet.to_csv(a.out / f"annotation_sheet_{who}.csv", index=False)

    # 3. QuPath round-trip: per-donor object_id -> review_id, for the importer groovy.
    qp = a.out / "qupath"
    qp.mkdir(exist_ok=True)
    for donor, sub in s.groupby("donor"):
        sub[["object_id", "review_id", "image"]].to_csv(qp / f"gold_standard_{donor}.csv", index=False)

    manifest = {
        "created": "2026-07-27",
        "seed": a.seed,
        "donors": donors,
        "n_per_donor": a.n_per_donor,
        "n_total": int(len(s)),
        "strata": STRATA,
        "islet_proximal_um": ISLET_PROXIMAL_UM,
        "vocabulary": VOCAB,
        "annotators": a.annotators,
        "immune_stratum_definition": (
            "any of CD3e/CD68/CD11b called positive under ANY completed divisor arm — the widest net, "
            "so immune precision is measured over every candidate the pipeline could produce"
        ),
        "weighting": "weight = stratum_N / stratum_n; sum(weight) reconstructs each donor's cell count",
        "blinding": "rows shuffled; sheets carry location only — no pipeline call, stratum or marker value",
    }
    (a.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(s.groupby(["donor", "stratum"]).size().unstack(fill_value=0).to_string())
    print(f"\nsampled {len(s):,} cells; sum(weight) = {s.weight.sum():,.0f} against "
          f"{s.groupby('donor').donor_N.first().sum():,} cells in the five donors")
    print(f"[done] {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
