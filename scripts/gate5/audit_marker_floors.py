#!/usr/bin/env python
"""Audit the RESTORE raw input floor for every viable marker, per donor, across candidate methods.

WHY THIS EXISTS
---------------
The raw input floor is not a diagnostic. It is load-bearing:

    raw input floor -> exclusive training arms (balanced_arm_fit_indices)
                    -> NNMF separator -> Step-1 anchor -> reference control set
                    -> Step-2 divisor -> the production positivity call

`estimate_marker_input_floor` applies ONE method per marker from the hardcoded
`MARKER_INPUT_FLOOR_METHODS` table. Only five markers (CD3e, CD20, CD31, E_cadherin, Vimentin) were
ever audited before that table was written (`data/restore_pair_validation/marker_floor_audit_v1`,
2026-07-21); the other nineteen gate markers got their method by analogy, and the script that produced
that audit no longer exists. Measured 2026-07-28 from `screen_v12_production`, four `linear_otsu`
markers place their floor in background on some donors -- SST calls 90.3% of cells positive on one donor
and 3.6% on another, INS 85.8% vs 2.3% with a 333x floor spread -- which collapses the arm separation and
therefore the divisor on exactly those donors.

This script rebuilds the lost instrument and extends it to every viable marker in the panel.

WHAT IT MEASURES
----------------
The floor is computed on the RAW (pre-REDSEA) QuPath whole-cell mean, over every canonical cell with a
finite non-negative value -- verified pair-invariant: across the 500 (donor, marker) combinations in
`screen_v12_production`, a marker's floor is bit-identical in every pair it appears in. So one floor per
(donor, marker), and this script reproduces it exactly by calling the production functions rather than
reimplementing them.

Per (donor, marker, method) it reports:

  threshold           the floor that method would set
  positive_fraction   share of cells above it
  class_fold          mean(above) / mean(below) -- a floor sitting in background gives ~1
  separation          Otsu's criterion at that threshold, between-class / total variance, in [0, 1];
                      frequency-free, and it answers "is there a real split here at all"
  membrane_margin     median(Membrane - Nucleus) among called-positive cells, for the twelve markers
                      that have compartment measurements. A floor in background has no membrane
                      preference. This is the diagnostic that does not depend on how abundant the
                      marker is.

NOT a frequency test. Nothing here asks whether a marker "should" be rare or abundant. The evidence a
floor has failed is the marker's OWN donor-to-donor instability -- a marker cannot be 3.6% positive on
one donor and 90.3% on another -- read against its own distribution, per the repository working rule.

USAGE
    ~/miniconda3/envs/scvi-env/bin/python scripts/gate5/audit_marker_floors.py
    ... --donors 6566 6380 --markers INS SST --no-figures
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# Each donor is already one process; letting BLAS/OpenMP fan out inside them oversubscribes badly
# (measured: 84 s wall for 45 min of CPU on a single donor). Must precede the numpy import.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np                                                             # noqa: E402
import pandas as pd                                                            # noqa: E402

# A marker with no compartment measurements has an all-NaN membrane margin; nanmedian says so loudly
# and repeatedly. The NaN is the intended answer, not a problem.
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encountered")

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))

from phenocycler import marker_taxonomy as mt                                  # noqa: E402
from phenocycler.config import load_config                                     # noqa: E402
from phenocycler.cohort import ensure_eligible_donors                          # noqa: E402
from phenocycler.restore_validation import (                                   # noqa: E402
    MARKER_INPUT_FLOOR_METHODS,
    raw_otsu_floor,
    raw_triangle_floor,
)

# The eight candidates recorded in marker_floor_audit_v1/manifest.json. `linear_otsu` and
# `linear_triangle` call the PRODUCTION functions so this audit reproduces the shipped floors exactly
# rather than approximating them -- that identity is asserted by --validate.
METHODS = (
    "linear_otsu",
    "linear_triangle",
    "linear_yen",
    "log_otsu",
    "log_triangle",
    "log_yen",
    "log_li",
    "log_gmm_crossover",
    # The composite production actually ships for 12 of the 24 gate markers (CD3e, CD20, CD79a, CD11c,
    # Iba1, CD206, MPO, CD4, CD8, CD34, Podoplanin, Keratin_5). Without it the audit cannot evaluate
    # half the gate set against the floor it really uses.
    "max_linear_otsu_triangle",
    # `raw_otsu_floor` clips its input at the 99.5th percentile before binning, to stop a single
    # extreme value from collapsing a 512-bin linear histogram into one bin. For a marker whose
    # positive population is SMALLER than 0.5% of cells, that clip deletes the entire positive
    # population before Otsu sees it, and the threshold necessarily lands inside the background.
    # Measured on donor 6566 INS: clipped floor 17.4 -> 23.2% of the tissue "positive"; unclipped
    # 864.2 -> 0.02%, about 356 cells, against the maintainer's independent image estimate of "a
    # couple hundred at most". Harmless where the positive population is large (115 EpCAM 616.5 vs
    # 604.6; 6461 E-cadherin 924.0 vs 926.2). These two arms measure the fix.
    "linear_otsu_p9999",
    "linear_otsu_unclipped",
    # No fixed clip quantile is right for every marker -- p99.99 is better for CD163 and CD11c, fully
    # unclipped is far better for INS (0.416 vs 0.057), because the correct clip depends on how rare
    # the positive population is, which is the thing being measured. The clip only exists to keep a
    # 512-bin LINEAR histogram from collapsing under a long tail, so the principled alternative is to
    # bin in log space and not discard anything. NOTE `log_otsu` above is contaminated by the same
    # clip (it applies raw_otsu_floor, defaults and all, to the logged values); this arm is the clean
    # comparison.
    "log_otsu_unclipped",
)

GMM_SUBSAMPLE = 200_000       # GMM is the only method too slow at 1.2M cells; deterministic subsample
GMM_SEED = 20260728


# --------------------------------------------------------------------------------------- thresholds

def _log_threshold(values: np.ndarray, fn) -> float:
    """Apply a skimage threshold in log1p space and map the answer back to raw intensity."""
    logged = np.log1p(values)
    if logged.min() == logged.max():
        raise ValueError("constant intensities")
    return float(np.expm1(fn(logged)))


def _gmm_crossover(values: np.ndarray) -> float:
    """Where a 2-component log-space Gaussian mixture's weighted densities cross.

    Uses the crossover rather than a component mean because the crossover is the decision boundary --
    the point at which a cell stops being better explained by the low component. Falls back to the
    midpoint between the means when the densities do not cross between them (which happens when one
    component is swallowed by the other, itself a signal that there is no second mode).
    """
    from sklearn.mixture import GaussianMixture

    logged = np.log1p(values)
    rng = np.random.default_rng(GMM_SEED)
    fit = logged if logged.size <= GMM_SUBSAMPLE else rng.choice(logged, GMM_SUBSAMPLE, replace=False)
    gm = GaussianMixture(n_components=2, random_state=GMM_SEED, n_init=2).fit(fit.reshape(-1, 1))
    order = np.argsort(gm.means_.ravel())
    lo, hi = gm.means_.ravel()[order]
    if not np.isfinite([lo, hi]).all() or hi <= lo:
        raise ValueError("GMM did not resolve two ordered components")
    grid = np.linspace(lo, hi, 2048).reshape(-1, 1)
    # log-density of each component, weighted; the boundary is where their difference changes sign
    from scipy.stats import norm
    w = gm.weights_[order]
    m = gm.means_.ravel()[order]
    s = np.sqrt(gm.covariances_.ravel()[order])
    diff = w[1] * norm.pdf(grid.ravel(), m[1], s[1]) - w[0] * norm.pdf(grid.ravel(), m[0], s[0])
    sign_change = np.flatnonzero(np.diff(np.signbit(diff)))
    x = grid.ravel()[sign_change[0]] if sign_change.size else 0.5 * (lo + hi)
    return float(np.expm1(x))


def candidate_thresholds(values: np.ndarray) -> dict[str, float]:
    """Every candidate floor for one donor-marker. A method that cannot fit returns NaN, never a guess."""
    from skimage.filters import threshold_li, threshold_triangle, threshold_yen

    out: dict[str, float] = {}
    for name in METHODS:
        try:
            if name == "linear_otsu":
                out[name] = raw_otsu_floor(values)          # production function, do not reimplement
            elif name == "linear_triangle":
                out[name] = raw_triangle_floor(values)      # production function, do not reimplement
            elif name == "linear_yen":
                out[name] = float(threshold_yen(values, nbins=256))
            elif name == "log_otsu":
                out[name] = _log_threshold(values, lambda v: raw_otsu_floor(v))
            elif name == "log_triangle":
                out[name] = _log_threshold(values, lambda v: threshold_triangle(v, nbins=256))
            elif name == "log_yen":
                out[name] = _log_threshold(values, lambda v: threshold_yen(v, nbins=256))
            elif name == "log_li":
                out[name] = _log_threshold(values, threshold_li)
            elif name == "log_gmm_crossover":
                out[name] = _gmm_crossover(values)
            elif name == "max_linear_otsu_triangle":
                out[name] = max(out["linear_otsu"], out["linear_triangle"])
            elif name == "linear_otsu_p9999":
                out[name] = raw_otsu_floor(values, upper_quantile=0.9999)
            elif name == "linear_otsu_unclipped":
                out[name] = raw_otsu_floor(values, upper_quantile=1.0)
            elif name == "log_otsu_unclipped":
                out[name] = _log_threshold(values, lambda v: raw_otsu_floor(v, upper_quantile=1.0))
        except Exception:                                   # a method that fails is reported as NaN
            out[name] = float("nan")
    return out


# -------------------------------------------------------------------------------------- diagnostics

def threshold_diagnostics(values: np.ndarray, threshold: float) -> dict[str, float]:
    """Frequency-free descriptions of what a threshold does to this donor's own distribution."""
    if not np.isfinite(threshold):
        return {"positive_n": np.nan, "positive_fraction": np.nan,
                "class_fold": np.nan, "separation": np.nan}
    above = values > threshold
    n_above, n = int(above.sum()), values.size
    if n_above == 0 or n_above == n:
        # A degenerate split has no fold and no separation; say so rather than emit 0 or inf.
        return {"positive_n": n_above, "positive_fraction": n_above / n,
                "class_fold": np.nan, "separation": 0.0}
    mu_hi, mu_lo, mu = values[above].mean(), values[~above].mean(), values.mean()
    var = values.var()
    w = n_above / n
    between = w * (mu_hi - mu) ** 2 + (1 - w) * (mu_lo - mu) ** 2
    return {
        "positive_n": n_above,
        "positive_fraction": n_above / n,
        "class_fold": float(mu_hi / mu_lo) if mu_lo > 0 else np.nan,
        "separation": float(between / var) if var > 0 else 0.0,
    }


def membrane_margin(compartments: pd.DataFrame, marker: str, positive: np.ndarray) -> float:
    """median(Membrane - Nucleus) among called-positive cells; NaN when the marker has no compartments.

    A real surface/structural signal concentrates at the membrane band, so its positives carry a
    positive margin. A floor sitting in background selects cells at random with respect to compartment,
    so the margin collapses toward zero. This does not care how many cells are called.
    """
    mem, nuc = f"Membrane__{marker}", f"Nucleus__{marker}"
    if compartments is None or mem not in compartments or nuc not in compartments:
        return np.nan
    d = (compartments[mem].to_numpy(float) - compartments[nuc].to_numpy(float))[positive]
    d = d[np.isfinite(d)]
    return float(np.median(d)) if d.size else np.nan


# ------------------------------------------------------------------------------------ per-donor run

def viable_markers(columns) -> list[str]:
    """Every marker in the panel except the taxonomy's EXCLUDED (DAPI, IAPP) -- maintainer, 2026-07-28."""
    known = set(mt.TYPE) | set(mt.PROCESS) | set(mt.MASK) | set(mt.PROLIFERATION)
    return [c for c in columns if c in known and c not in set(mt.EXCLUDED)]


HIST_BINS = 120
HIST_RANGE = (-1.0, 5.0)          # log10 raw intensity; QuPath means span roughly 10^-1..10^4


def _histogram(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fixed-grid log10 histogram so donors are directly comparable on one axis.

    Shared edges across donors are the point: the ridgeline exists to show that one donor's floor sits
    where another donor's background is, which is unreadable if every donor has its own binning.
    """
    v = values[values > 0]
    if v.size < 2:
        centers = 0.5 * (np.linspace(*HIST_RANGE, HIST_BINS + 1)[:-1]
                         + np.linspace(*HIST_RANGE, HIST_BINS + 1)[1:])
        return centers, np.zeros(HIST_BINS)
    counts, edges = np.histogram(np.log10(v), bins=HIST_BINS, range=HIST_RANGE)
    return 0.5 * (edges[:-1] + edges[1:]), counts / max(counts.sum(), 1)


def _density_at(centers: np.ndarray, density: np.ndarray, threshold: float) -> float:
    """Density at the floor as a fraction of this donor-marker's own peak density.

    THE test for "is this floor inside the background". A floor that has done its job sits where the
    background population has already decayed, so the density under it is a small fraction of the mode.
    A floor still on the flank of the background peak reads near 1.

    Entirely per-donor and per-marker, against that donor's own histogram: no frequency assumption, no
    cohort comparison, and no arbitrary positive-fraction cutoff -- which is why it catches 6566 INS
    (floor on the background flank at 23% positive) where a ">50% of cells" rule does not.
    """
    if not np.isfinite(threshold) or threshold <= 0 or density.max() <= 0:
        return float("nan")
    idx = int(np.argmin(np.abs(centers - np.log10(threshold))))
    return float(density[idx] / density.max())


def audit_donor(donor: str, cells_dir: Path, sample_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One donor: every viable marker x every candidate method.

    `cells_dir` selects the intensity the floor is estimated on:

      data/cells         RAW pre-REDSEA QuPath whole-cell mean -- what production uses today
      data/cells_redsea  the same cells after REDSEA lateral-spillover correction

    That production estimates its floors on RAW is worth stating plainly, because every MEMBERSHIP
    decision downstream inherits it: the target arm (`target_raw > target_floor & reference_raw <=
    reference_floor`), the negative control (`reference_raw > reference_floor`), and the NNMF
    separator features (`Cell__{marker}` is asserted identical to `raw__{marker}`). Only the Step-1
    separator, the Step-2 divisor and the final call read `target_redsea`. So a cell that is
    reference-positive only BECAUSE of spillover from a neighbour is admitted to the negative control,
    and its corrected value then helps set the divisor -- which is the same pan-positive-bright-cell
    failure REDSEA was introduced to remove. No rationale for the split is recorded anywhere in the
    code, the rebuild plan or the docs, so it reads as drift rather than a decision. This flag exists
    to measure it.
    """
    cells = pd.read_parquet(next((cells_dir / f"donor_id={donor}").glob("*.parquet")))
    markers = viable_markers(cells.columns)

    compartments = None
    sample_files = list((sample_root / f"donor_id={donor}").glob("*.parquet")) if sample_root else []
    if sample_files:
        comp = pd.read_parquet(sample_files[0])
        comp["object_id"] = comp["object_id"].astype(str)
        cells["object_id"] = cells["object_id"].astype(str)
        # left-join onto the canonical order so `positive` masks line up row-for-row
        compartments = cells[["object_id"]].merge(comp, on="object_id", how="left")

    rows, hist_rows = [], []
    for marker in markers:
        raw = cells[marker].to_numpy(float)
        finite = np.isfinite(raw) & (raw >= 0)          # exactly pair_input_qc's admission rule
        values = raw[finite]
        if values.size < 2:
            continue
        centers, density = _histogram(values)
        hist_rows.append(pd.DataFrame({"donor": donor, "marker": marker,
                                       "log_center": centers, "density": density}))
        for method, threshold in candidate_thresholds(values).items():
            diag = threshold_diagnostics(values, threshold)
            diag["density_at_floor_rel"] = _density_at(centers, density, threshold)
            positive_full = np.zeros(len(cells), bool)
            if np.isfinite(threshold):
                positive_full[np.flatnonzero(finite)[values > threshold]] = True
            rows.append({
                "donor": donor, "marker": marker, "method": method,
                "threshold": threshold, "available_n": int(values.size),
                **diag,
                "membrane_margin": membrane_margin(compartments, marker, positive_full),
                "production_method": MARKER_INPUT_FLOOR_METHODS.get(marker, ""),
                "is_gate_marker": marker in MARKER_INPUT_FLOOR_METHODS,
                "taxonomy": ("MASK" if marker in mt.MASK else
                             "PROLIFERATION" if marker in mt.PROLIFERATION else
                             "TYPE" if marker in mt.TYPE else "PROCESS"),
            })
    print(f"[floors] {donor}: {len(markers)} markers x {len(METHODS)} methods "
          f"on {len(cells):,} cells", flush=True)
    return pd.DataFrame(rows), pd.concat(hist_rows, ignore_index=True)


# ------------------------------------------------------------------------------------------ summary

def summarise(per_donor: pd.DataFrame) -> pd.DataFrame:
    """Per marker x method, across donors.

    THE HEADLINE STATISTIC IS `separation_consistency` = min(separation) / median(separation) across
    donors. It is the marker's own donor-to-donor stability, and it is what distinguishes the two
    things an absolute cutoff confuses:

      * a marker whose floor finds a real split on every donor -- EpCAM 0.92, E-cadherin 0.89,
        Vimentin 0.90 -- even where that split is weak or the marker is abundant;
      * a marker whose floor finds a real split on some donors and none on others -- INS 0.016,
        SST 0.16, CD163 0.26 -- which no biology explains and which is the failure this audit exists
        to find.

    Deliberately NOT a frequency test. Nothing here asks whether a marker should be rare or abundant;
    every quantity is that donor's own, and the comparison is a marker against itself.
    """
    def agg(g: pd.DataFrame) -> pd.Series:
        thr = g["threshold"].dropna()
        pf = g["positive_fraction"].dropna()
        sep = g["separation"].dropna()
        sep_med = sep.median() if sep.size else np.nan
        return pd.Series({
            "separation_consistency": (sep.min() / sep_med) if sep.size and sep_med > 0 else np.nan,
            # Donors where this method's floor finds a much weaker split than it finds on the same
            # marker's typical donor. The count, not a level.
            "n_donors_separation_collapsed": int((sep < 0.5 * sep_med).sum()) if sep_med > 0 else 0,
            "n_donors": int(len(g)),
            "n_donors_fitted": int(thr.size),
            "threshold_median": thr.median() if thr.size else np.nan,
            "floor_spread": (thr.max() / thr.min()) if thr.size and thr.min() > 0 else np.nan,
            "positive_fraction_median": pf.median() if pf.size else np.nan,
            "positive_fraction_min": pf.min() if pf.size else np.nan,
            "positive_fraction_max": pf.max() if pf.size else np.nan,
            "positive_fraction_spread": (pf.max() - pf.min()) if pf.size else np.nan,
            "n_donors_over_50pct": int((pf > 0.5).sum()),
            "class_fold_median": g["class_fold"].median(),
            "class_fold_min": g["class_fold"].min(),
            "separation_median": g["separation"].median(),
            "separation_min": g["separation"].min(),
            "density_at_floor_rel_median": g["density_at_floor_rel"].median(),
            "density_at_floor_rel_max": g["density_at_floor_rel"].max(),
            # The headline failure count: donors where this method's floor still sits on the flank of
            # the background peak (density under the floor > 10% of the donor's own modal density).
            "n_donors_floor_in_background": int((g["density_at_floor_rel"] > 0.10).sum()),
            "membrane_margin_median": g["membrane_margin"].median(),
            "production_method": g["production_method"].iloc[0],
            "is_gate_marker": bool(g["is_gate_marker"].iloc[0]),
            "taxonomy": g["taxonomy"].iloc[0],
        })
    return (per_donor.groupby(["marker", "method"], sort=False)
            .apply(agg, include_groups=False).reset_index())


# --------------------------------------------------------------------------------------- validation

def validate(per_donor: pd.DataFrame, screen: Path) -> int:
    """Two checks with known answers. If either fails, the instrument is wrong and nothing it says
    should be believed."""
    failures = 0

    # (1) Reproduce marker_floor_audit_v1's CD20 result, the case that drove the max(Otsu, Triangle)
    #     policy: linear_otsu ~49.4% positive, linear_triangle ~1.2%.
    #
    #     Those are MEDIANS OVER 20 DONORS, so the check is meaningless on a subset and is skipped
    #     rather than reported as a failure. It is also a soft check for a second reason: v1 ran
    #     2026-07-21, before `phenocycler.cell_qc` (2026-07-27) removed 6,303,013 objects that are not
    #     cells (31,454,041 -> 25,151,028, 20.04%). It therefore measured a different population, and an
    #     Otsu threshold is a property of the population. A moderate disagreement is expected and
    #     informative; an identical answer would be the surprise.
    cd20 = per_donor[per_donor.marker == "CD20"]
    n_donors = per_donor.donor.nunique()
    if cd20.empty or n_donors < 20:
        print(f"[validate] SKIP CD20 vs marker_floor_audit_v1 -- it is a 20-donor median and this run "
              f"has {n_donors}")
    else:
        got = {m: cd20[cd20.method == m]["positive_fraction"].median() for m in
               ("linear_otsu", "linear_triangle")}
        # Wide tolerance: the pre-QC/post-QC population differs by 20% of all objects.
        ok = (abs(got["linear_otsu"] - 0.494) < 0.15) and (abs(got["linear_triangle"] - 0.0118) < 0.05)
        print(f"[validate] CD20 vs marker_floor_audit_v1 (pre-cell-QC): linear_otsu "
              f"{got['linear_otsu']:.4f} (v1: 0.4940), linear_triangle {got['linear_triangle']:.4f} "
              f"(v1: 0.0118) -> {'PASS' if ok else 'FAIL'}")
        failures += 0 if ok else 1

    # (3) The maintainer's ground truth, 2026-07-28: donor 6566 HAS clear INS expression, and it is
    #     "a couple hundred at most" cells -- not the 0 the withdrawn CONFIRMED_ABSENT row asserted and
    #     not the 4,135 RESTORE called. Production's linear_otsu floor of 17.4 admits 23% of the whole
    #     tissue into the INS arm. Reported, never gating: the maintainer gave an order of magnitude,
    #     not a count, and the raw floor is not itself the positivity call.
    ins = per_donor[(per_donor.marker == "INS") & (per_donor.donor.astype(str) == "6566")]
    if not ins.empty:
        print("[validate] donor 6566 INS -- maintainer ground truth is a couple hundred cells:")
        for _, r in ins.sort_values("positive_n").iterrows():
            print(f"             {r.method:20s} floor {r.threshold:10.2f} -> {int(r.positive_n):>9,} "
                  f"cells  separation {r.separation:.3f}  fold {r.class_fold:8.1f}")

    # (2) Reproduce the floors production actually shipped, per donor, for every marker whose policy
    #     is linear_otsu. This is the identity that makes the audit comparable to the screen at all.
    ev = screen / "pair_evaluations.csv"
    if not ev.exists():
        print(f"[validate] SKIP screen comparison -- {ev} not found")
        return failures
    d = pd.read_csv(ev)
    shipped = pd.concat([
        d[["donor", "target", "target_input_floor", "target_input_floor_method"]]
        .rename(columns={"target": "marker", "target_input_floor": "floor",
                         "target_input_floor_method": "method"}),
        d[["donor", "reference", "reference_input_floor", "reference_input_floor_method"]]
        .rename(columns={"reference": "marker", "reference_input_floor": "floor",
                         "reference_input_floor_method": "method"}),
    ]).drop_duplicates(["donor", "marker"])
    shipped["donor"] = shipped["donor"].astype(str)

    mine = per_donor[per_donor.method == "linear_otsu"][["donor", "marker", "threshold"]]
    mine["donor"] = mine["donor"].astype(str)
    j = shipped[shipped.method == "linear_otsu"].merge(mine, on=["donor", "marker"], how="inner")
    if j.empty:
        print("[validate] SKIP screen comparison -- no overlapping linear_otsu donor-markers")
        return failures
    rel = ((j.threshold - j.floor).abs() / j.floor.replace(0, np.nan)).fillna(0)
    ok = bool((rel < 1e-6).all())
    print(f"[validate] linear_otsu floors vs screen_v12_production: {len(j)} donor-markers, "
          f"max relative delta {rel.max():.3e} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(j.assign(rel=rel).sort_values("rel", ascending=False).head(10).to_string(index=False))
    return failures + (0 if ok else 1)


# ------------------------------------------------------------------------------------------- driver

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", type=Path, default=REPO / "config.ini")
    p.add_argument("--outdir", type=Path, default=REPO / "data" / "floor_audit_v2")
    p.add_argument("--screen", type=Path,
                   default=REPO / "data/restore_pair_validation/screen_v12_production"
                                  "/ctrl-reference_only_div-p99")
    p.add_argument("--donors", nargs="*", default=None)
    p.add_argument("--markers", nargs="*", default=None, help="restrict for a quick check")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--figures", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--feature", choices=("raw", "redsea"), default="raw",
                   help="intensity the floor is estimated on. 'raw' is what production uses; "
                        "'redsea' measures the alternative (see audit_donor's docstring)")
    a = p.parse_args(argv)

    cfg = load_config(a.config)
    # Discovery goes through cfg.discover_donors, which filters 6457 and 6579 fail-closed; an explicit
    # list still goes through ensure_eligible_donors, which RAISES if either is named.
    donors = (ensure_eligible_donors(a.donors, context="marker floor audit") if a.donors
              else cfg.discover_donors(cfg.cells_dir))
    sample_root = cfg.restore_pair_validation_dir / "qupath_compartment_full_20donor"
    src = cfg.cells_dir if a.feature == "raw" else cfg.cells_redsea_dir
    a.outdir.mkdir(parents=True, exist_ok=True)

    print(f"[floors] {len(donors)} donors, feature={a.feature} ({src}) -> {a.outdir}", flush=True)
    if a.jobs > 1 and len(donors) > 1:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            results = list(ex.map(audit_donor, donors,
                                  [src] * len(donors), [sample_root] * len(donors)))
    else:
        results = [audit_donor(d, src, sample_root) for d in donors]

    per_donor = pd.concat([r[0] for r in results], ignore_index=True)
    hist = pd.concat([r[1] for r in results], ignore_index=True)
    if a.markers:
        per_donor = per_donor[per_donor.marker.isin(a.markers)]
        hist = hist[hist.marker.isin(a.markers)]
    per_donor = per_donor.sort_values(["marker", "method", "donor"]).reset_index(drop=True)
    per_donor.to_csv(a.outdir / "floor_audit_per_donor.csv", index=False)
    hist.to_csv(a.outdir / "floor_audit_histograms.csv", index=False)

    summary = summarise(per_donor)
    summary.to_csv(a.outdir / "floor_audit_summary.csv", index=False)

    if a.feature == "raw":
        failures = validate(per_donor, a.screen)
    else:
        # Check 2 asserts we reproduce the floors production SHIPPED, which are raw by construction;
        # it is meaningless here and is skipped rather than reported as a failure.
        print("[validate] feature=redsea -- skipping the shipped-floor identity check (raw-only)")
        failures = 0
        ins = per_donor[(per_donor.marker == "INS") & (per_donor.donor.astype(str) == "6566")]
        if not ins.empty:
            print("[validate] donor 6566 INS on REDSEA values (ground truth: a couple hundred cells):")
            for _, r in ins.sort_values("positive_n").iterrows():
                print(f"             {r.method:24s} floor {r.threshold:10.2f} -> {int(r.positive_n):>9,} "
                      f"cells  separation {r.separation:.3f}")

    (a.outdir / "manifest.json").write_text(json.dumps({
        "purpose": "Raw input floor audit for every viable marker; descriptive only, no method changed",
        "supersedes": "data/restore_pair_validation/marker_floor_audit_v1 (5 markers, 2026-07-21)",
        "source": {"raw_cells": str(cfg.cells_dir), "compartments": str(sample_root),
                   "screen": str(a.screen)},
        "feature": a.feature,
        "feature_source": str(src),
        "donors": list(donors),
        "excluded_donors": ["6457", "6579"],
        "excluded_markers": list(mt.EXCLUDED),
        "markers": sorted(per_donor.marker.unique().tolist()),
        "methods": list(METHODS),
        "validation_failures": failures,
    }, indent=2) + "\n")

    print(f"[floors] {per_donor.marker.nunique()} markers x {len(METHODS)} methods x "
          f"{len(donors)} donors = {len(per_donor):,} rows")
    if a.figures:
        from importlib import import_module
        import_module("floor_audit_figures").build(per_donor, summary, a.outdir, hist)
    if failures:
        print(f"[floors] {failures} VALIDATION FAILURE(S) -- do not trust this audit", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
