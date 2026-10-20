#!/usr/bin/env python
"""Per-cell genomic-DNA background for a Xenium bundle.

These runs are flagged for high genomic DNA control signal, and the flag is not cosmetic.
From the bundles' own metrics: `adjusted_genomic_control_probe_rate` 0.261 (0080870) and
0.306 (0081366) against a negative-control rate of 0.012, and
`estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts` of 126.0
and 128.6 against a median 433 / 403 transcripts per cell -- about 29% of a typical cell.

**The vendor's model, which this module reuses.** Every gene probe picks up gDNA at the same
rate a genomic control probe does. That reproduces exactly on 0080870:

    genomic_control_probe_counts_per_control_per_cell = 0.0252021   (vendor)
    0.0252021 x 5001 genes                            = 126.0       (vendor states 126.0)
    mean genomic_control_counts/cell 0.5292 / 21 probes = 0.0252     (same number)

So each cell carries its own background estimator, `genomic_control_counts_c / n_probes`.
`verify_against_vendor` re-derives the headline number and refuses on disagreement, because
the whole correction rests on this identity and it is checkable.

**Why the raw estimate is shrunk.** It is a Poisson draw with mean ~0.53 spread over 21
probes, so a single cell's estimate is far too noisy to use directly -- 59% of cells record
zero genomic counts, and those are not zero-background cells. Measured, the background scales
with cell size (mean gDNA/cell 0.212 -> 0.805 across area quintiles while gDNA per 100 um^2
stays ~1.0), so an area-conditioned expectation is the natural prior. The posterior is the
Gamma-Poisson conjugate mean; `shrink_a = n_probes` puts half the weight on each side.

**What this does NOT do.** It does not drop genes and it does not subtract counts. Genes whose
whole-section detection rate sits at or below background are flagged in
`var['below_gdna_background']` and kept -- a gene that is background-level across the section
can still be real in the one cluster that expresses it. Correction happens where decisions
happen: at the detection predicate and at group-level aggregation.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

#: Vendor key that states the per-genomic-control-probe rate.
VENDOR_RATE_KEY = "genomic_control_probe_counts_per_control_per_cell"
#: Vendor key for the headline false-positive estimate this module must reproduce.
VENDOR_FP_KEY = "estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts"
#: Relative tolerance when re-deriving the vendor's number. The identity is exact in
#: principle; anything above this means the panel or the column was misread.
VENDOR_TOL = 0.05
#: Area bins for the conditional prior. Deciles: enough to follow the trend, coarse enough
#: that every bin holds thousands of cells.
N_AREA_BINS = 10


class GdnaError(RuntimeError):
    """The background could not be established, or disagrees with the vendor."""


def panel_composition(bundle: pathlib.Path) -> dict[str, int]:
    """Count each target class in `gene_panel.json`.

    The probe count is read, never assumed: it sets the denominator of every rate below.
    """
    payload = json.loads((pathlib.Path(bundle) / "gene_panel.json").read_text())
    targets = payload.get("payload", {}).get("targets", payload.get("targets", []))
    counts: dict[str, int] = {}
    for t in targets:
        ty = t.get("type") or {}
        kind = ty.get("descriptor") if isinstance(ty, dict) else str(ty)
        counts[str(kind)] = counts.get(str(kind), 0) + 1
    if not counts.get("genomic_control"):
        raise GdnaError(f"{bundle}: gene_panel.json declares no genomic_control probes")
    return counts


def vendor_metrics(bundle: pathlib.Path) -> dict[str, float]:
    """`metrics_summary.csv` as floats, skipping the blank custom-panel columns."""
    df = pd.read_csv(pathlib.Path(bundle) / "metrics_summary.csv")
    row = df.iloc[0].to_dict()
    out: dict[str, float] = {}
    for k, v in row.items():
        try:
            out[str(k)] = float(v)
        except (TypeError, ValueError):
            continue
    return out


def area_prior(counts: np.ndarray, area: np.ndarray, n_bins: int = N_AREA_BINS
               ) -> tuple[np.ndarray, pd.DataFrame]:
    """E[genomic counts | cell area], evaluated per cell, plus the per-bin table.

    Ranked before binning so ties in area cannot collapse a bin, matching the convention
    `nullgate.detection_deciles` already uses for detection rate.
    """
    counts, area = np.asarray(counts, float), np.asarray(area, float)
    if counts.ndim != 1 or area.shape != counts.shape or not len(counts):
        raise GdnaError('counts and area must be nonempty matching vectors')
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any():
        raise GdnaError('genomic counts must be finite nonnegative integers')
    if n_bins < 1:
        raise GdnaError('n_bins must be positive')
    ok = np.isfinite(area) & (area > 0)
    if not ok.any():
        raise GdnaError('no positive finite cell areas for the area prior')
    # Small fixtures/cores cannot support ten independently estimated bins.
    n_bins = min(n_bins, max(1, int(ok.sum()) // 50))
    rank = pd.Series(np.where(ok, area, np.nan)).rank(method="average")
    try:
        b = pd.qcut(rank, n_bins, labels=False, duplicates="drop").to_numpy()
        if not np.isfinite(b[ok]).any():
            b[ok] = 0
    except ValueError as exc:                                          # pragma: no cover
        raise GdnaError(f"cannot bin cell_area into {n_bins} bins: {exc}") from exc
    b = np.where(ok, b, -1)
    mu = np.zeros(len(counts), float)
    global_mean = float(counts[ok].mean())
    rows = []
    for k in np.unique(b[b >= 0]):
        m = b == k
        mean_k = float(counts[m].mean())
        # A zero-count bin is insufficient evidence for exactly zero background.
        if mean_k == 0 and global_mean > 0:
            mean_k = global_mean
        mu[m] = mean_k
        rows.append({"bin": int(k), "n": int(m.sum()), "mean_area_um2": float(area[m].mean()),
                     "mean_genomic_counts": mean_k})
    if (b < 0).any():                       # cells with no area: fall back to the global mean
        mu[b < 0] = float(counts[b >= 0].mean()) if (b >= 0).any() else float(counts.mean())
    return mu, pd.DataFrame(rows)


def background_rates(genomic_counts: np.ndarray, area: np.ndarray, n_probes: int,
                     *, shrink_a: float | None = None
                     ) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, float]:
    """Per-cell expected gDNA counts **per gene**, raw and shrunk.

    Returns ``(lambda_shrunk, lambda_raw, area_table, shrink_a)``. The shrunk value is the
    Gamma-Poisson posterior mean with an area-conditioned prior; with ``shrink_a = n_probes``
    the observation and the prior carry equal weight, which is the default.
    """
    k = np.asarray(genomic_counts, dtype=float)
    a = float(n_probes if shrink_a is None else shrink_a)
    if not np.isfinite(n_probes) or n_probes <= 0 or int(n_probes) != n_probes:
        raise GdnaError('n_probes must be a positive integer')
    if not np.isfinite(a) or a <= 0:
        raise GdnaError('prior weight must be finite and positive')
    mu_counts, table = area_prior(k, np.asarray(area, dtype=float))
    r_prior = mu_counts / n_probes
    lam_raw = k / n_probes
    lam = (k + a * r_prior) / (n_probes + a)
    return lam, lam_raw, table, a


def verify_against_vendor(lam_mean: float, n_genes: int, metrics: dict[str, float],
                          *, tol: float = VENDOR_TOL) -> dict[str, float]:
    """Re-derive the vendor's headline false-positive number. Refuse on disagreement.

    This is the premise of the whole correction, and the vendor publishes the answer, so it
    is checked rather than trusted.
    """
    reproduced = lam_mean * n_genes
    stated = metrics.get(VENDOR_FP_KEY)
    out = {"reproduced_fp_per_cell": float(reproduced),
           "vendor_fp_per_cell": None if stated is None else float(stated),
           "vendor_rate": metrics.get(VENDOR_RATE_KEY)}
    if stated is None:
        raise GdnaError(f"metrics_summary.csv has no {VENDOR_FP_KEY!r}")
    if (not np.isfinite([lam_mean, n_genes, stated, tol]).all()
            or lam_mean < 0 or n_genes <= 0 or stated < 0 or tol <= 0):
        raise GdnaError('invalid background rate, panel size, vendor metric or tolerance')
    rel = abs(reproduced - stated) / max(abs(stated), 1e-12)
    out["relative_error"] = float(rel)
    if rel > tol:
        raise GdnaError(
            f"gDNA model does not reproduce the vendor: {reproduced:.2f} vs {stated:.2f} "
            f"false positives/cell ({100*rel:.1f}% off, tolerance {100*tol:.0f}%). "
            f"Check the genomic_control probe count and the genomic_control_counts column "
            f"before trusting any correction built on this rate.")
    return out


def below_background(detection_rate: np.ndarray, lam_mean: float) -> np.ndarray:
    """Genes whose section-wide detection rate is at or below pure background.

    ``1 - exp(-lambda)`` is P(at least one count from gDNA alone) for one gene in one cell.
    Flagged, never dropped.
    """
    return np.asarray(detection_rate, dtype=float) <= (1.0 - np.exp(-float(lam_mean)))


def summarise(bundle: pathlib.Path, *, shrink_a: float | None = None) -> dict:
    """Everything the pipeline stamps, computed from a bundle directory."""
    bundle = pathlib.Path(bundle)
    comp = panel_composition(bundle)
    metrics = vendor_metrics(bundle)
    n_probes, n_genes = comp["genomic_control"], comp.get("gene", 0)
    if not n_genes:
        raise GdnaError(f"{bundle}: gene_panel.json declares no gene targets")

    cells = pd.read_parquet(bundle / "cells.parquet",
                            columns=["cell_id", "genomic_control_counts", "cell_area"])
    lam, lam_raw, table, a = background_rates(
        cells["genomic_control_counts"].to_numpy(), cells["cell_area"].to_numpy(),
        n_probes, shrink_a=shrink_a)
    check = verify_against_vendor(float(lam_raw.mean()), n_genes, metrics)

    return {
        "bundle": str(bundle),
        "panel": comp,
        "n_probes": int(n_probes),
        "n_genes": int(n_genes),
        "n_cells": int(len(cells)),
        "shrink_a": float(a),
        "lambda_mean": float(lam.mean()),
        "lambda_raw_mean": float(lam_raw.mean()),
        "p_detect_by_gdna": float(1.0 - np.exp(-float(lam.mean()))),
        "frac_cells_zero_genomic": float((cells["genomic_control_counts"].to_numpy() == 0).mean()),
        "vendor_check": check,
        "area_table": table.to_dict("records"),
        "_lambda": lam,
        "_cell_id": cells["cell_id"].to_numpy(),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True, type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, default=None,
                    help="write the summary as JSON (the per-cell vector is not included)")
    ap.add_argument("--shrink-a", type=float, default=None,
                    help="prior weight; default = n genomic probes (50%% shrinkage)")
    args = ap.parse_args(argv)

    s = summarise(args.bundle, shrink_a=args.shrink_a)
    lam = s.pop("_lambda")
    s.pop("_cell_id")
    print(f"panel            : {s['panel']}")
    print(f"genomic probes   : {s['n_probes']}   genes: {s['n_genes']}")
    print(f"cells            : {s['n_cells']:,}  zero genomic counts: "
          f"{100*s['frac_cells_zero_genomic']:.1f}%")
    print(f"lambda (per gene): mean {s['lambda_mean']:.5f}  raw mean {s['lambda_raw_mean']:.5f}"
          f"  shrink_a {s['shrink_a']:.0f}")
    print(f"P(detect | gDNA) : {100*s['p_detect_by_gdna']:.2f}% of gene-cell pairs")
    v = s["vendor_check"]
    print(f"vendor check     : reproduced {v['reproduced_fp_per_cell']:.2f} vs stated "
          f"{v['vendor_fp_per_cell']:.2f} false positives/cell "
          f"({100*v['relative_error']:.2f}% off)  PASS")
    print(f"lambda quantiles : " + "  ".join(
        f"p{q}={np.percentile(lam, q):.4f}" for q in (5, 25, 50, 75, 95)))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(s, indent=2) + "\n")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
