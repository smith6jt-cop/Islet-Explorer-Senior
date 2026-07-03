#!/usr/bin/env python3
"""
RESTORE per-image intensity normalization of mutually-exclusive markers.

Applies RESTORE (Chang Lab / OHSU; github.com/smith6jt-cop/RESTORE, vendored at
external/RESTORE @ commit 38df59b) to the 59-plex protein single-cell data. **Pan_Cytokeratin**
(the abundant ~91% acinar majority) is the universal negative control for markers whose
background is brightness/autofluorescence-driven; the lone CD3e exception keeps a clean
exclusive immune reference (see below):
    exocrine:    Pan_Cytokeratin (acinar/ductal)         -- background from Vimentin
    mesenchyme:  Vimentin (stroma/fibroblast)            -- background from acinar (Pan_CK)
    endocrine:   INS (beta), GCG (alpha)                 -- background from acinar (Pan_CK)
    endothelium: CD31                                    -- background from acinar (Pan_CK)
    muscle:      SMA (smooth muscle/pericyte)            -- background from acinar (Pan_CK)
    immune:      CD20 (B), CD163 (macrophage)            -- background from acinar (Pan_CK)
                 CD3e (T)                                -- background from CD163 (clean exclusive)
(CD99 was dropped -- broadly expressed, not lineage-exclusive in this panel.)

RESTORE pairs each target marker with a mutually-exclusive *reference* marker. Cells
positive for the reference are, by biology, negative for the target, so the target
intensity within the reference-positive population is pure background/autofluorescence.
Per image, three 2-cluster models (KMeans, GMM, SSC) split background vs signal and the
threshold is set to  mean + 3*sigma  of the target intensity in the target-negative
cluster. This corrects per-image autofluorescence drift AND the per-cell background that
scales with cell brightness/autofluorescence, both of which bias fixed-cutoff gating.

WHY Pan_Cytokeratin for CD20/CD163 (the fix, Jun 2026): their background SCALES WITH CELL
BRIGHTNESS (CD20 MFI vs cytokeratin/DAPI Pearson r~0.5), so bright acinar cells carry elevated
CD20/CD163 background. The negative control must SPAN that bright range. Pan_Cytokeratin+
acinar cells are abundant, immune-negative, and bright (3.4x brighter in Ker8_18 than CD3e+
T cells) -> they represent the brightness-scaling background. The previous "immune-vs-immune"
references (CD20<-CD3e, CD163<-CD3e) used the dim, rare T-cell population as the negative
control, so RESTORE never saw the bright-acinar background -> threshold too low -> ~55%+ of the
Immune broad-lineage class in donor 6539 were bright acinar FALSE POSITIVES. Pan_Cytokeratin
removes them (CD20+ -87%) AND sharpens the disease signal (CD20/CD163 T1D/ND ratio 1.3x->~2.5x).

WHY CD3e is the EXCEPTION (keeps CD163, verified on the real cohort): CD3e was never a
bright-acinar false-positive (strong, specific; norm~3). Pan_Cytokeratin is POSITIVELY
correlated with CD3e (autofluorescence coupling, ref-qc pearson ~0.39), which drags the SSC
threshold too low and DESTROYS the validated T-cell infiltration biology (T-gain collapses
4.5x->0.9x in the real SSC cohort — a crude mean+3sigma approximation had wrongly suggested it
improved). CD3e keeps CD163 (macrophages — cleanly mutually exclusive, negatively correlated),
preserving ND 0.9% -> T1D 4.1% (4.5x). Principle: Pan_CK references the markers whose background
is brightness/autofluorescence-driven; a clean specific marker keeps a clean exclusive reference.
See `reference_exclusivity_qc()` for the data-driven audit (svd_ratio, pearson sign, abundance).

Deliberate deviation from the paper's literal red-flag: RESTORE (PMC7062831) recommends
rejecting POSITIVELY-correlated "mutually-exclusive" pairs as non-specific staining. Here
the positive CD20<->Pan_CK correlation IS the autofluorescence/brightness coupling we WANT
RESTORE to learn from (not antibody cross-reactivity), so Pan_Cytokeratin is kept as the
immune reference despite the positive correlation. The QC reports it; it does not auto-reject.

Operates on RAW QuPath mean intensities (RESTORE's idx_select uses an absolute >50 floor
and clusters on raw 2D intensities -- do NOT log-transform the input).

Run in the `scvi-env` conda env. Requires `spams` for the SSC model:
    pip install spams-bin        # prebuilt wheel (cp313 etc.)

Pipeline:
  data/cells/donor_id=*/data_0.parquet  (raw MFI, 23.3M cells)
    -> per-image subsample            (RESTORE input: scene=donor, batch=0)
    -> Normalization(save_figs=False).run()      -> data/restore/threshs.pkl
    -> flatten                                    -> data/restore_thresholds.csv
    -> apply per-image thresholds to ALL cells (chosen model, default SSC):
         <m>_pos   = raw >= thr
         <m>_norm  = raw / thr           (threshold -> 1.0, harmonized across images)
         <m>_log2r = log2((raw + 1) / thr)
       -> data/restore_gated/donor_id=<id>/data_0.parquet
    -> QC plots                                   -> data/restore/qc/*.png

Usage:
  python scripts/senior/restore_normalize.py                    # full run, all 15 images
  python scripts/senior/restore_normalize.py --limit-scenes 1 --skip-apply   # dry run
  python scripts/senior/restore_normalize.py --model GMM        # choose gating model
  python scripts/senior/restore_normalize.py --marker-pairs 'CD3e:CD163,CD99:CD163'
"""
from __future__ import annotations

import argparse
import glob
import pickle
import re
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BASE = Path(__file__).resolve().parents[2] / "data"
VENDOR = Path(__file__).resolve().parents[2] / "external" / "RESTORE" / "python_code"

# Default mutually-exclusive pairs: [target, reference].
# The reference's positive cells define the target's background. Pan_Cytokeratin (the abundant
# ~91% acinar majority) is the UNIVERSAL negative control: abundant, cleanly target-negative for
# every non-epithelial lineage, AND bright (so it spans the brightness-scaling background). Every
# lineage marker references it, except Pan_Cytokeratin itself (reciprocal Vimentin).
DEFAULT_MARKER_PAIRS = [
    # --- exocrine / mesenchyme reciprocal pair ---
    ["Pan_Cytokeratin", "Vimentin"],         # exocrine epithelium <- mesenchyme
    ["Vimentin",        "Pan_Cytokeratin"],  # mesenchyme/stroma   <- acinar
    # --- endocrine / endothelium / muscle <- acinar (validated; recovers beta-loss) ---
    ["INS",             "Pan_Cytokeratin"],  # endocrine beta      <- acinar
    ["GCG",             "Pan_Cytokeratin"],  # endocrine alpha     <- acinar
    ["SST",             "Pan_Cytokeratin"],  # endocrine delta     <- acinar (sparse, specific like INS/GCG)
    ["CD31",            "Pan_Cytokeratin"],  # endothelium         <- acinar
    ["SMA",             "Pan_Cytokeratin"],  # smooth muscle/pericyte <- acinar
    # --- immune (Jun 2026 fix). CD20/CD163 were the false-positive problem: their background is
    #     bright-acinar AUTOFLUORESCENCE (MFI scales with cell brightness), so the dim T-cell
    #     reference missed it -> bright acinar cells false-positive. Pan_Cytokeratin+ acinar is the
    #     bright, abundant negative control that captures it -> false positives removed AND sharper
    #     ND->T1D signal (CD20 T1D/ND 1.3x->2.5x, CD163 1.3x->2.4x). CD3e is the EXCEPTION: it was
    #     already clean/specific (never a bright-acinar false-positive), and Pan_Cytokeratin is
    #     POSITIVELY correlated with CD3e (autofluorescence coupling) which drags its SSC threshold
    #     too low and DESTROYS the validated T-cell infiltration biology (real cohort: T-gain
    #     4.5x->0.9x). So CD3e keeps a cleanly mutually-exclusive immune reference (CD163, negatively
    #     correlated) that preserves the T-gain. Principle: Pan_CK references the markers whose
    #     background is brightness/autofluorescence-driven; a clean specific marker keeps a clean
    #     exclusive reference. ---
    ["CD3e",            "CD163"],            # T cells             <- macrophages (unchanged; Pan_CK breaks T-gain)
    ["CD20",            "Pan_Cytokeratin"],  # B cells             <- acinar (was CD3e)
    ["CD163",           "Pan_Cytokeratin"],  # macrophages         <- acinar (was CD3e)
]

MODEL_COLORS = {"KMeans": "magenta", "GMM": "blue", "SSC": "green"}


# --------------------------------------------------------------------------- #
# Vendored RESTORE import (headless)
# --------------------------------------------------------------------------- #
def _install_headless_stubs():
    """RESTORE/RESTORE.py imports holoviews, selenium, and tqdm.notebook at module
    load, and builds hv.Curve/hv.Scatter/hv.Overlay objects during threshold
    computation (line ~272) even when save_figs=False. Inject lightweight stubs so the
    vendored source runs headless and *unmodified* -- we never render or save figures.
    """
    import tqdm

    class _HVObj:  # supports .opts(), .cols(), and the `*` overlay operator
        def opts(self, *a, **k):
            return self

        def cols(self, *a, **k):
            return self

        def __mul__(self, other):
            return self

        def __rmul__(self, other):
            return self

    def _factory(*a, **k):
        return _HVObj()

    hv = types.ModuleType("holoviews")
    for name in ("Curve", "Scatter", "VLine", "Overlay", "Layout"):
        setattr(hv, name, _factory)
    hv.save = lambda *a, **k: None
    hv.extension = lambda *a, **k: None
    sys.modules.setdefault("holoviews", hv)

    sel = types.ModuleType("selenium")
    swd = types.ModuleType("selenium.webdriver")
    sel.webdriver = swd
    sys.modules.setdefault("selenium", sel)
    sys.modules.setdefault("selenium.webdriver", swd)

    tn = types.ModuleType("tqdm.notebook")
    tn.tqdm = tqdm.tqdm
    sys.modules.setdefault("tqdm.notebook", tn)


def patch_restore_sigma(Normalization, sigma_map, default=3.0):
    """RESTORE hardcodes sigma_weight=3 (threshold = neg_mean + 3·neg_σ) — right for sparse
    markers but too CONSERVATIVE for abundant structural markers (Pan_Cytokeratin/Vimentin),
    which then flag only ~18% positive vs the true acinar/stromal majority. Monkeypatch the two
    threshold helpers so the target marker's sigma_weight comes from `sigma_map` (per target
    marker), leaving the vendored source untouched. marker_pair[0] is the target."""
    import functools
    if not sigma_map:
        return
    for fname in ("get_GMM_thresh", "get_clustering_thresh"):
        orig = getattr(Normalization, fname)

        def make(orig):
            @functools.wraps(orig)
            def wrapper(self, data, marker_pair, model, sigma_weight, color, batch):
                sw = sigma_map.get(marker_pair[0], default)
                return orig(self, data, marker_pair, model, sw, color, batch)
            return wrapper
        setattr(Normalization, fname, make(orig))
    print(f"[sigma] per-marker sigma_weight override: {sigma_map} (default {default})")


def patch_restore_cluster_cap(Normalization, cap, seed=0):
    """Subsample-AFTER-idx_select. RESTORE pre-subsampled each image to a flat 15k BEFORE
    idx_select, so a pair whose idx_select cloud is a small fraction of the image got a tiny,
    unstable negative-control and a wobbly threshold (e.g. CD3e's threshold swung ~10x and only
    converged near full data). We instead feed FULL per-image data (build_restore_input, no
    pre-subsample) so idx_select adapts to whatever the donor actually contains — making NO
    assumption about which markers are rare or abundant, and respecting real donor biology (e.g.
    a pancreatitis donor with very high T-cell content) — then cap the resulting idx_select cloud
    to `cap` cells HERE, right before clustering, so SSC/GMM/KMeans stay tractable. get_GMM_thresh
    and get_clustering_thresh receive the post-idx_select `marker_pair_data` as their first arg."""
    if not cap:
        return
    import functools
    rng = np.random.default_rng(seed)
    for fname in ("get_GMM_thresh", "get_clustering_thresh"):
        orig = getattr(Normalization, fname)

        def make(orig):
            @functools.wraps(orig)
            def wrapper(self, data, marker_pair, model, sigma_weight, color, batch):
                if data is not None and len(data) > cap:
                    data = data[rng.choice(len(data), size=cap, replace=False)]
                return orig(self, data, marker_pair, model, sigma_weight, color, batch)
            return wrapper
        setattr(Normalization, fname, make(orig))
    print(f"[cluster-cap] idx_select cloud capped to {cap:,} cells before clustering "
          "(subsample-AFTER-idx_select; full image feeds idx_select)")


def patch_restore_no_kde(Normalization):
    """RESTORE.py:221 computes a scipy gaussian_kde over the FULL post-idx_select cloud purely
    to color a scatter it only saves when save_figs=True. gaussian_kde evaluation is O(n^2); on
    full-image idx_select clouds (~1.3M cells) that stalls for HOURS per marker pair (holoviews
    is already stubbed, so the scatter itself is a no-op — the KDE was the sole remaining cost).
    With save_figs=False the density is discarded, so replace gaussian_kde with a length-
    preserving zero-density stub. Thresholds are unaffected (they come from the clustering, which
    runs on a representative sample via patch_restore_cluster_cap)."""
    import sys
    R = sys.modules[Normalization.__module__]        # the vendored RESTORE.RESTORE module
    R.gaussian_kde = lambda dataT: (lambda q: np.zeros(np.asarray(q).shape[-1]))
    print("[kde-stub] figure-only gaussian_kde neutralized (was O(n^2) on the full idx_select cloud)")


def import_restore():
    _install_headless_stubs()
    if str(VENDOR) not in sys.path:
        sys.path.insert(0, str(VENDOR))
    try:
        from RESTORE import Normalization, nested_dict  # noqa
    except ModuleNotFoundError as e:
        if "spams" in str(e):
            raise SystemExit(
                "[fatal] SSC model needs `spams`. Install with:  pip install spams-bin"
            )
        raise SystemExit(
            f"[fatal] cannot import RESTORE from {VENDOR} ({e}). "
            "Clone it:  git clone https://github.com/smith6jt-cop/RESTORE.git external/RESTORE"
        )
    return Normalization, nested_dict


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #
def donor_files(cells_dir, limit=None, donors=None):
    files = sorted(glob.glob(str(Path(cells_dir) / "donor_id=*" / "*.parquet")))
    pairs = [(re.search(r"donor_id=([^/]+)", f).group(1), f) for f in files]
    if donors:
        keep = {str(d) for d in donors}
        pairs = [p for p in pairs if p[0] in keep]
    return pairs[:limit] if limit else pairs


def build_restore_input(cells_dir, markers, subsample, seed, limit_scenes=None, donors=None):
    """Per-image (per-donor) subsample of raw MFI for threshold estimation.

    Returns a DataFrame with columns: scene (donor id), batch (0), + one column per
    marker. RESTORE only uses data.scene, data.batch, and data[marker_pair].
    """
    rng = np.random.default_rng(seed)
    frames = []
    for donor, f in donor_files(cells_dir, limit_scenes, donors):
        df = pd.read_parquet(f, columns=markers)
        n_raw = len(df)
        # REDSEA leaves a few "edge-cell" artifacts that are NaN across ALL markers (whole-cell
        # correction failed). Those must go, because RESTORE's idx_select (RESTORE.py
        # get_marker_pair_thresh) uses np.quantile, which is NOT NaN-aware — a single sampled NaN
        # makes the quantile NaN, empties the selection, and KMeans crashes on 0 samples.
        # We drop ONLY all-marker-NaN edge cells, NOT cells that merely lack one measurement:
        # raw QuPath marker means have no NaN, so a PARTIAL-NaN cell would be unexpected — warn
        # (loudly) rather than silently discard a partially-valid cell.
        allnan = df.isna().all(axis=1)
        partial = df.isna().any(axis=1) & ~allnan
        if partial.any():
            print(f"[input] WARN {donor}: {int(partial.sum())} PARTIAL-NaN cells (some markers "
                  "valid) — investigate upstream; they are dropped so the quantile stays finite")
        df = df[~df.isna().any(axis=1)]
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        df = df.astype("float64").reset_index(drop=True)
        df["scene"] = donor
        df["batch"] = 0
        frames.append(df)
        drop = n_raw - n
        print(f"[input] image {donor}: {n_raw:,} cells ({drop:,} NaN dropped) -> {len(df):,} sampled")
    data = pd.concat(frames, ignore_index=True)
    print(f"[input] total RESTORE input: {len(data):,} cells / {data['scene'].nunique()} images")
    return data


# --------------------------------------------------------------------------- #
# Reference exclusivity QC (data-driven audit of the mutually-exclusive pairs)
# --------------------------------------------------------------------------- #
def _idx_select(T, R, ratio_x=0.75, ratio_y=0.5):
    """Replicate RESTORE's idx_select prefilter (RESTORE.py get_marker_pair_thresh): keep cells
    with notable expression in at least one of the pair (>50 floor + quantile box)."""
    return (((T > 50) & (R > np.quantile(R, ratio_y))) |
            ((R > 50) & (T > np.quantile(T, ratio_x))))


def reference_exclusivity_qc(cells_dir, pairs, subsample, seed, out_dir, donors=None):
    """Data-driven audit of each [target, reference] pair (RESTORE paper, PMC7062831 §Methods).
    Over a pooled subsample, per pair reports:
      - svd_ratio r = σ2/σ1 on the log1p + z-scored, idx_select-prefiltered 2D [T,R] cloud
        (r→1 circular/mutually-exclusive; r→0 elongated/correlated),
      - pearson on the log1p cloud (sign: +corr = autofluorescence/brightness coupling OR antibody
        cross-reactivity; −corr = clean exclusivity),
      - reference abundance frac(R>50)/frac(R>200) — does the reference span the BRIGHT range that
        carries the brightness-scaling background?
      - target background in R-positive cells (median/p95) + the implied mean+3σ threshold,
      - median dist_islet of R-positive cells (peri-islet contamination check — flags e.g. why
        Vimentin is a poor INS reference: its R+ cells sit near islets and carry INS spillover).
    AUDITS the curated references; does NOT auto-select. The paper's literal "reject +correlation"
    rule would wrongly reject Pan_Cytokeratin for the immune markers — that +correlation IS the
    autofluorescence we WANT the negative control to capture (deliberate, documented deviation).
    Writes <out_dir>/reference_selection_qc.csv + <out_dir>/qc/reference_selection_qc.png."""
    rng = np.random.default_rng(seed)
    markers = sorted({m for pair in pairs for m in pair})
    frames = []
    for donor, f in donor_files(cells_dir, None, donors):
        df = pd.read_parquet(f, columns=["object_id"] + markers)
        df = df.dropna(subset=markers)          # drop REDSEA NaN edge-cells (see build_restore_input)
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        df = df.reset_index(drop=True)
        df["object_id"] = df["object_id"].astype(str)
        try:  # dist_islet lives in canonical data/cells (cells_redsea lacks it); join by object_id
            raw = pd.read_parquet(BASE / "cells" / f"donor_id={donor}" / "data_0.parquet",
                                  columns=["object_id", "dist_islet"])
            raw["object_id"] = raw["object_id"].astype(str)
            df = df.merge(raw, on="object_id", how="left")
        except Exception:
            df["dist_islet"] = np.nan
        frames.append(df)
    pool = pd.concat(frames, ignore_index=True)

    rows = []
    for target, ref in pairs:
        T = pool[target].to_numpy(np.float64); R = pool[ref].to_numpy(np.float64)
        sel = _idx_select(T, R)
        lt = np.log1p(T[sel]); lr = np.log1p(R[sel])
        M = np.vstack([(lt - lt.mean()) / (lt.std() or 1.0),
                       (lr - lr.mean()) / (lr.std() or 1.0)]).T
        s = np.linalg.svd(M, compute_uv=False)
        r = float(s[1] / s[0]) if s[0] > 0 else np.nan
        pear = float(np.corrcoef(np.log1p(T), np.log1p(R))[0, 1])
        refpos = R >= np.quantile(R, 0.75)          # reference-positive = the negative control
        bg = T[refpos]
        di = pool["dist_islet"].to_numpy(np.float64)[refpos]
        rows.append(dict(
            target=target, reference=ref,
            svd_ratio=round(r, 3), pearson_log=round(pear, 3),
            ref_frac_gt50=round(float((R > 50).mean()), 3),
            ref_frac_gt200=round(float((R > 200).mean()), 3),
            bg_median=round(float(np.median(bg)), 1),
            bg_p95=round(float(np.percentile(bg, 95)), 1),
            thr_mean3sd=round(float(bg.mean() + 3 * bg.std()), 1),
            refpos_dist_islet_med=round(float(np.nanmedian(di)), 1) if np.isfinite(di).any() else np.nan,
            flag_not_exclusive=bool(np.isfinite(r) and r < 0.4),
            flag_pos_corr=bool(pear > 0.3),
        ))
    qc = pd.DataFrame(rows)
    out_csv = Path(out_dir) / "reference_selection_qc.csv"
    qc.to_csv(out_csv, index=False)
    print(f"[ref-qc] wrote {out_csv}")
    print(qc.to_string(index=False))

    try:
        fig, ax = plt.subplots(1, 3, figsize=(18, 0.45 * len(qc) + 2))
        y = np.arange(len(qc)); lab = [f"{t}<-{rf}" for t, rf in zip(qc.target, qc.reference)]
        ax[0].barh(y, qc.svd_ratio, color="#3366CC"); ax[0].axvline(0.4, color="r", ls="--")
        ax[0].set_yticks(y); ax[0].set_yticklabels(lab); ax[0].invert_yaxis()
        ax[0].set_xlabel("SVD ratio σ2/σ1  (→1 = mutually exclusive)"); ax[0].set_title("Exclusivity")
        ax[1].barh(y, qc.pearson_log, color="#EE7733")
        ax[1].axvline(0.3, color="r", ls="--"); ax[1].axvline(0, color="k", lw=0.6)
        ax[1].set_yticks(y); ax[1].set_yticklabels([]); ax[1].invert_yaxis()
        ax[1].set_xlabel("Pearson (log cloud)"); ax[1].set_title("+corr = autofluor / cross-react")
        ax[2].barh(y, qc.ref_frac_gt200, color="#117733")
        ax[2].set_yticks(y); ax[2].set_yticklabels([]); ax[2].invert_yaxis()
        ax[2].set_xlabel("frac(reference > 200 MFI)"); ax[2].set_title("Reference spans bright range?")
        fig.suptitle("RESTORE reference audit — Pan_Cytokeratin is the abundant, bright, "
                     "mutually-exclusive universal negative control", fontsize=12)
        fig.tight_layout(); fig.savefig(Path(out_dir) / "qc" / "reference_selection_qc.png",
                                        dpi=130, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print("[ref-qc] figure failed:", e)
    return qc


# --------------------------------------------------------------------------- #
# Thresholds
# --------------------------------------------------------------------------- #
def flatten_threshs(threshs):
    """threshs[batch][scene][target_marker][model] = value  ->  tidy DataFrame.

    The 'global' batch (RESTORE's all-cells pass) duplicates the real per-batch pass
    when there is a single batch; we keep the numeric batch and drop 'global'.
    """
    batches = [b for b in threshs if b != "global"] or list(threshs)
    rows = []
    for batch in batches:
        for scene in threshs[batch]:
            for marker in threshs[batch][scene]:
                for model, val in threshs[batch][scene][marker].items():
                    rows.append(
                        {"image": scene, "marker": marker, "model": model,
                         "threshold": float(val)}
                    )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Apply to all cells
# --------------------------------------------------------------------------- #
def build_chosen_lut(thr_df, model, robust=True, factor=3.0):
    """Lookup (image, marker) -> threshold for the chosen model.

    RESTORE is fit per image, but on a few images SSC/GMM/KMeans return a degenerate
    threshold far in the tail (~0% positive), which makes that image's per-cell
    normalization meaningless. With `robust`, any per-image threshold that is an
    outlier (> factor x or < 1/factor x the cohort median for that marker, where the
    cohort median is taken over the non-outlier images) falls back to that cohort
    median. The raw per-model thresholds in restore_thresholds.csv are left untouched;
    only the thresholds actually applied are imputed, and every imputation is logged.
    """
    sub = thr_df[thr_df.model == model]
    lut = {(r.image, r.marker): float(r.threshold) for r in sub.itertuples()}
    imputed = {}
    if not robust:
        return lut, imputed
    for marker in sorted(sub.marker.unique()):
        vals = sub[sub.marker == marker].set_index("image").threshold.astype(float)
        med0 = vals.median()
        ok = vals[(vals <= factor * med0) & (vals >= med0 / factor)]
        med = ok.median() if len(ok) else med0
        for img, v in vals.items():
            if not np.isfinite(v) or v <= 0 or v > factor * med or v < med / factor:
                lut[(img, marker)] = float(med)
                imputed[(img, marker)] = (float(v), float(med))
    return lut, imputed


def apply_thresholds(cells_dir, markers, lut, imputed, gated_dir, limit_scenes=None, donors=None):
    """Stream every donor parquet, gate + normalize the 4 markers, write partitions."""
    frac_rows = []
    for donor, f in donor_files(cells_dir, limit_scenes, donors):
        df = pd.read_parquet(f, columns=["object_id"] + markers)
        out = pd.DataFrame({"object_id": df["object_id"].astype("string"),
                            "donor_id": donor})
        for m in markers:
            thr = lut.get((donor, m))
            if thr is None or not np.isfinite(thr) or thr <= 0:
                print(f"[apply] WARN image {donor} {m}: missing/invalid threshold ({thr}); skipped")
                continue
            raw = df[m].to_numpy(dtype="float64")
            pos = raw >= thr
            out[f"{m}_pos"] = pos
            out[f"{m}_norm"] = (raw / thr).astype("float32")
            out[f"{m}_log2r"] = np.log2((raw + 1.0) / thr).astype("float32")
            frac_rows.append({"image": donor, "marker": m, "threshold": float(thr),
                              "imputed": (donor, m) in imputed,
                              "n": int(len(raw)), "n_pos": int(pos.sum()),
                              "frac_pos": float(pos.mean())})
        dst = Path(gated_dir) / f"donor_id={donor}"
        dst.mkdir(parents=True, exist_ok=True)
        out.to_parquet(dst / "data_0.parquet", index=False)
        print(f"[apply] image {donor}: wrote {len(out):,} rows -> {dst}/data_0.parquet")
    return pd.DataFrame(frac_rows)


# --------------------------------------------------------------------------- #
# QC plots
# --------------------------------------------------------------------------- #
def qc_histograms(data, thr_df, markers, qc_dir):
    for scene in sorted(data["scene"].unique()):
        sd = data[data["scene"] == scene]
        ncol = 2
        nrow = int(np.ceil(len(markers) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 4 * nrow))
        axes = np.atleast_1d(axes).ravel()
        for i, m in enumerate(markers):
            ax = axes[i]
            x = sd[m].to_numpy()
            hi = np.quantile(x, 0.995) if len(x) else 1.0
            ax.hist(np.clip(x, 0, hi), bins=80, color="0.75")
            ax.set_yscale("log")
            for model, col in MODEL_COLORS.items():
                r = thr_df[(thr_df.image == scene) & (thr_df.marker == m)
                           & (thr_df.model == model)]
                if len(r):
                    t = float(r.threshold.iloc[0])
                    ax.axvline(t, color=col, lw=1.5, label=f"{model}={t:.0f}")
            ax.set_title(m)
            ax.set_xlabel("raw MFI")
            ax.legend(fontsize=8)
        for j in range(len(markers), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"RESTORE thresholds — image {scene}")
        fig.tight_layout()
        fig.savefig(Path(qc_dir) / f"{scene}_hist.png", dpi=120)
        plt.close(fig)


def _annotated_heatmap(piv, title, cbar_label, path, fmt="{:.0f}", cmap="viridis"):
    fig, ax = plt.subplots(figsize=(1.0 * piv.shape[1] + 3, 0.6 * piv.shape[0] + 2))
    im = ax.imshow(piv.values.astype(float), aspect="auto", cmap=cmap)
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=45, ha="right")
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center",
                        color="w", fontsize=7)
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def qc_threshold_heatmap(thr_df, markers, model, qc_dir):
    piv = (thr_df[thr_df.model == model]
           .pivot(index="marker", columns="image", values="threshold")
           .reindex(markers))
    _annotated_heatmap(piv, f"RESTORE per-image thresholds ({model})",
                       f"{model} threshold (raw MFI)",
                       Path(qc_dir) / "threshold_heatmap.png")


def qc_positive_fractions(frac_df, markers, qc_dir):
    piv = (frac_df.pivot(index="marker", columns="image", values="frac_pos")
           .reindex(markers) * 100.0)
    _annotated_heatmap(piv, "RESTORE positive fraction per image (chosen model)",
                       "% cells positive",
                       Path(qc_dir) / "positive_fractions.png",
                       fmt="{:.1f}", cmap="magma")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-dir", type=Path, default=BASE / "cells")
    ap.add_argument("--out-dir", type=Path, default=BASE / "restore",
                    help="threshs.pkl + qc/ + positive_fractions.csv")
    ap.add_argument("--thresh-csv", type=Path, default=BASE / "restore_thresholds.csv")
    ap.add_argument("--gated-dir", type=Path, default=BASE / "restore_gated",
                    help="partitioned per-cell gated + normalized output")
    ap.add_argument("--subsample", type=int, default=15000,
                    help="cap on the post-idx_select cloud fed to clustering per image×pair "
                         "(SSC-tractable). idx_select itself runs on the FULL image, so this is a "
                         "subsample-AFTER-idx_select cap, not a pre-subsample of the whole image")
    ap.add_argument("--model", default="SSC", choices=["GMM", "KMeans", "SSC"],
                    help="model whose threshold drives gating + normalization. "
                         "SSC is the default: on this panel GMM/KMeans place the "
                         "abundant lineage markers (CD3e/CD20/CD163) at/above the "
                         "99.9th intensity pct (near-zero positives); SSC lands at the "
                         "background->signal shoulder. All 3 are saved in the CSV.")
    ap.add_argument("--marker-pairs", default=None,
                    help="override pairs, e.g. 'CD3e:CD163,CD99:CD163'")
    ap.add_argument("--marker-sigma", default=None,
                    help="per-target-marker sigma_weight override (threshold = neg_mean + sigma·neg_σ), "
                         "e.g. 'Pan_Cytokeratin:1.0,Vimentin:1.5' to capture the abundant structural "
                         "majority that mean+3σ misses (default 3.0 for all)")
    ap.add_argument("--limit-scenes", type=int, default=None,
                    help="process only the first N images (dry run)")
    ap.add_argument("--donors", nargs="*", default=None,
                    help="restrict to specific donor ids (e.g. --donors 6539); fits per-image "
                         "thresholds on just those images (use with --no-robust for a single image)")
    ap.add_argument("--skip-apply", action="store_true",
                    help="thresholds + QC only; do not gate/normalize all cells")
    ap.add_argument("--robust", action=argparse.BooleanOptionalAction, default=True,
                    help="impute degenerate per-image thresholds (outliers vs cohort "
                         "median) with the cohort median before applying (default on; "
                         "use --no-robust for pure per-image RESTORE thresholds)")
    ap.add_argument("--robust-factor", type=float, default=3.0,
                    help="an image threshold is an outlier if >factor x or <1/factor x "
                         "the cohort median for that marker")
    ap.add_argument("--reuse-threshs", action="store_true",
                    help="reuse an existing threshs.pkl in --out-dir (skip RESTORE fit; "
                         "deterministic, so safe for re-applying with different options)")
    ap.add_argument("--ref-qc", action=argparse.BooleanOptionalAction, default=True,
                    help="write the data-driven reference-exclusivity audit "
                         "(reference_selection_qc.csv + figure); --no-ref-qc to skip")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    pairs = (DEFAULT_MARKER_PAIRS if not a.marker_pairs
             else [p.split(":") for p in a.marker_pairs.split(",")])
    markers = sorted({m for pair in pairs for m in pair})
    print(f"[config] pairs (target<-reference): {pairs}")
    print(f"[config] markers={markers}  chosen_model={a.model}  subsample={a.subsample}")

    Normalization, _ = import_restore()
    if a.marker_sigma:
        sigma_map = {k: float(v) for k, v in (p.split(":") for p in a.marker_sigma.split(","))}
        patch_restore_sigma(Normalization, sigma_map)
    # subsample-AFTER-idx_select: cap the post-idx_select cloud at --subsample cells, and feed
    # FULL per-image data to idx_select (pass subsample=None to build_restore_input) so the
    # negative-control is representative for every pair regardless of marker/cell frequencies.
    patch_restore_cluster_cap(Normalization, a.subsample, a.seed)
    patch_restore_no_kde(Normalization)   # kill the figure-only O(n^2) KDE on the full cloud
    qc_dir = a.out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    print("[step1] building RESTORE input (full images; idx_select runs on all cells) ...")
    data = build_restore_input(a.cells_dir, markers, None, a.seed, a.limit_scenes, a.donors)

    pkl = a.out_dir / "threshs.pkl"
    if a.reuse_threshs and pkl.exists():
        print(f"[step2] --reuse-threshs: loading existing {pkl}")
    else:
        print("[step2] running RESTORE (KMeans + GMM + SSC, save_figs=False) ...")
        norm = Normalization(data, pairs, save_dir=str(a.out_dir), save_figs=False)
        norm.run()
        print(f"[step2] wrote {pkl}")

    with open(pkl, "rb") as fh:
        threshs = pickle.load(fh)
    thr_df = flatten_threshs(threshs)
    thr_df["chosen"] = thr_df["model"] == a.model
    thr_df.sort_values(["marker", "image", "model"]).to_csv(a.thresh_csv, index=False)
    print(f"[step3] wrote {a.thresh_csv}  ({len(thr_df)} rows)")
    chosen = (thr_df[thr_df.model == a.model]
              .pivot(index="marker", columns="image", values="threshold"))
    print(f"[step3] {a.model} thresholds per image (raw MFI):\n{chosen.round(1).to_string()}")

    print("[step4] QC plots ...")
    qc_histograms(data, thr_df, markers, qc_dir)
    qc_threshold_heatmap(thr_df, markers, a.model, qc_dir)
    if a.ref_qc:
        print("[step4] reference-exclusivity audit ...")
        reference_exclusivity_qc(a.cells_dir, pairs, a.subsample, a.seed, a.out_dir, a.donors)
    print(f"[step4] wrote QC -> {qc_dir}")

    if a.skip_apply:
        print("[done] --skip-apply: thresholds + QC only")
        return 0

    chosen_lut, imputed = build_chosen_lut(thr_df, a.model, a.robust, a.robust_factor)
    if imputed:
        print(f"[robust] imputed {len(imputed)} degenerate {a.model} thresholds "
              f"(>{a.robust_factor}x cohort median) with the cohort median:")
        for (img, mk), (orig, med) in sorted(imputed.items()):
            print(f"  image {img} {mk}: {orig:.0f} -> {med:.0f}")
    elif a.robust:
        print(f"[robust] no outlier {a.model} thresholds to impute")

    print("[step5] applying thresholds to ALL cells ...")
    frac_df = apply_thresholds(a.cells_dir, markers, chosen_lut, imputed,
                               a.gated_dir, a.limit_scenes, a.donors)
    frac_df.sort_values(["marker", "image"]).to_csv(
        a.out_dir / "positive_fractions.csv", index=False)
    qc_positive_fractions(frac_df, markers, qc_dir)
    print(f"[step5] wrote {a.gated_dir} (partitioned) + {a.out_dir}/positive_fractions.csv")
    print("[step5] positive fraction per image (%):")
    print(frac_df.pivot(index="marker", columns="image", values="frac_pos")
          .mul(100).round(1).to_string())
    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
