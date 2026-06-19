#!/usr/bin/env python3
"""
RESTORE per-image intensity normalization of mutually-exclusive markers.

Applies RESTORE (Chang Lab / OHSU; github.com/smith6jt-cop/RESTORE, vendored at
external/RESTORE @ commit 38df59b) to the 59-plex protein single-cell data, for the
mutually-exclusive markers spanning the broad pancreas cell types:
    immune:      CD3e (T), CD20 (B), CD163 (macrophage)  -- background from each other
    exocrine:    Pan_Cytokeratin (acinar/ductal)         -- background from Vimentin
    mesenchyme:  Vimentin (stroma/fibroblast)            -- background from Pan_Cytokeratin
    endocrine:   INS (beta), GCG (alpha)                 -- background from acinar (Pan_CK)
    endothelium: CD31                                    -- background from acinar (Pan_CK)
    muscle:      SMA (smooth muscle/pericyte)            -- background from acinar (Pan_CK)
(CD99 was dropped -- broadly expressed, not lineage-exclusive in this panel. See
DEFAULT_MARKER_PAIRS; override with --marker-pairs.)

RESTORE pairs each target marker with a mutually-exclusive *reference* marker. Cells
positive for the reference are, by biology, negative for the target, so the target
intensity within the reference-positive population is pure background/autofluorescence.
Per image, three 2-cluster models (KMeans, GMM, SSC) split background vs signal and the
threshold is set to  mean + 3*sigma  of the target intensity in the target-negative
cluster. This corrects per-image autofluorescence drift that biases fixed-cutoff gating.

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
# The reference's positive cells define the target's background. References are chosen
# to be abundant, cleanly target-negative populations (T cells / macrophages).
DEFAULT_MARKER_PAIRS = [
    # --- immune lineages (background from other immune cells) ---
    ["CD3e",            "CD163"],            # T cells              <- macrophages
    ["CD20",            "CD3e"],             # B cells              <- T cells
    ["CD163",           "CD3e"],             # macrophages         <- T cells
    # --- broad pancreas compartments. Pan_Cytokeratin marks ~91% of cells (acinar/
    #     ductal exocrine) and is cleanly negative in every non-epithelial population,
    #     so it is the abundant universal reference; epithelium uses the reciprocal
    #     mesenchymal Vimentin (epithelial/stromal corr ~ -0.5, the broadest exclusivity)
    ["Pan_Cytokeratin", "Vimentin"],         # exocrine epithelium <- mesenchyme
    ["Vimentin",        "Pan_Cytokeratin"],  # mesenchyme/stroma   <- epithelium
    ["INS",             "Pan_Cytokeratin"],  # endocrine beta      <- acinar
    ["GCG",             "Pan_Cytokeratin"],  # endocrine alpha     <- acinar
    ["CD31",            "Pan_Cytokeratin"],  # endothelium         <- acinar
    ["SMA",             "Pan_Cytokeratin"],  # smooth muscle/pericyte <- acinar
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
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        df = df.astype("float64").reset_index(drop=True)
        df["scene"] = donor
        df["batch"] = 0
        frames.append(df)
        print(f"[input] image {donor}: {n:,} cells -> {len(df):,} sampled")
    data = pd.concat(frames, ignore_index=True)
    print(f"[input] total RESTORE input: {len(data):,} cells / {data['scene'].nunique()} images")
    return data


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
                    help="cells per image for threshold estimation (SSC-tractable)")
    ap.add_argument("--model", default="SSC", choices=["GMM", "KMeans", "SSC"],
                    help="model whose threshold drives gating + normalization. "
                         "SSC is the default: on this panel GMM/KMeans place the "
                         "abundant lineage markers (CD3e/CD20/CD163) at/above the "
                         "99.9th intensity pct (near-zero positives); SSC lands at the "
                         "background->signal shoulder. All 3 are saved in the CSV.")
    ap.add_argument("--marker-pairs", default=None,
                    help="override pairs, e.g. 'CD3e:CD163,CD99:CD163'")
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
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    pairs = (DEFAULT_MARKER_PAIRS if not a.marker_pairs
             else [p.split(":") for p in a.marker_pairs.split(",")])
    markers = sorted({m for pair in pairs for m in pair})
    print(f"[config] pairs (target<-reference): {pairs}")
    print(f"[config] markers={markers}  chosen_model={a.model}  subsample={a.subsample}")

    Normalization, _ = import_restore()
    qc_dir = a.out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    print("[step1] building RESTORE input ...")
    data = build_restore_input(a.cells_dir, markers, a.subsample, a.seed, a.limit_scenes, a.donors)

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
