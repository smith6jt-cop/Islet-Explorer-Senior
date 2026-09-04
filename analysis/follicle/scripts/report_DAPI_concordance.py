"""Summarise how well the DAPI U-Net reproduces the QuPath pixel classifier.

The training labels came from `SpleenRegions*`, which used 25-33 marker
channels. This model sees DAPI only, so it is distilling a multi-channel teacher
into a single-channel student and some error is irreducible. This script states
that plainly instead of reporting a single headline number:

  * per-donor Follicle count and mean area, predicted vs teacher (Spearman rho —
    the same success bar the H&E pipeline uses, rho > 0.9)
  * per-class region IoU
  * the Follicle <-> PALS confusion rate

The last one is the point. DAPI cannot see CD20 vs CD3e, and both zones are
dense-nuclei white pulp, so a model can post a respectable follicle count while
systematically swapping the two. A report that led with mean IoU would hide
exactly the failure this cohort is most prone to.

Reads the artifacts `train_DAPI_unet.py` already writes — it does not re-run
inference.

  python scripts/report_DAPI_concordance.py
  python scripts/report_DAPI_concordance.py --models-dir analysis/DAPI_unet_models
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import DAPI_unet_utils as D
import HE_unet_utils as U


def load_cv(models_dir: Path) -> list[dict]:
    p = models_dir / "cv_summary.json"
    if not p.exists():
        raise SystemExit(f"No cv_summary.json in {models_dir}. Run "
                         "train_DAPI_unet.py first.")
    return json.loads(p.read_text())


def per_donor_frame(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        inst = r["best"].get("instance_metrics", {})
        regs = r["best"].get("region_metrics", {})
        for name, m in inst.items():
            row = {
                "Image": name,
                # New runs carry the declared ID. The exact manifest fallback
                # keeps older real-data summaries readable without parsing.
                "Donor": m.get("sample_id") or U.donor_id_from_name(name),
                "Fold": r["fold"],
                "Teacher follicle count": m["gt_count"],
                "U-Net follicle count": m["pred_count"],
                "Teacher mean area µm^2": round(m["gt_mean_area_um2"], 2),
                "U-Net mean area µm^2": round(m["pred_mean_area_um2"], 2),
                "Follicle IoU": round(r["best"]["per_donor"].get(name, float("nan")), 4),
            }
            reg = regs.get(name)
            if reg:
                for cls, v in reg["per_class_iou"].items():
                    row[f"IoU {cls}"] = None if v is None else round(v, 4)
                row["Follicle->PALS"] = round(reg.get("follicle_as_pals", float("nan")), 4)
                row["PALS->Follicle"] = round(reg.get("pals_as_follicle", float("nan")), 4)
            rows.append(row)
    return pd.DataFrame(rows)


def spearman(a, b) -> float:
    if len(a) < 3:
        return float("nan")
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).correlation)


def make_figure(df: pd.DataFrame, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    iou_cols = [c for c in df.columns if c.startswith("IoU ")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    for ax, (gt, pred, label) in zip(axes[:2], [
        ("Teacher follicle count", "U-Net follicle count", "Follicle count"),
        ("Teacher mean area µm^2", "U-Net mean area µm^2", "Mean follicle area (µm²)"),
    ]):
        x, y = df[gt].to_numpy(float), df[pred].to_numpy(float)
        ax.scatter(x, y, s=48, alpha=0.8, edgecolor="black", linewidth=0.5)
        lim = [0, max(1.0, np.nanmax(np.concatenate([x, y])) * 1.1)]
        ax.plot(lim, lim, "k--", lw=1, alpha=0.5, label="identity")
        ax.set(xlim=lim, ylim=lim, xlabel=f"Pixel classifier (teacher)",
               ylabel="DAPI U-Net", title=f"{label}\nSpearman ρ={spearman(x, y):.3f}")
        ax.legend(fontsize=8)

    ax = axes[2]
    if iou_cols:
        means, labels = [], []
        for c in iou_cols:
            vals = pd.to_numeric(df[c], errors="coerce").dropna()
            if len(vals):
                means.append(vals.mean())
                labels.append(c.replace("IoU ", ""))
        ax.barh(labels, means, color="#4C72B0")
        ax.set(xlim=(0, 1), xlabel="IoU", title="Per-class region IoU")
        for i, v in enumerate(means):
            ax.text(min(v + 0.02, 0.95), i, f"{v:.2f}", va="center", fontsize=8)
    else:
        ax.text(0.5, 0.5, "no region metrics\n(trained with --no-aux)",
                ha="center", va="center")
        ax.axis("off")

    fig.suptitle("DAPI U-Net vs SpleenRegions pixel classifier "
                 "(DAPI-only student, multi-channel teacher)", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models-dir", default=str(D.MODELS_DIR))
    p.add_argument("--out-csv", default=str(U.OUTPUT_ROOT / "analysis" / "tables"
                                            / "DAPI_unet_concordance.csv"))
    p.add_argument("--out-png", default=str(U.OUTPUT_ROOT / "analysis" / "figures"
                                            / "DAPI_unet_concordance.png"))
    p.add_argument("--no-figure", action="store_true")
    args = p.parse_args()

    results = load_cv(Path(args.models_dir))
    df = per_donor_frame(results)
    if df.empty:
        raise SystemExit("cv_summary.json has no instance_metrics — the run "
                         "produced no best checkpoint.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    rho_c = spearman(df["Teacher follicle count"], df["U-Net follicle count"])
    rho_a = spearman(df["Teacher mean area µm^2"], df["U-Net mean area µm^2"])

    print(f"Donors: {len(df)}   (Spearman needs >=3)")
    print(f"Follicle count   Spearman rho = {rho_c:.3f}"
          f"   {'PASS' if rho_c > 0.9 else 'below the 0.9 bar'}")
    print(f"Mean area        Spearman rho = {rho_a:.3f}"
          f"   {'PASS' if rho_a > 0.9 else 'below the 0.9 bar'}")
    print(f"Mean Follicle IoU            = {df['Follicle IoU'].mean():.3f}")

    iou_cols = [c for c in df.columns if c.startswith("IoU ")]
    if iou_cols:
        print("\nPer-class region IoU (teacher as reference):")
        for c in iou_cols:
            vals = pd.to_numeric(df[c], errors="coerce").dropna()
            label = c.replace("IoU ", "")
            print(f"  {label:<14}"
                  + (f"{vals.mean():.3f}  (n={len(vals)} donors)" if len(vals)
                     else "absent from every donor"))
    if "Follicle->PALS" in df.columns:
        fap = pd.to_numeric(df["Follicle->PALS"], errors="coerce").mean()
        paf = pd.to_numeric(df["PALS->Follicle"], errors="coerce").mean()
        print(f"\nWhite-pulp confusion:  Follicle->PALS {fap:.3f}   "
              f"PALS->Follicle {paf:.3f}")
        print("  This is the DAPI-only ceiling, not a bug: the teacher used CD20 "
              "vs CD3e to\n  separate these zones and both are dense-nuclei white "
              "pulp in DAPI alone.\n  Re-run with --merge-white-pulp to see how "
              "much of the region error it accounts for.")

    print(f"\nWrote {out_csv}")
    if not args.no_figure:
        try:
            make_figure(df, Path(args.out_png))
            print(f"Wrote {args.out_png}")
        except ImportError as e:
            # The CSV and the printed summary are the deliverable; matplotlib is
            # only needed for the figure. Don't fail the whole report over it.
            print(f"Skipped figure ({e}). "
                  "pip install matplotlib, or pass --no-figure.")


if __name__ == "__main__":
    main()
