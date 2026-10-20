"""Per-donor CD8 ROC, for the maintainer's CD8 ruling (asked 2026-10-02; 2026-10-18: "The true CD8
determination requires checking an ROC plot for each case").

The question is one rule for all 25 donors: keep UniFORM's CD8 shift (the method of record) or align CD8
on no donor. UniFORM moves a donor by ONE multiplier, so a donor's ROC curve is the same under both;
what differs is where the cohort's positivity level (0.10 on the cohort scale) lands on that curve.

Reference cells are defined without CD8, from the same aligned values and the same 0.10 level:
CD8-type T cells are CD3e+ CD4-; CD8-negative references are CD4+ T cells and hormone-positive cells
without CD3e. CD3e+ CD4- also holds double-negative T cells, so the curve is a lower bound on what CD8
can do. The figure is decision support; nothing is decided here.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .common import Context, binding, require, threaded, write_json

GEN = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008/own_compartment")
LEVEL = 0.10
N_PER_DONOR = 200_000
HORMONES = ("INS", "GCG", "SST")


def donors() -> list[str]:
    return sorted(p.name.split("=", 1)[1] for p in (GEN / "uniform").glob("donor_id=*"))


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "roc_stage", "Unknown CD8 operation")
    return threaded(binding, [GEN / "shifts.json", GEN / "scale.parquet",
                              *(GEN / "uniform" / f"donor_id={d}" / "data.parquet" for d in donors())])


def scaled(values: np.ndarray, floor: float, span: float) -> np.ndarray:
    return np.clip((values - floor) / span, 0, 1)


def frame(donor: str, scale: pd.DataFrame, multiplier: float) -> pd.DataFrame:
    columns = ["object_id", "qc_analysis_eligible", "in_alignment_block", "CD8", "CD3e", "CD4", *HORMONES]
    t = pq.read_table(GEN / "uniform" / f"donor_id={donor}" / "data.parquet", columns=columns).to_pandas()
    t = t.loc[t.qc_analysis_eligible.fillna(False) & t.in_alignment_block.fillna(False)]
    s = {m: scaled(t[m].to_numpy(float), scale.at[m, "floor"], scale.at[m, "span"]) for m in ("CD3e", "CD4", *HORMONES)}
    hormone = np.nanmax(np.column_stack([s[m] for m in HORMONES]), axis=1)
    out = pd.DataFrame({"cd8_raw": t.CD8.to_numpy(float) / multiplier,
                        "positive": (s["CD3e"] >= LEVEL) & (s["CD4"] < LEVEL),
                        "negative_CD4_T": (s["CD3e"] >= LEVEL) & (s["CD4"] >= LEVEL),
                        "negative_hormone": (hormone >= LEVEL) & (s["CD3e"] < LEVEL)})
    return out.loc[np.isfinite(out.cd8_raw)]


def unshifted_threshold(frames: dict) -> tuple[float, float, float]:
    """The cohort level if CD8 were aligned on no donor: run_step3.build_scale's floor and top on raw CD8."""
    rng = np.random.default_rng(0)
    pool = []
    for d in sorted(frames):
        v = frames[d].cd8_raw.to_numpy(float)
        pool.append(v if len(v) <= N_PER_DONOR else v[rng.choice(len(v), N_PER_DONOR, replace=False)])
    v = np.concatenate(pool)
    top = float(np.percentile(v[np.isfinite(v)], 99.9))
    v = v[np.isfinite(v) & (v > 0)]
    h, e = np.histogram(np.log1p(v), bins=512)
    floor = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
    return floor + LEVEL * (top - floor), floor, top


def roc(scores_pos: np.ndarray, scores_neg: np.ndarray, points: int = 200) -> tuple[pd.DataFrame, float]:
    pooled = np.concatenate([scores_pos, scores_neg])
    thresholds = np.unique(np.quantile(pooled, np.linspace(0, 1, points)))
    curve = pd.DataFrame({"threshold": thresholds,
                          "TPR": [(scores_pos >= t).mean() for t in thresholds],
                          "FPR": [(scores_neg >= t).mean() for t in thresholds]})
    ranks = pd.Series(pooled).rank().to_numpy()
    auc = (ranks[:len(scores_pos)].sum() - len(scores_pos) * (len(scores_pos) + 1) / 2) / (len(scores_pos) * len(scores_neg))
    return curve, float(auc)


def roc_stage(context: Context, stage: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    shifts = json.loads((GEN / "shifts.json").read_text())["shifts"]["CD8"]
    scale = pd.read_parquet(GEN / "scale.parquet").set_index("marker")
    multipliers = {d: float(shifts["scale"][d]) for d in donors()}
    frames = {d: frame(d, scale, multipliers[d]) for d in donors()}
    t_b, floor_b, top_b = unshifted_threshold(frames)
    floor_a, span_a = float(scale.at["CD8", "floor"]), float(scale.at["CD8", "span"])
    rows, curves = [], []
    fig, axes = plt.subplots(5, 5, figsize=(13.33, 13.33), sharex=True, sharey=True, constrained_layout=True)
    for ax, d in zip(axes.flat, sorted(frames)):
        f = frames[d]
        t_a = (floor_a + LEVEL * span_a) / multipliers[d]
        pos = f.loc[f.positive, "cd8_raw"].to_numpy()
        row = {"donor_id": d, "CD8_multiplier": multipliers[d], "n_reference_CD8_type_T": len(pos),
               "threshold_shift_kept_raw": t_a, "threshold_no_shift_raw": t_b}
        for name in ("negative_CD4_T", "negative_hormone"):
            neg = f.loc[f[name], "cd8_raw"].to_numpy()
            row[f"n_{name}"] = len(neg)
            if len(pos) and len(neg):
                curve, auc = roc(pos, neg)
                curves.append(curve.assign(donor_id=d, negatives=name))
                row[f"AUC_vs_{name}"] = auc
                for tag, t in (("shift_kept", t_a), ("no_shift", t_b)):
                    row[f"TPR_{tag}"] = float((pos >= t).mean())
                    row[f"FPR_{tag}_vs_{name}"] = float((neg >= t).mean())
                ax.plot(curve.FPR, curve.TPR, color="#3a6e94" if name == "negative_CD4_T" else "#bd8a29", linewidth=1.2)
                ax.scatter([row[f"FPR_shift_kept_vs_{name}"]], [row["TPR_shift_kept"]], marker="o", s=28, color="#a74151", zorder=3)
                ax.scatter([row[f"FPR_no_shift_vs_{name}"]], [row["TPR_no_shift"]], marker="s", s=28, facecolor="none",
                           edgecolor="#222222", zorder=3)
        rows.append(row)
        ax.plot([0, 1], [0, 1], color="#cccccc", linewidth=0.8)
        ax.set_title(f"{d}  ×{multipliers[d]:.2f}", fontsize=10)
    for ax in axes[-1, :]:
        ax.set_xlabel("FPR", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("TPR", fontsize=9)
    fig.suptitle("CD8 ROC per donor: CD3e+CD4− T cells against CD4+ T cells (blue) and hormone+ CD3e− cells (gold)\n"
                 "red dot = cohort level with UniFORM's CD8 shift (as run); open square = cohort level with no CD8 shift",
                 fontsize=12)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stage / f"cd8_roc_per_donor.{ext}", dpi=170 if ext == "png" else None)
    plt.close(fig)
    table = pd.DataFrame(rows)
    table.to_csv(stage / "cd8_roc_summary.csv", index=False)
    pd.concat(curves, ignore_index=True).to_csv(stage / "cd8_roc_curves.csv", index=False)
    write_json(stage / "cd8_roc_definition.json", {
        "positives": "CD3e >= 0.10 and CD4 < 0.10 (cohort scale, aligned values); includes double-negative T cells",
        "negatives": {"negative_CD4_T": "CD3e >= 0.10 and CD4 >= 0.10",
                      "negative_hormone": "max(INS, GCG, SST) >= 0.10 and CD3e < 0.10"},
        "population": "QC-eligible cells in each donor's alignment block",
        "score": "the donor's unshifted CD8 (aligned CD8 / UniFORM multiplier)",
        "shift_kept": {"level": LEVEL, "floor": floor_a, "span": span_a, "raw_threshold": "(floor + 0.10*span) / multiplier"},
        "no_shift": {"level": LEVEL, "floor": floor_b, "top": top_b, "raw_threshold": t_b,
                     "method": "run_step3.build_scale floor/top on unshifted CD8 pooled over 25 donors (200k cells each, seed 0)"}})
    return {"donors": len(rows), "no_shift_raw_threshold": t_b,
            "median_AUC_vs_CD4_T": float(table["AUC_vs_negative_CD4_T"].median())}
