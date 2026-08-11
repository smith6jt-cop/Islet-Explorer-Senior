"""Review figures for the 3-donor broad-phenotyping pilot.

Windowing rule: every intensity panel is drawn on the workflow's own measured
quantities -- the marker's input floor and its candidate divisors. Nothing is
clipped at a percentile black point; where the axis ends before the data does,
the panel states what fraction of cells lie beyond it.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
PILOT = REPO / "data" / "pilots" / "baseline_p999"
SURVEY = REPO / "data" / "pilots" / "divisor_survey"
OUT = REPO / "data" / "pilots" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
DONORS = ["6436", "6450", "6566"]

RULES = {
    "max": ("#d62728", "manuscript max"),
    "p999": ("#1f77b4", "p999 (production)"),
    "p99": ("#ff7f0e", "p99"),
    "p95": ("#9467bd", "p95"),
    "mean3sd": ("#8c564b", "mean+3SD"),
}
LABEL_COLORS = {
    "Epithelial": "#4c72b0", "Immune": "#c44e52", "Endothelial": "#55a868",
    "Mesenchymal": "#8172b2", "Neural": "#ccb974", "Other": "#999999",
    "Ambiguous": "#dd8452", "Unavailable": "#cccccc",
}


def fig_divisor_landscape(donor: str) -> None:
    survey = pd.read_parquet(SURVEY / "divisor_survey.parquet")
    survey = survey[survey.donor == donor].sort_values("target")
    n = len(survey)
    ncol, nrow = 4, int(np.ceil(n / 4))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 2.9 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, r) in zip(axes, survey.iterrows()):
        path = SURVEY / "dist" / f"{donor}__{r.target}.npz"
        if not path.exists():
            ax.axis("off")
            continue
        z = np.load(path)
        pop, ctrl = z["population"].astype(float), z["control"].astype(float)
        cand = {k: float(z[f"div_{k}"]) for k in RULES}
        right = max(max(cand.values()), np.quantile(pop, 0.999))
        beyond = float((pop > right).mean())
        bins = np.linspace(0, np.log1p(right), 120)

        ax.hist(np.log1p(np.clip(pop, 0, right)), bins=bins, color="#bbbbbb",
                label=f"all callable cells (n={len(pop):,})")
        ax.hist(np.log1p(np.clip(ctrl, 0, right)), bins=bins, color="#2b6cb0",
                alpha=0.65, label=f"reference-negative control (n={len(ctrl):,})")
        ax.axvline(np.log1p(r.target_input_floor), color="k", ls=":", lw=1.4)
        for k, (colour, _) in RULES.items():
            ax.axvline(np.log1p(cand[k]), color=colour, lw=1.3)
        ax.set_yscale("log")
        pos = 100 * float(r.get("pos_p999", np.nan))
        ax.set_title(f"{r.target} <- {r.reference}\n"
                     f"p999 calls {pos:.1f}%   anchor {int(r.target_anchor_n):,}"
                     + (f"   [{beyond*100:.2f}% beyond axis]" if beyond > 0.0005 else ""),
                     fontsize=8.5)
        ax.set_xlabel("log1p(REDSEA-corrected intensity)", fontsize=7)
        ax.tick_params(labelsize=7)

    for ax in axes[n:]:
        ax.axis("off")
    handles = [Patch(color="#bbbbbb", label="all callable cells"),
               Patch(color="#2b6cb0", alpha=.65, label="reference-negative control"),
               plt.Line2D([], [], color="k", ls=":", label="input floor (raw scale)")]
    handles += [plt.Line2D([], [], color=c, label=lab) for c, lab in RULES.values()]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.suptitle(
        f"Donor {donor}: where every candidate divisor lands on each marker's own control "
        f"distribution\nlog-y so the sparse positive tail stays visible; no percentile black point",
        fontsize=12)
    fig.tight_layout(rect=[0, 0.045, 1, 0.965])
    fig.savefig(OUT / f"fig1_divisor_landscape_{donor}.png", dpi=115)
    plt.close(fig)


def fig_rule_vs_anchor() -> None:
    s = pd.read_parquet(SURVEY / "divisor_survey.parquet")
    fig, ax = plt.subplots(figsize=(13, 6))
    targets = sorted(s.target.unique())
    xpos = {t: i for i, t in enumerate(targets)}
    for k, (colour, label) in RULES.items():
        x = [xpos[t] for t in s.target]
        y = s["anchor_ratio_" + k].replace([np.inf, -np.inf], np.nan)
        ax.scatter(np.array(x) + (list(RULES).index(k) - 2) * 0.13, y,
                   s=22, color=colour, label=label, alpha=0.85)
    ax.axhline(1, color="k", lw=1)
    ax.axhspan(0.5, 2, color="k", alpha=0.07)
    ax.set_yscale("log")
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, rotation=90, fontsize=8)
    ax.set_ylabel("cells called positive / that pair's own NNMF anchor")
    ax.set_title("No single divisor rule fits every marker and donor\n"
                 "each point is one donor-marker pair (3 donors); grey band = within 2x of the anchor")
    ax.legend(fontsize=9, ncol=5)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_rule_vs_anchor.png", dpi=125)
    plt.close(fig)


def fig_composition() -> None:
    order = ["Epithelial", "Immune", "Endothelial", "Mesenchymal", "Neural",
             "Other", "Ambiguous", "Unavailable"]
    frames = {}
    for d in DONORS:
        a = pd.read_parquet(PILOT / f"assignments/donor_id={d}.parquet",
                            columns=["broad_label"])
        frames[d] = a.broad_label.value_counts(normalize=True) * 100
    comp = pd.DataFrame(frames).reindex(order).fillna(0)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 5.2),
                                  gridspec_kw={"width_ratios": [1, 1.25]})
    bottom = np.zeros(len(DONORS))
    for lab in order:
        v = comp.loc[lab].to_numpy()
        ax.bar(DONORS, v, bottom=bottom, color=LABEL_COLORS[lab], label=lab)
        for i, (b, h) in enumerate(zip(bottom, v)):
            if h > 3:
                ax.text(i, b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=8)
        bottom += v
    ax.set_ylabel("% of all cells")
    ax.set_title("Broad composition, production p999\nAmbiguous / Other / Unavailable shown, never folded away")
    ax.legend(fontsize=8, bbox_to_anchor=(1.0, 1.0), loc="upper left")

    rows = []
    for d in DONORS:
        a = pd.read_parquet(PILOT / f"assignments/donor_id={d}.parquet",
                            columns=["broad_label", "broad_supported_gates"])
        amb = a[a.broad_label == "Ambiguous"]
        share = amb.broad_supported_gates.value_counts(normalize=True) * 100
        rows.append(share.rename(d))
    att = pd.DataFrame(rows).T.fillna(0)
    att = att.loc[att.mean(axis=1).sort_values(ascending=False).index][:8][::-1]
    y = np.arange(len(att))
    for i, d in enumerate(DONORS):
        ax2.barh(y + (i - 1) * 0.27, att[d], height=0.26, label=f"donor {d}")
    ax2.set_yticks(y)
    ax2.set_yticklabels(att.index, fontsize=8)
    ax2.set_xlabel("% of that donor's Ambiguous cells")
    ax2.set_title("Which gates collide\nMesenchymal (Vimentin) drives most of it")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_composition_and_ambiguity.png", dpi=125)
    plt.close(fig)


if __name__ == "__main__":
    for donor in DONORS:
        fig_divisor_landscape(donor)
        print("wrote landscape", donor, flush=True)
    fig_rule_vs_anchor()
    fig_composition()
    print("figures ->", OUT)
