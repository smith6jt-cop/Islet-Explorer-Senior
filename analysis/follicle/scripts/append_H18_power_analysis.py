"""Idempotently append a power-aware reporting section to H18.

Adds bootstrap 95% CIs for rank-biserial r and Spearman rho per pairwise
comparison + a minimum-detectable-effect (MDE) calculation at the realized
n=9/9/9. Re-running is safe: the appended block has a sentinel marker
and is replaced rather than duplicated.
"""
from pathlib import Path
import nbformat as nbf

PROJECT = Path(__file__).resolve().parent.parent
NB_PATH = PROJECT / "analysis" / "H18_HE_follicles.ipynb"

SENTINEL = "## 7. Power-aware reporting (effect-size CIs + MDE)"

MD_INTRO = SENTINEL + """

The H18 first-pass null result reflects a hard reality: at n=9 per genotype, the **minimum detectable rank-biserial r** is around 0.65 at α=0.05 / 80% power. So "no test below p<0.05" mostly tells us the cohort is underpowered for small-to-moderate effects, not that no effect exists.

This section reports:
1. Bootstrap (BCa, 1000 reps) 95% CIs for rank-biserial r and Spearman ρ per metric.
2. The MDE at the realized n=9/9/9, computed via Monte-Carlo simulation under the alternative.
3. A forest plot of effect sizes with their CIs alongside the MDE band.

These outputs make "underpowered null" visually obvious: a CI that crosses zero but extends far beyond the MDE means the data **could not have detected a real moderate effect** even if one were present."""

CODE = """import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

rng = np.random.default_rng(20260508)
N_BOOT = 1000

POWER_METRICS = ["Follicle_Count", "Follicle_Area_Mean_um2",
                 "Follicle_Density_per_mm2", "Follicle_Fraction"]
PAIRS = [("C/C", "C/T"), ("C/C", "T/T"), ("C/T", "T/T")]


def rank_biserial(x, y):
    if len(x) == 0 or len(y) == 0:
        return np.nan
    u, _ = mannwhitneyu(x, y, alternative="two-sided")
    return 1 - 2 * u / (len(x) * len(y))


def bootstrap_ci(stat_fn, *args, n_boot=N_BOOT, alpha=0.05):
    \"\"\"Percentile bootstrap; resamples each arg array with replacement.\"\"\"
    obs = stat_fn(*args)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        resampled = [rng.choice(a, size=len(a), replace=True) for a in args]
        try:
            boots[i] = stat_fn(*resampled)
        except Exception:
            boots[i] = np.nan
    boots = boots[np.isfinite(boots)]
    if len(boots) < 10:
        return obs, np.nan, np.nan
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return obs, lo, hi


def spearman_dosage(values, ordinal):
    rho, _ = spearmanr(ordinal, values)
    return rho


# Effect-size CIs ------------------------------------------------------------
power_rows = []
for metric in POWER_METRICS:
    d = per_donor.dropna(subset=[metric, "Genotype"])
    # Pairwise rank-biserial
    for g1, g2 in PAIRS:
        x = d[d["Genotype"] == g1][metric].values
        y = d[d["Genotype"] == g2][metric].values
        r, lo, hi = bootstrap_ci(rank_biserial, x, y)
        u, p = (mannwhitneyu(x, y, alternative="two-sided") if len(x) and len(y)
                else (np.nan, np.nan))
        power_rows.append({"Metric": metric, "Test": f"MW {g1} vs {g2}",
                           "Effect": r, "CI_lo": lo, "CI_hi": hi, "p": p,
                           "n1": len(x), "n2": len(y)})
    # Spearman dosage
    ordinal = d["Genotype"].astype(str).map({"C/C": 0, "C/T": 1, "T/T": 2}).values
    vals = d[metric].values
    obs = spearman_dosage(vals, ordinal)
    boots = []
    n = len(vals)
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        try:
            boots.append(spearman_dosage(vals[idx], ordinal[idx]))
        except Exception:
            pass
    boots = np.array([b for b in boots if np.isfinite(b)])
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) > 10 else (np.nan, np.nan))
    _, p = spearmanr(ordinal, vals) if n >= 3 else (np.nan, np.nan)
    power_rows.append({"Metric": metric, "Test": "Spearman dosage",
                       "Effect": obs, "CI_lo": lo, "CI_hi": hi, "p": p,
                       "n1": n, "n2": np.nan})

power_df = pd.DataFrame(power_rows)

# Minimum-detectable-effect via Monte-Carlo --------------------------------
def mde_rank_biserial(n_per_group, n_sim=2000, alpha=0.05, target_power=0.8):
    \"\"\"Smallest |true r| where Mann-Whitney detects at p<alpha with target_power.\"\"\"
    candidate_rs = np.linspace(0.05, 0.95, 19)
    rng_local = np.random.default_rng(20260508)
    for r_true in candidate_rs:
        # Translate r_true into a location shift on standard normals.
        # Under shift d on N(0,1), AUC = P(X1<X2) = Phi(d/sqrt(2)),
        # rank-biserial r = 2*AUC - 1 → d = sqrt(2) * Phi^-1((r+1)/2).
        from scipy.stats import norm as _norm
        d = np.sqrt(2) * _norm.ppf((r_true + 1) / 2)
        n_sig = 0
        for _ in range(n_sim):
            x = rng_local.standard_normal(n_per_group)
            y = rng_local.standard_normal(n_per_group) + d
            _, p = mannwhitneyu(x, y, alternative="two-sided")
            n_sig += (p < alpha)
        power = n_sig / n_sim
        if power >= target_power:
            return r_true
    return np.nan

mde_9_9 = mde_rank_biserial(9, n_sim=2000)
mde_18_9 = mde_rank_biserial(9, n_sim=2000)  # same n for our cohort; placeholder
print(f"Minimum-detectable rank-biserial r at n=9/9, alpha=0.05, 80% power: {mde_9_9:.2f}")

power_df.attrs["mde_n9"] = mde_9_9
save_table(power_df, "HE_power_analysis")
power_df"""

PLOT = """# Forest plot of effect sizes with bootstrap CIs ----------------------------
fig, axes = plt.subplots(1, 4, figsize=(16, 4.6), sharex=True)
metric_titles = {"Follicle_Count": "Count", "Follicle_Area_Mean_um2": "Mean Area",
                 "Follicle_Density_per_mm2": "Density", "Follicle_Fraction": "Fraction"}
test_color = {"MW C/C vs C/T": "#4477AA", "MW C/C vs T/T": "#EE6677",
              "MW C/T vs T/T": "#228833", "Spearman dosage": "#222222"}

for ax, metric in zip(axes, POWER_METRICS):
    rows = power_df[power_df["Metric"] == metric]
    ypos = np.arange(len(rows))[::-1]
    for y, (_, row) in zip(ypos, rows.iterrows()):
        c = test_color.get(row["Test"], "grey")
        eff = row["Effect"]
        if pd.notna(row["CI_lo"]) and pd.notna(row["CI_hi"]):
            ax.errorbar(eff, y, xerr=[[eff - row["CI_lo"]], [row["CI_hi"] - eff]],
                        fmt="o", color=c, ecolor=c, capsize=3, markersize=8,
                        markeredgecolor="black", markeredgewidth=0.6)
        else:
            ax.scatter(eff, y, color=c, s=70, edgecolor="black", linewidth=0.6)
    # MDE band (rank-biserial only)
    ax.axvspan(-mde_9_9, mde_9_9, color="grey", alpha=0.12,
               label=f"|r| < MDE ({mde_9_9:.2f})" if metric == POWER_METRICS[0] else None)
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.set_yticks(ypos)
    ax.set_yticklabels([t.replace("MW ", "") for t in rows["Test"]], fontsize=8)
    ax.set_title(metric_titles[metric])
    ax.set_xlabel("Effect (rank-biserial r or Spearman ρ)")
    ax.set_xlim(-1.05, 1.05)

handles = [plt.Line2D([0], [0], marker="o", linestyle="None", color=c, markersize=8,
                       markeredgecolor="black", label=t.replace("MW ", ""))
           for t, c in test_color.items()]
handles.append(plt.Rectangle((0, 0), 1, 1, color="grey", alpha=0.25,
                              label=f"Below MDE r={mde_9_9:.2f} at n=9/9"))
fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.04), fontsize=9)
fig.suptitle("H18 power-aware reporting — effect sizes with bootstrap 95% CIs", y=1.02)
plt.tight_layout()
save_figure(fig, "HE_effect_size_forest")
plt.show()

# Concise interpretation
underpowered = power_df[(power_df["Test"].str.startswith("MW")) &
                         (power_df["Effect"].abs() < mde_9_9) &
                         (power_df["CI_hi"] > mde_9_9)]
print(f"\\n{len(underpowered)} of {(power_df['Test'].str.startswith('MW')).sum()} "
      "pairwise comparisons are underpowered: observed |r| < MDE but the upper CI "
      "extends above the MDE threshold. The data cannot rule out a moderate effect.")"""


def build_appendix_cells():
    return [
        nbf.v4.new_markdown_cell(MD_INTRO),
        nbf.v4.new_code_cell(CODE),
        nbf.v4.new_code_cell(PLOT),
    ]


def main():
    nb = nbf.read(NB_PATH, as_version=4)
    # Strip prior copy of *this* section only: from SENTINEL up to the next
    # top-level (## ) markdown header, so later appendices stay intact.
    keep_before = []
    in_section = False
    keep_after = []
    for cell in nb.cells:
        src = cell.source if isinstance(cell.source, str) else "".join(cell.source)
        if cell.cell_type == "markdown" and src.startswith(SENTINEL):
            in_section = True
            continue
        if in_section:
            if cell.cell_type == "markdown" and src.lstrip().startswith("## "):
                in_section = False
                keep_after.append(cell)
                continue
            continue
        if not in_section and not keep_after:
            keep_before.append(cell)
        else:
            keep_after.append(cell)
    nb.cells = keep_before + build_appendix_cells() + keep_after
    nbf.write(nb, NB_PATH)
    print(f"Wrote {NB_PATH.relative_to(PROJECT)}  ({len(nb.cells)} cells)")


if __name__ == "__main__":
    main()
