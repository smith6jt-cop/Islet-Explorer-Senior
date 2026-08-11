"""Why the Otsu floor is being retired, on synthetic data with a known answer.

The point of synthetic data here is that the truth is known exactly, so the two
methods can be scored rather than argued about. Three prevalences spanning the
regimes the cohort actually contains: a sparse immune marker, a mid-abundance
structural marker, and an abundant epithelial one.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from phenocycler.marker_calibration import empirical_upper_threshold  # noqa: E402
from phenocycler.restore_validation import (  # noqa: E402
    raw_otsu_floor,
    raw_triangle_floor,
)

OUT = Path(__file__).resolve().parents[2] / "data/pilots/figures"
N = 200_000
CASES = [(0.003, "0.3%  (CD3e, MPO, CD8 regime)"),
         (0.05, "5%  (CD31 regime)"),
         (0.30, "30%  (E-cadherin regime)")]

rng = np.random.default_rng(0)
control = np.log1p(rng.gamma(2.0, 50.0, 30_000))

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
for ax, (fraction, label) in zip(axes, CASES):
    n_pos = int(N * fraction)
    values = np.concatenate([
        rng.gamma(2.0, 50.0, N - n_pos),
        rng.gamma(20.0, 200.0, n_pos),
    ])
    otsu = raw_otsu_floor(values)
    triangle = raw_triangle_floor(values)
    empirical = float(np.expm1(empirical_upper_threshold(control, 0.001)))
    clip = float(np.quantile(values, 0.995))

    bins = np.linspace(0, np.log1p(np.quantile(values, 0.9999)), 200)
    ax.hist(np.log1p(values[: N - n_pos]), bins=bins, color="#bbbbbb", label="true negative")
    ax.hist(np.log1p(values[N - n_pos:]), bins=bins, color="#4c72b0", alpha=.85,
            label=f"true positive (n={n_pos:,})")
    for value, colour, name in (
        (otsu, "#c44e52", "Otsu floor (shipped)"),
        (triangle, "#dd8452", "Triangle floor"),
        (empirical, "#55a868", "empirical upper tail (successor)"),
        (clip, "#000000", "p99.5 clip point"),
    ):
        ax.axvline(np.log1p(value), color=colour, lw=1.8,
                   ls=":" if name.startswith("p99.5") else "-", label=name)
    ax.set_yscale("log")
    ax.set_xlabel("log1p(intensity)")
    ratio_o = (values > otsu).sum() / n_pos
    ratio_e = (values > empirical).sum() / n_pos
    ax.set_title(f"{label}\nOtsu calls {ratio_o:.2f}x truth   "
                 f"successor calls {ratio_e:.2f}x truth", fontsize=10)

axes[0].legend(fontsize=8, loc="upper right")
fig.suptitle(
    "Retiring the Otsu input floor: one rule that holds across prevalence, versus one that does not\n"
    "synthetic data with a known positive count; the successor uses each pair's own control at the "
    "registry's familywise alpha",
    fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT / "fig5_floor_vs_successor.png", dpi=125)
print("wrote", OUT / "fig5_floor_vs_successor.png")
