"""Figure: the two ends of the intensity scale are undefined.

Left  -- exposure time per gate marker across every slide. Points that do not
         line up are the same antigen measured in different units.
Right -- fraction of analysis-eligible cells clamped to exactly 0 by REDSEA,
         per donor per marker. Zero is a clamp, and how much it hides varies
         by donor.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "data/runs/3bf4fc8e391d007766d1"
OUT = REPO / "data/pilots/figures"
DONORS = ["6436", "6450", "6566"]
GATES = ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM", "Keratin_5", "INS",
         "GCG", "CD31", "CD34", "Podoplanin", "CD3e", "CD68", "CD79a", "CD11c",
         "Iba1", "CD206", "MPO", "CD4", "CD8", "Vimentin", "SMA", "B3TUBB"]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(16.5, 7.6),
                               gridspec_kw={"width_ratios": [1.15, 1]})

# ---- A: exposure per marker across slides -------------------------------
acq = pd.read_parquet(REPO / "data/pilots/acquisition_by_biomarker.parquet")
acq = acq[acq.biomarker.isin(GATES)].dropna(subset=["exposure"])
order = (acq.groupby("biomarker").exposure.agg(lambda s: s.max() / s.min())
         .sort_values(ascending=False).index.tolist())
epoch_colour = {d: c for d, c in zip(sorted(acq.calib.dropna().unique()),
                                     ["#4c72b0", "#dd8452", "#55a868"])}
for i, m in enumerate(order):
    sub = acq[acq.biomarker == m]
    jitter = np.linspace(-0.22, 0.22, len(sub))
    axA.scatter(sub.exposure / 1000, i + jitter, s=34,
                c=[epoch_colour.get(c, "#999") for c in sub.calib], zorder=3)
    lo, hi = sub.exposure.min() / 1000, sub.exposure.max() / 1000
    if hi > lo:
        axA.plot([lo, hi], [i, i], color="#c44e52", lw=2.2, alpha=.55, zorder=2)
        axA.text(hi + 3, i, f"{hi/lo:.2f}x", va="center", fontsize=8, color="#c44e52")
axA.set_yticks(range(len(order)))
axA.set_yticklabels(order, fontsize=9)
axA.set_xlabel("exposure time (ms)")
axA.invert_yaxis()
axA.grid(axis="x", alpha=.25)
axA.set_title("Same marker, different exposure between donors\n"
              "16 of 22 gate markers vary; colour = responsivity calibration epoch",
              fontsize=11)
axA.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c, label=d)
                    for d, c in epoch_colour.items()],
           title="calibrated", fontsize=8, title_fontsize=8, loc="lower right")

# ---- B: cells clamped to exactly zero by REDSEA --------------------------
rows = {}
for d in DONORS:
    cols = [f"Cell__{m}" for m in GATES]
    have = set(pq.ParquetFile(next((RUN / f"redsea/donor_id={d}").glob("*.parquet"))).schema_arrow.names)
    cols = [c for c in cols if c in have]
    r = pq.read_table(next((RUN / f"redsea/donor_id={d}").glob("*.parquet")),
                      columns=["object_id"] + cols).to_pandas()
    g = pq.read_table(next((RUN / f"geometry_qc/donor_id={d}").glob("*.parquet")),
                      columns=["object_id", "analysis_eligible"]).to_pandas()
    t = r.merge(g, on="object_id")
    t = t[t.analysis_eligible]
    rows[d] = {c.replace("Cell__", ""): 100 * float((t[c] == 0).mean()) for c in cols}
z = pd.DataFrame(rows).reindex(order)

y = np.arange(len(z))
for i, d in enumerate(DONORS):
    axB.barh(y + (i - 1) * 0.27, z[d], height=0.26, label=f"donor {d}")
axB.set_yticks(y)
axB.set_yticklabels(order, fontsize=9)
axB.invert_yaxis()
axB.set_xlabel("% of analysis-eligible cells clamped to exactly 0 by REDSEA")
axB.grid(axis="x", alpha=.25)
axB.legend(fontsize=9)
axB.set_title("Zero is a clamp, not a measurement\n"
              "np.clip(corrected_sum, 0, None) — and how much it hides varies by donor",
              fontsize=11)

fig.suptitle("The intensity scale has no defined zero or maximum "
             "(no pixel anywhere is saturated: 0 of 63M at the 16-bit ceiling)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.945])
fig.savefig(OUT / "fig4_intensity_scale.png", dpi=125)
print("wrote", OUT / "fig4_intensity_scale.png")
