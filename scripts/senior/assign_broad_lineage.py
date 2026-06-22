#!/usr/bin/env python3
"""
Step 1 — broad lineage assignment (every cell typed, ZERO "Unassigned").

On the REDSEA-corrected, re-RESTORE'd data, each cell is typed by a HIERARCHICAL rule (every
cell typed — no "Unassigned" bucket; the user mandate is to fix the upstream cause, not invent
a junk category):

  1. Endocrine    if INS or GCG `_pos`              (islet hormones — specific, sparse)
  2. Immune       elif CD3e/CD20/CD163 `_pos`        (leukocytes — specific, sparse)
  3. Endothelial  elif CD31 `_pos`                   (vessels — specific)
  4. else STRUCTURAL background → argmax of (Pan_Cytokeratin, Vimentin, SMA) `_norm`
        → Epithelial / Fibroblast / Muscle
     ...with cells below ALL three structural thresholds defaulted to Epithelial (see below).

Note on the "Fibroblast" label: the six broad classes must be MUTUALLY EXCLUSIVE siblings (this
is what lets "Immune", a broad functional system, still be a valid sibling — no immune cell is
also epithelial/endothelial/etc.). Two rejected names for this compartment both fail that test:
  • "Mesenchymal" — names a developmental ORIGIN/state (the M in EMT), not a resident compartment,
    and its leaves (fibroblasts, stellate cells) are not called "mesenchymal cells".
  • "Stromal" — an UMBRELLA (stroma = the parenchyma's complement), so it OVERLAPS its own
    siblings: endothelium, mural/smooth-muscle, and tissue-resident immune cells are all stromal.
"Fibroblast" is the exclusive terminal compartment for the Vimentin-anchored connective-tissue
cells; pancreatic stellate cells are its fibroblast-like resident members (B3TUBB⁺ neural cells
that fall here are pulled out in Step 2). Vimentin remains a pan-mesenchymal MARKER — unchanged.

Why hierarchical + an Epithelial default (not plain argmax-of-`_norm` over all six): RESTORE
sets each threshold at mean+3σ of the marker's NEGATIVE population — right for the sparse
specific markers, but it inherently flags only the brightest ~20% of any ABUNDANT marker
(verified: Pan_Cytokeratin AND Ker8_18 both ≈ 20% `_pos`, even at low sigma, because the
Vimentin⁺ reference carries high cytokeratin background). The true acinar majority therefore
falls below every threshold, and a six-way `_norm` argmax scatters those dim cells to whatever
low-threshold marker has the largest ratio (empirically ~90% of "Immune" calls were such weak
cells). So: specific sparse lineages (endocrine/immune/endothelial) use the calibrated `_pos`;
the structural background is argmax over (Pan_Cytokeratin, Vimentin, SMA) `_norm`; and cells
positive for NONE of the three structural markers default to Epithelial — NOT arbitrary, these
"weak" cells are confirmed epithelial (median Ker8_18 549 vs 113, Pan_CK 274 vs 55 in the
positively-typed cells; flagged `epi_default` for Step-2 review).

Fibroblast/Endothelial/Muscle overlap at the marker level (Vimentin is pan-mesenchymal, marks
endothelium and activated stellate cells too); Step 2 within-lineage subclustering resolves
fibroblast vs. stellate vs. vasculature.

Inputs : data/restore_gated_redsea/donor_id=*/data_0.parquet  ({m}_norm for the 9 markers)
         data/cells/donor_id=*/data_0.parquet                 (cell_region, islet_num for context)
         ~/IO60panc2nd/donor_metadata_panc.xlsx               (disease.status)
Outputs: data/phenotype/broad/donor_id=*/data_0.parquet
         data/phenotype/broad_lineage_composition.png
"""
from __future__ import annotations
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2] / "data"
META = Path("/home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx")
GATED = BASE / "restore_gated_redsea"
OUT = BASE / "phenotype" / "broad"

LINEAGES = {
    "Epithelial": ["Pan_Cytokeratin"],
    "Fibroblast": ["Vimentin"],
    "Immune": ["CD3e", "CD20", "CD163"],
    "Endocrine": ["INS", "GCG"],
    "Endothelial": ["CD31"],
    "Muscle": ["SMA"],
}
LNAMES = list(LINEAGES)
LCOLORS = {"Epithelial": "#4477AA", "Fibroblast": "#EE6677", "Immune": "#228833",
           "Endocrine": "#CCBB44", "Endothelial": "#66CCEE", "Muscle": "#AA3377"}
STATUS_ORDER = ["ND", "AAB", "T1D"]


def status_map():
    m = pd.read_excel(META)
    return dict(zip(m["donor_id"].astype(str), m["disease.status"].astype(str)))


STRUCT = ["Epithelial", "Fibroblast", "Muscle"]


def assign_donor(donor, gated_f):
    g = pd.read_parquet(gated_f)
    norm = {m: g[f"{m}_norm"].to_numpy(np.float64) for ms in LINEAGES.values() for m in ms}
    pos = {m: g[f"{m}_pos"].to_numpy(bool) for ms in LINEAGES.values() for m in ms}
    scores = {L: np.max([norm[m] for m in LINEAGES[L]], axis=0) for L in LNAMES}

    # structural background: argmax over Pan_Cytokeratin / Vimentin / SMA _norm
    struct = np.vstack([scores[L] for L in STRUCT]).T            # N x 3
    sorder = np.argsort(-struct, axis=1)
    # dtype "<U16" so later assignments of longer labels ("Endothelial", 11 chars) are not
    # silently truncated by the narrow fixed width of np.array(STRUCT) (longest = "Epithelial", 10).
    struct_call = np.array(STRUCT)[sorder[:, 0]].astype("<U16")
    ss = np.take_along_axis(struct, sorder, axis=1)
    struct_margin = ss[:, 0] - ss[:, 1]
    # cells below every structural threshold → Epithelial. NOT arbitrary: these "weak" cells are
    # CONFIRMED epithelial (median Ker8_18 549 vs 113 in positively-typed cells, Pan_CK 274 vs 55) —
    # RESTORE's mean+3σ inherently flags only the brightest ~20% of any abundant marker (verified:
    # Ker8_18 also ~22% _pos), so the dim acinar majority is sub-threshold but genuinely epithelial.
    # (Acinar/ductal/fibroblast/muscle resolution within the structural compartment = Step 2 clustering.)
    struct_pos = pos["Pan_Cytokeratin"] | pos["Vimentin"] | pos["SMA"]
    struct_call[~struct_pos] = "Epithelial"

    # hierarchy: specific _pos lineages take precedence over the structural background
    endocrine = pos["INS"] | pos["GCG"]
    immune = pos["CD3e"] | pos["CD20"] | pos["CD163"]
    endoth = pos["CD31"]
    broad = struct_call.copy()
    margin = struct_margin.astype(np.float64)
    broad[endoth] = "Endothelial"; margin[endoth] = scores["Endothelial"][endoth] - 1.0
    broad[immune] = "Immune";      margin[immune] = scores["Immune"][immune] - 1.0
    broad[endocrine] = "Endocrine"; margin[endocrine] = scores["Endocrine"][endocrine] - 1.0
    specific = endocrine | immune | endoth
    epi_default = (~specific) & (~struct_pos)   # epithelial by background (validated), flagged for Step-2

    out = pd.DataFrame({"object_id": g["object_id"].astype(str), "donor_id": donor,
                        "broad_lineage": broad, "assign_margin": margin.astype(np.float32),
                        "epi_default": epi_default})
    for L in LNAMES:
        out[f"score_{L}"] = scores[L].astype(np.float32)
    # attach cell_region / islet_num for downstream
    cf = sorted(glob.glob(str(BASE / "cells" / f"donor_id={donor}" / "*.parquet")))[0]
    ctx = pd.read_parquet(cf, columns=["object_id", "cell_region", "islet_num"])
    out = out.merge(ctx.astype({"object_id": str}), on="object_id", how="left")
    return out


def main():
    smap = status_map()
    OUT.mkdir(parents=True, exist_ok=True)
    comp = {}
    weak_total = 0
    n_total = 0
    for gf in sorted(glob.glob(str(GATED / "donor_id=*" / "*.parquet"))):
        donor = gf.split("donor_id=")[1].split("/")[0]
        out = assign_donor(donor, gf)
        (OUT / f"donor_id={donor}").mkdir(parents=True, exist_ok=True)
        out.to_parquet(OUT / f"donor_id={donor}" / "data_0.parquet", index=False)
        vc = out["broad_lineage"].value_counts()
        comp[donor] = (vc / len(out) * 100).reindex(LNAMES).fillna(0)
        weak_total += int(out["epi_default"].sum()); n_total += len(out)
        print(f"[{donor}] {len(out):,} cells | " + " ".join(f"{L[:3]}={comp[donor][L]:.0f}%" for L in LNAMES) +
              f" | epi_default={100*out['epi_default'].mean():.1f}% | Unassigned=0", flush=True)
    C = pd.DataFrame(comp).T  # donors x lineages
    C["status"] = [smap.get(d, "?") for d in C.index]
    print(f"\n[total] {n_total:,} cells, 0 Unassigned, {100*weak_total/n_total:.1f}% epithelial-by-background "
          "(sub-threshold but confirmed epithelial via Ker8_18; resolved in Step 2)")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(22, 6))
    # A: per-donor stacked composition (ordered by status)
    Cs = C.sort_values("status")
    bottom = np.zeros(len(Cs))
    for L in LNAMES:
        ax[0].bar(range(len(Cs)), Cs[L], bottom=bottom, color=LCOLORS[L], label=L)
        bottom += Cs[L].to_numpy()
    ax[0].set_xticks(range(len(Cs))); ax[0].set_xticklabels([f"{d}\n{Cs.loc[d,'status']}" for d in Cs.index], fontsize=7)
    ax[0].set_ylabel("% of cells"); ax[0].set_title("Broad-lineage composition per donor (0 Unassigned)")
    ax[0].legend(fontsize=8, ncol=2, loc="upper right")
    # B: mean composition by status
    xs = np.arange(len(LNAMES)); bw = 0.25
    for si, s in enumerate(STATUS_ORDER):
        m = C[C.status == s][LNAMES].mean()
        ax[1].bar(xs + (si - 1) * bw, m.values, bw, label=s)
    ax[1].set_xticks(xs); ax[1].set_xticklabels(LNAMES, rotation=45, ha="right", fontsize=9)
    ax[1].set_ylabel("mean % of cells"); ax[1].set_title("Composition by disease status"); ax[1].legend()
    # C: Endocrine (β-loss) + Immune (T-gain) by status — the biology check
    endo = [C[C.status == s]["Endocrine"].mean() for s in STATUS_ORDER]
    imm = [C[C.status == s]["Immune"].mean() for s in STATUS_ORDER]
    ax[2].plot(STATUS_ORDER, endo, "o-", color=LCOLORS["Endocrine"], lw=2, ms=9, label="Endocrine (expect ↓)")
    ax[2].plot(STATUS_ORDER, imm, "s-", color=LCOLORS["Immune"], lw=2, ms=9, label="Immune (expect ↑)")
    for s_i, s in enumerate(STATUS_ORDER):
        for d in C[C.status == s].index:
            ax[2].scatter(s_i, C.loc[d, "Endocrine"], color=LCOLORS["Endocrine"], s=14, alpha=0.5)
    ax[2].set_ylabel("mean % of cells"); ax[2].set_title("Disease trend: endocrine loss + immune gain"); ax[2].legend()
    fig.suptitle("Step 1 — broad lineage (argmax RESTORE _norm on REDSEA-corrected data)", fontsize=14)
    fig.tight_layout(); fig.savefig(BASE / "phenotype" / "broad_lineage_composition.png", dpi=150, bbox_inches="tight")
    print(f"[saved] {BASE/'phenotype'/'broad_lineage_composition.png'}")
    print("\ncomposition by status (mean %):")
    print(C.groupby("status")[LNAMES].mean().round(1).to_string())


if __name__ == "__main__":
    raise SystemExit(main())
