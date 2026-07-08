#!/usr/bin/env python3
"""
Step 1 — broad lineage assignment (every cell typed, ZERO "Unassigned").

On the REDSEA-corrected, re-RESTORE'd data, each cell is typed by a HIERARCHICAL rule (every
cell typed — no "Unassigned" bucket; the user mandate is to fix the upstream cause, not invent
a junk category):

  1. Endocrine    if INS/GCG/SST `_pos`             (islet hormones β/α/δ — specific, sparse)
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
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from marker_taxonomy import CD99_BRIGHT   # CD99+ := bright only (≥ CD99_BRIGHT × per-image RESTORE thr)

BASE = Path(__file__).resolve().parents[2] / "data"
# Prefer the in-repo copy (move-proof); fall back to the original ~/IO60panc2nd location.
META = next((p for p in (BASE / "donor_metadata_panc.xlsx",
                         Path("/home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx"))
             if p.exists()), BASE / "donor_metadata_panc.xlsx")
GATED = BASE / "restore_gated_redsea"
# B3TUBB/CD99/MPO were gated separately (Pan_Cytokeratin ref) into this dir so the validated
# 10-marker gates in restore_gated_redsea stay byte-identical; merged per donor in assign_donor().
EXTRA_GATED = BASE / "restore_gated_redsea_extra"
EXTRA_MARKERS = ["B3TUBB", "CD99", "MPO"]
OUT = BASE / "phenotype" / "broad"

LINEAGES = {
    "Epithelial": ["Pan_Cytokeratin"],
    "Fibroblast": ["Vimentin"],
    "Immune": ["CD3e", "CD20", "CD163"],
    "Endocrine": ["INS", "GCG", "SST", "CD99"],   # +CD99 (islet-only in this panel)
    "Endothelial": ["CD31"],
    "Muscle": ["SMA"],
    "Neural": ["B3TUBB"],                         # pulled out of Fibroblast (B3TUBB+ nerves/Schwann)
    "Neutrophil": ["MPO"],                        # granulocytes, carved out of the structural bucket
}
LNAMES = list(LINEAGES)
LCOLORS = {"Epithelial": "#4477AA", "Fibroblast": "#EE6677", "Immune": "#228833",
           "Endocrine": "#CCBB44", "Endothelial": "#66CCEE", "Muscle": "#AA3377",
           "Neural": "#B5838D", "Neutrophil": "#E69F00"}  # dusty rose / Okabe-Ito orange (CVD-safe)
# unique 3-char codes for the terse per-donor progress line (Endocrine/Endothelial and
# Neural/Neutrophil would otherwise both collide under a naive name[:3]).
ABBR = {"Epithelial": "Epi", "Fibroblast": "Fib", "Immune": "Imm", "Endocrine": "Enc",
        "Endothelial": "Eth", "Muscle": "Mus", "Neural": "Nrl", "Neutrophil": "Neu"}
STATUS_ORDER = ["ND", "AAB", "T1D"]


def status_map():
    m = pd.read_excel(META)
    return dict(zip(m["donor_id"].astype(str), m["disease.status"].astype(str)))


STRUCT = ["Epithelial", "Fibroblast", "Muscle"]


def assign_donor(donor, gated_f):
    g = pd.read_parquet(gated_f)
    # merge the separately-gated B3TUBB/CD99/MPO (Neural/Endocrine/Neutrophil markers) by object_id;
    # keeps the validated 10-marker gates untouched. NaN (a cell absent from the extra run) → not-pos.
    ex_cols = ["object_id"] + [f"{m}_{s}" for m in EXTRA_MARKERS for s in ("pos", "norm")]
    ex = pd.read_parquet(EXTRA_GATED / f"donor_id={donor}" / "data_0.parquet", columns=ex_cols)
    g = g.merge(ex.astype({"object_id": str}), on="object_id", how="left")
    for m in EXTRA_MARKERS:
        g[f"{m}_pos"] = g[f"{m}_pos"].fillna(False).astype(bool)
        g[f"{m}_norm"] = g[f"{m}_norm"].fillna(0.0)
    # CD99 is broad; the standard mean+3σ gate over-calls a non-islet population. Keep only BRIGHT CD99
    # (≥ CD99_BRIGHT× the per-image RESTORE threshold) — validated islet-coherent (see marker_taxonomy).
    g["CD99_pos"] = g["CD99_norm"] >= CD99_BRIGHT
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

    # hierarchy: specific _pos lineages take precedence over the structural background.
    # Applied sequentially (last write wins) → precedence high→low:
    #   Endocrine > Immune > Endothelial > Neutrophil > Neural > structural.
    # Neural/Neutrophil sit lowest of the specific lineages: they pull B3TUBB+/MPO+ cells OUT of the
    # structural bucket but yield to a genuine endocrine/immune/endothelial call (MPO+CD163+ → Immune
    # macrophage; INS+B3TUBB+ → Endocrine β-cell with spillover).
    endocrine = pos["INS"] | pos["GCG"] | pos["SST"] | pos["CD99"]  # +CD99 (islet-only in this panel)
    immune = pos["CD3e"] | pos["CD20"] | pos["CD163"]
    endoth = pos["CD31"]
    neural = pos["B3TUBB"]
    neutrophil = pos["MPO"]
    broad = struct_call.copy()
    margin = struct_margin.astype(np.float64)
    broad[neural] = "Neural";          margin[neural] = scores["Neural"][neural] - 1.0
    broad[neutrophil] = "Neutrophil";  margin[neutrophil] = scores["Neutrophil"][neutrophil] - 1.0
    broad[endoth] = "Endothelial";     margin[endoth] = scores["Endothelial"][endoth] - 1.0
    broad[immune] = "Immune";          margin[immune] = scores["Immune"][immune] - 1.0
    broad[endocrine] = "Endocrine";    margin[endocrine] = scores["Endocrine"][endocrine] - 1.0
    specific = endocrine | immune | endoth | neural | neutrophil
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
    # β/α marker positivity fractions (% of ALL cells) for the disease-trend panel — the Endocrine
    # lineage lumps INS|GCG|SST so it masks β-loss (INS↓) behind α-persistence (GCG↑); track them apart.
    mfrac = {"INS": 100 * pos["INS"].mean(), "GCG": 100 * pos["GCG"].mean()}
    return out, mfrac


def main():
    smap = status_map()
    OUT.mkdir(parents=True, exist_ok=True)
    comp = {}
    mfracs = {}
    weak_total = 0
    n_total = 0
    for gf in sorted(glob.glob(str(GATED / "donor_id=*" / "*.parquet"))):
        donor = gf.split("donor_id=")[1].split("/")[0]
        out, mfrac = assign_donor(donor, gf)
        (OUT / f"donor_id={donor}").mkdir(parents=True, exist_ok=True)
        out.to_parquet(OUT / f"donor_id={donor}" / "data_0.parquet", index=False)
        vc = out["broad_lineage"].value_counts()
        comp[donor] = (vc / len(out) * 100).reindex(LNAMES).fillna(0)
        mfracs[donor] = mfrac
        weak_total += int(out["epi_default"].sum()); n_total += len(out)
        print(f"[{donor}] {len(out):,} cells | " + " ".join(f"{ABBR[L]}={comp[donor][L]:.0f}%" for L in LNAMES) +
              f" | epi_default={100*out['epi_default'].mean():.1f}% | Unassigned=0", flush=True)
    C = pd.DataFrame(comp).T  # donors x lineages
    C["status"] = [smap.get(d, "?") for d in C.index]
    M = pd.DataFrame(mfracs).T  # donors x {INS, GCG} marker positivity (% of all cells)
    M["status"] = [smap.get(d, "?") for d in M.index]
    print(f"\n[total] {n_total:,} cells, 0 Unassigned, {100*weak_total/n_total:.1f}% epithelial-by-background "
          "(sub-threshold but confirmed epithelial via Ker8_18; resolved in Step 2)")

    # ---------------------------------------------------------------- figure
    plt.rcParams.update({"font.size": 15, "xtick.labelsize": 14, "ytick.labelsize": 14,
                         "axes.titlesize": 17, "axes.labelsize": 16, "legend.fontsize": 14})
    _rank = {s: i for i, s in enumerate(STATUS_ORDER)}
    Cs = C.sort_values("status", key=lambda col: col.map(lambda s: _rank.get(s, len(STATUS_ORDER))))
    SLABEL = {"ND": "ND", "AAB": "Aab+", "T1D": "T1D"}
    SCOL = {"ND": "#4477AA", "AAB": "#CC6633", "T1D": "#228833"}   # Paul-Tol status colours (00_globals DONOR_COLORS)
    present_status = [s for s in STATUS_ORDER if (C["status"] == s).any()]
    # per-donor % for each (status, lineage) — the units of replication are donors, not cells.
    gdata = {s: {L: C.loc[C["status"] == s, L].to_numpy(np.float64) for L in LNAMES} for s in present_status}
    PAIRS = [(present_status[i], present_status[j])
             for i in range(len(present_status)) for j in range(i + 1, len(present_status))]

    def _breakmarks(ax_up, ax_dn):   # diagonal cut marks at the shared broken-axis edge
        d = 0.008
        for ax_, ys in ((ax_up, (-d, +d)), (ax_dn, (1 - d, 1 + d))):
            kw = dict(transform=ax_.transAxes, color="k", clip_on=False, lw=1)
            ax_.plot((-d, +d), ys, **kw); ax_.plot((1 - d, 1 + d), ys, **kw)

    # ---- statistics: per-lineage Kruskal–Wallis omnibus (BH across the 8 lineages), then pairwise
    #      Mann–Whitney U (BH across all lineage×pair tests). Non-parametric: small, unequal n
    #      (ND 8 / Aab+ 3 / T1D 11) and genuinely non-normal per-donor composition (working rule:
    #      donor biology is variable — don't assume normality). Brackets are drawn only where the
    #      omnibus survives (gates the pairwise family against false positives).
    omni_p = {}
    for L in LNAMES:
        arrs = [gdata[s][L] for s in present_status if len(gdata[s][L]) > 0]
        try:
            omni_p[L] = float(kruskal(*arrs).pvalue) if len(arrs) >= 2 else np.nan
        except ValueError:      # all values identical across groups
            omni_p[L] = np.nan
    _Lok = [L for L in LNAMES if np.isfinite(omni_p[L])]
    omni_q = {L: np.nan for L in LNAMES}
    if _Lok:
        omni_q.update(dict(zip(_Lok, multipletests([omni_p[L] for L in _Lok], method="fdr_bh")[1])))

    pw = []   # (lineage, a, b, p)
    for L in LNAMES:
        for a, b in PAIRS:
            xa, xb = gdata[a][L], gdata[b][L]
            try:
                p = float(mannwhitneyu(xa, xb, alternative="two-sided").pvalue) \
                    if (len(xa) and len(xb) and (len(xa) + len(xb) >= 3)) else np.nan
            except ValueError:
                p = np.nan
            pw.append([L, a, b, p])
    pw = pd.DataFrame(pw, columns=["L", "a", "b", "p"])
    _fin = pw["p"].notna()
    pw["q"] = np.nan
    if _fin.any():
        pw.loc[_fin, "q"] = multipletests(pw.loc[_fin, "p"].to_numpy(), method="fdr_bh")[1]
    qlut = {(r.L, r.a, r.b): r.q for r in pw.itertuples()}

    def _stars(q):
        if not np.isfinite(q):
            return None
        return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else None

    fig = plt.figure(figsize=(23, 9))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.95], wspace=0.13)
    # Panel A stacked: the 7 non-Epithelial bands sit high in the 0–100 stack → give the TOP segment
    # more height so Fibroblast/Immune/Endocrine/… differences stay legible.
    gsA = outer[0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
    a_top, a_bot = fig.add_subplot(gsA[0]), fig.add_subplot(gsA[1])
    # Panel B: one boxplot subpanel per lineage (2×4) — each auto-scales to its own range (Epithelial
    # ~60–70% vs Immune/Endothelial < ~10%), which a single shared axis could never show.
    gsB = outer[1].subgridspec(2, 4, hspace=0.62, wspace=0.34)
    bax = {L: fig.add_subplot(gsB[i // 4, i % 4]) for i, L in enumerate(LNAMES)}

    # A: per-donor stacked composition — split y-axis (Epithelial dominates the low range; the split
    # expands the top so Fibroblast/Immune/Endocrine/… differences are legible). Ordered ND→Aab+→T1D.
    xA = np.arange(len(Cs)); bottom = np.zeros(len(Cs)); handles = []
    for L in LNAMES:
        h = a_top.bar(xA, Cs[L], bottom=bottom, color=LCOLORS[L], width=0.85)
        a_bot.bar(xA, Cs[L], bottom=bottom, color=LCOLORS[L], width=0.85)
        bottom += Cs[L].to_numpy(); handles.append(h)
    a_hi = max(15.0, float(np.floor(Cs["Epithelial"].min())) - 4)      # break below the lowest Epithelial
    a_top.set_ylim(a_hi, 100); a_bot.set_ylim(0, 8)
    a_top.spines["bottom"].set_visible(False); a_bot.spines["top"].set_visible(False)
    a_top.tick_params(bottom=False, labelbottom=False); _breakmarks(a_top, a_bot)
    a_bot.set_xticks(xA)   # rotate donor labels vertical so 22 large labels never overlap
    a_bot.set_xticklabels([f"{d}  {Cs.loc[d,'status']}" for d in Cs.index], rotation=90, fontsize=12)

    # B: per-lineage box plots by disease status — box = median/IQR, diamond = mean ± SEM, dots =
    # individual donors (deterministic jitter). Significant pairwise differences bracketed above.
    for L in LNAMES:
        ax = bax[L]
        data = [gdata[s][L] for s in present_status]
        bp = ax.boxplot(data, positions=list(range(len(present_status))), widths=0.6,
                        patch_artist=True, showfliers=False, showmeans=False, zorder=1,
                        medianprops=dict(color="black", lw=1.3))
        for patch, s in zip(bp["boxes"], present_status):
            patch.set_facecolor(SCOL[s]); patch.set_alpha(0.30); patch.set_edgecolor(SCOL[s]); patch.set_linewidth(1.2)
        for grp in ("whiskers", "caps"):
            for ln in bp[grp]:
                ln.set_color("#555555")
        for i, s in enumerate(present_status):
            y = data[i]
            if len(y) == 0:
                continue
            jit = np.linspace(-0.17, 0.17, len(y)) if len(y) > 1 else np.array([0.0])
            ax.scatter(i + jit, y, s=24, color=SCOL[s], edgecolor="black", lw=0.4, alpha=0.9, zorder=3)
            m = float(np.mean(y)); sem = float(np.std(y, ddof=1) / np.sqrt(len(y))) if len(y) > 1 else 0.0
            ax.errorbar(i, m, yerr=sem, fmt="D", color="black", ms=6.5, mfc="white", mec="black",
                        capsize=4, elinewidth=1.4, zorder=4)
        ax.set_xticks(range(len(present_status)))
        ax.set_xticklabels([SLABEL[s] for s in present_status], fontsize=12)
        if LNAMES.index(L) % 4 == 0:
            ax.set_ylabel("% of cells", fontsize=13)
        ax.set_title(L, fontsize=14, color=LCOLORS[L])

        # y-limits framed to THIS lineage's own data range (equal fractional padding) so Epithelial
        # (~55–72%) and Endothelial (~7–13%) are not compressed against a shared 0 baseline.
        allv = np.concatenate([d for d in data if len(d)]) if any(len(d) for d in data) else np.array([0.0, 1.0])
        lo, hi = float(allv.min()), float(allv.max())
        span = (hi - lo) if hi > lo else max(1.0, 0.1 * hi)
        pad = 0.15 * span
        # significance brackets — only when the omnibus survives BH; stack to avoid overlap.
        step = 0.16 * span; base = hi + 0.10 * span; level = 0
        if np.isfinite(omni_q.get(L, np.nan)) and omni_q[L] < 0.05:
            for a, b in PAIRS:
                st = _stars(qlut.get((L, a, b), np.nan))
                if st is None:
                    continue
                x1, x2 = present_status.index(a), present_status.index(b)
                yb = base + level * step
                ax.plot([x1, x1, x2, x2], [yb, yb + 0.03 * span, yb + 0.03 * span, yb], color="black", lw=1.1)
                ax.text((x1 + x2) / 2, yb + 0.035 * span, st, ha="center", va="bottom", fontsize=13)
                level += 1
        ax.set_ylim(lo - pad, (base + (level + 0.4) * step) if level else hi + pad)

    # centered y-axis title spanning BOTH broken segments of panel A, just left of the tick labels.
    ap_t, ap_b = a_top.get_position(), a_bot.get_position()
    fig.text(ap_t.x0 - 0.032, (ap_b.y0 + ap_t.y1) / 2, "% of cells",
             rotation="vertical", va="center", ha="center", fontsize=16)
    # 8-class lineage legend BELOW panel A (2 rows × 4), centered under panel A — anchored to a_bot so
    # it tracks panel A's horizontal extent.
    a_bot.legend(handles=handles, labels=LNAMES, loc="upper center", ncol=4,
                 fontsize=13, frameon=False, bbox_to_anchor=(0.5, -0.9))
    # Panel A + Panel B titles on the SAME baseline (top edge of both columns' axes), so they read as
    # one row beneath the suptitle; the columns share the gridspec bottom, so bottoms stay even.
    ap = a_top.get_position()
    bp0 = bax[LNAMES[0]].get_position(); bp3 = bax[LNAMES[3]].get_position()
    ty = ap.y1 + 0.05   # above the per-lineage subpanel titles (which sit at the top-row axes' edge)
    fig.text((ap.x0 + ap.x1) / 2, ty, "Broad-lineage composition per donor",
             ha="center", va="bottom", fontsize=17)
    fig.text((bp0.x0 + bp3.x1) / 2, ty, "Composition by disease status (per-donor)",
             ha="center", va="bottom", fontsize=17)
    fig.text((bp0.x0 + bp3.x1) / 2, 0.045,
             "Box = median/IQR · whiskers = range · diamond = mean ± SEM · dots = donors (ND 8 / Aab+ 3 / T1D 11).",
             ha="center", va="center", fontsize=11)
    fig.suptitle("Cell Phenotyping Stage 1: Broad Lineages", fontsize=20, y=1.0)
    fig.savefig(BASE / "phenotype" / "broad_lineage_composition.png", dpi=150, bbox_inches="tight")
    print(f"[saved] {BASE/'phenotype'/'broad_lineage_composition.png'}")
    print("\ncomposition by status (mean %):")
    print(C.groupby("status")[LNAMES].mean().reindex(STATUS_ORDER).round(1).to_string())
    print("\ndisease trend (mean % of cells positive, by status):")
    print(pd.concat([M.groupby("status")[["INS", "GCG"]].mean(),
                     C.groupby("status")["Immune"].mean()], axis=1)
            .reindex(STATUS_ORDER).round(2).to_string())


if __name__ == "__main__":
    raise SystemExit(main())
