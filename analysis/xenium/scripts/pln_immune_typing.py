#!/usr/bin/env python
"""SUPERSEDED 2026-09-18. This is v1 of the pLN Xenium immune typer. The current one is
``analysis/xenium/scripts/pln_immune_typing_v3.py`` (v2 is superseded too), and the release it
writes is what ``scripts/senior/xenium_release.py`` resolves. Recorded as a verified orphan in
``scripts/retire/dispositions.py`` (code tier 1). Kept as the record of the detection predicate
and the Y-probe control below; do not read it as the current method.

Two-stage immune typing of the pLN Xenium sections (maintainer's method, 2026-09-17).

    STAGE 1  separate immune from non-immune with CANONICAL genes.
    STAGE 2  within the immune compartment, resolve identity with the core/panels SUBPANELS.

Why the pancreas supplies nothing but the method here.  The subpanel classification draws
reference levels from "the ordinary tissue".  In a lymph node the ordinary tissue IS
lymphocytes, so a pLN-derived CD3E level sits at T-zone intensity and calls almost nothing
positive.  This lane therefore uses **no pLN-derived reference level at all**: stage 1 is a
per-cell detection predicate, not a population threshold, and stage 2 is a within-cell argmax
of absolute mean lognorm across the identity cores.  Nothing is pooled across donors.

DETECTION PREDICATE -- one transcript, above that cell's own measured background.

    detected(c, g)  <=>  raw(c, g) >= 1  AND  raw(c, g) > lambda_c
    lambda_c        =    (genomic_control_counts_c + a * mu_c) / (n_probes + a),  a = n_probes = 21
    mu_c            =    mean(genomic_control_counts | cell-area decile of c) / n_probes

i.e. the expected gDNA counts PER GENE in cell c, ``docs/XENIUM_PIPELINE.md`` section 1, shrunk
because the per-cell estimate is a Poisson draw over 21 probes.  Counts are integers, so the
effective floor is ``max(1, floor(lambda_c) + 1)``.

WHAT THIS REPLACED, AND WHY.  The first version of this lane detected at ``lognorm >= log1p(1)``
on a per-cell normalisation to the section MEDIAN depth.  That is algebraically
``raw(c,g) >= depth(c) / median_depth``: for the ~half of each section above median depth it
demands MORE THAN ONE transcript, and in the deepest section (6442, median 556) a cell at depth
2800 must carry 6.  That is a raw-count floor by another name, it scales with cell depth, and
it is the reason ``no_canonical_evidence`` was the largest outcome cohort-wide.  The floor is
now the cell's own measured background instead of its depth relative to its neighbours.

THE CONTROL, not the argument.  ``docs/XENIUM_PIPELINE.md`` requires any predicate to hold the
zero-template Y probes at or below ~0.3%.  The floor group ``TBL1Y``/``USP9Y``/``EIF1AY`` is
zero-template in a FEMALE donor, and this script measures it per donor under BOTH predicates and
writes both into ``release.json``.  Measured on the five female donors here, the predicate above
holds them at 0.0009%-0.0057% -- two orders of magnitude inside the requirement -- so the depth
scaling was not what was buying the background control on these sections.  ``lognorm`` is still
what stage 2 SCORES on; it is simply no longer what decides detection.

SIZE BIAS.  Stage 1 compares the FRACTION of each canonical set detected, not the count, so a
larger set carries no advantage under background.  Stage 2's score is the absolute mean lognorm
over a core (size-normalised by construction); only its >= 2-gene corroboration floor scales
with core size, and this script MEASURES that null per donor from the Y probes rather than
arguing about it.

WHICH GENES MAY NAME A LINEAGE.  An identity core carries two kinds of gene, and
``core/panels`` -- the one authority -- already separates them in its own ``evidence`` column:

    evidence in {canonical, both}   note "canonical identity marker"  -> MAY NAME the lineage
    evidence == cellxgene           note "marks lineage (CL)"         -> SUPPORTING ONLY

This is ``lineages.SUPPORTING_ONLY`` for the Xenium panel: a supporting gene still enters the
score and still describes the cell, it simply cannot be the evidence that admits the label.  The
first version of this lane counted BOTH kinds toward the >= 2 corroboration floor, and the
consequence was measured: 5,417 of 9,780 Granulocyte calls (55.4%) detected none of
MPO/ELANE/FCGR3B/CSF3R/CEACAM8/PRTN3, and 9,044 of 17,177 Immune_NK_ILC calls (52.7%) detected
none of the NK canonicals -- they were resting on DDX3X, SMARCA2, XPO6, FLOT2 and the like.  A
call resting on no detected canonical marker of its own lineage is not a call, so the floor now
counts NAMING genes only and such a cell becomes ``Immune_unresolved`` with the reason
``no_naming_evidence``.  It is NEVER reassigned to another core: the evidence is absent, which
is not evidence for something else.

CELL SET.  The vendor ``cells.parquet`` is the complete cell set; every one of its cells gets a
row.  Cells absent from the typed matrix, and cells with ``nucleus_count == 0`` (a cell with no
nucleus is excluded from every calculation), are ``Unavailable`` with the reason recorded.

WHAT THIS DOES NOT DO.  No spatial claim, no distance, no neighbourhood: pancreas and pLN are
different blocks.  No non-immune lineage is resolved -- the ruling is about immune.  Exocrine
is one label and is not claimed here at all.

    .venv/bin/python analysis/xenium/scripts/pln_immune_typing.py --out <dir> [--donors ...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import time

import h5py
import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[3]
PANELS = REPO / "core" / "panels" / "data" / "identity_core"
XENIUM_PATHS = REPO / "core" / "phenocycler" / "data" / "integration" / "xenium_paths.csv"
GDNA_TRIAGE = REPO / "data" / "development" / "xenium_gdna" / "gdna_triage_2026-09-17.parquet"
H5AD_DIR = pathlib.Path("/data/xenium/h5ad_out")
XENIUM_ROOT = pathlib.Path("/data/xenium")
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

SCHEMA_VERSION = 2
RAW_FLOOR = 1                                # one transcript; never scaled by cell depth
GDNA_N_PROBES = 21                           # genomic control probes on this panel
GDNA_AREA_DECILES = 10                       # gDNA per unit area is flat; per cell it is not
ANCHOR_MIN_GENES = 2                         # decision #12: no single-gene call
ANCHOR_MARGIN_FRAC = 0.25                    # a clean stage-2 call beats the runner-up by 25%
#: core/panels' own ``evidence`` values for a gene that MAY NAME its lineage.  Everything else
#: in a core (``cellxgene`` / "marks lineage (CL)") is SUPPORTING ONLY -- it scores, it does not
#: corroborate.  See the module docstring.
NAMING_EVIDENCE = ("canonical", "both")

#: h5ad stem -> (donor_id, roi).  The multi-section donors are resolved by matching the h5ad's
#: cell count to the per-bundle row of the gDNA triage table; see ``resolve_roi``.
STEM_DONOR = {"HDL112": "112", "HDL115": "115"}

#: Zero-template ground truth.  These five Y probes are on the panel; in a FEMALE donor they
#: have no template of any kind, so every detection is background (docs/XENIUM_PIPELINE.md).
Y_PROBES = ("DDX3Y", "KDM5D", "TBL1Y", "USP9Y", "EIF1AY")
Y_FLOOR_GROUP = ("TBL1Y", "USP9Y", "EIF1AY")   # the non-cross-hybridising three -- THE control
Y_CROSSHYB = ("DDX3Y", "KDM5D")                # zero-template, but cross-hybridise with DDX3X/KDM5C

#: A second, sex-independent background set: the panel-audit categories with no template in a
#: lymph node either.  Weaker than the Y probes (no zero-template guarantee, only "no known
#: template here"), but it is measurable in the MALE donors, where the Y probes are not.
PANEL_AUDIT_CSV = REPO / "analysis" / "xenium" / "panel_audit.csv"
NO_TEMPLATE_CATEGORIES = ("skin_cornification", "photoreceptor_retinal", "nephron_specific")

# --------------------------------------------------------------------------- canonical sets
#: STAGE 1 immune canonicals.  Pan-leukocyte plus the lineage canonicals that corroborate it.
#: CD45/PTPRC is on this panel, unlike the PhenoCycler pLN panel, so it can anchor -- but it
#: never calls a cell alone (ANCHOR_MIN_GENES).
IMMUNE_CANONICAL = (
    "PTPRC",
    # T
    "CD3E", "CD3G", "CD2", "CD247", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "LEF1", "THEMIS",
    # B / plasma
    "MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5", "MZB1", "DERL3", "TNFRSF17",
    # NK / ILC
    "NKG7", "KLRD1", "KLRB1", "PRF1", "GZMB", "KLRF1", "NCR1",
    # myeloid
    "CD68", "ITGAX", "ITGAM", "CD14", "CSF1R", "MRC1", "CD163", "FCN1", "CD86", "MSR1",
    # granulocyte
    "MPO", "ELANE", "FCGR3B", "CSF3R", "CEACAM8", "PRTN3",
)
#: STAGE 1 structural canonicals -- NOT immune, used to exclude.  A lymph node's non-immune
#: compartment is endothelium (blood and lymphatic), the fibroblastic stroma, and the capsule.
STRUCTURAL_CANONICAL = (
    "PECAM1", "CDH5", "CLDN5", "EGFL7", "PROX1", "LYVE1",          # endothelium
    "PDGFRB", "PDGFRA", "ACTA2", "RGS5", "NOTCH3", "PDPN", "POSTN",  # stroma / mural / FRC
    "EPCAM", "KRT19", "CDH1",                                       # epithelium (capsule)
    # KRT14 is deliberately NOT here: analysis/xenium/panel_audit.csv marks it
    # skin_cornification / exclude_recommended=yes, i.e. a measured clean negative.
    "SOX10", "MPZ",                                                 # neural
)
#: Reported as evidence only.  The maintainer's note that clustering works in pLN because of
#: the regional cell types; these are the regional genes.  NO regional call is made from them.
REGIONAL = ("CR2", "CXCL13", "CCL19", "CCL21", "BCL6", "AICDA", "SELL", "CD27", "MKI67")

#: STAGE 2 label -> identity-core file.  Immune classes only: the ruling is about immune.
IMMUNE_CORES = {
    "Immune_T": "03_immune_T_cell.csv",
    "Immune_B_Plasma": "04_immune_B_plasma.csv",
    "Immune_NK_ILC": "05_immune_NK_ILC.csv",
    "Myeloid": "06_immune_myeloid.csv",
    "Granulocyte": "07_immune_granulocyte.csv",
}
#: The common vocabulary (phenocycler.integration.vocab).  NK folds into T_NK there.
NATIVE_TO_COMMON = {
    "Immune_T": "T_NK", "Immune_NK_ILC": "T_NK", "Immune_B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid", "Granulocyte": "Granulocyte",
}
AMBIGUOUS, UNAVAILABLE = "Ambiguous", "Unavailable"
IMMUNE_UNRESOLVED, NON_IMMUNE = "Immune_unresolved", "Non_immune"


# --------------------------------------------------------------------------- inputs
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_cores() -> dict:
    """The five immune identity cores, each split into NAMING and SUPPORTING genes.

    The split is ``core/panels``' own ``evidence`` column, not a judgement made here.
    """
    out = {}
    for label, fname in IMMUNE_CORES.items():
        t = pd.read_csv(PANELS / fname)
        nm = t["evidence"].astype(str).isin(NAMING_EVIDENCE)
        out[label] = {
            "genes": t["gene"].astype(str).tolist(),
            "naming": t.loc[nm, "gene"].astype(str).tolist(),
            "supporting": t.loc[~nm, "gene"].astype(str).tolist(),
            "hints": {h: t.loc[t["subtype_hint"] == h, "gene"].astype(str).tolist()
                      for h in sorted(t["subtype_hint"].dropna().unique())},
            "hint_naming": {h: t.loc[(t["subtype_hint"] == h) & nm, "gene"].astype(str).tolist()
                            for h in sorted(t["subtype_hint"].dropna().unique())},
            "file": str(PANELS / fname), "sha256": sha256_file(PANELS / fname)}
    return out


def pln_bundles() -> pd.DataFrame:
    """Local pLN bundles with a cells.parquet, for the eligible cohort only."""
    from phenocycler.cohort import ensure_eligible_donors
    t = pd.read_csv(XENIUM_PATHS, dtype=str)
    t = t[t["roi"].str.startswith("pln") & t["source_path"].notna()
          & (t["source_path"].str.strip() != "NA")].copy()
    t["source_path"] = t["source_path"].str.strip().str.rstrip("/")
    t["donor_id"] = t["donor_id"].str.strip()
    t["bundle"] = t["source_path"].map(lambda p: XENIUM_ROOT / p)
    t["cells_parquet"] = t["bundle"].map(lambda b: b / "cells.parquet")
    t["local"] = t["cells_parquet"].map(lambda p: p.is_file())
    keep = []
    for d in t["donor_id"].unique():
        try:
            ensure_eligible_donors([d], context="pLN Xenium immune typing")
            keep.append(d)
        except Exception:
            pass
    t["eligible"] = t["donor_id"].isin(keep)
    return t.reset_index(drop=True)


def h5ad_for(donor: str, roi: str, n_bundle_cells: int) -> pathlib.Path | None:
    """The pLN h5ad whose cell set is a subset of THIS bundle's cells.

    Multi-section donors have several h5ads (6579a/b/c); the section is resolved by subset
    containment against the bundle's own cell_id set -- never by filename order.
    """
    stem = {v: k for k, v in STEM_DONOR.items()}.get(donor, donor)
    cands = sorted(H5AD_DIR.glob(f"{stem}*_pLNXenium.h5ad"))
    return cands or None


def gdna_lambda(cells: pd.DataFrame) -> pd.Series:
    """Expected gDNA counts PER GENE in each cell -- docs/XENIUM_PIPELINE.md section 1.

        lambda_c = (genomic_control_counts_c + a * mu_c) / (n_probes + a),   a = n_probes = 21
        mu_c     = mean(genomic_control_counts | cell-area decile of c) / n_probes

    Shrunk toward the area-decile mean because the per-cell estimate is a Poisson draw over 21
    probes and most cells record zero genomic counts -- which does not make them zero-background
    cells.  Area is the conditioning variable because gDNA per unit area is roughly flat while
    gDNA per cell is not.  Indexed by ``cell_id``, in counts per gene, so it is directly
    comparable to a raw count.  It is a BACKGROUND floor, not a depth floor: nothing here is a
    function of the cell's transcript depth or of the section's median depth.
    """
    g = cells["genomic_control_counts"].astype(float)
    area = cells["cell_area"].astype(float)
    dec = pd.qcut(area.rank(method="first"), GDNA_AREA_DECILES, labels=False, duplicates="drop")
    mu = g.groupby(dec).transform("mean") / GDNA_N_PROBES
    lam = (g + GDNA_N_PROBES * mu) / (GDNA_N_PROBES + GDNA_N_PROBES)
    return pd.Series(lam.to_numpy(), index=cells["cell_id"].astype(str).to_numpy())


def detection_floor(lam: np.ndarray) -> np.ndarray:
    """Per-cell raw-count floor: one transcript, and strictly above that cell's own background.

    Counts are integers, so ``raw > lambda_c`` is ``raw >= floor(lambda_c) + 1``.
    """
    return np.maximum(float(RAW_FLOOR), np.floor(lam) + 1.0)


def read_matrix_subset(path, genes_wanted, floor_by_cell: pd.Series):
    """Return (cell_ids, depth, floor, detected_bool[n_cells, n_wanted], lognorm, genes, obs).

    ``floor_by_cell`` is the per-cell raw-count detection floor from :func:`detection_floor`,
    indexed by ``cell_id``.  Every cell of the typed matrix must be in it -- the h5ad is a
    strict subset of the bundle -- and a miss raises rather than defaulting to anything.
    """
    with h5py.File(path, "r") as h:
        var = np.array(h["var"]["_index"][:]).astype(str)
        cell_ids = np.array(h["obs"]["_index"][:]).astype(str)
        depth = np.asarray(h["obs"]["nCount_Xenium"][:], dtype=np.float64)
        obs_extra = {}
        for c in ("disease.status", "sex", "age.yrs", "AAb.count", "race", "batch"):
            if c in h["obs"]:
                node = h["obs"][c]
                if isinstance(node, h5py.Group):
                    cats = np.array(node["categories"][:]).astype(str)
                    obs_extra[c] = cats[np.asarray(node["codes"][:])]
                else:
                    obs_extra[c] = np.asarray(node[:])
        col = {g: i for i, g in enumerate(var)}
        present = [g for g in genes_wanted if g in col]
        absent = [g for g in genes_wanted if g not in col]
        X = h["X"]
        if X.attrs["encoding-type"] != "csc_matrix":
            raise RuntimeError(f"{path}: expected csc_matrix, got {X.attrs['encoding-type']}")
        n_cells = int(X.attrs["shape"][0])
        indptr = X["indptr"][:]
        target = float(np.median(depth[depth > 0])) if (depth > 0).any() else 1.0
        # detected(c,g) <=> raw >= 1 AND raw > lambda_c.  NOT a function of depth.
        need = floor_by_cell.reindex(cell_ids)
        if need.isna().any():
            raise RuntimeError(f"{path}: {int(need.isna().sum())} typed cells are absent from the "
                               "bundle cells.parquet, so their background floor is unavailable")
        need = need.to_numpy(dtype=np.float64)
        det = np.zeros((n_cells, len(present)), dtype=bool)
        dset_d, dset_i = X["data"], X["indices"]
        for k, g in enumerate(present):
            j = col[g]
            a, b = int(indptr[j]), int(indptr[j + 1])
            if b <= a:
                continue
            rows = dset_i[a:b]
            vals = dset_d[a:b]
            det[rows, k] = vals >= need[rows]
        lg = np.zeros((n_cells, len(present)), dtype=np.float32)
        for k, g in enumerate(present):
            j = col[g]
            a, b = int(indptr[j]), int(indptr[j + 1])
            if b <= a:
                continue
            rows = dset_i[a:b]
            vals = dset_d[a:b]
            lg[rows, k] = np.log1p(vals * (target / np.maximum(depth[rows], 1e-12)))
    return cell_ids, depth, need, target, det, lg, present, absent, obs_extra


def p_at_least(k: int, n: int, p: float) -> float:
    """P(at least k of n independent Bernoulli(p)) -- the background null for a core floor."""
    from math import comb
    return float(sum(comb(n, i) * p**i * (1 - p)**(n - i) for i in range(k, n + 1)))


# --------------------------------------------------------------------------- the two stages
def background_rates(h5ad_path, genes, floor_by_cell: pd.Series):
    """False-detection rate over ``genes`` under the predicate in force AND under the withdrawn one.

    Returns ``(rate_this_lane, rate_withdrawn_depth_scaled, genes_used)``.  Both are reported
    because the control is what decides which predicate is defensible, and a reader must be able
    to see the one that was replaced.
    """
    with h5py.File(h5ad_path, "r") as h:
        var = np.array(h["var"]["_index"][:]).astype(str)
        col = {g: i for i, g in enumerate(var)}
        cell_ids = np.array(h["obs"]["_index"][:]).astype(str)
        depth = np.asarray(h["obs"]["nCount_Xenium"][:], dtype=np.float64)
        n = len(depth)
        target = float(np.median(depth[depth > 0])) if (depth > 0).any() else 1.0
        need_depth = np.where(depth > 0, depth / target, np.inf)
        need_bg = floor_by_cell.reindex(cell_ids).to_numpy(dtype=np.float64)
        ip = h["X"]["indptr"][:]
        used = [g for g in genes if g in col]
        here = withdrawn = 0
        for g in used:
            a, b = int(ip[col[g]]), int(ip[col[g] + 1])
            if b <= a:
                continue
            rows = h["X"]["indices"][a:b]
            vals = h["X"]["data"][a:b]
            here += int((vals >= need_bg[rows]).sum())
            withdrawn += int((vals >= need_depth[rows]).sum())
    if not used:
        return None, None, []
    return here / (len(used) * n), withdrawn / (len(used) * n), used


def type_section(donor, roi, h5ad_path, cells_parquet, cores, report):
    t0 = time.time()
    core_genes = sorted({g for c in cores.values() for g in c["genes"]})
    wanted = sorted(set(IMMUNE_CANONICAL) | set(STRUCTURAL_CANONICAL) | set(REGIONAL)
                    | set(core_genes) | set(Y_PROBES))

    # ---- the cell set and the per-cell background floor, both from the bundle
    bundle = pd.read_parquet(cells_parquet,
                             columns=["cell_id", "nucleus_count", "total_counts", "cell_area",
                                      "genomic_control_counts", "control_probe_counts"])
    bundle["cell_id"] = bundle["cell_id"].astype(str)
    lam = gdna_lambda(bundle)
    floor_by_cell = pd.Series(detection_floor(lam.to_numpy()), index=lam.index)

    cell_ids, depth, floor, target, det, lg, present, absent, obs_extra = read_matrix_subset(
        h5ad_path, wanted, floor_by_cell)
    gi = {g: k for k, g in enumerate(present)}
    n = len(cell_ids)

    # ---- background, measured on the zero-template Y probes (female donors only)
    sex = str(obs_extra.get("sex", np.array(["?"]))[0]) if "sex" in obs_extra else "?"
    yp = [g for g in Y_FLOOR_GROUP if g in gi]
    if sex.lower().startswith("f") and yp:
        y_rate = float(det[:, [gi[g] for g in yp]].mean())
        y_basis = f"{len(yp)} zero-template Y probes ({','.join(yp)}), donor sex={sex}"
    else:
        y_rate, y_basis = float("nan"), f"UNMEASURED: donor sex={sex} (Y probes are not zero-template)"

    au = pd.read_csv(PANEL_AUDIT_CSV)
    no_template = au.loc[au["category"].isin(NO_TEMPLATE_CATEGORIES), "gene"].astype(str).tolist()
    female = sex.lower().startswith("f")
    y_here, y_withdrawn, y_used = background_rates(h5ad_path, Y_FLOOR_GROUP, floor_by_cell)
    x_here, x_withdrawn, x_used = background_rates(h5ad_path, Y_CROSSHYB, floor_by_cell)
    n_here, n_withdrawn, n_used = background_rates(h5ad_path, no_template, floor_by_cell)
    bg = {
        "requirement": "docs/XENIUM_PIPELINE.md: any detection predicate must hold a "
                       "zero-template set at or below ~0.3%.",
        "y_probes_floor_group": {
            "genes": y_used, "zero_template": female,
            "rate_this_lane": (round(y_here, 8) if female else None),
            "rate_withdrawn_depth_scaled": (round(y_withdrawn, 8) if female else None),
            "note": ("zero-template in this female donor: every detection is background"
                     if female else
                     "UNMEASURED: donor is male, these probes have template")},
        "y_probes_crosshybridising": {
            "genes": x_used, "zero_template": female,
            "rate_this_lane": (round(x_here, 8) if female else None),
            "rate_withdrawn_depth_scaled": (round(x_withdrawn, 8) if female else None),
            "note": "DDX3Y/KDM5D are zero-template too but cross-hybridise with the highly "
                    "expressed X paralogs DDX3X/KDM5C, so their elevation over the floor group "
                    "is cross-hybridisation, not background. Reported, not used as the control."},
        "panel_audit_no_template": {
            "categories": list(NO_TEMPLATE_CATEGORIES), "n_genes": len(n_used),
            "rate_this_lane": round(n_here, 8), "rate_withdrawn_depth_scaled": round(n_withdrawn, 8),
            "note": "no KNOWN template in a lymph node; weaker than the Y probes, but "
                    "sex-independent, so it is the only background measurement available in the "
                    "male donors."},
    }

    # ---- STAGE 1
    imm = [g for g in IMMUNE_CANONICAL if g in gi]
    stc = [g for g in STRUCTURAL_CANONICAL if g in gi]
    n_imm = det[:, [gi[g] for g in imm]].sum(axis=1)
    n_str = det[:, [gi[g] for g in stc]].sum(axis=1)
    f_imm = n_imm / len(imm)
    f_str = n_str / len(stc)
    ptprc = det[:, gi["PTPRC"]] if "PTPRC" in gi else np.zeros(n, bool)

    stage1 = np.full(n, AMBIGUOUS, dtype=object)
    reason = np.full(n, "", dtype=object)
    is_imm = (n_imm >= ANCHOR_MIN_GENES) & (f_imm > f_str)
    is_str = (n_str >= ANCHOR_MIN_GENES) & (f_str > f_imm)
    stage1[is_imm] = "Immune"
    stage1[is_str] = "Non_immune"
    rest = ~(is_imm | is_str)
    reason[rest & (n_imm < ANCHOR_MIN_GENES) & (n_str < ANCHOR_MIN_GENES)] = "no_canonical_evidence"
    reason[rest & (n_imm >= ANCHOR_MIN_GENES) & (n_str >= ANCHOR_MIN_GENES)] = "compartment_tie"
    reason[rest & (reason == "")] = "below_corroboration"

    # ---- STAGE 2 (immune cells only)
    labels = list(cores)
    scores = np.zeros((n, len(labels)), dtype=np.float32)
    ndet = np.zeros((n, len(labels)), dtype=np.int32)
    ndet_naming = np.zeros((n, len(labels)), dtype=np.int32)
    core_on_panel, naming_on_panel = {}, {}
    for k, L in enumerate(labels):
        cols = [gi[g] for g in cores[L]["genes"] if g in gi]
        ncols = [gi[g] for g in cores[L]["naming"] if g in gi]
        core_on_panel[L] = len(cols)
        naming_on_panel[L] = len(ncols)
        if not cols:
            continue
        # the SCORE is over the whole core: a supporting gene still describes the cell.
        scores[:, k] = lg[:, cols].mean(axis=1)
        ndet[:, k] = det[:, cols].sum(axis=1)
        # the FLOOR is over the naming genes only: a supporting gene cannot admit the label.
        if ncols:
            ndet_naming[:, k] = det[:, ncols].sum(axis=1)
    order = np.argsort(-scores, axis=1)
    win = order[:, 0]
    run = order[:, 1]
    s1 = scores[np.arange(n), win]
    s2 = scores[np.arange(n), run]
    margin_ok = s1 > 0
    margin_ok &= (s1 - s2) >= ANCHOR_MARGIN_FRAC * np.where(s1 > 0, s1, 1.0)
    win_naming = ndet_naming[np.arange(n), win]
    win_core = ndet[np.arange(n), win]
    # The corroboration floor counts NAMING genes only (module docstring).  A cell whose winning
    # core is corroborated only by supporting genes is UNCALLED -- never handed to another core.
    floor_ok = win_naming >= ANCHOR_MIN_GENES

    native = np.array(stage1, dtype=object)
    immune_mask = stage1 == "Immune"
    resolved = immune_mask & margin_ok & floor_ok
    native[immune_mask] = IMMUNE_UNRESOLVED
    native[resolved] = np.array(labels, dtype=object)[win[resolved]]
    # two distinguishable failures, because they mean different things
    no_naming = immune_mask & ~floor_ok & (win_core >= ANCHOR_MIN_GENES)
    reason[no_naming] = "no_naming_evidence"
    reason[immune_mask & ~floor_ok & ~no_naming] = "core_below_corroboration"
    reason[immune_mask & floor_ok & ~margin_ok] = "core_margin"

    # ---- fine type from the subtype hints of the winning core
    fine = np.full(n, "", dtype=object)
    for k, L in enumerate(labels):
        sel = resolved & (win == k)
        if not sel.any():
            continue
        hints = {h: [gi[g] for g in gs if g in gi] for h, gs in cores[L]["hints"].items()}
        hnam = {h: [gi[g] for g in gs if g in gi] for h, gs in cores[L]["hint_naming"].items()}
        # a hint needs >= 3 genes on panel to be scoreable, and at least one that may NAME it
        hints = {h: c for h, c in hints.items() if len(c) >= 3 and hnam.get(h)}
        if not hints:
            continue
        hn = list(hints)
        hs = np.stack([lg[np.ix_(sel, hints[h])].mean(axis=1) for h in hn], axis=1)
        # corroboration for a fine type is also naming-only
        hd = np.stack([det[np.ix_(sel, hnam[h])].sum(axis=1) for h in hn], axis=1)
        bo = np.argsort(-hs, axis=1)
        b0, b1 = bo[:, 0], bo[:, 1] if hs.shape[1] > 1 else bo[:, 0]
        a = np.arange(hs.shape[0])
        ok = (hd[a, b0] >= ANCHOR_MIN_GENES) & (hs[a, b0] > 0)
        if hs.shape[1] > 1:
            ok &= (hs[a, b0] - hs[a, b1]) >= ANCHOR_MARGIN_FRAC * np.maximum(hs[a, b0], 1e-12)
        idx = np.flatnonzero(sel)
        fine[idx[ok]] = np.array(hn, dtype=object)[b0[ok]]

    typed = pd.DataFrame({
        "cell_id": cell_ids,
        "lineage_native": native.astype(str),
        "celltype_fine": fine.astype(str),
        "stage1": np.array(stage1, dtype=str),
        "reason": reason.astype(str),
        "confidence": np.where(immune_mask, (s1 - s2) / np.where(s1 > 0, s1, np.nan), np.nan),
        "n_immune_canonical": n_imm.astype(np.int16),
        "n_structural_canonical": n_str.astype(np.int16),
        "n_core_genes_detected": win_core.astype(np.int16),
        "n_core_naming_detected": win_naming.astype(np.int16),
        "ptprc_detected": ptprc,
        "depth": depth.astype(np.int32),
        "detection_floor": floor.astype(np.float32),
    })
    for g in REGIONAL:
        typed[f"det__{g}"] = det[:, gi[g]] if g in gi else pd.array([pd.NA] * n, dtype="boolean")

    # ---- complete cell set from the bundle; nucleus rule
    out = bundle[["cell_id"]].merge(typed, on="cell_id", how="left")
    out["nucleus_count"] = bundle["nucleus_count"].to_numpy()
    out["control_probe_frac"] = (bundle["control_probe_counts"]
                                 / bundle["total_counts"].replace(0, np.nan)).to_numpy()
    out["genomic_control_counts"] = bundle["genomic_control_counts"].to_numpy()
    not_typed = out["lineage_native"].isna().to_numpy()
    out.loc[not_typed, "lineage_native"] = UNAVAILABLE
    out.loc[not_typed, "reason"] = "not_in_typed_matrix"
    no_nuc = (out["nucleus_count"].to_numpy() == 0)
    out.loc[no_nuc, "lineage_native"] = UNAVAILABLE
    out.loc[no_nuc, "reason"] = "no_nucleus"
    out.loc[no_nuc | not_typed, ["celltype_fine", "stage1"]] = ""
    out.loc[no_nuc | not_typed, "confidence"] = np.nan
    out["common_lineage"] = out["lineage_native"].map(NATIVE_TO_COMMON).fillna("").to_numpy()
    out["evidence_positive"] = out["lineage_native"].isin(list(NATIVE_TO_COMMON)).to_numpy()
    out["celltype_fine"] = out["celltype_fine"].fillna("").astype(str)
    out["stage1"] = out["stage1"].fillna("").astype(str)
    out["reason"] = out["reason"].fillna("").astype(str)
    out["ptprc_detected"] = pd.array(out["ptprc_detected"], dtype="boolean")
    for c in [f"det__{g}" for g in REGIONAL]:
        out[c] = pd.array(out[c], dtype="boolean")
    # the winning-core counts describe the stage-2 argmax and mean nothing outside the immune
    # compartment, so they are UNAVAILABLE there rather than 0.
    not_immune = (out["stage1"].to_numpy() != "Immune")
    for c in ("n_core_genes_detected", "n_core_naming_detected"):
        out.loc[not_immune, c] = np.nan
    for c in ("n_immune_canonical", "n_structural_canonical",
              "n_core_genes_detected", "n_core_naming_detected"):
        out[c] = pd.array(out[c], dtype="Int16")
    out["depth"] = pd.array(out["depth"], dtype="Int32")

    nulls = {L: (None if not np.isfinite(y_rate)
                 else round(p_at_least(ANCHOR_MIN_GENES, naming_on_panel[L], y_rate), 8))
             for L in labels}
    stats = {
        "donor_id": donor, "roi": roi,
        "n_cells_bundle": int(len(out)), "n_cells_typed_matrix": int(n),
        "n_not_in_typed_matrix": int(not_typed.sum()), "n_no_nucleus": int(no_nuc.sum()),
        "median_depth": target,
        "detection_floor": {
            "predicate": "raw >= 1 AND raw > lambda_c (gDNA background per gene, per cell)",
            "lambda_median": round(float(np.median(lam.to_numpy())), 8),
            "lambda_p99": round(float(np.quantile(lam.to_numpy(), 0.99)), 8),
            "lambda_max": round(float(np.max(lam.to_numpy())), 8),
            "n_cells_floor_above_1": int((floor_by_cell.to_numpy() > 1).sum()),
            "note": "the floor is 1 transcript wherever lambda_c < 1. It is a function of the "
                    "cell's measured gDNA background and area decile, never of its depth."},
        "withdrawn_depth_scaled_floor": {
            "predicate": "raw >= depth(c) / median_depth  (== lognorm >= log1p(1))",
            "median_depth": target,
            "n_cells_would_need_more_than_1_transcript": int((depth > target).sum()),
            "frac_cells_would_need_more_than_1_transcript": round(float((depth > target).mean()), 6),
            "max_raw_required": int(np.ceil(np.max(depth) / target)) if n else None,
            "note": "reported so the size of what was removed is readable, not applied."},
        "sex": sex, "disease_status": (str(obs_extra["disease.status"][0])
                                       if "disease.status" in obs_extra else "UNAVAILABLE"),
        "y_probe_background_rate": (None if not np.isfinite(y_rate) else round(y_rate, 8)),
        "y_probe_basis": y_basis,
        "background_rates": bg,
        "stage1_counts": out["stage1"].value_counts().to_dict(),
        "labels": out["lineage_native"].value_counts().to_dict(),
        "reasons": out.loc[out["reason"] != "", "reason"].value_counts().to_dict(),
        "fine": out.loc[out["celltype_fine"] != "", "celltype_fine"].value_counts().to_dict(),
        "ptprc_detected_frac_of_immune": (
            float(out.loc[out["stage1"] == "Immune", "ptprc_detected"].mean())
            if (out["stage1"] == "Immune").any() else None),
        "core_genes_on_panel": core_on_panel,
        "core_naming_genes_on_panel": naming_on_panel,
        "naming_floor_null_from_y_background": nulls,
        "naming_floor_null_basis": (
            "P(>= 2 of the core's NAMING genes detected by background alone) at this donor's "
            "measured zero-template Y-probe rate. null where that rate is UNMEASURED (male)."),
        "canonical_absent_from_panel": absent,
        "n_immune_canonical_on_panel": len(imm),
        "n_structural_canonical_on_panel": len(stc),
        "seconds": round(time.time() - t0, 1),
    }
    report.append(stats)
    return out, stats


# --------------------------------------------------------------------------- release
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "data/development/xenium_typing/release_pln_immune"))
    ap.add_argument("--donors", nargs="*", default=None)
    a = ap.parse_args(argv)
    out_dir = pathlib.Path(a.out)
    (out_dir / "calls").mkdir(parents=True, exist_ok=True)

    cores = load_cores()
    tbl = pln_bundles()
    runnable = tbl[tbl["local"] & tbl["eligible"]].copy()
    if a.donors:
        runnable = runnable[runnable["donor_id"].isin([str(d) for d in a.donors])]

    triage = pd.read_parquet(GDNA_TRIAGE) if GDNA_TRIAGE.is_file() else None
    report, inputs, per_donor = [], [], {}
    for _, r in runnable.sort_values("donor_id").iterrows():
        donor, roi = r["donor_id"], r["roi"]
        n_bundle = len(pd.read_parquet(r["cells_parquet"], columns=["cell_id"]))
        cands = h5ad_for(donor, roi, n_bundle)
        if not cands:
            print(f"[skip] {donor} {roi}: no pLN h5ad", flush=True)
            continue
        bundle_ids = set(pd.read_parquet(r["cells_parquet"], columns=["cell_id"])["cell_id"].astype(str))
        chosen = None
        for c in cands:
            with h5py.File(c, "r") as h:
                ids = np.array(h["obs"]["_index"][:20000]).astype(str)
            if set(ids).issubset(bundle_ids):
                chosen = c
                break
        if chosen is None:
            print(f"[skip] {donor} {roi}: no h5ad whose cells are a subset of this bundle", flush=True)
            continue
        print(f"[run ] donor {donor} roi={roi} <- {chosen.name}", flush=True)
        df, stats = type_section(donor, roi, chosen, r["cells_parquet"], cores, report)
        if triage is not None:
            m = triage[np.isclose(triage["n_cells"], stats["n_cells_bundle"])]
            stats["gdna_adjusted_genomic_control_probe_rate"] = (
                float(m["rate"].iloc[0]) if len(m) == 1 else None)
            stats["gdna_rate_basis"] = ("matched to gdna_triage_2026-09-17.parquet by bundle n_cells"
                                        if len(m) == 1 else "UNMEASURED: no unique triage row")
        p = out_dir / "calls" / f"donor_id={donor}" / f"roi={roi}"
        p.mkdir(parents=True, exist_ok=True)
        df.to_parquet(p / "data.parquet", index=False)
        stats["calls_sha256"] = sha256_file(p / "data.parquet")
        per_donor[donor] = stats
        inputs.append({"kind": "typed_matrix", "path": str(chosen), "sha256": sha256_file(chosen),
                       "donor_id": donor, "roi": roi})
        inputs.append({"kind": "cell_set", "path": str(r["cells_parquet"]),
                       "sha256": sha256_file(r["cells_parquet"]), "donor_id": donor, "roi": roi})
        print(json.dumps({k: stats[k] for k in ("n_cells_bundle", "labels", "y_probe_background_rate")},
                         indent=1), flush=True)

    space = {
        "schema_version": SCHEMA_VERSION,
        "space_id": "xenium-pln-immune-two-stage-v1",
        "producer_lane": "xenium-pln-immune-two-stage",
        "labels": [{"label": L, "common_lineage": NATIVE_TO_COMMON[L],
                    "reason": "stage 2: identity-core argmax within the immune compartment "
                              "(core/panels/data/identity_core), corroborated by >= 2 detected "
                              "genes that MAY NAME this lineage (evidence in "
                              f"{list(NAMING_EVIDENCE)}), not merely >= 2 core genes"}
                   for L in IMMUNE_CORES]
        + [{"label": IMMUNE_UNRESOLVED, "common_resolution": "ambiguous",
            "reason": "stage 1 called immune; stage 2 did not resolve a unique core. Reasons: "
                      "no_naming_evidence (the winning core is corroborated only by SUPPORTING "
                      "genes -- the cell is NOT reassigned to another core), "
                      "core_below_corroboration (< 2 detected genes of the winning core at all), "
                      "core_margin (the winner does not beat the runner-up by 25%)"},
           {"label": NON_IMMUNE, "common_resolution": "ambiguous",
            "reason": "stage 1 called non-immune; no non-immune lineage is resolved here -- "
                      "the ruling is about immune, and exocrine is one label that is not claimed"},
           {"label": AMBIGUOUS, "common_resolution": "ambiguous",
            "reason": "stage 1 did not separate the compartments for this cell"},
           {"label": UNAVAILABLE, "common_resolution": "unavailable",
            "reason": "no nucleus, or the cell is absent from the typed matrix"}],
    }
    (out_dir / "label_space.json").write_text(json.dumps(space, indent=1) + "\n")

    try:
        head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip())
    except Exception:
        head, dirty = "", True
    inputs += [{"kind": "identity_core", "path": c["file"], "sha256": c["sha256"], "label": L}
               for L, c in cores.items()]
    rel = {
        "schema_version": SCHEMA_VERSION,
        "release_id": f"xenium-pln-immune-two-stage-dev-{time.strftime('%Y%m%d', time.gmtime())}",
        "modality": "xenium", "tissue": "pln",
        "producer": {"lane": "xenium-pln-immune-two-stage", "commit": head, "dirty": dirty,
                     "script": "analysis/xenium/scripts/pln_immune_typing.py"},
        "method": {
            "stage_1": "immune vs non-immune from canonical genes; a cell is Immune iff it "
                       "detects >= 2 immune canonicals AND a strictly larger FRACTION of the "
                       "immune canonical set than of the structural set (fractions, not counts, "
                       "so set size carries no advantage under background); Non_immune is the "
                       "mirror; everything else is Ambiguous with the reason recorded",
            "stage_2": "within Immune, argmax of the ABSOLUTE MEAN LOGNORM over each immune "
                       "identity core (the score is over the WHOLE core, supporting genes "
                       "included), admitted only when the winning core has >= 2 detected NAMING "
                       "genes and beats the runner-up by >= 25% of the winner; otherwise "
                       "Immune_unresolved with the reason recorded. A cell that fails only the "
                       "naming floor is NOT reassigned to another core.",
            "naming_vs_supporting": "core/panels' own evidence column decides which genes may "
                                    "name a lineage: evidence in {canonical, both} (note "
                                    "'canonical identity marker') MAY NAME; evidence == "
                                    "cellxgene (note 'marks lineage (CL)') is SUPPORTING ONLY -- "
                                    "it scores and describes, it cannot admit the label. This is "
                                    "lineages.SUPPORTING_ONLY for the Xenium panel. Per-cell "
                                    "columns n_core_genes_detected and n_core_naming_detected "
                                    "carry both counts for the winning core so the floor is "
                                    "auditable from the parquet alone.",
            "detection_predicate": "raw >= 1 AND raw > lambda_c, where lambda_c = "
                                   "(genomic_control_counts_c + a*mu_c)/(n_probes + a) with "
                                   "a = n_probes = 21 and mu_c the cell-area-decile mean per "
                                   "probe (docs/XENIUM_PIPELINE.md section 1): the expected gDNA "
                                   "counts per gene in that cell. Counts are integers, so the "
                                   "effective floor is max(1, floor(lambda_c)+1). There is NO "
                                   "depth term: the floor is not a function of the cell's "
                                   "transcript depth or of the section's median depth. Per-donor "
                                   "lambda quantiles and the number of cells where the floor "
                                   "exceeds 1 are in per_donor[*].detection_floor.",
            "withdrawn_predicate": "This release previously detected at lognorm >= log1p(1) on a "
                                   "per-cell normalisation to the section MEDIAN depth, which is "
                                   "algebraically raw >= depth(c)/median_depth -- a raw-count "
                                   "floor that rises with cell depth, demanding more than one "
                                   "transcript in every cell above median depth. It is "
                                   "withdrawn. per_donor[*].withdrawn_depth_scaled_floor records "
                                   "how many cells of each section it bound.",
            "detection_control": "docs/XENIUM_PIPELINE.md requires a predicate to hold a "
                                 "zero-template set at or below ~0.3%. The control is the "
                                 "non-cross-hybridising Y floor group TBL1Y/USP9Y/EIF1AY, which "
                                 "is zero-template only in a female donor; it is measured per "
                                 "donor under BOTH predicates and written to "
                                 "per_donor[*].background_rates, null (never 0) where the donor "
                                 "is male. DDX3Y/KDM5D are reported separately because their "
                                 "elevation is cross-hybridisation with DDX3X/KDM5C.",
            "reference_levels": "NONE are derived from pLN. In a lymph node the ordinary tissue "
                                "is lymphocytes, so a pLN-derived immune reference level inverts "
                                "the measurement. Stage 1 is a per-cell detection predicate and "
                                "stage 2 a within-cell argmax; no donor is pooled with another.",
            "controls": "The typed matrices carry Gene Expression only (5101 targets, verified "
                        "equal to the bundle's Gene Expression feature set), so no control probe "
                        "was present in them to exclude. The bundle's per-cell control counts are "
                        "carried through as columns; genomic_control_counts and cell_area are "
                        "ALSO an input to the detection predicate, as lambda_c above. "
                        "control_probe_counts is carried as a column only and gates nothing.",
            "nucleus_rule": "nucleus_count == 0 -> Unavailable, excluded from every calculation.",
            "no_spatial_claim": "Pancreas and pLN are different blocks. Nothing here crosses that "
                                "boundary; the comparison this feeds is compositional, on "
                                "donor_id and disease status only.",
        },
        "label_space": "label_space.json", "calls_dir": "calls", "cell_key": "cell_id",
        "roi": sorted({s["roi"] for s in per_donor.values()}),
        "donors": sorted(per_donor), "inputs": inputs,
        "columns": {"label": "lineage_native", "common_lineage": "common_lineage",
                    "evidence_positive": "evidence_positive", "confidence": "confidence",
                    "fine_type": "celltype_fine", "stage1": "stage1"},
        "stamp": "development",
        "latch": "no production typing is declared; this is a development release",
        "constants": {"raw_floor": RAW_FLOOR, "gdna_n_probes": GDNA_N_PROBES,
                      "gdna_area_deciles": GDNA_AREA_DECILES,
                      "anchor_min_genes": ANCHOR_MIN_GENES,
                      "anchor_margin_frac": ANCHOR_MARGIN_FRAC,
                      "naming_evidence": list(NAMING_EVIDENCE)},
        "identity_cores": {L: {"n_genes": len(c["genes"]), "naming": c["naming"],
                               "supporting": c["supporting"]} for L, c in cores.items()},
        "canonical_sets": {"immune": list(IMMUNE_CANONICAL), "structural": list(STRUCTURAL_CANONICAL),
                           "regional_reported_only": list(REGIONAL)},
        "per_donor": per_donor,
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (out_dir / "release.json").write_text(json.dumps(rel, indent=1, default=str) + "\n")
    print(f"[release] wrote {out_dir}  donors={sorted(per_donor)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
