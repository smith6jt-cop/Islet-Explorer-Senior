#!/usr/bin/env python3
"""SUPERSEDED 2026-09-18. ``analysis/xenium/scripts/pln_immune_typing_v3.py`` replaces this
file, and its ``release.json`` declares that ``release_pln_immune_v3`` supersedes the
``release_pln_immune_v2`` tree written here. ``scripts/senior/xenium_release.py`` honours that
declaration, so no builder reads the v2 tree any more -- the pins that named it in
``build_compartment_gradient.py`` and ``modality_concordance.py`` are gone. Recorded as a
verified orphan in ``scripts/retire/dispositions.py`` (code tier 1). Kept as the record of the
cluster-then-name method below; do not read it as the current one.

pLN Xenium immune typing, v2: cluster-then-name.

THE METHOD, AND WHY IT IS NOT THE PANCREAS METHOD
-------------------------------------------------
Maintainer, 2026-09-17: *"All immune cell calls from pancreas can be applied to pLN."*
*"Separate out the immune cells with canonical markers/genes and then use subpanels."*
*"Clustering works well in pLN because of the regional cell types."*  *"tissue does not mean
anything at all."*

In the pancreas a QuPath pixel classifier hands every cell a context -- islet, vessel,
lymphatic, nerve, ordinary tissue -- and a reference level is fitted against that context.
A lymph node section has no such partition: it is one block of node.  So the context machinery
has nothing to attach to, and no "ordinary tissue" reference population exists to fit against.

**Clustering is the substitute for context.**  An unsupervised partition of the section
discovers the regional structure that the pixel classifier would have supplied, and each
cluster is then NAMED by which canonical genes stand out *against the other clusters on that
same section*.  The calibration is entirely internal to one section: no level, no quantile and
no threshold crosses a slide boundary, and nothing computed on pancreas cells reaches a pLN
call.

WHAT IS WITHDRAWN AND NOT REINTRODUCED HERE
-------------------------------------------
(a) Any number derived from pancreas cells -- level, threshold, scale factor, quantile.  This
    file reads no pancreas data of any kind.
(b) `subpanel.fit_levels` / `reference_level` / `context == 'tissue'`.  There is no context
    column on the pLN and no reference population; none is synthesised.
(c) `IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor on pLN values.  It is defined against a
    reference level that does not exist here.

Also not reintroduced: the draft release's detection predicate `raw >= depth / target` with
`target` the section median depth.  That is a raw-count floor that rises with cell depth, in
the tissue where gDNA background is worst.  See DETECTION below.

WHAT THIS DOES NOT LICENSE
--------------------------
No spatial claim crosses the pancreas/pLN boundary -- they are different blocks on one slide,
and nothing here computes a distance, a neighbourhood or a shared embedding.  Non-immune
lineages appear only as far as stage 1 must exclude them.  Acinar and ductal are never claimed;
in a node they would mean nothing in any case.

PIPELINE
--------
1.  CLUSTERS ARE CONSUMED, NEVER RECOMPUTED.  The within-section Leiden partition comes from
    `data/development/typing/pln/clusters/xenium/assignments/` (produced by
    `core/phenocycler/scripts/pln/cluster_pln.py`).  This file re-derives the per-cluster gene
    profile from the h5ads because it needs its own detection predicate (below) -- the
    clustering track's `det__` columns are `lognorm > 0` under a section-median normalisation,
    which is not the predicate declared here.  The partition itself is taken verbatim; the
    cluster ids are the clustering track's own `<section>_cNN`.

2.  STAGE 1, on clusters, within section.  A gene STANDS OUT in cluster k when its
    cluster-level mean lognorm is >= 1 SD above the mean across that section's own clusters
    AND its detection fraction in k exceeds its detection fraction in the rest of the section.
    Both halves are section-internal contrasts.  A cluster is Immune when >= 2 of the canonical
    immune genes stand out and a strictly larger FRACTION of the immune set than of the
    structural set stands out; Non_immune is the mirror; anything else is Uncalled_ambiguous
    with the reason recorded.  Fractions rather than counts, because the two sets differ in
    size and a count comparison would hand the larger set an advantage that is not biology.

3.  STAGE 2, within immune clusters, from the five `core/panels` identity-core subpanels.
    Score of a core = mean standardised enrichment (z) over that core's genes present on the
    panel; the winner is the argmax.  THE SMALL-CORE DEFECT IS THE THING TO FIX: a core with
    few genes can win an argmax on near-zero evidence.  So a winning core is admitted only when
    at least one gene that MAY NAME that lineage (the panel's own `evidence` column, values
    `canonical` or `both`) actually stands out in that cluster.  A cluster that fails this is
    `Immune_uncalled`, named and never reassigned to the runner-up.

4.  CELLS INHERIT THEIR CLUSTER'S CALL, and per-cluster heterogeneity is reported so a mixed
    cluster is visible: the fraction of the cluster's cells whose own per-cell core argmax
    agrees with the cluster's call, the fraction detecting at least one naming gene of the
    called lineage, and the full per-cell argmax distribution.

DETECTION
---------
`lognorm >= log1p(1)` per cell, with lognorm built as `log1p(raw * 1e4 / depth)` -- counts per
10,000, the depth being that cell's own total over all 5,101 genes.  The target is FIXED at
1e4 and is not the section median, so for every cell whose depth is at or below 10,000 counts
the predicate is exactly `raw >= 1`: no depth-scaled floor.  The number of cells whose depth
exceeds 10,000 (where the predicate would demand raw >= 2) is measured and reported per
section, so the claim is auditable rather than asserted.

CONTROLS
--------
*   Y-probe false detection: the panel carries DDX3Y, USP9Y, EIF1AY and KDM5D.  In a female
    donor every detection of those is false, so their detection rate there is a direct read of
    the background this predicate lets through.  Reported per section and per gene.
*   gDNA: `data/development/xenium_gdna/gdna_triage_2026-09-17.parquet` carries a measured pLN
    rate for the serials it covers; donors it does not cover are recorded UNMEASURED, never 0.
*   Control probes: `spatialdata_io` loads Gene Expression only, so the h5ads may carry no
    control features at all.  Whether any `NegControl*` / `Intergenic*` / `Unassigned*` /
    `Deprecated*` feature exists in the var index is MEASURED per section and reported; no
    output claims controls were excluded unless they were there to exclude.

Run:
    .venv/bin/python analysis/xenium/scripts/pln_immune_typing_v2.py --stage call
    .venv/bin/python analysis/xenium/scripts/pln_immune_typing_v2.py --stage compare
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
H5AD_DIR = Path("/data/xenium/h5ad_out")
CLUSTERS = REPO / "data/development/typing/pln/clusters/xenium"
OUT_ROOT = REPO / "data/development/xenium_typing/release_pln_immune_v2"
DRAFT = REPO / "data/development/xenium_typing/release_pln_immune"
IDENTITY_CORE = REPO / "core/panels/data/identity_core"
GDNA = REPO / "data/development/xenium_gdna/gdna_triage_2026-09-17.parquet"
XENIUM_ROOT = Path("/data/xenium")
XENIUM_PATHS = REPO / "core/phenocycler/data/integration/xenium_paths.csv"

SCHEMA_VERSION = 1
RELEASE_ID = "xenium-pln-immune-cluster-named-dev-20260917"
SPACE_ID = "xenium-pln-immune-cluster-named-v2"
LANE = "xenium-pln-immune-cluster-named"
ROI = "pln_1"

TARGET_SUM = 1.0e4
DETECTION_MIN_LOGNORM = float(np.log1p(1.0))
Z_STANDOUT = 1.0
Z_SENSITIVITY = (0.75, 1.0, 1.5)
STAGE1_MIN_STANDOUT = 2
STAGE2_MIN_NAMING_STANDOUT = 1

#: Stage 1 canonicals, exactly the maintainer's list.  CD3D and LYZ were named in the brief and
#: are measured absent from all 24 pLN h5ads; nothing is substituted for them and their absence
#: is recorded in the release.
STAGE1_IMMUNE = [
    "PTPRC", "CD3E", "CD3G", "CD4", "CD8A", "MS4A1", "CD79A", "CD19",
    "CD68", "ITGAX", "ITGAM", "NKG7", "KLRD1",
]
NAMED_IN_BRIEF_CHECK = ["CD3D", "LYZ"]

#: The five immune identity cores, read from the one panel authority.  Label names are kept
#: identical to the draft release's label space so the two are directly comparable.
IMMUNE_CORES = {
    "Immune_T": "03_immune_T_cell.csv",
    "Immune_B_Plasma": "04_immune_B_plasma.csv",
    "Immune_NK_ILC": "05_immune_NK_ILC.csv",
    "Myeloid": "06_immune_myeloid.csv",
    "Granulocyte": "07_immune_granulocyte.csv",
}
NON_IMMUNE_CORES = {
    "endothelial_vascular": "11_endothelial_vascular.csv",
    "pericyte_smooth_muscle": "12_pericyte_smooth_muscle.csv",
    "neural_glial": "15_neural_glial.csv",
    "epithelial_general": "16_epithelial_general.csv",
}
COMMON_LINEAGE = {
    "Immune_T": "T_NK", "Immune_NK_ILC": "T_NK", "Immune_B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid", "Granulocyte": "Granulocyte",
}
#: `evidence` values that MAY NAME a lineage, from the panel's own column.  `cellxgene`
#: ("marks lineage (CL)") scores and describes but cannot admit a label.
NAMING_EVIDENCE = frozenset({"canonical", "both"})

#: Y probes.  TBL1Y / USP9Y / EIF1AY are the FLOOR GROUP: zero-template in a female donor and,
#: unlike DDX3Y and KDM5D, they do not cross-hybridise with a highly expressed X paralog
#: (DDX3X, KDM5C).  The floor group is the control; DDX3Y and KDM5D are reported beside it and
#: are NOT the control.
Y_FLOOR = ["TBL1Y", "USP9Y", "EIF1AY"]
Y_CROSSHYB = ["DDX3Y", "KDM5D"]
Y_GENES = Y_FLOOR + Y_CROSSHYB
REGIONAL = ["CR2", "CXCL13", "CCL19", "CCL21", "BCL6", "AICDA", "SELL", "CD27", "MKI67"]
CONTROL_PREFIXES = ("NegControl", "Intergenic", "Unassigned", "Deprecated", "BLANK", "antisense")


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256(path: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def read_core(name: str) -> tuple[list[str], list[str]]:
    frame = pd.read_csv(IDENTITY_CORE / name)
    genes = [str(g) for g in frame["gene"].dropna().unique()]
    naming = [
        str(r.gene) for r in frame.itertuples()
        if str(getattr(r, "evidence", "")).strip() in NAMING_EVIDENCE
    ]
    return genes, sorted(set(naming))


def load_panels() -> dict:
    immune = {k: read_core(v) for k, v in IMMUNE_CORES.items()}
    non_immune = {k: read_core(v) for k, v in NON_IMMUNE_CORES.items()}
    structural_naming = sorted({g for _, naming in non_immune.values() for g in naming})
    wanted = sorted(
        set(STAGE1_IMMUNE) | set(Y_GENES) | set(REGIONAL) | set(structural_naming)
        | {g for genes, _ in immune.values() for g in genes}
        | {g for genes, _ in non_immune.values() for g in genes}
    )
    return {
        "immune": immune,
        "non_immune": non_immune,
        "structural_naming": structural_naming,
        "wanted": wanted,
    }


def sections() -> list[tuple[str, str, Path]]:
    """(donor_id, section_key, h5ad) for every section the clustering track partitioned.

    The donor set is whatever `assignments/` holds; the central cohort exclusion was already
    applied there by `phenocycler.cohort.filter_eligible_donors` and is not restated here.
    """
    out = []
    for path in sorted(CLUSTERS.glob("assignments/donor_id=*/section=*/clusters.parquet")):
        donor = path.parent.parent.name.split("=", 1)[1]
        section = path.parent.name.split("=", 1)[1]
        h5ad = H5AD_DIR / f"{section}_pLNXenium.h5ad"
        if not h5ad.exists():
            raise SystemExit(f"missing h5ad for section {section}: {h5ad}")
        out.append((donor, section, h5ad))
    return out


def nucleus_for(donor: str, cell_ids: np.ndarray) -> tuple[np.ndarray | None, str]:
    """`nucleus_count` per h5ad cell, from the vendor bundle, or (None, why-not).

    A cell with no nucleus is excluded from every calculation -- CLAUDE.md, and the rule does
    not bend for a modality.  The h5ads carry no nucleus geometry, so the only handle is the
    bundle `cells.parquet`.  Where the bundle is not on this machine the rule is UNMEASURABLE
    and that is recorded as such; it is never recorded as "passed", and no cell is treated as
    nucleated by assumption without saying so.

    The bundle is matched by donor AND by subset containment of the h5ad's own cell ids --
    never by filename order.
    """
    table = pd.read_csv(XENIUM_PATHS, dtype=str)
    table = table[
        table["roi"].astype(str).str.startswith("pln")
        & (table["donor_id"].astype(str) == str(donor))
    ]
    tried = []
    for row in table.itertuples():
        src = str(row.source_path or "").strip()
        if not src or src == "NA":
            continue
        cells = XENIUM_ROOT / src / "cells.parquet"
        if not cells.is_file():
            tried.append(f"{row.roi}:absent")
            continue
        frame = pd.read_parquet(cells, columns=["cell_id", "nucleus_count"])
        frame["cell_id"] = frame["cell_id"].astype(str)
        lookup = dict(zip(frame["cell_id"], frame["nucleus_count"].to_numpy()))
        missing = int(sum(1 for c in cell_ids if c not in lookup))
        if missing:
            tried.append(f"{row.roi}:not_a_superset({missing} ids missing)")
            continue
        counts = np.array([lookup[c] for c in cell_ids], dtype=np.int32)
        return counts, f"bundle {row.roi} {cells}"
    return None, "UNMEASURED: no local bundle cells.parquet is a superset of this h5ad" + (
        f" (tried {'; '.join(tried)})" if tried else ""
    )


def stage1_decide(n_imm, f_imm, n_str, f_str):
    if n_imm >= STAGE1_MIN_STANDOUT and f_imm > f_str:
        return "Immune", ""
    if n_str >= STAGE1_MIN_STANDOUT and f_str > f_imm:
        return "Non_immune", ""
    if n_imm < STAGE1_MIN_STANDOUT and n_str < STAGE1_MIN_STANDOUT:
        return "Uncalled_ambiguous", "no_set_reaches_two_standouts"
    if f_imm == f_str:
        return "Uncalled_ambiguous", "equal_standout_fractions"
    return "Uncalled_ambiguous", "standout_count_and_fraction_disagree"


def process_section(donor: str, section: str, h5ad: Path, panels: dict) -> dict:
    import anndata as ad

    t0 = time.time()
    adata = ad.read_h5ad(h5ad)
    cell_ids = adata.obs_names.to_numpy().astype(str)
    var_names = adata.var_names.to_numpy().astype(str)
    n_cells = int(adata.n_obs)

    controls_present = sorted(
        g for g in var_names if any(g.startswith(p) for p in CONTROL_PREFIXES)
    )
    named_absent = [g for g in NAMED_IN_BRIEF_CHECK if g not in set(var_names)]

    X = adata.X
    X = X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
    depth = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    obs_depth = adata.obs["nCount_Xenium"].to_numpy(dtype=float)
    nfeature = adata.obs["nFeature_Xenium"].to_numpy(dtype=float)
    depth_matches_obs = int(np.sum(np.abs(depth - obs_depth) < 0.5))

    gene_idx = {g: i for i, g in enumerate(var_names)}
    present = [g for g in panels["wanted"] if g in gene_idx]
    cols = np.array([gene_idx[g] for g in present])
    sub = X[:, cols].tocsr()
    del X
    gc.collect()

    safe_depth = np.where(depth > 0, depth, np.nan)
    rows = np.repeat(np.arange(sub.shape[0]), np.diff(sub.indptr))
    scaled = sub.data.astype(np.float64) * TARGET_SUM / safe_depth[rows]
    scaled = np.nan_to_num(scaled, nan=0.0)
    lognorm = np.log1p(scaled).astype(np.float32)
    detected = (lognorm >= DETECTION_MIN_LOGNORM).astype(np.float32)

    # how often the fixed-target predicate is stricter than raw >= 1
    deep = int(np.sum(depth > TARGET_SUM))
    raw_one = int(sub.data.size)
    det_total = int(detected.sum())

    L = sp.csr_matrix((lognorm, sub.indices, sub.indptr), shape=sub.shape)
    D = sp.csr_matrix((detected, sub.indices, sub.indptr), shape=sub.shape)
    del sub, scaled, lognorm, detected, rows
    gc.collect()

    assign = pd.read_parquet(
        CLUSTERS / f"assignments/donor_id={donor}/section={section}/clusters.parquet"
    )
    amap = dict(zip(assign["cell_id"].astype(str), assign["cluster_id"].astype(str)))
    cluster_of = np.array([amap.get(c, "") for c in cell_ids])
    in_cluster = cluster_of != ""
    nucleus_count, nucleus_source = nucleus_for(donor, cell_ids)
    if nucleus_count is None:
        has_nucleus = np.ones(n_cells, dtype=bool)
        n_no_nucleus = 0
        nucleus_measured = False
    else:
        has_nucleus = nucleus_count > 0
        n_no_nucleus = int((~has_nucleus).sum())
        nucleus_measured = True
    # a cell with no nucleus is excluded from every calculation: it leaves the cluster
    # aggregates as well as the calls
    clustered = in_cluster & has_nucleus
    order = sorted(set(cluster_of[clustered]))
    lookup = {c: i for i, c in enumerate(order)}
    krow = np.array([lookup[c] for c in cluster_of[clustered]])
    cell_pos = np.flatnonzero(clustered)
    indic = sp.csr_matrix(
        (np.ones(krow.size, dtype=np.float32), (krow, cell_pos)),
        shape=(len(order), n_cells),
    )
    sizes = np.asarray(indic.sum(axis=1)).ravel()

    means = np.asarray((indic @ L).todense()) / sizes[:, None]
    dets = np.asarray((indic @ D).todense()) / sizes[:, None]
    n_clustered = int(sizes.sum())
    det_tot = np.asarray(D[cell_pos].sum(axis=0)).ravel()
    det_counts = dets * sizes[:, None]
    other_n = n_clustered - sizes
    with np.errstate(invalid="ignore", divide="ignore"):
        d_other = (det_tot[None, :] - det_counts) / np.where(other_n > 0, other_n, np.nan)[:, None]
    d_other = np.nan_to_num(d_other, nan=0.0)

    mu = means.mean(axis=0)
    sd = means.std(axis=0, ddof=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (means - mu[None, :]) / np.where(sd > 0, sd, np.nan)[None, :]
    z = np.nan_to_num(z, nan=0.0)
    standout = {thr: (z >= thr) & (dets > d_other) for thr in Z_SENSITIVITY}

    pidx = {g: i for i, g in enumerate(present)}
    imm_cols = np.array([pidx[g] for g in STAGE1_IMMUNE if g in pidx])
    str_cols = np.array([pidx[g] for g in panels["structural_naming"] if g in pidx])
    core_cols = {
        name: np.array([pidx[g] for g in genes if g in pidx])
        for name, (genes, _) in panels["immune"].items()
    }
    core_naming_cols = {
        name: np.array([pidx[g] for g in naming if g in pidx])
        for name, (_, naming) in panels["immune"].items()
    }
    core_names = list(IMMUNE_CORES)

    # ---- cluster-level calls, at every sensitivity threshold -------------------------------
    per_thr = {}
    for thr in Z_SENSITIVITY:
        so = standout[thr]
        n_imm = so[:, imm_cols].sum(axis=1)
        n_str = so[:, str_cols].sum(axis=1)
        f_imm = n_imm / max(len(imm_cols), 1)
        f_str = n_str / max(len(str_cols), 1)
        s1 = []
        s1r = []
        for k in range(len(order)):
            call, reason = stage1_decide(int(n_imm[k]), float(f_imm[k]), int(n_str[k]), float(f_str[k]))
            s1.append(call)
            s1r.append(reason)
        scores = np.column_stack([
            z[:, core_cols[c]].mean(axis=1) if core_cols[c].size else np.full(len(order), -np.inf)
            for c in core_names
        ])
        naming_so = np.column_stack([
            so[:, core_naming_cols[c]].sum(axis=1) if core_naming_cols[c].size
            else np.zeros(len(order), dtype=int)
            for c in core_names
        ])
        win = scores.argmax(axis=1)
        srt = np.sort(scores, axis=1)
        margin = srt[:, -1] - srt[:, -2]
        s2 = []
        s2r = []
        for k in range(len(order)):
            if s1[k] != "Immune":
                s2.append(s1[k])
                s2r.append(s1r[k])
                continue
            w = core_names[win[k]]
            if naming_so[k, win[k]] >= STAGE2_MIN_NAMING_STANDOUT:
                s2.append(w)
                s2r.append("")
            else:
                s2.append("Immune_uncalled")
                s2r.append(f"winning_core_{w}_has_zero_standout_naming_genes")
        per_thr[thr] = {
            "stage1": s1, "stage1_reason": s1r, "stage2": s2, "stage2_reason": s2r,
            "n_imm": n_imm, "n_str": n_str, "f_imm": f_imm, "f_str": f_str,
            "scores": scores, "naming_so": naming_so, "win": win, "margin": margin,
        }

    main = per_thr[Z_STANDOUT]

    # ---- per-cell descriptive scores, for heterogeneity ------------------------------------
    gsel = np.zeros((len(present), len(core_names)), dtype=np.float32)
    nsel = np.zeros((len(present), len(core_names)), dtype=np.float32)
    for j, c in enumerate(core_names):
        if core_cols[c].size:
            gsel[core_cols[c], j] = 1.0 / core_cols[c].size
        if core_naming_cols[c].size:
            nsel[core_naming_cols[c], j] = 1.0
    cell_core_mean = np.asarray(L @ gsel)
    cell_naming_det = np.asarray(D @ nsel)
    cell_argmax = cell_core_mean.argmax(axis=1)
    cell_argmax_flat = cell_core_mean.max(axis=1) <= 0
    cell_argmax_label = np.array(core_names, dtype=object)[cell_argmax]
    cell_argmax_label[cell_argmax_flat] = ""

    n_imm_cell = np.asarray(D[:, imm_cols].sum(axis=1)).ravel().astype(np.int16)
    n_str_cell = np.asarray(D[:, str_cols].sum(axis=1)).ravel().astype(np.int16)
    ptprc = (
        np.asarray(D[:, [pidx["PTPRC"]]].todense()).ravel() > 0 if "PTPRC" in pidx
        else np.zeros(n_cells, dtype=bool)
    )
    y_cols = np.array([pidx[g] for g in Y_GENES if g in pidx])
    y_det = np.asarray(D[:, y_cols].sum(axis=1)).ravel() if y_cols.size else np.zeros(n_cells)
    floor_cols = np.array([pidx[g] for g in Y_FLOOR if g in pidx])
    y_floor_det = (
        np.asarray(D[cell_pos][:, floor_cols].sum(axis=1)).ravel() if floor_cols.size
        else np.zeros(len(cell_pos))
    )
    y_per_gene = {
        g: float(np.asarray(D[cell_pos][:, [pidx[g]]].todense()).ravel().mean())
        for g in Y_GENES if g in pidx
    }

    # ---- per-cell frame ---------------------------------------------------------------------
    cl_stage1 = np.array(main["stage1"], dtype=object)
    cl_stage2 = np.array(main["stage2"], dtype=object)
    cl_reason = np.array(main["stage2_reason"], dtype=object)
    cidx = np.full(n_cells, -1, dtype=int)
    cidx[cell_pos] = krow
    fallback_label = np.where(has_nucleus, "Uncalled_not_clustered", "Unavailable_no_nucleus")
    fallback_reason = np.where(
        has_nucleus, "below_cluster_qc_nCount_lt_10_or_nFeature_lt_5", "no_nucleus"
    )
    stage1_cell = np.where(cidx >= 0, cl_stage1[np.clip(cidx, 0, None)], fallback_label)
    stage2_cell = np.where(cidx >= 0, cl_stage2[np.clip(cidx, 0, None)], fallback_label)
    reason_cell = np.where(cidx >= 0, cl_reason[np.clip(cidx, 0, None)], fallback_reason)
    call_core_j = np.array([core_names.index(x) if x in core_names else -1 for x in cl_stage2])
    per_cell_core_j = np.where(cidx >= 0, call_core_j[np.clip(cidx, 0, None)], -1)
    naming_of_call = np.where(
        per_cell_core_j >= 0,
        cell_naming_det[np.arange(n_cells), np.clip(per_cell_core_j, 0, None)],
        np.nan,
    )
    agrees = np.where(
        per_cell_core_j >= 0, cell_argmax == np.clip(per_cell_core_j, 0, None), False
    ) & (cell_argmax_label != "")

    calls = pd.DataFrame({
        "cell_id": cell_ids,
        "donor_id": donor,
        "section_key": section,
        "roi": ROI,
        "cluster_id": np.where(cidx >= 0, cluster_of, None),
        "stage1": stage1_cell.astype(str),
        "lineage": stage2_cell.astype(str),
        "reason": reason_cell.astype(str),
        "common_lineage": [COMMON_LINEAGE.get(x, "") for x in stage2_cell],
        "cluster_inherited": cidx >= 0,
        "n_stage1_immune_detected": n_imm_cell,
        "n_structural_naming_detected": n_str_cell,
        "ptprc_detected": ptprc,
        "n_naming_detected_for_cluster_call": naming_of_call.astype("float32"),
        "cell_core_argmax": cell_argmax_label.astype(str),
        "cell_agrees_with_cluster_call": agrees,
        "y_probe_detected": (y_det > 0),
        "nucleus_count": (nucleus_count if nucleus_count is not None
                          else np.full(n_cells, -1, dtype=np.int32)),
        "nucleus_measured": nucleus_measured,
        "depth": depth.astype(np.int32),
        "n_genes_detected_in_wanted_set": np.asarray(D.sum(axis=1)).ravel().astype(np.int32),
    })

    # ---- cluster-level frame -----------------------------------------------------------------
    het_rows = []
    for k, cid in enumerate(order):
        sel = cell_pos[krow == k]
        lab = cell_argmax_label[sel]
        dist = pd.Series(lab).value_counts(normalize=True)
        called = cl_stage2[k]
        het_rows.append({
            "donor_id": donor,
            "section_key": section,
            "cluster_id": cid,
            "n_cells": int(sizes[k]),
            "frac_of_section": float(sizes[k] / n_clustered),
            "median_depth": float(np.median(depth[sel])),
            "median_n_features": float(np.median(nfeature[sel])),
            "stage1": cl_stage1[k],
            "stage1_reason": main["stage1_reason"][k],
            "lineage": called,
            "stage2_reason": cl_reason[k],
            "n_immune_canonical_standout": int(main["n_imm"][k]),
            "n_structural_naming_standout": int(main["n_str"][k]),
            "frac_immune_canonical_standout": float(main["f_imm"][k]),
            "frac_structural_naming_standout": float(main["f_str"][k]),
            "stage1_evidence_margin": float(main["f_imm"][k] - main["f_str"][k]),
            "stage1_thin_two_standouts": bool(
                (cl_stage1[k] == "Immune" and int(main["n_imm"][k]) == STAGE1_MIN_STANDOUT)
                or (cl_stage1[k] == "Non_immune" and int(main["n_str"][k]) == STAGE1_MIN_STANDOUT)
            ),
            "winning_core": core_names[int(main["win"][k])],
            "winning_core_score": float(main["scores"][k, int(main["win"][k])]),
            "winning_core_naming_standouts": int(main["naming_so"][k, int(main["win"][k])]),
            "score_margin_over_runner_up": float(main["margin"][k]),
            **{f"score__{c}": float(main["scores"][k, j]) for j, c in enumerate(core_names)},
            **{f"naming_standouts__{c}": int(main["naming_so"][k, j])
               for j, c in enumerate(core_names)},
            "purity_cell_argmax_agrees": (
                float(np.mean(cell_argmax_label[sel] == called)) if called in core_names else np.nan
            ),
            "frac_cells_with_naming_gene_of_call": (
                float(np.mean(cell_naming_det[sel, core_names.index(called)] > 0))
                if called in core_names else np.nan
            ),
            "cell_argmax_top_label": dist.index[0] if len(dist) else "",
            "cell_argmax_top_frac": float(dist.iloc[0]) if len(dist) else np.nan,
            "cell_argmax_entropy": float(
                -np.sum(dist.values * np.log(dist.values + 1e-12)) if len(dist) else np.nan
            ),
            "standout_immune_genes": ";".join(
                g for g in STAGE1_IMMUNE if g in pidx and standout[Z_STANDOUT][k, pidx[g]]
            ),
            "standout_structural_genes": ";".join(
                g for g in panels["structural_naming"]
                if g in pidx and standout[Z_STANDOUT][k, pidx[g]]
            ),
        })

    # ---- per-cluster gene statistics, for audit ---------------------------------------------
    gene_stats = pd.DataFrame({
        "donor_id": donor,
        "section_key": section,
        "cluster_id": np.repeat(order, len(present)),
        "gene": np.tile(present, len(order)),
        "mean_lognorm": means.ravel().astype(np.float32),
        "z_across_clusters": z.ravel().astype(np.float32),
        "det_frac": dets.ravel().astype(np.float32),
        "det_frac_rest_of_section": d_other.ravel().astype(np.float32),
        "standout": standout[Z_STANDOUT].ravel(),
    })

    sens = {}
    for thr in Z_SENSITIVITY:
        s2 = pd.Series(per_thr[thr]["stage2"])
        sens[str(thr)] = {
            "clusters": {k: int(v) for k, v in s2.value_counts().items()},
            "cells": {
                k: int(sizes[[i for i, x in enumerate(per_thr[thr]["stage2"]) if x == k]].sum())
                for k in s2.unique()
            },
        }

    report = {
        "donor_id": donor,
        "section_key": section,
        "h5ad": str(h5ad),
        "n_cells_in_h5ad": n_cells,
        "n_cells_clustered": n_clustered,
        "n_cells_not_clustered": n_cells - n_clustered,
        "nucleus_source": nucleus_source,
        "nucleus_rule_measurable": nucleus_measured,
        "n_cells_no_nucleus_excluded": n_no_nucleus,
        "n_cells_below_cluster_qc": int((~in_cluster & has_nucleus).sum()),
        "n_clusters": len(order),
        "n_genes_on_panel": int(len(var_names)),
        "n_wanted_genes_present": len(present),
        "wanted_genes_absent": [g for g in panels["wanted"] if g not in gene_idx],
        "named_in_brief_but_absent": named_absent,
        "control_features_present_in_var": controls_present,
        "depth_equals_obs_nCount_for_cells": depth_matches_obs,
        "cells_with_depth_gt_target": deep,
        "nonzero_gene_entries_in_wanted_set": raw_one,
        "detected_entries_in_wanted_set": det_total,
        "entries_dropped_by_fixed_target": raw_one - det_total,
        "sex": str(adata.obs["sex"].iloc[0]),
        "disease_status": str(adata.obs["disease.status"].iloc[0]),
        "y_probe_detect_frac_any": float(np.mean(y_det[cell_pos] > 0)),
        "y_floor_group": Y_FLOOR,
        "y_floor_detect_frac_any": float(np.mean(y_floor_det > 0)),
        "y_probe_detect_frac_per_gene": y_per_gene,
        "sensitivity_z": sens,
        "seconds": round(time.time() - t0, 1),
    }
    del L, D, adata
    gc.collect()
    return {
        "calls": calls,
        "clusters": pd.DataFrame(het_rows),
        "gene_stats": gene_stats,
        "report": report,
    }


def run_calls(args) -> None:
    panels = load_panels()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUT_ROOT / "clusters").mkdir(exist_ok=True)
    secs = sections()
    if args.only:
        keep = set(args.only.split(","))
        secs = [s for s in secs if s[1] in keep or s[0] in keep]
    log(f"{len(secs)} sections, {len(set(d for d, _, _ in secs))} donors")

    by_donor: dict[str, list[pd.DataFrame]] = {}
    clusters, gene_stats, reports, inputs = [], [], [], {}
    for donor, section, h5ad in secs:
        log(f"section {section} (donor {donor})")
        res = process_section(donor, section, h5ad, panels)
        by_donor.setdefault(donor, []).append(res["calls"])
        clusters.append(res["clusters"])
        gene_stats.append(res["gene_stats"])
        reports.append(res["report"])
        inputs[str(h5ad)] = sha256(h5ad)
        inputs[str(CLUSTERS / f"assignments/donor_id={donor}/section={section}/clusters.parquet")] = (
            sha256(CLUSTERS / f"assignments/donor_id={donor}/section={section}/clusters.parquet")
        )
        r = res["report"]
        log(f"  {r['n_cells_clustered']:,} clustered / {r['n_cells_in_h5ad']:,}; "
            f"{r['n_clusters']} clusters; {r['seconds']}s")
        gc.collect()

    dup_report = {}
    for donor, frames in by_donor.items():
        frame = pd.concat(frames, ignore_index=True)
        dup = int(frame.duplicated(subset=["cell_id"]).sum())
        if dup:
            dup_report[donor] = dup
        dest = OUT_ROOT / f"calls/donor_id={donor}/roi={ROI}"
        dest.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(dest / "data.parquet", index=False)
        log(f"wrote {dest/'data.parquet'} ({len(frame):,} rows)")

    cl = pd.concat(clusters, ignore_index=True)
    cl.to_parquet(OUT_ROOT / "clusters/cluster_calls.parquet", index=False)
    cl.to_csv(OUT_ROOT / "clusters/cluster_calls.csv", index=False)
    pd.concat(gene_stats, ignore_index=True).to_parquet(
        OUT_ROOT / "clusters/cluster_gene_stats.parquet", index=False
    )

    # per-donor summary, from what was written
    rows = []
    for donor in sorted(by_donor):
        frame = pd.read_parquet(OUT_ROOT / f"calls/donor_id={donor}/roi={ROI}/data.parquet",
                                columns=["stage1", "lineage"])
        s1 = frame["stage1"].value_counts()
        s2 = frame["lineage"].value_counts()
        rows.append({
            "donor_id": donor,
            "cells": len(frame),
            "immune": int(s1.get("Immune", 0)),
            "non_immune": int(s1.get("Non_immune", 0)),
            "uncalled_ambiguous": int(s1.get("Uncalled_ambiguous", 0)),
            "not_clustered": int(s1.get("Uncalled_not_clustered", 0)),
            "no_nucleus_excluded": int(s1.get("Unavailable_no_nucleus", 0)),
            **{f"s2__{k}": int(s2.get(k, 0)) for k in
               list(IMMUNE_CORES) + ["Immune_uncalled", "Non_immune", "Uncalled_ambiguous",
                                     "Uncalled_not_clustered", "Unavailable_no_nucleus"]},
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_ROOT / "per_donor_summary.csv", index=False)

    gdna = {}
    if GDNA.exists():
        g = pd.read_parquet(GDNA)
        g = g[g["tissue"] == "pln"]
        gdna = {str(r.donor): float(r.rate) for r in g.itertuples()}
    donors = sorted(by_donor)
    gdna_out = {d: gdna.get(d, "UNMEASURED") for d in donors}

    label_space = {
        "schema_version": SCHEMA_VERSION,
        "space_id": SPACE_ID,
        "producer_lane": LANE,
        "labels": (
            [{"label": k, "common_lineage": COMMON_LINEAGE[k],
              "reason": "stage 2: identity-core argmax over the cluster's standardised gene "
                        "enrichment, admitted only because at least one gene that MAY NAME this "
                        "lineage (core/panels evidence in {canonical, both}) stands out in that "
                        "cluster against the section's other clusters"}
             for k in IMMUNE_CORES]
            + [
                {"label": "Immune_uncalled", "common_lineage": "",
                 "reason": "stage 1 called the cluster Immune; the winning identity core had zero "
                           "standout naming genes of its own lineage, so the call is withheld and "
                           "the cluster is NOT reassigned to the runner-up"},
                {"label": "Non_immune", "common_lineage": "",
                 "reason": "stage 1: >= 2 structural naming genes stand out and a larger fraction "
                           "of the structural set than of the immune canonical set"},
                {"label": "Uncalled_ambiguous", "common_lineage": "",
                 "reason": "stage 1 reached neither side; the specific reason is in the per-cell "
                           "and per-cluster `reason` column"},
                {"label": "Unavailable_no_nucleus", "common_lineage": "",
                 "reason": "the bundle records nucleus_count == 0 for this cell. A cell with no "
                           "nucleus is excluded from every calculation, so it is out of the "
                           "cluster aggregates as well as out of the calls. Where the bundle is "
                           "not on this machine the rule is UNMEASURABLE and recorded as such, "
                           "never as passed"},
                {"label": "Uncalled_not_clustered", "common_lineage": "",
                 "reason": "the cell is in the h5ad but below the clustering track's QC "
                           "(nCount_Xenium < 10 or nFeature_Xenium < 5), so it has no cluster to "
                           "inherit from; it is carried in the release as uncalled, never as 0"},
            ]
        ),
        "not_a_label": {
            "acinar_ductal": "never claimed anywhere; exocrine is one label in this project and "
                             "has no meaning in a lymph node",
        },
    }
    (OUT_ROOT / "label_space.json").write_text(json.dumps(label_space, indent=1))

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                            text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=REPO, capture_output=True,
                                text=True).stdout.strip())
    release = {
        "schema_version": SCHEMA_VERSION,
        "release_id": RELEASE_ID,
        "modality": "xenium",
        "tissue": "pln",
        "roi": ROI,
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "latch": "development; nothing here is a release of record",
        "producer": {
            "lane": LANE,
            "script": "analysis/xenium/scripts/pln_immune_typing_v2.py",
            "commit": commit,
            "dirty": dirty,
            "supersedes": "data/development/xenium_typing/release_pln_immune (draft; its stage-2 "
                          "small-core labels were refuted by adversarial verification and its "
                          "detection predicate was a depth-scaled raw floor). The draft is left "
                          "in place for comparison and is not overwritten.",
        },
        "method": {
            "context": "A pLN section carries no context partition -- no islet/vessel/tissue "
                       "objects, no reference population. Clustering is the substitute: the "
                       "within-section Leiden partition is CONSUMED from "
                       "data/development/typing/pln/clusters/xenium/assignments (produced by "
                       "core/phenocycler/scripts/pln/cluster_pln.py). Nothing here re-clusters.",
            "normalisation": f"log1p(raw * {TARGET_SUM:.0f} / depth), depth = that cell's total "
                             "raw count over all panel genes. The target is FIXED, not the "
                             "section median.",
            "detection": f"lognorm >= log1p(1), i.e. raw * {TARGET_SUM:.0f} / depth >= 1. For a "
                         "cell of depth <= 10000 this is exactly raw >= 1, so no depth-scaled raw "
                         "floor operates; the count of cells above that depth is reported per "
                         "section in section_reports.json as cells_with_depth_gt_target.",
            "standout": f"a gene stands out in cluster k iff its cluster-mean lognorm is >= "
                        f"{Z_STANDOUT} SD above the mean over that section's own clusters AND its "
                        "detection fraction in k exceeds its detection fraction in the rest of "
                        "the section. Both contrasts are internal to one section; no level, "
                        "quantile or threshold is imported from any other slide or tissue.",
            "stage_1": f"cluster is Immune iff >= {STAGE1_MIN_STANDOUT} of the canonical immune "
                       "genes stand out AND a strictly larger FRACTION of the immune canonical "
                       "set than of the structural naming set stands out; Non_immune is the "
                       "mirror; otherwise Uncalled_ambiguous with the reason recorded.",
            "stage_2": f"within Immune clusters, argmax over the five core/panels immune identity "
                       "cores of the mean standardised enrichment (z) over that core's genes, "
                       f"admitted only when the winning core has >= {STAGE2_MIN_NAMING_STANDOUT} "
                       "standout gene that MAY NAME its lineage. A cluster failing that is "
                       "Immune_uncalled and is NEVER reassigned to the runner-up. This is the "
                       "small-core fix: a core with few genes cannot win on near-zero evidence.",
            "naming_vs_supporting": "core/panels' own evidence column decides: evidence in "
                                    "{canonical, both} MAY NAME; evidence == cellxgene "
                                    "('marks lineage (CL)') scores and describes only.",
            "inheritance": "cells inherit their cluster's call verbatim. Per-cluster "
                           "heterogeneity (purity of the per-cell core argmax, fraction of cells "
                           "detecting a naming gene of the call, argmax entropy) is in "
                           "clusters/cluster_calls.parquet so a mixed cluster is visible rather "
                           "than hidden behind its majority label.",
            "withdrawn_and_not_reintroduced": [
                "any number derived from pancreas cells (level, threshold, scale factor, quantile)",
                "subpanel.fit_levels / reference_level / context == 'tissue' on the pLN",
                "IMMUNE_REQUIRE_MIN = 0.10 as an absolute floor on pLN values",
                "the draft's detection predicate raw >= depth / section_median_depth",
            ],
            "not_licensed": [
                "any spatial claim across the pancreas/pLN boundary",
                "non-immune lineages beyond what stage 1 must exclude",
                "acinar or ductal, ever",
            ],
        },
        "panel": {
            "stage1_immune_canonicals": STAGE1_IMMUNE,
            "named_in_brief_but_absent_from_every_section": sorted(
                {g for r in reports for g in r["named_in_brief_but_absent"]}
            ),
            "structural_naming_genes": panels["structural_naming"],
            "immune_cores": {k: {"file": f"core/panels/data/identity_core/{v}",
                                 "n_genes": len(panels["immune"][k][0]),
                                 "n_naming_genes": len(panels["immune"][k][1])}
                             for k, v in IMMUNE_CORES.items()},
            "control_features_present_in_h5ad_var": sorted(
                {g for r in reports for g in r["control_features_present_in_var"]}
            ) or "none: these h5ads carry Gene Expression only, so no control probe was there to "
                 "exclude and no output claims one was",
        },
        "controls": {
            "y_probe_genes": Y_GENES,
            "y_floor_group": Y_FLOOR,
            "y_cross_hybridising_reported_not_control": Y_CROSSHYB,
            "y_probe_note": "in a female donor every Y detection is false; that rate is the "
                            "false-detection rate this predicate admits. Per section in "
                            "section_reports.json and rolled up in controls.json.",
            "gdna_rate_per_donor": gdna_out,
            "gdna_source": str(GDNA.relative_to(REPO)) if GDNA.exists() else "UNMEASURED",
        },
        "nucleus_rule": {
            "rule": "a cell with no nucleus is excluded from every calculation",
            "handle": "the vendor bundle cells.parquet nucleus_count; the h5ads carry no nucleus "
                      "geometry, so where the bundle is not on this machine the rule is "
                      "UNMEASURABLE and is recorded as such, never as passed",
            "sections_measured": sorted(r["section_key"] for r in reports
                                        if r["nucleus_rule_measurable"]),
            "sections_unmeasured": sorted(r["section_key"] for r in reports
                                          if not r["nucleus_rule_measurable"]),
            "cells_excluded_no_nucleus": int(sum(r["n_cells_no_nucleus_excluded"]
                                                 for r in reports)),
        },
        "donors": donors,
        "sections": {r["section_key"]: r["donor_id"] for r in reports},
        "cells": int(summary["cells"].sum()),
        "cells_clustered": int(sum(r["n_cells_clustered"] for r in reports)),
        "duplicate_cell_ids_within_donor": dup_report or "none",
        "inputs_sha256": inputs,
    }
    (OUT_ROOT / "release.json").write_text(json.dumps(release, indent=1))
    (OUT_ROOT / "section_reports.json").write_text(json.dumps(reports, indent=1))

    ctrl = pd.DataFrame([{
        "section_key": r["section_key"], "donor_id": r["donor_id"], "sex": r["sex"],
        "disease_status": r["disease_status"],
        "n_cells_clustered": r["n_cells_clustered"],
        "y_any_detect_frac": r["y_probe_detect_frac_any"],
        "y_floor_group_detect_frac": r["y_floor_detect_frac_any"],
        **{f"y__{g}": v for g, v in r["y_probe_detect_frac_per_gene"].items()},
        "cells_with_depth_gt_target": r["cells_with_depth_gt_target"],
        "gdna_pln_rate": gdna.get(r["donor_id"], np.nan),
    } for r in reports])
    ctrl.to_csv(OUT_ROOT / "controls.csv", index=False)
    log(f"wrote {OUT_ROOT/'release.json'}")
    print(summary.to_string(index=False))


def run_compare(args) -> None:
    """Cross-tab this release against the draft, on the cells both cover."""
    out = {"draft": str(DRAFT.relative_to(REPO)), "v2": str(OUT_ROOT.relative_to(REPO))}
    rows = []
    for path in sorted(OUT_ROOT.glob("calls/donor_id=*/roi=*/data.parquet")):
        donor = path.parent.parent.name.split("=", 1)[1]
        dpath = DRAFT / f"calls/donor_id={donor}/roi={ROI}/data.parquet"
        if not dpath.exists():
            rows.append({"donor_id": donor, "in_draft": False})
            continue
        mine = pd.read_parquet(path, columns=["cell_id", "stage1", "lineage"])
        draft = pd.read_parquet(dpath, columns=["cell_id", "stage1", "lineage_native"])
        j = mine.merge(draft, on="cell_id", how="inner", suffixes=("_v2", "_draft"))
        rows.append({
            "donor_id": donor, "in_draft": True,
            "v2_cells": len(mine), "draft_cells": len(draft), "shared_cells": len(j),
            "stage1_agree": int((j["stage1_v2"] == j["stage1_draft"]).sum()),
            "lineage_agree": int((j["lineage"] == j["lineage_native"]).sum()),
        })
        ct = pd.crosstab(j["lineage_native"], j["lineage"])
        ct.to_csv(OUT_ROOT / f"compare_draft_donor_{donor}.csv")
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_ROOT / "compare_draft_summary.csv", index=False)
    out["per_donor"] = rows
    (OUT_ROOT / "compare_draft.json").write_text(json.dumps(out, indent=1))
    print(frame.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", choices=["call", "compare"], default="call")
    ap.add_argument("--only", default="", help="comma-separated section keys or donors")
    args = ap.parse_args()
    if args.stage == "call":
        run_calls(args)
    else:
        run_compare(args)


if __name__ == "__main__":
    main()
