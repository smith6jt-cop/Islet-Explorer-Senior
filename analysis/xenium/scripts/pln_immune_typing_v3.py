#!/usr/bin/env python3
"""pLN Xenium immune typing, v3: cluster, then name, then SPLIT T_NK BY THE GENE PANEL.

WHAT CHANGED FROM v2, AND WHY
-----------------------------
Maintainer, 2026-09-17: *"You say NK/ILC cannot be named when all you have to do is separate
the cells based on CD8 and then apply the gene panel to those cells."*

v2 named clusters and stopped.  A cluster whose winning identity core was `Immune_T` swallowed
every NK cell in it, and five donors came back with `Immune_NK_ILC = 0` -- reported as if the
population were absent rather than as if the unit of the call were too coarse to see it.

**Negative evidence is evidence.**  NK/ILC is CD3E/CD3G-NEGATIVE with the NK identity core
standing up; T is CD3E/CD3G-POSITIVE and then splits on CD4 / CD8A+CD8B / FOXP3.  No single
probe is "the NK probe" and none is needed: the split is a presence/absence PATTERN inside one
coarse class, exactly as the PhenoCycler panel carries `T_NK: (CD3e, CD8, CD4, FOXP3, CD56)` as
one class and splits it by pattern.  So v3 adds STAGE 3: within the T_NK compartment -- the
cells whose cluster was called `Immune_T` or `Immune_NK_ILC` -- every cell is resolved by the
genes it detects, using `core/panels/data/identity_core/05_immune_NK_ILC.csv` and
`03_immune_T_cell.csv` and THEIR OWN `evidence` column to decide which genes may name a lineage.

Nothing in stage 3 is tuned to a frequency, a cross-donor uniformity or a disease ordering, and
nothing is tuned to make the two modalities agree.  The subsets sum to the compartment, and the
residue -- cells with neither a CD3 detection nor enough NK evidence -- is NAMED and counted,
never redistributed.

ALSO FIXED HERE (measured defects, not opinions)
------------------------------------------------
D2  ONE detection predicate, everywhere in this lane: `lognorm >= log1p(1)` with lognorm built
    at a FIXED target of 1e4.  v2 declared it and used it; the clustering track it consumed
    reported `det__` at `lognorm > 0` under a section-median normalisation, so two different
    predicates wore the word "detection" in one release.  v3 derives every detection fraction
    it reports from its own predicate, asserts the constant at run time, and MEASURES the
    divergence: `entries_nonzero_but_below_predicate` per section is exactly the set on which
    `> 0` and `>= log1p(1)` disagree.  The Y-probe false-detection rate in the female donors is
    the control on that predicate and is reported per section and per gene.

D3  THE PARTITION IS RECOMPUTED WITHOUT THE NO-NUCLEUS CELLS.  v2 consumed the clustering
    track's partition and then dropped `nucleus_count == 0` cells from the per-cluster means --
    but those cells were inside the Leiden subsample and inside the k-NN vote, so they had
    already shaped the clusters the calls rest on.  "A cell with no nucleus is excluded from
    every calculation" is not satisfied by excluding it from the last calculation.  Stage
    `partition` re-runs the same algorithm (`cluster_pln.leiden_on_subsample`, same seed
    derivation, same QC) on the nucleated cells only, into this release's own tree.  The
    shared `data/development/typing/pln/clusters/xenium/` is NOT overwritten.

D4  A partial run cannot land on a cohort-wide name.  Every cohort-shaped write goes through
    `partial_run_guard.RunScope`, the same rule the pLN lanes already carry.

EXTENSION
---------
Sections are DISCOVERED from the h5ad directory and the donor eligibility comes from
`phenocycler.cohort.filter_eligible_donors`; no donor list is written down here.  A section
already partitioned is skipped unless `--force`, and the cohort-wide rollups are rebuilt from
whatever per-section artifacts are on disk at the time.  Adding the donors now transferring is
therefore: drop the h5ad in, add its bundle row to `xenium_paths.csv`, re-run.  Nothing in this
file has to change.

NOT LICENSED
------------
No spatial claim; no cell-to-cell matching across modalities; no acinar or ductal, ever.  The
cross-modality comparison is COMPOSITIONAL and joined on `donor_id` alone -- different sections
of different blocks.

Run:
    .venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage partition
    .venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage call
    .venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage cross
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
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
OUT_ROOT = REPO / "data/development/xenium_typing/release_pln_immune_v3"
V2 = REPO / "data/development/xenium_typing/release_pln_immune_v2"
SHARED_CLUSTERS = REPO / "data/development/typing/pln/clusters/xenium"
IDENTITY_CORE = REPO / "core/panels/data/identity_core"
GDNA = REPO / "data/development/xenium_gdna/gdna_triage_2026-09-17.parquet"
XENIUM_ROOT = Path("/data/xenium")
XENIUM_PATHS = REPO / "core/phenocycler/data/integration/xenium_paths.csv"
PHENO_CALLS = REPO / "data/development/typing/pln/pheno"

sys.path.insert(0, str(REPO / "core/phenocycler/scripts/pln"))
from partial_run_guard import RunScope, write_scope_note  # noqa: E402

SCHEMA_VERSION = 1
RELEASE_ID = "xenium-pln-immune-tnk-split-dev-20260917"
SPACE_ID = "xenium-pln-immune-tnk-split-v3"
LANE = "xenium-pln-immune-tnk-split"
ROI = "pln_1"

# ---- the one detection predicate in this lane -------------------------------------------
TARGET_SUM = 1.0e4
DETECTION_MIN_LOGNORM = float(np.log1p(1.0))
DETECTION_PREDICATE = "lognorm >= log1p(1), lognorm = log1p(raw * 1e4 / cell_depth)"

# ---- clustering, carried verbatim from the clustering track ------------------------------
CLUSTER_SUBSAMPLE = 60_000
CLUSTER_RESOLUTION = 1.0
CLUSTER_MIN_COUNTS = 10
CLUSTER_MIN_FEATURES = 5
CLUSTER_SEED_BASE = 20260917

Z_STANDOUT = 1.0
Z_SENSITIVITY = (0.75, 1.0, 1.5)
STAGE1_MIN_STANDOUT = 2
STAGE2_MIN_NAMING_STANDOUT = 1
#: How many NK-exclusive naming genes the NK/ILC call needs.  One constant, one meaning, at
#: both units it is applied: DETECTED on a cell for the per-cell gate, STANDING OUT in a
#: sub-cluster for the sub-cluster call.  Two rather than one because a single gene named a
#: sub-cluster NK on the strength of RORC alone in the first draft, and a sub-cluster whose
#: only NK evidence is one ILC-hint gene is not an NK sub-cluster.
NK_MIN_NAMING = 2
NK_SENSITIVITY = (1, 2, 3)

STAGE1_IMMUNE = [
    "PTPRC", "CD3E", "CD3G", "CD4", "CD8A", "MS4A1", "CD79A", "CD19",
    "CD68", "ITGAX", "ITGAM", "NKG7", "KLRD1",
]
NAMED_IN_BRIEF_CHECK = ["CD3D", "LYZ"]

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
NAMING_EVIDENCE = frozenset({"canonical", "both"})

# ---- stage 3, the T/NK split -------------------------------------------------------------
#: The gate, in the maintainer's own terms.  CD3D is absent from this panel and is recorded as
#: such; CD247 (CD3-zeta) is NOT used, because NK cells carry CD3-zeta and it would collapse
#: the very distinction being drawn.  TRAC is genuinely T-restricted and highly expressed but
#: is not in the panel authority's T core, so it is reported as a SENSITIVITY, never as the gate.
CD3_GATE = ["CD3E", "CD3G"]
CD3_NOT_USED = {
    "CD3D": "absent from the 5K panel; measured, never substituted",
    "CD247": "CD3-zeta is expressed by NK cells as well as T cells, so it cannot separate them",
}
T_SENSITIVITY_GENE = "TRAC"
T_SUBSET_GENES = {"CD4": ["CD4"], "CD8": ["CD8A", "CD8B"], "Treg": ["FOXP3"]}
TREG_SUPPORT = ["IKZF2", "IL2RA", "CTLA4"]
CYTOTOXIC_SUPPORT = ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "KLRD1", "FCGR3A"]
#: ILC hints come from the NK/ILC core's own `subtype_hint` column; they are DESCRIPTIVE here.
ILC_HINT_VALUES = frozenset({"ILC1", "ILC2", "ILC3"})

TNK_LABELS = [
    "T_CD4", "T_CD8", "T_Treg", "T_CD4_CD8_both", "T_subset_unresolved",
    "NK_ILC", "T_NK_unresolved",
]

Y_FLOOR = ["TBL1Y", "USP9Y", "EIF1AY"]
Y_CROSSHYB = ["DDX3Y", "KDM5D"]
Y_GENES = Y_FLOOR + Y_CROSSHYB
REGIONAL = ["CR2", "CXCL13", "CCL19", "CCL21", "BCL6", "AICDA", "SELL", "CD27", "MKI67"]
CONTROL_PREFIXES = ("NegControl", "Intergenic", "Unassigned", "Deprecated", "BLANK", "antisense")

SCOPE = RunScope.build("section", None, ())


def set_scope(scope: RunScope) -> RunScope:
    global SCOPE
    SCOPE = scope
    return SCOPE


def cohort_out(out_dir: Path, name: str) -> Path:
    return SCOPE.cohort_path(out_dir, name)


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


def cluster_module():
    """`cluster_pln` by path, for `leiden_on_subsample` and `stable_seed`.

    Imported rather than copied so the partition this release rests on is produced by the same
    code the clustering track uses -- the only difference being which cells go in.
    """
    spec = importlib.util.spec_from_file_location(
        "cluster_pln_for_v3", REPO / "core/phenocycler/scripts/pln/cluster_pln.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------
# panels
# ------------------------------------------------------------------------------------------

def read_core(name: str) -> dict:
    frame = pd.read_csv(IDENTITY_CORE / name)
    frame["evidence"] = frame["evidence"].astype(str).str.strip()
    genes = [str(g) for g in frame["gene"].dropna().unique()]
    naming = sorted({str(r.gene) for r in frame.itertuples()
                     if r.evidence in NAMING_EVIDENCE})
    hints = {str(r.gene): str(getattr(r, "subtype_hint", "") or "").strip()
             for r in frame.itertuples()}
    return {"genes": genes, "naming": naming, "hint": hints}


def load_panels() -> dict:
    immune = {k: read_core(v) for k, v in IMMUNE_CORES.items()}
    non_immune = {k: read_core(v) for k, v in NON_IMMUNE_CORES.items()}
    structural_naming = sorted({g for c in non_immune.values() for g in c["naming"]})

    t_core = immune["Immune_T"]
    nk_core = immune["Immune_NK_ILC"]
    #: A gene that appears ANYWHERE in the T core -- naming or cellxgene -- cannot be NK
    #: evidence in a T/NK split, because its presence would be explained by either side.
    nk_exclusive = sorted(set(nk_core["naming"]) - set(t_core["genes"]))
    nk_shared_with_t = sorted(set(nk_core["naming"]) & set(t_core["genes"]))
    ilc_hint = sorted(g for g in nk_exclusive if nk_core["hint"].get(g, "") in ILC_HINT_VALUES)
    nk_like = sorted(set(nk_exclusive) - set(ilc_hint))

    wanted = sorted(
        set(STAGE1_IMMUNE) | set(Y_GENES) | set(REGIONAL) | set(structural_naming)
        | set(CD3_GATE) | {T_SENSITIVITY_GENE} | set(TREG_SUPPORT) | set(CYTOTOXIC_SUPPORT)
        | {g for gs in T_SUBSET_GENES.values() for g in gs}
        | {g for c in immune.values() for g in c["genes"]}
        | {g for c in non_immune.values() for g in c["genes"]}
    )
    return {
        "immune": immune,
        "non_immune": non_immune,
        "structural_naming": structural_naming,
        "nk_exclusive_naming": nk_exclusive,
        "nk_naming_shared_with_t_core": nk_shared_with_t,
        "nk_ilc_hint_genes": ilc_hint,
        "nk_like_genes": nk_like,
        "t_naming": t_core["naming"],
        "wanted": wanted,
    }


# ------------------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------------------

def sections() -> list[tuple[str, str, Path]]:
    """(donor_id, section_key, h5ad) for every pLN h5ad an eligible donor has.

    Discovered from the directory, never from a list in this file, so a donor that lands
    tomorrow is picked up by re-running.  The central cohort exclusion is asked for, not
    restated.
    """
    sys.path.insert(0, str(REPO / "core/phenocycler"))
    from phenocycler.cohort import filter_eligible_donors

    found = []
    for path in sorted(H5AD_DIR.glob("*_pLNXenium.h5ad")):
        stem = path.name.replace("_pLNXenium.h5ad", "")
        core = stem[3:] if stem.startswith("HDL") else stem
        donor = core.rstrip("abc") if core[-1:].isalpha() else core
        found.append((donor, stem, path))
    eligible = set(filter_eligible_donors(sorted({d for d, _, _ in found})))
    return [row for row in found if row[0] in eligible]


def refused_donors() -> list[str]:
    sys.path.insert(0, str(REPO / "core/phenocycler"))
    from phenocycler.cohort import filter_eligible_donors

    found = set()
    for path in sorted(H5AD_DIR.glob("*_pLNXenium.h5ad")):
        stem = path.name.replace("_pLNXenium.h5ad", "")
        core = stem[3:] if stem.startswith("HDL") else stem
        found.add(core.rstrip("abc") if core[-1:].isalpha() else core)
    return sorted(found - set(filter_eligible_donors(sorted(found))))


def nucleus_for(donor: str, cell_ids: np.ndarray) -> tuple[np.ndarray | None, str]:
    """`nucleus_count` per h5ad cell from the vendor bundle, or (None, why-not).

    Matched by donor AND by subset containment of this section's own cell ids -- never by
    filename order, never by row order.  Where no local bundle contains them the rule is
    UNMEASURABLE and is recorded as such; it is never recorded as passed.
    """
    table = pd.read_csv(XENIUM_PATHS, dtype=str)
    table = table[
        table["roi"].astype(str).str.startswith("pln")
        & (table["donor_id"].astype(str) == str(donor))
    ]
    want = set(cell_ids.tolist())
    tried = []
    for row in table.itertuples():
        src = str(row.source_path or "").strip()
        if not src or src.lower() == "nan" or src == "NA":
            continue
        cells = XENIUM_ROOT / src / "cells.parquet"
        if not cells.is_file():
            tried.append(f"{row.roi}:bundle_absent")
            continue
        frame = pd.read_parquet(cells, columns=["cell_id", "nucleus_count"])
        frame["cell_id"] = frame["cell_id"].astype(str)
        lookup = dict(zip(frame["cell_id"], frame["nucleus_count"].to_numpy()))
        missing = len(want - lookup.keys())
        if missing:
            tried.append(f"{row.roi}:not_a_superset({missing}_ids_missing)")
            continue
        counts = np.array([lookup[c] for c in cell_ids], dtype=np.int32)
        return counts, f"bundle {row.roi} {cells}"
    return None, "UNMEASURED: no local bundle cells.parquet is a superset of this section" + (
        f" (tried {'; '.join(tried)})" if tried else ""
    )


# ------------------------------------------------------------------------------------------
# stage: partition
# ------------------------------------------------------------------------------------------

def partition_section(donor: str, section: str, path: Path, mod) -> dict:
    """The clustering track's algorithm, on the NUCLEATED cells of this section only."""
    import scanpy as sc
    import anndata as ad

    t0 = time.time()
    adata = ad.read_h5ad(path)
    n_in = int(adata.n_obs)
    cell_ids = adata.obs_names.to_numpy().astype(str)

    nucleus_count, nucleus_source = nucleus_for(donor, cell_ids)
    if nucleus_count is None:
        has_nucleus = np.ones(n_in, dtype=bool)
        nucleus_measured = False
    else:
        has_nucleus = nucleus_count > 0
        nucleus_measured = True

    counts = adata.obs["nCount_Xenium"].to_numpy(dtype=float)
    genes_det = adata.obs["nFeature_Xenium"].to_numpy(dtype=float)
    passes_qc = (counts >= CLUSTER_MIN_COUNTS) & (genes_det >= CLUSTER_MIN_FEATURES)

    n_qc_only = int(passes_qc.sum())                       # D3 "before": the v2 partition rule
    keep = passes_qc & has_nucleus
    n_keep = int(keep.sum())                               # D3 "after"
    n_no_nucleus_in_qc = int((passes_qc & ~has_nucleus).sum())

    adata = adata[keep].copy()
    kept_ids = cell_ids[keep]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    work = adata.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=2000)
    work = work[:, work.var["highly_variable"]].copy()
    sc.pp.scale(work, max_value=10)
    feats = np.asarray(work.X, dtype=np.float32)
    del work, adata
    gc.collect()

    seed = mod.stable_seed("pln-cluster", CLUSTER_SEED_BASE, "xenium", section)
    labels, meta = mod.leiden_on_subsample(
        feats, n_pcs=30, subsample=CLUSTER_SUBSAMPLE,
        resolution=CLUSTER_RESOLUTION, seed=seed,
    )
    del feats
    gc.collect()
    cluster_id = np.array([f"{section}_c{int(v):02d}" for v in labels])

    dest = OUT_ROOT / "partition/assignments" / f"donor_id={donor}" / f"section={section}"
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "cell_id": kept_ids,
        "donor_id": donor,
        "section_key": section,
        "modality": "xenium",
        "roi": ROI,
        "cluster_id": cluster_id,
    }).to_parquet(dest / "clusters.parquet", index=False)

    report = {
        "donor_id": donor,
        "section_key": section,
        "h5ad": str(path),
        "n_cells_in_h5ad": n_in,
        "qc_rule": f"nCount_Xenium >= {CLUSTER_MIN_COUNTS} and "
                   f"nFeature_Xenium >= {CLUSTER_MIN_FEATURES}",
        "n_cells_passing_qc_only": n_qc_only,
        "n_cells_no_nucleus_inside_qc": n_no_nucleus_in_qc,
        "n_cells_partitioned": n_keep,
        "nucleus_rule_measurable": nucleus_measured,
        "nucleus_source": nucleus_source,
        "n_clusters": int(len(set(cluster_id))),
        "seed": int(seed),
        **{k: v for k, v in meta.items() if k != "n_cells"},
        "seconds": round(time.time() - t0, 1),
    }
    gc.collect()
    return report


def run_partition(args) -> None:
    mod = cluster_module()
    secs = sections()
    scope = set_scope(SCOPE.refine([s for _, s, _ in secs]))
    if scope.missing:
        raise SystemExit(f"no pLN h5ad for: {', '.join(scope.missing)}")
    scope.announce(OUT_ROOT, log)
    todo = [row for row in secs if row[1] in set(scope.requested)]

    rdir = OUT_ROOT / "partition/section_reports"
    rdir.mkdir(parents=True, exist_ok=True)
    for donor, section, path in todo:
        dest = (OUT_ROOT / "partition/assignments" /
                f"donor_id={donor}/section={section}/clusters.parquet")
        if dest.exists() and not args.force:
            log(f"section {section}: already partitioned, skipped (--force to redo)")
            continue
        log(f"partition {section} (donor {donor})")
        rep = partition_section(donor, section, path, mod)
        (rdir / f"{section}.json").write_text(json.dumps(rep, indent=1))
        log(f"  {rep['n_cells_partitioned']:,} partitioned of "
            f"{rep['n_cells_passing_qc_only']:,} passing QC; "
            f"{rep['n_cells_no_nucleus_inside_qc']:,} no-nucleus removed; "
            f"{rep['n_clusters']} clusters; {rep['seconds']}s")

    reports = [json.loads(p.read_text()) for p in sorted(rdir.glob("*.json"))]
    frame = pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, (dict, list))}
                          for r in reports])
    frame.to_csv(cohort_out(OUT_ROOT / "partition", "partition_summary.csv"), index=False)
    write_scope_note(scope, OUT_ROOT / "partition")
    print(frame[[
        "section_key", "n_cells_in_h5ad", "n_cells_passing_qc_only",
        "n_cells_no_nucleus_inside_qc", "n_cells_partitioned", "nucleus_rule_measurable",
        "n_clusters", "knn_holdout_agreement",
    ]].to_string(index=False))


# ------------------------------------------------------------------------------------------
# stage: call
# ------------------------------------------------------------------------------------------

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


def call_section(donor: str, section: str, h5ad: Path, panels: dict) -> dict:
    import anndata as ad

    assert abs(DETECTION_MIN_LOGNORM - float(np.log1p(1.0))) < 1e-12, (
        "the detection predicate must be lognorm >= log1p(1); it is what the release documents"
    )

    t0 = time.time()
    adata = ad.read_h5ad(h5ad)
    cell_ids = adata.obs_names.to_numpy().astype(str)
    var_names = adata.var_names.to_numpy().astype(str)
    n_cells = int(adata.n_obs)
    controls_present = sorted(g for g in var_names
                              if any(g.startswith(p) for p in CONTROL_PREFIXES))
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

    raw_nonzero = int(sub.data.size)
    # float64 accumulator: 6593 carries 27.4M stored entries and a float32 sum over them is not
    # exact above 2**24, which made this count come back one HIGHER than the entries it counts
    det_total = int(detected.sum(dtype=np.float64))
    deep = int(np.sum(depth > TARGET_SUM))

    L = sp.csr_matrix((lognorm, sub.indices, sub.indptr), shape=sub.shape)
    D = sp.csr_matrix((detected, sub.indices, sub.indptr), shape=sub.shape)
    del sub, scaled, lognorm, detected, rows
    gc.collect()

    assign = pd.read_parquet(
        OUT_ROOT / f"partition/assignments/donor_id={donor}/section={section}/clusters.parquet"
    )
    amap = dict(zip(assign["cell_id"].astype(str), assign["cluster_id"].astype(str)))
    cluster_of = np.array([amap.get(c, "") for c in cell_ids])
    clustered = cluster_of != ""

    nucleus_count, nucleus_source = nucleus_for(donor, cell_ids)
    nucleus_measured = nucleus_count is not None
    has_nucleus = np.ones(n_cells, dtype=bool) if not nucleus_measured else (nucleus_count > 0)
    n_no_nucleus = 0 if not nucleus_measured else int((~has_nucleus).sum())
    # the partition already excluded them; this asserts it rather than assuming it
    n_no_nucleus_in_partition = int((clustered & ~has_nucleus).sum())
    if n_no_nucleus_in_partition:
        raise SystemExit(
            f"{section}: {n_no_nucleus_in_partition} nucleus_count==0 cells are inside the "
            "partition -- re-run --stage partition --force for this section"
        )

    order = sorted(set(cluster_of[clustered]))
    lookup = {c: i for i, c in enumerate(order)}
    krow = np.array([lookup[c] for c in cluster_of[clustered]])
    cell_pos = np.flatnonzero(clustered)
    indic = sp.csr_matrix(
        (np.ones(krow.size, dtype=np.float32), (krow, cell_pos)),
        shape=(len(order), n_cells),
    )
    sizes = np.asarray(indic.sum(axis=1)).ravel()
    n_clustered = int(sizes.sum())

    means = np.asarray((indic @ L).todense()) / sizes[:, None]
    dets = np.asarray((indic @ D).todense()) / sizes[:, None]
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
    def idx_of(genes):
        return np.array([pidx[g] for g in genes if g in pidx], dtype=int)

    imm_cols = idx_of(STAGE1_IMMUNE)
    str_cols = idx_of(panels["structural_naming"])
    core_names = list(IMMUNE_CORES)
    core_cols = {c: idx_of(panels["immune"][c]["genes"]) for c in core_names}
    core_naming_cols = {c: idx_of(panels["immune"][c]["naming"]) for c in core_names}

    per_thr = {}
    for thr in Z_SENSITIVITY:
        so = standout[thr]
        n_imm = so[:, imm_cols].sum(axis=1)
        n_str = so[:, str_cols].sum(axis=1)
        f_imm = n_imm / max(len(imm_cols), 1)
        f_str = n_str / max(len(str_cols), 1)
        s1, s1r = [], []
        for k in range(len(order)):
            call, reason = stage1_decide(int(n_imm[k]), float(f_imm[k]),
                                         int(n_str[k]), float(f_str[k]))
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
        s2, s2r = [], []
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
    cl_stage1 = np.array(main["stage1"], dtype=object)
    cl_stage2 = np.array(main["stage2"], dtype=object)
    cl_reason = np.array(main["stage2_reason"], dtype=object)

    # ---- per-cell descriptive core argmax ------------------------------------------------
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
    flat = cell_core_mean.max(axis=1) <= 0
    cell_argmax_label = np.array(core_names, dtype=object)[cell_argmax]
    cell_argmax_label[flat] = ""

    # =======================================================================================
    # STAGE 3 -- split the T_NK compartment by the gene panel, per cell
    # =======================================================================================
    def det_any(genes) -> np.ndarray:
        ii = idx_of(genes)
        if not ii.size:
            return np.zeros(n_cells, dtype=bool)
        return np.asarray(D[:, ii].sum(axis=1)).ravel() > 0

    def det_count(genes) -> np.ndarray:
        ii = idx_of(genes)
        if not ii.size:
            return np.zeros(n_cells, dtype=np.int16)
        return np.asarray(D[:, ii].sum(axis=1)).ravel().astype(np.int16)

    cd3_pos = det_any(CD3_GATE)
    trac_pos = det_any([T_SENSITIVITY_GENE])
    cd4_pos = det_any(T_SUBSET_GENES["CD4"])
    cd8_pos = det_any(T_SUBSET_GENES["CD8"])
    foxp3_pos = det_any(T_SUBSET_GENES["Treg"])
    treg_support = det_count(TREG_SUPPORT)
    nk_n = det_count(panels["nk_exclusive_naming"])
    nk_like_n = det_count(panels["nk_like_genes"])
    ilc_hint_n = det_count(panels["nk_ilc_hint_genes"])
    t_naming_n = det_count(panels["t_naming"])
    cyto_n = det_count(CYTOTOXIC_SUPPORT)

    cl_common = np.array([COMMON_LINEAGE.get(x, "") for x in cl_stage2], dtype=object)
    cidx = np.full(n_cells, -1, dtype=int)
    cidx[cell_pos] = krow
    in_tnk = np.zeros(n_cells, dtype=bool)
    in_tnk[cell_pos] = cl_common[krow] == "T_NK"

    def tnk_labels(nk_min: int) -> np.ndarray:
        lab = np.full(n_cells, "", dtype=object)
        t = in_tnk & cd3_pos
        lab[t & foxp3_pos] = "T_Treg"
        rest = t & ~foxp3_pos
        lab[rest & cd8_pos & ~cd4_pos] = "T_CD8"
        lab[rest & cd4_pos & ~cd8_pos] = "T_CD4"
        lab[rest & cd4_pos & cd8_pos] = "T_CD4_CD8_both"
        lab[rest & ~cd4_pos & ~cd8_pos] = "T_subset_unresolved"
        neg = in_tnk & ~cd3_pos
        lab[neg & (nk_n >= nk_min)] = "NK_ILC"
        lab[neg & (nk_n < nk_min)] = "T_NK_unresolved"
        return lab

    tnk = tnk_labels(NK_MIN_NAMING)
    nk_sens = {
        str(m): {k: int(v) for k, v in pd.Series(tnk_labels(m)[in_tnk]).value_counts().items()}
        for m in NK_SENSITIVITY
    }
    # sensitivity on the T side: CD3E/CD3G is the gate; TRAC is reported, never the gate
    trac_rescue = int((in_tnk & ~cd3_pos & trac_pos & (nk_n < NK_MIN_NAMING)).sum())

    # ---- the dropout control the negative gate needs --------------------------------------
    # A compartment that is definitionally CD3-negative gives the false-positive floor of the
    # CD3 gate; the T_NK compartment gives its sensitivity ceiling.  Both are measured on this
    # section's own cells, and neither is a target.
    def compartment(name: str) -> np.ndarray:
        m = np.zeros(n_cells, dtype=bool)
        m[cell_pos] = cl_stage2[krow] == name
        return m

    b_comp = compartment("Immune_B_Plasma")
    my_comp = compartment("Myeloid")
    def nk_rule_rate(mask):
        """The NK_ILC rule, applied verbatim to a compartment that is not NK.

        Not `>= NK_MIN naming genes` on its own -- the RULE, CD3-negative AND the evidence
        count -- because that is what would have to fire falsely for an NK call to be wrong.
        """
        m = mask & ~cd3_pos
        if not m.any():
            return None
        return float((nk_n[m] >= NK_MIN_NAMING).mean())

    # where in the compartment the NK calls come from: a T-called cluster or an NK-called one
    tnk_from_t = np.zeros(n_cells, dtype=bool)
    tnk_from_t[cell_pos] = cl_stage2[krow] == "Immune_T"
    tnk_from_nk = np.zeros(n_cells, dtype=bool)
    tnk_from_nk[cell_pos] = cl_stage2[krow] == "Immune_NK_ILC"
    gate_control = {
        "cd3_detect_frac_in_T_NK": float(cd3_pos[in_tnk].mean()) if in_tnk.any() else None,
        "cd3_detect_frac_in_B_Plasma": float(cd3_pos[b_comp].mean()) if b_comp.any() else None,
        "cd3_detect_frac_in_Myeloid": float(cd3_pos[my_comp].mean()) if my_comp.any() else None,
        "nk_rule_false_positive_frac_in_B_Plasma": nk_rule_rate(b_comp),
        "nk_rule_false_positive_frac_in_Myeloid": nk_rule_rate(my_comp),
        "nk_rule_frac_in_T_NK": nk_rule_rate(in_tnk),
        "median_depth_T_NK": float(np.median(depth[in_tnk])) if in_tnk.any() else None,
        "median_depth_B_Plasma": float(np.median(depth[b_comp])) if b_comp.any() else None,
        "nk_calls_from_Immune_T_clusters": int(((tnk == "NK_ILC") & tnk_from_t).sum()),
        "nk_calls_from_Immune_NK_ILC_clusters": int(((tnk == "NK_ILC") & tnk_from_nk).sum()),
        "cells_in_Immune_T_clusters": int(tnk_from_t.sum()),
        "cells_in_Immune_NK_ILC_clusters": int(tnk_from_nk.sum()),
    }

    # ---- Y-probe control ------------------------------------------------------------------
    y_cols = idx_of(Y_GENES)
    y_any = np.asarray(D[:, y_cols].sum(axis=1)).ravel() if y_cols.size else np.zeros(n_cells)
    floor_cols = idx_of(Y_FLOOR)
    y_floor = (np.asarray(D[cell_pos][:, floor_cols].sum(axis=1)).ravel()
               if floor_cols.size else np.zeros(len(cell_pos)))
    y_per_gene = {g: float(np.asarray(D[cell_pos][:, [pidx[g]]].todense()).ravel().mean())
                  for g in Y_GENES if g in pidx}

    # ---- per-cell frame --------------------------------------------------------------------
    fallback = np.where(has_nucleus, "Uncalled_not_partitioned", "Unavailable_no_nucleus")
    fb_reason = np.where(
        has_nucleus,
        f"below partition QC (nCount_Xenium < {CLUSTER_MIN_COUNTS} or "
        f"nFeature_Xenium < {CLUSTER_MIN_FEATURES})",
        "nucleus_count == 0: excluded from the partition and from every calculation",
    )
    stage1_cell = np.where(cidx >= 0, cl_stage1[np.clip(cidx, 0, None)], fallback)
    stage2_cell = np.where(cidx >= 0, cl_stage2[np.clip(cidx, 0, None)], fallback)
    reason_cell = np.where(cidx >= 0, cl_reason[np.clip(cidx, 0, None)], fb_reason)
    lineage = np.where(in_tnk, tnk.astype(str), stage2_cell)

    call_core_j = np.array([core_names.index(x) if x in core_names else -1 for x in cl_stage2])
    per_cell_core_j = np.where(cidx >= 0, call_core_j[np.clip(cidx, 0, None)], -1)
    agrees = np.where(per_cell_core_j >= 0,
                      cell_argmax == np.clip(per_cell_core_j, 0, None), False) \
        & (cell_argmax_label != "")

    calls = pd.DataFrame({
        "cell_id": cell_ids,
        "donor_id": donor,
        "section_key": section,
        "roi": ROI,
        "cluster_id": np.where(cidx >= 0, cluster_of, None),
        "stage1": stage1_cell.astype(str),
        "cluster_lineage": stage2_cell.astype(str),
        "lineage": lineage.astype(str),
        "common_lineage": np.where(
            in_tnk, "T_NK",
            [COMMON_LINEAGE.get(x, "") for x in stage2_cell]).astype(str),
        "reason": reason_cell.astype(str),
        "in_T_NK_compartment": in_tnk,
        "cd3_detected": cd3_pos,
        "trac_detected": trac_pos,
        "cd4_detected": cd4_pos,
        "cd8_detected": cd8_pos,
        "foxp3_detected": foxp3_pos,
        "n_treg_support_detected": treg_support,
        "n_nk_exclusive_naming_detected": nk_n,
        "n_nk_like_detected": nk_like_n,
        "n_ilc_hint_detected": ilc_hint_n,
        "n_t_naming_detected": t_naming_n,
        "n_cytotoxic_support_detected": cyto_n,
        "cell_core_argmax": cell_argmax_label.astype(str),
        "cell_agrees_with_cluster_call": agrees,
        "y_probe_detected": (y_any > 0),
        "nucleus_count": (nucleus_count if nucleus_measured
                          else np.full(n_cells, -1, dtype=np.int32)),
        "nucleus_measured": nucleus_measured,
        "depth": depth.astype(np.int32),
    })

    # ---- per-cluster frame -------------------------------------------------------------------
    cd3_so = np.zeros(len(order), dtype=bool)
    ii = idx_of(CD3_GATE)
    if ii.size:
        cd3_so = standout[Z_STANDOUT][:, ii].any(axis=1)
    het = []
    for k, cid in enumerate(order):
        sel = cell_pos[krow == k]
        dist = pd.Series(cell_argmax_label[sel]).value_counts(normalize=True)
        called = cl_stage2[k]
        het.append({
            "donor_id": donor, "section_key": section, "cluster_id": cid,
            "n_cells": int(sizes[k]),
            "frac_of_section": float(sizes[k] / n_clustered),
            "median_depth": float(np.median(depth[sel])),
            "median_n_features": float(np.median(nfeature[sel])),
            "stage1": cl_stage1[k], "stage1_reason": main["stage1_reason"][k],
            "cluster_lineage": called, "stage2_reason": cl_reason[k],
            "common_lineage": COMMON_LINEAGE.get(called, ""),
            "cd3_stands_out": bool(cd3_so[k]),
            "cd3_det_frac": float(dets[k, pidx["CD3E"]]) if "CD3E" in pidx else np.nan,
            "cd3_det_frac_rest_of_section": (float(d_other[k, pidx["CD3E"]])
                                             if "CD3E" in pidx else np.nan),
            "n_immune_canonical_standout": int(main["n_imm"][k]),
            "n_structural_naming_standout": int(main["n_str"][k]),
            "frac_immune_canonical_standout": float(main["f_imm"][k]),
            "frac_structural_naming_standout": float(main["f_str"][k]),
            "winning_core": core_names[int(main["win"][k])],
            "winning_core_score": float(main["scores"][k, int(main["win"][k])]),
            "winning_core_naming_standouts": int(main["naming_so"][k, int(main["win"][k])]),
            "score_margin_over_runner_up": float(main["margin"][k]),
            **{f"score__{c}": float(main["scores"][k, j]) for j, c in enumerate(core_names)},
            **{f"naming_standouts__{c}": int(main["naming_so"][k, j])
               for j, c in enumerate(core_names)},
            **({f"tnk__{lab}": int((tnk[sel] == lab).sum()) for lab in TNK_LABELS}
               if COMMON_LINEAGE.get(called, "") == "T_NK"
               else {f"tnk__{lab}": 0 for lab in TNK_LABELS}),
            "purity_cell_argmax_agrees": (float(np.mean(cell_argmax_label[sel] == called))
                                          if called in core_names else np.nan),
            "cell_argmax_top_label": dist.index[0] if len(dist) else "",
            "cell_argmax_top_frac": float(dist.iloc[0]) if len(dist) else np.nan,
            "standout_immune_genes": ";".join(
                g for g in STAGE1_IMMUNE if g in pidx and standout[Z_STANDOUT][k, pidx[g]]),
            "standout_nk_exclusive_genes": ";".join(
                g for g in panels["nk_exclusive_naming"]
                if g in pidx and standout[Z_STANDOUT][k, pidx[g]]),
            "standout_structural_genes": ";".join(
                g for g in panels["structural_naming"]
                if g in pidx and standout[Z_STANDOUT][k, pidx[g]]),
        })

    gene_stats = pd.DataFrame({
        "donor_id": donor, "section_key": section,
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
            "cells": {k: int(sizes[[i for i, x in enumerate(per_thr[thr]["stage2"]) if x == k]].sum())
                      for k in s2.unique()},
        }

    report = {
        "donor_id": donor, "section_key": section, "h5ad": str(h5ad),
        "n_cells_in_h5ad": n_cells,
        "n_cells_partitioned": n_clustered,
        "n_cells_not_partitioned": n_cells - n_clustered,
        "n_cells_no_nucleus": n_no_nucleus,
        "n_no_nucleus_cells_inside_partition": n_no_nucleus_in_partition,
        "nucleus_rule_measurable": nucleus_measured,
        "nucleus_source": nucleus_source,
        "n_clusters": len(order),
        "n_genes_on_panel": int(len(var_names)),
        "n_wanted_genes_present": len(present),
        "wanted_genes_absent": [g for g in panels["wanted"] if g not in gene_idx],
        "named_in_brief_but_absent": named_absent,
        "control_features_present_in_var": controls_present,
        "depth_equals_obs_nCount_for_cells": depth_matches_obs,
        "detection_predicate": DETECTION_PREDICATE,
        "cells_with_depth_gt_target": deep,
        "entries_raw_nonzero_in_wanted_set": raw_nonzero,
        "entries_detected_in_wanted_set": det_total,
        "entries_nonzero_but_below_predicate": raw_nonzero - det_total,
        "sex": str(adata.obs["sex"].iloc[0]),
        "disease_status": str(adata.obs["disease.status"].iloc[0]),
        "y_probe_detect_frac_any": float(np.mean(y_any[cell_pos] > 0)),
        "y_floor_group": Y_FLOOR,
        "y_floor_detect_frac_any": float(np.mean(y_floor > 0)),
        "y_probe_detect_frac_per_gene": y_per_gene,
        "t_nk_compartment_cells": int(in_tnk.sum()),
        "t_nk_split": {k: int(v) for k, v in pd.Series(tnk[in_tnk]).value_counts().items()},
        "t_nk_split_sensitivity_nk_min": nk_sens,
        "cd3_negative_trac_positive_not_called_nk": trac_rescue,
        "cd3_gate_control": gate_control,
        "sensitivity_z": sens,
        "seconds": round(time.time() - t0, 1),
    }
    del L, D, adata
    gc.collect()
    return {"calls": calls, "clusters": pd.DataFrame(het),
            "gene_stats": gene_stats, "report": report}


def run_calls(args) -> None:
    panels = load_panels()
    secs = sections()
    # `available` is what the LANE could process -- every section an eligible donor has an h5ad
    # for -- NOT what happens to be partitioned today.  Scoping against the partitioned subset
    # would let a run that covers one section call itself whole and land on the cohort name,
    # which is the defect this guard exists to stop.
    scope = set_scope(SCOPE.refine([s for _, s, _ in secs]))
    if scope.missing:
        raise SystemExit(f"no pLN h5ad for: {', '.join(scope.missing)}")
    todo = [row for row in secs if row[1] in set(scope.requested)]
    not_partitioned = [s for _, s, _ in todo
                       if not (OUT_ROOT / f"partition/assignments/donor_id="
                               f"{dict((k, d) for d, k, _ in todo)[s]}/section={s}/"
                               "clusters.parquet").exists()]
    if not_partitioned:
        raise SystemExit(
            "not partitioned, run --stage partition first: " + ", ".join(not_partitioned))
    scope.announce(OUT_ROOT, log)
    log(f"{len(todo)} sections, {len(set(d for d, _, _ in todo))} donors")

    calls_root = (OUT_ROOT / "calls" if not scope.is_partial
                  else scope.partial_dir(OUT_ROOT) / "calls")
    by_donor: dict[str, list[pd.DataFrame]] = {}
    clusters, gene_stats, reports, inputs = [], [], [], {}
    for donor, section, h5ad in todo:
        log(f"call {section} (donor {donor})")
        res = call_section(donor, section, h5ad, panels)
        by_donor.setdefault(donor, []).append(res["calls"])
        clusters.append(res["clusters"])
        gene_stats.append(res["gene_stats"])
        reports.append(res["report"])
        inputs[str(h5ad)] = sha256(h5ad)
        r = res["report"]
        log(f"  {r['n_cells_partitioned']:,} cells, {r['n_clusters']} clusters, "
            f"T_NK {r['t_nk_compartment_cells']:,} -> {r['t_nk_split']}; {r['seconds']}s")
        gc.collect()

    dup_report = {}
    for donor, frames in by_donor.items():
        frame = pd.concat(frames, ignore_index=True)
        dup = int(frame.duplicated(subset=["cell_id"]).sum())
        if dup:
            dup_report[donor] = dup
        dest = calls_root / f"donor_id={donor}/roi={ROI}"
        dest.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(dest / "data.parquet", index=False)
        log(f"wrote {dest/'data.parquet'} ({len(frame):,} rows)")

    cdir = cohort_out(OUT_ROOT, "clusters")
    cdir.mkdir(parents=True, exist_ok=True)
    cl = pd.concat(clusters, ignore_index=True)
    cl.to_parquet(cdir / "cluster_calls.parquet", index=False)
    cl.to_csv(cdir / "cluster_calls.csv", index=False)
    pd.concat(gene_stats, ignore_index=True).to_parquet(
        cdir / "cluster_gene_stats.parquet", index=False)

    rows = []
    for donor in sorted(by_donor):
        frame = pd.read_parquet(calls_root / f"donor_id={donor}/roi={ROI}/data.parquet",
                                columns=["stage1", "cluster_lineage", "lineage",
                                         "in_T_NK_compartment"])
        s1 = frame["stage1"].value_counts()
        lin = frame["lineage"].value_counts()
        tnk_cells = int(frame["in_T_NK_compartment"].sum())
        immune = int(s1.get("Immune", 0))
        row = {
            "donor_id": donor,
            "cells": len(frame),
            "immune": immune,
            "non_immune": int(s1.get("Non_immune", 0)),
            "uncalled_ambiguous": int(s1.get("Uncalled_ambiguous", 0)),
            "not_partitioned": int(s1.get("Uncalled_not_partitioned", 0)),
            "no_nucleus_excluded": int(s1.get("Unavailable_no_nucleus", 0)),
            "pct_immune": round(100.0 * immune / max(len(frame), 1), 2),
            "T_NK_compartment": tnk_cells,
        }
        for lab in TNK_LABELS:
            row[lab] = int(lin.get(lab, 0))
        for lab in ["Immune_B_Plasma", "Myeloid", "Granulocyte", "Immune_uncalled",
                    "Non_immune", "Uncalled_ambiguous", "Uncalled_not_partitioned",
                    "Unavailable_no_nucleus"]:
            row[lab] = int(lin.get(lab, 0))
        subset_sum = sum(row[lab] for lab in TNK_LABELS)
        row["T_NK_subsets_sum"] = subset_sum
        row["T_NK_coherent"] = subset_sum == tnk_cells
        row["T_total"] = sum(row[lab] for lab in TNK_LABELS if lab.startswith("T_")
                             and lab != "T_NK_unresolved")
        row["frac_NK_ILC_of_T_NK"] = (row["NK_ILC"] / tnk_cells) if tnk_cells else np.nan
        row["frac_NK_ILC_of_immune"] = (row["NK_ILC"] / immune) if immune else np.nan
        row["frac_T_of_immune"] = (row["T_total"] / immune) if immune else np.nan
        rows.append(row)
    summary = pd.DataFrame(rows)
    bad = summary.loc[~summary["T_NK_coherent"], "donor_id"].tolist()
    if bad:
        raise SystemExit(f"T_NK subsets do not sum to the compartment for: {bad}")
    summary.to_csv(cohort_out(OUT_ROOT, "per_donor_summary.csv"), index=False)

    gdna = {}
    if GDNA.exists():
        g = pd.read_parquet(GDNA)
        g = g[g["tissue"] == "pln"]
        gdna = {str(r.donor): float(r.rate) for r in g.itertuples()}
    donors = sorted(by_donor)

    ctrl = pd.DataFrame([{
        "section_key": r["section_key"], "donor_id": r["donor_id"], "sex": r["sex"],
        "disease_status": r["disease_status"],
        "n_cells_partitioned": r["n_cells_partitioned"],
        "y_any_detect_frac": r["y_probe_detect_frac_any"],
        "y_floor_group_detect_frac": r["y_floor_detect_frac_any"],
        **{f"y__{g}": v for g, v in r["y_probe_detect_frac_per_gene"].items()},
        "cells_with_depth_gt_target": r["cells_with_depth_gt_target"],
        "entries_nonzero_but_below_predicate": r["entries_nonzero_but_below_predicate"],
        "entries_raw_nonzero_in_wanted_set": r["entries_raw_nonzero_in_wanted_set"],
        **{k: v for k, v in r["cd3_gate_control"].items()},
        "gdna_pln_rate": gdna.get(r["donor_id"], np.nan),
    } for r in reports])
    ctrl.to_csv(cohort_out(OUT_ROOT, "controls.csv"), index=False)

    label_space = {
        "schema_version": SCHEMA_VERSION, "space_id": SPACE_ID, "producer_lane": LANE,
        "coarse_class": {
            "T_NK": "the combined class the panels carry (PhenoCycler lineages.py: "
                    "T_NK = CD3e, CD8, CD4, FOXP3, CD56). A cluster called Immune_T or "
                    "Immune_NK_ILC puts its cells in this compartment; the compartment is then "
                    "SPLIT per cell by the presence/absence pattern below, and the split labels "
                    "sum exactly to the compartment.",
        },
        "labels": (
            [{"label": "T_CD4", "common_lineage": "T_NK",
              "reason": "CD3E or CD3G detected; CD4 detected; neither CD8A nor CD8B; FOXP3 not "
                        "detected"},
             {"label": "T_CD8", "common_lineage": "T_NK",
              "reason": "CD3E or CD3G detected; CD8A or CD8B detected; CD4 not detected; FOXP3 "
                        "not detected"},
             {"label": "T_Treg", "common_lineage": "T_NK",
              "reason": "CD3E or CD3G detected and FOXP3 detected. FOXP3 takes precedence over "
                        "the CD4/CD8 split; IKZF2/IL2RA/CTLA4 detections are carried beside it "
                        "as support, not as the rule"},
             {"label": "T_CD4_CD8_both", "common_lineage": "T_NK",
              "reason": "CD3E or CD3G detected with BOTH CD4 and a CD8 chain detected. Named "
                        "rather than assigned to one side"},
             {"label": "T_subset_unresolved", "common_lineage": "T_NK",
              "reason": "CD3E or CD3G detected but none of CD4, CD8A, CD8B, FOXP3. The cell is T; "
                        "its subset is unavailable on this cell, and it is never spread across "
                        "the subsets"},
             {"label": "NK_ILC", "common_lineage": "T_NK",
              "reason": f"CD3E and CD3G BOTH undetected, and >= {NK_MIN_NAMING} genes that MAY "
                        "NAME the NK/ILC lineage (core/panels 05_immune_NK_ILC evidence in "
                        "{canonical, both}, minus every gene that also appears in the T core) "
                        "detected on that cell. Negative evidence is evidence: this is the "
                        "CD3-negative half of the combined T_NK class"},
             {"label": "T_NK_unresolved", "common_lineage": "T_NK",
              "reason": "in the T_NK compartment, no CD3E/CD3G detection and fewer than "
                        f"{NK_MIN_NAMING} NK naming genes. This is the residue the sparsity of "
                        "the assay leaves; it is counted and named, never redistributed"}]
            + [{"label": k, "common_lineage": COMMON_LINEAGE[k],
                "reason": "cluster-level identity-core call (stage 2)"}
               for k in ["Immune_B_Plasma", "Myeloid", "Granulocyte"]]
            + [
                {"label": "Immune_uncalled", "common_lineage": "",
                 "reason": "stage 1 called the cluster Immune; the winning identity core had "
                           "zero standout naming genes, so the call is withheld and the cluster "
                           "is NOT reassigned to the runner-up"},
                {"label": "Non_immune", "common_lineage": "",
                 "reason": "stage 1: >= 2 structural naming genes stand out and a larger "
                           "fraction of the structural set than of the immune canonical set"},
                {"label": "Uncalled_ambiguous", "common_lineage": "",
                 "reason": "stage 1 reached neither side; the reason is in the `reason` column"},
                {"label": "Unavailable_no_nucleus", "common_lineage": "",
                 "reason": "nucleus_count == 0. Excluded from the PARTITION as well as from "
                           "every aggregate, so it did not shape the clusters the calls rest on"},
                {"label": "Uncalled_not_partitioned", "common_lineage": "",
                 "reason": f"below the partition QC (nCount_Xenium < {CLUSTER_MIN_COUNTS} or "
                           f"nFeature_Xenium < {CLUSTER_MIN_FEATURES}); carried as uncalled, "
                           "never as 0"},
            ]
        ),
        "not_a_label": {
            "acinar_ductal": "never claimed; exocrine is one label in this project and has no "
                             "meaning in a lymph node",
            "NK_vs_ILC": "the identity core is one combined NK/ILC core, so the LABEL is "
                         "combined. n_nk_like_detected and n_ilc_hint_detected are carried per "
                         "cell as description, from the core's own subtype_hint column",
        },
    }
    (cohort_out(OUT_ROOT, "label_space.json")).write_text(json.dumps(label_space, indent=1))

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                            capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                                capture_output=True, text=True).stdout.strip())
    part_reports = [json.loads(p.read_text())
                    for p in sorted((OUT_ROOT / "partition/section_reports").glob("*.json"))
                    if json.loads(p.read_text())["section_key"] in set(scope.requested)]
    release = {
        "schema_version": SCHEMA_VERSION, "release_id": RELEASE_ID,
        "modality": "xenium", "tissue": "pln", "roi": ROI,
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "latch": "development; nothing here is a release of record",
        "scope": scope.describe(),
        "producer": {
            "lane": LANE,
            "script": "analysis/xenium/scripts/pln_immune_typing_v3.py",
            "commit": commit, "dirty": dirty,
            "supersedes": "data/development/xenium_typing/release_pln_immune_v2 (its partition "
                          "contained nucleus_count == 0 cells, and it stopped at the cluster "
                          "label so NK/ILC inside a T-called cluster was invisible). v2 is left "
                          "in place for comparison and is not overwritten.",
        },
        "method": {
            "partition": "the clustering track's own algorithm "
                         "(cluster_pln.leiden_on_subsample, same seed derivation, same QC) "
                         "re-run on the NUCLEATED cells of each section, into this release's "
                         "partition/ tree. data/development/typing/pln/clusters/xenium/ is not "
                         "written to.",
            "normalisation": f"log1p(raw * {TARGET_SUM:.0f} / depth); the target is FIXED, not "
                             "the section median",
            "detection": DETECTION_PREDICATE + ". This is the ONLY detection predicate in this "
                         "lane: every detection fraction in this release comes from it. The "
                         "number of raw-nonzero entries that fall below it is reported per "
                         "section as entries_nonzero_but_below_predicate, so the divergence "
                         "from a `> 0` predicate is measured rather than described.",
            "standout": f"gene g stands out in cluster k iff its cluster-mean lognorm is >= "
                        f"{Z_STANDOUT} SD above the mean over that section's own clusters AND "
                        "its detection fraction in k exceeds its detection fraction in the rest "
                        "of the section. Both contrasts are internal to one section.",
            "stage_1": f"cluster is Immune iff >= {STAGE1_MIN_STANDOUT} canonical immune genes "
                       "stand out AND a strictly larger FRACTION of the immune set than of the "
                       "structural naming set; Non_immune is the mirror.",
            "stage_2": "within Immune clusters, argmax over the five identity cores of the mean "
                       "standardised enrichment, admitted only when the winning core has >= "
                       f"{STAGE2_MIN_NAMING_STANDOUT} standout gene that MAY NAME its lineage.",
            "stage_3": "the T_NK compartment (cells of clusters called Immune_T or "
                       "Immune_NK_ILC) is split PER CELL by the gene panel: CD3E/CD3G positive "
                       "-> T, then FOXP3 / CD8A+CD8B / CD4; CD3E/CD3G negative with >= "
                       f"{NK_MIN_NAMING} NK-exclusive naming genes -> NK_ILC; the remainder is "
                       "T_NK_unresolved and is counted, not redistributed. The subsets sum to "
                       "the compartment and that is asserted at write time.",
            "stage_3b": "the T_NK compartment is ALSO sub-partitioned inside each section "
                        "(--stage tnk) and each sub-cluster named against the other T_NK "
                        "sub-clusters of that section: CD3E or CD3G STANDS DOWN (the mirror of "
                        "standout: >= 1 SD below the compartment mean AND detected less often "
                        f"than in the rest), >= {NK_MIN_NAMING} NK-exclusive naming genes stand "
                        "up, and the fraction of its cells detecting at least one of them "
                        "exceeds that fraction in the rest of the compartment -> NK_ILC. The "
                        "per-cell gate is a FLOOR (each cell individually corroborated); the "
                        "sub-cluster call is the POPULATION estimate and is the unit the "
                        "PhenoCycler pLN lane also calls at. Both are reported, in "
                        "per_donor_tnk.csv, and neither was tuned to the other. Outputs: "
                        "tnk/tnk_subcluster_profiles.parquet, tnk/assignments/, "
                        "tnk/tnk_section_reports.json.",
            "why_negative_evidence": "no probe on this panel is NK-specific, and none is needed. "
                                     "NK/ILC is the CD3-negative half of a class the panels "
                                     "already carry as combined; the identity core resolves what "
                                     "it is once CD3 has said what it is not.",
            "genes_not_used_for_the_gate": CD3_NOT_USED,
            "naming_vs_supporting": "core/panels' own evidence column decides: {canonical, both} "
                                    "MAY NAME; cellxgene ('marks lineage (CL)') describes only.",
            "nk_evidence_set": "the NK/ILC core's naming genes MINUS every gene that appears "
                               "anywhere in the T core, so no gene whose presence is equally "
                               "explained by a T cell can name NK.",
            "not_licensed": ["any spatial claim across the pancreas/pLN boundary",
                             "cell-to-cell matching across modalities",
                             "acinar or ductal, ever"],
        },
        "panel": {
            "stage1_immune_canonicals": STAGE1_IMMUNE,
            "cd3_gate": CD3_GATE,
            "t_subset_genes": T_SUBSET_GENES,
            "treg_support_genes": TREG_SUPPORT,
            "nk_exclusive_naming_genes": panels["nk_exclusive_naming"],
            "nk_naming_genes_dropped_because_also_in_T_core":
                panels["nk_naming_shared_with_t_core"],
            "nk_like_genes": panels["nk_like_genes"],
            "ilc_hint_genes": panels["nk_ilc_hint_genes"],
            "t_naming_genes": panels["t_naming"],
            "named_in_brief_but_absent_from_every_section": sorted(
                {g for r in reports for g in r["named_in_brief_but_absent"]}),
            "immune_cores": {k: {"file": f"core/panels/data/identity_core/{v}",
                                 "n_genes": len(panels["immune"][k]["genes"]),
                                 "n_naming_genes": len(panels["immune"][k]["naming"])}
                             for k, v in IMMUNE_CORES.items()},
            "control_features_present_in_h5ad_var": sorted(
                {g for r in reports for g in r["control_features_present_in_var"]}) or
                "none: these h5ads carry Gene Expression only, so no control probe was there "
                "to exclude and no output claims one was",
        },
        "controls": {
            "y_probe_genes": Y_GENES, "y_floor_group": Y_FLOOR,
            "y_cross_hybridising_reported_not_control": Y_CROSSHYB,
            "y_probe_note": "in a female donor every Y detection is false; that rate is the "
                            "false-detection rate this predicate admits",
            "female_sections": sorted(r["section_key"] for r in reports
                                      if str(r["sex"]).strip().lower().startswith("f")),
            "cd3_gate_control_note": "the CD3 gate's false-positive floor is its detection rate "
                                     "in a compartment that is definitionally CD3-negative "
                                     "(B/Plasma, Myeloid); its ceiling is its rate in T_NK. "
                                     "Both per section in controls.csv.",
            "gdna_rate_per_donor": {d: gdna.get(d, "UNMEASURED") for d in donors},
            "gdna_source": str(GDNA.relative_to(REPO)) if GDNA.exists() else "UNMEASURED",
        },
        "nucleus_rule": {
            "rule": "a cell with no nucleus is excluded from every calculation -- including the "
                    "partition",
            "handle": "the vendor bundle cells.parquet nucleus_count, matched by donor and by "
                      "subset containment of this section's own cell ids",
            "sections_measured": sorted(r["section_key"] for r in reports
                                        if r["nucleus_rule_measurable"]),
            "sections_unmeasured": sorted(r["section_key"] for r in reports
                                          if not r["nucleus_rule_measurable"]),
            "cells_no_nucleus": int(sum(r["n_cells_no_nucleus"] for r in reports)),
            "no_nucleus_cells_inside_partition": int(
                sum(r["n_no_nucleus_cells_inside_partition"] for r in reports)),
            "cells_removed_from_the_partition_vs_qc_only": int(
                sum(r["n_cells_no_nucleus_inside_qc"] for r in part_reports)),
            "partition_cells_qc_only": int(
                sum(r["n_cells_passing_qc_only"] for r in part_reports)),
            "partition_cells_after_nucleus_rule": int(
                sum(r["n_cells_partitioned"] for r in part_reports)),
        },
        "donors": donors,
        "donors_refused_by_cohort": refused_donors(),
        "sections": {r["section_key"]: r["donor_id"] for r in reports},
        "cells": int(summary["cells"].sum()),
        "cells_partitioned": int(sum(r["n_cells_partitioned"] for r in reports)),
        "t_nk_compartment_cells": int(sum(r["t_nk_compartment_cells"] for r in reports)),
        "duplicate_cell_ids_within_donor": dup_report or "none",
        "inputs_sha256": inputs,
    }
    (cohort_out(OUT_ROOT, "release.json")).write_text(json.dumps(release, indent=1))
    (cohort_out(OUT_ROOT, "section_reports.json")).write_text(json.dumps(reports, indent=1))
    write_scope_note(scope, OUT_ROOT)
    log(f"wrote {cohort_out(OUT_ROOT, 'release.json')}")
    cols = ["donor_id", "cells", "pct_immune", "T_NK_compartment", "T_CD4", "T_CD8", "T_Treg",
            "T_subset_unresolved", "NK_ILC", "T_NK_unresolved", "frac_NK_ILC_of_T_NK"]
    print(summary[cols].to_string(index=False))


# ------------------------------------------------------------------------------------------
# stage: cross -- compositional comparison against the PhenoCycler pLN calls
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# stage: tnk -- sub-partition the T_NK compartment and name the sub-clusters
# ------------------------------------------------------------------------------------------
#
# WHY A SECOND LEVEL EXISTS, AND WHAT EACH ONE IS FOR
# ---------------------------------------------------
# Stage 3's per-cell gate is a FLOOR: every NK_ILC cell it names detects, on that cell, two
# NK-exclusive naming genes and neither CD3E nor CD3G.  At a median of ~200 transcripts per
# cell most real NK cells will not clear that, so the per-cell count undercounts the population
# by construction -- and by how much is not something a per-cell rule can measure.
#
# The POPULATION estimate needs a unit with enough counts to carry a contrast, which is what
# the section-level partition already gives the coarse classes.  So the T_NK compartment is
# sub-partitioned inside each section and each sub-cluster is named against the OTHER T_NK
# sub-clusters of that same section.  This is also the unit the PhenoCycler pLN lane calls at,
# which is why it is the one the cross-modality comparison uses.
#
# Both estimates are reported.  Neither replaces the other and neither was tuned to the other.
# ------------------------------------------------------------------------------------------

def tnk_section(donor: str, section: str, h5ad: Path, panels: dict, mod) -> dict:
    import anndata as ad
    import scanpy as sc

    t0 = time.time()
    calls = pd.read_parquet(
        OUT_ROOT / f"calls/donor_id={donor}/roi={ROI}/data.parquet",
        columns=["cell_id", "section_key", "in_T_NK_compartment", "lineage", "cluster_id",
                 "cd3_detected", "n_nk_exclusive_naming_detected"],
    )
    calls = calls[calls["section_key"] == section]
    tnk_ids = set(calls.loc[calls["in_T_NK_compartment"], "cell_id"].astype(str))
    if len(tnk_ids) < 200:
        return {"donor_id": donor, "section_key": section, "status": "too_few_T_NK_cells",
                "n_T_NK": len(tnk_ids)}

    adata = ad.read_h5ad(h5ad)
    cell_ids = adata.obs_names.to_numpy().astype(str)
    keep = np.array([c in tnk_ids for c in cell_ids])
    adata = adata[keep].copy()
    ids = cell_ids[keep]
    n = int(adata.n_obs)

    raw = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    var_names = adata.var_names.to_numpy().astype(str)

    work = adata.copy()
    sc.pp.normalize_total(work)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=2000)
    work = work[:, work.var["highly_variable"]].copy()
    sc.pp.scale(work, max_value=10)
    feats = np.asarray(work.X, dtype=np.float32)
    del work
    gc.collect()
    seed = mod.stable_seed("pln-tnk-subcluster", CLUSTER_SEED_BASE, "xenium", section)
    labels, meta = mod.leiden_on_subsample(
        feats, n_pcs=30, subsample=CLUSTER_SUBSAMPLE,
        resolution=CLUSTER_RESOLUTION, seed=seed)
    del feats
    gc.collect()
    sub_id = np.array([f"{section}_t{int(v):02d}" for v in labels])

    # the one detection predicate, on these cells, over the genes this stage names with
    gene_idx = {g: i for i, g in enumerate(var_names)}
    wanted = sorted(set(CD3_GATE) | set(panels["t_naming"]) | set(panels["nk_exclusive_naming"])
                    | {g for gs in T_SUBSET_GENES.values() for g in gs}
                    | set(TREG_SUPPORT) | set(CYTOTOXIC_SUPPORT) | {T_SENSITIVITY_GENE})
    present = [g for g in wanted if g in gene_idx]
    cols = np.array([gene_idx[g] for g in present])
    depth = np.asarray(raw.sum(axis=1)).ravel().astype(np.float64)
    block = raw[:, cols].tocsr()
    rows_ix = np.repeat(np.arange(block.shape[0]), np.diff(block.indptr))
    safe = np.where(depth > 0, depth, np.nan)
    scaled = np.nan_to_num(block.data.astype(np.float64) * TARGET_SUM / safe[rows_ix], nan=0.0)
    lognorm = np.log1p(scaled).astype(np.float32)
    detected = (lognorm >= DETECTION_MIN_LOGNORM).astype(np.float32)
    L = sp.csr_matrix((lognorm, block.indices, block.indptr), shape=block.shape)
    D = sp.csr_matrix((detected, block.indices, block.indptr), shape=block.shape)
    del raw, block, scaled, lognorm, detected, rows_ix, adata
    gc.collect()

    order = sorted(set(sub_id))
    lookup = {c: i for i, c in enumerate(order)}
    krow = np.array([lookup[c] for c in sub_id])
    indic = sp.csr_matrix((np.ones(n, dtype=np.float32), (krow, np.arange(n))),
                          shape=(len(order), n))
    sizes = np.asarray(indic.sum(axis=1)).ravel()
    means = np.asarray((indic @ L).todense()) / sizes[:, None]
    dets = np.asarray((indic @ D).todense()) / sizes[:, None]
    det_tot = np.asarray(D.sum(axis=0)).ravel()
    other_n = n - sizes
    with np.errstate(invalid="ignore", divide="ignore"):
        d_other = ((det_tot[None, :] - dets * sizes[:, None])
                   / np.where(other_n > 0, other_n, np.nan)[:, None])
    d_other = np.nan_to_num(d_other, nan=0.0)
    mu = means.mean(axis=0)
    sd = means.std(axis=0, ddof=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.nan_to_num((means - mu[None, :]) / np.where(sd > 0, sd, np.nan)[None, :], nan=0.0)
    so = (z >= Z_STANDOUT) & (dets > d_other)
    #: The MIRROR of standout, and the thing a NEGATIVE call needs.  Inside a compartment that
    #: is mostly T, "CD3 does not stand out" is the default state of nine sub-clusters in ten
    #: and says nothing; a first draft of this stage used it and returned 35% NK on donor 6347,
    #: which is a rule failing, not a finding.  STANDING DOWN is the symmetric statement: the
    #: sub-cluster is >= 1 SD BELOW the compartment's own mean for that gene AND detects it
    #: less often than the rest of the compartment does.  Same contrast, same units, opposite
    #: direction, and it is what "CD3-negative" has to mean when every number is relative.
    sd_ = (z <= -Z_STANDOUT) & (dets < d_other)

    pidx = {g: i for i, g in enumerate(present)}

    def sel(genes, mat=None):
        mat = so if mat is None else mat
        ii = [pidx[g] for g in genes if g in pidx]
        return mat[:, ii].any(axis=1) if ii else np.zeros(len(order), dtype=bool)

    def cnt(genes):
        ii = [pidx[g] for g in genes if g in pidx]
        return so[:, ii].sum(axis=1) if ii else np.zeros(len(order), dtype=int)

    cd3_so = sel(CD3_GATE)
    cd3_sd = sel(CD3_GATE, sd_)
    nk_so_n = cnt(panels["nk_exclusive_naming"])
    t_so_n = cnt(panels["t_naming"])
    cd4_so = sel(T_SUBSET_GENES["CD4"])
    cd8_so = sel(T_SUBSET_GENES["CD8"])
    foxp3_so = sel(T_SUBSET_GENES["Treg"])

    #: The SET-level detection contrast, and the second thing a POSITIVE NK call needs.  The
    #: per-gene standout rule already asks that each named gene be detected more often in k
    #: than in the rest; asked of two genes out of nineteen it can pass in a sub-cluster whose
    #: cells carry LESS NK evidence overall than the compartment average -- which is not thin
    #: evidence, it is contradicted evidence.  A draft of this stage admitted five such
    #: sub-clusters (6579a_t01: 7.2% of its cells detect any NK naming gene against ~15% in the
    #: compartment, on KIR2DL1 and KIR3DL1 alone).  So the same contrast is asked of the SET:
    #: the fraction of k's cells detecting at least one NK-exclusive naming gene must exceed
    #: that fraction in the rest of the compartment.  No new threshold -- the same rule, one
    #: level up -- and the count it removes is reported beside the count it keeps.
    nk_ii = [pidx[g] for g in panels["nk_exclusive_naming"] if g in pidx]
    any_nk = ((np.asarray(D[:, nk_ii].sum(axis=1)).ravel() > 0).astype(np.float64)
              if nk_ii else np.zeros(n))
    nk_set_det = np.array([any_nk[krow == k].mean() for k in range(len(order))])
    tot_nk = float(any_nk.sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        nk_set_det_other = np.where(
            other_n > 0, (tot_nk - nk_set_det * sizes) / np.where(other_n > 0, other_n, 1),
            np.nan)
    nk_set_enriched = nk_set_det > nk_set_det_other

    call, why = [], []
    for k in range(len(order)):
        if (cd3_sd[k] and not cd3_so[k] and nk_so_n[k] >= NK_MIN_NAMING
                and nk_set_enriched[k]):
            call.append("NK_ILC")
            why.append("CD3E or CD3G STANDS DOWN against the rest of this section's T_NK "
                       f"compartment; {int(nk_so_n[k])} NK-exclusive naming gene(s) stand up; "
                       f"and {nk_set_det[k]:.1%} of its cells detect at least one of them "
                       f"against {nk_set_det_other[k]:.1%} in the rest of the compartment")
            continue
        if cd3_so[k] or t_so_n[k] >= 2:
            if foxp3_so[k]:
                call.append("T_Treg"); why.append("FOXP3 stands out")
            elif cd8_so[k] and not cd4_so[k]:
                call.append("T_CD8"); why.append("CD8A/CD8B stand out, CD4 does not")
            elif cd4_so[k] and not cd8_so[k]:
                call.append("T_CD4"); why.append("CD4 stands out, CD8A/CD8B do not")
            elif cd4_so[k] and cd8_so[k]:
                call.append("T_CD4_CD8_both"); why.append("CD4 and a CD8 chain both stand out")
            else:
                call.append("T_subset_unresolved")
                why.append(f"T evidence ({int(t_so_n[k])} naming standouts) with no "
                           "CD4/CD8/FOXP3 standout")
            continue
        call.append("T_NK_unresolved")
        why.append(f"neither side resolves (CD3 stands out {bool(cd3_so[k])}, stands down "
                   f"{bool(cd3_sd[k])}; NK naming standouts {int(nk_so_n[k])}; NK set "
                   f"detection {nk_set_det[k]:.1%} against {nk_set_det_other[k]:.1%} in the "
                   f"rest; T naming standouts {int(t_so_n[k])})")

    call = np.array(call, dtype=object)
    cell_call = call[krow]
    dest = OUT_ROOT / "tnk/assignments" / f"donor_id={donor}" / f"section={section}"
    dest.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": ids, "donor_id": donor, "section_key": section,
                  "tnk_subcluster": sub_id,
                  "tnk_subcluster_call": cell_call.astype(str)}).to_parquet(
        dest / "tnk.parquet", index=False)

    per_cell = calls.set_index("cell_id").loc[ids]
    prof = []
    for k, cid in enumerate(order):
        m = krow == k
        prof.append({
            "donor_id": donor, "section_key": section, "tnk_subcluster": cid,
            "n_cells": int(sizes[k]), "frac_of_T_NK": float(sizes[k] / n),
            "call": call[k], "reason": why[k],
            "cd3_stands_out": bool(cd3_so[k]), "cd3_stands_down": bool(cd3_sd[k]),
            "cd3_detect_frac": float(dets[k, pidx["CD3E"]]) if "CD3E" in pidx else np.nan,
            "cd3_detect_frac_rest_of_T_NK": (float(d_other[k, pidx["CD3E"]])
                                             if "CD3E" in pidx else np.nan),
            "n_nk_exclusive_naming_standouts": int(nk_so_n[k]),
            "nk_set_detection": float(nk_set_det[k]),
            "nk_set_detection_rest_of_T_NK": float(nk_set_det_other[k]),
            "nk_set_enriched": bool(nk_set_enriched[k]),
            "nk_call_thin_two_or_three_genes": bool(call[k] == "NK_ILC"
                                                    and int(nk_so_n[k]) <= 3),
            "n_t_naming_standouts": int(t_so_n[k]),
            "standout_nk_genes": ";".join(g for g in panels["nk_exclusive_naming"]
                                          if g in pidx and so[k, pidx[g]]),
            "standout_t_naming_genes": ";".join(g for g in panels["t_naming"]
                                                if g in pidx and so[k, pidx[g]]),
            "cell_cd3_detect_frac": float(per_cell["cd3_detected"].to_numpy()[m].mean()),
            "cell_frac_ge1_nk_gene": float(
                (per_cell["n_nk_exclusive_naming_detected"].to_numpy()[m] >= 1).mean()),
            "cell_frac_ge2_nk_gene": float(
                (per_cell["n_nk_exclusive_naming_detected"].to_numpy()[m] >= 2).mean()),
            "cell_gate_NK_ILC": int((per_cell["lineage"].to_numpy()[m] == "NK_ILC").sum()),
            "median_depth": float(np.median(depth[m])),
            "coarse_clusters": ";".join(sorted(set(
                str(x) for x in per_cell["cluster_id"].to_numpy()[m]))[:6]),
        })

    nk_sens = {}
    for m_ in NK_SENSITIVITY:
        base = cd3_sd & ~cd3_so & (nk_so_n >= m_)
        nk_sens[str(m_)] = {
            "subclusters": int((base & nk_set_enriched).sum()),
            "cells": int(sizes[base & nk_set_enriched].sum()),
            "without_set_detection_contrast_subclusters": int(base.sum()),
            "without_set_detection_contrast_cells": int(sizes[base].sum()),
        }
    report = {
        "donor_id": donor, "section_key": section, "status": "subclustered",
        "n_T_NK": n, "n_subclusters": len(order),
        "nk_rule_sensitivity_min_standouts": nk_sens,
        "seed": int(seed),
        **{k: v for k, v in meta.items() if k != "n_cells"},
        "calls": {k: int(v) for k, v in pd.Series(cell_call).value_counts().items()},
        "seconds": round(time.time() - t0, 1),
    }
    del L, D
    gc.collect()
    return {"report": report, "profiles": pd.DataFrame(prof)}


def run_tnk(args) -> None:
    panels = load_panels()
    mod = cluster_module()
    secs = sections()
    scope = set_scope(SCOPE.refine([s for _, s, _ in secs]))
    if scope.missing:
        raise SystemExit(f"no pLN h5ad for: {', '.join(scope.missing)}")
    scope.announce(OUT_ROOT, log)
    todo = [row for row in secs if row[1] in set(scope.requested)]

    reports, profiles = [], []
    for donor, section, h5ad in todo:
        if not (OUT_ROOT / f"calls/donor_id={donor}/roi={ROI}/data.parquet").exists():
            raise SystemExit(f"{donor}: no calls -- run --stage call first")
        log(f"T_NK sub-partition {section} (donor {donor})")
        res = tnk_section(donor, section, h5ad, panels, mod)
        if "profiles" not in res:
            reports.append(res)
            log(f"  {res['status']}")
            continue
        reports.append(res["report"])
        profiles.append(res["profiles"])
        r = res["report"]
        log(f"  {r['n_T_NK']:,} T_NK cells, {r['n_subclusters']} sub-clusters -> "
            f"{r['calls']}; {r['seconds']}s")
        gc.collect()

    prof = pd.concat(profiles, ignore_index=True)
    out = cohort_out(OUT_ROOT, "tnk")
    out.mkdir(parents=True, exist_ok=True)
    prof.to_csv(out / "tnk_subcluster_profiles.csv", index=False)
    prof.to_parquet(out / "tnk_subcluster_profiles.parquet", index=False)
    (out / "tnk_section_reports.json").write_text(json.dumps(reports, indent=1))

    rows = []
    for donor in sorted({d for d, _, _ in todo}):
        parts = sorted((OUT_ROOT / "tnk/assignments" / f"donor_id={donor}").glob(
            "section=*/tnk.parquet"))
        if not parts:
            continue
        sub = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
        gate = pd.read_parquet(OUT_ROOT / f"calls/donor_id={donor}/roi={ROI}/data.parquet",
                               columns=["cell_id", "section_key", "lineage",
                                        "in_T_NK_compartment"])
        gate = gate[gate["in_T_NK_compartment"]]
        j = gate.merge(sub, on=["cell_id", "section_key"], how="left")
        g = j["lineage"].value_counts()
        c = j["tnk_subcluster_call"].value_counts()
        row = {"donor_id": donor, "T_NK": int(len(j))}
        for lab in TNK_LABELS:
            row[f"gate__{lab}"] = int(g.get(lab, 0))
            row[f"subcluster__{lab}"] = int(c.get(lab, 0))
        row["gate_frac_NK_ILC"] = row["gate__NK_ILC"] / max(len(j), 1)
        row["subcluster_frac_NK_ILC"] = row["subcluster__NK_ILC"] / max(len(j), 1)
        row["gate_NK_inside_subcluster_NK"] = int(
            ((j["lineage"] == "NK_ILC") & (j["tnk_subcluster_call"] == "NK_ILC")).sum())
        row["subcluster_sum"] = sum(row[f"subcluster__{lab}"] for lab in TNK_LABELS)
        row["coherent"] = row["subcluster_sum"] == len(j)
        rows.append(row)
    summary = pd.DataFrame(rows)
    bad = summary.loc[~summary["coherent"], "donor_id"].tolist()
    if bad:
        raise SystemExit(f"T_NK sub-cluster calls do not sum to the compartment for: {bad}")
    summary.to_csv(cohort_out(OUT_ROOT, "per_donor_tnk.csv"), index=False)
    write_scope_note(scope, OUT_ROOT)
    cols = ["donor_id", "T_NK", "gate__NK_ILC", "subcluster__NK_ILC",
            "gate_NK_inside_subcluster_NK", "gate_frac_NK_ILC", "subcluster_frac_NK_ILC",
            "subcluster__T_CD4", "subcluster__T_CD8", "subcluster__T_Treg",
            "subcluster__T_subset_unresolved", "subcluster__T_NK_unresolved"]
    print(summary[cols].to_string(index=False))


def run_cross(args) -> None:
    """Donor-level COMPOSITIONAL comparison against the PhenoCycler pLN calls.

    Different sections of different blocks, joined on `donor_id` and nothing else.  No cell is
    matched to a cell, no coordinate is compared, and NEITHER SIDE WAS TUNED TO MOVE THE
    AGREEMENT.  What it can show is whether two modalities, typed by two different methods on
    two different molecule classes, put a similar share of the immune compartment in the same
    place -- and where they do not.  A disagreement here is a measurement, not a failure.

    The comparison is made on FRACTIONS, never on counts: the two assays sample different
    section areas and a count comparison would be a comparison of section sizes.  The
    PhenoCycler side's column names are read off its file as they are, and whatever it does not
    carry is reported as absent rather than filled in.
    """
    xen_path = OUT_ROOT / "per_donor_summary.csv"
    if not xen_path.exists():
        raise SystemExit(f"no {xen_path} -- run --stage call first")
    xen = pd.read_csv(xen_path, dtype={"donor_id": str})
    tnk_path = OUT_ROOT / "per_donor_tnk.csv"
    xtnk = (pd.read_csv(tnk_path, dtype={"donor_id": str}) if tnk_path.exists() else None)

    pheno_counts = PHENO_CALLS / "per_donor_counts.csv"
    if not pheno_counts.exists():
        raise SystemExit(f"no PhenoCycler pLN counts at {pheno_counts}")
    ph = pd.read_csv(pheno_counts, dtype={"donor_id": str})
    ph_mtime = datetime.fromtimestamp(pheno_counts.stat().st_mtime,
                                      timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    ph_sha = sha256(pheno_counts)

    PH_TNK = {"T_CD4": "tnk_T_CD4", "T_CD8": "tnk_T_CD8", "T_Treg": "tnk_T_reg",
              "NK_ILC": "tnk_NK_ILC", "T_NK_unresolved": "tnk_T_NK_unresolved"}
    ph_has_tnk = all(c in ph.columns for c in PH_TNK.values())

    rows = []
    for r in ph.itertuples():
        d = {"donor_id": str(r.donor_id),
             "pheno_eligible": int(r.n_eligible),
             "pheno_immune": int(r.n_immune),
             "pheno_pct_immune": float(r.pct_immune)}
        if ph_has_tnk:
            vals = {k: int(getattr(r, v)) for k, v in PH_TNK.items()}
            comp = sum(vals.values())
            d["pheno_T_NK"] = comp
            for k, v in vals.items():
                d[f"pheno_{k}"] = v
            t_res = vals["T_CD4"] + vals["T_CD8"] + vals["T_Treg"]
            d["pheno_T_resolved"] = t_res
            d["pheno_frac_NK_ILC_of_T_NK"] = vals["NK_ILC"] / comp if comp else np.nan
            d["pheno_frac_NK_ILC_of_immune"] = (vals["NK_ILC"] / r.n_immune
                                                if r.n_immune else np.nan)
            d["pheno_frac_T_NK_of_immune"] = comp / r.n_immune if r.n_immune else np.nan
            d["pheno_CD8_over_CD4_plus_CD8"] = (
                vals["T_CD8"] / (vals["T_CD4"] + vals["T_CD8"])
                if (vals["T_CD4"] + vals["T_CD8"]) else np.nan)
            d["pheno_frac_Treg_of_resolved_T"] = vals["T_Treg"] / t_res if t_res else np.nan
        rows.append(d)
    phf = pd.DataFrame(rows)

    xr = []
    for r in xen.itertuples():
        comp = int(r.T_NK_compartment)
        t_res = int(r.T_CD4) + int(r.T_CD8) + int(r.T_Treg)
        d = {"donor_id": str(r.donor_id), "xen_cells": int(r.cells),
             "xen_immune": int(r.immune), "xen_pct_immune": float(r.pct_immune),
             "xen_T_NK": comp, "xen_T_CD4": int(r.T_CD4), "xen_T_CD8": int(r.T_CD8),
             "xen_T_Treg": int(r.T_Treg), "xen_T_CD4_CD8_both": int(r.T_CD4_CD8_both),
             "xen_T_subset_unresolved": int(r.T_subset_unresolved),
             "xen_NK_ILC_cellgate": int(r.NK_ILC),
             "xen_T_NK_unresolved": int(r.T_NK_unresolved),
             "xen_T_resolved": t_res,
             "xen_frac_NK_ILC_of_T_NK_cellgate": r.NK_ILC / comp if comp else np.nan,
             "xen_frac_NK_ILC_of_immune_cellgate": (r.NK_ILC / r.immune
                                                    if r.immune else np.nan),
             "xen_frac_T_NK_of_immune": comp / r.immune if r.immune else np.nan,
             "xen_CD8_over_CD4_plus_CD8": (r.T_CD8 / (r.T_CD4 + r.T_CD8)
                                           if (r.T_CD4 + r.T_CD8) else np.nan),
             "xen_frac_Treg_of_resolved_T": r.T_Treg / t_res if t_res else np.nan}
        xr.append(d)
    xf = pd.DataFrame(xr)
    if xtnk is not None:
        sub = xtnk[["donor_id", "subcluster__NK_ILC", "T_NK"]].rename(
            columns={"subcluster__NK_ILC": "xen_NK_ILC_subcluster",
                     "T_NK": "xen_T_NK_check"})
        xf = xf.merge(sub, on="donor_id", how="left")
        xf["xen_frac_NK_ILC_of_T_NK_subcluster"] = (
            xf["xen_NK_ILC_subcluster"] / xf["xen_T_NK"])
        xf["xen_frac_NK_ILC_of_immune_subcluster"] = (
            xf["xen_NK_ILC_subcluster"] / xf["xen_immune"])

    both = xf.merge(phf, on="donor_id", how="inner").sort_values("donor_id")
    if ph_has_tnk:
        both["delta_frac_NK_ILC_of_T_NK_cellgate"] = (
            both["xen_frac_NK_ILC_of_T_NK_cellgate"] - both["pheno_frac_NK_ILC_of_T_NK"])
        if "xen_frac_NK_ILC_of_T_NK_subcluster" in both:
            both["delta_frac_NK_ILC_of_T_NK_subcluster"] = (
                both["xen_frac_NK_ILC_of_T_NK_subcluster"] - both["pheno_frac_NK_ILC_of_T_NK"])
        both["delta_CD8_over_CD4_plus_CD8"] = (
            both["xen_CD8_over_CD4_plus_CD8"] - both["pheno_CD8_over_CD4_plus_CD8"])
        both["delta_frac_Treg_of_resolved_T"] = (
            both["xen_frac_Treg_of_resolved_T"] - both["pheno_frac_Treg_of_resolved_T"])
        both["delta_frac_T_NK_of_immune"] = (
            both["xen_frac_T_NK_of_immune"] - both["pheno_frac_T_NK_of_immune"])
    both["delta_pct_immune"] = both["xen_pct_immune"] - both["pheno_pct_immune"]

    out_dir = OUT_ROOT / "cross_modality"
    out_dir.mkdir(parents=True, exist_ok=True)
    both.to_csv(out_dir / "per_donor_compositional.csv", index=False)

    def agree(a, b):
        if a not in both or b not in both:
            return "column_absent"
        m = both[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(m) < 3 or m[a].std() == 0 or m[b].std() == 0:
            return {"n": int(len(m)), "pearson": None, "spearman": None,
                    "why": "fewer than 3 usable donors, or one side is constant"}
        return {"n": int(len(m)),
                "pearson": round(float(m[a].corr(m[b])), 4),
                "spearman": round(float(m[a].corr(m[b], method="spearman")), 4),
                "median_abs_difference": round(float((m[a] - m[b]).abs().median()), 4),
                "xenium_median": round(float(m[a].median()), 4),
                "phenocycler_median": round(float(m[b].median()), 4)}

    note = {
        "what_this_is": "a COMPOSITIONAL comparison of two modalities at donor level. The "
                        "PhenoCycler and Xenium pLN sections are different sections of "
                        "different blocks; they are joined on donor_id and on nothing else. No "
                        "cell is matched to a cell and no spatial claim is made. Fractions are "
                        "compared, never counts.",
        "what_this_is_not": "not a validation of either side, and not a target. Neither side "
                            "was tuned to move the agreement, and no number here was accepted "
                            "because it agreed.",
        "xenium": {"release": str(OUT_ROOT.relative_to(REPO)),
                   "donors": int(len(xf)),
                   "nk_estimates": {
                       "cellgate": "per-cell: CD3E and CD3G both undetected on that cell and "
                                   f">= {NK_MIN_NAMING} NK-exclusive naming genes detected on "
                                   "it. A FLOOR -- every such cell is individually "
                                   "corroborated, and most real NK cells will not clear it at "
                                   "this depth.",
                       "subcluster": "per T_NK sub-cluster: CD3E/CD3G stand DOWN against the "
                                     f"rest of the section's T_NK compartment and >= "
                                     f"{NK_MIN_NAMING} NK-exclusive naming genes stand up. The "
                                     "POPULATION estimate, and the one comparable with the "
                                     "PhenoCycler side, which is also a cluster-level call."}},
        "phenocycler": {"source": str(pheno_counts.relative_to(REPO)),
                        "mtime_utc": ph_mtime, "sha256": ph_sha,
                        "donors": int(len(phf)),
                        "t_nk_split_present": bool(ph_has_tnk),
                        "columns_used": (list(PH_TNK.values()) if ph_has_tnk
                                         else "none: that lane carried no T_NK split when this "
                                              "ran, so the subset comparison is one-sided")},
        "donors_in_both": sorted(both.donor_id),
        "n_donors_in_both": int(len(both)),
        "xenium_only": sorted(set(xf.donor_id) - set(phf.donor_id)),
        "phenocycler_only": sorted(set(phf.donor_id) - set(xf.donor_id)),
        "agreement": {
            "pct_immune": agree("xen_pct_immune", "pheno_pct_immune"),
            "frac_T_NK_of_immune": agree("xen_frac_T_NK_of_immune",
                                         "pheno_frac_T_NK_of_immune"),
            "frac_NK_ILC_of_T_NK__xenium_cellgate": agree(
                "xen_frac_NK_ILC_of_T_NK_cellgate", "pheno_frac_NK_ILC_of_T_NK"),
            "frac_NK_ILC_of_T_NK__xenium_subcluster": agree(
                "xen_frac_NK_ILC_of_T_NK_subcluster", "pheno_frac_NK_ILC_of_T_NK"),
            "CD8_share_of_CD4_plus_CD8": agree("xen_CD8_over_CD4_plus_CD8",
                                               "pheno_CD8_over_CD4_plus_CD8"),
            "Treg_share_of_resolved_T": agree("xen_frac_Treg_of_resolved_T",
                                              "pheno_frac_Treg_of_resolved_T"),
        },
        "caveats": [],
        "disagreements": [],
    }
    note["caveats"].append(
        "the two sides resolve their T_NK compartments to different degrees. On the Xenium "
        "side a cell with no CD3E/CD3G transcript and too little NK evidence is "
        "T_NK_unresolved, and that residue is large because the assay is sparse; the "
        "PhenoCycler side's residue is small. Comparisons are therefore made on fractions "
        "WITHIN what each side resolved (CD8 share of CD4+CD8, Treg share of resolved T) as "
        "well as on fractions of the whole compartment, and the two answer different questions.")
    if ph_has_tnk and float(ph[PH_TNK["NK_ILC"]].sum()) == 0:
        note["caveats"].append(
            "the PhenoCycler side reports 0 NK/ILC on every donor in the file read here; the "
            "NK comparison is against a known zero and measures that lane's state.")
    if ph_has_tnk:
        zero = sorted(str(r.donor_id) for r in ph.itertuples()
                      if int(getattr(r, PH_TNK["NK_ILC"])) == 0)
        if zero:
            note["caveats"].append(
                f"the PhenoCycler side reports 0 NK/ILC on {len(zero)} of its donors "
                f"({', '.join(zero)}); a zero there is a lane state, not a measured absence, "
                "and those donors dominate the disagreement below.")
    note["caveats"].append(
        "the PhenoCycler pLN lane and this one share no threshold, no gene, no marker and no "
        "code. That is why the comparison is worth making and also why a difference cannot be "
        "attributed to one side without a further measurement.")

    for r in both.itertuples():
        items = []
        if pd.notna(r.delta_pct_immune) and abs(r.delta_pct_immune) >= 10:
            items.append(f"immune share of all cells differs by {r.delta_pct_immune:+.1f} "
                         f"points (xenium {r.xen_pct_immune:.1f}%, phenocycler "
                         f"{r.pheno_pct_immune:.1f}%)")
        for lab, col, xc, pc in [
            ("NK/ILC share of T_NK (xenium sub-cluster)",
             "delta_frac_NK_ILC_of_T_NK_subcluster",
             "xen_frac_NK_ILC_of_T_NK_subcluster", "pheno_frac_NK_ILC_of_T_NK"),
            ("CD8 share of CD4+CD8", "delta_CD8_over_CD4_plus_CD8",
             "xen_CD8_over_CD4_plus_CD8", "pheno_CD8_over_CD4_plus_CD8"),
            ("Treg share of resolved T", "delta_frac_Treg_of_resolved_T",
             "xen_frac_Treg_of_resolved_T", "pheno_frac_Treg_of_resolved_T"),
            ("T_NK share of immune", "delta_frac_T_NK_of_immune",
             "xen_frac_T_NK_of_immune", "pheno_frac_T_NK_of_immune"),
        ]:
            if col not in both.columns:
                continue
            v = getattr(r, col, np.nan)
            if pd.notna(v) and abs(v) >= 0.10:
                items.append(f"{lab} differs by {v:+.3f} (xenium {getattr(r, xc):.3f}, "
                             f"phenocycler {getattr(r, pc):.3f})")
        if items:
            note["disagreements"].append({"donor_id": r.donor_id, "items": items})
    (out_dir / "cross_modality.json").write_text(json.dumps(note, indent=1))

    cols = [c for c in ["donor_id", "xen_pct_immune", "pheno_pct_immune",
                        "xen_frac_T_NK_of_immune", "pheno_frac_T_NK_of_immune",
                        "xen_frac_NK_ILC_of_T_NK_cellgate",
                        "xen_frac_NK_ILC_of_T_NK_subcluster", "pheno_frac_NK_ILC_of_T_NK",
                        "xen_CD8_over_CD4_plus_CD8", "pheno_CD8_over_CD4_plus_CD8",
                        "xen_frac_Treg_of_resolved_T", "pheno_frac_Treg_of_resolved_T"]
            if c in both.columns]
    print(both[cols].round(4).to_string(index=False))
    print(f"\ndonors in both: {len(both)}  xenium only: {note['xenium_only']}  "
          f"phenocycler only: {note['phenocycler_only']}")
    print(json.dumps(note["agreement"], indent=1))
    for d in note["disagreements"]:
        print(d["donor_id"], "::", "; ".join(d["items"]))


def run_compare(args) -> None:
    """Cross-tab v3 against v2 on the cells both cover."""
    rows = []
    out_dir = OUT_ROOT / "compare_v2"
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(OUT_ROOT.glob("calls/donor_id=*/roi=*/data.parquet")):
        donor = path.parent.parent.name.split("=", 1)[1]
        old = V2 / f"calls/donor_id={donor}/roi={ROI}/data.parquet"
        if not old.exists():
            rows.append({"donor_id": donor, "in_v2": False})
            continue
        mine = pd.read_parquet(path, columns=["cell_id", "stage1", "cluster_lineage", "lineage"])
        prev = pd.read_parquet(old, columns=["cell_id", "stage1", "lineage"])
        j = mine.merge(prev, on="cell_id", how="inner", suffixes=("_v3", "_v2"))
        rows.append({
            "donor_id": donor, "in_v2": True,
            "v3_cells": len(mine), "v2_cells": len(prev), "shared_cells": len(j),
            "stage1_agree": int((j["stage1_v3"] == j["stage1_v2"]).sum()),
            "cluster_lineage_agree": int((j["cluster_lineage"] == j["lineage_v2"]).sum()),
            "v2_NK_ILC": int((j["lineage_v2"] == "Immune_NK_ILC").sum()),
            "v3_NK_ILC": int((j["lineage_v3"] == "NK_ILC").sum()),
        })
        pd.crosstab(j["lineage_v2"], j["lineage_v3"]).to_csv(out_dir / f"donor_{donor}.csv")
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "summary.csv", index=False)
    print(frame.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage",
                    choices=["partition", "call", "tnk", "cross", "compare"], default="call")
    ap.add_argument("--sections", default="",
                    help="comma-separated section keys; a partial run never writes a "
                         "cohort-wide path")
    ap.add_argument("--force", action="store_true",
                    help="redo a section that is already partitioned")
    args = ap.parse_args()
    requested = tuple(s for s in args.sections.split(",") if s.strip())
    set_scope(RunScope.undiscovered("section", requested))
    if args.stage == "partition":
        run_partition(args)
    elif args.stage == "call":
        run_calls(args)
    elif args.stage == "tnk":
        run_tnk(args)
    elif args.stage == "cross":
        run_cross(args)
    else:
        run_compare(args)


if __name__ == "__main__":
    main()
