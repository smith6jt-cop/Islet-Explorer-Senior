#!/usr/bin/env python3
"""Within-section clustering and cluster naming for the pancreas Xenium 5K sections.

This is the worker module behind ``panc_xenium_release.py``, which runs it once per section in
its own process and then assembles the release from what it wrote.  One section is loaded,
clustered and named in a single pass, and the per-section artifacts are written to a work
directory so the cohort run is resumable and can be fanned out across processes.

METHOD
------
1.  CLUSTER WITHIN DONOR, PER SECTION.  A pancreas section is one slide; nothing is pooled and
    no joint embedding is built.  Normalisation, HVG selection, scaling, PCA, the Leiden
    partition and every threshold below are computed from that section's own cells only.

2.  NAME EACH CLUSTER BY CANONICAL ENRICHMENT AGAINST THE OTHER CLUSTERS OF THE SAME SECTION.
    A gene STANDS OUT in cluster k iff its cluster-mean lognorm is >= Z_STANDOUT SD above the
    mean over that section's own clusters AND its detection fraction in k exceeds its detection
    fraction in the rest of the section.  Both halves are section-internal contrasts.

3.  RESOLVE WITH THE ``core/panels`` IDENTITY CORES.  Score of a core = mean standardised
    enrichment (z) over THE GENES OF THAT CORE THAT MAY NAME THE LINEAGE -- the panel's own
    ``evidence`` column, values ``canonical`` or ``both``.  The winner is the argmax, admitted
    only when one of those naming genes actually stands out.  A cluster that fails is
    ``Uncalled_no_naming_evidence`` -- never reassigned to the runner-up, never to a larger core.

    WHY THE SCORE IS NAMING-ONLY, AND HOW THAT WAS DECIDED.  Scoring the mean over the WHOLE core
    is biased by core size, because the cores are not the same size and the parts that differ in
    size are the ``cellxgene`` supporting genes, which are broadly expressed.  03_immune_T_cell
    carries 101 genes of which 20 may name; 04_immune_B_plasma carries 58 of which 26; 07_immune
    _granulocyte carries 21 of which 19.  A mean over the whole core is therefore pulled towards
    the section background in proportion to the supporting-gene count.  Measured on section 6539,
    at each core's own best cluster, mean z over naming genes minus mean z over the whole core:
    Immune_T (81 supporting genes) 1.033, Vascular (64) 0.532, Fibroblast_Stellate (29) 0.688,
    Immune_B_Plasma (32) 0.344, Immune_NK_ILC (9) 0.233, Pericyte_SMC (1) 0.190, Endocrine (2)
    0.187, Granulocyte (2) 0.086, Myeloid (59) -0.260.  The consequence was a wrong call: in
    6539's lymphocyte cluster c16 the whole-core score put B_plasma 2.875 over Immune_T 2.628 and
    the cluster went out as B, while ALL TWENTY T naming genes stood out in it (CD3E detected in
    37.2% of its cells against 1.1% of the rest of the section, CD8A 30.5% against 0.3%) and the
    naming-only score puts Immune_T 3.661 over B_plasma 3.219.

    WHY IT IS NOT PER-CORE NORMALISATION.  The alternative considered was standardising each
    core's score against its own distribution across that section's clusters, so the cores compete
    on a common scale.  Measured on 6539 (robust centre and MAD across the 19 clusters) it makes
    things worse, not better: the denominator is small for a core absent from almost every cluster
    of the section, so the winner becomes whichever core happens to have the tightest across-
    cluster spread -- an artifact of that section's composition, not evidence.  It moved 6539_c04,
    a 307-cell cluster carrying six standout granulocyte naming genes (CTSG, ENPP3, GATA2, IL5RA,
    PRTN3, RNASE3), to Immune_T, and a 767-cell exocrine cluster to B_plasma.  It replaces a bias
    towards small cores with a bias towards narrow ones.  Rejected on those measurements.

    The whole-core mean is still computed and written out per cluster as ``full_core_score__*``,
    so the difference between the two is visible on every cluster rather than asserted here.

4.  CELLS INHERIT THE CLUSTER CALL.  Per-cluster heterogeneity is written out so a mixed
    cluster is visible rather than hidden behind its majority label.

5.  NEGATIVE EVIDENCE IS EVIDENCE.  A cluster whose argmax is Immune_T but in which neither
    CD3E nor CD3G stands out, while at least one NK naming gene does, is called
    Immune_NK_ILC.  The move is recorded on every cluster of both partitions and counted in the
    report SEPARATELY FOR EACH, because the labels that ship come from the final subcluster
    partition and the parent partition is a diagnostic no released label is drawn from.  It is a
    combination of measurements on hand, not a fallback.

WHAT 16_epithelial_general IS AND IS NOT
----------------------------------------
It is an embedding extra: its four genes (KRT19, EPCAM, CDH1, TP63 as the panel carries them)
are measured, reported per cluster and used nowhere as a lineage.  KRT/EPCAM are shared by
ductal, acinar and endocrine cells, so scoring them as one lineage biases the argmax.  It is
never a broad mask and Exocrine is never split: 02_pancreas_exocrine is ONE label whose genes
happen to carry ductal and acinar subtype hints, and no output names acinar or ductal.

ENDOCRINE
---------
The 5K panel carries no INS, GCG, SST or PPY -- measured absent from the var index of every
section and recorded in the report.  Endocrine is resolved only through the curated surrogate
identity panel 01_pancreas_endocrine (18 genes, 16 of them naming).  That call supports "this
cluster is endocrine"; it does not support any alpha / beta / delta statement, and none is made.

DETECTION
---------
``lognorm >= log1p(1)`` with lognorm = ``log1p(raw * 1e4 / depth)``, depth being that cell's own
total over all panel genes.  The target is FIXED at 1e4, never the section median, so for every
cell of depth <= 10,000 the predicate is exactly ``raw >= 1``.  No raw-count floor operates and
none is depth-scaled.  The number of cells above that depth is measured per section.

CONTROLS
--------
Y probes.  TBL1Y / USP9Y / EIF1AY are the FLOOR GROUP: zero-template in a female donor and,
unlike DDX3Y and KDM5D, not cross-hybridising with a highly expressed X paralog.  In a female
donor every detection of those is false, so the rate is a direct read of the background this
predicate admits.  Transcript depth is reported per section as the median transcripts per cell,
because a 49-transcript cell is weak substrate for a 5K-panel call whatever its control rate.
"""
from __future__ import annotations

import gc
import hashlib
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
H5AD_DIR = Path("/data/xenium/h5ad_out")
IDENTITY_CORE_FLAT = REPO / "core/panels/data/identity_core"
IDENTITY_CORE_PANC = REPO / "core/panels/data/tissues/pancreas/identity_core"
XENIUM_ROOT = Path("/data/xenium")
XENIUM_PATHS = REPO / "core/phenocycler/data/integration/xenium_paths.csv"

ROI = "panc"
TARGET_SUM = 1.0e4
DETECTION_MIN_LOGNORM = float(np.log1p(1.0))
Z_STANDOUT = 1.0
Z_SENSITIVITY = (0.75, 1.0, 1.5)
MIN_NAMING_STANDOUT = 1
SHARED_SECRETORY_NAMING_GENES = frozenset({"XBP1", "FKBP11", "DERL3"})
NAMING_POLICIES = ("native", "shared_er_exact_v1")


def naming_policy_panels(panels: dict, policy: str) -> dict:
    """Keep shared secretory genes measured as support, never as B identity.

    The native policy is retained for provenance/reproduction. The selectable
    repair changes only B naming-score/admission roles; embedding features and
    full-core support diagnostics are unchanged.
    """
    if policy not in NAMING_POLICIES:
        raise ValueError(f"unknown naming policy: {policy}")
    if policy == "native":
        return panels
    out = dict(panels)
    out["cores"] = dict(panels["cores"])
    genes, naming = panels["cores"]["Immune_B_Plasma"]
    out["cores"]["Immune_B_Plasma"] = (genes, [g for g in naming if g not in SHARED_SECRETORY_NAMING_GENES])
    return out


def exact_detection_q(dets: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Section-wide one-sided cluster-versus-rest detection test, BH corrected.

    Every group × every measured naming-basis feature belongs to the declared
    family. Clusters were selected from expression, so this is an exploratory
    contrast after selection, not a guaranteed lineage false-discovery rate.
    Statistical significance remains separate from program coherence.
    """
    from scipy.stats import hypergeom
    sizes = np.asarray(sizes, dtype=np.int64)
    if sizes.shape != (dets.shape[0],) or (sizes <= 0).any():
        raise ValueError("Exact detection contrast requires positive exact group sizes")
    if not np.isfinite(dets).all() or (dets < 0).any() or (dets > 1).any():
        raise ValueError("Detection fractions must be finite and lie between zero and one")
    detected = np.rint(dets * sizes[:, None]).astype(np.int64)
    # Fresh native statistics use float32. Valid ratios in large parents can
    # differ from their exact integer by >.01 purely through float32 rounding.
    # The tolerance follows the input precision and stays strictly below half a
    # count, so nearest-integer recovery remains unambiguous.
    precision = np.finfo(dets.dtype).eps
    tolerance = np.minimum(.49, np.maximum(.01, precision * sizes[:, None]))
    if (np.abs(detected - dets * sizes[:, None]) > tolerance).any():
        raise ValueError("Stored fractions do not recover exact integer detection counts")
    total = detected.sum(axis=0)
    p = hypergeom.sf(detected - 1, sizes.sum(), total[None, :], sizes[:, None])
    p[sizes.sum() - sizes == 0] = 1
    flat = p.ravel()
    order = np.argsort(flat, kind="stable")
    ranked = flat[order]
    values = np.minimum.accumulate((ranked * len(flat) / np.arange(1, len(flat) + 1))[::-1])[::-1]
    q = np.empty_like(values)
    q[order] = np.minimum(values, 1)
    return q.reshape(dets.shape)

#: cluster QC, identical in spirit to the pLN clustering track: a cell with almost nothing
#: measured cannot be partitioned on its measurements.  It is carried in the release as
#: Uncalled_below_cluster_qc, never as 0 and never as a lineage.
QC_MIN_COUNTS = 10
QC_MIN_FEATURES = 5

N_HVG = 2000
N_PCS = 30
PCA_FIT_MAX = 250_000
LEIDEN_SUBSAMPLE = 60_000
LEIDEN_RESOLUTION = 2.0
KNN_K = 15

#: Every parent cluster of at least this many cells is partitioned a second time and named
#: against its OWN subclusters. It is applied to every lineage, not only the immune ones, so
#: no label is given a resolution the others do not get.
REFINE_MIN_CELLS = 2_000
REFINE_RESOLUTION = 1.0

#: The lineage cores.  Every one of these MAY produce a label.
LINEAGE_CORES = {
    "Endocrine":                 (IDENTITY_CORE_PANC / "01_pancreas_endocrine.csv",   "Endocrine"),
    "Exocrine":                  (IDENTITY_CORE_PANC / "02_pancreas_exocrine.csv",    "Exocrine"),
    "Immune_T":                  (IDENTITY_CORE_FLAT / "03_immune_T_cell.csv",        "T_NK"),
    "Immune_B_Plasma":           (IDENTITY_CORE_FLAT / "04_immune_B_plasma.csv",      "B_Plasma"),
    "Immune_NK_ILC":             (IDENTITY_CORE_FLAT / "05_immune_NK_ILC.csv",        "T_NK"),
    "Myeloid":                   (IDENTITY_CORE_FLAT / "06_immune_myeloid.csv",       "Myeloid"),
    "Granulocyte":               (IDENTITY_CORE_FLAT / "07_immune_granulocyte.csv",   "Granulocyte"),
    "Vascular":                  (IDENTITY_CORE_FLAT / "11_endothelial_vascular.csv", "Vascular"),
    "Stromal_Pericyte_SMC":      (IDENTITY_CORE_FLAT / "12_pericyte_smooth_muscle.csv", "Stromal"),
    "Stromal_Fibroblast_Stellate": (IDENTITY_CORE_PANC / "13_fibroblast_stellate.csv", "Stromal"),
    "Neural":                    (IDENTITY_CORE_FLAT / "15_neural_glial.csv",         "Neural"),
}
#: Measured, reported, and NOT a lineage.  See the module docstring.
EMBEDDING_EXTRA = {"epithelial_general": IDENTITY_CORE_FLAT / "16_epithelial_general.csv"}

NAMING_EVIDENCE = frozenset({"canonical", "both"})

#: Which genes of a core carry its score.  "naming" is the decided rule (see the module
#: docstring); "all" is the superseded whole-core mean, kept runnable ONLY so the before/after
#: can be produced from the same code on the same partition.
SCORE_ON_DEFAULT = "naming"

#: Internal coherence, not an expectation check.  If the genes that MAY NAME a lineage are
#: detected in a section and that lineage still receives no cells, the scoring is at fault and
#: the section is flagged for investigation.  This compares a section against itself; it is not
#: a check against any prior biological finding, and no target fraction exists for any of them.
COHERENCE_GENES = ["PTPRC", "CD3E", "CD3G", "CD4", "CD8A", "NKG7", "KLRD1", "MS4A1", "CD68"]

Y_FLOOR = ["TBL1Y", "USP9Y", "EIF1AY"]
Y_CROSSHYB = ["DDX3Y", "KDM5D"]
Y_GENES = Y_FLOOR + Y_CROSSHYB

#: absent from the 5K panel; checked per section and reported rather than assumed
ENDOCRINE_HORMONES = ["INS", "GCG", "SST", "PPY", "GHRL"]
EXOCRINE_ENZYMES = ["PRSS1", "CPA1", "CELA3A", "CTRB1", "CTRB2", "AMY2A", "REG1A"]

#: negative-evidence redirect, T -> NK
T_DEFINING = ["CD3E", "CD3G"]
NK_NAMING_HINT = ["NKG7", "KLRD1", "GNLY", "KLRF1", "NCR1", "PRF1"]
#: The reason string a cluster carries when the redirect produced its label.  One authority:
#: ``decide`` writes it, ``redirect_counts`` counts it and panc_xenium_release.py reads it off
#: the shipped calls.  A redirect that was APPLIED is not the same event as one that SHIPPED a
#: label -- the moved cluster can still fail the enrichment gate below and go out uncalled.
REDIRECT_REASON = "cd3_negative_nk_positive_redirect"

CONTROL_PREFIXES = ("NegControl", "Intergenic", "Unassigned", "Deprecated", "BLANK", "antisense")


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def stable_seed(*parts: object) -> int:
    import hashlib
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()
    return int(h[:8], 16)


def read_core(path: Path) -> tuple[list[str], list[str]]:
    frame = pd.read_csv(path)
    genes = [str(g) for g in frame["gene"].dropna().unique()]
    naming = sorted({
        str(r.gene) for r in frame.itertuples()
        if str(getattr(r, "evidence", "")).strip() in NAMING_EVIDENCE
    })
    return genes, naming


def load_panels() -> dict:
    cores = {name: read_core(path) for name, (path, _) in LINEAGE_CORES.items()}
    extra = {name: read_core(path) for name, path in EMBEDDING_EXTRA.items()}
    wanted = sorted(
        {g for genes, _ in cores.values() for g in genes}
        | {g for genes, _ in extra.values() for g in genes}
        | set(Y_GENES) | set(T_DEFINING) | set(NK_NAMING_HINT) | set(COHERENCE_GENES)
    )
    return {"cores": cores, "extra": extra, "wanted": wanted}


def sections() -> list[tuple[str, str, Path]]:
    """(donor_id, section_key, h5ad) for every pancreas h5ad of a donor whose Xenium is eligible.

    The exclusion is asked of ``multiplex5000.cohort.excluded_on(donor, "xenium")``, never
    restated: 6457's central exclusion covers its PhenoCycler data only (maintainer,
    2026-10-22), so its Xenium pancreas is typed here.
    """
    import sys
    sys.path.insert(0, str(REPO / "core/multiplex5000"))
    from multiplex5000.cohort import excluded_on

    out = []
    for path in sorted(H5AD_DIR.glob("*_PancXenium.h5ad")):
        stem = path.name.replace("_PancXenium.h5ad", "")
        core = stem[3:] if stem.startswith("HDL") else stem
        donor = core.rstrip("abc") if core[-1:].isalpha() else core
        if excluded_on(donor, "xenium"):
            continue
        out.append((donor, stem, path))
    return out


def nucleus_for(donor: str, cell_ids: np.ndarray) -> tuple[np.ndarray | None, str]:
    """``nucleus_count`` per h5ad cell from the vendor bundle, or (None, why-not).

    A cell with no nucleus is excluded from every calculation.  The h5ads carry no nucleus
    geometry, so the only handle is the bundle ``cells.parquet``, matched by donor AND by subset
    containment of the h5ad's own cell ids -- never by filename order.  Where no local bundle is
    a superset the rule is UNMEASURABLE and recorded as such, never as passed.
    """
    table = pd.read_csv(XENIUM_PATHS, dtype=str)
    table = table[
        table["roi"].astype(str).str.startswith("panc")
        & (table["donor_id"].astype(str) == str(donor))
    ]
    want = set(cell_ids.tolist())
    tried = []
    for row in table.itertuples():
        src = str(row.source_path or "").strip()
        if not src or src == "NA":
            continue
        cells = XENIUM_ROOT / src / "cells.parquet"
        if not cells.is_file():
            tried.append(f"{row.batch}:absent")
            continue
        frame = pd.read_parquet(cells, columns=["cell_id", "nucleus_count"])
        frame["cell_id"] = frame["cell_id"].astype(str)
        lookup = dict(zip(frame["cell_id"], frame["nucleus_count"].to_numpy()))
        missing = len(want - lookup.keys())
        if missing:
            tried.append(f"{row.batch}:not_a_superset({missing} ids missing)")
            continue
        counts = np.array([lookup[c] for c in cell_ids], dtype=np.int32)
        return counts, f"bundle {cells}"
    return None, "UNMEASURED: no local bundle cells.parquet is a superset of this h5ad" + (
        f" (tried {'; '.join(tried)})" if tried else ""
    )


def _standardised_chunks(Xh: sp.csr_matrix, mu: np.ndarray, sd: np.ndarray, step: int):
    for start in range(0, Xh.shape[0], step):
        block = Xh[start:start + step].toarray()
        block -= mu
        block /= sd
        np.clip(block, -10.0, 10.0, out=block)
        yield start, block


def cluster_section(lognorm: sp.csr_matrix, seed: int, resolution: float = LEIDEN_RESOLUTION) -> tuple[np.ndarray, dict]:
    """Leiden on a PCA subsample, carried to every cell by a 15-NN vote in the same space.

    Leiden needs a neighbour graph and a graph over a million cells is not affordable per
    section, so the partition is learned on a random subsample and assigned outward.  The cost
    of that shortcut is MEASURED, not assumed: the subsample is split 80/20, the classifier is
    refit on the 80 and its agreement with Leiden's own labels on the held-out 20 is reported.
    The PCA basis is likewise fitted on a subsample and the explained variance is reported.
    """
    import anndata as ad
    import scanpy as sc
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KNeighborsClassifier

    rng = np.random.default_rng(seed)
    n_cells = lognorm.shape[0]

    hv = ad.AnnData(X=lognorm)
    sc.pp.highly_variable_genes(hv, n_top_genes=min(N_HVG, lognorm.shape[1] - 1), flavor="seurat")
    hvg_idx = np.flatnonzero(hv.var["highly_variable"].to_numpy())
    del hv
    gc.collect()

    Xh = lognorm[:, hvg_idx].astype(np.float32)
    mu = np.asarray(Xh.mean(axis=0)).ravel().astype(np.float32)
    sq = np.asarray(Xh.multiply(Xh).mean(axis=0)).ravel().astype(np.float64)
    sd = np.sqrt(np.maximum(sq - mu.astype(np.float64) ** 2, 0.0)).astype(np.float32)
    sd[sd <= 0] = 1.0

    n_pcs = int(min(N_PCS, len(hvg_idx) - 1, n_cells - 1))
    fit_n = min(PCA_FIT_MAX, n_cells)
    fit_idx = np.sort(rng.choice(n_cells, size=fit_n, replace=False)) if fit_n < n_cells \
        else np.arange(n_cells)
    fit_block = Xh[fit_idx].toarray()
    fit_block -= mu
    fit_block /= sd
    np.clip(fit_block, -10.0, 10.0, out=fit_block)
    pca = PCA(n_components=n_pcs, random_state=seed, svd_solver="randomized")
    pca.fit(fit_block)
    del fit_block
    gc.collect()

    coords = np.empty((n_cells, n_pcs), dtype=np.float32)
    for start, block in _standardised_chunks(Xh, mu, sd, 100_000):
        coords[start:start + block.shape[0]] = pca.transform(block).astype(np.float32)
        del block
    del Xh
    gc.collect()

    take = min(LEIDEN_SUBSAMPLE, n_cells)
    idx = np.sort(rng.choice(n_cells, size=take, replace=False)) if take < n_cells \
        else np.arange(n_cells)
    sub = ad.AnnData(X=coords[idx].copy())
    sub.obsm["X_pca"] = coords[idx].copy()
    sc.pp.neighbors(sub, n_neighbors=KNN_K, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(sub, resolution=resolution, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False, random_state=seed)
    sub_labels = sub.obs["leiden"].to_numpy().astype(str)
    del sub
    gc.collect()

    perm = rng.permutation(take)
    cut = int(take * 0.8)
    probe = KNeighborsClassifier(n_neighbors=KNN_K, n_jobs=8)
    probe.fit(coords[idx][perm[:cut]], sub_labels[perm[:cut]])
    holdout = float((probe.predict(coords[idx][perm[cut:]]) == sub_labels[perm[cut:]]).mean())
    del probe

    labels = np.empty(n_cells, dtype=object)
    labels[idx] = sub_labels
    rest = np.setdiff1d(np.arange(n_cells), idx, assume_unique=True)
    if rest.size:
        clf = KNeighborsClassifier(n_neighbors=KNN_K, n_jobs=8)
        clf.fit(coords[idx], sub_labels)
        for start in range(0, rest.size, 200_000):
            chunk = rest[start:start + 200_000]
            labels[chunk] = clf.predict(coords[chunk])
        del clf
    del coords
    gc.collect()

    meta = {
        "n_cells_clustered": int(n_cells),
        "n_hvg": int(len(hvg_idx)),
        "n_pcs": n_pcs,
        "pca_fitted_on_cells": int(fit_n),
        "pca_variance_explained": float(pca.explained_variance_ratio_.sum()),
        "leiden_subsample": int(take),
        "leiden_resolution": resolution,
        "n_clusters": int(len(set(sub_labels))),
        "knn_assignment_k": KNN_K,
        "knn_holdout_agreement": holdout,
        "cells_assigned_by_knn": int(rest.size),
        "seed": int(seed),
    }
    return labels.astype(str), meta


def partition_stats(L, D, krow, n_groups):
    """Cluster means, detection fractions, the rest-of-partition detection fraction, and z.

    Every quantity here is a contrast INSIDE one partition of one section. `z` standardises a
    gene's cluster means across that partition's own clusters; `d_other` is its detection
    fraction in everything else in the same partition. No level, quantile or threshold crosses
    a slide, and none crosses a tissue.
    """
    n = L.shape[0]
    indic = sp.csr_matrix(
        (np.ones(n, dtype=np.float32), (krow, np.arange(n))), shape=(n_groups, n)
    )
    sizes = np.asarray(indic.sum(axis=1)).ravel()
    means = np.asarray((indic @ L).todense()) / sizes[:, None]
    dets = np.asarray((indic @ D).todense()) / sizes[:, None]
    det_tot = np.asarray(D.sum(axis=0)).ravel()
    other_n = n - sizes
    with np.errstate(invalid="ignore", divide="ignore"):
        d_other = (det_tot[None, :] - dets * sizes[:, None]) / \
            np.where(other_n > 0, other_n, np.nan)[:, None]
    d_other = np.nan_to_num(d_other, nan=0.0)
    mu = means.mean(axis=0)
    sd = means.std(axis=0, ddof=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (means - mu[None, :]) / np.where(sd > 0, sd, np.nan)[None, :]
    return means, dets, d_other, np.nan_to_num(z, nan=0.0), sizes


def decide(z, dets, d_other, thr, core_names, score_cols, core_naming_cols, t_cols, nk_hint_cols,
           *, exact_detection_fdr=None, cluster_sizes=None):
    """Name every cluster of one partition by canonical enrichment against its siblings.

    A gene STANDS OUT in cluster k iff its cluster-mean lognorm is >= `thr` SD above the mean
    over that partition's own clusters AND its detection fraction in k exceeds its detection
    fraction in the rest of the partition.

    A core's SCORE is the mean z over `score_cols[core]` -- the genes of that core that may name
    its lineage.  Scoring the whole core instead biases the argmax by core size; the measurement
    behind that is in the module docstring.

    A cluster takes the argmax core's label only when all three hold:
      * the winning score is above the partition mean at all (z > 0) -- otherwise nothing is
        enriched here and the argmax is picking the least-negative of eleven negatives;
      * at least one gene that MAY NAME the winning lineage stands out (core/panels' own
        `evidence` column, values canonical or both);
      * that naming gene is detected above the cluster's OWN background: `bg` is the MEDIAN
        of the cluster's detection fractions over the measured set -- the 485 identity-core
        genes carried in this section, not the whole 5K panel -- and the gene must clear it
        strictly, `dets > bg`.

    WHERE THIS GUARD IS DEGENERATE.  `bg` is zero whenever at most half of those 485 genes are
    detected in any cell of the cluster (measured: every zero-background cluster in the release
    has 242 or fewer of the 485 detected, median 177).  The test then reduces to `dets > 0`, so
    the naming gene need be seen in ONE cell -- for those clusters the guard admits exactly what
    a background is supposed to exclude, and the claim it once carried here, that a near-empty
    cluster cannot be named by a gene seen in eight of its cells, was false for them.
    Recomputed 2026-09-18 from `_work/*.gene_stats.parquet` and `_work/*.subclusters.parquet`
    of `/data/xenium_typing/panc_2026-09-17/` (21 sections; the reconstruction reproduces
    `n_naming_above_background_of_winning_core` exactly for all 4,914 labelled subclusters):
      * 544 of 4,914 labelled subclusters (11.1 %) have `bg == 0`, covering 53,333 called cells
        (0.35 % of the 15,318,868 cells in labelled subclusters).  Per section the share of
        labelled subclusters with a zero background runs 3.5 % (6593, 10 of 289) to 29.9 %
        (6523, 61 of 204); 6539 is 22.1 % and 6380 20.3 %.
      * they are small: median 61 cells against 1,924 for the labelled subclusters whose
        background is above zero.
      * 72 of the 544 (10,354 cells) were admitted on a SINGLE naming gene above that zero
        background, and for 37 of those (4,706 cells, largest 563) the one gene was detected in
        exactly one cell -- e.g. 6523_c19_s07, 563 cells, Neural on TH in 1 cell; 6523_c06_s01,
        410 cells, Granulocyte on CXCR2 in 1 cell; HDL112_c10_s08, 360 cells, Granulocyte on
        PRG2 in 1 cell.
      * across all 544, 185 carry an admitted naming gene detected in exactly one cell, and 86
        (11,549 cells) have no admitted naming gene detected in more than two.
    A non-degenerate floor is a METHOD DECISION and is deliberately not taken in this function;
    the guard as written is a no-op above zero for those 544.
    A cluster that fails is uncalled, with the reason recorded. It is NEVER reassigned to the
    runner-up and never to a larger core.

    NEGATIVE EVIDENCE IS EVIDENCE: a cluster whose argmax is Immune_T but in which neither CD3E
    nor CD3G stands out, while NK genes do and an NK naming gene clears the bar, is NK.
    """
    n_groups = z.shape[0]
    so = (z >= thr) & (dets > d_other)
    bg = np.median(dets, axis=1)
    above_bg = dets > bg[:, None]
    native_above_bg = above_bg.copy()
    detection_q = None
    if exact_detection_fdr is not None:
        if not 0 < exact_detection_fdr < 1 or cluster_sizes is None:
            raise ValueError("Exact admission requires declared FDR and exact cluster sizes")
        detection_q = exact_detection_q(dets, cluster_sizes)
        above_bg &= (bg[:, None] > 0) | (detection_q <= exact_detection_fdr)
    scores = np.column_stack([
        z[:, score_cols[c]].mean(axis=1) if score_cols[c].size else np.full(n_groups, -np.inf)
        for c in core_names
    ])
    naming_so = np.column_stack([
        so[:, core_naming_cols[c]].sum(axis=1) if core_naming_cols[c].size
        else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    naming_ok = np.column_stack([
        (so[:, core_naming_cols[c]] & above_bg[:, core_naming_cols[c]]).sum(axis=1)
        if core_naming_cols[c].size else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    naming_before_exact = np.column_stack([
        (so[:, core_naming_cols[c]] & native_above_bg[:, core_naming_cols[c]]).sum(axis=1)
        if core_naming_cols[c].size else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    win = scores.argmax(axis=1)
    srt = np.sort(scores, axis=1)
    margin = srt[:, -1] - srt[:, -2]
    nk_j = core_names.index("Immune_NK_ILC")
    calls, reasons, redirected = [], [], []
    for k in range(n_groups):
        w = core_names[win[k]]
        moved = False
        if w == "Immune_T":
            t_out = bool(so[k, t_cols].any()) if t_cols.size else False
            nk_out = bool(so[k, nk_hint_cols].any()) if nk_hint_cols.size else False
            if (not t_out) and nk_out and naming_before_exact[k, nk_j] >= MIN_NAMING_STANDOUT:
                w, moved = "Immune_NK_ILC", True
        j = core_names.index(w)
        if scores[k, j] <= 0:
            calls.append("Uncalled_no_enrichment")
            reasons.append(f"winning_core_{w}_score_{scores[k, j]:.3f}_not_above_partition_mean")
        elif naming_ok[k, j] >= MIN_NAMING_STANDOUT:
            calls.append(w)
            reasons.append(REDIRECT_REASON if moved else "")
        elif naming_so[k, j] >= MIN_NAMING_STANDOUT:
            calls.append("Uncalled_no_naming_evidence")
            suffix = ("at_or_below_cluster_own_background" if detection_q is None
                      else "at_or_below_cluster_own_background_or_failed_zero_background_exact_detection")
            reasons.append(f"winning_core_{w}_naming_standouts_{suffix}")
        else:
            calls.append("Uncalled_no_naming_evidence")
            reasons.append(f"winning_core_{w}_has_zero_standout_naming_genes")
        redirected.append(moved)
    return {"standout": so, "scores": scores, "naming_so": naming_so, "naming_ok": naming_ok,
            "win": win, "margin": margin, "calls": calls, "reasons": reasons,
            "redirected": redirected, "exact_detection_q": detection_q}


def redirect_counts(dec: dict, sizes) -> dict:
    """The negative-evidence redirect, counted over ONE partition -- the one it is given.

    Two different events, and they are not interchangeable:

    * APPLIED -- the argmax was Immune_T, CD3E/CD3G did not stand out, NK genes did, so the
      winner moved.  ``dec["redirected"]``.
    * SHIPPED -- that moved winner then cleared the enrichment gate and the admission rule, so
      the cluster went out labelled Immune_NK_ILC with :data:`REDIRECT_REASON` as its reason.
      A moved winner whose NK score is not above the partition mean is uncalled instead, and
      counting only ``redirected`` would report it as a label that never existed.

    ``sizes`` is that partition's own cell count per group, so the cell number belongs to the
    same partition as the cluster numbers beside it.  Pass the parent partition and you get the
    parent's counts; pass the final subcluster partition and you get the counts of the labels
    that ship.  The caller must say which in the field name.
    """
    shipped = [j for j, why in enumerate(dec["reasons"]) if why == REDIRECT_REASON]
    sizes = np.asarray(sizes)
    return {
        "clusters_applied": int(sum(bool(x) for x in dec["redirected"])),
        "clusters_shipped": len(shipped),
        "cells_shipped": int(sizes[shipped].sum()) if shipped else 0,
    }


def process_section(donor: str, section: str, h5ad: Path, panels: dict, work: Path,
                    score_on: str = SCORE_ON_DEFAULT, naming_policy: str = "native",
                    uncertainty_decisions: Path | None = None, source_serial: str | None = None) -> dict:
    import anndata as ad

    t0 = time.time()
    if naming_policy != "native" and score_on != "naming":
        raise ValueError("Repaired secretory roles require the established naming-only score")
    if naming_policy != "native" and not source_serial:
        raise ValueError("Repaired generation requires an explicit acquisition source serial")
    matrix_sha256 = None
    if naming_policy != "native":
        with h5ad.open("rb") as input_file:
            matrix_sha256 = hashlib.file_digest(input_file, "sha256").hexdigest()
    panels = naming_policy_panels(panels, naming_policy)
    adata = ad.read_h5ad(h5ad)
    cell_ids = adata.obs_names.to_numpy().astype(str)
    var_names = adata.var_names.to_numpy().astype(str)
    gene_idx = {g: i for i, g in enumerate(var_names)}
    n_cells = int(adata.n_obs)

    obs_counts = adata.obs["nCount_Xenium"].to_numpy(dtype=float)
    obs_features = adata.obs["nFeature_Xenium"].to_numpy(dtype=float)
    sex = str(adata.obs["sex"].iloc[0])
    status = str(adata.obs["disease.status"].iloc[0])
    batch = str(adata.obs["batch"].iloc[0])

    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    del adata
    gc.collect()

    depth = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    depth_matches_obs = int(np.sum(np.abs(depth - obs_counts) < 0.5))

    nucleus_count, nucleus_source = nucleus_for(donor, cell_ids)
    if naming_policy != "native":
        serial_token = "__" + str(source_serial).zfill(7) + "__"
        if not source_serial or serial_token not in nucleus_source:
            raise ValueError("Reviewed source serial is not verified by the exact-key native vendor bundle")
    if nucleus_count is None:
        has_nucleus = np.ones(n_cells, dtype=bool)
        nucleus_measured = False
        n_no_nucleus = 0
    else:
        has_nucleus = nucleus_count > 0
        nucleus_measured = True
        n_no_nucleus = int((~has_nucleus).sum())

    qc_ok = (obs_counts >= QC_MIN_COUNTS) & (obs_features >= QC_MIN_FEATURES)
    keep = qc_ok & has_nucleus & (depth > 0)
    kept = np.flatnonzero(keep)
    n_keep = kept.size

    Xk = X[kept]
    del X
    gc.collect()
    dk = depth[kept]
    rows = np.repeat(np.arange(n_keep), np.diff(Xk.indptr))
    data = np.log1p(Xk.data.astype(np.float64) * TARGET_SUM / dk[rows]).astype(np.float32)
    del rows
    lognorm = sp.csr_matrix((data, Xk.indices, Xk.indptr), shape=Xk.shape)
    n_raw_nonzero = int(Xk.data.size)
    del Xk, data
    gc.collect()

    seed = stable_seed("panc-xenium-cluster-name", section)
    labels, cmeta = cluster_section(lognorm, seed)
    parent_of_kept = np.array([f"{section}_c{int(v):02d}" for v in labels])

    present = [g for g in panels["wanted"] if g in gene_idx]
    pidx = {g: i for i, g in enumerate(present)}
    cols = np.array([gene_idx[g] for g in present])
    L = lognorm[:, cols].astype(np.float32)
    D = L.copy()
    D.data = (D.data >= DETECTION_MIN_LOGNORM).astype(np.float32)
    D.eliminate_zeros()

    core_names = list(LINEAGE_CORES)
    core_cols = {c: np.array([pidx[g] for g in panels["cores"][c][0] if g in pidx], dtype=int)
                 for c in core_names}
    core_naming_cols = {c: np.array([pidx[g] for g in panels["cores"][c][1] if g in pidx],
                                    dtype=int) for c in core_names}
    extra_cols = {c: np.array([pidx[g] for g in panels["extra"][c][0] if g in pidx], dtype=int)
                  for c in panels["extra"]}
    t_cols = np.array([pidx[g] for g in T_DEFINING if g in pidx], dtype=int)
    nk_hint_cols = np.array([pidx[g] for g in NK_NAMING_HINT if g in pidx], dtype=int)
    if score_on == "naming":
        score_cols = core_naming_cols
    elif score_on == "all":
        score_cols = core_cols
    else:
        raise ValueError(f"score_on must be 'naming' or 'all', not {score_on!r}")
    dec_args = (core_names, score_cols, core_naming_cols, t_cols, nk_hint_cols)

    order = sorted(set(parent_of_kept))
    lookup = {c: i for i, c in enumerate(order)}
    krow = np.array([lookup[c] for c in parent_of_kept])
    means, dets, d_other, z, sizes = partition_stats(L, D, krow, len(order))
    exact_parent = ({"exact_detection_fdr": .05, "cluster_sizes": sizes}
                    if naming_policy == "shared_er_exact_v1" else {})
    main = decide(z, dets, d_other, Z_STANDOUT, *dec_args, **exact_parent)
    per_thr = {thr: decide(z, dets, d_other, thr, *dec_args, **exact_parent) for thr in Z_SENSITIVITY}

    # ---- REFINEMENT: a second partition inside each large parent, NAMED AGAINST THE SECTION --
    # Why two levels at all: at one level a lymphocyte cluster is a single group, the B core
    # wins its argmax and the T cells inside it inherit a B label. Splitting each large parent
    # again separates them.
    #
    # Why the naming contrast stays SECTION-WIDE: naming a subcluster against its siblings alone
    # destroys the shared identity of a homogeneous parent -- inside an all-exocrine parent no
    # exocrine gene can stand out against other exocrine subclusters, and the argmax then picks
    # noise. Measured on section 6539: sibling-internal naming cut Exocrine from 175,032 cells to
    # 52,614 and pushed 38,499 cells to Uncalled_no_naming_evidence. So the refinement supplies
    # RESOLUTION only; every z, detection contrast and threshold below is computed across the
    # whole section's final partition, exactly as at the parent level.
    sub_label = np.array(parent_of_kept, dtype=object)
    refined = np.zeros(n_keep, dtype=bool)
    refine_meta = []
    for k, parent in enumerate(order):
        rows_k = np.flatnonzero(krow == k)
        if rows_k.size < REFINE_MIN_CELLS:
            refine_meta.append({"parent": parent, "refined": False,
                                "why": f"{rows_k.size} cells is below REFINE_MIN_CELLS"})
            continue
        sub_seed = stable_seed("panc-xenium-refine", section, parent)
        try:
            slab, smeta = cluster_section(lognorm[rows_k], sub_seed, resolution=REFINE_RESOLUTION)
        except Exception as exc:  # a parent too homogeneous to partition is left as it is
            refine_meta.append({"parent": parent, "refined": False, "why": repr(exc)})
            continue
        if len(set(slab)) < 2:
            refine_meta.append({"parent": parent, "refined": False, "why": "one subcluster only"})
            continue
        sub_label[rows_k] = np.array([f"{parent}_s{int(v):02d}" for v in slab], dtype=object)
        refined[rows_k] = True
        refine_meta.append({"parent": parent, "refined": True,
                            "n_subclusters": int(len(set(slab))),
                            **{f"meta_{a}": b for a, b in smeta.items()}})
        gc.collect()
    del lognorm
    gc.collect()

    forder = sorted(set(sub_label))
    flookup = {c: i for i, c in enumerate(forder)}
    frow = np.array([flookup[c] for c in sub_label])
    fmeans, fdets, fd_other, fz, fsizes = partition_stats(L, D, frow, len(forder))
    exact_final = ({"exact_detection_fdr": .05, "cluster_sizes": fsizes}
                   if naming_policy == "shared_er_exact_v1" else {})
    fdec = decide(fz, fdets, fd_other, Z_STANDOUT, *dec_args, **exact_final)
    fper_thr = {thr: decide(fz, fdets, fd_other, thr, *dec_args, **exact_final) for thr in Z_SENSITIVITY}
    if uncertainty_decisions is not None:
        if naming_policy != "shared_er_exact_v1" or not source_serial:
            raise ValueError("Keyed uncertainty decisions require repaired policy and explicit source serial")
        adjudicated = pd.read_csv(uncertainty_decisions, dtype={"section_id": str, "source_serial": str})
        required = {"section_id", "source_serial", "cluster_id", "n_cells", "candidate_uncertainty_state", "uncertainty_reason"}
        if not required.issubset(adjudicated.columns):
            raise ValueError("Uncertainty decisions must carry exact section/serial/partition/count keys")
        here = adjudicated[adjudicated.section_id.eq(section)]
        if not here.source_serial.eq(str(source_serial).zfill(7)).all() or here.cluster_id.duplicated().any():
            raise ValueError("Uncertainty decision serial or partition identity differs")
        for row in here.itertuples():
            if row.cluster_id not in flookup or int(fsizes[flookup[row.cluster_id]]) != int(row.n_cells):
                raise ValueError("Uncertainty partition absent or membership size differs")
            k = flookup[row.cluster_id]
            fdec["calls"][k] = row.candidate_uncertainty_state
            fdec["reasons"][k] = row.uncertainty_reason
            for sensitivity in fper_thr.values():
                sensitivity["calls"][k] = row.candidate_uncertainty_state
                sensitivity["reasons"][k] = row.uncertainty_reason
    final_call = np.array(fdec["calls"], dtype=object)[frow]
    final_reason = np.array(fdec["reasons"], dtype=object)[frow]
    parent_call = np.array([main["calls"][k] for k in krow], dtype=object)

    # ---- per-cell descriptive scores, for heterogeneity -------------------------------------
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

    # How much naming evidence the cluster's own call rests on, and on what substrate.  The
    # method admits a label on ONE standout naming gene, so a call resting on exactly one gene in
    # a cluster of 24-transcript cells is admissible AND thin.  Both facts travel per cell rather
    # than being argued about: no threshold is invented here and nothing is withheld.
    fcall_j = np.array([core_names.index(x) if x in core_names else -1 for x in fdec["calls"]])
    naming_n_of_cluster = np.array([
        int(fdec["naming_ok"][j, fcall_j[j]]) if fcall_j[j] >= 0 else
        (np.nan if naming_policy == "shared_er_exact_v1" else -1)
        for j in range(len(forder))], dtype=np.float64)[frow]
    cluster_median_depth = np.array([
        float(np.median(depth[kept][frow == j])) for j in range(len(forder))],
        dtype=np.float32)[frow]

    call_j = np.array([core_names.index(x) if x in core_names else -1 for x in final_call])
    naming_of_call = np.where(
        call_j >= 0, cell_naming_det[np.arange(n_keep), np.clip(call_j, 0, None)], np.nan)
    agrees = np.where(call_j >= 0, cell_argmax == np.clip(call_j, 0, None), False) \
        & (cell_argmax_label != "")

    # ---- internal coherence: is the naming evidence of a lineage present in this section? ----
    # Section-wide detection fraction of each coherence gene over the cells that were clustered
    # and called.  Nothing here is compared to an expected value; it exists so that a lineage
    # that receives zero cells while its own naming genes are detected is visible as a SCORING
    # failure to investigate rather than reported as a measurement.
    coherence = {}
    for g in COHERENCE_GENES:
        if g not in pidx:
            coherence[g] = {"on_panel": False}
            continue
        col = np.asarray(D[:, [pidx[g]]].todense()).ravel()
        kmax = int(np.argmax(fdets[:, pidx[g]]))
        coherence[g] = {
            "on_panel": True,
            "det_frac_section": float(col.mean()),
            "n_cells_detected": int(col.sum()),
            "max_cluster_det_frac": float(fdets[kmax, pidx[g]]),
            "max_cluster_id": forder[kmax],
            "max_cluster_lineage": fdec["calls"][kmax],
        }

    y_cols = np.array([pidx[g] for g in Y_GENES if g in pidx], dtype=int)
    floor_cols = np.array([pidx[g] for g in Y_FLOOR if g in pidx], dtype=int)
    y_any = np.asarray(D[:, y_cols].sum(axis=1)).ravel() if y_cols.size else np.zeros(n_keep)
    y_floor = np.asarray(D[:, floor_cols].sum(axis=1)).ravel() if floor_cols.size \
        else np.zeros(n_keep)
    y_per_gene = {g: float(np.asarray(D[:, [pidx[g]]].todense()).ravel().mean())
                  for g in Y_GENES if g in pidx}

    # ---- per-cell frame ----------------------------------------------------------------------
    lineage = np.empty(n_cells, dtype=object)
    reason = np.empty(n_cells, dtype=object)
    lineage[:] = "Uncalled_below_cluster_qc"
    reason[:] = f"nCount_Xenium < {QC_MIN_COUNTS} or nFeature_Xenium < {QC_MIN_FEATURES}"
    no_nuc = ~has_nucleus
    lineage[no_nuc] = "Unavailable_no_nucleus"
    reason[no_nuc] = "nucleus_count == 0; excluded from every calculation"
    lineage[kept] = final_call
    reason[kept] = final_reason
    common = np.array([LINEAGE_CORES[x][1] if x in LINEAGE_CORES else "" for x in lineage],
                      dtype=object)

    def spread_o(values):
        out = np.empty(n_cells, dtype=object)
        out[:] = None
        out[kept] = values
        return out

    def spread_f(values):
        out = np.full(n_cells, np.nan, dtype=np.float32)
        out[kept] = np.asarray(values, dtype=np.float32)
        return out

    def spread_b(values):
        out = np.zeros(n_cells, dtype=bool)
        out[kept] = np.asarray(values, dtype=bool)
        return out

    argmax_full = np.empty(n_cells, dtype=object)
    argmax_full[:] = ""
    argmax_full[kept] = cell_argmax_label

    calls = pd.DataFrame({
        "cell_id": cell_ids,
        "donor_id": donor,
        "section_key": section,
        "roi": ROI,
        "cluster_id": spread_o(parent_of_kept),
        "subcluster_id": spread_o(sub_label),
        "lineage": lineage.astype(str),
        "common_lineage": common.astype(str),
        "reason": reason.astype(str),
        "parent_cluster_lineage": spread_o(parent_call),
        "refined_by_subclustering": spread_b(refined),
        "cluster_inherited": keep,
        "cell_core_argmax": argmax_full.astype(str),
        "cell_agrees_with_cluster_call": spread_b(agrees),
        "n_naming_detected_for_cluster_call": spread_f(naming_of_call),
        "n_naming_standouts_behind_cluster_call": spread_f(naming_n_of_cluster),
        "cluster_median_depth": spread_f(cluster_median_depth),
        "y_probe_detected": spread_b(y_any > 0),
        # A missing measurement is UNAVAILABLE -- not zero and not a negative sentinel.  Where no
        # local bundle is a superset of this h5ad, nucleus_count is pd.NA for every cell and
        # nucleus_measured is False.  Nullable Int32 so the column means the same thing in every
        # partition of the release.
        "nucleus_count": pd.array(
            nucleus_count if nucleus_count is not None else [pd.NA] * n_cells, dtype="Int32"),
        "nucleus_measured": nucleus_measured,
        "depth": depth.astype(np.int32),
        "n_genes_detected": obs_features.astype(np.int32),
    })
    if naming_policy != "native":
        calls["source_serial"] = str(source_serial).zfill(7)
        calls["cell_key"] = section + "|" + calls.cell_id.astype(str)
        calls["serial_qualified_cell_key"] = calls.source_serial + "|" + calls.cell_key
        calls["naming_policy"] = naming_policy
        calls["cell_core_argmax_policy"] = "descriptive_full_core_support_including_shared_ER_not_cell_identity"
        calls["identity_uncertainty"] = calls.lineage.str.startswith("Unresolved_")
        calls["lineage_analysis_eligible"] = calls.cluster_inherited & calls.lineage.isin(LINEAGE_CORES) & ~calls.identity_uncertainty

    # ---- parent-cluster frame ------------------------------------------------------------------
    het = []
    for k, cid in enumerate(order):
        sel = np.flatnonzero(krow == k)
        called = main["calls"][k]
        lab = cell_argmax_label[sel]
        dist = pd.Series(lab).value_counts(normalize=True)
        final_here = pd.Series(final_call[sel]).value_counts()
        het.append({
            "donor_id": donor, "section_key": section, "cluster_id": cid,
            "n_cells": int(sizes[k]),
            "frac_of_section": float(sizes[k] / n_keep),
            "median_depth": float(np.median(depth[kept][sel])),
            "median_n_features": float(np.median(obs_features[kept][sel])),
            "lineage_before_refinement": called,
            "reason_before_refinement": main["reasons"][k],
            "was_refined": bool(refined[sel].all()),
            "n_subclusters": int(pd.Series(sub_label[sel]).nunique()),
            "final_label_majority": final_here.index[0] if len(final_here) else "",
            "final_label_majority_frac": float(final_here.iloc[0] / sel.size)
            if len(final_here) else np.nan,
            "n_final_labels_in_this_cluster": int(len(final_here)),
            "winning_core": core_names[int(main["win"][k])],
            "winning_core_score": float(main["scores"][k, int(main["win"][k])]),
            "score_margin_over_runner_up": float(main["margin"][k]),
            "cd3_negative_nk_redirect": bool(main["redirected"][k]),
            "n_naming_standouts_of_winning_core": int(main["naming_so"][k, int(main["win"][k])]),
            "n_naming_above_background_of_winning_core": int(
                main["naming_ok"][k, int(main["win"][k])]),
            **{f"score__{c}": float(main["scores"][k, j]) for j, c in enumerate(core_names)},
            **{f"full_core_score__{c}": (float(z[k, core_cols[c]].mean())
                                         if core_cols[c].size else np.nan) for c in core_names},
            **{f"naming_standouts__{c}": int(main["naming_so"][k, j])
               for j, c in enumerate(core_names)},
            **{f"extra_score__{c}": (float(z[k, extra_cols[c]].mean())
                                     if extra_cols[c].size else np.nan) for c in extra_cols},
            "purity_cell_argmax_agrees": (float(np.mean(lab == called))
                                          if called in core_names else np.nan),
            "frac_cells_with_naming_gene_of_call": (
                float(np.mean(cell_naming_det[sel, core_names.index(called)] > 0))
                if called in core_names else np.nan),
            "cell_argmax_top_label": dist.index[0] if len(dist) else "",
            "cell_argmax_top_frac": float(dist.iloc[0]) if len(dist) else np.nan,
            "cell_argmax_entropy": float(-np.sum(dist.values * np.log(dist.values + 1e-12)))
            if len(dist) else np.nan,
            "standout_naming_genes_of_call": ";".join(
                g for g in (panels["cores"][called][1] if called in core_names else [])
                if g in pidx and main["standout"][k, pidx[g]]),
        })

    gene_stats = pd.DataFrame({
        "donor_id": donor, "section_key": section,
        "cluster_id": np.repeat(forder, len(present)),
        "gene": np.tile(present, len(forder)),
        "mean_lognorm": fmeans.ravel().astype(np.float32),
        "z_across_clusters": fz.ravel().astype(np.float32),
        "det_frac": fdets.ravel().astype(np.float32),
        "det_frac_rest_of_section": fd_other.ravel().astype(np.float32),
        "standout": fdec["standout"].ravel(),
    })
    if naming_policy != "native":
        gene_stats["exploratory_detection_q"] = fdec["exact_detection_q"].ravel()
        gene_stats["B_shared_ER_support_only"] = gene_stats.gene.isin(SHARED_SECRETORY_NAMING_GENES)
        gene_stats["B_identity_naming_gene_after_repair"] = gene_stats.gene.isin(panels["cores"]["Immune_B_Plasma"][1])

    sub_rows = []
    for j, sid in enumerate(forder):
        sel = np.flatnonzero(frow == j)
        parent = sid.split("_s")[0]
        lab = cell_argmax_label[sel]
        dist = pd.Series(lab).value_counts(normalize=True)
        called = fdec["calls"][j]
        sub_rows.append({
            "donor_id": donor, "section_key": section,
            "parent_cluster_id": parent, "subcluster_id": sid,
            "n_cells": int(fsizes[j]),
            "frac_of_section": float(fsizes[j] / n_keep),
            "median_depth": float(np.median(depth[kept][sel])),
            "median_n_features": float(np.median(obs_features[kept][sel])),
            "parent_lineage": main["calls"][lookup[parent]],
            "lineage": called,
            "common_lineage": LINEAGE_CORES.get(called, ("", ""))[1],
            "reason": fdec["reasons"][j],
            "moved_from_parent": bool(called != main["calls"][lookup[parent]]),
            "winning_core": core_names[int(fdec["win"][j])],
            "winning_core_score": float(fdec["scores"][j, int(fdec["win"][j])]),
            "score_margin_over_runner_up": float(fdec["margin"][j]),
            "cd3_negative_nk_redirect": bool(fdec["redirected"][j]),
            "n_naming_standouts_of_winning_core": int(fdec["naming_so"][j, int(fdec["win"][j])]),
            "n_naming_above_background_of_winning_core": int(
                fdec["naming_ok"][j, int(fdec["win"][j])]),
            **{f"score__{c}": float(fdec["scores"][j, m]) for m, c in enumerate(core_names)},
            **{f"full_core_score__{c}": (float(fz[j, core_cols[c]].mean())
                                         if core_cols[c].size else np.nan) for c in core_names},
            **{f"naming_standouts__{c}": int(fdec["naming_so"][j, m])
               for m, c in enumerate(core_names)},
            **{f"extra_score__{c}": (float(fz[j, extra_cols[c]].mean())
                                     if extra_cols[c].size else np.nan) for c in extra_cols},
            "purity_cell_argmax_agrees": (float(np.mean(lab == called))
                                          if called in core_names else np.nan),
            "frac_cells_with_naming_gene_of_call": (
                float(np.mean(cell_naming_det[sel, core_names.index(called)] > 0))
                if called in core_names else np.nan),
            "cell_argmax_top_label": dist.index[0] if len(dist) else "",
            "cell_argmax_top_frac": float(dist.iloc[0]) if len(dist) else np.nan,
            "cell_argmax_entropy": float(-np.sum(dist.values * np.log(dist.values + 1e-12)))
            if len(dist) else np.nan,
            "standout_naming_genes_of_call": ";".join(
                g for g in (panels["cores"][called][1] if called in core_names else [])
                if g in pidx and fdec["standout"][j, pidx[g]]),
        })

    sens = {}
    for thr in Z_SENSITIVITY:
        sp_ = pd.Series(fper_thr[thr]["calls"])
        sens[str(thr)] = {
            "subclusters": {k: int(v) for k, v in sp_.value_counts().items()},
            "cells": {k: int(fsizes[[i for i, x in enumerate(fper_thr[thr]["calls"]) if x == k]]
                             .sum()) for k in sp_.unique()},
        }

    fin_redirect = redirect_counts(fdec, fsizes)
    par_redirect = redirect_counts(main, sizes)

    report = {
        "donor_id": donor, "section_key": section, "h5ad": str(h5ad),
        "batch": batch, "sex": sex, "disease_status": status,
        "n_cells_in_h5ad": n_cells,
        "n_cells_clustered_and_called": int(n_keep),
        "n_cells_no_nucleus_excluded": n_no_nucleus,
        "n_cells_below_cluster_qc": int((~qc_ok).sum()),
        "nucleus_source": nucleus_source,
        "nucleus_rule_measurable": nucleus_measured,
        "n_genes_on_panel": int(len(var_names)),
        "n_wanted_genes_present": len(present),
        "wanted_genes_absent": [g for g in panels["wanted"] if g not in gene_idx],
        "endocrine_hormones_absent_from_panel": [g for g in ENDOCRINE_HORMONES
                                                 if g not in gene_idx],
        "exocrine_enzymes_absent_from_panel": [g for g in EXOCRINE_ENZYMES if g not in gene_idx],
        "control_features_present_in_var": sorted(
            g for g in var_names if any(g.startswith(p) for p in CONTROL_PREFIXES)),
        "depth_equals_obs_nCount_for_cells": depth_matches_obs,
        "cells_with_depth_gt_target": int(np.sum(depth[kept] > TARGET_SUM)),
        "median_transcripts_per_cell": float(np.median(depth[kept])),
        "median_genes_per_cell": float(np.median(obs_features[kept])),
        "raw_nonzero_entries_kept_cells_all_genes": n_raw_nonzero,
        "y_probe_detect_frac_any": float(np.mean(y_any > 0)),
        "y_floor_group": Y_FLOOR,
        "y_floor_detect_frac_any": float(np.mean(y_floor > 0)),
        "y_probe_detect_frac_per_gene": y_per_gene,
        "n_parent_clusters": len(order),
        "n_parents_refined": int(sum(1 for m in refine_meta if m.get("refined"))),
        "n_subclusters": len(forder),
        "n_cells_moved_by_refinement": int(np.sum(final_call != parent_call)),
        # THE REDIRECT, BY THE PARTITION EACH NUMBER IS ABOUT.  The labels that ship come from
        # the final subcluster partition; the parent partition is a diagnostic and no released
        # label is drawn from it.  Until 2026-09-18 a single field named "n_clusters_..." carried
        # the PARENT count and panc_xenium_release.py published it as n_cd3_negative_nk_redirects,
        # so 6482's redirected subcluster -- 6482_c18_s07, 215 cells, Immune_NK_ILC -- was
        # reported as 0 for every donor in the cohort.
        "n_subclusters_cd3_negative_nk_redirect_applied": fin_redirect["clusters_applied"],
        "n_subclusters_labelled_by_cd3_negative_nk_redirect": fin_redirect["clusters_shipped"],
        "n_cells_labelled_by_cd3_negative_nk_redirect": fin_redirect["cells_shipped"],
        "n_parent_clusters_cd3_negative_nk_redirect_applied": par_redirect["clusters_applied"],
        # The parent partition carries no released label -- it reaches calls/ only as the
        # descriptive column parent_cluster_lineage -- so this one is admitted, never shipped.
        "n_parent_clusters_cd3_negative_nk_redirect_admitted":
            par_redirect["clusters_shipped"],
        "called_cells_by_naming_evidence_count": {
            str(k): int(v) for k, v in
            pd.Series(naming_n_of_cluster[fcall_j[frow] >= 0]).value_counts().sort_index().items()
        },
        "called_cells_on_one_naming_gene_in_low_depth_cluster": int(np.sum(
            (fcall_j[frow] >= 0) & (naming_n_of_cluster <= 1) & (cluster_median_depth < 50))),
        "score_on": score_on,
        "naming_policy": naming_policy,
        "shared_secretory_genes_support_only": (sorted(SHARED_SECRETORY_NAMING_GENES)
                                                 if naming_policy != "native" else []),
        "median_zero_detection_fdr": .05 if naming_policy != "native" else None,
        "uncertainty_decisions_path": str(uncertainty_decisions) if uncertainty_decisions else None,
        "uncertainty_decisions_sha256": (hashlib.sha256(uncertainty_decisions.read_bytes()).hexdigest()
                                           if uncertainty_decisions else None),
        "source_serial": source_serial,
        "input_matrix_sha256": matrix_sha256,
        "exact_detection_inference_limit": "Exploratory detection contrast after expression-based clustering; not guaranteed lineage FDR or independent validation",
        "cell_core_argmax_is_descriptive_full_core_support": True,
        "core_sizes": {c: {"n_genes_on_panel": int(core_cols[c].size),
                           "n_naming_on_panel": int(core_naming_cols[c].size)}
                       for c in core_names},
        "coherence_genes": coherence,
        "cells_by_lineage": {k: int(v) for k, v in
                             pd.Series(lineage.astype(str)).value_counts().items()},
        "subclusters_by_lineage": {k: int(v) for k, v in
                                   pd.Series(fdec["calls"]).value_counts().items()},
        "sensitivity_z": sens,
        "cluster_meta": cmeta,
        "refine_meta": refine_meta,
        "seconds": round(time.time() - t0, 1),
    }

    work.mkdir(parents=True, exist_ok=True)
    calls.to_parquet(work / f"{section}.calls.parquet", index=False)
    pd.DataFrame(het).to_parquet(work / f"{section}.clusters.parquet", index=False)
    pd.DataFrame(sub_rows).to_parquet(work / f"{section}.subclusters.parquet", index=False)
    gene_stats.to_parquet(work / f"{section}.gene_stats.parquet", index=False)
    (work / f"{section}.report.json").write_text(json.dumps(report, indent=1))
    del L, D, calls, gene_stats
    gc.collect()
    return report


def main() -> None:
    import argparse
    global H5AD_DIR
    ap = argparse.ArgumentParser(description="one pancreas Xenium section: cluster, then name")
    ap.add_argument("--section", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--h5ad-dir", type=Path, default=H5AD_DIR,
                    help="directory of raw-count section h5ads")
    ap.add_argument("--score-on", choices=("naming", "all"), default=SCORE_ON_DEFAULT,
                    help="genes carrying a core's score; 'all' is the superseded whole-core mean")
    ap.add_argument("--naming-policy", choices=NAMING_POLICIES, default="native",
                    help="select shared_er_exact_v1 for support-only ER3 and exact median-zero admission")
    ap.add_argument("--uncertainty-decisions", type=Path,
                    help="reviewed CSV keyed by section, source serial, original final partition and exact size")
    ap.add_argument("--source-serial", help="exact acquisition serial; required for shared_er_exact_v1")
    args = ap.parse_args()
    H5AD_DIR = args.h5ad_dir
    panels = load_panels()
    for donor, section, h5ad in sections():
        if section != args.section:
            continue
        r = process_section(donor, section, h5ad, panels, Path(args.work),
                            score_on=args.score_on, naming_policy=args.naming_policy,
                            uncertainty_decisions=args.uncertainty_decisions, source_serial=args.source_serial)
        log(f"{section}: {r['n_cells_clustered_and_called']:,}/{r['n_cells_in_h5ad']:,} cells, "
            f"{r['cluster_meta']['n_clusters']} clusters, {r['seconds']}s")
        return
    raise SystemExit(f"no such pancreas section: {args.section}")


if __name__ == "__main__":
    main()
