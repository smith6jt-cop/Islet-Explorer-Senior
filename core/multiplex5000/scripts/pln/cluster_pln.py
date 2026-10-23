#!/usr/bin/env python3
"""Unsupervised, within-donor clustering of the pLN, in each modality separately.

WHY CLUSTERING HERE AND NOWHERE ELSE
------------------------------------
The maintainer's ruling, 2026-09-17: *"clustering works well in pLN because of the regional
cell types."*  The reason is structural rather than methodological.  A lymph node has real
spatial compartments -- B-cell follicles and germinal centres, the T zone / paracortex, the
medulla, the subcapsular sinus -- and each carries a distinct composition, so an unsupervised
partition of the measurements has something to recover.  A pancreas does not organise that
way, which is why clustering is not the method there.

WHAT A CLUSTER IS ALLOWED TO BE
-------------------------------
A description.  Not a verdict.  CLAUDE.md: existing QuPath classifications, annotations,
phenotypes, clusters, disease labels and gold-standard reviews are **audit context only** --
never training, tuning or release-gating inputs.  The endorsement of clustering as a method
does not lift that, so:

*   The clustering reads measurements only.  `qupath_classification` is carried nowhere near
    the feature matrix; neither is disease status, batch, node membership or any existing call.
*   No cluster gates, tunes or validates anything.  Nothing downstream is accepted or rejected
    because a cluster agreed with it.
*   No cluster is **named** for an anatomy.  Clusters are `<donor>_c<NN>` and the enrichment
    table says what each one measurably expresses.  The anatomical reading is the reader's,
    stated separately and as an inference, in
    `docs/PLN_REGIONAL_CLUSTERING_2026-09-17.md`.

CLUSTERING IS THE SUBSTITUTE FOR CONTEXT, AND THAT IS MEASURED
--------------------------------------------------------------
In the pancreas the QuPath pixel-classifier objects hand you the context for free -- islet,
blood vessel, lymphatic, nerve, ordinary tissue -- and a reference level is fitted inside it.
A pLN section has no such partition, and that is **measured, not assumed**.  In
`core/phenocycler/data/integration/cells_pln`:

*   there is no `context` column at all (128 columns, none of them it);
*   `parent_raw` is the single constant string `'Annotation (Tissue)'` on all 8,961,168 rows,
    with no nulls;
*   `qupath_classification` is NULL on every one of those rows.

"Tissue" here names the whole section and therefore names nothing.  There is no ordinary
tissue to draw a reference population from, because there is no partition of the section at
all.  Clustering is what replaces it: the node's regional cell types are distinct enough that
an unsupervised partition *discovers* the structure the pancreas is handed.  The naming that
follows is internal to one section -- a cluster is named by how its canonical markers stand
against **the other clusters on that same section** -- and no level, quantile, threshold or
scale factor is imported from anywhere.  `naming_profiles.parquet` is that table; see
`naming_profile_pheno` / `naming_profile_xen` below.

WITHDRAWN, AND NOT TO BE REINTRODUCED
-------------------------------------
*   **Any number derived from pancreas cells reaching a pLN call** -- level, threshold, scale
    factor or quantile.  An earlier attempt transferred each donor's pancreas p99 onto their
    pLN slide; the transferred CD66 level left 1.9% of that slide beneath it.  Withdrawn.
    Those numbers are a property of the mistake and are not findings.
*   **`subpanel.fit_levels` / `reference_level` / `context == 'tissue'` applied to the pLN.**
    There is no context column and no reference population to fit against.  An earlier
    version of this docstring proposed fitting the levels on the pLN slide itself instead of
    on the pancreas; that is withdrawn too, because it synthesises the reference population
    the measurement above says is absent.
*   **`IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor on pLN values.**  It is defined against
    a reference level that does not exist here.

Nothing here computes a reference level or calls a cell positive.  The scaling below is a
*within-donor* robust standardisation used to give every marker comparable weight in a
distance; it is not a threshold and never becomes one.

WITHIN DONOR FIRST
------------------
Every embedding and every partition is built inside one donor.  Donors are related only
afterwards, by correlating cluster marker profiles.  The cohort's acquisition batch is aliased
with disease status at the coarse level, and this project has already been bitten by exactly
that: a deployed single-cell pseudotime turned out to be a batch axis at batch eta-squared
0.988.  So anything this script produces *across* donors is reported with a measured batch
eta-squared beside it, and with the disease eta-squared beside that, because the two cannot be
separated at the coarse level and presenting one alone would hide it.

THE BATCH VARIABLE IS A PROXY, AND THAT IS STATED IN THE OUTPUT
--------------------------------------------------------------
`data/donor_acquisition_runs.csv` and `data/donor_acquisition_batch.csv` record the
**pancreas** scans (`6442_Scan1.er.qptiff`).  For most donors the pLN is a separate image
(`6442pLN_Scan1.er.qptiff`) that no tracked manifest in this checkout records.  A few donors
are the exception -- their recorded image is a combined block (`HDL115pancLN_Scan1.er.qptiff`,
`0080870_6533_pLNpanc_Scan3.er.qptiff`) and for those the covariate *is* the pLN acquisition.
`covariate_provenance.csv` says which is which, measured from the manifests rather than
asserted.  Everywhere else `stain_group` / `acq_run` / `batch` are *donor-level* grouping
variables, not a verified pLN scan batch.  A high eta-squared on them is still a warning; a
low one is not a clearance.

MODALITIES ARE CLUSTERED SEPARATELY, NEVER JOINTLY
--------------------------------------------------
Different measurements of different sections.  There is no shared embedding, no co-plot in
space, and no pLN-to-pancreas spatial comparison anywhere in this script: pancreas and pLN are
different blocks.

REPRODUCIBLE, AND THAT IS CHECKED RATHER THAN CLAIMED
-----------------------------------------------------
Every random draw is seeded from `stable_seed()`, a BLAKE2b digest of the run's fixed seed and
the donor/section string.  The previous version used `SEED + hash(donor)`; Python salts `str`
hashing per process unless PYTHONHASHSEED is set, so the documented commands did not reproduce
their own output.  A second instance of the same defect was fixed with it: `eta_squared`
accumulated its levels in `set` order, so a saturated eta came back as 1.0 in one interpreter
and 1.0000000000000002 in the next, which moved a permutation p.  Levels are summed sorted and
the result is clipped to [0, 1].

`--stage determinism` proves it rather than asserting it.  It runs the cluster stage twice in
two fresh interpreters with **different** PYTHONHASHSEED values, compares the assignment
parquet bytes by sha256, and checks the tracked build reproduces from the same command; the
result goes to `determinism_check.json`.  `--out-root` is what makes that possible and never
changes what is computed.

    .venv/bin/python core/multiplex5000/scripts/pln/cluster_pln.py --modality phenocycler
    .venv/bin/python core/multiplex5000/scripts/pln/cluster_pln.py --modality xenium

    # the stages, independently
    ... --modality phenocycler --stage cluster       # partitions, per donor/section
    ... --modality phenocycler --stage profile       # naming_profiles.parquet + coherence
    ... --modality phenocycler --stage determinism --donors 115 6461 6534
    ... --modality phenocycler --stage summarize     # enrichment, nodes, batch eta-squared
    ... --modality phenocycler --stage link
    ... --modality phenocycler --stage prose
    ... --modality phenocycler --stage published-eta  # batch term for the PUBLISHED quantities

A PARTIAL RUN NEVER WRITES A COHORT-WIDE NAME
---------------------------------------------
`--donors` makes a run PARTIAL, and a partial run's cohort-shaped outputs -- everything that
aggregates across donors: `naming_profiles.parquet`, `cluster_report.json`, `summary.json`,
`batch_eta_squared.csv`, `cluster_immune_link.*`, `prose_facts.json` and the rest -- are
written under `_partial/donors=<slug>/` instead.  The cohort-wide files are not touched, not
read and not merged.  The rule, the two lanes that share it and why it is not a merge are in
`partial_run_guard.py`; it exists because `--donors X --stage standout` truncated three
cohort-wide files on 2026-09-17 and nothing said so.  The per-donor outputs -- `assignments/`,
`profiles/<donor>.parquet`, `naming_profiles/<donor>.parquet` -- are keyed by their unit and
are written where they always were.

THE BATCH TERM FOR THE PUBLISHED QUANTITIES
-------------------------------------------
`--stage published-eta` reads the per-donor table each pLN immune lane publishes and reports
eta-squared, its permutation null and p for EVERY cross-donor quantity in it, against
`stain_group`, `acq_run`, `coarse_batch`, `software`, `disease_status` and -- for Xenium --
the node section's own instrument `xenium_batch`.  Nothing is dropped, corrected or gated on
the answer; the batch term is published beside the quantity so it is not absent from the text.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[4]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from partial_run_guard import RunScope, write_scope_note  # noqa: E402

PHENO_CELLS = REPO / "core/phenocycler/data/integration/cells_pln"
XEN_H5AD_DIR = Path("/data/xenium/h5ad_out")
XEN_CELLS = REPO / "core/phenocycler/data/integration/cells_xen"
#: The node-membership layer build_pln_layer.py publishes, under the SAME root it writes to: this path
#: was pinned to /data/islet_layers until 2026-10-01, so a build into another root ($ISLET_LAYERS_DIR)
#: was clustered against the old root's membership -- other cells, joined by object id to nothing.
PLN_LAYER = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")) / "pln"
ACQ_RUNS = REPO / "data/donor_acquisition_runs.csv"
ACQ_BATCH = REPO / "data/donor_acquisition_batch.csv"
OUT_ROOT = REPO / "data/development/typing/pln/clusters"

SEED = 20260917

#: What this invocation covers.  Set once in `main()` from `--donors`, then refined by whichever
#: stage enumerates the units it could have processed.  Every cohort-wide write in this file goes
#: through `cohort_out()`, which asks this object where the file belongs -- so a partial run
#: cannot land on a cohort-wide name.  See `partial_run_guard.py` for the rule and why it is this
#: one rather than a merge.
SCOPE = RunScope.build("donor", None, ())


def set_scope(scope: RunScope) -> RunScope:
    global SCOPE
    SCOPE = scope
    return SCOPE


def cohort_out(out_dir: Path, name: str) -> Path:
    """Where a cohort-wide artifact belongs for THIS run.

    A whole run gets the cohort-wide path.  A partial run gets
    `<out_dir>/_partial/donors=<slug>/<name>`, and the cohort-wide file of that name is left
    untouched -- not read, not merged, not removed.
    """
    return SCOPE.cohort_path(out_dir, name)


def stable_seed(*parts: object) -> int:
    """A seed derived from strings that is the same in every interpreter.

    `hash('6461')` is salted per process unless PYTHONHASHSEED is set, so seeding with
    `SEED + hash(donor)` made every run of this script a different run wearing the same
    command.  BLAKE2b over the joined parts has no such salt.  The value is reduced into the
    32-bit range because scanpy, scikit-learn and igraph all take `random_state` as an int.
    """
    digest = hashlib.blake2b("|".join(str(p) for p in parts).encode("utf-8"),
                             digest_size=8).digest()
    return int.from_bytes(digest, "big") % (2 ** 31 - 1)

#: DAPI is the nuclear counterstain and Beta_actin the loading control.  Neither carries
#: regional identity and both track section thickness and staining, so they would dominate a
#: distance for a reason that is not biology.  Every other measured marker is kept -- including
#: the structural ones, because "this region is stroma/vessel/capsule" is a regional statement
#: in a node and dropping them would decide the answer in advance.
PHENO_FEATURE_EXCLUDE = frozenset({"DAPI", "Beta_actin"})

#: The regional sets named in the maintainer's ruling.  These are REPORTED, not used to build
#: or steer the partition -- the embedding sees every kept marker equally.
PHENO_REGIONAL: dict[str, list[str]] = {
    "fdc": ["CD21"],
    "b_follicle": ["CD20", "CD79a"],
    "germinal_centre_proliferation": ["Ki67", "PCNA", "Bcl_2"],
    "t_zone": ["CD3e", "CD4", "CD8", "TCF_1", "TOX", "PD_1", "ICOS"],
    "medulla_sinus_macrophage": ["CD68", "CD163", "CD206", "CD209"],
    "frc_lymphatic_sinus_lining": ["Podoplanin"],
}
#: Non-immune structural markers on this panel.  Reported so a cluster that is stroma or vessel
#: is visible as such rather than being read as an immune population.
#:
#: NO IMMUNE CANONICAL MAY APPEAR HERE (D1, measured 2026-09-17).  CD66 sat in this list while
#: `lineages.LINEAGE_MARKERS["Granulocyte"] = ("MPO", "CD66")` names it a granulocyte canonical,
#: and that module's own docstring records an earlier lane dropping CD66 and thereby "folding
#: granulocytes into the myeloid label".  A granulocyte cluster was therefore scored as structural
#: BY ITS OWN CANONICAL: removing CD66 moves `set__non_immune_structural` by up to 0.129 and flips
#: its sign on 3 of the 293 clusters (89,442 cells).  `lineages.py` is the authority on what a
#: marker means and this list is checked against it by
#: `tests/current/test_pln_immune_pheno.py::test_no_immune_canonical_is_listed_as_structural`.
PHENO_NON_IMMUNE = [
    "CD31", "SMA", "Vimentin", "Collagen_IV", "Caveolin", "Podoplanin", "CD34",
    "E_cadherin", "EpCAM", "Pan_Cytokeratin", "TP63", "B3TUBB", "CD141",
    "SOX2", "GP100",
]

XEN_REGIONAL: dict[str, list[str]] = {
    "fdc_follicle": ["CR2"],
    "follicular": ["CXCL13"],
    "t_zone_frc": ["CCL19", "CCL21"],
    "germinal_centre": ["BCL6", "AICDA"],
    "proliferation": ["MKI67"],
    "naive_memory": ["SELL", "CD27"],
}
#: Canonical lineage genes.  CD3D and LYZ were named in the brief but are measured absent from
#: all 24 pLN h5ads (see the report's `named_but_absent`), so they are not silently substituted.
XEN_CANONICAL: dict[str, list[str]] = {
    "pan_leukocyte": ["PTPRC"],
    "t": ["CD3E", "CD3G", "CD4", "CD8A", "FOXP3"],
    "b": ["MS4A1", "CD79A", "CD19"],
    "myeloid": ["CD68", "ITGAX", "ITGAM"],
    "nk": ["NKG7", "KLRD1"],
}
XEN_NAMED_IN_BRIEF = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "MS4A1", "CD79A", "CD19", "CD68",
    "LYZ", "ITGAX", "ITGAM", "NKG7", "KLRD1", "FOXP3", "CR2", "CXCL13", "CCL19", "CCL21",
    "BCL6", "AICDA", "SELL", "CD27", "MKI67",
]

IDENTITY_CORE = REPO / "core/panels/data/identity_core"

#: Measured 2026-09-17 by `derive_pln_node_membership`/`build_pln_layer` and re-measured here:
#: 6534's one surviving node outline is 10,742,163 um2 against a 9,863,458 um2 tissue outline
#: (ratio 1.089) and contains all 147,959 of its cells.  A per-node number for 6534 is a
#: per-slide number wearing a node's name.
WHOLE_SLIDE_NODE_DONORS = frozenset({"6534"})


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def pheno_partitioning() -> ds.Partitioning:
    """The explicit schema the pLN root must be read with.

    Bare ``partitioning='hive'`` infers ``donor_id`` as int32 here, because every pLN donor id
    happens to be all digits.  That silently turns '0112'-style ids into integers which will
    not join against the string ``donor_id`` used everywhere else.
    """
    return ds.partitioning(
        pa.schema([("donor_id", pa.string()), ("roi", pa.string())]), flavor="hive"
    )


def subpanel_genes() -> dict[str, list[str]]:
    """The `core/panels` identity-core subpanels, read from the one authority.

    Read, never copied.  Nothing here writes a second copy of these CSVs.
    """
    out: dict[str, list[str]] = {}
    for path in sorted(IDENTITY_CORE.glob("*.csv")):
        frame = pd.read_csv(path)
        out[path.stem] = [str(g) for g in frame["gene"].tolist()]
    return out


# ----------------------------------------------------------------------------------
# shared numerics
# ----------------------------------------------------------------------------------

def robust_z(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Column-wise (median, MAD) standardisation, computed inside one donor.

    This is a scaling, not a threshold and not a reference level.  A column whose MAD is zero
    -- a marker that did not vary on that slide -- is returned as all-zero rather than being
    divided by zero, and its zero MAD is returned so the caller can record it.
    """
    med = np.median(values, axis=0)
    mad = np.median(np.abs(values - med), axis=0) * 1.4826
    safe = np.where(mad > 0, mad, 1.0)
    z = (values - med) / safe
    z[:, mad <= 0] = 0.0
    return z.astype(np.float32), med, mad


def leiden_on_subsample(
    features: np.ndarray,
    *,
    n_pcs: int,
    subsample: int,
    resolution: float,
    seed: int,
) -> tuple[np.ndarray, dict]:
    """Leiden on a random subsample, then k-NN assignment of every remaining cell.

    Leiden needs a neighbour graph, and a graph over five million cells is not affordable per
    donor.  The partition is therefore learned on a random subsample and carried to the rest by
    a 15-nearest-neighbour vote in the same PCA space.  The cost of that shortcut is measured
    rather than assumed: the subsample is split 80/20, the classifier is refit on the 80 and
    its agreement with Leiden's own labels on the held-out 20 is reported as
    `holdout_agreement`.  A low number there means the assignment step, not the partition, is
    what a reader should distrust.
    """
    import scanpy as sc
    import anndata as ad
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KNeighborsClassifier

    rng = np.random.default_rng(seed)
    n_cells = features.shape[0]
    n_pcs = int(min(n_pcs, features.shape[1] - 1, n_cells - 1))
    pca = PCA(n_components=n_pcs, random_state=seed)
    coords = pca.fit_transform(features).astype(np.float32)

    take = min(subsample, n_cells)
    idx = rng.choice(n_cells, size=take, replace=False) if take < n_cells else np.arange(n_cells)
    idx = np.sort(idx)

    sub = ad.AnnData(X=coords[idx].copy())
    sub.obsm["X_pca"] = coords[idx].copy()
    sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        sub,
        resolution=resolution,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=seed,
    )
    sub_labels = sub.obs["leiden"].to_numpy().astype(str)

    # measured cost of the assignment shortcut
    perm = rng.permutation(take)
    cut = int(take * 0.8)
    tr, te = perm[:cut], perm[cut:]
    probe = KNeighborsClassifier(n_neighbors=15, n_jobs=-1)
    probe.fit(coords[idx][tr], sub_labels[tr])
    holdout = float((probe.predict(coords[idx][te]) == sub_labels[te]).mean())

    labels = np.empty(n_cells, dtype=object)
    labels[idx] = sub_labels
    rest = np.setdiff1d(np.arange(n_cells), idx, assume_unique=True)
    if rest.size:
        clf = KNeighborsClassifier(n_neighbors=15, n_jobs=-1)
        clf.fit(coords[idx], sub_labels)
        step = 200_000
        for start in range(0, rest.size, step):
            chunk = rest[start:start + step]
            labels[chunk] = clf.predict(coords[chunk])

    meta = {
        "n_cells": int(n_cells),
        "n_pcs": n_pcs,
        "pca_variance_explained": float(pca.explained_variance_ratio_.sum()),
        "subsample_used_for_graph": int(take),
        "leiden_resolution": resolution,
        "n_clusters": int(len(set(sub_labels))),
        "knn_assignment_k": 15,
        "knn_holdout_agreement": holdout,
        "cells_assigned_by_knn": int(rest.size),
        "seed": int(seed),
    }
    return labels.astype(str), meta


#: eta-squared at or above this is treated as saturated: the grouping has as many levels as it
#: needs to partition the values exactly, and the number carries no information about batch.
ETA_SATURATED = 1.0 - 1e-9


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """Fraction of the variance in `values` that lies between levels of `groups`.

    Returned as NaN when there is no variance to partition or fewer than two levels, rather
    than as a zero that would read as "no batch effect".

    The levels are summed in **sorted** order, not set order.  Iterating a Python `set` of
    strings follows the per-process string hash, so the accumulation order of `ss_between`
    changed between interpreters and a saturated eta came back as 1.0 in one run and
    1.0000000000000002 in the next -- which flipped a permutation p from 0.120 to 0.023 on a
    near-saturated quantity.  That is the same defect as the seeding bug, one level down.
    The result is clipped into [0, 1] for the same reason.
    """
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    values, groups = values[ok], np.asarray(groups)[ok]
    levels = sorted({str(g) for g in groups.tolist()})
    if values.size < 3 or len(levels) < 2:
        return float("nan")
    grand = values.mean()
    ss_total = float(((values - grand) ** 2).sum())
    if ss_total <= 0:
        return float("nan")
    as_str = np.asarray([str(g) for g in groups.tolist()])
    ss_between = 0.0
    for level in levels:
        sel = as_str == level
        ss_between += int(sel.sum()) * (values[sel].mean() - grand) ** 2
    return float(min(max(ss_between / ss_total, 0.0), 1.0))


def eta_squared_with_null(
    values: np.ndarray, groups: np.ndarray, *, n_perm: int = 2000, seed: int = SEED
) -> dict:
    """eta-squared beside the eta-squared a *random* grouping of the same shape would give.

    WHY THE NULL IS NOT OPTIONAL HERE.  eta-squared rises mechanically with the number of
    levels.  This cohort has 14 clusterable pLN donors spread over 8 stain groups and 11
    acquisition dates, so a grouping with one or two donors per level partitions essentially
    all of the variance no matter what the values are -- `stain_group` returns 1.000 on data
    with no batch structure whatsoever.  Reporting that number alone would manufacture a batch
    effect out of degrees of freedom, and reporting it as "no batch effect" when it is low
    would be the same error in the other direction.  So every eta-squared here carries the
    median eta-squared of 2000 random relabellings with identical level sizes, and the
    fraction of those relabellings that reached the observed value.
    """
    values = np.asarray(values, dtype=float)
    observed = eta_squared(values, groups)
    if not np.isfinite(observed):
        return {"eta_squared": float("nan"), "null_median": float("nan"),
                "p_permutation": float("nan"), "n_perm": 0}
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups)
    draws = np.empty(n_perm)
    for i in range(n_perm):
        draws[i] = eta_squared(values, rng.permutation(groups))
    good = np.isfinite(draws)
    # A relabelling that ties the observed value HAS reached it.  Comparing with `>=` on raw
    # floats made a tie turn on the last bit; the tolerance makes the tie count, which is both
    # deterministic and the conservative direction -- it will not report a batch effect that
    # exists only in rounding.
    tol = 1e-12
    return {
        "eta_squared": observed,
        "null_median": float(np.median(draws[good])) if good.any() else float("nan"),
        "p_permutation": float((draws[good] >= observed - tol).mean())
            if good.any() else float("nan"),
        "n_perm": int(good.sum()),
        "eta_squared_saturated": bool(observed >= ETA_SATURATED),
        "null_fraction_also_saturated": float((draws[good] >= ETA_SATURATED).mean())
            if good.any() else float("nan"),
    }


def covariate_provenance(donors: list[str]) -> pd.DataFrame:
    """For each donor, whether the acquisition record is of a slide that carries the pLN.

    Measured from the two manifests' `image_file`, not assumed: a few donors were scanned as a
    combined pancreas+pLN block, and for those the recorded `acq_run` / `stain_group` is the
    pLN acquisition.  For the rest the pLN is a separate image with no tracked record, and the
    covariate is a donor-level proxy.
    """
    frames = []
    for name, path in (("runs", ACQ_RUNS), ("batch", ACQ_BATCH)):
        frame = pd.read_csv(path, dtype=str)[["donor_id", "image_file"]]
        frame["manifest"] = name
        frames.append(frame)
    both = pd.concat(frames, ignore_index=True)
    both["image_covers_the_pln"] = both["image_file"].str.upper().str.contains("LN", na=False)
    rows = []
    for donor in donors:
        sub = both[both["donor_id"] == donor]
        rows.append({
            "donor_id": donor,
            "recorded_images": ";".join(sorted(set(sub["image_file"].dropna()))),
            "an_image_record_covers_the_pln": bool(sub["image_covers_the_pln"].any()),
            "n_records": int(len(sub)),
        })
    return pd.DataFrame(rows)


def donor_covariates() -> pd.DataFrame:
    """Donor-level batch and disease, from the tracked pancreas acquisition manifests.

    These files record the **pancreas** scan.  No tracked manifest in this checkout records the
    pLN scans, so these are donor-level grouping variables and not a verified pLN scan batch.
    """
    runs = pd.read_csv(ACQ_RUNS, dtype=str)
    batch = pd.read_csv(ACQ_BATCH, dtype=str)
    merged = runs[["donor_id", "acq_run", "stain_group", "software"]].merge(
        batch[["donor_id", "disease.status", "batch"]], on="donor_id", how="outer"
    )
    merged = merged.rename(columns={"disease.status": "disease_status", "batch": "coarse_batch"})
    return merged.set_index("donor_id")


# ----------------------------------------------------------------------------------
# PhenoCycler
# ----------------------------------------------------------------------------------

def pheno_marker_columns(dataset: ds.Dataset) -> list[str]:
    """Cell-mean marker columns: every column that also has a `Nucleus__` twin."""
    names = set(dataset.schema.names)
    return [
        n for n in dataset.schema.names
        if not n.startswith("Nucleus__") and f"Nucleus__{n}" in names
    ]


def cluster_pheno_donor(
    donor: str,
    markers: list[str],
    out_dir: Path,
    *,
    subsample: int,
    resolution: float,
    seed_base: int,
) -> dict:
    part = pheno_partitioning()
    dataset = ds.dataset(PHENO_CELLS, format="parquet", partitioning=part)
    cols = ["object_id", "image", "X_centroid", "Y_centroid", "qc_analysis_eligible"] + markers
    table = dataset.to_table(
        columns=cols, filter=(ds.field("donor_id") == donor) & ds.field("qc_analysis_eligible")
    )
    frame = table.to_pandas()
    del table
    gc.collect()
    n_eligible = len(frame)
    if n_eligible == 0:
        return {"donor_id": donor, "status": "skipped_no_eligible_cells", "n_eligible": 0}

    feats = [m for m in markers if m not in PHENO_FEATURE_EXCLUDE]
    values = frame[feats].to_numpy(dtype=np.float64)
    # A marker with no measurement on this slide is UNAVAILABLE, never zero.  It leaves the
    # feature set for this donor and is named in the report.
    all_null = np.isnan(values).all(axis=0)
    unavailable = [f for f, bad in zip(feats, all_null) if bad]
    keep = [f for f, bad in zip(feats, all_null) if not bad]
    values = values[:, ~all_null]

    row_ok = np.isfinite(values).all(axis=1)
    dropped_partial = int((~row_ok).sum())
    frame = frame.loc[row_ok].reset_index(drop=True)
    values = values[row_ok]

    # log1p compresses the intensity tail; the robust z then gives every marker comparable
    # weight in the distance.  Both are within-donor.
    z, med, mad = robust_z(np.log1p(np.clip(values, 0.0, None)))
    zero_mad = [k for k, m in zip(keep, mad) if m <= 0]
    z = np.clip(z, -5.0, 10.0)

    seed = stable_seed("pln-cluster", seed_base, "phenocycler", donor)
    labels, meta = leiden_on_subsample(
        z, n_pcs=20, subsample=subsample, resolution=resolution, seed=seed
    )
    cluster_id = np.array([f"{donor}_c{int(v):02d}" for v in labels])

    assign = pd.DataFrame(
        {
            "object_id": frame["object_id"].astype(str),
            "donor_id": donor,
            "modality": "phenocycler",
            "roi": "pln_1",
            "cluster_id": cluster_id,
            "X_centroid": frame["X_centroid"].to_numpy(),
            "Y_centroid": frame["Y_centroid"].to_numpy(),
        }
    )
    dest = out_dir / "assignments" / f"donor_id={donor}"
    dest.mkdir(parents=True, exist_ok=True)
    assign.to_parquet(dest / "clusters.parquet", index=False)

    # per-cluster profiles: mean within-donor z and mean raw intensity, per marker
    prof_rows = []
    for cid in sorted(set(cluster_id)):
        sel = cluster_id == cid
        row = {
            "donor_id": donor,
            "cluster_id": cid,
            "n_cells": int(sel.sum()),
            "frac_of_donor": float(sel.mean()),
        }
        for j, marker in enumerate(keep):
            row[f"z__{marker}"] = float(z[sel, j].mean())
            row[f"raw__{marker}"] = float(values[sel, j].mean())
        prof_rows.append(row)
    profiles = pd.DataFrame(prof_rows)
    profiles.to_parquet(out_dir / "profiles" / f"{donor}.parquet", index=False)

    result = {
        "donor_id": donor,
        "status": "clustered",
        "image": str(frame["image"].iloc[0]),
        "n_eligible": int(n_eligible),
        "n_clustered": int(len(frame)),
        "cells_dropped_partial_measurement": dropped_partial,
        "markers_used": keep,
        "markers_unavailable_on_this_slide": unavailable,
        "markers_zero_mad_forced_to_zero_z": zero_mad,
        **meta,
    }
    del frame, values, z
    gc.collect()
    return result


def run_pheno(args) -> None:
    out_dir = OUT_ROOT / "phenocycler"
    (out_dir / "profiles").mkdir(parents=True, exist_ok=True)
    part = pheno_partitioning()
    dataset = ds.dataset(PHENO_CELLS, format="parquet", partitioning=part)
    markers = pheno_marker_columns(dataset)

    census = (
        dataset.to_table(columns=["donor_id", "qc_analysis_eligible"])
        .to_pandas()
        .groupby("donor_id", observed=True)["qc_analysis_eligible"]
        .agg(["size", "sum"])
        .rename(columns={"size": "rows", "sum": "eligible"})
    )
    log(f"pLN PhenoCycler root: {int(census['rows'].sum())} rows, "
        f"{int(census['eligible'].sum())} analysis-eligible, {len(census)} donors")

    donors = [d for d in census.index.tolist() if census.loc[d, "eligible"] > 0]
    no_measurements = [d for d in census.index.tolist() if census.loc[d, "eligible"] == 0]
    clusterable = list(donors)
    if args.donors:
        donors = [d for d in donors if d in set(args.donors)]
    scope = set_scope(SCOPE.refine(clusterable))
    scope.announce(out_dir, log)
    write_scope_note(scope, out_dir, {"stage": "cluster", "modality": "phenocycler"})
    log(f"clusterable donors ({len(donors)}): {donors}")
    log(f"donors with no intensity measurements ({len(no_measurements)}): {no_measurements} "
        "-- these need re-running in QuPath, a maintainer action")

    results = []
    for donor in donors:
        started = time.time()
        log(f"  {donor}: clustering")
        res = cluster_pheno_donor(
            donor, markers, out_dir, subsample=args.subsample, resolution=args.resolution,
            seed_base=args.seed,
        )
        res["seconds"] = round(time.time() - started, 1)
        log(f"  {donor}: {res.get('n_clusters', 0)} clusters over {res.get('n_clustered', 0)} "
            f"cells, k-NN holdout agreement {res.get('knn_holdout_agreement', float('nan')):.3f} "
            f"({res['seconds']}s)")
        results.append(res)

    report = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": "phenocycler",
        "source": str(PHENO_CELLS),
        "read_with": "explicit hive schema donor_id:string, roi:string",
        "rows_in_root": int(census["rows"].sum()),
        "rows_analysis_eligible": int(census["eligible"].sum()),
        "donors_in_root": census.index.tolist(),
        "donors_with_no_intensity_measurements": no_measurements,
        "donors_clustered": donors,
        "seed_base": int(args.seed),
        "features": {
            "markers_in_export": markers,
            "excluded_from_features": sorted(PHENO_FEATURE_EXCLUDE),
            "transform": "log1p then within-donor robust z (median / 1.4826*MAD), clipped [-5, 10]",
            "note": "a scaling, not a threshold; no reference level and no positivity call here",
        },
        "labels_never_seen_by_the_clustering": [
            "qupath_classification", "parent_raw", "disease_status", "batch", "node membership",
        ],
        "per_donor": results,
    }
    cohort_out(out_dir, "cluster_report.json").write_text(json.dumps(report, indent=2))
    log(f"wrote {cohort_out(out_dir, 'cluster_report.json')}")


# ----------------------------------------------------------------------------------
# Xenium
# ----------------------------------------------------------------------------------

def xen_files() -> list[tuple[str, str, Path]]:
    """(donor_id, section_key, path) for each pLN h5ad.

    The donor comes from the filename stem, which is the only identifier these files carry --
    and the one the integration lane already uses for `cells_xen`.  Sections are kept apart:
    6579 contributes three (a/b/c), 6457 two.  `cell_id` is preserved verbatim from `obs`.
    """
    out = []
    for path in sorted(XEN_H5AD_DIR.glob("*_pLNXenium.h5ad")):
        stem = path.name.replace("_pLNXenium.h5ad", "")
        core = stem[3:] if stem.startswith("HDL") else stem
        donor = core.rstrip("abc") if core[-1:].isalpha() else core
        out.append((donor, stem, path))
    return out


def cluster_xen_section(
    donor: str,
    section: str,
    path: Path,
    out_dir: Path,
    report_genes: list[str],
    *,
    subsample: int,
    resolution: float,
    min_counts: int,
    seed_base: int,
) -> dict:
    import scanpy as sc
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.read_h5ad(path)
    n_in = int(adata.n_obs)
    cell_ids = adata.obs_names.to_numpy().astype(str)

    counts = adata.obs["nCount_Xenium"].to_numpy(dtype=float)
    genes_det = adata.obs["nFeature_Xenium"].to_numpy(dtype=float)
    keep = (counts >= min_counts) & (genes_det >= 5)
    n_dropped = int((~keep).sum())
    adata = adata[keep].copy()
    cell_ids = cell_ids[keep]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    lognorm = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    var_names = adata.var_names.to_numpy().astype(str)

    work = adata.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=2000)
    work = work[:, work.var["highly_variable"]].copy()
    sc.pp.scale(work, max_value=10)
    feats = np.asarray(work.X, dtype=np.float32)
    del work
    gc.collect()

    seed = stable_seed("pln-cluster", seed_base, "xenium", section)
    labels, meta = leiden_on_subsample(
        feats, n_pcs=30, subsample=subsample, resolution=resolution, seed=seed
    )
    del feats
    gc.collect()
    cluster_id = np.array([f"{section}_c{int(v):02d}" for v in labels])

    assign = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "donor_id": donor,
            "section_key": section,
            "modality": "xenium",
            "roi": "pln_1",
            "cluster_id": cluster_id,
        }
    )
    dest = out_dir / "assignments" / f"donor_id={donor}" / f"section={section}"
    dest.mkdir(parents=True, exist_ok=True)
    assign.to_parquet(dest / "clusters.parquet", index=False)

    # per-cluster pseudobulk on the reported gene set, via a sparse indicator product
    order = sorted(set(cluster_id))
    lookup = {c: i for i, c in enumerate(order)}
    rows = np.array([lookup[c] for c in cluster_id])
    indic = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.float32), (rows, np.arange(rows.size))),
        shape=(len(order), rows.size),
    )
    sizes = np.asarray(indic.sum(axis=1)).ravel()
    gene_idx = {g: i for i, g in enumerate(var_names)}
    present = [g for g in report_genes if g in gene_idx]
    cols = np.array([gene_idx[g] for g in present])
    sums = np.asarray((indic @ lognorm[:, cols]).todense())
    means = sums / sizes[:, None]
    donor_mean = np.asarray(lognorm[:, cols].mean(axis=0)).ravel()
    det = np.asarray((indic @ (lognorm[:, cols] > 0)).todense()) / sizes[:, None]

    wide = {}
    for j, gene in enumerate(present):
        wide[f"mean__{gene}"] = means[:, j]
        wide[f"enr__{gene}"] = means[:, j] - donor_mean[j]
        wide[f"det__{gene}"] = det[:, j]
    prof = pd.concat(
        [
            pd.DataFrame({"donor_id": donor, "section_key": section, "cluster_id": order,
                          "n_cells": sizes.astype(int),
                          "frac_of_section": sizes / sizes.sum()}),
            pd.DataFrame(wide),
        ],
        axis=1,
    )
    prof.to_parquet(out_dir / "profiles" / f"{section}.parquet", index=False)

    result = {
        "donor_id": donor,
        "section_key": section,
        "path": str(path),
        "status": "clustered",
        "n_cells_in_file": n_in,
        "n_cells_dropped_low_counts": n_dropped,
        "min_counts": min_counts,
        "n_clustered": int(len(cell_ids)),
        "n_genes": int(len(var_names)),
        "genes_reported_present": len(present),
        **meta,
    }
    del adata, lognorm
    gc.collect()
    return result


def xen_eligible(donors) -> tuple[set[str], list[str]]:
    """(eligible, refused) among the donors with a Xenium pLN h5ad.

    Asked of ``multiplex5000.cohort.excluded_on(donor, "xenium")``: 6457's central exclusion
    covers its PhenoCycler data only (maintainer, 2026-10-22), so its two Xenium sections are
    clustered.  The PhenoCycler side reads the ingest, which never wrote it.
    """
    sys.path.insert(0, str(REPO / "core/multiplex5000"))
    from multiplex5000.cohort import excluded_on

    donors = set(donors)
    refused = sorted(d for d in donors if excluded_on(d, "xenium"))
    return donors - set(refused), refused


def run_xen(args) -> None:
    out_dir = OUT_ROOT / "xenium"
    (out_dir / "profiles").mkdir(parents=True, exist_ok=True)

    files = xen_files()
    eligible_donors, refused = xen_eligible({d for d, _, _ in files})

    subpanels = subpanel_genes()
    report_genes = sorted(
        set(g for gs in XEN_REGIONAL.values() for g in gs)
        | set(g for gs in XEN_CANONICAL.values() for g in gs)
        | set(g for gs in subpanels.values() for g in gs)
    )

    import h5py
    panels = {}
    for _, section, path in files:
        with h5py.File(path, "r") as handle:
            panels[section] = set(x.decode() for x in handle["var/_index"][:])
    shared = set.intersection(*panels.values())
    named_absent = [g for g in XEN_NAMED_IN_BRIEF if g not in shared]
    log(f"xenium pLN: {len(files)} h5ads, shared panel {len(shared)} genes; "
        f"named in brief but absent from every file: {named_absent}")

    todo = [(d, s, p) for d, s, p in files if d in eligible_donors]
    clusterable = sorted({d for d, _, _ in todo})
    if args.donors:
        todo = [(d, s, p) for d, s, p in todo if d in set(args.donors)]
    scope = set_scope(SCOPE.refine(clusterable))
    scope.announce(out_dir, log)
    write_scope_note(scope, out_dir, {"stage": "cluster", "modality": "xenium",
                                      "sections_in_this_run": [s for _, s, _ in todo]})
    log(f"clustering {len(todo)} sections across {len({d for d,_,_ in todo})} donors; "
        f"refused by phenocycler.cohort: {refused}")

    results = []
    for donor, section, path in todo:
        started = time.time()
        log(f"  {section}: clustering")
        try:
            res = cluster_xen_section(
                donor, section, path, out_dir, report_genes,
                subsample=args.subsample, resolution=args.resolution,
                min_counts=args.min_counts, seed_base=args.seed,
            )
        except Exception as exc:  # a failing section is recorded, never silently skipped
            res = {"donor_id": donor, "section_key": section, "status": "failed",
                   "error": f"{type(exc).__name__}: {exc}"}
            log(f"  {section}: FAILED {res['error']}")
        res["seconds"] = round(time.time() - started, 1)
        if res.get("status") == "clustered":
            log(f"  {section}: {res['n_clusters']} clusters over {res['n_clustered']} cells, "
                f"k-NN holdout agreement {res['knn_holdout_agreement']:.3f} ({res['seconds']}s)")
        results.append(res)

    report = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": "xenium",
        "source": str(XEN_H5AD_DIR),
        "n_h5ads": len(files),
        "seed_base": int(args.seed),
        "donors_refused_by_cohort": refused,
        "panel": {
            "genes_per_file": 5101,
            "shared_across_all_files": len(shared),
            "named_in_brief_but_absent_from_every_file": named_absent,
            "note": "three genes swap between bundles (CPE/HADH/PCSK2 vs B2M/IAPP/MT2A)",
        },
        "features": {
            "transform": "normalize_total -> log1p -> top-2000 HVG within section -> scale(max 10)",
            "note": "HVGs are selected inside each section; sections are never pooled",
        },
        "qc": {
            "rule": f"nCount_Xenium >= {args.min_counts} and nFeature_Xenium >= 5",
            "nucleus_rule": (
                "unavailable: these h5ads carry no nucleus geometry, so the no-nucleus "
                "exclusion has no per-cell handle here. Recorded as unavailable, not as passed."
            ),
        },
        "labels_never_seen_by_the_clustering": [
            "disease.status", "batch", "AAb.count", "any existing celltype call",
        ],
        "per_section": results,
    }
    cohort_out(out_dir, "cluster_report.json").write_text(json.dumps(report, indent=2))
    log(f"wrote {cohort_out(out_dir, 'cluster_report.json')}")


# ----------------------------------------------------------------------------------
# summarize: enrichment, node relation, cross-donor structure, batch eta-squared
# ----------------------------------------------------------------------------------

def load_node_layer() -> pd.DataFrame:
    frames = []
    for path in sorted(PLN_LAYER.glob("donor_id=*/pln_node_membership.parquet")):
        frames.append(pq.read_table(path).to_pandas())
    layer = pd.concat(frames, ignore_index=True)
    return layer


def meta_cluster(profiles: pd.DataFrame, feature_cols: list[str], k: int) -> pd.Series:
    """Relate per-donor clusters across donors by correlating their marker profiles.

    Average-linkage on a correlation distance.  This is a *description* of which per-donor
    clusters resemble each other; it is not a joint embedding and no cell was ever placed in a
    shared space with a cell from another donor.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist

    mat = profiles[feature_cols].to_numpy(dtype=float)
    mat = np.nan_to_num(mat, nan=0.0)
    dist = pdist(mat, metric="correlation")
    dist = np.nan_to_num(dist, nan=1.0)
    link = linkage(dist, method="average")
    k = int(min(k, mat.shape[0]))
    return pd.Series(fcluster(link, t=k, criterion="maxclust"), index=profiles.index)


def summarize(modality: str, meta_k: int) -> None:
    out_dir = OUT_ROOT / modality
    prof_files = sorted((out_dir / "profiles").glob("*.parquet"))
    if not prof_files:
        raise SystemExit(f"no profiles under {out_dir/'profiles'} -- run the cluster stage first")
    profiles = pd.concat([pd.read_parquet(p) for p in prof_files], ignore_index=True)

    cov = donor_covariates()
    is_pheno = modality == "phenocycler"
    prefix = "z__" if is_pheno else "enr__"
    feature_cols = [c for c in profiles.columns if c.startswith(prefix)]
    # a feature is comparable across donors only where every donor measured it
    per_donor_ok = profiles.groupby("donor_id")[feature_cols].apply(
        lambda g: g.notna().all()
    ).all()
    shared_feats = [c for c in feature_cols if bool(per_donor_ok[c])]

    # ---- regional enrichment table -------------------------------------------------
    groups = PHENO_REGIONAL if is_pheno else {**XEN_REGIONAL, **XEN_CANONICAL}
    enr_rows = []
    for rec in profiles.to_dict("records"):
        entry = {
            "donor_id": rec["donor_id"],
            "section_key": rec.get("section_key", ""),
            "cluster_id": rec["cluster_id"],
            "n_cells": rec["n_cells"],
            "frac_of_section": rec.get("frac_of_donor", rec.get("frac_of_section")),
        }
        for label, members in groups.items():
            vals = [rec[f"{prefix}{m}"] for m in members if f"{prefix}{m}" in rec
                    and rec[f"{prefix}{m}"] == rec[f"{prefix}{m}"]]
            entry[f"set__{label}"] = float(np.mean(vals)) if vals else float("nan")
            entry[f"n_measured__{label}"] = len(vals)
        for member in sorted({m for ms in groups.values() for m in ms}):
            col = f"{prefix}{member}"
            entry[member] = rec.get(col, float("nan"))
        if is_pheno:
            vals = [rec[f"{prefix}{m}"] for m in PHENO_NON_IMMUNE if f"{prefix}{m}" in rec]
            entry["set__non_immune_structural"] = float(np.mean(vals)) if vals else float("nan")
        ranked = sorted(
            ((c[len(prefix):], rec[c]) for c in feature_cols
             if rec.get(c) == rec.get(c)),
            key=lambda kv: kv[1], reverse=True,
        )
        entry["top_features"] = "; ".join(f"{n}={v:+.2f}" for n, v in ranked[:6])
        enr_rows.append(entry)
    enrichment = pd.DataFrame(enr_rows)
    enrichment.to_csv(cohort_out(out_dir, "regional_enrichment.csv"), index=False)
    enrichment.to_parquet(cohort_out(out_dir, "regional_enrichment.parquet"), index=False)

    # ---- cross-donor structure -----------------------------------------------------
    profiles["meta_cluster"] = meta_cluster(profiles, shared_feats, meta_k)
    profiles.to_parquet(cohort_out(out_dir, "cluster_profiles.parquet"), index=False)

    frac_col = "frac_of_donor" if "frac_of_donor" in profiles.columns else "frac_of_section"
    meta_rows = []
    for mid, grp in profiles.groupby("meta_cluster"):
        donors = sorted(set(grp["donor_id"]))
        prof_means = grp[shared_feats].mean()
        top = prof_means.sort_values(ascending=False).head(8)
        meta_rows.append({
            "meta_cluster": int(mid),
            "n_donor_clusters": len(grp),
            "n_donors": len(donors),
            "donors": ";".join(donors),
            "total_cells": int(grp["n_cells"].sum()),
            f"median_{frac_col}": float(grp[frac_col].median()),
            "top_features": "; ".join(
                f"{n[len(prefix):]}={v:+.2f}" for n, v in top.items()
            ),
        })
    meta = pd.DataFrame(meta_rows).sort_values("n_donors", ascending=False)
    meta.to_csv(cohort_out(out_dir, "meta_clusters.csv"), index=False)

    # the same regional enrichment, rolled up to the recurring groups, cell-weighted
    keyed = enrichment.merge(
        profiles[["donor_id", "cluster_id", "meta_cluster"]],
        on=["donor_id", "cluster_id"], how="left",
    )
    set_cols = [c for c in keyed.columns if c.startswith("set__")]
    member_cols = sorted({m for ms in groups.values() for m in ms})
    member_cols = [c for c in member_cols if c in keyed.columns]
    rolled = []
    for mid, grp in keyed.groupby("meta_cluster"):
        w = grp["n_cells"].to_numpy(dtype=float)
        row = {"meta_cluster": int(mid), "n_donors": int(grp["donor_id"].nunique()),
               "n_donor_clusters": int(len(grp)), "total_cells": int(w.sum())}
        for col in set_cols + member_cols:
            v = grp[col].to_numpy(dtype=float)
            ok = np.isfinite(v)
            row[col] = float(np.average(v[ok], weights=w[ok])) if ok.any() else float("nan")
        rolled.append(row)
    pd.DataFrame(rolled).sort_values("total_cells", ascending=False).to_csv(
        cohort_out(out_dir, "meta_cluster_regional_enrichment.csv"), index=False
    )

    # ---- batch eta-squared on everything cross-donor -------------------------------
    comp = (
        profiles.assign(cells=profiles["n_cells"])
        .pivot_table(index="donor_id", columns="meta_cluster", values="cells",
                     aggfunc="sum", fill_value=0)
    )
    comp = comp.div(comp.sum(axis=1), axis=0)
    donors = comp.index.tolist()
    covs = cov.reindex(donors)
    eta_rows = []
    unrecorded: dict[str, list[str]] = {}
    for var in ("stain_group", "acq_run", "coarse_batch", "software", "disease_status"):
        # A donor the manifest does not record is DROPPED from that grouping, never pooled into
        # an "unrecorded" level.  Pooling them would invent a level made of donors that have
        # nothing in common except an absent record, and eta-squared would then partly measure
        # the gap in the manifest.
        have = covs[var].notna().to_numpy()
        unrecorded[var] = [d for d, ok in zip(donors, have) if not ok]
        levels = covs.loc[have, var].to_numpy()
        for mid in comp.columns:
            eta_rows.append({
                "quantity": f"meta_cluster_{int(mid)}_fraction_of_donor",
                "grouping": var,
                "n_donors_used": int(have.sum()),
                "n_donors_unrecorded": int((~have).sum()),
                "n_levels": int(len(set(levels.tolist()))),
                **eta_squared_with_null(comp[mid].to_numpy()[have], levels),
            })
    # and on the leading axis of the cross-donor composition, which is the thing most likely
    # to be read as an ordering
    from sklearn.decomposition import PCA
    if comp.shape[0] >= 3:
        pcs = PCA(n_components=min(3, comp.shape[0] - 1, comp.shape[1]),
                  random_state=stable_seed("pln-cluster", SEED, modality, "composition-pca")
                  ).fit_transform(
            comp.to_numpy() - comp.to_numpy().mean(axis=0))
        for var in ("stain_group", "acq_run", "coarse_batch", "software", "disease_status"):
            have = covs[var].notna().to_numpy()
            levels = covs.loc[have, var].to_numpy()
            for j in range(pcs.shape[1]):
                eta_rows.append({
                    "quantity": f"composition_PC{j+1}",
                    "grouping": var,
                    "n_donors_used": int(have.sum()),
                    "n_donors_unrecorded": int((~have).sum()),
                    "n_levels": int(len(set(levels.tolist()))),
                    **eta_squared_with_null(pcs[have, j], levels),
                })
    eta = pd.DataFrame(eta_rows)
    eta.to_csv(cohort_out(out_dir, "batch_eta_squared.csv"), index=False)
    comp.to_csv(cohort_out(out_dir, "donor_meta_cluster_composition.csv"))
    prov = covariate_provenance(donors)
    prov.to_csv(cohort_out(out_dir, "covariate_provenance.csv"), index=False)

    # ---- node relation (PhenoCycler only) ------------------------------------------
    node_summary = None
    if is_pheno:
        layer = load_node_layer()[["object_id", "donor_id", "node_name", "node_state"]]
        assigns = []
        for path in sorted((out_dir / "assignments").glob("donor_id=*/clusters.parquet")):
            assigns.append(pd.read_parquet(path, columns=["object_id", "donor_id", "cluster_id"]))
        assign = pd.concat(assigns, ignore_index=True)
        joined = assign.merge(layer, on=["object_id", "donor_id"], how="left")
        nodes_per_donor = joined.groupby("donor_id")["node_name"].nunique()
        # Expected share of a cluster that lands in a given node if the cluster were spread in
        # proportion to node size.  Without it, "confined to one node" is untestable: a cluster
        # of 40,000 cells is never literally confined to one node, so the question that can
        # actually be answered is whether it is *over-represented* in one beyond node size.
        node_counts = (
            joined.dropna(subset=["node_name"]).groupby(["donor_id", "node_name"]).size()
        )
        node_share = node_counts / node_counts.groupby(level=0).transform("sum")
        rows = []
        for (donor, cid), grp in joined.groupby(["donor_id", "cluster_id"]):
            inside = grp["node_name"].notna()
            counts = grp.loc[inside, "node_name"].value_counts()
            total = int(counts.sum())
            if total:
                obs = counts / total
                exp = node_share.loc[donor].reindex(obs.index)
                ratio = np.log2(obs.to_numpy() / exp.to_numpy())
                hi = int(np.argmax(ratio))
                max_ratio = float(ratio[hi])
                max_ratio_node = str(obs.index[hi])
                min_ratio = float(ratio.min())
                tvd = float(0.5 * np.abs(
                    obs.reindex(node_share.loc[donor].index).fillna(0).to_numpy()
                    - node_share.loc[donor].to_numpy()
                ).sum())
            else:
                max_ratio, max_ratio_node, tvd = float("nan"), "", float("nan")
                min_ratio = float("nan")
            probs = (counts / total).to_numpy() if total else np.array([])
            ent = float(-(probs * np.log(probs)).sum()) if probs.size else float("nan")
            max_ent = math.log(len(counts)) if len(counts) > 1 else float("nan")
            rows.append({
                "donor_id": donor,
                "cluster_id": cid,
                "n_cells": int(len(grp)),
                "n_cells_inside_a_node": total,
                "n_nodes_on_donor": int(nodes_per_donor.get(donor, 0)),
                "n_nodes_cluster_occupies": int(len(counts)),
                "top_node": counts.index[0] if len(counts) else "",
                "top_node_frac": float(counts.iloc[0] / total) if total else float("nan"),
                "node_entropy_normalized": float(ent / max_ent) if max_ent and max_ent > 0 else float("nan"),
                "max_node_log2_obs_over_expected": max_ratio,
                "max_node_log2_node": max_ratio_node,
                "min_node_log2_obs_over_expected": min_ratio,
                "node_total_variation_vs_node_size": tvd,
                "single_node_donor": bool(nodes_per_donor.get(donor, 0) <= 1),
                "node_is_whole_slide": donor in WHOLE_SLIDE_NODE_DONORS,
            })
        node_summary = pd.DataFrame(rows)
        node_summary = node_summary.merge(
            profiles[["donor_id", "cluster_id", "meta_cluster"]],
            on=["donor_id", "cluster_id"], how="left",
        )
        # both: the CSV for reading, the parquet because every pLN donor id is all digits and a
        # naive CSV reader turns '0112'-style ids into integers
        node_summary.to_csv(cohort_out(out_dir, "node_relation.csv"), index=False)
        node_summary.to_parquet(cohort_out(out_dir, "node_relation.parquet"), index=False)

        # recurrence: a meta-cluster that recurs across nodes and donors is structure; one
        # confined to a single node on a single donor is a slide effect
        rec = []
        comparable = node_summary[
            (~node_summary["node_is_whole_slide"]) & (~node_summary["single_node_donor"])
        ]
        for mid, grp in node_summary.groupby("meta_cluster"):
            sub = comparable[comparable["meta_cluster"] == mid]
            rec.append({
                "meta_cluster": int(mid),
                "n_donors": int(grp["donor_id"].nunique()),
                "n_distinct_nodes_occupied": int(
                    node_summary.loc[node_summary["meta_cluster"] == mid, "n_nodes_cluster_occupies"].sum()
                ),
                "donor_clusters_confined_to_one_node": int(
                    (sub["n_nodes_cluster_occupies"] == 1).sum()
                ),
                "donor_clusters_in_multi_node_donors": int(len(sub)),
                "median_top_node_frac_multi_node_donors": float(sub["top_node_frac"].median())
                if len(sub) else float("nan"),
                "median_max_node_log2_obs_over_expected": float(
                    sub["max_node_log2_obs_over_expected"].median()
                ) if len(sub) else float("nan"),
                "donor_clusters_over_2x_in_one_node": int(
                    (sub["max_node_log2_obs_over_expected"] >= 1.0).sum()
                ),
            })
        pd.DataFrame(rec).to_csv(cohort_out(out_dir, "meta_cluster_node_recurrence.csv"), index=False)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": modality,
        "n_donor_clusters": int(len(profiles)),
        "n_donors": int(profiles["donor_id"].nunique()),
        "total_cells_assigned": int(profiles["n_cells"].sum()),
        "meta_clusters_k": int(profiles["meta_cluster"].nunique()),
        "shared_features_used_for_cross_donor_matching": len(shared_feats),
        "donors_whose_recorded_image_covers_the_pln": sorted(
            prov.loc[prov["an_image_record_covers_the_pln"], "donor_id"].tolist()),
        "donors_whose_pln_acquisition_is_unrecorded": sorted(
            prov.loc[~prov["an_image_record_covers_the_pln"], "donor_id"].tolist()),
        "batch_variable_provenance": (
            "data/donor_acquisition_runs.csv and donor_acquisition_batch.csv record the "
            "PANCREAS scans. No tracked manifest here records the pLN scans, so stain_group / "
            "acq_run / batch are donor-level grouping variables, not a verified pLN scan batch."
        ),
        "eta_squared_by_grouping": {
            var: {
                "n_levels": int(eta.loc[eta["grouping"] == var, "n_levels"].iloc[0]),
                "n_donors_used": int(eta.loc[eta["grouping"] == var, "n_donors_used"].iloc[0]),
                "n_donors_unrecorded": int(eta.loc[eta["grouping"] == var, "n_donors_unrecorded"].iloc[0]),
                "max_observed": float(eta.loc[eta["grouping"] == var, "eta_squared"].max(skipna=True)),
                "median_observed": float(eta.loc[eta["grouping"] == var, "eta_squared"].median(skipna=True)),
                "median_of_random_relabelling_null": float(
                    eta.loc[eta["grouping"] == var, "null_median"].median(skipna=True)
                ),
                "min_permutation_p": float(eta.loc[eta["grouping"] == var, "p_permutation"].min(skipna=True)),
                "quantities_with_p_below_0.05": int(
                    (eta.loc[eta["grouping"] == var, "p_permutation"] < 0.05).sum()
                ),
                "quantities_with_a_saturated_eta_squared": int(
                    eta.loc[eta["grouping"] == var, "eta_squared_saturated"].sum()
                ),
            }
            for var in eta["grouping"].unique()
        },
        "donors_with_no_recorded_value_dropped_from_that_grouping": unrecorded,
        "how_to_read_eta_squared": (
            "eta-squared rises with the number of levels. With 14 donors over 8 stain groups "
            "and 11 acquisition dates, a near-donor-unique grouping partitions almost all the "
            "variance whatever the values are. Compare each observed value against "
            "median_of_random_relabelling_null, and read p_permutation, not the raw number."
        ),
        "node_layer": (
            {
                "source": str(PLN_LAYER),
                "clusters_summarized": int(len(node_summary)),
                "donors_with_more_than_one_node": int(
                    (~node_summary["single_node_donor"]).groupby(
                        node_summary["donor_id"]).first().sum()
                ),
                "excluded_from_per_node_comparison": sorted(WHOLE_SLIDE_NODE_DONORS),
            }
            if node_summary is not None
            else "not applicable: the node layer is a PhenoCycler-frame geometry "
                 "(frame phenocycler_um:<donor>:pln_1) and does not address Xenium cells"
        ),
    }
    cohort_out(out_dir, "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log(f"wrote {cohort_out(out_dir, 'summary.json')}")
    log(json.dumps(summary["eta_squared_by_grouping"], indent=2))



# ----------------------------------------------------------------------------------
# link: relate clusters to the immune calls the typing tracks produce
# ----------------------------------------------------------------------------------

PHENO_CALLS = REPO / "data/development/typing/pln/pheno"
XEN_CALLS = REPO / "data/development/xenium_typing/release_pln_immune"

#: Labels in the Xenium pLN release whose small-core stage-2 calls did not survive the other
#: track's adversarial verification (a majority of them detect none of their own lineage genes).
#: They are carried in the contingency table but named here so no reader treats a cluster that is
#: "overwhelmingly Granulocyte" as a granulocyte result.
XEN_LABELS_FAILED_ADVERSARIAL = ("Granulocyte", "Immune_NK_ILC")


def _cluster_assignments(out_dir: Path, key: str) -> pd.DataFrame:
    paths = sorted((out_dir / "assignments").rglob("clusters.parquet"))
    if not paths:
        raise SystemExit(f"no assignments under {out_dir/'assignments'} -- run --stage cluster")
    cols = [key, "donor_id", "cluster_id"]
    return pd.concat([pd.read_parquet(x, columns=cols) for x in paths], ignore_index=True)


#: What the typing track writes for a cell it did not place in a class.  Kept as an explicit
#: level rather than dropped: a cluster that is mostly uncalled is a fact about the join, and
#: silently dropping those cells would turn a 9%-called cluster into a "100% Myeloid" one.
BLANK_LABEL = "(blank: no class written by the typing track)"


def _as_label(series: pd.Series) -> pd.Series:
    """Any column -> a plain string label, with missing and empty as one explicit level.

    Booleans become "True"/"False"; a null or an empty string becomes `BLANK_LABEL`.  Nothing
    is dropped: a cluster whose cells the typing track left unresolved is a fact about the
    join, and dropping those cells would turn a cluster that is 9% called into a cluster that
    is "100% Myeloid".
    """
    def one(value: object) -> str:
        # pd.NA, np.nan and None all mean "the typing track wrote nothing here"; str(pd.NA)
        # would silently become the literal '<NA>' and read as a label.
        try:
            if value is None or bool(pd.isna(value)):
                return BLANK_LABEL
        except (TypeError, ValueError):
            pass
        return str(value)

    out = series.astype("object").map(one)
    return out.where(out.str.strip().ne(""), BLANK_LABEL)


def _contingency(joined: pd.DataFrame, label_col: str) -> pd.DataFrame:
    """Per cluster: how many cells carry each label, and the dominant label's share.

    A cluster is *described* by the calls that land in it.  Nothing here accepts or rejects a
    call, and no cluster is renamed after its dominant label.  The denominator is every cell in
    the cluster, including the ones with no label.
    """
    work = joined[["donor_id", "cluster_id"]].copy()
    work["label"] = _as_label(joined[label_col]).to_numpy()
    tab = (
        work.groupby(["donor_id", "cluster_id", "label"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    total = tab.groupby(["donor_id", "cluster_id"], observed=True)["n_cells"].transform("sum")
    tab["frac_of_cluster"] = tab["n_cells"] / total
    tab["n_cells_in_cluster"] = total
    tab["label_column"] = label_col
    return tab.sort_values(["donor_id", "cluster_id", "n_cells"], ascending=[True, True, False])


def link_immune_calls(modality: str) -> None:
    """Relate each cluster to whatever immune calls the typing tracks have written.

    This runs on whatever is on disk *now*.  The two typing tracks are separate builds, so the
    report records the source and its own timestamp: a linkage number is only about the release
    it was computed against, and re-running this stage after a typing rebuild is the only way to
    refresh it.  If no usable calls are present the stage says so with the measured reason and
    writes nothing else -- it never invents a label.
    """
    out_dir = OUT_ROOT / modality
    report: dict = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": modality,
        "what_this_is": (
            "a description of which calls land in which cluster. No cluster gates, tunes or "
            "validates any call, and no cluster is renamed after its dominant label."
        ),
    }

    if modality == "phenocycler":
        root, key = PHENO_CALLS, "object_id"
        files = sorted(root.rglob("*.parquet"))
        unavailable = sorted(x.relative_to(root).as_posix() for x in root.rglob("UNAVAILABLE.json"))
        report["source"] = str(root)
        report["parquet_files_found"] = len(files)
        report["donors_marked_unavailable_by_the_typing_track"] = unavailable
        man = root / "manifest.json"
        if man.exists():
            meta = json.loads(man.read_text())
            report["source_manifest_generated_at_utc"] = meta.get("generated_at_utc")
            report["source_method_as_the_typing_track_states_it"] = meta.get("method")
            report["source_manifest_mtime_utc"] = datetime.fromtimestamp(
                man.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    else:
        root, key = XEN_CALLS, "cell_id"
        files = sorted((root / "calls").rglob("*.parquet"))
        report["source"] = str(root)
        report["parquet_files_found"] = len(files)
        rel = root / "release.json"
        if rel.exists():
            meta = json.loads(rel.read_text())
            report["source_release_generated_at_utc"] = meta.get("generated_at_utc")
        report["labels_that_failed_the_other_track_adversarial_verification"] = list(
            XEN_LABELS_FAILED_ADVERSARIAL
        )

    if not files:
        report["status"] = "no_calls_available"
        report["note"] = (
            "The typing track has written no per-cell calls under this root, so there is "
            "nothing to relate. This is recorded as unavailable, not as zero agreement."
        )
        cohort_out(out_dir, "cluster_immune_link.json").write_text(json.dumps(report, indent=2))
        log(f"link: {modality}: no calls under {root} -- recorded as unavailable")
        return

    frames = []
    for path in files:
        try:
            frame = pd.read_parquet(path)
        except Exception as exc:
            report.setdefault("unreadable", []).append(f"{path}: {type(exc).__name__}: {exc}")
            continue
        if key not in frame.columns:
            report.setdefault("skipped_no_key", []).append(str(path))
            continue
        if "donor_id" not in frame.columns:
            # hive-partitioned: donor_id is in the path, not the file
            token = [q for q in path.parts if q.startswith("donor_id=")]
            if not token:
                report.setdefault("skipped_no_donor", []).append(str(path))
                continue
            frame["donor_id"] = token[0].split("=", 1)[1]
        frames.append(frame)
    if not frames:
        report["status"] = "no_usable_calls"
        cohort_out(out_dir, "cluster_immune_link.json").write_text(json.dumps(report, indent=2))
        log(f"link: {modality}: files present but none usable -- recorded as unavailable")
        return

    calls = pd.concat(frames, ignore_index=True)
    calls[key] = calls[key].astype(str)
    calls["donor_id"] = calls["donor_id"].astype(str)
    label_cols = [c for c in ("stage1", "is_immune", "lineage_native", "celltype_fine",
                              "cell_class", "immune_class", "immune_subtype", "class")
                  if c in calls.columns]
    report["label_columns_found"] = label_cols
    report["calls_rows"] = int(len(calls))

    assign = _cluster_assignments(out_dir, key)
    assign[key] = assign[key].astype(str)
    assign["donor_id"] = assign["donor_id"].astype(str)
    joined = assign.merge(calls[["donor_id", key] + label_cols], on=["donor_id", key], how="left")
    report["clustered_cells"] = int(len(assign))
    report["first_label_column"] = label_cols[0]
    report["clustered_cells_with_a_named_label"] = {
        col: int(_as_label(joined[col]).ne(BLANK_LABEL).sum()) for col in label_cols
    }
    report["clustered_cells_with_a_call"] = max(
        report["clustered_cells_with_a_named_label"].values()
    )
    report["donors_in_both"] = sorted(
        set(assign["donor_id"]) & set(calls["donor_id"])
    )
    report["donors_clustered_but_not_called"] = sorted(
        set(assign["donor_id"]) - set(calls["donor_id"])
    )

    if not label_cols or report["clustered_cells_with_a_call"] == 0:
        report["status"] = "no_key_overlap"
        report["note"] = (
            f"No clustered cell joined a call on ({key}, donor_id). The two builds do not "
            "share a key here; nothing is inferred from that."
        )
        cohort_out(out_dir, "cluster_immune_link.json").write_text(json.dumps(report, indent=2))
        log(f"link: {modality}: no key overlap -- recorded as unavailable")
        return

    tabs = pd.concat([_contingency(joined, c) for c in label_cols], ignore_index=True)
    prof = out_dir / "cluster_profiles.parquet"
    if prof.exists():
        tabs = tabs.merge(
            pd.read_parquet(prof, columns=["donor_id", "cluster_id", "meta_cluster"]),
            on=["donor_id", "cluster_id"], how="left",
        )
    tabs.to_csv(cohort_out(out_dir, "cluster_immune_link.csv"), index=False)

    # how concentrated each cluster is in one label, per label column
    conc = {}
    for col in label_cols:
        sub = tabs[tabs["label_column"] == col]
        top = sub.sort_values("n_cells", ascending=False).groupby(
            ["donor_id", "cluster_id"], observed=True).first()
        named = sub[sub["label"] != BLANK_LABEL]
        top_named = named.sort_values("n_cells", ascending=False).groupby(
            ["donor_id", "cluster_id"], observed=True).first()
        conc[col] = {
            "n_clusters": int(len(top)),
            "median_dominant_label_share": float(top["frac_of_cluster"].median()),
            "clusters_at_least_80pct_one_label": int((top["frac_of_cluster"] >= 0.80).sum()),
            "clusters_at_least_90pct_one_label": int((top["frac_of_cluster"] >= 0.90).sum()),
            "dominant_label_counts": {
                str(k): int(v) for k, v in top["label"].value_counts().items()},
            "clusters_whose_dominant_label_is_blank": int(
                (top["label"] == BLANK_LABEL).sum()),
            "ignoring_blanks": {
                "dominant_named_label_counts": {
                    str(k): int(v) for k, v in top_named["label"].value_counts().items()},
                "median_named_share_of_cluster": float(top_named["frac_of_cluster"].median())
                    if len(top_named) else float("nan"),
                "clusters_at_least_50pct_one_named_label": int(
                    (top_named["frac_of_cluster"] >= 0.50).sum()),
            },
        }
    report["status"] = "linked"
    report["concentration"] = conc
    cohort_out(out_dir, "cluster_immune_link.json").write_text(json.dumps(report, indent=2))
    log(f"link: {modality}: wrote {out_dir/'cluster_immune_link.csv'}")
    log(json.dumps(conc, indent=2))



# ----------------------------------------------------------------------------------
# prose: every number the write-up cites, emitted from the generated files
# ----------------------------------------------------------------------------------

def prose_facts(modality: str) -> None:
    """Emit `prose_facts.json`: the numbers `docs/PLN_REGIONAL_CLUSTERING_2026-09-17.md` quotes.

    The previous write-up carried numbers no command produced.  This stage exists so the prose
    can be written *from* the generated files rather than beside them: anything in the document
    that is not in this file, or in one of the CSVs it points at, does not belong there.
    """
    out_dir = OUT_ROOT / modality
    is_pheno = modality == "phenocycler"
    prefix = "z__" if is_pheno else "enr__"
    profiles = pd.read_parquet(out_dir / "cluster_profiles.parquet")
    report = json.loads((out_dir / "cluster_report.json").read_text())
    eta = pd.read_csv(out_dir / "batch_eta_squared.csv")
    per = report["per_donor"] if is_pheno else report["per_section"]
    done = [r for r in per if r.get("status") == "clustered"]

    facts: dict = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": modality,
        "cluster_report_generated_at_utc": report["generated_at_utc"],
        "run": {
            "n_units_clustered": len(done),
            "n_donors": int(profiles["donor_id"].nunique()),
            "n_clusters": int(len(profiles)),
            "n_cells_clustered": int(sum(r["n_clustered"] for r in done)),
            "clusters_per_unit_min": int(min(r["n_clusters"] for r in done)),
            "clusters_per_unit_max": int(max(r["n_clusters"] for r in done)),
            "knn_holdout_agreement_min": float(min(r["knn_holdout_agreement"] for r in done)),
            "knn_holdout_agreement_median": float(
                np.median([r["knn_holdout_agreement"] for r in done])),
            "knn_holdout_agreement_max": float(max(r["knn_holdout_agreement"] for r in done)),
            "seeds": {str(r.get("donor_id") if is_pheno else r.get("section_key")): int(r["seed"])
                      for r in done},
        },
    }

    if is_pheno:
        facts["source"] = {
            "root": report["source"],
            "rows_in_root": report["rows_in_root"],
            "rows_analysis_eligible": report["rows_analysis_eligible"],
            "donors_in_root": report["donors_in_root"],
            "donors_clustered": report["donors_clustered"],
            "donors_with_no_intensity_measurements": report["donors_with_no_intensity_measurements"],
        }
        facts["measurement_handling"] = {
            "cells_dropped_partial_measurement_total": int(
                sum(r["cells_dropped_partial_measurement"] for r in done)),
            "n_markers_in_export": len(report["features"]["markers_in_export"]),
            "n_markers_used_min": int(min(len(r["markers_used"]) for r in done)),
            "n_markers_used_max": int(max(len(r["markers_used"]) for r in done)),
            "markers_unavailable_anywhere": sorted(
                {m for r in done for m in r["markers_unavailable_on_this_slide"]}),
            "markers_zero_mad_anywhere": sorted(
                {m for r in done for m in r["markers_zero_mad_forced_to_zero_z"]}),
        }
        # how far each marker can separate a cluster, and how widely that recurs
        zcols = [c for c in profiles.columns if c.startswith("z__")]
        rows = []
        for col in zcols:
            marker = col[len("z__"):]
            vals = profiles[col]
            hits = profiles.loc[vals > 2.0]
            rows.append({
                "marker": marker,
                "max_cluster_z": float(vals.max()),
                "clusters_z_over_2": int(len(hits)),
                "donors_with_a_z_over_2_cluster": int(hits["donor_id"].nunique()),
            })
        sep = pd.DataFrame(rows).sort_values("max_cluster_z", ascending=False)
        sep.to_csv(cohort_out(out_dir, "marker_separation.csv"), index=False)
        facts["marker_separation_csv"] = str(out_dir / "marker_separation.csv")
        facts["marker_separation"] = sep.to_dict("records")
        # the same statement without any scaling: top cluster raw mean over donor pooled mean
        ratios = []
        for marker in sep["marker"]:
            col = f"raw__{marker}"
            if col not in profiles.columns:
                continue
            per_donor = []
            for donor, grp in profiles.groupby("donor_id"):
                w = grp["n_cells"].to_numpy(dtype=float)
                v = grp[col].to_numpy(dtype=float)
                ok = np.isfinite(v)
                if not ok.any():
                    continue
                pooled = float(np.average(v[ok], weights=w[ok]))
                if pooled > 0:
                    per_donor.append(float(np.nanmax(v[ok]) / pooled))
            if per_donor:
                ratios.append({"marker": marker,
                               "median_top_cluster_over_slide_mean": float(np.median(per_donor)),
                               "n_donors": len(per_donor)})
        rat = pd.DataFrame(ratios).sort_values(
            "median_top_cluster_over_slide_mean", ascending=False)
        rat.to_csv(cohort_out(out_dir, "marker_raw_dynamic_range.csv"), index=False)
        facts["marker_raw_dynamic_range"] = rat.to_dict("records")
    else:
        facts["source"] = {
            "root": report["source"],
            "n_h5ads": report["n_h5ads"],
            "donors_refused_by_cohort": report["donors_refused_by_cohort"],
            "shared_across_all_files": report["panel"]["shared_across_all_files"],
            "named_in_brief_but_absent_from_every_file":
                report["panel"]["named_in_brief_but_absent_from_every_file"],
        }
        facts["measurement_handling"] = {
            "cells_in_files_total": int(sum(r["n_cells_in_file"] for r in done)),
            "cells_dropped_low_counts_total": int(
                sum(r["n_cells_dropped_low_counts"] for r in done)),
            "qc_rule": report["qc"]["rule"],
            "nucleus_rule": report["qc"]["nucleus_rule"],
        }
        enrichment = pd.read_csv(out_dir / "regional_enrichment.csv")
        genes = sorted({g for gs in {**XEN_REGIONAL, **XEN_CANONICAL}.values() for g in gs})
        rows = []
        for gene in genes:
            mcol, dcol, ecol = f"mean__{gene}", f"det__{gene}", f"enr__{gene}"
            if dcol not in profiles.columns:
                continue
            w = profiles["n_cells"].to_numpy(dtype=float)
            d = profiles[dcol].to_numpy(dtype=float)
            ok = np.isfinite(d)
            hits = profiles.loc[profiles[ecol] > 0.25] if ecol in profiles.columns else profiles.iloc[:0]
            rows.append({
                "gene": gene,
                "detected_in_frac_of_all_clustered_cells": float(np.average(d[ok], weights=w[ok])),
                "best_cluster_detection": float(np.nanmax(d)),
                "max_cluster_enrichment_over_section_mean": float(profiles[ecol].max())
                    if ecol in profiles.columns else float("nan"),
                "clusters_enriched_over_0.25": int(len(hits)),
                "sections_with_such_a_cluster": int(hits["section_key"].nunique())
                    if len(hits) else 0,
            })
        gt = pd.DataFrame(rows).sort_values(
            "max_cluster_enrichment_over_section_mean", ascending=False)
        gt.to_csv(cohort_out(out_dir, "gene_detection_and_enrichment.csv"), index=False)
        facts["gene_detection_and_enrichment"] = gt.to_dict("records")
        facts["n_sections"] = int(profiles["section_key"].nunique())

    # ---- recurring groups ---------------------------------------------------------
    meta = pd.read_csv(out_dir / "meta_clusters.csv")
    n_units = int(profiles["donor_id"].nunique())
    recur = meta[meta["n_donors"] == n_units]
    facts["meta_clusters"] = {
        "k": int(len(meta)),
        "n_donors_in_build": n_units,
        "groups_present_in_every_donor": int(len(recur)),
        "cells_in_those_groups": int(recur["total_cells"].sum()),
        "frac_of_clustered_cells_in_those_groups": float(
            recur["total_cells"].sum() / meta["total_cells"].sum()),
        "table": meta.to_dict("records"),
    }

    # ---- node relation ------------------------------------------------------------
    node_path = out_dir / "node_relation.parquet"
    if node_path.exists():
        node = pd.read_parquet(node_path)
        comparable = node[(~node["node_is_whole_slide"]) & (~node["single_node_donor"])]
        facts["node_relation"] = {
            "clusters_summarized": int(len(node)),
            "cells_summarized": int(node["n_cells"].sum()),
            "cells_inside_a_node": int(node["n_cells_inside_a_node"].sum()),
            "donors_with_one_node": sorted(
                node.loc[node["single_node_donor"], "donor_id"].unique().tolist()),
            "donors_whose_node_is_the_slide": sorted(
                node.loc[node["node_is_whole_slide"], "donor_id"].unique().tolist()),
            "comparable_clusters": int(len(comparable)),
            "comparable_donors": int(comparable["donor_id"].nunique()),
            "clusters_confined_to_one_node": int(
                (comparable["n_nodes_cluster_occupies"] == 1).sum()),
            "median_max_node_log2_obs_over_expected": float(
                comparable["max_node_log2_obs_over_expected"].median()),
            "clusters_over_1.5x_in_one_node": int(
                (comparable["max_node_log2_obs_over_expected"] >= math.log2(1.5)).sum()),
            "clusters_over_2x_in_one_node": int(
                (comparable["max_node_log2_obs_over_expected"] >= 1.0).sum()),
            "largest_max_node_log2_obs_over_expected": float(
                comparable["max_node_log2_obs_over_expected"].max()),
            "median_total_variation_vs_node_size": float(
                comparable["node_total_variation_vs_node_size"].median()),
        }

    # ---- batch: every cross-donor quantity, every grouping -------------------------
    facts["batch"] = {
        "quantities": sorted(eta["quantity"].unique().tolist()),
        "groupings": sorted(eta["grouping"].unique().tolist()),
        "per_grouping": {
            var: {
                "n_levels": int(sub["n_levels"].iloc[0]),
                "n_donors_used": int(sub["n_donors_used"].iloc[0]),
                "n_donors_unrecorded": int(sub["n_donors_unrecorded"].iloc[0]),
                "max_eta_squared": float(sub["eta_squared"].max()),
                "quantity_at_max": str(sub.loc[sub["eta_squared"].idxmax(), "quantity"]),
                "median_null": float(sub["null_median"].median()),
                "min_p": float(sub["p_permutation"].min()),
                "n_quantities": int(len(sub)),
                "n_quantities_p_below_0.05": int((sub["p_permutation"] < 0.05).sum()),
                "n_quantities_saturated": int(sub["eta_squared_saturated"].sum()),
            }
            for var, sub in eta.groupby("grouping")
        },
        "significant": eta[eta["p_permutation"] < 0.05].sort_values(
            "p_permutation").to_dict("records"),
        "composition_pcs": eta[eta["quantity"].str.startswith("composition_PC")].to_dict("records"),
    }
    # the full quantity x grouping grid for the three groupings that must appear in the text
    grid = eta[eta["grouping"].isin(["acq_run", "stain_group", "coarse_batch"])].pivot_table(
        index="quantity", columns="grouping", values=["eta_squared", "p_permutation"])
    grid.columns = [f"{a}__{b}" for a, b in grid.columns]
    grid = grid.reset_index()
    grid.to_csv(cohort_out(out_dir, "batch_eta_squared_required_grid.csv"), index=False)
    facts["batch"]["required_grid"] = grid.to_dict("records")

    prov_path = out_dir / "covariate_provenance.csv"
    if prov_path.exists():
        prov = pd.read_csv(prov_path, dtype={"donor_id": str})
        facts["batch"]["donors_whose_recorded_image_covers_the_pln"] = sorted(
            prov.loc[prov["an_image_record_covers_the_pln"], "donor_id"].tolist())
        facts["batch"]["donors_whose_pln_acquisition_is_unrecorded"] = sorted(
            prov.loc[~prov["an_image_record_covers_the_pln"], "donor_id"].tolist())

    link_json = out_dir / "cluster_immune_link.json"
    if link_json.exists():
        facts["immune_link"] = json.loads(link_json.read_text())

    cohort_out(out_dir, "prose_facts.json").write_text(json.dumps(facts, indent=2, default=str))
    log(f"wrote {cohort_out(out_dir, 'prose_facts.json')}")


# ----------------------------------------------------------------------------------
# naming profiles: the table the naming track consumes
# ----------------------------------------------------------------------------------
#
# The naming step of the maintainer's method is: *name each cluster by canonical marker/gene
# enrichment relative to the other clusters on that same section*.  This table is built for
# exactly that consumer and for nothing else, so its contract is stated rather than implied:
#
#   med__<M>      median of the marker in the cluster, on the measured scale (PhenoCycler:
#                 raw cell-mean intensity; Xenium: log-normalized expression)
#   mad__<M>      1.4826 x median absolute deviation of the same, within the cluster
#   iqr__<M>      p75 - p25 of the same, within the cluster
#   logmed__<M>   (PhenoCycler) median of log1p(raw), the scale the contrast is taken on
#   mean__<M>     (Xenium) cluster mean log-normalized expression
#   det__<M>      (Xenium) fraction of the cluster's cells with a non-zero count
#   vsrest__<M>   THE NAMING CONTRAST: this cluster's central value minus the median of the
#                 OTHER clusters' central values on the same section
#   sepz__<M>     the same contrast divided by 1.4826 x MAD of the other clusters' central
#                 values -- a dimensionless "how far does this cluster stand out from its
#                 section's other clusters".  NaN, never zero, when that MAD is zero (fewer
#                 than a few distinct clusters), because an unavailable spread is unavailable.
#   othermad__<M> sepz's denominator, written beside it so the ratio can be read honestly: a
#                 large sepz on a near-zero denominator is not separation.  Measured: AICDA
#                 reaches sepz 49.4 on a vsrest of +0.011.  Read sepz with vsrest, never alone.
#
# What this table deliberately does not contain: a threshold, a positivity call, a reference
# level, a quantile computed on any other section, and any number of pancreas origin.  Every
# quantity is a comparison between clusters measured in the same channel on the same section,
# which is what makes it calibration-free across slides.
#
# `mad__`/`iqr__` are the heterogeneity handle.  A cluster whose CD20 median is high but whose
# CD20 IQR is as wide as the section's is a mixed cluster; the naming track is expected to say
# so rather than hide it behind a majority label.

def _cluster_vs_other_clusters(
    central: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each cluster, its central value against the other clusters' central values.

    `central` is (n_clusters, n_features).  Returns (vsrest, sepz, othermad).  The reference
    is the median over the OTHER clusters -- never over the cells, which would be dominated
    by whichever cluster happens to be largest, and never over another section.

    `othermad` is `sepz`'s denominator and is returned so it can be written beside it.  It has
    to be readable, because `sepz` is a ratio and a ratio with a near-zero denominator is
    large for a reason that is not separation: measured here, AICDA reaches sepz 49.4 on a
    `vsrest` of +0.011 in log-normalized expression, because a gene detected in 0.08% of cells
    has essentially no spread between clusters to divide by.  **`sepz` is never read alone.**
    """
    n_clusters = central.shape[0]
    vsrest = np.full_like(central, np.nan, dtype=np.float64)
    sepz = np.full_like(central, np.nan, dtype=np.float64)
    othermad = np.full_like(central, np.nan, dtype=np.float64)
    if n_clusters < 2:
        return vsrest, sepz, othermad
    for i in range(n_clusters):
        others = np.delete(central, i, axis=0)
        ref = np.nanmedian(others, axis=0)
        spread = np.nanmedian(np.abs(others - ref), axis=0) * 1.4826
        vsrest[i] = central[i] - ref
        othermad[i] = spread
        with np.errstate(invalid="ignore", divide="ignore"):
            sepz[i] = np.where(spread > 0, (central[i] - ref) / spread, np.nan)
    return vsrest, sepz, othermad


def naming_coherence(modality: str, profiles: pd.DataFrame) -> dict:
    """Internal-consistency checks on the assignments and the naming table.

    These are COHERENCE checks, not expectation checks.  None of them passes by reproducing a
    biological finding and none of them would fail because a lymph node contained the wrong
    thing.  They ask only whether the build contradicts itself: whether every clustered cell
    is counted exactly once, whether the per-cluster sizes add up to the cells that were
    clustered, whether an ineligible cell reached a cluster, and whether the naming contrast
    sits near zero across a section's clusters (it is a value minus the median of the others,
    so the median over clusters lands near zero unless that section's clusters are strongly
    asymmetric on that marker).  That last one is an indicator, not a law: the measured
    maximum is a real property of the data, not a tolerance to be met.

    `assigned_key_is_unique` is reported over the whole modality and
    `assigned_key_is_unique_within_unit` within each donor or section.  For Xenium the first
    is expected to be false -- the bundles number their cells independently, so `cell_id` is
    unique within a section and the join key is (section_key, cell_id).  Stating both keeps
    that from reading as a defect or from being missed.
    """
    key = "object_id" if modality == "phenocycler" else "cell_id"
    unit = "donor_id" if modality == "phenocycler" else "section_key"
    out_dir = OUT_ROOT / modality
    frames = [
        pd.read_parquet(path, columns=[key, unit, "cluster_id"])
        for path in sorted((out_dir / "assignments").rglob("clusters.parquet"))
    ]
    assign = pd.concat(frames, ignore_index=True)

    sizes = assign.groupby("cluster_id", observed=True).size()
    prof_sizes = profiles.set_index("cluster_id")["n_cells"]
    common = sizes.index.intersection(prof_sizes.index)

    frac = profiles.groupby(unit, observed=True)["frac_of_section"].sum()
    contrast_cols = [c for c in profiles.columns if c.startswith("vsrest__")]
    centred = {}
    for unit_value, block in profiles.groupby(unit, observed=True):
        if len(block) < 2:
            continue
        centred[str(unit_value)] = float(
            np.nanmax(np.abs(block[contrast_cols].median(axis=0).to_numpy()))
        )

    checks = {
        "assigned_rows": int(len(assign)),
        "assigned_key_is_unique": bool(assign[key].is_unique),
        "assigned_key_is_unique_within_unit": bool(
            not assign.duplicated(subset=[unit, key]).any()
        ),
        "join_key": [unit, key] if not assign[key].is_unique else [key],
        "clusters_in_assignments": int(assign["cluster_id"].nunique()),
        "clusters_in_naming_profiles": int(len(profiles)),
        "cluster_ids_match": bool(
            set(assign["cluster_id"].unique()) == set(profiles["cluster_id"])
        ),
        "max_abs_size_disagreement": int((sizes[common] - prof_sizes[common]).abs().max()),
        "profile_sizes_sum": int(profiles["n_cells"].sum()),
        "max_abs_frac_of_section_sum_minus_one": float((frac - 1.0).abs().max()),
        "max_abs_median_of_vsrest_over_clusters": (
            float(max(centred.values())) if centred else None
        ),
        "unit_at_that_maximum": (
            max(centred, key=centred.get) if centred else None
        ),
        "units_with_a_single_cluster_contrast_unavailable": sorted(
            str(u) for u, n in profiles.groupby(unit, observed=True).size().items() if n < 2
        ),
    }
    if modality == "phenocycler":
        part = pheno_partitioning()
        dataset = ds.dataset(PHENO_CELLS, format="parquet", partitioning=part)
        eligible = dataset.to_table(
            columns=["object_id"], filter=ds.field("qc_analysis_eligible")
        ).column("object_id").to_pandas().astype(str)
        checks["analysis_eligible_cells_in_root"] = int(len(eligible))
        checks["assigned_cells_not_analysis_eligible"] = int(
            len(set(assign[key].astype(str)) - set(eligible))
        )
    return checks


def naming_profile_pheno(args) -> None:
    out_dir = OUT_ROOT / "phenocycler"
    (out_dir / "naming_profiles").mkdir(parents=True, exist_ok=True)
    part = pheno_partitioning()
    dataset = ds.dataset(PHENO_CELLS, format="parquet", partitioning=part)
    markers = [m for m in pheno_marker_columns(dataset) if m not in PHENO_FEATURE_EXCLUDE]

    donors = sorted(p.name.split("=", 1)[1]
                    for p in (out_dir / "assignments").glob("donor_id=*"))
    assigned = list(donors)
    if args.donors:
        donors = [d for d in donors if d in set(args.donors)]
    scope = set_scope(SCOPE.refine(assigned))
    scope.announce(out_dir, log)
    write_scope_note(scope, out_dir, {"stage": "profile", "modality": "phenocycler"})
    log(f"naming profiles, phenocycler: {len(donors)} donors")

    frames = []
    for donor in donors:
        assign = pd.read_parquet(
            out_dir / "assignments" / f"donor_id={donor}" / "clusters.parquet",
            columns=["object_id", "cluster_id"],
        )
        assign["object_id"] = assign["object_id"].astype(str)
        table = dataset.to_table(
            columns=["object_id"] + markers,
            filter=(ds.field("donor_id") == donor) & ds.field("qc_analysis_eligible"),
        )
        frame = table.to_pandas()
        del table
        gc.collect()
        frame["object_id"] = frame["object_id"].astype(str)
        # joined on the QuPath UUID, never on row order, filename or centroid
        merged = frame.merge(assign, on="object_id", how="inner", validate="one_to_one")
        if len(merged) != len(assign):
            raise SystemExit(
                f"{donor}: {len(assign)} assigned cells but {len(merged)} joined by object_id"
            )
        del frame
        gc.collect()

        values = merged[markers].to_numpy(dtype=np.float64)
        # a marker with no measurement on this slide is UNAVAILABLE; it leaves the table for
        # this donor rather than becoming a zero
        avail = ~np.isnan(values).all(axis=0)
        keep = [m for m, a in zip(markers, avail) if a]
        values = values[:, avail]
        logv = np.log1p(np.clip(values, 0.0, None))
        cid = merged["cluster_id"].to_numpy().astype(str)
        del merged
        gc.collect()

        order = sorted(set(cid))
        n_c, n_f = len(order), len(keep)
        med = np.empty((n_c, n_f)); mad = np.empty((n_c, n_f))
        iqr = np.empty((n_c, n_f)); logmed = np.empty((n_c, n_f))
        sizes = np.empty(n_c, dtype=int)
        for i, c in enumerate(order):
            sel = cid == c
            sizes[i] = int(sel.sum())
            block = values[sel]
            q = np.nanpercentile(block, [25, 50, 75], axis=0)
            med[i] = q[1]
            iqr[i] = q[2] - q[0]
            mad[i] = np.nanmedian(np.abs(block - q[1]), axis=0) * 1.4826
            logmed[i] = np.nanmedian(logv[sel], axis=0)
        vsrest, sepz, othermad = _cluster_vs_other_clusters(logmed)

        wide = {}
        for j, m in enumerate(keep):
            wide[f"med__{m}"] = med[:, j]
            wide[f"mad__{m}"] = mad[:, j]
            wide[f"iqr__{m}"] = iqr[:, j]
            wide[f"logmed__{m}"] = logmed[:, j]
            wide[f"vsrest__{m}"] = vsrest[:, j]
            wide[f"sepz__{m}"] = sepz[:, j]
            wide[f"othermad__{m}"] = othermad[:, j]
        prof = pd.concat(
            [
                pd.DataFrame({
                    "donor_id": donor,
                    "section_key": f"{donor}:pln_1",
                    "cluster_id": order,
                    "n_cells": sizes,
                    "frac_of_section": sizes / sizes.sum(),
                    "markers_unavailable_on_this_slide": ";".join(
                        m for m, a in zip(markers, avail) if not a
                    ),
                }),
                pd.DataFrame(wide),
            ],
            axis=1,
        )
        prof.to_parquet(out_dir / "naming_profiles" / f"{donor}.parquet", index=False)
        frames.append(prof)
        log(f"  {donor}: {n_c} clusters x {n_f} markers")
        del values, logv
        gc.collect()

    combined = pd.concat(frames, ignore_index=True)
    combined.to_parquet(cohort_out(out_dir, "naming_profiles.parquet"), index=False)
    log(f"wrote {cohort_out(out_dir, 'naming_profiles.parquet')} ({len(combined)} clusters)")
    if not args.donors:
        checks = naming_coherence("phenocycler", combined)
        cohort_out(out_dir, "naming_profiles_coherence.json").write_text(json.dumps(checks, indent=2))
        log(f"coherence: {json.dumps(checks)}")


def naming_profile_xen(args) -> None:
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp

    out_dir = OUT_ROOT / "xenium"
    (out_dir / "naming_profiles").mkdir(parents=True, exist_ok=True)

    subpanels = subpanel_genes()
    report_genes = sorted(
        set(g for gs in XEN_REGIONAL.values() for g in gs)
        | set(g for gs in XEN_CANONICAL.values() for g in gs)
        | set(g for gs in subpanels.values() for g in gs)
    )

    files = {s: (d, p) for d, s, p in xen_files()}
    sections = sorted(
        p.name.split("=", 1)[1]
        for p in (out_dir / "assignments").glob("donor_id=*/section=*")
    )
    assigned_donors = sorted({files[s][0] for s in sections})
    if args.donors:
        wanted = set(args.donors)
        sections = [s for s in sections if files[s][0] in wanted]
    scope = set_scope(SCOPE.refine(assigned_donors))
    scope.announce(out_dir, log)
    write_scope_note(scope, out_dir, {"stage": "profile", "modality": "xenium",
                                      "sections_in_this_run": list(sections)})
    log(f"naming profiles, xenium: {len(sections)} sections")

    frames = []
    for section in sections:
        donor, path = files[section]
        assign = pd.read_parquet(
            out_dir / "assignments" / f"donor_id={donor}" / f"section={section}"
            / "clusters.parquet",
            columns=["cell_id", "cluster_id"],
        )
        assign["cell_id"] = assign["cell_id"].astype(str)

        adata = ad.read_h5ad(path)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        cell_ids = adata.obs_names.to_numpy().astype(str)
        lognorm = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
        var_names = adata.var_names.to_numpy().astype(str)
        del adata
        gc.collect()

        # joined on cell_id, preserved verbatim -- never on row order
        pos = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
        take = pos.reindex(assign["cell_id"])
        if take.isna().any():
            raise SystemExit(f"{section}: {int(take.isna().sum())} assigned cell_ids not in the h5ad")
        take = take.to_numpy(dtype=int)
        lognorm = lognorm[take]
        cid = assign["cluster_id"].to_numpy().astype(str)

        gene_idx = {g: i for i, g in enumerate(var_names)}
        present = [g for g in report_genes if g in gene_idx]
        cols = np.array([gene_idx[g] for g in present])
        sub = lognorm[:, cols]
        del lognorm
        gc.collect()

        order = sorted(set(cid))
        n_c, n_f = len(order), len(present)
        med = np.empty((n_c, n_f)); mad = np.empty((n_c, n_f)); iqr = np.empty((n_c, n_f))
        mean = np.empty((n_c, n_f)); det = np.empty((n_c, n_f))
        sizes = np.empty(n_c, dtype=int)
        for i, c in enumerate(order):
            sel = np.flatnonzero(cid == c)
            sizes[i] = sel.size
            block = np.asarray(sub[sel].todense(), dtype=np.float32)
            q = np.percentile(block, [25, 50, 75], axis=0)
            med[i] = q[1]
            iqr[i] = q[2] - q[0]
            mad[i] = np.median(np.abs(block - q[1]), axis=0) * 1.4826
            mean[i] = block.mean(axis=0)
            det[i] = (block > 0).mean(axis=0)
            del block
            gc.collect()
        # the contrast is taken on the cluster mean: Xenium counts are sparse enough that the
        # median of most genes is zero in every cluster, so a median-based contrast would be
        # identically zero and say nothing.  The medians and IQRs stay in the table as the
        # honest statement of that sparsity.
        vsrest, sepz, othermad = _cluster_vs_other_clusters(mean)
        det_vsrest, _, _ = _cluster_vs_other_clusters(det)

        wide = {}
        for j, g in enumerate(present):
            wide[f"med__{g}"] = med[:, j]
            wide[f"mad__{g}"] = mad[:, j]
            wide[f"iqr__{g}"] = iqr[:, j]
            wide[f"mean__{g}"] = mean[:, j]
            wide[f"det__{g}"] = det[:, j]
            wide[f"vsrest__{g}"] = vsrest[:, j]
            wide[f"sepz__{g}"] = sepz[:, j]
            wide[f"othermad__{g}"] = othermad[:, j]
            wide[f"detvsrest__{g}"] = det_vsrest[:, j]
        prof = pd.concat(
            [
                pd.DataFrame({
                    "donor_id": donor,
                    "section_key": section,
                    "cluster_id": order,
                    "n_cells": sizes,
                    "frac_of_section": sizes / sizes.sum(),
                    "genes_reported_present": n_f,
                }),
                pd.DataFrame(wide),
            ],
            axis=1,
        )
        prof.to_parquet(out_dir / "naming_profiles" / f"{section}.parquet", index=False)
        frames.append(prof)
        log(f"  {section}: {n_c} clusters x {n_f} genes")
        del sub
        gc.collect()

    combined = pd.concat(frames, ignore_index=True)
    combined.to_parquet(cohort_out(out_dir, "naming_profiles.parquet"), index=False)
    log(f"wrote {cohort_out(out_dir, 'naming_profiles.parquet')} ({len(combined)} clusters)")
    if not args.donors:
        checks = naming_coherence("xenium", combined)
        cohort_out(out_dir, "naming_profiles_coherence.json").write_text(json.dumps(checks, indent=2))
        log(f"coherence: {json.dumps(checks)}")


# ----------------------------------------------------------------------------------
# determinism: proved by running, never asserted
# ----------------------------------------------------------------------------------

def determinism_check(args) -> None:
    """Cluster one unit twice, in two fresh interpreters, and compare byte for byte.

    The two interpreters are given **different** `PYTHONHASHSEED` values on purpose.  That is
    the exact variable the previous `SEED + hash(donor)` seeding was silently riding on: with
    it unset, CPython salts `str` hashing per process, so two runs of the same documented
    command produced different partitions.  If the assignments are identical under two
    different salts, the seeding no longer depends on the salt.

    A pipeline whose clusters move between runs cannot carry a cell-type call, so this is run
    rather than claimed, and its result is written next to the clusters.
    """
    import subprocess
    import tempfile

    out_dir = OUT_ROOT / args.modality
    out_dir.mkdir(parents=True, exist_ok=True)
    units = list(args.donors) if args.donors else (
        ["6461"] if args.modality == "phenocycler" else ["6347"]
    )
    id_col = "object_id" if args.modality == "phenocycler" else "cell_id"

    import os
    runs = []
    with tempfile.TemporaryDirectory(prefix="pln-determinism-") as tmp:
        for tag, hashseed in (("a", "0"), ("b", "12345")):
            root = Path(tmp) / tag
            cmd = [
                sys.executable, str(Path(__file__).resolve()),
                "--modality", args.modality, "--stage", "cluster",
                "--out-root", str(root), "--seed", str(args.seed),
                "--subsample", str(args.subsample), "--resolution", str(args.resolution),
                "--min-counts", str(args.min_counts), "--donors", *units,
            ]
            env = dict(os.environ, PYTHONHASHSEED=hashseed)
            log(f"determinism run {tag}: PYTHONHASHSEED={hashseed}")
            subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL)
            runs.append((tag, hashseed, root))

        files = {}
        for tag, hashseed, root in runs:
            for path in sorted((root / args.modality / "assignments").rglob("clusters.parquet")):
                rel = str(path.relative_to(root / args.modality))
                files.setdefault(rel, {})[tag] = hashlib.sha256(path.read_bytes()).hexdigest()

        per_file = []
        for rel, digests in sorted(files.items()):
            entry = {
                "assignments_file": rel,
                "sha256_run_a_PYTHONHASHSEED_0": digests.get("a"),
                "sha256_run_b_PYTHONHASHSEED_12345": digests.get("b"),
                "identical": digests.get("a") == digests.get("b") and digests.get("a") is not None,
            }
            # and does the tracked build reproduce from the same command?
            tracked = OUT_ROOT / args.modality / rel
            fresh = runs[0][2] / args.modality / rel
            if tracked.exists():
                t = pd.read_parquet(tracked, columns=[id_col, "cluster_id"])
                f = pd.read_parquet(fresh, columns=[id_col, "cluster_id"])
                merged = t.merge(f, on=id_col, suffixes=("_tracked", "_fresh"))
                entry["tracked_rows"] = int(len(t))
                entry["joined_on_" + id_col] = int(len(merged))
                entry["tracked_build_reproduced"] = bool(
                    len(merged) == len(t) == len(f)
                    and (merged["cluster_id_tracked"] == merged["cluster_id_fresh"]).all()
                )
            per_file.append(entry)

    record = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "modality": args.modality,
        "units_checked": units,
        "seed_base": int(args.seed),
        "method": (
            "the same --stage cluster command run in two fresh interpreters with different "
            "PYTHONHASHSEED values, assignments compared by sha256 of the parquet bytes"
        ),
        "all_identical": all(e["identical"] for e in per_file) and bool(per_file),
        "tracked_build_reproduced": all(
            e.get("tracked_build_reproduced", False) for e in per_file
        ) and bool(per_file),
        "per_file": per_file,
    }
    cohort_out(out_dir, "determinism_check.json").write_text(json.dumps(record, indent=2))
    log(f"wrote {cohort_out(out_dir, 'determinism_check.json')}: identical={record['all_identical']} "
        f"tracked_reproduced={record['tracked_build_reproduced']}")


# ----------------------------------------------------------------------------------
# published-eta: the batch term for every cross-donor quantity either lane PUBLISHES
# ----------------------------------------------------------------------------------
#: The per-donor tables the two pLN immune lanes publish.  This stage reads them; it does not
#: write them and does not own them.
PHENO_PUBLISHED = REPO / "data/development/typing/pln/pheno/per_donor_counts.csv"
XEN_PUBLISHED = REPO / ("data/development/xenium_typing/release_pln_immune_v2/"
                        "per_donor_summary.csv")
XEN_PATHS = REPO / "core/phenocycler/data/integration/xenium_paths.csv"


def xenium_batch_by_donor(paths_csv: Path = XEN_PATHS) -> pd.Series:
    """The Xenium instrument batch per donor, from the pLN rows of xenium_paths.csv.

    Only `roi` beginning `pln` is read: a donor's pancreas section can sit in a different batch
    from its node, and the quantities here are node quantities.  A ``PhenoXeno`` row is never
    read: it is a same-slide acquisition (Xenium on the PhenoCycler slide), not one of the
    cohort bundles whose nodes this lane clusters -- xenium_typing/pins.py keeps them out for the
    same reason.  A donor whose pLN sections disagree on batch has no single batch, and gets NaN
    rather than one of the two -- it is dropped from that grouping, never pooled into a level
    invented for it.
    """
    frame = pd.read_csv(paths_csv, dtype=str)
    pln = frame[frame["roi"].fillna("").str.startswith("pln")
                & (frame["batch"].fillna("").str.strip() != "PhenoXeno")]
    out = {}
    for donor, grp in pln.groupby("donor_id"):
        batches = sorted(set(grp["batch"].dropna()))
        out[donor] = batches[0] if len(batches) == 1 else float("nan")
    return pd.Series(out, dtype=object)


def published_quantities(modality: str, path: Path) -> tuple[pd.DataFrame, dict]:
    """Every per-donor fraction the lane publishes, plus the coherence check on the totals.

    Each quantity is a FRACTION of a stated denominator, because eta-squared on a count would
    partly measure how many cells each node happened to yield.  The denominators are the ones
    the lane's own table is built on, and the subsets are checked to sum to their compartment
    before anything is computed on them -- internal coherence, which is a legitimate check;
    no expected value and no ordering is asserted anywhere here.
    """
    frame = pd.read_csv(path, dtype={"donor_id": str})
    frame["donor_id"] = frame["donor_id"].astype(str)
    coherence: dict = {"source": str(path), "n_donors": int(len(frame))}
    rows = {}
    if modality == "phenocycler":
        classes = [c for c in frame.columns if c.startswith("class_")]
        parts = frame[["n_immune", "n_non_immune", "n_uncalled"]].sum(axis=1)
        coherence["immune_plus_non_immune_plus_uncalled_equals_eligible"] = bool(
            (parts == frame["n_eligible"]).all())
        coherence["class_columns_sum_to_immune"] = bool(
            (frame[classes].sum(axis=1) == frame["n_immune"]).all())
        coherence["pct_immune_column_matches_the_ratio"] = bool(
            np.allclose(frame["pct_immune"], 100.0 * frame["n_immune"] / frame["n_eligible"]))
        den = frame["n_eligible"].to_numpy(dtype=float)
        rows["immune_fraction_of_eligible"] = frame["n_immune"].to_numpy(dtype=float) / den
        rows["non_immune_fraction_of_eligible"] = frame["n_non_immune"].to_numpy(dtype=float) / den
        rows["uncalled_fraction_of_eligible"] = frame["n_uncalled"].to_numpy(dtype=float) / den
        rows["eligible_fraction_of_rows"] = (frame["n_eligible"].to_numpy(dtype=float)
                                             / frame["n_rows"].to_numpy(dtype=float))
        imm = frame["n_immune"].to_numpy(dtype=float)
        for col in classes:
            v = frame[col].to_numpy(dtype=float)
            rows[f"{col}_fraction_of_eligible"] = v / den
            with np.errstate(invalid="ignore", divide="ignore"):
                rows[f"{col}_fraction_of_immune"] = np.where(imm > 0, v / imm, np.nan)
    else:
        subsets = ["immune", "non_immune", "uncalled_ambiguous", "not_clustered",
                   "no_nucleus_excluded"]
        coherence["subsets_sum_to_cells"] = bool(
            (frame[subsets].sum(axis=1) == frame["cells"]).all())
        s2 = [c for c in frame.columns if c.startswith("s2__")]
        coherence["s2_columns_sum_to_cells"] = bool(
            (frame[s2].sum(axis=1) == frame["cells"]).all())
        den = frame["cells"].to_numpy(dtype=float)
        for col in subsets + s2:
            rows[f"{col}_fraction_of_cells"] = frame[col].to_numpy(dtype=float) / den
        imm = frame["immune"].to_numpy(dtype=float)
        for col in s2:
            if col in ("s2__Non_immune", "s2__Uncalled_ambiguous", "s2__Uncalled_not_clustered",
                       "s2__Unavailable_no_nucleus"):
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                rows[f"{col}_fraction_of_immune"] = np.where(
                    imm > 0, frame[col].to_numpy(dtype=float) / imm, np.nan)
    quantities = pd.DataFrame(rows, index=frame["donor_id"].tolist())
    quantities.index.name = "donor_id"
    return quantities, coherence


def published_eta(modality: str, pheno_path: Path, xen_path: Path, n_perm: int,
                  required: bool = True) -> None:
    """eta-squared, with its permutation null, for every cross-donor quantity this lane publishes.

    WHY THIS EXISTS.  Both lanes report their immune quantities per donor, and the text carried
    the per-donor numbers with no batch term beside them.  A reader then has no way to tell
    whether a spread across donors is a spread across acquisitions.  This stage does not remove
    a confounded quantity and does not gate anything on the answer: the batch term is REPORTED,
    beside the quantity, with the null that says what that number would be on a grouping of the
    same shape carrying no batch structure at all.

    The grouping variables are donor-level records of the PANCREAS scan (`data/
    donor_acquisition_runs.csv`, `data/donor_acquisition_batch.csv`) except `xenium_batch`,
    which is the instrument batch of the pLN section itself.  `covariate_provenance.csv` says,
    per donor, whether any tracked image record covers the node; where it does not, the
    PhenoCycler grouping is a donor-level proxy and is labelled as one in the output.
    """
    out_dir = OUT_ROOT / modality
    out_dir.mkdir(parents=True, exist_ok=True)
    path = pheno_path if modality == "phenocycler" else xen_path
    if not path.exists():
        message = (f"{path} does not exist -- this stage reports the batch term for the table "
                   "the lane publishes, and cannot invent one")
        if required:
            raise SystemExit(f"[err] {message}")
        log(f"[warn] skipping published-eta: {message}")
        return
    quantities, coherence = published_quantities(modality, path)
    # A partial run reports the batch term for the donors it names, and lands under
    # `_partial/donors=<slug>/`.  The coherence check above ran on the whole published table,
    # which is what it is a check of; the eta-squared below is of this run's donors only, and
    # the count is in the record.
    scope = set_scope(SCOPE.refine(quantities.index.tolist()))
    scope.announce(out_dir, log)
    if scope.is_partial:
        quantities = quantities.loc[[d for d in quantities.index if d in set(scope.requested)]]
        write_scope_note(scope, out_dir, {"stage": "published-eta", "modality": modality,
                                          "source_table": str(path)})
    donors = quantities.index.tolist()
    if len(donors) < 3:
        raise SystemExit(f"[err] {len(donors)} donor(s) after scoping -- eta-squared over "
                         "fewer than three donors measures nothing; name more, or drop --donors")

    cov = donor_covariates().reindex(donors)
    groupings = ["stain_group", "acq_run", "coarse_batch", "software", "disease_status"]
    if modality == "xenium":
        cov["xenium_batch"] = xenium_batch_by_donor().reindex(donors)
        groupings.append("xenium_batch")

    rows = []
    for var in groupings:
        have = cov[var].notna().to_numpy()
        for name in quantities.columns:
            v = quantities[name].to_numpy(dtype=float)
            ok = have & np.isfinite(v)
            if ok.sum() < 3 or len(set(cov.loc[ok, var].to_numpy().tolist())) < 2:
                rows.append({
                    "modality": modality, "quantity": name, "grouping": var,
                    "n_donors_used": int(ok.sum()),
                    "n_donors_unrecorded": int((~have).sum()),
                    "n_levels": int(len(set(cov.loc[ok, var].to_numpy().tolist()))),
                    "eta_squared": float("nan"), "null_median": float("nan"),
                    "p_permutation": float("nan"), "n_perm": 0,
                    "state": "unavailable(too few donors or one level)",
                })
                continue
            result = eta_squared_with_null(v[ok], cov.loc[ok, var].to_numpy(), n_perm=n_perm)
            # A quantity that is the same on every donor has no variance to partition, so
            # eta-squared is undefined rather than zero.  `class_NK_ILC` is exactly that on
            # the PhenoCycler table read here -- 0 on all 14 donors -- and the state says so
            # rather than letting a NaN read as a measured absence of batch structure.
            constant = float(np.nanstd(v[ok])) == 0.0
            rows.append({
                "modality": modality, "quantity": name, "grouping": var,
                "n_donors_used": int(ok.sum()),
                "n_donors_unrecorded": int((~have).sum()),
                "n_levels": int(len(set(cov.loc[ok, var].to_numpy().tolist()))),
                **result,
                "state": ("unavailable(quantity_is_constant_across_donors)" if constant
                          else "measured" if np.isfinite(result["eta_squared"])
                          else "unavailable(eta_squared_undefined)"),
            })
    eta = pd.DataFrame(rows)
    dest = cohort_out(out_dir, "published_quantity_batch_eta_squared.csv")
    eta.to_csv(dest, index=False)

    record = {
        "what": ("eta-squared, and the eta-squared a random grouping of the same shape gives, "
                 "for every cross-donor quantity this lane publishes"),
        "source_table": str(path),
        "source_mtime_utc": datetime.fromtimestamp(path.stat().st_mtime,
                                                   tz=timezone.utc).isoformat(),
        "n_donors": len(donors),
        "donors": donors,
        "scope": scope.describe(),
        "n_quantities": int(len(quantities.columns)),
        "quantities": list(quantities.columns),
        "groupings": groupings,
        "n_perm": n_perm,
        "coherence_of_the_published_totals": coherence,
        "covariate_provenance": covariate_provenance(donors).to_dict("records"),
        "reading": ("eta_squared is not comparable across groupings with different level "
                    "counts. Read it against null_median, which is the median of n_perm "
                    "relabellings with identical level sizes, and p_permutation, the fraction "
                    "of those relabellings reaching the observed value. A high eta-squared "
                    "with a high p is a number the grouping's degrees of freedom can produce "
                    "on their own."),
        "not_a_gate": ("nothing is dropped or corrected on these numbers here. They are "
                       "published beside the quantity so a reader can see the batch term."),
    }
    cohort_out(out_dir, "published_quantity_batch_eta_squared.json").write_text(
        json.dumps(record, indent=2, default=str) + "\n")
    log(f"wrote {dest} ({len(eta)} rows: {len(quantities.columns)} quantities x "
        f"{len(groupings)} groupings)")
    top = eta[eta["state"] == "measured"].sort_values("eta_squared", ascending=False)
    for _, r in top.head(8).iterrows():
        log(f"  {r['quantity']} ~ {r['grouping']}: eta2={r['eta_squared']:.3f} "
            f"null_median={r['null_median']:.3f} p_perm={r['p_permutation']:.3f} "
            f"(n={r['n_donors_used']}, levels={r['n_levels']})")


def main() -> None:
    global OUT_ROOT
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modality", choices=("phenocycler", "xenium"), required=True)
    ap.add_argument("--stage",
                    choices=("cluster", "profile", "determinism", "summarize", "link",
                             "prose", "published-eta", "all"),
                    default="all")
    ap.add_argument("--pheno-published", default=str(PHENO_PUBLISHED),
                    help="per_donor_counts.csv of the PhenoCycler pLN immune lane, read by "
                         "--stage published-eta")
    ap.add_argument("--xen-published", default=str(XEN_PUBLISHED),
                    help="per_donor_summary.csv of the Xenium pLN immune release, read by "
                         "--stage published-eta")
    ap.add_argument("--n-perm", type=int, default=2000,
                    help="permutations behind every eta-squared null")
    ap.add_argument("--out-root", default=None,
                    help="write under this root instead of the tracked output root. Exists so "
                         "the same command can be run twice from two fresh interpreters and "
                         "the assignments compared; it never changes what is computed.")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--subsample", type=int, default=60_000,
                    help="cells used to build the neighbour graph, per donor/section")
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--meta-k", type=int, default=14,
                    help="how many groups the cross-donor profile matching is cut into")
    ap.add_argument("--min-counts", type=int, default=10,
                    help="xenium only: minimum transcripts per cell")
    ap.add_argument("--seed", type=int, default=SEED,
                    help="the fixed seed every per-donor seed is derived from. Exists so the "
                         "partition's sensitivity to the seed can be MEASURED -- which is what "
                         "the old `SEED + hash(donor)` was varying silently. Changing it "
                         "changes the output; the tracked build uses the default.")
    args = ap.parse_args()

    sys.path.insert(0, str(REPO / "core/multiplex5000"))
    if args.out_root:
        OUT_ROOT = Path(args.out_root).resolve()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # A named selection is PARTIAL until a stage that enumerates its units proves otherwise, and
    # a partial run's cohort-shaped outputs go to `_partial/donors=<slug>/`.  `summarize`,
    # `link` and `prose` never enumerate -- they read whatever is on disk -- so under
    # `--donors` they stay scoped.  `--out-root` does not make a run whole: a private root is
    # still a root, and a partial file copied off it must not arrive wearing a cohort name.
    set_scope(RunScope.undiscovered("donor", args.donors))
    if SCOPE.is_partial:
        log(f"[scope] --donors given: this is a partial run. Cohort-wide artifacts under "
            f"{OUT_ROOT / args.modality} will NOT be rewritten; this run's cohort-shaped "
            f"outputs go to {SCOPE.partial_dir(OUT_ROOT / args.modality)}")

    if args.stage in ("cluster", "all"):
        if args.modality == "phenocycler":
            run_pheno(args)
        else:
            run_xen(args)
    if args.stage == "determinism":
        determinism_check(args)
        return
    if args.stage in ("profile", "all"):
        if args.modality == "phenocycler":
            naming_profile_pheno(args)
        else:
            naming_profile_xen(args)
    if args.stage in ("summarize", "all"):
        summarize(args.modality, args.meta_k)
    if args.stage in ("link", "all"):
        link_immune_calls(args.modality)
    if args.stage in ("prose", "all"):
        prose_facts(args.modality)
    if args.stage in ("published-eta", "all"):
        # `--stage all` is the whole clustering lane; the published immune table belongs to the
        # OTHER lane and may not exist yet.  Asked for by name it is an error; reached as part
        # of `all` it is skipped with a warning rather than failing a clustering run.
        published_eta(args.modality, Path(args.pheno_published), Path(args.xen_published),
                      args.n_perm, required=args.stage == "published-eta")


if __name__ == "__main__":
    main()
