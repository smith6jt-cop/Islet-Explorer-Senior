#!/usr/bin/env python3
"""Islet-level trajectory on composition -- the Trajectory tab's H5AD.

The ordering is over **what an islet is made of**, not over marker averages.  For islet *i* and
class *t*, in region *r*:

    alr[i, r, t] = log( (n[i, r, t] + c) / (n_core_endocrine[i] + c) )        c = 0.5

The denominator is the islet's own **core endocrine count** in every region, because the hormone
question is only asked inside the islet.  That makes the expansion axis read "immune cells around
this islet per endocrine cell in it" rather than a fraction of a peri-zone whose size varies.  An
islet with no core endocrine cell has no defined ratio: it is ``Unavailable`` with its reason
recorded, and the count is reported rather than hidden.

**Two things this script deliberately does not do.**

*It does not gate on a biological result.*  The legacy ``compute_trajectory_senior.py`` requires
pseudotime to anti-correlate with INS, to correlate with GCG, and to put T1D islets later than ND
ones, and rejects the build otherwise.  A check that passes by reproducing a previous finding
measures the expectation, not the build.  Those three quantities are computed here and **reported**;
what gates is the method: a finite pseudotime for every included islet, a root that is a real ND
islet above the cell floor, and an ordering that survives re-rooting and a change of neighbourhood
size.  The donor and disease eta-squared are reported for the same reason -- they say how much of
the axis is donor rather than islet, which is a property of the build worth knowing either way.

*It does not batch-correct.*  Correction needs a batch variable carrying within-group disease
contrast.  Whether one exists is a property of the cohort in front of the script, so it is
**measured per run** and written to ``uns['batch_correction_reason']``: the staining groups that
hold a single donor within a disease status are named, because correcting on ``stain_group`` would
take those donors' disease contributions out with their group effects.  Harmony is not applied;
the cohort is stratified and the confound reported.

**The cohort limit is likewise derived per run**, not written here, and recorded in
``uns['cohort_limit']`` so the app can surface it: the donors per status, and the share of the
pseudotime donor explains within each.  This ordering is an islet-composition axis.  It is not a
disease-progression axis, and nothing downstream should read it as one.

Those three strings used to be fixed text describing the ten-donor build.  Run on twenty donors
they kept asserting "eight ND donors, two T1D, no Aab+" beside a ``within_status_donor_count`` of
``{ND 8, Aab+ 3, T1D 9}`` in the same dict.  A cohort fact written as a literal is a leak whether
it sits in the app or in the producer that feeds it.

    python core/phenocycler/scripts/typing/trajectory_islet.py
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"

ALR_PSEUDOCOUNT = 0.5          #: declared, recorded in uns -- never silently tuned
#: Euclidean, not cosine. An ALR vector's entries are log-ratios against a denominator that is
#: almost always larger, so every islet's vector points into the same orthant and cosine put the
#: median pair at 0.97 similarity -- the neighbourhood graph was then near-arbitrary and the
#: ordering did not survive a change of neighbourhood size (rho 0.31). Euclidean distance on
#: log-ratios is the Aitchison geometry the representation already lives in.
METRIC = "euclidean"
N_NEIGHBORS = 30
N_NEIGHBORS_ALT = 15           #: the re-run the stability check compares against
K_ROOT = 50                    #: the root is a medoid over this many islets, not one argmax
ROOT_MIN_CELLS = 20
#: Minimum core cells for an islet to receive a pseudotime.  Chosen by where the ordering becomes
#: reproducible, not by any result.  Measured, recorded in uns, and reported at run time:
#:
#:     floor   islets  rho(n_neighbors 30 vs 15)     [measured on the ten-donor build, 2026-09-14]
#:         0    7,302   -0.537   <- the two runs order the islets OPPOSITELY
#:         5    5,707    0.632
#:        10    4,371    0.931   <- chosen
#:        20    3,101    0.958
#:        50    1,687    0.958
#:
#: That table is the evidence for the VALUE and stays here as such; it is not carried into `uns`,
#: where each run counts its own islets per floor.  A composition over eight classes is estimated
#: from very few cells -- the ten-donor median islet held 14 -- and below ten the ordering is not
#: identified at all, while above ten the stability plateaus.  Islets below the floor stay in the file -- the app keys per islet and
#: pairs core with peri -- but their pseudotime is Unavailable with a reason, never zero.
COMPOSITION_FLOOR_CELLS = 10
ENDOCRINE = ("Beta_cell", "Alpha_cell", "Delta_cell")
#: Below this, the pseudotime is not reproducible enough to call an ordering.
STABILITY_MIN_RHO = 0.90
#: Leiden resolutions for the Spatial tab's islet-level UMAP, matching the legacy column names
#: (`leiden_0.3` ... `leiden_1.0`) that `mod_spatial_server.R` greps for.
LEIDEN_RESOLUTIONS = (0.3, 0.5, 0.8, 1.0)
#: Ceiling on the share of pseudotime variance explained by donor *within* a disease status.
#: The legacy lane used the same value for the same reason: above it, the axis is a donor axis.
DONOR_ETA2_CEILING = 0.15
#: ...but eta-squared over two donors is one pairwise contrast, not a variance decomposition: it
#: cannot separate "the build has a donor artefact" from "these two people differ".  A status
#: needs at least three donors before the number is a check rather than a description.  On this
#: cohort that gates ND (eight donors) and reports T1D (two) -- which is the cohort's limit
#: showing through, not an exemption granted to a failing number.
DONOR_ETA2_MIN_DONORS = 3
#: The coarse acquisition batch, derived from the donors whose QuPath export was re-done rather
#: than restated: `cohort.BATCH2_REEXPORTED` IS the second batch. Batch is aliased with disease
#: status -- every ND donor is in batch 1 -- so it is reported and never corrected on.
BATCH2_DONORS = frozenset(COHORT.BATCH2_REEXPORTED)


def load_islets(donors: list[str], proportions: pathlib.Path | None = None) -> pd.DataFrame:
    frame = pd.read_csv(proportions or (APP / "islet_celltype_proportions.csv"))
    frame = frame[frame["case_id"].astype(str).isin(donors)].reset_index(drop=True)
    meta = pd.read_csv(APP / "donor_metadata.csv")[
        ["case_id", "donor_status", "age", "sex"]].rename(columns={"sex": "gender"})
    runs = pd.read_csv(S.REPO / "data/donor_acquisition_runs.csv")
    runs = runs.rename(columns={"donor_id": "case_id"})[["case_id", "stain_group", "acq_run"]]
    runs = runs.drop_duplicates("case_id")
    return frame.merge(meta, on="case_id", how="left").merge(runs, on="case_id", how="left")


def alr_matrix(frame: pd.DataFrame, classes: list[str], prefix: str,
               denominator: np.ndarray) -> np.ndarray:
    counts = np.column_stack([frame[f"{prefix}count_{c}"].to_numpy(float) for c in classes])
    return np.log((counts + ALR_PSEUDOCOUNT)
                  / (denominator + ALR_PSEUDOCOUNT)[:, None]).astype(np.float32)


def diffusion_pseudotime(adata: ad.AnnData, rep: str, root: int,
                         n_neighbors: int) -> np.ndarray:
    work = adata.copy()
    for key in ("X_diffmap",):
        work.obsm.pop(key, None)
    for key in ("diffmap_evals", "neighbors", "iroot"):
        work.uns.pop(key, None)
    sc.pp.neighbors(work, n_neighbors=n_neighbors, use_rep=rep, metric=METRIC)
    work.uns["iroot"] = root
    sc.tl.diffmap(work, n_comps=10)
    sc.tl.dpt(work)
    return work.obs["dpt_pseudotime"].to_numpy().copy(), work


def embedding_eta2(coords: np.ndarray, donor: np.ndarray,
                   donor_to_group: dict[str, str], permutations: int = 500,
                   seed: int = 0) -> dict[str, float]:
    """How much of a 2-D embedding is one grouping of DONORS, against a null that keeps nesting.

    eta-squared here is the between-group sum of squares over the total, summed over both UMAP
    axes.  The null is what makes the number readable: islets are nested in donors, so shuffling
    islet labels asks a question nobody has -- it would call any donor-level grouping significant.
    The donor-to-group map is permuted instead, group sizes preserved, so the null says "how much
    would an ARBITRARY partition of these donors explain".  A grouping that is one donor per group
    is therefore exactly its own null (p = 1) rather than a finding.
    """
    rng = np.random.default_rng(seed)

    def eta2_2d(groups: np.ndarray) -> float:
        ok = np.isfinite(coords).all(axis=1)
        points, groups = coords[ok], groups[ok]
        grand = points.mean(axis=0)
        between = sum((groups == g).sum() * ((points[groups == g].mean(axis=0) - grand) ** 2).sum()
                      for g in np.unique(groups))
        total = ((points - grand) ** 2).sum()
        return float(between / total) if total > 0 else float("nan")

    donors = np.array(sorted(donor_to_group))
    labels = np.array([donor_to_group[d] for d in donors])
    observed = eta2_2d(np.array([donor_to_group[d] for d in donor]))
    # One permutation PER DRAW, hoisted out of the per-islet lookup. Built inline it was
    # re-drawn for every islet, which silently turns the donor-block null back into the
    # islet-level shuffle it exists to avoid -- and that null is ~0 for every grouping, so the
    # p-value looks decisive for a grouping that is its own null.
    null = np.empty(permutations)
    for draw in range(permutations):
        shuffled = dict(zip(donors, rng.permutation(labels)))
        null[draw] = eta2_2d(np.array([shuffled[d] for d in donor]))
    return {"eta2": observed, "null_mean": float(null.mean()),
            "null_p95": float(np.percentile(null, 95)), "k": int(len(set(labels))),
            "p": float((1 + (null >= observed).sum()) / (1 + len(null)))}


def pick_root(adata: ad.AnnData, rep: str) -> tuple[int, dict]:
    """Medoid of the top-K ND islets by endocrine fraction, above the cell floor.

    A single argmax once selected a saturated three-cell islet, which is why this is a medoid
    over a pool rather than the extreme point.
    """
    status = adata.obs["donor_status"].astype(str).to_numpy()
    cells = adata.obs["total_cells"].to_numpy()
    endocrine_fraction = adata.obs["endocrine_fraction"].to_numpy()
    pool = np.where((status == "ND") & (cells >= ROOT_MIN_CELLS))[0]
    if pool.size < 10:
        pool = np.where(status == "ND")[0]
    if pool.size == 0:
        raise ValueError("no ND islet available to root the ordering")
    top = pool[np.argsort(-endocrine_fraction[pool])[:min(K_ROOT, pool.size)]]
    representation = np.asarray(adata.obsm[rep])
    centroid = representation[top].mean(axis=0)
    root = int(top[np.argmin(np.linalg.norm(representation[top] - centroid, axis=1))])
    donors = adata.obs["imageid"].iloc[top].astype(str).value_counts().to_dict()
    return root, {"pool_size": int(top.size), "pool_donors": donors,
                  "pool": [int(i) for i in top]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    #: The per-islet proportions to fit on, and where to write.  Both default to the app's domain
    #: root, which is what a deploy wants; both are settable so a rebuild can be fitted and
    #: inspected without overwriting the deployed files while a question about it is still open.
    parser.add_argument("--islets", type=pathlib.Path, default=None,
                        help="islet_celltype_proportions.csv to fit on (default: the app's)")
    parser.add_argument("--out-dir", type=pathlib.Path, default=None,
                        help="where islet_leiden.parquet and the H5AD are written (default: app)")
    parser.add_argument("--leiden-name", default="islet_leiden.parquet",
                        help="basename for the islet-level clustering, so a rebuild can sit "
                             "beside the deployed file instead of replacing it")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    out_dir = args.out_dir or APP
    out_dir.mkdir(parents=True, exist_ok=True)

    frame = load_islets(donors, args.islets)
    classes = sorted(c[len("count_"):] for c in frame.columns if c.startswith("count_"))
    print(f"{len(frame):,} islets, {len(donors)} donors, {len(classes)} classes\n")

    endocrine_core = np.sum(
        [frame[f"count_{c}"].to_numpy(float) for c in ENDOCRINE if f"count_{c}" in frame],
        axis=0)
    defined = endocrine_core > 0
    print(f"islets with a core endocrine cell: {int(defined.sum()):,} of {len(frame):,}")
    print(f"  Unavailable (no core endocrine cell, ratio undefined): "
          f"{int((~defined).sum()):,}")
    undefined_by_donor = (frame.loc[~defined, "case_id"].astype(str)
                          .value_counts().sort_index().to_dict())
    if undefined_by_donor:
        print(f"  by donor: {undefined_by_donor}")

    frame = frame[defined].reset_index(drop=True)
    endocrine_core = endocrine_core[defined]

    core = alr_matrix(frame, classes, "", endocrine_core)
    peri = alr_matrix(frame, classes, "peri_", endocrine_core)
    counts = np.column_stack([frame[f"count_{c}"].to_numpy(float) for c in classes])
    # X is what the Trajectory tab plots against pseudotime, under the axis label
    # "normalized expression (mean per islet, 0-1)". Counts there run to ~1,500 and make that
    # label a lie; the per-islet proportion is both what the label says and the right per-islet
    # feature for a composition ordering. The counts stay available in obs as total_cells and in
    # the ALR representations.
    core_total = frame["total_cells_core"].to_numpy(float)
    proportions = np.divide(counts, core_total[:, None],
                            out=np.full_like(counts, np.nan),
                            where=core_total[:, None] > 0)

    obs = pd.DataFrame({
        "imageid": frame["case_id"].astype(int).astype(str).to_numpy(),
        "base_islet_id": frame["islet_key"].to_numpy(),
        "combined_islet_id": frame["combined_islet_id"].to_numpy(),
        "donor_status": frame["donor_status"].to_numpy(),
        "age": frame["age"].to_numpy(),
        "gender": frame["gender"].to_numpy(),
        "stain_group": frame["stain_group"].to_numpy(),
        "acq_run": frame["acq_run"].to_numpy(),
        "total_cells": frame["total_cells_core"].to_numpy(),
        "total_cells_peri": frame["total_cells_peri"].to_numpy(),
        "endocrine_core": endocrine_core,
        "endocrine_fraction": np.divide(endocrine_core, frame["total_cells_core"].to_numpy(float),
                                        out=np.zeros(len(frame)),
                                        where=frame["total_cells_core"].to_numpy(float) > 0),
    })
    obs.index = obs["combined_islet_id"].astype(str)

    # var needs a real column AND an unnamed index, matching the CODEX file. An index-only
    # DataFrame, or one whose index name duplicates a column name, round-trips through the app's
    # reticulate reader as a bare integer index, and the Trajectory feature selector then offers
    # "1, 2, 3 ..." in place of the class names.
    adata = ad.AnnData(X=proportions.astype(np.float32), obs=obs,
                       var=pd.DataFrame({"feature_type": ["composition"] * len(classes)},
                                        index=pd.Index(classes)))
    adata.obsm["X_composition_alr"] = core
    adata.obsm["X_composition_alr_combined"] = np.hstack([core, peri]).astype(np.float32)

    # The graph is built on the islets above the floor only; the rest carry an Unavailable
    # pseudotime.  Including them poisons the graph for everyone, which is what the floor-0 row
    # of the sweep above measures.
    eligible = adata.obs["total_cells"].to_numpy() >= COMPOSITION_FLOOR_CELLS
    adata.obs["pseudotime_unavailable_reason"] = np.where(
        eligible, "", f"below_composition_floor_{COMPOSITION_FLOOR_CELLS}_core_cells")
    print(f"\nislets above the {COMPOSITION_FLOOR_CELLS}-cell composition floor: "
          f"{int(eligible.sum()):,} of {adata.n_obs:,} ({100 * eligible.mean():.1f}%)")
    ordered = adata[eligible].copy()

    root, root_meta = pick_root(ordered, "X_composition_alr")
    print(f"root: {ordered.obs_names[root]} "
          f"(donor {ordered.obs['imageid'].iloc[root]}, "
          f"{int(ordered.obs['total_cells'].iloc[root])} core cells, "
          f"endocrine fraction {ordered.obs['endocrine_fraction'].iloc[root]:.3f})")

    pseudotime, fitted = diffusion_pseudotime(ordered, "X_composition_alr", root, N_NEIGHBORS)
    combined, _ = diffusion_pseudotime(ordered, "X_composition_alr_combined", root, N_NEIGHBORS)

    # --- Method checks. None of these is a biological result. -------------------------------
    finite = np.isfinite(pseudotime)
    alt_neighbors, _ = diffusion_pseudotime(ordered, "X_composition_alr", root, N_NEIGHBORS_ALT)
    both = finite & np.isfinite(alt_neighbors)
    rho_neighbors = float(spearmanr(pseudotime[both], alt_neighbors[both])[0])

    # Root stability over the WHOLE pool the root is defined from, not one point in it. The
    # question the check asks is "does it matter which islet of the defined pool we start at",
    # and one adversarially-chosen member is a poor estimate of that. Median gates; the minimum
    # is reported, so a stable-in-the-median / unstable-at-the-edge ordering is visible.
    rhos = []
    for candidate in root_meta["pool"]:
        if candidate == root:
            continue
        alt, _ = diffusion_pseudotime(ordered, "X_composition_alr", candidate, N_NEIGHBORS)
        ok = finite & np.isfinite(alt)
        if ok.sum() > 10:
            rhos.append(float(spearmanr(pseudotime[ok], alt[ok])[0]))
    rho_root = float(np.median(rhos)) if rhos else float("nan")
    rho_root_min = float(np.min(rhos)) if rhos else float("nan")

    status = ordered.obs["donor_status"].astype(str).to_numpy()
    donor = ordered.obs["imageid"].astype(str).to_numpy()

    def eta2(values: np.ndarray, groups: np.ndarray) -> float:
        ok = np.isfinite(values)
        values, groups = values[ok], groups[ok]
        grand = values.mean()
        between = sum(((groups == g).sum()) * (values[groups == g].mean() - grand) ** 2
                      for g in np.unique(groups))
        total = ((values - grand) ** 2).sum()
        return float(between / total) if total > 0 else float("nan")

    # Within each disease status, how much of the pseudotime is explained by which donor an
    # islet came from? This is a property of the build, not a biological expectation: a high
    # value means the axis is a donor axis wearing a disease label. Statuses with one donor
    # cannot be measured and are skipped rather than assumed fine.
    within_status_donor_eta2, within_status_donors = {}, {}
    for s_ in np.unique(status):
        mask = status == s_
        n_donors = len(np.unique(donor[mask]))
        within_status_donors[str(s_)] = n_donors
        within_status_donor_eta2[str(s_)] = (eta2(pseudotime[mask], donor[mask])
                                             if n_donors >= 2 else float("nan"))
    gated_statuses = [s_ for s_, n in within_status_donors.items()
                      if n >= DONOR_ETA2_MIN_DONORS
                      and np.isfinite(within_status_donor_eta2[s_])]

    gates = {
        "all_islets_finite": bool(finite.all()),
        "root_is_ND": bool(status[root] == "ND"),
        "root_above_cell_floor": bool(ordered.obs["total_cells"].iloc[root] >= ROOT_MIN_CELLS),
        "stable_under_neighbourhood_size": bool(rho_neighbors >= STABILITY_MIN_RHO),
        "stable_under_reroot": bool(rho_root >= STABILITY_MIN_RHO),
        "pseudotime_is_not_mostly_donor": bool(all(
            within_status_donor_eta2[s_] <= DONOR_ETA2_CEILING for s_ in gated_statuses)),
    }

    # --- Reported, never gated. -------------------------------------------------------------
    by_status = {s: float(np.nanmedian(pseudotime[status == s]))
                 for s in np.unique(status)}
    ins_like = np.asarray(ordered[:, "Beta_cell"].X).ravel() if "Beta_cell" in classes else None
    reported = {
        "median_pseudotime_by_status": by_status,
        "eta2_donor": eta2(pseudotime, donor),
        "eta2_donor_status": eta2(pseudotime, status),
        "eta2_stain_group": eta2(pseudotime, ordered.obs["stain_group"].astype(str).to_numpy()),
        "spearman_pseudotime_beta_proportion": (
            float(spearmanr(pseudotime[finite], ins_like[finite])[0])
            if ins_like is not None else float("nan")),
        "rho_vs_smaller_neighbourhood": rho_neighbors,
        "rho_vs_alternate_root_median": rho_root,
        "rho_vs_alternate_root_min": rho_root_min,
        "within_status_donor_eta2": within_status_donor_eta2,
        "within_status_donor_count": within_status_donors,
        "within_status_donor_eta2_gated": gated_statuses,
        "within_status_donor_eta2_reported_only": [
            s_ for s_ in within_status_donors if s_ not in gated_statuses],
    }
    if {"ND", "T1D"} <= set(np.unique(status)):
        reported["mwu_T1D_gt_ND_p"] = float(mannwhitneyu(
            pseudotime[status == "T1D"], pseudotime[status == "ND"],
            alternative="greater").pvalue)

    print("\nmethod checks (these gate the build)")
    for name, passed in gates.items():
        print(f"  {name:<34s} {'PASS' if passed else 'FAIL'}")
    print("\nreported, not gated (a result is never a verification target)")
    for name, value in reported.items():
        print(f"  {name:<34s} {value}")

    # Map back onto the full islet set: below the floor, the pseudotime is Unavailable.
    into = pd.Series(np.nan, index=adata.obs_names)
    into.loc[ordered.obs_names] = pseudotime
    adata.obs["dpt_pseudotime"] = into.to_numpy()
    into.loc[:] = np.nan
    into.loc[ordered.obs_names] = combined
    adata.obs["dpt_pseudotime_combined"] = into.to_numpy()

    sc.tl.umap(fitted)
    # Islet-level Leiden on the same graph the ordering was built from -- the Spatial tab's
    # "Islet-Level UMAPs" card reads `leiden_<res>` and `leiden_umap_1/2` and had been showing
    # "Leiden clustering not available in current H5AD" because nothing produced them here.
    for resolution in LEIDEN_RESOLUTIONS:
        sc.tl.leiden(fitted, resolution=resolution, key_added=f"leiden_{resolution}",
                     flavor="igraph", n_iterations=2, directed=False)
        print(f"  leiden {resolution}: {fitted.obs[f'leiden_{resolution}'].nunique()} clusters")
    # --- How much of the EMBEDDING is batch? Reported, never gated. --------------------------
    # The pseudotime already carries a donor eta-squared; the UMAP is a different object and the
    # app draws it, so it gets its own number. Acquisition batch is aliased with disease status
    # at the coarse level -- every ND donor is in batch 1 -- so the coarse figure is printed and
    # marked unattributable, and the finer variables (`stain_group`, `acq_run`) are what a reader
    # should weigh. The within-status pass is where coarse batch becomes attributable at all: a
    # status containing both batches is the only place the two come apart.
    coords = np.asarray(fitted.obsm["X_umap"])
    fit_donor = fitted.obs["imageid"].astype(str).to_numpy()
    fit_status = fitted.obs["donor_status"].astype(str).to_numpy()
    one_per = {d: str(d) for d in np.unique(fit_donor)}
    donor_maps = {
        "coarse_batch": {d: ("batch2" if d in BATCH2_DONORS else "batch1")
                         for d in np.unique(fit_donor)},
        "stain_group": dict(zip(fit_donor, fitted.obs["stain_group"].astype(str))),
        "acq_run": dict(zip(fit_donor, fitted.obs["acq_run"].astype(str))),
        "donor_status": dict(zip(fit_donor, fit_status)),
        "donor": one_per,
    }
    embedding_batch = {name: embedding_eta2(coords, fit_donor, mapping)
                       for name, mapping in donor_maps.items()}
    print("\nembedding (2-D UMAP): eta-squared, against a null that permutes the donor->group map")
    for name, stats in embedding_batch.items():
        print(f"  {name:<14s} k={stats['k']:<3d} eta2={stats['eta2']:.3f}   "
              f"null mean={stats['null_mean']:.3f} p95={stats['null_p95']:.3f}  "
              f"p={stats['p']:.4f}")
    statuses_with_both = [s_ for s_ in np.unique(fit_status)
                          if len({("batch2" if d in BATCH2_DONORS else "batch1")
                                  for d in np.unique(fit_donor[fit_status == s_])}) > 1]
    print("  coarse_batch is ALIASED with donor status across the cohort and is not attributable "
          "there.")
    within_status_embedding = {}
    for s_ in np.unique(fit_status):
        mask = fit_status == s_
        here = np.unique(fit_donor[mask])
        entry = {"n_donors": int(len(here)),
                 "donor": embedding_eta2(coords[mask], fit_donor[mask],
                                         {d: str(d) for d in here})}
        if s_ in statuses_with_both:
            entry["coarse_batch"] = embedding_eta2(
                coords[mask], fit_donor[mask],
                {d: ("batch2" if d in BATCH2_DONORS else "batch1") for d in here})
        within_status_embedding[str(s_)] = entry
        batch_note = (f"coarse_batch eta2={entry['coarse_batch']['eta2']:.3f} "
                      f"(null {entry['coarse_batch']['null_mean']:.3f}, "
                      f"p={entry['coarse_batch']['p']:.4f})"
                      if "coarse_batch" in entry
                      else "coarse batch Unavailable here -- one batch only")
        print(f"  within {s_:<5s} ({entry['n_donors']} donors): "
              f"donor eta2={entry['donor']['eta2']:.3f}; {batch_note}")

    diffmap = pd.DataFrame(np.full((adata.n_obs, fitted.obsm["X_diffmap"].shape[1]), np.nan),
                           index=adata.obs_names)
    diffmap.loc[ordered.obs_names] = fitted.obsm["X_diffmap"]
    adata.obsm["X_diffmap"] = diffmap.to_numpy().astype(np.float32)
    umap = pd.DataFrame(np.full((adata.n_obs, 2), np.nan), index=adata.obs_names)
    umap.loc[ordered.obs_names] = fitted.obsm["X_umap"]
    adata.obsm["X_umap"] = umap.to_numpy().astype(np.float32)

    # The islet-level clustering, mapped back onto the full islet set and written beside the
    # H5AD for build_parquet_views.py to merge into islets.parquet, which is where the app reads
    # it from. Islets below the composition floor are Unavailable here too.
    leiden = pd.DataFrame(index=adata.obs_names)
    for resolution in LEIDEN_RESOLUTIONS:
        column = f"leiden_{resolution}"
        leiden[column] = pd.Series(fitted.obs[column].astype(str).to_numpy(),
                                   index=ordered.obs_names).reindex(adata.obs_names)
    leiden["leiden_umap_1"] = adata.obsm["X_umap"][:, 0]
    leiden["leiden_umap_2"] = adata.obsm["X_umap"][:, 1]
    leiden.insert(0, "islet_key", adata.obs["base_islet_id"].to_numpy())
    leiden.insert(0, "case_id", adata.obs["imageid"].astype(int).to_numpy())
    leiden.reset_index(drop=True).to_parquet(out_dir / args.leiden_name, index=False)
    for column in [c for c in leiden.columns if c.startswith("leiden_") and "umap" not in c]:
        print(f"  {column}: {leiden[column].notna().sum():,} islets clustered")

    # --- Three uns strings that were COHORT FACTS WRITTEN AS LITERALS -------------------------
    #
    # `cohort_limit`, `batch_correction_reason` and `composition_floor_sweep` were fixed text and
    # a fixed table describing the ten-donor build. Run on twenty donors they went on asserting
    # "Eight ND donors, two T1D, no Aab+" and a floor table topping out at 4,371 islets, while
    # `within_status_donor_count` and `islets_ordered` in the SAME dict said {ND 8, Aab+ 3, T1D 9}
    # and 8,665. A file that contradicts itself is worse than one that says less, and this is the
    # leak the app's own rules name: a cohort's size, its donor ids and its method printed as a
    # measurement of whichever dataset was loaded. All three are derived from this fit now.
    floor_cells = adata.obs["total_cells"].to_numpy(float)
    floor_sweep = {
        "note": ("islets at or above each candidate floor, counted on THIS fit. `rho` is the "
                 f"Spearman between n_neighbors {N_NEIGHBORS} and {N_NEIGHBORS_ALT} and is "
                 "measured only at the floor actually used; at the others it is Unavailable "
                 "rather than carried over from another cohort's sweep."),
        "floor_used": COMPOSITION_FLOOR_CELLS,
    }
    for candidate in (0, 5, 10, 20, 30, 50):
        kept = int(np.sum(np.nan_to_num(floor_cells, nan=-1.0) >= candidate))
        floor_sweep[str(candidate)] = {
            "islets": kept,
            "rho": (float(rho_neighbors) if candidate == COMPOSITION_FLOOR_CELLS else None),
        }

    # Which staining groups hold a single donor *within* a disease status: correcting on
    # stain_group would take that donor's disease contribution out with its group effect. This is
    # the actual test the old sentence was asserting about S4/S5, run on whatever cohort is here.
    fit_stain = fitted.obs["stain_group"].astype(str).to_numpy()
    single_donor_groups, groups_with_contrast = {}, []
    for s_ in sorted(set(fit_status)):
        here = fit_status == s_
        for group in sorted(set(fit_stain[here])):
            donors_here = sorted(set(fit_donor[here & (fit_stain == group)]))
            if len(donors_here) == 1:
                single_donor_groups.setdefault(s_, {})[group] = donors_here[0]
    for group in sorted(set(fit_stain)):
        if len(set(fit_status[fit_stain == group])) > 1:
            groups_with_contrast.append(group)
    confounded = "; ".join(
        f"{s_}: " + ", ".join(f"{g} is donor {d} alone" for g, d in sorted(v.items()))
        for s_, v in sorted(single_donor_groups.items())) or "none on this cohort"
    batch_correction_reason = (
        "Correction needs a batch variable carrying within-group disease contrast. Measured on "
        f"this fit, the staining groups holding a single donor within a status are -- {confounded}"
        " -- so correcting on stain_group would remove those donors' contributions with their "
        "group effects. Groups that do carry within-group disease contrast: "
        + (", ".join(groups_with_contrast) or "none")
        + ". The cohort is stratified and the confound reported instead.")

    status_phrase = ", ".join(
        f"{n} {s_}" for s_, n in sorted(within_status_donors.items(), key=lambda kv: -kv[1]))
    eta_phrase = "; ".join(
        f"within {s_} ({within_status_donors[s_]} donors) donor explains "
        f"{within_status_donor_eta2[s_] * 100:.1f}% of the pseudotime"
        + (" -- gated" if s_ in gated_statuses
           else f" -- reported only, under {DONOR_ETA2_MIN_DONORS} donors")
        for s_ in sorted(within_status_donor_eta2))
    cohort_limit = (
        f"{status_phrase} donors ({len(set(fit_donor))} in the fit). This is an "
        "islet-composition ordering, not a disease-progression axis, and must not be read as "
        f"one. {eta_phrase}. A status below {DONOR_ETA2_MIN_DONORS} donors yields one pairwise "
        "contrast rather than a variance decomposition, so its number describes the build and "
        f"does not gate it; the ceiling on a gated status is {DONOR_ETA2_CEILING}.")

    adata.uns["iroot"] = int(np.where(adata.obs_names == ordered.obs_names[root])[0][0])
    adata.uns["trajectory_validation"] = {
        **{f"gate_{k}": v for k, v in gates.items()},
        **reported,
        "hard_pass": bool(all(gates.values())),
        "root": str(ordered.obs_names[root]),
        "root_donor": str(ordered.obs["imageid"].iloc[root]),
        "root_pool_size": root_meta["pool_size"],
        "alr_pseudocount": ALR_PSEUDOCOUNT,
        "n_neighbors": N_NEIGHBORS,
        "composition_floor_cells": COMPOSITION_FLOOR_CELLS,
        "islets_ordered": int(eligible.sum()),
        "islets_below_floor": int((~eligible).sum()),
        "composition_floor_sweep": floor_sweep,
        "batch_correction": "none",
        "batch_correction_reason": batch_correction_reason,
        "cohort_limit": cohort_limit,
        "islets_unavailable_no_core_endocrine": int(len(undefined_by_donor)
                                                    and sum(undefined_by_donor.values())),
        "islets_unavailable_by_donor": undefined_by_donor,
    }
    adata.uns["embedding_batch"] = {
        "note": ("eta-squared of the 2-D UMAP on each donor-level grouping, with a null that "
                 "permutes the donor->group map (islets are nested in donors, so an islet-level "
                 "shuffle is not a null). Reported, never gated."),
        "coarse_batch_is_aliased_with_status": True,
        "statuses_containing_both_batches": statuses_with_both,
        "pooled": embedding_batch,
        "within_status": within_status_embedding,
    }
    adata.uns["phenotype_names"] = classes
    adata.uns["run_id"] = S.RUN_ID

    # anndata 0.13 writes pandas' new string dtype as `nullable-string-array`, an encoding the
    # app's R `anndata` reader does not understand -- it falls back to an integer index without
    # an error, so the Trajectory feature selector offered "1, 2, 3 ..." instead of the class
    # names and `combined_islet_id` arrived empty. Object dtype writes the older `string-array`
    # the reader does understand. Check the encodings, not just that the file opens.
    def _string_array_encoding(frame: pd.DataFrame) -> pd.DataFrame:
        # pandas 3 infers `str` from an object-dtype array, so dtype=object must be explicit.
        frame.index = pd.Index(np.asarray(frame.index, dtype=object), dtype=object)
        for column in frame.columns:
            if pd.api.types.is_string_dtype(frame[column]) and frame[column].dtype != object:
                frame[column] = frame[column].astype(object)
        return frame

    adata.obs = _string_array_encoding(adata.obs)
    adata.var = _string_array_encoding(adata.var)

    target = out_dir / "adata_ins_root.h5ad"
    if not all(gates.values()):
        rejected = out_dir / "adata_ins_root.rejected.h5ad"
        adata.write_h5ad(rejected)
        print(f"\nMETHOD CHECK FAILED -> {rejected} (the deployed file is untouched)")
        return 1

    adata.write_h5ad(target)
    print(f"\n{adata.n_obs:,} islets x {adata.n_vars} classes "
          f"({int(eligible.sum()):,} with a pseudotime) -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
