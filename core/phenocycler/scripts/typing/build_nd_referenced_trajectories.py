#!/usr/bin/env python3
"""Fit an ND-only MFA reference, then donor-specific diffusion pseudotime paths.

The reference model and settings are learned from ND donors only. Objects remain in the output
when a trajectory coordinate is unavailable or unstable. Each connected graph component is rooted
separately; pseudotime values have meaning only within that donor and component.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing as mp
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kendalltau
from sklearn.neighbors import kneighbors_graph
from scipy.sparse.csgraph import connected_components


K_VALUES = (5, 10, 20, 40, 80)
DIMS_VALUES = (2, 5, 10, 20)
MASK_REPLICATES = 20
MASK_FRACTION = 0.20
OBJECT_SUBSAMPLES = 20
OBJECT_SUBSAMPLE_FRACTION = 0.80
SEED_BASE = "io60-nd-reference-2026-10-03"
MIN_ND_DONORS_FOR_VALIDATION = 3
# A one-coordinate representation still defines a graph with many diffusion modes.
# Input coordinate rank and graph eigenvector count are separate quantities.
MIN_TRAJECTORY_RANK = 1
FRACTION_ROUNDOFF_TOLERANCE = 1e-12
EXCLUDED_FEATURE_SUFFIXES = ("unresolved_fraction",)
BIOLOGICAL_GROUPS = ("endocrine", "vascular", "lymphatic", "neural", "stromal", "immune")
# User-directed completion policy: retain stromal composition as unavailable, but do not block
# an ND reference fit based on the five measured groups or replace stromal identity with a signal
# proxy. Every other biological group remains required for the primary trajectory.
NONBLOCKING_MISSING_GROUPS = frozenset({"stromal"})


def seed_for(*parts: object) -> int:
    value = "|".join([SEED_BASE, *map(str, parts)])
    return int(hashlib.sha256(value.encode()).hexdigest()[:8], 16)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def nd_frame_sha256(frame: pd.DataFrame, features: list[str]) -> str:
    availability = ["available__" + f.removeprefix("feature__") for f in features]
    keep = ["donor_id", "endocrine_object_id", "donor_status"] + features + availability
    keep = [c for c in keep if c in frame]
    nd = frame[frame["donor_status"].astype(str).str.upper().eq("ND")][keep].copy()
    nd = nd.sort_values(["donor_id", "endocrine_object_id"], kind="stable").reset_index(drop=True)
    return hashlib.sha256(pd.util.hash_pandas_object(nd, index=False).to_numpy().tobytes()).hexdigest()


def group_of(feature: str) -> str:
    return feature.removeprefix("feature__").split("__", 1)[0]


def feature_group_gate(groups: list[str]) -> dict:
    missing = sorted(set(BIOLOGICAL_GROUPS) - set(groups))
    return {"missing_groups": missing,
            "nonblocking_missing_groups": sorted(set(missing) & NONBLOCKING_MISSING_GROUPS),
            "blocking_missing_groups": sorted(set(missing) - NONBLOCKING_MISSING_GROUPS)}


def is_fraction(feature: str) -> bool:
    return feature.endswith("_fraction")


def select_features(objects: pd.DataFrame) -> tuple[list[str], dict]:
    nd = objects[objects["donor_status"].astype(str).str.upper().eq("ND")]
    donors = sorted(nd["donor_id"].astype(str).unique())
    candidates = [c for c in objects if c.startswith("feature__")
                  and not any(c.endswith(suffix) for suffix in EXCLUDED_FEATURE_SUFFIXES)]
    selected, rejected = [], {}
    for feature in candidates:
        availability = "available__" + feature.removeprefix("feature__")
        if availability not in objects:
            rejected[feature] = "availability_field_missing"
            continue
        donor_ok = {}
        for donor, frame in nd.groupby(nd["donor_id"].astype(str)):
            mask = frame[availability].fillna(False).astype(bool)
            donor_ok[donor] = bool(mask.any() and pd.to_numeric(frame.loc[mask, feature],
                                                                  errors="coerce").notna().any())
        missing = sorted(d for d in donors if not donor_ok.get(d, False))
        if missing:
            rejected[feature] = "not_measured_in_all_ND_donors:" + ",".join(missing)
        else:
            selected.append(feature)
    return selected, {"nd_donors": donors, "candidate_features": candidates,
                      "selected_features": selected, "rejected_features": rejected,
                      "groups_selected": {g: [f for f in selected if group_of(f) == g]
                                          for g in sorted({group_of(f) for f in candidates})}}


def transform_frame(frame: pd.DataFrame, features: list[str]) -> np.ndarray:
    # pandas may expose a read-only view of Arrow-backed Parquet columns; transformations below
    # are in-place, so always materialize a writable float array.
    values = frame[features].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float, copy=True)
    for j, feature in enumerate(features):
        finite = np.isfinite(values[:, j])
        if is_fraction(feature):
            observed = values[finite, j]
            if np.any((observed < -FRACTION_ROUNDOFF_TOLERANCE) |
                      (observed > 1 + FRACTION_ROUNDOFF_TOLERANCE)):
                raise ValueError(f"{feature}: fraction outside [0,1]")
            # Polygon intersections can produce a few ulps above one. Correct only these
            # floating-point excursions; materially invalid coverage remains a hard error.
            values[finite, j] = np.clip(observed, 0.0, 1.0)
        else:
            if np.any(values[finite, j] < 0):
                raise ValueError(f"{feature}: abundance ratio/density must be nonnegative")
            values[finite, j] = np.log1p(values[finite, j])
    return values


def fraction_roundoff_audit(frame: pd.DataFrame, features: list[str]) -> dict:
    audit = {}
    for feature in features:
        if not is_fraction(feature):
            continue
        values = pd.to_numeric(frame[feature], errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        outside = finite[(finite < 0) | (finite > 1)]
        if len(outside) and np.any((outside < -FRACTION_ROUNDOFF_TOLERANCE) |
                                   (outside > 1 + FRACTION_ROUNDOFF_TOLERANCE)):
            raise ValueError(f"{feature}: fraction outside [0,1] beyond roundoff tolerance")
        if len(outside):
            audit[feature] = {
                "count_clipped": int(len(outside)),
                "minimum_raw_value": float(outside.min()),
                "maximum_raw_value": float(outside.max()),
                "absolute_tolerance": FRACTION_ROUNDOFF_TOLERANCE,
            }
    return audit


def donor_weights(donors: np.ndarray) -> np.ndarray:
    weights = np.zeros(len(donors), dtype=float)
    for donor in np.unique(donors):
        idx = np.flatnonzero(donors == donor)
        weights[idx] = 1.0 / len(idx)
    return weights


def fit_scaling(X: np.ndarray, donor: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    means = np.full(X.shape[1], np.nan)
    sds = np.full(X.shape[1], np.nan)
    audit = {}
    for j in range(X.shape[1]):
        donor_means, donor_variances = [], []
        for name in sorted(np.unique(donor)):
            vals = X[(donor == name) & np.isfinite(X[:, j]), j]
            if len(vals):
                donor_means.append(float(vals.mean()))
                donor_variances.append(float(np.mean((vals - vals.mean()) ** 2)))
        if not donor_means:
            continue
        means[j] = np.mean(donor_means)
        variance = np.mean([v + (m - means[j]) ** 2
                            for m, v in zip(donor_means, donor_variances)])
        sds[j] = np.sqrt(variance) if variance > 1e-12 else 1.0
        audit[str(j)] = {"donors_observed": len(donor_means), "equal_donor_mean": means[j],
                         "equal_donor_sd": sds[j]}
    return means, sds, audit


def iterative_low_rank_impute(X: np.ndarray, mask: np.ndarray, weights: np.ndarray,
                              rank: int, iterations: int = 8) -> np.ndarray:
    """Low-rank completion with missing values initialized at the ND-centered mean, not zero abundance."""
    fill = np.where(mask, X, 0.0)
    root_weight = np.sqrt(weights / np.mean(weights[weights > 0]))[:, None]
    for _ in range(iterations):
        weighted = fill * root_weight
        u, s, vt = np.linalg.svd(weighted, full_matrices=False)
        r = max(1, min(rank, len(s)))
        reconstructed = ((u[:, :r] * s[:r]) @ vt[:r]) / root_weight
        fill[~mask] = reconstructed[~mask]
        fill[mask] = X[mask]
    return fill


def fit_mfa(frame: pd.DataFrame, features: list[str], rank: int) -> dict:
    X = transform_frame(frame, features)
    donors = frame["donor_id"].astype(str).to_numpy()
    means, sds, scale_audit = fit_scaling(X, donors)
    if np.isnan(means).any() or np.isnan(sds).any():
        raise ValueError("selected feature has no observations in this ND training fold")
    Z = (X - means) / sds
    mask = np.isfinite(Z)
    weights = donor_weights(donors)
    Z_complete = iterative_low_rank_impute(Z, mask, weights,
                                           min(rank, max(1, min(Z.shape) - 1)))
    groups = np.array([group_of(f) for f in features])
    group_weights = {}
    group_first_singular = {}
    for group in sorted(set(groups)):
        cols = np.flatnonzero(groups == group)
        A = Z_complete[:, cols] * np.sqrt(weights / np.mean(weights))[:, None]
        singular = np.linalg.svd(A, compute_uv=False)
        first = float(singular[0]) if len(singular) and singular[0] > 1e-10 else 1.0
        group_weights[group] = 1.0 / first
        group_first_singular[group] = first
    gw = np.array([group_weights[g] for g in groups])
    balanced = Z_complete * gw
    weighted = balanced * np.sqrt(weights / np.mean(weights))[:, None]
    _, singular, vt = np.linalg.svd(weighted, full_matrices=False)
    max_rank = max(1, min(rank, vt.shape[0], vt.shape[1]))
    loadings = vt[:max_rank].T
    return {"features": features, "feature_groups": groups.tolist(), "means": means,
            "sds": sds, "group_weights": group_weights, "loadings": loadings,
            "group_first_singular": group_first_singular,
            "singular_values": singular[:max_rank], "scaling_audit": scale_audit,
            "training_donors": sorted(set(donors)), "training_rows": len(frame),
            "rank": max_rank}


def project_mfa(frame: pd.DataFrame, model: dict, rank: int, require_groups: bool = True):
    features = model["features"]
    X = transform_frame(frame, features)
    Z = (X - model["means"]) / model["sds"]
    groups = np.array(model["feature_groups"])
    gw = np.array([model["group_weights"][g] for g in groups])
    loadings = model["loadings"][:, :rank]
    scores = np.full((len(frame), rank), np.nan)
    valid = np.zeros(len(frame), dtype=bool)
    n_observed = np.isfinite(Z).sum(axis=1)
    available_groups = np.zeros(len(frame), dtype=int)
    reference_distance = np.full(len(frame), np.nan)
    reference_distance_complete = np.full(len(frame), np.nan)
    for i in range(len(frame)):
        observed = np.isfinite(Z[i])
        present_groups = sorted(set(groups[observed]))
        available_groups[i] = len(present_groups)
        if not observed.any():
            continue
        x = Z[i, observed] * gw[observed]
        reference_distance[i] = float(np.sqrt(np.mean(x ** 2) * len(features)))
        if observed.all():
            reference_distance_complete[i] = float(np.linalg.norm(x))
        if n_observed[i] < rank or (require_groups and len(present_groups) < len(model["group_weights"])):
            continue
        local_basis = loadings[observed]
        score, _, _, _ = np.linalg.lstsq(local_basis, x, rcond=None)
        if len(score) == rank and np.isfinite(score).all():
            scores[i] = score
            valid[i] = True
    return {"scores": scores, "projectable": valid, "n_features_observed": n_observed,
            "n_groups_observed": available_groups, "reference_distance_observed": reference_distance,
            "reference_distance_complete": reference_distance_complete}


def masked_reconstruction_fold(job):
    nd, features, rank_max, heldout = job
    train = nd[nd["donor_id"].astype(str) != heldout].reset_index(drop=True)
    valid = nd[nd["donor_id"].astype(str) == heldout].reset_index(drop=True)
    feasible = min(rank_max, max(1, len(train) - 1), len(features))
    model = fit_mfa(train, features, feasible)
    values = transform_frame(valid, features)
    Z = (values - model["means"]) / model["sds"]
    mask_observed = np.isfinite(Z)
    groups = np.array(model["feature_groups"])
    gw = np.array([model["group_weights"][g] for g in groups])
    errors_by_rank = {r: [] for r in range(1, rank_max + 1)}
    for rep in range(MASK_REPLICATES):
        rng = np.random.default_rng(seed_for("rank-mask", heldout, rep))
        observed_positions = np.flatnonzero(mask_observed)
        take = rng.choice(observed_positions, size=max(1, int(round(
            MASK_FRACTION * len(observed_positions)))), replace=False)
        hidden = np.zeros(mask_observed.size, dtype=bool)
        hidden[take] = True
        hidden = hidden.reshape(mask_observed.shape)
        keep = mask_observed & ~hidden
        for rank in range(1, rank_max + 1):
            basis = model["loadings"][:, :rank]
            errors, truth = [], []
            for i in range(len(valid)):
                visible = keep[i]
                if visible.sum() < rank:
                    continue
                local_basis = basis[visible]
                x = Z[i, visible] * gw[visible]
                score, _, _, _ = np.linalg.lstsq(local_basis, x, rcond=None)
                predict_mask = hidden[i]
                if not predict_mask.any():
                    continue
                recon = score @ basis.T
                errors.extend((recon[predict_mask] - Z[i, predict_mask] * gw[predict_mask]).tolist())
                truth.extend((Z[i, predict_mask] * gw[predict_mask]).tolist())
            if errors:
                errors_by_rank[rank].append(float(np.sqrt(np.mean(np.square(errors)))))
    return heldout, errors_by_rank


def masked_reconstruction_cv(nd: pd.DataFrame, features: list[str], rank_max: int, workers: int = 1):
    donors = sorted(nd["donor_id"].astype(str).unique())
    max_training_rows = min(len(nd[nd["donor_id"].astype(str) != d]) for d in donors)
    rank_max = max(1, min(rank_max, len(features), max_training_rows - 1))
    donor_errors: dict[int, dict[str, list[float]]] = {r: {d: [] for d in donors}
                                                       for r in range(1, rank_max + 1)}
    slim = nd[["donor_id", *features]].copy()
    jobs = [(slim, features, rank_max, heldout) for heldout in donors]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            outputs = list(pool.map(masked_reconstruction_fold, jobs))
    else:
        outputs = list(map(masked_reconstruction_fold, jobs))
    for heldout, fold_errors in outputs:
        for rank, values in fold_errors.items():
            donor_errors[rank][heldout] = values
    rows = []
    for rank, by_donor in donor_errors.items():
        donor_scores = {d: float(np.mean(v)) for d, v in by_donor.items() if v}
        values = np.array(list(donor_scores.values()), dtype=float)
        mean = float(values.mean()) if len(values) else np.nan
        se = float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan
        rows.append({"rank": rank, "mean_donor_balanced_rmse": mean,
                     "se_across_heldout_donors": se,
                     "heldout_donors_validated": len(values),
                     "per_donor_rmse": donor_scores})
    table = pd.DataFrame([{k: row[k] for k in ("rank", "mean_donor_balanced_rmse",
                                                "se_across_heldout_donors",
                                                "heldout_donors_validated")} for row in rows])
    chosen = select_rank_one_se(table, MIN_TRAJECTORY_RANK)
    if chosen is None:
        return None, rows
    return chosen, rows


def select_rank_one_se(table: pd.DataFrame, minimum_rank: int = MIN_TRAJECTORY_RANK) -> int | None:
    """Choose the smallest validated positive input rank within one SE."""
    feasible = table[(pd.to_numeric(table["rank"], errors="coerce") >= minimum_rank)
                     & pd.to_numeric(table["mean_donor_balanced_rmse"], errors="coerce").notna()]
    if feasible.empty:
        return None
    best = feasible.loc[feasible["mean_donor_balanced_rmse"].idxmin()]
    threshold = float(best["mean_donor_balanced_rmse"] + best["se_across_heldout_donors"])
    if not np.isfinite(threshold):
        threshold = float(best["mean_donor_balanced_rmse"])
    selected = feasible.loc[feasible["mean_donor_balanced_rmse"] <= threshold, "rank"]
    return int(selected.min()) if len(selected) else int(best["rank"])


def graph_dpt(scores: np.ndarray, object_ids: np.ndarray, reference_center: np.ndarray,
              n_neighbors: int, dimensions: int) -> dict:
    n = len(scores)
    result = {"pseudotime": np.full(n, np.nan), "component": np.full(n, -1, dtype=int),
              "root_object_id": np.full(n, "", dtype=object), "component_size": np.zeros(n, dtype=int),
              "status": np.full(n, "unresolved", dtype=object)}
    if n == 0:
        return result
    k = min(n_neighbors, n - 1)
    if k < 1:
        result["component"][:] = 0
        result["component_size"][:] = 1
        result["pseudotime"][:] = 0.0
        result["root_object_id"][:] = object_ids.astype(str)
        result["status"][:] = "singleton_component"
        return result
    graph = kneighbors_graph(scores, n_neighbors=k, mode="connectivity", include_self=False)
    graph = graph.maximum(graph.T)
    n_comp, labels = connected_components(graph, directed=False)
    ids_text = object_ids.astype(str)
    comp_sets = []
    for old_label in range(n_comp):
        idx = np.flatnonzero(labels == old_label)
        idx = idx[np.argsort(ids_text[idx], kind="stable")]
        comp_sets.append(idx)
    comp_sets.sort(key=lambda ix: ids_text[ix[0]])
    for comp_id, indices in enumerate(comp_sets):
        result["component"][indices] = comp_id
        result["component_size"][indices] = len(indices)
        if len(indices) == 1:
            result["pseudotime"][indices] = 0.0
            result["root_object_id"][indices] = ids_text[indices]
            result["status"][indices] = "singleton_component"
            continue
        if len(indices) < dimensions + 2:
            result["status"][indices] = "component_too_small_for_diffusion_dimensions"
            continue
        # Root in the frozen low-dimensional reference representation. Distances in the
        # incompletely observed input-feature space would favor objects with fewer measured
        # features, so missingness must not determine the root.
        root_dist = np.linalg.norm(scores[indices] - reference_center[None, :], axis=1)
        root_dist[~np.isfinite(root_dist)] = np.inf
        min_dist = float(np.nanmin(root_dist))
        tied = indices[np.isclose(root_dist, min_dist, rtol=1e-12, atol=1e-12)]
        root_idx = min(tied, key=lambda i: ids_text[i])
        local_root = int(np.flatnonzero(indices == root_idx)[0])
        subset = ad.AnnData(X=np.asarray(scores[indices], dtype=np.float32))
        subset.obs_names = pd.Index(ids_text[indices])
        subset.uns["iroot"] = local_root
        local_k = min(k, len(indices) - 1)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sc.pp.neighbors(subset, n_neighbors=local_k, use_rep="X", metric="euclidean",
                                transformer="sklearn", random_state=0)
                sc.tl.diffmap(subset, n_comps=min(dimensions + 1, len(indices) - 1))
                subset.uns["iroot"] = local_root
                sc.tl.dpt(subset, n_dcs=dimensions)
            values = pd.to_numeric(subset.obs["dpt_pseudotime"], errors="coerce").to_numpy(float)
            finite = np.isfinite(values)
            result["pseudotime"][indices[finite]] = values[finite]
            result["root_object_id"][indices] = ids_text[root_idx]
            result["status"][indices[finite]] = "ordered"
            result["status"][indices[~finite]] = "dpt_nonfinite"
        except Exception as exc:  # retain the component and reason instead of dropping objects
            result["root_object_id"][indices] = ids_text[root_idx]
            result["status"][indices] = "dpt_error:" + type(exc).__name__
    return result


def pairwise_within_component_tau(a: dict, b: dict, ids: list[str]) -> tuple[float, float]:
    index_a = {str(k): i for i, k in enumerate(a["ids"])}
    index_b = {str(k): i for i, k in enumerate(b["ids"])}
    common = [x for x in ids if x in index_a and x in index_b]
    finite_common = [x for x in common if np.isfinite(a["pseudo"][index_a[x]])
                     and np.isfinite(b["pseudo"][index_b[x]])]
    if not common:
        return np.nan, 0.0
    coverage = len(finite_common) / len(ids) if ids else 0.0

    # Compare ranks only inside intersections of components from the two resamples. The direct
    # nested-pair implementation scales quadratically in thousands of objects and tens of
    # thousands of resample comparisons. Kendall tau-b gives C-D in O(n log n); correcting by
    # the exact comparable-pair count recovers the original (C-D)/(C+D) estimand while ignoring
    # ties in either ordering exactly as the previous implementation did.
    ia = np.asarray([index_a[x] for x in finite_common], dtype=int)
    ib = np.asarray([index_b[x] for x in finite_common], dtype=int)
    comp_pairs = np.column_stack([a["component"][ia], b["component"][ib]])
    _, component_codes = np.unique(comp_pairs, axis=0, return_inverse=True)
    sorted_rows = np.argsort(component_codes, kind="stable")
    split_points = np.flatnonzero(np.diff(component_codes[sorted_rows])) + 1
    numerator = 0.0
    denominator = 0
    for group_rows in np.split(sorted_rows, split_points):
        x = np.asarray(a["pseudo"][ia[group_rows]], dtype=float)
        y = np.asarray(b["pseudo"][ib[group_rows]], dtype=float)
        n = len(x)
        if n < 2:
            continue
        total_pairs = n * (n - 1) // 2

        def tied_pairs(values: np.ndarray) -> int:
            _, counts = np.unique(values, return_counts=True)
            return int(np.sum(counts * (counts - 1) // 2))

        tied_x = tied_pairs(x)
        tied_y = tied_pairs(y)
        _, joint_counts = np.unique(np.column_stack([x, y]), axis=0, return_counts=True)
        tied_both = int(np.sum(joint_counts * (joint_counts - 1) // 2))
        comparable = total_pairs - tied_x - tied_y + tied_both
        if comparable <= 0:
            continue
        tau_b = kendalltau(x, y, variant="b").statistic
        if not np.isfinite(tau_b):
            continue
        c_minus_d = float(tau_b * np.sqrt((total_pairs - tied_x) * (total_pairs - tied_y)))
        numerator += c_minus_d
        denominator += comparable
    tau = numerator / denominator if denominator else np.nan
    return float(tau) if np.isfinite(tau) else np.nan, coverage


def graph_validation_config(job):
    donor, data, cfg = job
    ids = data["ids"]
    samples = []
    for rep in range(OBJECT_SUBSAMPLES):
        rng = np.random.default_rng(seed_for("graph-validation", donor, rep))
        size = max(2, int(np.floor(OBJECT_SUBSAMPLE_FRACTION * len(ids))))
        chosen = np.sort(rng.choice(len(ids), size=min(size, len(ids)), replace=False))
        samples.append(chosen)
    k, dimensions = cfg
    results = []
    for chosen in samples:
        run = graph_dpt(data["scores"][chosen], ids[chosen],
                        data["reference_center"], k, dimensions)
        results.append({"ids": ids[chosen].tolist(), "pseudo": run["pseudotime"],
                        "component": run["component"]})
    stability = []
    for i in range(len(results)):
        for j in range(i + 1, len(results)):
            tau, coverage = pairwise_within_component_tau(results[i], results[j], ids.tolist())
            if np.isfinite(tau):
                stability.append(tau * coverage)
            else:
                stability.append(0.0)
    donor_score = float(np.mean(stability)) if stability else 0.0
    return donor, cfg, donor_score


def choose_graph_settings(nd: pd.DataFrame, features: list[str], rank: int,
                          reference_model: dict, workers: int = 1) -> tuple[dict | None, list[dict]]:
    donor_names = sorted(nd["donor_id"].astype(str).unique())
    fold_inputs = {}
    for heldout in donor_names:
        train = nd[nd["donor_id"].astype(str) != heldout].reset_index(drop=True)
        held = nd[nd["donor_id"].astype(str) == heldout].reset_index(drop=True)
        model = fit_mfa(train, features, rank)
        projection = project_mfa(held, model, rank)
        train_projection = project_mfa(train, model, rank)
        train_scores = train_projection["scores"][train_projection["projectable"]]
        reference_center = np.nanmean(train_scores, axis=0)
        keep = projection["projectable"]
        frame = held.loc[keep].reset_index(drop=True)
        scores = projection["scores"][keep]
        fold_inputs[heldout] = {"frame": frame, "scores": scores,
                                "reference_center": reference_center,
                                "ids": frame["endocrine_object_id"].astype(str).to_numpy()}
    candidates = []
    for k in K_VALUES:
        for dimensions in DIMS_VALUES:
            feasible = all(len(v["ids"]) > k and len(v["ids"]) >= dimensions + 2
                           for v in fold_inputs.values())
            if feasible:
                candidates.append((k, dimensions))
    if not candidates:
        return None, []
    scores_by_config = {cfg: [] for cfg in candidates}
    per_donor_config = {cfg: {} for cfg in candidates}
    jobs = [(donor, data, cfg) for donor, data in fold_inputs.items()
            if len(data["ids"]) >= 2 for cfg in candidates]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            outputs = list(pool.map(graph_validation_config, jobs))
    else:
        outputs = list(map(graph_validation_config, jobs))
    for donor, cfg, donor_score in outputs:
        per_donor_config[cfg][donor] = donor_score
        scores_by_config[cfg].append(donor_score)
    rows = []
    for cfg, scores in scores_by_config.items():
        values = np.asarray(scores, dtype=float)
        rows.append({"n_neighbors": cfg[0], "diffusion_dimensions": cfg[1],
                     "mean_donor_balanced_rank_reproducibility": float(values.mean()),
                     "se_across_heldout_donors": float(values.std(ddof=1) / np.sqrt(len(values)))
                     if len(values) > 1 else np.nan,
                     "heldout_nd_donors": len(values), "per_donor_scores": per_donor_config[cfg]})
    best = max(rows, key=lambda x: x["mean_donor_balanced_rank_reproducibility"])
    se = best["se_across_heldout_donors"]
    threshold = best["mean_donor_balanced_rank_reproducibility"] - (se if np.isfinite(se) else 0.0)
    eligible = [r for r in rows
                if r["mean_donor_balanced_rank_reproducibility"] >= threshold]
    chosen = min(eligible, key=lambda r: (r["diffusion_dimensions"], r["n_neighbors"]))
    return {"n_neighbors": chosen["n_neighbors"],
            "diffusion_dimensions": chosen["diffusion_dimensions"],
            "selection_rule": "maximum donor-balanced within-component Kendall ordering stability times finite-object coverage; choose lower dimensions, then lower neighbors, within 1 SE",
            "selected_validation_score": chosen["mean_donor_balanced_rank_reproducibility"],
            "best_validation_score": best["mean_donor_balanced_rank_reproducibility"],
            "one_se_threshold": threshold}, rows


def donor_trajectory_worker(job):
    donor, frame, model, rank, settings, center, groups, group_gate = job
    frame = frame.reset_index(drop=True)
    projection = project_mfa(frame, model, rank)
    ids = frame['endocrine_object_id'].astype(str).to_numpy()
    valid_idx = np.flatnonzero(projection['projectable'])
    base_run = graph_dpt(projection['scores'][valid_idx], ids[valid_idx], center, settings['n_neighbors'], settings['diffusion_dimensions'])
    full = pd.DataFrame({'endocrine_object_id': ids[valid_idx], 'component_id': base_run['component'], 'dpt_pseudotime': base_run['pseudotime'], 'component_root_object_id': base_run['root_object_id'], 'component_size': base_run['component_size'], 'trajectory_status': base_run['status']})
    full['pseudotime_rank_fraction'] = np.nan
    for _, idx in full.groupby('component_id').groups.items():
        take = np.asarray(list(idx), dtype=int)
        finite = np.isfinite(full.loc[take, 'dpt_pseudotime'].to_numpy())
        if finite.any():
            ranks = pd.Series(full.loc[take, 'dpt_pseudotime'].to_numpy()).rank(method='average', pct=True).to_numpy()
            full.loc[take, 'pseudotime_rank_fraction'] = ranks
    sampled_pseudo: dict[str, list[float]] = {x: [] for x in ids}
    sampled_tau, sampled_coverage = ([], [])
    if len(valid_idx) >= 2:
        samples = []
        for rep in range(OBJECT_SUBSAMPLES):
            rng = np.random.default_rng(seed_for('production-resample', donor, rep))
            size = max(2, int(np.floor(OBJECT_SUBSAMPLE_FRACTION * len(valid_idx))))
            chosen = np.sort(rng.choice(valid_idx, size=min(size, len(valid_idx)), replace=False))
            run = graph_dpt(projection['scores'][chosen], ids[chosen], center, settings['n_neighbors'], settings['diffusion_dimensions'])
            pct = np.full(len(chosen), np.nan)
            for component_id in np.unique(run['component']):
                if component_id < 0:
                    continue
                local = np.flatnonzero(run['component'] == component_id)
                local_time = run['pseudotime'][local]
                finite = np.isfinite(local_time)
                if finite.any():
                    pct[local[finite]] = pd.Series(local_time[finite]).rank(method='average', pct=True).to_numpy()
            samples.append({'ids': ids[chosen].tolist(), 'pseudo': run['pseudotime'], 'component': run['component'], 'rank': pct, 'status': run['status']})
            for local_i, object_id in enumerate(ids[chosen]):
                if np.isfinite(run['pseudotime'][local_i]) and run['component'][local_i] >= 0:
                    sampled_pseudo[object_id].append(float(pct[local_i]))
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                tau, coverage = pairwise_within_component_tau(samples[i], samples[j], ids[valid_idx].tolist())
                if np.isfinite(tau):
                    sampled_tau.append(tau * coverage)
                else:
                    sampled_tau.append(0.0)
                sampled_coverage.append(coverage)
    by_id = {str(row.endocrine_object_id): i for i, row in full.iterrows()}
    full['resampling_finite_fraction'] = 0.0
    full['resampling_rank_sd_within_component'] = np.nan
    for i, row in full.iterrows():
        object_id = str(row['endocrine_object_id'])
        full['resampling_finite_fraction'] = full['resampling_finite_fraction'].astype(float)
        vals = sampled_pseudo.get(object_id, [])
        full.loc[i, 'resampling_finite_fraction'] = len(vals) / OBJECT_SUBSAMPLES
        if len(vals) >= 2:
            full.loc[i, 'resampling_rank_sd_within_component'] = float(np.std(vals, ddof=1))
    stability = {'donor_id': donor, 'objects_total': len(frame), 'objects_projectable': len(valid_idx), 'objects_finite_full_graph': int(np.isfinite(base_run['pseudotime']).sum()), 'finite_fraction_full_graph': float(np.isfinite(base_run['pseudotime']).mean()) if len(base_run['pseudotime']) else 0.0, 'within_component_pairwise_rank_stability': float(np.mean(sampled_tau)) if sampled_tau else None, 'mean_resampling_finite_fraction': float(np.mean(sampled_coverage)) if sampled_coverage else None, 'component_count': int(pd.Series(base_run['component']).nunique())}
    full['donor_id'] = donor
    output = frame[['donor_id', 'donor_status', 'endocrine_object_id', 'islet_key', 'original_object_class', 'n_cells_in_object', 'endocrine_denominator']].merge(full, on=['donor_id', 'endocrine_object_id'], how='left', validate='one_to_one')
    output['trajectory_status'] = output['trajectory_status'].fillna('missing_feature_group_or_rank')
    output['mfa_projection_available'] = projection['projectable']
    output['n_selected_features_observed'] = projection['n_features_observed']
    output['n_feature_groups_observed'] = projection['n_groups_observed']
    output['nd_reference_distance_mfa'] = np.where(projection['projectable'], np.linalg.norm(projection['scores'] - center[None, :], axis=1), np.nan)
    output['nd_reference_distance_observed_features'] = projection['reference_distance_observed']
    output['nd_reference_distance_complete_features'] = projection['reference_distance_complete']
    output['within_donor_component_normalized_time'] = output['pseudotime_rank_fraction']
    output['trajectory_feature_groups_used'] = ';'.join(groups)
    output['trajectory_unavailable_feature_groups'] = ';'.join(group_gate['nonblocking_missing_groups'])
    for component in range(rank):
        output[f'mfa_coordinate_{component + 1}'] = projection['scores'][:, component]
    return (output, stability)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composition", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-nd-donors", type=int, default=MIN_ND_DONORS_FOR_VALIDATION)
    parser.add_argument("--rank-max", type=int, default=None)
    parser.add_argument("--workers", type=int, default=1,
                        help="parallel independent ND folds/configurations and donor jobs; validation/seeds unchanged")
    parser.add_argument("--rank-validation-from", type=Path,
                        help="reuse exact rank reconstruction errors from a manifest with identical ND inputs, features, and candidate ranks")
    args = parser.parse_args()
    if args.composition.is_dir():
        source = args.composition / "endocrine_objects.parquet"
        if not source.is_file():
            parts = sorted(args.composition.glob("donor_id=*/data.parquet"))
            if not parts:
                raise FileNotFoundError(f"no composition parquet under {args.composition}")
            objects = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
        else:
            objects = pd.read_parquet(source)
    else:
        source = args.composition
        objects = pd.read_parquet(source)
    args.out.mkdir(parents=True, exist_ok=True)
    if objects[["donor_id", "endocrine_object_id", "donor_status"]].isna().any().any():
        raise ValueError("composition requires donor, object, and disease status keys")
    if objects.duplicated(["donor_id", "endocrine_object_id"]).any():
        raise ValueError("composition object key is not unique")
    if "composition_analysis_eligible" in objects:
        if objects.composition_analysis_eligible.isna().any():
            raise ValueError("explicit composition analysis eligibility cannot be null")
        analysis_eligible = objects.composition_analysis_eligible.astype(bool)
    else:
        analysis_eligible = pd.Series(True, index=objects.index)
    analysis_objects = objects.loc[analysis_eligible].copy()
    nd = analysis_objects[analysis_objects["donor_status"].astype(str).str.upper().eq("ND")].copy()
    nd_donors = sorted(nd["donor_id"].astype(str).unique())
    features, feature_audit = select_features(analysis_objects)
    feature_audit["fraction_roundoff_corrections"] = fraction_roundoff_audit(analysis_objects, features)
    groups = sorted({group_of(f) for f in features})
    group_gate = feature_group_gate(groups)
    missing_groups = group_gate["missing_groups"]
    base = {"title": "ND-referenced endocrine-object trajectories",
            "composition_input": str(source.resolve()), "composition_sha256": sha256(source)
            if source.is_file() else None,
            "trajectory_builder": str(Path(__file__).resolve()),
            "trajectory_builder_sha256": sha256(Path(__file__).resolve()),
            "nd_donors": nd_donors, "nd_object_count": len(nd),
            "nd_reference_input_sha256": nd_frame_sha256(analysis_objects, features),
            "original_object_count": len(objects), "analysis_eligible_object_count": len(analysis_objects),
            "explicitly_excluded_object_count": int((~analysis_eligible).sum()),
            "explicitly_excluded_nd_object_count": int((~analysis_eligible & objects.donor_status.astype(str).str.upper().eq("ND")).sum()),
            "excluded_objects": objects.loc[~analysis_eligible, [c for c in ["donor_id", "endocrine_object_id", "composition_exclusion_reason", "excluded_cell_ids"] if c in objects]].to_dict("records"),
            "execution_workers": args.workers,
            "feature_selection": feature_audit,
            "feature_group_gate": group_gate,
            "rank_selection": {
                "criterion": "smallest rank within one SE of minimum leave-one-ND-donor-out masked-feature reconstruction RMSE",
                "minimum_feasible_rank": MIN_TRAJECTORY_RANK,
                "minimum_rank_reason": "one input coordinate can define a graph; diffusion modes are constrained by graph component size",
            },
            "rank_mask_protocol": {"replicates": MASK_REPLICATES, "fraction": MASK_FRACTION,
                                   "seed_base": SEED_BASE,
                                   "seed_namespace": "rank-mask"},
            "graph_dimension_feasibility": {
                "rule": "all declared diffusion dimensions are independent of MFA coordinate rank; each component requires at least diffusion_dimensions + 2 objects",
                "declared_diffusion_dimensions": list(DIMS_VALUES),
                "declared_neighbors": list(K_VALUES),
                "small_components": "retain component and report component_too_small_for_diffusion_dimensions",
            },
            "trajectory_feature_policy": {
                "mode": "reduced_measured_feature_set" if group_gate["nonblocking_missing_groups"] else "complete_six_group_set",
                "used_groups": groups,
                "unavailable_groups": group_gate["nonblocking_missing_groups"],
                "unavailable_group_handling": "stromal cell-body identity remains unavailable and is excluded from fitting; no signal proxy or inferred stromal call is substituted",
                "group_balance": "MFA group scaling is fit from observed ND variation across available groups only",
            },
            "fixed_guardrails": [
                "ND donors alone fit centering, scaling, group weights, and representation loadings.",
                "Non-ND labels/measurements do not select features, ranks, or graph settings.",
                "Ratios and densities use log1p; fractions remain fractions; only excursions <=1e-12 from [0,1] are clipped as floating-point roundoff; missing values are not zeros.",
                "No object cell-count floor is used; EndocrineSingle and EndocrineCluster remain eligible.",
                "Objects containing only typing-excluded cells remain in release tables with original counts/geometry, but do not enter ND feature selection, fitting, validation, or donor graphs. All cell-derived features omit quarantined cells.",
                "If stromal cell-body identity is unavailable, fit the ND reference on the five measured groups and retain stromal unavailability as an explicit limitation; do not infer a stromal call from marker intensity.",
                "Representation rank is selected by one-SE masked reconstruction among all validated positive ranks; graph diffusion dimension is independent of MFA input rank.",
                "Diffusion components are rooted separately at the object closest to the frozen ND center; IDs break exact ties.",
                "Pseudotime is within-donor and within-component; it is not a shared cross-donor clock.",
            ]}
    if len(nd_donors) < args.min_nd_donors:
        base.update({"state": "insufficient_nd_donors_for_validation",
                     "required_nd_donors": args.min_nd_donors})
        (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
        print(base["state"], len(nd_donors), "ND donors")
        return 0
    if group_gate["blocking_missing_groups"]:
        base.update({"state": "insufficient_shared_nd_feature_groups",
                     "required_groups": BIOLOGICAL_GROUPS,
                     "missing_groups": missing_groups,
                     "blocking_missing_groups": group_gate["blocking_missing_groups"],
                     "reason": "A required measured biological group other than the explicitly nonblocking stromal group has no shared primary feature across ND donors."})
        (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
        print(base["state"], "missing", ", ".join(missing_groups))
        return 0

    if not features:
        base.update({"state": "no_shared_features"})
        (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
        return 0
    max_training_rows = min(len(nd[nd["donor_id"].astype(str) != d]) for d in nd_donors)
    rank_max = args.rank_max or min(len(features), max(1, max_training_rows - 1))
    print(f"ND rank validation: {len(nd):,} admitted objects / {len(nd_donors)} donors; {MASK_REPLICATES} masks; {rank_max} candidate ranks; workers={args.workers}", flush=True)
    if args.rank_validation_from:
        cached = json.loads(args.rank_validation_from.read_text())
        expected = {"composition_sha256": base["composition_sha256"],
                    "nd_reference_input_sha256": base["nd_reference_input_sha256"],
                    "selected_features": features, "nd_donors": nd_donors,
                    "nd_object_count": len(nd)}
        for field, value in expected.items():
            if cached.get(field) != value:
                raise ValueError(f"cached rank validation has different {field}")
        if "rank_mask_protocol" in cached and cached["rank_mask_protocol"] != base["rank_mask_protocol"]:
            raise ValueError("cached rank validation used a different mask protocol")
        rank_cv = cached["rank_validation"]
        if [row["rank"] for row in rank_cv] != list(range(1, rank_max + 1)):
            raise ValueError("cached rank validation has different candidate ranks")
        if any(set(row["per_donor_rmse"]) - set(nd_donors) for row in rank_cv):
            raise ValueError("cached rank validation includes donors outside the ND reference")
        rank = select_rank_one_se(pd.DataFrame(rank_cv), MIN_TRAJECTORY_RANK)
        base["rank_validation_reuse"] = {
            "source_manifest": str(args.rank_validation_from.resolve()),
            "source_manifest_sha256": sha256(args.rank_validation_from),
            "validation": "identical composition hash, ND feature-input hash, selected features, ND donors/object count, and exhaustive candidate ranks; reconstruction algorithm/mask seeds unchanged",
            "selection_recomputed": True,
        }
        print(f"Reused exact cached rank errors; one-SE selection includes rank1: {rank}", flush=True)
    else:
        rank, rank_cv = masked_reconstruction_cv(nd.reset_index(drop=True), features, rank_max, args.workers)
    if rank is None:
        base.update({"state": "rank_validation_unavailable", "rank_validation": rank_cv})
        (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
        return 0
    model = fit_mfa(nd.reset_index(drop=True), features, rank)
    print(f"ND rank selected {rank}; graph validation: {len(nd_donors)} heldout donors x {OBJECT_SUBSAMPLES} subsamples x feasible declared configurations", flush=True)
    settings, graph_cv = choose_graph_settings(nd.reset_index(drop=True), features, rank, model, args.workers)
    if settings is None:
        base.update({"state": "graph_validation_unavailable", "representation_rank": rank,
                     "rank_validation": rank_cv, "graph_validation": graph_cv})
        (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
        return 0

    reference_rows = nd.reset_index(drop=True)
    reference_projection = project_mfa(reference_rows, model, rank)
    center = np.nanmean(reference_projection["scores"], axis=0)
    # Model is centered on the ND training set; store measured score mean to make numeric residuals explicit.
    model["reference_score_center"] = center.tolist()
    base["state"] = ("fitted_reduced_feature_set" if group_gate["nonblocking_missing_groups"]
                     else "fitted")
    base["selected_features"] = features
    base["feature_groups"] = groups
    base["unavailable_feature_groups"] = group_gate["nonblocking_missing_groups"]
    base["mfa_group_weights"] = model["group_weights"]
    contributions = {}
    for group in groups:
        idx = [i for i, feature in enumerate(features) if group_of(feature) == group]
        contributions[group] = (model["loadings"][idx, :rank] ** 2).sum(axis=0).tolist()
    total_contribution = np.sum(np.asarray(list(contributions.values())), axis=0)
    base["mfa_group_contributions_by_component"] = {
        group: (np.asarray(values) / total_contribution).tolist()
        for group, values in contributions.items()}
    base["mfa_group_first_singular_value"] = model["group_first_singular"]
    base["mfa_loadings"] = {f: model["loadings"][i, :rank].tolist()
                            for i, f in enumerate(features)}
    base["representation_rank"] = rank
    base["rank_validation"] = rank_cv
    base["graph_settings"] = settings
    base["graph_validation"] = graph_cv
    base["reference_model"] = {
        "training_donors": nd_donors, "means_transformed_features":
            {f: float(v) for f, v in zip(features, model["means"])},
        "sd_transformed_features": {f: float(v) for f, v in zip(features, model["sds"])},
        "equal_donor_standardization": True,
        "missing_entry_strategy": "iterative low-rank completion initialized at equal-donor standardized ND mean; original missingness is retained and never interpreted as a measured zero",
        "reference_score_center": center.tolist(),
        "composition_sha256": sha256(source) if source.is_file() else None,
    }

    output_frames = []
    donor_stability = []
    print(f"Graph selected {settings}; fitting {analysis_objects.donor_id.nunique()} donor paths with {OBJECT_SUBSAMPLES} resamples each", flush=True)
    frame_columns = ["donor_id", "donor_status", "endocrine_object_id", "islet_key",
                     "original_object_class", "n_cells_in_object", "endocrine_denominator", *features]
    jobs = [(donor, frame[frame_columns].copy(), model, rank, settings, center, groups, group_gate)
            for donor, frame in analysis_objects.groupby(analysis_objects["donor_id"].astype(str), sort=True)]
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
            donor_results = list(pool.map(donor_trajectory_worker, jobs))
    else:
        donor_results = list(map(donor_trajectory_worker, jobs))
    original_totals = objects.groupby(objects["donor_id"].astype(str)).size().to_dict()
    for output, stability in donor_results:
        stability["objects_analysis_eligible"] = stability["objects_total"]
        stability["objects_total"] = original_totals[stability["donor_id"]]
        stability["objects_explicitly_excluded"] = stability["objects_total"] - stability["objects_analysis_eligible"]
        output_frames.append(output)
        donor_stability.append(stability)
    excluded = objects.loc[~analysis_eligible, frame_columns[:7]].copy()
    if len(excluded):
        excluded["component_id"] = -1
        excluded["component_size"] = 0
        excluded["component_root_object_id"] = ""
        excluded["trajectory_status"] = "excluded_explicit_cell_identity_decision"
        excluded["mfa_projection_available"] = False
        excluded["n_selected_features_observed"] = 0
        excluded["n_feature_groups_observed"] = 0
        excluded["resampling_finite_fraction"] = 0.0
        excluded["trajectory_feature_groups_used"] = ";".join(groups)
        excluded["trajectory_unavailable_feature_groups"] = ";".join(group_gate["nonblocking_missing_groups"])
        output_frames.append(excluded)
    trajectories = pd.concat(output_frames, ignore_index=True)
    order = objects[["donor_id", "endocrine_object_id"]].copy()
    trajectories = order.merge(trajectories, on=["donor_id", "endocrine_object_id"], how="left", validate="one_to_one")
    metadata = [c for c in ["composition_analysis_eligible", "composition_exclusion_reason", "excluded_cell_ids", "n_original_cells_in_object", "n_analysis_admitted_cells_in_object", "n_analysis_excluded_cells_in_object", "original_endocrine_denominator", "endocrine_ratio_availability_reason"] if c in objects]
    if metadata:
        trajectories = trajectories.merge(objects[["donor_id", "endocrine_object_id", *metadata]],
                                          on=["donor_id", "endocrine_object_id"], how="left", validate="one_to_one")
    trajectories["component_id"] = trajectories.component_id.astype("Int64")
    trajectories["component_size"] = trajectories.component_size.astype("Int64")
    if len(trajectories) != len(objects) or trajectories.trajectory_status.isna().any():
        raise ValueError("Every original object must retain an explicit trajectory result or exclusion reason")
    trajectories.to_parquet(args.out / "donor_endocrine_object_trajectories.parquet", index=False)
    pd.DataFrame(donor_stability).to_csv(args.out / "donor_trajectory_stability.csv", index=False)
    (args.out / "trajectory_manifest.json").write_text(json.dumps(base, indent=2) + "\n")
    print(f"rank {rank}; graph k={settings['n_neighbors']}, dims={settings['diffusion_dimensions']}; "
          f"{len(trajectories):,} objects across {trajectories.donor_id.nunique()} donors -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
