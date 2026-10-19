"""H5 and H6: RNA spatial endpoints under the frozen spatial_v2 protocol as amended.

Ported from the archived evaluate_spatial_RNA_H5_H6.py (2026-10-06) with the 2026-10-17 rulings:
no support cutoff, only computability. Identity is the current approved RNA annotation, read as
released (named-lineage eligible); raw marker counts are measurements, never an identity gate.
Post-Xenium IF supplies islet localization only. Sections are separate spatial frames.
"""
from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from .common import Context, binding, parallel_map, peak_rss_gib, read_json, require, threaded, write_json

POSTIF_TABLES = "results/v3/islet_immune_reassessment_20261006/tables"
FOOTPRINT_CODE = ("core/phenocycler/scripts/restain/pheno_annotations.py",
                  "core/phenocycler/scripts/restain/islet_annotations.py",
                  "core/phenocycler/phenocycler/integration/rasterize.py")
ORGANIZER_GENES = ("CCL19", "CCL21")
MEMORY_GENES = ("TCF7", "LEF1")
CYTOTOXIC_GENES = ("TOX", "GZMB", "PRF1", "NKG7")
ORGANIZER_LINEAGES = ("Stromal_Fibroblast_Stellate", "Vascular")
# The approved pLN release's T identities under their current names: T_FOXP3_detected_pattern is the
# release's spelling of the protocol's T_Treg, and Immune_T is a T call the release left unsplit.
# T_NK_unresolved is not named-lineage eligible; the NK/ILC marker-pattern label is not a T identity.
PLN_T_LINEAGES = ("Immune_T", "T_CD4", "T_CD8", "T_FOXP3_detected_pattern", "T_CD4_CD8_both",
                  "T_subset_unresolved")
# (name, minimum raw library depth, neighbourhood radius um). The approved release QC (>=10 counts,
# >=5 genes, nucleus) is primary; deeper libraries are measurement sensitivities, never typing rules.
CONFIGS = (("primary", 10, 50.0), ("depth50", 50, 50.0), ("depth100", 100, 50.0), ("depth200", 200, 50.0),
           ("radius25", 10, 25.0), ("radius75", 10, 75.0))
H6_NEIGHBOURS = 8
# Largest measured peaks on 2026-10-18: 1.4 GiB (H5) and 1.2 GiB (H6).
GIB_PER_TASK = {"h5": 2, "h6": 2}


def frozen(context: Context) -> dict:
    return read_json(context.task("protocol") / "frozen_protocol.json")


def scopes(context: Context, lane: str) -> list[dict]:
    qualified = context.task("qualify")
    table = pd.read_csv(qualified / "rna_scope_qualification.csv", dtype={"donor_id": str, "source_serial": str})
    disease = pd.read_csv(qualified / "protein_scope_qualification.csv", dtype={"donor_id": str}).set_index(
        "donor_id").disease_status
    table = table.loc[table.lane.eq(lane)]
    require(table.donor_id.isin(disease.index).all(), "An RNA donor has no disease authority")
    return [{"scope_id": row.scope_id, "donor": row.donor_id, "serial": row.source_serial,
             "section": row.section_id, "disease": disease[row.donor_id], "rows": int(row.biological_rows),
             "cells": str(qualified / row.relative_data_path)} for row in table.itertuples()]


def postif(context: Context) -> Path:
    return context.original_repo / POSTIF_TABLES


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation in GIB_PER_TASK, f"Unknown RNA endpoint {operation}")
    records = [binding(context.task("qualify") / "completion.json")]
    if operation == "h5":
        tables = postif(context)
        records.extend(threaded(binding, [tables / "postIF_islet_objects.parquet",
                                          tables / "qc_postIF_islet_manifest.json",
                                          tables / "qc_postIF_islet_coverage.csv",
                                          *(context.original_repo / p for p in FOOTPRINT_CODE)]))
    return records


def seed(*parts) -> int:
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:8], 16)


def program(frame: pd.DataFrame, genes) -> np.ndarray:
    depth = frame.raw_library_counts.to_numpy(float)
    return np.mean(np.column_stack([np.log1p(1e4 * frame[f"RNA_{g}"].to_numpy(float) / depth) for g in genes]), axis=1)


def zscore(x):
    x = np.asarray(x, float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 1e-12 else np.zeros(len(x))


def numerical_columns(frame):
    columns, names = [], []
    for name in frame:
        x = frame[name].to_numpy(float)
        if np.isfinite(x).all() and x.std() > 1e-12:
            columns.append(zscore(x))
            names.append(name)
    return (np.column_stack(columns) if columns else np.empty((len(frame), 0))), names


def weighted_center(x, blocks, w):
    out = np.array(x, float)
    one = out.ndim == 1
    out = out[:, None] if one else out
    for b in np.unique(blocks):
        k = blocks == b
        out[k] -= np.average(out[k], axis=0, weights=w[k])
    return out[:, 0] if one else out


def block_wild(X0, interaction, y, blocks, w, null_n, rng):
    """Interaction coefficient beyond the pairwise baseline; whole-islet Rademacher residual signs."""
    X = weighted_center(X0, blocks, w)
    interaction = weighted_center(interaction, blocks, w)
    y = weighted_center(y, blocks, w)
    root = np.sqrt(w)
    A = X * root[:, None]
    r_i = interaction - X @ np.linalg.lstsq(A, interaction * root, rcond=None)[0]
    r_y = y - X @ np.linalg.lstsq(A, y * root, rcond=None)[0]
    denominator = float(np.sum(w * r_i * r_i))
    if denominator < 1e-10:
        return None
    beta = float(np.sum(w * r_i * r_y) / denominator)
    unique = np.unique(blocks)
    contribution = np.array([np.sum(w[blocks == b] * r_i[blocks == b] * r_y[blocks == b]) / denominator for b in unique])
    centred = np.array([np.sum(w[blocks == b] * r_i[blocks == b] * (r_y[blocks == b] - beta * r_i[blocks == b])) / denominator
                        for b in unique])
    exact = len(unique) <= 12
    if exact:
        numbers = np.arange(2 ** len(unique), dtype=np.uint32)
        signs = 2 * ((numbers[:, None] >> np.arange(len(unique))) & 1).astype(float) - 1
    else:
        signs = rng.choice([-1.0, 1.0], size=(null_n, len(unique)))
    null = signs @ contribution
    interval = beta + signs @ centred
    p = float(np.mean(null >= beta - 1e-12)) if exact else float((1 + (null >= beta).sum()) / (len(null) + 1))
    low, high = np.quantile(interval, [0.025, 0.975])
    return {"donor_effect": beta, "within_donor_p_positive": p, "CI95_low": float(low), "CI95_high": float(high),
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)) if len(null) > 1 else np.nan,
            "n_null": len(null), "exact_block_sign_enumeration": exact, "n_islet_blocks": len(unique),
            "baseline_rank_after_islet_centering": int(np.linalg.matrix_rank(A)),
            "interaction_residual_sum_squares": denominator}


# ---------------------------------------------------------------- H5

def footprint(xy, original_repo: str):
    """The section's tissue silhouette from its admitted cell centroids (v3 pheno_annotations)."""
    for path in (Path(original_repo) / "core/phenocycler/scripts/restain", Path(original_repo) / "core/phenocycler"):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    import pheno_annotations
    require(Path(pheno_annotations.__file__).resolve() == (Path(original_repo) / FOOTPRINT_CODE[0]).resolve(),
            "Executing a foreign footprint helper")
    geometry, _ = pheno_annotations.xenium_tissue_outline(xy)
    return geometry


def h5_config(cells, inside, islet, outline, config, null_n, rng):
    name, depth, radius = config
    out = {"configuration": name, "depth_threshold": depth, "neighbour_radius_um": radius,
           "donor_effect": np.nan, "status": "unavailable", "reason": ""}
    deep = cells.raw_library_counts.ge(depth).to_numpy()
    lineage = cells.candidate_lineage.astype(str).to_numpy()
    named = cells.handoff_named_lineage_use_eligible.fillna(False).to_numpy(bool)
    organizer = deep & named & inside & np.isin(lineage, ORGANIZER_LINEAGES)
    b_cell = deep & named & (lineage == "Immune_B_Plasma")
    t_cell = deep & named & (lineage == "Immune_T")
    out.update(n_organizers=int(organizer.sum()), n_B=int(b_cell.sum()), n_T=int(t_cell.sum()))
    if not organizer.any():
        return dict(out, reason="no released stromal/vascular organizer cell inside an IF islet"), None
    if not b_cell.any():
        return dict(out, reason="no released B/plasma cell in the section"), None
    if not t_cell.any():
        return dict(out, reason="no released T cell in the section"), None
    xy = cells[["x_um", "y_um"]].to_numpy(float)
    org = cells.loc[organizer]
    oxy = xy[organizer]
    b_index = np.flatnonzero(b_cell)
    b_ltb = program(cells.iloc[b_index], ["LTB"])
    b_near = cKDTree(xy[b_cell]).query_ball_point(oxy, radius)
    t_near = cKDTree(xy[t_cell]).query_ball_point(oxy, radius, return_length=True)
    everyone = cKDTree(xy[deep])
    degree = everyone.query_ball_point(oxy, radius, return_length=True) - 1
    wide = everyone.query_ball_point(oxy, 100.0, return_length=True) - 1
    area = shapely.area(shapely.intersection(outline, shapely.buffer(shapely.points(oxy), 100.0, quad_segs=8)))
    density = np.divide(wide, area, out=np.full(len(wide), np.nan), where=area > 0)
    b_count = np.array([len(v) for v in b_near])
    exposure = np.array([b_ltb[v].sum() for v in b_near]) / np.maximum(degree, 1)
    response = t_near / np.maximum(degree, 1)
    score = program(org, ORGANIZER_GENES)
    out.update(n_B_neighbour_context=int(len(np.unique(np.concatenate([np.asarray(v, int) for v in b_near])))),
               n_T_neighbour_organizers=int((t_near > 0).sum()), n_islet_blocks_with_organizer=int(len(np.unique(islet[organizer]))))
    for label, values in (("organizer CCL19/CCL21 state", score), ("neighbouring B-cell LTB exposure", exposure),
                          ("T adjacency", response)):
        if np.std(values) < 1e-12:
            return dict(out, reason=f"{label} does not vary"), None
    frame = pd.DataFrame({
        "CCL_state": zscore(score), "B_LTB_exposure": zscore(exposure), "B_neighbour_count": np.log1p(b_count),
        "local_degree": np.log1p(np.maximum(degree, 0)), "local_density_100um": np.log1p(density * 1e6),
        "library_depth": np.log1p(org.raw_library_counts.to_numpy(float)),
        "nucleus_area": np.log1p(org.nucleus_area_um2.to_numpy(float)), "cell_area": np.log1p(org.cell_area_um2.to_numpy(float)),
        "organizer_vascular": (org.candidate_lineage.astype(str) == "Vascular").to_numpy(float)})
    controls = [f for f in ("genomic_control_counts", "control_probe_counts", "control_codeword_counts")
                if f in org and np.isfinite(org[f].to_numpy(float)).all()]
    if "genomic_control_counts" not in controls:
        return dict(out, reason="genomic DNA control counts are not measured"), None
    for field in controls:
        frame[f"{field}_per1000"] = org[field].to_numpy(float) / org.raw_library_counts.to_numpy(float) * 1000
    fit = np.isfinite(frame.to_numpy(float)).all(axis=1)
    X, names = numerical_columns(frame.loc[fit])
    blocks = islet[organizer][fit]
    interaction = (frame.CCL_state * frame.B_LTB_exposure).to_numpy()[fit]
    result = block_wild(X, interaction, response[fit], blocks, np.maximum(degree[fit], 1), null_n, rng)
    if result is None:
        return dict(out, reason="the interaction is not identifiable beyond the pairwise baseline"), None
    out.update(result, status="available", n_organizers_model=int(fit.sum()), baseline_columns=",".join(names))
    nodes = pd.DataFrame({"cell_id": org.cell_id.to_numpy()[fit], "islet_id": blocks, "x_um": oxy[fit, 0],
                          "y_um": oxy[fit, 1], "CCL_state": score[fit], "B_LTB_exposure": exposure[fit],
                          "T_neighbour_fraction": response[fit], "interaction_predictor": interaction})
    return out, nodes


def h5_task(task: dict) -> dict:
    started = time.perf_counter()
    base = {"hypothesis": "H5", "scope_id": task["scope_id"], "donor_id": task["donor"],
            "disease_status": task["disease"], "source_serial": task["serial"]}
    if not task["masks"]:
        rows = [dict(base, configuration=c[0], depth_threshold=c[1], neighbour_radius_um=c[2], donor_effect=np.nan,
                     status="unavailable", reason=task["mask_reason"]) for c in CONFIGS]
        return {"rows": rows, "resource": {"scope_id": task["scope_id"], "peak_rss_gib": peak_rss_gib(),
                                           "seconds": time.perf_counter() - started}}
    cells = pd.read_parquet(task["cells"])
    require(cells.serial_qualified_cell_key.is_unique, f"Duplicate RNA keys: {task['scope_id']}")
    xy = cells[["x_um", "y_um"]].to_numpy(float)
    masks = pd.read_parquet(task["masks"])
    masks = masks.loc[masks.donor_id.astype(str).eq(task["donor"]) & masks.xenium_serial.astype(str).eq(task["serial"])]
    polygons = shapely.from_wkb(masks.geometry_wkb.to_numpy())
    hits = shapely.STRtree(polygons).query(shapely.points(xy), predicate="covered_by")
    membership = np.bincount(hits[0], minlength=len(xy))
    islet = np.full(len(xy), "", dtype=object)
    single = membership[hits[0]] == 1
    islet[hits[0][single]] = masks.islet_id.to_numpy()[hits[1][single]]
    inside = membership == 1
    outline = footprint(xy, task["original_repo"])
    rows, nodes = [], None
    for config in CONFIGS:
        rng = np.random.default_rng(seed(task["seed"], task["scope_id"], "H5", config[0]))
        row, example = h5_config(cells, inside, islet.astype(str), outline, config, task["null"], rng)
        rows.append(dict(base, n_IF_islets=len(masks), cells_in_more_than_one_mask=int((membership > 1).sum()),
                         cells_inside_one_mask=int(inside.sum()), **row))
        if config[0] == "primary" and example is not None:
            nodes = example.assign(donor_id=task["donor"], scope_id=task["scope_id"])
    if nodes is not None:
        nodes.to_parquet(Path(task["stage"]) / "organizer_nodes" / f"donor_id={task['donor']}.parquet", index=False)
    print(f"H5 {task['donor']}: " + "; ".join(f"{r['configuration']} {r['status']} {r.get('reason', '')}" for r in rows), flush=True)
    return {"rows": rows, "resource": {"scope_id": task["scope_id"], "rows": task["rows"], "peak_rss_gib": peak_rss_gib(),
                                       "seconds": time.perf_counter() - started}}


def h5(context: Context, stage: Path) -> dict:
    statistics = frozen(context)["statistics"]
    tables = postif(context)
    coverage = pd.read_csv(tables / "qc_postIF_islet_coverage.csv", dtype={"donor_id": str, "xenium_serial": str})
    (stage / "organizer_nodes").mkdir()
    tasks = []
    for scope in scopes(context, "pancreas_RNA"):
        record = coverage.loc[coverage.donor_id.eq(scope["donor"]) & coverage.xenium_serial.eq(scope["serial"])]
        available = bool(record.mask_available.any())
        reason = "" if available else ("no exact-serial post-IF islet mask: " + "; ".join(
            record.exclusion_reason.dropna().astype(str)) if len(record) else "no post-IF acquisition recorded for this section")
        tasks.append(dict(scope, stage=str(stage), masks=str(tables / "postIF_islet_objects.parquet") if available else "",
                          mask_reason=reason, original_repo=str(context.original_repo),
                          seed=statistics["seed"], null=statistics["conditional_null_replicates"]))
    results = parallel_map(h5_task, tasks, GIB_PER_TASK["h5"], "h5", cost=[t["rows"] * bool(t["masks"]) for t in tasks])
    rows = pd.DataFrame([row for result in results for row in result["rows"]])
    rows.to_csv(stage / "donor_endpoints.csv", index=False)
    pd.DataFrame([r["resource"] for r in results]).to_csv(stage / "resource_usage.csv", index=False)
    primary = rows.loc[rows.configuration.eq("primary")]
    return {"hypothesis": "H5", "donors": len(primary), "donors_available": int(primary.status.eq("available").sum()),
            "unavailable": {r.donor_id: r.reason for r in primary.itertuples() if r.status != "available"}}


# ---------------------------------------------------------------- H6

def bins(x, n):
    return np.minimum((rankdata(x, method="average") - 1) * n / max(len(x), 1), n - 1).astype(int)


def residual_state(score, covariates, labels):
    X, names = numerical_columns(covariates)
    dummies = pd.get_dummies(pd.Series(labels).astype(str), drop_first=True, dtype=float).to_numpy()
    X = np.column_stack([np.ones(len(score)), X, dummies])
    residual = score - X @ np.linalg.lstsq(X, score, rcond=None)[0]
    return zscore(residual), int(np.linalg.matrix_rank(X)), names


def h6_config(cells, config, null_n, rng):
    name, depth, radius = config
    out = {"configuration": name, "depth_threshold": depth, "neighbour_radius_um": radius,
           "architecture_excess": np.nan, "status": "unavailable", "reason": ""}
    lineage = cells.candidate_lineage.astype(str).to_numpy()
    keep = (np.isin(lineage, PLN_T_LINEAGES) & cells.handoff_named_lineage_use_eligible.fillna(False).to_numpy(bool)
            & cells.raw_library_counts.ge(depth).to_numpy())
    z = cells.loc[keep].reset_index(drop=True)
    n = len(z)
    out["n_T_cells"] = n
    if n < 3:
        return dict(out, reason="fewer than three released T cells: no neighbour contrast"), None
    if "genomic_control_counts" not in z or not np.isfinite(z.genomic_control_counts.to_numpy(float)).all():
        return dict(out, reason="genomic DNA control counts are not measured"), None
    memory, cyto = program(z, MEMORY_GENES), program(z, CYTOTOXIC_GENES)
    out.update(mean_memory_state=float(memory.mean()), mean_cytotoxic_state=float(cyto.mean()),
               mean_library_depth=float(z.raw_library_counts.mean()))
    library = z.raw_library_counts.to_numpy(float)
    covariates = pd.DataFrame({"log_depth": np.log1p(library), "log_nucleus_area": np.log1p(z.nucleus_area_um2.to_numpy(float)),
                               "log_cell_area": np.log1p(z.cell_area_um2.to_numpy(float)),
                               "gDNA_per1000": z.genomic_control_counts.to_numpy(float) / library * 1000})
    for field in ("control_probe_counts", "control_codeword_counts"):
        if field in z and np.isfinite(z[field].to_numpy(float)).all():
            covariates[f"{field}_per1000"] = z[field].to_numpy(float) / library * 1000
    labels = z.candidate_lineage.astype(str).to_numpy()
    contrast, rank, controls = residual_state(memory - cyto, covariates, labels)
    out.update(state_nuisance_rank=rank, state_control_fields=",".join(controls))
    if contrast.std() < 1e-10:
        return dict(out, reason="the memory-minus-cytotoxic contrast has no residual variance"), None
    xy = z[["x_um", "y_um"]].to_numpy(float)
    k = min(H6_NEIGHBOURS, n - 1)
    distance, index = cKDTree(xy).query(xy, k=k + 1, distance_upper_bound=radius)
    ii = np.repeat(np.arange(n), k)
    jj = index[:, 1:].ravel()
    valid = (jj < n) & (distance[:, 1:].ravel() <= radius) & (jj != ii)
    ii, jj = ii[valid], jj[valid]
    out.update(n_directed_edges=len(ii), n_cells_with_T_neighbour=int(np.unique(ii).size))
    if not len(ii):
        return dict(out, reason="no T cell has a T neighbour within the radius"), None
    code = pd.factorize(labels, sort=True)[0]
    strata = (code * 10000 + bins(library, 5) * 1000 + bins(z.nucleus_area_um2.to_numpy(float), 4) * 100
              + (z.genomic_control_counts.to_numpy() > 0).astype(int))
    _, inverse, counts = np.unique(strata, return_inverse=True, return_counts=True)
    out["fraction_cells_exchangeable"] = float(np.mean(counts[inverse] >= 2))
    observed = float(np.mean(contrast[ii] * contrast[jj]))
    order = np.argsort(strata, kind="stable")
    null = np.empty(null_n)
    permuted = np.empty(n)
    for r in range(null_n):
        permuted[order] = contrast[np.lexsort((rng.random(n), strata))]
        null[r] = np.mean(permuted[ii] * permuted[jj])
    if null.std() < 1e-12:
        return dict(out, reason="no stratum has two exchangeable cells: the label null is the identity"), None
    out.update(status="available", observed_assortativity=observed, null_mean_assortativity=float(null.mean()),
               null_sd_assortativity=float(null.std(ddof=1)), architecture_excess=observed - float(null.mean()),
               within_section_p_positive=float((1 + (null >= observed).sum()) / (null_n + 1)), n_null=null_n)
    example = z[["cell_id", "x_um", "y_um", "candidate_lineage"]].assign(memory_state=memory, cytotoxic_state=cyto,
                                                                          residual_contrast=contrast)
    return out, example


def h6_task(task: dict) -> dict:
    started = time.perf_counter()
    cells = pd.read_parquet(task["cells"])
    require(cells.serial_qualified_cell_key.is_unique, f"Duplicate RNA keys: {task['scope_id']}")
    rows = []
    for config in CONFIGS:
        rng = np.random.default_rng(seed(task["seed"], task["scope_id"], "H6", config[0]))
        row, example = h6_config(cells, config, task["null"], rng)
        rows.append({"hypothesis": "H6", "scope_id": task["scope_id"], "donor_id": task["donor"],
                     "section_id": task["section"], "disease_status": task["disease"], **row})
        if config[0] == "primary" and example is not None:
            example.to_parquet(Path(task["stage"]) / "pln_states" / f"{task['section']}.parquet", index=False)
    print(f"H6 {task['scope_id']}: " + "; ".join(f"{r['configuration']} {r['status']}" for r in rows), flush=True)
    return {"rows": rows, "resource": {"scope_id": task["scope_id"], "rows": task["rows"], "peak_rss_gib": peak_rss_gib(),
                                       "seconds": time.perf_counter() - started}}


def h6(context: Context, stage: Path) -> dict:
    statistics = frozen(context)["statistics"]
    (stage / "pln_states").mkdir()
    tasks = [dict(scope, stage=str(stage), seed=statistics["seed"], null=statistics["conditional_null_replicates"])
             for scope in scopes(context, "pLN_RNA")]
    results = parallel_map(h6_task, tasks, GIB_PER_TASK["h6"], "h6", cost=[t["rows"] for t in tasks])
    sections = pd.DataFrame([row for result in results for row in result["rows"]])
    sections.to_csv(stage / "pln_section_architecture.csv", index=False)
    donors = []
    for (donor, configuration), z in sections.groupby(["donor_id", "configuration"], sort=True):
        usable = z.loc[z.status.eq("available")]
        row = {"donor_id": donor, "disease_status": z.disease_status.iloc[0], "configuration": configuration,
               "n_sections": len(z), "n_sections_available": len(usable),
               "section_aggregation": "unweighted mean of section metrics; coordinates never merged",
               "architecture_excess": np.nan, "status": "unavailable",
               "reason": "; ".join(sorted(set(z.reason))) if not len(usable) else ""}
        if len(usable):
            weights = usable.n_T_cells.to_numpy(float)
            row.update(status="available", architecture_excess=float(usable.architecture_excess.mean()),
                       mean_memory_state=float(np.average(usable.mean_memory_state, weights=weights)),
                       mean_cytotoxic_state=float(np.average(usable.mean_cytotoxic_state, weights=weights)),
                       n_T_cells=int(usable.n_T_cells.sum()))
        donors.append(row)
    donors = pd.DataFrame(donors)
    donors.to_csv(stage / "pln_donor_architecture.csv", index=False)
    h3 = pd.read_csv(context.task("h3") / "donor_endpoints.csv", dtype={"donor_id": str})
    h3 = h3.loc[h3.configuration.eq("primary"), ["donor_id", "disease_status", "donor_effect", "status", "reason"]].rename(
        columns={"disease_status": "H3_disease_status", "donor_effect": "H3", "status": "H3_status", "reason": "H3_reason"})
    paired = donors.loc[donors.configuration.eq("primary")].merge(h3, on="donor_id", how="outer", validate="one_to_one")
    no_pln = paired.status.isna()
    paired.loc[no_pln, ["status", "reason"]] = ["unavailable", "no pLN RNA section for this donor"]
    paired["disease_status"] = paired.disease_status.fillna(paired.H3_disease_status)
    require(paired.disease_status.notna().all(), "A paired donor has no disease status")
    paired = paired.drop(columns="H3_disease_status")
    paired["paired_status"] = np.where(paired.status.eq("available") & paired.H3_status.eq("available"),
                                       "available", "unavailable")
    paired.to_csv(stage / "paired_donors.csv", index=False)
    pd.DataFrame([r["resource"] for r in results]).to_csv(stage / "resource_usage.csv", index=False)
    write_json(stage / "pln_T_lineages.json", {"included": list(PLN_T_LINEAGES),
                                                "excluded": ["T_NK_unresolved", "T_NK_subtype_unresolved_NK_ILC_marker_pattern"]})
    return {"hypothesis": "H6", "pLN_sections": len(tasks),
            "pLN_donors_available": int(donors.loc[donors.configuration.eq("primary")].status.eq("available").sum()),
            "paired_donors_available": int(paired.paired_status.eq("available").sum())}
