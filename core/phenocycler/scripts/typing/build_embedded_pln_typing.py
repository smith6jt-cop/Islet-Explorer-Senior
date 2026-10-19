#!/usr/bin/env python3
"""Type the current embedded pLN cell populations with the dedicated pLN method.

Current v3 pancreas partitions intentionally omit cells whose Parent is Annotation (pLN).
This producer reads those cells from the current full-cell QuPath measurements, uses the
authoritative annotation/block sidecar, and runs the established within-section PhenoCycler pLN
cluster-and-name chain separately for each pLN annotation. It never imports a pancreas level or
uses pancreas-vs-node distance.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/pln"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/celltype"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import alignment_block as AB  # noqa: E402
import cluster_pln as CP  # noqa: E402
import build_pln_immune_pheno as BPI  # noqa: E402
from phenocycler.cells_parquet import sanitize  # noqa: E402


DEFAULT_RELEASE = REPO / "results/v3/release_20261003"
DEFAULT_BLOCK_MANIFEST = REPO / "data/development/celltype/alignment_block/manifest.json"
DEFAULT_INPUT_MANIFEST = DEFAULT_RELEASE / "current_input_manifest.json"
DEFAULT_QUPATH_MANIFEST = REPO / "core/phenocycler/data/qupath_manifest.json"
DEFAULT_GEOMETRY_QC = REPO / "core/phenocycler/data/geometry_qc"
DEFAULT_DEDICATED_CELLS = REPO / "core/phenocycler/data/integration/cells_pln"
DEFAULT_DEDICATED_CLUSTERS = REPO / "data/development/typing/pln/clusters/phenocycler"
SEED = 20260917


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalized_name(value: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "", value).casefold()


def current_channel_columns(names: list[str], dedicated_channels: list[str]) -> dict[str, str]:
    """Map raw QuPath Cell means onto this project's canonical PhenoCycler channel names."""
    dedicated_by_key: dict[str, list[str]] = {}
    for marker in dedicated_channels:
        dedicated_by_key.setdefault(normalized_name(marker), []).append(marker)
    nuclei = set(names)
    mapping: dict[str, str] = {}
    seen: dict[str, str] = {}
    for cell_col in names:
        match = re.fullmatch(r"Cell: (.+): Mean", cell_col)
        if not match:
            continue
        raw_marker = match.group(1)
        nucleus_col = f"Nucleus: {raw_marker}: Mean"
        if nucleus_col not in nuclei:
            continue
        key = normalized_name(sanitize(raw_marker))
        options = dedicated_by_key.get(key, [])
        if len(options) > 1:
            raise ValueError(f"ambiguous channel mapping for {raw_marker!r}: {options}")
        marker = options[0] if options else sanitize(raw_marker)
        if marker in seen:
            raise ValueError(f"two current QuPath channels map to {marker!r}: {seen[marker]!r} and {cell_col!r}")
        seen[marker] = cell_col
        mapping[marker] = cell_col
    return mapping


def load_marker_dictionary(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path)
    if frame.duplicated(["donor_id", "object_id"]).any():
        raise ValueError(f"duplicate geometry QC keys in {path}")
    return frame.set_index("object_id").to_dict("index")


def marker_tests(st: dict, cluster_ids: list[str], markers: list[str]) -> pd.DataFrame:
    rows = []
    for i, cid in enumerate(cluster_ids):
        for j, marker in enumerate(markers):
            rows.append({
                "cluster_id": cid,
                "marker": marker,
                "auc_vs_other_cells_same_node": float(st["auc"][i, j]),
                "offset_adjusted_auc": float(st["rel"][i, j]),
                "rank_test_p": float(st["p"][i, j]),
                "bonferroni_alpha": float(st["alpha"]),
                "enriched_by_native_rule": bool(st["enriched"][i, j]),
                "median_over_section_median": float(st["magnitude"][i, j]),
                "n_finite": int(st["n_finite"][i, j]),
            })
    return pd.DataFrame(rows)


def build_annotation(
    donor: str,
    annotation_id: str,
    block: str,
    frame: pd.DataFrame,
    markers: list[str],
    channels: list[str],
    out: Path,
    *,
    subsample: int,
    resolution: float,
    seed_base: int,
) -> tuple[pd.DataFrame, dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Cluster and call one physical pLN annotation independently."""
    feats = [m for m in markers if m not in CP.PHENO_FEATURE_EXCLUDE]
    values = frame[feats].to_numpy(dtype=np.float64)
    all_null = np.isnan(values).all(axis=0)
    unavailable = [m for m, bad in zip(feats, all_null) if bad]
    kept_features = [m for m, bad in zip(feats, all_null) if not bad]
    values = values[:, ~all_null]
    row_ok = np.isfinite(values).all(axis=1)
    n_partial = int((~row_ok).sum())
    typed = frame.loc[row_ok].reset_index(drop=True)
    values = values[row_ok]
    if len(typed) < 2 or values.shape[1] < 2:
        unavailable_result = {
            "state": "unavailable",
            "reason": "fewer than two analysis-eligible cells with at least two measured pLN features",
            "n_analysis_eligible": int(len(frame)),
            "n_cells_partial_measurement": n_partial,
            "markers_unavailable": unavailable,
        }
        assignments = pd.DataFrame({
            "object_id": frame["object_id"].astype(str),
            "donor_id": donor,
            "pln_annotation_id": annotation_id,
            "block_relationship": block,
            "typing_status": "unavailable",
            "typing_reason": unavailable_result["reason"],
        })
        return assignments, unavailable_result, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    z, med, mad = CP.robust_z(np.log1p(np.clip(values, 0.0, None)))
    zero_mad = [m for m, v in zip(kept_features, mad) if v <= 0]
    z = np.clip(z, -5.0, 10.0)
    seed = CP.stable_seed("pln-cluster", seed_base, "phenocycler", f"{donor}:{annotation_id}")
    labels, cluster_meta = CP.leiden_on_subsample(
        z, n_pcs=20, subsample=subsample, resolution=resolution, seed=seed
    )
    cluster_ids = [f"{donor}_embedded_{annotation_id[:8]}_c{int(v):02d}"
                   for v in sorted(set(labels.astype(int).tolist()))]
    labels = np.array([f"{donor}_embedded_{annotation_id[:8]}_c{int(v):02d}" for v in labels])
    codes = pd.Categorical(labels, categories=cluster_ids).codes.astype(np.int64)
    sizes = np.bincount(codes, minlength=len(cluster_ids))

    # The same established pLN cell typing rules are applied on this node's own measurements.
    panel_markers = sorted(set(channels) | set(markers))
    members = BPI.subpanel_membership(panel_markers)
    named = BPI.naming_markers(members)
    nk_all, nk_any, nk_none = BPI.nk_licence()
    decision_markers = sorted(
        {m for lin in BPI.IMMUNE_LINEAGES for m in BPI.PRIMARY_IMMUNE[lin]}
        | set(BPI.STRUCTURAL)
        | set(BPI.SUPPORTING_IMMUNE)
        | {m for group in named.values() for m in group}
        | set(nk_all)
        | {m for group in nk_any for m in group}
        | set(nk_none)
    )
    evidence_markers = sorted(set(panel_markers) | set(decision_markers))
    marker_frame = typed.copy()
    for marker in evidence_markers:
        if marker not in marker_frame:
            marker_frame[marker] = np.nan
    marker_frame = marker_frame[evidence_markers]
    st = BPI.standout(marker_frame, codes, len(cluster_ids), evidence_markers,
                      offset_over=[m for m in channels if m in evidence_markers])
    st["magnitude"] = BPI.magnitude(marker_frame, codes, len(cluster_ids), evidence_markers)
    dyn_markers = [m for m in decision_markers if m in marker_frame]
    dyn = BPI.dynamic_range(marker_frame[dyn_markers], codes, len(cluster_ids), dyn_markers)
    ev = BPI.cell_evidence(st, members, len(typed))
    calls = BPI.call_clusters(donor, cluster_ids, sizes, st, members)
    subpanel_argmax = BPI.per_cell_subpanel(calls, cluster_ids, codes, st)
    het = BPI.heterogeneity(calls, cluster_ids, codes, ev, subpanel_argmax)
    resolution_table, tnk = BPI.resolve_t_nk(donor, calls, cluster_ids, codes, st, channels)

    by_cluster = calls.set_index("cluster_id")
    take = lambda col: by_cluster[col].reindex(cluster_ids).to_numpy()[codes]  # noqa: E731
    assignments = pd.DataFrame({
        "object_id": typed["object_id"].astype(str).to_numpy(),
        "donor_id": donor,
        "pln_annotation_id": annotation_id,
        "block_relationship": block,
        "cluster_id": labels,
        "typing_status": "typed",
        "is_immune": pd.array(take("is_immune"), dtype="boolean"),
        "immune_class": pd.Categorical(take("immune_class").astype(str)),
        "immune_subtype": pd.Categorical(take("immune_subtype").astype(str)),
        "t_nk_class": pd.Categorical(tnk["t_nk_class"].astype(str)),
        "t_nk_resolution_mode": pd.Categorical(tnk["t_nk_resolution_mode"].astype(str)),
        "stage1_rule": pd.Categorical(take("stage1_rule").astype(str)),
        "stage2_reason": pd.Categorical(take("stage2_reason").astype(str)),
        "primary_evidence_marker": pd.Categorical(take("top_immune_marker").astype(str)),
        "immune_marker_enriched": take("immune_marker_enriched").astype(bool),
        "structural_marker_enriched": take("structural_marker_enriched").astype(bool),
        "cell_argmax_subpanel": pd.Categorical(subpanel_argmax.astype(str)),
        "cell_immune_magnitude": ev["immune_magnitude"].astype(float),
        "cell_structural_magnitude": ev["structural_magnitude"].astype(float),
        "cell_agrees_with_cluster_class": pd.array(
            np.where(subpanel_argmax == "", None,
                     subpanel_argmax == take("immune_class").astype(str)), dtype="boolean"),
        "typing_reason": "inherited from pLN annotation-specific cluster",
    })
    cluster_calls = calls.copy()
    cluster_calls["pln_annotation_id"] = annotation_id
    cluster_calls["block_relationship"] = block
    cluster_calls = cluster_calls.merge(het, on=["donor_id", "cluster_id"], suffixes=("", "_heterogeneity"),
                                        validate="one_to_one")
    cluster_calls["knn_holdout_agreement"] = cluster_meta["knn_holdout_agreement"]
    cluster_calls["pca_variance_explained"] = cluster_meta["pca_variance_explained"]
    cluster_calls["markers_unavailable"] = ";".join(unavailable)
    cluster_calls["markers_zero_mad"] = ";".join(zero_mad)
    profile_rows = []
    for i, cid in enumerate(cluster_ids):
        mask = codes == i
        row = {"donor_id": donor, "pln_annotation_id": annotation_id,
               "block_relationship": block, "cluster_id": cid, "n_cells": int(mask.sum()),
               "fraction_of_annotation": float(mask.mean())}
        for j, marker in enumerate(kept_features):
            row[f"z__{marker}"] = float(z[mask, j].mean())
            row[f"raw__{marker}"] = float(values[mask, j].mean())
        profile_rows.append(row)
    profiles = pd.DataFrame(profile_rows)
    tests = marker_tests(st, cluster_ids, evidence_markers)
    tests["donor_id"] = donor
    tests["pln_annotation_id"] = annotation_id
    tests["block_relationship"] = block
    result = {
        "state": "typed",
        "n_analysis_eligible": int(len(frame)),
        "n_clustered": int(len(typed)),
        "n_cells_partial_measurement": n_partial,
        "n_clusters": int(len(cluster_ids)),
        "markers_used": kept_features,
        "markers_unavailable": unavailable,
        "markers_zero_mad": zero_mad,
        "cluster_algorithm": cluster_meta,
        "n_cells_called_immune": int((assignments.is_immune == True).sum()),  # noqa: E712
        "n_cells_called_nonimmune": int((assignments.is_immune == False).sum()),  # noqa: E712
        "n_cells_unresolved": int(assignments.is_immune.isna().sum()),
        "immune_class_counts": assignments.loc[assignments.is_immune == True, "immune_class"]
                                       .astype(str).value_counts().to_dict(),  # noqa: E712
        "t_nk_summary": BPI.t_nk_summary(donor, resolution_table,
                                          int(calls.loc[calls.immune_class == BPI.T_NK_CLASS, "n_cells"].sum())),
        "section_channel_dynamic_range": dyn,
    }
    # Keep a cell row for partial measurements with a reason instead of dropping it.
    partial = frame.loc[~frame["object_id"].astype(str).isin(typed["object_id"].astype(str)),
                        ["object_id"]].copy()
    if len(partial):
        partial_assignments = pd.DataFrame({
            "object_id": partial.object_id.astype(str), "donor_id": donor,
            "pln_annotation_id": annotation_id, "block_relationship": block,
            "typing_status": "unavailable", "typing_reason": "partial_marker_measurement",
        })
        assignments = pd.concat([assignments, partial_assignments], ignore_index=True)
    return assignments, result, cluster_calls, profiles, tests


def compare_dedicated_profiles(
    embedded_profiles: pd.DataFrame,
    dedicated_root: Path,
) -> pd.DataFrame:
    rows = []
    if embedded_profiles.empty:
        return pd.DataFrame()
    for (donor, annotation), emb in embedded_profiles.groupby(["donor_id", "pln_annotation_id"]):
        dedicated_path = dedicated_root / f"profiles/{donor}.parquet"
        if not dedicated_path.is_file():
            rows.append({"donor_id": donor, "pln_annotation_id": annotation,
                         "state": "unavailable", "reason": "no same-donor dedicated pLN cluster profile"})
            continue
        ded = pq.read_table(dedicated_path).to_pandas()
        shared = sorted(set(c.removeprefix("z__") for c in emb if c.startswith("z__"))
                        & set(c.removeprefix("z__") for c in ded if c.startswith("z__")))
        if len(shared) < 3:
            rows.append({"donor_id": donor, "pln_annotation_id": annotation,
                         "state": "unavailable", "n_shared_markers": len(shared),
                         "reason": "fewer than three shared measured protein markers"})
            continue
        best_records = []
        for _, erow in emb.iterrows():
            best = None
            for _, drow in ded.iterrows():
                x = np.array([erow.get(f"z__{m}", np.nan) for m in shared], dtype=float)
                y = np.array([drow.get(f"z__{m}", np.nan) for m in shared], dtype=float)
                ok = np.isfinite(x) & np.isfinite(y)
                if int(ok.sum()) < 3:
                    continue
                rho = float(spearmanr(x[ok], y[ok]).statistic)
                rec = {"donor_id": donor, "pln_annotation_id": annotation,
                       "embedded_cluster_id": erow.cluster_id,
                       "dedicated_cluster_id": drow.cluster_id,
                       "state": "descriptive_profile_match",
                       "n_shared_markers": int(ok.sum()),
                       "spearman_profile_correlation": rho,
                       "n_embedded_cells": int(erow.n_cells),
                       "n_dedicated_cells": int(drow.n_cells),
                       "comparison_scope": "same donor, distinct slides; cluster profiles only; no cell matching"}
                rows.append(rec)
                if best is None or (np.isfinite(rho) and rho > best["spearman_profile_correlation"]):
                    best = rec
            if best:
                best_records.append(best)
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-manifest", type=Path, default=DEFAULT_INPUT_MANIFEST)
    ap.add_argument("--block-manifest", type=Path, default=DEFAULT_BLOCK_MANIFEST)
    ap.add_argument("--qupath-manifest", type=Path, default=DEFAULT_QUPATH_MANIFEST)
    ap.add_argument("--geometry-qc", type=Path, default=DEFAULT_GEOMETRY_QC)
    ap.add_argument("--dedicated-cells", type=Path, default=DEFAULT_DEDICATED_CELLS)
    ap.add_argument("--dedicated-clusters", type=Path, default=DEFAULT_DEDICATED_CLUSTERS)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--subsample", type=int, default=60_000)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    input_manifest = json.loads(args.input_manifest.read_text())
    block_manifest = json.loads(args.block_manifest.read_text())
    qmanifest = json.loads(args.qupath_manifest.read_text())
    images = {str(r["donor_id"]): r for r in qmanifest["images"]}
    annotations_dir = Path(block_manifest["annotations"])
    sidecars = {str(r["donor_id"]): r for r in block_manifest.get("donors", [])}
    if args.donors:
        requested = {str(x) for x in args.donors}
        sidecars = {d: r for d, r in sidecars.items() if d in requested}

    out = args.out
    calls_root = out / "calls"
    for p in (calls_root, out / "cluster_calls", out / "cluster_profiles", out / "marker_tests"):
        p.mkdir(parents=True, exist_ok=True)
    dedicated_channels = BPI.panel_channels()
    all_assignments, cluster_rows, profile_rows, test_rows, node_rows = [], [], [], [], []

    for donor in sorted(sidecars):
        image_row = images.get(donor)
        if image_row is None:
            # Some old manifests spell 0112 as 112; compare numerically only for this ID lookup.
            alt = next((r for k, r in images.items() if k.lstrip("0") == donor.lstrip("0")), None)
            image_row = alt
        if image_row is None:
            raise ValueError(f"{donor}: no QuPath image row in current manifest")
        image = str(image_row["image_id"])
        measurement = Path(image_row["measurement_csv"]["path"])
        if not measurement.is_file():
            raise FileNotFoundError(f"{donor}: current full-cell measurement missing: {measurement}")
        side_path = args.block_manifest.parent / f"donor_id={donor}.parquet"
        if not side_path.is_file():
            # block manifest paths can have the legacy unpadded spelling
            alt_side = args.block_manifest.parent / f"donor_id={donor.lstrip('0')}.parquet"
            if not alt_side.is_file():
                raise FileNotFoundError(f"{donor}: embedded pLN sidecar missing")
            side_path = alt_side
        side = pq.read_table(side_path).to_pandas()
        if side.object_id.duplicated().any():
            raise ValueError(f"{donor}: duplicate object IDs in block sidecar")
        ann_map = sidecars[donor]["annotations"]
        geojson = annotations_dir / f"{AB.image_tag(image)}.geojson"
        annotations = AB.pln_annotations(geojson)
        on_image = {annotation for annotation, _ in annotations}
        if on_image != set(ann_map):
            raise ValueError(f"{donor}: current pLN annotation geometry and confirmed block manifest disagree")

        file = pq.ParquetFile(measurement)
        source_names = file.schema_arrow.names
        channels_map = current_channel_columns(source_names, dedicated_channels)
        channel_markers = sorted(channels_map)
        columns = ["Object ID", "Image", "Parent", "Centroid X µm", "Centroid Y µm",
                   "Cell: Area µm^2", "Cell: Solidity", "Nucleus: Area µm^2",
                   "Nucleus/Cell area ratio"] + [channels_map[m] for m in channel_markers]
        columns = [c for c in columns if c in source_names]
        raw = file.read(columns=columns).to_pandas()
        raw["Object ID"] = raw["Object ID"].astype(str)
        is_pln_parent = raw["Parent"].fillna("").astype(str).str.fullmatch(
            r"Annotation\s*\(\s*pLN\s*\)", case=False)
        raw = raw.loc[is_pln_parent].copy()
        if raw["Object ID"].duplicated().any():
            raise ValueError(f"{donor}: duplicate Object IDs in current full-cell QuPath source")

        # Determine every annotated pLN cell from its current outline and centroid. For already
        # admitted cells, require exact agreement with the confirmed block sidecar.
        assigned_id, found_by = AB.assign(raw["Centroid X µm"].to_numpy(float),
                                          raw["Centroid Y µm"].to_numpy(float), annotations)
        raw["pln_annotation_id"] = assigned_id
        raw["membership_found_by"] = found_by
        raw["block_relationship"] = [ann_map.get(str(a), "unassigned") if a else "unassigned"
                                     for a in assigned_id]
        current_assign = raw.set_index("Object ID")["pln_annotation_id"]
        current_block = raw.set_index("Object ID")["block_relationship"]
        side = side.set_index("object_id", drop=False)
        side_ann = current_assign.reindex(side.index)
        side_block = current_block.reindex(side.index)
        if side_ann.isna().any() or side_block.isna().any():
            raise ValueError(f"{donor}: sidecar contains IDs missing in the current full-cell export")
        if not np.array_equal(side_ann.astype(str).to_numpy(), side["pln_annotation_id"].astype(str).to_numpy()):
            raise ValueError(f"{donor}: current annotation assignment disagrees with the existing sidecar")
        if not np.array_equal(side_block.astype(str).to_numpy(), side["block"].astype(str).to_numpy()):
            raise ValueError(f"{donor}: confirmed block assignment disagrees with current annotation")

        qc_path = args.geometry_qc / f"donor_id={donor}/data.parquet"
        qc = pq.read_table(qc_path, columns=["object_id", "analysis_eligible", "analysis_reason"]).to_pandas()
        if qc.object_id.duplicated().any():
            raise ValueError(f"{donor}: duplicate IDs in current geometry QC")
        qc = qc.set_index("object_id")
        qc_aligned = qc.reindex(raw["Object ID"].astype(str))
        raw["qc_analysis_eligible"] = qc_aligned["analysis_eligible"].fillna(False).to_numpy(bool)
        raw["analysis_reason"] = qc_aligned["analysis_reason"].fillna(
            "unavailable in current geometry QC").astype(str).to_numpy()
        # The block sidecar is a geometric membership layer and can include cells later rejected
        # by morphology QC. Typing eligibility comes only from the current geometry QC table.

        for marker, col in channels_map.items():
            raw[marker] = pd.to_numeric(raw[col], errors="coerce")
        raw["donor_id"] = donor
        raw["image"] = image
        raw["object_id"] = raw["Object ID"].astype(str)
        raw["X_centroid"] = raw["Centroid X µm"]
        raw["Y_centroid"] = raw["Centroid Y µm"]

        donor_calls = []
        donor_profiles = []
        donor_cluster_calls = []
        donor_tests = []
        for annotation_id, ann_frame in raw.groupby("pln_annotation_id", dropna=False, sort=True):
            annotation_id = str(annotation_id)
            ann_frame = ann_frame.copy()
            block = ann_map.get(annotation_id, "unassigned") if annotation_id else "unassigned"
            eligible = ann_frame.loc[ann_frame.qc_analysis_eligible].copy()
            geometry_scope = "whole_slide_annotation" if donor == "6534" else "embedded_node_annotation"
            base_cols = ["object_id", "donor_id", "image", "pln_annotation_id", "block_relationship",
                         "X_centroid", "Y_centroid", "qc_analysis_eligible", "analysis_reason"]
            assignments = pd.DataFrame({
                "object_id": ann_frame.object_id,
                "donor_id": donor,
                "image": image,
                "pln_annotation_id": annotation_id,
                "block_relationship": block,
                "qc_analysis_eligible": ann_frame.qc_analysis_eligible.astype(bool),
                "analysis_reason": ann_frame.analysis_reason.astype(str),
                "typing_status": "unavailable",
                "typing_reason": np.where(ann_frame.qc_analysis_eligible,
                                           "not_clustered_or_no_measured_identity_evidence",
                                           ann_frame.analysis_reason.astype(str)),
            })
            result = {"state": "unavailable", "reason": "not run"}
            cluster_calls = pd.DataFrame()
            profiles = pd.DataFrame()
            tests = pd.DataFrame()
            if block != "unassigned" and not eligible.empty:
                value_frame = eligible[["object_id"] + channel_markers].copy()
                calls, result, cluster_calls, profiles, tests = build_annotation(
                    donor, annotation_id, block, value_frame, channel_markers, channel_markers,
                    out, subsample=args.subsample, resolution=args.resolution, seed_base=args.seed)
                call_fields = [c for c in calls if c not in {"donor_id", "pln_annotation_id",
                                                               "block_relationship"}]
                assignments = assignments.drop(columns=["typing_status", "typing_reason"]).merge(
                    calls[call_fields], on="object_id", how="left", validate="one_to_one")
                typed_mask = assignments["typing_status"].eq("typed")
                assignments.loc[assignments["typing_status"].isna(), "typing_status"] = "unavailable"
                assignments.loc[typed_mask, "typing_reason"] = assignments.loc[typed_mask, "typing_reason"].fillna(
                    "inherited from pLN annotation-specific cluster")
                missing_reason = assignments["typing_reason"].isna()
                assignments.loc[missing_reason, "typing_reason"] = assignments.loc[missing_reason,
                                                                                     "analysis_reason"]
            else:
                result = {"state": "unavailable", "reason": "unassigned annotation" if block == "unassigned"
                          else "no current analysis-eligible cells"}
            result.update({
                "donor_id": donor, "pln_annotation_id": annotation_id,
                "block_relationship": block,
                "geometry_scope": geometry_scope,
                "n_source_pLN_parent_cells": int(len(ann_frame)),
                "n_analysis_eligible": int(ann_frame.qc_analysis_eligible.sum()),
                "n_analysis_ineligible": int((~ann_frame.qc_analysis_eligible).sum()),
                "analysis_ineligible_reasons": ann_frame.loc[~ann_frame.qc_analysis_eligible,
                                                               "analysis_reason"].value_counts().to_dict(),
                "membership_methods": ann_frame.membership_found_by.value_counts().to_dict(),
                "source_measurement": str(measurement),
                "geometry_qc": str(qc_path),
                "block_sidecar": str(side_path),
            })
            dest = calls_root / f"donor_id={donor}" / f"annotation_id={annotation_id}"
            dest.mkdir(parents=True, exist_ok=True)
            assignments.to_parquet(dest / "data.parquet", index=False)
            if len(cluster_calls):
                cluster_calls.to_parquet(out / "cluster_calls" / f"donor_id={donor}__annotation={annotation_id}.parquet",
                                         index=False)
                profiles.to_parquet(out / "cluster_profiles" / f"donor_id={donor}__annotation={annotation_id}.parquet",
                                    index=False)
                # Marker tests are a separate table so profile comparisons remain easy to consume.
                tests.to_parquet(out / "marker_tests" / f"donor_id={donor}__annotation={annotation_id}.parquet",
                                 index=False)
                donor_cluster_calls.append(cluster_calls)
                donor_profiles.append(profiles)
                donor_tests.append(tests)
            donor_calls.append(assignments)
            node_rows.append(result)
            print(f"{donor} {annotation_id[:8]} {block}: {result['state']} | "+
                  f"cells {result.get('n_analysis_eligible',0):,}; typed {result.get('n_clustered',0):,}; "+
                  f"clusters {result.get('n_clusters','')}", flush=True)
        if donor_calls:
            all_assignments.extend(donor_calls)
        if donor_cluster_calls:
            cluster_rows.extend(donor_cluster_calls)
        if donor_profiles:
            profile_rows.extend(donor_profiles)
        if donor_tests:
            test_rows.extend(donor_tests)

    node_manifest = pd.DataFrame(node_rows)
    node_manifest.to_csv(out / "node_coverage.csv", index=False)
    if profile_rows:
        profiles = pd.concat(profile_rows, ignore_index=True, sort=False)
        profiles.to_parquet(out / "cluster_profiles.parquet", index=False)
        comparison = compare_dedicated_profiles(profiles, args.dedicated_clusters)
        comparison.to_csv(out / "dedicated_pLN_profile_comparison.csv", index=False)
    else:
        profiles = pd.DataFrame()
        pd.DataFrame().to_csv(out / "dedicated_pLN_profile_comparison.csv", index=False)
    if all_assignments:
        # A single cell-level partition per donor/annotation is the native release interface.
        all_cells = pd.concat(all_assignments, ignore_index=True, sort=False)
        all_cells.to_parquet(out / "embedded_pln_cell_typing.parquet", index=False)
    report = {
        "title": "Current embedded pLN cell typing",
        "input_manifest": str(args.input_manifest.resolve()),
        "block_manifest": str(args.block_manifest.resolve()),
        "block_manifest_sha256": sha256(args.block_manifest),
        "qupath_manifest": str(args.qupath_manifest.resolve()),
        "source": "current full-cell QuPath measurement exports; pLN-parented cells only",
        "methods": {
            "membership": "current pLN annotation outline and centroid; checked cell-for-cell against the confirmed block sidecar",
            "block_assignment": "same-block/different-block from pln_blocks.csv via the signed sidecar; never inferred",
            "qc": "current geometry_qc analysis_eligible; block membership is independent of QC eligibility",
            "typing": "established pLN within-node PhenoCycler robust-z PCA/Leiden plus cluster-relative enrichment, magnitude arbitration, identity-core subpanels, and T/NK resolution",
            "scale": "raw within-section measurements; no pancreas value, scale, reference, or positivity threshold transferred",
            "feature_transform": "log1p then within-node robust z (median / 1.4826*MAD), clipped [-5, 10]; DAPI and Beta_actin excluded from clustering as in dedicated pLN",
            "cluster_seed": "stable BLAKE2b seed from donor and annotation ID",
            "cross_assay": "within-donor cluster-profile Spearman comparison to dedicated pLN PhenoCycler profiles on shared channels; distinct slides and no cell matching",
        },
        "analysis_guardrails": [
            "Cell calls remain cluster-inherited and unresolved when the established pLN chain cannot name them.",
            "Same-block and different-block cells remain separate annotation records and retain the physical block relationship.",
            "The 6534 annotation is flagged as whole-slide scope; its pLN-region labeling is not interpreted as an isolated node cluster.",
            "Profile correlation is descriptive and does not transfer labels between sections.",
            "There is no islet-to-node distance analysis.",
        ],
        "n_cells_in_release": int(sum(len(x) for x in all_assignments)),
        "n_annotations": len(node_rows),
        "nodes": node_rows,
        "outputs": {
            "cell_level": str((out / "embedded_pln_cell_typing.parquet").resolve()),
            "node_coverage": str((out / "node_coverage.csv").resolve()),
            "cluster_profiles": str((out / "cluster_profiles.parquet").resolve()) if profile_rows else None,
            "dedicated_profile_comparison": str((out / "dedicated_pLN_profile_comparison.csv").resolve()),
        },
    }
    (out / "embedded_pln_typing.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(f"wrote embedded pLN typing: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
