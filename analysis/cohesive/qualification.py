"""Reconcile current analytical populations; never refit or promote identity."""
from __future__ import annotations

from collections import Counter
import hashlib
import importlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .common import Context, binding, read_json, require, unique_bindings, verify_binding, write_json
from .rna_measurements import exact_index, project_expression
from .sources import audit_binding, original_reader, recorded_bindings


MARKERS = ("INS", "GCG", "SST", "CD99", "HLA_A", "PD_L1", "PD_1",
           "Granzyme_B", "CD107a", "CD3e", "CD31", "B3TUBB", "Collagen_IV")
MEMBER_FLAGS = ("qualified_broad_endocrine", "qualified_broad_immune",
                "qualified_structure_vascular", "qualified_structure_neural",
                "qualified_structure_lymphatic", "qualified_identity_unresolved")
# The semantic labels exactly as the handoff's build_qualified_broad_composition.py reads them.
ENDOCRINE_LABELS = ("Alpha", "Beta", "Delta", "Endocrine (ambiguous subtype)",
                    "Endocrine (unresolved subtype)")
STRUCTURE_LABELS = {"Blood vessel-associated allocation": "vascular",
                    "Nerve-associated allocation": "neural",
                    "Lymphatic-associated allocation": "lymphatic"}
# The refined_class spelling of each broad label; "Identity unresolved" has none, so refinement
# can never promote an unresolved cell.
REFINED_NAME = {**{name: name for name in ENDOCRINE_LABELS + ("Immune",)},
                "Blood vessel-associated allocation": "Blood vessel",
                "Nerve-associated allocation": "Nerve",
                "Lymphatic-associated allocation": "Lymphatic"}


def selection(context):
    return read_json(context.task("reconcile") / "selected_inputs.json")


def data_bindings(context):
    chosen = selection(context)
    records = [binding(Path(__file__).with_name("rna_measurements.py")),
               binding(Path(__file__).with_name("sources.py"))]
    for row in chosen["rna_exports"]:
        records.extend(row[key] for key in ("output", "annotation_sidecar_binding", "native_call_nucleus_binding"))
        matrix = row["native_expression_matrix_binding"]
        records.append({"path": matrix["path"], "sha256": matrix["saved_upstream_sha256"]})
        vendor = row["native_geometry_bindings"].get("vendor_cell_table")
        if vendor is not None:
            records.append({"path": vendor["path"], "sha256": vendor["saved_sha256"]})
    for scope in chosen["scopes"]:
        if scope["lane"] == "pancreas_protein":
            records.extend(({"path": scope["path"], "sha256": scope["sha256"]},
                            {"path": scope["native_source_path"], "sha256": scope["native_source_sha256"]}))
    # The selected release manifest names these files but does not hash its outputs.
    # Bind their current full bytes, then independently join the already pinned assignments.
    records.extend(binding(row["full_cell_release_path"]) for row in chosen["protein_records"])
    for row in chosen["protein_records"]:
        records.append(audit_binding(context, "protein_composition/qualified_broad_composition/"
                                     f"donor_id={row['donor_id']}/exact_member_evidence.parquet"))
    records.append(audit_binding(context, "protein_composition/qualified_broad_composition/qualified_broad_objects.parquet"))
    generation = read_json(chosen["contracts"]["corrected_protein_measurement_generation"]["path"])
    block_manifest = Path(generation["current_block_sidecars"]) / "manifest.json"
    source_manifest = Path(chosen["contracts"]["corrected_protein_measurement_generation"]["path"]).parent / "source_manifest.json"
    source = read_json(source_manifest)
    records.extend((binding(source_manifest),
                    {"path": str(block_manifest), "sha256": source["current_block_manifest_sha256"]},
                    binding(Path(generation["geometry_root"]).parent / "repaired_primary_stratum_evidence_QC.csv"),
                    binding(Path(generation["geometry_root"]).parent / "repaired_primary_evidence_QC.json")))
    records.extend({"path": row["snapshot"], "sha256": row["sha256"]}
                   for row in generation["source_scripts"]
                   if Path(row["path"]).name in ("reprocess_protein_own_compartment_candidate.py",
                                                 "run_step2.py", "alignment_block.py"))
    return unique_bindings(records)


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "qualify", "Unknown qualification operation")
    return data_bindings(context)


def block_dependency_audit(context, stage):
    chosen = selection(context)
    generation = read_json(chosen["contracts"]["corrected_protein_measurement_generation"]["path"])
    original = read_json(Path(generation["current_block_sidecars"]) / "manifest.json")
    corrected = read_json(chosen["contracts"]["corrected_block_metadata"]["path"])
    verify_binding({"path": original["blocks_file"], "sha256": original["blocks_sha256"]})
    require(original["blocks_sha256"] == corrected["authority"]["sha256"],
            "Calibration used a different block authority; targeted recalibration is required")
    relations = {(row["image"], uid): relation for row in original["donors"]
                 for uid, relation in row["annotations"].items()}
    rows = []
    for row in corrected["records"]:
        key = (row["acquisition_identifier"], row["annotation_UUID"])
        require(key in relations, f"Corrected block not found in fitted calibration: {key}")
        agree = relations[key] == row["effective_block_relationship"]
        rows.append({"scope_id": row["scope_id"], "fitted_block": relations[key],
                     "effective_block": row["effective_block_relationship"], "agrees": agree})
    require(all(row["agrees"] for row in rows),
            "Some fitted block selections differ; do not reuse affected calibration")
    pd.DataFrame(rows).to_csv(stage / "calibration_block_relationship_replay.csv", index=False)
    write_json(stage / "calibration_block_dependency.json", {
        "status": "all_exact_annotation_relationships_already_correct_in_calibration",
        "fitted_block_manifest": binding(Path(generation["current_block_sidecars"]) / "manifest.json"),
        "corrected_relationships": chosen["contracts"]["corrected_block_metadata"],
        "matched_annotations": len(rows), "recalibration_required_for_this_metadata_change": False,
        "reason": "The fitted calibration already used the same exact block authority and all corrected annotation relationships.",
        "does_not_establish_biological_accuracy": True,
    })


def protein_donor(context, stage, record, scope, view):
    donor = record["donor_id"]
    print(f"qualify protein {donor}: current semantic/nucleus/measurement join", flush=True)
    full = pq.ParquetFile(record["full_cell_release_path"])
    columns = ["object_id", "X_centroid", "Y_centroid", "nucleus_area_um2",
               "fundamentals_cell_area_um2", "qc_analysis_eligible", "endocrine_object_id",
               "endocrine_object_class", "cell_class", "immune_class", "hormone_measurement_state",
               "endocrine_subtype_status", "islet_route", "islet_stratum",
               "islet_cd99_level", "islet_elongation_level", "acq_run", "stain_group",
               "refined_class", "refinement_state"]
    require(set(columns).issubset(full.schema_arrow.names), f"Missing current source fields: {donor}")
    measured = [name for name in (["redsea_aligned_" + marker for marker in MARKERS] +
                                 ["typing_scaled_" + marker for marker in ("INS", "GCG", "SST", "CD99")])
                if name in full.schema_arrow.names]
    native = full.read(columns=columns + measured).to_pandas()
    source_index = pd.Index(native.object_id)
    require(source_index.is_unique and not source_index.hasnans, f"Duplicate native UUIDs: {donor}")
    assignments = pq.ParquetFile(scope["native_source_path"]).read(
        columns=["object_id", "X_centroid", "Y_centroid", "qc_analysis_eligible",
                 "analysis_eligible", "endocrine_object_id", "endocrine_object_class"]).to_pandas()
    positions = exact_index(source_index, assignments.object_id, donor + " pinned assignments")
    require(len(assignments) == len(native), f"Native/assignment universe differs: {donor}")
    for field in ("endocrine_object_id", "endocrine_object_class", "qc_analysis_eligible"):
        pd.testing.assert_series_equal(native.iloc[positions][field].reset_index(drop=True),
                                       assignments[field].reset_index(drop=True), check_names=False,
                                       check_dtype=False)
    require(np.allclose(native.iloc[positions][["X_centroid", "Y_centroid"]].to_numpy(float),
                        assignments[["X_centroid", "Y_centroid"]].to_numpy(float), equal_nan=True),
            f"Current native/assignment coordinates differ: {donor}")
    admission = pd.Series(assignments.analysis_eligible.to_numpy(), index=assignments.object_id)
    folder = stage / "protein" / f"donor_id={donor}"
    folder.mkdir(parents=True)
    output = folder / "cells.parquet"
    schema = pa.schema(
        [full.schema_arrow.field(name) for name in columns + measured] +
        [pa.field("annotation_name", pa.string()), pa.field("evidence_level", pa.string()),
         pa.field("broad_endocrine", pa.bool_()), pa.field("broad_immune", pa.bool_()),
         pa.field("structure_associated", pa.bool_()), pa.field("refinement_keeps_broad_label", pa.bool_()),
         pa.field("donor_id", pa.string()),
         pa.field("disease_status", pa.string()), pa.field("coordinates_available", pa.bool_())]
    )
    writer = None
    totals = Counter()
    marker_rows, evidence_parts = [], []
    seen = np.zeros(len(native), dtype=bool)
    diameter_sum = 0.0
    diameter_n = 0
    subtype_counts = Counter()
    try:
        for batch in view.iter_annotations(scope["scope_id"], columns=["object_id", "candidate_cell_class"],
                                           verify_file=False, batch_size=131072):
            sem = batch.to_pandas()
            indices = exact_index(source_index, sem.object_id, donor + " current semantics")
            require(not seen[indices].any(), f"Repeated semantic UUID: {donor}")
            seen[indices] = True
            source = native.iloc[indices].reset_index(drop=True)
            qc = sem.qc_analysis_eligible.fillna(False).to_numpy(bool)
            adm = sem.analysis_eligible.fillna(False).to_numpy(bool)
            require(np.array_equal(qc, source.qc_analysis_eligible.fillna(False).to_numpy(bool)),
                    f"Semantic/native QC differs: {donor}")
            require(np.array_equal(adm, admission.reindex(sem.object_id).fillna(False).to_numpy(bool)),
                    f"Semantic/pinned admission differs: {donor}")
            nucleus = source.nucleus_area_um2.to_numpy(float)
            positive = np.isfinite(nucleus) & (nucleus > 0)
            biological = qc & adm & positive
            totals.update(source_rows=len(sem), QC_rows=int(qc.sum()), admission_rows=int(adm.sum()),
                          QC_nonpositive_nucleus=int((qc & ~positive).sum()),
                          admitted_nonpositive_nucleus=int((qc & adm & ~positive).sum()),
                          biological_rows=int(biological.sum()))
            broad = sem.eligible_for_broad_lineage.fillna(False).to_numpy(bool) & biological
            structure = sem.eligible_for_structure.fillna(False).to_numpy(bool) & biological
            label = sem.annotation_name.astype("string")
            endocrine = broad & label.isin(ENDOCRINE_LABELS).fillna(False).to_numpy(bool)
            immune = broad & label.eq("Immune").fillna(False).to_numpy(bool)
            # Broad markers first, then refinement to the true identity (maintainer, 2026-10-17):
            # a cell keeps its lineage only where the refined class is that same lineage.
            kept = label.map(REFINED_NAME).eq(source.refined_class.astype("string")).fillna(False).to_numpy(bool)
            eligible = source.loc[biological].copy().reset_index(drop=True)
            eligible["annotation_name"] = sem.loc[biological, "annotation_name"].to_numpy()
            eligible["evidence_level"] = sem.loc[biological, "evidence_level"].to_numpy()
            eligible["broad_endocrine"] = endocrine[biological]
            eligible["broad_immune"] = immune[biological]
            eligible["structure_associated"] = structure[biological]
            eligible["refinement_keeps_broad_label"] = kept[biological]
            eligible["donor_id"] = donor
            eligible["disease_status"] = record["disease_status"]
            eligible["coordinates_available"] = np.isfinite(
                eligible[["X_centroid", "Y_centroid"]].to_numpy(float)).all(axis=1)
            totals["invalid_coordinates"] += int((~eligible.coordinates_available).sum())
            areas = eligible.fundamentals_cell_area_um2.to_numpy(float)
            available_area = np.isfinite(areas) & (areas > 0)
            diameter_sum += float((2 * np.sqrt(areas[available_area] / np.pi)).sum())
            diameter_n += int(available_area.sum())
            totals["broad_endocrine"] += int(endocrine.sum())
            totals["broad_immune"] += int(immune.sum())
            totals["structure_associated"] += int(structure.sum())
            totals["refined_endocrine"] += int((endocrine & kept).sum())
            totals["refined_immune"] += int((immune & kept).sum())

            endocrine_core = endocrine & source.endocrine_object_class.eq("Islet").fillna(False).to_numpy(bool)
            hormone_columns = ["typing_scaled_" + marker for marker in ("INS", "GCG", "SST")]
            if set(hormone_columns).issubset(source):
                hormones = source.loc[endocrine_core, hormone_columns].to_numpy(float)
                finite = np.isfinite(hormones)
                safe = np.where(finite, hormones, -np.inf)
                maxima = safe.max(axis=1)
                unique = finite.any(axis=1) & ((safe == maxima[:, None]).sum(axis=1) == 1)
                names = np.asarray(["Beta", "Alpha", "Delta"])[safe.argmax(axis=1)]
                labels = source.loc[endocrine_core, "cell_class"].to_numpy()
                usable = source.loc[endocrine_core, "hormone_measurement_state"].eq("usable_for_typing").to_numpy()
                subtype_counts.update(true_islet_broad_endocrine=int(endocrine_core.sum()),
                                      any_hormone_measured=int(finite.any(axis=1).sum()),
                                      all_hormones_measured=int(finite.all(axis=1).sum()),
                                      unique_usable_hormone_argmax=int((unique & usable).sum()),
                                      unique_argmax_disagrees_with_native=int((unique & (names != labels)).sum()),
                                      unmeasured_or_tied_hormones=int((~unique).sum()))
            route = source.islet_route.astype("string").fillna("")
            cd99 = route.eq("cd99").to_numpy(bool)
            elongation = route.eq("elongation").to_numpy(bool)
            cd99_level = source.islet_cd99_level.to_numpy(float)
            shape_level = source.islet_elongation_level.to_numpy(float)
            subtype_counts["cd99_routes_without_positive_level"] += int(
                (biological & cd99 & (~np.isfinite(cd99_level) | (cd99_level <= 0))).sum())
            subtype_counts["shape_routes_without_ratio_above_one"] += int(
                (biological & elongation & (~np.isfinite(shape_level) | (shape_level <= 1))).sum())
            subtype_counts["shape_routes_below_large_stratum"] += int(
                (biological & elongation & ~source.islet_stratum.eq(">10").fillna(False).to_numpy(bool)).sum())

            evidence_parts.append(pd.DataFrame({
                "object_id": sem.object_id, "biological": biological,
                "qualified_broad_endocrine": endocrine, "qualified_broad_immune": immune,
                **{"qualified_structure_" + kind: structure & label.eq(name).fillna(False).to_numpy(bool)
                   for name, kind in STRUCTURE_LABELS.items()},
            }))
            table = pa.Table.from_pandas(eligible, schema=schema, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output, table.schema, compression="zstd")
            writer.write_table(table)
        require(seen.all(), f"Semantic reader omitted native rows: {donor}")
    finally:
        if writer is not None:
            writer.close()
    require(diameter_n > 0, f"No measured cell area for donor scale: {donor}")
    members = pq.ParquetFile(context.inputs.parent / "protein_composition/qualified_broad_composition" /
                             f"donor_id={donor}/exact_member_evidence.parquet").read().to_pandas()
    evidence = pd.concat(evidence_parts, ignore_index=True).set_index("object_id", verify_integrity=True)
    actual = evidence.reindex(members.object_id)
    require(not actual.biological.isna().any(), f"Original member UUID missing: {donor}")
    require(np.array_equal(actual.biological, members.qualified_member), f"Member admission differs: {donor}")
    for flag in MEMBER_FLAGS[:-1]:
        require(np.array_equal(actual[flag], members[flag]), f"Current member semantic mismatch: {donor}/{flag}")
    require(totals["source_rows"] == scope["expected_source_rows"] and
            totals["QC_rows"] == scope["expected_QC_rows"] and
            totals["admission_rows"] == scope["expected_additional_admission_rows"],
            f"Current protein denominator mismatch: {donor}")
    mask = evidence.reindex(native.object_id).biological.to_numpy(bool)
    for marker in MARKERS:
        column = "redsea_aligned_" + marker
        values = native[column].to_numpy(float) if column in native else np.full(len(native), np.nan)
        finite = np.isfinite(values) & mask
        marker_rows.append({"donor_id": donor, "marker": marker, "population": "admitted_positive_nucleus",
                            "denominator": int(mask.sum()), "measured": int(finite.sum()),
                            "mean_aligned_intensity": float(values[finite].mean()) if finite.any() else np.nan,
                            "availability_reason": "measured" if finite.any() else "no_finite_aligned_measurement"})
    write_json(folder / "qualification.json", {
        "donor_id": donor, "scope_id": scope["scope_id"], **dict(totals),
        "mean_cell_diameter_um": diameter_sum / diameter_n, "cell_area_support": diameter_n,
        "exact_original_members_reconciled": len(members),
        "fitted_route_and_hormone_accounting": dict(subtype_counts),
        "new_biological_accuracy_claim": False,
    })
    return {"donor_id": donor, "disease_status": record["disease_status"], **dict(totals),
            **dict(subtype_counts), "mean_cell_diameter_um": diameter_sum / diameter_n}, marker_rows


def qualify_rna(context, stage, chosen):
    reader_path = context.inputs.parent / "identity_review"
    sys.path.insert(0, str(reader_path))
    reader = importlib.import_module("read_current_rna")
    require(Path(reader.__file__).resolve() == reader_path / "read_current_rna.py", "Foreign RNA reader")
    inventory = reader.load_inventory()
    totals, counts, support = [], [], []
    for export in chosen["rna_exports"]:
        sid = export["scope_id"]
        print(f"qualify RNA {sid}: exact native replay and raw measurements", flush=True)
        original = reader.read_current_rna(inventory[sid], biological_only=True,
                                           columns=["candidate_lineage", "candidate_common_lineage"],
                                           verify_hashes=False)
        saved = pq.ParquetFile(export["output"]["path"]).read(columns=original.column_names)
        require(original.equals(saved), f"Current RNA export differs from exact native replay: {sid}")
        frame = saved.to_pandas()
        require(len(frame) == export["output_rows"], f"RNA biological denominator differs: {sid}")
        named = frame.handoff_named_lineage_use_eligible
        require(int(named.sum()) == export["named_lineage_eligible_rows"], f"RNA named denominator differs: {sid}")
        for (label, admitted), part in frame.groupby(["candidate_lineage", "handoff_named_lineage_use_eligible"],
                                                    dropna=False, observed=True):
            counts.append({"scope_id": sid, "lane": export["lane"], "donor_id": export["donor_id"],
                           "source_serial": export["source_serial"], "section_id": export["section_id"],
                           "annotation": label, "named_lineage_eligible": bool(admitted),
                           "cells": len(part), "biological_denominator": len(frame)})
        measurement, audit = project_expression(export, frame)
        destination = stage / "rna" / hashlib.sha256(sid.encode()).hexdigest()[:16]
        destination.mkdir(parents=True)
        measurement.to_parquet(destination / "cells.parquet", index=False, compression="zstd")
        write_json(destination / "qualification.json", audit)
        totals.append({"scope_id": sid, "lane": export["lane"], "donor_id": export["donor_id"],
                       "section_id": export["section_id"], "source_serial": export["source_serial"],
                       "source_rows": inventory[sid]["expected_source_rows"],
                       "biological_rows": len(frame), "named_lineage_rows": int(named.sum()),
                       "unresolved_or_not_named_rows": int((~named).sum()),
                       "relative_data_path": str((destination / "cells.parquet").relative_to(stage)),
                       "coordinates_available": audit["rows_with_finite_coordinates"]})
        support.append(audit)
    pd.DataFrame(totals).to_csv(stage / "rna_scope_qualification.csv", index=False)
    pd.DataFrame(counts).to_csv(stage / "rna_annotation_counts.csv", index=False)
    write_json(stage / "rna_measurement_support.json", support)
    return totals


def qualify(context: Context, stage: Path) -> dict:
    pa.set_cpu_count(4)
    pa.set_io_thread_count(4)
    chosen = selection(context)
    block_dependency_audit(context, stage)
    view = original_reader(context)
    scopes = {row["donor_id"]: row for row in chosen["scopes"] if row["lane"] == "pancreas_protein"}
    protein, markers = [], []
    for row in chosen["protein_records"]:
        result, measured = protein_donor(context, stage, row, scopes[row["donor_id"]], view)
        protein.append(result)
        markers.extend(measured)
    pd.DataFrame(protein).to_csv(stage / "protein_scope_qualification.csv", index=False)
    pd.DataFrame(markers).to_csv(stage / "protein_marker_support.csv", index=False)
    rna = qualify_rna(context, stage, chosen)
    write_json(stage / "qualification_summary.json", {
        "protein_pancreas_scopes": len(protein), "RNA_scopes": len(rna),
        "protein_admitted_positive_nucleus_observations": sum(row["biological_rows"] for row in protein),
        "RNA_nativeQC_positive_nucleus_observations_nonadditive": sum(row["biological_rows"] for row in rna),
        "source_cells_labels_and_latches_unchanged": True,
        "protein_pLN_counts_not_recomputed_without_a_separate_nucleus_binding": True,
        "subtype_endpoint_qualification_not_inferred_from_integrity": True,
        "calibration_already_used_corrected_block_relationships": True,
    })
    return read_json(stage / "qualification_summary.json")
