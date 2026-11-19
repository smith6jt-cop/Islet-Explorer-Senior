#!/usr/bin/env python3
"""Record dated, evidence-backed evaluation observations; never fit or approve science."""
import json
import math
import shutil
import sqlite3
from pathlib import Path

import pandas as pd
from multiplex5000.parameter_observations import (
    CONTEXT_KEYS,
    SCHEMA,
    ObservationHistory,
    canonical,
    file_hash,
    record_or_queue,
)

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"
DOC = BASE / "docs/evaluations/dynamic_portability_20261102"
PROJECT = RUN / "parameter-history"


def read(path):
    return json.loads(Path(path).read_text())


def ref(path):
    path = Path(path)
    return {"id": str(path.relative_to(BASE)) if path.is_relative_to(BASE) else str(path),
            "path": str(path), "sha256": file_hash(path)}


def setting(name, value, source, kind="fixed", units=None):
    return {"name": name, "value": value, "source": str(source), "kind": kind, "units": units}


def available(value):
    """Represent producer NaN/Infinity as unavailable while retaining its original artifact."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: available(v) for k, v in value.items()}
    if isinstance(value, list):
        return [available(v) for v in value]
    return value


def context(modality, assay, tissue, panel, relationship, ids, producer, units, receipt=None):
    c = {**dict.fromkeys(CONTEXT_KEYS), "modality": modality, "assay": assay,
         "tissue": tissue, "panel": panel, "relationship_type": relationship,
         "acquisition_ids": ids, "section_ids": [], "segmentation_ids": [],
         "measurement_units": units, "producer_version": producer,
         "registry_hashes": {}, "panel_hashes": {}, "input_receipts": []}
    if receipt:
        c["input_receipts"] = [ref(receipt)]
        c["panel_hashes"] = {x["name"]: x["sha256"] for x in read(receipt)["inputs"]
                             if ("panel" in x["name"] or "manifest" in x["name"] or "rule_sheet" == x["name"])
                             and x.get("sha256")}
    return c


def observation(name, status, operation, ctx, subject, settings, outcome, evidence):
    return {"schema": SCHEMA, "observation_id": "portability-20261102:" + name,
            "revision": 1, "supersedes": None, "status": status, "operation": operation,
            "context": ctx, "subject": subject, "settings": settings,
            "quality": {"before": None, "after": None, "biological_acceptance": None},
            "outcome": outcome, "evidence": [ref(p) for p in evidence]}


def table_records(path):
    return json.loads(pd.read_csv(path, keep_default_na=True).to_json(orient="records"))


def build():
    pins = read(DOC / "pins.json")
    ledger = read(DOC / "acquisition_dispositions.json")
    by_id = {x["acquisition_id"]: x for x in ledger}
    records = []
    embedding_path = RUN / "remote/rna/pln_embedding_normalization.json"
    embedding_targets = {r["acquisition_id"]: r for r in read(embedding_path)["acquisitions"]}
    for row in read(RUN / "remote/rna/per_acquisition.json"):
        aid = row["acquisition_id"]
        inv = by_id[aid]
        rp = RUN / "remote/rna/receipt.json"
        c = context("rna", "Xenium", row["tissue"], inv["panel"], "independent_acquisition",
                    [aid], pins["remote_RNA_producer"], "raw transcript counts; log1p(raw*10000/cell_depth)", rp)
        c["segmentation_ids"] = [inv["segmentation_id"]] if inv["segmentation_id"] else []
        evidence = [rp, RUN / "remote/rna/per_acquisition.json", RUN / "remote/rna/state_reconciliation.json"]
        settings = [setting("detection_normalization_target_sum", 10000, "tissues/*_xenium_worker.TARGET_SUM", units="transcripts/cell"),
                    setting("detection_predicate", "lognorm >= log1p(1)", "tissues/*_xenium_worker.DETECTION_MIN_LOGNORM"),
                    setting("region_min_median_detected_genes", 20, str(rp) + "#settings.rna_min_median_detected_genes", units="genes/cell"),
                    setting("nucleus_required", True, "stages/io60_rna_typing_pilot._worker_result"),
                    setting("frozen_partition_overrides", "not_applied_to_new_clusters", rp, "selected_method"),
                    setting("state_counts", {k: row[k] for k in ("cells", "called_cells", "broad_cells", "unresolved_cells", "uncalled_cells", "unavailable_cells")},
                            RUN / "remote/rna/state_reconciliation.json", "fitted", "cells")]
        if row["tissue"] == "pancreatic_lymph_node":
            settings += [setting("clustering_fixed_settings", {"min_counts": 10, "min_features": 5,
                                       "subsample": 60000, "resolution": 1.0, "seed_base": 20260917,
                                       "n_hvg": 2000, "n_pcs": 30, "knn_k": 15, "scale_max_value": 10},
                                 "tissues/lymph_node_xenium_worker.partition_section and pln_clustering"),
                         setting("embedding_normalization_method", "normalize_total(target_sum=None): section median after nucleus and cluster QC", "tissues/lymph_node_xenium_worker.partition_section", "selected_method"),
                         setting("reconstructed_embedding_target_sum", embedding_targets[aid]["normalization_target"], embedding_path, "fitted", "gene counts/cell")]
            evidence.append(embedding_path)
        else:
            settings += [setting("clustering_fixed_settings", {"min_counts": 10, "min_features": 5,
                                       "n_hvg": 2000, "n_pcs": 30, "pca_fit_max": 250000,
                                       "subsample": 60000, "resolution": 2.0, "knn_k": 15,
                                       "refine_min_cells": 2000, "refine_resolution": 1.0},
                                 "tissues/pancreas_xenium_worker constants"),
                         setting("embedding_input_normalization", "fixed target 10000 log1p matrix", "tissues/pancreas_xenium_worker.process_section", "selected_method")]
        report = RUN / f"remote/rna/reports/{row['sample_key']}.json"
        if report.exists():
            data = read(report)
            evidence.append(report)
            for key in ("cluster_meta", "refine_meta", "partition", "n_genes_on_panel", "n_wanted_genes_present",
                        "wanted_genes_absent", "named_in_brief_but_absent", "endocrine_hormones_absent_from_panel",
                        "median_transcripts_per_cell", "median_genes_per_cell", "cells_with_depth_gt_target",
                        "y_floor_detect_frac_any", "y_probe_detect_frac_any", "cd3_gate_control", "t_nk_split"):
                if key in data:
                    settings.append(setting(key, data[key], report, "fitted"))
            if "naming_policy" in data:
                settings.append(setting("naming_policy", data["naming_policy"], report, "selected_method"))
        records.append(observation("rna:" + aid, "blocked" if row["status"] == "unavailable" else "evaluated",
                                   "rna_typing", c, {"kind": "genes", "identifiers": []}, settings,
                                   row["reason"] + "; computational reconciliation only; identity-core gene coverage and fitted diagnostics are in the hashed report; unavailable call resolution remains unavailable", evidence))
    rp = RUN / "remote/spleen/receipt.json"
    registry = RUN / "remote/spleen-code/core/multiplex5000/multiplex5000/positivity/registry"
    for row in read(RUN / "remote/spleen/per_acquisition.json"):
        inv = next(x for x in ledger if x["cohort"] == "spleen" and x["modality"] == "phenocycler" and x["donor_id"] == row["section_id"])
        c = context("protein", "PhenoCycler", "spleen", row["panel"], "independent_acquisition", [inv["acquisition_id"]],
                    pins["remote_spleen_producer"], "raw REDSEA intensity; levels in log1p and raw units", rp)
        c["section_ids"] = [inv["physical_section_id"]]
        c["segmentation_ids"] = [inv["segmentation_id"]]
        c["registry_hashes"] = {p.name: file_hash(p) for p in sorted(registry.glob("phenocycler-pc67__spleen*.csv"))}
        settings = [setting("method", read(rp)["settings"]["method"], rp, "selected_method")]
        evidence = [rp, RUN / "remote/spleen/per_acquisition.json"]
        if row["identity_check"] == "pass":
            for table in ("levels", "levels_differentiation", "canary_levels", "selection", "differentiation_selection"):
                path = RUN / f"remote/spleen/tables/{table}.csv"
                rows = [r for r in table_records(path) if "image" not in r or r["image"] == row["section_id"]]
                settings.append(setting(table, rows, path, "selected_method" if "selection" in table else "fitted",
                                        "per-row measurement and raw/log1p levels"))
                evidence.append(path)
        else:
            path = RUN / "remote/spleen_input_geometry/per_acquisition.json"
            source = next(r for r in read(path) if r.get("section_id") == row["section_id"])
            settings.append(setting("source_cell_and_geometry_audit", source, path, "fitted"))
            evidence.append(path)
        records.append(observation("spleen:" + row["section_id"], "evaluated" if row["identity_check"] == "pass" else "blocked",
                                   "protein_typing", c, {"kind": "markers", "identifiers": []}, settings,
                                   row.get("reason", "All eligible native IDs reconciled; unavailable raw image hashes and biological transfer review remain; not original CODEX-port parity"), evidence))
    proposals = RUN / "typing_port/core/multiplex5000/multiplex5000/positivity/registry/phenocycler-pc67__spleen__proposals.csv"
    sp = RUN / "remote/spleen_rna_inputs.json"
    for row in read(sp)["acquisitions"]:
        c = context("rna", "Xenium", "spleen", "gene_set_sha256:" + row["gene_panel_sha256"],
                    "independent_acquisition", [row["acquisition_id"]], "audit_spleen_rna_inputs.py/1", "raw transcript counts")
        c["section_ids"] = [row["physical_section_id"]]
        c["panel_hashes"] = {"observed_gene_names": row["gene_panel_sha256"], **read(sp)["core_hashes"]}
        records.append(observation("spleen-rna-inputs:" + row["acquisition_id"], "evaluated", "rna_input_validation", c,
                                   {"kind": "feature_set", "identifiers": ["gene_set_sha256:" + row["gene_panel_sha256"]]},
                                   [setting(k, row[k], sp, "fitted") for k in ("cells", "matrix_shape", "gene_count", "nucleus_eligible", "no_nucleus", "core_coverage")],
                                   "All matrix/native IDs reconcile; actual gene coverage recorded. Full independent spleen RNA typing remains blocked; no borrowed protein labels", [sp]))
    for row in table_records(proposals):
        c = context("protein", "PhenoCycler", "spleen", "phenocycler-pc67", "independent_acquisition", [], pins["typing_port"], None)
        c["registry_hashes"] = {proposals.name: file_hash(proposals)}
        records.append(observation("proposal:" + row["id"], "proposed", "protein_typing_specification", c,
                                   {"kind": "markers", "identifiers": []}, [setting(row["id"], row, proposals, "review_criterion")],
                                   "Unapproved at the pinned typing port; later producer changes are a separate method/version", [proposals]))
    for row in read(RUN / "codex_parity/summary.json"):
        c = context("protein", "CODEX", "pancreas", "codex-uniform32", "independent_acquisition", [row["image"]], pins["typing_port"], "whole-cell/membrane raw intensity")
        path = RUN / f"codex_parity/{row['image']}.log"
        records.append(observation("codex:" + row["image"], "blocked", "protein_typing_parity", c,
                                   {"kind": "acquisition", "identifiers": [row["image"]]},
                                   [setting("comparison", "exact tables, check_exact=True", "pinned test_positivity_parity.py", "review_criterion")],
                                   "Provenance fixture rejected pre-cutoff table mtimes; no numerical parity comparison occurred", [path]))
    recovered = DOC / "codex_numerical_parity.json"
    if recovered.exists():
        numerical = read(recovered)
        assert numerical["status"] == "complete" and numerical["images_completed"] == 20
        direct_path = DOC / "codex_source_to_port_diagnostics.json"
        direct = read(direct_path)
        assert direct["status"] == "complete" and direct["images_completed"] == 20
        for result in numerical["results"]:
            record = next(r for r in records if r["observation_id"] == "portability-20261102:codex:" + result["image"])
            record["status"] = "failed"
            record["settings"] += [setting("journal_validated_numerical_results", result, recovered, "fitted"),
                                   setting("exact_original_source_vs_port", next(r["direct_source_to_port"] for r in direct["results"] if r["image"] == result["image"]), direct_path, "fitted")]
            record["outcome"] = ("Recovered journal binds the original producer and inputs with stated historical hash limits. "
                                 "Unchanged identity assertion passes exactly; levels/canary exact CSV assertions fail. "
                                 "Direct unchanged original-source versus port computation passes exactly in the same environment. "
                                 "Producer-thread and CSV-parser diagnostics are separate evidence; no criterion was weakened, "
                                 "and the original mtime fixture remains unchanged.")
            record["evidence"] += [ref(p) for p in (DOC / "codex_journal_proof.json", recovered, direct_path,
                                                    DOC / "codex_producer_single_thread_tables.json")]
            records.append(observation("codex-source-port:" + result["image"], "evaluated", "original_source_port_comparison",
                                       dict(record["context"]), dict(record["subject"]),
                                       [setting("exact_original_source_vs_port", next(r["direct_source_to_port"] for r in direct["results"] if r["image"] == result["image"]), direct_path, "fitted"),
                                        setting("eligible_identity_cells", result["identity_cells"], recovered, "fitted", "cells"),
                                        setting("criterion", "exact equality; no tolerance", direct_path, "review_criterion")],
                                       "Unchanged original and port levels/canary/identity computations pass exactly in one environment; "
                                       "eligible IDs reconcile. This is computational portability evidence, not biological or v20 approval. "
                                       "Separate artifact-check observation retains the original guard/CSV failures.",
                                       [DOC / "codex_completion_summary.json", DOC / "codex_journal_proof.json", recovered, direct_path]))
    vp = RUN / "codex_v20/summary.json"
    v20 = read(vp)
    c = context("protein", "CODEX", "pancreas", "codex-uniform32", "independent_acquisition", [],
                v20["producer"], "per-image levels and identity tables")
    records.append(observation("codex-v20-provenance", "failed", "provenance_validation", c,
                               {"kind": "acquisition", "identifiers": []},
                               [setting("bound_input_hashes", v20["checks"], vp, "fitted"),
                                setting("required_binding", "exact recorded SHA256", vp, "review_criterion")],
                               "22 of 24 v20 manifest bindings mismatch current files, including all 20 identities; historical pass is not current compatibility evidence", [vp]))
    for row in table_records(RUN / "maldi/sections.csv"):
        aid = row["acquisition_id"]
        c = context("maldi", "MALDI", None, "38 N-glycans + 10 glycogen", "serial_section", [aid], pins["maldi_branch"], row["source_units"])
        c["section_ids"] = [row["tissue_id"]]
        for op in ("maldi_masking", "registration", "maldi_object_summary"):
            settings = [setting("spot_pitch", 50, "profiles/platforms/maldi.yaml", units="um")]
            if op == "maldi_masking":
                settings += [setting("on_tissue", "sum of finite feature intensities > 0", "multiplex5000.maldi.on_tissue", "selected_method"),
                             setting("zero_semantics", "not detected; preserved measurement", "multiplex5000.maldi"),
                             setting("mask_counts", {k: row[k] for k in ("spots", "on_tissue", "all_zero", "all_zero_on_tissue", "missing_values")},
                                     RUN / "maldi/sections.csv", "fitted", "spots/values")]
                outcome = "Mask and paired spots audited; all-zero raster spots excluded from masked detection denominator; calibrated intensity units unavailable"
            elif op == "registration":
                settings += [setting("accepted_transform", None, "blank inspection inventory", "fitted", "um"),
                             setting("physical_scale", "fixed 50 um sampling; no same-cell interpretation", "profiles/platforms/maldi.yaml")]
                outcome = "Blocked: inspected tissue and exact neighbour absent; Xenium-neighbour registration not built"
            else:
                settings += [setting("eligibility_margin", "registration_RMSE + spot_pitch/2", "stages/maldi_serial_integration.py", "selected_method", "um"),
                             setting("distance", "continuous signed distance", "stages/maldi_serial_integration.py", "selected_method", "um"),
                             setting("object_values", None, "no accepted transform or object geometry", "fitted")]
                outcome = "Blocked: no accepted registration; pancreas islet/structure export dependency; dedicated pLN object contract absent; no cell types or distance bands"
            records.append(observation(op + ":" + aid, "evaluated" if op == "maldi_masking" else "blocked", op, c,
                                       {"kind": "feature_set" if op == "maldi_masking" else "geometry", "identifiers": []},
                                       settings, outcome, [RUN / "maldi/summary.json", RUN / "maldi/sections.csv", RUN / "maldi/detection_denominators.csv"]))
    for pair in read(DOC / "pair_dispositions.json")["explicit_pairs"]:
        ids = ([pair["protein_entry"], pair["xenium_acquisition_id"]] if "protein_entry" in pair
               else pair["protein_acquisition_ids"] + pair["rna_acquisition_ids"])
        c = context("protein/rna", "PhenoCycler/Xenium", pair["tissue"], "PC67 / RNA panel requires verification" if pair["tissue"] == "spleen" else None,
                    pair["relationship"], ids, pins["maldi_branch"], "um")
        for op in ("registration", "same_slide_matching"):
            settings = ([setting("accepted_transform", None, "no receipt-bound accepted Stage 07 transform", "fitted")]
                        if op == "registration" else [setting("max_distance", 5, "stages/same_slide_integration.py MAX_DIST_UM", units="um"),
                                                      setting("matching", "mutual nearest neighbours", "integration/paired_typing.py", "selected_method"),
                                                      setting("registration_uncertainty", None, "not quantified on these pairs", "fitted", "um")])
            records.append(observation(op + ":" + pair["pair_id"], "blocked", op, c,
                                       {"kind": "acquisition_pair", "identifiers": ids}, settings,
                                       "Sources/declarations retained; missing accepted geometry and uncertainty model; no cell links, hierarchy concordance or identity transfer asserted",
                                       [DOC / "pair_dispositions.json"]))
    for row in read(RUN / "pln/per_acquisition.json"):
        rp = RUN / "pln/receipt.json"
        c = context("protein", "PhenoCycler", "pancreatic_lymph_node", "pilot 55 measured markers", "independent_acquisition",
                    [row["acquisition_id"]], read(rp)["code"]["commit"], "section-fitted normalized background-to-bright units", rp)
        c["segmentation_ids"] = [row["segmentation_id"]]
        tp = RUN / "pln_recheck/protein_thresholds.csv"
        thresholds = [r for r in table_records(tp) if r.get("protein_entry") == row["protein_entry"]]
        records.append(observation("pln:" + row["protein_entry"], "evaluated", "protein_typing", c,
                                   {"kind": "markers", "identifiers": []},
                                   [setting("candidate_cutoff", 0.10, "stages/pln_typing_pilot.py; protein_evidence.py", units="background-to-bright unit"),
                                    setting("eligible_cells", int(row["eligible"]), RUN / "pln/per_acquisition.json", "fitted", "cells"),
                                    setting("fitted_thresholds", thresholds, tp, "fitted", "per-row normalized intensity units")],
                                   "Pilot candidate reconciliation only; original manifest recovered exactly; no final cell classes or full 19-section support",
                                   [rp, RUN / "pln/per_acquisition.json", RUN / "pln/provenance_recovery.json", tp]))
    for inv in ledger:
        if (inv["cohort"] == "io60" and inv["modality"] == "phenocycler"
                and inv["tissue"] == "pancreatic_lymph_node" and inv["disposition"] == "blocked"):
            c = context("protein", "PhenoCycler", inv["tissue"], inv["panel"], "independent_acquisition",
                        [inv["acquisition_id"]], None, None)
            c["segmentation_ids"] = [inv["segmentation_id"]] if inv["segmentation_id"] else []
            records.append(observation("pln-blocked:" + inv["acquisition_id"], "blocked", "protein_typing", c,
                                       {"kind": "acquisition", "identifiers": [inv["acquisition_id"]]}, [],
                                       "No resolved or fitted settings: this declared acquisition lacks the bound upstream artifacts available for the two-section pilot. The configurable candidate route is not intrinsically limited to two sections; no final cell classes or production acceptance.",
                                       [DOC / "inventories/io60_acquisitions.csv", RUN / "pln/per_acquisition.json",
                                        DOC / "dependencies.csv"]))
    failed_root = next((RUN / "releases/protein-evaluation").glob(".06_typing.incomplete-*"))
    primary_audit = failed_root / "pancreas/subpanel_audit.json"
    recovery = RUN / "acquisition_runs_recovery.json"
    retry = RUN / "protein_primary"
    source_audit = retry / "pancreas/subpanel_audit.json"
    full = retry / "per_acquisition.json"
    full_rows = {r["donor_id"]: r for r in read(full)}
    nucleus = RUN / "protein_nucleus.json"
    coordinates = RUN / "protein_primary_coordinates.json"
    nucleus_rows = {r["donor_id"]: r for r in read(nucleus)["rows"]}
    coordinate_rows = {r["donor_id"]: r for r in read(coordinates)["rows"]}
    rep = retry / "pancreas/refinement/representation_params.parquet"
    representation = json.loads(pd.read_parquet(rep).reset_index().to_json(orient="records"))
    for row in read(source_audit):
        inv = next(x for x in ledger if x["cohort"] == "io60" and x["modality"] == "phenocycler"
                   and x["tissue"] == "pancreas" and x["donor_id"] == row["donor_id"])
        c = context("protein", "PhenoCycler", "pancreas_with_embedded_pLN", inv["panel"], "independent_acquisition",
                    [inv["acquisition_id"]], pins["maldi_branch"], "aligned protein units with cohort scale; section-fitted geometry rules")
        c["segmentation_ids"] = [inv["segmentation_id"]] if inv["segmentation_id"] else []
        inner = read(failed_root / "pancreas/run.json")
        c["registry_hashes"] = {"typing_registry_log": inner["registry_criteria_sha256"]}
        evidence = [primary_audit, failed_root / "pancreas/run.json", failed_root / "failure.json", recovery]
        evidence += [source_audit, full, retry / "pancreas/run.json", retry / "failure.json",
                     retry / "pancreas/refine.log", rep, retry / "pancreas/refinement/representation_provenance.json",
                     RUN / "protein_primary_expected_universe.json", nucleus, coordinates]
        c["input_receipts"] = [ref(retry / "upstream" / s / "receipt.json")
                               for s in ("00_inventory", "03_protein_measurement", "05_normalization")]
        scale = retry / "pancreas_scale.parquet"
        scale_values = json.loads(pd.read_parquet(scale).to_json(orient="records"))
        evidence.append(scale)
        record = observation("protein-mixed:" + row["donor_id"], "failed", "protein_typing", c,
                             {"kind": "markers", "identifiers": full_rows[row["donor_id"]]["markers_carried"]},
                             [setting("primary_fitted_diagnostics", available(row), source_audit, "fitted"),
                              setting("cohort_scale_method", {"scope": "QC eligible and in_alignment_block", "sample_max_per_donor": 200000,
                                          "seed": 0, "histogram_bins": 512, "min_positive_for_background": 1000,
                                          "top_percentile": 99.9, "missing_floor_fallback": 0, "invalid_span_fallback": 1}, "cohort_scale.py"),
                              setting("fitted_cohort_scale", scale_values, scale, "fitted", "aligned intensity"),
                              setting("refinement_fixed_settings", {"K": 30, "LOW": 0.40, "BOUND_PCT": 5, "MIN_CLEAR": 50,
                                          "seed": 0, "representation_sample_per_donor": 200000, "cofactor_median_divisor": 5,
                                          "invalid_cofactor_fallback": 1, "invalid_spread_fallback": 1}, "scripts/typing/refine.py"),
                              setting("refinement_representation_scope", "all 25 configured donors; geometry-QC eligible inner join; no additional alignment-block filter", "scripts/typing/refine.py:representation_params", "selected_method"),
                              setting("fitted_representation", representation, rep, "fitted", "arcsinh intensity and standardized SD"),
                              setting("refinement_admission_diagnostics", available(full_rows[row["donor_id"]]), full, "fitted"),
                              setting("nucleus_admission_audit", nucleus_rows[row["donor_id"]], nucleus, "fitted"),
                              setting("coordinate_convention_audit", coordinate_rows[row["donor_id"]], coordinates, "fitted", "um; calibration um/pixel"),
                              setting("pln_candidate_cutoff", 0.10, "stages/io60_mixed_protein_pilot.py", units="background-to-bright units")],
                             "Corrected metadata retry built all 25 primary tables, then failed at 0112 refinement: an unacquired Ker8_18 column is treated as carried and leaves zero measurable cells. Fixed settings and fitted representation were recorded, but no refinement bounds or mixed output exist. Original metadata failure remains in history.", evidence)
        records.append(record)
    records.extend(followup_records(pins, ledger))
    return records


def followup_records(pins, ledger):
    """Read-only source/parameter follow-up; never promote observations to approval."""
    records = []
    visual_path = DOC / "maldi_source_review.json"
    if visual_path.exists():
        for row in read(visual_path)["sections"]:
            c = context("maldi", "MALDI", None, "38 N-glycans + 10 glycogen", "serial_section",
                        [row["acquisition_id"]], pins["maldi_branch"], "exported intensity; physical grid um")
            c["section_ids"] = [row["tissue_id"]]
            records.append(observation("maldi-source-review:" + row["tissue_id"], "proposed", "serial_source_inspection", c,
                                       {"kind": "geometry", "identifiers": []},
                                       [setting("visual_review_observation", row, visual_path, "fitted"),
                                        setting("accepted_transform", None, visual_path, "fitted", "um")],
                                       row["observation"] + " No neighbour, tissue or transform approval; exact serial distance is unrecorded.",
                                       [visual_path, DOC / "maldi_source_evidence.json"]))
    objects_path = DOC / "pln_object_availability.json"
    if objects_path.exists():
        for row in read(objects_path)["rows"]:
            inv = next(r for r in ledger if r["acquisition_id"] == row["acquisition_id"])
            c = context("protein", "PhenoCycler", "pancreatic_lymph_node", inv["panel"], "serial_section",
                        [row["acquisition_id"]], "existing QuPath annotation export; source hashes retained", "level-zero image pixels; calibration um/pixel")
            records.append(observation("pln-object-source:" + row["qupath_entry"], "evaluated", "object_source_audit", c,
                                       {"kind": "geometry", "identifiers": row["native_object_ids"]},
                                       [setting("source_geometry_and_availability", row, objects_path, "fitted")],
                                       "Exact entry/image/acquisition and native object IDs reconciled. Rectangles are annotation ROIs, "
                                       "not inferred tissue or follicles; non-area geometry unavailable for area summaries. "
                                       "No MALDI pairing or accepted registration.", [objects_path]))
    trace_path = DOC / "parameter_trace_summary.json"
    if trace_path.exists():
        c = context("multimodal", None, None, None, "evaluation", [], pins["maldi_branch"], None)
        records.append(observation("parameter-execution-audit", "evaluated", "parameter_source_audit", c,
                                   {"kind": "feature_set", "identifiers": []},
                                   [setting("trace_coverage_and_limits", read(trace_path), trace_path, "fitted")],
                                   "All catalog declarations and registry rows dispositioned; literal candidates separately enumerated. "
                                   "Static reads do not prove runtime branches or values; source/report/receipt bindings are explicit. "
                                   "No automatic application or scientific instrumentation change.",
                                   [trace_path, DOC / "parameter_traces.json.gz", DOC / "parameter_literal_candidates.json.gz",
                                    DOC / "parameter_source_manifest.json"]))
    return records


def main():
    records = build()
    global_db = Path.home() / ".kintsugi/parameter_learning.db"
    before = file_hash(global_db) if global_db.exists() else None
    store = ObservationHistory(PROJECT, BASE / "external/KINTSUGI")
    # Rebuilding the evidence may correct an attribution or add a completed route.
    # Keep the previous payload; append only if the actual content changed.
    for record in records:
        previous = store.history(record["observation_id"])
        if previous:
            last = previous[-1]["observation"]
            record.update(revision=last["revision"], supersedes=last["supersedes"])
            if canonical(record) != canonical(last):
                record.update(revision=last["revision"] + 1, supersedes=last["revision"])
    results = [record_or_queue(PROJECT, BASE / "external/KINTSUGI", r) for r in records]
    if any(r["status"] != "recorded" for r in results):
        raise RuntimeError("retryable recording failure; inspect durable outbox")
    retried = [record_or_queue(PROJECT, BASE / "external/KINTSUGI", r) for r in records]
    assert all(r["duplicate"] for r in retried)
    exported = store.export()
    after = file_hash(global_db) if global_db.exists() else None
    assert before == after
    shutil.copyfile(exported["path"], DOC / "parameter_history.json")
    with sqlite3.connect(store.db) as connection:
        approved = connection.execute("SELECT COUNT(*) FROM parameter_records WHERE user_approved = 1").fetchone()[0]
    assert approved == 0
    latest = {r["observation"]["observation_id"]: r["observation"] for r in store.history()}
    summary = {"observations_submitted": len(records), "records": exported["records"],
               "latest_status_counts": {status: sum(r["status"] == status for r in latest.values())
                                        for status in ("proposed", "evaluated", "failed", "blocked")},
               "approved_records": approved,
               "all_retries_duplicate": True, "global_before": before, "global_after": after,
               "automatic_application": False, "recommendations_enabled": False,
               "export": exported, "database_sha256": file_hash(store.db)}
    (PROJECT / "recording_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
