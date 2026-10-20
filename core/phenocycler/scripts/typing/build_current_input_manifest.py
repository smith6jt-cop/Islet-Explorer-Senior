#!/usr/bin/env python3
"""Inventory current slide, tissue, modality, block evidence and producer coverage.

The manifest is intentionally evidence-led: a missing block or acquisition record is
written as unresolved instead of being inferred from donor identity or slide proximity.
"""
from __future__ import annotations

import argparse
import csv
import configparser
import hashlib
import json
import re
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
from phenocycler.cohort import DONOR_EXCLUSIONS

REPO = Path(__file__).resolve().parents[4]
PROJECT = REPO.parent
INTEGRATION = REPO / "core/phenocycler/data/integration"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(line for line in f if not line.lstrip().startswith("#")))


def key(value: object) -> str:
    return str(value or "").strip().lstrip("0") or "0"


def count_value(value: object) -> int | None:
    if value is None or str(value).strip().lower() in {"", "nan", "none"}:
        return None
    return int(float(value))


def parquet_info(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {"exists": False, "rows": None, "columns": []}
    pf = pq.ParquetFile(path)
    return {"exists": True, "rows": pf.metadata.num_rows, "columns": pf.schema_arrow.names,
            "bytes": path.stat().st_size}


_XENIUM_BUNDLE_INDEX: dict[str, Path] | None = None


def local_xenium(source: str) -> str | None:
    if not source or source.strip().upper() == "NA":
        return None
    path = Path(source.strip())
    if path.is_absolute():
        if path.exists():
            return str(path)
    # The project table stores paths relative to Xenium_Analysis/data, while the local
    # acquisition mirror preserves the full HuPrime5kPanTissue hierarchy under /data/xenium.
    # Try that exact relative path before falling back to a bundle-name search (needed for
    # source paths that point to the remote /blue mount).
    if not path.is_absolute():
        candidate = Path("/data/xenium") / path
        if candidate.is_dir():
            return str(candidate)
    match = re.search(r"(output-XETG\d+__[^/]+)$", source.strip().rstrip("/"))
    if match:
        bundle = match.group(1)
        roots = [Path("/data/xenium"), Path("/data/xenium_bundles"),
                 Path("/data/xenium_pln"), Path("/data/QP_Xenium"),
                 Path("/data/PancLN_Xenium"), Path("/data/Phenocycler_Xenium")]
        for root in roots:
            candidate = root / bundle
            if candidate.is_dir():
                return str(candidate)
        global _XENIUM_BUNDLE_INDEX
        if _XENIUM_BUNDLE_INDEX is None:
            _XENIUM_BUNDLE_INDEX = {}
            for root in roots:
                if root.is_dir():
                    for candidate in root.rglob("output-XETG*"):
                        if candidate.is_dir():
                            _XENIUM_BUNDLE_INDEX.setdefault(candidate.name, candidate)
        candidate = _XENIUM_BUNDLE_INDEX.get(bundle)
        if candidate is not None:
            return str(candidate)
    if not path.is_absolute():
        candidate = Path("/home/smith6jt/IO60panc2nd") / path
        if candidate.exists():
            return str(candidate)
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, type=Path,
                    help="output directory for current_input_manifest.json and .csv")
    ap.add_argument("--embedded-block-ledger", type=Path,
                    default=Path("/data/io60-calibration/alignment_block_current/manifest.json"))
    ap.add_argument("--block-annotations", type=Path,
                    default=Path("/home/smith6jt/IO60panc2nd/results/v3/pln_blocks.csv"))
    ap.add_argument("--rna-inventory", type=Path,
                    default=Path("/data/xenium/io60_cohort_20261003/input_inventory.json"))
    ap.add_argument("--pancreas-rna-release", type=Path,
                    default=Path("/data/xenium_typing/panc_20261003_bundle/release/release.json"))
    ap.add_argument("--pln-rna-release", type=Path,
                    default=Path("/data/xenium_typing/pln_20261003/release.json"))
    ap.add_argument("--postxenium-export", type=Path,
                    default=Path("/data/io60-restain/valis/export_panc.json"))
    ap.add_argument("--postxenium-coverage", type=Path,
                    default=Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/release_20261003/same_slide_postxenium/coverage.csv"))
    ap.add_argument("--postxenium-measurements", type=Path,
                    default=Path("/data/io60-restain/protein_rna_panc.json"))
    ap.add_argument("--postxenium-channel-assignment", type=Path,
                    default=Path("/data/io60-restain/protein_rna/channel_assignment.json"))
    args = ap.parse_args()

    qpath = REPO / "core/phenocycler/data/qupath_manifest.json"
    qman = json.loads(qpath.read_text())
    images = {key(r["donor_id"]): r for r in qman["images"]}
    declared_expected = qman.get("expected_donors", [])
    expected = (set(map(key, declared_expected)) if isinstance(declared_expected, list)
                else set(images)) or set(images)
    if set(images) != expected:
        raise ValueError(f"QuPath manifest donor mismatch: absent {sorted(expected-set(images))}, "
                         f"unexpected {sorted(set(images)-expected)}")

    runs = {key(r["donor_id"]): r for r in read_csv(REPO / "data/donor_acquisition_runs.csv")}
    batch = {key(r["donor_id"]): r for r in read_csv(REPO / "data/donor_acquisition_batch.csv")}
    node_rows = read_csv(INTEGRATION / "pln_node_annotations.csv")
    xenium_rows = read_csv(INTEGRATION / "xenium_paths.csv")
    section_modes = read_csv(INTEGRATION / "section_modes.csv")
    mode_by_key = {(key(r["donor_id"]), r["roi"].strip(), r["xenium_serial"].strip()): r
                   for r in section_modes if r.get("mode")}
    integration_config = REPO / "core/phenocycler/config.ini"
    config = configparser.ConfigParser()
    config.read(integration_config)
    fallback_relationship = config.get("integration", "mode", fallback="unresolved").strip()

    block_manifest = json.loads(args.embedded_block_ledger.read_text()) \
        if args.embedded_block_ledger.exists() else {"donors": [], "donors_without_pln_cells": []}
    embedded = {key(r["donor_id"]): r for r in block_manifest.get("donors", [])}
    block_annotations = read_csv(args.block_annotations) if args.block_annotations.is_file() else []
    block_annotations_by_donor: dict[str, list[dict]] = {}
    for annotation in block_annotations:
        # The historical annotation ledger records 0112/0115 as HDL112/HDL115.
        donor_token = re.sub(r"^HDL", "", str(annotation.get("donor", "")), flags=re.IGNORECASE)
        donor = key(donor_token)
        relationship = str(annotation.get("block_from_notes", "")).strip()
        ledger_relationship = embedded.get(donor, {}).get("annotations", {}).get(
            annotation.get("object_id", ""))
        if ledger_relationship and relationship and ledger_relationship != relationship:
            raise ValueError(
                f"{donor}/{annotation.get('object_id')}: block annotation says "
                f"{relationship!r}, cell-assignment ledger says {ledger_relationship!r}")
        note = str(annotation.get("notes", ""))
        region = note.split(" and ", 1)[0].strip() if note else "unresolved"
        block_annotations_by_donor.setdefault(donor, []).append({
            "annotation_id": annotation.get("object_id", ""),
            "annotation_index": annotation.get("pln_annotation", ""),
            "pancreas_block_region_from_notes": region,
            "block_relationship_to_pancreas": relationship or "unresolved",
            "measured_relation_to_tissue": annotation.get("measured_relation", "unresolved"),
            "gap_to_tissue_um": annotation.get("gap_to_tissue_um", ""),
            "area_mm2": annotation.get("area_mm2", ""),
            "source_notes": note,
        })
    for donor, annotations in block_annotations_by_donor.items():
        ledger_ids = set(embedded.get(donor, {}).get("annotations", {}))
        source_ids = {r["annotation_id"] for r in annotations}
        if ledger_ids != source_ids:
            raise ValueError(
                f"{donor}: block ledger and authoritative annotation table differ; "
                f"missing={sorted(source_ids-ledger_ids)}, stale={sorted(ledger_ids-source_ids)}")
    pprotein_root = REPO / "core/phenocycler/data/integration/cells_pln"
    ptyped_root = REPO / "data/development/typing/pln/pheno"
    pingest_report = pprotein_root.parent / "cells_pln_ingest_report.json"
    pprotein_source = ""
    if pingest_report.is_file():
        ingest = json.loads(pingest_report.read_text())
        pprotein_source = ingest.get("source", {}).get("path", "")
    rna_sections = []
    if args.rna_inventory.is_file():
        rna_sections = json.loads(args.rna_inventory.read_text()).get("sections", [])
    disease_by_donor: dict[str, str] = {}
    excluded_donors = {key(d) for d in DONOR_EXCLUSIONS}
    for section in rna_sections:
        donor = key(re.sub(r"^HDL", "", str(section.get("donor_id", "")),
                            flags=re.IGNORECASE))
        disease = str(section.get("disease_status", "")).strip()
        if not donor or not disease:
            continue
        if donor in disease_by_donor and disease_by_donor[donor] != disease:
            raise ValueError(f"{donor}: inconsistent disease status in Xenium input inventory")
        disease_by_donor[donor] = disease
    panc_rna_release = json.loads(args.pancreas_rna_release.read_text()) \
        if args.pancreas_rna_release.is_file() else {}
    pln_rna_release = json.loads(args.pln_rna_release.read_text()) \
        if args.pln_rna_release.is_file() else {}
    panc_rna_qc = {(key(r.get("donor_id")), str(r.get("section_key", ""))): r
                   for r in panc_rna_release.get("per_donor", [])}
    pln_reports_path = args.pln_rna_release.parent / "section_reports.json"
    pln_rna_qc = {}
    if pln_reports_path.is_file():
        pln_rna_qc = {(key(r.get("donor_id")), str(r.get("section_key", ""))): r
                      for r in json.loads(pln_reports_path.read_text())}

    def qc_summary(report: dict, tissue: str) -> dict:
        if not report:
            return {}
        wanted = ("donor_id", "section_key", "disease_status", "batch", "n_cells",
                  "n_nucleus_excluded", "called_fraction_of_nucleated_cells",
                  "median_transcripts_per_cell", "median_genes_per_cell", "low_depth_flag",
                  "n_clusters", "n_subclusters", "knn_holdout_agreement", "y_probe_detect_frac_any",
                  "y_floor_detect_frac_any", "n_cells_in_h5ad", "n_cells_partitioned",
                  "n_cells_not_partitioned", "n_cells_no_nucleus", "n_genes_on_panel",
                  "n_wanted_genes_present", "wanted_genes_absent", "detection_predicate",
                  "entries_raw_nonzero_in_wanted_set", "entries_detected_in_wanted_set",
                  "entries_nonzero_but_below_predicate")
        return {"tissue": tissue, "qc_source": "producer section report",
                **{k: report[k] for k in wanted if k in report}}

    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                              check=True, capture_output=True, text=True).stdout.strip()
    status = subprocess.run(["git", "status", "--short"], cwd=REPO,
                            check=True, capture_output=True, text=True).stdout.splitlines()
    rows: list[dict] = []

    for donor, image in sorted(images.items()):
        q = image
        runs_row, batch_row = runs.get(donor, {}), batch.get(donor, {})
        current = REPO / f"data/development/typing/revised_cells/donor_id={q['donor_id']}/data.parquet"
        fundamentals = REPO / f"data/development/typing/qupath_fundamentals/donor_id={q['donor_id']}/data.parquet"
        expression = Path(q.get("measurement_csv", {}).get("path", ""))
        geom = q.get("cell_geometry", {})
        rows.append({
            "donor_id": str(q["donor_id"]), "tissue": "pancreas", "slide_id": q["image_id"],
            "physical_block_id": "unresolved", "block_relationship_to_pancreas": "index tissue",
            "modality": "protein_if", "acquisition": runs_row.get("acq_run", "unresolved"),
            "software": runs_row.get("software", batch_row.get("acquisition_software", "unresolved")),
            "panel": q.get("panel_id", "unresolved"),
            "geometry_revision": q.get("segmentation_version", "unresolved"),
            "pixel_size_um_x": q.get("pixel_size_um_x"), "pixel_size_um_y": q.get("pixel_size_um_y"),
            "source_path": q.get("qptiff", {}).get("path", ""),
            "source_sha256": q.get("qptiff", {}).get("sha256", ""),
            "current_measurement": str(expression), "current_cell_geometry": parquet_info(current),
            "current_qupath_fundamentals": parquet_info(fundamentals),
            "current_geometry_feature_count": geom.get("feature_count"),
            "current_geometry_object_id_sha256": geom.get("object_id_sha256"),
            "coverage_status": "available; expansion-ring analysis removed; current geometry/cell IDs reconciled; 8-arm calibration dossier complete and final production policy pending",
            "evidence_source": str(qpath), "code_revision": revision,
        })

        if donor in embedded:
            b = embedded[donor]
            rows.append({
                "donor_id": str(q["donor_id"]), "tissue": "embedded_pLN", "slide_id": q["image_id"],
                "physical_block_id": "unresolved", "block_relationship_to_pancreas": {
                    "same_block_cells": b.get("same-block", 0),
                    "different_block_cells": b.get("different-block", 0),
                    "unassigned_cells": b.get("unassigned", 0)},
                "embedded_node_annotations": block_annotations_by_donor.get(donor, []),
                "modality": "protein_if", "acquisition": runs_row.get("acq_run", "unresolved"),
                "software": runs_row.get("software", batch_row.get("acquisition_software", "unresolved")),
                "panel": q.get("panel_id", "unresolved"),
                "geometry_revision": q.get("segmentation_version", "unresolved"),
                "coverage_status": "embedded pLN annotations, same/different-block assignment, native cell typing, and dedicated-pLN shared-marker comparison available; panel-specific analysis retained separately",
                "evidence_source": str(args.embedded_block_ledger), "code_revision": revision,
            })

    for r in node_rows:
        donor = key(r["donor_id"])
        image_file = r["image_file"]
        ppath = pprotein_root / f"donor_id={r['donor_id']}" / f"roi={r['roi']}" / "data.parquet"
        typed = ptyped_root / f"donor_id={r['donor_id']}" / f"roi={r['roi']}" / "data.parquet"
        raw_info = parquet_info(ppath)
        excluded = (str(r.get("excluded", "")).strip().lower() == "true"
                    or donor in {key(d) for d in DONOR_EXCLUSIONS})
        if raw_info["exists"]:
            raw_info["observed_marker_columns"] = [c for c in raw_info["columns"]
                                                   if c not in {"donor_id", "roi", "object_id"}]
        rows.append({
            "donor_id": str(r["donor_id"]), "tissue": "dedicated_pLN", "slide_id": image_file,
            "excluded": excluded, "exclusion_reason": DONOR_EXCLUSIONS.get(donor, ""),
            "physical_block_id": "unresolved", "block_relationship_to_pancreas": "unresolved",
            "modality": "protein_if", "acquisition": "date/run unavailable; exact QuPath image ID recorded",
            "software": "QuPath; acquisition software/version unavailable",
            "panel": "IO60 pLN panel; exact slide-specific measured marker columns listed in current_pLN_cells",
            "node_annotation_classes": r.get("node_classes", "unresolved"),
            "node_annotation_count": r.get("n_nodes", "unresolved"),
            "geometry_revision": f"IO60pLN export; QuPath entry {r.get('qupath_entry','unresolved')}",
            "source_path": pprotein_source, "source_sha256": "",
            "current_measurement": str(ppath),
            "current_pLN_cells": raw_info, "current_pLN_typing": parquet_info(typed),
            "coverage_status": ("explicitly excluded; source annotation retained for inventory" if excluded
                                else "current dedicated pLN cell and native phenotype partitions available"),
            "evidence_source": str(INTEGRATION / "pln_node_annotations.csv"), "code_revision": revision,
        })

    for r in xenium_rows:
        donor, roi, source = key(r["donor_id"]), r["roi"].strip(), r.get("source_path", "").strip()
        match = re.search(r"output-XETG\d+__(\d+)__", source)
        serial = match.group(1) if match else ""
        mode = mode_by_key.get((donor, roi, serial), {})
        section_relationship = mode.get("mode", "").strip() or fallback_relationship or "unresolved"
        section_relationship_evidence = mode.get("evidence", "").strip() or (
            f"{integration_config} [integration] mode={fallback_relationship}; "
            "no exact donor/ROI/serial override")
        local_path = local_xenium(source)
        local = local_path is not None
        tissue = "panc" if roi == "panc" else "pln"
        candidates = [s for s in rna_sections
                      if key(re.sub(r"^HDL", "", str(s.get("donor_id", "")),
                                    flags=re.IGNORECASE)) == donor
                      and str(s.get("tissue", "")).lower() == tissue
                      and (str(s.get("source_roi", "")).lower() == roi.lower()
                           or (not s.get("source_roi") and
                               key(re.sub(r"^HDL", "", str(s.get("section_key", "")),
                                          flags=re.IGNORECASE)) == donor))]
        donor_acquisitions = [x for x in xenium_rows if key(x.get("donor_id")) == donor
                              and str(x.get("roi", "")).strip().lower() == roi.lower()]
        prepared_matches = [s for s in candidates
                            if s.get("source_serial") and s.get("source_serial") == serial]
        if prepared_matches:
            prepared = prepared_matches[0] if len(prepared_matches) == 1 else None
        elif (len(candidates) == 1 and len(donor_acquisitions) == 1
              and (not serial or not candidates[0].get("source_serial"))):
            # A source with no serial can be associated only when the manifest declares one
            # acquisition for that donor/ROI. Never copy one prepared h5ad onto both a serial
            # section and a separate same-slide acquisition, or use a known-different serial.
            prepared = candidates[0]
        else:
            prepared = None
        ambiguous_candidates = candidates if candidates and prepared is None else []
        exact_prepared = bool(prepared and serial and prepared.get("source_serial") == serial)
        if prepared and prepared.get("h5ad_source_kind") == "vendor_bundle_gene_expression_matrix":
            rna_coverage = ("exact-serial local bundle matrix and cells prepared for tissue-specific QC and typing"
                           if exact_prepared else
                           "prepared local bundle matrix selected by a single-candidate/single-acquisition fallback; no exact-serial claim")
        elif prepared:
            rna_coverage = ("exact-serial prepared expression h5ad available; no complete local bundle matrix"
                            if exact_prepared else
                            "prepared expression h5ad available by a unique candidate; exact acquisition is unresolved")
        elif ambiguous_candidates:
            local_note = "; bundle directory is local" if local else ""
            rna_coverage = ("prepared h5ad exists, but its acquisition is ambiguous" + local_note +
                            "; no same-cell or serial match is claimed")
        elif local:
            rna_coverage = "bundle directory local but no complete expression matrix was prepared"
        else:
            rna_coverage = "recorded remotely or unavailable locally; no RNA call inferred"
        excluded = donor in excluded_donors
        if excluded:
            rna_coverage = "explicitly excluded by the central cohort rule; no RNA typing or analysis admitted"
        qc_candidates = []
        for section in candidates:
            qc_key = (donor, str(section.get("section_key", "")))
            report = ((panc_rna_qc if tissue == "panc" else pln_rna_qc).get(qc_key, {}))
            if report:
                qc_candidates.append({
                    "section_key": qc_key[1],
                    "h5ad": section.get("h5ad", ""),
                    "source_kind": section.get("h5ad_source_kind", ""),
                    "cell_geometry_sha256": section.get("cells_parquet_sha256", ""),
                    "panel_gene_count": section.get("n_genes"),
                    "qc": qc_summary(report, tissue),
                })
        if qc_candidates:
            panel_description = (f"Xenium 5K; candidate prepared matrix has "
                                 f"{qc_candidates[0].get('panel_gene_count')} genes")
        elif local:
            panel_description = "local Xenium bundle; exact prepared panel/QC unavailable"
        else:
            panel_description = "remote or unavailable; panel not locally verified"
        geometry_description = (
            "bundle cells.parquet cell geometry; source hash recorded per prepared candidate"
            if qc_candidates else "bundle-specific cell geometry not locally verified")
        rows.append({
            "donor_id": str(r["donor_id"]), "tissue": "pancreas" if roi == "panc" else "pLN",
            "disease_status": disease_by_donor.get(donor, "unavailable"),
            "excluded": excluded, "exclusion_reason": DONOR_EXCLUSIONS.get(donor, ""),
            "slide_id": f"Xenium serial {serial or 'unresolved'}", "physical_block_id": "unresolved",
            "block_relationship_to_pancreas": "unresolved_from_local_block_id_inventory",
            "modality": "rna_xenium", "acquisition": serial or "unresolved",
            "software": "Xenium", "panel": panel_description,
            "geometry_revision": geometry_description,
            "source_path": source, "local_source_available": local,
            "local_source_path": local_path or "",
            "section_relationship": section_relationship,
            "section_relationship_evidence": section_relationship_evidence,
            "prepared_h5ad": prepared.get("h5ad", "") if prepared else "",
            "prepared_h5ad_source_kind": prepared.get("h5ad_source_kind", "unavailable")
                                         if prepared else "unavailable_or_ambiguous",
            "prepared_h5ad_cells": prepared.get("n_cells") if prepared else None,
            "prepared_h5ad_genes": prepared.get("n_genes") if prepared else None,
            "prepared_match_is_exact": exact_prepared,
            "prepared_resolution_method": (
                "exact_serial_match" if exact_prepared else
                "unique_candidate_fallback" if prepared else
                "ambiguous_candidates" if ambiguous_candidates else "no_local_match"),
            "prepared_section_qc_candidates": qc_candidates,
            "nucleus_count_measured": prepared.get("nucleus_count_measured") if prepared else False,
            "nucleus_count_zero_cells": prepared.get("nucleus_count_zero_cells") if prepared else None,
            "median_transcripts_per_cell": prepared.get("median_nCount_Xenium") if prepared else None,
            "median_genes_per_cell": prepared.get("median_nFeature_Xenium") if prepared else None,
            "bundle_matrix_sha256": prepared.get("cell_feature_matrix_sha256", "") if prepared else "",
            "bundle_cells_sha256": prepared.get("cells_parquet_sha256", "") if prepared else "",
            "ambiguous_prepared_h5ad_candidates": [
                {"h5ad": s.get("h5ad", ""), "section_key": s.get("section_key", ""),
                 "source_kind": s.get("h5ad_source_kind", "")}
                for s in ambiguous_candidates],
            "rna_input_inventory": str(args.rna_inventory),
            "coverage_status": rna_coverage,
            "evidence_source": str(INTEGRATION / "xenium_paths.csv"), "code_revision": revision,
        })

    # The post-Xenium IF lane is a distinct exact same-section pair: it uses Xenium's cell
    # polygons and an independently gated transform. Keep it separate from the 70-plex
    # PhenoCycler section relationship declared in section_modes.csv.
    postx_export = json.loads(args.postxenium_export.read_text()) \
        if args.postxenium_export.is_file() else {"donors": []}
    postx_coverage = read_csv(args.postxenium_coverage) \
        if args.postxenium_coverage.is_file() else []
    postx_coverage_by_serial = {str(r.get("xenium_serial", "")): r for r in postx_coverage}
    postx_measurements = json.loads(args.postxenium_measurements.read_text()) \
        if args.postxenium_measurements.is_file() else []
    postx_measurement_by_serial = {str(r.get("xenium_serial", "")): r for r in postx_measurements}
    channel_report = json.loads(args.postxenium_channel_assignment.read_text()) \
        if args.postxenium_channel_assignment.is_file() else {}
    channel_by_serial = {str(r.get("xenium_serial", "")): r
                         for r in channel_report.get("resolved", [])}
    pancreas_donor_by_serial: dict[str, str] = {}
    for r in xenium_rows:
        if str(r.get("roi", "")).strip().lower() != "panc":
            continue
        source = str(r.get("source_path", ""))
        match = re.search(r"output-XETG\d+__(\d+)__", source)
        if not match:
            continue
        serial = match.group(1)
        donor = str(r.get("donor_id", ""))
        if serial in pancreas_donor_by_serial and pancreas_donor_by_serial[serial] != donor:
            raise ValueError(f"Xenium serial {serial} maps to multiple pancreas donors")
        pancreas_donor_by_serial[serial] = donor
    disease_by_exact_donor = {key(r.get("donor_id")): str(r.get("disease_status", "unavailable"))
                              for r in rows if r.get("tissue") == "pancreas"
                              and r.get("modality") == "rna_xenium"}
    for section in postx_export.get("donors", []):
        serial = str(section.get("serial", ""))
        registration_donor = str(section.get("donor_id", ""))
        if key(registration_donor) in {key(d) for d in DONOR_EXCLUSIONS}:
            # Keep the user exclusion intact even if an archived image/export reappears.
            continue
        donor = pancreas_donor_by_serial.get(serial, registration_donor)
        coverage = postx_coverage_by_serial.get(serial, {})
        measured_pair = postx_measurement_by_serial.get(serial, {})
        channels = channel_by_serial.get(serial, {})
        if coverage:
            state = coverage.get("state", "unresolved")
        elif measured_pair:
            state = "cell_measurement_written; exact pair coverage manifest pending"
        elif str(section.get("state", "")).lower() == "complete":
            state = "source_available; registration_or_measurement_unresolved"
        else:
            state = "source_unavailable"
        diagnostics_raw = coverage.get("registration_diagnostics", "")
        if isinstance(diagnostics_raw, str) and diagnostics_raw:
            try:
                diagnostics = json.loads(diagnostics_raw)
            except json.JSONDecodeError:
                diagnostics = [{"state": "unparsed_csv_diagnostics", "raw": diagnostics_raw}]
        else:
            diagnostics = diagnostics_raw or []
        rows.append({
            "donor_id": donor, "disease_status": disease_by_exact_donor.get(key(donor), "unavailable"),
            "tissue": "pancreas", "slide_id": f"Xenium serial {serial}",
            "physical_block_id": "unresolved",
            "block_relationship_to_pancreas": "same physical Xenium section; not the separate PhenoCycler serial-section pair",
            "modality": "postxenium_if", "acquisition": "post-Xenium IF of exact Xenium section",
            "software": "QuPath channel names + polygon-integrated IF; exact registration-gated producer",
            "panel": "INS, GCG; source channel names recorded per donor",
            "geometry_revision": "Xenium cell_boundaries.parquet with exact PASS transform; see registration fields",
            "source_path": section.get("qptiff", ""), "source_serial": serial,
            "xenium_bundle": section.get("xenium_bundle", ""),
            "qptiff_pixel_size_um": section.get("if_pixel_um"),
            "xenium_pixel_size_um": section.get("xenium_um_per_px"),
            "if_image_origin_um": section.get("if_origin_um"),
            "xenium_image_origin_um": section.get("xenium_origin_um"),
            "pixel_size_source": section.get("pixel_size_source", "unresolved"),
            "section_relationship": "same_slide",
            "section_relationship_evidence": f"{args.postxenium_export}: exact Xenium serial {serial}",
            "channel_source": measured_pair.get("channel_source", "unresolved"),
            "qupath_entry": measured_pair.get("qupath_entry", "unresolved"),
            "channel_assignment": channels.get("channel_names", "unresolved"),
            "channel_indices": channels.get("channels", "unresolved"),
            "qptiff_channel_names": channels.get("qptiff_channel_names", []),
            "channel_assignment_conflict": channels.get("csv_conflict", []),
            "channel_assignment_resolution": ("QuPath project channel names take precedence" if channels else "unresolved"),
            "registration_transform": coverage.get("registration_transform", measured_pair.get("transform_source", "unresolved")),
            "registration_gate": coverage.get("registration_gate", measured_pair.get("alignment_gate", "unresolved")),
            "registration_transform_sha256": coverage.get("registration_transform_sha256", "unresolved"),
            # This is the row count with a verified same-section cell identity (including cells
            # whose protein vector is unavailable). Keep it distinct from how many cells have a
            # finite INS measurement.
            "paired_cell_count": (count_value(coverage.get("n_xenium_cells"))
                                  if state == "paired_cell_features_measured" else None),
            "protein_INS_measured_cell_count": count_value(
                coverage.get("n_cells_with_INS_measurement")),
            "source_cell_count": count_value(coverage.get(
                "n_xenium_cells", measured_pair.get("n_cells", None))),
            "registration_diagnostics": diagnostics,
            "coverage_status": state,
            "reason": coverage.get("reason", ""),
            "provenance_source": str(args.postxenium_export),
            "measurement_sidecar": str(args.postxenium_measurements),
            "coverage_evidence": str(args.postxenium_coverage),
            "code_revision": revision,
        })

    # The exclusion is recorded because coverage reporting must distinguish it from a missing run.
    for excluded_donor, reason in DONOR_EXCLUSIONS.items():
        rows.append({
            "donor_id": excluded_donor, "tissue": "pancreas_and_pLN",
            "slide_id": "historical acquisition; retained for exclusion audit",
            "physical_block_id": "unresolved", "block_relationship_to_pancreas": "unresolved",
            "modality": "all", "acquisition": "not run", "software": "not run",
            "panel": "not run", "geometry_revision": "not admitted to active typing exports",
            "coverage_status": "explicitly excluded by the central cohort rule; not missing or negative",
            "exclusion_reason": reason, "evidence_source": "phenocycler.cohort.DONOR_EXCLUSIONS",
            "code_revision": revision,
        })

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "title": "Current input and coverage manifest",
        "generated_local": "2026-10-03",
        "active_pancreas_donors": sorted(str(r["donor_id"]) for r in images.values()),
        "active_pancreas_donor_count": len(images),
        "dedicated_pln_slide_count": len({key(r["donor_id"]) for r in node_rows
                                           if str(r.get("excluded", "")).strip().lower() != "true"}),
        "excluded_donors": DONOR_EXCLUSIONS,
        "qupath_manifest": {"path": str(qpath), "content_id": qman.get("content_id"),
                            "manifest_version": qman.get("manifest_version")},
        "embedded_pln_block_manifest": {"path": str(args.embedded_block_ledger),
                                        "annotations_path": str(args.block_annotations),
                                        "annotations_sha256": (hashlib.sha256(args.block_annotations.read_bytes()).hexdigest()
                                                               if args.block_annotations.is_file() else None),
                                        "donors_with_nodes": sorted(embedded),
                                        "donors_without_nodes": block_manifest.get("donors_without_pln_cells", [])},
        "xenium_analysis_releases": {
            "pancreas": {"path": str(args.pancreas_rna_release),
                         "release_id": panc_rna_release.get("release_id"),
                         "n_donors": panc_rna_release.get("n_donors"),
                         "n_sections": panc_rna_release.get("n_sections"),
                         "n_cells_total": panc_rna_release.get("n_cells_total"),
                         "donor_qc_sections": len(panc_rna_qc)},
            "pLN": {"path": str(args.pln_rna_release),
                    "release_id": pln_rna_release.get("release_id"),
                    "sections_requested": pln_rna_release.get("scope", {}).get("n_requested"),
                    "sections_available": pln_rna_release.get("scope", {}).get("n_available"),
                    "section_reports": len(pln_rna_qc),
                    "section_reports_path": str(pln_reports_path)},
        },
        "postxenium_same_section_coverage": {
            "source_sections": sum(r.get("modality") == "postxenium_if" for r in rows),
            "paired_cell_sections": sum(r.get("modality") == "postxenium_if"
                                         and r.get("coverage_status") == "paired_cell_features_measured"
                                         for r in rows),
            "states": {str(s): sum(r.get("modality") == "postxenium_if"
                                   and r.get("coverage_status") == s for r in rows)
                       for s in sorted({str(r.get("coverage_status")) for r in rows
                                        if r.get("modality") == "postxenium_if"})},
            "coverage_manifest": str(args.postxenium_coverage),
            "producer_measurements": str(args.postxenium_measurements),
            "channel_assignment": str(args.postxenium_channel_assignment),
            "relationship_note": "These are exact same-section post-Xenium IF measurements; they do not replace separate PhenoCycler 70-plex sequential-slide or declared 70-plex same-slide pairs.",
        },
        "section_modes_path": str(INTEGRATION / "section_modes.csv"),
        "redsea_review_status": "8-arm paired comparison and per-marker dossier complete; provisional recommendation written; production policy awaits user ruling",
        "app_development": "halted per user instruction",
        "code_revision": revision,
        "worktree_modified_paths": status,
        "unresolved_metadata": [
            "exact physical block IDs are unavailable for dedicated pLN slides; embedded nodes carry annotation-level same/different-block relationships and pancreas-region notes",
            "exact dedicated pLN acquisition date/run and physical block IDs are unavailable; slide IDs and observed marker columns are recorded",
            "Xenium source bundles are remote/unavailable locally unless local_source_available is true",
            "Xenium panel, geometry, and QC provenance is attached only to matched prepared sections or explicitly listed ambiguous candidates; remote/unprepared bundles remain unavailable",
        ],
        "records": rows,
    }
    json_path = out / "current_input_manifest.json"
    json_path.write_text(json.dumps(payload, indent=2))
    flat = []
    for row in rows:
        flat.append({k: (json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v)
                     for k, v in row.items()})
    columns = sorted(set().union(*(r.keys() for r in flat)))
    with (out / "current_input_manifest.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat)
    print(f"{len(images)} active pancreas donors; {len(node_rows)} dedicated pLN slides; "
          f"{len(embedded)} embedded-node donors; {len(xenium_rows)} recorded Xenium rows")
    print(f"-> {json_path}")
    print(f"-> {out / 'current_input_manifest.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
