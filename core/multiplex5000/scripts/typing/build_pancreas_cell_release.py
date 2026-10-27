#!/usr/bin/env python3
"""Assemble a provenance-complete per-cell pancreas protein release from one frozen arm.

The source geometry table is authoritative for row coverage. Native calls are left-joined so
untyped/QC-excluded objects remain visible with their reason. The REDSEA measurements are kept
beside the original QuPath measurements; hormone values used for typing are masked by the same
fragment-overlap sidecar as the typing producer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from phenocycler.cohort import DONOR_EXCLUSIONS

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(HERE))
import subpanel_io as S  # noqa: E402


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def donor_key(value: object) -> str:
    text = str(value or "").strip()
    return text.lstrip("0") or "0"


def check_unique(frame: pd.DataFrame, path: Path, name: str) -> None:
    if "object_id" not in frame:
        raise ValueError(f"{path}: {name} has no object_id")
    if frame.object_id.isna().any() or frame.object_id.duplicated().any():
        raise ValueError(f"{path}: {name} has missing or duplicate object_id values")


def merge_one_to_one(left: pd.DataFrame, right: pd.DataFrame, donor: str,
                     source: str) -> pd.DataFrame:
    check_unique(right, Path(source), source)
    overlap = (set(left.columns) & set(right.columns)) - {"object_id"}
    if overlap:
        raise ValueError(f"{donor}: {source} would overwrite columns {sorted(overlap)}")
    return left.merge(right, on="object_id", how="left", validate="one_to_one")


def merge_reconciled(left: pd.DataFrame, right: pd.DataFrame, donor: str,
                     source: str) -> pd.DataFrame:
    """Join a producer table, verifying duplicate answers before keeping one copy."""
    check_unique(right, Path(source), source)
    overlap = sorted((set(left.columns) & set(right.columns)) - {"object_id"})
    if overlap:
        index = left.set_index("object_id")
        side = right.set_index("object_id")
        common_ids = side.index.intersection(index.index)
        for column in overlap:
            a, b = index.loc[common_ids, column], side.loc[common_ids, column]
            both = a.notna() & b.notna()
            if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                equal = np.isclose(a[both].to_numpy(dtype=float), b[both].to_numpy(dtype=float),
                                   rtol=0, atol=0, equal_nan=True).all()
            else:
                equal = a[both].astype(str).to_numpy().tolist() == b[both].astype(str).to_numpy().tolist()
            if not equal:
                raise ValueError(f"{donor}: {source} disagrees with existing {column} values")
    keep = ["object_id"] + [c for c in right.columns if c not in overlap and c != "object_id"]
    return left.merge(right[keep], on="object_id", how="left", validate="one_to_one")


def marker_frame(aligned: pd.DataFrame, donor: str, scale_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    check_unique(aligned, Path("aligned REDSEA inputs"), "aligned REDSEA inputs")
    markers = [c for c in aligned.columns if c not in {
        "donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
        "cell_region", "islet_num", "in_alignment_block"}]
    values = aligned[["object_id"] + markers].copy()
    masked = S.mask_fragment_hormones(values, donor)
    scale = pd.read_parquet(scale_path).set_index("marker")
    available = [m for m in markers if m in scale.index]
    columns: dict[str, object] = {"object_id": aligned.object_id.astype(str).to_numpy()}
    if available:
        scaled_values = S.apply_scale(masked[available].to_numpy(dtype=float), available, scale)
        for j, marker in enumerate(available):
            if marker not in {"DAPI", "IAPP"}:
                columns[f"typing_scaled_{marker}"] = scaled_values[:, j]
    for marker in markers:
        # Keep the pre-marker-scaling REDSEA output; this is the value used by the structural
        # highest-marker rule and it retains the calibration scale's interpretable units.
        columns[f"redsea_aligned_{marker}"] = aligned[marker].to_numpy()
        if marker not in scale.index and marker not in {"DAPI", "IAPP"}:
            columns[f"typing_scaled_{marker}"] = np.full(len(aligned), np.nan, dtype=np.float32)
    frac = S.fragment_coverage(aligned.object_id, donor).to_numpy(dtype=float)
    columns["hormone_fragment_overlap_fraction"] = frac
    columns["hormone_measurement_state"] = np.select(
        [~np.isfinite(frac), frac > 0],
        ["overlap_measurement_unavailable", "masked_fragment_overlap"],
        default="usable_for_typing",
    )
    scaled = pd.DataFrame(columns)
    scale_table = scale.reset_index()[[c for c in ("marker", "floor", "top", "span")
                                       if c in scale.reset_index().columns]].copy()
    scale_table["donor_id"] = donor
    return scaled, scale_table


def build_donor(donor: str, arm_root: Path, run: dict, disease: str,
                out_dir: Path) -> dict:
    geometry_root = Path(run["geometry_root"])
    aligned_root = Path(run["aligned_root"])
    scale_path = Path(run["scale"])
    raw_path = geometry_root / "revised_cells" / f"donor_id={donor}" / "data.parquet"
    aligned_path = aligned_root / f"donor_id={donor}" / "data.parquet"
    call_path = arm_root / "cells" / f"donor_id={donor}" / "data.parquet"
    immune_path = arm_root / "immune" / f"donor_id={donor}" / "data.parquet"
    fundamentals_path = geometry_root / "qupath_fundamentals" / f"donor_id={donor}" / "data.parquet"
    crossed_path = geometry_root / "crossed" / f"crossed_{donor}.parquet"
    for path in (raw_path, aligned_path, call_path, immune_path, fundamentals_path, crossed_path):
        if not path.is_file():
            raise FileNotFoundError(f"{donor}: required release input missing: {path}")

    raw = pd.read_parquet(raw_path)
    check_unique(raw, raw_path, "revised QuPath cells")
    if "islet_num" in raw:
        # The revised QuPath export's source islet number and the typing chain's geometric
        # islet membership are distinct answers. Preserve both; the typed membership remains
        # `islet_num`, while this original export value is explicitly named.
        raw = raw.rename(columns={"islet_num": "source_parent_islet_num"})
    raw_markers = [c for c in raw.columns if c not in {
        "donor_id", "object_id", "image", "parent_raw", "qupath_classification",
        "X_centroid", "Y_centroid", "nucleus_area_um2", "cell_area", "cell_solidity",
        "nucleus_cell_area_ratio", "signed_dist_islet_um", "islet_num", "cell_region",
        "qc_analysis_eligible", "qc_estimation_eligible", "analysis_reason",
        "estimation_reason", "in_alignment_block"}]
    raw = raw.rename(columns={marker: f"raw_{marker}" for marker in raw_markers})
    raw["tissue"] = "pancreas"
    raw["modality"] = "protein_if"
    raw["disease_status"] = disease
    raw["donor_id"] = donor

    # A corrected generation can add mandatory admission exclusions to its
    # geometry QC. Keep the authoritative original QC answers beside those
    # current decisions rather than losing the original failure reasons.
    original_qc_provenance = None
    if run.get("original_geometry_root"):
        original_path = Path(run["original_geometry_root"]) / "revised_cells" / f"donor_id={donor}" / "data.parquet"
        fields = ["qc_analysis_eligible", "qc_estimation_eligible", "analysis_reason", "estimation_reason"]
        original = pd.read_parquet(original_path, columns=["object_id", *fields])
        check_unique(original, original_path, "original geometry QC provenance")
        if len(original) != len(raw) or set(original.object_id) != set(raw.object_id):
            raise ValueError(f"{donor}: original/current geometry QC universes differ")
        original = original.rename(columns={c: "original_geometry_" + c for c in fields})
        raw = merge_one_to_one(raw, original, donor, "original geometry QC provenance")
        previous = raw.original_geometry_qc_analysis_eligible.fillna(False).astype(bool)
        current = raw.qc_analysis_eligible.fillna(False).astype(bool)
        original_qc_provenance = {
            "path": str(original_path.resolve()), "sha256": digest(original_path), "rows": len(original),
            "fields": {field: "original_geometry_" + field for field in fields},
            "original_analysis_eligible_rows": int(previous.sum()), "current_analysis_eligible_rows": int(current.sum()),
            "old_eligible_now_ineligible_rows": int((previous & ~current).sum()),
            "old_ineligible_now_eligible_rows": int((~previous & current).sum()),
            "already_ineligible_current_reason_changed_rows": int((~previous &
                raw.original_geometry_analysis_reason.fillna("").ne(raw.analysis_reason.fillna(""))).sum()),
            "join": "Full exact object_id universe; original and current flags/reasons remain distinct",
        }

    calls = pd.read_parquet(call_path)
    check_unique(calls, call_path, "primary and refinement calls")
    called_ids = set(calls.object_id.astype(str))
    eligible_ids = set(raw.loc[raw.qc_analysis_eligible.fillna(False), "object_id"].astype(str))
    if called_ids != eligible_ids:
        raise ValueError(f"{donor}: primary call keys differ from analysis-eligible geometry; "
                         f"missing={len(eligible_ids-called_ids)}, extra={len(called_ids-eligible_ids)}")
    frame = raw
    frame = merge_reconciled(frame, calls, donor, "primary/refined calls")
    immune = pd.read_parquet(immune_path)
    check_unique(immune, immune_path, "immune subtype calls")
    immune_columns = [c for c in immune.columns if c not in {"donor_id", "roi"}]
    frame = merge_one_to_one(frame, immune[immune_columns], donor, "immune subtype calls")
    fundamentals = pd.read_parquet(fundamentals_path)
    check_unique(fundamentals, fundamentals_path, "QuPath cell fundamentals")
    fundamentals = fundamentals.rename(columns={c: f"fundamentals_{c}"
                                                for c in fundamentals.columns if c != "object_id"})
    frame = merge_one_to_one(frame, fundamentals, donor, "QuPath cell fundamentals")
    crossed = pd.read_parquet(crossed_path)
    check_unique(crossed, crossed_path, "measured cell/object overlap fractions")
    overlap_columns = [c for c in crossed.columns if c not in {"object_id", "donor_id"}]
    frame = merge_one_to_one(frame, crossed[["object_id"] + overlap_columns], donor,
                             "measured cell/object overlap fractions")

    aligned = pd.read_parquet(aligned_path)
    calibrated, scale_table = marker_frame(aligned, donor, scale_path)
    # Raw source covers every geometry row; calibrated inputs can miss true no-measurement cells.
    frame = merge_one_to_one(frame, calibrated, donor, "masked aligned REDSEA inputs")
    frame = frame.copy()  # compact blocks once before adding per-row release status columns
    frame["primary_call_state"] = np.where(frame["cell_class"].notna(), "called",
                                            np.where(frame["qc_analysis_eligible"].fillna(False),
                                                     "eligible_without_call", "not_analysis_eligible"))
    frame["primary_call_unavailability_reason"] = np.where(
        frame["primary_call_state"].eq("not_analysis_eligible"),
        frame["analysis_reason"].fillna("analysis_eligibility_unavailable"),
        np.where(frame["primary_call_state"].eq("eligible_without_call"),
                 "no_primary_call", ""))
    frame["refinement_state"] = np.where(
        frame["refine_pass"].eq(True), "refinement_passed",
        np.where(frame["refine_pass"].eq(False), "refinement_not_passed",
                 "refinement_not_evaluated"))
    frame["immune_subtyping_state"] = np.where(
        frame["immune_class"].notna() | frame["immune_subtype"].notna(),
        "immune_subtype_call_available", "no_immune_subtype_call")
    frame["cell_id"] = frame.object_id.astype(str)
    if frame.duplicated(["donor_id", "cell_id"]).any():
        raise ValueError(f"{donor}: release contains duplicate donor/cell keys")
    if "cell_region_20um" in frame or any("20um" in c.lower() for c in frame.columns):
        raise ValueError(f"{donor}: obsolete expansion-derived feature entered the release")

    out = out_dir / f"donor_id={donor}"
    out.mkdir(parents=True, exist_ok=True)
    destination = out / "data.parquet"
    frame.to_parquet(destination, index=False, compression="zstd")
    scale_table.to_csv(out / "marker_scale.csv", index=False)
    return {
        "donor_id": donor, "disease_status": disease,
        "source_rows": len(raw), "release_rows": len(frame),
        "analysis_eligible_rows": int(frame.qc_analysis_eligible.fillna(False).sum()),
        "primary_called_rows": int(frame.primary_call_state.eq("called").sum()),
        "no_primary_call_rows": int(frame.primary_call_state.eq("eligible_without_call").sum()),
        "not_analysis_eligible_rows": int(frame.primary_call_state.eq("not_analysis_eligible").sum()),
        "fragment_overlap_masked_rows": int(frame.hormone_measurement_state.eq("masked_fragment_overlap").sum()),
        "fragment_overlap_unavailable_rows": int(frame.hormone_measurement_state.eq("overlap_measurement_unavailable").sum()),
        "cell_table": str(destination.resolve()), "cell_table_rows": pq.ParquetFile(destination).metadata.num_rows,
        "cell_table_bytes": destination.stat().st_size,
        "raw_input_sha256": digest(raw_path), "aligned_input_sha256": digest(aligned_path),
        "calls_input_sha256": digest(call_path), "immune_input_sha256": digest(immune_path),
        "overlap_input_sha256": digest(crossed_path), "fundamentals_input_sha256": digest(fundamentals_path),
        "original_geometry_QC_provenance": original_qc_provenance,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--input-manifest", type=Path,
                    default=REPO / "results/v3/release_20261003/current_input_manifest.json")
    ap.add_argument("--ruling-json", type=Path, default=None,
                    help="Required with --production; records the user-selected policy.")
    ap.add_argument("--production", action="store_true")
    ap.add_argument("--donors", nargs="*", default=None,
                    help="Optional donor subset for a validation build; cannot be labelled production.")
    args = ap.parse_args()
    arm_root, out = args.arm_root.resolve(), args.out.resolve()
    run_path = arm_root / "run.json"
    run = json.loads(run_path.read_text())
    if run.get("status") != "complete":
        raise ValueError(f"{arm_root.name}: calibration arm is not complete")
    if args.production and (not args.ruling_json or not args.ruling_json.is_file()):
        raise ValueError("production release requires the recorded user production ruling")
    if args.production and args.donors:
        raise ValueError("a partial donor selection cannot be labelled a production release")
    ruling = json.loads(args.ruling_json.read_text()) if args.ruling_json else None
    manifest = json.loads(args.input_manifest.read_text())
    disease = {donor_key(r.get("donor_id")): str(r.get("disease_status", "unavailable"))
               for r in manifest.get("records", [])
               if r.get("tissue") == "pancreas" and r.get("modality") == "rna_xenium"}
    donors = list(map(str, run["donors"]))
    expected_donors = {donor_key(d) for d in manifest.get("active_pancreas_donors", [])}
    observed_donors = {donor_key(d) for d in donors}
    excluded_donors = {donor_key(d) for d in DONOR_EXCLUSIONS}
    if observed_donors & excluded_donors or observed_donors != expected_donors:
        raise ValueError(f"unexpected pancreas donor universe in arm {arm_root.name}: {donors}")
    selected = [d for d in donors if not args.donors or donor_key(d) in
                {donor_key(x) for x in args.donors}]
    if args.donors and {donor_key(x) for x in args.donors} != {donor_key(x) for x in selected}:
        raise ValueError(f"requested donor subset is not present in {arm_root.name}")
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for donor in selected:
        if donor_key(donor) not in disease:
            raise ValueError(f"{donor}: disease status is absent from the current input manifest")
        summaries.append(build_donor(donor, arm_root, run, disease[donor_key(donor)], out))

    keys = {(r["donor_id"], r["cell_table_rows"]) for r in summaries}
    if len(keys) != len(summaries):
        raise ValueError("duplicate donor/row-count key in assembled partitions")
    production = bool(args.production)
    revision = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              check=True, capture_output=True, text=True).stdout.strip()
    payload = {
        "title": "Pancreas PhenoCycler cell-level release",
        "status": "production" if production else "validation_build",
        "production_policy": ruling if ruling else "not_promoted; production ruling pending",
        "calibration_arm": arm_root.name,
        "calibration_run": str(run_path), "calibration_run_sha256": digest(run_path),
        "calibration_source_scripts": run.get("source_scripts", {}),
        "input_manifest": str(args.input_manifest.resolve()),
        "input_manifest_sha256": digest(args.input_manifest),
        "release_builder": str(Path(__file__).resolve()),
        "release_builder_sha256": digest(Path(__file__).resolve()),
        "production_ruling_sha256": digest(args.ruling_json) if args.ruling_json else None,
        "number_of_donors": len(summaries), "donors": summaries,
        "excluded_donors": DONOR_EXCLUSIONS,
        "row_key": ["donor_id", "cell_id"],
        "source_rows": "Every row in the authoritative revised QuPath cell geometry export is retained; calls are left-joined.",
        "columns": {
            "raw_*": "Original QuPath measurements; original object class, geometry, signed nearest-islet distance, and QC reasons are retained.",
            "original_geometry_*": "Authoritative geometry QC flags/reasons before the current generation adds keyed admission exclusions; current qc_* and reason fields remain separate.",
            "redsea_aligned_*": "Selected-arm aligned marker values before per-marker scale; marker units are the producer's aligned values.",
            "typing_scaled_*": "Selected-arm marker-scaled values. Hormones are null for any positive measured hormone-fragment overlap; missing overlap evidence remains unavailable.",
            "cell_class": "Primary identity call.", "refined_class": "Separate refinement result; does not overwrite cell_class.",
            "immune_class/immune_subtype/macrophage_profile": "Immune typing sidecar from the same hierarchy; missing profile markers remain unavailable.",
            "crossed_frac_*": "Measured fraction of the cell body overlapped by the named structural object; no expansion-ring feature is included.",
            "signed_dist_islet_um": "Continuous signed distance to the nearest true islet object in microns.",
        },
        "spatial_units": "microns for coordinates and signed distances; square microns for measured areas",
        "no_expansion_ring": True,
        "diagnostic_field_validity": {
            "top_identity_marker/argmax_scaled_no_supporting_only": "Historical cohort-floor/span-scaled marker argmax diagnostic. The carried FOXP3 dominance is a scale artifact; this is not cell identity, Treg biology, marker function or biological accuracy evidence.",
            "top_identity_score": "If present, the scaled-marker diagnostic score is not identity confidence or biological accuracy evidence.",
            "identity_fields": ["cell_class", "refined_class", "phenotype", "refined_phenotype", "immune_class", "immune_subtype", "macrophage_profile"],
            "top_expressed_marker": "Expressed-value comparison used in the declared primary structural rule; separately retained from the scaled argmax diagnostic and not independent accuracy evidence.",
        },
        "code_revision": revision,
        "inputs": {"geometry_root": run["geometry_root"], "aligned_root": run["aligned_root"],
                   "scale": run["scale"], "ruling_json": str(args.ruling_json) if args.ruling_json else None},
    }
    (out / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    pd.DataFrame(summaries).to_csv(out / "donor_coverage.csv", index=False)
    print(f"wrote {len(summaries)} donor partitions to {out} ({payload['status']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
