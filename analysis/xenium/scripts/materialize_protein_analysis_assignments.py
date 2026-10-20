#!/usr/bin/env python3
"""Write bounded-memory current identity assignments with explicit exclusions.

Measurements remain in the indexed native release. These alternative assignment
partitions retain every source cell, native identities, object membership and QC.
Historical mode applies its 25-cell overlay; current mode preserves the mandatory
exact-key registry independently of the new source's QC eligibility.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))
from protein_analysis_admission import DEFAULT_EXCLUSIONS, load_protein_exclusions

FIELDS = ["cell_id", "object_id", "donor_id", "cell_class", "refined_class", "phenotype",
          "refined_phenotype", "immune_class", "immune_subtype", "immune_leaf",
          "qc_analysis_eligible", "analysis_reason", "excluded_reason", "is_endocrine",
          "endocrine_route", "islet_route", "resolved_by", "endocrine_object_class",
          "endocrine_object_id", "islet_key", "X_centroid", "Y_centroid", "disease_status",
          "primary_call_state", "hormone_measurement_state", "analysis_eligible",
          "analysis_exclusion_reason", "endocrine_subtype_status", "immune_rules_unavailable_count"]


def validate_assignment_census(records: list[dict], selected_keys: list[tuple],
                                expected_keys: set[tuple], historical: bool) -> dict:
    """Verify exact exclusion coverage; historical totals apply only to their frozen input."""
    if len(selected_keys) != len(set(selected_keys)) or set(selected_keys) != expected_keys:
        raise ValueError("mandatory exclusion keys did not match exactly once in assignment sources")
    census = {name: sum(r[name] for r in records) for name in
              ("source_cells", "QC_cells", "typing_analysis_cells", "typing_excluded_cells")}
    census["mandatory_exclusion_keys"] = len(expected_keys)
    if not historical:
        for name in ("mandatory_exclusions_QC_eligible", "mandatory_exclusions_already_QC_ineligible"):
            census[name] = sum(r[name] for r in records)
        if census["mandatory_exclusion_keys"] != (census["mandatory_exclusions_QC_eligible"] +
                                                  census["mandatory_exclusions_already_QC_ineligible"]):
            raise ValueError("mandatory exclusion QC status census does not cover every exact key")
    if historical:
        assert len(selected_keys) == 25
        assert census["source_cells"] == 37_642_209
        assert census["QC_cells"] == 31_397_179
        assert census["typing_analysis_cells"] == 31_385_695
    return census


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--package", type=Path,
                   help="Required for historical_20261005; current_generation does not modify a package")
    p.add_argument("--source", type=Path, default=Path("/data/io60-calibration/pancreas_release_production"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--policy", choices=["historical_20261005", "current_generation"],
                   default="historical_20261005",
                   help="current_generation writes a fresh view and computes its source/key census")
    p.add_argument("--mandatory-exclusions", type=Path, default=DEFAULT_EXCLUSIONS)
    args = p.parse_args()
    historical = args.policy == "historical_20261005"
    if historical and args.package is None:
        p.error("--package is required for --policy historical_20261005")
    package = args.package.resolve() if args.package is not None else None
    out = args.out.resolve()
    if not historical and out.exists() and any(out.iterdir()):
        raise FileExistsError("current_generation requires a fresh output directory")
    overlay_path = (package / "source_tables/protein_analysis_exclusions_20261005.parquet"
                    if historical else args.mandatory_exclusions.resolve())
    overlay = pd.read_parquet(overlay_path) if historical else None
    if historical:
        overlay["donor_id"] = overlay.donor_id.astype(str).map(lambda x: str(int(x)))
    records, selected_keys, expected_keys = [], [], set()
    sources = sorted(args.source.glob("donor_id=*/data.parquet"))
    if not sources:
        raise ValueError("no protein assignment source partitions")
    for source in sources:
        donor_native = source.parent.name.split("=", 1)[1]
        donor = str(int(donor_native))
        parquet = pq.ParquetFile(source)
        excluded = (overlay.loc[overlay.donor_id.eq(donor)] if historical
                    else load_protein_exclusions(overlay_path, donor))
        key = next((c for c in ("cell_id", "object_id")
                    if c in excluded and c in parquet.schema.names), None)
        if key is None:
            raise ValueError(f"{donor}: exclusions and assignment source share no exact cell key")
        changes = excluded.set_index(key)
        expected_keys.update((donor, str(c)) for c in changes.index)
        columns = [c for c in FIELDS if c in parquet.schema.names]
        boolean_columns = {field.name for field in parquet.schema_arrow if pa.types.is_boolean(field.type)}
        boolean_columns.add("native_is_endocrine")
        path = out / f"donor_id={donor}" / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        temp = path.with_suffix(".parquet.tmp")
        writer = None
        n_rows = n_qc = n_identity = n_new = n_excluded = 0
        n_mandatory_qc = n_mandatory_qc_ineligible = 0
        for batch in parquet.iter_batches(columns=columns, batch_size=100000):
            frame = batch.to_pandas()
            frame["donor_id_native"] = donor_native
            frame["donor_id"] = donor
            frame["native_cell_class"] = frame.cell_class.copy()
            frame["native_refined_class"] = frame.refined_class.copy()
            frame["native_is_endocrine"] = frame.is_endocrine.copy()
            source_analysis = (frame.analysis_eligible.fillna(False).astype(bool).copy()
                               if "analysis_eligible" in frame else pd.Series(True, index=frame.index))
            source_reason = (frame.analysis_exclusion_reason.copy() if "analysis_exclusion_reason" in frame
                             else pd.Series("source_analysis_ineligible", index=frame.index))
            hit = frame[key].astype(str).isin(changes.index.astype(str))
            if hit.any():
                selected_keys.extend((donor, str(c)) for c in frame.loc[hit, key])
                reason = frame.loc[hit, key].astype(str).map(
                    pd.Series(changes.analysis_exclusion_reason.to_numpy(), index=changes.index.astype(str)))
                for field in ["cell_class", "refined_class", "phenotype", "refined_phenotype"]:
                    if field in frame:
                        frame.loc[hit, field] = "Excluded"
                frame.loc[hit, "is_endocrine"] = False
                for field in ["endocrine_route", "islet_route"]:
                    if field in frame:
                        frame.loc[hit, field] = ""
                frame.loc[hit, "resolved_by"] = ("excluded_single_hormone" if historical
                                                  else "mandatory_keyed_analysis_exclusion")
                frame.loc[hit, "excluded_reason"] = reason.to_numpy()
            qc = frame.qc_analysis_eligible.fillna(False).astype(bool)
            excluded = frame.cell_class.eq("Excluded").fillna(False)
            frame["analysis_eligible"] = qc & ~excluded & source_analysis
            frame["analysis_exclusion_reason"] = ""
            frame.loc[~qc, "analysis_exclusion_reason"] = frame.loc[~qc, "analysis_reason"].fillna("geometry_QC_ineligible")
            frame.loc[qc & excluded, "analysis_exclusion_reason"] = frame.loc[qc & excluded, "excluded_reason"]
            frame.loc[qc & ~excluded & ~source_analysis, "analysis_exclusion_reason"] = source_reason
            frame["analysis_assignment_state"] = "admitted"
            frame.loc[~qc, "analysis_assignment_state"] = "QC_ineligible"
            frame.loc[qc & excluded, "analysis_assignment_state"] = "typing_excluded"
            frame.loc[qc & ~excluded & ~source_analysis, "analysis_assignment_state"] = "source_analysis_ineligible"
            if not historical:
                # A registry key remains a mandatory exclusion even if the
                # corrected geometry has already removed it from the QC mask.
                frame["mandatory_analysis_exclusion"] = hit
                if hit.any():
                    frame.loc[hit, "analysis_exclusion_reason"] = reason.to_numpy()
                frame.loc[hit, "analysis_assignment_state"] = "mandatory_typing_excluded"
                n_mandatory_qc += int((hit & qc).sum())
                n_mandatory_qc_ineligible += int((hit & ~qc).sum())
            frame["measurement_source_path"] = str(source)
            # Nullable Arrow booleans arrive as object dtype when they include
            # missing entries. Preserve boolean meaning before normalizing text.
            for column in boolean_columns.intersection(frame.columns):
                frame[column] = frame[column].astype("boolean")
            for column in frame.columns:
                if frame[column].dtype == object or isinstance(frame[column].dtype, pd.CategoricalDtype):
                    frame[column] = frame[column].astype("string")
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(temp, table.schema, compression="zstd", use_dictionary=True)
            writer.write_table(table.cast(writer.schema))
            n_rows += len(frame)
            n_qc += int(qc.sum())
            n_identity += int(frame.analysis_eligible.sum())
            n_new += int(hit.sum())
            n_excluded += int((qc & excluded).sum())
        assert writer is not None
        if n_rows != parquet.metadata.num_rows:
            raise ValueError(f"{donor}: source row census changed during assignment materialization")
        writer.close()
        temp.replace(path)
        record = {"donor_id": donor, "donor_id_native": donor_native, "path": str(path),
                        "native_source_path": str(source), "source_cells": n_rows, "QC_cells": n_qc,
                        "typing_analysis_cells": n_identity, "typing_excluded_cells": n_excluded,
                        ("new_single_exclusions" if historical else "mandatory_typing_exclusions"): n_new,
                        "source_key": key, "bytes": path.stat().st_size, "sha256": digest(path)}
        if not historical:
            record.update(mandatory_exclusions_QC_eligible=n_mandatory_qc,
                          mandatory_exclusions_already_QC_ineligible=n_mandatory_qc_ineligible)
        records.append(record)
        print(f"{donor}: {n_rows:,} retained, {n_identity:,} admitted, {n_new} new exclusions", flush=True)
    if historical:
        expected_keys = set(zip(overlay.donor_id, overlay.cell_id.astype(str)))
    census = validate_assignment_census(records, selected_keys, expected_keys, historical)
    source_keys = sorted({r["source_key"] for r in records})
    manifest = {"state": "complete", "analysis_revision": ("2026-10-05 human-approved Single exclusions"
                    if historical else "current generation with pinned mandatory exact-key exclusions"),
                "policy": args.policy, "census": census,
                "overlay": str(overlay_path), "stable_keys": ["donor_id", *source_keys],
                "stable_keys_by_donor": {r["donor_id"]: ["donor_id", r["source_key"]] for r in records},
                "source_measurements_copied": False, "native_source_files_modified": False,
                "alternative_view_warning": "These partitions contain the same cells as native sources; never sum both views",
                "native_evidence_preserved": ["native_cell_class", "native_refined_class", "native_is_endocrine"],
                "current_identity_fields": ["cell_class", "refined_class", "is_endocrine"],
                "eligibility": "analysis_eligible; QC eligibility and typing exclusions are separate",
                "donors": records}
    if not historical:
        manifest["mandatory_exclusion_status"] = ("Exact registry membership is retained in mandatory_analysis_exclusion, "
            "including keys already QC-ineligible in this generation; none are re-admitted. The independent QC-status "
            "census covers every registry key; typing_excluded_cells counts only QC-eligible excluded identities.")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not historical:
        print(json.dumps({"state": "complete", "output": str(out), "census": census,
                          "historical_package_index_modified": False}, indent=2))
        return
    index_path = package / "cell_level_release_index.csv"
    index = pd.read_csv(index_path, dtype=str, keep_default_na=False)
    role = "current_analysis_revision_assignments"
    index = index.loc[~index.role.eq(role)].copy()
    extra = []
    for r in records:
        row = {c: "" for c in index.columns}
        row.update({"lane": "pancreas_protein", "tissue": "pancreas", "modality": "protein", "role": role,
                    "source_path": r["path"], "exists": "True", "retained_rows": str(r["source_cells"]),
                    "bytes": str(r["bytes"]), "donor_id": r["donor_id"], "donor_id_native": r["donor_id_native"],
                    "stable_keys": json.dumps(["donor_id", "cell_id"]),
                    "native_identity_fields": json.dumps(manifest["native_evidence_preserved"]),
                    "qualified_identity_fields": json.dumps(manifest["current_identity_fields"]),
                    "availability_and_QC_fields": json.dumps(["analysis_eligible", "analysis_exclusion_reason", "qc_analysis_eligible"]),
                    "cell_scope": "All native source cells, current analysis assignments; measurements join to native source by stable keys",
                    "native_release_status": "current analysis revision; native evidence retained",
                    "duplicates_other_indexed_cells": "True", "excluded_donor": "6457", "source_copied": "False",
                    "audit_source_manifest": str(out / "manifest.json")})
        extra.append(row)
    pd.concat([index, pd.DataFrame(extra)], ignore_index=True).to_csv(index_path, index=False)
    print("Complete:25 donors,25 exact new exclusions; all rows retained", flush=True)


if __name__ == "__main__":
    main()
