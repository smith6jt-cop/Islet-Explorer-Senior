#!/usr/bin/env python3
"""Build one no-expansion composition row per typed endocrine object.

Primary features are cell-body counts per endocrine cell in each object. The measured outline
coverage of vessels, lymphatics, and nerves is kept beside (not substituted for) those cell-body
counts. Hormone intensities and refined labels are supplementary. Missing denominators or missing
panel evidence stay null; no pseudocount is added.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from shapely.ops import unary_union
from shapely.strtree import STRtree
from multiplex5000.tissues.macrophage import labels as macrophage_labels

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_islet_geometry as GEOM
import context as CX
import islet_membership as MEMBERSHIP
import subpanel_io as S
from protein_analysis_admission import ENDOCRINE_UNCERTAINTY, protein_analysis_mask

OUTLINE_STRUCTURES = {
    "vascular": ("Blood vessel", "blood_vessel"),
    "lymphatic": ("Lymphatic", "lymphatic"),
    "neural": ("Nerve", "nerve"),
}
GROUP_MARKERS = {
    "endocrine": ("INS", "GCG", "SST"),
    "vascular": ("CD31", "Caveolin", "CD34"),
    "lymphatic": ("Podoplanin",),
    "neural": ("B3TUBB", "CD56"),
    "stromal": ("SMA", "Vimentin", "Collagen_IV"),
    "immune": ("CD3e", "CD68", "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11c"),
}
IMMUNE_LABELS = (
    "Imm-T-cell", "Imm-B-cell", "Imm-macrophage", "Imm-dendritic-cell",
    "Imm-NK-cell", "Imm-neutrophil", "Imm-plasma-cell",
)
IMMUNE_SUBTYPES = (
    "T-regulatory", "T-exhausted", "T-progenitor-memory", "T-CD8-cytotoxic", "T-CD8",
    "T-CD4-helper",
)
ENDOCRINE_TYPES = ("Beta", "Alpha", "Delta")
STROMAL_LABELS = {"Stromal", "Fibroblast", "Stellate", "Pericyte", "Smooth muscle"}


def read_cell_exclusions(path: Path | None) -> pd.DataFrame:
    """Read explicit scientific exclusions without changing original object membership."""
    if path is None:
        return pd.DataFrame(columns=["donor_id", "object_id", "endocrine_object_id", "exclusion_reason"])
    frame = pd.read_parquet(path) if Path(path).suffix == ".parquet" else pd.read_csv(path, dtype=str)
    required = {"donor_id", "object_id", "endocrine_object_id"}
    if not required.issubset(frame):
        raise ValueError(f"cell exclusion overlay requires {sorted(required)}")
    frame = frame.copy()
    for column in required:
        frame[column] = frame[column].astype(str)
    if frame.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("explicit cell exclusion overlay repeats a donor/cell key")
    if "exclusion_reason" not in frame:
        frame["exclusion_reason"] = frame.get("analysis_exclusion_reason", "unresolved_or_ambiguous_EndocrineSingle_without_supported_hormone_identity")
    return frame


def native_single_exclusions(typed: pd.DataFrame, donor: str) -> pd.DataFrame:
    """Carry the declared native hormone-Single exclusion into future object builds."""
    fields = {"object_id", "endocrine_object_id", "endocrine_object_class", "resolved_by"}
    if not fields.issubset(typed):
        return read_cell_exclusions(None)
    selected = typed[typed["endocrine_object_class"].astype(str).eq("EndocrineSingle")
                     & typed["resolved_by"].astype(str).eq("excluded_single_hormone")]
    if selected.empty:
        return read_cell_exclusions(None)
    if "excluded_reason" not in selected or selected["excluded_reason"].isna().any():
        raise ValueError("native excluded hormone Singles require their explicit exclusion reason")
    allowed = {"endocrine_single_no_usable_hormone", "endocrine_single_tied_highest_hormones"}
    if not selected["excluded_reason"].astype(str).isin(allowed).all():
        raise ValueError("native hormone Single has an undeclared exclusion reason")
    result = selected[["object_id", "endocrine_object_id", "excluded_reason"]].copy()
    result["donor_id"] = donor
    return result.rename(columns={"excluded_reason": "exclusion_reason"})


def merge_cell_exclusions(explicit: pd.DataFrame, native: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate matching explicit/native records and reject conflicting decisions."""
    fields = ["donor_id", "object_id", "endocrine_object_id", "exclusion_reason"]
    result = pd.concat([explicit[fields], native[fields]], ignore_index=True)
    if result.empty:
        return result
    result["donor_id"] = result["donor_id"].map(donor_id_key)
    result["object_id"] = result["object_id"].astype(str)
    for _, block in result.groupby(["donor_id", "object_id"]):
        if block["endocrine_object_id"].astype(str).nunique() != 1 or block["exclusion_reason"].astype(str).nunique() != 1:
            raise ValueError("explicit/native Single exclusions disagree on object or reason")
    return result.drop_duplicates(["donor_id", "object_id"])


def apply_composition_cell_exclusions(objects: pd.DataFrame, exclusions: pd.DataFrame) -> pd.DataFrame:
    """Preserve single-object source counts/geometry; exclude its sole cell from analysis.

    This helper intentionally accepts only explicitly enumerated, one-cell endocrine
    Singles. An exclusion of a member of a larger object requires cell-level rebuilding
    and cannot be approximated by subtracting a summary count.
    """
    result = objects.copy()
    if result.duplicated(["donor_id", "endocrine_object_id"]).any():
        raise ValueError("composition requires unique donor/object keys")
    raw_count_columns = [c for c in result if c.startswith(("primary_count__", "refined_count__", "count_"))]
    source_measurement_columns = [c for c in result if c.startswith(("process_overlap_", "supplementary_"))]
    original_fields = ["endocrine_denominator", "identity_unresolved_count", *raw_count_columns, *source_measurement_columns]
    for column in original_fields:
        name = "original_" + column
        if name not in result:
            result[name] = result[column]
    result["n_original_cells_in_object"] = result.get("n_original_cells_in_object", result["n_cells_in_object"])
    result["n_analysis_excluded_cells_in_object"] = result.get("n_analysis_excluded_cells_in_object", 0)
    result["n_analysis_admitted_cells_in_object"] = result.get("n_analysis_admitted_cells_in_object", result["n_cells_in_object"])
    result["composition_analysis_eligible"] = result.get("composition_analysis_eligible", True)
    result["composition_exclusion_reason"] = result.get("composition_exclusion_reason", "")
    result["excluded_cell_ids"] = result.get("excluded_cell_ids", "")
    result["endocrine_ratio_availability_reason"] = result.get("endocrine_ratio_availability_reason", "")
    if exclusions.empty:
        return result
    keys = result.set_index(["donor_id", "endocrine_object_id"], verify_integrity=True).index
    donor_keys = {donor_id_key(d): str(d) for d in result.donor_id.astype(str).unique()}
    for exclusion in exclusions.itertuples(index=False):
        donor = donor_keys.get(donor_id_key(exclusion.donor_id))
        if donor is None:
            continue  # A per-donor producer receives the complete cohort overlay.
        key = (donor, str(exclusion.endocrine_object_id))
        if key not in keys:
            raise ValueError(f"explicit excluded cell has no preserved composition object: {key}")
        index = keys.get_loc(key)
        row = result.iloc[index]
        if row.original_object_class != "EndocrineSingle" or int(row.n_cells_in_object) != 1:
            raise ValueError(f"summary exclusion is licensed only for a sole-cell Single: {key}")
        if int(row.endocrine_denominator) not in (0, 1):
            raise ValueError(f"sole-cell Single has an invalid endocrine denominator: {key}")
        result.loc[result.index[index], "original_endocrine_denominator"] = 1
        result.loc[result.index[index], "original_count_endocrine__unresolved_subtype"] = 1
        result.loc[result.index[index], "n_analysis_excluded_cells_in_object"] = 1
        result.loc[result.index[index], "n_analysis_admitted_cells_in_object"] = 0
        result.loc[result.index[index], "composition_analysis_eligible"] = False
        result.loc[result.index[index], "composition_exclusion_reason"] = str(exclusion.exclusion_reason)
        result.loc[result.index[index], "excluded_cell_ids"] = str(exclusion.object_id)
        result.loc[result.index[index], "endocrine_denominator"] = 0
        result.loc[result.index[index], "identity_unresolved_count"] = 0
        result.loc[result.index[index], "refined_typing_present"] = False
        result.loc[result.index[index], "endocrine_ratio_availability_reason"] = "zero_admitted_endocrine_denominator_after_explicit_single_exclusion"
        for column in raw_count_columns:
            # A missing unmeasured subtype remains unavailable, not a measured zero.
            if pd.notna(row[column]):
                result.loc[result.index[index], column] = 0
        for column in result.columns:
            if column.startswith("feature__") and not column.endswith("outline_coverage_fraction"):
                result.loc[result.index[index], column] = np.nan
                availability = "available__" + column.removeprefix("feature__")
                if availability in result:
                    result.loc[result.index[index], availability] = False
        for column in source_measurement_columns:
            if column.startswith("supplementary_available__"):
                result.loc[result.index[index], column] = False
            elif column.startswith(("supplementary_n_measured__", "process_overlap_measured_cells__")):
                result.loc[result.index[index], column] = 0
            else:
                result.loc[result.index[index], column] = np.nan
    return result


def donor_id_key(value: object) -> str:
    """Normalize numeric donor IDs so padded typing partitions match inventory IDs."""
    text = str(value).strip()
    if text.isdecimal():
        return str(int(text))
    return text


def composition_cell_analysis_mask(cells: pd.DataFrame, donor: str,
                                   exclusions: pd.DataFrame) -> pd.Series:
    """Use eligibility and exact exclusions; ordinary endocrine uncertainty stays admitted."""
    return protein_analysis_mask(cells, donor, exclusions)


def apply_and_rebuild_composition_cell_exclusions(objects: pd.DataFrame, cells: pd.DataFrame,
                                                donor: str, exclusions: pd.DataFrame) -> pd.DataFrame:
    """Validate exact cell/object keys and rebuild larger objects from cell membership.

    An annotation named EndocrineSingle can contain multiple segmented cells.
    Only actual one-cell Singles use the summary-only helper. Every object is
    then recomputed with the exact cell mask, retaining admitted co-members.
    """
    if cells.object_id.isna().any() or cells.object_id.astype(str).duplicated().any():
        raise ValueError("composition source requires unique nonmissing cell keys")
    if objects.duplicated(["donor_id", "endocrine_object_id"]).any():
        raise ValueError("composition requires unique donor/object keys")
    if not objects.donor_id.map(donor_id_key).eq(donor_id_key(donor)).all():
        raise ValueError("cell-level composition routing requires one exact donor scope")
    sole = exclusions.iloc[:0].copy()
    if len(exclusions):
        required = {"donor_id", "object_id", "endocrine_object_id", "exclusion_reason"}
        if not required.issubset(exclusions):
            raise ValueError(f"cell-level exclusions require {sorted(required)}")
        selected = exclusions.loc[exclusions.donor_id.map(donor_id_key).eq(donor_id_key(donor))].copy()
        if selected[["object_id", "endocrine_object_id"]].isna().any().any() or selected.object_id.astype(str).duplicated().any():
            raise ValueError("cell-level exclusions repeat or omit a donor/cell/object key")
        membership = cells.set_index(cells.object_id.astype(str)).endocrine_object_id.astype(str)
        object_rows = objects.set_index(objects.endocrine_object_id.astype(str))
        for exclusion in selected.itertuples(index=False):
            cell_key, object_key = str(exclusion.object_id), str(exclusion.endocrine_object_id)
            if cell_key not in membership.index:
                raise ValueError(f"explicit excluded cell has no exact source membership: {(donor, cell_key)}")
            if membership.loc[cell_key] != object_key:
                raise ValueError(f"explicit excluded cell/object membership disagrees: {(donor, cell_key, object_key)}")
            if object_key not in object_rows.index:
                raise ValueError(f"explicit excluded cell has no preserved composition object: {(donor, object_key)}")
        sole_objects = set(object_rows.index[object_rows.original_object_class.eq("EndocrineSingle")
                                           & object_rows.n_cells_in_object.eq(1)])
        sole = selected.loc[selected.endocrine_object_id.astype(str).isin(sole_objects)]
    initialized = apply_composition_cell_exclusions(objects, sole)
    return rebuild_analytic_object_composition(initialized, cells, donor, exclusions)


def rebuild_analytic_object_composition(objects: pd.DataFrame, cells: pd.DataFrame,
                                       donor: str, exclusions: pd.DataFrame) -> pd.DataFrame:
    """Recompute every cell contribution; preserve raw object counts and geometry.

    Candidate/reviewed identities never override the primary analysis mask. A
    retained object with no admitted cells is unavailable for composition or ND
    analysis, even when its independent outline geometry remains measured.
    """
    result = objects.copy()
    if result.duplicated(['donor_id', 'endocrine_object_id']).any():
        raise ValueError('composition has duplicate donor/object keys')
    raw_fields = [c for c in result if c.startswith(('primary_count__', 'refined_count__', 'count_', 'process_overlap_', 'supplementary_'))]
    for c in ['endocrine_denominator', 'identity_unresolved_count', *raw_fields]:
        if 'original_' + c not in result:
            result['original_' + c] = result[c]
    admitted = composition_cell_analysis_mask(cells, donor, exclusions)
    parts = {str(k): g for k, g in cells.groupby('endocrine_object_id', sort=False)}
    rows = []
    for original in result.to_dict('records'):
        row = original.copy()
        object_id = str(row['endocrine_object_id'])
        if object_id not in parts:
            raise ValueError('retained object has no original cell membership: ' + object_id)
        raw = parts[object_id]
        block = raw.loc[admitted.loc[raw.index]]
        omitted = raw.loc[~admitted.loc[raw.index]]
        if int(row['n_cells_in_object']) != len(raw):
            raise ValueError('original cell count disagrees with exact source membership: ' + object_id)
        row['n_original_cells_in_object'] = len(raw)
        row['n_analysis_admitted_cells_in_object'] = len(block)
        row['n_analysis_excluded_cells_in_object'] = len(omitted)
        row['composition_analysis_eligible'] = bool(len(block))
        row['composition_exclusion_reason'] = '' if len(block) else 'no_analysis_admitted_cells_after_typing_exclusions'
        row['excluded_cell_ids'] = '|'.join(omitted.object_id.astype(str))
        row['refined_typing_present'] = bool(block.refined_class.notna().any())
        counts = block.cell_class.astype(str).value_counts()
        endocrine = block.is_endocrine.fillna(False).astype(bool)
        denominator = int(endocrine.sum())
        row['endocrine_denominator'] = denominator
        row['endocrine_ratio_availability_reason'] = '' if denominator else 'zero_admitted_endocrine_denominator_after_typing_exclusions'
        row['identity_unresolved_count'] = int(block.cell_class.astype(str).isin(
            {'Ambiguous', 'Unresolved', 'Unresolved (no lineage evidence)'} | ENDOCRINE_UNCERTAINTY).sum())
        for label in [*ENDOCRINE_TYPES, 'Immune', 'Blood vessel', 'Lymphatic', 'Nerve', *STROMAL_LABELS, 'Exocrine']:
            suffix = label.lower().replace(' ', '_')
            row['primary_count__' + suffix] = int(counts.get(label, 0))
            row['refined_count__' + suffix] = int(block.refined_class.astype(str).eq(label).sum())
        unresolved = int((endocrine & ~block.cell_class.astype(str).isin(ENDOCRINE_TYPES)).sum())
        row['endocrine_subtype_composition_complete'] = bool(denominator and not unresolved)
        row['endocrine_subtype_availability_reason'] = (
            'contains_admitted_endocrine_cells_without_supported_subtype' if unresolved
            else '' if denominator else 'zero_admitted_endocrine_denominator_after_typing_exclusions')
        for cell_type in ENDOCRINE_TYPES:
            name = cell_type.lower()
            hormone = {'Beta': 'INS', 'Alpha': 'GCG', 'Delta': 'SST'}[cell_type]
            n = int((endocrine & block.cell_class.astype(str).eq(cell_type)).sum())
            available = bool(denominator > 0 and not unresolved and hormone in block
                             and np.isfinite(pd.to_numeric(block[hormone], errors='coerce')).any())
            row['count_endocrine__' + name] = n
            row['feature__endocrine__' + name + '_fraction'] = n / denominator if available else np.nan
            row['available__endocrine__' + name + '_fraction'] = available
        row['count_endocrine__unresolved_subtype'] = unresolved
        row['feature__endocrine__unresolved_fraction'] = unresolved / denominator if denominator else np.nan
        row['available__endocrine__unresolved_fraction'] = bool(denominator)
        for group, (_, suffix) in OUTLINE_STRUCTURES.items():
            label = {'blood_vessel': 'Blood vessel', 'lymphatic': 'Lymphatic', 'nerve': 'Nerve'}[suffix]
            n = int(counts.get(label, 0))
            row['count_' + suffix + '_cell_bodies'] = n
            feature = 'feature__' + group + '__cell_bodies_per_endocrine'
            availability = 'available__' + group + '__cell_bodies_per_endocrine'
            available = bool(original[availability]) and denominator > 0
            row[feature] = n / denominator if available else np.nan
            row[availability] = available
        stromal = int(sum(counts.get(label, 0) for label in STROMAL_LABELS))
        row['count_stromal_cell_bodies'] = stromal
        available = bool(row['stromal_identity_available_in_donor']) and denominator > 0
        row['feature__stromal__cell_bodies_per_endocrine'] = stromal / denominator if available else np.nan
        row['available__stromal__cell_bodies_per_endocrine'] = available
        immune_n = int(counts.get('Immune', 0))
        row['count_immune_cell_bodies'] = immune_n
        available = bool(original['available__immune__all_immune_per_endocrine']) and denominator > 0
        row['feature__immune__all_immune_per_endocrine'] = immune_n / denominator if available else np.nan
        row['available__immune__all_immune_per_endocrine'] = available
        for labels, field, prefix, feature_prefix in [(IMMUNE_LABELS, 'immune_class', 'count_immune__', 'feature__immune__'),
                                                      (IMMUNE_SUBTYPES, 'immune_subtype', 'count_immune_subtype__', 'feature__immune__subtype_')]:
            for label in labels:
                n = int(block[field].eq(label).sum())
                measured = pd.notna(original[prefix + label])
                row[prefix + label] = n if measured else np.nan
                feature = feature_prefix + label + '_per_endocrine'
                availability = feature.replace('feature__', 'available__', 1)
                row[feature] = n / denominator if measured and denominator else np.nan
                row[availability] = bool(measured and denominator)
        for label in macrophage_labels():
            key = 'count_macrophage_profile__' + label
            measured = pd.notna(original.get(key, np.nan))
            n = int(block['macrophage_profile'].eq(label).sum()) if 'macrophage_profile' in block else 0
            row[key] = n if measured else np.nan
            feature = 'feature__immune__macrophage_profile_' + label + '_per_endocrine'
            row[feature] = n / denominator if measured and denominator else np.nan
            row[feature.replace('feature__', 'available__', 1)] = bool(measured and denominator)
        for structure, col in [('blood_vessel', 'passed_through_Blood_vessel'), ('lymphatic', 'passed_through_Lymphatic'), ('nerve', 'passed_through_Nerve')]:
            values = pd.to_numeric(block[col], errors='coerce') if col in block else pd.Series(np.nan, index=block.index)
            finite = values.notna()
            row['process_overlap_fraction__' + structure] = float(values[finite].mean()) if finite.any() else np.nan
            row['process_overlap_cell_fraction__' + structure] = float(values[finite].gt(0).mean()) if finite.any() else np.nan
            row['process_overlap_measured_cells__' + structure] = int(finite.sum())
        for hormone in ['INS', 'GCG', 'SST']:
            values = pd.to_numeric(block[hormone], errors='coerce') if hormone in block else pd.Series(np.nan, index=block.index)
            finite = values.notna()
            row['supplementary_mean_aligned_intensity__' + hormone] = float(values[finite].mean()) if finite.any() else np.nan
            row['supplementary_n_measured__' + hormone] = int(finite.sum())
            row['supplementary_available__' + hormone] = bool(original['original_supplementary_available__' + hormone] and len(block))
        rows.append(row)
    return pd.DataFrame(rows)


def disease_status_map(sections: list[dict]) -> dict[str, str]:
    """Build an unambiguous disease-status lookup keyed by normalized donor ID."""
    statuses: dict[str, str] = {}
    for section in sections:
        donor = donor_id_key(section.get("donor_id", ""))
        status = str(section.get("disease_status", "")).strip()
        if donor and status:
            if donor in statuses and statuses[donor] != status:
                raise ValueError(f"{donor}: inconsistent disease status in RNA inventory")
            statuses[donor] = status
    return statuses


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def marker_availability(path: Path, wanted: tuple[str, ...]) -> dict[str, bool]:
    parquet = pq.ParquetFile(path)
    md = parquet.metadata
    names = {md.schema.column(i).name: i for i in range(md.num_columns)}
    out = {}
    for marker in wanted:
        index = names.get(marker)
        if index is None:
            out[marker] = False
            continue
        stats = [md.row_group(g).column(index).statistics for g in range(md.num_row_groups)]
        if any(s is None for s in stats):
            out[marker] = True  # conservative; a later cell-level measurement check refines this
        else:
            out[marker] = sum(s.null_count for s in stats) < md.num_rows
    return out


def structure_coverage(donor: str, stem: str, pixels: float,
                       outlines: dict[str, object]) -> dict[str, dict[str, float | None]]:
    """Measure unioned structure area inside each true core outline, in physical units."""
    output = {name: {key: {"area_um2": None, "fraction_of_object_area": None}
                     for key in outlines}
              for name in OUTLINE_STRUCTURES}
    scale_area = float(pixels) ** 2
    for group, (structure_name, _) in OUTLINE_STRUCTURES.items():
        try:
            objects = CX.load_objects(stem, structure_name)
        except FileNotFoundError:
            continue
        if not objects:
            continue
        tree = STRtree(objects)
        for key, outline in outlines.items():
            area = float(outline.area)
            candidates = tree.query(outline)
            if len(candidates) == 0:
                covered = 0.0
            else:
                merged = unary_union([objects[int(i)] for i in candidates])
                covered = float(merged.intersection(outline).area)
            output[group][key] = {
                "area_um2": covered * scale_area,
                "fraction_of_object_area": covered / area if area > 0 else None,
            }
    return output


def read_immune(typing_root: Path, donor: str) -> tuple[pd.DataFrame, dict]:
    path = typing_root / "immune" / f"donor_id={donor}" / "data.parquet"
    if not path.is_file():
        return pd.DataFrame(columns=["object_id", "immune_class", "immune_subtype",
                                     "macrophage_profile", "macrophage_profile_status"]), {}
    keep = ["object_id", "immune_class", "immune_subtype"]
    schema = pq.ParquetFile(path).schema_arrow.names
    keep += [name for name in ("macrophage_profile", "macrophage_profile_status")
             if name in schema]
    frame = pq.read_table(path, columns=keep).to_pandas()
    coverage_path = typing_root / "immune" / "coverage.json"
    coverage = json.loads(coverage_path.read_text()).get("per_donor", {}).get(donor, {}) \
        if coverage_path.is_file() else {}
    return frame, coverage


def build_donor(typing_root: Path, aligned_root: Path, donor: str,
                disease_status: str, cell_exclusions: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    typed_path = typing_root / "cells" / f"donor_id={donor}" / "data.parquet"
    calls_path = typing_root / "subpanel_rebuilt" / f"donor_id={donor}" / "data.parquet"
    if not typed_path.is_file():
        raise FileNotFoundError(f"{donor}: missing final cell typing {typed_path}")
    typed = pq.read_table(typed_path).to_pandas()
    if typed["object_id"].duplicated().any():
        raise ValueError(f"{donor}: final cell table repeats object_id")
    required = ["object_id", "cell_class", "refined_class", "is_endocrine",
                "endocrine_object_class", "endocrine_object_id", "islet_key", "cell_region"]
    missing = sorted(set(required) - set(typed.columns))
    if missing:
        raise ValueError(f"{donor}: required typing columns missing: {missing}")
    typed["donor_id"] = donor
    if "signed_dist_islet_um" not in typed:
        typed["signed_dist_islet_um"] = np.nan

    # Only geometry-backed true-core object membership is used. Tissue/background rows have no
    # endocrine_object_id and therefore cannot leak into any islet composition.
    inside = typed[typed["endocrine_object_id"].notna()].copy()
    if inside.empty:
        return pd.DataFrame(), {"donor_id": donor, "state": "no_endocrine_objects",
                                "typed_rows": len(typed), "object_rows": 0}

    immune, immune_coverage = read_immune(typing_root, donor)
    if len(immune):
        inside = inside.merge(immune, on="object_id", how="left", validate="one_to_one")
    else:
        inside["immune_class"] = pd.NA
        inside["immune_subtype"] = pd.NA
    if "macrophage_profile" not in inside:
        inside["macrophage_profile"] = pd.NA
    if "macrophage_profile_status" not in inside:
        inside["macrophage_profile_status"] = pd.NA

    calls = None
    if calls_path.is_file():
        calls = pq.read_table(calls_path).to_pandas()
        passed = [c for c in calls.columns if c.startswith("passed_through_") and c not in inside]
        if passed:
            inside = inside.merge(calls[["object_id", *passed]], on="object_id", how="left",
                                  validate="one_to_one")

    aligned_path = aligned_root / f"donor_id={donor}/data.parquet"
    hormones = ["INS", "GCG", "SST"]
    hormones_measured = {m: False for m in hormones}
    if aligned_path.is_file():
        available = marker_availability(aligned_path, tuple(hormones))
        aligned = pq.read_table(aligned_path, columns=["object_id"] +
                                [h for h in hormones if available[h]]).to_pandas()
        aligned = S.mask_fragment_hormones(aligned, donor)
        hormones_measured = {m: bool(m in aligned.columns and aligned[m].notna().any())
                             for m in hormones}
        inside = inside.merge(aligned, on="object_id", how="left", validate="one_to_one")
    else:
        for hormone in hormones:
            inside[hormone] = np.nan

    pixels, stems = GEOM._manifest()
    stem = stems[donor]
    outlines = MEMBERSHIP.outlines(stem)
    object_area = {f"Islet_{k}": float(shape.area) * float(pixels[donor]) ** 2
                   for k, shape in outlines.items()}
    geometry_coverage = structure_coverage(donor, stem, pixels[donor], outlines)

    if "cell_class" not in inside:
        raise AssertionError("cell class disappeared before composition build")
    if inside[["endocrine_object_id", "islet_key", "endocrine_object_class"]].isna().any().any():
        raise ValueError(f"{donor}: endocrine object is missing an authoritative ID/class/key")

    group_available = {
        group: any(marker_availability(aligned_path, markers).values())
        if aligned_path.is_file() else False
        for group, markers in GROUP_MARKERS.items()
    }
    immune_rule_unavailable = set(immune_coverage.get("unavailable_rules", {}))
    immune_class_counts = immune_coverage.get("immune_class", {})
    immune_parent_available = {
        label: bool(immune_coverage.get("state") == "built")
        and label not in immune_rule_unavailable for label in IMMUNE_LABELS
    }
    subtype_parent = {"T-regulatory": "Imm-T-cell", "T-exhausted": "Imm-T-cell",
                      "T-progenitor-memory": "Imm-T-cell", "T-CD8-cytotoxic": "Imm-T-cell",
                      "T-CD8": "Imm-T-cell", "T-CD4-helper": "Imm-T-cell"}
    immune_subtype_available = {
        label: bool(immune_coverage.get("state") == "built")
        and label not in immune_rule_unavailable
        and int(immune_class_counts.get(subtype_parent[label], 0)) > 0
        for label in IMMUNE_SUBTYPES
    }
    stromal_identity_available = bool(
        typed["cell_class"].astype(str).isin(STROMAL_LABELS).any())
    # Marker presence is not a stromal identity call. Until the native typing
    # chain emits a stromal class, the cell-body profile remains unavailable.
    group_available["stromal"] = stromal_identity_available
    immune_profile_state = "available" if any(immune_parent_available.values()) else \
        "unavailable_rule_markers"

    rows = []
    for object_id, block in inside.groupby("endocrine_object_id", sort=False, dropna=False):
        key_values = block["islet_key"].astype(str).unique()
        class_values = block["endocrine_object_class"].astype(str).unique()
        if len(key_values) != 1 or len(class_values) != 1:
            raise ValueError(f"{donor}/{object_id}: object ID maps to multiple keys/classes")
        key = key_values[0]
        object_class = class_values[0]
        counts = block["cell_class"].astype(str).value_counts(dropna=False)
        endocrine = block["is_endocrine"].fillna(False).astype(bool)
        denominator = int(endocrine.sum())
        row = {
            "donor_id": donor, "donor_status": disease_status,
            "endocrine_object_id": str(object_id), "islet_key": key,
            "original_object_class": object_class,
            "n_cells_in_object": int(len(block)), "endocrine_denominator": denominator,
            "endocrine_denominator_available": True,
            "outline_area_um2": object_area.get(key),
            "outline_area_available": key in object_area,
            "refined_typing_present": bool(block["refined_class"].notna().any()),
            "immune_typing_state": immune_coverage.get("state", "unavailable"),
            "immune_profile_state": immune_profile_state,
            "stromal_identity_available_in_donor": stromal_identity_available,
            "signed_distance_to_nearest_islet_available": "signed_dist_islet_um" in block.columns,
            "identity_unresolved_count": int(block["cell_class"].astype(str).isin(
                {"Ambiguous", "Unresolved", "Unresolved (no lineage evidence)", "Excluded"}
                | ENDOCRINE_UNCERTAINTY).sum()),
        }
        for label in [*ENDOCRINE_TYPES, "Immune", "Blood vessel", "Lymphatic", "Nerve",
                      *STROMAL_LABELS, "Exocrine"]:
            row[f"primary_count__{label.lower().replace(' ', '_')}"] = int(counts.get(label, 0))
            row[f"refined_count__{label.lower().replace(' ', '_')}"] = int(
                block["refined_class"].astype(str).eq(label).sum())
        unresolved_endo = int((endocrine & ~block["cell_class"].astype(str).isin(ENDOCRINE_TYPES)).sum())
        row["endocrine_subtype_composition_complete"] = bool(denominator and not unresolved_endo)
        row["endocrine_subtype_availability_reason"] = (
            "contains_admitted_endocrine_cells_without_supported_subtype" if unresolved_endo
            else "" if denominator else "zero_endocrine_denominator")
        for cell_type in ENDOCRINE_TYPES:
            n = int((endocrine & block["cell_class"].astype(str).eq(cell_type)).sum())
            hormone = {"Beta": "INS", "Alpha": "GCG", "Delta": "SST"}[cell_type]
            hormone_available_object = (hormone in block
                                        and np.isfinite(pd.to_numeric(block[hormone], errors='coerce')).any())
            row[f"count_endocrine__{cell_type.lower()}"] = n
            row[f"feature__endocrine__{cell_type.lower()}_fraction"] = (
                n / denominator if denominator > 0 and not unresolved_endo and hormone_available_object else np.nan)
            row[f"available__endocrine__{cell_type.lower()}_fraction"] = (
                bool(hormone_available_object)
                and denominator > 0 and not unresolved_endo)
        row["count_endocrine__unresolved_subtype"] = unresolved_endo
        row["feature__endocrine__unresolved_fraction"] = unresolved_endo / denominator \
            if denominator > 0 else np.nan
        row["available__endocrine__unresolved_fraction"] = denominator > 0

        for group, (_, suffix) in OUTLINE_STRUCTURES.items():
            phenotype = {"blood_vessel": "Blood vessel", "lymphatic": "Lymphatic", "nerve": "Nerve"}[suffix]
            n = int(counts.get(phenotype, 0))
            row[f"count_{suffix}_cell_bodies"] = n
            row[f"feature__{group}__cell_bodies_per_endocrine"] = (
                n / denominator if denominator > 0 and group_available[group] else np.nan)
            row[f"available__{group}__cell_bodies_per_endocrine"] = (
                bool(group_available[group]) and denominator > 0)
            outline_key = key.removeprefix("Islet_")
            geom = geometry_coverage[group].get(outline_key, {})
            row[f"structure_area_um2__{suffix}"] = geom.get("area_um2")
            row[f"feature__{group}__outline_coverage_fraction"] = geom.get(
                "fraction_of_object_area")
            row[f"available__{group}__outline_coverage_fraction"] = (
                geom.get("fraction_of_object_area") is not None)

        stroma_n = int(sum(counts.get(label, 0) for label in STROMAL_LABELS))
        row["count_stromal_cell_bodies"] = stroma_n
        row["feature__stromal__cell_bodies_per_endocrine"] = (
            stroma_n / denominator if denominator > 0 and stromal_identity_available else np.nan)
        row["available__stromal__cell_bodies_per_endocrine"] = (
            stromal_identity_available and denominator > 0)
        row["stromal_feature_availability_reason"] = (
            "explicit_stromal_identity_measured" if stroma_n > 0 else
            ("no_stromal_identity_cells_called_in_donor" if stromal_identity_available else
             "no_explicit_stromal_identity_call_in_primary_typing"))

        immune_n = int(counts.get("Immune", 0))
        row["count_immune_cell_bodies"] = immune_n
        row["feature__immune__all_immune_per_endocrine"] = (
            immune_n / denominator if denominator > 0 and group_available["immune"] else np.nan)
        row["available__immune__all_immune_per_endocrine"] = (
            bool(group_available["immune"]) and denominator > 0)
        for label in IMMUNE_LABELS:
            n = int(block["immune_class"].eq(label).sum())
            available = immune_parent_available[label]
            row[f"count_immune__{label}"] = n if available else np.nan
            row[f"feature__immune__{label}_per_endocrine"] = (
                n / denominator if available and denominator > 0 else np.nan)
            row[f"available__immune__{label}_per_endocrine"] = available and denominator > 0
        for label in IMMUNE_SUBTYPES:
            n = int(block["immune_subtype"].eq(label).sum())
            available = immune_subtype_available[label]
            row[f"count_immune_subtype__{label}"] = n if available else np.nan
            row[f"feature__immune__subtype_{label}_per_endocrine"] = (
                n / denominator if available and denominator > 0 else np.nan)
            row[f"available__immune__subtype_{label}_per_endocrine"] = (
                available and denominator > 0)
        macrophages = block["immune_class"].eq("Imm-macrophage")
        profile_available = (bool(immune_coverage.get("state") == "built")
                             and int(immune_class_counts.get("Imm-macrophage", 0)) > 0
                             and immune_coverage.get("macrophage_profile") is not None
                             and not block.loc[macrophages, "macrophage_profile"].isna().any())
        for label in macrophage_labels():
            n = int(block["macrophage_profile"].eq(label).sum())
            key = f"count_macrophage_profile__{label}"
            row[key] = n if profile_available else np.nan
            feature = f"feature__immune__macrophage_profile_{label}_per_endocrine"
            row[feature] = n / denominator if profile_available and denominator > 0 else np.nan
            row[feature.replace("feature__", "available__", 1)] = (
                profile_available and denominator > 0)

        # The intersection values are an independent morphology measurement of passing processes,
        # reported separately from structure-body counts and outline coverage.
        for structure, col in [("blood_vessel", "passed_through_Blood_vessel"),
                               ("lymphatic", "passed_through_Lymphatic"),
                               ("nerve", "passed_through_Nerve")]:
            if col in block:
                values = pd.to_numeric(block[col], errors="coerce")
                finite = values.notna()
                row[f"process_overlap_fraction__{structure}"] = (
                    float(values[finite].mean()) if finite.any() else np.nan)
                row[f"process_overlap_cell_fraction__{structure}"] = (
                    float((values[finite] > 0).mean()) if finite.any() else np.nan)
                row[f"process_overlap_measured_cells__{structure}"] = int(finite.sum())
            else:
                row[f"process_overlap_fraction__{structure}"] = np.nan
                row[f"process_overlap_cell_fraction__{structure}"] = np.nan
                row[f"process_overlap_measured_cells__{structure}"] = 0

        for hormone in hormones:
            values = pd.to_numeric(block[hormone], errors="coerce") if hormone in block else pd.Series(
                np.nan, index=block.index)
            finite = values.notna()
            row[f"supplementary_mean_aligned_intensity__{hormone}"] = (
                float(values[finite].mean()) if finite.any() else np.nan)
            row[f"supplementary_n_measured__{hormone}"] = int(finite.sum())
            row[f"supplementary_available__{hormone}"] = bool(hormones_measured.get(hormone))
        rows.append(row)

    result = pd.DataFrame(rows)
    exclusion_decisions = merge_cell_exclusions(
        cell_exclusions if cell_exclusions is not None else read_cell_exclusions(None),
        native_single_exclusions(inside, donor))
    result = apply_and_rebuild_composition_cell_exclusions(result, inside, donor, exclusion_decisions)
    if result["endocrine_object_id"].duplicated().any():
        raise ValueError(f"{donor}: composition repeats endocrine_object_id")
    counts_manifest = {
        "donor_id": donor, "state": "built", "typed_rows": len(typed),
        "object_rows": len(result), "core_membership_rows": len(inside),
        "native_or_explicit_single_exclusion_count": int((~result["composition_analysis_eligible"]).sum()),
        "unmatched_outline_keys": sorted(set(result.islet_key) - set(object_area)),
        "class_counts": result["original_object_class"].value_counts().to_dict(),
        "available_groups": group_available,
        "immune_rule_unavailable": sorted(immune_rule_unavailable),
        "hormone_channels_measured": hormones_measured,
        "explicitly_excluded_single_objects": int((~result.composition_analysis_eligible).sum()),
    }
    return result, counts_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typing-root", required=True, type=Path)
    parser.add_argument("--aligned-root", type=Path, default=None)
    parser.add_argument("--inventory", type=Path,
                        default=Path("/data/xenium/io60_cohort_20261003/input_inventory.json"))
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--analysis-status", default="candidate_preview",
                        help="candidate_preview until the production calibration ruling is recorded")
    parser.add_argument("--cell-exclusions", type=Path, default=None,
                        help="explicit donor/cell/object exclusion overlay; preserves source object records and counts")
    args = parser.parse_args()
    typing_root = args.typing_root.resolve()
    run_path = typing_root / "run.json"
    run = json.loads(run_path.read_text()) if run_path.is_file() else {}
    aligned_root = args.aligned_root.resolve() if args.aligned_root else Path(
        run.get("aligned_root", S.ALIGNED)).resolve()
    inventory = json.loads(args.inventory.read_text()) if args.inventory.is_file() else {}
    statuses = disease_status_map(inventory.get("sections", []))
    donors = args.donors or run.get("donors") or sorted(
        p.name.split("=", 1)[1] for p in (typing_root / "cells").glob("donor_id=*")
        if (p / "data.parquet").is_file())
    if not donors:
        raise ValueError("no typing donor partitions found")
    args.out.mkdir(parents=True, exist_ok=True)
    manifests = []
    cell_exclusions = read_cell_exclusions(args.cell_exclusions)
    for donor in donors:
        frame, manifest = build_donor(typing_root, aligned_root, str(donor),
                                      statuses.get(donor_id_key(donor), "unavailable"), cell_exclusions)
        if len(frame):
            target = args.out / f"donor_id={donor}"
            target.mkdir(parents=True, exist_ok=True)
            frame.to_parquet(target / "data.parquet", index=False)
        manifests.append(manifest)
        print(f"{donor}: {manifest['object_rows']} endocrine objects; "
              f"{manifest['class_counts']}", flush=True)
    all_objects = pd.concat([
        pd.read_parquet(args.out / f"donor_id={d}" / "data.parquet")
        for d in donors if (args.out / f"donor_id={d}" / "data.parquet").is_file()
    ], ignore_index=True) if any((args.out / f"donor_id={d}" / "data.parquet").is_file()
                                 for d in donors) else pd.DataFrame()
    if len(all_objects) and all_objects["endocrine_object_id"].duplicated().any():
        raise ValueError("cohort composition repeats endocrine_object_id")
    if len(all_objects):
        all_objects.to_parquet(args.out / "endocrine_objects.parquet", index=False)
    data_dictionary = {
        "title": "Endocrine-object composition (no spatial expansion rings)",
        "analysis_status": args.analysis_status,
        "typing_root": str(typing_root), "aligned_root": str(aligned_root),
        "typing_run": str(run_path.resolve()) if run_path.is_file() else None,
        "typing_run_sha256": file_sha256(run_path) if run_path.is_file() else None,
        "input_inventory": str(args.inventory.resolve()) if args.inventory.is_file() else None,
        "input_inventory_sha256": file_sha256(args.inventory) if args.inventory.is_file() else None,
        "composition_builder": str(Path(__file__).resolve()),
        "composition_builder_sha256": file_sha256(Path(__file__).resolve()),
        "cell_exclusion_overlay": str(args.cell_exclusions.resolve()) if args.cell_exclusions else None,
        "cell_exclusion_overlay_sha256": file_sha256(args.cell_exclusions) if args.cell_exclusions else None,
        "explicitly_excluded_single_objects": int((~all_objects.composition_analysis_eligible).sum()) if len(all_objects) else 0,
        "object_key": ["donor_id", "endocrine_object_id"],
        "membership": "true QuPath islet-core polygon membership; no expansion or distance band",
        "classes": "original_object_class preserves Islet, EndocrineCluster, EndocrineSingle",
        "denominator": "all primary-typed endocrine cells in the object; zero yields null ratios",
        "primary_basis": "phenotyped cell bodies per endocrine cell, in six groups",
        "groups": {"endocrine": "Beta/Alpha/Delta and unresolved subtype fractions",
                   "vascular": "Blood vessel cell bodies per endocrine + measured outline coverage",
                   "lymphatic": "Lymphatic cell bodies per endocrine + measured outline coverage",
                   "neural": "Nerve cell bodies per endocrine + measured outline coverage",
                   "stromal": "explicit primary stromal identity only; no expression proxy is recoded",
                   "immune": "broad Immune cell bodies plus rule-based available subtype calls"},
        "supplementary": "fragment-masked aligned mean INS/GCG/SST intensity; refined classes; cell-body overlap with structural object outlines",
        "missingness": "unmeasured channels, unavailable rules, and zero endocrine denominators remain null with availability fields/reasons",
        "per_donor": manifests,
    }
    (args.out / "manifest.json").write_text(json.dumps(data_dictionary, indent=2, default=str) + "\n")
    print(f"wrote {len(all_objects)} endocrine objects -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
