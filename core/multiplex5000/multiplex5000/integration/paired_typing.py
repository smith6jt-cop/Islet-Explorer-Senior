"""Candidate links between independently segmented cells on one physical section.

This is the shared cell-linking function for tissue-specific same-slide typing.
Protein and RNA calls remain separate evidence; a spatial candidate is not a
verified cell identity or a final type. The caller supplies a transform whose
metadata names the exact source pair and coordinate frames. The metadata's
``reviewed`` flag is a declaration, not proof of gate acceptance or registration
quality; the stage using this function must check the registration receipt.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..same_slide_sources import SameSlidePair
from .manifest import parse_bundle_path
from .transform import Transform


class PairedTypingError(ValueError):
    """The two assay tables or their registration cannot be paired as declared."""


@dataclass(frozen=True)
class CandidateLinks:
    links: pd.DataFrame
    protein_cells: int
    rna_cells: int
    protein_eligible: int
    rna_eligible: int
    unmatched_protein_eligible: int
    unmatched_rna_eligible: int
    candidates: pd.DataFrame = field(default_factory=pd.DataFrame)
    cell_dispositions: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict = field(default_factory=dict)


IDENTITY_COLUMNS = {"candidate_phenotype", "cell_class", "refined_class", "cell_type",
                    "label", "identity", "lineage", "subtype", "identity_core"}
SOURCE_COLUMNS = {"object_id", "donor_id", "tissue", "acquisition_id", "x_um", "y_um",
                  "typing_eligible", "protein_eligible", "coordinate_frame", "coordinate_units",
                  "segmentation_id", "segmentation_version"}


def _cells(table: pd.DataFrame, *, name: str, eligible_column: str) -> pd.DataFrame:
    if table.columns.duplicated().any():
        raise PairedTypingError(f"{name}: duplicate cell columns")
    needed = {"cell_uid", "donor_id", "x_um", "y_um", eligible_column}
    missing = sorted(needed - set(table.columns))
    if missing:
        raise PairedTypingError(f"{name}: missing columns {missing}")
    if (table.empty or table["cell_uid"].isna().any()
            or not table["cell_uid"].map(lambda value: isinstance(value, str) and bool(value)).all()
            or table["cell_uid"].duplicated().any()):
        raise PairedTypingError(f"{name}: cell_uid must be present and unique")
    if table[eligible_column].isna().any():
        raise PairedTypingError(f"{name}: missing eligibility is not false")
    if not table[eligible_column].isin([True, False]).all():
        raise PairedTypingError(f"{name}: eligibility must be boolean")
    xy = table[["x_um", "y_um"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise PairedTypingError(f"{name}: non-finite coordinates")
    return table


def _calls(cells: pd.DataFrame, calls: pd.DataFrame, *, name: str) -> pd.DataFrame:
    if calls.columns.duplicated().any():
        raise PairedTypingError(f"{name}: duplicate call columns")
    if "cell_uid" not in calls or calls["cell_uid"].isna().any():
        raise PairedTypingError(f"{name}: calls need cell_uid")
    if calls["cell_uid"].duplicated().any() or set(calls["cell_uid"]) != set(cells["cell_uid"]):
        raise PairedTypingError(f"{name}: calls must cover exactly the ingested cells once")
    joined = cells.merge(calls, on="cell_uid", how="left", validate="one_to_one",
                         suffixes=("", "_call"))
    for column in SOURCE_COLUMNS & set(cells) & set(calls):
        check = f"{column}_call"
        equal = (joined[column] == joined[check]) | (joined[column].isna() & joined[check].isna())
        if not equal.fillna(False).all():
            raise PairedTypingError(f"{name}: {column} changed between ingest and typing")
    return joined


def validate_pair_transform(pair: SameSlidePair, transform: Transform) -> None:
    """Check declared same-section identity and the transform's finite geometry."""
    if pair.relationship != "same_section_after_phenocycler":
        raise PairedTypingError("cell links require a recorded same-section relationship")
    if (pair.protein_frame, pair.xenium_frame) != ("crop_local_um", "region_local_um"):
        raise PairedTypingError("unsupported or missing coordinate frames")
    expected = {"protein_entry": pair.protein_entry,
                "xenium_acquisition_id": pair.xenium_acquisition_id,
                "moving_frame": pair.xenium_frame, "fixed_frame": pair.protein_frame,
                "reviewed": True}
    if any(transform.meta.get(key) != value for key, value in expected.items()):
        raise PairedTypingError("registration does not name this exact reviewed source pair")
    if transform.meta.get("reviewed") is not True:
        raise PairedTypingError("registration review declaration must be boolean")
    if (not np.isfinite(transform.affine).all()
            or (transform.displacement is not None
                and not np.isfinite(transform.displacement).all())):
        raise PairedTypingError("registration has non-finite geometry")
    if np.linalg.det(transform.affine[:, :2]) == 0:
        raise PairedTypingError("registration has singular geometry")


def validate_pair_cells(pair: SameSlidePair, protein_cells: pd.DataFrame,
                        rna_cells: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Require both assay inventories to name this tissue and exact physical section."""
    protein = _cells(protein_cells, name="protein", eligible_column="protein_eligible")
    rna = _cells(rna_cells, name="RNA", eligible_column="typing_eligible")
    slide, _ = parse_bundle_path(pair.xenium_output.name)
    if (set(protein["donor_id"]) != {pair.donor_id}
            or set(rna["donor_id"]) != {pair.donor_id}
            or "tissue" not in protein or set(protein["tissue"]) != {pair.tissue}
            or "tissue" not in rna or set(rna["tissue"]) != {pair.tissue}
            or "acquisition_id" not in protein
            or set(protein["acquisition_id"]) != {f"phenocycler:{slide}"}
            or "acquisition_id" not in rna
            or set(rna["acquisition_id"]) != {pair.xenium_acquisition_id}):
        raise PairedTypingError("cell tables do not name the recorded donor, tissue or acquisition")
    if (not protein["cell_uid"].str.startswith(
            f"phenocycler:{slide}:{pair.protein_entry}:").all()
            or not rna["cell_uid"].str.startswith(
                f"{pair.xenium_acquisition_id}:cells:").all()):
        raise PairedTypingError("cell keys do not name the recorded physical section")
    for name, cells, frame in [("protein", protein, pair.protein_frame), ("RNA", rna, pair.xenium_frame)]:
        if "coordinate_frame" in cells and set(cells["coordinate_frame"]) != {frame}:
            raise PairedTypingError(f"{name}: mixed or incorrect source coordinate frame")
        if "coordinate_units" in cells and set(cells["coordinate_units"]) != {"um"}:
            raise PairedTypingError(f"{name}: source coordinate units must be um")
        for column in ["segmentation_id", "segmentation_version"]:
            if column in cells and (cells[column].isna().any() or cells[column].nunique() != 1):
                raise PairedTypingError(f"{name}: mixed or missing {column}")
    return protein, rna


def _uncertainty(value: dict | None, pair: SameSlidePair) -> dict:
    if value is None:
        return {"status": "unavailable", "bound_um": None, "evidence_ref": None,
                "reason": "No registration uncertainty bound was provided."}
    if not isinstance(value, dict):
        raise PairedTypingError("registration uncertainty must be an explicit bound record")
    bound = value.get("bound_um")
    if (isinstance(bound, bool) or not isinstance(bound, (int, float))
            or not np.isfinite(bound) or bound < 0 or value.get("frame") != pair.protein_frame
            or value.get("kind") != "radial_displacement_bound"
            or not isinstance(value.get("evidence_ref"), str) or not value["evidence_ref"].strip()):
        raise PairedTypingError("uncertainty needs a finite nonnegative radial bound, fixed frame and evidence reference")
    return {**value, "bound_um": float(bound), "status": "provided",
            "method": "conservative distance intervals and nearest-neighbour margins",
            "scope": "Each transformed RNA centroid may move within the supplied bound; segmentation error is not modeled."}


def _geometry(moving: np.ndarray, fixed: np.ndarray, cap: float, uncertainty: dict) -> pd.DataFrame:
    from scipy.spatial import cKDTree

    cols = ["moving_idx", "fixed_idx", "distance_um", "within_nominal_cap", "nominal_state",
            "rna_candidate_count", "protein_candidate_count", "rna_nearest_tie_count",
            "protein_nearest_tie_count", "distance_lower_um", "distance_upper_um",
            "stable_under_provided_bound", "uncertainty_state"]
    if not len(moving) or not len(fixed):
        return pd.DataFrame(columns=cols)
    bound = uncertainty["bound_um"]
    tm, tf = cKDTree(moving), cKDTree(fixed)
    pairs = tm.sparse_distance_matrix(tf, cap + (bound or 0), output_type="ndarray")
    rows = pd.DataFrame({"moving_idx": pairs["i"], "fixed_idx": pairs["j"], "distance_um": pairs["v"]})
    if rows.empty:
        return pd.DataFrame(columns=cols)
    rows = rows.sort_values(["moving_idx", "fixed_idx"]).reset_index(drop=True)
    rows["within_nominal_cap"] = rows["distance_um"] <= cap
    nominal = rows.loc[rows["within_nominal_cap"]].copy()
    for index, prefix in [("moving_idx", "rna"), ("fixed_idx", "protein")]:
        minimum = nominal.groupby(index)["distance_um"].transform("min")
        ties = nominal.loc[nominal["distance_um"] == minimum].groupby(index).size()
        rows[f"{prefix}_candidate_count"] = rows[index].map(nominal.groupby(index).size()).fillna(0).astype(int)
        rows[f"{prefix}_nearest_tie_count"] = rows[index].map(ties).fillna(0).astype(int)
        rows[f"_{prefix}_nearest"] = rows["distance_um"] == rows[index].map(nominal.groupby(index)["distance_um"].min())
    both_nearest = rows["_rna_nearest"] & rows["_protein_nearest"]
    unique = both_nearest & rows["rna_nearest_tie_count"].eq(1) & rows["protein_nearest_tie_count"].eq(1)
    rows["nominal_state"] = np.select(
        [unique, both_nearest, rows["within_nominal_cap"]],
        ["unique_mutual_candidate", "ambiguous_nearest", "nonmutual_alternative"],
        default="outside_nominal_cap")
    rows["distance_lower_um"] = None
    rows["distance_upper_um"] = None
    rows["stable_under_provided_bound"] = pd.Series(pd.NA, index=rows.index, dtype="boolean")
    rows["uncertainty_state"] = "unavailable"
    if bound is not None:
        rows["distance_lower_um"] = np.maximum(0, rows["distance_um"] - bound)
        rows["distance_upper_um"] = rows["distance_um"] + bound
        dm, _ = tf.query(moving, k=2)
        df, _ = tm.query(fixed, k=2)
        margin_m = dm[rows["moving_idx"], 1] - rows["distance_um"].to_numpy()
        margin_f = df[rows["fixed_idx"], 1] - rows["distance_um"].to_numpy()
        stable = unique & (rows["distance_upper_um"] <= cap) & (margin_m > 2 * bound) & (margin_f > 2 * bound)
        rows["stable_under_provided_bound"] = pd.array(stable, dtype="boolean")
        rows["uncertainty_state"] = np.where(stable, "stable_under_provided_bound", "sensitive_or_alternative")
    return rows[cols]


def _evidence(calls: pd.DataFrame, ids: pd.Series, prefix: str) -> dict:
    selected = calls.set_index("cell_uid").loc[ids]
    return {(col if prefix == "rna" and col.startswith("gene:") else f"{prefix}_{col}"): selected[col].to_numpy()
            for col in calls if col != "cell_uid"}


def _compare(links: pd.DataFrame, mapping: dict | None, pair: SameSlidePair,
             protein_calls: pd.DataFrame, rna_calls: pd.DataFrame) -> dict:
    links["identity_concordance"] = "unavailable"
    if mapping is None:
        return {"status": "unavailable", "reason": "No explicit hierarchy mapping was supplied."}
    if (not isinstance(mapping, dict) or mapping.get("tissue") != pair.tissue
            or not isinstance(mapping.get("level"), str) or not mapping["level"].strip()):
        raise PairedTypingError("hierarchy mapping needs the exact tissue and an explicit level")
    for assay, calls in [("protein", protein_calls), ("rna", rna_calls)]:
        column, vocabulary = mapping.get(f"{assay}_column"), mapping.get(assay)
        if column not in IDENTITY_COLUMNS or column not in calls or not isinstance(vocabulary, dict) or not vocabulary:
            raise PairedTypingError(f"hierarchy mapping lacks the {assay} call column or vocabulary")
        if any(not isinstance(v, str) or not v.strip() for item in vocabulary.items() for v in item):
            raise PairedTypingError("hierarchy vocabulary keys and values must be nonempty strings")
        values = links[f"{assay}_{column}"]
        available = values.notna() & values.ne("")
        state = mapping.get(f"{assay}_state_column")
        if state is not None:
            accepted = mapping.get(f"{assay}_available_states")
            if state not in calls or not isinstance(accepted, list) or not accepted or any(not isinstance(v, str) for v in accepted):
                raise PairedTypingError("hierarchy state filtering needs a declared column and available states")
            available &= links[f"{assay}_{state}"].isin(accepted)
        links[f"{assay}_comparison_label"] = values.map(vocabulary).where(available)
        links[f"_{assay}_comparison_available"] = available
    both = links["protein_comparison_label"].notna() & links["rna_comparison_label"].notna()
    comparable = links["_protein_comparison_available"] & links["_rna_comparison_available"]
    links.loc[comparable & ~both, "identity_concordance"] = "unmapped"
    links.loc[both, "identity_concordance"] = np.where(
        links.loc[both, "protein_comparison_label"] == links.loc[both, "rna_comparison_label"], "agree", "disagree")
    links["comparison_level"] = mapping["level"]
    links.drop(columns=["_protein_comparison_available", "_rna_comparison_available"], inplace=True)
    return {"status": "declared", "mapping": mapping, "scope": "Comparison evidence only; original calls are retained."}


def link_candidate_cells(
    pair: SameSlidePair,
    protein_cells: pd.DataFrame,
    protein_calls: pd.DataFrame,
    rna_cells: pd.DataFrame,
    rna_calls: pd.DataFrame,
    transform: Transform,
    *,
    max_dist_um: float = 5.0,
    registration_uncertainty: dict | None = None,
    hierarchy_map: dict | None = None,
) -> CandidateLinks:
    """Return mutual-nearest candidate links in the protein crop frame.

    The transform maps the Xenium region-local micron frame to the protein
    crop-local micron frame. It must carry exact pair IDs and ``reviewed: true``
    in its metadata. The flag alone does not authenticate review. The result
    links contain only unique mutual-nearest eligible cells. All nearby edges,
    including ambiguous ties, remain in ``candidates``; ``cell_dispositions``
    accounts for every ingested cell. Optional uncertainty is a provided radial
    bound in the fixed frame, never a fit to candidate agreement. An explicit
    hierarchy map permits comparison while retaining both original calls.
    """
    validate_pair_transform(pair, transform)
    if not 0 < max_dist_um <= 20:
        raise PairedTypingError("invalid cell-distance limit")
    protein, rna = validate_pair_cells(pair, protein_cells, rna_cells)
    if not IDENTITY_COLUMNS & set(protein_calls):
        raise PairedTypingError("protein hierarchy candidates are missing")
    if not any(name.startswith("gene:") for name in rna_calls) and not IDENTITY_COLUMNS & set(rna_calls):
        raise PairedTypingError("RNA identity-gene evidence is missing")
    protein = _calls(protein, protein_calls, name="protein")
    rna = _calls(rna, rna_calls, name="RNA")
    # Optional source metadata can first appear in a typing table. Validate it
    # too, rather than silently ignoring a call-only coordinate frame.
    protein, rna = validate_pair_cells(pair, protein, rna)
    p = protein.loc[protein["protein_eligible"].astype(bool)].sort_values("cell_uid").reset_index(drop=True)
    x = rna.loc[rna["typing_eligible"].astype(bool)].sort_values("cell_uid").reset_index(drop=True)
    moving = transform.apply(x[["x_um", "y_um"]].to_numpy(dtype=float))
    if not np.isfinite(moving).all():
        raise PairedTypingError("registration produced non-finite cell coordinates")
    uncertainty = _uncertainty(registration_uncertainty if registration_uncertainty is not None
                               else transform.meta.get("registration_uncertainty"), pair)
    geometry = _geometry(moving, p[["x_um", "y_um"]].to_numpy(dtype=float), max_dist_um, uncertainty)
    pi = geometry["fixed_idx"].to_numpy(dtype=int)
    xi = geometry["moving_idx"].to_numpy(dtype=int)
    candidates = pd.DataFrame({
        "protein_cell_uid": p.iloc[pi]["cell_uid"].to_numpy(),
        "rna_cell_uid": x.iloc[xi]["cell_uid"].to_numpy(),
        **_evidence(protein_calls, p.iloc[pi]["cell_uid"], "protein"),
        **_evidence(rna_calls, x.iloc[xi]["cell_uid"], "rna"),
    })
    reserved = (set(geometry) - {"moving_idx", "fixed_idx"}) | {
        "donor_id", "physical_section", "xenium_acquisition_id", "segmentation_relation_state",
        "registration_uncertainty_um", "protein_identity_call_present", "rna_identity_call_present",
        "protein_comparison_label", "rna_comparison_label"}
    collisions = reserved & set(candidates)
    if collisions:
        raise PairedTypingError(f"independent call columns collide with link metadata: {sorted(collisions)}")
    for column in geometry:
        if column not in {"moving_idx", "fixed_idx"}:
            candidates[column] = geometry[column].reset_index(drop=True)
    candidates["donor_id"] = pair.donor_id
    candidates["physical_section"] = pair.protein_entry
    candidates["xenium_acquisition_id"] = pair.xenium_acquisition_id
    candidates["segmentation_relation_state"] = "unavailable_without_overlap_evidence"
    candidates["registration_uncertainty_um"] = uncertainty["bound_um"]
    for assay, calls in [("protein", protein_calls), ("rna", rna_calls)]:
        available = pd.Series(False, index=candidates.index)
        for column in IDENTITY_COLUMNS & set(calls):
            available |= candidates[f"{assay}_{column}"].notna() & candidates[f"{assay}_{column}"].ne("")
        candidates[f"{assay}_identity_call_present"] = available
    links = candidates.loc[candidates["nominal_state"] == "unique_mutual_candidate"].copy().reset_index(drop=True)
    comparison = _compare(links, hierarchy_map, pair, protein_calls, rna_calls)
    dispositions = []
    for assay, cells, eligible, calls, index, frame in [
            ("protein", protein, p, protein_calls, "fixed_idx", pair.protein_frame),
            ("rna", rna, x, rna_calls, "moving_idx", pair.xenium_frame)]:
        cells = cells.sort_values("cell_uid").reset_index(drop=True)
        uid_to_index = pd.Series(eligible.index, index=eligible["cell_uid"])
        indices = cells["cell_uid"].map(uid_to_index)
        is_eligible = indices.notna()
        nominal_count = geometry.groupby(index)[f"{assay}_candidate_count"].first()
        tie_count = geometry.groupby(index)[f"{assay}_nearest_tie_count"].first()
        possible_count = geometry.groupby(index).size()
        count = indices.map(nominal_count).fillna(0).astype(int)
        possible = indices.map(possible_count).fillna(0).astype(int)
        ties = indices.map(tie_count).fillna(0).astype(int)
        ambiguous = indices.isin(geometry.loc[geometry["nominal_state"] == "ambiguous_nearest", index]) | ties.gt(1)
        linked = cells["cell_uid"].isin(links[f"{assay}_cell_uid"])
        reason = np.select([~is_eligible, linked, ambiguous, count.eq(0) & possible.gt(0), count.eq(0)],
                           ["excluded_by_assay_eligibility", "unique_mutual_candidate", "ambiguous_nearest",
                            "only_uncertainty_candidates", "no_candidate_within_cap"], default="not_mutual_nearest")
        dispositions.append(pd.DataFrame({"assay": assay, "cell_uid": cells["cell_uid"],
            "eligible": is_eligible, "link_state": reason,
            "candidate_count": count.where(is_eligible).astype("Int64"),
            "possible_candidate_count": possible.where(is_eligible).astype("Int64"),
            "nearest_tie_count": ties.where(is_eligible).astype("Int64"),
            "coordinate_frame": frame,
            "segmentation_id": cells.get("segmentation_id", pd.Series(None, index=cells.index)),
            "segmentation_version": cells.get("segmentation_version", pd.Series(None, index=cells.index)),
            **_evidence(calls, cells["cell_uid"], assay)}))
    metadata = {"registration_uncertainty": uncertainty, "hierarchy_comparison": comparison,
                "coordinate_frame_evidence": {
                    "protein": "cell_metadata_and_pair" if "coordinate_frame" in protein else "declared_pair_only",
                    "rna": "cell_metadata_and_pair" if "coordinate_frame" in rna else "declared_pair_only"},
                "segmentation_provenance": {
                    assay: {column: cells[column].tolist()[0] if column in cells else None
                            for column in ["segmentation_id", "segmentation_version"]}
                    for assay, cells in [("protein", protein), ("rna", rna)]},
                "tie_rule": "Equal nearest distances in either direction remain ambiguous; no row-order tie break.",
                "segmentation_relation": "Unavailable without overlap evidence; multiple candidates are not asserted splits or merges.",
                "call_presence": "A nonempty supplied label is present; its independent state/confidence fields remain authoritative, not an inferred resolved class.",
                "method": "unique_mutual_nearest_centroid_candidate", "max_dist_um": max_dist_um}
    return CandidateLinks(links, len(protein), len(rna), len(p), len(x),
                          len(p) - len(links), len(x) - len(links), candidates,
                          pd.concat(dispositions, ignore_index=True), metadata)
