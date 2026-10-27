"""Candidate links between independently segmented cells on one physical section.

This is the shared cell-linking function for tissue-specific same-slide typing.
Protein and RNA calls remain separate evidence; a spatial candidate is not a
verified cell identity or a final type. The caller supplies a transform whose
metadata names the exact source pair and coordinate frames. The metadata's
``reviewed`` flag is a declaration, not proof of gate acceptance or registration
quality; the stage using this function must check the registration receipt.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..same_slide_sources import SameSlidePair
from .manifest import parse_bundle_path
from .sameslide import match_by_centroid
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


def _cells(table: pd.DataFrame, *, name: str, eligible_column: str) -> pd.DataFrame:
    needed = {"cell_uid", "donor_id", "x_um", "y_um", eligible_column}
    missing = sorted(needed - set(table.columns))
    if missing:
        raise PairedTypingError(f"{name}: missing columns {missing}")
    if table.empty or table["cell_uid"].isna().any() or table["cell_uid"].duplicated().any():
        raise PairedTypingError(f"{name}: cell_uid must be present and unique")
    if table[eligible_column].isna().any():
        raise PairedTypingError(f"{name}: missing eligibility is not false")
    xy = table[["x_um", "y_um"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise PairedTypingError(f"{name}: non-finite coordinates")
    return table


def _calls(cells: pd.DataFrame, calls: pd.DataFrame, *, name: str) -> pd.DataFrame:
    if "cell_uid" not in calls or calls["cell_uid"].isna().any():
        raise PairedTypingError(f"{name}: calls need cell_uid")
    if calls["cell_uid"].duplicated().any() or set(calls["cell_uid"]) != set(cells["cell_uid"]):
        raise PairedTypingError(f"{name}: calls must cover exactly the ingested cells once")
    joined = cells.merge(calls, on="cell_uid", how="left", validate="one_to_one",
                         suffixes=("", "_call"))
    for column in ("object_id", "typing_eligible"):
        check = f"{column}_call"
        if check in joined and not joined[column].equals(joined[check]):
            raise PairedTypingError(f"{name}: {column} changed between ingest and typing")
    return joined


def link_candidate_cells(
    pair: SameSlidePair,
    protein_cells: pd.DataFrame,
    protein_calls: pd.DataFrame,
    rna_cells: pd.DataFrame,
    rna_calls: pd.DataFrame,
    transform: Transform,
    *,
    max_dist_um: float = 5.0,
) -> CandidateLinks:
    """Return mutual-nearest candidate links in the protein crop frame.

    The transform maps the Xenium region-local micron frame to the protein
    crop-local micron frame. It must carry exact pair IDs and ``reviewed: true``
    in its metadata. The flag alone does not authenticate review. The result
    contains only cells eligible in both assays;
    the counts report excluded and unmatched cells separately. No candidate
    link is called ground truth or used as a final cell type here.
    """
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
    if not np.isfinite(transform.affine).all() or not 0 < max_dist_um <= 20:
        raise PairedTypingError("invalid transform or cell-distance limit")
    protein = _cells(protein_cells, name="protein", eligible_column="protein_eligible")
    rna = _cells(rna_cells, name="RNA", eligible_column="typing_eligible")
    slide, _ = parse_bundle_path(pair.xenium_output.name)
    if (set(protein["donor_id"]) != {pair.donor_id}
            or set(rna["donor_id"]) != {pair.donor_id}
            or "acquisition_id" not in protein
            or set(protein["acquisition_id"]) != {f"phenocycler:{slide}"}
            or "acquisition_id" not in rna
            or set(rna["acquisition_id"]) != {pair.xenium_acquisition_id}):
        raise PairedTypingError("cell tables do not name the recorded donor and RNA acquisition")
    if "candidate_state" not in protein_calls or "candidate_phenotype" not in protein_calls:
        raise PairedTypingError("protein hierarchy candidates are missing")
    if not any(name.startswith("gene:") for name in rna_calls):
        raise PairedTypingError("RNA identity-gene evidence is missing")
    protein = _calls(protein, protein_calls, name="protein")
    rna = _calls(rna, rna_calls, name="RNA")
    p = protein.loc[protein["protein_eligible"].astype(bool)].reset_index(drop=True)
    x = rna.loc[rna["typing_eligible"].astype(bool)].reset_index(drop=True)
    moving = transform.apply(x[["x_um", "y_um"]].to_numpy(dtype=float))
    if not np.isfinite(moving).all():
        raise PairedTypingError("registration produced non-finite cell coordinates")
    matched = match_by_centroid(moving, p[["x_um", "y_um"]].to_numpy(dtype=float),
                                max_dist_um=max_dist_um)
    columns = ["donor_id", "physical_section", "xenium_acquisition_id",
               "protein_cell_uid", "rna_cell_uid", "distance_um",
               "protein_candidate_phenotype", "protein_candidate_state"]
    genes = sorted(name for name in rna_calls if name.startswith("gene:"))
    columns += genes
    if matched.empty:
        links = pd.DataFrame(columns=columns)
    else:
        pi, xi = matched["fixed_idx"].to_numpy(), matched["moving_idx"].to_numpy()
        links = pd.DataFrame({
            "donor_id": pair.donor_id, "physical_section": pair.protein_entry,
            "xenium_acquisition_id": pair.xenium_acquisition_id,
            "protein_cell_uid": p.iloc[pi]["cell_uid"].to_numpy(),
            "rna_cell_uid": x.iloc[xi]["cell_uid"].to_numpy(),
            "distance_um": matched["distance_um"].to_numpy(),
            "protein_candidate_phenotype": p.iloc[pi]["candidate_phenotype"].to_numpy(),
            "protein_candidate_state": p.iloc[pi]["candidate_state"].to_numpy(),
            **{gene: x.iloc[xi][gene].to_numpy() for gene in genes},
        })
    return CandidateLinks(links, len(protein), len(rna), len(p), len(x),
                          len(p) - len(links), len(x) - len(links))
