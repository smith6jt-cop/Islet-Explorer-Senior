"""Join one mixed pancreas slide's separate tissue typing evidence by native cell ID.

The pancreas hierarchy covers eligible pancreas cells. The lymph-node sheet
covers pLN cells on the same slide. Neither is allowed to make a call for the
other tissue, and every segmented cell remains in the candidate index.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..identity import cell_uid


class MixedTypingError(ValueError):
    """A tissue route lost cells, crossed tissues or collided on native IDs."""


def assemble(aligned: pd.DataFrame, pancreas: pd.DataFrame, pln: pd.DataFrame,
             *, donor_id: str, acquisition_id: str, segmentation_id: str,
             pln_blocks: pd.DataFrame | None = None) -> pd.DataFrame:
    """Make one candidate index, retaining untyped cells and unavailable pLN calls.

    ``pancreas`` is the final step of one historical hierarchical method,
    carrying both ``cell_class`` and ``refined_class`` for audit. ``pln`` is the
    reviewed-sheet candidate, with ``candidate_phenotype`` and
    ``candidate_state``. It does not become a final class here.
    """
    required = {"object_id", "qc_analysis_eligible", "cell_region", "in_alignment_block"}
    if required - set(aligned):
        raise MixedTypingError(f"aligned table lacks {sorted(required - set(aligned))}")
    if aligned.object_id.isna().any() or not aligned.object_id.is_unique:
        raise MixedTypingError("aligned native object_id is missing or repeated")
    if not aligned.object_id.map(lambda value: isinstance(value, str) and bool(value)).all():
        raise MixedTypingError("aligned native object_id must be nonempty text")
    known_regions = {"tissue", "core", "cluster", "single", "pLN"}
    if not set(aligned.cell_region.dropna()) <= known_regions or aligned.cell_region.isna().any():
        raise MixedTypingError("aligned cell_region is missing or has an unknown tissue route")
    for name, frame, columns in (
            ("pancreas", pancreas, {"object_id", "cell_class", "refined_class"}),
            ("pLN", pln, {"object_id", "candidate_phenotype", "candidate_state"})):
        if columns - set(frame):
            raise MixedTypingError(f"{name} calls lack {sorted(columns - set(frame))}")
        if frame.object_id.isna().any() or not frame.object_id.is_unique:
            raise MixedTypingError(f"{name} calls have missing or repeated object_id")
    node = aligned.cell_region.eq("pLN")
    eligible = aligned.qc_analysis_eligible.fillna(False).astype(bool)
    expected_pancreas = set(aligned.loc[eligible & ~node, "object_id"])
    if set(pancreas.object_id) != expected_pancreas:
        raise MixedTypingError("pancreas hierarchy must cover exactly the eligible non-pLN cells")
    if set(pln.object_id) != set(aligned.loc[node, "object_id"]):
        raise MixedTypingError("pLN candidate must cover every recorded pLN cell")
    if not set(pln.candidate_state) <= {"unique_candidate", "ambiguous", "unavailable", "no_match"}:
        raise MixedTypingError("pLN candidate has an unknown state")
    if (pln.loc[~pln.object_id.isin(aligned.loc[eligible & node, "object_id"]),
                "candidate_state"] != "unavailable").any():
        raise MixedTypingError("ineligible pLN cell received a candidate state")

    out = aligned.loc[:, ["object_id", "qc_analysis_eligible", "cell_region",
                          "in_alignment_block"]].copy()
    out.insert(0, "donor_id", donor_id)
    out.insert(1, "cell_uid", [cell_uid(acquisition_id, segmentation_id, str(oid))
                               for oid in out.object_id])
    out.insert(2, "acquisition_id", acquisition_id)
    out.insert(3, "segmentation_id", segmentation_id)
    out["native_cell_id"] = out.object_id
    out["tissue"] = np.where(node, "pancreatic_lymph_node", "pancreas")
    out = out.merge(pancreas[["object_id", "cell_class", "refined_class"]],
                    on="object_id", how="left", validate="one_to_one", sort=False)
    out = out.merge(pln[["object_id", "candidate_phenotype", "candidate_state"]],
                    on="object_id", how="left", validate="one_to_one", sort=False)
    if pln_blocks is not None:
        if ({"object_id", "block"} - set(pln_blocks) or pln_blocks.object_id.isna().any()
                or not pln_blocks.object_id.is_unique
                or set(pln_blocks.object_id) != set(aligned.loc[node, "object_id"])):
            raise MixedTypingError("pLN block sidecar must cover exactly the pLN cells")
        out = out.merge(pln_blocks[["object_id", "block"]].rename(
            columns={"block": "pln_block"}), on="object_id", how="left",
            validate="one_to_one", sort=False)
    else:
        out["pln_block"] = pd.NA
    out["typing_state"] = np.where(~eligible, "unavailable_qc",
                                   np.where(node, out.candidate_state,
                                            "pancreas_hierarchy"))
    if len(out) != len(aligned) or not out.cell_uid.is_unique:
        raise MixedTypingError("candidate index changed cell count or repeated cell_uid")
    return out
