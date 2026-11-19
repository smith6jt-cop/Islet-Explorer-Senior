"""Carry protein cell classes onto Xenium cells by measured polygon overlap.

The two assays segment the same physical section independently and disagree by design, so there
is no one-to-one cell match to find. Xenium polygons are warped forward into the protein frame
and every positive-area overlap is kept as an edge. Three things are carried over those edges:

* The class, by largest shared area. A weighted mean of labels has no meaning, so the Xenium cell
  takes the class of the protein cell it overlaps most. ``call_share`` and ``call_purity`` record
  what that is worth.
* The lineage composition, as each lineage's share of the overlapped area. Region calling uses
  this, because a cell straddling a B/T boundary should count partly toward each.
* Per-marker positivity, as an area-weighted share of the overlap with a measured call.

Missing stays missing. An empty final class means Cell Typing found no supported class, not that
the lineage markers were negative, so such a cell is ``empty_final_class``, never a negative. A
``-1`` marker call is unavailable and contributes neither numerator nor denominator. A Xenium cell
that overlaps no protein cell is ``no_overlap``; one that overlaps only protein cells outside the
typed table (excluded upstream) is ``untyped_only``. Both have NaN measures.

``call_share`` and ``overlap_fraction`` are over typed overlap. ``call_share_all``,
``overlap_fraction_all`` and ``untyped_overlap_share`` count every protein cell, so a large untyped
neighbour lowers confidence instead of vanishing from the denominator.

This moves ``overlap_edges`` and ``transfer_calls`` from ``scripts/sameslide`` into the package and
keeps their operators. The transform is forward only, because ``Transform.invert`` drops any
displacement field. Transferred names keep their own columns and never replace an RNA-derived name.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import shapely

from .transform import Transform

UNAVAILABLE_CALL = -1
STATUS_AWARDED = "awarded"
STATUS_EMPTY = "empty_final_class"
STATUS_NO_OVERLAP = "no_overlap"
STATUS_UNTYPED = "untyped_only"


class CallTransferError(ValueError):
    """The inputs cannot support a defined overlap transfer."""


def _require_unique(values: np.ndarray, name: str) -> np.ndarray:
    values = np.asarray(values).astype(str)
    if (values == "").any() or (values == "None").any() or (values == "nan").any():
        raise CallTransferError(f"{name} has a null or empty identifier")
    if len(np.unique(values)) != len(values):
        raise CallTransferError(f"{name} has a duplicate identifier")
    return values


def _valid(polygons: np.ndarray) -> np.ndarray:
    polygons = np.array(polygons, dtype=object)
    bad = ~shapely.is_valid(polygons)
    if bad.any():
        polygons[bad] = shapely.make_valid(polygons[bad])
    return polygons


def overlap_edges(xenium_ids, xenium_polygons, protein_ids, protein_polygons,
                  transform: Transform) -> tuple[pd.DataFrame, np.ndarray]:
    """Positive-area edges in the protein frame, plus each Xenium cell's warped area.

    Sorted by Xenium ID then protein ID, so the result never depends on input row order.
    """
    xenium_ids = _require_unique(xenium_ids, "xenium_ids")
    protein_ids = _require_unique(protein_ids, "protein_ids")
    warped = _valid(shapely.transform(np.asarray(xenium_polygons, dtype=object), transform.apply))
    protein = _valid(protein_polygons)
    xi, pi = shapely.STRtree(protein).query(warped, predicate="intersects")
    area = shapely.area(shapely.intersection(warped[xi], protein[pi]))
    keep = np.isfinite(area) & (area > 0)
    edges = pd.DataFrame({"cell_id": xenium_ids[xi[keep]], "object_id": protein_ids[pi[keep]],
                          "overlap_area_um2": area[keep]})
    return (edges.sort_values(["cell_id", "object_id"], kind="stable").reset_index(drop=True),
            shapely.area(warped))


def _final(values: pd.Series) -> pd.Series:
    """An empty string in a typing table means no class. Normalize it to missing."""
    return values.astype(object).where(values.notna() & (values.astype(str) != ""), other=None)


def transfer_classes(xenium_ids, edges: pd.DataFrame, typed: pd.DataFrame, warped_area,
                     *, class_column: str = "refined_class", lineage_column: str = "cell_class",
                     marker_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    """One row per Xenium cell, in ``xenium_ids`` order. See the module docstring for meanings."""
    xenium_ids = _require_unique(xenium_ids, "xenium_ids")
    area = np.asarray(warped_area, dtype=float)
    if area.shape != (len(xenium_ids),) or not np.isfinite(area).all() or (area <= 0).any():
        raise CallTransferError("warped Xenium areas must be finite, positive, and one per cell")
    if not np.isfinite(edges.overlap_area_um2).all() or (edges.overlap_area_um2 <= 0).any():
        raise CallTransferError("overlap areas must be finite and positive")
    if edges.duplicated(["cell_id", "object_id"]).any():
        raise CallTransferError("duplicate overlap edge")
    typed = typed.copy()
    typed["object_id"] = _require_unique(typed["object_id"], "typed object_id")
    typed["_class"] = _final(typed[class_column])
    typed["_lineage"] = _final(typed[lineage_column])

    index = pd.Index(xenium_ids)
    n = len(xenium_ids)
    all_xi = index.get_indexer(edges.cell_id)
    if (all_xi < 0).any():
        raise CallTransferError("overlap edges name an unknown Xenium cell")
    total_all = np.bincount(all_xi, weights=edges.overlap_area_um2.to_numpy(float), minlength=n)
    # Protein cells outside the typed table were excluded upstream: they never name a class, but
    # their area is kept in the *_all measures.
    kept = edges[edges.object_id.isin(typed.object_id)]
    xi = index.get_indexer(kept.cell_id)
    weights = kept.overlap_area_um2.to_numpy(float)
    total = np.bincount(xi, weights=weights, minlength=n)

    joined = kept.merge(typed[["object_id", "_class", "_lineage", *marker_columns]], on="object_id",
                        validate="many_to_one")
    # Equal areas resolve by the immutable protein ID, never by input row order.
    winner = (joined.sort_values(["overlap_area_um2", "object_id"], ascending=[False, True], kind="stable")
              .drop_duplicates("cell_id").set_index("cell_id").reindex(xenium_ids))

    def share(numerator: np.ndarray) -> np.ndarray:
        return np.divide(numerator, total, out=np.full(n, np.nan), where=total > 0)

    joined_xi = index.get_indexer(joined.cell_id)
    joined_area = joined.overlap_area_um2.to_numpy(float)
    winning_class = winner["_class"].to_numpy(dtype=object)
    same = joined["_class"].notna().to_numpy() & (joined["_class"].to_numpy(dtype=object)
                                                  == winning_class[joined_xi])
    has_overlap = total > 0
    has_class = pd.notna(winning_class)
    status = np.where(total_all <= 0, STATUS_NO_OVERLAP,
                      np.where(~has_overlap, STATUS_UNTYPED,
                               np.where(has_class, STATUS_AWARDED, STATUS_EMPTY)))
    winner_area = np.nan_to_num(winner["overlap_area_um2"].to_numpy(float), nan=0.0)

    out = pd.DataFrame({
        "cell_id": xenium_ids,
        "protein_transferred_class": winning_class,
        "protein_transferred_lineage": winner["_lineage"].to_numpy(dtype=object),
        "protein_winning_object_id": winner["object_id"].to_numpy(dtype=object),
        "transfer_status": status,
        "call_share": share(winner_area),
        "call_purity": share(np.bincount(joined_xi[same], weights=joined_area[same], minlength=n)),
        "n_protein_cells": np.bincount(xi, minlength=n),
        "overlap_fraction": total / area,
        "call_share_all": np.divide(winner_area, total_all, out=np.full(n, np.nan), where=has_overlap),
        "overlap_fraction_all": total_all / area,
        "untyped_overlap_share": np.divide(total_all - total, total_all, out=np.full(n, np.nan),
                                           where=total_all > 0),
        "n_protein_cells_all": np.bincount(all_xi, minlength=n),
    })
    out["call_purity"] = out["call_purity"].where(has_class)

    for lineage in sorted(joined["_lineage"].dropna().unique()):
        chosen = joined["_lineage"].to_numpy(dtype=object) == lineage
        out[f"lineage_share__{lineage}"] = share(
            np.bincount(joined_xi[chosen], weights=joined_area[chosen], minlength=n))

    for marker in marker_columns:
        calls = joined[marker].to_numpy()
        measured = (calls != UNAVAILABLE_CALL) & pd.notna(calls)
        denominator = np.bincount(joined_xi[measured], weights=joined_area[measured], minlength=n)
        positive = measured & (calls == 1)
        numerator = np.bincount(joined_xi[positive], weights=joined_area[positive], minlength=n)
        name = marker.removeprefix("call__")
        out[f"positive_share__{name}"] = np.divide(numerator, denominator, out=np.full(n, np.nan),
                                                   where=denominator > 0)
    return out
