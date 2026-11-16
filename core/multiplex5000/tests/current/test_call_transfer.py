"""Class transfer by overlap: the semantics Region annotation and the gene comparison rely on."""

import numpy as np
import pandas as pd
import pytest
import shapely

from multiplex5000.integration.call_transfer import (
    STATUS_AWARDED, STATUS_EMPTY, STATUS_NO_OVERLAP, STATUS_UNTYPED, CallTransferError, overlap_edges,
    transfer_classes,
)
from multiplex5000.integration.transform import Transform

IDENTITY = Transform()


def box(x0, y0, x1, y1):
    return shapely.box(x0, y0, x1, y1)


def typed(rows, **markers):
    frame = pd.DataFrame(rows, columns=["object_id", "refined_class", "cell_class"])
    for name, values in markers.items():
        frame[name] = values
    return frame


def run(xenium, protein, table, transform=IDENTITY, markers=()):
    ids = [f"x{i}" for i in range(len(xenium))]
    edges, area = overlap_edges(ids, np.array(xenium, dtype=object), [p[0] for p in protein],
                                np.array([p[1] for p in protein], dtype=object), transform)
    return transfer_classes(ids, edges, table, area, marker_columns=markers).set_index("cell_id")


def test_largest_shared_area_wins_and_is_scored():
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 7, 10)), ("p2", box(7, 0, 10, 10))],
              typed([("p1", "Macrophage", "Myeloid"), ("p2", "B_cell", "B_cell")]))
    row = out.loc["x0"]
    assert row.protein_transferred_class == "Macrophage"
    assert row.protein_transferred_lineage == "Myeloid"
    assert row.transfer_status == STATUS_AWARDED
    assert row.call_share == pytest.approx(0.7)
    assert row.call_purity == pytest.approx(0.7)
    assert row.overlap_fraction == pytest.approx(1.0)
    assert row.n_protein_cells == 2


def test_purity_counts_every_protein_cell_with_the_winning_class():
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 5, 10)), ("p2", box(5, 0, 8, 10)), ("p3", box(8, 0, 10, 10))],
              typed([("p1", "B_cell", "B_cell"), ("p2", "B_cell", "B_cell"), ("p3", "Treg", "T_cell")]))
    assert out.loc["x0"].call_share == pytest.approx(0.5)
    assert out.loc["x0"].call_purity == pytest.approx(0.8)


def test_lineage_shares_keep_a_straddling_cell_in_both():
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 6, 10)), ("p2", box(6, 0, 10, 10))],
              typed([("p1", "B_cell", "B_cell"), ("p2", "Treg", "T_cell")]))
    assert out.loc["x0"].lineage_share__B_cell == pytest.approx(0.6)
    assert out.loc["x0"].lineage_share__T_cell == pytest.approx(0.4)


def test_an_empty_final_class_is_unavailable_never_negative():
    out = run([box(0, 0, 10, 10)], [("p1", box(0, 0, 10, 10))], typed([("p1", "", "")]))
    row = out.loc["x0"]
    assert row.transfer_status == STATUS_EMPTY
    assert pd.isna(row.protein_transferred_class)
    assert pd.isna(row.call_purity)
    assert row.overlap_fraction == pytest.approx(1.0)


def test_empty_winner_still_records_the_lineage_the_cell_overlaps():
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 7, 10)), ("p2", box(7, 0, 10, 10))],
              typed([("p1", "", ""), ("p2", "B_cell", "B_cell")]))
    row = out.loc["x0"]
    assert row.transfer_status == STATUS_EMPTY
    assert row.lineage_share__B_cell == pytest.approx(0.3)


def test_no_overlap_is_missing_not_zero():
    out = run([box(100, 100, 110, 110)], [("p1", box(0, 0, 10, 10))], typed([("p1", "B_cell", "B_cell")]))
    row = out.loc["x0"]
    assert row.transfer_status == STATUS_NO_OVERLAP
    assert np.isnan(row.call_share)
    assert pd.isna(row.protein_transferred_class)
    assert pd.isna(row.call_purity)
    assert row.overlap_fraction == 0
    assert row.n_protein_cells == 0


def test_cells_outside_the_typed_table_take_no_part():
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 7, 10)), ("p2", box(7, 0, 10, 10))],
              typed([("p2", "B_cell", "B_cell")]))
    row = out.loc["x0"]
    assert row.protein_transferred_class == "B_cell"
    assert row.call_share == pytest.approx(1.0)
    assert row.overlap_fraction == pytest.approx(0.3)
    # The untyped neighbour still counts against confidence in the *_all measures.
    assert row.call_share_all == pytest.approx(0.3)
    assert row.overlap_fraction_all == pytest.approx(1.0)
    assert row.untyped_overlap_share == pytest.approx(0.7)
    assert row.n_protein_cells == 1 and row.n_protein_cells_all == 2


def test_overlap_with_only_untyped_cells_is_not_called_no_overlap():
    out = run([box(0, 0, 10, 10), box(100, 100, 110, 110)],
              [("p1", box(0, 0, 10, 10)), ("p2", box(50, 50, 60, 60))],
              typed([("p2", "B_cell", "B_cell")]))
    row = out.loc["x0"]
    assert row.transfer_status == STATUS_UNTYPED
    assert pd.isna(row.protein_transferred_class) and np.isnan(row.call_share) and np.isnan(row.call_share_all)
    assert row.overlap_fraction == 0 and row.overlap_fraction_all == pytest.approx(1.0)
    assert row.untyped_overlap_share == pytest.approx(1.0)
    far = out.loc["x1"]
    assert far.transfer_status == STATUS_NO_OVERLAP and np.isnan(far.untyped_overlap_share)


def test_unavailable_marker_calls_leave_the_denominator():
    table = typed([("p1", "B_cell", "B_cell"), ("p2", "Treg", "T_cell"), ("p3", "Treg", "T_cell")],
                  call__CD20=[1, 0, -1])
    out = run([box(0, 0, 10, 10)],
              [("p1", box(0, 0, 2, 10)), ("p2", box(2, 0, 4, 10)), ("p3", box(4, 0, 10, 10))],
              table, markers=("call__CD20",))
    assert out.loc["x0"].positive_share__CD20 == pytest.approx(0.5)


def test_a_marker_with_no_measured_overlap_is_missing():
    table = typed([("p1", "B_cell", "B_cell")], call__CD20=[-1])
    out = run([box(0, 0, 10, 10)], [("p1", box(0, 0, 10, 10))], table, markers=("call__CD20",))
    assert np.isnan(out.loc["x0"].positive_share__CD20)


def test_ties_resolve_by_protein_id_not_row_order():
    table = typed([("p2", "Treg", "T_cell"), ("p1", "B_cell", "B_cell")])
    out = run([box(0, 0, 10, 10)], [("p2", box(5, 0, 10, 10)), ("p1", box(0, 0, 5, 10))], table)
    assert out.loc["x0"].protein_transferred_class == "B_cell"


def test_the_transform_warps_xenium_into_the_protein_frame():
    moved = Transform(affine=np.array([[1.0, 0.0, 100.0], [0.0, 1.0, 0.0]]))
    out = run([box(0, 0, 10, 10)], [("p1", box(100, 0, 110, 10))],
              typed([("p1", "B_cell", "B_cell")]), transform=moved)
    assert out.loc["x0"].overlap_fraction == pytest.approx(1.0)


def test_row_order_does_not_change_the_result():
    xenium = [box(0, 0, 10, 10), box(20, 0, 30, 10)]
    protein = [("p1", box(0, 0, 10, 10)), ("p2", box(20, 0, 26, 10)), ("p3", box(26, 0, 30, 10))]
    table = typed([("p1", "B_cell", "B_cell"), ("p2", "Treg", "T_cell"), ("p3", "B_cell", "B_cell")])
    forward = run(xenium, protein, table)
    backward = run(xenium, protein[::-1], table.iloc[::-1])
    pd.testing.assert_frame_equal(forward, backward[forward.columns])


def test_duplicate_identifiers_are_refused():
    with pytest.raises(CallTransferError):
        overlap_edges(["x", "x"], np.array([box(0, 0, 1, 1), box(1, 1, 2, 2)], dtype=object),
                      ["p"], np.array([box(0, 0, 1, 1)], dtype=object), IDENTITY)
