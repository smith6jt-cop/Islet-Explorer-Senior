"""The Xenium Explorer export, on geometry whose answer is arithmetic.

The bundle's whole job is a change of frame, and the failure mode of a change of frame is that
it produces a plausible number in the wrong place -- an inverted direction, a dropped half
pixel, a pyramid factor taken as a shape ratio. None of those raise. So every one of them is
planted here and checked against the value it must produce, with no export, bundle or qptiff
on the machine.
"""
import csv
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import explorer as EX                                          # noqa: E402
from phenocycler.integration.transform import Transform        # noqa: E402

XEN = EX.XENIUM_UM_PER_PX          # 0.2125


def _affine(scale=1.0, rot_deg=0.0, tx=0.0, ty=0.0):
    t = np.radians(rot_deg)
    R = scale * np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return np.column_stack([R, [tx, ty]])


def _through_matrix(M, if_px):
    """IF pixel coordinates pushed through the alignment matrix, in Xenium pixels."""
    if_px = np.asarray(if_px, dtype=float).reshape(-1, 2)
    return (np.column_stack([if_px, np.ones(len(if_px))]) @ np.asarray(M).T)[:, :2]


def _one(M, if_px):
    """The single point form, flattened -- ``pytest.approx`` will not take a nested list."""
    return list(_through_matrix(M, [if_px])[0])


# --------------------------------------------------------------------------- #
# The direction, which no amount of care in the caller can recover
# --------------------------------------------------------------------------- #

def test_matrix_maps_if_pixels_to_xenium_pixels_not_the_reverse():
    """A planted translation must come out INVERTED and in Xenium pixels.

    The lane's affine maps Xenium microns to IF microns. Explorer wants the other direction.
    A matrix built the forward way is the same size and shape and lands 2x the offset away.
    """
    A = _affine(tx=700.0, ty=280.0)
    M = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)

    # The Xenium origin (0, 0) um sits at IF micron (700, 280), i.e. IF pixel
    # ((700 - 0.125)/0.25, (280 - 0.125)/0.25). Pushing that pixel through M must land on
    # Xenium pixel (0, 0).
    if_px = [(700.0 - 0.125) / 0.25, (280.0 - 0.125) / 0.25]
    assert _one(M, if_px) == pytest.approx([0.0, 0.0], abs=1e-9)

    # And the forward-built matrix would put it at -2x the offset in Xenium pixels.
    wrong = _one(EX.alignment_matrix(_affine(tx=-700.0, ty=-280.0),
                                     if_um_per_px=0.25, if_origin_um=0.125), if_px)
    assert abs(wrong[0]) > 1000.0


def test_matrix_target_unit_is_xenium_pixels_at_0_2125_um():
    A = _affine()
    M = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)
    # IF pixel 0 is IF micron 0.125 = Xenium micron 0.125 = Xenium pixel 0.125/0.2125.
    assert _one(M, [0.0, 0.0]) == pytest.approx([0.125 / XEN, 0.125 / XEN])
    # One IF pixel is 0.25 um is 0.25/0.2125 Xenium pixels.
    step = _one(M, [1.0, 0.0])[0] - _one(M, [0.0, 0.0])[0]
    assert step == pytest.approx(0.25 / XEN)


@pytest.mark.parametrize("scale,rot,tx,ty", [(1.0, 0.0, 0.0, 0.0),
                                             (0.9976, -0.1217, 714.3188, 281.7308),
                                             (1.0031, 2.5, -120.0, 3000.0)])
def test_round_trip_is_exact_for_any_affine(scale, rot, tx, ty):
    A = _affine(scale, rot, tx, ty)
    M = EX.alignment_matrix(A, if_um_per_px=1.0, if_origin_um=0.5)
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 15_000, size=(500, 2))
    err = EX.matrix_round_trip_um(M, A, if_um_per_px=1.0, if_origin_um=0.5, xy_xenium_um=xy)
    assert err < 1e-8


def test_dropping_the_pixel_origin_displaces_by_exactly_that_origin():
    """The half-pixel is worth its own value in every axis, and it is not noise."""
    A = _affine()
    right = EX.alignment_matrix(A, if_um_per_px=1.0, if_origin_um=0.5)
    wrong = EX.alignment_matrix(A, if_um_per_px=1.0, if_origin_um=0.0)
    delta_px = np.array(_one(right, [10.0, 10.0])) - np.array(_one(wrong, [10.0, 10.0]))
    assert list(delta_px * XEN) == pytest.approx([0.5, 0.5])


def test_pyramid_level_changes_the_matrix_by_exactly_the_factor():
    """Shipping level 2 instead of level 0 must rescale the matrix, not shift it."""
    A = _affine(tx=700.0, ty=280.0)
    lvl0 = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)
    lvl2 = EX.alignment_matrix(A, if_um_per_px=1.0, if_origin_um=0.5)
    # level-2 pixel k covers level-0 pixels 4k..4k+3, so the same tissue point is at
    # level-0 pixel (4k + 1.5); both must land on the same Xenium pixel.
    k = np.array([[123.0, 456.0]])
    assert _through_matrix(lvl2, k) == pytest.approx(_through_matrix(lvl0, 4 * k + 1.5))


def test_alignment_matrix_refuses_singular_and_malformed():
    with pytest.raises(EX.ExplorerExportError, match="2x3"):
        EX.alignment_matrix(np.eye(3), if_um_per_px=1.0, if_origin_um=0.5)
    with pytest.raises(EX.ExplorerExportError, match="singular"):
        EX.alignment_matrix(np.array([[1.0, 1.0, 0.0], [2.0, 2.0, 0.0]]),
                            if_um_per_px=1.0, if_origin_um=0.5)
    with pytest.raises(EX.ExplorerExportError, match="non-finite"):
        EX.alignment_matrix(np.array([[1.0, 0.0, np.nan], [0.0, 1.0, 0.0]]),
                            if_um_per_px=1.0, if_origin_um=0.5)
    with pytest.raises(EX.ExplorerExportError, match="positive"):
        EX.alignment_matrix(_affine(), if_um_per_px=0.0, if_origin_um=0.5)


# --------------------------------------------------------------------------- #
# The file 10x specifies
# --------------------------------------------------------------------------- #

def test_alignment_csv_is_three_rows_no_header_and_survives_the_write(tmp_path):
    A = _affine(0.9976023372822801, -0.1217, 714.3188440038771, 281.7308222673026)
    M = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)
    p = EX.write_alignment_csv(tmp_path / "a.csv", M)

    rows = list(csv.reader(p.open()))
    assert len(rows) == 3 and all(len(r) == 3 for r in rows)
    assert [float(v) for v in rows[2]] == [0.0, 0.0, 1.0]
    # Read back the way spatialdata-io reads it: pd.read_csv(..., header=None).values
    back = np.array([[float(v) for v in r] for r in rows])
    assert np.array_equal(back, M), "the written text must not round the matrix at all"


def test_alignment_csv_refuses_a_matrix_that_is_not_an_alignment(tmp_path):
    with pytest.raises(EX.ExplorerExportError, match="3x3"):
        EX.write_alignment_csv(tmp_path / "a.csv", np.eye(2))
    bad = np.eye(3)
    bad[2] = [0.0, 0.1, 1.0]
    with pytest.raises(EX.ExplorerExportError, match=r"\[0,0,1\]"):
        EX.write_alignment_csv(tmp_path / "a.csv", bad)


# --------------------------------------------------------------------------- #
# The cell groups
# --------------------------------------------------------------------------- #

def test_cell_groups_header_and_colour_are_exactly_what_explorer_reads(tmp_path):
    rec = EX.write_cell_groups(tmp_path / "g.csv", ["a-1", "b-1", "c-1"],
                               ["hi", "lo", "hi"], {"hi": "#ff0055", "lo": "#adb5bd"})
    rows = list(csv.reader((tmp_path / "g.csv").open()))
    assert rows[0] == ["cell_id", "group", "color"]
    assert rows[1] == ["a-1", "hi", "#ff0055"]
    assert rec["groups"] == {"hi": 2, "lo": 1}


def test_cell_groups_refuse_a_group_with_no_colour(tmp_path):
    with pytest.raises(EX.ExplorerExportError, match="no colour"):
        EX.write_cell_groups(tmp_path / "g.csv", ["a-1"], ["hi"], {"lo": "#adb5bd"})


def test_cell_groups_refuse_a_length_mismatch(tmp_path):
    with pytest.raises(EX.ExplorerExportError, match="against"):
        EX.write_cell_groups(tmp_path / "g.csv", ["a-1", "b-1"], ["hi"], {"hi": "#ff0055"})


# --------------------------------------------------------------------------- #
# The pyramid, and the gate the labels come from
# --------------------------------------------------------------------------- #

def test_integer_factors_refuse_a_pyramid_that_is_not_exactly_halved():
    ok = [(68119, 37187), (34059, 18593), (17029, 9296)]
    assert EX.integer_factors(ok) == [1, 2, 4]
    with pytest.raises(ValueError):
        EX.integer_factors([(1000, 1000), (300, 300)])


def test_tile_stream_pads_the_last_partial_tile_with_zeros():
    a = np.ones((1, 1500, 1200), dtype=np.uint8)
    tiles = list(EX._tile_stream(a, [0], 1500, 1200, a.dtype, tile=1024))
    assert len(tiles) == 4
    assert all(t.shape == (1024, 1024) for t in tiles)
    assert tiles[0].sum() == 1024 * 1024                     # wholly inside
    assert tiles[3].sum() == (1500 - 1024) * (1200 - 1024)   # the corner offcut


def test_reset_clears_the_bundle_but_skip_image_spares_the_ome_tiff(tmp_path):
    """`--skip-image` must not delete the 2 GB file it exists to avoid rewriting."""
    import run_explorer_export as RX

    dest = tmp_path / "donor_6436_0059907"
    (dest / "images").mkdir(parents=True)
    (dest / "cell_groups").mkdir()
    (dest / "images" / "0059907_IF.ome.tif").write_bytes(b"pretend-pyramid")
    (dest / "images" / "0059907_IF_imagealignment.csv").write_text("1,0,0\n0,1,0\n0,0,1\n")
    (dest / "cell_groups" / "stale.csv").write_text("cell_id,group,color\n")

    kept = RX._reset(dest, keep_images=True)
    assert kept == ["0059907_IF.ome.tif"]
    assert (dest / "images" / "0059907_IF.ome.tif").read_bytes() == b"pretend-pyramid"
    assert not (dest / "cell_groups" / "stale.csv").exists(), "a stale group file survived"
    assert not (dest / "images" / "0059907_IF_imagealignment.csv").exists()

    assert RX._reset(dest, keep_images=False) == []
    assert not (dest / "images" / "0059907_IF.ome.tif").exists()
    assert (dest / "images").is_dir()


def test_reset_on_a_directory_that_does_not_exist_yet(tmp_path):
    import run_explorer_export as RX

    dest = tmp_path / "donor_9999_0000000"
    assert RX._reset(dest, keep_images=True) == []
    assert (dest / "images").is_dir()
