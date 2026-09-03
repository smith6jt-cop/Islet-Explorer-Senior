"""The Xenium Explorer export, on geometry whose answer is arithmetic.

The bundle's whole job is a change of frame, and the failure mode of a change of frame is that
it produces a plausible number in the wrong place -- an inverted direction, a dropped half
pixel, a pyramid factor taken as a shape ratio. None of those raise. So every one of them is
planted here and checked against the value it must produce, with no export, bundle or qptiff
on the machine.
"""
import csv
import json
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

#: 10x's published Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_imagealignment.csv, all 73 bytes,
#: hexdumped 2026-09-02: CRLF between rows, no newline after the last, literal ``0,0,1``,
#: numbers of 6-9 significant digits. This is what Explorer's own export looks like.
EXEMPLAR_BYTES = (b"0.00184638,0.64417243,1433.27392\r\n"
                  b"-0.6441724,0.00184638,12184.7121\r\n"
                  b"0,0,1")


def _sig_digits(token: str) -> int:
    return len(token.replace("-", "").replace(".", "").lstrip("0").rstrip("0") or "0")


def test_alignment_csv_reproduces_explorers_own_export_byte_for_byte(tmp_path):
    """Parse the exemplar to floats, write it back with this writer: the bytes must not move."""
    M = np.array([[float(v) for v in r.split(",")] for r in EXEMPLAR_BYTES.decode().split("\r\n")])
    p = EX.write_alignment_csv(tmp_path / "a.csv", M)
    assert p.read_bytes() == EXEMPLAR_BYTES
    assert len(p.read_bytes()) == 73


def test_alignment_csv_bytes_are_shaped_like_the_exemplar_on_a_real_matrix(tmp_path):
    """Two CRLF, no trailing newline, literal last row, <= 9 significant digits, ASCII, no header.

    The previous writer produced a trailing CRLF, ``0.0,0.0,1.0`` and 17-digit floats, and
    Explorer 4.x refused the file. Every reader still agrees on the numbers.
    """
    import pandas as pd

    A = _affine(0.9976023372822801, -0.1217, 714.3188440038771, 281.7308222673026)
    M = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)
    p = EX.write_alignment_csv(tmp_path / "a.csv", M)
    b = p.read_bytes()
    assert b.count(b"\r\n") == 2 and b.count(b"\n") == 2
    assert not b.endswith((b"\n", b"\r"))
    assert b.split(b"\r\n")[-1] == b"0,0,1"
    b.decode("ascii")
    for tok in b.decode().replace("\r\n", ",").split(","):
        assert "e" not in tok.lower()
        assert _sig_digits(tok) <= EX.ALIGNMENT_SIGNIFICANT_DIGITS, tok
        float(tok)

    via_csv = np.array([[float(v) for v in r] for r in csv.reader(p.open(newline=""))])
    via_np = np.loadtxt(p, delimiter=",")
    via_pd = pd.read_csv(p, header=None).values
    via_ex = EX.read_alignment_csv(p)
    assert via_csv.shape == via_np.shape == via_pd.shape == via_ex.shape == (3, 3)
    assert np.array_equal(via_csv, via_ex) and np.array_equal(via_np, via_ex)
    assert np.array_equal(via_pd, via_ex)
    assert np.array_equal(via_ex, EX.rounded_alignment_matrix(M))
    assert np.allclose(via_ex, M, rtol=1e-8, atol=0.0)

    # The rounding is worth nothing: the far corner of a 37187 x 68119 Xenium field moves by
    # less than a thousandth of a pixel.
    fov = np.array([[0.0, 0.0, 1.0], [37187.0, 0.0, 1.0], [0.0, 68119.0, 1.0],
                    [37187.0, 68119.0, 1.0]])
    if_px = fov @ np.linalg.inv(M).T
    moved = (if_px @ via_ex.T - if_px @ M.T)[:, :2]
    assert np.abs(moved).max() < 1e-3


def test_alignment_number_format_never_uses_an_exponent_and_normalises_zero():
    assert EX.format_alignment_number(5e-05) == "0.00005"
    assert EX.format_alignment_number(-0.0) == "0"
    assert EX.format_alignment_number(0.0) == "0"
    assert EX.format_alignment_number(1433.27392) == "1433.27392"
    assert EX.format_alignment_number(1.1791626834145463) == "1.17916268"
    with pytest.raises(EX.ExplorerExportError, match="non-finite"):
        EX.format_alignment_number(float("nan"))


def test_read_alignment_csv_tolerates_lf_and_refuses_the_wrong_shape(tmp_path):
    (tmp_path / "lf.csv").write_bytes(b"1,0,5\n0,1,6\n0,0,1\n")
    assert np.array_equal(EX.read_alignment_csv(tmp_path / "lf.csv"),
                          [[1, 0, 5], [0, 1, 6], [0, 0, 1]])
    (tmp_path / "bad.csv").write_bytes(b"1,0,5\r\n0,1,6")
    with pytest.raises(EX.ExplorerExportError, match="3 rows"):
        EX.read_alignment_csv(tmp_path / "bad.csv")


def test_alignment_csv_refuses_a_matrix_that_is_not_an_alignment(tmp_path):
    with pytest.raises(EX.ExplorerExportError, match="3x3"):
        EX.write_alignment_csv(tmp_path / "a.csv", np.eye(2))
    bad = np.eye(3)
    bad[2] = [0.0, 0.1, 1.0]
    with pytest.raises(EX.ExplorerExportError, match=r"\[0,0,1\]"):
        EX.write_alignment_csv(tmp_path / "a.csv", bad)


# --------------------------------------------------------------------------- #
# The annotations
# --------------------------------------------------------------------------- #

def _square(x0, y0, side):
    import shapely
    return shapely.geometry.box(x0, y0, x0 + side, y0 + side)


def test_annotations_are_warped_backwards_into_the_xenium_frame():
    """Drawn on the IF, so they need the INVERSE transform. Forward lands 2x away."""
    t = Transform(affine=_affine(tx=700.0, ty=280.0))
    out = EX.warp_geoms_to_xenium([_square(1000.0, 1000.0, 50.0)], t)
    assert out[0].bounds == pytest.approx((300.0, 720.0, 350.0, 770.0))


def test_multipolygon_tissue_becomes_one_feature_per_part(tmp_path):
    import shapely
    multi = shapely.geometry.MultiPolygon([_square(0, 0, 10), _square(100, 100, 10)])
    rec = EX.annotations_geojson({"Tissue": [multi]}, tmp_path / "a.geojson")
    assert rec["n_features"] == 2 and rec["per_class"] == {"Tissue": 2}
    doc = json.loads((tmp_path / "a.geojson").read_text())
    assert {f["properties"]["classification"]["name"] for f in doc["features"]} == {"Tissue"}
    assert all(f["geometry"]["type"] == "Polygon" for f in doc["features"])


def test_geojson_keeps_holes_and_the_csv_says_it_dropped_them(tmp_path):
    import shapely
    ring = shapely.geometry.Polygon(
        [(0, 0), (100, 0), (100, 100), (0, 100)],
        holes=[[(40, 40), (60, 40), (60, 60), (40, 60)]])
    gj = EX.annotations_geojson({"Islet": [ring]}, tmp_path / "a.geojson")
    doc = json.loads((tmp_path / "a.geojson").read_text())
    assert len(doc["features"][0]["geometry"]["coordinates"]) == 2, "hole lost from GeoJSON"
    assert gj["n_features"] == 1

    csv_rec = EX.annotations_csv({"Islet": [ring]}, tmp_path / "a.csv")
    assert csv_rec["n_holes_dropped"] == 1
    header = next(csv.reader((tmp_path / "a.csv").open()))
    assert header == ["Selection", "X", "Y"]


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
# The threshold, and the fact that it is never asked for a frequency
# --------------------------------------------------------------------------- #

def test_two_component_threshold_lands_between_two_planted_modes():
    rng = np.random.default_rng(0)
    lo = rng.normal(2.0, 0.30, 8_000)
    hi = rng.normal(5.0, 0.30, 2_000)
    v = np.expm1(np.concatenate([lo, hi]))
    fit = EX.two_component_threshold(v)
    assert 2.0 < fit["threshold_log1p"] < 5.0
    assert fit["means_log1p"][0] == pytest.approx(2.0, abs=0.1)
    assert fit["means_log1p"][1] == pytest.approx(5.0, abs=0.1)
    # The weight is ESTIMATED and reported. It is not an input, and this asserts it landed
    # near the planted 0.2 -- which is a property of the data, not a target it was given.
    assert fit["weights"][1] == pytest.approx(0.2, abs=0.03)


def test_two_component_threshold_is_recovered_whatever_the_mixing_weight():
    """The same two modes at three different abundances give the same cut.

    This is the property the repository's no-frequency rule actually needs: the fit is not
    steered by how many cells are positive, so it cannot be used to reach a chosen fraction.
    """
    rng = np.random.default_rng(1)
    cuts = []
    for w in (0.05, 0.25, 0.50):
        n = 20_000
        k = int(n * w)
        v = np.expm1(np.concatenate([rng.normal(2.0, 0.3, n - k), rng.normal(5.0, 0.3, k)]))
        cuts.append(EX.two_component_threshold(v)["threshold_log1p"])
    assert max(cuts) - min(cuts) < 0.5


def test_two_component_threshold_refuses_too_few_values():
    with pytest.raises(EX.ExplorerExportError, match="refusing to fit"):
        EX.two_component_threshold(np.arange(10.0))


def test_decile_bands_are_ten_bands_and_nan_is_unavailable():
    v = np.arange(1000.0)
    v[:7] = np.nan
    labels, meta = EX.decile_bands(v)
    assert meta["n_available"] == 993
    assert list(labels[:7]) == ["Unavailable"] * 7
    assert set(np.unique(labels[7:])) == {f"d{i:02d}" for i in range(1, 11)}
    counts = {b: int((labels == b).sum()) for b in np.unique(labels[7:])}
    assert all(90 <= c <= 110 for c in counts.values()), counts


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


# --------------------------------------------------------------------------- #
# The morphology grid
# --------------------------------------------------------------------------- #

def _write_morphology(path, um: float, shape=(70, 50), n_levels: int = 3):
    import tifffile

    rng = np.random.default_rng(3)
    levels = [rng.integers(0, 4000, size=shape, dtype=np.uint16)]
    for L in range(1, n_levels):
        prev = levels[-1]
        levels.append(prev[: (prev.shape[0] // 2) * 2, : (prev.shape[1] // 2) * 2]
                      .reshape(prev.shape[0] // 2, 2, prev.shape[1] // 2, 2).mean(axis=(1, 3))
                      .astype(np.uint16))
    path.parent.mkdir(parents=True, exist_ok=True)
    with tifffile.TiffWriter(str(path), bigtiff=True, ome=True) as tw:
        tw.write(levels[0], subifds=n_levels - 1, tile=(16, 16), photometric="minisblack",
                 metadata={"axes": "YX", "PhysicalSizeX": um, "PhysicalSizeY": um,
                           "PhysicalSizeXUnit": "µm", "PhysicalSizeYUnit": "µm"})
        for arr in levels[1:]:
            tw.write(arr, subfiletype=1, tile=(16, 16), photometric="minisblack")
    return levels


def test_morphology_grid_reads_the_pyramid_and_refuses_a_foreign_pixel_size(tmp_path):
    good = tmp_path / "good"
    _write_morphology(good / "morphology_focus" / "morphology_focus_0000.ome.tif", 0.2125)
    g = EX.morphology_grid(good)
    assert (g["height"], g["width"], g["n_levels"]) == (70, 50, 3)
    assert g["level_shapes"] == [[70, 50], [35, 25], [17, 12]]
    assert g["um_per_px"] == 0.2125

    bad = tmp_path / "bad"
    _write_morphology(bad / "morphology_focus" / "morphology_focus_0000.ome.tif", 0.5)
    with pytest.raises(EX.ExplorerExportError, match="not Xenium"):
        EX.morphology_grid(bad)
    with pytest.raises(EX.ExplorerExportError, match="absent"):
        EX.morphology_grid(tmp_path / "nowhere")


# --------------------------------------------------------------------------- #
# The keypoints, and the similarity Explorer would fit through them
# --------------------------------------------------------------------------- #

MORPH_HW = (68119, 37187)
IF_HW = (83520, 34560)


def _cells(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    # Xenium microns, inside the 6436 field: 37187 x 68119 px at 0.2125 um
    return rng.uniform([500.0, 500.0], [7400.0, 13900.0], size=(n, 2))


def _real_matrix():
    A = _affine(0.9976023372822801, -0.1217, 714.3188440038771, 281.7308222673026)
    return EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125), A


def test_keypoints_are_exact_preimages_of_the_matrix_and_sit_on_nuclei():
    M, _ = _real_matrix()
    cells = _cells()
    rec = EX.keypoint_pairs(cells, M, morph_hw=MORPH_HW, if_hw=IF_HW, grid=(4, 5),
                            cell_ids=[f"c{i}" for i in range(len(cells))])
    kp = rec["keypoints"]
    assert 3 <= len(kp) <= 20 and kp.shape[1] == 4
    assert EX.keypoints_residual_px(kp, M) < 1e-6
    # every fixed point IS a cell centroid in morphology pixels
    assert np.allclose(kp[:, 0:2], cells[rec["index"]] / XEN)
    assert rec["cell_ids"] == [f"c{i}" for i in rec["index"]]
    assert _through_matrix(M, kp[:, 2:4]) == pytest.approx(kp[:, 0:2], abs=1e-6)
    assert (kp[:, 2] >= 0).all() and (kp[:, 2] < IF_HW[1]).all()
    assert (kp[:, 3] >= 0).all() and (kp[:, 3] < IF_HW[0]).all()


def test_keypoints_refuse_a_pair_off_either_image():
    M, _ = _real_matrix()
    with pytest.raises(EX.ExplorerExportError, match="outside"):
        EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=(1000, 1000))
    with pytest.raises(EX.ExplorerExportError, match="outside"):
        EX.keypoint_pairs(_cells(), M, morph_hw=(1000, 1000), if_hw=IF_HW)
    with pytest.raises(EX.ExplorerExportError, match="at least 3"):
        EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW, grid=(1, 2))


def test_keypoints_csv_header_bytes_and_which_column_is_which_frame(tmp_path):
    # A pure translation: Xenium px = IF px + (100, 200). fixed - alignment must be (100, 200).
    M = np.array([[1.0, 0.0, 100.0], [0.0, 1.0, 200.0], [0.0, 0.0, 1.0]])
    rec = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)
    kp = rec["keypoints"]
    assert np.allclose(kp[:, 0:2] - kp[:, 2:4], [100.0, 200.0])

    p = EX.write_keypoints_csv(tmp_path / "k.csv", kp)
    b = p.read_bytes()
    assert b.split(b"\r\n")[0] == b"fixedX,fixedY,alignmentX,alignmentY"
    assert b.count(b"\r\n") == len(kp) and not b.endswith((b"\n", b"\r"))
    b.decode("ascii")
    rows = list(csv.reader(p.open(newline="")))
    assert rows[0] == list(EX.KEYPOINTS_HEADER) and len(rows) == len(kp) + 1
    assert np.array_equal(EX.read_keypoints_csv(p), kp), "repr floats must round-trip exactly"
    with pytest.raises(EX.ExplorerExportError, match=r"\(n>=3, 4\)"):
        EX.write_keypoints_csv(tmp_path / "bad.csv", kp[:2])


def test_similarity_through_keypoints_is_exact_for_a_similarity_and_measures_a_planted_shear():
    A = _affine(0.9976, -0.1217, 714.3, 281.7)
    M = EX.alignment_matrix(A, if_um_per_px=0.25, if_origin_um=0.125)
    kp = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"]
    S = EX.similarity_from_keypoints(kp)
    assert np.allclose(S, M, atol=1e-9)
    assert EX.similarity_deviation_um(M, S, kp[:, 2:4]).max() < 1e-6
    assert EX.is_exact_similarity(S)

    # Plant 1e-3 of anisotropy: the deviation is nonzero, grows away from the centre, and is
    # exactly |(M - S) p| * 0.2125 at every point.
    A2 = A.copy()
    A2[0, 0] *= 1.001
    M2 = EX.alignment_matrix(A2, if_um_per_px=0.25, if_origin_um=0.125)
    kp2 = EX.keypoint_pairs(_cells(), M2, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"]
    S2 = EX.similarity_from_keypoints(kp2)
    dev = EX.similarity_deviation_um(M2, S2, kp2[:, 2:4])
    assert dev.max() > 0.1
    h = np.column_stack([kp2[:, 2:4], np.ones(len(kp2))])
    analytic = np.linalg.norm((h @ M2.T - h @ S2.T)[:, :2], axis=1) * XEN
    assert dev == pytest.approx(analytic)
    centre = kp2[:, 2:4].mean(axis=0)
    r = np.linalg.norm(kp2[:, 2:4] - centre, axis=1)
    assert dev[np.argmax(r)] > dev[np.argmin(r)]
    an = EX.matrix_anisotropy(M2)
    assert an["a_minus_d"] != 0.0 and an["axis_ratio"] > 1.0005 and not an["is_exact_similarity"]
    assert EX.matrix_anisotropy(M)["axis_ratio"] == pytest.approx(1.0, abs=1e-12)


def test_similarity_is_bit_exact_after_the_csv_round_trip_and_equals_umeyama(tmp_path):
    from phenocycler.integration.transform import fit_rigid_umeyama

    M, _ = _real_matrix()
    M[0, 0] *= 1.0003          # make it a real affine, like the shipped ones
    kp = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"]
    S = EX.similarity_from_keypoints(kp)
    U = fit_rigid_umeyama(kp[:, 2:4], kp[:, 0:2], allow_scale=True, allow_reflection=False)
    assert np.allclose(S[:2], U.affine, atol=1e-9)
    # The raw fit is a similarity only to floating-point noise; the canonicalised one is one
    # BEFORE any rounding -- which is the property the writer must not have to supply.
    assert U.affine[0, 0] != U.affine[1, 1] or U.affine[0, 1] != -U.affine[1, 0]
    assert S[0, 0] == S[1, 1] and S[0, 1] == -S[1, 0]
    back = EX.read_alignment_csv(EX.write_alignment_csv(tmp_path / "s.csv", S))
    assert back[0, 0] == back[1, 1] and back[0, 1] == -back[1, 0]
    assert np.linalg.det(back[:2, :2]) > 0
    assert EX.is_exact_similarity(back)


def test_similarity_refuses_a_correspondence_it_cannot_fit():
    """A reflected correspondence still yields a proper-rotation matrix; it must be refused.

    The determinant and the a==d / b==-c equalities hold for ANY output of the canonical
    rebuild, so the only guard that can fire is the fit's own residual against the field.
    """
    M, _ = _real_matrix()
    kp = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"].copy()
    kp[:, 1] = -kp[:, 1] + 70000.0                      # flip the fixed side in y
    with pytest.raises(EX.ExplorerExportError, match="misses the keypoints"):
        EX.similarity_from_keypoints(kp)
    swapped = kp[:, [2, 3, 0, 1]]                       # columns in the wrong frames
    kp_ok = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"]
    EX.similarity_from_keypoints(kp_ok)                 # the real thing still fits
    assert EX.keypoints_residual_px(kp_ok, EX.similarity_from_keypoints(kp_ok)) < \
        EX.SIMILARITY_MAX_RESIDUAL_FRACTION * np.ptp(kp_ok[:, :2], axis=0).max()
    del swapped


def test_readers_name_the_file_on_ragged_or_non_numeric_rows(tmp_path):
    (tmp_path / "ragged.csv").write_bytes(b"1,0,5\r\n0,1\r\n0,0,1")
    with pytest.raises(EX.ExplorerExportError, match="ragged.csv.*3 comma-separated"):
        EX.read_alignment_csv(tmp_path / "ragged.csv")
    (tmp_path / "text.csv").write_bytes(b"1,0,5\r\n0,x,6\r\n0,0,1")
    with pytest.raises(EX.ExplorerExportError, match="text.csv.*not a number"):
        EX.read_alignment_csv(tmp_path / "text.csv")
    (tmp_path / "k.csv").write_bytes(b"fixedX,fixedY,alignmentX,alignmentY\r\n1,2,3,4\r\n5,6,7")
    with pytest.raises(EX.ExplorerExportError, match="k.csv.*4 comma-separated"):
        EX.read_keypoints_csv(tmp_path / "k.csv")


def test_similarity_consistency_pin_catches_a_direction_error():
    M, A = _real_matrix()
    kp = EX.keypoint_pairs(_cells(), M, morph_hw=MORPH_HW, if_hw=IF_HW)["keypoints"]
    S = EX.similarity_from_keypoints(kp)
    c = EX.similarity_consistency(S, A, if_um_per_px=0.25)
    assert c["consistent"]
    assert abs(c["rotation_diff_deg"]) < 1e-2 and abs(c["scale_rel_diff"]) < 1e-5
    # The inverse direction (Xenium px -> IF px) is the same shape of number and must FAIL.
    wrong = EX.similarity_consistency(np.linalg.inv(S), A, if_um_per_px=0.25)
    assert not wrong["consistent"]


# --------------------------------------------------------------------------- #
# The image on the Xenium grid
# --------------------------------------------------------------------------- #

def _T(tx, ty):
    return np.array([[1.0, 0.0, tx], [0.0, 1.0, ty], [0.0, 0.0, 1.0]])


def test_warp_samples_the_if_at_minv_of_the_pixel_centre_exactly():
    """Output pixel (k, l) is Xenium pixel [k, k+1): its centre (k+0.5, l+0.5) maps to the IF.

    With M = T(10.5, 20.5) the centre of output (60, 80) lands exactly on IF pixel (50, 60),
    so a single planted 1.0 must come out as a single 1.0 there and nowhere else, and every
    output pixel whose source is off the image must be exactly 0.
    """
    plane = np.zeros((400, 400), np.float32)
    plane[60, 50] = 1.0
    out = EX.warp_plane_to_xenium_grid(plane, _T(10.5, 20.5), 200, 200, tile=64)
    nz = np.argwhere(out != 0)
    assert nz.tolist() == [[80, 60]] and out[80, 60] == 1.0
    assert (out[:20, :] == 0).all() and (out[:, :10] == 0).all()
    assert np.array_equal(out, EX.warp_plane_to_xenium_grid(plane, _T(10.5, 20.5), 200, 200,
                                                              tile=400))


def test_warp_identity_is_the_2x2_mean_and_a_half_pixel_shift_is_a_bit_exact_copy():
    rng = np.random.default_rng(0)
    plane = rng.uniform(0, 255, size=(64, 64)).astype(np.float32)
    # T(0.5, 0.5): output (k, l) samples IF (k, l) exactly -> a copy.
    out = EX.warp_plane_to_xenium_grid(plane, _T(0.5, 0.5), 64, 64, tile=32)
    assert np.array_equal(out, plane)
    # identity: output (k, l) samples IF (k + 0.5, l + 0.5) -> the mean of its 2x2 neighbours.
    out = EX.warp_plane_to_xenium_grid(plane, np.eye(3), 63, 63, tile=32)
    expect = 0.25 * (plane[:-1, :-1] + plane[:-1, 1:] + plane[1:, :-1] + plane[1:, 1:])
    assert np.allclose(out, expect, atol=1e-3)


def test_warp_recovers_a_blob_centroid_under_a_general_affine_to_0_05_px():
    yy, xx = np.mgrid[0:600, 0:600]
    cx, cy = 123.3, 210.7
    plane = np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 3.0 ** 2)).astype(np.float32)
    s = 0.25 / 0.2125
    t = np.radians(-0.12)
    M = np.array([[s * np.cos(t), -s * np.sin(t) + 2.6e-4, -33.7],
                  [s * np.sin(t), s * np.cos(t), 57.2],
                  [0.0, 0.0, 1.0]])
    out64 = EX.warp_plane_to_xenium_grid(plane, M, 500, 500, tile=64)
    out1000 = EX.warp_plane_to_xenium_grid(plane, M, 500, 500, tile=1000)
    assert np.allclose(out64, out1000, atol=1e-4)
    yy2, xx2 = np.mgrid[0:500, 0:500]
    w = out64.astype(np.float64)
    got = np.array([(w * xx2).sum() / w.sum(), (w * yy2).sum() / w.sum()])
    expect = (M @ np.array([cx, cy, 1.0]))[:2] - 0.5     # Xenium coordinate -> pixel index
    assert np.abs(got - expect).max() < 0.05


def test_warp_refuses_a_mismatched_out_buffer_and_leaves_off_image_tiles_zero():
    plane = np.full((32, 32), 7, np.uint8)
    with pytest.raises(EX.ExplorerExportError, match="out"):
        EX.warp_plane_to_xenium_grid(plane, np.eye(3), 32, 32, out=np.zeros((10, 10), np.uint8))
    out = EX.warp_plane_to_xenium_grid(plane, _T(1000.0, 1000.0), 64, 64, tile=16)
    assert (out == 0).all()
    # a reused, dirty buffer must come out zero where the source is off the image too
    dirty = np.full((64, 64), 9, np.uint8)
    out = EX.warp_plane_to_xenium_grid(plane, _T(1000.0, 1000.0), 64, 64, tile=16, out=dirty)
    assert out is dirty and (out == 0).all()


def test_keypoint_bounds_use_the_centre_indexed_if_frame():
    """IF pixel 0 is centred at coordinate 0, so -0.3 is on the image and W - 0.2 is off it."""
    # M = T(t): alignment = fixed - t. Put one fixed point so its alignment x is exactly -0.3.
    cells = _cells(500)
    cells[0] = [100.3 * XEN, 5000.0 * XEN]          # fixed x = 100.3 px
    M = np.array([[1.0, 0.0, 100.6], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])   # alignment x = -0.3
    rec = EX.keypoint_pairs(cells, M, morph_hw=MORPH_HW, if_hw=IF_HW, grid=(1, 3))
    assert rec["keypoints"][:, 2].min() >= -0.5
    M = np.array([[1.0, 0.0, -(IF_HW[1] - 0.2 - 100.3)], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(EX.ExplorerExportError, match="outside"):
        EX.keypoint_pairs(cells[:1].repeat(3, axis=0) + [[0, 0], [0, 1], [0, 2]], M,
                          morph_hw=MORPH_HW, if_hw=IF_HW, grid=(1, 3))


def test_xenium_pyramid_is_a_floor_halved_box_mean():
    rng = np.random.default_rng(1)
    lvl0 = rng.integers(0, 256, size=(3, 37, 23), dtype=np.uint8)
    levels = EX.xenium_pyramid(lvl0, 4)
    assert [lv.shape for lv in levels] == [(3, 37, 23), (3, 18, 11), (3, 9, 5), (3, 4, 2)]
    for prev, cur in zip(levels, levels[1:]):
        C, h1, w1 = cur.shape
        blocks = prev[:, :2 * h1, :2 * w1].astype(np.float64).reshape(C, h1, 2, w1, 2)
        expect = np.floor(blocks.mean(axis=(2, 4)) + 0.5).astype(np.uint8)
        assert np.array_equal(cur, expect)
    with pytest.raises(EX.ExplorerExportError, match="halve"):
        EX.xenium_pyramid(lvl0, 7)


def test_warped_ome_tiff_carries_the_xenium_pixel_size_and_the_morphology_dims(tmp_path):
    import re
    import tifffile

    rng = np.random.default_rng(2)
    levels = EX.xenium_pyramid(rng.integers(0, 256, size=(3, 70, 50), dtype=np.uint8), 3)
    out = tmp_path / "x_IF_on_xenium_grid.ome.tif"
    rec = EX.write_warped_ome_tiff(levels, out, ["DAPI", "INS", "GCG"], tile=16, workers=1)
    assert [(lv["height"], lv["width"]) for lv in rec["levels"]] == [(70, 50), (35, 25), (17, 12)]
    with tifffile.TiffFile(str(out)) as t:
        assert t.is_ome and t.is_bigtiff
        s = t.series[0]
        assert s.shape == (3, 70, 50) and len(s.levels) == 3
        assert [lv.shape for lv in s.levels] == [(3, 70, 50), (3, 35, 25), (3, 17, 12)]
        for L, lv in enumerate(s.levels):
            assert np.array_equal(lv.asarray(), levels[L])
        pg = t.pages[0]
        assert pg.tilewidth == 16 and pg.compression == 8
        assert len(pg.subifds) == 2
        num, den = pg.tags["XResolution"].value
        assert num / den == pytest.approx(1e4 / 0.2125, rel=1e-6)
        assert pg.tags["ResolutionUnit"].value == 3
        assert re.search(r'PhysicalSizeX="0.2125"', t.ome_metadata)
        assert re.findall(r'Channel[^>]*Name="([^"]+)"', t.ome_metadata) == ["DAPI", "INS", "GCG"]
    with pytest.raises(EX.ExplorerExportError, match="names"):
        EX.write_warped_ome_tiff(levels, tmp_path / "y.ome.tif", ["DAPI"], tile=16, workers=1)


# --------------------------------------------------------------------------- #
# The driver: reset, carry-over, and the README the maintainer reads
# --------------------------------------------------------------------------- #

def test_reset_spares_both_ome_tiffs_and_deletes_every_csv(tmp_path):
    """`--skip-image` / `--skip-warped` must not delete the multi-GB files they exist to keep,
    and no CSV may survive a reset -- a stale matrix beside a new image is the trap."""
    import run_explorer_export as RX

    dest = tmp_path / "donor_6436_0059907"
    (dest / "images").mkdir(parents=True)
    (dest / "cell_groups").mkdir()
    (dest / "images" / "0059907_IF.ome.tif").write_bytes(b"pretend-pyramid")
    (dest / "images" / "0059907_IF_on_xenium_grid.ome.tif").write_bytes(b"pretend-warped")
    for name in ("0059907_IF_imagealignment.csv", "0059907_IF_similarity_imagealignment.csv",
                 "0059907_IF_keypoints.csv"):
        (dest / "images" / name).write_text("1,0,0\n0,1,0\n0,0,1")
    (dest / "cell_groups" / "stale.csv").write_text("cell_id,group,color\n")

    kept = RX._reset(dest, keep_images=True)
    assert kept == ["0059907_IF.ome.tif", "0059907_IF_on_xenium_grid.ome.tif"]
    assert (dest / "images" / "0059907_IF.ome.tif").read_bytes() == b"pretend-pyramid"
    assert (dest / "images" / "0059907_IF_on_xenium_grid.ome.tif").read_bytes() == b"pretend-warped"
    assert sorted(p.name for p in (dest / "images").iterdir()) == kept
    assert not (dest / "cell_groups").exists()
    assert not (dest.parent / ".donor_6436_0059907.keep").exists()


def test_reset_refuses_a_stash_left_by_an_interrupted_run_and_never_resurrects_it(tmp_path):
    import run_explorer_export as RX

    dest = tmp_path / "donor_6436_0059907"
    (dest / "images").mkdir(parents=True)
    (dest / "images" / "0059907_IF.ome.tif").write_bytes(b"current")
    stash = tmp_path / ".donor_6436_0059907.keep"
    stash.mkdir()
    (stash / "0059907_IF_OLDNAME.ome.tif").write_bytes(b"orphan")
    with pytest.raises(RX.ExportError, match="interrupted run"):
        RX._reset(dest, keep_images=False)
    # nothing was touched by the refusal
    assert (dest / "images" / "0059907_IF.ome.tif").read_bytes() == b"current"
    assert (stash / "0059907_IF_OLDNAME.ome.tif").exists()


def test_a_kept_image_without_a_pixel_size_is_refused_not_described_as_nan(tmp_path):
    import tifffile

    import run_explorer_export as RX

    p = tmp_path / "0059907_IF.ome.tif"
    with tifffile.TiffWriter(str(p), ome=True) as tw:
        tw.write(np.zeros((2, 16, 16), np.uint8), tile=(16, 16), photometric="minisblack",
                 metadata={"axes": "CYX", "Channel": {"Name": ["DAPI", "INS"]}})
    with pytest.raises(RX.ExportError, match="PhysicalSizeX"):
        RX._describe_ome_tiff(p)


def test_a_kept_image_is_described_from_its_own_bytes(tmp_path):
    import run_explorer_export as RX

    rng = np.random.default_rng(5)
    levels = EX.xenium_pyramid(rng.integers(0, 256, size=(3, 70, 50), dtype=np.uint8), 3)
    out = tmp_path / "0059907_IF_on_xenium_grid.ome.tif"
    EX.write_warped_ome_tiff(levels, out, ["DAPI", "INS", "GCG"], tile=16, workers=1)
    rec = RX._describe_ome_tiff(out)
    assert rec["skipped"] is True and rec["preserved_from_a_previous_run"] == [out.name]
    assert rec["channels"] == ["DAPI", "INS", "GCG"] and rec["dtype"] == "uint8"
    assert [(lv["height"], lv["width"], lv["factor"]) for lv in rec["levels"]] == \
        [(70, 50, 1), (35, 25, 2), (17, 12, 4)]
    assert rec["levels"][1]["um_per_px"] == pytest.approx(0.425)
    assert rec["bytes"] == out.stat().st_size and rec["tile"] == 16
    assert RX._describe_ome_tiff(tmp_path / "missing.ome.tif") is None


def _manifest_fixture() -> dict:
    return {
        "generated_utc": "2026-09-02T00:00:00+00:00", "donor_id": "6436",
        "xenium_serial": "0059907",
        "xenium_bundle": "/x/output-XETG00298__0059907__Panc__20250328__202338",
        "transform": {"translation_um": [714.3, 281.7]},
        "morphology": {"width": 37187, "height": 68119, "if_over_morphology_pixel_ratio": 1.14},
        "alignment": {"file": "images/0059907_IF_imagealignment.csv",
                      "applies_to": "images/0059907_IF.ome.tif",
                      "footprint_xenium_px": {"x_range": [-3591.0, 37386.0],
                                              "y_range": [-1336.0, 97257.0]},
                      "if_um_per_px": 0.25, "round_trip_n_cells": 695113,
                      "round_trip_max_um": 4e-12, "round_trip_as_written_max_um": 6e-5,
                      "unresolved_convention_um": 0.125},
        "keypoints": {"file": "images/0059907_IF_keypoints.csv", "n": 20,
                      "residual_px_through_shipped_matrix": 1e-10},
        "similarity_diagnostic": {"file": "images/0059907_IF_similarity_imagealignment.csv",
                                  "deviation_from_affine_um": {
                                      "over_cells": {"p50": 0.42, "max": 0.87},
                                      "at_fov_corners_max": 1.0}},
        "image_on_xenium_grid": {"file": "images/0059907_IF_on_xenium_grid.ome.tif"},
        "positivity": {"counts": {"Negative": 10, "INS_positive": 1}, "measure": "mean",
                       "threshold_counts": 20.0, "sampled_um_per_px": 0.5,
                       "channels": {"INS": "Opal 690", "GCG": "Opal 520"},
                       "n_without_polygon": 0, "n_cells": 11},
    }


def test_readme_names_every_file_and_every_explorer_choice_and_forbids_the_wrong_ones():
    import run_explorer_export as RX

    text = RX._readme(_manifest_fixture())
    for needed in ("Yes, use an existing alignment file", "No, image is not aligned",
                   "Import a New Alignment File", "Yes, skip alignment", "Help → About",
                   "__Panc__", "__pLN__", "Remove any IF image layer",
                   "0059907_IF.ome.tif", "0059907_IF_imagealignment.csv",
                   "0059907_IF_similarity_imagealignment.csv", "0059907_IF_keypoints.csv",
                   "0059907_IF_on_xenium_grid.ome.tif", "37187 × 68119", "Never apply a CSV",
                   "never choose \"Yes, skip", "Report back", "hormone_positivity.csv"):
        assert needed in text, needed
    assert "not keypoints" not in text and "annotations/" not in text


# --------------------------------------------------------------------------- #
# The bundle QC: the identity sampler and its sweep
# --------------------------------------------------------------------------- #

def test_xenium_grid_sampler_uses_floor_not_rint():
    import explorer_figures as FIG

    plane = np.arange(100, dtype=np.float32).reshape(10, 10)
    # a point at um = 0.2125 * (k + 0.9) lies inside pixel k; rint would give k + 1.
    xy = np.array([[XEN * 3.9, XEN * 6.9]])
    assert FIG._sample_on_xenium_grid(plane, xy, 0, 0).tolist() == [plane[6, 3]]
    assert FIG._sample_on_xenium_grid(plane, xy, 1, 2).tolist() == [plane[4, 2]]
    assert FIG._sample_on_xenium_grid(plane, np.array([[XEN * 50, XEN * 50]]), 0, 0).size == 0


def test_identity_sweep_peaks_at_zero_on_synthetic_blobs():
    import explorer_figures as FIG

    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:200, 0:200]
    plane = np.zeros((200, 200), np.float32)
    ks = rng.integers(10, 190, size=(60, 2))
    for k, l in ks:
        plane += np.exp(-((xx - k) ** 2 + (yy - l) ** 2) / (2 * 2.0 ** 2)).astype(np.float32)
    # nuclei recorded in microns at their pixel CENTRES (corner-origin grid)
    xy_um = (ks + 0.5) * XEN
    offsets = (0.0, 0.5, 1.0, 2.0, 4.0)
    sw = FIG.sweep_on_xenium_grid(plane, xy_um, 0, 0, offsets)
    assert int(np.argmax(sw)) == 0
    assert sw[0] > sw[1] > sw[2]                  # falls off the nucleus
    assert (sw[2:] < 0.2 * sw[0]).all()          # and stays on the floor past 1 um
