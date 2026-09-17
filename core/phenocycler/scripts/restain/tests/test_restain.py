"""The restain lane, on geometry whose answer is known.

Nothing here needs an export or a bundle: the point is that the machinery is checked against
synthetic shapes where the right answer is arithmetic, so it still runs when the data is not
on this machine.
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import geo                                                     # noqa: E402
import validate as VAL                                         # noqa: E402
from phenocycler.integration import rasterize as RZ            # noqa: E402


def _square_geojson(tmp, x0, y0, side_px, n_islets=0, islet_px=10):
    """A Tissue square plus n small Islet squares, in qptiff PIXELS."""
    feats = [{"type": "Feature", "properties": {"objectType": "annotation",
                                                "classification": {"name": "Tissue"}},
              "geometry": {"type": "Polygon", "coordinates": [[
                  [x0, y0], [x0 + side_px, y0], [x0 + side_px, y0 + side_px],
                  [x0, y0 + side_px], [x0, y0]]]}}]
    for i in range(n_islets):
        cx, cy = x0 + 100 + 200 * i, y0 + 100
        h = islet_px / 2
        feats.append({"type": "Feature",
                      "properties": {"objectType": "annotation",
                                     "classification": {"name": "Islet"}},
                      "geometry": {"type": "Polygon", "coordinates": [[
                          [cx - h, cy - h], [cx + h, cy - h], [cx + h, cy + h],
                          [cx - h, cy + h], [cx - h, cy - h]]]}})
    p = tmp / "0059999_Scan1.qptiff - resolution #1.geojson"
    p.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    return p


def test_annotations_arrive_in_microns_not_pixels(tmp_path):
    """The conversion is a constant, and getting it wrong does not error -- it silently puts
    a section's annotations a whole tissue away from its own cells."""
    p = _square_geojson(tmp_path, 1000, 2000, 4000)
    ann = geo.load_annotations(p)
    tis = geo.tissue_polygon(ann)
    x0, y0, x1, y1 = tis.bounds
    assert (x0, y0) == pytest.approx((1000 * geo.PIXEL_UM, 2000 * geo.PIXEL_UM))
    assert (x1 - x0) == pytest.approx(4000 * geo.PIXEL_UM)


def test_tissue_presence_requires_a_nonempty_tissue_class():
    assert geo.has_tissue({geo.TISSUE: [object()], geo.ISLET: [object()]})
    assert not geo.has_tissue({geo.ISLET: [object()]})
    assert not geo.has_tissue({geo.TISSUE: []})


def test_the_islet_size_filter_keeps_islets_and_drops_fragments(tmp_path):
    """95% of the classifier's `Islet` objects are sub-cell fragments -- median 23-27 um^2,
    about one cell across.  Without the filter a centroid set is mostly noise, so the filter
    is part of the contract, not a convenience."""
    # 3 objects at 40 px (10 um across -> below the bar) and 3 at 400 px (100 um -> above)
    small = _square_geojson(tmp_path, 0, 0, 8000, n_islets=3, islet_px=40)
    ann = geo.load_annotations(small)
    assert len(ann[geo.ISLET]) == 3
    assert geo.islets(ann) == []
    big = _square_geojson(tmp_path, 0, 0, 8000, n_islets=3, islet_px=400)
    ann = geo.load_annotations(big)
    kept = geo.islets(ann)
    assert len(kept) == 3
    assert geo.centroids(kept).shape == (3, 2)


def test_a_wrong_pixel_size_is_refused_rather_than_registered(tmp_path):
    """The failure this guards is silent: a wrong scale produces a registration that cannot
    converge and an overlay that merely looks misaligned."""
    p = _square_geojson(tmp_path, 0, 0, 4000)             # 1000 x 1000 um at 0.25 um/px
    tis = geo.tissue_polygon(geo.load_annotations(p))
    matching = np.random.default_rng(0).uniform([0, 0], [1000, 1000], size=(500, 2))
    geo.assert_scale_is_plausible(tis, matching)          # must not raise
    with pytest.raises(geo.RestainGeoError, match="outside"):
        geo.assert_scale_is_plausible(tis, matching * 2.0)


def test_rasterise_puts_a_polygon_where_the_grid_says_it_is(tmp_path):
    """A grid silently disagreeing with a transform is how a plausible, completely wrong
    overlay gets produced, so the mask is checked against the arithmetic."""
    p = _square_geojson(tmp_path, 0, 0, 4000)             # 0..1000 um square
    tis = geo.tissue_polygon(geo.load_annotations(p))
    grid = RZ.GridSpec(0.0, 0.0, 10.0, 150, 150)
    mask = geo.rasterize([tis], grid)
    assert mask.shape == (150, 150)
    assert mask[:100, :100].all(), "the square must fill its own 1000 um"
    assert not mask[120:, 120:].any(), "and nothing outside it"


def test_the_gene_set_score_is_not_carried_by_one_bright_gene():
    """The same correction the Xenium lane applies, for the same measured reason: without
    per-gene SD scaling one gene carries the set and the score stops being about the set."""
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    n = 400
    X = np.zeros((n, 3))
    X[:, 0] = rng.poisson(1.0, n)                 # ordinary
    X[:, 1] = rng.poisson(1.0, n)                 # ordinary
    X[:, 2] = rng.poisson(1.0, n) * 100.0         # one very bright gene
    idx = {"a": 0, "b": 1, "loud": 2}
    s = VAL.gene_set_score(sp.csr_matrix(X), ["a", "b", "loud"], idx)
    only_loud = VAL.gene_set_score(sp.csr_matrix(X), ["loud"], idx)
    # If `loud` carried the set, the 3-gene score would track it almost perfectly.
    assert np.corrcoef(s, only_loud)[0, 1] < 0.95


def test_enrichment_finds_a_planted_signal_and_refuses_a_flat_one():
    """Both directions, because a statistic that cannot come back negative is not a test."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 1000, size=(4000, 2))
    landmarks = np.array([[250.0, 250.0], [750.0, 750.0]])
    d = np.minimum(*[np.hypot(*(xy - p).T) for p in landmarks])
    planted = np.where(d < 25.0, 5.0, 0.1)
    hit = VAL.enrichment(planted, xy, landmarks, draws=10)
    assert hit["ratio"] > 5 and hit["clears_null"]
    flat = rng.normal(1.0, 0.1, size=len(xy))
    miss = VAL.enrichment(flat, xy, landmarks, draws=10)
    assert 0.8 < miss["ratio"] < 1.25
    assert not miss["clears_null"]


# --------------------------------------------------------------------- the hormone channels
import hormone as HZ                                             # noqa: E402


def test_a_channel_recorded_failed_is_refused_not_read(tmp_path):
    """`Opal 570` is recorded FAILED. A caller who does not know that must not be able to
    read it by asking for it -- which is why the refusal lives in the accessor rather than in
    a comment."""
    row = pd.Series({"image_file": "x.qptiff", "channel_opal570": "FAILED",
                     "channel_opal690": "INS", "channel_sampleaf": "AF"})
    with pytest.raises(HZ.HormoneError, match="FAILED"):
        HZ.channel_index(tmp_path / "x.qptiff", "FAILED", row)
    with pytest.raises(HZ.HormoneError, match="must not be read"):
        HZ.channel_index(tmp_path / "x.qptiff", "AF", row)


def test_the_section_map_refuses_a_blank_rather_than_guessing(tmp_path):
    """Both facts in this file are in no image. A blank is the one thing that must fail."""
    good = ("image_file,donor_id,roi,xenium_serial,pixel_size_um\n"
            "a.qptiff,6442,panc,0059947,0.25\n")
    p = tmp_path / "ok.csv"; p.write_text(good)
    assert len(HZ.load_section_map(p)) == 1
    bad = tmp_path / "bad.csv"
    bad.write_text("image_file,donor_id,roi,xenium_serial,pixel_size_um\na.qptiff,,panc,0059947,0.25\n")
    with pytest.raises(HZ.HormoneError, match="blank"):
        HZ.load_section_map(bad)


def test_sampling_returns_nan_off_the_image_not_zero():
    """A cell that fell off the edge is unavailable, not dark. A 0 would sort it to the
    bottom of every rank and quietly bias the AUROC."""
    img = np.full((50, 50), 7.0, dtype=np.float32)
    vals = HZ.sample_cells(img, 1.0, np.array([[25.0, 25.0], [-100.0, -100.0], [999.0, 999.0]]),
                           radius_um=2.0)
    assert vals[0] == pytest.approx(7.0)
    assert np.isnan(vals[1]) and np.isnan(vals[2])


def test_auroc_is_rank_based_and_has_no_threshold():
    perfect = HZ.auroc(np.array([1.0, 2, 3, 4]), np.array([False, False, True, True]))
    assert perfect == pytest.approx(1.0)
    inverted = HZ.auroc(np.array([4.0, 3, 2, 1]), np.array([False, False, True, True]))
    assert inverted == pytest.approx(0.0)
    # a channel that is flat scores 0.5 rather than being rewarded for its tie structure
    flat = HZ.auroc(np.ones(100), np.arange(100) < 50)
    assert flat == pytest.approx(0.5)


def test_the_permutation_null_keeps_both_distributions_and_destroys_only_the_pairing():
    rng = np.random.default_rng(0)
    score = rng.normal(size=2000)
    pos = score > 0.5                                   # perfectly paired
    real = HZ.auroc_with_null(score, pos, draws=50)
    assert real["auroc"] > 0.9 and real["clears_null"]
    assert real["null_mean"] == pytest.approx(0.5, abs=0.05)
    shuffled = HZ.auroc_with_null(score, rng.permutation(pos), draws=50)
    assert not shuffled["clears_null"]


def test_agreement_by_scale_separates_an_offset_from_structure_level_agreement():
    """The diagnostic exists to tell those two apart, so it is checked on both.

    A fixed OFFSET plateaus once the radius covers it; a purely large-scale correspondence
    keeps climbing past the size of the structures.
    """
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 400, size=(1500, 2))
    block = ((xy[:, 0] // 100) + (xy[:, 1] // 100)) % 2      # 100 um chequerboard
    rna = block + rng.normal(0, 2.0, len(xy))                # buried in per-cell noise
    prot = block + rng.normal(0, 2.0, len(xy))
    prof = HZ.agreement_by_scale(xy, rna, prot, radii_um=(0, 10, 50))
    assert prof[0]["spearman"] < prof[-1]["spearman"], (
        "averaging must recover a structure-level correspondence buried in per-cell noise")
    assert prof[0]["radius_um"] == 0 and prof[0]["mean_neighbours"] == 1.0


def _fake_qptiff(path, names, *, base=(64, 96), n_levels=3, px_per_cm=39924.958):
    """A minimal pyramid qptiff: `n_levels` levels, one page per channel, names in the XML.

    `px_per_cm` defaults to what donor 6442's file records, which works out to 0.2504699
    um/px -- deliberately NOT the 0.25 the same file's XML and the section registry record,
    so a test can tell which of the two a caller used.
    """
    import tifffile
    with tifffile.TiffWriter(str(path)) as tw:
        for lvl in range(n_levels):
            h, w = base[0] >> lvl, base[1] >> lvl
            for ci, nm in enumerate(names):
                tw.write(
                    np.full((h, w), ci + 1, dtype=np.uint8),
                    description=f"<Root><Name>{nm}</Name></Root>",
                    photometric="minisblack",
                    subfiletype=0 if lvl == 0 else 1,
                    resolution=(px_per_cm, px_per_cm),
                    resolutionunit="CENTIMETER",
                    contiguous=False,
                )
    return path


def test_the_annotation_scale_and_the_recorded_scale_are_the_same_number():
    """The failure this catches is not a scale error, so it does not look like one.

    `geo.PIXEL_UM` puts the annotations into the microns the registration is fitted in;
    `pixel_size_um` is what turns those microns back into pixels of the image being sampled.
    If they differ the round trip is not the identity, and the residual is proportional to the
    coordinate -- ~0 at the image origin and growing across the slide. A fitted affine cannot
    absorb it because it happens after the transform.
    """
    sections = HERE.parents[2] / "data" / "integration" / "restain_sections.csv"
    df = HZ.load_section_map(sections)
    assert len(df) >= 1
    recorded = sorted({float(v) for v in df["pixel_size_um"]})
    assert recorded == [geo.PIXEL_UM], (
        f"restain_sections.csv records pixel_size_um {recorded} but annotations are scaled by "
        f"geo.PIXEL_UM={geo.PIXEL_UM}; these are the two ends of one round trip.")


def test_a_recorded_pixel_size_overrides_the_resolution_tag(tmp_path):
    """The tag is not wrong, it just disagrees, and the caller knows which frame it is in."""
    p = _fake_qptiff(tmp_path / "s.qptiff", ["DAPI", "Opal 690"])

    _, um_tag = RZ.read_qptiff_channel(p, 0, target_um=1.0)
    assert um_tag == pytest.approx(1.0018796, rel=1e-6)          # 4 x 0.2504699

    _, um_rec = RZ.read_qptiff_channel(p, 0, target_um=1.0, base_um_per_px=geo.PIXEL_UM)
    assert um_rec == pytest.approx(1.0, rel=1e-12)               # 4 x 0.25, exactly

    # The size of the thing: 0.188% of a coordinate, not 0.188% of a length.
    drift_um = 15_402.0 * (1.0 - um_rec / um_tag)
    assert drift_um == pytest.approx(28.9, abs=0.5)              # ~3 cell diameters

    with pytest.raises(ValueError, match="must be positive"):
        RZ.read_qptiff_channel(p, 0, target_um=1.0, base_um_per_px=0.0)


def test_the_af_channel_is_found_by_its_recorded_role_not_by_page_order(tmp_path, monkeypatch):
    """The previous rule zipped the file's page order against a hardcoded column order. It was
    right for the five sections recorded so far, and the first section acquired in a different
    order would have regressed the wrong covariate out -- with no error and no test.

    The channel NAMES are stubbed rather than written into a fixture qptiff, because what
    changed is the role -> column -> name lookup, not the file parsing.
    """
    names = ["Sample AF", "Opal 520", "DAPI", "Opal 690", "Opal 570"]   # NOT the recorded order
    monkeypatch.setattr(HZ, "channel_names", lambda _p: names)
    p = tmp_path / "shuffled.qptiff"
    row = pd.Series({"image_file": "shuffled.qptiff", "channel_dapi": "NUCLEUS",
                     "channel_opal520": "GCG", "channel_opal570": "FAILED",
                     "channel_opal690": "INS", "channel_sampleaf": "AF"})

    assert names[HZ.channel_index_for_role(p, row, "AF")] == "Sample AF"
    assert names[HZ.channel_index_for_role(p, row, "FAILED")] == "Opal 570"
    assert HZ.channel_index_for_role(p, row, "SST") is None
    # the marker accessor agrees with itself on the same shuffled file
    assert names[HZ.channel_index(p, "INS", row)] == "Opal 690"
    assert names[HZ.channel_index(p, "GCG", row)] == "Opal 520"
    # and the positional rule this replaced would have picked the DAPI page for AF
    assert names[0] != "DAPI"


def test_an_ambiguous_role_is_refused_rather_than_picked(tmp_path, monkeypatch):
    monkeypatch.setattr(HZ, "channel_names", lambda _p: ["DAPI", "Opal 690", "Sample AF"])
    row = pd.Series({"image_file": "amb.qptiff", "channel_opal690": "AF",
                     "channel_sampleaf": "AF"})
    with pytest.raises(HZ.HormoneError, match="ambiguous role is refused"):
        HZ.channel_index_for_role(tmp_path / "amb.qptiff", row, "AF")


# --------------------------------------------------------------------------- #
# Windowed crops. Both failure modes here are silent: a clamped read is the WRONG
# SIZE, and a chunk the file never stored decodes to a legal-looking black.
# --------------------------------------------------------------------------- #

import tiles as TL                                                # noqa: E402


def _tiled_tif(path, h=64, w=96, tile=32, sparse=()):
    """A tiled uint8 TIFF whose pixel value encodes its position, with optional empty tiles."""
    import tifffile
    ny, nx = -(-h // tile), -(-w // tile)

    def chunks():
        for i in range(ny):
            for j in range(nx):
                if (i, j) in sparse:
                    yield None                       # written as a chunk of zero byte count
                else:
                    yield np.full((tile, tile), 1 + i * nx + j, dtype=np.uint8)

    tifffile.imwrite(str(path), chunks(), shape=(h, w), dtype=np.uint8,
                     tile=(tile, tile), photometric="minisblack")
    return path


def test_a_window_running_off_the_page_is_placed_not_stretched(tmp_path):
    """`read_window` returns the INTERSECTION, so it comes back smaller than requested. The
    grid-shaped buffer must place that block at its true offset; assigning it at [0, 0] would
    shift every pixel in the crop and nothing would raise."""
    import tifffile
    p = _tiled_tif(tmp_path / "t.tif", h=64, w=96, tile=32)
    with tifffile.TiffFile(str(p)) as tf:
        page = tf.series[0].pages[0]
        # a grid starting 4 px above and 8 px left of the page
        grid = RZ.GridSpec(origin_x_um=-8 * 0.25, origin_y_um=-4 * 0.25,
                           pixel_um=0.25, height=20, width=24)
        img, valid = TL.read_window_on_grid(page, grid, 0.25)

    assert img.shape == grid.shape and valid.shape == grid.shape
    assert not valid[:4, :].any() and not valid[:, :8].any()      # off-page rows/cols
    assert valid[4:, 8:].all()
    # the on-page part is the page's own top-left corner, placed at (4, 8)
    assert img[4, 8] == 1
    assert (img[:4, :] == 0).all() and (img[:, :8] == 0).all()


def test_a_chunk_the_file_never_stored_is_unavailable_not_black(tmp_path):
    """Zero is a legal intensity. Without the mask nothing downstream -- statistic or figure --
    can tell an unwritten tile from measured background."""
    import tifffile
    p = _tiled_tif(tmp_path / "s.tif", h=64, w=96, tile=32, sparse={(0, 1)})
    with tifffile.TiffFile(str(p)) as tf:
        page = tf.series[0].pages[0]
        grid, _px = TL.window_grid((0.0, 0.0, 96 * 0.25, 64 * 0.25), 0.25)
        img, valid = TL.read_window_on_grid(page, grid, 0.25)

    assert img.shape == (64, 96)
    assert (img[0:32, 32:64] == 0).all()          # the sparse tile decoded to zeros ...
    assert not valid[0:32, 32:64].any()           # ... and is marked unavailable
    assert valid[0:32, 0:32].all() and valid[32:64, 0:32].all()
    assert (img[0:32, 0:32] == 1).all()


def test_the_grid_shares_the_images_own_pixel_lattice(tmp_path):
    """Snapped outward to whole pixels, so geometry and pixels are never resampled and a
    half-pixel cannot creep in between them."""
    grid, (y0, y1, x0, x1) = TL.window_grid((10.3, 20.7, 30.1, 44.9), 0.25)
    assert (x0, y0) == (41, 82)
    assert grid.origin_x_um == pytest.approx(41 * 0.25)
    assert grid.origin_y_um == pytest.approx(82 * 0.25)
    assert float(grid.origin_x_um / 0.25).is_integer()
    assert grid.width == x1 - x0 and grid.height == y1 - y0
    # outward: the grid covers at least the requested bounds
    assert grid.origin_x_um <= 10.3 and grid.origin_x_um + grid.width * 0.25 >= 30.1
    with pytest.raises(TL.TileError, match="empty bounds"):
        TL.window_grid((5.0, 5.0, 5.0, 9.0), 0.25)


def test_the_grid_refuses_a_read_at_a_different_pixel_size(tmp_path):
    import tifffile
    p = _tiled_tif(tmp_path / "m.tif")
    with tifffile.TiffFile(str(p)) as tf:
        grid, _ = TL.window_grid((0.0, 0.0, 8.0, 8.0), 0.25)
        with pytest.raises(TL.TileError, match="pixel_um"):
            TL.read_window_on_grid(tf.series[0].pages[0], grid, 0.5)


def test_pyramid_level_choice_uses_the_recorded_scale_not_a_tag():
    shapes = [(82080, 40320), (41040, 20160), (20520, 10080), (10260, 5040)]
    assert TL.level_for(shapes, 0.25, 0.25) == 0
    assert TL.level_for(shapes, 0.25, 0.5) == 1
    assert TL.level_for(shapes, 0.25, 1.0) == 2


# --------------------------------------------------------------------------- #
# The alignment gate. Everything downstream is scoped by its verdict, so the sign
# convention and the thresholds are pinned by tests rather than by comments.
# --------------------------------------------------------------------------- #

import align as AL                                                # noqa: E402
import rasterize_labels as RL                                     # noqa: E402


def _nuclei_field(shape, centres_px, radius_px=6.0):
    """Discs at `centres_px`, softened, standing in for a DAPI field."""
    from scipy.ndimage import gaussian_filter
    img = np.zeros(shape, dtype=np.float64)
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    for cy, cx in centres_px:
        img[((yy - cy) ** 2 + (xx - cx) ** 2) <= radius_px ** 2] = 1.0
    return gaussian_filter(img, 1.5)


def _grid_of_centres(n=9, step=28, first=30):
    return [(first + step * i, first + step * j) for i in range(n) for j in range(n)]


def test_a_known_shift_is_recovered_with_the_documented_sign(tmp_path):
    """`dx_um`/`dy_um` is the shift to ADD to the warped mask to land it on the DAPI.

    Constructed so an inverted sign fails: the mask sits 3 px LEFT and 2 px UP of the DAPI, so
    the answer must be positive in both axes.
    """
    shape, px = (300, 300), 0.5
    centres = _grid_of_centres()
    dapi = _nuclei_field(shape, centres)
    mask = np.zeros(shape, dtype=bool)
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    for cy, cx in centres:                       # mask discs displaced by (-2, -3) px
        mask[((yy - (cy - 2)) ** 2 + (xx - (cx - 3)) ** 2) <= 6.0 ** 2] = True

    off = AL.offset_by_ncc(dapi, np.ones(shape, bool), mask, pixel_um=px,
                           n_nuclei=len(centres))
    assert off.reason == "measured", off.reason
    assert off.dx_um == pytest.approx(3 * px, abs=0.3)
    assert off.dy_um == pytest.approx(2 * px, abs=0.3)
    assert off.magnitude_um == pytest.approx(np.hypot(3, 2) * px, abs=0.3)
    assert off.ncc_peak > AL.MIN_NCC_PEAK


def test_a_perfectly_aligned_pair_returns_zero_not_a_small_bias():
    shape = (300, 300)
    centres = _grid_of_centres()
    dapi = _nuclei_field(shape, centres)
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    mask = np.zeros(shape, dtype=bool)
    for cy, cx in centres:
        mask[((yy - cy) ** 2 + (xx - cx) ** 2) <= 6.0 ** 2] = True
    off = AL.offset_by_ncc(dapi, np.ones(shape, bool), mask, pixel_um=0.5, n_nuclei=len(centres))
    assert off.reason == "measured"
    assert off.magnitude_um == pytest.approx(0.0, abs=0.15)


def test_no_correspondence_is_unavailable_not_zero():
    """Random noise against a real mask must not return a confident 0 um offset."""
    rng = np.random.default_rng(0)
    shape = (300, 300)
    noise = rng.normal(size=shape)
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    mask = np.zeros(shape, dtype=bool)
    for cy, cx in _grid_of_centres():
        mask[((yy - cy) ** 2 + (xx - cx) ** 2) <= 6.0 ** 2] = True
    off = AL.offset_by_ncc(noise, np.ones(shape, bool), mask, pixel_um=0.5, n_nuclei=81)
    assert off.magnitude_um is None
    assert off.reason in {"no nucleus-scale peak", "peak is ambiguous"}


def test_a_blank_window_is_named_rather_than_silently_excluded():
    shape = (120, 120)
    mask = np.zeros(shape, dtype=bool); mask[10:20, 10:20] = True
    off = AL.offset_by_ncc(np.zeros(shape), np.ones(shape, bool), mask,
                           pixel_um=0.5, n_nuclei=50)
    assert off.reason == "DAPI window blank or constant" and off.dx_um is None


def test_too_few_nuclei_is_refused_before_a_peak_is_believed():
    shape = (120, 120)
    off = AL.offset_by_ncc(np.zeros(shape), np.ones(shape, bool),
                           np.zeros(shape, bool), pixel_um=0.5, n_nuclei=3)
    assert off.reason == "too few nuclei to localise"


def test_the_self_overlap_thresholds_are_geometry_not_preference():
    """Two circles of radius r with centres t apart share
    (2/pi)(acos(u) - u*sqrt(1-u^2)) of one circle, u = t/2r."""
    for d in (9.083, 9.513):
        t50 = AL.self_overlap_offset_um(d, 0.50)
        t75 = AL.self_overlap_offset_um(d, 0.75)
        assert t75 < t50
        for t, want in ((t50, 0.50), (t75, 0.75)):
            u = t / d
            got = (2 / np.pi) * (np.arccos(u) - u * np.sqrt(1 - u * u))
            assert got == pytest.approx(want, abs=1e-6)
    assert AL.self_overlap_offset_um(9.083, 0.50) == pytest.approx(3.669, abs=0.002)
    assert AL.self_overlap_offset_um(9.513, 0.50) == pytest.approx(3.843, abs=0.002)


def test_the_verdict_is_a_pure_function_of_the_declared_thresholds():
    ok = dict(diameter_um=9.083, n_measured=100, n_total=125, contrast_fraction=0.9)
    assert AL.verdict(1.20, **ok)[0] == "PASS"
    assert AL.verdict(2.50, **ok)[0] == "INCONCLUSIVE"
    assert AL.verdict(4.50, **ok)[0] == "FAIL"
    # a minority of islets cannot carry a verdict, however good the median
    assert AL.verdict(0.1, diameter_um=9.083, n_measured=30, n_total=125,
                      contrast_fraction=0.9)[0] == "INCONCLUSIVE"
    # nor can a section where DAPI is not brighter inside the nuclei
    assert AL.verdict(0.1, diameter_um=9.083, n_measured=100, n_total=125,
                      contrast_fraction=0.1)[0] == "FAIL"
    # an estimator that cannot resolve the threshold cannot certify it
    assert AL.verdict(0.1, **ok, precision_um=5.0)[0] == "INCONCLUSIVE"
    # and the cell-mask fallback can never PASS
    assert AL.verdict(0.1, **ok, mode="cell")[0] == "INCONCLUSIVE"
    assert AL.verdict(None, **ok)[0] == "INCONCLUSIVE"


def test_contrast_returns_none_rather_than_inf_so_the_record_stays_writable():
    import json
    shape = (40, 40)
    mask = np.zeros(shape, bool); mask[10:20, 10:20] = True
    dapi = np.zeros(shape); dapi[mask] = 5.0            # gaps are exactly 0
    assert AL.contrast_inside_vs_gaps(dapi, np.ones(shape, bool), mask) is None
    dapi2 = np.full(shape, 1.0); dapi2[mask] = 4.0
    assert AL.contrast_inside_vs_gaps(dapi2, np.ones(shape, bool), mask) == pytest.approx(4.0)
    json.dumps({"contrast": AL.contrast_inside_vs_gaps(dapi, np.ones(shape, bool), mask)},
               allow_nan=False)


# --------------------------------------------------------------------------- #
# The honest-background rule, which until now was prose in another lane's docstring.
# --------------------------------------------------------------------------- #

import dataclasses                                                # noqa: E402

import display as DP                                              # noqa: E402


def test_the_black_point_is_zero_on_every_distribution():
    """Including one where a percentile black point would look most defensible: a channel
    whose 1st percentile is ~200 still starts its ramp at 0, because the zeros elsewhere in
    the section are the thing the panel exists to show."""
    for values in (np.random.default_rng(0).uniform(200, 255, 5000),
                   np.zeros(100),
                   np.r_[np.zeros(900), np.arange(1, 101)],
                   np.array([3.0])):
        r = DP.unsigned_display_range(values)
        assert r.vmin == 0.0
        assert r.vmax > 0.0


def test_there_is_no_black_point_knob_to_turn():
    """A `vmin=0` DEFAULT is one keyword away from being set. The repo's one non-compliant
    display path got there exactly that way."""
    import inspect
    params = inspect.signature(DP.unsigned_display_range).parameters
    assert not any("min" in p or "low" in p or "black" in p for p in params)
    r = DP.unsigned_display_range(np.arange(10.0))
    with pytest.raises(DP.DisplayError, match="measured floor"):
        dataclasses.replace(r, vmin=12.0)


def test_the_zero_inflation_survives_the_window_and_is_counted():
    """Pinning the black point is only worth something if what it protects is reported."""
    r = DP.unsigned_display_range(np.r_[np.zeros(90), np.arange(1, 11.0)])
    assert r.frac_at_zero == pytest.approx(0.9)
    assert "=0 90.0%" in r.caption() and "gamma 1" in r.caption()
    assert r.as_dict()["gamma"] is None


def test_an_all_zero_channel_falls_back_with_a_reason_not_a_zero_range():
    r = DP.unsigned_display_range(np.zeros(500))
    assert r.vmax == 255.0 and r.fallback_reason is not None
    assert "at the floor" in r.fallback_reason
    r2 = DP.unsigned_display_range(np.array([np.nan, np.nan]))
    assert r2.n_finite == 0 and r2.n_nonfinite == 2 and r2.fallback_reason is not None


def test_a_signed_score_is_windowed_symmetrically_about_zero():
    """The RNA contrast is a difference, so 0 is its neutral value and must sit at the centre
    of the colour map. The two arms are counted separately because they are not the same size."""
    from matplotlib.colors import Normalize
    s = DP.symmetric_display_range(np.r_[-np.arange(1, 50.0), np.arange(1, 200.0)])
    assert s.vmin == -s.vmax
    assert Normalize(s.vmin, s.vmax)(0.0) == pytest.approx(0.5)
    assert s.n_below != s.n_above          # reported apart, never pooled
    assert "0 at centre" in s.caption()


def test_unavailable_has_its_own_colour_rather_than_the_bottom_of_the_ramp():
    """A NaN drawn at the colormap's minimum is indistinguishable from a real minimum -- the
    same failure as writing 0 for a missing measurement."""
    from matplotlib.colors import to_rgba
    cmap = DP.mapped_colormap("inferno")
    assert to_rgba(cmap.get_bad()) == to_rgba(DP.UNAVAILABLE)
    assert to_rgba(cmap.get_bad()) != to_rgba(cmap(0.0))


def test_the_white_point_does_not_move_with_how_much_is_drawn():
    """The window is computed over the selection-independent population, so choosing to draw
    four islets instead of six cannot change the ramp."""
    rng = np.random.default_rng(1)
    population = rng.gamma(2.0, 20.0, 20000)
    a = DP.unsigned_display_range(population)
    b = DP.unsigned_display_range(population)
    assert a.vmax == b.vmax
    subset = DP.unsigned_display_range(population[:2000])
    assert subset.vmax != a.vmax        # which is exactly why the population must be fixed


# --------------------------------------------------------------------------- #
# Turning VALIS's point mapping into the lane's transform. What matters is not that
# an affine can be fitted -- one always can -- but that the fit REPORTS when it does
# not describe the mapping, instead of flattening a deformation into six parameters.
# --------------------------------------------------------------------------- #

import valis_transform as VT                                      # noqa: E402


def _valis_dir(tmp_path, src, dst, donor="6442"):
    d = tmp_path / f"donor_id={donor}"
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "xenium_centroids_um.npy", src)
    np.save(d / "warped_if_um.npy", dst)
    (d / "export.json").write_text(json.dumps(
        {"donor_id": donor, "serial": "0059947", "roi": "panc",
         "xenium_um_per_px": 0.85, "if_um_per_px": 1.0}))
    return d


def test_an_affine_mapping_is_recovered_and_declared_affine(tmp_path):
    rng = np.random.default_rng(0)
    src = rng.uniform(0, 15000, size=(5000, 2))
    A = np.array([[1.006, -0.009, 120.0], [0.009, 1.004, -75.0]])
    rec = VT.transform_one(_valis_dir(tmp_path, src, src @ A[:, :2].T + A[:, 2]))
    assert rec["is_affine"]
    assert np.allclose(np.asarray(rec["affine"]), A, atol=1e-6)
    assert rec["affine_residual_um"]["median"] < 1e-6
    assert "displacement" not in rec


def test_a_deformation_is_not_flattened_into_an_affine(tmp_path):
    """A smooth non-affine warp must be REPORTED, not absorbed. Six parameters will always
    fit something; the residual is what says whether that something is the mapping."""
    rng = np.random.default_rng(1)
    src = rng.uniform(0, 15000, size=(20000, 2))
    A = np.array([[1.0, 0.0, 10.0], [0.0, 1.0, -5.0]])
    dst = src @ A[:, :2].T + A[:, 2]
    dst[:, 0] += 8.0 * np.sin(src[:, 1] / 2500.0)      # a bend an affine cannot express
    rec = VT.transform_one(_valis_dir(tmp_path, src, dst))
    assert not rec["is_affine"]
    assert rec["affine_residual_um"]["median"] > VT.AFFINE_RESIDUAL_MAX_UM
    assert "displacement" in rec
    z = np.load(rec["displacement"]["file"])
    assert z["field"].shape[0] == 2
    assert rec["displacement"]["max_um"] > 1.0


def test_the_displacement_field_actually_removes_the_residual(tmp_path):
    """Keeping a field is only worth doing if applying it helps."""
    from phenocycler.integration.rasterize import GridSpec
    from phenocycler.integration.transform import Transform
    rng = np.random.default_rng(2)
    src = rng.uniform(0, 15000, size=(20000, 2))
    dst = src.copy()
    dst[:, 0] += 8.0 * np.sin(src[:, 1] / 2500.0)
    rec = VT.transform_one(_valis_dir(tmp_path, src, dst))
    z = np.load(rec["displacement"]["file"])
    tf = Transform(affine=np.asarray(rec["affine"]), displacement=z["field"],
                   grid=GridSpec.from_dict(rec["displacement"]["grid"]))
    before = np.median(np.linalg.norm(
        dst - (src @ np.asarray(rec["affine"])[:, :2].T + np.asarray(rec["affine"])[:, 2]),
        axis=1))
    after = np.median(np.linalg.norm(dst - tf.apply(src), axis=1))
    assert after < before / 2.0


def test_non_finite_points_are_dropped_and_counted_not_fitted(tmp_path):
    """VALIS returns NaN for points it cannot map. Feeding those to a least-squares fit
    poisons every parameter, so they are removed and the count is recorded."""
    rng = np.random.default_rng(3)
    src = rng.uniform(0, 15000, size=(2000, 2))
    A = np.array([[1.002, 0.0, 5.0], [0.0, 1.001, -3.0]])
    dst = src @ A[:, :2].T + A[:, 2]
    dst[:17] = np.nan
    rec = VT.transform_one(_valis_dir(tmp_path, src, dst))
    assert rec["n_dropped_nonfinite"] == 17
    assert rec["n_points"] == len(src) - 17
    assert np.allclose(np.asarray(rec["affine"]), A, atol=1e-6)


# --------------------------------------------------------------------------- #
# The pixel-to-micron handoff: two pyramids, built differently
# --------------------------------------------------------------------------- #

import export_dapi as ED                                          # noqa: E402
import registration_inputs as RI                                 # noqa: E402


def _registration_inventory(tmp_path, *, donor="6436", serial="0059907"):
    images = tmp_path / "images"
    images.mkdir()
    (images / f"{serial}_Scan1.qptiff").touch()
    inventory = tmp_path / "inventory.csv"
    inventory.write_text(
        "image_file,donor_id,roi,xenium_serial,provenance_source\n"
        f"{serial}_Scan1.qptiff,{donor},panc,{serial},transfer-ledger.csv\n"
    )
    xenium = tmp_path / "xenium_paths.csv"
    source = f"batch/run/output-XETG__{serial}__Panc__20260101__000000"
    xenium.write_text(
        "donor_id,roi,source_path\n"
        f"{donor},panc,{source}\n"
    )
    return inventory, xenium, images, source


def test_registration_inventory_cross_checks_serial_and_reads_image_facts(tmp_path, monkeypatch):
    inventory, xenium, images, source = _registration_inventory(tmp_path)
    monkeypatch.setattr(RI, "qptiff_registration_facts", lambda _path: {
        "pixel_um": 0.251,
        "dapi_page": 2,
        "channel_names": ("Opal 480", "Opal 520", "DAPI"),
        "pyramid_shapes": ((100, 80), (50, 40)),
        "dtype": "uint8",
    })

    records = RI.load_registration_inputs(
        inventory, xenium, images, tmp_path / "geojson", tmp_path / "xenium-root"
    )
    record = records["6436"]
    assert record.serial == "0059907"
    assert record.pixel_um == pytest.approx(0.251)
    assert record.dapi_page == 2
    assert record.bundle == tmp_path / "xenium-root" / source
    assert record.provenance_source == "transfer-ledger.csv"


def test_registration_inventory_refuses_donor_disagreement(tmp_path, monkeypatch):
    inventory, xenium, images, _source = _registration_inventory(tmp_path)
    xenium.write_text(xenium.read_text().replace("6436,panc", "6442,panc"))
    monkeypatch.setattr(RI, "qptiff_registration_facts", lambda _path: {})

    with pytest.raises(RI.RegistrationInputError, match="says donor 6436.*says 6442"):
        RI.load_registration_inputs(
            inventory, xenium, images, tmp_path / "geojson", tmp_path / "xenium-root"
        )


def test_registration_selection_enforces_central_exclusions():
    def record(donor):
        return RI.RegistrationInput(
            donor_id=donor,
            serial=f"00{donor}",
            roi="panc",
            qptiff=pathlib.Path(f"{donor}.qptiff"),
            bundle=pathlib.Path(donor),
            geojson=None,
            pixel_um=0.25,
            dapi_page=0,
            channel_names=("DAPI",),
            pyramid_shapes=((100, 100),),
            provenance_source="test",
        )

    records = {donor: record(donor) for donor in ("6436", "6457", "6579")}
    assert [item.donor_id for item in RI.select_registration_inputs(records, None)] == ["6436"]
    with pytest.raises(ValueError, match="6457"):
        RI.select_registration_inputs(records, ["6457"])


def test_the_pyramid_factor_is_the_integer_not_the_shape_ratio():
    """Donor 6436's own level shapes: 68119/17029 is 4.000176, and the factor is 4.

    That 4.4e-5 is not cosmetic. It is a scale error, so it is a displacement proportional to
    the coordinate -- zero at the origin and 0.61 um at the far corner of a 14.5 mm section,
    which is a third of the residual this lane is trying to remove.
    """
    shapes = [(68119, 37187), (34059, 18593), (17029, 9296), (8514, 4648)]
    assert ED.pyramid_factor(shapes, 2) == 4
    ratio = shapes[0][0] / shapes[2][0]
    assert ratio != 4.0
    drift = abs(ratio - 4) / 4 * (shapes[2][0] * 0.85)
    assert 0.5 < drift < 0.7                       # microns at the far corner

    # 6442's pyramid divides exactly and must give the same answer.
    assert ED.pyramid_factor([(74968, 45350), (37484, 22675), (18742, 11337)], 2) == 4


def test_a_pyramid_that_is_not_a_clean_downsample_is_refused():
    with pytest.raises(ValueError):
        ED.pyramid_factor([(1000, 1000), (300, 300)], 1)     # 3.33x, not an integer step


def test_the_two_pyramids_get_different_origins_because_they_are_built_differently():
    """Measured on donor 6436, three tiles each, level 0 against level 2:

        Xenium morphology  decimation phase (0,0) r=0.9936   box mean r=0.9766
        IF qptiff (DAPI)   box mean               r=0.99992  best decimation r=0.9891

    A decimated level KEEPS level-0 pixel ``f*k``, centred at ``(f*k + 0.5) * base``; an
    averaged one covers ``f*k .. f*k+f-1``, centred at ``(f*k + f/2) * base``. Left uncorrected
    the difference is a constant, and a constant is exactly what this lane is chasing.
    """
    assert ED.origin_um(0.2125, 4, averaged=False) == pytest.approx(0.10625)
    assert ED.origin_um(0.25, 4, averaged=True) == pytest.approx(0.5)
    # Level 0 is the identity under both rules, whichever way the pyramid was built.
    assert ED.origin_um(0.2125, 1, averaged=False) == ED.origin_um(0.2125, 1, averaged=True)


# --------------------------------------------------------------------------- #
# Nucleus detection and the point-set refinement
# --------------------------------------------------------------------------- #

import nuclei as NU                                               # noqa: E402
import refine as RF                                               # noqa: E402
from phenocycler.integration.transform import Transform           # noqa: E402


def _dapi(centres_xy, shape=(240, 240), sigma=2.0, amplitude=200.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    img = np.zeros(shape, np.float32)
    for cx, cy in centres_xy:
        img += amplitude * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2)))
    return img


def test_a_planted_nucleus_is_found_at_the_micron_position_the_origin_implies():
    """Pixel k is at ``k * um_per_px + origin_um`` -- the detector must not reinvent that."""
    img = _dapi([(60, 50), (140, 120), (170, 33)])
    xy, z = NU.detect(img, 1.0, 0.5, min_response_z=3.0)
    assert len(xy) == 3
    got = xy[np.lexsort((xy[:, 1], xy[:, 0]))]
    want = np.array([[60.5, 50.5], [140.5, 120.5], [170.5, 33.5]])
    assert np.allclose(got, want, atol=0.1)
    assert (z > 10).all()

    # The same image read at a different declared origin moves by exactly that difference.
    xy2, _ = NU.detect(img, 1.0, 0.10625, min_response_z=3.0)
    shift = np.median(xy2[np.lexsort((xy2[:, 1], xy2[:, 0]))] - got, axis=0)
    assert np.allclose(shift, 0.10625 - 0.5, atol=1e-6)


def test_noise_alone_yields_no_nuclei_rather_than_a_field_of_them():
    rng = np.random.default_rng(0)
    img = rng.normal(50.0, 1.0, size=(240, 240)).astype(np.float32)
    xy, _ = NU.detect(img, 1.0, 0.5, min_response_z=8.0)
    assert len(xy) == 0


def test_mutual_matching_leaves_an_orphan_unpaired_rather_than_pairing_it_with_a_stranger():
    """Plain nearest-neighbour would give the orphan a partner and call the distance a residual."""
    moving = np.array([[0.0, 0.0], [10.0, 0.0], [200.0, 200.0]])
    fixed = np.array([[0.3, 0.0], [10.3, 0.0]])
    i, j, d = RF.mutual_pairs(moving, fixed, max_dist_um=4.5)
    assert sorted(i.tolist()) == [0, 1]                # the far point is not paired
    assert np.allclose(sorted(d.tolist()), [0.3, 0.3])


def test_a_planted_translation_is_recovered_by_the_smallest_candidate():
    rng = np.random.default_rng(1)
    xen = rng.uniform(0, 4000, size=(4000, 2))
    truth = np.array([1.75, 2.18])
    ifn = xen + truth
    tf = RF.refine_translation(Transform(), xen, ifn)
    assert np.allclose(tf.translation_um, truth, atol=0.01)
    assert RF.score_transform(tf, xen, ifn)["pair_distance_um"]["p50"] < 0.01


def test_icp_recovers_a_similarity_and_survives_points_that_have_no_partner():
    """At the real nucleus density and the real seed error: 4,000 points over 600 um is ~9.5 um
    apart, and the planted similarity displaces by at most ~2.4 um -- inside ``PAIR_MAX_UM``,
    which is the condition the whole point-set route rests on. 500 nuclei the IF never detects
    and 400 spurious detections are added because both are true of the real fields.
    """
    rng = np.random.default_rng(2)
    xen = rng.uniform(0, 600, size=(4000, 2))
    th = np.radians(0.1)
    lin = 1.0003 * np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    ifn = xen @ lin.T + np.array([1.0, -0.6])
    assert np.linalg.norm(ifn - xen, axis=1).max() < RF.PAIR_MAX_UM
    ifn = np.vstack([ifn[:3500], rng.uniform(0, 600, size=(400, 2))])
    tf = RF.refine_icp(Transform(), xen, ifn)
    s = RF.score_transform(tf, xen, ifn)
    assert s["pair_distance_um"]["p50"] < 0.05
    assert tf.rotation_deg == pytest.approx(0.1, abs=0.02)
    assert tf.scale == pytest.approx(1.0003, abs=1e-4)


def test_icp_cannot_bootstrap_from_a_seed_outside_the_matching_radius():
    """The precondition, stated as a test rather than assumed. Nearest-neighbour matching is
    only meaningful while the seed is inside roughly half the spacing between nuclei; past that
    every point is matched to a stranger and the fit is confidently wrong. VALIS's seed is at
    1.4-2.8 um against a ~9 um spacing, which is why this route is available at all.
    """
    rng = np.random.default_rng(6)
    xen = rng.uniform(0, 600, size=(4000, 2))
    ifn = xen + np.array([30.0, 0.0])                 # three nuclei out
    tf = RF.refine_icp(Transform(), xen, ifn)
    assert np.linalg.norm(tf.translation_um) < 5.0    # it does not find the 30 um shift
    assert RF.score_transform(tf, xen, ifn)["pair_distance_um"].get("p50", 99) > 1.0


def test_the_score_carries_a_null_that_says_what_chance_would_give():
    """Nuclei are ~9 um apart, so a displaced field still has neighbours. The null is what
    stops a matched fraction being read as agreement when it is only density."""
    rng = np.random.default_rng(3)
    xen = rng.uniform(0, 3000, size=(5000, 2))
    ifn = xen + 0.2
    s = RF.score_transform(Transform(), xen, ifn, null_offset_um=25.0)
    assert s["matched_within"]["2um"] > 0.99
    assert s["matched_within_null"]["2um"] < 0.25
    assert s["nn_distance_um_null"]["p50"] > 10 * s["nn_distance_um"]["p50"]


def test_the_local_field_is_sampled_at_the_node_it_was_binned_at():
    """A field whose bins are half a bin from the sampler's nodes cannot be corrected later:
    the displacement is added AFTER the affine, so nothing upstream sees it."""
    rng = np.random.default_rng(4)
    xen = rng.uniform(0, 6000, size=(20000, 2))
    bend = np.column_stack([3.0 * np.sin(xen[:, 1] / 1200.0), np.zeros(len(xen))])
    ifn = xen + bend
    before = RF.score_transform(Transform(), xen, ifn)["pair_distance_um"]["p50"]
    tf = RF.refine_local(Transform(), xen, ifn, pixel_um=250.0, max_dist_um=6.0)
    after = RF.score_transform(tf, xen, ifn, max_dist_um=6.0)["pair_distance_um"]["p50"]
    assert before > 1.5
    assert after < 0.25 * before
    # A half-bin error would show up as a residual with the bend's own period still in it.
    assert tf.grid.pixel_um == 250.0


def test_a_flat_residual_is_not_given_a_field_shaped_like_its_noise():
    rng = np.random.default_rng(5)
    xen = rng.uniform(0, 4000, size=(20000, 2))
    ifn = xen + rng.normal(0.0, 0.05, size=xen.shape)
    tf = RF.refine_local(Transform(), xen, ifn, pixel_um=250.0, max_dist_um=6.0)
    assert tf.max_displacement_um < 0.1
