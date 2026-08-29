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
