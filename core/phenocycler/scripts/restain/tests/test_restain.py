"""The restain lane, on geometry whose answer is known.

Nothing here needs an export or a bundle: the point is that the machinery is checked against
synthetic shapes where the right answer is arithmetic, so it still runs when the data is not
on this machine.
"""
import json
import pathlib
import sys

import numpy as np
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
