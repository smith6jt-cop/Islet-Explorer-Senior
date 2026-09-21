"""The AnnData attachment for TissUUmaps, and the scale error it is built to refuse.

`build_tissuumaps_h5ad.py` writes one AnnData per donor; `build_tissuumaps_project.py --h5ad`
attaches it to the .tmap. The thing worth testing is the coordinate factor, because getting it
wrong is silent: cells still land on tissue, just not on the right cells, and only a
cell-level zoom shows it.

Measured 2026-09-21: `typing/cells` X_centroid/Y_centroid are at QuPath resolution #1
(0.5078 um/px), while the DZI tiles are at the full-res base (0.2539). The factor is 2.0. The
AnnData's own docstring claimed the coordinates were already in tile pixels, which would have
placed the whole cohort in the top-left quadrant.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
BUILDER = REPO / "scripts/senior/build_tissuumaps_project.py"


@pytest.fixture(scope="module")
def btp():
    spec = importlib.util.spec_from_file_location("btp_under_test", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_h5ad(path: Path, n=64, markers=("CD3", "INS", "GCG"), extent=(1000.0, 2000.0)):
    """A minimal AnnData in the shape the real builder writes."""
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({
        "cell_class": pd.Categorical(rng.choice(["Beta", "Alpha", "Immune"], n)),
        "resolved_by": pd.Categorical(rng.choice(["hormone", "cd99"], n)),
        "cell_um2": rng.uniform(20, 90, n),          # numeric: must NOT be groupable
    })
    obs.index = [f"uuid-{i}" for i in range(n)]
    a = ad.AnnData(X=rng.random((n, len(markers)), dtype=np.float32), obs=obs)
    a.var_names = list(markers)
    xy = rng.random((n, 2))
    xy[0] = [1.0, 1.0]
    xy[-1] = list(extent)                              # pin the extent exactly
    a.obsm["spatial"] = xy * np.array([extent[0], extent[1]])
    a.obsm["spatial"][-1] = list(extent)
    a.write_h5ad(path)
    return path


def test_facts_report_markers_and_only_categorical_obs(btp, tmp_path):
    f = btp.h5ad_facts(_write_h5ad(tmp_path / "d.h5ad"))
    assert f is not None
    assert f["markers"] == ["CD3", "INS", "GCG"]
    assert "cell_class" in f["obs_categoricals"]
    assert "cell_um2" not in f["obs_categoricals"], (
        "a float column would make one colour group per distinct value")
    assert f["extent"] == pytest.approx((1000.0, 2000.0))


def test_the_coordinate_factor_is_two_not_one(btp, tmp_path):
    """The whole point. 0.5078 / 0.2539 == 2.0, and it must reach expectedHeader."""
    assert btp.DEFAULT_REGION_UM_PER_PX / btp.DEFAULT_UM_PER_PX == pytest.approx(2.0, abs=1e-3)
    f = btp.h5ad_facts(_write_h5ad(tmp_path / "d.h5ad"))
    entries = btp.h5ad_marker_files("d.h5ad", f, 2.0)
    for e in entries:
        assert e["expectedHeader"]["coord_factor"] == 2.0, (
            "coord_factor 1 puts every cell at half its true position")


def test_the_entries_address_the_anndata_the_way_tissuumaps_reads_it(btp, tmp_path):
    f = btp.h5ad_facts(_write_h5ad(tmp_path / "d.h5ad"))
    calls, expression = btp.h5ad_marker_files("d.h5ad", f, 2.0)

    for e in (calls, expression):
        assert e["expectedHeader"]["X"] == "/obsm/spatial;0"
        assert e["expectedHeader"]["Y"] == "/obsm/spatial;1"
        assert e["path"] == "d.h5ad"

    # the call layer groups by a categorical; the marker layer colours by an X column
    assert calls["expectedHeader"]["gb_col"].startswith("/obs/")
    assert calls["expectedRadios"]["cb_gr"] is True
    assert expression["expectedHeader"]["cb_col"] == "/X;0"
    assert expression["expectedRadios"]["cb_col"] is True

    # one dropdown entry per groupable column, and per marker, in X order
    assert len(calls["dropdownOptions"]) == len(f["obs_categoricals"])
    assert [o["optionName"] for o in expression["dropdownOptions"]] == f["markers"]
    assert [o["expectedHeader.cb_col"] for o in expression["dropdownOptions"]] == [
        "/X;0", "/X;1", "/X;2"]


def test_a_file_without_spatial_is_declined_rather_than_attached(btp, tmp_path):
    ad = pytest.importorskip("anndata")
    a = ad.AnnData(X=np.zeros((4, 2), dtype=np.float32))
    p = tmp_path / "nospatial.h5ad"
    a.write_h5ad(p)
    assert btp.h5ad_facts(p) is None, "a project that opens empty is worse than no project"


def test_a_missing_file_is_declined_rather_than_raising(btp, tmp_path):
    assert btp.h5ad_facts(tmp_path / "absent.h5ad") is None
