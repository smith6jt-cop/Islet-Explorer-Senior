"""Nucleus-detection behavior for GeoJSON files without a Tissue annotation."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import nuclei as NU  # noqa: E402


def test_islet_only_geojson_uses_the_whole_slide(tmp_path, monkeypatch):
    import tifffile

    donor_dir = tmp_path / "donor"
    donor_dir.mkdir()
    meta = {
        "donor_id": "6539", "serial": "0059865", "roi": "panc",
        "if_um_per_px": 1.0, "if_origin_um": 0.5, "xenium_origin_um": 0.5,
        "if": {"level": 1}, "xenium_bundle": str(tmp_path / "bundle"),
    }
    (donor_dir / "export.json").write_text(json.dumps(meta))
    monkeypatch.setattr(tifffile, "imread", lambda _path: np.zeros((4, 4)))
    monkeypatch.setattr(
        NU, "detect",
        lambda *_args, **_kwargs: (np.array([[1.0, 2.0], [3.0, 4.0]]), np.array([5.0, 6.0])))
    monkeypatch.setattr(
        NU, "xenium_nucleus_centroids",
        lambda _bundle: (np.array([[1.0, 2.0]]), np.array(["cell-1"])))
    monkeypatch.setattr(NU.geo, "load_annotations",
                        lambda _path: {NU.geo.ISLET: [object()]})

    geojson = tmp_path / "islets.geojson"
    rec = NU.detect_one(donor_dir, gj=geojson)

    assert rec["if"]["spatial_filter"] == "whole_slide_no_tissue_annotation"
    assert rec["if"]["geojson"] == str(geojson)
    assert "n_in_tissue" not in rec["if"]
    assert np.load(donor_dir / "if_nuclei_um.npy").shape == (2, 2)
