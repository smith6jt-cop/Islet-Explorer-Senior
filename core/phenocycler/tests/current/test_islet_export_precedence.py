"""Which islet export the geometry builder picks, when more than one matches.

Re-exports arrive in dated subdirectories beside the originals and under a different spelling --
``2026SeptNew/6579_Scan2.er.qptiff - resolution #1.geojson`` against June's
``Islets__6579_Scan2.er.qptiff_-_resolution__1.geojson``.  Both match the same stem, so the
builder has to choose, and the choice has to survive the maintainer adding more of them.

It also has to survive a *detection* export landing in that directory under an islet export's
name.  One such file is there now: 66 GiB of per-cell measurements written to
``2026SeptNew/6436_Scan1.er.qptiff - resolution #1.geojson``.  Reading it as islet outlines would
not fail -- it carries annotations too -- it would just be catastrophically slow and wrong about
what an islet is.  The discriminator is ``objectType``: islet exports are 100% ``annotation``,
that file carries ``cell``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

TYPING = Path(__file__).resolve().parents[2] / "scripts/typing"
sys.path.insert(0, str(TYPING))


def _feature(object_type: str, name: str = "Islet_1") -> dict:
    return {
        "type": "Feature",
        "geometry": {"type": "Polygon",
                     "coordinates": [[[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]]]},
        "properties": {"objectType": object_type, "name": name,
                       "classification": {"name": "Islet"}},
    }


def _write(path: Path, object_types, mtime: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"type": "FeatureCollection",
                                "features": [_feature(t) for t in object_types]}))
    import os
    os.utime(path, (mtime, mtime))
    return path


@pytest.fixture()
def builder(tmp_path, monkeypatch):
    import build_islet_geometry as B
    monkeypatch.setattr(B.S, "ISLETS", tmp_path)
    return B


def test_newest_matching_export_wins_and_subdirectories_are_searched(builder, tmp_path):
    _write(tmp_path / "Islets__6579_Scan2.er.qptiff_-_resolution__1.geojson",
           ["annotation"], mtime=1_000_000)
    newer = _write(tmp_path / "2026SeptNew/6579_Scan2.er.qptiff - resolution #1.geojson",
                   ["annotation"], mtime=2_000_000)
    assert builder.islet_export("6579_Scan2") == newer


def test_older_export_wins_when_it_is_the_only_one(builder, tmp_path):
    only = _write(tmp_path / "Islets__6380_Scan1.er.qptiff_-_resolution__1.geojson",
                  ["annotation"], mtime=1_000_000)
    assert builder.islet_export("6380_Scan1") == only


def test_detection_export_is_refused_and_the_search_falls_through(builder, tmp_path, capsys):
    islets = _write(tmp_path / "Islets__6436_Scan1.er.qptiff_-_resolution__1.geojson",
                    ["annotation"], mtime=1_000_000)
    _write(tmp_path / "2026SeptNew/6436_Scan1.er.qptiff - resolution #1.geojson",
           ["annotation", "cell", "cell"], mtime=2_000_000)

    assert builder.islet_export("6436_Scan1") == islets, (
        "the newer file is a detection export; the builder must fall through to the islet export")
    refusal = capsys.readouterr().err
    assert "[refuse]" in refusal and "objectType=cell" in refusal, (
        "a refused export must say so with its measured reason, never be skipped silently")


def test_all_candidates_refused_raises_rather_than_returning_nothing(builder, tmp_path):
    _write(tmp_path / "2026SeptNew/6436_Scan1.er.qptiff - resolution #1.geojson",
           ["cell"], mtime=2_000_000)
    with pytest.raises(FileNotFoundError, match="refused as a detection export"):
        builder.islet_export("6436_Scan1")


def test_no_match_at_all_raises(builder):
    with pytest.raises(FileNotFoundError, match="found 0 files"):
        builder.islet_export("6999_Scan1")
