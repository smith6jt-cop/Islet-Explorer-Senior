"""Dedicated pLN exports use the recorded roster and an explicit site root."""

import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import roots
from multiplex5000.pln_sources import PLNSourceError, load_exports


def _csv(path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _source(path, image, *, cells=2):
    pq.write_table(pa.table({
        "Image": [image + " - resolution #1"] * cells,
        "Object ID": [f"native-{i}" for i in range(cells)],
        "Centroid X µm": [1.0] * cells,
        "Centroid Y µm": [2.0] * cells,
        "Nucleus: Area µm^2": [3.0] * cells,
        "Cell: CD20: Mean": [4.0] * cells,
    }), path)


def _fixture(tmp_path):
    image = "HDL115pLN_Scan1.er.qptiff"
    root = tmp_path / "site"
    root.mkdir()
    export = root / "misleading-donor-999.parquet"
    _source(export, image)
    roster, index = tmp_path / "roster.csv", tmp_path / "index.csv"
    _csv(roster, ["image_file", "donor_id", "roi", "n_cells", "excluded",
                  "excluded_reason"], [{"image_file": image, "donor_id": "115",
                                        "roi": "pln_1", "n_cells": "251349",
                                        "excluded": "False", "excluded_reason": ""}])
    _csv(index, ["image_file", "source_path", "segmentation_id"], [
        {"image_file": image, "source_path": export.name,
         "segmentation_id": "recorded-v3"}])
    return roster, index, root, image


def test_exact_source_uses_roster_identity_and_keeps_old_count_as_comparison(tmp_path):
    roster, index, root, image = _fixture(tmp_path)
    (found,) = load_exports(roster, index, source_root=root, require_sources=True)
    assert (found.image_file, found.donor_id, found.roi) == (image, "115", "pln_1")
    assert found.source_path == root / "misleading-donor-999.parquet"
    assert (found.roster_cells, found.export_cells) == (251349, 2)
    assert found.segmentation_id == "recorded-v3"


@pytest.mark.parametrize("change,match", [
    (lambda rows: rows.append(dict(rows[0])), "repeat"),
    (lambda rows: rows[0].update(source_path="../another.parquet"), "relative"),
    (lambda rows: rows[0].update(source_path="/tmp/another.parquet"), "relative"),
    (lambda rows: rows[0].update(image_file="unrecorded.qptiff"), "differ"),
])
def test_index_refuses_ambiguous_or_outside_sources(tmp_path, change, match):
    roster, index, root, _ = _fixture(tmp_path)
    with index.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    change(rows)
    _csv(index, ["image_file", "source_path", "segmentation_id"], rows)
    with pytest.raises(PLNSourceError, match=match):
        load_exports(roster, index, source_root=root)


def test_selected_source_must_match_export_image_and_cohort_state(tmp_path):
    roster, index, root, image = _fixture(tmp_path)
    _source(root / "misleading-donor-999.parquet", "wrong.qptiff")
    with pytest.raises(PLNSourceError, match="source Image column"):
        load_exports(roster, index, source_root=root, require_sources=True)
    with roster.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0].update(excluded="True", excluded_reason="damaged image")
    _csv(roster, ["image_file", "donor_id", "roi", "n_cells", "excluded",
                  "excluded_reason"], rows)
    with pytest.raises(PLNSourceError, match="damaged image"):
        load_exports(roster, index, source_root=root, selected_images={image})


def test_tracked_index_covers_the_recorded_nineteen_images_without_local_files():
    repo = roots.REPO
    roster = repo / "core/phenocycler/data/integration/pln_node_annotations.csv"
    index = repo / "core/multiplex5000/sites/pln_qupath_exports_v3.csv"
    exports = load_exports(roster, index, source_root=Path("/tmp"))
    assert len(exports) == 19
    assert len({entry.image_file for entry in exports}) == 19
    assert len({entry.source_path for entry in exports}) == 19
