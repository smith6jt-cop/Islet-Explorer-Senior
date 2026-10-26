"""Physical same-slide pairs are explicit; an unresolved frame cannot be used."""

import csv
import gzip
import hashlib
from pathlib import Path

import pytest
import yaml

from multiplex5000.same_slide_sources import SameSlideSourceError, load_pairs
from multiplex5000.xenium_sources import CELL_FILES, INSTRUMENT_FILES

REPO = Path(__file__).resolve().parents[4]
HPG = REPO / "core/multiplex5000/sites/hipergator/spleen_same_slide_20261026.yaml"


def test_recorded_spleen_pilot_has_eight_distinct_explicit_pairs():
    pairs = load_pairs(HPG)
    assert len(pairs) == 8
    assert {p.xenium_acquisition_id for p in pairs} == {
        "xenium:0110618:090", "xenium:0110618:082", "xenium:0110618:104",
        "xenium:0110618:014", "xenium:0110979:059", "xenium:0110979:009",
        "xenium:0110979:034", "xenium:0110979:080"}
    assert all(p.crop_origin_px is None for p in pairs)  # sources not verified on this host


def _synthetic(tmp_path):
    section = tmp_path / "sections.csv"
    with section.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["slide", "section", "bbox_x", "bbox_y", "donor_id", "entry_name"])
        writer.writerow(["0110618", "S1", "720", "416", "HDL090", "HDL090_PC67_0110618-S1"])
    header = b"object_id\tx_um\ty_um\tCell: CD20: Mean\n"
    export = tmp_path / "HDL090_PC67_0110618-S1.csv.gz"
    with gzip.open(export, "wb") as handle:
        handle.write(header)
    image = tmp_path / "0110618_Scan1.er.qptiff"
    image.write_bytes(b"image")
    output = tmp_path / "output-XETG00298__0110618__090__20260826__224731"
    output.mkdir()
    for name in (*INSTRUMENT_FILES, *CELL_FILES):
        (output / name).write_bytes(b"source")
    manifest = {
        "schema": "multiplex5000-same-slide-sources/1",
        "tissue": "spleen", "relationship": "same_section_after_phenocycler",
        "protein_frame": "crop_local_um", "xenium_frame": "region_local_um",
        "expected_pairs": 1,
        "section_table": {"path": str(section),
                          "sha256": hashlib.sha256(section.read_bytes()).hexdigest()},
        "protein_header_sha256": hashlib.sha256(header).hexdigest(),
        "pairs": [{"donor_id": "HDL090", "protein_entry": "HDL090_PC67_0110618-S1",
                   "protein_export": str(export), "protein_image": str(image),
                   "xenium_acquisition_id": "xenium:0110618:090",
                   "xenium_output": str(output)}],
    }
    path = tmp_path / "pairs.yaml"
    path.write_text(yaml.safe_dump(manifest))
    return path, manifest, section, output


def test_verified_pair_uses_section_table_for_crop_origin(tmp_path):
    path, _, _, _ = _synthetic(tmp_path)
    pair, = load_pairs(path, require_sources=True)
    assert pair.crop_origin_px == (720, 416)
    assert pair.tissue == "spleen"
    assert pair.relationship == "same_section_after_phenocycler"


def test_changed_section_table_or_missing_cell_file_blocks(tmp_path):
    path, _, section, output = _synthetic(tmp_path)
    section.write_text(section.read_text() + "\n")
    with pytest.raises(SameSlideSourceError, match="section table changed"):
        load_pairs(path, require_sources=True)
    second = tmp_path / "second"
    second.mkdir()
    path, _, _, output = _synthetic(second)
    (output / "cells.parquet").unlink()
    with pytest.raises(SameSlideSourceError, match="missing sources"):
        load_pairs(path, require_sources=True)


def test_a_pair_cannot_swap_xenium_region_for_the_same_donor(tmp_path):
    path, manifest, _, _ = _synthetic(tmp_path)
    manifest["pairs"][0]["xenium_acquisition_id"] = "xenium:0110618:082"
    path.write_text(yaml.safe_dump(manifest))
    with pytest.raises(SameSlideSourceError, match="acquisition differs"):
        load_pairs(path)
