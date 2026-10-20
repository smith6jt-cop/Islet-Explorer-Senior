"""Fail-closed tests for the manifest-backed post-Xenium registration cohort."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

from phenocycler.cohort import DONOR_EXCLUSIONS


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import registration_inputs as RI  # noqa: E402


def _facts(_path: Path, **_kw) -> dict:
    # **_kw because qptiff_registration_facts() takes keyword-only options -- allow_dtypes
    # today -- and a stub that pins the old signature turns any addition into five failures
    # that look like a data problem. These five had been failing for exactly that reason,
    # unnoticed, because scripts/*/tests is run by no CI workflow.
    return {
        "pixel_um": 0.25,
        # qptiff_registration_facts() records WHERE the size came from -- the XML tag or a
        # fitted inventory. Production reads it; the stub predated it.
        "pixel_size_source": "qptiff XML PixelSizeMicrons",
        "dapi_page": 0,
        "channel_names": ("DAPI", "Opal 570"),
        "pyramid_shapes": ((64, 48), (32, 24)),
        "dtype": "uint8",
    }


def _inventory(tmp_path: Path, monkeypatch, rows: list[dict]) -> dict[str, RI.RegistrationInput]:
    images = tmp_path / "images"
    images.mkdir()
    for row in rows:
        (images / row["image_file"]).touch()
    inventory = tmp_path / "registration.csv"
    pd.DataFrame(rows).to_csv(inventory, index=False)
    xenium = tmp_path / "xenium.csv"
    pd.DataFrame([
        {
            "donor_id": row["donor_id"],
            "roi": row["roi"],
            "source_path": f"run/output-XET__{row['xenium_serial']}__Panc__date/",
        }
        for row in rows
    ]).to_csv(xenium, index=False)
    monkeypatch.setattr(RI, "qptiff_registration_facts", _facts)
    return RI.load_registration_inputs(
        inventory, xenium, images, tmp_path / "geojson", tmp_path / "bundles")


def test_missing_geojson_is_valid_for_registration(tmp_path, monkeypatch):
    rows = [{
        "image_file": "slide.qptiff",
        "donor_id": "6436",
        "roi": "panc",
        "xenium_serial": "0059907",
        "provenance_source": "audited fixture",
    }]
    records = _inventory(tmp_path, monkeypatch, rows)
    assert records["6436"].geojson is None
    assert records["6436"].channel_names == ("DAPI", "Opal 570")


def test_automatic_selection_filters_central_exclusions(tmp_path, monkeypatch):
    # Derived from the authority, never named here. This test asserted the pair {6457, 6579} and
    # so failed when 6579's exclusion was lifted on 2026-09-17 -- a correct change breaking a test
    # that had restated the list, which is the failure mode "no lane may restate it" prevents.
    excluded = sorted(DONOR_EXCLUSIONS)
    assert excluded, "the cohort must carry at least one exclusion for this test to mean anything"
    keep = "6436"
    rows = [{"image_file": f"{keep}.qptiff", "donor_id": keep, "roi": "panc",
             "xenium_serial": "0059907", "provenance_source": "audited fixture"}]
    rows += [{"image_file": f"{d}.qptiff", "donor_id": d, "roi": "panc",
              "xenium_serial": f"00000{i:02d}", "provenance_source": "audited fixture"}
             for i, d in enumerate(excluded)]
    records = _inventory(tmp_path, monkeypatch, rows)
    assert [record.donor_id for record in RI.select_registration_inputs(records, None)] == [keep]


def test_explicit_excluded_donor_is_refused(tmp_path, monkeypatch):
    for donor in sorted(DONOR_EXCLUSIONS):
        rows = [{"image_file": f"{donor}.qptiff", "donor_id": donor, "roi": "panc",
                 "xenium_serial": "0037637", "provenance_source": "audited fixture"}]
        records = _inventory(tmp_path, monkeypatch, rows)
        with pytest.raises(ValueError, match="excluded donor"):
            RI.select_registration_inputs(records, [donor], context="test registration")


def test_xenium_donor_disagreement_is_refused(tmp_path, monkeypatch):
    rows = [{"image_file": "slide.qptiff", "donor_id": "6436", "roi": "panc",
             "xenium_serial": "0059907", "provenance_source": "audited fixture"}]
    images = tmp_path / "images"
    images.mkdir()
    (images / "slide.qptiff").touch()
    inventory = tmp_path / "registration.csv"
    pd.DataFrame(rows).to_csv(inventory, index=False)
    xenium = tmp_path / "xenium.csv"
    pd.DataFrame([{"donor_id": "6442", "roi": "panc",
                   "source_path": "run/output-XET__0059907__Panc__date/"}]).to_csv(
                       xenium, index=False)
    monkeypatch.setattr(RI, "qptiff_registration_facts", _facts)
    with pytest.raises(RI.RegistrationInputError, match="says donor 6436"):
        RI.load_registration_inputs(
            inventory, xenium, images, tmp_path / "geojson", tmp_path / "bundles")


def test_unmanifested_local_qptiff_is_refused(tmp_path, monkeypatch):
    rows = [{"image_file": "slide.qptiff", "donor_id": "6436", "roi": "panc",
             "xenium_serial": "0059907", "provenance_source": "audited fixture"}]
    images = tmp_path / "images"
    records = _inventory(tmp_path, monkeypatch, rows)
    assert records
    (images / "unrecorded.qptiff").touch()
    # Reuse the generated manifests, now against the enlarged local image set.
    with pytest.raises(RI.RegistrationInputError, match="unrecorded"):
        RI.load_registration_inputs(
            tmp_path / "registration.csv", tmp_path / "xenium.csv", images,
            tmp_path / "geojson", tmp_path / "bundles")
