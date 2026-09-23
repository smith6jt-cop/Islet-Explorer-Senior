"""Cohort readback and synthetic checks for Xenium typing and bundle selection.

One test carries ``pytest.mark.data``: the cohort readback opens the acquired cohort
(``core/phenocycler/data/qupath_manifest.json``, through ``pins.donors``) and the donor ->
bundle manifest (``core/phenocycler/data/integration/xenium_paths.csv``, through
``pins.serial_map`` / ``pins.run_map``).  Neither is in version control, so on a checkout
without them it skips naming the file it wanted.  Every other test here is synthetic, or
writes the file it reads under ``tmp_path``.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import gate  # noqa: E402
import pins  # noqa: E402
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402


@pytest.mark.data
def test_pins_serial_and_run_maps_cover_the_cohort():
    for pinned in (pins.COHORT_MANIFEST, pins.XENIUM_PATHS):
        if not pinned.is_file():
            pytest.skip(f"{pinned} not present (data is gitignored)")
    donors = pins.donors()
    assert len(donors) == len(set(donors)) == 22
    assert {"112", "6533"} <= set(donors)
    assert not set(donors) & DONOR_EXCLUSIONS.keys()
    serials, runs = pins.serial_map(), pins.run_map()
    assert set(donors) <= serials.keys() & runs.keys()
    for donor in donors:
        assert pins.serial(donor) == serials[donor]
        assert pins.xenium_run(donor) == runs[donor]
    # New donors use the cohort MM036 bundles, not same-slide 0080870/0081366.
    for donor, expected in {
        "6436": ("0059907", "MM044"),
        "6442": ("0059947", "MM045"),
        "6533": ("0041393", "MM036"),
        "112": ("0041396", "MM036"),
    }.items():
        assert (serials[donor], runs[donor]) == expected
    assert len({runs[d] for d in donors}) == 12


def test_donors_follow_the_acquired_manifest_and_normalize_aliases(tmp_path, monkeypatch):
    manifest = tmp_path / "qupath_manifest.json"
    manifest.write_text(json.dumps({"expected_donors": ["0112", "112", "6533"]}))
    monkeypatch.setattr(pins, "COHORT_MANIFEST", manifest)
    assert pins.donors() == ("112", "6533")


@pytest.mark.parametrize("donor", sorted(DONOR_EXCLUSIONS))
def test_donors_cannot_bypass_central_exclusions(tmp_path, monkeypatch, donor):
    manifest = tmp_path / "qupath_manifest.json"
    manifest.write_text(json.dumps({"expected_donors": ["0" + donor]}))
    monkeypatch.setattr(pins, "COHORT_MANIFEST", manifest)
    with pytest.raises(ValueError, match="excluded donor"):
        pins.donors()


def test_missing_cohort_is_unavailable_instead_of_falling_back(tmp_path, monkeypatch):
    monkeypatch.setattr(pins, "COHORT_MANIFEST", tmp_path / "absent.json")
    with pytest.raises(pins.PinError, match="acquired cohort is unavailable"):
        pins.donors()


@pytest.fixture
def bundle_rows():
    rows = []
    for donor, cohort, same_slide in [("6533", "0041393", "0080870"),
                                     ("112", "0041396", "0081366")]:
        rows.extend([
            {"donor_id": donor, "roi": "panc", "batch": "6KTRAF",
             "source_path": f"HuPrime5kPanTissue/6KTRAF/MM036_test/output-XETG00298__{cohort}__panc/ "},
            {"donor_id": donor, "roi": "panc", "batch": "PhenoXeno",
             "source_path": f"/same-slide/output-XETG00298__{same_slide}__panc"},
        ])
    return pd.DataFrame(rows)


@pytest.mark.parametrize("reverse", [False, True])
def test_cohort_bundles_are_selected_even_with_same_slide_rows(tmp_path, bundle_rows, reverse):
    csv = tmp_path / "xenium_paths.csv"
    bundle_rows.iloc[::-1 if reverse else 1].to_csv(csv, index=False)
    assert pins.serial_map(csv) == {"6533": "0041393", "112": "0041396"}
    assert pins.run_map(csv) == {"6533": "MM036", "112": "MM036"}
    for donor, serial in pins.serial_map(csv).items():
        path = pins.source_path(donor, csv)
        assert path.startswith("HuPrime5kPanTissue/")
        assert path.endswith(f"__{serial}__panc")


def test_same_slide_cannot_replace_an_absent_cohort_bundle(tmp_path, bundle_rows):
    csv = tmp_path / "xenium_paths.csv"
    bundle_rows[bundle_rows.batch == "PhenoXeno"].to_csv(csv, index=False)
    with pytest.raises(pins.PinError, match="0 bundle rows"):
        pins.source_path("6533", csv)


@pytest.mark.parametrize("reader", [pins.serial_map, pins.run_map,
                                  lambda csv: pins.source_path("6533", csv)])
def test_duplicate_cohort_bundles_are_refused_by_every_map(tmp_path, bundle_rows, reader):
    csv = tmp_path / "xenium_paths.csv"
    pd.concat([bundle_rows, bundle_rows.iloc[[0]]], ignore_index=True).to_csv(csv, index=False)
    with pytest.raises(pins.PinError, match="duplicate panc cohort bundle rows"):
        reader(csv)


def test_validate_calls_refuses_a_missing_reason_and_an_unknown_label():
    df = pd.DataFrame({"cell_id": ["a", "b"], "lineage_native": ["Endocrine", gate.UNAVAILABLE], "celltype_fine": ["Beta", ""],
                       "confidence": [0.9, np.nan], "reason": ["", ""], "cycling": [False, False]})
    with pytest.raises(ValueError):
        gate.validate_calls(df)
    df.loc[1, "reason"] = "qc"
    gate.validate_calls(df)
    bad = df.copy(); bad.loc[0, "lineage_native"] = "Neural_X"
    with pytest.raises(ValueError):
        gate.validate_calls(bad)
