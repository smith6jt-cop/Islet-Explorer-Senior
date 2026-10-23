"""Fail-closed tests for the manifest-backed post-Xenium registration cohort."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

from phenocycler import cohort
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


@pytest.fixture
def fully_excluded(monkeypatch):
    """Every exclusion with no declared platform, as 6457's was before 2026-10-22."""
    for donor in list(cohort.EXCLUDED_PLATFORMS):
        monkeypatch.delitem(cohort.EXCLUDED_PLATFORMS, donor)


def test_automatic_selection_filters_central_exclusions(tmp_path, monkeypatch, fully_excluded):
    # Derived from the authority, never named here. This test asserted the pair {6457, 6579} and
    # so failed when 6579's exclusion was lifted on 2026-09-17 -- a correct change breaking a test
    # that had restated the list, which is the failure mode "no lane may restate it" prevents.
    # An exclusion with no declared platform covers the Xenium and the IF (fail closed).
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


def test_explicit_excluded_donor_is_refused(tmp_path, monkeypatch, fully_excluded):
    for donor in sorted(DONOR_EXCLUSIONS):
        rows = [{"image_file": f"{donor}.qptiff", "donor_id": donor, "roi": "panc",
                 "xenium_serial": "0037637", "provenance_source": "audited fixture"}]
        records = _inventory(tmp_path, monkeypatch, rows)
        with pytest.raises(ValueError, match="excluded donor"):
            RI.select_registration_inputs(records, [donor], context="test registration")


def test_a_phenocycler_only_exclusion_does_not_refuse_the_registration(tmp_path, monkeypatch):
    """6457's central exclusion covers its PhenoCycler data only (maintainer, 2026-10-22).

    A registration reads the Xenium DAPI and the post-Xenium IF DAPI, so 6457 is selected
    automatically and by name, while the PhenoCycler side still refuses it.
    """
    assert cohort.excluded_on("6457", "phenocycler")
    assert RI.registration_excluded_on("6457") == []
    rows = [{"image_file": f"{d}.qptiff", "donor_id": d, "roi": "panc", "xenium_serial": serial,
             "provenance_source": "audited fixture"}
            for d, serial in (("6436", "0059907"), ("6457", "0037637"))]
    records = _inventory(tmp_path, monkeypatch, rows)
    assert [r.donor_id for r in RI.select_registration_inputs(records, None)] == ["6436", "6457"]
    assert [r.donor_id for r in RI.select_registration_inputs(records, ["6457"])] == ["6457"]
    with pytest.raises(ValueError, match="excluded donor"):
        cohort.ensure_eligible_donors(["6457"], context="PhenoCycler")


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


def _valis_main(tmp_path, monkeypatch, argv):
    """Run valis_register.main over two donor directories, recording what it would register."""
    import valis_register as V

    for donor in ("6436", "6457"):
        (tmp_path / f"donor_id={donor}").mkdir(exist_ok=True)
    registered = []
    monkeypatch.setattr(V, "register_one", lambda d, **_kw: registered.append(d.name))
    assert V.main(["--root", str(tmp_path), *argv]) == 0
    return registered


def test_valis_registers_a_donor_whose_exclusion_covers_only_its_phenocycler(
        tmp_path, monkeypatch):
    """6457 (maintainer, 2026-10-22): its Xenium and IF DAPI are registered, by name or not."""
    assert _valis_main(tmp_path, monkeypatch, []) == ["donor_id=6436", "donor_id=6457"]
    assert _valis_main(tmp_path, monkeypatch, ["--donors", "6457"]) == ["donor_id=6457"]


def test_valis_still_refuses_an_exclusion_without_a_platform(tmp_path, monkeypatch,
                                                             fully_excluded):
    assert _valis_main(tmp_path, monkeypatch, []) == ["donor_id=6436"]
    with pytest.raises(ValueError, match="excluded donor.*6457 \\(on xenium, if\\)"):
        _valis_main(tmp_path, monkeypatch, ["--donors", "6457"])
