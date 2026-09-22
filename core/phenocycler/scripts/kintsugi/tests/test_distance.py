"""Step 0's claims, pinned so a re-run cannot quietly change them.

These read the written parquet rather than re-extracting: the extraction is a 212 GB scan and
its result is the thing under test.

Everything but the synthetic rendering test carries ``pytest.mark.data``: the report tests read
``data/development/kintsugi/distance_validation.json`` and the two ring-width tests re-measure
the islet annotation GeoJSON under ``data/Islets`` through the cohort manifest.  None of those
is in version control, so each group takes a fixture (``report``, ``ring_geometry``) that skips,
naming the file it wanted, when that file is absent: a data-lane test that cannot find its data
skips with the reason, never fails.
"""
import json
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import distance                                     # noqa: E402
import ring_width                                   # noqa: E402
from phenocycler.cohort import DONOR_EXCLUSIONS     # noqa: E402

REPORT = distance.LANE_OUT / "distance_validation.json"
#: What ``ring_width.measure`` opens besides the donor's GeoJSON (``ring_width._pixel_size_um``).
MANIFEST = ring_width.panel.PKG / "data" / "qupath_manifest.json"


def _report_is_current() -> tuple[bool, str]:
    """Is the written report readable AND on the current schema?

    Regenerating it needs ``expression`` from the pinned run, so when that run is not on disk
    the report cannot be refreshed and a stale one must not be asserted against as though it
    were.  Skipping with the cause named is the honest state; silently accepting either
    schema would hide exactly the drift these tests exist to catch.
    """
    if not REPORT.exists():
        return False, "run scripts/kintsugi/run_step0.py first"
    donors = json.loads(REPORT.read_text()).get("donors") or []
    if donors and "ring_width_um_median" not in donors[0]:
        return False, (
            "distance_validation.json predates the ring_width_um_median rename. Re-run "
            "`run_step0.py --validate-only` once the run pin is resolved -- it needs "
            f"expression from {distance.panel.RUN_ID}, which is not currently on disk.")
    return True, ""


@pytest.fixture(scope="module")
def report():
    """The written Step 0 report -- skips (with the cause) if it is absent or stale.

    Scoped to the fixture rather than the module so the tests that read only the annotation
    GeoJSON keep running: those do not depend on a run being on disk, and take
    ``ring_geometry`` instead.
    """
    current, why = _report_is_current()
    if not current:
        pytest.skip(why)
    return json.loads(REPORT.read_text())


@pytest.fixture
def ring_geometry():
    """``ring_width.measure`` for one donor -- skips, naming the file, when its input is absent.

    The measurement opens the donor's islet annotation GeoJSON under ``data/Islets`` and the
    QuPath manifest for the pixel size; neither is in version control.  Only ABSENCE skips: a
    GeoJSON that is present and measures wrong still fails, which is the finding these tests
    pin.
    """
    def _measure(donor_id: str) -> dict:
        geojson = ring_width._geojson_by_donor().get(str(donor_id))
        if geojson is None:
            pytest.skip(f"no islet GeoJSON for donor {donor_id} under {ring_width.ISLETS} "
                        "(data is gitignored)")
        if not MANIFEST.is_file():
            pytest.skip(f"{MANIFEST} not present (data is gitignored)")
        return ring_width.measure(donor_id)
    return _measure


@pytest.mark.data
def test_every_object_id_resolves(report):
    """The join is on the QuPath UUID, so an unmatched cell is a contract failure."""
    unmatched = {r["donor_id"]: r["unmatched_object_id"]
                 for r in report["donors"] if r["unmatched_object_id"]}
    assert not unmatched, f"unmatched object_id: {unmatched}"


@pytest.mark.data
def test_all_twenty_donors_present(report):
    """The count is the written report's; the exclusion list is ``phenocycler.cohort``'s.

    This used to restate ``("6457", "6579")`` as a literal.  6579's exclusion was lifted on
    2026-09-17 and only ``cohort.py`` may say who is excluded, so the check is derived from
    ``DONOR_EXCLUSIONS`` -- a re-run of Step 0 that admits 6579 changes the count, and the
    pin at 20 is what makes that a visible decision rather than a silent one.
    """
    donors = {r["donor_id"] for r in report["donors"]}
    assert len(donors) == 20
    for excluded in DONOR_EXCLUSIONS:
        assert excluded not in donors, f"{excluded}: {DONOR_EXCLUSIONS[excluded]}"


@pytest.mark.data
def test_sign_convention_holds_against_cell_region(report):
    """cell_region was derived independently at ingest, so this is a real cross-check."""
    bar = distance.MAX_VIOLATION_RATE
    bad = {r["donor_id"]: r["core_positive_rate"]
           for r in report["donors"] if r["core_positive_rate"] > bar}
    assert not bad, f"core cells reading a positive distance above {bar}: {bad}"


@pytest.mark.data
def test_no_cell_lies_beyond_its_donors_own_ring(report):
    """Checked against the MEASURED width, never the nominal 20 um the name asserts."""
    bar = distance.MAX_VIOLATION_RATE
    bad = {r["donor_id"]: r["peri_beyond_ring_rate"]
           for r in report["donors"] if r["peri_beyond_ring_rate"] > bar}
    assert not bad, f"peri cells beyond the measured ring: {bad}"


@pytest.mark.data
def test_donor_6482_ring_is_25um_and_is_recorded_as_such(report):
    """The one donor whose expansion is not the width its annotation is named for.

    Pinned because it is invisible in ``expression`` -- ``cell_region`` says 'peri' either
    way -- and because any cross-donor peri statistic that assumes 20 um is wrong for 6482.
    """
    by_donor = {r["donor_id"]: r for r in report["donors"]}
    assert by_donor["6482"]["ring_width_um_median"] == pytest.approx(25.06, abs=0.2)
    assert by_donor["6482"]["ring_matches_nominal"] is False
    others = [r for r in report["donors"] if r["donor_id"] != "6482"]
    assert all(r["ring_matches_nominal"] for r in others)
    for r in others:
        assert r["ring_width_um_median"] == pytest.approx(20.1, abs=0.2)


@pytest.mark.data
def test_6539_tissue_distance_is_available(report):
    """has_tissue_annotation:false is about parenting, not about the annotation existing."""
    by_donor = {r["donor_id"]: r for r in report["donors"]}
    assert by_donor["6539"]["tissue_available"] is True
    assert not by_donor["6539"]["tissue_unavailable_reason"]
    assert all(r["tissue_available"] for r in report["donors"])


@pytest.mark.data
def test_ring_width_is_measured_from_geometry_not_assumed(ring_geometry):
    """Re-measure one donor from the polygons; it must not read the nominal constant."""
    measured = ring_geometry("6482")
    assert measured["ring_width_um_median"] > ring_width.NOMINAL_RING_UM + 1.0
    assert measured["n_rings_measured"] > 0
    assert measured["n_rings_in_geojson"] > 400


@pytest.mark.data
def test_6539_rings_were_drawn_the_cells_were_not_assigned_to_them(ring_geometry):
    """Corrects the inherited "expansion run on 3 of 539" account.

    6539's GeoJSON carries 569 IsletExpanded polygons at a correct ~20.1 um median. What is
    missing is the cell-to-ring parent assignment, not the annotations -- which changes the
    remedy from "redraw 566 rings" to "re-resolve one donor's hierarchy".
    """
    measured = ring_geometry("6539")
    assert measured["n_rings_in_geojson"] > 500
    assert measured["ring_width_um_median"] == pytest.approx(20.1, abs=0.2)
    assert measured["matches_nominal"] is True


def test_the_report_line_renders_for_a_donor_whose_ring_is_not_nominal():
    """The non-nominal branch of Step 0's report line, which only donor 6482 reaches.

    ``run_step0`` formatted the flag from ``r['ring_width_um']`` while ``distance.validate``
    returns ``ring_width_um_median``.  A wrong key there is invisible on 19 of 20 donors --
    ``ring_matches_nominal`` is true, so the branch never runs -- and raises ``KeyError`` on
    the one donor the entire ring-width finding is about, which is also the one donor a
    reader would most want the flag for.

    This asserts the rendering against a synthetic row rather than the cohort, so it holds
    with no run on disk and cannot be satisfied by the branch simply not being taken.
    """
    import run_step0                                  # noqa: E402

    row = {"donor_id": "6482", "n_expression": 1_234_567, "unmatched_object_id": 0,
           "core_positive_rate": 0.01, "peri_nonpositive_rate": 0.01,
           "peri_beyond_ring_rate": 0.0, "ring_matches_nominal": False,
           "ring_width_um_median": 25.06}
    assert run_step0.format_validation_line(row).endswith("ring=25.1um!")

    row_ok = dict(row, donor_id="6374", ring_matches_nominal=True,
                  ring_width_um_median=20.1)
    assert "ring=" not in run_step0.format_validation_line(row_ok)
