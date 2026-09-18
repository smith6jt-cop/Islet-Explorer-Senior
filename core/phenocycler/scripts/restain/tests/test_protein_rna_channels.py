"""Which channel is which marker, and at what pixel size -- resolved by NAME, never by index.

The post-Xenium IF lane took both facts from ``restain_sections.csv`` through
``sections[sections.xenium_serial == serial].iloc[0]``. That file has five data rows and
``registered_serials()`` returns thirteen donors, so eleven of thirteen raised ``IndexError`` on
that one line -- and the two that did not were read at a pixel size the CSV declares identically
for every row while the images themselves differ by a factor of two.

Most of what follows runs on synthetic QuPath entries, so it checks the machinery rather than the
state of a QuPath renaming pass that is happening while this is being written. The tests that do
read the real project assert PROPERTIES that stay true as more channels are renamed, never a
count that a rename would move.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import protein_rna as PR                                       # noqa: E402

#: The instrument's own channel names, which is all the qptiff XML ever carries. Read 2026-09-18,
#: all twenty post-Xenium IF images spell their channels exactly this way in their ``<Name>``
#: tags; the marker assignment exists only in the QuPath project.
INSTRUMENT_CHANNELS = ("DAPI", "Opal 570", "Opal 690", "Opal 520", "Sample AF")


# --------------------------------------------------------------------------------------- fixtures

def _entry(project: pathlib.Path, entry: str, image_file: str, names):
    """One QuPath project entry, exactly as QuPath writes it after a channel rename."""
    d = project / "data" / entry
    d.mkdir(parents=True, exist_ok=True)
    (d / "server.json").write_text(json.dumps({
        "builderType": "uri",
        "uri": f"file:/somewhere/Images/{image_file}",
        "metadata": {"name": f"{image_file} - resolution #1",
                     "channels": [{"name": n, "color": -1} for n in names]}}))
    return d


@pytest.fixture
def inputs_csv(tmp_path):
    path = tmp_path / "restain_registration_inputs.csv"
    path.write_text(
        "# synthetic\n"
        "image_file,donor_id,roi,xenium_serial,provenance_source\n"
        "renamed_a.qptiff,6493,panc,0037191,if_transfer_provenance.csv\n"
        "renamed_b.qptiff,6461,panc,0059847,if_transfer_provenance.csv\n"
        "untouched.qptiff,6482,panc,0037354,if_transfer_provenance.csv\n"
        "typo.qptiff,6380,panc,0037502,if_transfer_provenance.csv\n"
        "partial.qptiff,6457,panc,0037637,if_transfer_provenance.csv\n"
        "tenx.qptiff,6591,panc,0041116,restain_sections.csv\n")
    return path


@pytest.fixture
def sections_csv():
    """The one CSV row, with the blanket assignment the real file actually carries.

    One row rather than six, because that asymmetry is the defect: the CSV covers five of the
    twenty images, and every donor beyond them had no row at all for ``.iloc[0]`` to find.
    """
    return pd.DataFrame([
        {"image_file": "tenx.qptiff", "donor_id": "6591", "roi": "panc",
         "xenium_serial": "0041116", "pixel_size_um": "0.2505", "channel_dapi": "NUCLEUS",
         "channel_opal520": "GCG", "channel_opal570": "FAILED", "channel_opal690": "INS",
         "channel_sampleaf": "AF"},
    ])


@pytest.fixture
def project(tmp_path):
    p = tmp_path / "PostXeniumIHC"
    # Two renamed entries that put INS and GCG at DIFFERENT pages, which is the production case.
    _entry(p, "20", "renamed_a.qptiff", ["DAPI", "INS", "GCG", "Opal 520", "Sample AF"])
    _entry(p, "8", "renamed_b.qptiff", ["DAPI", "Opal 570", "INS", "GCG", "Sample AF"])
    _entry(p, "18", "untouched.qptiff", ["DAPI", "Opal 570", "Opal 690", "Opal 520", "Sample AF"])
    _entry(p, "15", "typo.qptiff", ["DAPI", "Opal 570", "CGC", "Opal 520", "Sample AF"])
    _entry(p, "13", "partial.qptiff", ["DAPI", "Opal 570", "GCG", "Opal 520", "Sample AF"])
    _entry(p, "5", "tenx.qptiff", ["DAPI", "INS", "GCG", "Opal 520", "Sample AF"])
    return p


@pytest.fixture
def images(tmp_path):
    d = tmp_path / "Images"
    d.mkdir()
    for name in ("renamed_a.qptiff", "renamed_b.qptiff", "untouched.qptiff",
                 "typo.qptiff", "partial.qptiff", "tenx.qptiff"):
        (d / name).write_bytes(b"")
    return d


@pytest.fixture
def image_facts():
    """What each image's XML says: the instrument's names, and the 20x/10x scale split.

    Stubbed rather than written into a fixture qptiff, for the reason ``test_restain.py`` stubs
    ``HZ.channel_names``: what is under test here is the name -> page lookup and the source
    precedence, not tifffile's series grouping. The parser itself is covered by
    :func:`test_qptiff_facts_reads_the_pixel_size_out_of_the_xml` and, against the real file
    format, by :func:`test_qptiff_facts_on_the_real_images`.
    """
    def facts(path):
        name = pathlib.Path(path).name
        return {"channel_names": INSTRUMENT_CHANNELS,
                "pixel_size_um": 0.5 if name == "tenx.qptiff" else 0.25,
                "pixel_size_detail": "qptiff_xml_PixelSizeMicrons"}
    return facts


@pytest.fixture
def resolve(project, inputs_csv, images, sections_csv, image_facts):
    def _resolve(serial, **kw):
        kw.setdefault("sections", sections_csv)
        kw.setdefault("facts_reader", image_facts)
        return PR.channel_roles(serial, project_dir=project, inputs_csv=inputs_csv,
                                images_dir=images, **kw)
    return _resolve


# ------------------------------------------------------------------------------- name, not index

def test_the_same_markers_resolve_to_different_pages_on_different_images(resolve):
    """The defence against index-based resolution, stated as an assertion.

    Measured on the real project 2026-09-18: INS/GCG sit at pages (2,3) on six entries, (1,2) on
    six, (1,3) on two and (3,1) on donor 6579. A reader that took a fixed index would return
    another antibody's channel on most of the cohort and raise nothing at all.
    """
    a = resolve("0037191")
    b = resolve("0059847")
    assert a.ok and b.ok
    assert a.channels == {"INS": 1, "GCG": 2}
    assert b.channels == {"INS": 2, "GCG": 3}
    assert a.channels != b.channels
    # The page carries the INSTRUMENT's name in the image; only QuPath knows the marker.
    assert a.qptiff_names[a.channels["INS"]] == "Opal 570"
    assert b.qptiff_names[b.channels["INS"]] == "Opal 690"
    assert a.names == {"INS": "INS", "GCG": "GCG"}
    assert a.source == b.source == "qupath_project"
    assert a.qupath_entry == "20"


def test_two_channels_claiming_one_marker_are_refused_never_picked(
        project, inputs_csv, images, sections_csv, image_facts):
    _entry(project, "20", "renamed_a.qptiff", ["DAPI", "INS", "GCG", "INS", "Sample AF"])
    with pytest.raises(PR.ProteinRnaError, match="ambiguous"):
        PR.channel_roles("0037191", project_dir=project, inputs_csv=inputs_csv,
                         images_dir=images, sections=sections_csv, facts_reader=image_facts)


# ------------------------------------------------------------------------- the named absences

def test_an_unrenamed_image_is_a_named_unavailable_not_an_exception(resolve):
    r = resolve("0037354")
    assert not r.ok
    assert r.unavailable == PR.CHANNELS_NOT_YET_ASSIGNED
    assert r.missing_markers == ("INS", "GCG")
    assert r.unrecognised_names == ()
    assert "untouched.qptiff" in r.detail and "Opal 690" in r.detail
    assert r.donor_id == "6482"


def test_a_transposed_name_is_its_own_state_and_is_not_silently_missed(resolve):
    """Entry 15 spells glucagon ``CGC``. 6380 passes its registration gate.

    ``channels_not_yet_assigned`` means nobody has been through that image;
    ``channel_name_unrecognised`` means somebody has and the spelling is not one this reader
    knows. The two need different actions, so they are different states.
    """
    r = resolve("0037502")
    assert not r.ok
    assert r.unavailable == PR.CHANNEL_NAME_UNRECOGNISED
    assert r.unrecognised_names == ("CGC",)
    assert "CGC" in r.detail
    assert r.missing_markers == ("INS", "GCG")


def test_a_partly_renamed_image_names_only_the_marker_that_is_missing(resolve):
    r = resolve("0037637")
    assert not r.ok
    assert r.unavailable == PR.CHANNELS_NOT_YET_ASSIGNED
    assert r.missing_markers == ("INS",)
    assert r.unrecognised_names == ()


def test_a_serial_the_inventory_never_heard_of_is_named_not_an_indexerror(resolve):
    r = resolve("9999999")
    assert not r.ok
    assert r.unavailable == "serial_not_in_registration_inputs"


def test_a_donor_disagreement_between_two_ledgers_is_refused(resolve):
    r = resolve("0037191", expect_donor="6442")
    assert not r.ok
    assert r.unavailable == "donor_disagrees_with_registration_inputs"
    assert "6442" in r.detail and "6493" in r.detail


def test_one_image_in_two_qupath_entries_is_refused_rather_than_picked(
        project, inputs_csv, images, sections_csv, image_facts):
    _entry(project, "21", "renamed_a.qptiff", ["DAPI", "INS", "GCG", "Opal 520", "Sample AF"])
    r = PR.channel_roles("0037191", project_dir=project, inputs_csv=inputs_csv,
                         images_dir=images, sections=sections_csv, facts_reader=image_facts)
    assert not r.ok
    assert r.unavailable == "duplicate_qupath_entries"


# ---------------------------------------------------------------------------- the CSV fallback

def test_the_csv_still_resolves_an_image_qupath_has_not_reached(resolve, project):
    """With the QuPath entry removed, the five-row CSV is what is left -- and it still works.

    It resolves through the image's OWN ``<Name>`` XML (``channel_opal690`` -> the page spelled
    ``Opal 690``), not through the column's position in the header.
    """
    (project / "data" / "5" / "server.json").unlink()
    r = resolve("0041116")
    assert r.ok
    assert r.source == "restain_sections.csv"
    assert r.channels == {"INS": 2, "GCG": 3}
    assert r.names == {"INS": "Opal 690", "GCG": "Opal 520"}
    # The image still supplies the scale: the CSV is a fallback for the ROLES, not for the grid.
    assert r.pixel_size_um == 0.5
    assert r.csv_conflict == ()


def test_qupath_wins_a_disagreement_and_the_disagreement_is_recorded(resolve):
    """The two 10x images are where the two sources genuinely contradict each other.

    Measured on the real project 2026-09-18: QuPath names ``Opal 570`` INS and ``Opal 690`` GCG
    on ``HDL115_Scan1.qptiff`` and ``npod6591_Scan1.qptiff``, while the CSV records ``Opal 690``
    INS, ``Opal 520`` GCG and ``Opal 570`` FAILED for all five of its rows. 6591 passes its
    registration gate. The newer per-image assignment is used; the contradiction is carried out
    by name rather than averaged away or quietly dropped.
    """
    r = resolve("0041116")
    assert r.ok
    assert r.source == "qupath_project"
    assert r.channels == {"INS": 1, "GCG": 2}
    assert len(r.csv_conflict) == 2
    assert any("INS" in c and "Opal 690" in c for c in r.csv_conflict)


# --------------------------------------------------------------------------------- the pixel size

def test_pixel_size_comes_from_the_image_not_the_csv(resolve):
    """The CSV declares 0.2505 for every row; the two 10x images are 0.5 um/px.

    Reading 6591's coordinates at 0.2505 selects the wrong pyramid level and halves every
    distance. The value used is recorded per row so nobody has to ask which it was.
    """
    twenty = resolve("0037191")
    ten = resolve("0041116")
    assert twenty.pixel_size_um == 0.25
    assert ten.pixel_size_um == 0.5
    assert twenty.pixel_size_source == ten.pixel_size_source == "qptiff_xml_PixelSizeMicrons"


def test_no_image_and_no_declared_size_is_unavailable_never_assumed(resolve, images):
    (images / "renamed_a.qptiff").unlink()
    r = resolve("0037191")
    assert not r.ok
    assert r.unavailable == "pixel_size_unavailable"


def test_a_declared_size_is_used_only_when_the_image_cannot_be_read(resolve, images):
    (images / "tenx.qptiff").unlink()
    r = resolve("0041116")
    assert r.ok
    assert r.pixel_size_um == 0.2505
    assert "restain_sections.csv" in r.pixel_size_source


@pytest.mark.parametrize("bad", ["", "0", "-0.25", "nan"])
def test_a_blank_or_impossible_declared_size_is_never_substituted(resolve, images, sections_csv,
                                                                  bad):
    (images / "tenx.qptiff").unlink()
    sections_csv.loc[sections_csv.xenium_serial == "0041116", "pixel_size_um"] = bad
    r = resolve("0041116", sections=sections_csv)
    assert not r.ok
    assert r.unavailable == "pixel_size_unavailable"


def test_qptiff_facts_reads_the_pixel_size_out_of_the_xml(tmp_path):
    """The parser itself, on a file it can be handed.

    The TIFF resolution rational is deliberately not consulted: both real values carry numerator
    4294967295, so the trailing digits of 0.25047 / 0.50139 are a uint32 clamp rather than a
    measurement.
    """
    import tifffile

    path = tmp_path / "one.qptiff"
    tifffile.imwrite(path, np.zeros((4, 4), dtype=np.uint8),
                     description=("<PerkinElmer-QPI-ImageDescription><Name>Opal 690</Name>"
                                  "<PixelSizeMicrons>0.5</PixelSizeMicrons>"
                                  "</PerkinElmer-QPI-ImageDescription>"),
                     resolution=(39924.958, 39924.958), resolutionunit="CENTIMETER")
    facts = PR.qptiff_facts(path)
    assert facts["pixel_size_um"] == 0.5
    assert facts["channel_names"] == ("Opal 690",)
    assert facts["pixel_size_detail"] == "qptiff_xml_PixelSizeMicrons"


def test_two_disagreeing_pixel_sizes_in_one_file_are_not_averaged(tmp_path, monkeypatch):
    class _Page:
        def __init__(self, description):
            self.description = description

    class _Level:
        pages = [_Page("<Name>DAPI</Name><PixelSizeMicrons>0.25</PixelSizeMicrons>"),
                 _Page("<Name>Opal 690</Name><PixelSizeMicrons>0.5</PixelSizeMicrons>")]

    class _Series:
        levels = [_Level()]

    class _Tif:
        series = [_Series()]

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    import tifffile
    monkeypatch.setattr(tifffile, "TiffFile", lambda *a, **k: _Tif())
    facts = PR.qptiff_facts(tmp_path / "x.qptiff")
    assert facts["pixel_size_um"] is None
    assert "found [0.25, 0.5]" in facts["pixel_size_detail"]


# ----------------------------------------------------------------------------------- the artifact

def test_resolve_all_never_raises_and_the_report_groups_by_named_reason(
        project, inputs_csv, images, sections_csv, image_facts, tmp_path):
    serials = {"0037191": "6493", "0059847": "6461", "0037354": "6482",
               "0037502": "6380", "0037637": "6457", "9999999": "6499"}
    resolved = PR.resolve_all(serials, project_dir=project, inputs_csv=inputs_csv,
                              images_dir=images, sections=sections_csv,
                              facts_reader=image_facts)
    assert set(resolved) == set(serials)
    report = PR.channel_assignment_report(resolved, project_dir=project)
    assert report["n_registered"] == 6
    assert report["n_resolved"] == 2
    assert report["donors_by_reason"][PR.CHANNELS_NOT_YET_ASSIGNED] == ["6457", "6482"]
    assert report["donors_by_reason"][PR.CHANNEL_NAME_UNRECOGNISED] == ["6380"]
    path = PR.write_channel_assignment(report, tmp_path)
    assert path == tmp_path / "protein_rna" / "channel_assignment.json"
    written = json.loads(path.read_text())
    assert {r["donor_id"] for r in written["unavailable"]} == {"6482", "6380", "6457", "6499"}
    assert all(r["pixel_size_source"] for r in written["resolved"])
    # Every unavailable row names the channels it saw, so the artifact is actionable on its own.
    for record in written["unavailable"]:
        assert record["unavailable"] and record["detail"]


def test_measure_one_refuses_by_name_instead_of_raising_indexerror(monkeypatch):
    monkeypatch.setattr(PR, "donor_of", lambda serial, paths: "6482")
    roles = PR.unavailable(PR.CHANNELS_NOT_YET_ASSIGNED, "0037354_Scan1.qptiff: not assigned",
                           serial="0037354", donor_id="6482",
                           image_file="0037354_Scan1.qptiff")
    with pytest.raises(PR.ProteinRnaError, match=PR.CHANNELS_NOT_YET_ASSIGNED):
        PR.measure_one("0037354", pd.DataFrame(), roles=roles)


# --------------------------------------------------------------- the real project, properties only

_REAL = pytest.mark.skipif(not (PR.QUPATH_PROJECT / "data").is_dir(),
                           reason="the QuPath project is not on this machine")


@_REAL
def test_the_real_project_resolves_by_name_and_never_by_position():
    """Properties that stay true as the maintainer keeps renaming, not a count that moves.

    Asserting "eleven of thirteen resolve today" would make this test a record of one afternoon
    of QuPath work. What is asserted instead is that every resolution came from a NAME, that no
    two markers ever land on one page, and that every absence carries one of the named states.
    """
    assert PR.qupath_entries(), "the project holds no readable server.json"
    serials = PR.registered_serials("panc")
    assert serials, "no donor passes the registration gate"
    resolved = PR.resolve_all(serials, sections=None)
    named = {PR.CHANNELS_NOT_YET_ASSIGNED, PR.CHANNEL_NAME_UNRECOGNISED,
             "serial_not_in_registration_inputs", "donor_disagrees_with_registration_inputs",
             "duplicate_qupath_entries", "pixel_size_unavailable", "channel_resolution_failed"}
    for serial, roles in resolved.items():
        if not roles.ok:
            assert roles.unavailable in named, f"{serial}: {roles.unavailable}"
            assert roles.detail
            continue
        assert set(roles.channels) == set(PR.MARKERS)
        assert len(set(roles.channels.values())) == 2, f"{serial}: two markers on one page"
        for marker, page in roles.channels.items():
            assert 0 <= page < len(roles.qptiff_names)
            assert PR._norm(roles.names[marker]) in PR.MARKER_ALIASES[marker]
        assert roles.pixel_size_um in (0.25, 0.5), roles.pixel_size_um
        assert roles.pixel_size_source == "qptiff_xml_PixelSizeMicrons"


@_REAL
def test_an_unrecognised_name_in_the_real_project_is_reported_as_such():
    """If any entry carries a rename this reader does not know, it is named, not missed."""
    for roles in PR.resolve_all(PR.registered_serials("panc"), sections=None).values():
        if roles.ok or not roles.unrecognised_names:
            continue
        assert roles.unavailable == PR.CHANNEL_NAME_UNRECOGNISED
        assert all(name in roles.detail for name in roles.unrecognised_names)


@pytest.mark.skipif(not PR.IMAGES.is_dir(), reason="the IF images are not on this machine")
def test_qptiff_facts_on_the_real_images():
    """The parser against the real file format, on one image of each objective.

    Both magnifications are read because the two groups are the whole point: the 20x images are
    0.25 um/px and the 10x ones 0.5, while ``restain_sections.csv`` declares 0.2505 for both.
    """
    for name, expected in (("0059865_Scan1.qptiff", 0.25), ("npod6591_Scan1.qptiff", 0.5)):
        path = PR.IMAGES / name
        if not path.is_file():
            pytest.skip(f"{name} is not on this machine")
        facts = PR.qptiff_facts(path)
        assert facts["pixel_size_um"] == expected, name
        # The image's XML carries the INSTRUMENT's names and never the marker's.
        assert all(PR.INSTRUMENT_NAME.match(PR._norm(n)) for n in facts["channel_names"]), \
            facts["channel_names"]
