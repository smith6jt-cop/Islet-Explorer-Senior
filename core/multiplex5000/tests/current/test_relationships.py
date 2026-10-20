"""Recorded relationships: read as written, stopped where the record stops, never inferred.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_relationships.py -q

Synthetic sources hold the reader to the notes' grammar, to stopping after the post-Xenium IF
section (a planted line after it must never surface), to refusing what it cannot read, and to
reporting -- not resolving -- every disagreement between the sources. The one data test reads the
real sources this site names.
"""

from __future__ import annotations

import hashlib
import textwrap

import pytest

from multiplex5000 import relationships as R

ALIASES = {"HDL112": "0112", "112": "0112", "HDL115": "115"}
PLANTED = "PLANTED-LINE-AFTER-THE-IF-SECTION do not read"

NOTES = textwrap.dedent(f"""\
    IO60panc2nd Notes
    6374: pancTail no pLN, sequential Xenium
    6539: pancHead and 1 pLN same-block, huge vessel between two pieces of panc, sequential Xenium
    HDL115: pancTail and 1 pLN different-block, sequential Xenium
    6457: NEEDS TO BE RERUN OR REPROCESSED; HAS 3 pLNs same-block and Tissue-embedded TLS?, sequential Xenium
    6579: pancHead and 2 pLN same-block; bulk T and B aggregates, but kept in Tissue, sequential Xenium
    HDL112: pancTail and 1 pLN different-block, area of pLN is about the same as panc, same-slide post-phenocycler Xenium
    6347: pancHead and 1 pLN different-block and 1 pLN same-block, same-slide post-phenocycler Xenium


    Post-Xenium IF Notes
    6457: post-Xenium IF too much punctate background, not exported, skipped

    {PLANTED}
    this line is not a donor line and must not be parsed
    """)

PLN_BLOCKS = textwrap.dedent("""\
    image,entry_id,donor,pln_annotation,object_id,pieces,area_mm2,overlap_with_tissue_um2,gap_to_tissue_um,measured_relation,block_from_notes,notes
    img6539,1,6539,1,uuid-6539,1,1.0,0.0,0.0,abuts Tissue (gap 0),same-block,x
    imgHDL115,2,HDL115,1,uuid-115,1,1.0,0.0,627.7,separate,different-block,x
    img6579,3,6579,1,uuid-6579,2,1.0,0.0,0.0,abuts Tissue (gap 0),same-block,x
    imgHDL112,4,HDL112,1,uuid-0112,1,1.0,0.0,664.7,separate,different-block,x
    img6347,5,6347,1,uuid-6347a,1,1.0,0.0,4409.2,separate,different-block,x
    img6347,5,6347,2,uuid-6347b,1,1.0,0.0,0.0,abuts Tissue (gap 0),same-block,x
    """)

MODES = textwrap.dedent("""\
    # comments are allowed
    donor_id,roi,xenium_serial,mode,evidence,note
    112,panc,,sequential,unqualified callers get the serial pair,
    112,panc,0041396,sequential,serial section from the same block,
    112,panc,0081366,same_slide,same slide after PhenoCycler,
    6347,panc,0080579,same_slide,same slide after PhenoCycler,
    """)


@pytest.fixture
def sources(tmp_path, monkeypatch):
    """Synthetic sources only: the site's real ones are hidden, so a source left out is unnamed."""
    monkeypatch.setattr(R.roots, "source", lambda name: None)
    paths = {}
    for name, text in (("notes.txt", NOTES), ("pln_blocks.csv", PLN_BLOCKS),
                       ("section_modes.csv", MODES)):
        paths[name] = tmp_path / name
        paths[name].write_text(text)
    return paths


def _load(sources, **over):
    kw = dict(notes=sources["notes.txt"], notes_copy=None, pln_blocks=sources["pln_blocks.csv"],
              section_modes=sources["section_modes.csv"], aliases=ALIASES,
              cohort_donors=["6374", "6539", "115", "0112", "6347", "6579"])
    kw.update(over)
    return R.load(**kw)


def test_each_donor_line_is_read_as_written(sources):
    donors, _ = R.read_donor_notes(sources["notes.txt"], ALIASES)
    assert set(donors) == {"6374", "6539", "115", "6457", "6579", "0112", "6347"}
    assert (donors["6374"].pancreas_region, donors["6374"].same_block_pln,
            donors["6374"].different_block_pln) == ("tail", 0, 0)
    assert donors["6539"].cautions == ("huge vessel between two pieces of panc",)
    assert (donors["115"].label, donors["115"].different_block_pln) == ("HDL115", 1)
    assert (donors["6347"].same_block_pln, donors["6347"].different_block_pln) == (1, 1)
    assert donors["6347"].xenium_design == R.SAME_SECTION_AFTER_PHENOCYCLER
    assert donors["6579"].cautions == ("bulk T and B aggregates, but kept in Tissue",)
    assert donors["0112"].cautions == ("area of pLN is about the same as panc",)


def test_a_count_in_free_text_is_kept_as_stated_and_flagged(sources):
    donors, _ = R.read_donor_notes(sources["notes.txt"], ALIASES)
    d = donors["6457"]
    assert d.pancreas_region is None and d.pln_from_caution
    assert (d.same_block_pln, d.different_block_pln) == (3, None)      # unstated stays unknown
    assert d.cautions == ("NEEDS TO BE RERUN OR REPROCESSED",
                          "HAS 3 pLNs same-block and Tissue-embedded TLS?")
    assert d.post_xenium_if == "too much punctate background, not exported, skipped"


def test_reading_stops_after_the_if_section(sources):
    donors, evidence = R.read_donor_notes(sources["notes.txt"], ALIASES)
    # the header, seven donor lines, two blank lines, the IF header and its one line
    consumed = "".join(NOTES.splitlines(keepends=True)[:12])
    assert evidence.where == "lines 1-12"
    assert evidence.sha256 == hashlib.sha256(consumed.encode()).hexdigest()
    report = R.report(_load(sources))
    assert PLANTED not in report and "must not be parsed" not in report
    assert all(PLANTED not in repr(d) for d in donors.values())


@pytest.mark.parametrize("bad, message", [
    ("6374 pancTail no pLN, sequential Xenium", "not a donor line"),
    ("6374: pancTail no pLN, Xenium later", "does not end with a Xenium design"),
    ("6374: pancTail no pLN-ish, sequential Xenium", "cannot read"),
])
def test_a_line_it_cannot_read_is_an_error(tmp_path, bad, message):
    path = tmp_path / "notes.txt"
    path.write_text(f"IO60panc2nd Notes\n{bad}\n")
    with pytest.raises(R.RelationshipError, match=message):
        R.read_donor_notes(path, ALIASES)


def test_two_labels_for_one_donor_or_a_stray_if_note_are_refused(tmp_path):
    path = tmp_path / "notes.txt"
    path.write_text("IO60panc2nd Notes\nHDL112: pancTail no pLN, sequential Xenium\n"
                    "112: pancTail no pLN, sequential Xenium\n")
    with pytest.raises(R.RelationshipError, match="same donor"):
        R.read_donor_notes(path, ALIASES)
    path.write_text("IO60panc2nd Notes\n6374: pancTail no pLN, sequential Xenium\n\n"
                    "Post-Xenium IF Notes\n9999: post-Xenium IF skipped\n")
    with pytest.raises(R.RelationshipError, match="no donor line"):
        R.read_donor_notes(path, ALIASES)
    path.write_text("Some other file\n")
    with pytest.raises(R.RelationshipError, match="header"):
        R.read_donor_notes(path, ALIASES)


@pytest.mark.parametrize("edit, message", [
    (lambda t: t.replace("block_from_notes", "block"), "expected the columns"),
    (lambda t: t.replace("same-block,x\n", "adjacent,x\n", 1), "not same-block or different-block"),
    (lambda t: t.replace("uuid-115", "uuid-6539"), "empty or repeated"),
    (lambda t: t.replace(",uuid-6539,1,", ",uuid-6539,0,"), "at least 1"),
])
def test_pln_blocks_is_read_strictly(sources, edit, message):
    path = sources["pln_blocks.csv"]
    path.write_text(edit(PLN_BLOCKS))
    with pytest.raises(R.RelationshipError, match=message):
        R.read_pln_blocks(path, ALIASES)


def test_the_relationships_come_from_the_records(sources):
    reading = _load(sources)
    kinds = {(r.kind, r.donor_id, r.b) for r in reading.relationships}
    assert (R.SAME_SECTION, "0112", "Xenium section 0081366") in kinds         # section_modes
    assert (R.SERIAL_SECTION, "0112", "Xenium section 0041396") in kinds
    assert (R.SERIAL_SECTION, "6374", "Xenium section (the notes name no serial)") in kinds
    blocks = sorted((r.kind, r.donor_id) for r in reading.relationships
                    if r.kind in (R.SAME_BLOCK, R.DIFFERENT_BLOCK))
    assert blocks == [(R.DIFFERENT_BLOCK, "0112"), (R.DIFFERENT_BLOCK, "115"),
                      (R.DIFFERENT_BLOCK, "6347"), (R.SAME_BLOCK, "6347"),
                      (R.SAME_BLOCK, "6539"), (R.SAME_BLOCK, "6579")]
    assert sum(r.kind == R.SHARED_SLIDE for r in reading.relationships) == 6
    assert all(r.evidence for r in reading.relationships if r.kind != R.UNKNOWN)


def test_disagreements_are_reported_and_nothing_is_chosen(sources):
    reading = _load(sources)
    conflicts = {(c.donor_id, c.what) for c in reading.conflicts}
    # the notes give 6457 three same-block nodes; pln_blocks has none
    assert ("6457", "same_block pLN: the notes and pln_blocks.csv") in conflicts
    assert len(reading.conflicts) == 1
    sources["section_modes.csv"].write_text(MODES + "6374,panc,0041000,same_slide,said so,\n")
    reading = _load(sources)
    assert ("6374", "Xenium design: the notes and section_modes.csv") in {
        (c.donor_id, c.what) for c in reading.conflicts}


def test_the_two_copies_of_the_notes_must_agree(sources, tmp_path):
    copy = tmp_path / "copy.md"
    copy.write_text(NOTES.replace("huge vessel", "large vessel"))
    reading = _load(sources, notes_copy=copy)
    assert ("6539", "the two copies of the donor notes differ") in {
        (c.donor_id, c.what) for c in reading.conflicts}


def test_what_nobody_recorded_is_unknown_with_the_reason(sources):
    reading = _load(sources, cohort_donors=["6374", "7000"])
    unknown = [r for r in reading.relationships if r.kind == R.UNKNOWN]
    assert [(r.donor_id, r.reason) for r in unknown] == [
        ("7000", "the cohort lists this donor; the donor notes have no line for it")]
    reading = _load(sources, notes=None, pln_blocks=None)
    assert {c.what for c in reading.conflicts} >= {"donor_notes: not named by this site",
                                                   "pln_blocks: not named by this site"}


#: The known disagreements are kept once, beside the reader; a new one fails here until someone
#: has looked at it.
KNOWN = R.KNOWN_DISAGREEMENTS


@pytest.mark.data
def test_the_real_sources_agree_except_where_known():
    reading = R.load()
    assert len(reading.donors) == 26, "the donor notes describe 26 donors"
    assert {e.source for e in reading.sources} >= {"donor_notes", "donor_notes_copy", "pln_blocks",
                                                   "section_modes"}
    found = {(c.donor_id, c.what) for c in reading.conflicts}
    assert found - set(KNOWN) == set(), sorted(found - set(KNOWN))
