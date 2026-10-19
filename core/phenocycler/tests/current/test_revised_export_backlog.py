"""Which revised QuPath exports still have to be rebuilt, and who each one belongs to.

`data/RevisedQuPathExports` is the first place a donor's cell matrix is looked for, and the
maintainer drops new exports into it while the lane runs.  Two things then have to hold.

The backlog has to be *discovered*.  A tuple of donor ids in the source goes stale the moment an
export lands -- that is exactly how `cohort._rebuilt()`'s literal list went on refusing donors
whose re-export had already arrived.

And the donor has to come from the export's own `Image` value.  The filenames in that directory
do not agree with the donor ids: `0112.tsv` and `HDL112_with_pLN.tsv` are one image, `HDL115_with_pLN.tsv`
is donor `115`, and `6533.tsv` and `6533_with_pLN.tsv` are one image too.  Parsing any of those
spellings is donor-name parsing, which this repo forbids because it is how a right-looking number
ends up describing the wrong donor.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

TYPING = Path(__file__).resolve().parents[2] / "scripts/typing"
sys.path.insert(0, str(TYPING))

import rebuild_donor_from_revised_export as R  # noqa: E402
import subpanel_io as S  # noqa: E402

DIST = "Signed distance to annotation with Islet µm"
HEADER = ["Image", "Object ID", "Object type", "Name", "Classification", "Parent",
          "Centroid X µm", "Centroid Y µm", "Nucleus: Area µm^2", "Cell: Area µm^2",
          "Cell: Solidity", "Nucleus/Cell area ratio", DIST]


def write_export(path: Path, image: str, parents: list[str], *, mtime: float) -> Path:
    """A minimal export: one image, one row per parent, all geometry inside the thresholds."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(HEADER)]
    for index, parent in enumerate(parents):
        lines.append("\t".join([
            image, f"uuid-{index:03d}", "Cell", "", "", parent,
            str(10.0 + index), str(20.0 + index), "35.0", "120.0", "0.95", "0.29",
            str(-5.0 + index * 12.0),
        ]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    import os
    os.utime(path, (mtime, mtime))
    return path


@pytest.fixture()
def lane(tmp_path, monkeypatch):
    revised = tmp_path / "RevisedQuPathExports"
    out = tmp_path / "revised_cells"
    aligned = tmp_path / "uniform"
    revised.mkdir()

    manifest = tmp_path / "qupath_manifest.json"
    manifest.write_text(json.dumps({"images": [
        {"donor_id": "0112", "image_id": "0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1"},
        {"donor_id": "115", "image_id": "HDL115pancLN_Scan1.er.qptiff - resolution #1"},
        {"donor_id": "6374", "image_id": "6374_Scan1.er.qptiff - resolution #1"},
        {"donor_id": "6457", "image_id": "6457_Scan1.er.qptiff - resolution #1"},
    ]}))
    acquisition = tmp_path / "donor_acquisition_runs.csv"
    acquisition.write_text("donor_id,image_file,excluded\n"
                           "6579,6579_Scan2.er.qptiff,True\n"
                           "6374,6374_Scan1.er.qptiff,False\n")

    geometry_qc = tmp_path / "geometry_qc"
    for donor in ("0112", "115", "6374", "6457", "6579"):
        target = aligned / f"donor_id={donor}"
        target.mkdir(parents=True)
        ids = [f"uuid-{i:03d}" for i in range(12)]
        pd.DataFrame({"object_id": ids, "CD3e": [1.0] * 12}).to_parquet(target / "data.parquet", index=False)
        # The source workflow's QC decisions a rebuild keeps: every synthetic cell passes.
        qc = geometry_qc / f"donor_id={donor}"
        qc.mkdir(parents=True)
        pd.DataFrame({"object_id": ids, "analysis_eligible": [True] * 12, "estimation_eligible": [True] * 12,
                      "analysis_reason": [""] * 12, "estimation_reason": [""] * 12}).to_parquet(
            qc / "data.parquet", index=False)

    monkeypatch.setattr(R, "GEOMETRY_QC", geometry_qc)
    monkeypatch.setattr(R, "REVISED", revised)
    monkeypatch.setattr(R, "OUT", out)
    monkeypatch.setattr(R, "SECTION_PARENTS", out / "section_parents.json")
    monkeypatch.setattr(R, "MANIFEST", manifest)
    monkeypatch.setattr(R, "ACQUISITION", acquisition)
    monkeypatch.setattr(S, "ALIGNED", aligned)
    return revised, out, aligned


def state_of(records: list[dict], tsv: str) -> dict:
    return next(r for r in records if r["tsv"] == tsv)


def test_the_donor_comes_from_the_image_not_the_file_name(lane):
    """A file whose name says nothing still resolves, and one whose name lies does not mislead."""
    revised, _, _ = lane
    write_export(revised / "anything_at_all.tsv",
                 "6374_Scan1.er.qptiff - resolution #1", ["Annotation (Tissue)"], mtime=2_000_000)
    write_export(revised / "HDL112_with_pLN.tsv",
                 "0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    records = R.survey()
    assert state_of(records, "anything_at_all.tsv")["donor"] == "6374"
    assert state_of(records, "HDL112_with_pLN.tsv")["donor"] == "0112"


def test_both_manifests_are_read(lane):
    """Neither covers the directory alone: 0112 is only in one, 6579 only in the other."""
    revised, _, _ = lane
    write_export(revised / "6579.tsv", "6579_Scan2.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    assert state_of(R.survey(), "6579.tsv")["donor"] == "6579"


def test_an_image_neither_manifest_names_is_refused_not_guessed(lane):
    revised, _, _ = lane
    write_export(revised / "6999.tsv", "6999_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    record = state_of(R.survey(), "6999.tsv")
    assert record["state"] == "refused"
    assert record["donor"] is None
    assert "6999_Scan1.er.qptiff" in record["reason"]


def test_two_manifests_disagreeing_about_one_image_raises(lane, tmp_path):
    acquisition = tmp_path / "donor_acquisition_runs.csv"
    acquisition.write_text("donor_id,image_file,excluded\n"
                           "9999,6374_Scan1.er.qptiff,False\n")
    with pytest.raises(ValueError, match="refusing to pick"):
        R.image_to_donor()


def test_the_newest_export_for_a_donor_supersedes_the_older_one(lane):
    revised, _, _ = lane
    write_export(revised / "0112.tsv",
                 "0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=1_000_000)
    write_export(revised / "HDL112_with_pLN.tsv",
                 "0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1",
                 ["Annotation (Tissue)", "Annotation (pLN)"], mtime=2_000_000)
    records = R.survey()
    assert state_of(records, "0112.tsv")["state"] == "superseded"
    assert "HDL112_with_pLN.tsv" in state_of(records, "0112.tsv")["reason"]
    assert state_of(records, "HDL112_with_pLN.tsv")["state"] == "backlog"


def test_a_centrally_excluded_donor_is_refused(lane):
    """The exclusion list is asked, never restated -- 6457's export is not rebuilt."""
    revised, _, _ = lane
    write_export(revised / "6457.tsv", "6457_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    record = state_of(R.survey(), "6457.tsv")
    assert record["state"] == "refused"
    assert "donor excluded" in record["reason"]


def test_a_donor_with_no_aligned_partition_is_refused_with_that_named(lane):
    """Joining on object_id needs rows on both sides; there is nothing to impute."""
    revised, _, aligned = lane
    import shutil
    shutil.rmtree(aligned / "donor_id=6374")
    write_export(revised / "6374.tsv", "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    record = state_of(R.survey(), "6374.tsv")
    assert record["state"] == "refused"
    assert "no aligned partition" in record["reason"]


def test_freshness_is_the_export_against_the_partition_not_a_list(lane):
    """Both directions, because the whole point is that neither is remembered anywhere."""
    import os
    revised, out, _ = lane
    write_export(revised / "6374.tsv", "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=1_000_000)
    target = out / "donor_id=6374"
    target.mkdir(parents=True)
    pd.DataFrame({"object_id": ["uuid-000"]}).to_parquet(target / "data.parquet", index=False)
    os.utime(target / "data.parquet", (2_000_000, 2_000_000))
    assert state_of(R.survey(), "6374.tsv")["state"] == "fresh"

    os.utime(revised / "6374.tsv", (3_000_000, 3_000_000))
    assert state_of(R.survey(), "6374.tsv")["state"] == "backlog"


def test_list_mode_writes_nothing(lane, monkeypatch, capsys):
    revised, out, _ = lane
    write_export(revised / "6374.tsv", "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    monkeypatch.setattr(sys, "argv", ["rebuild", "--list"])
    assert R.main() == 0
    printed = capsys.readouterr().out
    assert "6374.tsv" in printed and "backlog" in printed
    assert not (out / "donor_id=6374").exists()
    assert not (out / "section_parents.json").exists()


def test_naming_a_refused_donor_stops_the_run(lane, monkeypatch):
    revised, _, _ = lane
    write_export(revised / "6457.tsv", "6457_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    monkeypatch.setattr(sys, "argv", ["rebuild", "6457"])
    with pytest.raises(SystemExit, match="6457"):
        R.main()


def test_the_second_tissue_does_not_enter_the_pancreas_partition(lane):
    """A cell parented to a pLN annotation has no islet to be at a distance from."""
    revised, out, _ = lane
    parents = ["Annotation (Tissue)", "Islet_7", "Annotation (pLN)", "Annotation (pLN)",
               "Root object (Image)"]
    write_export(revised / "HDL115_with_pLN.tsv",
                 "HDL115pancLN_Scan1.er.qptiff - resolution #1", parents, mtime=2_000_000)
    record = state_of(R.survey(), "HDL115_with_pLN.tsv")
    R.rebuild(record, R.load_thresholds())

    written = pd.read_parquet(out / "donor_id=115/data.parquet")
    assert len(written) == 3
    assert not written["parent_raw"].str.contains("pLN").any()

    census = json.loads((out / "section_parents.json").read_text())["donors"]["115"]
    assert census["annotation_parent_counts"] == {"Tissue": 1, "pLN": 2}
    assert census["pln_annotation_classes"] == ["pLN"]
    assert census["islet_parent_rows"] == 1
    assert census["root_object_rows"] == 1
    assert census["rows_in_export"] == 5
    assert census["rows_written_to_pancreas_partition"] == 3
    assert census["image"] == "HDL115pancLN_Scan1.er.qptiff - resolution #1"


def test_an_image_with_no_pln_annotation_records_a_measured_zero(lane):
    """Absence of evidence is recorded as an empty list, not as a missing key."""
    revised, out, _ = lane
    write_export(revised / "6374.tsv", "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)", "Islet_1"], mtime=2_000_000)
    R.rebuild(state_of(R.survey(), "6374.tsv"), R.load_thresholds())
    census = json.loads((out / "section_parents.json").read_text())["donors"]["6374"]
    assert census["pln_annotation_classes"] == []
    assert census["annotation_parent_counts"] == {"Tissue": 1}


def test_the_census_accumulates_across_donors(lane):
    revised, out, _ = lane
    write_export(revised / "6374.tsv", "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    write_export(revised / "HDL115_with_pLN.tsv",
                 "HDL115pancLN_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)", "Annotation (pLN)"], mtime=2_000_000)
    thresholds = R.load_thresholds()
    for record in R.survey():
        if record["state"] == "backlog":
            R.rebuild(record, thresholds)
    assert set(json.loads((out / "section_parents.json").read_text())["donors"]) == {"6374", "115"}


def test_an_export_carrying_two_images_is_refused(lane):
    """One image rebuilds into one donor partition; splitting them is not this lane's call."""
    revised, _, _ = lane
    path = revised / "mixed.tsv"
    write_export(path, "6374_Scan1.er.qptiff - resolution #1",
                 ["Annotation (Tissue)"], mtime=2_000_000)
    with path.open("a", encoding="utf-8") as handle:
        handle.write("\t".join([
            "HDL115pancLN_Scan1.er.qptiff - resolution #1", "uuid-900", "Cell", "", "",
            "Annotation (Tissue)", "1.0", "2.0", "35.0", "120.0", "0.95", "0.29", "1.0",
        ]) + "\n")
    record = state_of(R.survey(), "mixed.tsv")
    with pytest.raises(ValueError, match="distinct Image values"):
        R.rebuild(record, R.load_thresholds())


def test_classify_parents_names_every_class_it_saw(lane):
    """A class it does not recognise is counted under its own name, never folded into pancreas."""
    is_pln, census = R.classify_parents(pd.Series([
        "Annotation (Tissue)", "Annotation (pLN)", "Annotation (Spleen)", "Islet_12",
        "Root object (Image)", "",
    ]))
    assert list(is_pln) == [False, True, False, False, False, False]
    assert census["annotation_parent_counts"] == {"Spleen": 1, "Tissue": 1, "pLN": 1}
    assert census["pln_annotation_classes"] == ["pLN"]
    assert census["islet_parent_rows"] == 1
    assert census["root_object_rows"] == 1
    assert census["other_parent_rows"] == 1
