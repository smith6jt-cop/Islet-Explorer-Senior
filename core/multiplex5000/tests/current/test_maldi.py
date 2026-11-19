"""MALDI imaging as a sequential-slide platform (2026-10-31; docs/PLAN_2026-10-31_MALDI_SEQUENTIAL.md).

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_maldi.py -q

On a small synthetic export written the way the real one is (two files over the same spots, an
unnamed row-number column, ``NA`` for a missing value): the platform and its panel name, the
Source Inventory's MALDI rows and inspection gate, MALDI Ingest's keys and zeros, the serial
registration recovering a known rotation, mirror and shift, and the object rules of Serial
Integration.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pytest
import yaml

from multiplex5000 import gates, maldi, panel_names, profiles, stages, tasks
from multiplex5000 import project as P
from multiplex5000.cohort import PLATFORMS, excluded_on

GLYCANS = ("Hex5HexNAc2..Na", "Hex6HexNAc2..Na")
GLYCOGEN = ("X689", "X851")
INSPECTION = "tissue_id,tissue,neighbour_acquisition,inspected_by,note\n"


def _disk(cx: float, cy: float, radius_spots: int) -> list[tuple[float, float]]:
    out = []
    for i in range(-radius_spots, radius_spots + 1):
        for j in range(-radius_spots, radius_spots + 1):
            if i * i + j * j <= radius_spots * radius_spots:
                out.append((cx + 50.0 * i, cy + 50.0 * j))
    return out


def write_export(folder: Path, tissues: dict[str, list[tuple[float, float]]], *,
                 zero_spot: bool = True, glycogen_shift: float = 0.0) -> Path:
    """The two export files, as R writes them, and their SHA256SUMS."""
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    spot = 0
    for raster, (tissue, points) in enumerate(tissues.items()):
        for k, (x, y) in enumerate(points):
            rows.append((spot, raster, x, y, tissue, k))
            spot += 1
    for name, features, shift in (("syn_glycan_pbp.csv", GLYCANS, 0.0),
                                  ("syn_glycogen_pbp.csv", GLYCOGEN, glycogen_shift)):
        with (folder / name).open("w", newline="") as handle:
            writer = csv.writer(handle, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
            writer.writerow(["", "spotId", "raster", "x", "y", "z", "Tissue_id", "Date", *features])
            for n, (spot_id, raster, x, y, tissue, k) in enumerate(rows, start=1):
                values = [0.0 if (zero_spot and k == 0) else float(10 + (k % 7) * i)
                          for i in range(len(features))]
                if features == GLYCANS and k == 1:
                    values[0] = 0.0          # one feature not detected in one spot: kept as 0
                handle.write(",".join([f'"{n}"', str(spot_id), str(raster), repr(x + shift),
                                       repr(y), "0", f'"{tissue}"', "NA",
                                       *("NA" if (features == GLYCOGEN and k == 2 and i == 0)
                                         else repr(v) for i, v in enumerate(values))]) + "\n")
    sums = "".join(f"{hashlib.sha256((folder / n).read_bytes()).hexdigest()}  {n}\n"
                   for n in ("syn_glycan_pbp.csv", "syn_glycogen_pbp.csv"))
    (folder / maldi.CHECKSUMS).write_text(sums)
    return folder


# ------------------------------------------------------------------------------ the platform ----
def test_maldi_is_a_platform_and_6457s_exclusion_does_not_cover_it():
    assert "maldi" in PLATFORMS
    assert excluded_on("6457", "phenocycler") and not excluded_on("6457", "maldi")
    with pytest.raises(ValueError):
        excluded_on("6457", "MALDI")


def test_the_panel_is_named_from_the_feature_sets():
    rule = profiles.load_all()[("platform", "maldi")].panel_names
    assert panel_names.maldi({"glycan": 38, "glycogen": 10}, rule) == "38 N-glycans + 10 glycogen"
    with pytest.raises(panel_names.PanelNameError):
        panel_names.maldi({"glycan": 38, "lipid": 4, "glycogen": 10}, rule)
    with pytest.raises(panel_names.PanelNameError):
        panel_names.maldi({"glycan": 38, "glycogen": 0}, rule)


def test_a_tissue_id_gives_its_donor_and_nothing_else():
    aliases = {"HDL115": "115"}
    assert maldi.tissue_donor("HDL115", aliases) == ("HDL115", "115", "")
    assert maldi.tissue_donor("6551_03B", aliases) == ("6551", "6551", "03B")
    assert maldi.tissue_donor("6461_02_PLN", aliases) == ("6461", "6461", "02_PLN")


def test_the_two_files_must_hold_the_same_spots(tmp_path):
    good = write_export(tmp_path / "good", {"6374_04": _disk(0, 0, 3)})
    con = maldi.connect()
    files = maldi.source_files(good)
    assert not any(v for k, v in maldi.pairing_problems(con, files).items()
                   if not k.endswith("_rows"))
    bad = write_export(tmp_path / "bad", {"6374_04": _disk(0, 0, 3)}, glycogen_shift=1.0)
    assert maldi.pairing_problems(con, maldi.source_files(bad))["spots_that_differ"] > 0


def test_a_source_folder_without_its_checksums_is_refused(tmp_path):
    folder = write_export(tmp_path / "m", {"6374_04": _disk(0, 0, 2)})
    (folder / maldi.CHECKSUMS).unlink()
    with pytest.raises(maldi.MaldiSourceError):
        maldi.source_files(folder)


def test_stage_plan_and_routes():
    order = [s for s, _ in stages.PLAN]
    assert order.index("04_rna_ingest") < order.index("04b_maldi_ingest") < order.index("05_normalization")
    project = P.load("io60_maldi_serial")
    assert project.stages == ("00_inventory", "04b_maldi_ingest", "07_registration", "08_integration")
    registration, integration = stages.get("07_registration"), stages.get("08_integration")
    assert registration._route(project).__name__.endswith("maldi_serial_registration")
    assert integration._route(project).__name__.endswith("maldi_serial_integration")
    assert registration._route(P.load("spleen_pilot")).__name__.endswith("same_slide_registration")
    assert integration._route(P.load("spleen_pilot")).__name__.endswith("same_slide_integration")


# ------------------------------------------------------------------------- the inventory ----
def _with_maldi(path: Path, export: Path, inspection: str) -> Path:
    data = yaml.safe_load(path.read_text())
    (path.parent / "maldi_tissue_inspection.csv").write_text(inspection)
    data["settings"]["maldi"] = True
    data["sources"]["maldi_spots"] = str(export)
    data["sources"]["maldi_tissue_inspection"] = str(path.parent / "maldi_tissue_inspection.csv")
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


TISSUES = {"6374_04": _disk(10_000, 10_000, 4), "6374_01_pln": _disk(20_000, 10_000, 2),
           "HPAP_001": _disk(30_000, 10_000, 2)}
FILLED = (INSPECTION + "6374_04,pancreas,phenocycler:6374_Scan1 - resolution #1,tester,\n"
          "6374_01_pln,pLN,phenocycler:6374pLN_Scan1.qptiff,tester,\n")


def _inventory(path):
    project = P.load(path)
    result = tasks.run_stage(stages.get("00_inventory"), project)
    checks = {c["name"]: c for c in gates.load(result.path)["checks"]}
    return project, result.path, checks


def _rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def test_the_inventory_lists_io60_maldi_sections_and_blocks_until_inspected(
        synthetic_project, tmp_path):
    export = write_export(tmp_path / "maldi", TISSUES)
    _, folder, checks = _inventory(_with_maldi(synthetic_project(), export, INSPECTION))
    rows = {r["acquisition_id"]: r for r in _rows(folder / "maldi_acquisitions.csv")}
    assert set(rows) == {"maldi:6374_04", "maldi:6374_01_pln"}, "HPAP is out of cohort"
    assert rows["maldi:6374_04"]["panel"] == "2 N-glycans + 2 glycogen"
    assert rows["maldi:6374_04"]["id_block_token"] == "04"
    assert [r["tissue_id"] for r in _rows(folder / "maldi_out_of_cohort.csv")] == ["HPAP_001"]
    acquisitions = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}
    assert acquisitions["maldi:6374_04"]["platform"] == "maldi"
    not_inspected = checks["every IO60 MALDI section has been inspected: its tissue and the "
                           "slide it neighbours"]
    assert not not_inspected["passed"] and not_inspected["blocking"]
    assert (folder / "gate" / "maldi_sections.png").stat().st_size > 1000


def test_an_inspected_inventory_passes_and_records_serial_sections(synthetic_project, tmp_path):
    export = write_export(tmp_path / "maldi", TISSUES)
    _, folder, checks = _inventory(_with_maldi(synthetic_project(), export, FILLED))
    assert all(c["passed"] for c in checks.values() if c["blocking"])
    serial = [r for r in _rows(folder / "relationships.csv")
              if r["a"].startswith("MALDI section")]
    assert {r["b"] for r in serial} == {"phenocycler:6374_Scan1 - resolution #1",
                                        "phenocycler:6374pLN_Scan1.qptiff"}
    assert {r["kind"] for r in serial} == {"serial_section"}
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert coverage["6374"]["MALDI section: MALDI"] == "2"


def test_an_inspection_row_must_name_the_same_donors_acquisition(synthetic_project, tmp_path):
    export = write_export(tmp_path / "maldi", TISSUES)
    wrong = (INSPECTION + "6374_04,pancreas,phenocycler:HDL112_Scan3 - resolution #1,tester,\n"
             "6374_01_pln,lymphoid,phenocycler:6374pLN_Scan1.qptiff,tester,\n"
             "HPAP_001,pancreas,phenocycler:6374_Scan1 - resolution #1,,\n")
    _, _, checks = _inventory(_with_maldi(synthetic_project(), export, wrong))
    unusable = checks["every MALDI inspection row names an IO60 MALDI section, a known tissue, an "
                      "acquisition of the same donor and who inspected it"]
    assert not unusable["passed"]
    text = json.dumps(unusable["value"])
    assert "is donor 0112, not 6374" in text and "unknown tissue" in text
    assert "no IO60 MALDI tissue" in text and "nobody named" in text


# -------------------------------------------------------------------------------- the ingest ----
def test_ingest_keys_every_spot_and_keeps_zeros(synthetic_project, tmp_path, monkeypatch):
    monkeypatch.setattr(tasks.roots, "tmp_dir", lambda: None)
    export = write_export(tmp_path / "maldi", TISSUES)
    path = _with_maldi(synthetic_project(stages=("00_inventory", "04b_maldi_ingest"),
                                         automatic=("04b_maldi_ingest",)), export, FILLED)
    project, folder, _ = _inventory(path)
    gates.accept(folder, by="test")
    result = tasks.run_stage(stages.get("04b_maldi_ingest"), project)
    assert result.gate_passed and result.accepted_automatically
    spots = pq.read_table(result.path / "data/spots/tissue_id=6374_04/data.parquet").to_pandas()
    assert len(spots) == len(TISSUES["6374_04"])
    assert spots["spot_uid"].is_unique
    assert spots["spot_uid"].iloc[0].startswith("maldi%3A6374_04:spot50:")
    first = spots.iloc[0]
    assert first["on_tissue"] == False and first[GLYCANS[0]] == 0.0, "nothing detected: off tissue"  # noqa: E712
    assert spots[GLYCANS[0]].iloc[1] == 0.0 and spots["on_tissue"].iloc[1], "a 0 is kept"
    assert np.isnan(spots[GLYCOGEN[0]].iloc[2]), "a missing value stays missing"
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["summary"]["missing_values"] == 2
    assert receipt["settings"]["normalization"] == "none"


# -------------------------------------------------------------------------- the registration ----
def test_serial_registration_recovers_a_rotation_mirror_and_shift():
    from scipy import ndimage

    from multiplex5000.integration.rasterize import GridSpec
    from multiplex5000.integration.transform import from_rigid
    from multiplex5000.stages.maldi_serial_registration import register_section

    # an asymmetric section: an L of spots in the MALDI stage frame
    xy = np.array([(x, y) for x in np.arange(0, 3000, 50) for y in np.arange(0, 1200, 50)]
                  + [(x, y) for x in np.arange(0, 1000, 50) for y in np.arange(1200, 4000, 50)],
                  float) + 50_000
    truth = from_rigid(37.0, (-41_000.0, 7_500.0), mirror=True, center_um=(51_500.0, 52_000.0))
    fixed_xy = truth.apply(xy)
    grid = GridSpec(float(fixed_xy[:, 0].min() - 2000), float(fixed_xy[:, 1].min() - 2000), 25.0,
                    int(np.ptp(fixed_xy[:, 1]) / 25) + 160, int(np.ptp(fixed_xy[:, 0]) / 25) + 160)
    mask = np.zeros(grid.shape, bool)
    px = grid.to_pixels(fixed_xy)
    mask[np.rint(px[:, 1]).astype(int), np.rint(px[:, 0]).astype(int)] = True
    mask = ndimage.binary_closing(ndimage.binary_dilation(mask, iterations=2), iterations=2)
    rows, cols = np.nonzero(mask)
    centre = grid.to_microns(np.column_stack([cols.mean(), rows.mean()]))[0]
    fixed = dict(mask=mask, mask_grid=grid, origin=(grid.origin_x_um, grid.origin_y_um),
                 centre=centre, extent=float(np.hypot(*np.array(mask.shape[::-1]) * 25.0) / 2))
    moving = dict(xy=xy, centre=xy.mean(axis=0),
                  extent=float(np.max(np.hypot(*(xy - xy.mean(axis=0)).T))))
    transform, qc = register_section(moving, fixed)
    error = np.hypot(*(transform.apply(xy) - fixed_xy).T)
    assert qc["mirrored"] and transform.is_mirrored
    assert qc["tissue_dice"] > 0.85
    assert np.median(error) < 50.0, f"median spot error {np.median(error):.0f} um (one spot is 50)"


# --------------------------------------------------------------------------- the integration ----
def test_objects_are_evaluable_only_beyond_the_registration_error():
    import shapely

    from multiplex5000.stages.maldi_serial_integration import (_containing, _object_table,
                                                               _signed_distance)

    big = dict(object_id="a", object_name="Islet_1", object_class="Islet", size_class="Islet",
               geometry=shapely.box(0, 0, 400, 400))
    small = dict(object_id="b", object_name="Islet_2", object_class="Islet", size_class="Single",
                 geometry=shapely.box(1000, 1000, 1060, 1060))
    xy = np.array([(x, y) for x in range(25, 1200, 50) for y in range(25, 1200, 50)], float)
    points = shapely.points(xy)
    distance = _signed_distance(points, [big, small])
    inside = (xy[:, 0] < 400) & (xy[:, 1] < 400)
    in_small = (xy[:, 0] > 1000) & (xy[:, 0] < 1060) & (xy[:, 1] > 1000) & (xy[:, 1] < 1060)
    assert (distance[inside | in_small] < 0).all()
    assert (distance[~inside & ~in_small] > 0).all()
    assert (_containing(points, [big, small])[inside] == 0).all()
    values = np.ones((len(xy), 2))
    values[::2, 0] = 0.0
    rows = {r["object_id"]: r for r in _object_table([big, small], xy, values, ["f1", "f2"],
                                                     shrink_um=60 + 25, base={"tissue_id": "t"})}
    assert rows["a"]["evaluable"] and 0 < rows["a"]["detected_fraction__f1"] < 1
    assert rows["a"]["detected_fraction__f2"] == 1.0
    assert not rows["b"]["evaluable"] and rows["b"]["reason"]
    assert np.isnan(rows["b"]["detected_fraction__f1"]), "an unevaluable object carries no values"
