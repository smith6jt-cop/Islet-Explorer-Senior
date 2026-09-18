"""Availability and section contracts, including both CLIs against the real trees.

Run with PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider
core/phenocycler/tests/current/test_availability_layer.py. Artifacts go to tmp_path.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts.senior import xenium_release

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior import build_availability as ledger  # noqa: E402
from scripts.senior import build_sections_layer as sections  # noqa: E402

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

AVAILABILITY_COLUMNS = ["donor_id", "capability", "state", "reason", "n", "source"]
SECTION_COLUMNS = ["donor_id", "roi", "modality", "image_file", "xenium_serial",
                   "pixel_size_um", "frame",
                   # Derived from the per-image annotation census; see
                   # test_section_relationship.py for what each value has to be evidenced by.
                   "section_relationship", "section_relationship_basis"]



def _release_called_cells(donor: str) -> int | None:
    """Cells the current pancreas Xenium release carries a LINEAGE for, from its own summary.

    Read from `lineage_counts.csv` rather than from `calls/`, so that this expectation and the
    builder do not both come from the same file. Returns None when no release is resolvable, so
    the caller skips rather than asserting against an absence.
    """
    release = xenium_release.resolve("pancreas")
    if not getattr(release, "is_measured", lambda: bool(getattr(release, "root", None)))():
        return None
    counts = Path(release.root) / "lineage_counts.csv"
    if not counts.exists():
        return None
    with counts.open() as handle:
        for row in csv.DictReader(handle):
            if str(row.get("donor_id")) != str(donor):
                continue
            return sum(int(float(v)) for k, v in row.items()
                       if k != "donor_id" and v not in (None, "")
                       and not k.startswith(("Uncalled_", "Unavailable_")))
    return None

def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_parquet(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.table(values), path)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def source_tree(tmp_path):
    repo = tmp_path / "repo"
    inventory = repo / "core/phenocycler/data/integration"
    excluded = next(iter(DONOR_EXCLUSIONS))
    paths = []
    for donor, serial in (("0112", "0041396"), ("115", "0041117"),
                          ("6436", "0059907"), ("6442", "0059947"),
                          ("6591", "0041116"), (excluded, "0099999")):
        paths.append(dict(donor_id=donor, roi="panc", batch="fixture",
                          source_path=f"batch/output-XETG00298__{serial}__Panc__20240919__203848"))
    paths += [
        dict(donor_id="112", roi="pln_1", batch="fixture",
             source_path="batch/output-XETG00298__0041396__pLN__20240919__203848"),
        dict(donor_id="6591", roi="pln_1", batch="fixture", source_path="NA"),
        dict(donor_id="112", roi="panc", batch="second_section",
             source_path="batch/output-XETG00298__0080870__panc__20260910__211917"),
    ]
    write_csv(inventory / "xenium_paths.csv", paths)
    write_csv(inventory / "restain_sections.csv", [
        dict(image_file=f"{r['donor_id']}.qptiff", donor_id=r["donor_id"], roi="panc",
             xenium_serial=sections.parse_bundle_path(r["source_path"])[0], pixel_size_um="0.2712")
        for r in paths if r["donor_id"] in {"115", "6436", "6442", excluded}
    ])
    write_csv(inventory / "pln_node_annotations.csv", [
        dict(donor_id="0112", roi="pln_1", image_file="a.qptiff", n_nodes="2"),
        dict(donor_id=excluded, roi="pln_1", image_file="excluded.qptiff", n_nodes="999"),
    ])
    typing = repo / "data/development/typing"
    write_parquet(typing / "cells/donor_id=0112/data.parquet", {
        "object_id": ["a", "b", "c"], "qc_analysis_eligible": [True, False, True],
    })
    write_parquet(typing / f"cells/donor_id={excluded}/data.parquet", {
        "object_id": ["excluded"], "qc_analysis_eligible": [True],
    })
    # Physical donor columns win over a deliberately misleading directory. The
    # large_string type also reproduces the hive donor_id collision in real data.
    write_parquet(inventory / "cells_xen/donor_id=9999/roi=panc/data.parquet", {
        "donor_id": pa.array(["0112"] * 4 + ["115", "6436", "6442", excluded],
                             type=pa.large_string()),
        "roi": ["panc"] * 8,
        "lineage_native": ["", "Unavailable", "Immune_T", None, "", "Unavailable",
                           "Endocrine", "Immune_T"],
        "lineage_common": [""] * 8,
        "feat__if_INS": [None] * 5 + [1.0, 2.0, None],
        "feat__if_GCG": [None] * 5 + [3.0, 4.0, None],
    })
    write_parquet(typing / "geometry/outlines.parquet", {
        "donor_id": ["0112", "112", "115", excluded],
        "islet_key": ["Islet_1", "Islet_2", "Islet_1", "Islet_1"],
        "wkb": [b"exported_wkb", None, b"", b"excluded"],
    })
    write_parquet(typing / "trajectory/singlecell_islet_summary.parquet", {
        "case_id": [112] * 4 + [115],
        "islet_key": ["Islet_1", "Islet_2", "Islet_3", "Islet_4", "Islet_1"],
        "sc_pseudotime_core": [1.0, None, float("nan"), 3.0, None],
        "sc_pseudotime_peri": [None, 2.0, float("inf"), None, None],
    })
    restain = repo / "data/development/restain"
    write_json(restain / "hormone_panc_refined.json", [
        dict(donor_id=donor, roi="panc", hormone={
            "preregistered_auroc_bar": 0.75,
            "INS": dict(auroc=value, clears_preregistered_bar=passed),
            "GCG": dict(auroc=value, clears_preregistered_bar=passed),
        }) for donor, value, passed in (("6436", 0.52, False), ("6442", 0.85, True))
    ])
    write_json(restain / "valis/donor_id=115/alignment_refined.json", {
        "donor_id": "115", "roi": "panc", "alignment": {
            "verdict": "FAIL", "median_offset_um": 7.198, "verdict_why": "recorded failure",
        },
    })
    return repo


def assert_ledger_contract(table):
    assert table.column_names == AVAILABILITY_COLUMNS
    assert table.schema == ledger.AVAILABILITY_SCHEMA
    donors = set(table["donor_id"].to_pylist())
    assert not donors.intersection(DONOR_EXCLUSIONS)
    assert all(d == sections.norm_donor(d) for d in donors)
    rows = table.to_pylist()
    assert len(rows) == len(donors) * len(ledger.CAPABILITIES)
    assert len({(r["donor_id"], r["capability"]) for r in rows}) == len(rows)
    assert {r["capability"] for r in rows} == set(ledger.CAPABILITIES)
    for row in rows:
        assert row["state"] in ledger.STATES
        assert row["source"]
        if row["state"] != "measured":
            assert row["n"] is None, row  # Null, never a missing measurement disguised as zero.
            assert row["reason"].strip(), row
        else:
            assert isinstance(row["n"], int) and row["n"] >= 0, row


def test_five_states_and_null_counts_survive_parquet(source_tree, tmp_path):
    table = ledger.build_availability(source_tree)
    path = sections.write_layer(table, tmp_path / "layers", "availability")
    table = pq.ParquetFile(path).read()
    assert_ledger_contract(table)
    assert set(table["state"].to_pylist()) == set(ledger.STATES)
    rows = {(r["donor_id"], r["capability"]): r for r in table.to_pylist()}
    assert rows["112", "phenocycler_panc"]["n"] == 2
    assert rows["112", "xenium_panc"]["n"] == 1
    assert "blank=2; qc_unavailable=1" in rows["112", "xenium_panc"]["reason"]
    assert rows["115", "xenium_panc"]["state"] == "blank"
    assert rows["6436", "xenium_panc"]["state"] == "unavailable"
    assert rows["6436", "restain_if"]["state"] == "withheld"
    assert "INS_auroc=0.52" in rows["6436", "restain_if"]["reason"]
    assert rows["6442", "restain_if"]["n"] == 1
    assert rows["115", "restain_if"]["state"] == "failed"
    assert "7.198" in rows["115", "restain_if"]["reason"]
    assert rows["112", "islet_outline"]["n"] == 1
    assert rows["115", "islet_outline"]["state"] == "unavailable"
    assert rows["112", "trajectory"]["n"] == 3
    assert rows["115", "trajectory"]["state"] == "unavailable"
    assert rows["112", "pln_nodes"]["n"] == 2
    assert rows["112", "phenocycler_pln"]["n"] == 1
    assert rows["6591", "xenium_pln"]["reason"] == "xenium_source_path_recorded_unavailable"


@pytest.mark.parametrize("state", ["blank", "unavailable", "withheld", "failed"])
def test_nonmeasured_cannot_be_written_as_zero_or_without_reason(state):
    with pytest.raises(ValueError, match="n=None"):
        ledger.availability("112", "trajectory", state, Path("source"), reason="absent", n=0)
    with pytest.raises(ValueError, match="requires a reason"):
        ledger.availability("112", "trajectory", state, Path("source"))
    assert ledger.availability("112", "trajectory", state, Path("source"), reason="absent")["n"] is None


def test_a_measured_zero_is_distinct_from_absence():
    row = ledger.availability("0112", "pln_nodes", "measured", Path("source"), n=0)
    assert row["n"] == 0 and row["donor_id"] == "112"


def test_sections_preserve_na_scales_serials_and_tissue_frames(source_tree):
    table = sections.build_sections(source_tree / "core/phenocycler/data/integration")
    assert table.column_names == SECTION_COLUMNS
    assert table.schema == sections.SECTIONS_SCHEMA
    rows = table.to_pylist()
    assert not set(table["donor_id"].to_pylist()).intersection(DONOR_EXCLUSIONS)
    if_rows = [r for r in rows if r["modality"] == "restain_if"]
    assert {r["pixel_size_um"] for r in if_rows} == {0.2712}
    xen = [r for r in rows if r["modality"] == "xenium"]
    assert all(r["pixel_size_um"] is None for r in xen)
    missing = next(r for r in xen if r["donor_id"] == "6591" and r["roi"] == "pln_1")
    assert missing["image_file"] is None and missing["xenium_serial"] is None
    assert missing["frame"] == "xenium_um:6591:pln_1"
    same_serial = [r for r in xen if r["xenium_serial"] == "0041396"]
    assert {r["roi"] for r in same_serial} == {"panc", "pln_1"}
    assert len({r["frame"] for r in same_serial}) == 2
    assert len([r for r in xen if r["donor_id"] == "112" and r["roi"] == "panc"]) == 2


def test_blank_required_inventory_fields_are_refused(source_tree):
    inventory = source_tree / "core/phenocycler/data/integration"
    path = inventory / "xenium_paths.csv"
    path.write_text(path.read_text().replace(",NA", ","))
    with pytest.raises(ValueError, match="required field source_path"):
        sections.build_sections(inventory)


def test_an_inconclusive_gate_does_not_release_paired_restain(source_tree):
    gate = source_tree / "data/development/restain/valis/donor_id=6442/alignment_refined.json"
    write_json(gate, {"donor_id": "6442", "roi": "panc", "alignment": {
        "verdict": "INCONCLUSIVE", "verdict_why": "insufficient measurable windows",
    }})
    rows = ledger.build_availability(source_tree).to_pylist()
    row = next(r for r in rows if r["donor_id"] == "6442" and r["capability"] == "restain_if")
    assert row["state"] == "unavailable" and row["n"] is None
    assert "inconclusive" in row["reason"]


@pytest.mark.data
def test_both_builders_run_against_real_trees(tmp_path):
    core_data = ledger.default_phenocycler_data(REPO)
    if not (core_data / "integration/cells_xen").is_dir():
        pytest.skip("Real PhenoCycler data tree is not mounted")
    if not (REPO / "data/development/typing/cells").is_dir():
        pytest.skip("Real typing data tree is not mounted")
    layers = tmp_path / "islet_layers"
    environment = dict(os.environ, ISLET_LAYERS_DIR=str(layers), PYTHONDONTWRITEBYTECODE="1")
    for builder in ("build_sections_layer.py", "build_availability.py"):
        subprocess.run([sys.executable, str(REPO / "scripts/senior" / builder)],
                       cwd=tmp_path, env=environment, check=True, capture_output=True, text=True)
    table = pq.ParquetFile(layers / "availability/data.parquet").read()
    assert_ledger_contract(table)
    # SUBSET, not equality. Every state the builder emits must be in the vocabulary, but the
    # vocabulary does not have to be exercised on the real trees. `blank` was only ever produced
    # by reading the UNTYPED cells_xen ingest, whose lineage column is empty on every row. Now
    # that the Xenium arms resolve through xenium_release, a donor the release carries a call for
    # is `measured`, and requiring `blank` to appear here would assert that the defect is still
    # present. `blank` stays in ledger.STATES because a never-typed section is a real state and
    # the app renders it -- test_five_states_and_null_counts_survive_parquet still pins all five
    # against a fixture that constructs them.
    assert set(table["state"].to_pylist()) <= set(ledger.STATES)
    assert table.num_rows > 0, "the builder emitted no rows at all"
    row_map = {(r["donor_id"], r["capability"]): r for r in table.to_pylist()}
    for donor in ("6436", "6442"):
        assert row_map[donor, "restain_if"]["state"] == "withheld"
        # The xenium_panc count is pinned to the RELEASE, not to cells_xen. The untyped ingest
        # carries no call on any row, so counting its non-empty lineage_native counted cells that
        # were never typed -- which is what the builder used to report as `blank`. The expectation
        # here is read from the release's own published per-donor summary, a different file from
        # the calls/ tree the builder reads, so the check is not circular. A call is a lineage that
        # is neither Uncalled_* (no naming evidence, no enrichment, below cluster QC) nor
        # Unavailable_no_nucleus, which is the nucleus rule this repo applies everywhere.
        expected = _release_called_cells(donor)
        if expected is None:
            pytest.skip(f"no current pancreas Xenium release to pin {donor} against")
        assert row_map[donor, "xenium_panc"]["n"] == expected
    assert row_map["115", "restain_if"]["state"] == "failed"
    assert row_map["6380", "pln_nodes"]["n"] == 1  # Its Tissue annotation is the one node.
    output = pq.ParquetFile(layers / "sections/data.parquet").read()
    assert output.column_names == SECTION_COLUMNS
    assert output.schema == sections.SECTIONS_SCHEMA
    inventory = REPO / "core/phenocycler/data/integration"
    restain, xenium = sections.section_inventories(inventory)
    assert output.num_rows == len(restain) + len(xenium)
    rows = output.to_pylist()
    expected_scales = {(r["donor_id"], r["roi"]): float(r["pixel_size_um"]) for r in restain}
    for row in rows:
        assert row["donor_id"] not in DONOR_EXCLUSIONS
        if row["modality"] == "restain_if":
            assert row["pixel_size_um"] == expected_scales[row["donor_id"], row["roi"]]
        else:
            assert row["pixel_size_um"] is None
    missing = next(r for r in rows if r["donor_id"] == "6591" and r["roi"] == "pln_1")
    assert missing["image_file"] is None and missing["xenium_serial"] is None
    assert row_map["6591", "xenium_pln"]["state"] == "unavailable"
