"""Exercise the maintained command through all four stages on native fixtures."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from scipy import sparse

from multiplex5000.gammadelta.cli import sha
import test_gammadelta_selection as fixtures

CRITERIA = fixtures.CRITERIA


def complete_fixture(case, parquet):
    args = case.fixture(parquet=parquet)
    inventory = json.loads(args.inventory.read_text())
    row = inventory["runs"][0]
    bundle = Path(row["bundle"])
    with h5py.File(bundle / "cell_feature_matrix.h5") as handle:
        matrix = handle["matrix"]
        dense = sparse.csc_matrix((matrix["data"][:], matrix["indices"][:],
                                   matrix["indptr"][:]), shape=matrix["shape"][:]).toarray()
        features = [v.decode() for v in matrix["features/name"][:]]
        cells = [v.decode() for v in matrix["barcodes"][:]]
    meta_file = bundle / ("cells.parquet" if parquet else "cells.csv.gz")
    meta = pd.read_parquet(meta_file) if parquet else pd.read_csv(meta_file, float_precision="round_trip")
    rows = []
    for c, cell in enumerate(cells):
        for g, gene in enumerate(features):
            for count in range(int(dense[g, c])):
                rows.append((cell, gene, float(meta.x_centroid[c]), float(meta.y_centroid[c]), 30.))
    # Preserve a low-QV assigned transcript and unassigned native background.
    rows.extend([(cells[0], features[0], float(meta.x_centroid[0]) + 1000,
                  float(meta.y_centroid[0]) + 1000, 2.),
                 ("UNASSIGNED", "background", float(meta.x_centroid[0]),
                  float(meta.y_centroid[0]), 5.)])
    frame = pd.DataFrame(rows, columns=["cell_id", "feature_name", "x_location", "y_location", "qv"])
    frame["transcript_id"] = np.arange(len(frame), dtype=np.uint64)
    frame["z_location"] = 1.
    frame["overlaps_nucleus"] = False
    frame["native_extra"] = [None] + ["preserve"] * (len(frame) - 1)
    table = pa.Table.from_pandas(frame, preserve_index=False).replace_schema_metadata({b"fixture": b"native"})
    transcript = bundle / "transcripts.parquet"
    pq.write_table(table, transcript, row_group_size=7)
    row["sources"][transcript.name] = {
        "path": str(transcript), "sha256": sha(transcript), "bytes": transcript.stat().st_size,
        "mtime_ns": transcript.stat().st_mtime_ns, "source_rows": len(table),
    }
    inputs = case.root / "inputs"
    inputs.mkdir()
    (inputs / "native_sources.json").write_text(json.dumps(inventory))
    config = json.loads(args.config.read_text())
    config.update(source_inventory_sha256=sha(inputs / "native_sources.json"),
                  roi_half_width_um=30., qv_reconciliation_minimum=20., preserve_native_transcript_columns=True)
    (inputs / "approved_runtime_config.json").write_text(json.dumps(config))
    (inputs / "historical_and").mkdir()
    shutil.copyfile(args.historical_and_cells, inputs / "historical_and" / (case.region + ".csv"))
    shutil.copyfile(args.historical_or2_cells, inputs / "historical_or2_cells.csv")
    shutil.copyfile(CRITERIA, inputs / "reviewed_criteria.json")
    return inputs


@pytest.mark.parametrize("parquet", [False, True])
def test_normal_cli_runs_all_stages_preserves_flags_and_separate_eligibility(parquet):
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet)
        output = case.root / "result"
        # Use the normal public CLI from a neutral working directory. The same
        # module layout is exercised in the separately recorded wheel check.
        command = [sys.executable, "-m", "multiplex5000.cli", "gammadelta",
                   "--inputs-root", str(inputs), "--region", case.region,
                   "--output-dir", str(output), "--chunk-cells", "3", "--batch-size", "5"]
        env = dict(os.environ, PYTHONPATH=str(Path(CRITERIA).parents[2]), PYTHONDONTWRITEBYTECODE="1")
        result = subprocess.run(command, cwd=case.root, env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        receipt = json.loads((output / "pipeline.json").read_text())
        assert receipt["status"] == "COMPLETE_TECHNICAL_VALIDATION"
        assert len(receipt["steps"]) == 4
        assert len(receipt["validation_receipts"]) == 2
        eligible = pd.read_csv(output / "selector/eligible_raw_cells.csv")
        assert len(eligible) == 20
        assert eligible.foreign_flag_union.sum() == 16
        assert {0, 2} <= set(eligible.nucleus_count)
        assert not receipt["nucleus_filter_applied"] and not receipt["tissue_filter_applied"]
        assert not receipt["biological_identity_accepted"] and not receipt["release_activated"]
        before = sha(output / "pipeline.json")
        repeat = subprocess.run(command, cwd=case.root, env=env, capture_output=True, text=True)
        assert repeat.returncode != 0
        assert "Refusing existing" in repeat.stderr
        assert sha(output / "pipeline.json") == before
    finally:
        case.tearDown()


def test_main_help_lists_gammadelta(capsys):
    from multiplex5000 import cli
    with pytest.raises(SystemExit) as result:
        cli.main(["--help"])
    assert result.value.code == 0
    assert "gammadelta" in capsys.readouterr().out


def _gamma_args(inputs, case, output, *extra):
    return ["--inputs-root", str(inputs), "--region", case.region,
            "--output-dir", str(output), *extra]


def test_dry_run_reports_all_steps_without_creating_output(capsys):
    from multiplex5000.gammadelta.cli import main
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet=True)
        output = case.root / "planned all-stage result"
        before = {str(p): sha(p) for p in inputs.rglob("*") if p.is_file()}
        assert main(_gamma_args(inputs, case, output, "--dry-run")) == 0
        report = json.loads(capsys.readouterr().out)
        assert report["status"] == "READY_FOR_EXECUTION"
        assert report["side_effects"] is False
        assert report["native_content_verified"] is False
        assert report["scientific_acceptance"] is False
        assert report["inventory_row_key"] == case.region
        assert report["output_dir"] == str(output)
        assert [c[2].rsplit(".", 1)[-1] for c in report["commands"]] == [
            "select_corrected", "validate_selection", "extract_corrected", "validate_coordinates"]
        assert set(report["required_sources"]) == {
            "cell_feature_matrix.h5", "cells.parquet", "transcripts.parquet"}
        assert all(check["matches_pinned_stat"] for check in report["required_sources"].values())
        assert not output.exists()
        assert {str(p): sha(p) for p in inputs.rglob("*") if p.is_file()} == before
    finally:
        case.tearDown()


def test_missing_transcripts_block_coordinate_work_before_output_but_not_selection(capsys):
    from multiplex5000.gammadelta.cli import main
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet=False)
        row = json.loads((inputs / "native_sources.json").read_text())["runs"][0]
        Path(row["sources"]["transcripts.parquet"]["path"]).unlink()  # synthetic fixture only
        output = case.root / "must not be created"
        argv = _gamma_args(inputs, case, output)
        assert main([*argv, "--dry-run"]) == 2
        report = json.loads(capsys.readouterr().out)
        assert report["status"] == "BLOCKED"
        assert report["required_sources"]["transcripts.parquet"]["status"] == "not_found_on_this_host"
        assert report["blockers"]
        assert not output.exists()
        assert main(argv) == 1
        assert "not_found_on_this_host" in capsys.readouterr().err
        assert not output.exists()
        assert main([*argv, "--stage", "selection", "--dry-run"]) == 0
        selection_report = json.loads(capsys.readouterr().out)
        assert "transcripts.parquet" not in selection_report["required_sources"]
        assert len(selection_report["commands"]) == 2
        assert not output.exists()
    finally:
        case.tearDown()


def test_dry_run_preserves_independently_pinned_transcript_location(capsys):
    from multiplex5000.gammadelta.cli import main
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet=True)
        inventory_path = inputs / "native_sources.json"
        inventory = json.loads(inventory_path.read_text())
        pin = inventory["runs"][0]["sources"]["transcripts.parquet"]
        replica = case.root / "separate transcript replica" / "transcripts.parquet"
        replica.parent.mkdir()
        shutil.move(pin["path"], replica)  # preserve bytes/mtime within a synthetic fixture
        pin["path"] = str(replica)
        inventory_path.write_text(json.dumps(inventory))
        config_path = inputs / "approved_runtime_config.json"
        config = json.loads(config_path.read_text())
        config["source_inventory_sha256"] = sha(inventory_path)
        config_path.write_text(json.dumps(config))
        output = case.root / "mixed-location result"
        assert main(_gamma_args(inputs, case, output, "--dry-run")) == 0
        report = json.loads(capsys.readouterr().out)
        assert report["required_sources"]["transcripts.parquet"]["path"] == str(replica)
        assert report["required_sources"]["transcripts.parquet"]["matches_pinned_stat"]
        assert report["native_bundle"] == inventory["runs"][0]["bundle"]
        assert not output.exists()
    finally:
        case.tearDown()


@pytest.mark.parametrize("change, message", [
    ("threshold", "threshold differs"),
    ("criteria", "Reviewed criteria hash mismatch"),
    ("historical", "Historical membership hash mismatch"),
])
def test_dry_run_rejects_altered_approved_contract_before_output(capsys, change, message):
    from multiplex5000.gammadelta.cli import main
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet=True)
        if change == "threshold":
            path = inputs / "approved_runtime_config.json"
            config = json.loads(path.read_text())
            config["per_gene_minimum"] = 3
            path.write_text(json.dumps(config))
        elif change == "criteria":
            with (inputs / "reviewed_criteria.json").open("a") as handle:
                handle.write("\n")
        else:
            with (inputs / "historical_or2_cells.csv").open("a") as handle:
                handle.write("\n")
        output = case.root / "invalid-contract result"
        assert main(_gamma_args(inputs, case, output, "--dry-run")) == 1
        assert message in capsys.readouterr().err
        assert not output.exists()
    finally:
        case.tearDown()


def test_coordinate_preflight_requires_unchanged_independently_validated_selection(capsys):
    from multiplex5000.gammadelta.cli import main
    case = fixtures.SelectionTests()
    case.setUp()
    try:
        inputs = complete_fixture(case, parquet=True)
        selected = case.root / "existing selection"
        assert main(_gamma_args(inputs, case, selected, "--stage", "selection")) == 0
        capsys.readouterr()
        output = case.root / "new coordinates"
        selection_dir = selected / "selector"
        argv = _gamma_args(inputs, case, output, "--stage", "coordinates",
                           "--selection-dir", str(selection_dir), "--dry-run")
        assert main(argv) == 0
        report = json.loads(capsys.readouterr().out)
        assert report["selection_dir"] == str(selection_dir)
        assert len(report["commands"]) == 2
        assert not output.exists()
        with (selection_dir / "eligible_raw_cells.csv").open("a") as handle:
            handle.write("\n")
        assert main(argv) == 1
        assert "Eligible table differs" in capsys.readouterr().err
        assert not output.exists()
    finally:
        case.tearDown()
