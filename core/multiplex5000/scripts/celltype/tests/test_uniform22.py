"""Per-donor inputs and isolated output for a changed UniFORM cohort."""
from __future__ import annotations

import json
import pathlib
import sys
import types
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import panel
import run_step2 as S2


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    flat, retired = tmp_path / "expression", tmp_path / "retired"
    monkeypatch.setattr(panel, "EXPRESSION_FLAT", flat)
    monkeypatch.setattr(panel, "EXPRESSION", retired)
    availability = {"112": ["CD38", "b_Catenin1"], "115": ["CD38"], "6591": []}
    monkeypatch.setattr(panel, "donor_markers", lambda d: availability[d])
    monkeypatch.setattr(panel, "cohort", lambda: pd.DataFrame({"donor_id": list(availability)}))

    def write(root, donor, markers, values=(2.0, 8.0, 16.0)):
        path = root / f"donor_id={donor}" / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame({
            "donor_id": [donor] * 3, "object_id": [10, 11, 12],
            "qc_analysis_eligible": [True, False, True],
            "X_centroid": [1.0, 2.0, 3.0], "Y_centroid": [4.0, 5.0, 6.0],
            "cell_region": ["tissue"] * 3, "islet_num": [0, 0, 1],
            **{m: values for m in markers},
        })
        df.to_parquet(path, index=False)
        return path

    write(flat, "0112", availability["112"])
    write(retired, "115", availability["115"])
    write(retired, "6591", availability["6591"])
    return flat, retired, write


@pytest.mark.parametrize("donor", ["112", "0112"])
def test_flat_partition_wins_even_when_retired_has_an_exact_alias(inputs, donor):
    flat, retired, write = inputs
    write(retired, "112", ["CD38"], (100.0, 200.0, 300.0))
    with warnings.catch_warnings(record=True) as caught:
        assert panel.expression_path(donor) == flat / "donor_id=0112" / "data.parquet"
    assert not caught


def test_retired_fallback_is_per_donor_and_warns(inputs):
    flat, retired, _ = inputs
    (flat / "donor_id=0115").mkdir()  # An empty export is not a usable partition.
    with pytest.warns(RuntimeWarning, match="Donor 115.*retired run tree"):
        assert panel.expression_path("0115") == retired / "donor_id=115" / "data.parquet"


def test_missing_expression_is_unavailable_not_an_empty_measurement(inputs):
    with pytest.raises(FileNotFoundError, match="unavailable for donor 999"):
        panel.expression_path("0999")


def test_duplicate_populated_aliases_are_ambiguous(inputs):
    flat, _, write = inputs
    write(flat, "112", ["CD38"])
    with pytest.raises(ValueError, match="Ambiguous expression partitions"):
        panel.expression_path("112")


def test_marker_values_use_both_roots_qc_and_actual_panel_membership(inputs):
    with pytest.warns(RuntimeWarning, match="Donor 115"):
        values = S2.marker_values("CD38", ["112", "115", "6591"])
    assert set(values) == {"112", "115"}
    for value in values.values():
        np.testing.assert_array_equal(value, [2.0, 16.0])
    assert set(S2.marker_values("b_Catenin1", ["112", "115", "6591"])) == {"112"}


@pytest.mark.parametrize("values", [(float("nan"),) * 3, (0.1, 0.2, 0.3)])
def test_unmeasurable_acquired_marker_has_no_synthetic_histogram(inputs, values):
    flat, _, write = inputs
    write(flat, "0112", ["CD38"], values)
    unavailable = {}
    assert S2.marker_values("CD38", ["112"], unavailable) == {}
    assert "No finite" in unavailable["112"]


def test_derive_excludes_missing_measurements_and_apply_preserves_them(inputs, tmp_path, monkeypatch):
    flat, _, write = inputs
    write(flat, "0112", ["CD38"], (float("nan"),) * 3)
    write(flat, "6533", ["CD38"])
    monkeypatch.setattr(panel, "donor_markers", lambda d: ["CD38"])

    def registration(hist, markers, reference_markers, donors, *args):
        assert donors == ["115", "6533"]
        assert all(sum(counts) == 2 for counts in hist["CD38"]["counts"])
        return {"CD38": [0, 0]}, ["115"], None

    reg = types.SimpleNamespace(automatic_registration=registration)
    norm = types.SimpleNamespace(calculate_shift_in_log_pixels=lambda *args: {"CD38": [0.0, 0.0]})
    monkeypatch.setattr(S2.UA, "_import_uniform", lambda: (reg, norm, None))
    monkeypatch.setattr(S2.UA, "gmm_cutoff", lambda *args: 2.0)
    with pytest.warns(RuntimeWarning, match="retired"):
        result = S2.derive(["CD38"], ["112", "115", "6533"],
                           pd.DataFrame({"marker": ["CD38"], "acq_ratio": [0.1]}),
                           save_figs=False, figs=tmp_path / "unused_figures")
    shift = result["CD38"]
    assert shift["n_donors"] == 2 and shift["n_acquired_donors"] == 3
    assert shift["scale"]["112"] == "unavailable"
    assert "112" not in shift["shift_bins"] and "112" not in shift["shift_log"]
    assert "No finite values" in shift["unavailable_donors"]["112"]
    assert not (tmp_path / "unused_figures").exists()
    S2.apply_shift(result, ["CD38"], ["112"], out=tmp_path / "uniform22")
    table = pq.ParquetFile(tmp_path / "uniform22/donor_id=112/data.parquet").read()
    assert table["CD38"].null_count == 3


def test_all_missing_measurements_are_recorded_as_unavailable(inputs, tmp_path, monkeypatch, capsys):
    flat, _, write = inputs
    write(flat, "0112", ["CD38"], (float("nan"),) * 3)
    monkeypatch.setattr(panel, "OUT", tmp_path)
    monkeypatch.setattr(S2.UA, "_import_uniform", lambda: (None, None, None))
    pd.DataFrame({"marker": ["CD38"], "acq_ratio": [0.1]}).to_parquet(
        tmp_path / "marker_donor_usability.parquet")
    out = tmp_path / "uniform22"
    assert S2.main(["--donors", "112", "--out", str(out), "--dry-run"]) == 0
    raw = (out / "shifts.json").read_text()
    shift = json.loads(raw)["shifts"]["CD38"]
    assert shift["provenance"] == "too_few_donors_to_align"
    assert shift["n_donors"] == 0 and shift["n_acquired_donors"] == 1
    assert shift["mode_sd_after"] == "unavailable"
    assert shift["reference_donor"] == "unavailable"
    assert shift["unavailable_reasons"]["reference_donor"]
    assert "NaN" not in raw
    assert "median unavailable" in capsys.readouterr().out


def test_apply_keeps_cells_keys_and_unacquired_markers_absent(inputs, tmp_path):
    out = tmp_path / "uniform22"
    shifts = {m: {"scale": {"112": 2.5, "115": 0.5, "6591": 1.0}}
              for m in ["CD38", "b_Catenin1"]}
    with pytest.warns(RuntimeWarning):
        S2.apply_shift(shifts, list(shifts), ["112", "115", "6591"], out=out)
    for donor in ["112", "115", "6591"]:
        table = pq.ParquetFile(out / f"donor_id={donor}" / "data.parquet").read().to_pandas()
        assert table.donor_id.tolist() == [donor] * 3
        assert table.object_id.tolist() == [10, 11, 12]
        assert table.qc_analysis_eligible.tolist() == [True, False, True]
        assert set(table.columns) & set(shifts) == set(panel.donor_markers(donor))
        for marker in panel.donor_markers(donor):
            np.testing.assert_array_equal(table[marker], np.array([2, 8, 16]) * shifts[marker]["scale"][donor])


@pytest.mark.parametrize("custom,explicit_shift,dry_run", [
    (False, False, False), (True, False, False), (True, True, False), (True, True, True),
])
def test_cli_routes_all_outputs_and_deduplicates_donor_aliases(
        inputs, tmp_path, monkeypatch, custom, explicit_shift, dry_run):
    historical = tmp_path / "celltype"
    historical.mkdir()
    monkeypatch.setattr(panel, "OUT", historical)
    monkeypatch.setattr(S2, "OUT", historical / "uniform")
    monkeypatch.setattr(S2, "FIGS", historical / "figures_step2")
    pd.DataFrame({"marker": ["CD38"], "acq_ratio": [0.1]}).to_parquet(
        historical / "marker_donor_usability.parquet")
    sentinels = [historical / "shifts.json", S2.OUT / "sentinel", S2.FIGS / "CD38.png"]
    for p in sentinels:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("existing baseline")

    def derive(markers, donors, usab, *, figs):
        assert donors == ["112", "115"]
        assert markers == ["CD38"]
        figs.mkdir(parents=True, exist_ok=True)
        (figs / "CD38.png").write_text("new figure")
        return {"CD38": {"scale": {d: 2.0 for d in donors}, "mode_sd_after": 0.01,
                         "acq_ratio": float("nan"),
                         "provenance": "automatic_registration", "over_bar": False}}

    monkeypatch.setattr(S2, "derive", derive)
    out = historical / "uniform22" if custom else S2.OUT
    args = ["--donors", "0112", "112", "115"]
    if custom:
        args += ["--out", str(out)]
    shifts_path = (out / "shifts.json") if custom else historical / "shifts.json"
    if explicit_shift:
        shifts_path = historical / "shifts_22donor.json"
        args += ["--shifts-out", str(shifts_path)]
    if dry_run:
        args += ["--dry-run"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert S2.main(args) == 0
    meta = json.loads(shifts_path.read_text())
    assert meta["donors"] == ["112", "115"]
    assert meta["expression_inputs"]["112"].endswith("donor_id=0112/data.parquet")
    assert meta["shifts"]["CD38"]["acq_ratio"] == "unavailable"
    assert "Step-1 usability input" in meta["shifts"]["CD38"]["unavailable_reasons"]["acq_ratio"]
    assert "NaN" not in shifts_path.read_text()
    assert (out / "donor_id=112" / "data.parquet").exists() is not dry_run
    assert (out / "donor_id=115" / "data.parquet").exists() is not dry_run
    assert not (out / "donor_id=0112").exists()
    if custom:
        assert all(p.read_text() == "existing baseline" for p in sentinels)
        assert (out / "figures_step2" / "CD38.png").read_text() == "new figure"


@pytest.mark.parametrize("donor", ["6457", "6579", "unknown"])
def test_cli_rejects_ineligible_or_unknown_donors_before_io(inputs, donor):
    with pytest.raises(SystemExit) as exc:
        S2.main(["--donors", donor])
    assert exc.value.code == 2


def test_cli_defaults_to_the_entire_eligible_cohort(inputs, tmp_path, monkeypatch):
    monkeypatch.setattr(panel, "OUT", tmp_path)
    pd.DataFrame({"marker": ["CD38"], "acq_ratio": [0.1]}).to_parquet(
        tmp_path / "marker_donor_usability.parquet")

    def derive(markers, donors, usab, *, figs):
        assert donors == ["112", "115", "6591"]
        return {"CD38": {"scale": {d: 1.0 for d in donors}, "mode_sd_after": 0.01,
                         "provenance": "automatic_registration", "over_bar": False}}

    monkeypatch.setattr(S2, "derive", derive)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert S2.main(["--out", str(tmp_path / "uniform22"), "--dry-run"]) == 0
    meta = json.loads((tmp_path / "uniform22" / "shifts.json").read_text())
    assert meta["donors"] == ["112", "115", "6591"]
