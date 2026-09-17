"""Restain acceptance: measured direction, CSV calibration, and explicit absence.

Run with PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider
core/phenocycler/tests/current/test_restain_layer.py
"""

from __future__ import annotations

import builtins
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "scripts/senior/build_restain_layer.py"
SPEC = importlib.util.spec_from_file_location("restain_layer", SCRIPT)
layer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(layer)


@pytest.fixture
def declared_sections(tmp_path):
    """A real TIFF disagrees with its CSV, as the two 10x images do in production."""
    import tifffile

    image = tmp_path / "source.qptiff"
    tifffile.imwrite(image, np.zeros((4, 4), dtype=np.uint8),
                     description="<PixelSizeMicrons>0.50139</PixelSizeMicrons>")
    path = tmp_path / "restain_sections.csv"
    path.write_text(
        "# Maintainer declares 0.2505; the image records 0.50139.\n"
        "donor_id,roi,image_file,xenium_serial,pixel_size_um\n"
        f"06436,panc,{image},0059907,0.2505\n"
        f"6442,panc,{image},0059947,0.2505\n"
        f"0112,panc,{image},0041396,0.2505\n"
    )
    return path


@pytest.fixture
def forbid_qptiff(monkeypatch):
    """Trip on file access as well as on either TIFF-reader entry point."""
    import tifffile

    attempts = []
    for owner in (builtins, io):
        original = owner.open

        def checked(file, *args, _open=original, **kwargs):
            if isinstance(file, (str, bytes, os.PathLike)) and str(file).lower().endswith(".qptiff"):
                attempts.append(str(file))
                raise AssertionError("builder attempted to read a qptiff")
            return _open(file, *args, **kwargs)

        monkeypatch.setattr(owner, "open", checked)

    def forbidden(*args, **kwargs):
        attempts.append("TIFF reader")
        raise AssertionError("builder attempted to read a qptiff")

    monkeypatch.setattr(tifffile, "TiffFile", forbidden)
    monkeypatch.setattr(tifffile, "imread", forbidden)
    return attempts


def test_csv_pixel_size_overrules_the_qptiff(declared_sections, forbid_qptiff):
    sections = layer.read_sections(declared_sections)
    assert sections.pixel_size_um.tolist() == [0.2505] * 3
    assert sections.donor_id.tolist() == ["6436", "6442", "112"]
    assert not forbid_qptiff


@pytest.fixture
def builder_inputs(tmp_path, declared_sections):
    """Geometry order differs from the bundle; one cell has no recorded nucleus."""
    source = tmp_path / "cells_xen"
    source.mkdir()
    pq.write_table(pa.table({"donor_id": ["06436"] * 4, "roi": ["panc"] * 4,
                            "cell_id": ["cell-b", "cell-d", "cell-a", "cell-c"],
                            "feat__if_INS": [20.0, 999.0, 0.0, np.nan],
                            "feat__if_GCG": [40.0, 999.0, 30.0, np.nan]}), source / "cells.parquet")
    pq.write_table(pa.table({"donor_id": ["0112"], "roi": ["panc"],
                            "cell_id": ["unpaired"]}), source / "unpaired.parquet")
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    order = np.array(["cell-a", "cell-b", "cell-c", "cell-d"])
    pq.write_table(pa.table({"cell_id": order}), bundle / "cells.parquet")
    restain = tmp_path / "inputs"
    warp = restain / "valis" / "export-one"
    warp.mkdir(parents=True)
    (warp / "export.json").write_text(json.dumps({
        "donor_id": "6436", "roi": "panc", "serial": "0059907",
        "xenium_bundle": str(bundle), "qptiff": str(tmp_path / "source.qptiff"),
        "if_um_per_px": 1.002, "if_pixel_um": 0.2505,
        "pixel_size_source": "fixture producer calibration", "if": {"factor": 4}
    }))
    np.save(warp / "xenium_centroids_um.npy", np.array([[1., 1.], [2., 2.], [3., 3.], [4., 4.]]))
    from phenocycler.integration.rasterize import GridSpec

    layer.Transform(affine=[[2, 0, 100], [0, 3, 300]],
                    displacement=np.array([np.full((2, 2), 5), np.full((2, 2), -7)]),
                    grid=GridSpec(origin_x_um=0, origin_y_um=0, pixel_um=1000,
                                  height=2, width=2)).save(warp, stem="refined")
    np.save(warp / "warped_if_um.npy", np.array([[101.0, 301.0], [102.0, 302.0],
                                               [103.0, 303.0], [104.0, 304.0]]))
    np.save(warp / "xenium_cell_ids.npy", order.astype(object))
    np.save(warp / "xenium_nuclei_cell_ids.npy", np.array(["cell-b", "cell-a", "cell-c", "cell-b"]))
    (restain / "protein_rna_panc.json").write_text(json.dumps([{
        "donor_id": "6436", "roi": "panc", "xenium_serial": "0059907",
        "xenium_bundle": str(bundle), "qptiff": str(tmp_path / "source.qptiff"),
        "n_cells": 4, "transform_source": "refined",
        "protein_um_per_px": 0.501, "pyramid_factor": 2
    }]))
    (restain / "hormone_panc_refined.json").write_text(json.dumps([{
        "donor_id": "6436", "roi": "panc", "hormone": {
            "preregistered_auroc_bar": 0.75,
            "INS": {"auroc": 0.5143, "clears_preregistered_bar": False},
            "GCG": {"auroc": 0.5220, "clears_preregistered_bar": False}
        }
    }]))
    paths = tmp_path / "xenium_paths.csv"
    paths.write_text("donor_id,roi\n6436,panc\n0112,panc\n6347,panc\n6457,panc\n6579,panc\n")
    return dict(cells_xen_root=source, restain_root=restain,
                sections_csv=declared_sections, xenium_paths_csv=paths,
                layers_dir=tmp_path / "layers")


def test_builder_csv_calibration_and_if_overlay_acceptance(builder_inputs, forbid_qptiff):
    """The full builder uses CSV calibration and joins bundle order by cell_id."""
    result = layer.build_restain_layer(**builder_inputs).set_index("donor_id")
    assert result.loc["6436", "pixel_size_um"] == 0.2505
    assert not forbid_qptiff
    assert set(result.index) == {"112", "6347", "6436", "6442"}
    assert result.loc["6436", "state"] == "measured"
    assert result.loc["6436", "n_xenium_cells"] == result.loc["6436", "n_xenium_warped_to_if"] == 4
    assert result.loc["6436", "n"] == 3
    assert result.loc["6436", "n_without_nucleus"] == 1
    assert result.loc["6436", "n_nucleated_cells"] == 3
    assert result.loc["6436", "nucleus_filter"] == "cell_id_in_xenium_nuclei_cell_ids"
    for marker in layer.MARKERS:
        assert result.loc["6436", f"n_measured_if_{marker}_before_nucleus_filter"] == 3
        assert result.loc["6436", f"n_measured_if_{marker}_excluded_without_nucleus"] == 1
    assert result.loc["6436", "n_measured_if_INS"] == result.loc["6436", "n_measured_if_GCG"] == 2
    assert result.loc[["112", "6442"], "state"].eq("unavailable").all()
    assert result.loc[["112", "6442"], "reason"].eq("pending_qupath_classifier_and_geojson_export").all()
    assert result.loc["6347", "reason"] == "no_restain_acquisition_record"
    assert result.loc[["112", "6442"], "n"].isna().all()
    path = builder_inputs["layers_dir"] / "restain/cells/6436.parquet"
    overlay = pq.ParquetFile(path).read().to_pandas().set_index("cell_id")
    assert set(overlay.frame) == {"if_um:6436:panc"}
    assert set(overlay.donor_id) == {"6436"}
    assert set(overlay.index) == {"cell-a", "cell-b", "cell-c"}
    # Refined affine AND displacement, not the seed's [101, 301] / [102, 302].
    assert overlay.loc["cell-a", ["x_um", "y_um", "feat__if_INS", "feat__if_GCG"]].tolist() == [107, 296, 0, 30]
    assert overlay.loc["cell-b", ["x_um", "y_um", "feat__if_INS", "feat__if_GCG"]].tolist() == [109, 299, 20, 40]
    assert overlay.loc["cell-c", ["x_um", "y_um"]].tolist() == [111, 302]
    assert overlay.transform_source.eq("refined").all()
    assert overlay.transform_direction.eq("xenium_um_to_if_um").all()
    assert overlay.transform_sha256.eq(result.loc["6436", "transform_sha256"]).all()
    assert overlay.transform_displacement_sha256.str.len().eq(64).all()
    assert overlay.pixel_size_um.eq(0.2505).all()
    assert overlay.producing_pixel_size_um.eq(0.2505).all()
    assert overlay.object_id.isna().all()
    assert overlay.object_id_state.eq("unavailable").all()
    assert overlay.object_id_reason.eq("no_qupath_object_id").all()
    for marker in layer.MARKERS:
        field = f"feat__if_{marker}"
        assert pd.isna(overlay.loc["cell-c", field])
        assert overlay.loc["cell-c", f"{field}_state"] == "unavailable"
        assert overlay.loc["cell-c", f"{field}_reason"] == "no_finite_if_measurement"
        assert overlay.loc["cell-a", f"{field}_state"] == "measured"
        assert overlay[f"hormone_positivity_{marker}"].isna().all()
        assert overlay[f"hormone_positivity_{marker}_state"].eq("withheld").all()
        assert overlay[f"hormone_positivity_{marker}_reason"].eq("auroc_below_preregistered_bar").all()


def test_calibration_mismatch_records_both_scales_and_retires_stale_overlay(builder_inputs, forbid_qptiff):
    layer.build_restain_layer(**builder_inputs)
    path = builder_inputs["restain_root"] / "valis/export-one/export.json"
    export = json.loads(path.read_text())
    export.update(if_pixel_um=0.25, if_um_per_px=1.0)
    path.write_text(json.dumps(export))
    path = builder_inputs["restain_root"] / "protein_rna_panc.json"
    protein = json.loads(path.read_text())
    protein[0]["protein_um_per_px"] = 0.5
    path.write_text(json.dumps(protein))
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6436"]
    assert row.state == row.calibration_state == "failed"
    assert "declared=0.2505" in row.reason and "producing=0.25" in row.reason
    assert "tolerance=5e-05" in row.reason and "coordinates_not_rescaled" in row.reason
    assert row.pixel_size_um == 0.2505
    assert row.producing_pixel_size_um == row.export_if_pixel_um == 0.25
    assert row.export_if_um_per_px == 1.0
    assert row.export_pixel_size_source == "fixture producer calibration"
    assert row.pixel_size_delta_um == pytest.approx(0.0005)
    assert row.transform_source == "refined"
    assert row.n_nucleated_cells == 3 and pd.isna(row.n)
    assert not (builder_inputs["layers_dir"] / "restain/cells/6436.parquet").exists()
    assert not forbid_qptiff


def test_legacy_export_scale_is_derived_from_producer_not_dapi_resolution(builder_inputs):
    path = builder_inputs["restain_root"] / "valis/export-one/export.json"
    export = json.loads(path.read_text())
    for field in ("if_pixel_um", "pixel_size_source", "if"):
        del export[field]
    path.write_text(json.dumps(export))
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6436"]
    assert row.state == "measured"
    assert row.export_if_um_per_px == 1.002
    assert pd.isna(row.export_if_pixel_um)
    assert row.export_pixel_size_source == "not_recorded_in_export"
    assert row.producing_pixel_size_um == 0.2505
    assert row.producing_pixel_size_source == "protein_rna_panc.json:protein_um_per_px/pyramid_factor"


@pytest.mark.parametrize("producing, state", [(0.25047, "measured"), (0.25045, "measured"),
                                             (0.25044, "failed"), (0.25, "failed"), (0.5, "failed")])
def test_calibration_tolerance_allows_only_declared_rounding(producing, state):
    result = layer.calibration_record({"if_pixel_um": producing, "if_um_per_px": 4 * producing}, 0.2505)
    assert result["calibration_state"] == state
    assert result["pixel_size_tolerance_um"] == 0.00005


def test_export_and_intensity_scales_cannot_disagree(builder_inputs):
    path = builder_inputs["restain_root"] / "protein_rna_panc.json"
    protein = json.loads(path.read_text())
    protein[0]["protein_um_per_px"] = 0.5
    path.write_text(json.dumps(protein))
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6436"]
    assert row.state == "failed"
    assert row.reason.startswith("inconsistent_producing_pixel_sizes:")
    assert not (builder_inputs["layers_dir"] / "restain/cells/6436.parquet").exists()


def test_cli_reports_failed_calibration_with_nonzero_exit(builder_inputs, monkeypatch, tmp_path, capsys):
    # Exercise the real CLI with its normal cells_xen location in a fixture repo.
    data = tmp_path / "core/phenocycler/data/integration"
    data.mkdir(parents=True)
    (data / "cells_xen").symlink_to(builder_inputs["cells_xen_root"], target_is_directory=True)
    monkeypatch.setattr(layer, "REPO", tmp_path)
    path = builder_inputs["restain_root"] / "valis/export-one/export.json"
    export = json.loads(path.read_text())
    export.update(if_pixel_um=0.25)
    path.write_text(json.dumps(export))
    args = [arg for name in ("restain_root", "sections_csv", "xenium_paths_csv", "layers_dir")
            for arg in ("--" + name.replace("_", "-"), str(builder_inputs[name]))]
    assert layer.main(args) == 1
    assert "FAILED donor 6436: inconsistent_producing_pixel_sizes" in capsys.readouterr().err


def add_registration_gate(inputs, donor="115", verdict="FAIL"):
    warp = inputs["restain_root"] / "valis" / ("export-one" if donor == "6436" else "failed-export")
    warp.mkdir(exist_ok=True)
    template = inputs["restain_root"] / "valis/export-one/export.json"
    export = json.loads(template.read_text())
    export.update(donor_id=donor)
    (warp / "export.json").write_text(json.dumps(export))
    gate = {"donor_id": donor, "roi": "panc", "serial": export["serial"], "alignment": {
        "transform_source": "refined", "verdict": verdict,
        "median_offset_um": 7.198144117833294, "verdict_why": "measured residual exceeds gate"
    }}
    (warp / "alignment_refined.json").write_text(json.dumps(gate))
    return warp


@pytest.mark.parametrize("donor", ["115", "6436"])
def test_registration_failure_precedes_missing_pairing_or_measured_values(builder_inputs, donor):
    warp = add_registration_gate(builder_inputs, donor)
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc[donor]
    assert row.state == "failed"
    assert row.reason.startswith("registration_gate_fail:") and "7.198 um" in row.reason
    assert row.registration_median_offset_um == 7.198144117833294
    assert row.source == str(warp / "alignment_refined.json")
    assert pd.isna(row.n)
    assert not (builder_inputs["layers_dir"] / f"restain/cells/{donor}.parquet").exists()


def test_registration_pass_does_not_claim_missing_intensity_is_measured(builder_inputs):
    add_registration_gate(builder_inputs, verdict="PASS")
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["115"]
    assert row.state == "unavailable"
    assert row.reason == "pending_qupath_classifier_and_geojson_export"


def test_acquired_donor_without_warp_is_pending_work(builder_inputs):
    path = builder_inputs["sections_csv"].with_name("restain_registration_inputs.csv")
    path.write_text("donor_id,roi,image_file,xenium_serial\n6347,panc,recorded.qptiff,recorded-serial\n")
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6347"]
    assert row.state == "unavailable"
    assert row.reason == "pending_qupath_classifier_and_geojson_export"
    assert row.source == str(path)
    assert row.image_file == "recorded.qptiff"


@pytest.mark.parametrize("defect, message", [
    ("warp_count", "n_xenium_warped_to_if"),
    ("bundle_count", "bundle cells"),
    ("bundle_duplicate", "duplicate bundle cell_id"),
    ("bundle_null", "missing or duplicate bundle cell_id"),
    ("geometry_unknown", "geometry cell_id keys disagree"),
    ("geometry_duplicate", "duplicate Xenium cell_id"),
    ("warp_order", "warp cell_id order disagrees"),
    ("nucleus_unknown", "invalid Xenium nucleus cell_id"),
])
def test_identifier_mismatches_refuse_every_write(builder_inputs, defect, message):
    warp = builder_inputs["restain_root"] / "valis/export-one"
    bundle = Path(json.loads((warp / "export.json").read_text())["xenium_bundle"])
    if defect == "warp_count":
        np.save(warp / "xenium_centroids_um.npy", np.ones((3, 2)))
    elif defect.startswith("bundle_"):
        ids = ["cell-a", "cell-b", "cell-c", "cell-d"]
        if defect == "bundle_count":
            ids.pop()
        elif defect == "bundle_duplicate":
            ids[-1] = "cell-a"
        else:
            ids[-1] = None
        pq.write_table(pa.table({"cell_id": ids}), bundle / "cells.parquet")
    elif defect.startswith("geometry_"):
        path = builder_inputs["cells_xen_root"] / "cells.parquet"
        cells = pq.ParquetFile(path).read().to_pandas()
        cells.loc[0, "cell_id"] = "unknown" if defect == "geometry_unknown" else "cell-a"
        pq.write_table(pa.Table.from_pandas(cells, preserve_index=False), path)
    elif defect == "warp_order":
        np.save(warp / "xenium_cell_ids.npy", np.array(["cell-b", "cell-a", "cell-c", "cell-d"]))
    else:
        np.save(warp / "xenium_nuclei_cell_ids.npy", np.array(["unknown"]))
    with pytest.raises((ValueError, AssertionError), match=message):
        layer.build_restain_layer(**builder_inputs)
    assert not builder_inputs["layers_dir"].exists()


def test_later_donor_mismatch_prevents_writing_earlier_donor(builder_inputs):
    restain = builder_inputs["restain_root"]
    warp = restain / "valis/export-two"
    shutil.copytree(restain / "valis/export-one", warp)
    export = json.loads((warp / "export.json").read_text())
    export.update(donor_id="6442", serial="0059947")
    (warp / "export.json").write_text(json.dumps(export))
    path = restain / "protein_rna_panc.json"
    records = json.loads(path.read_text())
    records.append({**records[0], "donor_id": "6442", "xenium_serial": "0059947"})
    path.write_text(json.dumps(records))
    source = builder_inputs["cells_xen_root"]
    cells = pq.ParquetFile(source / "cells.parquet").read().to_pandas()
    cells["donor_id"] = "6442"
    pq.write_table(pa.Table.from_pandas(cells, preserve_index=False), source / "second.parquet")
    np.save(warp / "xenium_centroids_um.npy", np.ones((3, 2)))
    with pytest.raises(AssertionError, match="6442: n_xenium_warped_to_if"):
        layer.build_restain_layer(**builder_inputs)
    assert not builder_inputs["layers_dir"].exists()


@pytest.mark.parametrize("missing", ["nuclei", "bundle", "refined.json", "refined_displacement.npz"])
def test_missing_input_is_explicitly_unavailable(builder_inputs, missing):
    warp = builder_inputs["restain_root"] / "valis/export-one"
    path = warp / "xenium_nuclei_cell_ids.npy"
    if missing == "bundle":
        path = Path(json.loads((warp / "export.json").read_text())["xenium_bundle"]) / "cells.parquet"
    elif missing.startswith("refined"):
        path = warp / missing
    path.unlink()
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6436"]
    assert row.state == "unavailable"
    assert row.reason == f"missing_restain_input:{path.name}"
    assert pd.isna(row.n)
    assert not (builder_inputs["layers_dir"] / "restain/cells/6436.parquet").exists()


@pytest.mark.parametrize("missing, reason", [
    ("nuclei", "no_recorded_xenium_nuclei"),
    ("intensities", "no_finite_if_measurements"),
])
def test_absent_measurements_never_produce_a_measured_empty_overlay(builder_inputs, missing, reason):
    if missing == "nuclei":
        path = builder_inputs["restain_root"] / "valis/export-one/xenium_nuclei_cell_ids.npy"
        np.save(path, np.array([], dtype=str))
    else:
        path = builder_inputs["cells_xen_root"] / "cells.parquet"
        cells = pq.ParquetFile(path).read().to_pandas()
        cells[["feat__if_INS", "feat__if_GCG"]] = np.nan
        pq.write_table(pa.Table.from_pandas(cells, preserve_index=False), path)
    row = layer.build_restain_layer(**builder_inputs).set_index("donor_id").loc["6436"]
    assert row.state == "unavailable"
    assert row.reason == reason
    assert pd.isna(row.n)
    assert not (builder_inputs["layers_dir"] / "restain/cells/6436.parquet").exists()


@pytest.mark.parametrize("bad", ["", "NaN", "0", "-0.1", "inf"])
def test_bad_csv_calibration_never_falls_back_to_image(declared_sections, forbid_qptiff, bad):
    declared_sections.write_text(declared_sections.read_text().replace(",0.2505", f",{bad}"))
    with pytest.raises(ValueError):
        layer.read_sections(declared_sections)
    assert not forbid_qptiff


def test_all_five_real_sections_declare_maintainer_pixel_size():
    path = REPO / "core/phenocycler/data/integration/restain_sections.csv"
    d = layer.read_sections(path)
    assert len(d) == 5
    assert d.pixel_size_um.eq(0.2505).all()


@pytest.mark.parametrize("shape", [(2, 2), (4, 2), (3, 3)])
def test_warp_count_and_shape_are_asserted(shape):
    with pytest.raises((AssertionError, ValueError)):
        layer.validate_warp_count(np.ones(shape), 3, "6436")


def test_nonfinite_warp_is_rejected():
    with pytest.raises(ValueError, match="nonfinite"):
        layer.validate_warp_count(np.array([[np.nan, 1.0]]), 1, "6436")


def test_recorded_hormone_validation_is_withheld_not_negative():
    record = {"preregistered_auroc_bar": 0.75,
              "INS": {"auroc": 0.5143, "clears_preregistered_bar": False}}
    result = layer.hormone_availability(record, "INS")
    assert result == {"state": "withheld", "reason": "auroc_below_preregistered_bar",
                      "auroc": 0.5143, "bar": 0.75}
    assert layer.hormone_availability(record, "GCG")["state"] == "unavailable"


def test_validation_does_not_invent_hormone_calls():
    record = {"preregistered_auroc_bar": 0.75,
              "INS": {"auroc": 0.9, "clears_preregistered_bar": True}}
    result = layer.hormone_availability(record, "INS")
    assert result["state"] == "unavailable"
    assert result["reason"] == "no_released_hormone_calls"


def test_source_inventory_uses_columns_and_handles_large_string(tmp_path):
    for name, donor in [("misleading-9999", "0112"), ("included", "6436"),
                        ("excluded-1", "06457"), ("excluded-2", "6579")]:
        p = tmp_path / name / "donor_id=9999" / "data.parquet"
        p.parent.mkdir(parents=True)
        pq.write_table(pa.table({"donor_id": pa.array([donor], type=pa.large_string()),
                                 "roi": ["panc"]}), p)
    inventory = layer.source_inventory(tmp_path)
    assert set(inventory) == {"112", "6436"}


def test_r_cards_parse_and_missing_artifact_is_a_named_note(tmp_path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    code = r'''
    source("apps/R/ieCore/R/ie_frames.R")
    source("apps/R/ieCore/R/ie_availability.R")
    source("apps/islet/R/card_x2_overlay.R")
    source("apps/islet/R/card_x3_identity.R")
    index <- ie_restain_index()
    stopifnot(index$availability$state == "unavailable")
    text <- as.character(ie_restain_note(index$availability, "x2_overlay"))
    stopifnot(grepl("x2_overlay", text), grepl("restain_artifact_not_built", text))
    a <- ie_attach_frame(data.frame(x_um = 1, y_um = 2), ie_frame("if_um", "6436", "panc"))
    b <- ie_attach_frame(data.frame(x_um = 1, y_um = 2), ie_frame("xenium_um", "6436", "panc"))
    stopifnot(inherits(try(ie_assert_same_frame(a, b), silent = TRUE), "try-error"))
    '''
    result = subprocess.run([rscript, "--vanilla", "-e", code], cwd=REPO, text=True,
                            capture_output=True,
                            env={**os.environ, "ISLET_LAYERS_DIR": str(tmp_path)})
    assert result.returncode == 0, result.stdout + result.stderr


def test_real_paired_warp_counts_and_nucleus_exclusions():
    root = REPO / "core/phenocycler/data/integration/cells_xen"
    if not root.is_dir():
        pytest.skip("external Xenium source data are not mounted")
    sections = layer.read_sections(REPO / "core/phenocycler/data/integration/restain_sections.csv")
    sections = sections.set_index("donor_id")
    for donor in ("6436", "6442"):
        paths = sorted(root.glob(f"donor_id={donor}/roi=panc/*.parquet"))
        warp_dir = REPO / f"data/development/restain/valis/donor_id={donor}"
        protein = layer.protein_records(REPO / "data/development/restain/protein_rna_panc.json")[donor]
        cells, warped, provenance = layer.inspect_donor(paths, warp_dir, donor, sections.loc[donor], protein)
        assert len(cells) == len(warped) > 0
        assert sections.loc[donor, "pixel_size_um"] == 0.2505
        assert "object_id" not in cells  # Do not silently fabricate QuPath UUIDs.
        ids = pd.Index(np.load(warp_dir / "xenium_cell_ids.npy", allow_pickle=True))
        nuclei = pd.Index(np.load(warp_dir / "xenium_nuclei_cell_ids.npy", allow_pickle=False))
        assert set(ids) == set(cells.cell_id)
        # Assembly, not just the input counts, must work on the real Xenium IDs.
        status = {**layer._status(donor, root), **provenance}
        overlay, details = layer.prepare_restain_cells(cells, warped, warp_dir, status)
        actual = overlay.to_pandas().set_index("cell_id")
        assert details["state"] == "measured"
        assert len(actual) == nuclei.nunique()
        assert details["n_without_nucleus"] == len(cells) - len(actual)
        np.testing.assert_array_equal(actual[["x_um", "y_um"]], warped[ids.get_indexer(actual.index)])
        bundle = Path(protein["xenium_bundle"])
        native = layer.read_source(bundle / "cells.parquet", ["x_centroid", "y_centroid"]).to_numpy()
        np.testing.assert_array_equal(warped, layer.Transform.load(warp_dir, stem="refined").apply(native))
        assert not np.allclose(warped, np.load(warp_dir / "warped_if_um.npy"))
        expected = cells.set_index("cell_id").loc[actual.index]
        for marker in layer.MARKERS:
            np.testing.assert_array_equal(actual[f"feat__if_{marker}"], expected[f"feat__if_{marker}"])
        assert actual.object_id.isna().all()
        assert actual.object_id_reason.eq("no_qupath_object_id").all()
        print(f"{donor}: n_xenium_warped_to_if = n_xenium_cells = {len(cells):,}; "
              f"{len(actual):,} nucleated cells retained; "
              f"{details['n_without_nucleus']:,} without recorded nuclei; "
              "CSV pixel size = 0.2505 um/px")


def test_r_cards_render_artifacts_and_refuse_wrong_frames(builder_inputs):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    add_registration_gate(builder_inputs)
    layer.build_restain_layer(**builder_inputs)
    layers = builder_inputs["layers_dir"]
    (layers / "availability").mkdir()
    pq.write_table(pa.table({"donor_id": ["6436", "112"], "capability": ["restain_if"] * 2,
                            "state": ["withheld", "failed"],
                            "reason": ["auroc_below_preregistered_bar", "registration failed: 7.198 um median offset"],
                            "n": pa.array([None, None], type=pa.int64())}),
                   layers / "availability/data.parquet")
    code = r'''
    source("apps/R/ieCore/R/ie_frames.R")
    source("apps/R/ieCore/R/ie_availability.R")
    source("apps/islet/R/card_x2_overlay.R")
    source("apps/islet/R/card_x3_identity.R")
    theme_islet <- function(...) ggplot2::theme_minimal()
    shiny::testServer(function(input, output, session) {
      card <- card_x2_overlay_server(input, output, session)
    }, {
      session$setInputs(x2_overlay_donor = "6436", x2_overlay_marker = "INS")
      stopifnot(ie_is_measured(card$drawing()$availability))
      stopifnot(inherits(card$drawing()$plot, "ggplot"))
      stopifnot(grepl("0.2505", card$drawing()$caption))
      stopifnot(grepl("2 measured; 1 unavailable", card$drawing()$caption))
      stopifnot(grepl("4 total Xenium cells; 3 retained; 1 excluded", card$drawing()$caption))
      stopifnot(grepl("3 IF measurements before filtering; 1 excluded", card$drawing()$caption))
      stopifnot(grepl("Nucleus filter: retain cells with a recorded Xenium nucleus", card$drawing()$caption))
      stopifnot(grepl("mapped with refined", card$drawing()$caption))
      stopifnot(length(ggplot2::ggplot_build(card$drawing()$plot)$data) == 2L)
      stopifnot(nrow(ggplot2::ggplot_build(card$drawing()$plot)$data[[2]]) == 2L)
      stopifnot(nzchar(output$x2_overlay_plot$src))
      session$setInputs(x2_overlay_donor = "112")
      stopifnot(card$drawing()$availability$state == "failed")
      stopifnot(grepl("7.198", output$x2_overlay_content$html))
      session$setInputs(x2_overlay_donor = "115")
      stopifnot(card$drawing()$availability$state == "failed")
      stopifnot(grepl("7.198", output$x2_overlay_content$html))
      session$setInputs(x2_overlay_donor = "6442")
      stopifnot(grepl("pending_qupath_classifier_and_geojson_export", output$x2_overlay_content$html))
      stopifnot(!grepl("dropped at QC", output$x2_overlay_content$html))
    })
    shiny::testServer(function(input, output, session) {
      card <- card_x3_identity_server(input, output, session)
    }, {
      session$setInputs(x3_identity_donor = "6436")
      d <- card$identity()$data
      stopifnot(all(grepl("Withheld", d$State)), identical(d$AUROC, c("0.5143", "0.5220")))
      stopifnot(grepl("x3_identity_table", output$x3_identity_content$html))
      session$setInputs(x3_identity_donor = "112")
      stopifnot(card$identity()$availability$state == "failed")
      stopifnot(grepl("7.198", output$x3_identity_content$html))
    })
    unlink(file.path(Sys.getenv("ISLET_LAYERS_DIR"), "availability"), recursive = TRUE)
    stopifnot(ie_is_measured(ie_restain_selection(ie_restain_index(), "6436")$availability))
    stopifnot(ie_restain_selection(ie_restain_index(), "115")$availability$state == "failed")
    row <- ie_restain_index()$data
    row <- row[row$donor_id == "6436", , drop = FALSE]
    row$hormone_positivity_INS_state <- "unavailable"
    row$hormone_positivity_INS_reason <- "no_hormone_validation"
    stopifnot(grepl("Unavailable", ie_restain_identity(row)$State[[1]]))
    stopifnot(!grepl("dropped at QC", ie_restain_identity(row)$State[[1]]))
    row <- ie_restain_index()$data
    row <- row[row$donor_id == "6436", , drop = FALSE]
    file <- file.path(ie_restain_root(), "cells", "6436.parquet")
    d <- as.data.frame(arrow::read_parquet(file))
    d$frame <- "xenium_um:6436:panc"
    arrow::write_parquet(d, file)
    stopifnot(inherits(try(ie_restain_cells(row), silent = TRUE), "try-error"))
    d$frame <- "if_um:6436:panc"
    d$feat__if_INS_state[is.na(d$feat__if_INS)] <- "measured"
    arrow::write_parquet(d, file)
    stopifnot(inherits(try(ie_restain_cells(row), silent = TRUE), "try-error"))
    '''
    result = subprocess.run([rscript, "--vanilla", "-e", code], cwd=REPO, text=True,
                            capture_output=True,
                            env={**os.environ, "ISLET_LAYERS_DIR": str(layers)})
    assert result.returncode == 0, result.stdout + result.stderr
