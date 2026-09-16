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


def test_builder_csv_calibration_and_if_overlay_acceptance(
    tmp_path, declared_sections, forbid_qptiff
):
    """Exercise the builder, including the measured absence of QuPath UUIDs.

    This acceptance case intentionally does not fabricate an object_id column
    absent from the real source. It cannot pass the brief's UUID-only join rule
    without a producer crosswalk or a clarified native-Xenium-ID exception.
    """
    source = tmp_path / "cells_xen"
    source.mkdir()
    pq.write_table(pa.table({"donor_id": ["06436"] * 2, "roi": ["panc"] * 2,
                            "cell_id": ["cell-b", "cell-a"],
                            "feat__if_INS": [20.0, 10.0],
                            "feat__if_GCG": [40.0, 30.0]}), source / "cells.parquet")
    restain = tmp_path / "inputs"
    warp = restain / "valis" / "export-one"
    warp.mkdir(parents=True)
    (warp / "export.json").write_text(json.dumps({
        "donor_id": "6436", "roi": "panc", "serial": "0059907"
    }))
    np.save(warp / "warped_if_um.npy", np.array([[101.0, 301.0], [102.0, 302.0]]))
    np.save(warp / "xenium_cell_ids.npy", np.array(["cell-a", "cell-b"]))
    np.save(warp / "xenium_nuclei_cell_ids.npy", np.array(["cell-a", "cell-b"]))
    paths = tmp_path / "xenium_paths.csv"
    paths.write_text("donor_id,roi\n6436,panc\n0112,panc\n6457,panc\n6579,panc\n")
    result = layer.build_restain_layer(cells_xen_root=source, restain_root=restain,
        sections_csv=declared_sections, xenium_paths_csv=paths, layers_dir=tmp_path / "layers")
    assert result.loc[result.donor_id.eq("6436"), "pixel_size_um"].item() == 0.2505
    assert not forbid_qptiff
    overlay = pq.ParquetFile(tmp_path / "layers/restain/cells/6436.parquet").read().to_pandas()
    assert set(overlay.frame) == {"if_um:6436:panc"}
    assert set(overlay.x_um) == {101.0, 102.0}


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
        # This worktree keeps tracked CSVs in data/, with the external tree
        # linked at data/data. Do not modify that shared setup.
        root = REPO / "core/phenocycler/data/data/integration/cells_xen"
    if not root.is_dir():
        pytest.skip("external Xenium source data are not mounted")
    sections = layer.read_sections(REPO / "core/phenocycler/data/integration/restain_sections.csv")
    sections = sections.set_index("donor_id")
    for donor in ("6436", "6442"):
        paths = sorted(root.glob(f"donor_id={donor}/roi=panc/*.parquet"))
        warp_dir = REPO / f"data/development/restain/valis/donor_id={donor}"
        cells, warped = layer.inspect_donor(paths, warp_dir, donor, sections.loc[donor])
        assert len(cells) == len(warped) > 0
        assert sections.loc[donor, "pixel_size_um"] == 0.2505
        assert "object_id" not in cells  # Do not silently fabricate QuPath UUIDs.
        ids = np.load(warp_dir / "xenium_cell_ids.npy", allow_pickle=True).astype(str)
        nuclei = np.load(warp_dir / "xenium_nuclei_cell_ids.npy", allow_pickle=False)
        assert set(ids) == set(cells.cell_id)
        print(f"{donor}: n_xenium_warped_to_if = n_xenium_cells = {len(cells):,}; "
              f"{np.count_nonzero(~np.isin(ids, nuclei)):,} without recorded nuclei; "
              "CSV pixel size = 0.2505 um/px")


def test_r_cards_render_artifacts_and_refuse_wrong_frames(tmp_path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    root = tmp_path / "restain"
    (root / "cells").mkdir(parents=True)
    row = layer._status("6436", tmp_path)
    row.update(state="measured", reason="", n=2, pixel_size_um=0.2505)
    for marker in layer.MARKERS:
        row[f"hormone_positivity_{marker}_state"] = "withheld"
        row[f"hormone_positivity_{marker}_reason"] = "auroc_below_preregistered_bar"
        row[f"hormone_positivity_{marker}_auroc"] = 0.5143
        row[f"hormone_positivity_{marker}_bar"] = 0.75
    pq.write_table(pa.Table.from_pylist([row], schema=layer._schema()), root / "availability.parquet")
    pq.write_table(pa.table({"donor_id": ["6436"] * 2, "cell_id": ["a", "b"],
                            "object_id": ["uuid-a", "uuid-b"], "x_um": [100.0, 200.0],
                            "y_um": [300.0, 400.0], "frame": ["if_um:6436:panc"] * 2,
                            "feat__if_INS": [1.0, 2.0], "feat__if_GCG": [3.0, 4.0]}),
                   root / "cells/6436.parquet")
    (tmp_path / "availability").mkdir()
    pq.write_table(pa.table({"donor_id": ["6436"], "capability": ["restain_if"],
                            "state": ["measured"], "reason": [""], "n": [2]}),
                   tmp_path / "availability/data.parquet")
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
      stopifnot(length(ggplot2::ggplot_build(card$drawing()$plot)$data) == 2L)
      session$setInputs(x2_overlay_donor = "112")
      stopifnot(!ie_is_measured(card$drawing()$availability))
    })
    shiny::testServer(function(input, output, session) {
      card <- card_x3_identity_server(input, output, session)
    }, {
      session$setInputs(x3_identity_donor = "6436")
      d <- card$identity()$data
      stopifnot(all(grepl("Withheld", d$State)), all(d$AUROC == "0.5143"))
    })
    row <- ie_restain_index()$data
    file <- file.path(ie_restain_root(), "cells", "6436.parquet")
    d <- as.data.frame(arrow::read_parquet(file))
    d$frame <- "xenium_um:6436:panc"
    arrow::write_parquet(d, file)
    stopifnot(inherits(try(ie_restain_cells(row), silent = TRUE), "try-error"))
    '''
    result = subprocess.run([rscript, "--vanilla", "-e", code], cwd=REPO, text=True,
                            capture_output=True,
                            env={**os.environ, "ISLET_LAYERS_DIR": str(tmp_path)})
    assert result.returncode == 0, result.stdout + result.stderr
