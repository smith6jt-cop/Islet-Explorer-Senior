"""Acceptance and failure cases for the small islet size/composition layer.

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
        core/phenocycler/tests/current/test_islet_summary.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pytest
from shapely.geometry import box

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "core" / "phenocycler"))
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

from scripts.senior import build_islet_summary as builder  # noqa: E402


@pytest.fixture
def source(tmp_path):
    root = tmp_path / "typing"
    (root / "app").mkdir(parents=True)
    (root / "geometry").mkdir()
    geometry, composition, outlines = [], [], []
    cases = [
        ("0112", "Islet_1", 1, 1),
        ("0112", "Islet_2", 12, 10),
        ("6539", "Islet_1", 24, 24),
        ("6539", "Islet_2", 0, None),
    ]
    cases += [(f"0{donor}", "Islet_1", 9, 9) for donor in DONOR_EXCLUSIONS]
    for donor, key, n, denominator in cases:
        for region in ("core", "expansion"):
            geometry.append(
                dict(
                    donor_id=donor,
                    islet_key=key,
                    combined_islet_id=f"opaque-{donor}-{key}",
                    region=region,
                    region_um2=64.0,
                    n_cells=n if region == "core" else 100,
                    islet_diam_um=9.0,
                    centroid_x_um=20.0,
                    centroid_y_um=30.0,
                )
            )
        if denominator is not None:
            # Different donor spellings and deliberately opaque combined IDs force
            # the read normalization and explicit aggregate-key contract to matter.
            counts = [6, 4] if denominator == 10 else [denominator]
            for phenotype, count in zip(["Beta_cell", "Alpha_cell"], counts, strict=False):
                composition.append(
                    dict(
                        donor_id=builder.norm_donor(donor),
                        islet_key=key,
                        combined_islet_id=f"another-opaque-id-{donor}-{key}",
                        region="core",
                        phenotype=phenotype,
                        n_cells=count,
                        region_total=denominator,
                        prop=count / denominator,
                    )
                )
        if key == "Islet_1":
            outlines.append(
                dict(
                    donor_id=builder.norm_donor(donor),
                    islet_key=key,
                    endocrine_um2_from_outlines=16.0 if donor == "6539" else 1.0,
                    n_endocrine_outlines=1,
                    wkb=box(0, 0, 2, 2).wkb,
                )
            )
    pd.DataFrame(geometry).to_parquet(root / "app/islet_geometry.parquet", index=False)
    pd.DataFrame(composition).to_parquet(root / "app/islet_composition_long.parquet", index=False)
    pd.DataFrame(outlines).to_parquet(root / "geometry/outlines.parquet", index=False)
    # A directory scan that tries to parse every file as parquet would fail here.
    (root / "geometry/nonendocrine.geojson").write_text(
        '{"type":"FeatureCollection","features":[]}'
    )
    return root


def test_emitted_per_islet_counts_and_donor_groups_survive_serialization(source, tmp_path):
    output = builder.write_summary(builder.build_islet_summary(source), tmp_path / "layers")
    islets = pq.ParquetFile(output / "islets.parquet").read().to_pandas()
    donors = pq.ParquetFile(output / "donors.parquet").read().to_pandas()
    counts = islets.set_index(builder.KEY)["n_cells"].to_dict()
    assert counts == {
        ("112", "Islet_1"): 1,
        ("112", "Islet_2"): 12,
        ("6539", "Islet_1"): 24,
        ("6539", "Islet_2"): 0,
    }
    assert islets["region"].eq("core").all()
    assert islets["state"].eq("measured").all()  # measured zero is retained
    assert set(donors["donor_id"]) == {"112", "6539"}
    assert donors.groupby("donor_id")["n_islets"].sum().to_dict() == {"112": 2, "6539": 2}
    assert donors.groupby("donor_id")["n_cells"].sum().to_dict() == {"112": 13, "6539": 24}
    assert set(donors["size_stratum"]) == set(builder.STRATA)
    assert islets.set_index(builder.KEY).loc[("112", "Islet_2"), "size_stratum"] == ">10 cells"
    assert not set(islets["donor_id"]) & DONOR_EXCLUSIONS.keys()
    schema = pq.ParquetFile(output / "islets.parquet").schema_arrow
    assert pa.types.is_string(schema.field("donor_id").type) or pa.types.is_large_string(
        schema.field("donor_id").type
    )
    assert schema.field("n_cells").type == pa.int64()
    # The converted outline is measured and then dropped: 147.3 of the file's 148.7 MB, no reader.
    assert "wkb" not in schema.names


def test_composition_keeps_its_denominator_and_unavailable_is_not_zero(source):
    tables = builder.build_islet_summary(source)
    islets = tables["islets"].set_index(builder.KEY)
    comp = tables["composition"]
    different = comp[(comp.donor_id == "112") & (comp.islet_key == "Islet_2")]
    assert different["n_cells"].tolist() == [12, 12]
    assert different["region_total"].tolist() == [10, 10]
    assert different["phenotype_n_cells"].sum() == 10
    assert different["prop"].sum() == pytest.approx(1)
    assert set(comp["composition_kind"]) == {"phenotype_counts"}
    no_composition = islets.loc[("6539", "Islet_2")]
    assert no_composition["n_cells"] == 0
    assert no_composition["composition_state"] == "unavailable"
    assert no_composition["composition_reason"] == "no_core_composition_exported"
    assert pd.isna(no_composition["region_total"])
    assert islets["positivity_state"].eq("unavailable").all()
    assert islets["positivity_reason"].str.contains("not_thresholded").all()


def test_real_outline_units_and_kind_are_preserved_without_fabricating_a_border(source):
    islets = builder.build_islet_summary(source)["islets"].set_index(builder.KEY)
    for donor, expected_area, expected_scale in [("112", 1.0, 0.5), ("6539", 16.0, 2.0)]:
        row = islets.loc[(donor, "Islet_1")]
        assert row["outline_state"] == "measured"
        assert row["outline_kind"] == "endocrine_cell_union"
        assert row["frame"] == f"phenocycler_um:{donor}:panc"
        assert row["outline_um_per_px"] == expected_scale
        assert row["endocrine_um2_from_outlines"] == pytest.approx(expected_area)
        assert row["endocrine_solidity"] == pytest.approx(1)
        assert row["endocrine_circularity"] == pytest.approx(np.pi / 4)
    missing = islets.loc[("112", "Islet_2")]
    assert missing["outline_state"] == "unavailable"
    assert missing["outline_reason"] == "no_endocrine_outline_exported"
    assert pd.isna(missing["outline_um_per_px"])
    assert pd.isna(missing["endocrine_circularity"])
    # The circularity and solidity above are measured on the converted outline, so they still
    # prove the conversion; the outline itself is not carried out of the builder.
    assert "wkb" not in islets.columns


def test_all_reads_disable_partition_discovery_and_stay_in_the_small_sources(source, monkeypatch):
    original = builder.ds.dataset
    calls = []

    def recording(paths, **kwargs):
        assert "partitioning" in kwargs and kwargs["partitioning"] is None
        calls.extend(Path(p).relative_to(source) for p in paths)
        return original(paths, **kwargs)

    monkeypatch.setattr(builder.ds, "dataset", recording)
    builder.build_islet_summary(source)
    assert set(calls) == {
        Path("app/islet_geometry.parquet"),
        Path("app/islet_composition_long.parquet"),
        Path("geometry/outlines.parquet"),
    }


@pytest.mark.parametrize("filename", ["islet_hormone_composition.parquet", "renamed.parquet"])
def test_refuses_degenerate_hormone_composition_by_path_and_by_content(source, tmp_path, filename):
    bad = source / "app" / filename
    pd.DataFrame(
        {"Ins_any": [4], "Ins_single": [4], "Multi_Pos": [0], "Triple_Neg": [0]}
    ).to_parquet(bad)
    with pytest.raises(ValueError, match="Refusing degenerate.*7,634.*Argmax"):
        builder.main(
            [
                "--source-root",
                str(source),
                "--composition-source",
                str(bad),
                "--layers-dir",
                str(tmp_path / "layers"),
            ]
        )
    assert not (tmp_path / "layers").exists()


def test_missing_sources_emit_reasons_and_remove_stale_measurements(source, tmp_path):
    layers = tmp_path / "layers"
    builder.write_summary(builder.build_islet_summary(source), layers)
    (source / "app/islet_composition_long.parquet").unlink()
    (source / "geometry/outlines.parquet").unlink()
    tables = builder.build_islet_summary(source)
    target = builder.write_summary(tables, layers)
    assert "composition" not in tables
    assert not (target / "composition.parquet").exists()
    islets = pq.ParquetFile(target / "islets.parquet").read().to_pandas()
    assert islets["composition_state"].eq("unavailable").all()
    assert islets["composition_reason"].str.startswith("source_missing:").all()
    assert islets["outline_state"].eq("unavailable").all()
    assert islets["outline_reason"].str.startswith("no_parquet_files:").all()
    assert "wkb" not in islets.columns
    assert islets["region_total"].isna().all()
    (source / "app/islet_geometry.parquet").unlink()
    tables = builder.build_islet_summary(source)
    builder.write_summary(tables, layers)
    assert set(tables) == {"availability"}
    assert tables["availability"]["state"].eq("unavailable").all()
    assert tables["availability"]["n"].isna().all()
    assert not (target / "islets.parquet").exists()
    assert not (target / "donors.parquet").exists()


def test_missing_cell_count_has_a_separate_state_and_is_not_a_zero_stratum(source):
    path = source / "app/islet_geometry.parquet"
    geometry = pq.ParquetFile(path).read().to_pandas()
    geometry.loc[(geometry.donor_id == "0112") & (geometry.islet_key == "Islet_1"), "n_cells"] = (
        None
    )
    geometry.to_parquet(path, index=False)
    tables = builder.build_islet_summary(source)
    row = tables["islets"].set_index(builder.KEY).loc[("112", "Islet_1")]
    assert row["state"] == "unavailable"
    assert row["reason"] == "missing_or_invalid_core_n_cells"
    assert row["size_stratum"] == "unavailable"
    assert pd.isna(row["n_cells"])
    composition = tables["composition"]
    absent = composition[(composition.donor_id == "112") & (composition.islet_key == "Islet_1")]
    assert absent["n_cells_state"].eq("unavailable").all()
    assert absent["n_cells_reason"].eq("missing_or_invalid_core_n_cells").all()
    donors = tables["donors"]
    unknown = donors[(donors.donor_id == "112") & (donors.size_stratum == "unavailable")].iloc[0]
    assert unknown["n_islets"] == 1
    assert pd.isna(unknown["n_cells"])


@pytest.mark.parametrize("problem", ["duplicate", "bad_fraction", "bad_total", "orphan"])
def test_invalid_aggregate_keys_and_denominators_fail_closed(source, problem):
    path = source / "app/islet_composition_long.parquet"
    comp = pq.ParquetFile(path).read().to_pandas()
    if problem == "duplicate":
        comp = pd.concat([comp, comp.iloc[[0]]], ignore_index=True)
    elif problem == "bad_fraction":
        comp.loc[0, "prop"] = 0.3
    elif problem == "bad_total":
        comp.loc[0, "region_total"] = 2
        comp.loc[0, "prop"] = 0.5
    else:
        comp.loc[0, "islet_key"] = "Islet_absent"
    comp.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="core composition"):
        builder.build_islet_summary(source)


def test_cli_honors_layer_root_environment(source, tmp_path, monkeypatch):
    monkeypatch.setenv("ISLET_LAYERS_DIR", str(tmp_path / "env_layers"))
    assert builder.main(["--source-root", str(source)]) == 0
    assert (tmp_path / "env_layers/islet_summary/islets.parquet").is_file()


def test_card_renders_sizes_before_composition_and_keeps_donors_separate(source, tmp_path):
    rscript = Path("/opt/R/4.5.1/bin/Rscript")
    if not rscript.exists():
        found = shutil.which("Rscript")
        if not found:
            pytest.skip("Rscript is not installed")
        rscript = Path(found)
    layers = tmp_path / "layers"
    builder.write_summary(builder.build_islet_summary(source), layers)
    script = r"""
    args <- commandArgs(TRUE)
    source("apps/R/ieCore/R/ie_availability.R")
    source("apps/R/ieCore/R/ie_frames.R")
    source("apps/islet/R/card_s1_size_composition.R")
    for (p in c("shiny", "ggplot2", "arrow", "scales")) {
      if (!requireNamespace(p, quietly = TRUE)) stop("Missing test dependency: ", p)
    }
    # Read the real navigation keys and the real card registry without loading the
    # app's data and packages.
    navigation <- new.env(parent = baseenv())
    navigation$SHOW_INFLAMMATION_TAB <- FALSE
    take <- function(file, names) {
      for (expr in parse(file)) {
        if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
            is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names) {
          eval(expr, navigation)
        }
      }
    }
    take("apps/islet/R/00_globals.R", c("IE_PIPELINE", "IE_REFERENCE", "IE_TAB_KEYS"))
    take("apps/islet/R/mod_crossmodal_ui.R", c("IE_CROSSMODAL_TAB", "IE_CROSSMODAL_CARDS"))
    # The card does not guess its own tab key -- the registry that mounts it supplies one.
    # Take the key from that registry instead of restating a tab name here, which is how
    # this card previously shipped pointing at a tab it was not mounted on.
    mount <- navigation$IE_CROSSMODAL_TAB
    others <- setdiff(navigation$IE_TAB_KEYS, mount)
    stopifnot("s1_size_composition" %in%
                vapply(navigation$IE_CROSSMODAL_CARDS, `[[`, "", "id"),
              mount %in% navigation$IE_TAB_KEYS,
              length(others) > 0L)
    missing <- card_s1_size_composition_read(file.path(args[1], "absent"), args[2])
    stopifnot(inherits(missing$note, "ie_guard_note"),
              grepl("islet_summary_not_built", format(missing$note)),
              grepl("s1_size_composition", format(missing$note)))
    wrong_domain <- card_s1_size_composition_read(args[1], file.path(args[2], "other"))
    stopifnot(grepl("different_domain", format(wrong_domain$note)))
    shiny::testServer(function(input, output, session) {
      card <- card_s1_size_composition_server(input, output, session,
                    active_tab = shiny::reactive(input$tabs), tab_key = mount,
                    layers_dir = args[1], active_domain_root = args[2])
    }, {
      session$setInputs(tabs = mount, s1_stratum = "large")
      stopifnot(card$active())
      tryCatch(card$layer(), shiny.silent.error = function(e) {
        stop("S1 must render when its own tab is active: ", mount)
      })
      stopifnot(nrow(card$sizes()) == 2L, all(card$sizes()$donor_id == "112"),
                identical(ie_frame_of(card$sizes()), "none"),
                all(card$composition()$n_cells == 12))
      before <- card$sizes()
      html <- as.character(output$s1_body$html)
      stopifnot(grepl('<option value="112" selected>', html, fixed = TRUE))
      stopifnot(regexpr("Size distribution", html)[1] < regexpr("Composition within", html)[1])
      stopifnot(inherits(card$size_plot(), "ggplot"), inherits(card$composition_plot(), "ggplot"))
      ggplot2::ggplot_build(card$size_plot())
      ggplot2::ggplot_build(card$composition_plot())
      session$setInputs(s1_stratum = "small")
      stopifnot(identical(card$sizes(), before), all(card$composition()$n_cells == 1))
      session$setInputs(s1_donor = "6539", s1_stratum = "all")
      # A donor change updates plots without rebuilding the controls.
      stopifnot(identical(as.character(output$s1_body$html), html))
      stopifnot(all(card$sizes()$donor_id == "6539"),
                identical(sort(card$sizes()$n_cells), c(0L, 24L)),
                all(card$composition()$donor_id == "6539"))
      session$setInputs(s1_stratum = "small")
      stopifnot(nrow(card$composition()) == 0L,
                grepl("no_composition_in_selected_size_stratum", output$s1_composition_view$html))
      for (other_tab in others) {
        session$setInputs(tabs = other_tab)
        stopifnot(!card$active(),
                  inherits(tryCatch(card$layer(), error = identity), "shiny.silent.error"))
        session$setInputs(tabs = mount)
        stopifnot(grepl('<option value="6539" selected>', output$s1_body$html, fixed = TRUE),
                  all(card$sizes()$donor_id == "6539"),
                  identical(card$size_plot()$labels$title, "Donor 6539"))
      }
      # Stale input from an artifact that no longer contains a donor falls back
      # to the same available donor in both the rebuilt menu and the plot.
      session$setInputs(s1_donor = "not-an-exported-donor", tabs = others[[1L]])
      session$setInputs(tabs = mount)
      stopifnot(grepl('<option value="112" selected>', output$s1_body$html, fixed = TRUE),
                all(card$sizes()$donor_id == "112"))
    })
    cat("S1 Shiny acceptance passed\n")
    """
    completed = subprocess.run(
        [str(rscript), "--vanilla", "-e", script, str(layers), str(source / "app")],
        cwd=REPO,
        env={**os.environ, "R_DEFAULT_DEVICE": "pdf"},
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "S1 Shiny acceptance passed" in completed.stdout


@pytest.mark.data
def test_real_sources_emit_every_core_islet_without_pooling_and_refuse_real_hormones(tmp_path):
    root = builder.DEFAULT_SOURCE
    if not (root / "app/islet_geometry.parquet").is_file():
        pytest.skip("real typing aggregates not present")
    tables = builder.build_islet_summary(root)
    target = builder.write_summary(tables, tmp_path / "real_layers")
    actual = pq.ParquetFile(target / "islets.parquet").read().to_pandas()
    source = ds.dataset(root / "app/islet_geometry.parquet", format="parquet", partitioning=None)
    expected = source.to_table().to_pandas().query("region == 'core'").copy()
    expected["donor_id"] = expected["donor_id"].map(builder.norm_donor)
    expected = expected[
        expected.donor_id.isin(builder.filter_eligible_donors(expected.donor_id.unique()))
    ]
    expected = expected.sort_values(builder.KEY).reset_index(drop=True)
    pd.testing.assert_frame_equal(actual[builder.KEY], expected[builder.KEY])
    assert actual.n_cells.tolist() == expected.n_cells.tolist()
    assert tables["donors"].groupby("donor_id").n_islets.sum().to_dict() == (
        expected.groupby("donor_id").size().to_dict()
    )
    assert actual.loc[actual.outline_state == "unavailable", "outline_reason"].str.len().gt(0).all()
    assert actual.loc[actual.composition_state == "unavailable", "region_total"].isna().all()
    hormone_path = root / "app/islet_hormone_composition.parquet"
    hormones = ds.dataset(hormone_path, format="parquet", partitioning=None).to_table().to_pandas()
    assert hormones.Ins_any.eq(hormones.Ins_single).all()
    assert hormones[["Multi_Pos", "Triple_Neg"]].eq(0).all().all()
    with pytest.raises(ValueError, match="Refusing degenerate"):
        builder.build_islet_summary(root, composition_source=hormone_path)
