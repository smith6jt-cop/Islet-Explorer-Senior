"""Geometry export: a contract table from ``data/expression``, carrying no identity.

Synthetic tree throughout — these tests must not depend on the real run data.
"""

from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest

from phenocycler.config import PipelineConfig
from phenocycler.integration import export_pheno
from phenocycler.integration.contract import feature_columns, has_identity

MARKERS = ("INS", "GCG", "SST", "CD3e", "CD31")


def _make_run(tmp_path, run_id=None, donor="6539", n=6, *, with_area=True,
              eligible=None, drop_area_for=()):
    """Write a stage tree. Flat by default; ``run_id`` builds a retired run directory."""
    root = tmp_path if run_id is None else tmp_path / "runs" / run_id
    ids = [f"uuid-{i}" for i in range(n)]
    if eligible is None:
        eligible = [True] * n

    expr = pd.DataFrame({
        "donor_id": donor,
        "object_id": ids,
        "image": f"{donor}_Scan1.er.qptiff - resolution #1",
        "X_centroid": np.arange(n, dtype="float64") * 10.0,
        "Y_centroid": np.arange(n, dtype="float64") * 20.0,
        "cell_region": ["core"] * (n // 2) + ["tissue"] * (n - n // 2),
        "islet_num": ["7"] * (n // 2) + [""] * (n - n // 2),
        "qc_analysis_eligible": eligible,
        "qc_estimation_eligible": [True] * n,
    })
    for j, m in enumerate(MARKERS):
        expr[m] = np.arange(n, dtype="float64") + j
    d = root / "expression" / f"donor_id={donor}"
    d.mkdir(parents=True)
    expr.to_parquet(d / "data.parquet", index=False)

    if with_area:
        cells = pd.DataFrame({
            "object_id": ids,
            "cell_area": np.arange(n, dtype="float64") + 100.0,
        })
        cells = cells[~cells["object_id"].isin(drop_area_for)]
        c = root / "cells" / f"donor_id={donor}"
        c.mkdir(parents=True)
        cells.to_parquet(c / "data.parquet", index=False)

    if run_id is not None:
        (tmp_path / "runs" / "LATEST").write_text(run_id)
    (tmp_path / "qupath_manifest.json").write_text(json.dumps({
        "images": [{
            "donor_id": donor,
            "image_id": f"{donor}_Scan1.er.qptiff - resolution #1",
        }]
    }))
    return root


@pytest.fixture
def cfg(tmp_path):
    overrides = tmp_path / "section_overrides.csv"
    overrides.write_text("image_id,donor_id,roi,slide,reason\n")
    return PipelineConfig(
        data_dir=tmp_path,
        qupath_manifest=tmp_path / "qupath_manifest.json",
        section_overrides_csv=overrides,
    )


def test_export_geometry_builds_a_contract_table_with_no_identity(cfg, tmp_path):
    _make_run(tmp_path)
    out = export_pheno.export_geometry(cfg, "6539")

    assert len(out) == 6
    assert not has_identity(out)
    assert (out["lineage_native"] == "").all()
    assert set(feature_columns(out)) == {"feat__CD3e", "feat__CD31"}


def test_object_id_is_the_join_key_and_survives_as_cell_id(cfg, tmp_path):
    _make_run(tmp_path)
    out = export_pheno.export_geometry(cfg, "6539")
    assert list(out["cell_id"]) == [f"uuid-{i}" for i in range(6)]


def test_area_comes_from_the_cells_stage_joined_on_object_id(cfg, tmp_path):
    _make_run(tmp_path)
    out = export_pheno.export_geometry(cfg, "6539").set_index("cell_id")
    for i in range(6):
        assert out.loc[f"uuid-{i}", "area_um2"] == pytest.approx(100.0 + i)


def test_a_cell_with_no_area_is_null_not_zero(cfg, tmp_path):
    """A missing measurement is unavailable, not zero — the rule the whole contract rests on."""
    _make_run(tmp_path, drop_area_for=("uuid-2",))
    out = export_pheno.export_geometry(cfg, "6539").set_index("cell_id")
    assert pd.isna(out.loc["uuid-2", "area_um2"])
    assert out.loc["uuid-3", "area_um2"] == pytest.approx(103.0)


def test_islet_geometry_is_carried_but_never_as_identity(cfg, tmp_path):
    _make_run(tmp_path)
    out = export_pheno.export_geometry(cfg, "6539")
    assert set(out["region"]) == {"core", "tissue"}
    assert set(out.loc[out["region"] == "core", "qupath_islet_num"]) == {"7"}
    assert (out["structure_id"] == "").all()
    # islet_num is registration/audit context; S3 alone may create structure_id.
    assert not has_identity(out)


def test_qc_ineligible_cells_are_dropped_by_default(cfg, tmp_path):
    _make_run(tmp_path, eligible=[True, False, True, True, False, True])
    out = export_pheno.export_geometry(cfg, "6539")
    assert len(out) == 4
    assert "uuid-1" not in set(out["cell_id"])

    keep_all = export_pheno.export_geometry(cfg, "6539", analysis_eligible_only=False)
    assert len(keep_all) == 6


def test_no_run_id_is_stamped_on_any_row(cfg, tmp_path):
    """Rows used to carry ``phenocycler_run_id``, set to the run directory's name.

    There is no run directory now, so the column would have said "data" on every row.
    Nothing outside these tests ever read it.
    """

    _make_run(tmp_path)
    out = export_pheno.export_geometry(cfg, "6539")
    assert "phenocycler_run_id" not in out.columns


def test_the_flat_tree_is_used_when_no_run_id_is_given(cfg, tmp_path):
    _make_run(tmp_path)
    assert len(export_pheno.export_geometry(cfg, "6539")) == 6


def test_a_retired_run_directory_is_still_readable_by_name(cfg, tmp_path):
    """The pre-2026-09-16 cohort output is 76 GB and stays where it is."""

    _make_run(tmp_path, run_id="retired-run")
    assert len(export_pheno.export_geometry(cfg, "6539", run_id="retired-run")) == 6


def test_an_unknown_run_or_donor_fails_loudly(cfg, tmp_path):
    _make_run(tmp_path)
    with pytest.raises(export_pheno.ExportError, match="not found"):
        export_pheno.export_geometry(cfg, "6539", run_id="nope")
    with pytest.raises(export_pheno.ExportError, match="no expression parquet"):
        export_pheno.export_geometry(cfg, "9999")


def test_identity_export_is_still_refused(cfg, tmp_path):
    """Adding a geometry path must not open an identity path."""
    _make_run(tmp_path)
    with pytest.raises(export_pheno.ExportError, match="no marker-pair"):
        export_pheno.export_donor("6539", cfg)


def test_run_export_geometry_writes_one_partition_per_donor(cfg, tmp_path):
    _make_run(tmp_path)
    written = export_pheno.run_export_geometry(cfg, ["6539"])
    assert len(written) == 1 and written[0].exists()

    back = pd.read_parquet(written[0])
    assert len(back) == 6
    assert not has_identity(back)


def test_run_export_geometry_refuses_a_donor_the_run_does_not_have(cfg, tmp_path):
    _make_run(tmp_path)
    with pytest.raises(export_pheno.ExportError, match="6457"):
        export_pheno.run_export_geometry(cfg, ["6457"])
