"""Public phenotype export stays closed while pair design is unverified."""

from __future__ import annotations

import importlib

import pytest

from phenocycler.config import PipelineConfig
from phenocycler.integration import export_pheno


@pytest.fixture
def cfg(tmp_path):
    return PipelineConfig(data_dir=tmp_path)


def test_donor_export_has_no_historical_opt_in(cfg):
    with pytest.raises(TypeError, match="unexpected keyword"):
        export_pheno.export_donor(
            "6436",
            cfg,
            allow_legacy_restore_gated_redsea_for_historical_comparison=True,
        )

    with pytest.raises(export_pheno.ExportError, match="no marker-pair") as exc:
        export_pheno.export_donor("6436", cfg)
    assert "no compatibility opt-in" in str(exc.value)


def test_cohort_export_refuses_before_section_discovery(cfg, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("section discovery must not run")

    monkeypatch.setattr(export_pheno, "_sections_in", forbidden)
    with pytest.raises(export_pheno.ExportError, match="unavailable"):
        export_pheno.run_export_pheno(cfg)


def test_cli_exposes_no_legacy_flag():
    with pytest.raises(SystemExit) as exc:
        export_pheno.main(
            ["--allow-legacy-restore-gated-redsea-for-historical-comparison"]
        )
    assert exc.value.code == 2


def test_normal_integration_pipeline_remains_fail_closed(cfg):
    from phenocycler.integration import pipeline

    with pytest.raises(export_pheno.ExportError, match="no marker-pair"):
        pipeline.run_pipeline(
            cfg,
            only=["export_pheno"],
            force=True,
            tissues=["pancreas"],
            use_images=False,
        )


def test_cached_historical_export_cannot_bypass_the_gate(cfg):
    from phenocycler.integration import pipeline

    cached = cfg.data_dir / "integration/cells_pheno/donor_id=6436"
    cached.mkdir(parents=True)
    (cached / "part.parquet").write_bytes(b"retired")

    with pytest.raises(export_pheno.ExportError, match="cached files ignored|never reusable"):
        pipeline.run_pipeline(
            cfg,
            only=["export_pheno"],
            force=False,
            tissues=["pancreas"],
            use_images=False,
        )


def test_downstream_only_cannot_bypass_the_missing_current_export(cfg):
    from phenocycler.integration import pipeline

    with pytest.raises(export_pheno.ExportError, match="donor"):
        pipeline.run_pipeline(
            cfg,
            only=["donor"],
            force=False,
            tissues=["pancreas"],
            use_images=False,
        )


def test_status_marks_cached_historical_export_blocked(cfg, capsys):
    from phenocycler.integration import pipeline

    (cfg.data_dir / "integration/cells_pheno/donor_id=6436").mkdir(parents=True)
    pipeline.status(cfg, tissues=["pancreas"])

    export_line = next(
        line for line in capsys.readouterr().out.splitlines() if "export_pheno" in line
    )
    assert "BLOCKED" in export_line
    assert "cached files ignored" in export_line


@pytest.mark.parametrize(
    ("module_name", "runner_name"),
    [
        ("postxen", "run_postxen"),
        ("vocab", "run_vocab"),
        ("structures", "run_structures"),
        ("register", "run_register"),
        ("match", "run_match"),
        ("cellmatch", "run_cellmatch"),
        ("sameslide", "run_sameslide"),
        ("grid", "run_grid"),
        ("donor", "run_donor"),
        ("crossmodal", "run_crossmodal"),
        ("qc", "run_qc"),
        ("figures", "run_figures"),
    ],
)
def test_public_pheno_dependent_runners_refuse_before_reading_cache(
    cfg, module_name, runner_name
):
    module = importlib.import_module(f"phenocycler.integration.{module_name}")
    runner = getattr(module, runner_name)

    with pytest.raises(export_pheno.ExportError, match="Cached files are ignored"):
        runner(cfg)
