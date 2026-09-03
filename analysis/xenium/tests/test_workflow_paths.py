"""Fast contracts for local/nested Xenium workflow inputs."""

from pathlib import Path

import pytest

from scripts.workflow_paths import (
    require_existing_files,
    resolve_insulitis_jobs,
    resolve_phenotyping_jobs,
    sample_from_phenotyped_path,
)


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def test_explicit_nested_phenotyped_path_preserves_roi_directory(tmp_path):
    path = _touch(
        tmp_path / "processed" / "0059907" / "panc" / "0059907_phenotyped.h5ad"
    )
    assert resolve_phenotyping_jobs(
        samples=["ignored"], paths=[path], data_root=tmp_path
    ) == [("0059907", path)]


def test_legacy_sample_resolution_remains_flat(tmp_path):
    path = _touch(tmp_path / "processed" / "0041323" / "0041323_phenotyped.h5ad")
    assert resolve_phenotyping_jobs(
        samples=["0041323"], paths=None, data_root=tmp_path
    ) == [("0041323", path)]


def test_all_missing_inputs_are_reported_before_work_starts(tmp_path):
    one = tmp_path / "0059907_phenotyped.h5ad"
    two = tmp_path / "0059947_phenotyped.h5ad"
    with pytest.raises(FileNotFoundError) as excinfo:
        resolve_phenotyping_jobs(
            samples=[], paths=[one, two], data_root=tmp_path
        )
    message = str(excinfo.value)
    assert str(one) in message
    assert str(two) in message


def test_explicit_path_requires_exact_phenotyped_filename_contract(tmp_path):
    path = _touch(tmp_path / "0059907.h5ad")
    with pytest.raises(ValueError, match=r"\*_phenotyped\.h5ad"):
        sample_from_phenotyped_path(path)


def test_duplicate_serials_are_rejected_even_in_different_roi_directories(tmp_path):
    panc = _touch(tmp_path / "0059907" / "panc" / "0059907_phenotyped.h5ad")
    pln = _touch(tmp_path / "0059907" / "pln_1" / "0059907_phenotyped.h5ad")
    with pytest.raises(ValueError, match="duplicate serial.*0059907"):
        resolve_phenotyping_jobs(
            samples=[], paths=[panc, pln], data_root=tmp_path
        )


def test_insulitis_resolves_immune_sidecar_next_to_nested_base(tmp_path):
    base = _touch(tmp_path / "0059947" / "panc" / "0059947_phenotyped.h5ad")
    immune = _touch(
        tmp_path / "0059947" / "panc" / "0059947_immune_phenotyped.h5ad"
    )
    assert resolve_insulitis_jobs(
        samples=[], paths=[base], data_root=tmp_path
    ) == [("0059947", base, immune)]


def test_insulitis_missing_immune_sidecar_is_an_error(tmp_path):
    base = _touch(tmp_path / "0059947" / "panc" / "0059947_phenotyped.h5ad")
    with pytest.raises(FileNotFoundError, match="0059947_immune_phenotyped.h5ad"):
        resolve_insulitis_jobs(samples=[], paths=[base], data_root=tmp_path)


def test_generic_preflight_accepts_existing_files(tmp_path):
    path = _touch(tmp_path / "arbitrary_input.h5ad")
    assert require_existing_files([path], context="test") == [path]
