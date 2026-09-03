"""Shared integration-contract and declared-vocabulary fixtures."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from phenocycler.config import PipelineConfig
from phenocycler.integration.contract import coerce_cell_table
from phenocycler.integration.phenospace import load_phenotype_space, load_release
from phenocycler.integration.vocab import load_legacy_space


@pytest.fixture
def make_contract_table():
    """One factory for synthetic geometry and identity tables.

    Identity fixtures declare complete evidence by default.  Tests exercising unavailable
    evidence opt into ``evidence=pd.NA`` explicitly, so a forgotten field cannot quietly
    change the donor-composition denominator.
    """

    def _make(
        *,
        n: int = 20,
        modality: str = "phenocycler",
        donor_id: str = "6436",
        roi: str = "panc",
        profile: str = "identity",
        release_id: str = "synthetic-test-release",
        lineage_native: str = "Endocrine",
        lineage_common: str = "Endocrine",
        evidence=True,
        **columns,
    ) -> pd.DataFrame:
        raw = pd.DataFrame({
            "cell_id": [f"{modality}-{donor_id}-{roi}-{i}" for i in range(n)],
            "x_um": np.linspace(0.0, 500.0, n),
            "y_um": np.linspace(0.0, 400.0, n),
        })
        if profile == "identity":
            raw["lineage_native"] = lineage_native
            raw["lineage_common"] = lineage_common
            raw["evidence_positive"] = pd.Series([evidence] * n, dtype="boolean")
        for name, value in columns.items():
            raw[name] = value
        return coerce_cell_table(
            raw,
            modality=modality,
            donor_id=donor_id,
            roi=roi,
            profile=profile,
            identity_release_id=release_id if profile == "identity" else "",
        )

    return _make


@pytest.fixture
def legacy_space():
    """The retired seven-class declaration, loaded explicitly for audit tests only."""

    return load_legacy_space()


@pytest.fixture
def v4_space(tmp_path):
    """A second producer space proving mappings are per-release, not module state."""

    path = tmp_path / "v4_space.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "space_id": "synthetic-v4-space",
        "producer_lane": "synthetic-v4",
        "labels": [
            {"label": "Endo-beta", "common_lineage": "Endocrine"},
            {"label": "Endo-vascular", "common_lineage": "Vascular"},
            {"label": "Unavailable", "common_resolution": "unavailable", "reason": "missing"},
        ],
    }))
    return load_phenotype_space(path)


@pytest.fixture
def synthetic_run(tmp_path):
    """A content-pinned expression run suitable for release-loader tests."""

    data = tmp_path / "data"
    manifest = data / "runs" / "run-synthetic" / "manifests" / "expression.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({"content_id": "expression-content-synthetic"}))
    return PipelineConfig(data_dir=data)


@pytest.fixture
def synthetic_release(synthetic_run, tmp_path, v4_space):
    """A minimal valid PhenoCycler release bound to ``synthetic_run``."""

    release_dir = tmp_path / "identity_release"
    release_dir.mkdir()
    (release_dir / "label_space.json").write_text(v4_space.source.read_text())
    (release_dir / "release.json").write_text(json.dumps({
        "schema_version": 1,
        "release_id": "synthetic-release-v1",
        "modality": "phenocycler",
        "producer": {"lane": "synthetic-v4", "commit": "abc123"},
        "label_space": "label_space.json",
        "calls_dir": "calls",
        "cell_key": "object_id",
        "donors": ["6436", "6442"],
        "inputs": [{
            "kind": "expression",
            "run_id": "run-synthetic",
            "content_id": "expression-content-synthetic",
            "path": "expression",
        }],
        "columns": {"label": "phenotype"},
    }))
    synthetic_run.phenocycler_identity_release = release_dir
    return load_release(synthetic_run, "phenocycler")
