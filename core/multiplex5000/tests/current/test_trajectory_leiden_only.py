"""`trajectory_islet.py --leiden-only` writes the clustering without fitting the ordering.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_trajectory_leiden_only.py -q

The layer driver's islet_leiden node runs it this way since 2026-09-24. The full run's pooled
ordering fails its own `pseudotime_is_not_mostly_donor` check on this cohort and exits 1, which
stopped every node after it; the Leiden table it writes first is identical to this one (measured
cell for cell on the real inputs on 2026-09-24).  Here a synthetic islet table stands in for
`load_islets`, and the root pick and the ordering raise if they are reached at all.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "core/multiplex5000/scripts/typing/trajectory_islet.py"


@pytest.fixture
def traj():
    spec = importlib.util.spec_from_file_location("trajectory_islet_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _islets(n: int = 90) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    classes = ("Beta_cell", "Alpha_cell", "Delta_cell", "Exocrine")
    frame = pd.DataFrame({
        "case_id": np.repeat([6436, 6442, 6534], n // 3),
        "islet_key": [f"Islet_{i}" for i in range(n)],
        "combined_islet_id": [f"islet-{i}" for i in range(n)],
        "donor_status": np.repeat(["ND", "T1D", "T1D"], n // 3),
        "age": 20.0, "gender": "F", "stain_group": "S0", "acq_run": "run",
    })
    for c in classes:
        frame[f"count_{c}"] = rng.integers(1, 40, n)
        frame[f"peri_count_{c}"] = rng.integers(0, 20, n)
    frame["total_cells_core"] = frame[[f"count_{c}" for c in classes]].sum(axis=1)
    frame["total_cells_peri"] = frame[[f"peri_count_{c}" for c in classes]].sum(axis=1)
    return frame


def test_leiden_only_writes_the_table_and_never_fits_the_ordering(traj, tmp_path, monkeypatch):
    def refuse(*_args, **_kwargs):
        raise AssertionError("--leiden-only reached the pooled ordering")

    monkeypatch.setattr(traj, "load_islets", lambda donors, proportions=None: _islets())
    monkeypatch.setattr(traj, "pick_root", refuse)
    monkeypatch.setattr(traj, "diffusion_pseudotime", refuse)
    monkeypatch.setattr(sys, "argv", ["trajectory_islet.py", "--leiden-only",
                                      "--out-dir", str(tmp_path), "--leiden-name", "leiden.parquet"])
    assert traj.main() == 0
    table = pd.read_parquet(tmp_path / "leiden.parquet")
    resolutions = [f"leiden_{r}" for r in traj.LEIDEN_RESOLUTIONS]
    assert list(table.columns) == ["case_id", "islet_key", *resolutions,
                                   "leiden_umap_1", "leiden_umap_2"]
    assert len(table) == 90
    assert table[resolutions].notna().all().all()      # every islet is above the 10-cell floor
    assert not list(tmp_path.glob("*.h5ad")), "--leiden-only wrote an H5AD"


def test_graph_policy_selects_regions_and_pools_gal3_before_alr(traj):
    classes = [
        "Beta_cell", "Alpha_cell", "Exocrine", "Blood_Vessel", "Unresolved",
        "Tissue", "Doublet_suspected", "Conflict_B_T",
        "CD68+CD163+Macrophage_LGALS3+CD4+",
        "CD68+CD163-Macrophage_LGALS3+CD4-",
        "CD68-CD163+Macrophage_LGALS3-CD4+",
    ]
    frame = pd.DataFrame({prefix + c: [i + 1., i + 2.]
                          for prefix in ("count_", "peri_count_")
                          for i, c in enumerate(classes)})
    denominator = np.array([10., 20.])
    core_names, core = traj.trajectory_composition_inputs(frame, classes, denominator, "core")
    peri_names, peri = traj.trajectory_composition_inputs(frame, classes, denominator, "peri")
    assert core_names == ["core_Beta_cell", "core_Alpha_cell", "core_Blood_Vessel",
                          "core_Macrophage_LGALS3+", "core_Macrophage_LGALS3-"]
    assert peri_names == ["peri_Exocrine", "peri_Blood_Vessel",
                          "peri_Macrophage_LGALS3+", "peri_Macrophage_LGALS3-"]
    expected = np.log((frame["count_" + classes[8]] + frame["count_" + classes[9]]
                       + traj.ALR_PSEUDOCOUNT) / (denominator + traj.ALR_PSEUDOCOUNT))
    np.testing.assert_allclose(core[:, core_names.index("core_Macrophage_LGALS3+")],
                               expected, rtol=1e-6)
    assert peri.shape == (2, len(peri_names))


def test_graph_policy_refuses_missing_measurements(traj):
    frame = pd.DataFrame({"count_Beta_cell": [4.], "count_Blood_Vessel": [np.nan]})
    with pytest.raises(ValueError, match="unavailable"):
        traj.trajectory_composition_inputs(frame, ["Beta_cell", "Blood_Vessel"],
                                           np.array([4.]), "core")
    with pytest.raises(ValueError, match="missing measured"):
        traj.trajectory_composition_inputs(frame, ["Exocrine"], np.array([4.]), "peri")


def test_graph_policy_refuses_overlapping_prepooled_macrophages(traj):
    labels = ["Macrophage_Gal3+", "CD68+CD163+Macrophage_LGALS3+CD4+"]
    frame = pd.DataFrame({"count_" + label: [1.] for label in labels})
    with pytest.raises(ValueError, match="overlapping"):
        traj.trajectory_composition_inputs(frame, labels, np.array([4.]), "core")


def test_supported_coarse_immune_and_ambiguous_senior_labels(traj):
    denominator = np.array([5.])
    labels = ["Beta_cell", "Immune", "Unresolved", "Endocrine_unresolved",
              "Endocrine_ambiguous", "Excluded", "class_unavailable"]
    frame = pd.DataFrame({"count_" + label: [2.] for label in labels})
    names, values = traj.trajectory_composition_inputs(frame, labels, denominator, "core")
    assert names == ["core_Beta_cell", "core_Immune"]
    np.testing.assert_allclose(values[0, 1], np.log(2.5 / 5.5))


def test_unknown_coarse_immune_subtype_overlap_is_refused(traj):
    labels = ["Beta_cell", "Immune", "B_cell"]
    frame = pd.DataFrame({"count_" + label: [0.] for label in labels})
    with pytest.raises(ValueError, match="unknown coverage"):
        traj.trajectory_composition_inputs(frame, labels, np.array([5.]), "core")
