"""The v3 galectin layer: rebuilt from typed cells, on the UniFORM cohort scale.

Exercises ``derive_counts_from_typed_cells`` and the two helpers it composes
(``donor_m2gal3_reference_level``, ``build_donor_galectin``) against synthetic typed-cells and
UniFORM-aligned fixtures -- no real ``data/`` tree is read.  ``M2Gal3``'s cohort ``scale`` is
built with floor 0 / span 1, so the "scaled" value equals the raw value exactly and every
expectation below can be hand-computed.

Covers: per-islet/per-band counts (including the core vs. expansion/"peri" split), the
Metabolically_active definition (galectin-positive and NOT Myeloid), null propagation at both the
donor level (M2Gal3 not acquired) and the islet-region level (no admitted cell in the region
carries a measured M2Gal3 value) while ``count_Macrophage`` stays measured either way, and
column-name stability against the same ``count_<class>`` / ``prop_<class>`` names correction J
used.

Run with
    cd core/phenocycler && PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest \\
      -p no:cacheprovider -q tests/current/test_galectin_typed_cells.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))

import build_galectin_layer as G  # noqa: E402

COUNT_NAMES = ("count_Galectin_positive", "count_Galectin_macrophage",
               "count_Macrophage", "count_Metabolically_active")
PROP_NAMES = tuple(c.replace("count_", "prop_") for c in COUNT_NAMES)


def _write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def _typed_cells_row(object_id, cell_class, context="islet", immune_definition="",
                     islet_key=None, cell_region_20um=None):
    return dict(object_id=object_id, cell_class=cell_class, context=context,
                immune_definition=immune_definition, islet_key=islet_key,
                cell_region_20um=cell_region_20um)


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Redirects ``typing/cells`` and the UniFORM ``ALIGNED`` tree to ``tmp_path``, and hands back
    a floor-0/span-1 scale for M2Gal3 so scaled == raw exactly."""
    monkeypatch.setattr(G.S, "OUT", tmp_path / "typing")
    monkeypatch.setattr(G.S, "ALIGNED", tmp_path / "aligned")
    scale = pd.DataFrame({"floor": [0.0], "span": [1.0]}, index=pd.Index(["M2Gal3"], name="marker"))

    def write_donor(donor: str, cells_rows: list[dict], aligned: dict[str, float] | None,
                    marker_present: bool = True) -> None:
        """``cells_rows`` -- typed-cells rows (no ``m2gal3`` column, that lives in ALIGNED).
        ``aligned`` -- object_id -> M2Gal3 raw value; a NaN value writes ``NaN`` on that object.
        ``marker_present=False`` writes an ALIGNED file with NO ``M2Gal3`` column at all, the
        "not acquired in this donor's panel" case."""
        cells = pd.DataFrame(cells_rows)
        _write_parquet(tmp_path / f"typing/cells/donor_id={donor}/data.parquet", cells)

        object_ids = cells["object_id"].tolist()
        if marker_present:
            values = [aligned.get(o, np.nan) for o in object_ids] if aligned else \
                [np.nan] * len(object_ids)
            frame = pd.DataFrame({"object_id": object_ids, "M2Gal3": values})
        else:
            frame = pd.DataFrame({"object_id": object_ids, "DAPI": [1.0] * len(object_ids)})
        _write_parquet(tmp_path / f"aligned/donor_id={donor}/data.parquet", frame)

    return scale, write_donor


# --------------------------------------------------------------------- reference level fitting


def test_reference_level_is_p99_of_the_tissue_context_population(env):
    scale, write_donor = env
    tissue_values = np.linspace(0.0, 0.5, 50)
    rows = [_typed_cells_row(f"t{i}", "Exocrine", context="tissue")
            for i in range(len(tissue_values))]
    write_donor("D1", rows, {f"t{i}": v for i, v in enumerate(tissue_values)})

    cells = G.load_typed_donor("D1")
    m2gal3, level = G.donor_m2gal3_reference_level(cells, scale)
    assert level == pytest.approx(np.percentile(tissue_values, 99), abs=1e-6)
    # `apply_scale` casts through float32 (correction J's own convention); compare loosely.
    expected = cells.merge(
        pd.DataFrame({"object_id": [f"t{i}" for i in range(50)], "M2Gal3": tissue_values}),
        on="object_id")["M2Gal3"].to_numpy()
    assert np.allclose(m2gal3, expected, atol=1e-6)


def test_no_tissue_context_cells_leaves_the_level_unavailable(env):
    scale, write_donor = env
    rows = [_typed_cells_row("a", "Beta", context="islet", islet_key="Islet_1",
                             cell_region_20um="core")]
    write_donor("D1", rows, {"a": 0.5})

    cells = G.load_typed_donor("D1")
    _, level = G.donor_m2gal3_reference_level(cells, scale)
    assert level is None


def test_marker_not_acquired_leaves_the_level_unavailable_without_raising(env):
    scale, write_donor = env
    rows = [_typed_cells_row(f"t{i}", "Exocrine", context="tissue") for i in range(10)]
    write_donor("D1", rows, aligned=None, marker_present=False)

    cells = G.load_typed_donor("D1")
    m2gal3, level = G.donor_m2gal3_reference_level(cells, scale)
    assert level is None
    assert np.isnan(m2gal3).all()


# ---------------------------------------------------------------------- per-islet/band counts


def _two_islet_fixture(write_donor, donor="D1"):
    """Islet_1 core (6 admitted + 1 excluded + 1 non-islet), Islet_1 peri (2 admitted),
    Islet_2 core (2 admitted, both with a missing M2Gal3 value -- region-level unavailable)."""
    rows = [
        # Islet_1, core -- c1..c6 admitted, c7 Excluded, c8 outside every islet
        _typed_cells_row("c1", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c2", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c3", "Immune", immune_definition="Imm-T-cell",
                         islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c4", "Beta", islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c5", "Immune", immune_definition="Imm-dendritic-cell",
                         islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c6", "Exocrine", islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c7", "Excluded", islet_key="Islet_1", cell_region_20um="core"),
        _typed_cells_row("c8", "Beta", islet_key=None, cell_region_20um="tissue"),
        # Islet_1, expansion/"peri"
        _typed_cells_row("p1", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_1", cell_region_20um="expansion"),
        _typed_cells_row("p2", "Beta", islet_key="Islet_1", cell_region_20um="expansion"),
        # Islet_2, core -- every admitted cell has a missing M2Gal3 value
        _typed_cells_row("e1", "Beta", islet_key="Islet_2", cell_region_20um="core"),
        _typed_cells_row("e2", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_2", cell_region_20um="core"),
        # tissue reference population -- p99 lands at 0.495, well between 0.01 and 0.99 below
        *[_typed_cells_row(f"t{i}", "Exocrine", context="tissue") for i in range(50)],
    ]
    aligned = {
        "c1": 0.99, "c2": 0.01, "c3": 0.99, "c4": 0.99, "c5": 0.99, "c6": np.nan,
        "c7": 0.99, "c8": 0.99,
        "p1": 0.99, "p2": 0.01,
        "e1": np.nan, "e2": np.nan,
        **{f"t{i}": v for i, v in enumerate(np.linspace(0.0, 0.5, 50))},
    }
    write_donor(donor, rows, aligned)


def test_counts_and_proportions_per_islet_and_band(env):
    scale, write_donor = env
    _two_islet_fixture(write_donor)

    out = G.derive_counts_from_typed_cells(["D1"], scale)
    out = out.set_index(["islet_key", "region"])

    core = out.loc[("Islet_1", "core")]
    assert core["n_cells_eligible"] == 6, "c7 (Excluded) and c8 (no islet) are not admitted"
    assert core["n_cells_galectin_evaluable"] == 5, "c6's M2Gal3 is missing"
    assert core["count_Galectin_positive"] == 4        # c1, c3, c4, c5
    assert core["count_Galectin_macrophage"] == 2       # c1, c5 -- positive AND Myeloid
    assert core["count_Metabolically_active"] == 2      # c3, c4 -- positive AND NOT Myeloid
    assert core["count_Macrophage"] == 3                # c1, c2, c5 -- Myeloid regardless of M2Gal3
    assert (core["count_Galectin_macrophage"] + core["count_Metabolically_active"]
           == core["count_Galectin_positive"]), "the split is exhaustive and disjoint"
    assert core["prop_Galectin_positive"] == pytest.approx(4 / 6)
    assert core["prop_Macrophage"] == pytest.approx(3 / 6)
    assert core[G.COUNT_STATE] == G.MEASURED
    assert pd.isna(core[G.COUNT_REASON])

    peri = out.loc[("Islet_1", "peri")]
    assert peri["n_cells_eligible"] == 2
    assert peri["count_Galectin_positive"] == 1          # p1
    assert peri["count_Galectin_macrophage"] == 1         # p1
    assert peri["count_Metabolically_active"] == 0
    assert peri["count_Macrophage"] == 1
    assert peri["prop_Metabolically_active"] == 0.0


def test_a_region_with_no_measured_m2gal3_is_unavailable_but_macrophage_is_not(env):
    scale, write_donor = env
    _two_islet_fixture(write_donor)

    out = G.derive_counts_from_typed_cells(["D1"], scale)
    row = out[(out["islet_key"] == "Islet_2") & (out["region"] == "core")].iloc[0]

    assert row["n_cells_eligible"] == 2
    assert row["n_cells_galectin_evaluable"] == 0
    for c in ("count_Galectin_positive", "count_Galectin_macrophage",
              "count_Metabolically_active"):
        assert pd.isna(row[c]), c
        assert pd.isna(row[c.replace("count_", "prop_")]), c
    # e2 is Myeloid; the count needs no M2Gal3 at all and stays a real, measured 1 -- not null.
    assert row["count_Macrophage"] == 1
    assert row["prop_Macrophage"] == pytest.approx(0.5)
    assert row[G.COUNT_STATE] == G.UNAVAILABLE
    reason = row[G.COUNT_REASON]
    assert "count_Galectin_positive" in reason and "count_Macrophage" not in reason
    assert "no admitted cell" in reason


def test_a_donor_whose_m2gal3_reference_cannot_be_fit_is_unavailable_throughout(env):
    """M2Gal3 was not acquired for this donor's panel: every islet-region's galectin counts are
    null, but a Myeloid cell still makes ``count_Macrophage`` a real, measured count -- computing
    it never touches M2Gal3."""
    scale, write_donor = env
    rows = [
        _typed_cells_row("m1", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_5", cell_region_20um="core"),
        _typed_cells_row("m2", "Beta", islet_key="Islet_5", cell_region_20um="core"),
        *[_typed_cells_row(f"t{i}", "Exocrine", context="tissue") for i in range(10)],
    ]
    write_donor("D2", rows, aligned=None, marker_present=False)

    out = G.derive_counts_from_typed_cells(["D2"], scale)
    row = out.iloc[0]
    assert row["n_cells_eligible"] == 2
    assert row["n_cells_galectin_evaluable"] == 0
    for c in ("count_Galectin_positive", "count_Galectin_macrophage",
              "count_Metabolically_active"):
        assert pd.isna(row[c]), c
    assert row["count_Macrophage"] == 1
    assert pd.isna(row["m2gal3_reference_level"])
    assert row[G.COUNT_STATE] == G.UNAVAILABLE
    assert "D2" in row[G.COUNT_REASON] and "reference level is unavailable" in row[G.COUNT_REASON]


def test_multiple_donors_do_not_leak_reference_levels_or_regions(env):
    """D1's fitted level and D2's unavailable one must each apply only to their own cells."""
    scale, write_donor = env
    _two_islet_fixture(write_donor, donor="D1")
    write_donor("D2", [
        _typed_cells_row("m1", "Immune", immune_definition="Imm-macrophage",
                         islet_key="Islet_5", cell_region_20um="core"),
        *[_typed_cells_row(f"t{i}", "Exocrine", context="tissue") for i in range(10)],
    ], aligned=None, marker_present=False)

    out = G.derive_counts_from_typed_cells(["D1", "D2"], scale)
    assert set(out["case_id"]) == {"D1", "D2"}
    d1_core = out[(out["case_id"] == "D1") & (out["islet_key"] == "Islet_1")
                 & (out["region"] == "core")].iloc[0]
    assert d1_core["count_Galectin_positive"] == 4
    d2_row = out[out["case_id"] == "D2"].iloc[0]
    assert pd.isna(d2_row["count_Galectin_positive"])
    assert d2_row["count_Macrophage"] == 1


# ---------------------------------------------------------------------------- column stability


def test_column_names_and_nullable_dtypes_match_correction_j(env):
    scale, write_donor = env
    _two_islet_fixture(write_donor)
    out = G.derive_counts_from_typed_cells(["D1"], scale)

    assert {c for c in out.columns if c.startswith("count_")} == set(COUNT_NAMES)
    assert {c for c in out.columns if c.startswith("prop_")} == set(PROP_NAMES)
    assert not G.COUNT_STATE.startswith(("count_", "prop_"))
    assert not G.COUNT_REASON.startswith(("count_", "prop_"))
    for c in COUNT_NAMES:
        assert str(out[c].dtype) == "Int64", c
    for p in PROP_NAMES:
        assert out[p].dtype == np.dtype("float64"), p
    assert set(out["region"].unique()) <= {"core", "peri"}
    assert (out["method"] == G.METHOD).all()
    assert G.METHOD == "typed-cells-uniform-v3"
