"""Contracts for the four-cell grid layer.

The layer presents per-donor composition for pancreas and pLN in protein and
RNA. Everything asserted here is a rule the layer would break silently if it
regressed: a hole printing as a zero, a never-measured nucleus counted as an
absent one, a cohort fact typed into a string, a distance crossing the tissue
boundary, or a donor the cohort module excludes reaching the roster.

Nothing here checks that a number matches a previous biological result. The one
quantitative check is internal coherence -- the buckets of a grid cell sum to
the cells that grid cell assessed -- which is a statement about the build, not
about the biology.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
BUILDER = REPO / "scripts" / "senior" / "build_grid_composition.py"
CARD = REPO / "apps" / "islet" / "R" / "card_g1_grid.R"


def _load():
    spec = importlib.util.spec_from_file_location("build_grid_composition", BUILDER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def grid():
    return _load()


@pytest.fixture(scope="module")
def layer(grid):
    """The built layer, when this checkout has one. Never built by the test."""
    out = Path(grid.OUT_DEFAULT)
    if not (out / "manifest.json").exists():
        pytest.skip(f"no grid_composition layer at {out}; run {BUILDER.name} first")
    return {
        "manifest": json.loads((out / "manifest.json").read_text()),
        "composition": pd.read_parquet(out / "composition.parquet"),
        "presence": pd.read_parquet(out / "presence.parquet"),
        "fine": pd.read_parquet(out / "fine.parquet"),
        "crossmodal": pd.read_csv(out / "crossmodal_tnk.csv", dtype={"donor_id": str}),
        "eta": pd.read_csv(out / "batch_eta_squared.csv"),
        "dir": out,
    }


# ---------------------------------------------------------------------------
# The exclusion list lives once.
# ---------------------------------------------------------------------------
def test_excluded_donor_cannot_reach_the_roster(grid):
    excluded = sorted(grid.DONOR_EXCLUSIONS)
    assert excluded, "the cohort module declares no exclusion; this test is meaningless"
    assert grid.eligible(excluded + ["9999"]) == ["9999"]


def test_builder_does_not_restate_the_exclusion_list(grid):
    """A lane that retyped the excluded ids would go stale the day cohort.py moved.

    6579 is exactly that: it was excluded, cohort.py reinstated it on 2026-09-17,
    and a lane holding its own copy would still be dropping it.
    """
    source = BUILDER.read_text()
    body = "\n".join(
        line for line in source.splitlines()
        if not line.lstrip().startswith("#") and "REINSTATED" not in line
    )
    for donor in grid.DONOR_EXCLUSIONS:
        assert f'"{donor}"' not in body and f"'{donor}'" not in body


# ---------------------------------------------------------------------------
# A missing measurement is unavailable -- never zero, never negative.
# ---------------------------------------------------------------------------
def test_absence_never_carries_a_count(layer):
    holes = layer["composition"][layer["composition"]["state"].ne("measured")]
    assert len(holes), "no absence in the layer; the rule would go untested"
    assert holes["n_cells"].isna().all()
    assert holes["frac_of_assessed"].isna().all()


def test_every_absence_names_a_reason(layer):
    for table in ("composition", "presence"):
        rows = layer[table][layer[table]["state"].ne("measured")]
        reasons = rows["reason"].fillna("")
        assert (reasons.str.len() > 0).all(), f"{table} carries an absence with no reason"


def test_states_are_the_five(grid, layer):
    for table in ("composition", "presence"):
        assert set(layer[table]["state"]) <= set(grid.STATES)


def test_lineage_absence_is_never_measured_and_always_explains(grid):
    crosswalk = grid.load_crosswalk()
    for _, row in crosswalk.iterrows():
        for tissue in ("pancreas", "pln"):
            for modality in ("protein", "RNA"):
                state, reason = grid.lineage_absence(row, tissue, modality)
                assert state in grid.STATES and state != grid.STATE_MEASURED
                assert reason and row["common"] in reason


# ---------------------------------------------------------------------------
# The nucleus rule, and the absence that is not it.
# ---------------------------------------------------------------------------
def _xenium_fixture(tmp_path: Path) -> Path:
    """Three cells: one called, one measured to have no nucleus, one never measured."""
    frame = pd.DataFrame(
        {
            "lineage": ["Myeloid", "Unavailable_no_nucleus", "Myeloid"],
            "common_lineage": ["Myeloid", "", "Myeloid"],
            "nucleus_count": [3, 0, None],
            "nucleus_measured": [True, True, False],
        }
    )
    path = tmp_path / "data.parquet"
    frame.to_parquet(path, index=False)
    return path


def test_nucleus_never_measured_is_not_a_nucleus_that_is_absent(grid, tmp_path):
    broad, _fine, no_nucleus, extras = grid.read_xenium([_xenium_fixture(tmp_path)])
    assert no_nucleus == 1, "the measured no-nucleus cell must be excluded"
    assert extras["nucleus_not_measured"] == 1, (
        "a cell whose nucleus was never measured is a missing measurement, not a "
        "measured absence of a nucleus"
    )
    assert sum(broad.values()) == 1, "both excluded cells must leave the denominator"


def test_a_partition_with_no_nucleus_measurement_is_unavailable_not_empty(layer):
    """Measured in this release: three pancreas RNA donors carry nucleus_measured
    False on every cell. That is unavailable with a named reason, not a donor with
    zero cells."""
    presence = layer["presence"]
    rows = presence[presence["n_nucleus_not_measured"].fillna(0) > 0]
    if not len(rows):
        pytest.skip("no partition in this checkout has an unmeasured nucleus")
    assert (rows["state"] == "unavailable").all()
    assert rows["reason"].str.contains("never measured").all()
    assert (rows["n_cells_assessed"].fillna(0) == 0).all()


# ---------------------------------------------------------------------------
# Internal coherence: the subsets sum to their parent.
# ---------------------------------------------------------------------------
def test_buckets_sum_to_the_cells_assessed(layer):
    comp = layer["composition"]
    measured = comp[comp["state"].eq("measured")]
    summed = measured.groupby(["tissue", "modality", "donor_id"])["n_cells"].sum()
    declared = (
        layer["presence"][layer["presence"]["state"].eq("measured")]
        .set_index(["tissue", "modality", "donor_id"])["n_cells_assessed"]
    )
    joined = pd.concat([summed.rename("s"), declared.rename("d")], axis=1).dropna()
    assert len(joined) == len(declared)
    assert (joined["s"] == joined["d"]).all()
    assert layer["manifest"]["coherence"]["n_mismatched"] == 0


def test_every_measured_row_has_a_denominator_it_divides(layer):
    measured = layer["composition"][layer["composition"]["state"].eq("measured")]
    frac = measured["n_cells"] / measured["denom_assessed"]
    assert np.allclose(frac, measured["frac_of_assessed"])
    assert (measured["n_cells"] >= 0).all()
    assert (measured["denom_assessed"] > 0).all()


def test_the_fine_level_never_exceeds_the_lineage_it_sits_in(layer):
    """t_nk_class sits inside the T_NK class and must never be added to it."""
    fine = layer["fine"]
    comp = layer["composition"]
    lineage = (
        comp[comp["state"].eq("measured")]
        .set_index(["tissue", "modality", "donor_id", "lineage_common"])["n_cells"]
    )
    grouped = fine.groupby(["tissue", "modality", "donor_id", "lineage_common", "fine_level"])[
        "n_cells"
    ].sum()
    for key, value in grouped.items():
        parent = lineage.get(key[:4])
        if parent is None or isinstance(parent, pd.Series):
            continue
        assert value <= parent, f"fine level {key} exceeds its lineage"


# ---------------------------------------------------------------------------
# The two tissues are different blocks.
# ---------------------------------------------------------------------------
def test_no_spatial_column_reaches_the_layer(layer):
    forbidden = {"signed_dist_islet_um", "distance_bin_um", "compartment",
                 "neighbourhood", "x_centroid", "y_centroid"}
    for table in ("composition", "presence", "fine"):
        assert not (forbidden & set(layer[table].columns))


def test_a_cross_tissue_quantity_is_refused(grid):
    with pytest.raises(grid.GridError):
        grid.refuse_cross_tissue("peri-islet expansion", ["pancreas", "pln"])
    grid.refuse_cross_tissue("peri-islet expansion", ["pancreas", "pancreas"])


def test_the_cross_tissue_join_carries_only_donor_and_status(layer):
    comp = layer["composition"]
    assert {"donor_id", "disease_status"} <= set(comp.columns)
    per_donor = comp.groupby("donor_id")["disease_status"].nunique(dropna=True)
    assert (per_donor <= 1).all(), "a donor carries two disease statuses across the grid"


# ---------------------------------------------------------------------------
# The crosswalk is read, never narrowed.
# ---------------------------------------------------------------------------
def test_every_crosswalk_lineage_appears_in_every_grid_cell(grid, layer):
    crosswalk = grid.load_crosswalk()
    wanted = set(crosswalk["common"]) - {"Unassigned"}
    comp = layer["composition"]
    for (tissue, modality, donor), rows in comp.groupby(["tissue", "modality", "donor_id"]):
        present = set(rows[rows["bucket_kind"].eq(grid.BK_LINEAGE)]["lineage_common"])
        assert wanted <= present, (
            f"{tissue}/{modality}/{donor} drops {sorted(wanted - present)}; a lineage is "
            "carried with a reason, never dropped for disagreeing"
        )


def test_resolution_columns_travel_with_every_row(layer):
    comp = layer["composition"]
    lineage_rows = comp[comp["in_crosswalk"]]
    for column in ("pheno_resolution", "xenium_resolution", "is_immune"):
        assert (lineage_rows[column].fillna("").str.len() > 0).all()


def test_a_disagreement_with_the_crosswalk_is_named_not_dropped(layer):
    comp = layer["composition"]
    conflicts = comp[comp["resolution_conflict"].fillna("").str.len() > 0]
    if not len(conflicts):
        pytest.skip("this build found no crosswalk disagreement")
    assert (conflicts["bucket_kind"] == "lineage").all()
    assert conflicts["resolution_conflict"].str.contains("crosswalk").all()


def test_map_targets_must_exist_in_the_crosswalk(grid, monkeypatch):
    crosswalk = grid.load_crosswalk()
    monkeypatch.setitem(grid.PANC_CLASS_TO_COMMON, "Beta", "Acinar")
    with pytest.raises(grid.GridError):
        grid.validate_maps(crosswalk)


def test_no_exocrine_subdivision_is_claimed(grid, layer):
    """The 60-plex carries no acinar-specific and no ductal-specific marker."""
    forbidden = re.compile(r"acinar|ductal|amylase|trypsin|CPA1|elastase|CK19|SOX9|HNF1B|CFTR|MUC1",
                           re.IGNORECASE)
    labels = set(layer["composition"]["lineage_common"]) | set(layer["fine"]["lineage_fine"])
    assert not [label for label in labels if forbidden.search(str(label))]


# ---------------------------------------------------------------------------
# Eta-squared never travels without its null.
# ---------------------------------------------------------------------------
def test_eta_squared_is_reported_beside_a_permutation_null(grid, layer):
    eta = layer["eta"]
    measured = eta[eta["state"].eq("measured")]
    assert len(measured)
    assert measured["null_median"].notna().all()
    assert measured["p_permutation"].notna().all()
    assert (measured["n_perm"] > 0).all()
    blank = eta[eta["state"].ne("measured")]
    assert blank["eta_squared"].isna().all()
    assert (blank["reason"].str.len() > 0).all()


def test_eta_with_null_is_bounded_and_honest(grid):
    rng = np.random.default_rng(0)
    values = rng.normal(size=12)
    groups = np.array(list("aaaabbbbcccc"), dtype=object)
    stat = grid.eta_with_null(values, groups, n_perm=200)
    assert 0.0 <= stat["eta_squared"] <= 1.0
    assert 0.0 < stat["p_permutation"] <= 1.0
    assert stat["n_donors"] == 12


def test_the_fusion_version_is_never_a_batch_variable(grid, layer):
    """Correction is only ever defensible on the finer acquisition variables."""
    assert "software" not in set(layer["eta"]["variable"])
    code = "\n".join(
        line for line in BUILDER.read_text().splitlines()
        if not line.lstrip().startswith("#")
    )
    for forbidden in ('"software"', "'software'", '"Fusion"', "'Fusion'"):
        assert forbidden not in code


# ---------------------------------------------------------------------------
# The layer is written where its reader looks, and belongs to one domain.
# ---------------------------------------------------------------------------
def test_the_layer_root_is_the_one_the_card_resolves(grid, monkeypatch, tmp_path):
    monkeypatch.setenv("ISLET_LAYERS_DIR", str(tmp_path))
    module = _load()
    assert module.OUT_DEFAULT == tmp_path / "grid_composition"
    monkeypatch.delenv("ISLET_LAYERS_DIR")
    module = _load()
    assert module.OUT_DEFAULT == Path("/data/islet_layers") / "grid_composition"
    card = CARD.read_text()
    assert 'Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers")' in card
    assert "grid_composition" in card


def test_the_layer_records_the_domain_it_belongs_to(layer):
    recorded = Path(layer["manifest"]["source_domain_root"])
    assert recorded.is_absolute() and recorded.exists()
    assert CARD.read_text().count("domain_data_root()") >= 1


def test_the_card_refuses_a_foreign_domain(layer):
    card = CARD.read_text()
    assert "grid_composition_belongs_to_a_different_domain" in card


# ---------------------------------------------------------------------------
# No cohort fact is written as a literal in the card.
# ---------------------------------------------------------------------------
def test_the_card_names_no_donor_no_cohort_size_and_no_method(layer):
    card = CARD.read_text()
    code = "\n".join(line for line in card.splitlines() if not line.lstrip().startswith("#"))
    for donor in sorted(set(layer["presence"]["donor_id"])):
        assert donor not in code, f"card names donor {donor} as a literal"
    for method in sorted({str(m) for m in layer["presence"]["method"] if str(m) != "nan"}):
        assert method not in code, f"card names the typing method {method!r} as a literal"
    for status in sorted({str(s) for s in layer["presence"]["disease_status"].dropna()}):
        assert f'"{status}"' not in code, f"card names the status level {status!r} as a literal"


def test_the_card_renders_a_named_note_for_every_state(grid):
    card = CARD.read_text()
    assert "ie_state_label" in card
    assert ".g1_note_ui" in card
    # A card that returned NULL for a missing cell would vanish without a trace.
    assert "ie_require_measured" in card


# ---------------------------------------------------------------------------
# Presence: every donor, every cell of the grid, with a reason.
# ---------------------------------------------------------------------------
def test_presence_is_complete_over_the_grid(layer):
    presence = layer["presence"]
    donors = set(presence["donor_id"])
    cells = set(map(tuple, presence[["tissue", "modality"]].drop_duplicates().to_numpy()))
    assert len(presence) == len(donors) * len(cells)
    assert not presence.duplicated(["tissue", "modality", "donor_id"]).any()


def test_a_measured_cell_has_cells_and_an_absent_one_has_none(layer):
    presence = layer["presence"]
    measured = presence[presence["state"].eq("measured")]
    assert (measured["n_cells_assessed"] > 0).all()
    blank = presence[presence["state"].eq("blank")]
    assert blank["n_cells_assessed"].isna().all()


def test_a_producer_reason_that_outlived_its_citation_is_flagged(layer):
    presence = layer["presence"]
    flagged = presence[presence["citation_conflict"].fillna("").str.len() > 0]
    if not len(flagged):
        pytest.skip("no stale citation in this build")
    assert (flagged["state"] != "measured").all()
    assert flagged["citation_conflict"].str.contains("no longer|now records").any()
