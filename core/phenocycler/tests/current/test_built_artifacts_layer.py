"""Acceptance for the five persisted artifacts, after running their real builders.

Uses ISLET_LAYERS_DIR (default /data/islet_layers), independently reconciles
distance denominators/numerators against physical typing columns, and checks the
current acquisition roster rather than assuming that every source covers it.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior.build_sections_layer import filter_eligible_donors, norm_donor

pytestmark = pytest.mark.data


@pytest.fixture(scope="module")
def artifacts():
    root = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers"))
    roster = REPO / "data/donor_acquisition_batch.csv"
    if not root.is_dir() or not roster.is_file():
        pytest.skip("Real layer artifacts and acquisition roster are not mounted")
    paths = {
        "availability": "availability/data.parquet",
        "sections": "sections/data.parquet",
        "islets": "islet_summary/islets.parquet",
        "composition": "islet_summary/composition.parquet",
        "donors": "islet_summary/donors.parquet",
        "summary_availability": "islet_summary/availability.parquet",
        "trajectory": "trajectory/data.parquet",
        "distance": "distance_profiles/data.parquet",
    }
    missing = [str(root / p) for p in paths.values() if not (root / p).is_file()]
    assert not missing, f"Run the five layer builders first: {missing}"
    with roster.open(newline="") as handle:
        cohort = set(filter_eligible_donors(
            norm_donor(row["donor_id"]) for row in csv.DictReader(handle)
            if row["eligible"].lower() == "true"
        ))
    tables = {name: pq.ParquetFile(root / path).read() for name, path in paths.items()}
    return cohort, tables


def test_persisted_layers_preserve_cohort_exclusions_frames_and_absence(artifacts):
    cohort, tables = artifacts
    for name, table in tables.items():
        assert table.num_rows > 0, name
        if "donor_id" in table.column_names:
            donors = set(table["donor_id"].to_pylist())
            assert None not in donors
            assert donors == set(filter_eligible_donors(donors)), name
            assert all(d == norm_donor(d) for d in donors), name
        if "frame" in table.column_names:
            assert table["frame"].null_count == 0, name
        frame = table.to_pandas()
        if "state" in frame:
            absent = frame.loc[frame.state != "measured"]
            assert absent.reason.str.strip().str.len().gt(0).all(), name
            for column in ("n", "n_cells", "n_immune", "immune_fraction", "within_donor_rank"):
                if column in absent:
                    assert absent[column].isna().all(), (name, column)
    ledger = tables["availability"].to_pandas()
    assert cohort <= set(ledger.donor_id)
    assert cohort <= set(tables["sections"]["donor_id"].to_pylist())
    assert not ledger.duplicated(["donor_id", "capability"]).any()
    assert ledger.groupby("donor_id").capability.nunique().eq(8).all()


def test_persisted_islets_and_trajectory_reconcile_with_their_sources(artifacts):
    _, tables = artifacts
    typing = REPO / "data/development/typing"
    expected = pq.ParquetFile(typing / "app/islet_geometry.parquet").read().to_pandas()
    expected["donor_id"] = expected.donor_id.map(norm_donor)
    expected = expected.loc[
        expected.region.eq("core") & expected.donor_id.isin(filter_eligible_donors(expected.donor_id))
    ].set_index(["donor_id", "islet_key"]).sort_index()
    actual = tables["islets"].to_pandas().set_index(["donor_id", "islet_key"]).sort_index()
    assert actual.index.equals(expected.index)
    np.testing.assert_array_equal(actual.n_cells, expected.n_cells)
    donors = tables["donors"].to_pandas()
    assert donors.groupby("donor_id").n_islets.sum().to_dict() == (
        expected.groupby("donor_id").size().to_dict()
    )

    source = pq.ParquetFile(typing / "trajectory/singlecell_islet_summary.parquet").read().to_pandas()
    expected_keys = {
        (norm_donor(str(row.case_id)), row.islet_key) for row in source.itertuples()
        if filter_eligible_donors([norm_donor(str(row.case_id))])
    }
    trajectory = tables["trajectory"].to_pandas()
    assert set(zip(trajectory.donor_id, trajectory.islet_key)) == expected_keys
    measured = trajectory.loc[trajectory.state == "measured"]
    for donor, group in measured.groupby("donor_id"):
        expected_rank = group.sc_pseudotime_core.rank(method="average")
        np.testing.assert_array_equal(group.within_donor_rank, expected_rank, err_msg=donor)
        assert group.n_islets_ordered.eq(len(group)).all()
    paired = trajectory.loc[trajectory.delta_state == "measured"]
    np.testing.assert_allclose(paired.sc_pseudotime_peri_minus_core,
                               paired.sc_pseudotime_peri - paired.sc_pseudotime_core)
    provenance = json.loads((typing / "trajectory/singlecell_trajectory.json").read_text())
    assert trajectory.n_ordered_cells.eq(provenance["n_ordered"]).all()
    assert trajectory.n_population_cells.eq(provenance["n_population"]).all()


def test_persisted_distance_counts_reconcile_with_every_typing_partition(artifacts):
    _, tables = artifacts
    table = tables["distance"]
    profile = table.to_pandas()
    assert profile.frame.eq("signed_distance_um").all()
    audit = json.loads(table.schema.metadata[b"audit"])
    cells_root = REPO / "data/development/typing/cells"
    eligible_counts, immune_counts, source_counts = {}, {}, {}
    for path in sorted(cells_root.glob("donor_id=*/*.parquet")):
        donor = norm_donor(path.parent.name.removeprefix("donor_id="))
        if not filter_eligible_donors([donor]):
            continue
        eligible_counts.setdefault(donor, 0)
        immune_counts.setdefault(donor, 0)
        source_counts.setdefault(donor, 0)
        parquet = pq.ParquetFile(path)
        assert "donor_id" not in parquet.schema_arrow.names
        source_counts[donor] += parquet.metadata.num_rows
        for batch in parquet.iter_batches(columns=[
            "qc_analysis_eligible", "cell_class", "excluded_reason",
        ]):
            eligible = pc.and_(pc.fill_null(batch["qc_analysis_eligible"], False), pc.and_(
                pc.fill_null(pc.not_equal(batch["cell_class"], "Excluded"), True),
                pc.equal(pc.fill_null(batch["excluded_reason"], ""), ""),
            ))
            immune = pc.and_(eligible, pc.fill_null(pc.equal(batch["cell_class"], "Immune"), False))
            eligible_counts[donor] += pc.sum(eligible).as_py() or 0
            immune_counts[donor] += pc.sum(immune).as_py() or 0
    assert set(profile.donor_id) == set(eligible_counts)
    for donor, group in profile.groupby("donor_id"):
        assert audit[donor]["source_cells"] == source_counts[donor]
        assert audit[donor]["eligible_cells"] == eligible_counts[donor]
        if group.state.eq("measured").all():
            totals = group.loc[group.resolved_by == "all"]
            assert totals.n_cells.sum() == eligible_counts[donor]
            assert totals.n_immune.sum() == immune_counts[donor]
            np.testing.assert_allclose(group.immune_fraction, group.n_immune / group.n_cells)
            assert audit[donor]["missing_membership"] == audit[donor]["invalid_membership"] == 0
        else:
            assert len(group) == 1
            assert group.state.eq("unavailable").all()
            assert group[["n_cells", "n_immune", "immune_fraction"]].isna().all().all()
