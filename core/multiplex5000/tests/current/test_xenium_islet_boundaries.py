import json

import numpy as np
import pandas as pd
import pytest
import shapely

from multiplex5000.integration.xenium_islets import (
    BoundaryParams, _digest, cluster_boundaries, read_section, run_manifest,
)


def cells():
    xy = np.array([[0, 0], [0, 1], [1, 0], [1, 1],
                   [10, 0], [10, 1], [11, 0], [11, 1], [50, 50]], float)
    return pd.DataFrame(dict(cell_id=[f"c{i}" for i in range(len(xy))],
        x_um=xy[:, 0], y_um=xy[:, 1], nucleus_count=1, depth=50, seed_evidence=True))


def test_boundaries_preserve_seed_keys_area_and_order_invariance():
    frame = cells()
    a, members, qc = cluster_boundaries(frame, acquisition_id="0000001:panc", params=BoundaryParams(2, 3))
    b, reordered, _ = cluster_boundaries(frame.sample(frac=1, random_state=7),
        acquisition_id="0000001:panc", params=BoundaryParams(2, 3))
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(members, reordered)
    assert a.area_um2.tolist() == [1., 1.]
    assert qc["n_noise_seeds"] == 1 and qc["n_estimated_boundaries"] == 2
    for row in a.itertuples():
        selected = members[members.structure_id.eq(row.structure_id)]
        assert shapely.covers(shapely.from_wkb(row.geometry_wkb),
                              shapely.points(selected.x_um, selected.y_um)).all()
    other, _, _ = cluster_boundaries(frame, acquisition_id="0000002:panc", params=BoundaryParams(2, 3))
    assert set(a.structure_id).isdisjoint(other.structure_id)


def test_no_nucleus_or_depth_cannot_seed_even_when_requested():
    frame = cells()
    frame.loc[:3, "nucleus_count"] = 0
    frame.loc[4:7, "depth"] = 0
    objects, members, qc = cluster_boundaries(frame, acquisition_id="s", params=BoundaryParams(2, 3))
    assert objects.empty and len(members) == 1
    assert qc["n_seed_evidence_excluded_by_nucleus_or_depth"] == 8


def test_degenerate_objects_and_empty_seeds_remain_explicit():
    frame = cells().iloc[:4].copy()
    frame.y_um = 0
    objects, members, _ = cluster_boundaries(frame, acquisition_id="s", params=BoundaryParams(2, 3))
    assert objects.boundary_status.tolist() == ["unavailable_degenerate_cluster"]
    assert objects.geometry_wkb.isna().all()
    frame.seed_evidence = False
    objects, members, qc = cluster_boundaries(frame, acquisition_id="s", params=BoundaryParams(2, 3))
    assert objects.empty and members.empty and qc["n_seeds"] == 0


def test_bad_keys_coordinates_or_missing_seed_evidence_are_rejected():
    frame = cells()
    with pytest.raises(ValueError, match="unique"):
        cluster_boundaries(pd.concat([frame, frame]), acquisition_id="s")
    frame.loc[0, "x_um"] = np.nan
    with pytest.raises(ValueError, match="Finite"):
        cluster_boundaries(frame, acquisition_id="s")
    frame = cells(); frame.seed_evidence = "False"
    with pytest.raises(ValueError, match="Boolean"):
        cluster_boundaries(frame, acquisition_id="s")
    with pytest.raises(ValueError):
        BoundaryParams(float("nan"), 10)


def fixture_manifest(tmp_path):
    frame = cells()
    raw = frame.rename(columns={"x_um": "x_centroid", "y_um": "y_centroid"})
    calls = frame.assign(source_serial="0000001", donor_id="test", roi="panc", lineage="Endocrine",
        lineage_analysis_eligible=True, n_naming_detected_for_cluster_call=1, nucleus_measured=True)
    # A zero-evidence call is not a seed even though the broad label says endocrine.
    calls.loc[8, "n_naming_detected_for_cluster_call"] = 0
    raw.to_parquet(tmp_path / "cells.parquet", index=False)
    calls.to_parquet(tmp_path / "calls.parquet", index=False)
    spec = dict(donor_id="test", source_serial="0000001", roi="panc", acquisition_id="0000001:panc",
        coordinate_units="micrometres", source_release_id="test-candidate", source_release_status="for_review",
        cells_path=str(tmp_path / "cells.parquet"), cells_sha256=_digest(tmp_path / "cells.parquet"),
        calls_path=str(tmp_path / "calls.parquet"), calls_sha256=_digest(tmp_path / "calls.parquet"))
    return dict(schema="xenium_islet_boundary_inputs_v1", sections=[spec],
                profiles={"primary": dict(eps_um=2., min_samples=3)})


def test_pinned_adapter_and_export_bundle(tmp_path):
    manifest = fixture_manifest(tmp_path)
    path = tmp_path / "manifest.json"; path.write_text(json.dumps(manifest))
    out = tmp_path / "out"
    result = run_manifest(path, out)
    assert result.n_seeds.tolist() == [8]
    geo = json.loads((out / "primary/0000001/boundaries_um.geojson").read_text())
    assert len(geo["features"]) == 2 and geo["coordinate_units"] == "micrometres"
    assert not geo["individual_islet_accuracy_validated"]
    with pytest.raises(FileExistsError):
        run_manifest(path, out)
    bad = dict(manifest["sections"][0], source_serial="0000002")
    with pytest.raises(ValueError, match="source_serial"):
        read_section(bad)
    bad = dict(manifest["sections"][0], calls_sha256="bad")
    with pytest.raises(ValueError, match="SHA-256"):
        read_section(bad)


def test_declared_subset_keeps_missing_call_evidence_explicit(tmp_path):
    spec = fixture_manifest(tmp_path)["sections"][0]
    path = tmp_path / "calls.parquet"
    pd.read_parquet(path).iloc[:8].to_parquet(path, index=False)
    spec["calls_sha256"] = _digest(path)
    with pytest.raises(ValueError, match="subset"):
        read_section(spec)
    spec["calls_scope"] = "declared_subset_of_vendor_cells"
    frame = read_section(spec)
    assert len(frame) == 9 and frame.attrs["input_accounting"]["n_vendor_cells_without_calls"] == 1
    assert frame.loc[frame.cell_id.eq("c8"), "depth"].isna().all()
    assert not frame.loc[frame.cell_id.eq("c8"), "seed_evidence"].any()
