from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import tifffile

from phenocycler.artifacts import ContractError
from scripts.review.boundary_contact_audit import (
    DEFAULT_DEVELOPMENT_ROOT,
    DEFAULT_LABEL,
    SCRATCH_AVAILABLE,
    SCRATCH_INVALID,
    MarkerPair,
    TiledPageReader,
    _candidate_frame,
    audit_scratch_mapping,
    contact_metrics,
    display_marker,
    extract_tiled_window,
    load_live_images,
    pairs_from_roster,
    parse_pairs,
    prepare_development_workspace,
    select_crop_manifest,
    selected_contact_edges,
    selection_source_kind,
    selection_source_measurement,
)
from scripts.review.qupath_fundamental_inputs import (
    METHOD_VERSION as PROJECTION_METHOD_VERSION,
)
from scripts.review.qupath_fundamental_inputs import (
    SOURCE_MEASUREMENTS,
)
from scripts.review.qupath_fundamental_inputs import (
    WORKSPACE_MARKER as PROJECTION_WORKSPACE_MARKER,
)


def _declare_projection(root, donor: str = "D", row_count: int = 2) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / PROJECTION_WORKSPACE_MARKER).write_text(
        PROJECTION_METHOD_VERSION + "\n", encoding="utf-8"
    )
    (root / "projection_summary.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "method_version": PROJECTION_METHOD_VERSION,
                "source_kind": "upstream QuPath full-cell measurement export",
                "source_measurements": SOURCE_MEASUREMENTS,
                "restore_outputs_read": False,
                "donors": {donor: row_count},
                "donor_images": {donor: "image"},
            }
        ),
        encoding="utf-8",
    )


def _minimal_scratch(*, include_extra: bool = True) -> dict[str, np.ndarray]:
    object_ids = ["cell-a", "cell-b", "cell-c"]
    if include_extra:
        object_ids.append("scratch-only")
    dense = np.zeros((len(object_ids), len(object_ids)), dtype=np.float32)
    dense[0, 1] = dense[1, 0] = 3
    dense[1, 2] = dense[2, 1] = 2
    if include_extra:
        dense[0, 3] = dense[3, 0] = 1
    contact = sp.csr_matrix(dense)
    return {
        "oids": np.asarray(object_ids, dtype=object),
        "m_data": contact.data,
        "m_indices": contact.indices,
        "m_indptr": contact.indptr,
        "m_shape": np.asarray(contact.shape),
    }


def test_default_pair_roster_and_explicit_iba1_display_mapping():
    assert parse_pairs(None) == (
        MarkerPair("CD3e", "Ker8_18"),
        MarkerPair("Iba1_norm", "Ker8_18"),
    )
    assert parse_pairs(["CD4:Iba1_norm"]) == (MarkerPair("CD4", "Iba1_norm"),)
    assert display_marker("Iba1_norm") == "Iba1"
    assert display_marker("CD3e") == "CD3e"
    assert selection_source_measurement("Iba1_norm") == "Cell: Iba1: Mean (norm)"
    assert selection_source_kind("Iba1_norm") == ("qupath_acquisition_normalized_cell_mean")
    assert selection_source_measurement("Ker8_18") == "Cell: Ker8-18: Mean"
    assert DEFAULT_LABEL == "boundary_contact_diagnostics"
    assert DEFAULT_DEVELOPMENT_ROOT.as_posix().endswith("/data/development")


def test_boundary_defaults_follow_a_reduced_shared_roster(tmp_path):
    from scripts.review.fundamental_pair_survey import ROSTER_PATH

    payload = json.loads(ROSTER_PATH.read_text(encoding="utf-8"))
    payload["pairs"] = payload["pairs"][:1]
    roster = tmp_path / "roster.json"
    roster.write_text(json.dumps(payload), encoding="utf-8")

    expected = (MarkerPair("CD3e", "Ker8_18"),)
    assert pairs_from_roster(roster) == expected
    assert parse_pairs(None, roster_path=roster) == expected


def test_live_image_registration_ignores_content_ids_and_file_hashes(tmp_path):
    path = tmp_path / "qupath.json"
    path.write_text(
        json.dumps(
            {
                "content_id": "intentionally-wrong-and-ignored",
                "images": [
                    {
                        "donor_id": "D",
                        "image_id": "image",
                        "qptiff": {"path": "/current/image.qptiff", "sha256": "wrong"},
                        "cell_geojson": {"path": "/current/cells.geojson", "sha256": "wrong"},
                        "nucleus_geojson": {
                            "path": "/current/nuclei.geojson",
                            "sha256": "wrong",
                        },
                        "combined_geometry_file": False,
                        "pixel_size_um_x": 0.5,
                        "pixel_size_um_y": 0.5,
                        "channel_map": [
                            {"channel_index": 0, "channel_name": "DAPI", "marker": "DAPI"}
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    images = load_live_images(path)

    assert len(images) == 1
    assert images[0].qptiff.path == "/current/image.qptiff"


def test_minimal_contact_only_scratch_is_available_and_emits_metrics():
    result = audit_scratch_mapping(
        _minimal_scratch(),
        ["cell-a", "cell-b", "cell-c"],
    )

    assert result.state == SCRATCH_AVAILABLE
    assert result.validation is not None
    assert result.validation.scratch_only_n == 1
    assert result.validation.current_only_n == 0
    metrics = contact_metrics(result.validation).set_index("object_id")
    assert metrics.loc["cell-a", "contact_neighbor_n"] == 2
    assert metrics.loc["cell-c", "contact_neighbor_n"] == 1


def test_partial_current_overlap_is_reported_without_suppressing_contacts():
    result = audit_scratch_mapping(
        _minimal_scratch(include_extra=False),
        ["cell-a", "cell-b", "current-only"],
    )

    assert result.state == SCRATCH_AVAILABLE
    validation = result.validation
    assert validation is not None
    assert validation.current_only_n == 1
    assert validation.current_only_examples == ("current-only",)
    assert validation.scratch_only_n == 1
    assert validation.scratch_only_examples == ("cell-c",)
    assert contact_metrics(validation)["object_id"].tolist() == ["cell-a", "cell-b"]

    edges = selected_contact_edges(validation, ["cell-a", "current-only"])
    assert set(edges["object_id"]) == {"cell-a"}
    assert "neighbor_is_current" in edges


def test_noncontact_cell_metadata_are_not_scratch_acceptance_inputs():
    values = _minimal_scratch(include_extra=False)
    live = pd.DataFrame({"object_id": ["cell-a", "cell-b", "cell-c"], "cell_area_px": [1, 2, 3]})
    first = audit_scratch_mapping(values, live["object_id"])
    live["cell_area_px"] = [1000, 2000, 3000]
    second = audit_scratch_mapping(values, live["object_id"])

    assert first.state == second.state == SCRATCH_AVAILABLE
    pd.testing.assert_frame_equal(
        contact_metrics(first.validation),
        contact_metrics(second.validation),
    )


def test_old_compartment_arrays_are_optional_and_ignored():
    values = _minimal_scratch(include_extra=False)
    values["data"] = np.array([[np.nan]])
    values["sizes"] = np.array([-99])

    result = audit_scratch_mapping(values, ["cell-a", "cell-b", "cell-c"])

    assert result.state == SCRATCH_AVAILABLE


def test_iba1_norm_selection_uses_exact_upstream_qupath_projection(tmp_path):
    run_root = tmp_path / "run"
    projection_root = tmp_path / "qupath-inputs"
    donor = "D"
    _declare_projection(projection_root, donor)
    (run_root / "geometry_qc" / f"donor_id={donor}").mkdir(parents=True)
    (projection_root / f"donor_id={donor}").mkdir(parents=True)
    pd.DataFrame(
        {
            "donor_id": [donor, donor],
            "object_id": ["a", "b"],
            "image": ["image", "image"],
            "X_centroid": [1.0, 2.0],
            "Y_centroid": [3.0, 4.0],
            "cell_area": [10.0, 11.0],
            "cell_solidity": [0.8, 0.9],
            "DAPI": [5.0, 6.0],
            "CD3e": [7.0, 8.0],
            # This alias is created only from exact Cell: Iba1: Mean (norm)
            # by qupath_fundamental_inputs; there is no raw-Iba1 fallback.
            "Iba1_norm": [-0.25, 0.75],
        }
    ).to_parquet(projection_root / f"donor_id={donor}" / "data.parquet")
    pd.DataFrame(
        {
            "object_id": ["a", "b"],
            "analysis_eligible": [True, True],
            "raster_area_ratio": [1.0, 1.0],
            "nucleus_cell_area_ratio": [0.5, 0.5],
        }
    ).to_parquet(run_root / "geometry_qc" / f"donor_id={donor}" / "data.parquet")
    observed, coverage = _candidate_frame(
        run_root,
        projection_root,
        donor,
        (MarkerPair("CD3e", "Iba1_norm"),),
    )

    assert observed["Iba1_norm"].tolist() == [-0.25, 0.75]
    assert coverage.projection_cell_n == coverage.current_cell_n == 2
    assert coverage.matched_cell_n == 2
    assert coverage.projection_only_n == coverage.current_only_n == 0


def test_candidate_frame_inner_joins_partial_projection_and_geometry_coverage(tmp_path):
    run_root = tmp_path / "run"
    projection_root = tmp_path / "qupath-inputs"
    donor = "D"
    _declare_projection(projection_root, donor, row_count=4)
    (run_root / "geometry_qc" / f"donor_id={donor}").mkdir(parents=True)
    (projection_root / f"donor_id={donor}").mkdir(parents=True)
    pd.DataFrame(
        {
            "donor_id": [donor] * 4,
            "object_id": ["matched-b", "projection-z", "matched-a", "projection-a"],
            "image": ["image"] * 4,
            "X_centroid": [1.0, 2.0, 3.0, 4.0],
            "Y_centroid": [5.0, 6.0, 7.0, 8.0],
            "cell_area": [10.0, 11.0, 12.0, 13.0],
            "cell_solidity": [0.8, 0.9, 0.85, 0.95],
            "DAPI": [5.0, 6.0, 7.0, 8.0],
            "CD3e": [0.1, 0.2, 0.3, 0.4],
            "Iba1_norm": [-0.25, 1.25, 0.75, -1.0],
        }
    ).to_parquet(projection_root / f"donor_id={donor}" / "data.parquet")
    pd.DataFrame(
        {
            "object_id": ["current-z", "matched-a", "current-a", "matched-b"],
            "analysis_eligible": [True, True, False, True],
            "raster_area_ratio": [1.0, 1.0, 1.0, 1.0],
            "nucleus_cell_area_ratio": [0.5, 0.5, 0.5, 0.5],
        }
    ).to_parquet(run_root / "geometry_qc" / f"donor_id={donor}" / "data.parquet")

    observed, coverage = _candidate_frame(
        run_root,
        projection_root,
        donor,
        (MarkerPair("CD3e", "Iba1_norm"),),
    )

    observed_by_id = observed.assign(object_id=observed["object_id"].astype(str)).set_index(
        "object_id"
    )
    assert set(observed_by_id.index) == {"matched-a", "matched-b"}
    assert observed_by_id["Iba1_norm"].to_dict() == {
        "matched-a": 0.75,
        "matched-b": -0.25,
    }
    assert coverage.projection_cell_n == 4
    assert coverage.current_cell_n == 4
    assert coverage.matched_cell_n == 2
    assert coverage.projection_only_n == 2
    assert coverage.current_only_n == 2
    assert coverage.projection_only_examples == ("projection-a", "projection-z")
    assert coverage.current_only_examples == ("current-a", "current-z")


@pytest.mark.parametrize("source", ["projection", "geometry"])
@pytest.mark.parametrize("fault", ["null", "duplicate"])
def test_candidate_frame_rejects_invalid_object_ids(tmp_path, source: str, fault: str):
    run_root = tmp_path / "run"
    projection_root = tmp_path / "qupath-inputs"
    donor = "D"
    _declare_projection(projection_root, donor)
    (run_root / "geometry_qc" / f"donor_id={donor}").mkdir(parents=True)
    (projection_root / f"donor_id={donor}").mkdir(parents=True)
    projection_ids: list[str | None] = ["a", "b"]
    geometry_ids: list[str | None] = ["a", "b"]
    target = projection_ids if source == "projection" else geometry_ids
    target[1] = None if fault == "null" else "a"
    pd.DataFrame(
        {
            "donor_id": [donor, donor],
            "object_id": projection_ids,
            "image": ["image", "image"],
            "X_centroid": [1.0, 2.0],
            "Y_centroid": [3.0, 4.0],
            "cell_area": [10.0, 11.0],
            "cell_solidity": [0.8, 0.9],
            "DAPI": [5.0, 6.0],
            "CD3e": [7.0, 8.0],
            "Iba1_norm": [-0.25, 0.75],
        }
    ).to_parquet(projection_root / f"donor_id={donor}" / "data.parquet")
    pd.DataFrame(
        {
            "object_id": geometry_ids,
            "analysis_eligible": [True, True],
            "raster_area_ratio": [1.0, 1.0],
            "nucleus_cell_area_ratio": [0.5, 0.5],
        }
    ).to_parquet(run_root / "geometry_qc" / f"donor_id={donor}" / "data.parquet")

    with pytest.raises(ContractError, match="object IDs must be"):
        _candidate_frame(
            run_root,
            projection_root,
            donor,
            (MarkerPair("CD3e", "Iba1_norm"),),
        )


@pytest.mark.parametrize("fault", ["missing", "duplicate_ids", "diagonal", "asymmetric"])
def test_invalid_contact_structure_is_rejected(fault: str):
    values = _minimal_scratch(include_extra=False)
    if fault == "missing":
        del values["m_indptr"]
    elif fault == "duplicate_ids":
        values["oids"][1] = values["oids"][0]
    else:
        dense = np.zeros((3, 3), dtype=np.float32)
        if fault == "diagonal":
            dense[0, 0] = 1
        else:
            dense[0, 1] = 1
        contact = sp.csr_matrix(dense)
        values.update(
            m_data=contact.data,
            m_indices=contact.indices,
            m_indptr=contact.indptr,
            m_shape=np.asarray(contact.shape),
        )

    result = audit_scratch_mapping(values, ["cell-a", "cell-b", "cell-c"])

    assert result.state == SCRATCH_INVALID
    assert result.validation is None


def test_same_label_refreshes_only_an_owned_development_workspace(tmp_path):
    pilots = tmp_path / "data" / "pilots"
    first = prepare_development_workspace(pilots, "boundary-debug")
    (first / "old-output.txt").write_text("replace me", encoding="utf-8")

    second = prepare_development_workspace(pilots, "boundary-debug")

    assert second == first
    assert not (second / "old-output.txt").exists()
    assert (second / ".boundary_contact_development").is_file()

    unrelated = pilots / "other-tool"
    unrelated.mkdir()
    (unrelated / "result.txt").write_text("preserve me", encoding="utf-8")
    with pytest.raises(ContractError, match="not owned"):
        prepare_development_workspace(pilots, "other-tool")
    assert (unrelated / "result.txt").read_text(encoding="utf-8") == "preserve me"


@pytest.mark.parametrize("label", ["", ".", "..", "nested/name"])
def test_development_workspace_rejects_non_component_labels(tmp_path, label: str):
    with pytest.raises(ContractError, match="one nonempty path component"):
        prepare_development_workspace(tmp_path / "pilots", label)


def test_mixed_dtype_crop_selection_is_rank_stable_and_outcome_blind():
    cells = pd.DataFrame(
        {
            "object_id": [f"cell-{index:02d}" for index in range(18)],
            "X_centroid": np.linspace(10, 180, 18),
            "Y_centroid": np.linspace(20, 190, 18),
            "DAPI": np.arange(18) + 1,
            "Iba1_norm": np.array(
                [-0.4, 1.9, -0.2, 1.7, 0.0, 1.5, 0.2, 1.3, 0.4] * 2,
                dtype=np.float64,
            ),
            "Ker8_18": np.array(
                [900, 100, 800, 200, 700, 300, 600, 400, 500] * 2,
                dtype=np.uint16,
            ),
            "analysis_eligible": True,
            "cell_area": np.arange(18) + 20,
            "cell_solidity": 0.9,
            "raster_area_ratio": 1.0,
            "nucleus_cell_area_ratio": 0.5,
            "disease_status": ["never-read"] * 18,
            "phenotype": ["never-read"] * 18,
            "called_positive": np.arange(18) % 2 == 0,
        }
    )
    kwargs = {
        "donor_id": "donor",
        "image_id": "image",
        "pairs": (MarkerPair("Iba1_norm", "Ker8_18"),),
        "image_shape": (256, 256),
        "pixel_size_um_x": 1.0,
        "pixel_size_um_y": 1.0,
        "crops_per_stratum": 3,
        "crop_size_px": 32,
    }

    first = select_crop_manifest(cells, **kwargs)
    second = select_crop_manifest(cells.sample(frac=1, random_state=11), **kwargs)
    rescaled = cells.copy()
    rescaled["Ker8_18"] = rescaled["Ker8_18"].astype(np.float64) * 1e-3 + 7.0
    third = select_crop_manifest(rescaled, **kwargs)

    pd.testing.assert_frame_equal(first, second)
    invariant = [column for column in first.columns if column != "selection_reference_cell_value"]
    pd.testing.assert_frame_equal(first[invariant], third[invariant])
    np.testing.assert_allclose(
        third["selection_reference_cell_value"],
        first["selection_reference_cell_value"] * 1e-3 + 7.0,
    )
    assert len(first) == 9
    assert set(first["stratum"]) == {
        "dual_bright",
        "target_bright",
        "reference_bright",
    }
    assert set(first["target_display_marker"]) == {"Iba1"}
    assert set(first["reference_display_marker"]) == {"Ker8_18"}
    assert set(first["target_selection_source_measurement"]) == {"Cell: Iba1: Mean (norm)"}
    assert set(first["reference_selection_source_measurement"]) == {"Cell: Ker8-18: Mean"}
    assert set(first["target_display_channel"]) == {"Iba1"}
    assert set(first["reference_display_channel"]) == {"Ker8_18"}
    assert set(first["pair_finite_n"]) == {18}
    assert set(first["stratum_requested_n"]) == {3}
    assert set(first["stratum_selected_n"]) == {3}
    assert not {"disease_status", "phenotype", "called_positive"}.intersection(first.columns)


def test_crop_selection_fails_on_bad_or_insufficient_pair_values():
    base = pd.DataFrame(
        {
            "object_id": [f"cell-{index}" for index in range(4)],
            "X_centroid": np.arange(4, dtype=float),
            "Y_centroid": np.arange(4, dtype=float),
            "DAPI": np.ones(4),
            "Iba1_norm": np.linspace(-0.2, 0.5, 4),
            "Ker8_18": np.arange(4, dtype=float),
            "analysis_eligible": True,
        }
    )
    kwargs = {
        "donor_id": "D",
        "image_id": "image",
        "pairs": (MarkerPair("Iba1_norm", "Ker8_18"),),
        "image_shape": (32, 32),
        "pixel_size_um_x": 1.0,
        "pixel_size_um_y": 1.0,
        "crops_per_stratum": 1,
        "crop_size_px": 8,
    }

    bad = base.copy()
    bad["Iba1_norm"] = "not-a-number"
    with pytest.raises(ContractError, match="not numeric"):
        select_crop_manifest(bad, **kwargs)

    no_finite = base.copy()
    no_finite["Iba1_norm"] = np.nan
    with pytest.raises(ContractError, match="no finite crop candidates"):
        select_crop_manifest(no_finite, **kwargs)

    with pytest.raises(ContractError, match="requested 2 crops"):
        select_crop_manifest(base.iloc[:5], **{**kwargs, "crops_per_stratum": 2})


def test_tiled_window_uses_pixel_values_api_not_whole_page_read():
    class ValuesOnlyReader:
        def __init__(self):
            self.calls = 0

        def values(self, rows, columns):
            self.calls += 1
            return np.asarray(rows) * 100 + np.asarray(columns)

        def asarray(self):  # pragma: no cover - must never be reached
            raise AssertionError("whole-page reads are forbidden")

    reader = ValuesOnlyReader()
    observed = extract_tiled_window(reader, x0=7, y0=5, width=3, height=2)

    np.testing.assert_array_equal(observed, [[507, 508, 509], [607, 608, 609]])
    assert reader.calls == 1


def test_tiled_page_reader_returns_exact_requested_pixels(tmp_path):
    values = np.arange(32 * 32, dtype=np.uint16).reshape(32, 32)
    path = tmp_path / "image.tiff"
    tifffile.imwrite(path, values, tile=(16, 16))

    with tifffile.TiffFile(path) as image:
        reader = TiledPageReader(image.filehandle, image.pages[0], cache_tiles=2)
        rows = np.array([0, 15, 16, 31, 3, 20])
        columns = np.array([0, 15, 16, 31, 20, 3])
        observed = reader.values(rows, columns)

    np.testing.assert_array_equal(observed, values[rows, columns])
