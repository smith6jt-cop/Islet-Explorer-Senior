from pathlib import Path
import sys
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[2] / "scripts" / "typing"
sys.path.insert(0, str(SCRIPT_DIR))
import build_nd_referenced_trajectories as trajectories  # noqa: E402


def test_unavailable_stromal_identity_does_not_block_measured_group_trajectory():
    gate = trajectories.feature_group_gate(
        ["endocrine", "vascular", "lymphatic", "neural", "immune"])

    assert gate["missing_groups"] == ["stromal"]
    assert gate["nonblocking_missing_groups"] == ["stromal"]
    assert gate["blocking_missing_groups"] == []


def test_other_missing_biological_group_still_blocks_primary_trajectory():
    gate = trajectories.feature_group_gate(
        ["endocrine", "lymphatic", "neural", "immune"])

    assert gate["missing_groups"] == ["stromal", "vascular"]
    assert gate["nonblocking_missing_groups"] == ["stromal"]
    assert gate["blocking_missing_groups"] == ["vascular"]


def test_complete_feature_groups_remain_complete_mode():
    gate = trajectories.feature_group_gate(list(trajectories.BIOLOGICAL_GROUPS))

    assert gate["missing_groups"] == []
    assert gate["blocking_missing_groups"] == []
    assert gate["nonblocking_missing_groups"] == []


def test_transform_frame_materializes_writable_numeric_values():
    source = pd.DataFrame({"feature__vascular__cell_bodies_per_endocrine": [0.0, 1.0]})
    transformed = trajectories.transform_frame(
        source, ["feature__vascular__cell_bodies_per_endocrine"])

    np.testing.assert_allclose(transformed[:, 0], [0.0, np.log(2.0)])
    np.testing.assert_allclose(source.iloc[:, 0], [0.0, 1.0])


def test_transform_frame_clips_only_fraction_roundoff():
    feature = "feature__vascular__outline_coverage_fraction"
    source = pd.DataFrame({feature: [0.0, 1.0 + 1.8e-15]})

    transformed = trajectories.transform_frame(source, [feature])

    np.testing.assert_allclose(transformed[:, 0], [0.0, 1.0])
    audit = trajectories.fraction_roundoff_audit(source, [feature])
    assert audit[feature]["count_clipped"] == 1


def test_transform_frame_rejects_materially_invalid_fraction():
    feature = "feature__vascular__outline_coverage_fraction"
    source = pd.DataFrame({feature: [1.0001]})

    try:
        trajectories.transform_frame(source, [feature])
    except ValueError as exc:
        assert "outside [0,1]" in str(exc)
    else:
        raise AssertionError("materially invalid fraction was accepted")


def test_rank_selection_includes_one_coordinate_independently_of_diffusion_modes():
    table = pd.DataFrame({
        "rank": [1, 2, 3],
        "mean_donor_balanced_rmse": [0.5, 0.8, 0.75],
        "se_across_heldout_donors": [0.1, 0.2, 0.1],
    })

    assert trajectories.select_rank_one_se(table) == 1


def test_diffusion_modes_are_not_limited_by_input_coordinate_count():
    scores = np.linspace(-1.0, 1.0, 120).reshape(-1, 1)
    ids = np.asarray([f"object-{i:03d}" for i in range(len(scores))])
    result = trajectories.graph_dpt(scores, ids, np.asarray([0.0]),
                                    n_neighbors=20, dimensions=20)

    assert np.isfinite(result["pseudotime"]).all()
    assert set(result["status"]) == {"ordered"}
    assert set(result["component_size"]) == {120}


def test_fast_pairwise_tau_matches_brute_force_with_ties_and_disconnected_components():
    rng = np.random.default_rng(20261003)
    ids = [f"object-{i}" for i in range(80)]
    left = {"ids": ids,
            "pseudo": rng.integers(0, 9, size=len(ids)).astype(float),
            "component": rng.integers(0, 5, size=len(ids))}
    right = {"ids": ids,
             "pseudo": rng.integers(0, 11, size=len(ids)).astype(float),
             "component": rng.integers(0, 6, size=len(ids))}
    fast_tau, fast_coverage = trajectories.pairwise_within_component_tau(left, right, ids)

    concordant = discordant = 0
    for i, left_id in enumerate(ids):
        for right_id in ids[i + 1:]:
            ia, ja = ids.index(left_id), ids.index(right_id)
            if left["component"][ia] != left["component"][ja] or \
                    right["component"][ia] != right["component"][ja]:
                continue
            dx = np.sign(left["pseudo"][ia] - left["pseudo"][ja])
            dy = np.sign(right["pseudo"][ia] - right["pseudo"][ja])
            if not dx or not dy:
                continue
            if dx == dy:
                concordant += 1
            else:
                discordant += 1
    expected = ((concordant - discordant) / (concordant + discordant)
                if concordant + discordant else np.nan)
    assert fast_coverage == 1.0
    assert np.isclose(fast_tau, expected)
