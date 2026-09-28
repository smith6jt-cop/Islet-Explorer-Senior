"""Synthetic tests for centroid-based neighbour features: k-NN, max-intensity, fraction, and
distance to the nearest cell exceeding a threshold (the bleed-through distance)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.neighbour_features import (
    distance_to_nearest_exceeding,
    max_neighbour_intensity,
    neighbour_fraction,
    query_neighbours,
)

# A line of four evenly-spaced cells plus one isolated cell far away.
LINE_COORDS = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [100.0, 0.0]])


def test_k_nearest_within_radius_matches_hand_computed_neighbours():
    query = query_neighbours(LINE_COORDS, k=2, radius_um=1.5)
    assert query.indices.shape == (5, 2)

    def neighbours_of(row):
        valid = query.indices[row] >= 0
        return sorted(query.indices[row][valid].tolist())

    assert neighbours_of(0) == [1]
    assert neighbours_of(1) == [0, 2]
    assert neighbours_of(2) == [1, 3]
    assert neighbours_of(3) == [2]
    assert neighbours_of(4) == []
    assert np.isnan(query.distances[4]).all()


def test_self_is_never_returned_even_with_a_coincident_duplicate():
    coords = np.array([[0.0, 0.0], [0.0, 0.0], [5.0, 0.0]])
    query = query_neighbours(coords, k=1, radius_um=1.0)
    assert query.indices[0, 0] == 1
    assert query.distances[0, 0] == pytest.approx(0.0)
    assert query.indices[1, 0] == 0
    assert query.indices[2, 0] == -1  # cell 2 is 5um from the pair, outside radius


def test_max_neighbour_intensity_picks_the_highest_neighbour_value():
    query = query_neighbours(LINE_COORDS, k=2, radius_um=1.5)
    values = np.array([10.0, 20.0, 30.0, 5.0, 999.0])
    out = max_neighbour_intensity(query, values)
    assert out[0] == pytest.approx(20.0)  # cell 0's only neighbour is cell 1
    assert out[1] == pytest.approx(30.0)  # neighbours 0 (10) and 2 (30)
    assert out[2] == pytest.approx(20.0)  # neighbours 1 (20) and 3 (5)
    assert out[3] == pytest.approx(30.0)  # cell 3's only neighbour is cell 2
    assert np.isnan(out[4])  # isolated cell has no neighbours


def test_max_neighbour_intensity_ignores_a_neighbours_missing_value():
    query = query_neighbours(LINE_COORDS, k=2, radius_um=1.5)
    values = np.array([10.0, np.nan, 30.0, 5.0, 999.0])
    out = max_neighbour_intensity(query, values)
    assert out[0] == pytest.approx(np.nan) or np.isnan(out[0])  # only neighbour (1) is NaN


def test_neighbour_fraction_excludes_unavailable_neighbours_from_denominator():
    query = query_neighbours(LINE_COORDS, k=2, radius_um=1.5)
    positive = pd.array([True, False, pd.NA, True, False], dtype="boolean")
    out = neighbour_fraction(query, positive)
    assert out[0] == pytest.approx(0.0)  # neighbour 1 is False
    assert out[1] == pytest.approx(1.0)  # neighbour 0 is True; neighbour 2 (NA) excluded
    assert out[2] == pytest.approx(0.5)  # neighbours 1 (False) and 3 (True)
    assert np.isnan(out[3])  # cell 3's only neighbour (2) is NA -> no known neighbour
    assert np.isnan(out[4])  # no neighbours at all


def test_distance_to_nearest_exceeding_finds_the_planted_bleed_through_neighbour():
    coords = np.array([[float(i), 0.0] for i in range(6)])
    values = np.array([1.0, 1.0, 1.0, 10.0, 1.0, 1.0])  # cell 3 is the lone GCG+ cell
    out = distance_to_nearest_exceeding(coords, values, threshold=5.0)
    assert out[0] == pytest.approx(3.0)
    assert out[1] == pytest.approx(2.0)
    assert out[2] == pytest.approx(1.0)
    assert np.isnan(out[3])  # excludes self by default, and no other cell exceeds threshold
    assert out[4] == pytest.approx(1.0)
    assert out[5] == pytest.approx(2.0)


def test_distance_to_nearest_exceeding_can_include_self():
    coords = np.array([[float(i), 0.0] for i in range(6)])
    values = np.array([1.0, 1.0, 1.0, 10.0, 1.0, 1.0])
    out = distance_to_nearest_exceeding(coords, values, threshold=5.0, exclude_self=False)
    assert out[3] == pytest.approx(0.0)


def test_distance_to_nearest_exceeding_is_unavailable_when_nothing_exceeds_threshold():
    coords = np.array([[float(i), 0.0] for i in range(6)])
    values = np.full(6, 1.0)
    out = distance_to_nearest_exceeding(coords, values, threshold=5.0)
    assert np.isnan(out).all()


def test_distance_to_nearest_exceeding_treats_missing_values_as_not_exceeding():
    coords = np.array([[float(i), 0.0] for i in range(4)])
    values = np.array([np.nan, 1.0, 1.0, 1.0])
    out = distance_to_nearest_exceeding(coords, values, threshold=0.5)
    assert out[0] == pytest.approx(
        1.0
    )  # nearest OTHER qualifying cell, not itself (NaN doesn't qualify)


def test_query_neighbours_rejects_bad_inputs():
    with pytest.raises(ContractError, match="k must be"):
        query_neighbours(LINE_COORDS, k=0, radius_um=1.0)
    with pytest.raises(ContractError, match="radius_um"):
        query_neighbours(LINE_COORDS, k=1, radius_um=0.0)
    with pytest.raises(ContractError, match="coords"):
        query_neighbours(np.zeros((0, 2)), k=1, radius_um=1.0)


def test_distance_to_nearest_exceeding_rejects_a_non_finite_threshold():
    coords = np.array([[0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ContractError, match="threshold"):
        distance_to_nearest_exceeding(coords, np.array([1.0, 2.0]), threshold=float("nan"))


def test_chunking_never_changes_the_result():
    rng = np.random.default_rng(42)
    n = 613
    coords = rng.uniform(0.0, 200.0, size=(n, 2))
    values = rng.normal(1.0, 1.0, n)

    whole = query_neighbours(coords, k=5, radius_um=15.0, chunk_size=n)
    chunked = query_neighbours(coords, k=5, radius_um=15.0, chunk_size=17)
    np.testing.assert_array_equal(whole.indices, chunked.indices)
    np.testing.assert_allclose(whole.distances, chunked.distances, equal_nan=True)

    whole_dist = distance_to_nearest_exceeding(coords, values, threshold=1.5, chunk_size=n)
    chunked_dist = distance_to_nearest_exceeding(coords, values, threshold=1.5, chunk_size=23)
    np.testing.assert_allclose(whole_dist, chunked_dist, equal_nan=True)


def test_result_is_deterministic_across_repeated_calls():
    rng = np.random.default_rng(7)
    n = 240
    coords = rng.uniform(0.0, 100.0, size=(n, 2))
    values = rng.normal(0.0, 1.0, n)

    first = query_neighbours(coords, k=4, radius_um=10.0)
    second = query_neighbours(coords, k=4, radius_um=10.0)
    np.testing.assert_array_equal(first.indices, second.indices)
    np.testing.assert_allclose(first.distances, second.distances, equal_nan=True)

    first_max = max_neighbour_intensity(first, values)
    second_max = max_neighbour_intensity(second, values)
    np.testing.assert_allclose(first_max, second_max, equal_nan=True)

    first_dist = distance_to_nearest_exceeding(coords, values, threshold=0.5)
    second_dist = distance_to_nearest_exceeding(coords, values, threshold=0.5)
    np.testing.assert_allclose(first_dist, second_dist, equal_nan=True)
