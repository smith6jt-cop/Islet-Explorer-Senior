"""The SOM port must reproduce KINTSUGI's own function, and recover known structure."""
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import som_adapter                                  # noqa: E402

pyflowsom = pytest.importorskip("pyFlowSOM")


def _three_populations(n=900, seed=0):
    """Three well-separated populations over six markers; the answer is known."""
    rng = np.random.default_rng(seed)
    markers = [f"M{i}" for i in range(6)]
    blocks, truth = [], []
    for k in range(3):
        centre = np.zeros(6)
        centre[2 * k:2 * k + 2] = 10.0
        blocks.append(rng.normal(centre, 0.5, size=(n // 3, 6)))
        truth += [k] * (n // 3)
    return np.vstack(blocks), markers, np.array(truth)


def test_som_recovers_known_populations():
    values, markers, truth = _three_populations()
    clusters, profile = som_adapter.som_cluster(values, markers, num_nodes=9, num_passes=5)
    # Every SOM node must be pure in the true population it came from.
    for node in np.unique(clusters):
        assert len(set(truth[clusters == node])) == 1, f"node {node} mixes populations"
    assert len(profile) == len(np.unique(clusters))


def test_profile_is_on_the_input_scale_not_the_scaled_one():
    """The annotation gate is stated on the aligned scale, so the profile must be too."""
    values, markers, truth = _three_populations()
    clusters, profile = som_adapter.som_cluster(values, markers, num_nodes=9, num_passes=5)
    assert profile.to_numpy().max() > 5.0, "profile looks z-scored, not on the input scale"
    for node in np.unique(clusters):
        expected = values[clusters == node].mean(axis=0)
        got = profile.loc[f"Cluster_{node + 1}"].to_numpy()
        assert np.allclose(got, expected, atol=1e-6)


def test_matches_the_last_frame_of_visualize_soms_progressive_loop():
    """visualize_som refits per iteration; one fit at rlen=num_passes is the same answer."""
    from sklearn.preprocessing import StandardScaler

    values, markers, _ = _three_populations()
    passes = 4
    ours, _ = som_adapter.som_cluster(values, markers, num_nodes=9, num_passes=passes)

    train = StandardScaler().fit_transform(values.astype(np.float64))
    node_map = None
    for rlen in range(1, passes + 1):               # the notebook's loop, verbatim
        som = pyflowsom.som(train, 3, 3, rlen, alpha_range=(
            som_adapter.NOTEBOOK_SOM["lr_start"], som_adapter.NOTEBOOK_SOM["lr_end"]),
            seed=som_adapter.NOTEBOOK_SOM["seed"])
        node_map, _ = pyflowsom.map_data_to_nodes(som, train)
    unique = np.unique(node_map)
    remap = {old: new for new, old in enumerate(unique)}
    theirs = np.array([remap[c] for c in node_map])
    assert np.array_equal(ours, theirs)


def test_non_finite_input_is_refused_not_imputed():
    values, markers, _ = _three_populations()
    values[0, 0] = np.nan
    with pytest.raises(ValueError, match="unavailable, not zero"):
        som_adapter.som_cluster(values, markers, num_nodes=9, num_passes=2)


def test_non_square_node_count_is_refused():
    values, markers, _ = _three_populations()
    with pytest.raises(ValueError, match="not square"):
        som_adapter.som_cluster(values, markers, num_nodes=10, num_passes=2)


def test_meta_cluster_collapses_by_correlation_and_refuses_bad_counts():
    import pandas as pd
    profile = pd.DataFrame(np.array([[10.0, 0, 0], [10.1, 0, 0], [0, 0, 10.0], [0, 0, 9.9]]),
                           index=[f"Cluster_{i}" for i in range(1, 5)], columns=list("abc"))
    mapping = som_adapter.meta_cluster(profile, 2)
    assert mapping["Cluster_1"] == mapping["Cluster_2"]
    assert mapping["Cluster_3"] == mapping["Cluster_4"]
    assert mapping["Cluster_1"] != mapping["Cluster_3"]
    with pytest.raises(ValueError):
        som_adapter.meta_cluster(profile, 99)
