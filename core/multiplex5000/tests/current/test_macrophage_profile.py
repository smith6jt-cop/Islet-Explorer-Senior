"""Owner-approved macrophage labels use measured signs and one hierarchy."""

import numpy as np
import pandas as pd
import pytest

from multiplex5000.tissues.macrophage import classify_profile, labels


def test_twelve_labels_exclude_double_negative_parent():
    all_labels = labels()
    assert len(all_labels) == len(set(all_labels)) == 12
    assert all("Macrophage_LGALS3" in label for label in all_labels)
    assert all(not label.startswith("CD68-CD163-") for label in all_labels)


def test_profile_never_turns_missing_marker_into_negative():
    scaled = pd.DataFrame({
        "CD68": [0.11, 0.0, 0.11, 0.0],
        "CD163": [0.0, 0.2, 0.0, 0.0],
        "M2Gal3": [0.2, 0.0, np.nan, 0.4],
        "CD4": [0.0, 0.3, 0.3, 0.0],
    })
    assert classify_profile(scaled).tolist() == [
        "CD68+CD163-Macrophage_LGALS3+CD4-",
        "CD68-CD163+Macrophage_LGALS3-CD4+",
        None, None,
    ]


def test_profile_uses_image_specific_marker_levels():
    values = pd.DataFrame({
        "CD68": [0.08, 0.06], "CD163": [0.22, 0.26],
        "M2Gal3": [0.05, 0.03], "CD4": [0.17, 0.21],
    })
    levels = {"CD68": 0.07, "CD163": 0.25, "M2Gal3": 0.04, "CD4": 0.20}
    assert classify_profile(values, levels).tolist() == [
        "CD68+CD163-Macrophage_LGALS3+CD4-",
        "CD68-CD163+Macrophage_LGALS3-CD4+",
    ]
    with pytest.raises(ValueError, match="missing markers"):
        classify_profile(values, {"CD68": 0.07})
    with pytest.raises(ValueError, match="must be finite"):
        classify_profile(values, {**levels, "CD4": np.nan})
