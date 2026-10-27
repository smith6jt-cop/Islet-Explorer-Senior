"""Owner-approved macrophage labels use measured signs and one hierarchy."""

import numpy as np
import pandas as pd

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
