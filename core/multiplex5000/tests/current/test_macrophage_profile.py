"""Regressions for the tissue-specific macrophage gate and two Gal3 classes."""

import numpy as np
import pandas as pd

from multiplex5000.tissues.macrophage import classify_profile, current_profile_names, labels


def test_positive_anchor_determines_membership_and_only_gal3_splits_the_parent():
    pos, neg = labels()
    assert [pos, neg] == ["Macrophage_Gal3+", "Macrophage_Gal3-"]
    frame = pd.DataFrame({
        "CD68": [1, 0, 1, 1, np.nan, 0, 0, np.nan, 1],
        "CD163": [0, 1, 1, np.nan, 1, 0, np.nan, np.nan, 1],
        "M2Gal3": [1, 1, 0, 0, 1, 1, 0, 1, np.nan],
    })
    expected = [pos, pos, neg, neg, pos, None, None, None, None]
    assert classify_profile(frame).tolist() == expected
    assert classify_profile(frame.assign(CD4=np.arange(len(frame)))).tolist() == expected
    assert classify_profile(frame.drop(columns="CD163")).tolist() == [
        pos, None, neg, neg, None, None, None, None, None]


def test_image_specific_thresholds_on_available_evidence():
    frame = pd.DataFrame({"CD68": [.08, .06], "CD163": [.22, .26], "M2Gal3": [.05, .03]})
    thresholds = {"CD68": .07, "CD163": .25, "M2Gal3": .04}
    pos, neg = labels()
    assert classify_profile(frame, thresholds).tolist() == [pos, neg]
    assert classify_profile(frame, {"CD163": .25, "M2Gal3": .04}).tolist() == [None, neg]
    assert classify_profile(frame, {"CD68": .07}).tolist() == [None, None]
    assert classify_profile(frame, {**thresholds, "M2Gal3": np.nan}).tolist() == [None, None]


def test_archived_names_pool_by_gal3_only():
    prior = pd.Series(["CD68+CD163-Macrophage_LGALS3+CD4-", "CD68-CD163+Macrophage_LGALS3+",
                       "CD68+CD163+Macrophage_LGALS3-CD4+", "Macrophage_Gal3-", None])
    result = current_profile_names(prior)
    assert result.iloc[:4].tolist() == [labels()[0], labels()[0], labels()[1], labels()[1]]
    assert pd.isna(result.iloc[4])
