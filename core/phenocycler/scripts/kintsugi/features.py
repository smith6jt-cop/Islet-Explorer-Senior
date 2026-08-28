"""Step 1: what goes into the clustering, on what scale, and the check that the scale is honest.

The SOM's distance is Euclidean over the marker block, so two decisions determine everything
downstream: which markers are in it, and what units they are in.

**Which markers.**  Derived, never hand-listed.  ``panel.markers()`` is the cohort
intersection; ``panel.REGISTRY_EXCLUDED`` drops IAPP because the registry records it
``kind: excluded, uses: failed_marker``; DAPI is dropped because the registry records it
``kind: structural`` and it is on every nucleus, so it cannot separate anything; and
``panel.PANEL_VARIANT_MARKERS`` (CD38 on 13 donors, b_Catenin1 on 4) are dropped because a
cohort-wide feature block has no honest value for a donor that never acquired them -- a
missing measurement is unavailable, not zero.

**CD11b stays.**  The registry flags it ``assay_qc``, but ``scripts/celltype/panel.py`` keeps
it and says why: it was set aside for the retired RESTORE lane, and the rule table's Myeloid
gate is ``CD68 | CD11c | CD11b``, so dropping it silently removes a third of that gate's
evidence.  This lane answers against the same rule table, so it makes the same choice -- and
it has to, or the concordance measurement compares two lanes that were not given the same
evidence.

**What units.**  Raw expression is not comparable between donors: donor-median fold-range is
4.2x across the panel.  So the block is built from Step 2's UniFORM-aligned parquet, which is
the repo's sanctioned cross-donor alignment -- a per-(donor, marker) multiplicative rescale,
verified pure to 2-9e-16 over all 31.5M rows.

**The scaling tension, and how it is settled.**  ``Kanalysis.visualize_som`` calls
``StandardScaler().fit_transform`` on the pooled block, and ``docs/INTENSITY_SCALE.md``
permanently forbids "z-scoring across donors".  These are not the same operation.  UniFORM has
already done the cross-donor alignment; the scaler that follows it equalises **marker**
variance so a bright marker does not dominate Euclidean distance.  That is feature weighting,
not donor normalisation.

The distinction is checkable rather than rhetorical, and :func:`scaler_is_feature_weighting`
checks it: a per-marker affine transform applied to the pooled data moves every donor's median
by the same amount, so the *spread of donor medians* -- measured in units of the marker's own
pooled SD -- is invariant.  If the scaler were correcting donors, that spread would shrink.
The check fails the arm rather than warning.

**The twelve that alignment cannot fix.**  ``shifts.json`` marks 12 markers ``over_bar``: they
differ between donors in distribution *shape*, not location, so a rigid shift cannot align
them and ``celltype/README.md`` records that only within-donor rank statistics are safe on
them.  Six of the twelve -- INS, E_cadherin, Pan_Cytokeratin, Ker8_18, CD31, Vimentin -- are
the epithelial, endocrine and endothelial backbone of the rule table, so dropping them is not
a neutral act; it would remove the axis the annotation step needs most.  Three arms are
therefore declared in advance and decided by measurement, not by preference:

===  ===========================================================================
arm  the twelve
===  ===========================================================================
A    aligned value as-is -- accepts the residual, keeps the axis
B    within-donor rank -- comparable by construction, at the cost of shape
C    excluded -- the conservative reading of "only rank statistics are safe"
===  ===========================================================================

CD57 is excluded from all three: ``shifts.json`` records it
``provenance: refused_reagent_change`` with an acquisition-group ratio of 35.1, and a rigid
shift models a per-donor gain, not a different antibody clone.

The arms are judged on two measurements declared before they are run -- cluster donor-purity
(:func:`donor_purity`) and concordance with the cosine lane -- and the bars are not relaxed
after the answer is seen.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
import panel                                        # noqa: E402

UNIFORM = panel.OUT / "uniform"
SHIFTS = panel.OUT / "shifts.json"

#: Registry ``kind: structural``, on every nucleus.  Notebook 4 drops DAPI from the cell
#: compartment block for the same reason ("Remove DAPI and any markers not expressed in the
#: cell ... compartment").
STRUCTURAL = ("DAPI",)

ARMS = ("A", "B", "C")


def shift_classes() -> dict[str, list[str]]:
    """Derive the marker classes from Step 2's own record, so they cannot drift.

    Hand-listing the twelve would let this lane and ``celltype/README.md`` disagree silently
    after a re-run of Step 2.  ``shifts.json`` is the thing that decided them.
    """
    shifts = json.loads(SHIFTS.read_text())["shifts"]
    return {
        "over_bar": sorted(m for m, v in shifts.items() if v.get("over_bar")),
        "refused": sorted(m for m, v in shifts.items()
                          if v.get("provenance") == "refused_reagent_change"),
        "aligned": sorted(m for m, v in shifts.items()
                          if v.get("provenance") == "automatic_registration"),
    }


def feature_markers(arm: str) -> list[str]:
    """The marker block for one arm, derived from the registry, the panel and shifts.json."""
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, got {arm!r}")
    classes = shift_classes()
    drop = set(panel.REGISTRY_EXCLUDED) | set(STRUCTURAL) \
        | set(panel.PANEL_VARIANT_MARKERS) | set(classes["refused"])
    if arm == "C":
        drop |= set(classes["over_bar"])
    return [m for m in panel.markers(include_variant=False) if m not in drop]


def within_donor_rank(values: np.ndarray, donors: np.ndarray) -> np.ndarray:
    """Rank each donor's cells among themselves, scaled to [0, 1].

    A per-donor monotone transform, so it preserves within-donor order exactly -- which is
    the statistic ``celltype/README.md`` says is safe on the shape-mismatched markers.  It
    does make each donor's marginal uniform, which is a real cost: it discards how far apart
    the positive and negative populations sit.  That is why it is arm B and not the default.
    """
    out = np.empty(len(values), dtype=np.float64)
    for d in np.unique(donors):
        m = donors == d
        v = values[m]
        order = v.argsort(kind="mergesort").argsort()
        out[m] = order / max(len(v) - 1, 1)
    return out


def scaler_is_feature_weighting(values: np.ndarray, donors: np.ndarray,
                                scaled: np.ndarray, tol: float = 1e-6) -> dict:
    """Prove the scaler weighted markers rather than correcting donors.

    A per-marker affine map shifts every donor's median identically, so the spread of donor
    medians expressed in units of the marker's pooled SD is unchanged.  Donor correction would
    shrink it.  Reports the largest relative change over the markers; the caller refuses the
    arm if it exceeds ``tol``.
    """
    worst, worst_marker = 0.0, None
    for j in range(values.shape[1]):
        before = np.array([np.median(values[donors == d, j]) for d in np.unique(donors)])
        after = np.array([np.median(scaled[donors == d, j]) for d in np.unique(donors)])
        sd_before, sd_after = values[:, j].std(), scaled[:, j].std()
        if sd_before <= 0 or sd_after <= 0:
            continue
        rel = abs(before.std() / sd_before - after.std() / sd_after)
        if rel > worst:
            worst, worst_marker = float(rel), j
    return {"max_relative_change": worst, "worst_marker_index": worst_marker,
            "tolerance": tol, "passes": worst <= tol}


def donor_purity(clusters: np.ndarray, donors: np.ndarray) -> pd.DataFrame:
    """Per-cluster donor-composition entropy against the cohort's.

    The cohort's islet-core immune gradient is aliased with acquisition batch -- batch 2 is
    every high-immune donor and no ND donor -- so a clustering that separates donors rather
    than cell types would reproduce that alias and look like biology.  ``purity`` near 1 means
    a cluster is one donor; near 0 means it mirrors the cohort.  This is reported per arm and
    is one of the two pre-declared measurements the arms are judged on.
    """
    frame = pd.DataFrame({"cluster": clusters, "donor": donors})
    cohort_share = frame["donor"].value_counts(normalize=True)
    max_entropy = -np.sum(cohort_share * np.log(cohort_share))
    rows = []
    for cluster, part in frame.groupby("cluster"):
        share = part["donor"].value_counts(normalize=True)
        entropy = -np.sum(share * np.log(share))
        rows.append({"cluster": cluster, "n_cells": len(part),
                     "n_donors": part["donor"].nunique(),
                     "entropy": float(entropy),
                     "purity": float(1.0 - entropy / max_entropy) if max_entropy else np.nan,
                     "top_donor": share.index[0],
                     "top_donor_share": float(share.iloc[0])})
    return pd.DataFrame(rows).sort_values("purity", ascending=False).reset_index(drop=True)


CONTEXT = ["donor_id", "object_id", "qc_analysis_eligible", "cell_region", "islet_num"]


def load_block(donors: list[str], markers: list[str], arm: str) -> tuple[pd.DataFrame, dict]:
    """The aligned block for these donors, QC-eligible cells only, with what was dropped.

    ``qc_analysis_eligible`` is honoured rather than ignored: a cell the geometry stage
    refused is not a cell this lane should let steer a SOM node.  Cells carrying a
    non-finite value in any feature marker are dropped and counted -- never imputed, because
    a missing measurement is unavailable, not zero.
    """
    frames, dropped_qc, dropped_nonfinite = [], 0, 0
    for donor in donors:
        table = pq.read_table(UNIFORM / f"donor_id={donor}" / "data.parquet",
                              columns=CONTEXT + markers).to_pandas()
        before = len(table)
        table = table[table["qc_analysis_eligible"]]
        dropped_qc += before - len(table)
        finite = np.isfinite(table[markers].to_numpy()).all(axis=1)
        dropped_nonfinite += int((~finite).sum())
        frames.append(table[finite])
    block = pd.concat(frames, ignore_index=True)
    if arm == "B":
        donors_arr = block["donor_id"].to_numpy()
        for marker in shift_classes()["over_bar"]:
            if marker in markers:
                block[marker] = within_donor_rank(
                    block[marker].to_numpy(), donors_arr)
    return block, {"dropped_qc_ineligible": dropped_qc,
                   "dropped_nonfinite": dropped_nonfinite,
                   "n_cells": len(block)}
