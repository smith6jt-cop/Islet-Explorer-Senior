"""Step 3: pSAM's scoring, at cell level, on what this workflow already measured.

pSAM (`smith6jt-cop/pSAM`) computes the cosine of every **pixel**'s marker vector against a
set of reference pseudospectra, then thresholds the maps and matches them to segmentation.
What we adopt is the *scoring*, not the pipeline: the cosine is taken per **cell**, against
the REDSEA-corrected, UniFORM-aligned values already on disk.  No images, no Cellpose, no
re-segmentation, and the ``object_id`` contract is untouched.

**Why cosine.** It is magnitude-invariant by construction: scaling a cell's whole vector
leaves the angle unchanged.  That matters here because the largest axis of variation between
cells is not identity but overall brightness -- average pairwise marker correlation
r = 0.16-0.30, forcing PC1 to 21.8-36.3% of variance with 40-47 of 47 loadings sharing a
sign.  And unlike dividing by a chosen marker, it cannot manufacture a pairwise association:
the DAPI test showed division turning a negative control (TOX & E-cadherin) from 0.1x to
14.8x, while a rescaling of the whole vector cannot do that to any pair.  The invariance is
to scaling the whole vector, and the measured axis is "bright in everything except stroma",
so it removes most of that axis rather than all of it.

**Two things the cosine needs that the aligned values do not yet provide.**

* *A common per-marker scale.*  Cosine over a 0-1 reference weights every marker by its own
  numeric size, so a marker with a wide dynamic range would dominate the angle whatever the
  reference says.  pSAM gets this from a per-image ``clip(background, max) / max``.  Here the
  same shape is used with a measured background and a robust top, computed **cohort-wide on
  the aligned values** -- never per donor, which would put back the distribution matching
  Step 2 exists to avoid.
* *Per-donor marker availability.*  The cohort is three panel variants, so a reference naming
  CD38 is scored on the 13 donors that carry it.  The cosine is taken over the markers that
  donor actually has, and the number of donors a phenotype was scored on is recorded.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402

def nuclear_ratio(df, markers: list, denominator: str = "DAPI",
                  nuclear: tuple = ()) -> tuple:
    """Divide each raw-nucleus marker by DAPI, in place, returning the columns changed.

    Eight markers reach ``expression`` as UNCOMPENSATED *nucleus* means
    (``qupath:Nucleus__<m>``) while the other 52 are REDSEA spillover-corrected whole-cell
    means.  Mixing the two in one cosine compares quantities that were not measured the same
    way.

    REDSEA is compartment-aware and computed Nucleus / Cytoplasm / Membrane / Cell for all 60
    markers.  What it did not do for these is apply the spillover compensation:
    ``marker_taxonomy.NO_SPILLOVER`` lists Collagen_IV, DAPI, FOXP3, Ki67, PCNA, SOX2, TCF_1,
    TOX and TP63, and ``config.ini`` sets ``exclude_no_spillover = true``.  So the REDSEA nucleus
    value for them is uncompensated too -- measured on donor 115 it differs from the QuPath
    nucleus mean by a median 2%, re-rasterisation rather than correction, and the two carry the
    same tail evidence (partner medians 55.6 against 53.6, controls 0.2 against 0.2).  Which is
    why ``expression`` selects the QuPath value and why substituting REDSEA's would change
    nothing measurable.
      DAPI is the natural internal reference for that group: same compartment, same
    optical path, present on every nucleus, and measured with the same exposure on every
    slide in this cohort.

    The ratio removes what the group shares -- nuclear content in the optical section,
    staining penetration, how much nucleus fell inside the mask.  What it also removes is
    real: DAPI carries nuclear size and chromatin density, which are cell-type traits.

    **Measured on the aligned values, 5 donors, top 0.1% tails.**  It strips the magnitude
    coupling as intended (TOX's Spearman against the per-cell magnitude axis 0.63 -> 0.09)
    and it takes the lineage evidence with it: partner tails fall from a median 53.6x to
    0.7x -- TOX & CD3e 105.1 -> 0.7, TCF_1 & CD3e 53.6 -> 0.2, FOXP3 & CD3e 48.6 -> 6.0 --
    while the TOX & E-cadherin negative control goes 0.0 -> 10.5, an association the ratio
    invents rather than reveals.  Ki67 & PCNA survives best at 375.7 -> 147.0 because BOTH
    are divided by the same quantity, which is the shape of the one case where this is safe:
    nuclear against nuclear.  Against a REDSEA whole-cell marker it is not.

    Regressing on ``log(DAPI)`` at its fitted slope (0.46-0.72 across the group) is milder
    and still leaves that control at 13.8x, so the problem is not the exponent.
    """
    if denominator not in df.columns:
        return df, []
    d = df[denominator].to_numpy("float64")
    d = np.where(np.isfinite(d) & (d > 0), d, np.nan)
    changed = [m for m in markers if m in nuclear and m != denominator and m in df.columns]
    for m in changed:
        df[m] = df[m].to_numpy("float64") / d
    return df, changed


#: Upper point of the per-marker scale.  A high quantile rather than the maximum: pSAM
#: divides by the image max, which one hot pixel can set.
SCALE_Q = 99.9


def marker_scale(values_by_marker: dict, floor_by_marker: dict | None = None,
                 q: float = SCALE_Q) -> pd.DataFrame:
    """Cohort-wide floor and top per marker, so every marker enters the angle comparably."""
    rows = []
    for m, v in values_by_marker.items():
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v)]
        floor = float(floor_by_marker.get(m, 0.0)) if floor_by_marker else 0.0
        top = float(np.percentile(v, q)) if len(v) else np.nan
        rows.append({"marker": m, "floor": floor, "top": top,
                     "span": top - floor if np.isfinite(top) else np.nan})
    return pd.DataFrame(rows).set_index("marker")


def apply_scale(X: np.ndarray, markers: list, scale: pd.DataFrame) -> np.ndarray:
    """(cells x markers) -> [0, 1], clipped at the measured floor and the robust top."""
    lo = scale.loc[markers, "floor"].to_numpy(float)
    span = scale.loc[markers, "span"].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    return np.clip((np.asarray(X, dtype=np.float32) - lo) / span, 0.0, 1.0)


def reference_ceiling(R: np.ndarray) -> np.ndarray:
    """The largest cosine each reference can attain against a NON-NEGATIVE cell vector.

    Cell values are scaled to [0, 1], so the best a cell can do against reference ``r`` is to
    align with its positive part alone: ``max cos = ||r_positive|| / ||r||``.  A reference
    with no negatives can reach 1; one with many cannot get close.  Measured on donor 6450,
    that ceiling ran from 1.000 down to 0.147 across the 37 references, and it PREDICTED the
    call share -- Spearman +0.48.  The argmax was partly choosing on reference geometry
    rather than on the cell.
    """
    R = np.asarray(R, dtype=np.float32)
    rn = np.linalg.norm(R, axis=0)
    return np.divide(np.linalg.norm(np.clip(R, 0, None), axis=0), rn,
                     out=np.zeros_like(rn), where=rn > 0)


def cosine_scores(X: np.ndarray, R: np.ndarray, normalise_ceiling: bool = True) -> np.ndarray:
    """(cells x markers) against (markers x references) -> (cells x references).

    With ``normalise_ceiling`` each column is divided by what that reference could attain at
    best, so the references compete on equal terms and a score of 1 means "as close to this
    reference as any non-negative cell could be".  Without it, references carrying negative
    evidence are structurally unable to win.

    A cell with no signal at all has a zero norm and no angle; those rows come back NaN
    rather than 0, so "nothing to compare" is distinguishable from "compared and unlike".
    """
    X = np.asarray(X, dtype=np.float32)
    R = np.asarray(R, dtype=np.float32)
    xn = np.linalg.norm(X, axis=1)
    rn = np.linalg.norm(R, axis=0)
    rn = np.where(rn > 0, rn, np.nan)
    S = (X @ R) / (xn[:, None] * rn[None, :])
    if normalise_ceiling:
        c = reference_ceiling(R)
        S = S / np.where(c > 0, c, np.nan)[None, :]
    S[xn == 0] = np.nan
    return S


def call_from_scores(S: np.ndarray, names: list) -> pd.DataFrame:
    """Best reference per cell, with the margin to the runner-up.

    The margin is carried because an argmax alone hides how close the decision was, and a
    cell that ties two references is a different object from one that clearly matches one.
    """
    S = np.asarray(S, dtype=np.float32)
    ok = np.isfinite(S).any(axis=1)
    order = np.argsort(-np.nan_to_num(S, nan=-np.inf), axis=1)
    best, second = order[:, 0], order[:, 1] if S.shape[1] > 1 else order[:, 0]
    top = S[np.arange(len(S)), best]
    run = S[np.arange(len(S)), second] if S.shape[1] > 1 else np.zeros(len(S), np.float32)
    return pd.DataFrame({
        "call": np.where(ok, np.asarray(names, dtype=object)[best], "Unscored"),
        "score": np.where(ok, top, np.nan),
        "runner_up": np.where(ok, np.asarray(names, dtype=object)[second], ""),
        "margin": np.where(ok, top - run, np.nan),
    })
