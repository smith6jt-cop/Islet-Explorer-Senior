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

import numpy as np
import pandas as pd

from multiplex5000.cohort_scale import SCALE_Q as SCALE_Q
from multiplex5000.cohort_scale import apply_scale as apply_scale
from multiplex5000.cohort_scale import background_mode_floor as background_mode_floor
from multiplex5000.cohort_scale import marker_scale as marker_scale


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


def nuclear_only_references(R: pd.DataFrame, nucleus: tuple) -> np.ndarray:
    """Which references constrain ONLY nucleus-measured markers.

    The DAPI ratio is safe exactly when both sides of a comparison were measured the same
    way.  Ki67 & PCNA survives it (375.7x -> 147.0x) because both are divided by the same
    quantity; every cross-compartment pair does not, and the measured cost is not a
    degradation but an inversion -- partner tails 53.6x -> 0.7x, and the TOX & E-cadherin
    negative control 0.0x -> 10.5x, an association the ratio manufactures.

    So the ratio is applied per REFERENCE, not per cell: a reference whose every constraint
    is on a ``NUCLEUS_MEASURED`` marker is scored on the ratio, and every other reference is
    scored on the aligned values untouched.  On the current 37-reference table that is the
    proliferation split alone.
    """
    nuc = set(nucleus)
    constrained = (R != 0)
    return np.array([
        bool(constrained[c].any()) and set(R.index[constrained[c]]).issubset(nuc)
        for c in R.columns
    ], dtype=bool)


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


def meets_requirements(Xs: np.ndarray, markers: list, reqs: dict, names: list,
                       min_scaled: float) -> np.ndarray:
    """(cells x references) boolean: does each cell show what each rule demands?

    ``Xs`` is already scaled to [0, 1] over [aligned background mode, robust top], so
    "elevated" means the value sits ``min_scaled`` of the way above this marker's own measured
    background.  That is a statement about background, not about prevalence -- deliberately,
    because a fixed quantile would assume how common each cell type is, which is the trap that
    made the DeepCellTypes decile contrast unreadable.  One number for the whole panel, stated
    and swept, rather than a cut per marker.

    ``reqs`` values are ``(all_of, any_groups)`` or ``(all_of, any_groups, none_of)`` -- the
    shapes ``rules.requirements`` returns without and with ``with_negatives``.  A ``none_of``
    marker at or above ``min_scaled`` fails the gate: the rule says the cell must NOT show it,
    and until this was enforced the negatives entered only the signed cosine as a direction,
    so a cell bright in a forbidden marker could still win the reference on angle.

    A marker a donor does not carry cannot be shown, so a reference needing it fails here and
    is scored on the donors that have it.  The same holds for a negative: absence cannot be
    verified on an unmeasured marker, so the reference fails closed rather than passing on a
    constraint nobody checked (a missing measurement is unavailable, not negative).  For the
    same reason a NaN value on a ``none_of`` marker fails the cell -- ``NaN < min_scaled`` is
    False -- exactly as it fails an ``all_of`` marker.  A caller that means to waive a
    constraint on an unmeasured marker removes it from ``reqs`` and records that it did.
    """
    idx = {m: i for i, m in enumerate(markers)}
    hot = Xs >= min_scaled
    cold = Xs < min_scaled                     # NaN is neither hot nor cold
    ok = np.ones((len(Xs), len(names)), dtype=bool)
    for j, ph in enumerate(names):
        spec = reqs.get(ph, ([], []))
        all_of, any_groups = list(spec[0]), list(spec[1])
        none_of = list(spec[2]) if len(spec) > 2 else []
        col = np.ones(len(Xs), dtype=bool)
        unverifiable = False
        for m in all_of:
            if m not in idx:
                unverifiable = True
                break
            col &= hot[:, idx[m]]
        if not unverifiable:
            for grp in any_groups:
                have = [idx[m] for m in grp if m in idx]
                if not have:
                    unverifiable = True
                    break
                col &= hot[:, have].any(axis=1)
        if not unverifiable:
            for m in none_of:
                if m not in idx:
                    unverifiable = True
                    break
                col &= cold[:, idx[m]]
        ok[:, j] = False if unverifiable else col
    return ok


def score_and_gate(Xs, Xs_nuc, R, markers: list, reqs: dict, names: list,
                   min_scaled: float, nuc_ref, ceiling: bool = True):
    """Cosine + requirement gate, with the DAPI ratio used ONLY where it is safe.

    ``Xs`` is the aligned matrix untouched and scores every reference.  ``Xs_nuc`` is the
    same cells with the nucleus-measured group divided by DAPI, and it replaces the score and
    the gate for exactly those references :func:`nuclear_only_references` admits -- i.e. those
    comparing nucleus against nucleus, where both sides carry the same denominator.

    On the current rule table no reference qualifies, because ``rules.flatten`` composes each
    phenotype with its ancestors and every leaf therefore inherits pan-epithelial or
    pan-immune constraints.  That is not a special case to code around: it is the measurement
    saying the ratio has no safe application here, and the ratio is correspondingly inert.
    Add a rule that constrains only nucleus markers and it switches on for that rule alone.
    """
    S = cosine_scores(Xs, R, normalise_ceiling=ceiling)
    allowed = meets_requirements(Xs, markers, reqs, names, min_scaled) if reqs is not None else None
    idx = np.flatnonzero(np.asarray(nuc_ref, dtype=bool)) if nuc_ref is not None else np.array([], int)
    if Xs_nuc is not None and idx.size:
        S[:, idx] = cosine_scores(Xs_nuc, np.asarray(R)[:, idx], normalise_ceiling=ceiling)
        if allowed is not None:
            allowed[:, idx] = meets_requirements(
                Xs_nuc, markers, reqs, [names[i] for i in idx], min_scaled)
    return S, allowed


def call_from_scores(S: np.ndarray, names: list, allowed: np.ndarray | None = None,
                     ) -> pd.DataFrame:
    """Best reference per cell, with the margin to the runner-up.

    With ``allowed`` the gate sets the candidate set and the cosine only ranks within it, so a
    reference cannot win a cell that fails what its rule demands.  A cell meeting nothing is
    ``Ambiguous`` with a recorded reason rather than given its least-bad match.

    The margin is carried because an argmax alone hides how close the decision was, and a
    cell that ties two references is a different object from one that clearly matches one.
    """
    S = np.asarray(S, dtype=np.float32).copy()
    n_allowed = None
    if allowed is not None:
        n_allowed = np.asarray(allowed).sum(axis=1)
        S[~np.asarray(allowed)] = np.nan       # the gate sets the candidates, cosine ranks
    ok = np.isfinite(S).any(axis=1)
    order = np.argsort(-np.nan_to_num(S, nan=-np.inf), axis=1)
    best, second = order[:, 0], order[:, 1] if S.shape[1] > 1 else order[:, 0]
    top = S[np.arange(len(S)), best]
    run = S[np.arange(len(S)), second] if S.shape[1] > 1 else np.zeros(len(S), np.float32)
    out = pd.DataFrame({
        "call": np.where(ok, np.asarray(names, dtype=object)[best], "Unscored"),
        "score": np.where(ok, top, np.nan),
        "runner_up": np.where(ok, np.asarray(names, dtype=object)[second], ""),
        "margin": np.where(ok, top - run, np.nan),
    })
    if n_allowed is not None:
        out["n_candidates"] = n_allowed
        out.loc[(n_allowed == 0), ["call", "runner_up"]] = ["Ambiguous", ""]
        out["reason"] = np.where(n_allowed == 0, "no_reference_requirements_met", "")
    return out
