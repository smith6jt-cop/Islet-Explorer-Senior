"""Shared fixtures for the broadtype spine tests.  Synthetic data only -- no cohort data.

Two things live here because more than one test module needs them:

* the import path.  The spine modules import each other as bare names (``import pins``) and
  ``pins`` puts ``scripts/common``, ``scripts/celltype`` and ``scripts/dct`` on ``sys.path``,
  so a test only has to see ``scripts/broadtype`` itself.
* a synthetic donor with PLANTED labels and PLANTED marker tails, and a toy arm that reads
  them back from expression alone.  The gate is a set of measurements against nulls; the only
  way to know the bars measure the ARM and not the gate is to hand the gate a case whose
  answer is known in advance and demand it be verified.  The planted structure is the answer.

The planted donor is deliberately artificial (unit-variance noise plus a fixed shift on a
label's declared markers) -- a positive control, not a tissue model.  Nothing in it, and
nothing in the toy arm, consults a frequency: shares are an INPUT the tests scramble on
purpose, and the arm is an argmax over declared-marker signatures.
"""
from __future__ import annotations

import pathlib
import sys
import uuid
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

BROADTYPE = pathlib.Path(__file__).resolve().parents[1]
if str(BROADTYPE) not in sys.path:
    sys.path.insert(0, str(BROADTYPE))

import declarations  # noqa: E402
from records import UNAVAILABLE  # noqa: E402


@pytest.fixture(scope="session")
def dec() -> dict:
    """The frozen declarations, as every arm and the gate read them."""
    return declarations.load()


# --------------------------------------------------------------------------- planted donor
@dataclass
class PlantedDonor:
    ids: np.ndarray            # object_id (uuid4 strings), the frozen cell set
    X: pd.DataFrame            # aligned values, columns = the declared identity markers
    labels: np.ndarray         # the planted truth, one declared label per cell
    xy: np.ndarray             # level-0 pixel centroids
    groups: np.ndarray         # region index (the frame's ``region_idx`` analogue)
    cell_region: np.ndarray    # a synthetic QuPath parentage column -- AUDIT ONLY
    shares: dict               # the requested share per label (an input, not a target)

    @property
    def n(self) -> int:
        return len(self.ids)


def _uuid4(rng: np.random.Generator) -> str:
    return str(uuid.UUID(bytes=rng.bytes(16), version=4))


def _patch_labels(rng, xy, groups, labels, shares, patches_per_region):
    """Spatially coherent labels: nearest-patch-centre membership, patches filled to quota.

    Every label gets whole patches until its quota (share x n) is met, in a random patch
    order, so a small share still gets at least one patch and the achieved share tracks the
    requested one.  Leftover patches go to whichever label is furthest below quota.
    """
    n = len(xy)
    patch = np.empty(n, dtype=int)
    offset = 0
    for g in np.unique(groups):
        sel = np.flatnonzero(groups == g)
        centres = xy[rng.choice(sel, patches_per_region, replace=False)]
        d = ((xy[sel, None, :] - centres[None, :, :]) ** 2).sum(-1)
        patch[sel] = d.argmin(1) + offset
        offset += patches_per_region
    sizes = np.bincount(patch, minlength=offset)
    order = rng.permutation(offset)
    quota = {L: shares[L] * n for L in labels}
    got = {L: 0.0 for L in labels}
    patch_label = np.empty(offset, dtype=object)
    li = 0
    for p in order:
        while li < len(labels) - 1 and got[labels[li]] >= quota[labels[li]]:
            li += 1
        L = labels[li]
        if got[L] >= quota[L]:                      # every quota met: fill the largest deficit
            L = min(labels, key=lambda k: got[k] - quota[k])
        patch_label[p] = L
        got[L] += sizes[p]
    return patch_label[patch].astype(str)


def make_planted_donor(dec: dict, *, n: int = 20_000, n_regions: int = 4, shares: dict | None = None,
                       seed: int = 0, spatial: bool = True, signal: float = 3.0,
                       patches_per_region: int = 40) -> PlantedDonor:
    """A synthetic donor whose labels and marker tails are known in advance.

    ``shares`` maps declared label -> share of cells (defaults to equal).  With ``spatial``
    the labels form patches inside each region; without it they are drawn cell by cell with
    no structure.  ``signal`` is the shift added to a cell's OWN label's declared markers
    (0 plants no tails at all).
    """
    rng = np.random.default_rng(seed)
    labels_all = list(dec["labels"])
    markers = list(dec["identity_markers"])
    if shares is None:
        shares = {L: 1.0 / len(labels_all) for L in labels_all}
    if abs(sum(shares.values()) - 1.0) > 1e-9:
        raise ValueError("shares must sum to 1")
    per = n // n_regions
    groups = np.repeat(np.arange(n_regions), per)
    n = len(groups)
    # regions are 1000-px squares on a 2 x 2 grid, far enough apart that no kNN edge crosses
    xy = np.column_stack([rng.uniform(0, 1000, n) + (groups % 2) * 5000,
                          rng.uniform(0, 1000, n) + (groups // 2) * 5000])
    if spatial:
        labels = _patch_labels(rng, xy, groups, labels_all, shares, patches_per_region)
    else:
        labels = rng.choice(labels_all, n, p=[shares[L] for L in labels_all]).astype(str)
    X = pd.DataFrame(rng.normal(size=(n, len(markers))), columns=markers)
    if signal:
        for L in labels_all:
            sel = labels == L
            X.loc[sel, dec["labels"][L]["markers"]] += signal
    ids = np.array([_uuid4(rng) for _ in range(n)], dtype=object)
    # a synthetic QuPath parentage column, so the audit-only path is exercised; it is built
    # from the planted truth only to have SOMETHING enriched to audit, and the tests assert
    # the verdict is identical with and without it
    cell_region = np.where((labels == "Endocrine") & (rng.random(n) < 0.8), "core",
                           np.where(rng.random(n) < 0.05, "core", "tissue")).astype(str)
    return PlantedDonor(ids=ids, X=X, labels=labels, xy=xy, groups=groups, cell_region=cell_region,
                        shares=dict(shares))


def signature_arm(X: pd.DataFrame, dec: dict) -> np.ndarray:
    """The toy arm: argmax over per-label mean z-score of the declared markers PRESENT in X.

    It reads nothing but the aligned values it is handed, so an ablation (a column removed
    from X) genuinely changes its input -- which is what the held-out criterion measures.
    """
    Z = (X - X.mean()) / X.std(ddof=0).replace(0, 1.0)
    Z = Z.fillna(0.0)
    scores = {}
    for L, spec in dec["labels"].items():
        cols = [m for m in spec["markers"] if m in Z.columns]
        if cols:
            scores[L] = Z[cols].mean(axis=1).to_numpy()
    S = pd.DataFrame(scores)
    return S.columns.to_numpy()[S.to_numpy().argmax(axis=1)].astype(str)


def subsample_replicate(labels: np.ndarray, ids: np.ndarray, frac: float, salt: str) -> np.ndarray:
    """A 90 %-subsample replicate as the screen driver aligns it: cells outside the draw are
    Unavailable so every replicate stays aligned to the frozen ``ids``."""
    from gate import hash_subsample
    keep = hash_subsample(ids, int(round(frac * len(ids))), salt=salt)
    out = np.full(len(labels), UNAVAILABLE, dtype=object)
    out[keep] = labels[keep]
    return out.astype(str)


@pytest.fixture
def planted_donor(dec):
    """Factory: ``planted_donor(**kwargs) -> PlantedDonor``."""
    def _make(**kw):
        return make_planted_donor(dec, **kw)
    return _make
