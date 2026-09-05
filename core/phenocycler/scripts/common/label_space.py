"""The label space: which of the model's 51 classes may be assigned at all.

DeepCellTypes emits a 51-way softmax over a vocabulary built for 17 tissue types.
Twenty of those classes cannot occur in a human pancreas section and five more have
no supporting marker anywhere in this panel, so 25 of the 51 are classes the argmax
can reach on no evidence -- and it does: on donor 6539 the 2026-06-15 checkpoint puts
a mean 0.239 of its probability mass there, and for 34.1 % of cells that mass exceeds
whatever survives.

Restricting the label space is **exactly** a renormalization of the returned softmax
over the kept columns, which is algebraically identical to masking the excluded logits
before the softmax (softmax is a ratio, so the shared denominator cancels).  There is
no model surgery here and no new inference; ``restrict`` runs on a probability matrix
``predict(..., return_probabilities=True)`` already produced.

The keep/drop decision is NOT made here.  It lives in ``broad_class_crosswalk.json``
as ``in_label_space``, derived from ``expected_in_pancreas`` (frozen 2026-08-21, three
days before any inference output existed) and ``panel_supported`` (judged against the
51-marker cross-donor intersection).  Both are decidable by someone who has never seen
the model's output; that is what keeps this a configuration choice and not tuning.
Read ``policy`` in that file before changing anything.

Nothing is discarded silently.  Every cell keeps ``impossible_mass``, its
``restricted_call`` and its ``unrestricted_call``, and a cell whose removed mass
exceeds its surviving top class is ``ambiguous`` -- a comparison between two quantities
the model produced, not a tuned cut.
"""
from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass

import numpy as np

__all__ = ["LabelSpace", "Restricted", "load_label_space", "AMBIGUOUS"]

HERE = pathlib.Path(__file__).resolve().parent
CROSSWALK = HERE.parent / "dct" / "broad_class_crosswalk.json"  # the DCT lane owns this artifact (module moved to scripts/common, Track E.3)

#: Sentinel for a cell the restriction cannot call.  Never "Unassigned": every use
#: carries an enumerated reason alongside it.
AMBIGUOUS = "Ambiguous"

_TOL = 1e-4


@dataclass(frozen=True)
class Restricted:
    """One arm's calls under the restricted label space, plus what was removed."""

    probabilities: np.ndarray   # (n, n_kept) renormalized; NaN row where no mass survived
    kept: tuple                 # column names of `probabilities`
    restricted_call: np.ndarray # (n,) str -- AMBIGUOUS where unavailable
    unrestricted_call: np.ndarray
    impossible_mass: np.ndarray # (n,) probability removed
    ambiguous_reason: np.ndarray  # (n,) str, "" where the cell was called

    def __len__(self) -> int:
        return len(self.impossible_mass)

    @property
    def ambiguous(self) -> np.ndarray:
        return self.restricted_call == AMBIGUOUS


@dataclass(frozen=True)
class LabelSpace:
    classes: tuple              # all 51, in the order the caller's columns are in
    kept: tuple
    excluded: tuple
    exclusion_reason: dict      # excluded class -> code
    broad_class: dict           # class -> broad compartment
    method_version: str

    def broad_groups(self, *, restricted: bool = True) -> dict:
        """compartment -> the classes in it, from the kept set or the full 51."""
        pool = self.kept if restricted else self.classes
        out: dict = {}
        for c in pool:
            out.setdefault(self.broad_class[c], []).append(c)
        return {k: tuple(v) for k, v in out.items()}

    def restrict(self, probabilities, classes=None) -> Restricted:
        """Renormalize a (n, 51) softmax over the kept columns.

        Fails closed: a non-finite entry, a row that is not a probability vector, or a
        column list that is not this vocabulary is an error, never a silent pass.
        """
        P = np.asarray(probabilities, dtype=np.float64)
        cols = tuple(classes) if classes is not None else self.classes
        if P.ndim != 2 or P.shape[1] != len(cols):
            raise ValueError(f"probabilities {P.shape} do not match {len(cols)} class columns")
        if set(cols) != set(self.classes):
            missing = sorted(set(self.classes) - set(cols))
            extra = sorted(set(cols) - set(self.classes))
            raise ValueError(f"column set is not the vocabulary; missing={missing} extra={extra}")
        if not np.isfinite(P).all():
            raise ValueError("probability matrix carries non-finite values -- refusing to renormalize")
        total = P.sum(axis=1)
        if not np.allclose(total, 1.0, atol=_TOL):
            bad = int((~np.isclose(total, 1.0, atol=_TOL)).sum())
            raise ValueError(f"{bad} of {len(P)} rows do not sum to 1 within {_TOL}")

        idx = {c: i for i, c in enumerate(cols)}
        keep_ix = np.array([idx[c] for c in self.kept])
        drop_ix = np.array([idx[c] for c in self.excluded])

        kept_mass = P[:, keep_ix].sum(axis=1)
        impossible = P[:, drop_ix].sum(axis=1)

        alive = kept_mass > 0
        R = np.full((len(P), len(self.kept)), np.nan)
        R[alive] = P[np.ix_(alive, keep_ix)] / kept_mass[alive, None]

        kept_arr = np.asarray(self.kept)
        top = np.full(len(P), AMBIGUOUS, dtype=object)
        top[alive] = kept_arr[np.nanargmax(R[alive], axis=1)]

        reason = np.full(len(P), "", dtype=object)
        reason[~alive] = "no_mass_in_label_space"
        # A cell whose removed mass outweighs its surviving best is not one the model
        # placed in this label space; calling it anyway would invent the answer.  Both
        # sides are read off the ORIGINAL softmax, so this compares like with like and
        # introduces no cut of its own.
        outweighed = alive & (impossible > P[:, keep_ix].max(axis=1))
        top[outweighed] = AMBIGUOUS
        reason[outweighed] = "mass_on_impossible_classes"

        full = np.asarray(cols)
        return Restricted(
            probabilities=R,
            kept=self.kept,
            restricted_call=top.astype(str),
            unrestricted_call=full[P.argmax(axis=1)],
            impossible_mass=impossible,
            ambiguous_reason=reason.astype(str),
        )


def load_label_space(classes, path: pathlib.Path = CROSSWALK) -> LabelSpace:
    """Build the label space for `classes`, which must be the model's own column order."""
    d = json.loads(pathlib.Path(path).read_text())
    E = d["entries"]
    cols = tuple(classes)
    if set(cols) != set(E):
        raise ValueError(
            f"crosswalk covers {len(E)} classes, caller passed {len(cols)}; "
            f"symmetric difference {sorted(set(cols) ^ set(E))}"
        )
    for c, v in E.items():
        if "in_label_space" not in v:
            raise ValueError(f"{c} has no in_label_space decision -- an undecided class is not a pass")
        if not v["in_label_space"] and not (v.get("reason") or v.get("panel_reason")):
            raise ValueError(f"{c} is excluded with no written reason")
    kept = tuple(c for c in cols if E[c]["in_label_space"])
    excluded = tuple(c for c in cols if not E[c]["in_label_space"])
    if not kept:
        raise ValueError("the label space is empty -- an empty vocabulary is not a configuration")
    return LabelSpace(
        classes=cols,
        kept=kept,
        excluded=excluded,
        exclusion_reason={c: E[c]["exclusion_reason"] for c in excluded},
        broad_class={c: E[c]["broad_class"] for c in cols},
        method_version=d["method_version"],
    )
