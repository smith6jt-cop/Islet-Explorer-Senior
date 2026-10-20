"""Verdict helpers that fail closed.

Every gate in this pilot answers "did the check pass?".  The dangerous default is
a comparison that is *False* for missing data and therefore reports PASS: on the
first run of the plane-identity check every correlation was ``nan`` (the run's
expression table records a missing measurement as NaN, not zero -- see CLAUDE.md
"a missing measurement is unavailable, not zero or negative"), and because
``nan <= nan`` is False the gate reported PASS on 52 markers it had not actually
compared.

The rule these helpers enforce: **absent, non-finite or under-powered evidence is
a FAILURE, never a pass.**  A gate may only pass on evidence it actually has.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

__all__ = ["Verdict", "assess_self_match", "assess_anchor_concordance", "require_evidence"]


@dataclass(frozen=True)
class Verdict:
    passed: bool
    rows: list           # assess_self_match:        (name, self, best_other, who, n, ok, why)
                         # assess_anchor_concordance: (..., why, group, rank) -- 9 fields
    reason: str

    def __bool__(self) -> bool:            # so `if verdict:` cannot silently be truthy
        return self.passed


def require_evidence(n_obs: int, min_obs: int, what: str) -> str | None:
    """Return a failure reason when there is too little evidence to judge."""
    if n_obs < min_obs:
        return f"{what}: only {n_obs} usable observations, need >= {min_obs}"
    return None


def assess_self_match(labels, scores, counts, *, min_obs: int = 50,
                      min_self: float = 0.0, margin: float = 0.0) -> Verdict:
    """Each row of ``scores`` must score highest against its own column.

    ``scores[i, j]`` is how well measurement *i* matches reference *j*.
    ``counts[i]`` is how many observations backed the diagonal comparison.

    Fails when any of these is true -- none of them may be treated as a pass:
      * no labels at all;
      * the self score or the best rival score is not finite;
      * fewer than ``min_obs`` observations backed the comparison;
      * the self score does not beat every rival by more than ``margin``;
      * the self score is below ``min_self``.
    """
    S = np.asarray(scores, dtype=float)
    counts = np.asarray(counts, dtype=int)
    if S.ndim != 2 or S.shape[0] != S.shape[1] or S.shape[0] != len(labels):
        return Verdict(False, [], f"score matrix {S.shape} does not match {len(labels)} labels")
    if len(labels) == 0:
        return Verdict(False, [], "no markers to assess -- an empty check is not a passing check")

    diag = np.diag(S)
    rival = S.copy()
    np.fill_diagonal(rival, -np.inf)
    best = rival.max(axis=1)
    who = [labels[k] for k in rival.argmax(axis=1)]

    rows, failures = [], []
    for i, name in enumerate(labels):
        why = require_evidence(int(counts[i]), min_obs, str(name))
        if why is None:
            if not np.isfinite(diag[i]):
                why = "self score is not finite"
            elif not np.isfinite(best[i]):
                why = "best rival score is not finite"
            elif diag[i] < min_self:
                why = f"self score {diag[i]:.3f} below floor {min_self:.3f}"
            elif not (diag[i] > best[i] + margin):
                why = f"matches {who[i]} at least as well as itself"
        ok = why is None
        if not ok:
            failures.append(f"{name} ({why})")
        rows.append((name, diag[i], best[i], who[i], int(counts[i]), ok, why))

    return Verdict(not failures, rows,
                   "all markers match themselves" if not failures
                   else f"{len(failures)} of {len(labels)} failed: " + "; ".join(failures[:6])
                        + (" ..." if len(failures) > 6 else ""))


def assess_anchor_concordance(anchors, groups, scores, counts, *, min_obs: int = 50,
                              min_self: float = 0.0, margin: float = 0.0,
                              max_rank: int | None = None) -> Verdict:
    """Each anchor marker must agree best with the group it is a marker *for*.

    The rectangular sibling of :func:`assess_self_match`.  ``assess_self_match`` compares
    a thing to itself, so its matrix is square; here the rows are marker measurements and
    the columns are cell-type groups, which is a different shape and a weaker claim -- a
    marker is *evidence for* a compartment, not the same quantity as it.

    ``anchors[i]`` is ``(marker, group)``; ``scores[i, j]`` is how well marker *i* agrees
    with group *j*; ``groups[j]`` names the columns; ``counts[i]`` is how many observations
    backed anchor *i*.

    Fails -- never passes -- when any of these is true:
      * there are no anchors, or an anchor names a group that is not in ``groups``;
      * the intended score or the best rival score is not finite;
      * fewer than ``min_obs`` observations backed the anchor;
      * the intended score is below ``min_self`` (with 0.0 that is the sign test: a marker
        whose compartment probability moves *against* it is a failure, not a weak pass);
      * the intended score does not beat every rival by more than ``margin``, when
        ``max_rank`` is None;
      * the intended group ranks worse than ``max_rank``, when that is given.
    """
    S = np.asarray(scores, dtype=float)
    counts = np.asarray(counts, dtype=int)
    groups = list(groups)
    if not len(anchors):
        return Verdict(False, [], "no anchors to assess -- an empty check is not a passing check")
    if S.ndim != 2 or S.shape != (len(anchors), len(groups)):
        return Verdict(False, [], f"score matrix {S.shape} does not match "
                                  f"{len(anchors)} anchors x {len(groups)} groups")
    if len(counts) != len(anchors):
        return Verdict(False, [], f"{len(counts)} counts for {len(anchors)} anchors")

    rows, failures = [], []
    for i, (marker, group) in enumerate(anchors):
        if group not in groups:
            rows.append((marker, np.nan, np.nan, "-", int(counts[i]), False,
                         f"group {group!r} is not one of the {len(groups)} assessed",
                         group, len(groups)))
            failures.append(f"{marker} (unknown group {group!r})")
            continue
        j = groups.index(group)
        self_s = S[i, j]
        rival = np.delete(S[i], j)
        rival_names = [g for k, g in enumerate(groups) if k != j]
        best = rival.max() if rival.size else -np.inf
        who = rival_names[int(rival.argmax())] if rival.size else "-"
        # 1 = the intended group agrees best of all
        rank = 1 + int((rival > self_s).sum()) if np.isfinite(self_s) else len(groups)

        why = require_evidence(int(counts[i]), min_obs, str(marker))
        if why is None:
            if not np.isfinite(self_s):
                why = "score against its own group is not finite"
            elif rival.size and not np.isfinite(best):
                why = "best rival score is not finite"
            elif self_s < min_self:
                why = f"{self_s:.3f} against {group} is below floor {min_self:.3f}"
            elif max_rank is not None:
                if rank > max_rank:
                    why = f"{group} ranks {rank} of {len(groups)}; {who} agrees better"
            elif not (self_s > best + margin):
                why = f"agrees with {who} at least as well as with {group}"
        ok = why is None
        if not ok:
            failures.append(f"{marker} ({why})")
        rows.append((marker, self_s, best, who, int(counts[i]), ok, why, group, rank))

    return Verdict(not failures, rows,
                   f"all {len(anchors)} anchors agree best with their own compartment" if not failures
                   else f"{len(failures)} of {len(anchors)} failed: " + "; ".join(failures[:6])
                        + (" ..." if len(failures) > 6 else ""))
