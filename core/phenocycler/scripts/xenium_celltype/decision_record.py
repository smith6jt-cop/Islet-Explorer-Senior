"""Stage 6 — one decision record for three typing lanes, and one vocabulary for absence.

Three lanes in this repository emit the same KIND of statement — "this unit is (or is not)
this type, and here is what that rests on" — in three shapes that share no field:

    ``nullgate.GroupCall``                group, n_cells, lineage, reason, real, null_mean,
                                          null_p95, z, runners_up
    ``hierarchical_phenotype`` tier 1      obs: celltype_broad, broad_confidence,
                                          broad_n_identity_genes, broad_lowsignal,
                                          broad_marker_margin
    ``celltype/score.call_from_scores``    call, score, runner_up, margin, n_candidates, reason

A consumer that wants "what was called, on what evidence, and where was nothing called"
currently needs three readers and has to know three private conventions: that ``lineage == ""``
means no award, that ``call == "Unscored"`` means the cell had no vector to compare, that
``broad_n_identity_genes == 0`` means the cell expressed no identity gene of any lineage. This
module gives those one shape, and nothing else — the three lanes are untouched.

**Uniform shape must not manufacture uniform evidence.** The numbers are not interchangeable
and the record refuses to pretend they are:

* ``nullgate``'s ``real`` is a CONTRAST against a size- and detection-matched null — the only
  one of the three with a null at all, and the only one whose award means "beat something".
  Its ``z`` is documented in the lane as "for ranking only", so it is NOT copied into
  ``confidence``; a ranking device rendered as a confidence is a fabricated probability.
* ``score.py``'s ``score`` is a ceiling-normalised COSINE — a similarity, not a probability,
  so it also leaves ``confidence`` empty.
* ``hierarchical_phenotype``'s ``broad_confidence`` is a scANVI posterior, genuinely on
  [0, 1] — and it is the lane with NO null, so a high confidence there and a cleared null in
  ``nullgate`` are not comparable quantities however similar the columns look.

So ``statistic`` carries its own NAME, and ``null`` is ``None`` for the two lanes that compute
none. A reader comparing 0.9 to 0.9 across lanes has to pass the names first.

**The vocabulary and the NaN rule come from ``run_donor.composition``**, which is where they
were paid for. Rows there used to be built only from lineages some group won, so a lineage
that won nothing had no row — and a reader turns a missing row into a zero. That misreported
donor 6539's spatial arm as ``Immune_T`` 0.0% when Immune_T had in fact cleared its null and
been outranked in six separate groups. Three distinct states — awarded, cleared-but-outranked,
never-cleared — collapsed into one number naming none of them. Hence four statuses:

    ``awarded``                this unit carries this label, on the evidence recorded here
    ``unavailable``            nothing could be awarded for this unit; ``reason`` says what
                               stopped it, and ``label`` is None rather than the nearest type
    ``cleared_but_outranked``  this label was in contention here and lost to another
    ``no_group_cleared``       this label was offered here and never cleared

and one rule, which is the repository rule (``CLAUDE.md``: *a missing measurement is
unavailable — not zero, not negative*): **an unmeasured number is NaN and never 0.** The
record enforces it rather than documenting it — constructing a non-awarded record with
``confidence=0.0``, or a label-level record with a measured ``frac``, raises. That is the one
conversion this module exists to prevent, and :func:`assert_no_absence_as_zero` re-checks it
on a frame that has been through CSV or someone else's ``fillna``.

A measured zero stays a zero: ``Statistic.support`` of 0 ("zero identity genes were detected",
"zero references met their requirements") is a measurement and is written as 0, in the same
record whose ``value`` and ``confidence`` are NaN. Those two columns disagreeing is the point.

**This is an adapter layer, not a refactor.** Each adapter reads one lane's row and returns
one record; no lane imports this module and none of their shapes change.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd

NA = float("nan")

#: This unit carries this label, on the evidence in the record.
AWARDED = "awarded"
#: Nothing could be awarded for this unit. ``label`` is None; ``reason`` says what stopped it.
UNAVAILABLE = "unavailable"
#: This label was in contention for this unit and lost to another.
CLEARED_BUT_OUTRANKED = "cleared_but_outranked"
#: This label was offered for this unit and never cleared.
NO_GROUP_CLEARED = "no_group_cleared"

STATUSES = (AWARDED, UNAVAILABLE, CLEARED_BUT_OUTRANKED, NO_GROUP_CLEARED)

#: The two statuses that speak ABOUT a label rather than about a call. They exist because a
#: label that won nothing must still get a record, or its absence reads as a zero. Nothing was
#: measured for the label here, so the statistic and the confidence must both be NaN.
LABEL_LEVEL = (CLEARED_BUT_OUTRANKED, NO_GROUP_CLEARED)

#: ``score.call_from_scores`` sentinels. A cell with no signal has a zero norm and no angle;
#: a cell meeting no reference's requirements is not given its least-bad match.
UNSCORED = "Unscored"
AMBIGUOUS = "Ambiguous"

#: Reasons this module supplies where a lane encodes its refusal in the call string instead of
#: the reason column. Kept as constants so a consumer can match on them.
REASON_NOTHING_TO_COMPARE = "nothing_to_compare"
REASON_NO_REQUIREMENTS_MET = "no_reference_requirements_met"
REASON_NO_IDENTITY_GENE = "no_identity_gene_detected"
REASON_NO_BROAD_CALL = "no_broad_call"

#: ``nullgate`` matches its null on BOTH detection decile and magnitude decile. Named in the
#: record so "cleared its null" cannot be read as "cleared some null".
MATCHED_NULL = "size_and_detection_matched"

COLUMNS = ("unit", "n", "label", "reason", "statistic", "statistic_value",
           "statistic_margin", "statistic_lo", "statistic_hi", "statistic_support",
           "null_kind", "null_mean", "null_p95", "null_z", "confidence",
           "runners_up", "status")


class DecisionRecordError(ValueError):
    """A record that would render an absence as a number, or name a status that is not one."""


# ---------------------------------------------------------------------------
# coercion — every one of these maps "not measured" to NaN/None and never to 0
# ---------------------------------------------------------------------------
def _missing(v: Any) -> bool:
    if v is None:
        return True
    try:
        r = pd.isna(v)
    except (TypeError, ValueError):
        return False
    return bool(r) if isinstance(r, (bool, np.bool_)) else False


def _as_float(v: Any) -> float:
    """Absent -> NaN. Never 0: that is the conversion the module exists to refuse."""
    if _missing(v) or (isinstance(v, str) and not v.strip()):
        return NA
    return float(v)


def _as_label(v: Any) -> str | None:
    """Absent, empty, or NaN -> None. ``nullgate`` writes ``""`` for "nothing was awarded"."""
    if _missing(v):
        return None
    s = str(v).strip()
    return s or None


def _as_reason(v: Any) -> str:
    return "" if _missing(v) else str(v).strip()


def _as_support(v: Any) -> int | None:
    """None means the lane reports no support count — distinct from a measured 0."""
    return None if _missing(v) else int(v)


def _as_runners(v: Any) -> tuple[str, ...]:
    if _missing(v):
        return ()
    if isinstance(v, str):
        return tuple(x for x in v.split("|") if x)
    if isinstance(v, (list, tuple, set, np.ndarray, pd.Index, pd.Series)):
        return tuple(str(x) for x in v if not _missing(x) and str(x))
    return (str(v),)


def _get(row: Any, key: str, default: Any = None) -> Any:
    """One reader for a dataclass, a Mapping, a Series or an itertuples row."""
    if isinstance(row, (Mapping, pd.Series)):
        return row[key] if key in row else default
    return getattr(row, key, default)


# ---------------------------------------------------------------------------
# the schema
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Statistic:
    """The number the call rests on, carrying WHAT was measured.

    ``name`` is not decoration. A contrast against a matched null, a ceiling-normalised
    cosine and a scANVI posterior all land in ``value``, and only the name stops them being
    compared to each other.

    ``margin`` is the distance to the runner-up ON THE SAME SCALE as ``value``; ``lo``/``hi``
    are an interval where the lane computed one (``run_donor.composition`` is the only lane
    that does). ``support`` is how many measurements back the value — identity-core genes
    detected, references that met their requirements. ``support`` is the one field here where
    0 is a legitimate answer, and it is what lets "measured, and the answer was none" stand
    beside a NaN ``value`` in the same record.
    """
    name: str = ""
    value: float = NA
    margin: float = NA
    lo: float = NA
    hi: float = NA
    support: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", "" if _missing(self.name) else str(self.name))
        for f in ("value", "margin", "lo", "hi"):
            object.__setattr__(self, f, _as_float(getattr(self, f)))
        object.__setattr__(self, "support", _as_support(self.support))

    @property
    def measured(self) -> bool:
        return bool(np.isfinite(self.value))


@dataclass(frozen=True)
class Null:
    """What the statistic was compared against, where anything was.

    Only ``nullgate`` builds one. ``None`` in :attr:`DecisionRecord.null` means the lane
    computes no null at all; a ``Null`` whose fields are NaN means a null belonged here and
    was not computed (a group too small to pool, for instance). Those are different states and
    the record keeps them apart.
    """
    kind: str = MATCHED_NULL
    mean: float = NA
    p95: float = NA
    z: float = NA

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", "" if _missing(self.kind) else str(self.kind))
        for f in ("mean", "p95", "z"):
            object.__setattr__(self, f, _as_float(getattr(self, f)))


@dataclass(frozen=True)
class DecisionRecord:
    """One typing statement, whichever lane made it.

    ``unit`` is what the statement is about — a group id, a cell id, a (donor, method, k)
    tuple. The caller owns its uniqueness within a frame; the adapters default it to the
    lane's own identifier, which is only unique within one run of one lane.

    ``n`` is how many cells the statement pools over: the group size for ``nullgate``, 1 for a
    per-cell call, the cells carrying a share for a composition row. It is a count and 0 is a
    real answer for a label that won no cells.

    ``label`` is ``None`` whenever nothing was awarded — never ``""``, never the nearest type.
    For the two :data:`LABEL_LEVEL` statuses it names the label the statement is ABOUT, which
    is how a label that won nothing keeps a row.

    ``confidence`` is populated only where a lane produced a probability. Two of the three do
    not, and NaN there is a statement about the lane, not about the call.
    """
    unit: Any
    n: int = 0
    label: str | None = None
    reason: str = ""
    statistic: Statistic = field(default_factory=Statistic)
    null: Null | None = None
    confidence: float = NA
    runners_up: tuple[str, ...] = ()
    status: str = AWARDED

    def __post_init__(self) -> None:
        object.__setattr__(self, "n", 0 if _missing(self.n) else int(self.n))
        object.__setattr__(self, "label", _as_label(self.label))
        object.__setattr__(self, "reason", _as_reason(self.reason))
        object.__setattr__(self, "confidence", _as_float(self.confidence))
        object.__setattr__(self, "runners_up", _as_runners(self.runners_up))
        if not isinstance(self.statistic, Statistic):
            object.__setattr__(self, "statistic", Statistic(**dict(self.statistic or {})))
        if self.null is not None and not isinstance(self.null, Null):
            object.__setattr__(self, "null", Null(**dict(self.null)))

        if self.status not in STATUSES:
            raise DecisionRecordError(
                f"status {self.status!r} is not one of {STATUSES}; a lane's private word for "
                "a refusal belongs in `reason`, which is carried verbatim"
            )
        if self.n < 0:
            raise DecisionRecordError(f"n must be a count, got {self.n}")

        if self.status == AWARDED:
            if self.label is None:
                raise DecisionRecordError(
                    f"unit {self.unit!r} is awarded with no label; an award without a label is "
                    f"the absence this record exists to keep visible — use {UNAVAILABLE!r}"
                )
            if np.isfinite(self.confidence) and not 0.0 <= self.confidence <= 1.0:
                raise DecisionRecordError(
                    f"confidence {self.confidence} is outside [0, 1]; only a probability "
                    "belongs in this field — a cosine or a z goes in `statistic`"
                )
        else:
            if np.isfinite(self.confidence):
                raise DecisionRecordError(
                    f"unit {self.unit!r} has status {self.status!r} and confidence "
                    f"{self.confidence!r}. Nothing was awarded, so there is no call to be "
                    "confident in; the value must be NaN. A 0 here is an absence rendered as "
                    "a number, which is the one thing this module refuses."
                )
            if not self.reason:
                raise DecisionRecordError(
                    f"unit {self.unit!r} has status {self.status!r} with no reason; a refusal "
                    "keeps its reason, or a reader cannot tell pooling from evidence"
                )

        if self.status == UNAVAILABLE and self.label is not None:
            raise DecisionRecordError(
                f"unit {self.unit!r} is unavailable but labelled {self.label!r}; the lane's "
                "unawarded pick belongs in `runners_up`, where it reads as a contender "
                "rather than a call"
            )
        if self.status in LABEL_LEVEL:
            if self.label is None:
                raise DecisionRecordError(
                    f"status {self.status!r} is a statement about a label and needs one"
                )
            if self.statistic.measured:
                raise DecisionRecordError(
                    f"label {self.label!r} has status {self.status!r} and a measured "
                    f"{self.statistic.name or 'statistic'} of {self.statistic.value}. It won "
                    "nothing here, so its share is unmeasured and must be NaN — a number "
                    "there is exactly how Immune_T was reported as 0.0%."
                )

    # -- flat form -------------------------------------------------------
    def to_row(self) -> dict:
        """A flat dict for a DataFrame or a CSV. Unmeasured fields stay NaN."""
        s, z = self.statistic, self.null
        return {
            "unit": self.unit,
            "n": self.n,
            "label": self.label,
            "reason": self.reason,
            "statistic": s.name,
            "statistic_value": s.value,
            "statistic_margin": s.margin,
            "statistic_lo": s.lo,
            "statistic_hi": s.hi,
            "statistic_support": s.support,
            "null_kind": "" if z is None else z.kind,
            "null_mean": NA if z is None else z.mean,
            "null_p95": NA if z is None else z.p95,
            "null_z": NA if z is None else z.z,
            "confidence": self.confidence,
            "runners_up": "|".join(self.runners_up),
            "status": self.status,
        }

    @classmethod
    def from_row(cls, row: Any) -> "DecisionRecord":
        """Rebuild from :meth:`to_row`, including after a CSV round trip.

        An empty ``null_kind`` restores ``null=None`` — the lane computed no null — while a
        named kind with NaN fields restores a null that belonged and was not computed.
        """
        kind = _as_reason(_get(row, "null_kind", ""))
        return cls(
            unit=_get(row, "unit"),
            n=_get(row, "n", 0),
            label=_get(row, "label"),
            reason=_get(row, "reason", ""),
            statistic=Statistic(
                name=_get(row, "statistic", ""),
                value=_get(row, "statistic_value"),
                margin=_get(row, "statistic_margin"),
                lo=_get(row, "statistic_lo"),
                hi=_get(row, "statistic_hi"),
                support=_get(row, "statistic_support"),
            ),
            null=None if not kind else Null(
                kind=kind,
                mean=_get(row, "null_mean"),
                p95=_get(row, "null_p95"),
                z=_get(row, "null_z"),
            ),
            confidence=_get(row, "confidence"),
            runners_up=_get(row, "runners_up", ()),
            status=_as_reason(_get(row, "status", AWARDED)),
        )


def records_to_frame(records: Iterable[DecisionRecord]) -> pd.DataFrame:
    """Records -> one frame with the fixed column order, empty frame included.

    Nothing here fills, and nothing downstream should: ``frame.fillna(0)`` on this table
    converts every absence into a measured zero, which :func:`assert_no_absence_as_zero`
    exists to catch.
    """
    rows = [r.to_row() for r in records]
    return pd.DataFrame(rows, columns=list(COLUMNS)) if rows else pd.DataFrame(
        {c: pd.Series(dtype="object") for c in COLUMNS})


def frame_to_records(frame: pd.DataFrame) -> list[DecisionRecord]:
    return [DecisionRecord.from_row(r) for _, r in frame.iterrows()]


def assert_no_absence_as_zero(frame: pd.DataFrame) -> None:
    """Re-check the rule on a frame that has been somewhere — a CSV, a merge, a ``fillna``.

    The constructor enforces this, so a frame built here is already clean; a frame read back
    from disk or passed through pandas is not, and this is the cheapest place to find out.
    """
    missing = [c for c in ("status", "confidence", "statistic_value") if c not in frame]
    if missing:
        raise DecisionRecordError(f"not a decision-record frame; missing {missing}")
    bad = []
    st = frame["status"].astype(str)
    conf = pd.to_numeric(frame["confidence"], errors="coerce")
    hit = st.ne(AWARDED) & conf.notna()
    if hit.any():
        bad.append(f"confidence is a number on {int(hit.sum())} non-awarded row(s): "
                   f"{sorted(frame.loc[hit, 'status'].unique())}")
    val = pd.to_numeric(frame["statistic_value"], errors="coerce")
    hit2 = st.isin(LABEL_LEVEL) & val.notna()
    if hit2.any():
        bad.append(f"a statistic is measured on {int(hit2.sum())} row(s) whose label won "
                   "nothing here")
    if bad:
        raise DecisionRecordError(
            "an absence has been rendered as a number: " + "; ".join(bad)
            + ". A missing measurement is unavailable — not zero, not negative."
        )


# ---------------------------------------------------------------------------
# adapter 1 — nullgate.GroupCall / annotate_groups rows
# ---------------------------------------------------------------------------
def from_group_call(call: Any, *, unit: Any = None) -> DecisionRecord:
    """A ``GroupCall`` (or a row of ``annotate_groups``) -> one record.

    ``confidence`` stays NaN on purpose. The lane's ``z`` is annotated "for ranking only" and
    is a standardised distance from a 50-draw null, not a probability of anything; copying it
    into ``confidence`` would put it alongside a scANVI posterior as though the two could be
    read the same way. It is carried in ``null.z``, where it says what it is.

    ``runners_up`` changes meaning with the reason, exactly as the lane writes it: on an award
    it is the OTHER cores that also cleared, and on ``cores_cleared_but_are_not_separable`` it
    is every core that cleared, the top one included — because none of them won.
    """
    lineage = _as_label(_get(call, "lineage"))
    group = _get(call, "group")
    return DecisionRecord(
        unit=group if unit is None else unit,
        n=_get(call, "n_cells", 0),
        label=lineage,
        reason=_as_reason(_get(call, "reason")),
        statistic=Statistic(name="contrast", value=_get(call, "real")),
        null=Null(kind=MATCHED_NULL, mean=_get(call, "null_mean"),
                  p95=_get(call, "null_p95"), z=_get(call, "z")),
        confidence=NA,
        runners_up=_get(call, "runners_up"),
        status=AWARDED if lineage is not None else UNAVAILABLE,
    )


def from_group_calls(calls: pd.DataFrame, *, unit: Any = None) -> Iterator[DecisionRecord]:
    """``annotate_groups`` output -> records. ``unit`` may be a callable taking the row."""
    for _, row in calls.iterrows():
        yield from_group_call(row, unit=unit(row) if callable(unit) else unit)


# ---------------------------------------------------------------------------
# adapter 2 — celltype/score.call_from_scores rows
# ---------------------------------------------------------------------------
def from_cell_call(row: Any, *, unit: Any) -> DecisionRecord:
    """A row of ``call_from_scores`` -> one record.

    Two refusals arrive as call strings rather than in ``reason``, and both become
    ``unavailable`` with the reason spelled out: ``Unscored`` is a cell whose vector had zero
    norm, so there was no angle to take — "nothing to compare", which the lane already keeps
    apart from "compared and unlike" by returning NaN instead of 0. ``Ambiguous`` is a cell
    that met no reference's requirements and was not handed its least-bad match.

    The cosine goes in ``statistic``, not ``confidence``: it is ceiling-normalised so 1 means
    "as close to this reference as any non-negative cell could be", which is a similarity, not
    a probability. ``n_candidates`` becomes ``statistic.support`` — and its measured 0 on an
    ``Ambiguous`` cell sits beside a NaN score in the same record.
    """
    call = _as_label(_get(row, "call"))
    reason = _as_reason(_get(row, "reason"))
    support = _get(row, "n_candidates")
    if call == UNSCORED:
        label, status = None, UNAVAILABLE
        reason = reason or REASON_NOTHING_TO_COMPARE
    elif call == AMBIGUOUS:
        label, status = None, UNAVAILABLE
        reason = reason or REASON_NO_REQUIREMENTS_MET
    elif call is None:
        label, status = None, UNAVAILABLE
        reason = reason or REASON_NOTHING_TO_COMPARE
    else:
        label, status = call, AWARDED
    return DecisionRecord(
        unit=unit,
        n=1,
        label=label,
        reason=reason,
        statistic=Statistic(name="cosine", value=_get(row, "score"),
                            margin=_get(row, "margin"), support=support),
        null=None,                      # this lane compares against no null
        confidence=NA,                  # a cosine is a similarity, not a probability
        runners_up=_get(row, "runner_up"),
        status=status,
    )


def from_cell_calls(calls: pd.DataFrame, *, units: Sequence | None = None
                    ) -> Iterator[DecisionRecord]:
    """``call_from_scores`` output -> records, one per cell.

    ``call_from_scores`` returns positional rows, so pass the ``object_id`` column the caller
    holds as ``units``; the frame's index is only a fallback.
    """
    idx = list(calls.index) if units is None else list(units)
    if len(idx) != len(calls):
        raise DecisionRecordError(
            f"{len(idx)} units for {len(calls)} calls; the identifier must not be positional")
    for u, (_, row) in zip(idx, calls.iterrows()):
        yield from_cell_call(row, unit=u)


# ---------------------------------------------------------------------------
# adapter 3 — hierarchical_phenotype tier-1 obs columns
# ---------------------------------------------------------------------------
_BROAD_COLS = ("celltype_broad", "broad_confidence", "broad_n_identity_genes",
               "broad_lowsignal", "broad_marker_margin")


def from_broad_call(row: Any, *, unit: Any) -> DecisionRecord:
    """One cell's tier-1 obs fields -> one record.

    This lane has no unavailable state of its own: ``celltype_broad`` is Indeterminate-free by
    contract, and a cell expressing no identity-core gene of ANY lineage
    (``broad_lowsignal``) still receives the absolute-score argmax. So it is the lane that
    renders an absence as a call, and the adapter is where that becomes visible: a cell whose
    ``broad_n_identity_genes`` is 0 becomes ``unavailable`` with the lane's pick moved into
    ``runners_up``, ``support`` recorded as the measured 0 it is, and the posterior dropped —
    a posterior over a cell with no identity evidence is a number about the model, not about
    the cell. The lane's own label is unchanged on disk; nothing here writes back.

    ``broad_marker_margin`` is carried as the lane computes it — the gap between the top two
    posteriors BEFORE the expression mask. For the cells where the mask overrode the free
    argmax it is therefore not the margin of the decision that was taken, and it is reported
    rather than gated for that reason.
    """
    label = _as_label(_get(row, "celltype_broad"))
    n_genes = _as_support(_get(row, "broad_n_identity_genes"))
    low = _get(row, "broad_lowsignal")
    lowsignal = bool(low) if not _missing(low) else (n_genes == 0)

    if label is None:
        return DecisionRecord(unit=unit, n=1, label=None, reason=REASON_NO_BROAD_CALL,
                              statistic=Statistic(name="scanvi_posterior", support=n_genes),
                              null=None, confidence=NA, status=UNAVAILABLE)
    if n_genes == 0 or (n_genes is None and lowsignal):
        return DecisionRecord(
            unit=unit, n=1, label=None, reason=REASON_NO_IDENTITY_GENE,
            statistic=Statistic(name="scanvi_posterior", support=0),
            null=None, confidence=NA, runners_up=(label,), status=UNAVAILABLE)
    conf = _as_float(_get(row, "broad_confidence"))
    return DecisionRecord(
        unit=unit, n=1, label=label, reason="",
        statistic=Statistic(name="scanvi_posterior", value=conf,
                            margin=_get(row, "broad_marker_margin"), support=n_genes),
        null=None,                      # this lane compares against no null
        confidence=conf,                # the one lane of the three with a probability
        status=AWARDED,
    )


def from_broad_obs(obs: pd.DataFrame, *, units: Sequence | None = None
                   ) -> Iterator[DecisionRecord]:
    """``adata.obs`` -> records, one per cell.

    A generator, not a list: a donor is 1.2M cells and a per-cell dataclass for each is an
    inspection or export tool, not a table to materialise.
    """
    cols = {c: obs[c].to_numpy() for c in _BROAD_COLS if c in obs.columns}
    if "celltype_broad" not in cols:
        raise DecisionRecordError("obs carries no `celltype_broad`; this is not tier-1 output")
    idx = list(obs.index) if units is None else list(units)
    for i, u in enumerate(idx):
        yield from_broad_call({c: v[i] for c, v in cols.items()}, unit=u)


# ---------------------------------------------------------------------------
# adapter 4 — run_donor.composition rows (where the vocabulary comes from)
# ---------------------------------------------------------------------------
def from_composition_row(row: Any, *, unit: Any) -> DecisionRecord:
    """A ``composition`` row -> one record, keeping all four statuses.

    ``composition`` is the lane the vocabulary is taken from, so it is adapted too: it is the
    only shape that emits ``cleared_but_outranked`` and ``no_group_cleared``, and without it
    those two words would be in the schema and in no output.

    An unavailable row's lineage is ``Unavailable:<reason>`` — the reason is kept, because
    ``too_few_cells`` is a statement about pooling and ``no_core_beats_its_matched_null`` is a
    statement about evidence, and folding them together loses which one happened. Its share IS
    measured (those cells were counted), so ``statistic.value`` is finite while ``confidence``
    stays NaN; a label-level row's share is not measured and stays NaN, which is what the
    constructor checks.
    """
    name = str(_get(row, "lineage", ""))
    status = _as_reason(_get(row, "status", AWARDED))
    label: str | None = _as_label(name)
    reason = ""
    if status == UNAVAILABLE or name.startswith("Unavailable:"):
        status = UNAVAILABLE
        reason = name.split(":", 1)[1] if ":" in name else "unavailable"
        label = None
    elif status in LABEL_LEVEL:
        reason = status
    return DecisionRecord(
        unit=unit,
        n=_get(row, "n_cells", 0),
        label=label,
        reason=reason,
        statistic=Statistic(name="frac", value=_get(row, "frac"),
                            lo=_get(row, "ci_lo"), hi=_get(row, "ci_hi")),
        null=None,                      # the null was applied per group, not to the share
        confidence=NA,
        status=status,
    )


def from_composition(comp: pd.DataFrame, *, unit: Any) -> Iterator[DecisionRecord]:
    for _, row in comp.iterrows():
        yield from_composition_row(row, unit=unit)
