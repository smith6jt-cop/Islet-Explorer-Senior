#!/usr/bin/env python3
"""One rule for partial runs, shared by the pLN lanes: a partial run never writes a
cohort-wide path.

WHAT BIT US, TWICE IN ONE DAY (2026-09-17)
------------------------------------------
`build_pln_immune_pheno.py --donors 6539 --stage standout` rewrote three cohort-wide files
with the one donor it had just processed: `cluster_calls.csv` 293 rows -> 19,
`cluster_heterogeneity.csv` 293 -> 19, `per_donor_separation.csv` 14 -> 1.  Earlier the same
day `protein_rna --donors 6539` dropped 6436 and 6442 from its summary index the same way.
Nothing was corrupted and nothing errored -- the file simply became a smaller file with the
same name, and the next reader had no way to tell.  `cluster_pln.py --stage profile --donors
115` had the identical shape: it concatenates the per-donor naming profiles it just wrote and
lands them on the cohort-wide `naming_profiles.parquet`.

THE RULE, AND WHY THIS ONE
--------------------------
    A run whose unit set is a strict subset of the units the lane could process is PARTIAL,
    and a partial run writes its cohort-shaped outputs under `_partial/<kind>s=<slug>/`.
    The cohort-wide path is left exactly as it was.

Three candidates were on the table.

*Merge* -- re-read every per-unit shard on disk and concatenate -- was rejected.  The shards on
disk are not all from the same generation: a shard written by yesterday's parameterisation
merges just as silently as today's, and this repository has already been bitten by a typing
tree that stacks generations side by side with nothing marking which is which.  A merge would
turn a visible truncation into an invisible mixture, which is worse.

*Refuse* is right where the output is a published artifact and a scoped copy of it would be
meaningless -- there is no such thing as a one-donor node-membership LAYER.  `build_pln_layer.py`
uses `refuse_shrink` for exactly that case.

*Scope* is right where the partial run is a legitimate thing to do and only its claim to cohort
coverage is wrong.  A one-donor naming profile is a useful object; it is just not
`naming_profiles.parquet`.  So it is written, in full, one directory down, and the reader can
see from the path what it covers.

Both are the same rule -- the cohort path belongs to the cohort run -- so `RunScope` carries it
and each lane picks the consequence its artifact deserves.

WHAT COUNTS AS PARTIAL
----------------------
Partial is measured against what the lane could have processed on this invocation, not against
a remembered cohort size: `requested` narrower than `available`.  A run that names every
available unit explicitly is NOT partial, because its output does cover them all.  A run that
names a unit which is not available at all is an error the lane raises, not a partial run.

`--out-root` does not make a run whole: a partial run into a private root still writes to the
scoped path there, so that copying such a tree onto the tracked root cannot smuggle a partial
file onto a cohort name.
"""
from __future__ import annotations

import dataclasses
import hashlib
import pathlib
import re

#: Everything a partial run writes lands under this directory inside the lane's output root.
PARTIAL_DIR = "_partial"

#: Longer than this and the slug is a hash, so a 20-donor subset cannot make a path no
#: filesystem will take.  The full member list is always written beside it in SCOPE.json.
_SLUG_MAX = 64


def scope_slug(units) -> str:
    """A deterministic, filesystem-safe name for a set of units.

    Sorted, so the same subset named in a different order is the same directory rather than a
    second one holding the same thing under another name.
    """
    members = sorted({str(u) for u in units})
    if not members:
        return "none"
    clean = [re.sub(r"[^A-Za-z0-9._-]", "_", m) for m in members]
    joined = "+".join(clean)
    if len(joined) <= _SLUG_MAX:
        return joined
    digest = hashlib.sha256("\n".join(members).encode()).hexdigest()[:12]
    return f"{len(members)}units-{digest}"


@dataclasses.dataclass(frozen=True)
class RunScope:
    """What this run covers, and therefore where its cohort-shaped outputs may be written."""

    unit_kind: str                      # "donor", "section", ...
    requested: tuple[str, ...]          # what the run processes
    available: tuple[str, ...]          # what it could have processed

    @classmethod
    def build(cls, unit_kind: str, requested, available) -> "RunScope":
        av = tuple(sorted({str(u) for u in available}))
        req = tuple(sorted({str(u) for u in (requested if requested else av)}))
        return cls(unit_kind=unit_kind, requested=req, available=av)

    @classmethod
    def undiscovered(cls, unit_kind: str, requested) -> "RunScope":
        """A selection made before the lane has listed what it holds.

        A stage that never enumerates its units -- `summarize` reads whatever assignments are
        on disk -- cannot prove that a named selection covers them all, so a named selection is
        treated as partial until a stage that does enumerate calls `refine`.  Conservative in
        the one direction that matters: the cohort file is never written by a run that cannot
        show it covers the cohort.
        """
        return cls(unit_kind=unit_kind,
                   requested=tuple(sorted({str(u) for u in (requested or ())})),
                   available=())

    def refine(self, available) -> "RunScope":
        """The same selection, once the lane has enumerated what it could have processed."""
        return RunScope.build(self.unit_kind, self.requested, available)

    @property
    def is_partial(self) -> bool:
        if not self.available:
            return bool(self.requested)
        return set(self.requested) < set(self.available)

    @property
    def slug(self) -> str:
        return scope_slug(self.requested)

    @property
    def missing(self) -> tuple[str, ...]:
        """Units the run was asked for that the lane does not have. The lane decides what to do."""
        return tuple(sorted(set(self.requested) - set(self.available)))

    def partial_dir(self, out_dir) -> pathlib.Path:
        return pathlib.Path(out_dir) / PARTIAL_DIR / f"{self.unit_kind}s={self.slug}"

    def cohort_path(self, out_dir, name: str) -> pathlib.Path:
        """Where a COHORT-WIDE artifact called `name` belongs for THIS run.

        A whole run gets the cohort path.  A partial run gets the scoped path, and the
        cohort-wide file of the same name is not touched, not read and not deleted.
        """
        out_dir = pathlib.Path(out_dir)
        if not self.is_partial:
            return out_dir / name
        dest = self.partial_dir(out_dir) / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        return dest

    def describe(self) -> dict:
        return {
            "unit_kind": self.unit_kind,
            "n_requested": len(self.requested),
            "n_available": len(self.available) if self.available else None,
            "requested": list(self.requested),
            "available": list(self.available) if self.available else "not_enumerated_by_this_stage",
            "is_partial": self.is_partial,
            "slug": self.slug if self.is_partial else None,
            "rule": ("a partial run writes cohort-shaped outputs under "
                     f"{PARTIAL_DIR}/<kind>s=<slug>/ and leaves the cohort-wide file alone"),
        }

    def announce(self, out_dir, log=print) -> None:
        """Say, once, where the cohort-shaped outputs of this run are going."""
        if not self.is_partial:
            return
        kept = sorted(set(self.available) - set(self.requested))
        of = f"of {len(self.available)} " if self.available else "of an unenumerated cohort; "
        tail = (f"({len(kept)} {self.unit_kind}s not in this run)" if self.available
                else "(this stage does not enumerate its units)")
        log(f"[scope] PARTIAL run: {len(self.requested)} {of}"
            f"{self.unit_kind}s ({', '.join(self.requested)}). Cohort-wide outputs go to "
            f"{self.partial_dir(out_dir)} -- the cohort files are left as they are {tail}.")


def write_scope_note(scope: RunScope, out_dir, extra: dict | None = None) -> pathlib.Path | None:
    """Drop a SCOPE.json beside a partial run's outputs saying what it does and does not cover."""
    import json

    if not scope.is_partial:
        return None
    dest = scope.partial_dir(out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    note = scope.describe()
    note["not_in_this_run"] = (sorted(set(scope.available) - set(scope.requested))
                               if scope.available else "not_enumerated_by_this_stage")
    if extra:
        note.update(extra)
    path = dest / "SCOPE.json"
    path.write_text(json.dumps(note, indent=2) + "\n")
    return path


def refuse_shrink(*, artifact: str, unit_kind: str, losing: dict, hint: str) -> None:
    """Stop a run that would take units OUT of a published artifact for want of an input.

    `losing` maps unit -> why this run does not write it.  Only pass the units whose absence is
    an ABSENT INPUT.  A unit dropped by a decision -- a cohort exclusion, a tracked manifest
    withholding it -- is not a partial run and must not be passed here: the artifact is supposed
    to lose it, and the reconciliation is supposed to say so.
    """
    if not losing:
        return
    lines = "\n".join(f"    {unit}: {reason}" for unit, reason in sorted(losing.items()))
    raise SystemExit(
        f"[err] refusing to write {artifact}: this run does not have an input for "
        f"{len(losing)} {unit_kind}(s) that the populated output root already holds, so "
        "publishing it would silently remove them:\n"
        f"{lines}\n"
        f"  {hint}\n"
        "  (A unit withheld by a DECISION -- phenocycler.cohort, or the tracked roster's "
        "excluded flag -- is removed as before; this refusal is only for a missing input.)"
    )
