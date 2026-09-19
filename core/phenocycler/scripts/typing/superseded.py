"""The one table of superseded artifacts, and the refusal every reader gets from it.

WHY THIS FILE EXISTS.  Being retired in prose is not the same as being unreachable.  Every
superseded artifact in this tree was already described as superseded -- in ``CLAUDE.md``, in
``scripts/retire/dispositions.py``, in a generated ``SUPERSEDED.md`` beside the data itself --
and every one of them was still a plain path that ``pq.read_table`` would open without a word.
On 2026-09-11 and again later, ``typing/subpanel/subpanel_calls.parquet`` was read as the
current endocrine typing.  Nothing stopped it, because nothing was asked.

So the knowledge lives in ONE place -- :data:`REGISTRY` below -- and the refusal is attached to
the paths themselves, through :class:`GuardedPath`, which ``subpanel_io`` hands out as ``OUT``
and ``RUN``.  Every module in this lane already builds its paths from those two, so every one of
them is guarded without calling anything.

    from superseded import refuse_if_superseded        # explicit, for a path from elsewhere
    S.OUT / "umap/refinement_log.json"                 # guarded for free: raises on open

TWO SEVERITIES, AND THE DIFFERENCE IS DELIBERATE.

``REFUSE`` raises.  It is for an artifact whose replacement is on disk right now: reading the old
one is a mistake with a named fix.

``ANNOUNCE`` prints the same notice to stderr, once per process, and lets the read proceed.  It is
for an artifact that is historical but is still the ONLY copy of something the current work needs
-- ``core/phenocycler/data/runs/<run_id>/`` is the whole of it.  The run framework was removed on
2026-09-16, but the 20-donor cohort expression output lives nowhere else, and ``subpanel_io``'s
own ``cell_source`` reads it for every donor that has not been re-exported.  Refusing it would
stop the typing lane, which is a method decision, and this module does not make method decisions.
It makes reads loud.  "Must not be SILENTLY readable" is the requirement, and an announcement
satisfies it where a refusal would break the build.

THE OTHER ROUTE TO THE SAME FILE.  A ``GuardedPath`` catches everything BUILT from ``S.OUT`` or
``S.RUN``.  It cannot catch a literal -- ``Path("/home/.../umap/refinement_log.json").read_text()``
or a path arriving from a config -- because such a path never passes through one.  For the one
artifact where that mattered, ``refine.py`` installs a process-wide ``open`` audit hook
(``refine.refuse_superseded_refinement_log``) that refuses the 2026-09-11 orphan log however it is
spelled.  The two are deliberately separate and cover different routes; this table is where the
artifact is RECORDED, and that hook is a second door on one of them.

WHAT IS NOT REGISTERED IS NOT A JUDGEMENT.  A path this table does not carry is simply a path
nobody has measured a superseding replacement for.  It is not thereby blessed.

ONE THING DELIBERATELY LEFT OUT, so it is not "fixed" later.  The superseded pLN typers --
``data/development/xenium_typing/release_pln_immune`` and ``_v2``, both 2026-09-17, sitting beside
the current ``_v3`` -- are NOT registered here.  They already carry their own supersession, in
each release's ``release.json`` under ``producer.supersedes``, and ``scripts/senior/xenium_release``
reads and honours that declaration; the artifact says which generation it is, which is stronger
than a path table asserting it from outside.  A path refusal would also break ``_v3``'s own
producer, which reads ``_v2`` to compare against it.  Registering them would duplicate a better
mechanism and take a live comparison down with it.

WHAT IS EXEMPT, AND WHY.  Two live inputs sit INSIDE retired trees.  They are named in
:data:`EXEMPT` rather than left to a lucky prefix miss:

  * ``data/development/broadtype/scratch/`` -- the withdrawn producer comparison's directory, but
    five live modules use it as their cache root and ``subpanel_io.load_scale()`` reads
    ``cosine_scale.parquet`` out of it as the cohort scale the classification was made on.
  * ``core/phenocycler/data/runs/<run_id>/expression`` is NOT exempt; it is ANNOUNCE, above.

THE ONE ALLOWED READER.  ``subpanel.py::gate`` reads ``subpanel_calls.parquet`` on purpose: it is
the row-for-row reconstruction gate, and ``CLAUDE.md`` keeps the four pooled reference levels
alive for exactly that ("survive only so ``subpanel.py --gate`` can reproduce it").  That caller
is recorded in the record itself, as ``allowed_readers``, and it STILL gets the full notice on
stderr.  Allowed is not silent.  No other caller may read that file, and nothing else in this
table has an allowed reader at all.

WHAT IT COSTS.  Measured 2026-09-19 on this machine: ``os.fspath`` on a guarded path that is NOT
registered takes 0.23 us more than on a plain one (0.022 s -> 0.068 s over 200,000 calls).  There
is no syscall on that route -- the string test is cached per distinct path and the symlink resolve
happens only in the explicit :func:`refuse_if_superseded` a builder calls once -- so a bulk read
pays it once per file open, against milliseconds of parquet I/O.

MEASURED 2026-09-19 in this checkout, with the numbers this module quotes:

    data/development/typing/subpanel/subpanel_calls.parquet  25,150,537 rows  2026-09-11 17:14
    data/development/typing/cells/ (20 donor partitions)     25,230,109 rows  2026-09-15 13:16
                                                             ---------
                                                                79,572 cells apart
"""
from __future__ import annotations

import functools
import inspect
import os
import pathlib
import sys
from dataclasses import dataclass, field

#: Repository root.  ``core/phenocycler/scripts/typing/superseded.py`` -> four parents up.
REPO = pathlib.Path(__file__).resolve().parents[4]

REFUSE = "refuse"
ANNOUNCE = "announce"


class SupersededArtifactError(ValueError):
    """Raised instead of reading an artifact a measured replacement has superseded.

    A ``ValueError`` because that is what the readers in this tree already raise and catch for a
    refused input; the subclass exists so a caller that wants to distinguish this refusal from a
    bad argument can.
    """


@dataclass(frozen=True)
class Record:
    """One superseded artifact: what it is, when it was superseded, by what, and by how much."""

    key: str
    #: Repo-relative path fragments that name this artifact.  A file record's fragment is the
    #: whole path; a tree record's fragment is the directory, and every child matches.
    fragments: tuple[str, ...]
    kind: str  # "file" or "tree"
    superseded_on: str
    superseded_by: str
    #: The measurement that makes the supersession a fact rather than an opinion.  UNAVAILABLE
    #: where nothing was counted -- never 0, never an estimate.
    counts: str
    why: str
    read_instead: str
    severity: str = REFUSE
    #: "<repo-relative file>::<function>" of the one caller permitted to read this anyway.
    allowed_readers: tuple[str, ...] = ()
    #: Basenames that identify this artifact wherever it has been copied or symlinked to.
    aliases: tuple[str, ...] = field(default_factory=tuple)


CELLS_CURRENT = "data/development/typing/cells/ (25,230,109 cells, 20 donors, written 2026-09-15)"

REGISTRY: tuple[Record, ...] = (
    Record(
        key="subpanel_calls",
        fragments=("data/development/typing/subpanel/subpanel_calls.parquet",),
        aliases=("subpanel_calls.parquet",),
        kind="file",
        superseded_on="2026-09-15",
        superseded_by=CELLS_CURRENT,
        counts=(
            "25,150,537 rows written 2026-09-11 against typing/cells 25,230,109 rows written "
            "2026-09-15 -- 79,572 cells apart; endocrine "
            "Beta 191,485 / Alpha 246,055 / Delta 62,797 "
            "against Beta 226,183 / Alpha 279,686 / Delta 90,658"
        ),
        why=(
            "the PREVIOUS generation of the calls. Its endocrine counts differ materially from "
            "the build on disk, so a number taken from it is a number from a method that was "
            "replaced four days later. The typing tree stacks generations side by side with "
            "nothing in the filename saying which is which"
        ),
        read_instead=(
            "data/development/typing/cells/ for the calls and "
            "data/development/typing/islet_membership/ for islet membership"
        ),
        allowed_readers=("core/phenocycler/scripts/typing/subpanel.py::gate",),
    ),
    Record(
        key="umap_refinement_log",
        fragments=("data/development/typing/umap/refinement_log.json",),
        aliases=("refinement_log.json",),
        kind="file",
        superseded_on="2026-09-17",
        superseded_by="nothing on disk -- the current own-marker table is UNAVAILABLE",
        counts=(
            "6,156 bytes written 2026-09-11 13:38, an 8-class own_markers table over 21 markers; "
            "NO PRODUCER for it exists in the repo (grep: zero writers). The refinement that was "
            "written on 2026-09-17, scripts/typing/refine.py, writes refine_log.json instead and "
            "its marker set is 54, not 21"
        ),
        why=(
            "it is the 2026-09-11 generation's log and it has no producer, so it cannot be "
            "regenerated, checked, or dated against the chain it is being used to describe. "
            "refine.py is WRITTEN AND NOT WIRED and its 5th-percentile bound demotes about 5% of "
            "every class by construction -- that is an open question, not a decision, so its "
            "refine_log.json is not a replacement either. A missing own-marker table is "
            "UNAVAILABLE, not the 2026-09-11 one"
        ),
        read_instead=(
            "data/development/typing/subpanel/subpanel_log.json for the markers the "
            "classification actually used (that file is current and is NOT refused); "
            "core/phenocycler/scripts/typing/refine.py's OWN_MARKERS for the committed "
            "defining-marker table, and REFINEMENT.md for what is still open; for the immune "
            "definitions' markers, scripts/typing/RECONSTRUCTION.md, reading the "
            "[ON DISK] / [CHAIN NOW] mark before quoting it. refine.py also installs a "
            "process-wide open() hook that refuses this file when it is named by a literal "
            "path, which no GuardedPath can see"
        ),
    ),
    Record(
        key="content_addressed_runs",
        fragments=("core/phenocycler/data/runs",),
        kind="tree",
        superseded_on="2026-09-16",
        superseded_by=(
            "the stable stage roots core/phenocycler/data/{cells,geometry_qc,redsea,expression}"
        ),
        counts=(
            "5 run directories plus LATEST; the 20-donor cohort expression output exists in "
            "30f6689e4d27d295f890 and NOWHERE ELSE, which is why this is announced and not "
            "refused"
        ),
        why=(
            "the content-addressed run framework was removed on 2026-09-16. There is no run_id, "
            "no stage manifest, no donor receipt and no fingerprint validation any more, so "
            "nothing here is validated by anything and nothing writes to it"
        ),
        read_instead=(
            "core/phenocycler/data/<stage>/ for anything re-run since 2026-09-16; for the "
            "re-exported batch-2 donors, data/development/typing/revised_cells/ via "
            "subpanel_io.cell_source(), which prefers it over this tree"
        ),
        severity=ANNOUNCE,
    ),
    Record(
        key="cells_redsea",
        fragments=("data/cells_redsea",),
        kind="tree",
        superseded_on="2026-08-14",
        superseded_by="core/phenocycler/data/redsea/ (the REDSEA stage of the source workflow)",
        counts="ABSENT from this checkout, measured 2026-09-19: no data/cells_redsea on disk",
        why=(
            "a retired flat root from the lane RESTORE was retired with on 2026-08-14. It is on "
            "the never-revive list"
        ),
        read_instead="core/phenocycler/data/redsea/ and core/phenocycler/data/expression/",
    ),
    Record(
        key="restore_roots",
        fragments=(
            "data/restore_gated_redsea",
            "data/restore_gated_redsea_extra",
            "data/restore_pair_validation",
        ),
        kind="tree",
        superseded_on="2026-08-14",
        superseded_by=CELLS_CURRENT,
        counts=(
            "data/restore_gated_redsea is still on disk (directory mtime 2026-07-27); "
            "restore_pair_validation is still on disk; the rest are ABSENT. Row counts "
            "UNAVAILABLE -- nothing in this lane has counted them since RESTORE was retired"
        ),
        why=(
            "RESTORE was retired on 2026-08-14 and the pair webs are on the never-revive list. "
            "The ten <marker>_pos booleans are not a positivity convention this project holds: "
            "the current typing asks one question per context against reference levels drawn "
            "from the ordinary tissue, and never asks whether a cell is <marker>-positive as a "
            "standalone question. You cannot indict a marker, only a method"
        ),
        read_instead=(
            "data/development/typing/cells/ for the calls; for the run-cohort per-cell distance "
            "artifact, core/phenocycler/scripts/typing/build_app_data.py::write_cell_distance_"
            "parquet, which omits the _pos columns for the reason above"
        ),
    ),
    Record(
        key="retired_typing_chain",
        fragments=(
            "data/development/typing/astir",
            "data/development/typing/tacit",
            "data/development/typing/calls",
            "data/development/typing/positivity",
            "data/development/typing/qupath_export",
            "data/development/typing/substrate",
        ),
        kind="tree",
        superseded_on="2026-09-08",
        superseded_by=CELLS_CURRENT,
        counts=(
            "the withdrawn two-donor comparison arms; astir/tacit/calls carry per-class "
            "posteriors and positivity carries 39 per-marker pixel-positivity columns, none of "
            "which has a counterpart in typing/cells -- retired by ruling, not replaced, so "
            "their counts are not comparable to anything current and are UNAVAILABLE as such"
        ),
        why=(
            "the retired typing chain (substrate -> positivity -> call -> validate -> "
            "export_qupath) and the broad-typing producer comparison withdrawn on 2026-09-08. "
            "Hierarchical, simultaneous and probabilistic typing are on the never-revive list"
        ),
        read_instead=CELLS_CURRENT,
    ),
    Record(
        key="dct",
        fragments=("data/development/dct", "core/phenocycler/data/development/dct"),
        kind="tree",
        superseded_on="2026-09-08",
        superseded_by=CELLS_CURRENT,
        counts="ABSENT from this checkout, measured 2026-09-19: no data/development/dct on disk",
        why="the retired DeepCell-types comparison lane; withdrawn with the producer comparison",
        read_instead=CELLS_CURRENT,
    ),
    Record(
        key="broadtype",
        fragments=(
            "data/development/broadtype/release",
            "data/development/broadtype/arms",
            "data/development/broadtype/frame",
            "data/development/broadtype/verdict",
            "data/development/broadtype/logs",
        ),
        kind="tree",
        superseded_on="2026-09-08",
        superseded_by=CELLS_CURRENT,
        counts=(
            "the withdrawn broad-typing producer comparison's output. Cell counts UNAVAILABLE as "
            "a comparison: it is a different method over a different class set, so its numbers "
            "do not sit beside typing/cells' 25,230,109"
        ),
        why=(
            "the producer comparison was WITHDRAWN on 2026-09-08 and release/ must not be read "
            "as a verified typing. It is announced rather than refused because config.ini binds "
            "release/ and the --verify manifest declares frame/ files as inputs, so a refusal "
            "here would fail a lane that is not this one's to change"
        ),
        read_instead=CELLS_CURRENT,
        severity=ANNOUNCE,
    ),
)

#: Live inputs that sit inside a retired tree.  Checked BEFORE the registry, so the prefix cannot
#: sweep them up.  Each one is a path the current method reads today.
EXEMPT: tuple[str, ...] = (
    # subpanel_io.load_scale() reads cosine_scale.parquet out of here -- the cohort floor/top per
    # marker that the classification on disk was made on -- and five live modules cache here.
    "data/development/broadtype/scratch",
    # The classification's own log.  Current, and the replacement this module points at for the
    # 2026-09-11 refinement_log.json.  Named so no future subpanel/ prefix can catch it.
    "data/development/typing/subpanel/subpanel_log.json",
)

_BY_KEY = {r.key: r for r in REGISTRY}
#: Cheap pre-filter.  ``__fspath__`` runs on every open in this lane, so the common case must be
#: a handful of substring tests and nothing else -- no resolve(), no stat, no I/O.
_NEEDLES: tuple[str, ...] = tuple(sorted({f for r in REGISTRY for f in r.fragments}
                                         | {a for r in REGISTRY for a in r.aliases}))

_announced: set[str] = set()


def _norm(path) -> str:
    text = str(path).replace("\\", "/")
    while "//" in text:
        text = text.replace("//", "/")
    return text.rstrip("/")


def _matches(text: str, fragment: str) -> bool:
    """True when ``text`` names ``fragment`` itself or something inside it.

    Fragment boundaries are respected on both sides, so ``data/restore_gated_redsea`` does not
    match ``data/restore_gated_redsea_v2`` and ``typing/calls`` does not match
    ``celltype/calls`` -- a path segment is matched, never a bare substring.
    """
    index = text.find(fragment)
    while index != -1:
        before_ok = index == 0 or text[index - 1] == "/"
        end = index + len(fragment)
        after_ok = end == len(text) or text[end] == "/"
        if before_ok and after_ok:
            return True
        index = text.find(fragment, index + 1)
    return False


def _is_exempt(text: str) -> bool:
    """A live input that happens to sit inside a retired tree, or to share a basename."""
    return any(_matches(text, fragment) for fragment in EXEMPT)


@functools.lru_cache(maxsize=4096)
def _lookup(text: str) -> Record | None:
    # EXEMPT first, and before the needle pre-filter, so an exemption holds no matter what a
    # later registry entry's prefix happens to reach. It is cheap because this whole function is
    # cached per distinct path string.
    if _is_exempt(text):
        return None
    if not any(needle in text for needle in _NEEDLES):
        return None
    for record in REGISTRY:
        for fragment in record.fragments:
            if _matches(text, fragment):
                return record
    return None


_FILE_ALIASES: dict[str, Record] = {
    alias: record for record in REGISTRY if record.kind == "file" for alias in record.aliases
}


@functools.lru_cache(maxsize=4096)
def _match_text(text: str) -> Record | None:
    """String-only registry lookup.  No syscall, so it is safe on every ``__fspath__``."""
    found = _lookup(text)
    if found is not None:
        return found
    if _is_exempt(text):          # an exemption also beats the basename fallback
        return None
    return _FILE_ALIASES.get(text.rsplit("/", 1)[-1])


@functools.lru_cache(maxsize=4096)
def _match_target(text: str) -> Record | None:
    """What the path resolves to.  One ``resolve()``, cached per distinct string."""
    try:
        resolved = _norm(pathlib.Path(text).resolve())
    except (OSError, RuntimeError, ValueError):
        return None
    if resolved == text:
        return None
    found = _lookup(resolved)
    if found is not None and found.kind == "file":
        return found
    if _is_exempt(resolved):
        return None
    return _FILE_ALIASES.get(resolved.rsplit("/", 1)[-1])


def match(path, follow_symlinks: bool = True) -> Record | None:
    """The registry record naming ``path``, or None.

    A file record matches by its full path, by basename, and -- with ``follow_symlinks`` --
    by the basename of what the path resolves to.  An alias pointing at
    ``subpanel_calls.parquet`` IS that file, whatever the link is called; copying or linking a
    superseded artifact under another name is exactly how one gets read as a current one, so
    the spelling alone is not what is being trusted.

    ``follow_symlinks=False`` is the string-only test, used by :class:`GuardedPath` so that
    guarding costs no syscall on an ordinary read.  The explicit :func:`refuse_if_superseded`
    that a builder calls once per run does follow them.
    """
    text = _norm(path)
    found = _match_text(text)
    if found is not None or not follow_symlinks:
        return found
    return _match_target(text)


def notice(record: Record, path, verb: str = "Refusing superseded") -> str:
    """The whole of what a reader needs: what was refused, why, the counts, what to read."""
    return "\n".join((
        f"{verb} {_norm(path)}",
        f"  what it is         : {record.key} ({record.kind})",
        f"  superseded on      : {record.superseded_on}",
        f"  superseded by      : {record.superseded_by}",
        f"  measured           : {record.counts}",
        f"  why                : {record.why}",
        f"  read instead       : {record.read_instead}",
        "  the one table      : core/phenocycler/scripts/typing/superseded.py"
        f" (record {record.key!r})",
    ))


def _caller_is_allowed(record: Record) -> str | None:
    """The recorded allowed reader currently on the stack, or None.

    Only the callers written into the record match, by file tail AND function name together.
    Twenty frames is far deeper than any read path in this lane, and the walk only happens on a
    path that has already matched the registry, so it never runs on an ordinary read.
    """
    if not record.allowed_readers:
        return None
    wanted = {}
    for entry in record.allowed_readers:
        tail, _, function = entry.partition("::")
        wanted[(tail, function)] = entry
    frame = inspect.currentframe()
    try:
        for _ in range(20):
            frame = frame.f_back if frame is not None else None
            if frame is None:
                return None
            filename = _norm(frame.f_code.co_filename)
            function = frame.f_code.co_name
            for (tail, wanted_function), entry in wanted.items():
                if wanted_function == function and filename.endswith(tail):
                    return entry
    finally:
        del frame
    return None


def announce(record: Record, path, verb: str, once_key: str | None = None) -> None:
    """Print the full notice to stderr, at most once per key per process."""
    key = once_key if once_key is not None else record.key
    if key in _announced:
        return
    _announced.add(key)
    print(notice(record, path, verb), file=sys.stderr, flush=True)


def refuse_if_superseded(path, follow_symlinks: bool = True):
    """Refuse, by name and with the counts, to read a superseded artifact.

    Returns ``path`` unchanged when it is not registered, so it can be used inline::

        table = pq.read_table(refuse_if_superseded(candidate))

    An ``ANNOUNCE`` record prints the same notice to stderr and returns; a ``REFUSE`` record
    raises :class:`SupersededArtifactError`, unless the calling frame is one this record's
    ``allowed_readers`` names -- and that read is announced too.  Nothing here is silent.
    """
    record = match(path, follow_symlinks=follow_symlinks)
    if record is None:
        return path
    if record.severity == ANNOUNCE:
        announce(record, path, "Reading historical (not current)")
        return path
    allowed = _caller_is_allowed(record)
    if allowed is not None:
        announce(record, path, f"Reading superseded, allowed for {allowed}",
                 once_key=f"{record.key}:{allowed}")
        return path
    raise SupersededArtifactError(notice(record, path))


def is_superseded(path) -> bool:
    """True when the registry names this path at all, at either severity."""
    return match(path) is not None


_GuardedBase = type(pathlib.Path())


class GuardedPath(_GuardedBase):  # type: ignore[misc,valid-type]
    """A path that refuses to be opened when it names a superseded artifact.

    ``pathlib`` builds derived paths through ``with_segments``, so ``GuardedPath(root) / "a" / "b"``
    is a ``GuardedPath`` too.  That is the whole mechanism: ``subpanel_io`` publishes ``OUT`` and
    ``RUN`` as ``GuardedPath``, every module in this lane composes its paths from them, and the
    guard is therefore attached to the artifact rather than to the thirty-three readers.

    Only the reading doors are guarded -- ``__fspath__``, ``open``, ``read_text``, ``read_bytes``.
    ``str()``, ``exists()``, ``mkdir()`` and the rest are untouched, so a caller can still name,
    test for and report a superseded path without tripping anything.  Refusing ``str()`` would
    make the refusal message itself impossible to write.
    """

    def __fspath__(self) -> str:
        refuse_if_superseded(_GuardedBase.__str__(self), follow_symlinks=False)
        return _GuardedBase.__fspath__(self)

    def open(self, mode="r", *args, **kwargs):
        if "r" in mode and "+" not in mode:
            refuse_if_superseded(_GuardedBase.__str__(self), follow_symlinks=False)
        return _GuardedBase.open(self, mode, *args, **kwargs)

    def read_text(self, *args, **kwargs) -> str:
        refuse_if_superseded(_GuardedBase.__str__(self), follow_symlinks=False)
        return _GuardedBase.read_text(self, *args, **kwargs)

    def read_bytes(self) -> bytes:
        refuse_if_superseded(_GuardedBase.__str__(self), follow_symlinks=False)
        return _GuardedBase.read_bytes(self)


def guarded(path) -> GuardedPath:
    """``path`` as a :class:`GuardedPath`, so everything derived from it is guarded too."""
    return GuardedPath(os.fspath(path) if not isinstance(path, str) else path)


def table() -> str:
    """The registry as text -- what this module knows, for a reader or a report."""
    lines = []
    for record in REGISTRY:
        lines.append(f"[{record.severity.upper():8}] {record.key}")
        for fragment in record.fragments:
            on_disk = (REPO / fragment).exists()
            lines.append(f"             {fragment}  ({'on disk' if on_disk else 'absent'})")
        lines.append(f"             superseded {record.superseded_on} by {record.superseded_by}")
        lines.append(f"             measured: {record.counts}")
        if record.allowed_readers:
            lines.append(f"             allowed reader: {', '.join(record.allowed_readers)}")
        lines.append("")
    lines.append("exempt (live inputs inside a retired tree):")
    lines.extend(f"             {fragment}" for fragment in EXEMPT)
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", help="paths to check against the registry")
    args = parser.parse_args(argv)
    if not args.paths:
        print(table())
        return 0
    status = 0
    for candidate in args.paths:
        record = match(candidate)
        if record is None:
            print(f"ok        {candidate}")
            continue
        verb = "REFUSED  " if record.severity == REFUSE else "ANNOUNCED"
        print(f"{verb} {candidate}")
        print(notice(record, candidate,
                     "Refusing superseded" if record.severity == REFUSE
                     else "Reading historical (not current)"))
        if record.severity == REFUSE:
            status = 1
    return status


if __name__ == "__main__":
    raise SystemExit(main())
