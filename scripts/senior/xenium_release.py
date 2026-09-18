#!/usr/bin/env python3
"""One place that answers "which Xenium release is current for this tissue".

WHY THIS EXISTS.  Measured 2026-09-18, five builders named five different Xenium
sources and only one of them read the current pair:

    build_availability.py            core/phenocycler/data/integration/cells_xen
                                     -- the UNTYPED ingest
    build_distance_profiles.py       the newest release.json by DIRECTORY MTIME
    build_compartment_gradient.py    release_pln_immune_v2, superseded by _v3
    modality_concordance.py          release_pln_immune_v2, release_panc_v1, release
                                     -- release_panc_v1 holds no calls/ at all, so the
                                        pancreas comparison fell through to a 2-donor tree
    build_grid_composition.py        /data/xenium_typing/panc_2026-09-17/release
                                     and release_pln_immune_v3 -- the only current pair

The consequence reached the app: ``availability/data.parquet`` reported
``xenium_panc`` blank on 18 of 20 donors and ``xenium_pln`` blank on 13 of 13
while both releases carried calls for them.

WHAT THIS MODULE PROMISES.

* It searches BOTH roots.  The pancreas release is written outside the repo, so a
  repo-relative search alone can never see it.  Roots are an argument, then an
  environment variable, then a default -- never a literal in a builder.
* It selects on the release's OWN DECLARED METADATA: ``release_id``, the declared
  write time (``generated_at_utc`` / ``written_utc`` / ``generated_at``), and the
  declared ``tissue``/``roi``.  **Nothing here reads a filesystem mtime.**  A
  ``touch`` must not be able to change which release the cohort is built from,
  and the previous mtime scan could be reordered by one.
* It honours a release that declares it supersedes another.  ``release_pln_immune_v3``
  records ``producer.supersedes: data/development/xenium_typing/release_pln_immune_v2``,
  and v2 in turn supersedes the draft.  A superseded release is never current even
  if its timestamp says otherwise.
* A release is only a candidate once a REAL ``calls/`` tree is behind it: the
  directory exists and holds at least one ``donor_id=`` partition with at least one
  parquet file, and a candidate that fails that is carried with the reason it
  failed rather than dropped.  ``release_panc_v1`` never becomes a candidate at
  all -- it holds only ``_work/6539.*``, no ``calls/`` and no ``release.json``, so
  it declares nothing -- and the comparison that used to fall through it to a
  two-donor tree now says which release it read.
* When nothing is current it returns a NAMED :class:`Unavailable` carrying the
  reason and every candidate it rejected and why.  It never returns ``None``,
  because a ``None`` that reaches a count reads as a zero, and a missing
  measurement is unavailable -- not zero.
* It carries the releases' own reading instruction.  The pancreas release spells
  it out in ``partition_key_note`` / ``read_calls_with``: read ``calls/`` as one
  hive dataset with ``donor_id`` as a **STRING**, because every donor id in this
  cohort is all digits and an inferring reader turns ``"112"`` into ``112`` and it
  stops matching the donor it names.  Measured on this checkout: inferred hive
  partitioning types ``donor_id`` as ``int32``.  :meth:`Release.partitioning`
  never infers.

WHAT IT DOES NOT DO.  It does not filter the cohort.  ``DONOR_EXCLUSIONS`` lives
once, in ``core/phenocycler/phenocycler/cohort.py``, and the builders apply it;
restating it here would be a second copy.  It does not type anything, does not
move the ``NO_VERIFIED_TYPING`` latch, and reads no result as a target.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]

#: The environment variable that overrides the search roots entirely, as an
#: ``os.pathsep``-separated list.  Set it and NOTHING else is searched.
ROOTS_ENV = "XENIUM_RELEASE_ROOTS"

#: The root outside the repository.  The pancreas release is written here and a
#: repo-relative search alone cannot see it.  Overridable so a test never reaches
#: real data through a fixture repo root.
EXTERNAL_ROOT_ENV = "XENIUM_RELEASE_EXTERNAL_ROOT"
EXTERNAL_ROOT_DEFAULT = Path("/data/xenium_typing")

#: The repository-relative root, relative to the repo root handed in.
REPO_RELATIVE_ROOT = Path("data") / "development" / "xenium_typing"

#: Canonical tissue tokens.  Everything else is an alias of one of these.
PANC = "panc"
PLN = "pln"
TISSUES = (PANC, PLN)

#: How a release may spell its tissue.  ``build_distance_profiles`` uses the same
#: two spellings for the pLN, so the tokens agree across the lanes.
TISSUE_ALIASES = {
    "panc": PANC,
    "pancreas": PANC,
    "pancreatic": PANC,
    "pln": PLN,
    "pln_1": PLN,
    "lymph_node": PLN,
    "lymphnode": PLN,
    "pancreatic_lymph_node": PLN,
}

#: The keys a release may use for its own declared write time, in the order they
#: are trusted.  There is no mtime fallback: a release that declares no time is
#: rejected by name.
TIME_KEYS = ("generated_at_utc", "written_utc", "generated_at", "written_at")

#: Column roles, and the spellings a release may use for each.  A role is
#: resolved through the release's declared ``columns`` map first; the candidates
#: are the fallback for a release that declares none (the pLN releases do not).
ROLE_CANDIDATES: Mapping[str, tuple[str, ...]] = {
    "cell_key": ("cell_id",),
    "label": ("lineage", "lineage_native"),
    "common_lineage": ("common_lineage", "lineage_common"),
    "reason": ("reason",),
    "section_key": ("section_key",),
    "nucleus_count": ("nucleus_count",),
    "nucleus_measured": ("nucleus_measured",),
    "depth": ("depth",),
}

#: How a label is read.  These are STATES OF THE CALL, not lineages.
CALL = "called"            # a real determination, including a measured negative
UNCALLED = "uncalled"      # looked at, no call resolved
LABEL_UNAVAILABLE = "unavailable"  # the release itself says the answer is unavailable
BLANK = "blank"            # no label written at all
CALL_STATES = (CALL, UNCALLED, LABEL_UNAVAILABLE, BLANK)


class XeniumReleaseError(RuntimeError):
    """A release contract this module will not paper over."""


def classify_label(value: Any) -> str:
    """One label -> one of :data:`CALL_STATES`.

    ``Non_immune`` is a CALL: the pLN lanes resolve immune identity only, and a
    cell measured as not immune has been answered.  Collapsing it into "uncalled"
    would turn a measured negative into an absence.

    ``Unavailable_no_nucleus`` and any other ``Unavailable*`` label is the
    release saying the answer is unavailable, which is not the same as blank and
    is never a zero.  ``Uncalled_*`` and ``*_uncalled`` are "looked at, nothing
    resolved" -- both spellings are in use (``Uncalled_ambiguous`` in the pLN
    release, ``Immune_uncalled`` in the same file).
    """
    if value is None:
        return BLANK
    text = str(value).strip()
    if not text:
        return BLANK
    if text.startswith("Unavailable"):
        return LABEL_UNAVAILABLE
    if text.startswith("Uncalled") or text.endswith("_uncalled"):
        return UNCALLED
    return CALL


# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------
def default_roots(repo: Path | str = REPO) -> tuple[Path, ...]:
    """Both roots, unless the environment or the caller names its own.

    ``$XENIUM_RELEASE_ROOTS`` REPLACES the list rather than extending it, so a
    caller that names roots gets exactly those.

    THE EXTERNAL ROOT IS MACHINE CONFIGURATION, not something derivable from a
    repository path.  It is searched for THIS checkout, or wherever
    ``$XENIUM_RELEASE_EXTERNAL_ROOT`` is set explicitly.  A caller that hands in a
    different repository root -- a synthetic fixture tree, a second checkout --
    gets that tree and nothing else, so a build against a tree someone made for a
    test can never silently reach a real release it did not create.  That is not a
    convenience: a fixture that reads real data is measuring the machine.
    """
    override = os.environ.get(ROOTS_ENV, "").strip()
    if override:
        return tuple(Path(p).expanduser() for p in override.split(os.pathsep) if p.strip())
    roots = [Path(repo) / REPO_RELATIVE_ROOT]
    named = os.environ.get(EXTERNAL_ROOT_ENV, "").strip()
    if named:
        roots.append(Path(named).expanduser())
    else:
        try:
            is_this_checkout = Path(repo).resolve() == REPO.resolve()
        except OSError:
            is_this_checkout = False
        if is_this_checkout:
            roots.append(EXTERNAL_ROOT_DEFAULT)
    return tuple(roots)


# ---------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------
def _parse_time(declared: Mapping[str, Any]) -> tuple[datetime | None, str]:
    for key in TIME_KEYS:
        raw = declared.get(key)
        if not raw:
            continue
        text = str(raw).strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        try:
            stamp = datetime.fromisoformat(text)
        except ValueError:
            continue
        if stamp.tzinfo is None:
            stamp = stamp.replace(tzinfo=timezone.utc)
        return stamp.astimezone(timezone.utc), key
    return None, ""


def _declared_tissue(declared: Mapping[str, Any]) -> str | None:
    """The tissue the release declares, from ``tissue`` then ``roi``.

    The pancreas release declares no ``tissue`` key at all and only ``roi:
    "panc"``; the pLN draft declares ``roi`` as a LIST.  Both are the release's
    own words, which is the point -- nothing here reads a directory name.
    """
    token = str(declared.get("tissue", "")).strip().lower()
    if token in TISSUE_ALIASES:
        return TISSUE_ALIASES[token]
    for value in _declared_roi(declared):
        token = value.strip().lower()
        if token in TISSUE_ALIASES:
            return TISSUE_ALIASES[token]
        if token.startswith("pln"):
            return PLN
    return None


def _declared_roi(declared: Mapping[str, Any]) -> tuple[str, ...]:
    raw = declared.get("roi")
    if raw is None:
        return ()
    if isinstance(raw, str):
        return (raw,) if raw.strip() else ()
    if isinstance(raw, Sequence):
        return tuple(str(v) for v in raw if str(v).strip())
    return (str(raw),)


def _supersedes_text(declared: Mapping[str, Any]) -> str:
    """Whatever the release says it supersedes, flattened to text.

    ``release_pln_immune_v3`` writes it as free prose beginning with a path:
    ``"data/development/xenium_typing/release_pln_immune_v2 (its partition
    contained ...)"``.  A future release may write a list or a mapping; all three
    flatten here, and the MATCHING below is on resolved paths and release ids, not
    on substrings, so ``release_pln_immune`` is not mistaken for
    ``release_pln_immune_v2``.
    """
    parts: list[str] = []
    for holder in (declared.get("producer") or {}, declared):
        if not isinstance(holder, Mapping):
            continue
        value = holder.get("supersedes")
        if value is None:
            continue
        if isinstance(value, str):
            parts.append(value)
        elif isinstance(value, Mapping):
            parts.extend(str(v) for v in value.values())
        elif isinstance(value, Sequence):
            parts.extend(str(v) for v in value)
        else:
            parts.append(str(value))
    return " ".join(parts)


_TOKEN = re.compile(r"[^\s,;()\[\]{}'\"]+")


def _path_tokens(text: str) -> tuple[str, ...]:
    return tuple(m.group(0) for m in _TOKEN.finditer(text))


def _calls_state(root: Path, declared: Mapping[str, Any]) -> tuple[Path, str]:
    """(calls dir, "" if usable else the reason it is not).

    "A real ``calls/`` directory" is not "a directory called calls": it must hold
    at least one ``donor_id=`` partition with at least one parquet file in it.
    """
    calls = root / str(declared.get("calls_dir") or "calls")
    if not calls.is_dir():
        return calls, f"declares no usable calls tree: {calls} is not a directory"
    partitions = sorted(calls.glob("donor_id=*"))
    if not partitions:
        return calls, f"{calls} holds no donor_id= partition"
    for part in partitions:
        if next(iter(part.rglob("*.parquet")), None) is not None:
            return calls, ""
    return calls, f"{calls} holds donor_id= partitions but no parquet file in any of them"


@dataclass(frozen=True)
class Candidate:
    """One ``release.json`` found under a root, and what it declares."""

    root: Path
    manifest: Path
    declared: Mapping[str, Any]
    release_id: str
    tissue: str | None
    generated_at: datetime | None
    generated_at_field: str
    roi: tuple[str, ...]
    calls_dir: Path
    supersedes: str
    rejected: str = ""

    def describe(self) -> dict:
        return {
            "root": str(self.root),
            "release_id": self.release_id,
            "tissue": self.tissue,
            "generated_at": self.generated_at.isoformat() if self.generated_at else None,
            "generated_at_field": self.generated_at_field or None,
            "roi": list(self.roi),
            "rejected": self.rejected or None,
        }


def _iter_manifests(root: Path) -> Iterator[Path]:
    """Every ``release.json`` under ``root``, never descending INTO a release.

    ``release_pln_immune_v3/_partial/sections=6347/release.json`` is a real
    manifest with the same ``release_id`` and a real ``calls/`` tree: it is one
    section of a partial run, which the release's own scope rule says lives under
    ``_partial/``.  Stopping the walk at the first manifest keeps it out, and the
    ``scope.is_partial`` check below keeps out any that another layout lets past.
    """
    if not root.is_dir():
        return
    stack = [root]
    while stack:
        here = stack.pop()
        manifest = here / "release.json"
        if manifest.is_file():
            yield manifest
            continue
        try:
            children = sorted(p for p in here.iterdir() if p.is_dir())
        except OSError:
            continue
        stack.extend(children)


def scan(roots: Iterable[Path]) -> list[Candidate]:
    """Every candidate under every root, each carrying its rejection or "".

    Rejections are recorded rather than dropped, so an Unavailable can name what
    it looked at and why each one lost.
    """
    seen: set[Path] = set()
    out: list[Candidate] = []
    for root in roots:
        root = Path(root)
        for manifest in _iter_manifests(root):
            here = manifest.parent.resolve()
            if here in seen:
                continue
            seen.add(here)
            try:
                declared = json.loads(manifest.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                out.append(Candidate(here, manifest, {}, "", None, None, "", (),
                                     here / "calls", "",
                                     rejected=f"release.json is unreadable: {exc}"))
                continue
            if not isinstance(declared, Mapping):
                out.append(Candidate(here, manifest, {}, "", None, None, "", (),
                                     here / "calls", "",
                                     rejected="release.json is not an object"))
                continue
            release_id = str(declared.get("release_id", "")).strip()
            tissue = _declared_tissue(declared)
            stamp, stamp_field = _parse_time(declared)
            calls, calls_reason = _calls_state(here, declared)
            scope = declared.get("scope")
            reason = ""
            if not release_id:
                reason = "declares no release_id"
            elif tissue is None:
                reason = ("declares neither a tissue nor an roi this module recognises "
                          f"(tissue={declared.get('tissue')!r}, roi={declared.get('roi')!r})")
            elif stamp is None:
                reason = (f"declares no parsable write time under any of {list(TIME_KEYS)}; "
                          "a directory mtime is not a substitute and is never read here")
            elif isinstance(scope, Mapping) and scope.get("is_partial"):
                reason = (f"declares scope.is_partial with slug {scope.get('slug')!r}: "
                          "a partial run is not the cohort release")
            elif calls_reason:
                reason = calls_reason
            out.append(Candidate(
                root=here, manifest=manifest, declared=declared, release_id=release_id,
                tissue=tissue, generated_at=stamp, generated_at_field=stamp_field,
                roi=_declared_roi(declared), calls_dir=calls,
                supersedes=_supersedes_text(declared), rejected=reason,
            ))
    return out


def _superseded_by(candidate: Candidate, others: Sequence[Candidate]) -> Candidate | None:
    """The candidate that declares it supersedes ``candidate``, if any.

    Matched on a RESOLVED PATH or an exact ``release_id``, never on a substring:
    ``release_pln_immune`` is a prefix of ``release_pln_immune_v2`` and a
    substring test would retire the wrong one.
    """
    target = candidate.root.resolve()
    for other in others:
        if other is candidate or not other.supersedes:
            continue
        tokens = _path_tokens(other.supersedes)
        if candidate.release_id and candidate.release_id in tokens:
            return other
        for token in tokens:
            for base in (Path(token), REPO / token, other.root.parent / token):
                try:
                    if base.resolve() == target:
                        return other
                except OSError:
                    continue
    return None


# ---------------------------------------------------------------------------
# The two return types
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Unavailable:
    """No current release for this tissue, and the reason, by name.

    It deliberately carries NO ``root`` and NO ``calls_dir``.  Code that assumed a
    release raises ``AttributeError`` here rather than counting an absence as a
    zero.  Check ``.available`` -- truthiness is not overloaded, because a
    ``if not release:`` that quietly succeeds is the failure this type exists to
    prevent.
    """

    tissue: str
    reason: str
    roots: tuple[Path, ...] = ()
    rejected: tuple[dict, ...] = ()

    available: bool = field(default=False, init=False)

    def describe(self) -> dict:
        return {
            "tissue": self.tissue,
            "available": False,
            "reason": self.reason,
            "roots_searched": [str(r) for r in self.roots],
            "candidates_rejected": list(self.rejected),
        }

    def __str__(self) -> str:  # pragma: no cover - diagnostics
        return f"no current Xenium {self.tissue} release: {self.reason}"


@dataclass(frozen=True)
class Release:
    """A current release, and every way the builders need to read it."""

    tissue: str
    root: Path
    manifest: Path
    release_id: str
    generated_at: datetime
    generated_at_field: str
    declared_roi: tuple[str, ...]
    calls_dir: Path
    supersedes: str
    metadata: Mapping[str, Any]

    available: bool = field(default=True, init=False)

    # -- declared facts ----------------------------------------------------
    @property
    def producer(self) -> Mapping[str, Any]:
        value = self.metadata.get("producer")
        return value if isinstance(value, Mapping) else {}

    @property
    def lane(self) -> str:
        return str(self.producer.get("lane", ""))

    @property
    def declared_donors(self) -> tuple[str, ...]:
        """The donors the release says it covers, as STRINGS, in its own order."""
        raw = self.metadata.get("donors")
        if isinstance(raw, Sequence) and not isinstance(raw, str):
            return tuple(str(d) for d in raw)
        return ()

    @property
    def read_calls_note(self) -> str:
        """The release's own instruction for reading ``calls/``.

        The pancreas release writes both ``partition_key_note`` and
        ``read_calls_with``; the pLN releases write neither, so this is "" there
        and :meth:`partitioning` supplies the same rule from the data.
        """
        parts = [str(self.metadata.get(k, "")).strip()
                 for k in ("partition_key_note", "read_calls_with")]
        return "\n".join(p for p in parts if p)

    # -- the calls tree ----------------------------------------------------
    def partition_files(self) -> tuple[tuple[str, str | None, Path], ...]:
        """``(donor_id, roi, parquet)`` for every file, keys parsed as STRINGS.

        The hive segment is read as text.  Not inferred, not cast: ``"112"`` is a
        donor id, and the moment it becomes ``112`` it stops matching the donor it
        names.
        """
        out: list[tuple[str, str | None, Path]] = []
        for part in sorted(self.calls_dir.glob("donor_id=*")):
            donor = part.name.split("=", 1)[1]
            for path in sorted(part.rglob("*.parquet")):
                roi = None
                for segment in path.relative_to(part).parts:
                    if segment.startswith("roi="):
                        roi = segment.split("=", 1)[1]
                        break
                out.append((donor, roi, path))
        return tuple(out)

    def donor_partitions(self) -> dict[str, list[Path]]:
        """``donor_id`` (string) -> its parquet files."""
        out: dict[str, list[Path]] = {}
        for donor, _roi, path in self.partition_files():
            out.setdefault(donor, []).append(path)
        return out

    def observed_donors(self) -> tuple[str, ...]:
        return tuple(sorted(self.donor_partitions(), key=lambda d: (len(d), d)))

    def observed_roi(self) -> tuple[str, ...]:
        return tuple(sorted({roi for _d, roi, _p in self.partition_files() if roi}))

    def partition_keys(self) -> tuple[str, ...]:
        """Declared partition keys, else the ones the tree actually uses."""
        raw = self.metadata.get("partition_keys")
        if isinstance(raw, Sequence) and not isinstance(raw, str) and raw:
            return tuple(str(k) for k in raw)
        files = self.partition_files()
        if not files:
            return ("donor_id",)
        _donor, roi, _path = files[0]
        return ("donor_id", "roi") if roi is not None else ("donor_id",)

    def calls_schema(self) -> pa.Schema:
        files = self.partition_files()
        if not files:
            raise XeniumReleaseError(f"{self.root}: calls tree holds no parquet file")
        return pq.read_schema(files[0][2])

    def partitioning(self) -> ds.Partitioning:
        """Hive partitioning with ``donor_id`` as a STRING, never inferred.

        Two measured facts shape this:

        * Inferred hive partitioning types ``donor_id`` as ``int32`` on the
          pancreas release, which is exactly what the release's own
          ``partition_key_note`` warns about.
        * The pLN releases ALSO carry ``donor_id`` and ``roi`` as physical
          ``large_string`` columns, so declaring them as ``string`` raises
          ``ArrowTypeError: Unable to merge: Field donor_id has incompatible
          types: large_string vs string``.  The key therefore takes the file's own
          string type where the file has one, and ``string`` where it does not.
          Both are strings; neither is inferred.
        """
        names = set(self.calls_schema().names)
        fields = []
        for key in self.partition_keys():
            if key in names:
                declared = self.calls_schema().field(key).type
                if not (pa.types.is_string(declared) or pa.types.is_large_string(declared)):
                    raise XeniumReleaseError(
                        f"{self.root}: partition key {key!r} is a physical column of type "
                        f"{declared}, which is not a string. A donor id that is not a string "
                        "stops matching the donor it names."
                    )
                fields.append(pa.field(key, declared))
            else:
                fields.append(pa.field(key, pa.string()))
        return ds.partitioning(pa.schema(fields), flavor="hive")

    def dataset(self) -> ds.Dataset:
        return ds.dataset(self.calls_dir, format="parquet", partitioning=self.partitioning())

    def read_calls(self, columns: Sequence[str] | None = None,
                   donors: Iterable[str] | None = None) -> pa.Table:
        """The whole calls tree, or the named donors of it, as one table."""
        data = self.dataset()
        expression = None
        if donors is not None:
            wanted = [str(d) for d in donors]
            expression = ds.field("donor_id").isin(wanted)
        return data.to_table(columns=list(columns) if columns else None, filter=expression)

    # -- columns -----------------------------------------------------------
    def column(self, role: str) -> str | None:
        """The column that fills ``role`` in THIS release, or None.

        The declared ``columns`` map wins; the pancreas release has one and the
        pLN releases have none, so the fallback candidates carry them.  Guessing
        wrong reads every cell as uncalled, so a missing role comes back as None
        and the caller names the absence.
        """
        names = set(self.calls_schema().names)
        declared = self.metadata.get("columns")
        if isinstance(declared, Mapping):
            value = declared.get(role)
            if value and str(value) in names:
                return str(value)
        for candidate in ROLE_CANDIDATES.get(role, ()):
            if candidate in names:
                return candidate
        return None

    def require_column(self, role: str) -> str:
        name = self.column(role)
        if name is None:
            raise XeniumReleaseError(
                f"{self.release_id}: no column fills the role {role!r}; the calls schema is "
                f"{sorted(self.calls_schema().names)}"
            )
        return name

    # -- sidecars ----------------------------------------------------------
    def sidecar(self, name: str) -> Path | None:
        """A file beside ``release.json`` (``per_donor_summary.csv`` and friends)."""
        path = self.root / name
        return path if path.is_file() else None

    def require_sidecar(self, name: str) -> Path:
        path = self.sidecar(name)
        if path is None:
            raise XeniumReleaseError(f"{self.release_id}: {self.root / name} is not on disk")
        return path

    # -- provenance --------------------------------------------------------
    def describe(self) -> dict:
        return {
            "tissue": self.tissue,
            "available": True,
            "root": str(self.root),
            "release_id": self.release_id,
            "generated_at": self.generated_at.isoformat(),
            "generated_at_field": self.generated_at_field,
            "declared_roi": list(self.declared_roi),
            "lane": self.lane,
            "supersedes": self.supersedes,
            "n_donors_declared": len(self.declared_donors),
            "n_donors_on_disk": len(self.observed_donors()),
            "selected_by": ("the release's own declared release_id, write time and tissue; "
                            "never a directory mtime"),
        }

    def __str__(self) -> str:  # pragma: no cover - diagnostics
        return f"{self.release_id} ({self.tissue}) at {self.root}"


# ---------------------------------------------------------------------------
# resolve
# ---------------------------------------------------------------------------
def canonical_tissue(tissue: str) -> str:
    token = str(tissue).strip().lower()
    if token in TISSUE_ALIASES:
        return TISSUE_ALIASES[token]
    raise XeniumReleaseError(
        f"unknown tissue {tissue!r}; this cohort has {list(TISSUES)} and nothing else"
    )


def resolve(tissue: str, *, roots: Iterable[Path] | None = None,
            repo: Path | str = REPO) -> Release | Unavailable:
    """The current release for ``tissue``, or a named :class:`Unavailable`.

    Selection, in order, all of it from the release's own declaration:

    1. the tissue it declares (``tissue``, else ``roi``);
    2. a real ``calls/`` tree behind it;
    3. not declared superseded by another candidate, and not a partial run;
    4. the latest DECLARED write time; ties broken on ``release_id`` then path so
       two runs of this function agree.

    No step reads a filesystem mtime.
    """
    want = canonical_tissue(tissue)
    search = tuple(Path(r) for r in roots) if roots is not None else default_roots(repo)
    candidates = scan(search)
    mine = [c for c in candidates if c.tissue == want]
    if not mine:
        return Unavailable(
            tissue=want,
            reason=("no release under any searched root declares this tissue"
                    if candidates else "no release.json found under any searched root"),
            roots=search,
            rejected=tuple(c.describe() for c in candidates),
        )

    rejected = [c for c in mine if c.rejected]
    usable = [c for c in mine if not c.rejected]
    for c in list(usable):
        winner = _superseded_by(c, mine)
        if winner is not None:
            usable.remove(c)
            rejected.append(Candidate(
                **{**c.__dict__,
                   "rejected": (f"declared superseded by {winner.release_id} "
                                f"({winner.root})")}))
    if not usable:
        return Unavailable(
            tissue=want,
            reason=(f"{len(mine)} release(s) declare tissue {want!r} and none is current; "
                    "see candidates_rejected for why each one lost"),
            roots=search,
            rejected=tuple(c.describe() for c in rejected),
        )

    best = max(usable, key=lambda c: (c.generated_at, c.release_id, str(c.root)))
    return Release(
        tissue=want, root=best.root, manifest=best.manifest, release_id=best.release_id,
        generated_at=best.generated_at, generated_at_field=best.generated_at_field,
        declared_roi=best.roi, calls_dir=best.calls_dir, supersedes=best.supersedes,
        metadata=best.declared,
    )


def require(tissue: str, *, roots: Iterable[Path] | None = None,
            repo: Path | str = REPO) -> Release:
    """:func:`resolve`, but an absence raises instead of travelling onward."""
    found = resolve(tissue, roots=roots, repo=repo)
    if not found.available:
        raise XeniumReleaseError(str(found))
    return found  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# One tally, shared, because two builders need exactly the same one
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DonorTally:
    """What one donor's partition holds, with the two absences kept apart.

    ``n_no_nucleus_excluded`` is the rule: a cell whose nucleus WAS measured and
    is absent is not a cell we can be sure of, so it is excluded from every
    calculation.  ``n_nucleus_not_measured`` is a MISSING MEASUREMENT -- the
    nucleus was never measured on that cell, the rule never applied to it, and
    counting it as an absent nucleus would silently drop whole donors.  Measured
    2026-09-17 in the pancreas release: every cell of 6347, 6579 and 6580.  The
    two are counted apart and reported apart, and neither is ever a zero.
    """

    donor_id: str
    n_rows: int
    n_called: int
    n_uncalled: int
    n_label_unavailable: int
    n_blank: int
    n_no_nucleus_excluded: int
    n_nucleus_not_measured: int
    nucleus_rule_applied: bool

    def reason(self) -> str:
        nucleus = (f"no_nucleus_excluded={self.n_no_nucleus_excluded}; "
                   f"nucleus_not_measured={self.n_nucleus_not_measured}"
                   if self.nucleus_rule_applied
                   else "this release carries no nucleus column, so the nucleus rule "
                        "could NOT be applied here and no count from it is reported")
        return (f"called={self.n_called}; uncalled={self.n_uncalled}; "
                f"label_unavailable={self.n_label_unavailable}; blank={self.n_blank}; "
                f"{nucleus}")


def donor_tallies(release: Release) -> dict[str, DonorTally]:
    """Per donor, how the calls read -- the one place this counting lives.

    Nothing is inferred from a filename, a row order or a centroid: the donor is
    the hive segment, read as text, and the label is the release's own column.
    """
    label_column = release.require_column("label")
    nucleus_count = release.column("nucleus_count")
    nucleus_measured = release.column("nucleus_measured")
    columns = [label_column] + [c for c in (nucleus_count, nucleus_measured) if c]

    out: dict[str, DonorTally] = {}
    for donor, files in release.donor_partitions().items():
        counts = dict.fromkeys(CALL_STATES, 0)
        rows = no_nucleus = not_measured = 0
        for path in files:
            table = pq.read_table(path, columns=columns)
            rows += table.num_rows
            if nucleus_count is not None:
                # A null flag is NOT a measurement: it is treated as never
                # measured, the same way the grid layer reads it.
                measured = (pc.fill_null(table[nucleus_measured], False)
                            if nucleus_measured is not None
                            else pa.chunked_array([pa.array([True] * table.num_rows)]))
                not_measured += pc.sum(pc.invert(measured)).as_py() or 0
                count = pc.cast(table[nucleus_count], pa.int64(), safe=False)
                # A null count on a cell whose nucleus WAS measured is an absent
                # nucleus, not a licence to keep the cell.
                absent = pc.and_(measured, pc.fill_null(pc.less_equal(count, 0), True))
                no_nucleus += pc.sum(absent).as_py() or 0
                labels = table.filter(pc.and_(measured, pc.invert(absent)))[label_column]
            else:
                labels = table[label_column]
            for entry in labels.value_counts().to_pylist():
                counts[classify_label(entry["values"])] += int(entry["counts"])
        out[donor] = DonorTally(
            donor_id=donor, n_rows=rows,
            n_called=counts[CALL], n_uncalled=counts[UNCALLED],
            n_label_unavailable=counts[LABEL_UNAVAILABLE], n_blank=counts[BLANK],
            n_no_nucleus_excluded=no_nucleus, n_nucleus_not_measured=not_measured,
            nucleus_rule_applied=nucleus_count is not None,
        )
    return out


# ---------------------------------------------------------------------------
# CLI: say what is current, and why, without building anything
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tissue", choices=list(TISSUES), action="append",
                        help="repeatable; defaults to both")
    parser.add_argument("--root", type=Path, action="append", dest="roots",
                        help=f"repeatable; overrides ${ROOTS_ENV} and the defaults")
    args = parser.parse_args(argv)
    report = {}
    for tissue in (args.tissue or list(TISSUES)):
        found = resolve(tissue, roots=args.roots)
        report[tissue] = found.describe()
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
