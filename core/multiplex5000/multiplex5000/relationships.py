"""Recorded relationships between sections, blocks and slides -- read, never inferred.

Which pieces of tissue share a section, a block or a slide decides which comparisons are allowed:
a same-section PhenoCycler/Xenium pair supports per-cell measurement, a serial pair only regional
correspondence, a same-block lymph node a block-level pairing with its pancreas, and a
different-block node on the same slide nothing anatomical at all. None of it can be read off two
annotations appearing in one image. It is what the maintainer recorded, and nothing else.

Three recorded sources, read as they are:

* **The donor notes** (site source ``donor_notes``, ``ProcessingNotes.txt``): one line per donor --
  the pancreas region, its lymph nodes and their block, and the Xenium design -- then the
  post-Xenium IF notes. Only those two sections are parsed. The reader stops at the first line
  after them, so nothing else in the file is parsed, kept or copied. The repo-root copy
  (``donor_notes_copy``) must say the same thing, donor by donor.
* **pln_blocks.csv** (site source ``pln_blocks``, ``results/v3/pln_blocks.csv``): one row per pLN
  annotation on a pancreas-project image, with its QuPath UUID, its block as the notes record it
  and the measured gap to the Tissue annotation, which stays a measurement. The maintainer
  confirmed it on 2026-10-01.
* **section_modes.csv** (``core/phenocycler/data/integration``): the PhenoCycler/Xenium section
  pairs that are not the sequential default, read by ``integration.modes``.

A line the reader cannot read is an error, never a guess; free text is kept verbatim as a
caution. Where the sources disagree, the disagreement is reported and nothing is chosen. A
relationship nobody recorded is ``unknown``, with the reason.

MultiPlex5000 D2 (2026-10-20). Read-only: ``python -m multiplex5000.relationships`` prints the
transcription for review; stage 00 (D3) writes it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Optional

from . import roots

SAME_SECTION = "same_section"
SERIAL_SECTION = "serial_section"
SAME_BLOCK = "same_block"
DIFFERENT_BLOCK = "different_block"
SHARED_SLIDE = "shared_slide"
UNKNOWN = "unknown"
KINDS = (SAME_SECTION, SERIAL_SECTION, SAME_BLOCK, DIFFERENT_BLOCK, SHARED_SLIDE, UNKNOWN)

#: The Xenium design a donor line ends with.
SEQUENTIAL = "sequential"
SAME_SECTION_AFTER_PHENOCYCLER = "same_section_after_phenocycler"

IF_HEADER = "Post-Xenium IF Notes"
_DONOR_LINE = re.compile(r"^(?P<label>[A-Za-z]*\d+): (?P<body>.+)$")
_UNDELIMITED_DONOR_LINE = re.compile(r"^[A-Za-z]*\d+\s")
_IF_LINE = re.compile(r"^(?P<label>[A-Za-z]*\d+): post-Xenium IF (?P<note>.+)$")
_DESIGN_SUFFIX = {", sequential Xenium": SEQUENTIAL,
                  ", same-slide post-phenocycler Xenium": SAME_SECTION_AFTER_PHENOCYCLER}
_STRUCTURED = re.compile(
    r"^panc(?P<region>Head|Body|Tail|Other) "
    r"(?P<pln>no pLN|and \d+ pLNs? (?:same|different)-block(?: and \d+ pLNs? (?:same|different)-block)*)")
_PLN = re.compile(r"(?P<n>\d+) pLNs? (?P<block>same|different)-block")

PLN_BLOCKS_COLUMNS = ("image", "entry_id", "donor", "pln_annotation", "object_id", "pieces",
                      "area_mm2", "overlap_with_tissue_um2", "gap_to_tissue_um",
                      "measured_relation", "block_from_notes", "notes")
_BLOCK = {"same-block": SAME_BLOCK, "different-block": DIFFERENT_BLOCK}
_SECTION_MODE = {"sequential": SERIAL_SECTION, "same_slide": SAME_SECTION}


class RelationshipError(ValueError):
    """A recorded source this reader cannot read as written."""


@dataclass(frozen=True)
class Evidence:
    """Where a statement was read: the source, the file, a digest of what was read, and where."""

    source: str
    path: str
    sha256: str
    where: str


@dataclass(frozen=True)
class DonorNote:
    label: str                          # as the notes write it
    donor_id: str                       # through the cohort's declared donor aliases
    pancreas_region: Optional[str]      # head / body / tail / other; None = not stated
    same_block_pln: Optional[int]       # None = not stated
    different_block_pln: Optional[int]
    xenium_design: str
    pln_from_caution: bool              # the counts came from free text, not the structured clause
    cautions: tuple[str, ...]
    post_xenium_if: Optional[str]
    evidence: Evidence
    body: str                           # the line as written, for comparing the two copies


@dataclass(frozen=True)
class PlnAnnotation:
    donor_id: str
    label: str
    image_id: str
    entry_id: str
    annotation_index: int
    annotation_uuid: str
    pieces: int
    block: str                          # SAME_BLOCK or DIFFERENT_BLOCK, as the notes record it
    measured_relation: str              # a measurement, kept verbatim and never used as the block
    gap_to_tissue_um: float
    evidence: Evidence


@dataclass(frozen=True)
class Relationship:
    kind: str
    donor_id: str
    a: str
    b: str
    evidence: tuple[Evidence, ...]
    reason: str = ""


@dataclass(frozen=True)
class Conflict:
    donor_id: str
    what: str
    detail: str
    evidence: tuple[Evidence, ...] = ()


@dataclass
class Reading:
    sources: list[Evidence] = field(default_factory=list)
    donors: dict[str, DonorNote] = field(default_factory=dict)
    pln: list[PlnAnnotation] = field(default_factory=list)
    modes: list = field(default_factory=list)
    relationships: list[Relationship] = field(default_factory=list)
    conflicts: list[Conflict] = field(default_factory=list)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _alias(label: str, aliases: dict[str, str]) -> str:
    return aliases.get(label, label)


def read_donor_notes(path: Path, aliases: dict[str, str], *, source: str = "donor_notes"
                     ) -> tuple[dict[str, DonorNote], Evidence]:
    """The donor and post-Xenium IF sections of the notes, with any introductory text.

    The file is read line by line and reading stops at the first line that belongs to neither
    section; the digest is of the lines read, so it changes only when they do.
    """
    path = Path(path)
    read: list[str] = []
    donors: dict[str, dict] = {}
    if_notes: dict[str, tuple[str, int]] = {}
    with path.open(encoding="utf-8") as handle:
        section = "before_donors"
        for number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip()
            if section == "before_donors":
                if _DONOR_LINE.match(line):
                    section = "donors"
                elif _UNDELIMITED_DONOR_LINE.match(line):
                    raise RelationshipError(f"{path}:{number}: not a donor line: {line!r}")
                elif line == IF_HEADER:
                    raise RelationshipError(f"{path}:{number}: post-Xenium IF notes precede donors")
                else:
                    read.append(raw)
                    continue
            if section == "donors":
                if not line:
                    section = "between"
                    read.append(raw)
                    continue
                m = _DONOR_LINE.match(line)
                if not m:
                    raise RelationshipError(f"{path}:{number}: not a donor line: {line!r}")
                if m["label"] in donors:
                    raise RelationshipError(f"{path}:{number}: {m['label']} has a second line")
                donors[m["label"]] = {"body": m["body"].rstrip(), "line": number}
                read.append(raw)
            elif section == "between":
                if not line:
                    read.append(raw)
                    continue
                if line == IF_HEADER:
                    section = "if"
                    read.append(raw)
                    continue
                break                                   # neither section: stop reading
            else:                                       # the post-Xenium IF section
                if not line:
                    break                               # its end: stop reading
                m = _IF_LINE.match(line)
                if not m:
                    raise RelationshipError(f"{path}:{number}: not a post-Xenium IF line: {line!r}")
                if m["label"] not in donors:
                    raise RelationshipError(f"{path}:{number}: {m['label']} has an IF note but "
                                            f"no donor line")
                if_notes[m["label"]] = (m["note"], number)
                read.append(raw)
    if not donors:
        raise RelationshipError(f"{path}: no donor lines found")
    digest = _sha("".join(read).encode("utf-8"))
    evidence = Evidence(source, str(path), digest, f"lines 1-{len(read)}")
    out = {}
    for label, entry in donors.items():
        where = Evidence(source, str(path), digest, f"line {entry['line']}")
        note, _ = if_notes.get(label, (None, None))
        donor = _parse_donor(label, entry["body"], _alias(label, aliases), where, note, path)
        if donor.donor_id in out:
            raise RelationshipError(f"{path}: {label} and {out[donor.donor_id].label} are the "
                                    f"same donor ({donor.donor_id})")
        out[donor.donor_id] = donor
    return out, evidence


def _parse_donor(label, body, donor_id, evidence, if_note, path) -> DonorNote:
    design = head = None
    for suffix, name in _DESIGN_SUFFIX.items():
        if body.endswith(suffix):
            design, head = name, body[:-len(suffix)]
            break
    if design is None:
        raise RelationshipError(f"{path} {evidence.where}: {label} does not end with a Xenium "
                                f"design ({' or '.join(s[2:] for s in _DESIGN_SUFFIX)})")
    clauses = head.split("; ")
    first, rest = clauses[0], clauses[1:]
    region = same = different = None
    cautions: list[str] = []
    from_caution = False
    m = _STRUCTURED.match(first)
    if m:
        region = m["region"].lower()
        same = different = 0
        for pm in _PLN.finditer(m["pln"]):
            if pm["block"] == "same":
                same += int(pm["n"])
            else:
                different += int(pm["n"])
        tail = first[m.end():]
        if tail:
            if not tail.startswith(", "):
                raise RelationshipError(f"{path} {evidence.where}: {label}: cannot read "
                                        f"{tail!r} after {m.group(0)!r}")
            cautions.append(tail[2:])
    else:
        cautions.append(first)
    cautions.extend(rest)
    if m is None:
        stated = [(int(pm["n"]), pm["block"]) for c in cautions for pm in _PLN.finditer(c)]
        if stated:
            from_caution = True
            if any(b == "same" for _, b in stated):
                same = sum(n for n, b in stated if b == "same")
            if any(b == "different" for _, b in stated):
                different = sum(n for n, b in stated if b == "different")
    return DonorNote(label, donor_id, region, same, different, design, from_caution,
                     tuple(c for c in cautions if c), if_note, evidence, body)


def read_pln_blocks(path: Path, aliases: dict[str, str], *, source: str = "pln_blocks"
                    ) -> tuple[list[PlnAnnotation], Evidence]:
    path = Path(path)
    data = path.read_bytes()
    digest = _sha(data)
    rows = list(csv.reader(data.decode("utf-8").splitlines()))
    if not rows or tuple(rows[0]) != PLN_BLOCKS_COLUMNS:
        raise RelationshipError(f"{path}: expected the columns {', '.join(PLN_BLOCKS_COLUMNS)}")
    out: list[PlnAnnotation] = []
    seen: set[str] = set()
    for number, values in enumerate(rows[1:], start=2):
        if len(values) != len(PLN_BLOCKS_COLUMNS):
            raise RelationshipError(f"{path} row {number}: {len(values)} fields, expected "
                                    f"{len(PLN_BLOCKS_COLUMNS)}")
        row = dict(zip(PLN_BLOCKS_COLUMNS, values, strict=True))
        block = _BLOCK.get(row["block_from_notes"])
        if block is None:
            raise RelationshipError(f"{path} row {number}: block_from_notes is "
                                    f"{row['block_from_notes']!r}, not same-block or different-block")
        uuid = row["object_id"].strip()
        if not uuid or uuid in seen:
            raise RelationshipError(f"{path} row {number}: object_id {uuid!r} is empty or repeated")
        seen.add(uuid)
        try:
            pieces, index, gap = int(row["pieces"]), int(row["pln_annotation"]), float(
                row["gap_to_tissue_um"])
        except ValueError as err:
            raise RelationshipError(f"{path} row {number}: {err}") from err
        if pieces < 1:
            raise RelationshipError(f"{path} row {number}: pieces must be at least 1")
        out.append(PlnAnnotation(_alias(row["donor"], aliases), row["donor"], row["image"],
                                 row["entry_id"], index, uuid, pieces, block,
                                 row["measured_relation"], gap,
                                 Evidence(source, str(path), digest, f"row {number}")))
    return out, Evidence(source, str(path), digest, f"rows 2-{len(rows)}")


def read_section_modes(path: Path, aliases: dict[str, str]) -> tuple[list, Evidence]:
    """``section_modes.csv`` through the reader that already validates it, donors aliased."""
    from .integration import modes

    path = Path(path)
    loaded = modes.load_manifest(SimpleNamespace(section_modes_csv=path))
    entries = sorted(loaded.values(), key=lambda e: e.key)
    digest = _sha(path.read_bytes())
    out = [SimpleNamespace(donor_id=_alias(e.donor_id, aliases), label=e.donor_id, roi=e.roi,
                           xenium_serial=e.xenium_serial, mode=e.mode, evidence_text=e.evidence,
                           evidence=Evidence("section_modes", str(path), digest,
                                             f"({e.donor_id}, {e.roi}, {e.xenium_serial or '*'})"))
           for e in entries]
    return out, Evidence("section_modes", str(path), digest, f"{len(out)} declarations")


def _relationships(reading: Reading) -> list[Relationship]:
    out: list[Relationship] = []
    for donor in reading.donors.values():
        kind = SERIAL_SECTION if donor.xenium_design == SEQUENTIAL else SAME_SECTION
        out.append(Relationship(kind, donor.donor_id, "PhenoCycler pancreas section",
                                "Xenium section (the notes name no serial)", (donor.evidence,)))
    for ann in reading.pln:
        name = f"pLN annotation {ann.annotation_uuid} ({ann.pieces} piece(s))"
        out.append(Relationship(ann.block, ann.donor_id, "pancreas tissue", name, (ann.evidence,)))
        out.append(Relationship(SHARED_SLIDE, ann.donor_id, f"image {ann.image_id}", name,
                                (ann.evidence,), "the annotation is drawn on the pancreas image"))
    for mode in reading.modes:
        b = (f"Xenium section {mode.xenium_serial}" if mode.xenium_serial
             else "every Xenium section of this donor and ROI not declared on its own")
        out.append(Relationship(_SECTION_MODE[mode.mode], mode.donor_id,
                                f"PhenoCycler {mode.roi} section", b, (mode.evidence,)))
    return out


def _conflicts(reading: Reading, copy: Optional[dict[str, DonorNote]]) -> list[Conflict]:
    out: list[Conflict] = []
    if copy is not None:
        for donor_id in sorted(set(reading.donors) | set(copy)):
            a, b = reading.donors.get(donor_id), copy.get(donor_id)
            if a is None or b is None or a.body != b.body:
                out.append(Conflict(donor_id, "the two copies of the donor notes differ",
                                    f"{(a.body if a else '(no line)')!r} / "
                                    f"{(b.body if b else '(no line)')!r}",
                                    tuple(x.evidence for x in (a, b) if x)))
    counted: dict[tuple[str, str], int] = {}
    for ann in reading.pln:
        counted[(ann.donor_id, ann.block)] = counted.get((ann.donor_id, ann.block), 0) + ann.pieces
    for donor in reading.donors.values():
        for block, stated in ((SAME_BLOCK, donor.same_block_pln),
                              (DIFFERENT_BLOCK, donor.different_block_pln)):
            pieces = counted.pop((donor.donor_id, block), 0)
            if stated is not None and stated != pieces:
                out.append(Conflict(donor.donor_id, f"{block} pLN: the notes and pln_blocks.csv",
                                    f"the notes say {stated}; pln_blocks.csv has {pieces} piece(s)",
                                    (donor.evidence,)))
    for (donor_id, block), pieces in sorted(counted.items()):
        out.append(Conflict(donor_id, f"{block} pLN: pln_blocks.csv and the notes",
                            f"pln_blocks.csv has {pieces} piece(s); the notes have no count"))
    same_slide = {m.donor_id for m in reading.modes if m.mode == "same_slide"}
    for donor in reading.donors.values():
        declared = donor.donor_id in same_slide
        if (donor.xenium_design == SAME_SECTION_AFTER_PHENOCYCLER) != declared:
            out.append(Conflict(donor.donor_id, "Xenium design: the notes and section_modes.csv",
                                f"the notes say {donor.xenium_design}; section_modes.csv "
                                f"{'declares' if declared else 'does not declare'} a same_slide "
                                f"section", (donor.evidence,)))
    for mode in reading.modes:
        if mode.donor_id not in reading.donors:
            out.append(Conflict(mode.donor_id, "section_modes.csv names a donor the notes do not",
                                mode.label, (mode.evidence,)))
    return out


def _given(value, name: str):
    """``None`` means the site's source; ``False`` means do not read one; else the path given."""
    if value is None:
        return roots.source(name)
    return value or None


#: Disagreements between the real sources that are known and explained. Stage 00's gate blocks
#: on any other until someone has looked at it and added it here with the reason.
KNOWN_DISAGREEMENTS: dict[tuple[str, str], str] = {
    ("6457", "same_block pLN: the notes and pln_blocks.csv"):
        "6457 is centrally excluded and was never processed, so no pln_blocks row was made",
}


def load(*, notes: Optional[Path] = None, notes_copy: Optional[Path] = None,
         pln_blocks: Optional[Path] = None, section_modes: Optional[Path] = None,
         aliases: Optional[dict[str, str]] = None, cohort_donors=None) -> Reading:
    """Read every recorded source this site names, relate them, and list where they disagree.

    A source the site does not name, and no argument gives, is reported as unknown -- the
    reading goes on without it rather than inventing what it would have said. Pass ``False``
    for a source that must not be read at all.
    """
    if aliases is None or cohort_donors is None:
        from . import profiles
        cohort = profiles.load_all()[("cohort", "io60_phenocycler")]
        aliases = cohort.donor_aliases if aliases is None else aliases
        if cohort_donors is None:
            import json
            cohort_donors = json.loads((profiles.REPO / cohort.membership).read_text()).get(
                "expected_donors", [])
    notes = _given(notes, "donor_notes")
    notes_copy = _given(notes_copy, "donor_notes_copy")
    pln_blocks = _given(pln_blocks, "pln_blocks")
    section_modes = section_modes or (roots.INTEGRATION / "section_modes.csv")
    reading = Reading()
    copy = None
    if notes is not None:
        reading.donors, ev = read_donor_notes(notes, aliases)
        reading.sources.append(ev)
    if notes_copy is not None:
        copy, ev = read_donor_notes(notes_copy, aliases, source="donor_notes_copy")
        reading.sources.append(ev)
    if pln_blocks is not None:
        reading.pln, ev = read_pln_blocks(pln_blocks, aliases)
        reading.sources.append(ev)
    reading.modes, ev = read_section_modes(section_modes, aliases)
    reading.sources.append(ev)
    reading.relationships = _relationships(reading)
    reading.conflicts = _conflicts(reading, copy)
    for name, value in (("donor_notes", notes), ("pln_blocks", pln_blocks)):
        if value is None:
            reading.conflicts.append(Conflict("*", f"{name}: not named by this site",
                                              "every relationship it records is unknown here"))
    if notes is not None:
        for donor_id in sorted(set(map(str, cohort_donors)) - set(reading.donors)):
            reading.relationships.append(Relationship(
                UNKNOWN, donor_id, "PhenoCycler pancreas section",
                "its Xenium section and any lymph node", (),
                "the cohort lists this donor; the donor notes have no line for it"))
    return reading


def report(reading: Reading) -> str:
    lines = ["Recorded relationships (read-only; nothing is written)", "", "Sources read:"]
    for ev in reading.sources:
        lines.append(f"  {ev.source:<17} {ev.where:<16} sha256 {ev.sha256[:16]}  {ev.path}")
    lines += ["", f"Donors ({len(reading.donors)}): label -> id | region | pLN same/different "
                  f"block | Xenium | post-Xenium IF | cautions"]
    for d in reading.donors.values():
        counts = f"{d.same_block_pln if d.same_block_pln is not None else '?'}/" \
                 f"{d.different_block_pln if d.different_block_pln is not None else '?'}"
        flag = " (from free text)" if d.pln_from_caution else ""
        lines.append(f"  {d.label:>6} -> {d.donor_id:<5} | {d.pancreas_region or 'not stated':<10}"
                     f" | {counts}{flag} | {d.xenium_design} | {d.post_xenium_if or '-'}"
                     f" | {'; '.join(d.cautions) or '-'}")
    kinds = {k: sum(r.kind == k for r in reading.relationships) for k in KINDS}
    lines += ["", "Relationships: " + ", ".join(f"{k} {n}" for k, n in kinds.items())]
    for r in reading.relationships:
        why = f"  [{r.reason}]" if r.reason else ""
        where = ", ".join(f"{e.source} {e.where}" for e in r.evidence) or "no record"
        lines.append(f"  {r.kind:<15} {r.donor_id:<5} {r.a} <-> {r.b}  ({where}){why}")
    lines += ["", f"Disagreements and gaps ({len(reading.conflicts)}), reported, not resolved:"]
    for c in reading.conflicts:
        lines.append(f"  {c.donor_id:<5} {c.what}: {c.detail}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m multiplex5000.relationships",
                                     description=__doc__.split("\n\n")[0])
    parser.add_argument("--notes", type=Path, help="the donor notes (default: the site's)")
    parser.add_argument("--notes-copy", type=Path, help="their repo-root copy (default: the site's)")
    parser.add_argument("--pln-blocks", type=Path, help="pln_blocks.csv (default: the site's)")
    parser.add_argument("--section-modes", type=Path, help="section_modes.csv (default: the repo's)")
    args = parser.parse_args(argv)
    reading = load(notes=args.notes, notes_copy=args.notes_copy, pln_blocks=args.pln_blocks,
                   section_modes=args.section_modes)
    sys.stdout.write(report(reading))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
