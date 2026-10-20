#!/usr/bin/env python
"""Protein INS and GCG inside each Xenium cell's OWN polygon, against that cell's own RNA.

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/protein_rna.py --donors 6436 6442

**What changed, and why.** The earlier arm sampled a 4 um disc at a warped centroid and then
restricted to cells inside a QuPath ``Islet`` polygon. The maintainer's instruction of
2026-09-01 removes both from THIS comparison: it runs on the Xenium cell and nucleus
segmentations against the IF insulin and glucagon channels only.

    CORRECTED 2026-09-21 -- this paragraph used to end "the QuPath islet annotations are **not
    to be used on Xenium or on the post-Xenium IF**", and a later session read that as a
    repo-wide prohibition and backed out working code on the strength of it.  It is not one,
    and nothing outside this docstring ever said it was: the string appears in no other file,
    no doc and no commit message; ``run_register.py`` and ``run_hormone.py`` -- live, in this
    same directory -- read ``geo.islets`` and ``islet_residents`` throughout and produced
    ``hormone_panc_refined.json`` on the very same day; and the maintainer exported fresh
    ``Islet`` geojsons for ten donors on 2026-09-18, seventeen days later.  What the
    instruction scoped is the sentence above: what THIS script measures.  Describe what this
    script does; do not write a rule for the rest of the repo into its docstring.

That is a stronger design, not merely a different one. Every object here is Xenium's own --
``cell_boundaries.parquet`` and ``nucleus_boundaries.parquet`` are segmentations of the same
DAPI the registration was fitted on -- so a cell's protein and a cell's RNA are the same cell
by construction, and no annotation, no islet call and no third party's polygon enters the path.
A disc also straddles neighbours at 9 um spacing; a polygon does not.

**The QuPath project is read for CHANNEL NAMES and nothing else.** Which Opal carried INS and
which carried GCG is recorded in no image -- the qptiff XML names only the instrument's channels
-- so the assignment lives in ``PostXeniumIHC/data/<entry>/server.json``, where the maintainer's
rename writes it, with ``restain_sections.csv``'s five hand-recorded rows as the fallback. That
is a name, not an object: no annotation, no classification and no QuPath-derived geometry enters
this script. See :func:`channel_roles`.

**Scored against the subpanel, not the whole panel.** ``core/panels``' 18-gene pancreas
endocrine identity core is the authority, and most of the 5,101-gene panel is not expressed in
pancreas at all. An earlier version ranked every gene and reported the top 25; that list was
mostly noise, and reading a cell type off the top of it is what produced a wrong conclusion.

**No cell is labelled here, and no AUROC is computed.** An AUROC needs a confirmed positive
class, and this lane has none: nobody has verified that a cell with one detected transcript is a
cell of that type. An earlier version defined "positive = count > 0", scored an AUROC against
it, and reported the result as though the label were truth -- an invented ground truth, which
this repository's rules forbid. What is reported is Spearman rho between the per-cell IF
intensity and the per-cell transcript count, plus the median IF intensity at each count level.
Those are statements about two measurements agreeing; they assert nothing about what a cell is.

**Defining a positive is the maintainer's call, not this file's.** If a truth set is ever
supplied, a classification statistic can be computed against it here.

**A whole-section AUROC is the wrong instrument for a focal signal, and this file cannot fix
that on its own.** Donor 6436 is T1D: insulin survives in a minority of lobules and is absent
from the rest, which is what the disease is. Averaged over the section its INS channel looks
inert -- and it is not, it is bright and cellular in the islets that still have beta cells. Any
number in this file that is computed over all cells inherits that, so read it beside
``protein_rna_figures.py --window INS``, which shows one islet at a time.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Mapping

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import hormone as HZ                                             # noqa: E402
import xenio                                                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                         # noqa: E402
from run_register import OUT, XENIUM_PATHS, donor_of              # noqa: E402
# `discover` and `GEOJSON_DIR` are deliberately NOT imported. Enumerating from the QuPath GeoJSON
# directory is what skipped ten donors with a PASS gate in silence (f54e5ce); the import outlived
# the call and is removed so nothing reaches for it again.

#: Resolution the protein is read at, in microns per pixel. 0.5 um/px puts ~280 pixels inside a
#: 9 um cell and ~110 inside a 6 um nucleus, which is a mean over a cell rather than over a
#: handful of pixels. It is a real pyramid level of these qptiffs (level 1), so no resampling
#: happens and the box-mean origin rule stays exact.
PROTEIN_UM_PER_PX = 0.5

#: The subpanel this lane scores against: ``core/panels``' pancreas endocrine identity core,
#: read through the regenerable snapshot. Not the whole 5,101-gene panel -- most of it is not
#: expressed in pancreas, and ranking all of it produced a list that invited over-reading.
SUBPANEL = "identity_core/pancreas/01_pancreas_endocrine.csv"


class ProteinRnaError(RuntimeError):
    """The per-cell comparison cannot be made as asked -- never approximated."""


#: The QuPath project whose ``server.json`` files carry the CURRENT channel assignment.
#: ``restain_sections.csv`` hand-records five of the twenty images. The maintainer is renaming
#: the Opal channels in QuPath, and that rename is written into
#: ``PostXeniumIHC/data/<entry>/server.json`` -- so for every image beyond those five the
#: assignment exists there and nowhere else. The qptiff's own XML never learns it: read today,
#: all twenty images still name their channels ``DAPI / Opal 520 / Opal 570 / Opal 690 /
#: Sample AF`` in their ``<Name>`` tags.
QUPATH_PROJECT = IMAGES.parent

#: Image -> donor / serial for ALL twenty images. ``restain_sections.csv`` has five data rows,
#: which is why ``.iloc[0]`` on it raised ``IndexError`` for eleven of the thirteen donors whose
#: registration passes its gate.
REGISTRATION_INPUTS = HERE.parents[1] / "data" / "integration" / "restain_registration_inputs.csv"

#: The two markers this lane measures.
MARKERS = ("INS", "GCG")

#: Spellings accepted as a marker, after lowercasing and dropping non-alphanumerics. A rename
#: this file does not recognise is REPORTED under its own name rather than read as absence:
#: entry 15 (``0037502_Scan1.qptiff``, donor 6380, a PASS donor) spells glucagon ``CGC``, and a
#: name-matching reader that merely failed to match it would drop a registered donor for a
#: transposition nobody had been told about.
MARKER_ALIASES: dict[str, tuple[str, ...]] = {
    "INS": ("ins", "insulin"),
    "GCG": ("gcg", "glucagon"),
}

#: A channel still carrying the instrument's own name has not been assigned a marker yet.
INSTRUMENT_NAME = re.compile(r"^(dapi|opal\d+|sampleaf|autofluorescence|af)$")

#: Named states a resolution can be unavailable in. Each is a fact about the QuPath pass or the
#: inventory, never about the donor's biology, and each is emitted by name so ``measure_one``
#: skips that donor out loud instead of raising.
CHANNELS_NOT_YET_ASSIGNED = "channels_not_yet_assigned"
CHANNEL_NAME_UNRECOGNISED = "channel_name_unrecognised"


def _norm(name) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


@dataclass(frozen=True)
class ChannelRoles:
    """Which page of one qptiff carries which marker, and at what pixel size.

    ``channels`` maps marker -> PAGE INDEX, resolved by NAME. Never by position: the INS/GCG
    index differs between images in this cohort -- measured 2026-09-18 over the twenty QuPath
    entries, INS/GCG sit at pages (2,3) on six entries, (1,2) on six, (1,3) on two and (3,1) on
    donor 6579 -- so a reader that took a fixed index would return another antibody's channel
    on most of the cohort with no error at all.
    """

    serial: str
    donor_id: str
    image_file: str
    source: str
    channels: dict[str, int] = field(default_factory=dict)
    names: dict[str, str] = field(default_factory=dict)
    qptiff_names: tuple[str, ...] = ()
    qupath_entry: str | None = None
    pixel_size_um: float | None = None
    pixel_size_source: str = ""
    unavailable: str | None = None
    detail: str = ""
    missing_markers: tuple[str, ...] = ()
    unrecognised_names: tuple[str, ...] = ()
    csv_conflict: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.unavailable is None

    def as_record(self) -> dict:
        return {"donor_id": self.donor_id, "xenium_serial": self.serial,
                "image_file": self.image_file, "qupath_entry": self.qupath_entry,
                "channel_source": self.source, "channels": dict(self.channels),
                "channel_names": dict(self.names),
                "qptiff_channel_names": list(self.qptiff_names),
                "pixel_size_um": self.pixel_size_um,
                "pixel_size_source": self.pixel_size_source,
                "unavailable": self.unavailable, "detail": self.detail,
                "missing_markers": list(self.missing_markers),
                "unrecognised_names": list(self.unrecognised_names),
                "csv_conflict": list(self.csv_conflict)}


def unavailable(reason: str, detail: str, **kw) -> ChannelRoles:
    """A NAMED absence. ``measure_one`` skips on it; nothing downstream guesses past it."""
    kw.setdefault("serial", "")
    kw.setdefault("donor_id", "")
    kw.setdefault("image_file", "")
    kw.setdefault("source", "")
    return ChannelRoles(unavailable=reason, detail=detail, **kw)


def qptiff_facts(qptiff: Path) -> dict:
    """Channel names in page order and ``<PixelSizeMicrons>``, from the image's own XML.

    **The pixel size comes from the image, not from the CSV.** Measured 2026-09-18 over all
    twenty post-Xenium IF images: every 20x image's XML reads 0.25 um/px and both 10x images
    (``HDL115_Scan1.qptiff``, ``npod6591_Scan1.qptiff``) read 0.5, with no exceptions.
    ``restain_sections.csv`` declares 0.2505 for all five of its rows -- within 0.2 % of the XML
    for the three 20x ones and HALF the XML's value for the two 10x ones. 6591 passes its
    registration gate, so reading its image at 0.2505 selects the wrong pyramid level and puts
    every measurement on a grid at twice the intended scale.

    The TIFF resolution rational is deliberately not read: both carry numerator 4294967295
    (uint32 max, a saturated value), so the trailing digits of 0.25047 / 0.50139 are an artifact
    of that clamp rather than a measurement -- the same reading ``registration_inputs.py`` makes.
    """
    import tifffile

    qptiff = Path(qptiff)
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        level = list(getattr(series, "levels", [series]))[0]
        names, pixels = [], []
        for page in level.pages:
            description = page.description or ""
            match = re.search(r"<Name>([^<]*)</Name>", description)
            names.append(match.group(1).strip() if match else "")
            pixels.extend(float(v) for v in re.findall(
                r"<PixelSizeMicrons>\s*([^<]+?)\s*</PixelSizeMicrons>", description))
    distinct = sorted({p for p in pixels if np.isfinite(p) and p > 0})
    return {"channel_names": tuple(names),
            "pixel_size_um": distinct[0] if len(distinct) == 1 else None,
            "pixel_size_detail": ("qptiff_xml_PixelSizeMicrons" if len(distinct) == 1
                                  else f"expected one positive PixelSizeMicrons, found {distinct}")}


def qupath_entries(project_dir: Path = QUPATH_PROJECT) -> dict[str, list[dict]]:
    """image file name -> the QuPath entries that hold it, with their CURRENT channel names.

    Keyed by the basename of each ``server.json``'s ``uri`` rather than by its ``metadata.name``:
    the metadata name carries a donor token and a ``- resolution #1`` suffix that the inventory
    CSV does not, and matching on a decorated string is how a rename becomes a silent miss.
    """
    project_dir = Path(project_dir)
    out: dict[str, list[dict]] = {}
    data_dir = project_dir / "data"
    if not data_dir.is_dir():
        return out
    for entry_dir in sorted(data_dir.iterdir(), key=lambda p: p.name):
        server = entry_dir / "server.json"
        if not server.is_file():
            continue
        try:
            j = json.loads(server.read_text())
        except (OSError, ValueError):
            continue
        uri = str(j.get("uri") or "")
        image_file = uri.rsplit("/", 1)[-1]
        if not image_file:
            continue
        metadata = j.get("metadata") or {}
        names = tuple(str(c.get("name") or "") for c in (metadata.get("channels") or []))
        out.setdefault(image_file, []).append(
            {"entry": entry_dir.name, "server_json": str(server), "channel_names": names})
    return out


def _match_by_name(names: Iterable[str]) -> tuple[dict[str, int], dict[str, str], list[str], list[str]]:
    """``(marker -> page, marker -> name as spelled, missing markers, unrecognised names)``."""
    names = list(names)
    found: dict[str, list[int]] = {}
    for i, name in enumerate(names):
        key = _norm(name)
        for marker, aliases in MARKER_ALIASES.items():
            if key in aliases:
                found.setdefault(marker, []).append(i)
    ambiguous = sorted(m for m, pages in found.items() if len(pages) > 1)
    if ambiguous:
        raise ProteinRnaError(
            f"channel roles are ambiguous: {ambiguous} each name more than one channel in "
            f"{names}. A hand-made assignment with two channels for one marker is refused, "
            "never picked.")
    channels = {m: pages[0] for m, pages in found.items()}
    spelled = {m: names[p] for m, p in channels.items()}
    missing = [m for m in MARKERS if m not in channels]
    known = {a for aliases in MARKER_ALIASES.values() for a in aliases}
    unrecognised = [name for name in names
                    if name and _norm(name) not in known
                    and not INSTRUMENT_NAME.match(_norm(name))]
    return channels, spelled, missing, unrecognised


def _csv_channels(row: pd.Series, qptiff_names: Iterable[str]) -> tuple[dict[str, int], dict[str, str]]:
    """The five-row CSV's assignment, resolved against the image's own ``<Name>`` XML.

    Same rule as the QuPath path and as ``hormone.channel_index``: the CSV names an INSTRUMENT
    channel (``channel_opal690``) and the page is found by matching that name in the file, never
    by the column's position in the header.
    """
    qptiff_names = list(qptiff_names)
    pages = {_norm(name): i for i, name in enumerate(qptiff_names)}
    channels: dict[str, int] = {}
    spelled: dict[str, str] = {}
    for column, value in row.items():
        if not (isinstance(column, str) and column.startswith("channel_")):
            continue
        marker = str(value).strip().upper()
        if marker not in MARKERS or marker in channels:
            continue
        page = pages.get(column[len("channel_"):])
        if page is None:
            continue
        channels[marker] = page
        spelled[marker] = qptiff_names[page]
    return channels, spelled


def channel_roles(serial: str, *, sections: pd.DataFrame | None = None,
                  project_dir: Path = QUPATH_PROJECT,
                  inputs_csv: Path = REGISTRATION_INPUTS,
                  images_dir: Path = IMAGES,
                  facts_reader: Callable[[Path], dict] = qptiff_facts,
                  expect_donor: str | None = None) -> ChannelRoles:
    """Which page carries INS and which carries GCG, for one Xenium serial.

    Resolution order, and why it is that way round:

    1. **The QuPath project.** The maintainer is assigning the roles there now, one image at a
       time, and that is the only place fifteen of the twenty images carry them. Matched BY
       NAME against the entry's ``server.json`` channel list.
    2. **``restain_sections.csv``**, for an image QuPath has not been through. It covers five
       images and is hand-recorded; it is a fallback, not the authority, and it is the file that
       ``.iloc[0]`` used to index into for donors it has no row for.
    3. Otherwise a NAMED unavailable, so the donor is skipped with a printed reason rather than
       raising ``IndexError`` eleven times out of thirteen.

    Where both sources resolve and DISAGREE, the QuPath assignment is used and the disagreement
    is recorded in ``csv_conflict`` and printed -- it is not averaged, not silently preferred,
    and not swallowed. Measured 2026-09-18, the two 10x images disagree: QuPath names ``Opal
    570`` INS and ``Opal 690`` GCG on both, while the CSV records ``Opal 690`` INS, ``Opal 520``
    GCG and ``Opal 570`` FAILED for every one of its five rows.
    """
    inventory = pd.read_csv(inputs_csv, comment="#", dtype=str)
    hits = inventory[inventory["xenium_serial"].astype(str).str.strip() == str(serial).strip()]
    if len(hits) != 1:
        return unavailable(
            "serial_not_in_registration_inputs",
            f"{serial} has {len(hits)} rows in {inputs_csv}; expected one",
            serial=str(serial), donor_id=str(expect_donor or ""))
    row = hits.iloc[0]
    image_file = str(row["image_file"]).strip()
    donor = str(row["donor_id"]).strip()
    base = dict(serial=str(serial), donor_id=donor, image_file=image_file)
    if expect_donor is not None and str(expect_donor).strip() != donor:
        return unavailable(
            "donor_disagrees_with_registration_inputs",
            f"serial {serial}: the registration export says donor {expect_donor}, "
            f"{inputs_csv.name} says {donor}", **base)

    qptiff = Path(images_dir) / image_file
    qptiff_names: tuple[str, ...] = ()
    pixel_um: float | None = None
    pixel_source = ""
    if qptiff.is_file():
        facts = facts_reader(qptiff)
        qptiff_names = tuple(facts["channel_names"])
        pixel_um = facts["pixel_size_um"]
        pixel_source = facts["pixel_size_detail"]

    csv_row = None
    if sections is not None and len(sections):
        csv_hits = sections[sections["xenium_serial"].astype(str).str.strip() == str(serial).strip()]
        if len(csv_hits) == 1:
            csv_row = csv_hits.iloc[0]
    if pixel_um is None and csv_row is not None:
        declared = str(csv_row.get("pixel_size_um", "")).strip()
        try:
            candidate = float(declared)
        except ValueError:
            candidate = float("nan")
        if np.isfinite(candidate) and candidate > 0:
            pixel_um, pixel_source = candidate, "restain_sections.csv (image unreadable)"

    entries = qupath_entries(project_dir).get(image_file, [])
    if len(entries) > 1:
        return unavailable(
            "duplicate_qupath_entries",
            f"{image_file} is held by {len(entries)} QuPath entries "
            f"({[e['entry'] for e in entries]}); an ambiguous project is refused, never picked",
            qptiff_names=qptiff_names, **base)

    entry = entries[0] if entries else None
    q_channels: dict[str, int] = {}
    q_names: dict[str, str] = {}
    missing: list[str] = list(MARKERS)
    unrecognised: list[str] = []
    if entry is not None:
        q_channels, q_names, missing, unrecognised = _match_by_name(entry["channel_names"])

    csv_channels: dict[str, int] = {}
    csv_names: dict[str, str] = {}
    if csv_row is not None and qptiff_names:
        csv_channels, csv_names = _csv_channels(csv_row, qptiff_names)

    # Recorded, never reconciled. Two hand-made authorities that disagree about which
    # antibody a channel carries is a fact for the maintainer, not an average to take.
    conflict = tuple(
        f"{m}: QuPath names page {q_channels[m]} ({q_names[m]}), "
        f"restain_sections.csv names page {csv_channels[m]} ({csv_names[m]})"
        for m in MARKERS if m in q_channels and m in csv_channels
        and q_channels[m] != csv_channels[m])

    if not missing:
        channels, names, source = q_channels, q_names, "qupath_project"
    elif all(m in csv_channels for m in MARKERS):
        channels, names, source = csv_channels, csv_names, "restain_sections.csv"
    else:
        reason = CHANNEL_NAME_UNRECOGNISED if unrecognised else CHANNELS_NOT_YET_ASSIGNED
        detail = (f"{image_file}: QuPath entry "
                  f"{entry['entry'] if entry else '(none)'} names "
                  f"{list(entry['channel_names']) if entry else []}; "
                  f"{', '.join(missing)} not assigned")
        if unrecognised:
            detail += (f"; unrecognised channel name(s) {unrecognised} -- a rename this reader "
                       f"does not know, not an absence")
        return unavailable(reason, detail, qupath_entry=entry["entry"] if entry else None,
                           qptiff_names=qptiff_names, missing_markers=tuple(missing),
                           unrecognised_names=tuple(unrecognised),
                           pixel_size_um=pixel_um, pixel_size_source=pixel_source, **base)

    if pixel_um is None:
        return unavailable(
            "pixel_size_unavailable",
            f"{image_file}: {pixel_source or 'the image is not on this machine'} and no usable "
            f"pixel_size_um is recorded for {serial}. A pixel size is never guessed: it chooses "
            f"the pyramid level every measurement is read from.",
            qupath_entry=entry["entry"] if entry else None, qptiff_names=qptiff_names,
            csv_conflict=conflict, **base)

    return ChannelRoles(source=source, channels=channels, names=names,
                        qptiff_names=qptiff_names,
                        qupath_entry=entry["entry"] if entry else None,
                        pixel_size_um=float(pixel_um), pixel_size_source=pixel_source,
                        csv_conflict=conflict,
                        unrecognised_names=tuple(unrecognised),
                        serial=str(serial), donor_id=donor, image_file=image_file)


def resolve_all(serials: Mapping[str, str], **kw) -> dict[str, ChannelRoles]:
    """``serial -> ChannelRoles`` for every registered serial. Never raises on one bad row."""
    out: dict[str, ChannelRoles] = {}
    for serial, donor in sorted(serials.items()):
        try:
            out[serial] = channel_roles(serial, expect_donor=donor, **kw)
        except (ProteinRnaError, OSError, ValueError) as exc:
            out[serial] = unavailable("channel_resolution_failed", f"{type(exc).__name__}: {exc}",
                                      serial=str(serial), donor_id=str(donor))
    return out


def channel_assignment_report(resolved: Mapping[str, ChannelRoles], *, roi: str = "panc",
                              project_dir: Path = QUPATH_PROJECT) -> dict:
    """What the QuPath renaming pass has and has not reached, as a build artifact.

    Emitted every run so the not-yet-assigned list is a file somebody can act on rather than a
    line that scrolled past. The unavailable donors are grouped BY NAMED REASON, because
    ``channels_not_yet_assigned`` (nobody has been through that image) and
    ``channel_name_unrecognised`` (somebody has, and spelled it ``CGC``) need different actions.
    """
    resolved_rows, blocked = [], []
    for serial in sorted(resolved):
        roles = resolved[serial]
        (resolved_rows if roles.ok else blocked).append(roles.as_record())
    by_reason: dict[str, list[str]] = {}
    for record in blocked:
        by_reason.setdefault(str(record["unavailable"]), []).append(str(record["donor_id"]))
    return {"roi": roi, "qupath_project": str(project_dir),
            "registration_inputs": str(REGISTRATION_INPUTS),
            "n_registered": len(resolved), "n_resolved": len(resolved_rows),
            "n_unavailable": len(blocked),
            "resolved": resolved_rows, "unavailable": blocked,
            "donors_by_reason": {k: sorted(v) for k, v in sorted(by_reason.items())},
            "note": ("channel roles are matched BY NAME, never by index -- the INS/GCG page "
                     "differs between images in this cohort; pixel size is read from each "
                     "image's own XML, never from the sections CSV, which declares one value "
                     "for images that differ by a factor of two")}


def write_channel_assignment(report: dict, out_dir: Path = OUT) -> Path:
    path = Path(out_dir) / "protein_rna" / "channel_assignment.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False))
    return path


def warp_polygons(polys, transform):
    """Xenium-micron polygons -> IF-micron polygons.

    ``transform.apply`` is the forward direction: this lane fits Xenium -> IF, and the cells
    being warped are Xenium's. Getting this backwards does not error, it puts every polygon
    ~2 x (714, 282) um away on 6436 and every protein value becomes a sample of empty slide.
    """
    import shapely

    return [shapely.transform(g, lambda a: transform.apply(a)) for g in polys]


def pyramid_level_for(qptiff: Path, um_per_px: float, base_um_per_px: float) -> tuple[int, int]:
    """The level whose scale is exactly ``um_per_px``, and its integer factor.

    Exact rather than nearest. A level that is merely close would put the label raster and the
    image on different grids, and the resulting error grows across the slide instead of showing
    up as an offset.
    """
    import export_dapi as ED

    import tifffile
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        shapes = [(int(lv.shape[-2]), int(lv.shape[-1]))
                  for lv in getattr(series, "levels", [series])]
    want = int(round(um_per_px / base_um_per_px))
    for lvl in range(len(shapes)):
        f = ED.pyramid_factor(shapes, lvl)
        if f == want:
            return lvl, f
    raise ProteinRnaError(
        f"{qptiff.name} has no pyramid level at {um_per_px} um/px "
        f"(base {base_um_per_px}, factors {[ED.pyramid_factor(shapes, l) for l in range(len(shapes))]})")


def mean_in_labels(image: np.ndarray, labels: np.ndarray, n: int) -> np.ndarray:
    """Mean of ``image`` over each label 1..n, NaN where a label drew no pixel.

    NaN rather than 0: a cell that fell outside the imaged area has an unavailable protein
    value, not a dark one, and a 0 would rank it at the bottom of every statistic below.
    """
    from scipy import ndimage

    counts = np.bincount(labels.ravel(), minlength=n + 1)[1:]
    sums = ndimage.sum_labels(image, labels, index=np.arange(1, n + 1))
    out = np.full(n, np.nan, dtype=np.float64)
    hit = counts > 0
    out[hit] = sums[hit] / counts[hit]
    return out


def gene_rank_association(counts, values: np.ndarray, gene_cols: list[int]) -> dict:
    """Association between two MEASUREMENTS on the same cell. No labels, no truth set.

    **Why this is not an AUROC.** An AUROC scores a ranking against a positive class, and this
    lane has no confirmed positive class: nobody has verified that a cell with one detected
    transcript of a gene is a cell of that type. An earlier version of this file defined
    "positive = count > 0", computed an AUROC against it, and reported the result as if the
    label were ground truth. It was a definition invented here, and the repository's own rule
    is that existing classifications are audit context and never a gate.

    So what is reported instead needs no class at all: Spearman rho between the per-cell IF
    intensity and the per-cell transcript count for one gene, plus the median IF intensity at
    each count level. Both are statements about two measurements agreeing or not. Neither
    asserts that any cell is anything, and the count-stratified medians are there so the shape
    can be read directly rather than through a coefficient.

    Cells with no IF value are dropped, never zero-filled: nothing was measured for them.
    """
    import scipy.sparse as sp
    from scipy.stats import rankdata

    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    if ok.sum() < 1000:
        raise ProteinRnaError(f"only {int(ok.sum())} cells carry a protein value")
    r_if = rankdata(v[ok])
    X = sp.csc_matrix(counts)[ok][:, gene_cols].toarray()

    out = {}
    for j, col in enumerate(gene_cols):
        c = X[:, j].astype(np.float64)
        r_c = rankdata(c)
        rho = float(np.corrcoef(r_if, r_c)[0, 1]) if c.std() > 0 else float("nan")
        strata = {}
        for k in (0, 1, 2, 3):
            m = (c >= k) if k == 3 else (c == k)
            key = "3+" if k == 3 else str(k)
            strata[key] = {"n": int(m.sum()),
                           "median_if": float(np.median(v[ok][m])) if m.any() else None}
        out[col] = {"spearman_rho": rho, "n_cells": int(ok.sum()),
                    "n_cells_detected": int((c > 0).sum()), "median_if_by_count": strata}
    return out





def rank_contrast(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``rank(a) - rank(b)`` per cell, on the cells where both were measured.

    **This is the statistic that separates alpha from beta, and neither channel alone can.**
    Insulin and glucagon are both confined to islets, so each on its own ranks *any* islet gene
    near the top -- which is why 6442's GCG channel returns a beta gene list. Subtracting the
    two within one cell removes everything the two channels share, which is islet membership,
    and leaves only what differs between them, which is subtype.

    Percentile ranks rather than raw values because the two channels are different fluorophores
    at different exposures on an 8-bit scanner: their units are not comparable and a difference
    of raw intensities would be dominated by whichever was brighter. A cell missing either side
    is NaN, never zero -- zero is the middle of this scale and would read as "equally both".
    """
    from scipy.stats import rankdata

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    out = np.full(a.shape, np.nan)
    if both.sum() < 2:
        return out
    n = float(both.sum())
    out[both] = rankdata(a[both]) / n - rankdata(b[both]) / n
    return out


def score_against_genes(values: np.ndarray, dd, by_hint: dict) -> dict:
    """One per-cell IF vector against the SUBPANEL's endocrine genes, grouped by subtype hint.

    Restricted to ``core/panels``' 18-gene pancreas endocrine identity core, which is the panel
    taxonomy's one authority. Ranking all 5,101 genes was a mistake and is gone: most of the
    panel is not expressed in pancreas, so a ranked list of everything is mostly noise dressed
    as a result.

    ``subtype_hint`` is carried through from the subpanel and is a HINT, in its own file's
    words. It is not used to label any cell here, and the per-hint means below are means of
    association statistics -- not a classification score, and not evidence that a cell is a
    beta cell.
    """
    gi = {g: i for i, g in enumerate(dd.genes)}
    wanted = [(g, hint or "pan") for hint, genes in by_hint.items() for g in genes if g in gi]
    cols = [gi[g] for g, _ in wanted]
    assoc = gene_rank_association(dd.counts, values, cols)

    out = {"n_cells_measured": int(np.isfinite(values).sum()),
           "statistic": "spearman_rho(IF intensity, transcript count) per cell",
           "no_ground_truth": ("no cell is labelled here; there is no confirmed positive class "
                               "for cell type in this lane and none is invented"),
           "genes": {}}
    for (g, hint), col in zip(wanted, cols):
        out["genes"][g] = {"subtype_hint": hint, **assoc[col]}
    for hint in ("beta", "alpha"):
        vals = [v["spearman_rho"] for v in out["genes"].values()
                if v["subtype_hint"] == hint and np.isfinite(v["spearman_rho"])]
        out[f"mean_rho_{hint}"] = float(np.mean(vals)) if vals else None
    return out




def polygons_per_cell(bundle: Path, compartment: str, cell_ids) -> list:
    """One geometry per Xenium cell, in Xenium microns, or ``None`` where there is none.

    Neither boundary table is one row per cell, and the two are not irregular in the same way:

    * ``cell_boundaries.parquet`` can carry more than one ring for a cell;
    * ``nucleus_boundaries.parquet`` is per **nucleus** -- a cell may have two, or none.

    Parts belonging to one cell are combined and drawn at a single label, so a two-nucleus cell
    contributes both nuclei to that cell's mean rather than whichever one a naive reindex
    happened to keep. A cell with no nucleus is ``None`` and ends up NaN, which is the honest
    value: nothing was measured for it, and it must not rank as dark.
    """
    import shapely

    from phenocycler.integration.boundaries import (load_xenium_nucleus_polygons,
                                                    load_xenium_polygons)

    if compartment == "cell":
        series = load_xenium_polygons(bundle / "cell_boundaries.parquet", cell_ids=cell_ids)
        owner = np.asarray(series.index, dtype=object)
        geoms = list(series.to_numpy())
    elif compartment == "nucleus":
        series, owner = load_xenium_nucleus_polygons(bundle / "nucleus_boundaries.parquet")
        geoms = list(series.to_numpy())
        owner = np.asarray(owner, dtype=object)
        want = set(map(str, cell_ids))
        keep = np.array([str(o) in want for o in owner], dtype=bool)
        geoms = [g for g, k in zip(geoms, keep) if k]
        owner = owner[keep]
    else:
        raise ProteinRnaError(f"unknown compartment {compartment!r}")

    by: dict[str, list] = {}
    for o, g in zip(owner, geoms):
        if g is None or g.is_empty:
            continue
        by.setdefault(str(o), []).append(g)

    out = []
    for cid in map(str, cell_ids):
        parts = by.get(cid)
        if not parts:
            out.append(None)
        elif len(parts) == 1:
            out.append(parts[0])
        else:
            # A MultiPolygon, deliberately not a union. `rasterize_labels` draws every part of
            # a multipart geometry at the same label, so the union buys nothing -- and
            # `union_all` raises TopologyException on real Xenium rings ("unable to assign free
            # hole to a shell"), which would drop a cell for a reason that does not matter here.
            flat = []
            for g in parts:
                flat.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
            out.append(shapely.MultiPolygon(flat))
    return out


def measure_one(serial: str, sections: pd.DataFrame, *, roi: str = "panc",
                out_dir: Path = OUT, um_per_px: float = PROTEIN_UM_PER_PX,
                roles: ChannelRoles | None = None,
                transform_stem: str = "refined") -> dict:
    """One section measured. ``roles`` is resolved here when the caller has not resolved it.

    The channel assignment and the pixel size used to come from ``restain_sections.csv`` through
    ``sections[sections.xenium_serial == serial].iloc[0]``. That CSV has five data rows and
    ``registered_serials()`` returns thirteen, so eleven of thirteen donors raised ``IndexError``
    on that one line. Both facts now come from :func:`channel_roles`.
    """
    import tifffile
    import zarr

    from phenocycler.integration.boundaries import load_xenium_polygons
    from phenocycler.integration.rasterize import GridSpec
    from rasterize_labels import rasterize_labels
    import validate as VAL
    from run_align import _saved_transform
    from run_register import SUBPANELS

    donor = donor_of(serial, XENIUM_PATHS)
    if roles is None:
        roles = channel_roles(serial, sections=sections, expect_donor=donor)
    if not roles.ok:
        # Refused before any work is done. A named absence is a fact about the QuPath pass, and
        # it belongs at the top of this function rather than halfway through a section.
        raise ProteinRnaError(f"{donor}/{serial}: {roles.unavailable}: {roles.detail}")
    _core, by_hint = VAL.endocrine_core(SUBPANELS)
    qptiff = IMAGES / roles.image_file
    if not qptiff.exists():
        raise ProteinRnaError(f"{qptiff} is absent")
    base_um = float(roles.pixel_size_um)
    transform = _saved_transform(donor, transform_stem, out_dir)

    dd = xenio.load_donor(donor, XENIUM_PATHS, roi=roi, prefer="bundle")
    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    print(f"[protein-rna] {donor}: {len(dd.cell_id):,} cells, {len(dd.genes):,} genes "
          f"({dd.source})", flush=True)

    lvl, factor = pyramid_level_for(qptiff, um_per_px, base_um)
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        H = int(series.levels[lvl].shape[-2]); W = int(series.levels[lvl].shape[-1])
    # The grid IS the image: same origin, same pixel size, same extent, so a label index and a
    # pixel index are the same thing and no second alignment can creep in between them.
    grid = GridSpec(origin_x_um=0.5 * um_per_px, origin_y_um=0.5 * um_per_px,
                    pixel_um=um_per_px, height=H, width=W)

    rec = {"donor_id": donor, "xenium_serial": serial, "roi": roi,
           "qptiff": str(qptiff), "xenium_bundle": str(bundle),
           "transform_source": transform_stem,
           "alignment_gate": str(out_dir / "valis" / f"donor_id={donor}"
                                 / f"alignment_{transform_stem}.json"),
           "protein_um_per_px": um_per_px,
           "pyramid_level": lvl, "pyramid_factor": factor,
           # Recorded PER ROW, because it is not one value across the cohort: the 20x images
           # read 0.25 um/px and the two 10x ones read 0.5, and the CSV declares 0.2505 for all.
           "image_pixel_size_um": base_um, "pixel_size_source": roles.pixel_size_source,
           "channel_source": roles.source, "qupath_entry": roles.qupath_entry,
           "n_cells": int(len(dd.cell_id)), "n_genes": int(len(dd.genes)),
           "uses_qupath_annotations": False, "subpanel": SUBPANEL,
           "note": ("no QuPath annotation of any class is read by THIS lane; every object is "
                    "Xenium's own segmentation, per the maintainer's instruction of "
                    "2026-09-01 scoping this comparison. The QuPath project is read for "
                    "CHANNEL NAMES only. This records what this lane does; it is not a "
                    "statement about what other lanes may read"),
           "compartments": {}}
    if roles.csv_conflict:
        rec["channel_assignment_conflict"] = list(roles.csv_conflict)
        print(f"[protein-rna] {donor}: channel assignment DISAGREES between sources -- "
              + "; ".join(roles.csv_conflict)
              + f". Used {roles.source}.", flush=True)
    names = roles.qptiff_names or tuple(HZ.channel_names(qptiff))
    channels = {marker: (roles.channels[marker],
                         names[roles.channels[marker]] if roles.channels[marker] < len(names)
                         else roles.names[marker])
                for marker in MARKERS}
    print(f"[protein-rna] {donor}: channels from {roles.source} -- "
          + ", ".join(f"{m}=page {p} ({roles.names[m]!r} on {names[p] if p < len(names) else '?'})"
                      for m, (p, _n) in sorted(channels.items()))
          + f"; pixel size {base_um} um/px ({roles.pixel_size_source})", flush=True)

    for compartment in ("cell", "nucleus"):
        geoms = polygons_per_cell(bundle, compartment, dd.cell_id)
        keep = np.array([g is not None for g in geoms])
        missing = int((~keep).sum())
        idx = np.nonzero(keep)[0]
        warped = warp_polygons([geoms[i] for i in idx], transform)
        lab = rasterize_labels(warped, grid, order="input")
        print(f"[protein-rna] {donor}: {compartment} {len(idx):,} polygons rasterised "
              f"({missing:,} absent)", flush=True)

        comp = {"n_polygons": int(len(idx)), "n_absent": missing, "markers": {}}
        measured: dict[str, np.ndarray] = {}
        with tifffile.TiffFile(str(qptiff)) as tif:
            store = tif.series[0].aszarr(level=lvl)
            try:
                z = zarr.open(store, mode="r")
                for marker, (page, _chname) in channels.items():
                    img = np.asarray(z[page]).astype(np.float32)
                    vals_sub = mean_in_labels(img, lab.labels, len(idx))
                    del img
                    vals = np.full(len(dd.cell_id), np.nan)
                    vals[idx] = vals_sub
                    measured[marker] = vals
            finally:
                store.close()

        scored = dict(measured)
        scored["INS_minus_GCG"] = rank_contrast(measured["INS"], measured["GCG"])

        for name, vals in scored.items():
            # The per-cell vector is kept because everything downstream -- cell groups, figures,
            # any re-scoring -- must use the SAME numbers these AUROCs were computed from.
            vp = out_dir / "protein_rna" / f"{donor}_{serial}_{compartment}_{name}.npy"
            vp.parent.mkdir(parents=True, exist_ok=True)
            np.save(vp, vals)
            entry = score_against_genes(vals, dd, by_hint)
            entry["values_npy"] = str(vp)
            if name in channels:
                entry.update({"qptiff_channel": channels[name][1],
                              "assigned_name": roles.names[name],
                              "channel_source": roles.source,
                              "channel_page": int(channels[name][0]),
                              "recorded_role": name,
                              "median": float(np.nanmedian(vals)),
                              "p99": float(np.nanpercentile(vals, 99))})
            else:
                entry["definition"] = ("percentile rank of INS minus percentile rank of GCG, "
                                       "within one cell; removes what the two channels share "
                                       "(islet membership) and keeps what differs (subtype)")
            comp["markers"][name] = entry
            b, a = entry.get("mean_rho_beta"), entry.get("mean_rho_alpha")
            print(f"[protein-rna]   {compartment}/{name}: mean rho beta {b:+.3f} / "
                  f"alpha {a:+.3f}   "
                  + ", ".join(f"{g} {v['spearman_rho']:+.3f}"
                              for g, v in entry["genes"].items()
                              if v["subtype_hint"] in ("beta", "alpha")), flush=True)
        del lab
        rec["compartments"][compartment] = comp

    return rec


def registered_transform_records(roi: str = "panc", *, out_dir: Path = OUT) -> dict[str, dict]:
    """serial -> exact donor transform whose named all-window gate passed.

    This is what the module actually needs: a refined transform that cleared the release gate.
    It used to enumerate from the QuPath GeoJSON directory -- three files -- while stating in three
    places that it reads no annotation, so ten donors with a PASS gate were skipped in silence and
    the lane looked like a registration failure. A gate admits only the exact saved transform it
    names and hashes; a PASS for another candidate cannot certify the file this consumer measures.
    A passing ``refined`` transform is preferred, followed by a separately passing ``image_seed``.
    """
    import json as _json
    out: dict[str, dict] = {}
    root = Path(out_dir) / "valis"
    if not root.is_dir():
        return out
    for donor_dir in sorted(root.glob("donor_id=*")):
        donor = donor_dir.name.split("=", 1)[1]
        try:
            export = _json.loads((donor_dir / "export.json").read_text())
        except (OSError, ValueError):
            continue
        if export.get("roi") != roi:
            continue
        serial = str(export.get("serial") or "")
        if serial:
            # Fail closed when different transforms passed; the extractor's choice must be
            # deterministic, recorded, and bound to its own exact PASS gate.
            for stem in ("refined", "image_seed"):
                transform = donor_dir / f"{stem}.json"
                gate = donor_dir / f"alignment_{stem}.json"
                if not transform.is_file() or not gate.is_file():
                    continue
                try:
                    record = _json.loads(gate.read_text())
                    alignment = record.get("alignment", record)
                except (OSError, ValueError):
                    continue
                if str(alignment.get("verdict", "")).upper() != "PASS":
                    continue
                certified = alignment.get("transform_file")
                if not certified or Path(certified).resolve() != transform.resolve():
                    continue
                expected_sha = alignment.get("transform_sha256")
                if not expected_sha:
                    continue
                from phenocycler.artifacts import sha256_file
                if sha256_file(transform) != expected_sha:
                    continue
                out[serial] = {"donor_id": donor, "transform_stem": stem,
                               "gate": str(gate.resolve()),
                               "transform_sha256": expected_sha}
                break
    return out


def registered_serials(roi: str = "panc", *, out_dir: Path = OUT) -> dict[str, str]:
    """serial -> donor for cells admitted by their exact passing transform gate."""
    return {serial: record["donor_id"]
            for serial, record in registered_transform_records(roi, out_dir=out_dir).items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--um-per-px", type=float, default=PROTEIN_UM_PER_PX)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    (args.out / "protein_rna").mkdir(parents=True, exist_ok=True)

    registrations = registered_transform_records(args.roi, out_dir=args.out)
    serials = {serial: record["donor_id"] for serial, record in registrations.items()}
    resolved = resolve_all(serials, sections=sections)
    report = channel_assignment_report(resolved, roi=args.roi)
    report_path = write_channel_assignment(report, args.out)
    print(f"[protein-rna] channel assignment: {report['n_resolved']} of "
          f"{report['n_registered']} registered sections resolve; wrote {report_path}",
          flush=True)
    for reason, donors in report["donors_by_reason"].items():
        print(f"[protein-rna]   {reason}: {', '.join(donors)}", flush=True)

    recs = {}
    skipped: list[str] = []
    for serial, donor in sorted(serials.items()):
        if args.donors and donor not in args.donors:
            continue
        roles = resolved[serial]
        if not roles.ok:
            # SKIP BY NAME, out loud. This used to be an IndexError on `.iloc[0]` against a
            # five-row CSV, which stopped the whole run on the first donor the CSV had never
            # heard of -- eleven of the thirteen that pass their registration gate.
            print(f"[protein-rna] {donor}: SKIPPED -- {roles.unavailable}: {roles.detail}",
                  flush=True)
            skipped.append(donor)
            continue
        rec = measure_one(serial, sections, roi=args.roi, out_dir=args.out,
                          um_per_px=args.um_per_px, roles=roles,
                          transform_stem=registrations[serial]["transform_stem"])
        rec.setdefault("serial", serial)
        rec.setdefault("donor", donor)
        rec["registration_gate"] = registrations[serial]["gate"]
        rec["registration_transform_sha256"] = registrations[serial]["transform_sha256"]
        recs[serial] = rec
    if not recs:
        print(f"no donor measured; {len(skipped)} skipped for a named reason "
              f"({', '.join(skipped) or 'none'}) -- see {report_path}", file=sys.stderr)
        return 1
    p = args.out / f"protein_rna_{args.roi}.json"
    # MERGE, never replace. `--donors 6539` used to rewrite this file with that one record and
    # drop 6436 and 6442 from it -- the per-donor .npy arrays survived, so nothing was
    # measured twice, but the index silently lost every donor not named on the command line.
    # Keyed by serial so a re-run of a donor replaces its own record and no other.
    if p.exists():
        try:
            prior = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            prior = []
        # Records written before `serial` existed carry None, so they cannot be matched by
        # key. Identify those by what they measured: a prior record whose n_cells equals one
        # just written is the same section, and the fresh record wins.
        measured_now = {r.get("n_cells") for r in recs.values()}
        for rec in prior if isinstance(prior, list) else []:
            key = rec.get("serial")
            if key is None:
                if rec.get("n_cells") in measured_now:
                    continue
                key = f"n_cells={rec.get('n_cells')}"
            recs.setdefault(key, rec)
    merged = sorted(recs.values(), key=lambda r: str(r.get("serial") or ""))
    p.write_text(json.dumps(merged, indent=2, allow_nan=False))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
