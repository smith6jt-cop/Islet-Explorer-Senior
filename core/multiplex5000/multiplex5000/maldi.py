"""MALDI imaging spots: the recorded export, read as it is (2026-10-31).

The MALDI acquisition is exported as two CSV files over the same spots -- one of N-glycan
compositions, one of glycogen m/z features -- each row a spot on a 50 um grid in the instrument's
stage frame, keyed by ``spotId`` and carrying the ``Tissue_id`` the MALDI lab gave the section.
Both files must hold exactly the same spots, ids, positions, rasters and tissues; a difference is
refused, never reconciled.

What the maintainer decided (2026-10-31), and this module keeps:

* a value of 0 is **not detected** and is kept as a measured value; a value missing from the file
  is unavailable and counted, never filled;
* the MALDI sections are serial sections from the same block series as the PhenoCycler and
  Xenium slides -- different cells, never a same-cell claim. Which tissue a Tissue_id is, and which
  slide it neighbours, is recorded from what is seen (``maldi_tissue_inspection.csv``), never
  parsed from the id;
* only donors with PhenoCycler or Xenium data are in scope (maintainer, 2026-11-03); any other
  Tissue_id (the HPAP sections) is counted as out of cohort and not ingested.

A ``Tissue_id``'s donor is the text before its first underscore (the whole id where it has none),
through the cohort's declared donor aliases (HDL115 is donor 115), unless the maintainer has
recorded that the id names the wrong donor (``TISSUE_DONOR_CORRECTIONS``); the rest of the id is
kept as text (``id_block_token``) and never interpreted. The Tissue_id itself is never rewritten:
it stays the section's key, as the export writes it.

Each MALDI raster is drawn to its tissue's outline, so the tissue mask is the raster footprint:
a spot is on the tissue unless no ion was detected in it at all (``ON_TISSUE_METHOD``). That is a
property of the spot, never a cutoff on a feature's value.

Each N-glycan feature is described by its composition, read from its name
(``glycan_annotation``): the monosaccharide counts, the ion written after them, and the
composition class and modification flags of the CRO N-glycan report's rules (``MALDI/
T1D_Nglycan_CRO_Report.pdf``, 2026-07-09, section 9.1). These are descriptions of a feature,
never of a spot's or a cell's identity, and composition only (Schymanski level 3-4): linkage,
core versus antennary fucose and branching are not resolved by the mass. The report covers the
same acquisition through a different glycan table (49 columns against the 38 here, the shared
ones peak-picked slightly differently), so a feature list is always read from the file, never
assumed.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .integration.rasterize import GridSpec

#: The columns of both export files that are not features, as the files write them. The first,
#: unnamed column is R's row number and is not read.
META = ("spotId", "raster", "x", "y", "z", "Tissue_id", "Date")
#: The feature sets, by the file they are exported in: (role, file-name suffix).
ROLES = (("glycan", "_glycan_pbp.csv"), ("glycogen", "_glycogen_pbp.csv"))
#: The grid the spots sit on, in um (the platform profile's pixel_size_um).
SPOT_PITCH_UM = 50.0
#: How a spot is placed on or off the tissue: inside the raster drawn to the tissue's outline, and
#: at least one ion detected in it.
ON_TISSUE_METHOD = "raster_footprint_any_ion"
#: The checksum file that names, and binds, the export files in the source folder.
CHECKSUMS = "SHA256SUMS"


class MaldiSourceError(ValueError):
    """A MALDI export that cannot be read as recorded."""


@dataclass(frozen=True)
class SourceFile:
    role: str
    path: Path
    sha256: str


def source_files(folder: Path) -> list[SourceFile]:
    """The export files of a MALDI source folder, each with the digest its ``SHA256SUMS`` records:
    exactly one file per feature set."""
    folder = Path(folder)
    sums = folder / CHECKSUMS
    if not sums.is_file():
        raise MaldiSourceError(f"{folder}: no {CHECKSUMS}; the export files must be bound by digest")
    recorded = {}
    for line in sums.read_text().splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            recorded[name.lstrip("*").strip()] = digest
    out = []
    for role, suffix in ROLES:
        names = [n for n in recorded if n.endswith(suffix)]
        if len(names) != 1:
            raise MaldiSourceError(f"{sums}: expected one *{suffix}, found {names}")
        out.append(SourceFile(role, folder / names[0], recorded[names[0]]))
    return out


#: The monosaccharide residues an N-glycan composition is written in, as the export names them.
RESIDUES = ("Hex", "HexNAc", "dHex", "NeuAc")
_RESIDUE = re.compile(r"(dHex|HexNAc|NeuAc|Hex)(\d+)")
_COMPOSITION = re.compile(r"(?:(?:dHex|HexNAc|NeuAc|Hex)\d+)+")
#: Where the composition classes come from, recorded with them.
COMPOSITION_RULE = ("CRO N-glycan report (2026-07-09) section 9.1, by composition only: "
                    "high_mannose = HexNAc 2, Hex >= 5, no dHex or NeuAc; paucimannose = "
                    "HexNAc 2, Hex <= 4; hybrid = HexNAc 3, Hex >= 5; complex = any other "
                    "HexNAc >= 3; otherwise unclassified")


def composition_class(hex_: int, hexnac: int, dhex: int, neuac: int) -> str:
    """An N-glycan composition's class by ``COMPOSITION_RULE``.

    The report's text bounds high-mannose at Hex 5-9, but its six high-mannose species include
    Hex10HexNAc2 (figures A15, A28); the rule here follows the species it classified (Hex >= 5).
    A composition no rule reaches (HexNAc 2 with Hex >= 5 and a dHex or NeuAc; HexNAc under 2) is
    ``unclassified``, never assigned the nearest class."""
    if hexnac == 2:
        if hex_ <= 4:
            return "paucimannose"
        return "high_mannose" if not dhex and not neuac else "unclassified"
    if hexnac == 3 and hex_ >= 5:
        return "hybrid"
    if hexnac >= 3:
        return "complex"
    return "unclassified"


def glycan_annotation(feature: str) -> dict:
    """What an N-glycan feature's name says: its composition (the residue counts written before
    the first ``.``), the ion as the export writes it after them (``..Na``, ``...2Na``; R's
    rewriting of the ion label, kept as text), the composition class and the fucosylated
    (dHex >= 1) and sialylated (NeuAc >= 1) flags. A name that is not a composition is
    ``unparsed`` with every count missing, never guessed."""
    head, dot, tail = str(feature).partition(".")
    counts = {r: 0 for r in RESIDUES}
    if not _COMPOSITION.fullmatch(head):
        return dict(feature=feature, composition="", ion_label=dot + tail, parsed=False,
                    **{r: None for r in RESIDUES}, composition_class="unparsed",
                    fucosylated=None, sialylated=None)
    seen = []
    for residue, n in _RESIDUE.findall(head):
        seen.append(residue)
        counts[residue] = int(n)
    if len(seen) != len(set(seen)):
        return dict(feature=feature, composition=head, ion_label=dot + tail, parsed=False,
                    **{r: None for r in RESIDUES}, composition_class="unparsed",
                    fucosylated=None, sialylated=None)
    return dict(feature=feature, composition=head, ion_label=dot + tail, parsed=True, **counts,
                composition_class=composition_class(counts["Hex"], counts["HexNAc"],
                                                    counts["dHex"], counts["NeuAc"]),
                fucosylated=counts["dHex"] >= 1, sialylated=counts["NeuAc"] >= 1)


def feature_table(feats: list[tuple[str, str]]) -> list[dict]:
    """One row per feature, in the export's order: its feature set, its name and, for an N-glycan,
    what its composition says (``glycan_annotation``), with ``ions_of_composition`` -- how many
    features measure the same composition as different ions, so a composition is never counted
    twice. A glycogen feature is named by its m/z (``X689``) and carries no composition."""
    rows = []
    for role, name in feats:
        row = dict(role=role, feature=name)
        if role == "glycan":
            row.update({k: v for k, v in glycan_annotation(name).items() if k != "feature"})
        rows.append(row)
    compositions = Counter(r["composition"] for r in rows if r.get("parsed"))
    for row in rows:
        row["ions_of_composition"] = compositions[row["composition"]] if row.get("parsed") else None
        if row["role"] == "glycan":
            row["composition_rule"] = COMPOSITION_RULE
    return rows


#: Tissue_ids whose donor label names the wrong donor: the donor the section belongs to, and the
#: maintainer's reason. Only the donor changes; the Tissue_id stays the section's key.
TISSUE_DONOR_CORRECTIONS: dict[str, tuple[str, str]] = {
    "6426_02_pln": ("6436", "there is no donor 6426 in the project, so the section is 6436's "
                            "(maintainer, 2026-11-03); the CRO N-glycan report (2026-07-09, figure "
                            "A8) labels the same section 6436_02_pln"),
}


def tissue_donor(tissue_id: str, aliases: dict) -> tuple[str, str, str]:
    """(the donor label the id carries, the donor id -- through the cohort aliases, or the
    recorded correction for this exact id -- and the rest of the id as text). Nothing in the rest
    is interpreted."""
    label, _, rest = str(tissue_id).partition("_")
    corrected = TISSUE_DONOR_CORRECTIONS.get(str(tissue_id))
    return label, corrected[0] if corrected else str(aliases.get(label, label)), rest


def _quoted(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _relation(path: Path) -> str:
    # The format is stated, never sniffed: comma-separated, double-quoted names, R's NA.
    return (f"read_csv('{Path(path)}', delim=',', quote='\"', escape='\"', header=true, "
            f"nullstr='NA', sample_size=-1)")


def feature_columns(con, path: Path) -> list[str]:
    """A file's feature columns, in the order it writes them."""
    names = [row[0] for row in con.execute(f"DESCRIBE SELECT * FROM {_relation(path)}").fetchall()]
    missing = [m for m in META if m not in names]
    if missing:
        raise MaldiSourceError(f"{path}: missing column(s) {missing}")
    first = names.index(META[-1]) + 1
    return names[first:]


def connect():
    import duckdb

    return duckdb.connect()


def features(con, files: list[SourceFile]) -> list[tuple[str, str]]:
    """(role, feature) for every feature, glycans first then glycogen; a name in two sets is
    refused, because the spot table holds each feature once."""
    out = [(f.role, name) for f in files for name in feature_columns(con, f.path)]
    names = [n for _, n in out]
    twice = sorted({n for n in names if names.count(n) > 1})
    if twice:
        raise MaldiSourceError(f"feature(s) named in more than one set: {twice}")
    return out


def spots_sql(files: list[SourceFile], feats: list[tuple[str, str]],
              tissue_ids: list[str] | None = None) -> str:
    """One row per spot, both feature sets joined on ``spotId``; ``mismatch`` is true where the two
    files disagree on the spot's position, raster or tissue."""
    by_role = {f.role: f for f in files}
    glycan, glycogen = by_role["glycan"], by_role["glycogen"]
    cols = ", ".join(f"{'g' if role == 'glycan' else 'k'}.{_quoted(n)}" for role, n in feats)
    where = ""
    if tissue_ids is not None:
        listed = ", ".join("'" + t.replace("'", "''") + "'" for t in tissue_ids)
        where = f"WHERE g.Tissue_id IN ({listed})"
    return (f"SELECT g.Tissue_id AS tissue_id, g.spotId AS spot_id, g.raster AS raster, "
            f"g.x AS x_um, g.y AS y_um, "
            f"(k.spotId IS NULL OR g.x IS DISTINCT FROM k.x OR g.y IS DISTINCT FROM k.y "
            f"OR g.raster IS DISTINCT FROM k.raster OR g.Tissue_id IS DISTINCT FROM k.Tissue_id) "
            f"AS mismatch, {cols} "
            f"FROM {_relation(glycan.path)} g LEFT JOIN {_relation(glycogen.path)} k "
            f"ON g.spotId = k.spotId {where} ORDER BY g.spotId")


def pairing_problems(con, files: list[SourceFile]) -> dict:
    """How far the two files are from holding the same spots: rows each has, spot ids only one
    has, ids repeated, and spots whose position, raster or tissue differ."""
    g, k = (_relation(f.path) for f in files)
    row = con.execute(
        f"WITH g AS (SELECT spotId, x, y, raster, Tissue_id FROM {g}), "
        f"k AS (SELECT spotId, x, y, raster, Tissue_id FROM {k}) "
        "SELECT (SELECT count(*) FROM g), (SELECT count(*) FROM k), "
        "(SELECT count(*) - count(DISTINCT spotId) FROM g), "
        "(SELECT count(*) - count(DISTINCT spotId) FROM k), "
        "(SELECT count(*) FROM g ANTI JOIN k USING (spotId)), "
        "(SELECT count(*) FROM k ANTI JOIN g USING (spotId)), "
        "(SELECT count(*) FROM g JOIN k USING (spotId) WHERE g.x IS DISTINCT FROM k.x "
        " OR g.y IS DISTINCT FROM k.y OR g.raster IS DISTINCT FROM k.raster "
        " OR g.Tissue_id IS DISTINCT FROM k.Tissue_id)").fetchall()[0]
    keys = ("glycan_rows", "glycogen_rows", "glycan_repeated_ids", "glycogen_repeated_ids",
            "only_in_glycan", "only_in_glycogen", "spots_that_differ")
    return dict(zip(keys, map(int, row), strict=True))


def tissue_summary(con, files: list[SourceFile]) -> list[dict]:
    """Per Tissue_id: spots, rasters and the stage-frame bounding box (um)."""
    g = _relation(next(f.path for f in files if f.role == "glycan"))
    rows = con.execute(
        f"SELECT Tissue_id, count(*), count(DISTINCT raster), min(x), max(x), min(y), max(y) "
        f"FROM {g} GROUP BY Tissue_id ORDER BY Tissue_id").fetchall()
    return [dict(tissue_id=r[0], n_spots=int(r[1]), rasters=int(r[2]), x_min_um=float(r[3]),
                 x_max_um=float(r[4]), y_min_um=float(r[5]), y_max_um=float(r[6])) for r in rows]


def summed_ions(con, files: list[SourceFile], feats: list[tuple[str, str]],
                tissue_ids: list[str]) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """{tissue: (x_um, y_um, summed ions)} for the named tissues, a missing value counted as
    nothing detected only in this sum (it is a picture, never a measurement)."""
    if not tissue_ids:
        return {}
    total = " + ".join(f"coalesce({'g' if r == 'glycan' else 'k'}.{_quoted(n)}, 0)"
                       for r, n in feats)
    sql = spots_sql(files, feats, tissue_ids)
    sql = sql.replace(" AS mismatch, ", f" AS mismatch, {total} AS summed, ", 1)
    data = con.execute(f"SELECT tissue_id, x_um, y_um, summed FROM ({sql})").fetchnumpy()
    out = {}
    for tissue in tissue_ids:
        mine = data["tissue_id"] == tissue
        out[tissue] = (np.asarray(data["x_um"][mine], float), np.asarray(data["y_um"][mine], float),
                       np.asarray(data["summed"][mine], float))
    return out


def on_tissue(summed: np.ndarray) -> np.ndarray:
    """``ON_TISSUE_METHOD``: a spot of the raster with at least one ion detected."""
    return np.asarray(summed, float) > 0


def spot_grid(x_um: np.ndarray, y_um: np.ndarray, *, margin_spots: int = 2) -> GridSpec:
    """The 50 um grid the spots sit on, one pixel per spot centred on it, with a margin.

    A pixel's centre is ``GridSpec.to_microns`` of its index (pixel i is centred on
    ``origin + i * pitch``), the convention the registration's rasters and transforms share: a
    half-pixel disagreement turns into a 70-110 um error once a rotation and mirror act on it."""
    margin = margin_spots * SPOT_PITCH_UM
    x0 = float(np.min(x_um)) - margin
    y0 = float(np.min(y_um)) - margin
    width = int(np.rint((float(np.max(x_um)) + margin - x0) / SPOT_PITCH_UM)) + 1
    height = int(np.rint((float(np.max(y_um)) + margin - y0) / SPOT_PITCH_UM)) + 1
    return GridSpec(x0, y0, SPOT_PITCH_UM, height, width)


def spot_image(x_um: np.ndarray, y_um: np.ndarray, values: np.ndarray, grid: GridSpec,
               *, fill: float = np.nan) -> np.ndarray:
    """The spots' values as an image on ``grid``, each in the pixel whose centre is nearest it; a
    pixel no spot falls in holds ``fill``."""
    image = np.full(grid.shape, fill, dtype=np.float32)
    px = grid.to_pixels(np.column_stack([x_um, y_um]))
    col = np.rint(px[:, 0]).astype(int)
    row = np.rint(px[:, 1]).astype(int)
    keep = (row >= 0) & (row < grid.height) & (col >= 0) & (col < grid.width)
    image[row[keep], col[keep]] = np.asarray(values, np.float32)[keep]
    return image


def signal_range(image: np.ndarray, tissue: np.ndarray | None = None) -> tuple[float, float]:
    """A section's actual biological low and high, in measured units: the minimum and maximum of
    its signal after a 3 x 3 median filter, over the on-tissue spots (``tissue``: at least one ion
    detected, ``ON_TISSUE_METHOD``) whose whole 3 x 3 neighbourhood is on the tissue. The filter
    takes out single-spot spikes; the neighbourhood rule leaves out the rim's partly ablated spots
    and the edges of tears and lumens, none of which is tissue signal. Nothing is a percentile or
    relative to another section (maintainer, 2026-11-01: "each actual biological max and min is
    the high and low of the colorscale")."""
    from scipy import ndimage

    image = np.asarray(image, float)
    spots = np.isfinite(image) & (np.asarray(tissue, bool) if tissue is not None else True)
    if not spots.any():
        return 0.0, 0.0
    filled = np.where(spots, np.maximum(image, 0), 0.0)
    smooth = ndimage.median_filter(filled, size=3, mode="constant", cval=0.0)
    interior = ndimage.binary_erosion(spots, np.ones((3, 3), bool))
    values = smooth[interior] if interior.any() else smooth[spots]
    return float(values.min()), float(values.max())


def display(image: np.ndarray, low: float, high: float) -> np.ndarray:
    """A MALDI image for display only: the measured value placed linearly between the section's
    own biological ``low`` and ``high`` (``signal_range``), as 0 to 1; beyond them, the end
    colours. A pixel with no spot stays NaN (drawn blank, never as 0). A section whose range is
    empty (high equal to low, e.g. a feature set it does not detect) is drawn at the low colour."""
    image = np.asarray(image, float)
    out = np.full(image.shape, np.nan)
    spots = np.isfinite(image)
    span = high - low
    out[spots] = (np.clip((np.maximum(image[spots], 0) - low) / span, 0, 1) if span > 0
                  else 0.0)
    return out


def range_text(low: float, high: float) -> str:
    return f"{low:,.0f} to {high:,.0f}" if high > low else f"no signal range ({high:,.0f})"


def safe_key(tissue_id: str) -> str:
    """A Tissue_id as a partition or file-name key: unchanged where it is already safe."""
    if not re.fullmatch(r"[A-Za-z0-9._-]+", tissue_id):
        raise MaldiSourceError(f"Tissue_id {tissue_id!r} is not a safe partition key")
    return tissue_id
