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
* only the IO60 donors are in scope; a Tissue_id whose donor the IO60 sources do not know (6426,
  the HPAP sections) is counted as out of cohort and not ingested.

A ``Tissue_id``'s donor is the text before its first underscore (the whole id where it has none),
through the cohort's declared donor aliases (HDL115 is donor 115); the rest of the id is kept as
text (``id_block_token``) and never interpreted.

Each MALDI raster is drawn to its tissue's outline, so the tissue mask is the raster footprint:
a spot is on the tissue unless no ion was detected in it at all (``ON_TISSUE_METHOD``). That is a
property of the spot, never a cutoff on a feature's value.
"""

from __future__ import annotations

import re
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


def tissue_donor(tissue_id: str, aliases: dict) -> tuple[str, str, str]:
    """(the donor label the id carries, the donor id through the cohort aliases, the rest of the
    id as text). Nothing in the rest is interpreted."""
    label, _, rest = str(tissue_id).partition("_")
    return label, str(aliases.get(label, label)), rest


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


def safe_key(tissue_id: str) -> str:
    """A Tissue_id as a partition or file-name key: unchanged where it is already safe."""
    if not re.fullmatch(r"[A-Za-z0-9._-]+", tissue_id):
        raise MaldiSourceError(f"Tissue_id {tissue_id!r} is not a safe partition key")
    return tissue_id
