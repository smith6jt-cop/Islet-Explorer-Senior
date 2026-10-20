#!/usr/bin/env python3
"""Recover what the ingest dropped: the object-classifier calls and the cell morphometry.

The pipeline's per-cell table keeps intensities, a parent string and two shape numbers
(``cell_area``, ``cell_solidity``).  The raw QuPath exports it was built from carry far more, and
two of those columns are load-bearing for work the plan asks for:

  ``Classification``   the QuPath OBJECT-classifier call per cell -- the galectin apparatus
                       (``M2Galectin3``, ``Iba1_CD68`` and the rest).  Correction J had to read
                       inflated per-annotation counts instead because this column never reached
                       the run.  The maintainer confirms these calls are unchanged by the
                       re-export, so they come from the FIRST export for every donor.
  13 morphometry cols  nucleus and cell area, length, circularity, solidity, max/min diameter,
                       and the nucleus/cell area ratio.  The islet chain decides endocrine-vs-not
                       on CD99 and, in the large stratum only, on **NUCLEUS** shape --
                       ``nucleus_elongation``, never a cell-boundary feature, because cell
                       circularity runs 0.514-0.523 across every class and is a tessellation
                       artefact.  With only area and solidity there is no shape basis at all.
                       (Corrected 2026-09-18: this said "cell shape".)

Everything is keyed by ``Object ID`` -- the QuPath UUID -- so it joins to the run by ``object_id``
with no centroid matching anywhere.  Verified on 6505: 2,863,277 of the revised export's
2,870,418 ids (99.8%) are present in the run.

A v3 image's source is its own cells Parquet (``extract_parquet``); the three CSV exports below
are the v2 cohort's, ~215 GB, streamed once each with nothing loaded whole:

    results/Cellmeasurements.csv                  147.9 GB   12 batch-1 donors
    results/CellmeasurementsBatch2.csv             64.2 GB    7 batch-2 donors
    results/Cellmeasurements_6539parentFixed.csv    3.4 GB    6539

    python core/multiplex5000/scripts/typing/extract_qupath_fundamentals.py [--donors 6505]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import subpanel_io as S

OUT = S.FUNDAMENTALS

#: Carried verbatim from the export header. The exports differ in width (992 vs 994 columns) and
#: in marker order, so every column is resolved BY NAME per file, never by position.
IDENTITY = ["Object ID", "Image", "Classification", "Parent"]
MORPHOMETRY = [
    "Nucleus: Area µm^2", "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
    "Nucleus: Max diameter µm", "Nucleus: Min diameter µm",
    "Cell: Area µm^2", "Cell: Length µm", "Cell: Circularity", "Cell: Solidity",
    "Cell: Max diameter µm", "Cell: Min diameter µm",
    "Nucleus/Cell area ratio",
]
#: Output names: the source headers carry units, colons and a micro sign.
RENAME = {
    "Object ID": "object_id", "Image": "image", "Classification": "qupath_classification",
    "Parent": "parent_raw",
    "Nucleus: Area µm^2": "nucleus_area_um2", "Nucleus: Length µm": "nucleus_length_um",
    "Nucleus: Circularity": "nucleus_circularity", "Nucleus: Solidity": "nucleus_solidity",
    "Nucleus: Max diameter µm": "nucleus_max_diam_um",
    "Nucleus: Min diameter µm": "nucleus_min_diam_um",
    "Cell: Area µm^2": "cell_area_um2", "Cell: Length µm": "cell_length_um",
    "Cell: Circularity": "cell_circularity", "Cell: Solidity": "cell_solidity",
    "Cell: Max diameter µm": "cell_max_diam_um", "Cell: Min diameter µm": "cell_min_diam_um",
    "Nucleus/Cell area ratio": "nucleus_cell_area_ratio",
}
FLOAT_COLUMNS = {RENAME[c] for c in MORPHOMETRY}
#: Rows per parquet row-group; keeps peak memory flat regardless of donor size.
CHUNK = 500_000


def source_files() -> dict[str, list[str]]:
    """raw export path -> the donors it contains, from the manifest rather than a guess."""
    manifest = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    mapping: dict[str, list[str]] = {}
    stems: dict[str, str] = {}
    for record in manifest["images"]:
        path = record["measurement_csv"]["path"]
        mapping.setdefault(path, []).append(record["donor_id"])
        stems[record["image_id"]] = record["donor_id"]
    return mapping, stems


class DonorWriter:
    """One parquet writer per donor, flushed in row groups."""

    def __init__(self, donor: str, schema: pa.Schema):
        target = OUT / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        self.path = target / "data.parquet"
        self.writer = pq.ParquetWriter(self.path, schema, compression="zstd")
        self.schema = schema
        self.rows: list[list] = []
        self.n = 0

    def add(self, row: list) -> None:
        self.rows.append(row)
        if len(self.rows) >= CHUNK:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return
        columns = list(zip(*self.rows))
        arrays = []
        for name, values in zip(self.schema.names, columns):
            if name in FLOAT_COLUMNS:
                arrays.append(pa.array([_f(v) for v in values], type=pa.float32()))
            else:
                arrays.append(pa.array([v if v != "" else None for v in values], type=pa.string()))
        self.writer.write_table(pa.Table.from_arrays(arrays, schema=self.schema))
        self.n += len(self.rows)
        self.rows = []

    def close(self) -> int:
        self.flush()
        self.writer.close()
        return self.n


def _f(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def extract(path: str, donors: list[str], stems: dict[str, str], wanted: set[str] | None) -> None:
    keep = [d for d in donors if wanted is None or d in wanted]
    if not keep:
        print(f"  skip {path} (no requested donor in it)")
        return
    print(f"\n{path}\n  donors: {' '.join(sorted(keep))}", flush=True)

    fields = IDENTITY + MORPHOMETRY
    schema = pa.schema([
        pa.field(RENAME[c], pa.float32() if RENAME[c] in FLOAT_COLUMNS else pa.string())
        for c in fields
    ])

    started = time.time()
    writers: dict[str, DonorWriter] = {}
    unmatched: dict[str, int] = {}
    seen = 0
    # csv.reader on a plain file object streams; the sniffed delimiter is taken from the suffix
    # because a 148 GB sniff is not free and both exports are unambiguous.
    delimiter = "\t" if path.endswith(".tsv") else ","
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        missing = [c for c in fields if c not in index]
        if missing:
            raise KeyError(f"{path}: missing {missing}")
        picks = [index[c] for c in fields]
        image_at = index["Image"]

        for row in reader:
            seen += 1
            donor = stems.get(row[image_at])
            if donor is None:
                unmatched[row[image_at]] = unmatched.get(row[image_at], 0) + 1
                continue
            if donor not in keep:
                continue
            writer = writers.get(donor)
            if writer is None:
                writer = writers[donor] = DonorWriter(donor, schema)
            writer.add([row[i] for i in picks])
            if seen % 5_000_000 == 0:
                print(f"    {seen:,} rows  {time.time() - started:.0f}s", flush=True)

    for donor, writer in sorted(writers.items()):
        n = writer.close()
        print(f"    {donor}: {n:>10,} cells -> {writer.path}", flush=True)
    if unmatched:
        print(f"    images with no manifest donor: {dict(list(unmatched.items())[:4])}")
    print(f"  {seen:,} rows in {time.time() - started:.0f}s", flush=True)


def extract_parquet(path: str, donors: list[str], stems: dict[str, str],
                    wanted: set[str] | None) -> None:
    """The same columns from a v3 cells Parquet -- typed, one image per file, read by name.

    A v3 image's measurement source is a Parquet written from the QuPath export with its header
    unchanged, so the seventeen columns are asked for under exactly the names the CSV carries and
    a missing one is the same refusal.  The output is the same per-donor table the CSV path writes.
    """
    keep = [d for d in donors if wanted is None or d in wanted]
    if not keep:
        print(f"  skip {path} (no requested donor in it)")
        return
    print(f"\n{path}\n  donors: {' '.join(sorted(keep))}", flush=True)
    fields = IDENTITY + MORPHOMETRY
    missing = [c for c in fields if c not in set(pq.read_schema(path).names)]
    if missing:
        raise KeyError(f"{path}: missing {missing}")
    schema = pa.schema([
        pa.field(RENAME[c], pa.float32() if RENAME[c] in FLOAT_COLUMNS else pa.string())
        for c in fields
    ])
    started = time.time()
    table = pq.read_table(path, columns=fields)
    image = table.column("Image").to_pandas()
    donor_of = image.map(stems)
    unmatched = image[donor_of.isna()].value_counts().to_dict()
    for donor in sorted(keep):
        rows = table.filter(pa.array((donor_of == donor).to_numpy()))
        if rows.num_rows == 0:
            continue
        arrays = [rows.column(c).cast(schema.field(RENAME[c]).type, safe=False) for c in fields]
        target = OUT / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_arrays(arrays, schema=schema), target / "data.parquet",
                       compression="zstd", row_group_size=CHUNK)
        print(f"    {donor}: {rows.num_rows:>10,} cells -> {target / 'data.parquet'}", flush=True)
    if unmatched:
        print(f"    images with no manifest donor: {dict(list(unmatched.items())[:4])}")
    print(f"  {table.num_rows:,} rows in {time.time() - started:.0f}s", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    args = parser.parse_args()
    wanted = set(args.donors) if args.donors else None

    mapping, stems = source_files()
    OUT.mkdir(parents=True, exist_ok=True)
    for path, donors in mapping.items():
        if str(path).endswith(".parquet"):
            extract_parquet(path, donors, stems, wanted)
        else:
            extract(path, donors, stems, wanted)
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
