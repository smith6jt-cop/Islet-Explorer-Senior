#!/usr/bin/env python3
"""Recover the real batch structure from the qptiff headers.

`data/donor_acquisition_batch.csv` carries the panel variant and the Fusion version, and the
Fusion version is confounded with disease status -- no ND donor was acquired on 2.4.1 -- so on
its own it cannot serve as a batch variable: correcting for it would remove the effect being
measured.

Two finer factors are recoverable, and they are not confounded the same way:

  * slides **run together** share a scan day;
  * slides **stained together** fall on consecutive days with no gap.

Both come from the TIFF ``DateTime`` tag, which is per slide.  The qptiff XML ``<Date>`` is the
export date and is useless here -- it collapses to exactly the three Fusion versions.

Four of the nine staining groups this yields contain more than one disease status, which is what
makes correction possible at all:

    S0  6521 Aab+ · 6476 ND          S4  6493 ND · 6539 ND · 6593 T1D · 6442 ND · 6461 ND
    S6  6505 Aab+ · 6623 T1D         S7  6450 Aab+ · 6534 T1D · 6566 T1D · 6591 T1D

Reads only each qptiff's first page (no pixel data), so it is safe to run against /data.

THE SHIPPED CSV IS STALE ON ONE DONOR, measured 2026-09-18.  ``excluded`` is derived here from
``phenocycler.cohort.DONOR_EXCLUSIONS`` and is never a literal, but ``data/donor_acquisition_runs.csv``
was last written 2026-09-12 and records ``6579,...,True``.  6579's exclusion was LIFTED
2026-09-17 -- ``cohort.py`` now carries it in ``REINSTATED`` -- so the file disagrees with the
authority until this generator reruns against the images.  The pancreatitis it was excluded for
is a disease context, not a reason to drop it.  Only 6457 is centrally excluded today.  A reader
that needs the exclusion state should ask ``phenocycler.cohort``, not this column.

    python scripts/build_acquisition_runs.py --images ~/IO60panc2nd/Images
    -> data/donor_acquisition_runs.csv
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import pathlib
import re
import sys

import pandas as pd
import tifffile

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "core" / "multiplex5000"))
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

#: `<donor>[panc|pacn]_Scan<n>.er.qptiff`, plus the HDL-prefixed lymph-node scan.  Used ONLY for a
#: slide the cohort specification does not declare (6457, which is listed nowhere because it is
#: excluded); a declared slide takes its donor from the specification -- see `declared_donors`.
NAME = re.compile(r"(?:HDL)?(\d+)(?:panc|pacn)?(?:LN)?_Scan\d+\.er\.qptiff$")

#: The tracked cohort specification: the one place an image is assigned to a donor.
COHORT_SPEC = REPO / "core" / "multiplex5000" / "cohort-spec.json"


def declared_donors(spec: pathlib.Path = COHORT_SPEC) -> dict[str, str]:
    """qptiff file name -> donor id, for every image the cohort specification declares.

    The same-slide images carry a slide number before the donor (`0080870_6533_pLNpanc_Scan3`),
    which the name pattern above does not read -- so they had no row here, and a lane that needs a
    donor's acquisition run (`build_cells.py`) refused them.  Their donor is declared, not parsed.
    """
    import json
    if not spec.is_file():
        return {}
    return {os.path.basename(row["qptiff"]): str(row["donor_id"])
            for row in json.loads(spec.read_text())["images"]}


def _read_header(path: str) -> dict:
    with tifffile.TiffFile(path) as handle:
        page = handle.pages[0]
        tags = {tag.name: tag.value for tag in page.tags}
        description = page.description
    def field(name: str) -> str | None:
        match = re.search(rf"<{name}[^>]*>([^<]*)</{name}>", description)
        return match.group(1) if match else None
    return {
        "scan_time": tags.get("DateTime"),
        "software": field("AcquisitionSoftware"),
        "operator": field("OperatorName"),
        "computer": field("ComputerName"),
        "validation_code": field("ValidationCode"),
    }


def build(images: pathlib.Path, declared: dict[str, str] | None = None) -> pd.DataFrame:
    declared = declared_donors() if declared is None else declared
    rows = []
    for path in sorted(glob.glob(str(images / "*.qptiff"))):
        name = os.path.basename(path)
        donor = declared.get(name)
        if donor is None:
            match = NAME.match(name)
            if not match:
                continue
            donor = match.group(1)
        header = _read_header(path)
        if not header["scan_time"]:
            raise ValueError(f"{path}: no TIFF DateTime tag -- cannot place it in a run")
        rows.append({"donor_id": donor, "image_file": name, **header})
    frame = pd.DataFrame(rows)
    frame["scan_time"] = pd.to_datetime(frame["scan_time"], format="%Y:%m:%d %H:%M:%S")
    frame = frame.sort_values("scan_time").reset_index(drop=True)

    # A run is a scan day.  A staining group is a maximal span of consecutive days: the gap to
    # the previous scan day is what separates one staining session from the next.
    frame["acq_run"] = frame["scan_time"].dt.date.astype(str)
    days = sorted(frame["scan_time"].dt.date.unique())
    group_of, group = {}, 0
    for index, day in enumerate(days):
        if index and (day - days[index - 1]).days > 1:
            group += 1
        group_of[day] = group
    frame["stain_group"] = frame["scan_time"].dt.date.map(group_of).map("S{}".format)
    frame["excluded"] = frame["donor_id"].isin(DONOR_EXCLUSIONS)
    return frame[["donor_id", "image_file", "scan_time", "acq_run", "stain_group",
                  "software", "operator", "computer", "validation_code", "excluded"]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--images", type=pathlib.Path,
                        default=pathlib.Path.home() / "IO60panc2nd/Images")
    parser.add_argument("--cohort-spec", type=pathlib.Path, default=COHORT_SPEC,
                        help="the cohort specification whose images take their donor from it")
    parser.add_argument("--out", type=pathlib.Path,
                        default=REPO / "data/donor_acquisition_runs.csv")
    args = parser.parse_args()

    frame = build(args.images, declared_donors(args.cohort_spec))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)

    runs = frame["acq_run"].nunique()
    groups = frame["stain_group"].nunique()
    print(f"{len(frame)} slides -> {runs} acquisition runs, {groups} staining groups")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
