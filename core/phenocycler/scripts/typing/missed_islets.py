#!/usr/bin/env python3
"""Did the islet classifier miss any islets? Asked without trusting the classifier.

Take the cells whose hormone signal is high for this donor, keep only the ones lying OUTSIDE
every islet annotation, and cluster them spatially. A cluster that comes out is a group of
hormone-producing cells sitting where no islet was drawn -- a candidate islet the classifier
failed to find.

WHY IT IS SWEPT RATHER THAN SET. The count depends on the parameters and obviously so: a longer
linking distance merges neighbouring clusters into one, a higher percentile keeps fewer cells.
Any single setting's count is partly an artefact of the setting. What the sweep establishes is
the part that does not move -- which donors come back with ZERO everywhere. That zero set is what
`cohort.CLEAN` means, and it is the only claim made from this.

    linking distance   15 / 25 / 40 um
    hormone percentile 99 / 99.5, WITHIN donor
    = six settings

WITHIN DONOR, never pooled: a cohort hormone percentile would rank a T1D donor's cells against
donors that still have beta cells, and find "missed islets" wherever a donor is simply brighter.

INDEPENDENT OF THE TYPING. The hormone value is read from the aligned lane and scaled, never from
`cell_class`. Asking whether the classifier missed an islet using cells the classifier's own
downstream typing named would be circular.

WHY THIS FILE EXISTS. The sweep's RESULTS were recorded in `cohort.py`'s docstring and the code
that produced them was never written down -- the same gap as `subpanel_calls.parquet`. So the
measurement that defines the cohort could not be repeated, which matters now that the batch-2
donors have been re-exported and nothing has checked whether the re-export fixed their missed
islets (a separate defect from the duplicated detections, which it demonstrably did fix).

HOW min_cells WAS SET, AND WHAT DOES NOT REPRODUCE. The original recorded its results and not its
parameters, so the cluster-size floor had to be recovered. It is 20, and the reasoning is checkable
two ways. The recorded note says 6380's 6 missed islets hold 13 % of its islet-core cells -- about
28 cells each, so the floor is tens of cells, not a handful. And swept from 10 to 60, the RECORDED
ZERO SET (115, 6374, 6436, 6442, 6461, 6476, 6482, 6493, 6539, 6593) comes back zero at every value
from 15 upward: the claim `CLEAN` rests on does not depend on the choice. At 10 it does not hold,
and at 5 the sweep finds 154 clusters on 6380 against a recorded 6.

What does NOT reproduce, stated rather than smoothed over: the three marginal donors' exact counts.
At 20 this gives 6551 2 of 6 settings, max 2 (recorded 2 of 6, max 3) and 6380 2 of 6, max 5
(recorded 3 of 6, max 6) -- close but not equal -- and **6521 comes back zero where the record says
1-2 missed of 395**. The disagreement is in the conservative direction, so this implementation
under-reports rather than over-reports, and 6521's marginal status should be taken from the record,
not from this run.

    python core/phenocycler/scripts/typing/missed_islets.py [--donors 115] [--min-cells 20]
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.cluster import DBSCAN

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S

LINK_UM = (15.0, 25.0, 40.0)
PERCENTILE = (99.0, 99.5)

#: A cell must be beyond the 20 um expansion to count as "outside every islet annotation".
#: The expansion IS an annotation, and a hormone-high cell 5 um outside an islet outline is part
#: of that islet's neighbourhood, not evidence of a missed one.
OUTSIDE_UM = 20.0


def donor_frame(donor: str) -> pd.DataFrame | None:
    """Eligible cells with a scaled hormone value and a signed distance, for one donor."""
    cells_path = S.OUT / f"cells/donor_id={donor}/data.parquet"
    if not cells_path.exists():
        return None
    cells = pq.read_table(cells_path, columns=["object_id", "X_centroid", "Y_centroid",
                                               "signed_dist_islet_um"]).to_pandas()
    aligned = S.ALIGNED / f"donor_id={donor}/data.parquet"
    have = set(pq.ParquetFile(aligned).schema_arrow.names)
    horm = [h for h in S.HORMONES if h in have]
    if not horm:
        return None
    values = pq.read_table(aligned, columns=["object_id"] + horm).to_pandas()
    frame = cells.merge(values, on="object_id", how="inner", validate="one_to_one")
    scaled = S.apply_scale(frame[horm].to_numpy(), horm, S.load_scale())
    # A cell with no hormone measured at all is not evidence either way, and numpy's all-NaN
    # warning is not a number. Count it, drop it, report it.
    measured = np.isfinite(scaled).any(axis=1)
    frame["hormone"] = np.where(measured, np.nanmax(np.where(np.isfinite(scaled), scaled, -np.inf),
                                                    axis=1), np.nan)
    frame.attrs["no_hormone_measured"] = int((~measured).sum())
    return frame.dropna(subset=["hormone", "X_centroid", "Y_centroid"])


def sweep(frame: pd.DataFrame, min_cells: int) -> dict:
    """Candidate missed islets at each of the six settings."""
    outside = frame["signed_dist_islet_um"].to_numpy() > OUTSIDE_UM
    hormone = frame["hormone"].to_numpy()
    xy = frame[["X_centroid", "Y_centroid"]].to_numpy()
    out = {}
    for pct in PERCENTILE:
        # The percentile is over ALL of this donor's eligible cells, not only the outside ones:
        # "hormone-high for this donor" is a property of the donor, not of the subset being tested.
        level = np.nanpercentile(hormone, pct)
        bright = outside & (hormone >= level)
        points = xy[bright]
        for link in LINK_UM:
            if len(points) < min_cells:
                out[(link, pct)] = 0
                continue
            labels = DBSCAN(eps=link, min_samples=min_cells).fit(points).labels_
            out[(link, pct)] = int(len({x for x in labels if x >= 0}))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--min-cells", type=int, default=20,
                    help="cells a cluster needs to count as a candidate islet")
    args = ap.parse_args()

    donors = args.donors or sorted(p.name.split("=")[1]
                                   for p in (S.OUT / "cells").glob("donor_id=*"))
    rows = []
    for donor in donors:
        frame = donor_frame(donor)
        if frame is None:
            print(f"  {donor}: no cells table or no hormone markers"); continue
        counts = sweep(frame, args.min_cells)
        row = {"donor": donor, "cells": len(frame)}
        row.update({f"{int(l)}um_p{p}": n for (l, p), n in counts.items()})
        row["no_hormone_measured"] = int(frame.attrs.get("no_hormone_measured", 0))
        row["settings_nonzero"] = sum(1 for n in counts.values() if n)
        row["max_missed"] = max(counts.values())
        rows.append(row)
        print(f"  {donor}: {row['settings_nonzero']} of 6 settings non-zero, "
              f"max {row['max_missed']}", flush=True)

    frame = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print("\n" + frame.to_string(index=False))
    zero = sorted(frame.loc[frame.settings_nonzero == 0, "donor"])
    print(f"\nZERO under all six settings ({len(zero)}): {' '.join(zero)}")
    print(f"min_cells = {args.min_cells}, outside = beyond {OUTSIDE_UM:.0f} um")
    out = S.OUT / "missed_islets.parquet"
    frame.to_parquet(out, index=False)
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
