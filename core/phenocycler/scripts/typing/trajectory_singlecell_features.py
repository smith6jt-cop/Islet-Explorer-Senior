#!/usr/bin/env python3
"""Feature matrix for the single-cell trajectory over islet and peri-islet cells.

**The algorithm, pinned before implementing, from the paper's own Methods.**
Steenbuck, Damond, Engler, Kusmartseva, Posgai, Drotar, Williams, de Souza, Brusko, Brusko,
Wasserfall, Atkinson, Bodenmiller, *"Imaging mass cytometry reveals early beta-cell dysfunction and
changes in immune signatures during type 1 diabetes progression in human pancreata"*, bioRxiv
2025.03.05.641526 (v1, 11 March 2025), section **Pseudotime analysis**:

    "Pseudotime analysis of beta-cell expression profiles was performed using slingshot (v. 2.10.0).
     We used default parameters to identify the global lineage structure, and fit smooth branching
     curves to infer the pseudotime variables for each beta-cell.  Notably, we did not provide start
     and leaf nodes or cluster information and thus, slingshot inferred the global trajectory in an
     unbiased manner.  As input, we used diffusion map embedded single-cell beta-cell expression
     profiles, which were computed using destiny (v. 3.16.0).  We considered the expression of 28
     beta-cell lineage and functional markers for computation of the diffusion map and pseudotime."

So: **destiny diffusion map -> slingshot principal curve, default parameters, no start/leaf nodes
and no cluster labels.**  Their values are arcsinh counts, cofactor 1.  Both packages are installed
here, so the ordering is run as published rather than approximated by a different algorithm.

**Two deliberate departures, and why each is not a silent change:**

1. **The cells are not beta-cells.**  The paper orders beta-cells to describe beta-cell dysfunction.
   The question here is the immune environment of the islet as disease progresses, and the
   population named for it is *islet and peri-islet cells* -- all of them, endocrine and not.  That
   is a different question answered with the same machinery, not a failed reproduction of theirs.

2. **The markers are ours, and there are not 28 beta-cell functional markers in this panel.**  The
   feature set is exactly the markers the subpanel classification itself uses to decide what a cell
   is, read from ``subpanel_log.json`` and ``umap/refinement_log.json`` rather than chosen here:
   the islet subpanel, each structure's own markers, the exocrine trio, the immune definitions and
   the supporting tier.  Ordering cells on the evidence that named them is the defensible choice
   when the published feature list does not exist in this panel.

Values are the UniFORM-aligned ones on the classification's own floor and span **with no lower
clip**.  The clip pins 55.1 % of values to exactly zero, which is right for *is this marker present*
and wrong for *how similar are two cells* -- the same distinction the coarse embedding was built on.

**COVERAGE, 2026-09-17.**  The first build of this file ran on ``COHORT.resolve(None)``, which on
2026-09-14 was the ten batch-1 donors -- the seven re-exported batch-2 donors had not been rebuilt
yet and so were not in ``CLEAN``, and the three ``MARGINAL`` donors need ``--include-marginal``.
That is the whole reason the deployed trajectory covers ten donors, and it is why it carries no
Aab+ case at all: all three Aab+ donors (6450, 6505, 6521) were outside that ten.  The default is
now **every eligible donor** -- ``CLEAN`` plus ``MARGINAL``, with the centrally excluded donors
removed by ``phenocycler.cohort`` rather than by any list restated here.  ``--clean-only`` asks for
the narrower set.

    python core/phenocycler/scripts/typing/trajectory_singlecell_features.py
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

sys.path.insert(0, str(S.REPO / "core/phenocycler"))
from phenocycler.cohort import filter_eligible_donors  # noqa: E402  -- the one exclusion authority

OUT = S.OUT / "trajectory"

#: Markers that encode which panel a donor was on rather than which cell it is. Held out of any
#: similarity representation, exactly as the coarse embedding holds them out.
PANEL_IDENTITY = ("CD38", "b_Catenin1", "bCatenin1", "CD57")


def feature_markers() -> list[str]:
    """The markers the classification uses, read from its own logs."""
    log = json.loads((S.OUT / "subpanel/subpanel_log.json").read_text())
    markers: list[str] = list(log["islet_subpanel"])
    for panel in log["subpanels"].values():
        markers += list(panel)
    markers += [log["hormones"][k] for k in log["hormones"]]
    markers += list(log["bottom_tier_immune"])
    markers += list(S.SUBPANEL["exocrine"])
    # The immune definitions' markers, from the refinement log's own-marker table.
    refinement = json.loads((S.OUT / "umap/refinement_log.json").read_text())
    for step in refinement:
        for panel in step.get("own_markers", {}).values():
            markers += list(panel)
    markers += ["CD11b", "IDO1", "HLA_DR", "Granzyme_B"]
    ordered = sorted(set(markers) - set(PANEL_IDENTITY))
    return ordered


def donor_frame(donor: str, markers: list[str], scale: pd.DataFrame) -> pd.DataFrame:
    """Islet-core and peri-islet cells of one donor, scaled with no lower clip."""
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["object_id", "cell_class", "cell_region_20um", "islet_key", "phenotype",
                 "X_centroid", "Y_centroid", "signed_dist_islet_um"],
    ).to_pandas()
    cells = cells[(cells["cell_class"] != "Excluded")
                  & cells["cell_region_20um"].isin(["core", "expansion"])]
    if cells.empty:
        return cells

    available = pq.ParquetFile(S.ALIGNED / f"donor_id={donor}/data.parquet").schema_arrow.names
    have = [m for m in markers if m in available]
    values = pq.read_table(S.ALIGNED / f"donor_id={donor}/data.parquet",
                           columns=["object_id"] + have).to_pandas()
    cells = cells.merge(values, on="object_id", how="left", validate="one_to_one")

    # apply_scale's floor and span, WITHOUT the lower clip: the clip answers "is this marker
    # present", and this representation answers "how similar are two cells".
    floor = scale.loc[have, "floor"].to_numpy(float)
    span = scale.loc[have, "span"].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    scaled = (cells[have].to_numpy(dtype=np.float32) - floor) / span
    cells[have] = np.minimum(scaled, 1.0).astype(np.float32)

    cells["donor_id"] = donor
    cells["region"] = cells["cell_region_20um"].map({"core": "core", "expansion": "peri"})
    for missing in [m for m in markers if m not in have]:
        cells[missing] = np.nan          # unavailable on this donor's panel, never zero
    return cells


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--clean-only", action="store_true",
                        help="CLEAN only -- the 2026-09-14 default, which is how the deployed "
                             "build ended up with ten donors and no Aab+ case")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=not args.clean_only)
    # The exclusion list has exactly one home. This lane derives from it and never restates it.
    donors = filter_eligible_donors(donors)

    markers = feature_markers()
    scale = S.load_scale()
    markers = [m for m in markers if m in scale.index]
    print(f"{len(markers)} features: {', '.join(markers)}\n")

    meta = pd.read_csv(S.OUT / "app/donor_metadata.csv")[["case_id", "donor_status"]]
    meta["case_id"] = meta["case_id"].astype(str)
    status = meta.set_index("case_id")["donor_status"]
    runs = pd.read_csv(S.REPO / "data/donor_acquisition_runs.csv").rename(
        columns={"donor_id": "case_id"}).drop_duplicates("case_id")
    runs["case_id"] = runs["case_id"].astype(str)
    stain = runs.set_index("case_id")["stain_group"]

    frames, per_donor = [], {}
    print(f"{len(donors)} donors: {', '.join(donors)}\n")
    print(f"  {'donor':>6} {'status':>5} {'cells':>9} {'core':>8} {'peri':>8} {'islets':>7} "
          f"{'markers':>7}  source generation")
    for donor in donors:
        frame = donor_frame(donor, markers, scale)
        frame["donor_status"] = status[donor]
        frame["stain_group"] = stain.get(donor, "unknown")
        counts = frame["region"].value_counts().to_dict()
        have = sorted(m for m in markers if frame[m].notna().any())
        source = S.OUT / f"cells/donor_id={donor}/data.parquet"
        written = _dt.datetime.fromtimestamp(source.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        per_donor[donor] = {
            "donor_status": str(status[donor]), "stain_group": str(stain.get(donor, "unknown")),
            "cells": int(len(frame)), "core": int(counts.get("core", 0)),
            "peri": int(counts.get("peri", 0)),
            "islets": int(frame["islet_key"].nunique()),
            "markers_available": len(have),
            "markers_absent_on_this_panel": sorted(set(markers) - set(have)),
            "source": str(source.relative_to(S.REPO)), "source_written": written,
        }
        print(f"  {donor:>6} {str(status[donor]):>5} {len(frame):>9,} "
              f"{counts.get('core', 0):>8,} {counts.get('peri', 0):>8,} "
              f"{frame['islet_key'].nunique():>7,} {len(have):>7}  {written}", flush=True)
        frames.append(frame)

    out = pd.concat(frames, ignore_index=True)
    keep = ["donor_id", "object_id", "islet_key", "region", "phenotype", "donor_status",
            "stain_group", "X_centroid", "Y_centroid", "signed_dist_islet_um"] + markers
    out = out[keep]

    absent = {m: int(out[m].isna().sum()) for m in markers if out[m].isna().any()}
    if absent:
        print(f"\nmarkers unavailable on some donor's panel (left Unavailable, not zero): {absent}")

    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "singlecell_features.parquet"
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), target, compression="zstd")
    (OUT / "singlecell_features.json").write_text(json.dumps({
        "markers": markers,
        "n_cells": int(len(out)),
        "donors": list(donors),
        "per_donor": per_donor,
        "typing_generation": "data/development/typing/cells (current); NOT typing/subpanel",
        "scaling": "apply_scale floor/span, upper clip only -- no lower clip",
        "population": "islet core and peri-islet (20 um) cells, Excluded removed",
        "held_out": list(PANEL_IDENTITY),
        "built": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "run_id": S.RUN_ID,
    }, indent=2))
    print(f"\n{len(out):,} cells x {len(markers)} features -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
