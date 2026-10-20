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
   feature set is the markers the subpanel classification itself uses to decide what a cell is:
   the islet subpanel, each structure's own markers, the exocrine trio, the immune definitions and
   the supporting tier.  Ordering cells on the evidence that named them is the defensible choice
   when the published feature list does not exist in this panel.  Thirty-four of the thirty-eight
   are read -- from ``subpanel/subpanel_log.json``, ``subpanel_io.SUBPANEL["exocrine"]`` and
   ``refine.OWN_MARKERS`` -- and **four are chosen here** by name in :func:`feature_markers`
   (CD11b, IDO1, HLA_DR, Granzyme_B).  That distinction was recorded in
   ``docs/TRAJECTORY_PARAMETER_AUDIT_2026-09-15.md`` and is stated here rather than left to a
   docstring that said "rather than chosen here".

**WHERE THE DEFINING-MARKER TABLE COMES FROM, 2026-09-19.**  It comes from ``refine.OWN_MARKERS``,
which is committed code.  Until 2026-09-19 it was parsed out of
``data/development/typing/umap/refinement_log.json`` -- a 2026-09-11 file that **nothing in this
repository produces**, whose second step is recorded ``"accepted": false`` and was being walked
anyway, and whose counts are a 320,000-cell subsample against a 25,230,109-cell cohort.  That was
superseded work serving as a live input to a current build.  The two tables were compared marker
for marker before the change and are **identical**, so the feature set is unchanged at 38 markers;
that is why the read could be closed while closing it still costs nothing.  The orphan now
**refuses to be read at all** -- importing ``refine`` installs the guard, so a read of it anywhere
in this process raises ``refine.SupersededRefinementLog``.  The file is not deleted.

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

    python core/multiplex5000/scripts/typing/trajectory_singlecell_features.py
"""
from __future__ import annotations

import argparse
import datetime as _dt
import importlib.util
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cohort as COHORT
import subpanel_io as S

sys.path.insert(0, str(S.REPO / "core/multiplex5000"))
from phenocycler.cohort import filter_eligible_donors  # noqa: E402  -- the one exclusion authority


def _load_refine():
    """``scripts/typing/refine.py``, loaded BY PATH and never by name.

    ``import refine`` is ambiguous in this repository -- ``scripts/restain/refine.py`` is a
    different module that the restain lane imports under exactly that name -- so a plain import
    here would depend on whose directory reached ``sys.path`` first and on what is already in
    ``sys.modules``.  Loading by file location cannot resolve to the wrong file.

    There is no fallback.  If this fails the build stops: the alternative it used to have was
    the 2026-09-11 orphan log, and falling back to superseded work is the defect, not the
    recovery.  Importing it also installs ``refine``'s process-wide refusal of that orphan.
    """
    spec = importlib.util.spec_from_file_location("typing_refine", HERE / "refine.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("typing_refine", module)
    spec.loader.exec_module(module)
    return module


REFINE = _load_refine()

#: Each class's own defining markers.  AUTHORITY: ``refine.OWN_MARKERS`` in
#: ``core/multiplex5000/scripts/typing/refine.py`` -- committed code, not an artifact on disk.
#: Taken from the module rather than restated, so this file cannot hold a stale copy of it.
OWN_MARKERS = REFINE.OWN_MARKERS

OUT = S.OUT / "trajectory"

#: Markers that encode which panel a donor was on rather than which cell it is. Held out of any
#: similarity representation, exactly as the coarse embedding holds them out.
PANEL_IDENTITY = ("CD38", "b_Catenin1", "bCatenin1", "CD57")


#: Markers named here and read from nothing.  Four of the thirty-eight features are a choice
#: made in this file, not evidence the classification reports using; they are listed separately
#: so a reader can see which is which, and they are recorded as such in the sidecar.
CHOSEN_HERE = ("CD11b", "IDO1", "HLA_DR", "Granzyme_B")


def feature_markers() -> list[str]:
    """The markers the classification uses: read from its own parameter sidecar and from code.

    ``subpanel/subpanel_log.json`` HAS NO PRODUCER EITHER, and saying otherwise here was wrong.
    Measured 2026-09-19: nothing in the repository writes it.  ``subpanel.py:318`` is
    ``def _log(): return json.loads((OUT / "subpanel_log.json").read_text())`` -- a READ.  It is
    688 bytes, written 2026-09-11 14:10, the same generation as the ``umap/refinement_log.json``
    orphan this function stopped reading and 32 minutes later than it, and its top-level key is
    still ``"bars"``, the vocabulary the root CLAUDE.md rejects in favour of *reference level*.

    So 13 of the 38 markers below -- CD107a, CD163, CD206, CD39, CD4, CD8, CD99, ICOS, M2Gal3,
    PD_1, TCF_1, TOX, iNOS -- still come only from an unrebuildable Sep-11 artifact.  15 come
    from committed code.  That is REPORTED, not fixed here: giving this file a producer, or
    moving its contents into code, changes what the trajectory features are and is the
    maintainer's, not this function's.

    The defining-marker table is NOT read from disk: it comes from :data:`OWN_MARKERS`, i.e.
    ``refine.OWN_MARKERS``.  See the module docstring for the 2026-09-11 orphan this replaced
    and the comparison that showed the two tables identical.
    """
    log = json.loads((S.OUT / "subpanel/subpanel_log.json").read_text())
    markers: list[str] = list(log["islet_subpanel"])
    for panel in log["subpanels"].values():
        markers += list(panel)
    markers += [log["hormones"][k] for k in log["hormones"]]
    markers += list(log["bottom_tier_immune"])
    markers += list(S.SUBPANEL["exocrine"])
    # Each class's OWN defining markers -- the immune definitions among them -- from refine.py,
    # the committed authority.  The old reader walked EVERY step of umap/refinement_log.json,
    # including step 2, which that file records as "accepted": false.  There are no steps here.
    for panel in OWN_MARKERS.values():
        markers += list(panel)
    markers += list(CHOSEN_HERE)
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
        "own_markers_authority": (
            "core/multiplex5000/scripts/typing/refine.py OWN_MARKERS (committed code). "
            "NOT data/development/typing/umap/refinement_log.json, which was the source until "
            "2026-09-19, has no producer in the repository, and is now refused on read."),
        "own_markers": {k: list(v) for k, v in OWN_MARKERS.items()},
        "markers_chosen_in_this_file": list(CHOSEN_HERE),
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
