"""Inputs the subpanel classification reads, and the scale it decides on.

Every number here was reproduced against the classification already on disk before this
module was written, so it is a reconstruction rather than a reimplementation:

    celltype/uniform/donor_id=<d>/data.parquet      UniFORM-aligned per-cell values
      -> apply_scale(floor, span)                   broadtype/scratch/cosine_scale.parquet
      -> reference level = p99 of max(subpanel) over cells in NO classified object
                           and outside every islet

POOLED IS NOT THE METHOD.  Those four levels reproduce to four decimals when the reference
population is pooled over all 20 donors, which is how the OLD table was built, and pooling
survives only so ``subpanel.py --gate`` can reproduce it.  Checked 2026-09-12 against
``subpanel/subpanel_log.json``:

    hormone 0.1010 · Nerve 0.2272 · Lymphatic 0.1776 · Blood vessel 0.3073

Fitted PER DONOR -- the default, and the method -- the same reference population moves the
hormone level 19.5x across the cohort and separates almost perfectly by acquisition batch
(batch 1 median 0.0418, batch 2 0.3167), which is aliased with disease status.  See
``subpanel.fit_levels``.

The scale is cohort-wide by the celltype lane's own design (``run_step3.build_scale``):
one floor per marker at the aligned background mode, one top at the cohort p99.9.  The
floor is the half that UniFORM did not flatten for twelve markers -- see the plan's
"Raised, not resolved here" -- and it is read, not recomputed, so this lane cannot drift
from the values the classification was made on.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

# The one table of superseded artifacts, and the refusal that comes with it.  It is imported
# here, rather than by each of the thirty-three modules that read through this one, because the
# paths below are what those modules compose everything else from: publishing OUT and RUN as
# guarded paths attaches the refusal to the ARTIFACT instead of to the readers, and a reader
# added tomorrow is covered without knowing this module exists.  See superseded.py's own
# docstring for the two severities and the single allowed reader.
_HERE = str(pathlib.Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import superseded as SUPERSEDED  # noqa: E402
from superseded import (  # noqa: E402  -- re-exported so `S.refuse_if_superseded(p)` works
    SupersededArtifactError,
    guarded,
    is_superseded,
    refuse_if_superseded,
)

REPO = pathlib.Path(__file__).resolve().parents[4]

#: The HISTORICAL run directory, not a live run key.  The content-addressed run framework was
#: removed on 2026-09-16 -- there is no ``run_id`` any more, and the pipeline writes to stable
#: roots.  This particular tree is kept because the 20-donor cohort output still lives ONLY
#: there; nothing writes to it, and a donor re-exported since is served by :data:`REVISED`
#: through :func:`cell_source` instead.
RUN_ID = "30f6689e4d27d295f890"

#: GUARDED.  Reading anything under the historical run tree prints the full superseded notice
#: to stderr once per process -- what it is, when the run framework was removed, what replaced
#: it, and what to read instead.  It ANNOUNCES rather than refuses because the 20-donor cohort
#: expression output exists nowhere else; refusing it would stop the typing lane, which is a
#: method decision and not this module's to make.  Loud, not silent, is the requirement.
RUN = guarded(REPO / "core/phenocycler/data/runs" / RUN_ID)
ALIGNED = REPO / "data/development/celltype/uniform"
#: NOT superseded, though it sits inside the withdrawn producer comparison's tree: this is the
#: cohort floor/top per marker that the classification on disk was made on, and
#: `superseded.EXEMPT` names it so no prefix can sweep it up.
SCALE_CSV = REPO / "data/development/broadtype/scratch/cosine_scale.parquet"

#: GUARDED, and this is the one that matters.  Every module in this lane composes its paths from
#: `S.OUT`, so every one of them now refuses the superseded artifacts under it without calling
#: anything: `S.OUT / "umap/refinement_log.json"` raises on open, naming the 2026-09-11 file, the
#: fact that nothing in the repo produces it, and what to read instead.  `S.OUT / "cells"`,
#: `S.OUT / "subpanel/subpanel_log.json"` and the rest are untouched -- the registry names
#: artifacts, not directories it happens to be able to reach.
OUT = guarded(REPO / "data/development/typing")

#: The v3 QuPath export root: what `IO60panc2nd/scripts/v3/run_v3.sh --out` writes for the 25
#: images re-segmented on the Tissue annotations redrawn 2026-09-30.  ONE root for every object
#: this lane reads, and it is NOT `results/Geojson`: that folder holds the v2 objects, drawn
#: against the old annotations, and the old resolver took "the newest file that matches" -- so a
#: v3 file dropped beside them would have been mixed into a v2 build, or the reverse, by mtime.
#: Nothing v3 is ever written into `results/Geojson`.
EXPORT = pathlib.Path.home() / "IO60panc2nd/results/v3/export"

#: Structure objects, one file per image and class: `<Class>/<Class>__<image tag>.geojson`, written
#: by the driver's step 7 (the maintainer's tiled pixel classifiers, unchanged).
GEOJSON = EXPORT / "structures"

#: Islet outlines for this lane, one file per image, written by the export adapter
#: (`IO60panc2nd/scripts/v3/islet_export_adapter.py`) from the reviewed QuPath objects: features
#: named `Islet_<N>` with classification `Islet`, and `Islet_<N>_exp20um` where an expansion is
#: declared -- the shape `islet_membership.outlines` and `build_islet_geometry.islet_outlines`
#: have always read.
#:
#: Read from the export where it LANDS, never a hand-copied folder inside the repo: `REPO/data/Islets`
#: was that copy and went stale on 2026-09-15 (the re-exported outlines for 6523 and 6566 sat in
#: the export directory for three hours while the repo held the July versions).
ISLETS = EXPORT / "typing_islets"

#: The rebuild of the re-exported batch-2 donors. Authoritative for them: the historical run tree
#: still holds the defective export whose islet detections were duplicated.
REVISED = OUT / "revised_cells"

#: Where the corrected typing is written.
TYPING = OUT / "subpanel_rebuilt"


def cell_source(donor: str) -> tuple[pathlib.Path, str]:
    """(path, which) for a donor's cell geometry -- the revised rebuild if it exists, else the run.

    Every consumer must go through this. `build_cells.py` read the run unconditionally until
    2026-09-15, which on the seven re-exported donors meant taking coordinates and islet numbers
    from the very export they were re-exported to replace.
    """
    revised = REVISED / f"donor_id={donor}/data.parquet"
    if revised.exists():
        return revised, "revised"
    return RUN / f"expression/donor_id={donor}/data.parquet", "run"


def cell_area_source(donor: str) -> pathlib.Path:
    """The file a donor's QuPath `cell_area` is read from -- the SAME source `cell_source` names.

    The area readers took it from the historical run tree unconditionally until 2026-10-01, which
    for a donor whose cells come from a revised export meant an area from the superseded
    segmentation, joined by `object_id` and silently NaN for every cell the run never had.
    """
    path, which = cell_source(donor)
    if which == "revised":
        return path
    return RUN / f"cells/donor_id={donor}/data_0.parquet"


def cell_areas(donor: str) -> pd.DataFrame:
    """(object_id, cell_area), one row per cell, from :func:`cell_area_source`."""
    import pyarrow.parquet as pq
    return pq.read_table(cell_area_source(donor), columns=["object_id", "cell_area"]).to_pandas()


def cell_outlines(donor: str) -> pathlib.Path:
    """The per-cell outline GeoJSON of one donor, as the QuPath manifest records it.

    The manifest is where an image's geometry file is declared (`cell_geojson.path`), so it is
    asked rather than a folder being globbed for a name: the glob was
    `data/redsea_scratch/geojson/cells__*<stem>*.geojson`, which is where 21 of the 23 v2 images
    keep theirs and where none of the v3 images do.
    """
    import json
    raw = json.loads((REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    found = [pathlib.Path(r["cell_geojson"]["path"]) for r in raw["images"]
             if str(r["donor_id"]) == str(donor)]
    if len(found) != 1:
        raise FileNotFoundError(f"{donor}: {len(found)} images in the QuPath manifest")
    return found[0]

#: context -> the markers that answer that context's question.  Vessel is Caveolin/CD34 and
#: NOT CD31: the objects were drawn by `Cavdil1200` / `CD34dil1430`, so CD31 was never what
#: made them.  Exocrine is three markers because EpCAM (0.9-3.1x) and Keratin-5 (1.0-1.3x)
#: were measured not to separate the bulk from any removed class.
SUBPANEL: dict[str, tuple[str, ...]] = {
    "islet": ("INS", "GCG", "SST", "CD99"),
    "Nerve": ("B3TUBB", "CD56"),
    "Lymphatic": ("Podoplanin",),
    "Blood vessel": ("Caveolin", "CD34"),
    "exocrine": ("Ker8_18", "E_cadherin", "Pan_Cytokeratin"),
}

#: The three hormones, read together against the non-islet reference level.  CD99 is in the
#: islet subpanel but cannot hold a cell endocrine: when it wins the argmax the cell is not
#: endocrine and continues down the ordinary chain.
HORMONES: tuple[str, ...] = ("INS", "GCG", "SST")

#: The quantile of the reference population that sets every structure's level.
REFERENCE_QUANTILE = 99.0

#: Islets of STRICTLY MORE than this many cells are "the larger islets" the maintainer's
#: CD99/size/shape rule is scoped to.  Not a second definition: it is the same boundary
#: ``islet_identity_mixture.py`` takes as ``--min-islet-cells`` (default 10, strictly greater),
#: and that module's own reason for it is quoted there -- a two-cell classifier fragment has no
#: architecture to impose a shape on its cells.
LARGE_ISLET_MIN_CELLS = 10

#: The size strata, named exactly as the 20-donor measurement of 2026-09-18 reported them.  The
#: measurement is the reason they exist: nucleus elongation's mixture signature runs
#: 1.519 / 1.139 / 1.657 / 2.089 across these four and is unanimous across donors only in the
#: last, while CD99 carries no mixture signature in any of them and a 2.07 / 2.63 / 2.34 / 3.85x
#: median ratio instead.  So the strata are what tells a shape criterion where it may fire.
SIZE_STRATA: tuple[str, ...] = ("1", "2-3", "4-10", ">10")
LARGE_STRATUM = ">10"


def size_stratum(islet_cells) -> np.ndarray:
    """Stratum label per cell, from the number of cells its islet contributes to the analysis.

    ``islet_cells`` must be counted over the cells the analysis ADMITS -- islet core,
    non-structural, not ``Excluded`` -- which is the count ``islet_identity.load`` already
    derives and says why.  A cell outside every islet has no stratum and is given ``""``;
    an islet that contributes nothing lands in the smallest stratum rather than a fifth one,
    because a size of zero is a size, not a separate kind of thing.
    """
    n = np.asarray(islet_cells, dtype=float)
    out = np.full(n.shape, "", dtype=object)
    known = np.isfinite(n)
    out[known] = SIZE_STRATA[0]
    out[known & (n >= 2) & (n <= 3)] = SIZE_STRATA[1]
    out[known & (n >= 4) & (n <= LARGE_ISLET_MIN_CELLS)] = SIZE_STRATA[2]
    out[known & (n > LARGE_ISLET_MIN_CELLS)] = LARGE_STRATUM
    return out


#: The mirror of :data:`REFERENCE_QUANTILE`, and derived from it rather than chosen again.
#:
#: A structure's level is the 99th percentile of a population that CANNOT contain the structure:
#: "1 % of ordinary tissue reaches this".  Endocrine identity inside an islet is the same
#: statement upside down -- the reference population IS the thing being detected, so the level
#: is the 1st percentile of it: "1 % of confirmed endocrine cells fall this low".  Same 1 %
#: convention, same direction of error, no second threshold to defend.
ENDOCRINE_QUANTILE = 100.0 - REFERENCE_QUANTILE

#: Fewest hormone-confirmed cells a donor x stratum must contribute before a level is fitted
#: from them.  Taken from ``islet_identity.main``'s existing reference guard (``len(ref) < 500``)
#: rather than invented here.  Below it the level is UNAVAILABLE -- not defaulted, not borrowed
#: from another stratum, and not silently replaced by the cohort's.
MIN_REFERENCE_CELLS = 500


def load_scale() -> pd.DataFrame:
    """The cohort floor/top/span per marker, as the classification decided on it."""
    return pd.read_parquet(SCALE_CSV).set_index("marker")


def apply_scale(values: np.ndarray, markers: list[str], scale: pd.DataFrame) -> np.ndarray:
    """(cells x markers) -> [0, 1], clipped at the measured floor and the cohort top.

    Identical to ``celltype/score.apply_scale``; duplicated rather than imported so this
    lane's inputs cannot change under it when that lane is edited.
    """
    lo = scale.loc[markers, "floor"].to_numpy(float)
    span = scale.loc[markers, "span"].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    return np.clip((np.asarray(values, dtype=np.float32) - lo) / span, 0.0, 1.0)


def reference_level(pooled: np.ndarray, quantile: float = REFERENCE_QUANTILE) -> float:
    """The level a structure's own markers must clear, from the tissue it sits in.

    A non-finite level would silently classify nothing -- ``value >= nan`` is False for every
    cell -- which is exactly how eight NaN cells in 10.4 million disabled a whole tier and
    left every downstream number looking plausible.  It fails loudly instead.
    """
    level = float(np.nanpercentile(pooled, quantile))
    if not np.isfinite(level):
        raise ValueError(f"reference level is not finite ({level!r}) -- refusing to classify")
    return level


def endocrine_level(values: np.ndarray, quantile: float = ENDOCRINE_QUANTILE,
                    min_cells: int = MIN_REFERENCE_CELLS,
                    floor: float = 0.0) -> float | None:
    """The level a cell must CLEAR to be endocrine, fitted from cells that already are.

    The argument is the hormone-confirmed cells of ONE donor and ONE size stratum -- never the
    cohort, never a donor pooled across strata.  Pooling was measured to be the wrong thing on
    both axes: the same reference population fitted per donor moves the hormone level 19.5x
    across the cohort and separates almost perfectly by acquisition batch (see
    :func:`subpanel.fit_levels`), and the 2026-09-18 measurement puts CD99's confirmed/undecided
    median ratio at 2.07x in one-cell islets against 3.85x above ten, so one level for both
    strata is a blend of two populations exactly as the pooled hormone level was.

    Returns ``None`` -- unavailable -- in three cases, and an unfittable stratum has to be
    visible as one, because a level that does not exist must not be invented:

      * fewer than ``min_cells`` confirmed cells are on hand;
      * none of them carries the measurement;
      * the level lands AT the feature's floor.

    The third is not hypothetical and it is the more dangerous.  ``apply_scale`` clips at the
    measured background, so a marker at or below background reads exactly 0.0, and on donor 6534
    20.3 % of islet cells sit there -- which puts the 1st percentile of its confirmed endocrine
    cells at 0.0 too.  ``value >= 0.0`` is then true of every cell that was measured at all, so
    the "level" admits the entire islet while looking like a fitted number.  That is the same
    failure as the eight NaN cells that disabled a whole tier and left every downstream number
    plausible, so it fails the same way: refused, named, and handed back to the caller.

    ``floor`` is the feature's own lower bound -- 0.0 for a scaled marker, 1.0 for an elongation
    that is a ratio of diameters.
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < min_cells:
        return None
    level = float(np.percentile(finite, quantile))
    if not np.isfinite(level):
        raise ValueError(f"endocrine level is not finite ({level!r}) -- refusing to classify")
    if level <= floor:
        return None
    return level
