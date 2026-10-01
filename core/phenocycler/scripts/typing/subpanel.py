#!/usr/bin/env python3
"""The subpanel classification: context from geometry, then one question per context.

WHY THIS FILE EXISTS. `subpanel/subpanel_calls.parquet` -- 25,150,537 rows, written 2026-09-11,
now SUPERSEDED -- had **no producer**. It was made in a session and the script was never saved: an
exhaustive search of the working tree, every branch, the full git history and the other checkouts
finds only `subpanel_io.py` (its inputs), `build_cells.py` (a reader) and prose. So the typing
could not be re-run, which meant no correction to it could ever be applied.

This is that producer, written from `RECONSTRUCTION.md` -- which was itself recovered by
reproducing the calls from their inputs and checked against them, not inferred from prose.

THE ISLET CHAIN BELOW IS ON DISK. Measured 2026-09-26: the stage roots `typing/subpanel_rebuilt/`
and `typing/cells/` hold 25,230,109 rows over 20 donors and were both rewritten on 2026-09-19 with
it, and 6374 alone was retyped 2026-09-26 on its re-exported structure objects -- every partition
carries `islet_route` (cd99 408,934, elongation 29,047, hormone 190,571) and
`islet_stratum`. The row count equals the 2026-09-15 pre-chain build's, so read `islet_route`, not
a count, to tell which build is on disk.

THE METHOD, unchanged from what was recovered:

    context comes from geometry; each context is asked only what its own markers can answer;
    every level comes from a REFERENCE POPULATION, never from a constant.

  reference population   cells in no classified object and outside every islet -- the ordinary
                         tissue a structure sits in. It excludes the objects and the islets so it
                         cannot contain the thing being detected. A structure's level is its p99:
                         "1 % of ordinary tissue reaches this".
  levels                 fitted PER DONOR (`--levels per-donor`, the default). The POOLED four --
                         hormone 0.1010 · Nerve 0.2272 · Lymphatic 0.1776 · Blood vessel 0.3073 --
                         belong to the OLD table and survive only so `--gate` can reproduce it.
                         Pooling was measured wrong on 2026-09-15: see `fit_levels`.
  scale                  READ from the celltype lane, never recomputed -- `build_scale` draws a
                         seeded 200,000-cell sample per donor, so recomputing it would drift.

THE ISLET CHAIN, as the maintainer ruled it (2026-09-18). Inside an islet outline:

    1. STRUCTURE     in a classified object, clearing that object's level, and the object's own
                     marker is the cell's `top_identity_marker`. The classified objects separate
                     vessels and nerves, inside islets as well as outside.
                     THAT COLUMN IS NOT "THE TOP EXPRESSED MARKER" -- it is the argmax of a
                     per-marker-normalised value, and a broadly expressed marker wins it. The
                     measurement, the two alternative readings and the one constant that says
                     which of them correction C reads are on `IDENTITY_MARKERS` below.
    2. EXOCRINE      inside an islet is an ANNOTATION ERROR by definition -- exocrine markers get
                     into an islet annotation only by mistake. `Excluded`, with the reason
                     written on the row.
    3. CD99          the key differentiator between endocrine and non-endocrine, at a level
                     fitted PER DONOR **AND PER SIZE STRATUM** from that donor-and-stratum's own
                     hormone-confirmed cells -- the 1 % convention upside down, since here the
                     reference population IS the thing being detected. A FLOOR GUARD refuses a
                     level at or below the feature's own floor and returns it UNAVAILABLE. The
                     donor's ordinary-tissue CD99 level is carried alongside as
                     `donor_reference_level` and REPORTED, never substituted for the per-stratum
                     level. See `fit_islet_levels` and `S.endocrine_level`.
    4. ELONGATION    a second ROUTE that can only ADD, in the large stratum only (islets of
                     STRICTLY more than ten admitted cells). ITS LEVEL IS NOT SETTLED -- see the
                     note below.
    5. ENDOCRINE     -> the hormone argmax names Beta / Alpha / Delta. At least one hormone must
                     have been READ first: argmax over three NaNs returns index 0, "Beta", every
                     time.
    6. NOT ENDOCRINE -> the ordinary tissue chain, unchanged. That population is the analysis.

WHY CD99 SITS AT A LEVEL AND ELONGATION AT A SHAPE. Measured over 20 donors, 2026-09-18, with
spread(undecided)/spread(confirmed) by islet size 1 / 2-3 / 4-10 / >10:

    nucleus_elongation    1.519 (1 of 20 donors >1) | 1.139 (5/20) | 1.657 (16/20) | 2.089 (20/20)
    nucleus_circularity   0.988 (0/20) | 1.078 (3/20) | 1.329 (14/20) | 1.529 (20/20)
    CD99                  1.021 (1/20) | 0.483 (0/20) | 0.665 (4/20) | 0.460 (1/20)
    CD99 confirmed/undecided MEDIAN ratio   2.07x | 2.63x | 2.34x | 3.85x

Nucleus shape carries a MIXTURE signature -- the undecided cells are WIDER than the confirmed
ones, which is what two populations in one group look like -- and it is unanimous only above ten
cells per islet and noise below four. CD99 carries no mixture signature anywhere; the undecided
cells are NARROWER. What it carries is 2.1-3.9x LESS of the marker. So CD99 is a LOCATION
discriminator and belongs at a level, and elongation is a SHAPE discriminator.

AND THAT IS WHY ELONGATION'S LEVEL IS AN OPEN QUESTION, STATED HERE RATHER THAN HIDDEN IN THE
CODE. A level can only separate populations that differ in LOCATION, and on elongation the
confirmed and undecided medians are the same to within 5 % (ratio 0.95-1.02 at every islet
size). The threshold currently fitted lands a few percent above 1.0 -- the floor a max/min
diameter ratio cannot go below -- and on 6534 it calls 22,022 of 88,736 large-islet cells
endocrine on shape alone, more than CD99 itself moves. That inverts the rule it implements:
CD99 is the key differentiator and elongation is the additional one. `islet_identity_mixture.py`
exists for exactly this -- a dispersion signal needs a mixture test, not a threshold -- and it was
run for the first time on 2026-09-18 over 20 donors (output at `typing/islet_identity/`): in the
large stratum, 14 of 20 donor rows are testable, and CD99 separates far more strongly than
nucleus elongation does. THAT DOES NOT CLOSE IT. The test asks one threshold against a
cohort-pooled null per stratum, and a single stratification is not expected to hold across donors
whose islets differ in shape by disease. So the elongation route stays provisional, its count is
quoted separately (`n_by_elongation_only` in the audit), never pooled into "endocrine", and its
WEIGHT IS OPEN FOR THE MAINTAINER.

TWO THINGS THIS CHANGES FROM THE TABLE ON DISK, both deliberate, both behind a flag:

  --membership geometry  Islet membership from the islet OUTLINES, not QuPath's `Parent`. The
                         Parent labels cells `peri` that sit up to 74.6 um inside a core outline;
                         because a false `in_islet` means the cell is never asked the islet
                         question, 11,556 islet-core cells fell through to Exocrine, 5,250 of them
                         hormone-positive. `--membership parent` reproduces the old behaviour and
                         is the DEFAULT, because stage 1 is a reconstruction gate.

                         MIND THE MIXED DEFAULTS. `--levels` defaults to `per-donor` (the method)
                         while `--membership` defaults to `parent` (the gate), so a bare
                         `subpanel.py` writes the stage root with the membership `context.py`
                         says was measured wrong and must not be reinstated. A real rebuild is
                         `--membership geometry --levels per-donor`; `--gate --levels pooled
                         --membership parent` is the reconstruction gate. Neither is what you
                         get by typing nothing.
  --emit-detail          Carries `immune_definition` and `resolved_by` as columns. The table on
                         disk collapsed both into the parent label and threw them away.

THE GATE COMES FIRST. `--gate` rebuilds with `--membership parent` and compares row-for-row on
(donor_id, object_id, cell_class) against the table on disk. Until that matches, a difference
cannot be attributed: it is either a reconstruction bug or a correction working, and there is no
way to tell which. Two runs were lost to exactly that confusion before.

WHAT IS NOT RECOVERED, and is reported rather than hidden: which cells the supporting tier is
eligible to claim. Applying it to every unclaimed cell reproduces 99.998 % on 6442 and 115 but
99.59 % on 6436, so an eligibility condition remains unrecovered. This runs the permissive form
and prints the per-donor residual every time, because a reconstruction that quietly papers over
its own known gap is worth less than one that names it.

    python core/phenocycler/scripts/typing/subpanel.py --gate --donors 6442 6436 115
    python core/phenocycler/scripts/typing/subpanel.py --membership geometry --emit-detail
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "celltype"))
import context as CX
import islet_identity as ID
import islet_membership as MEMBERSHIP
import lineages as LIN
import rules as RULES
import subpanel_io as S

OUT = S.OUT / "subpanel"
WORK = S.OUT / "scratch/subpanel_context"

#: Strongest first. Recovered from the calls rather than assumed -- see RECONSTRUCTION.md.
STRUCTURES = ("Blood vessel", "Nerve", "Lymphatic")

#: The marker each pixel classifier was actually BUILT ON -- what drew the object.
#:
#: Deliberately not the same as `S.SUBPANEL`, which is what ANSWERS the question once a cell is
#: in the object. Nerve's subpanel reads CD56 as well, but the object came from
#: `B3Tub_Ero0.5_Thr600` -- B3TUBB alone. Vessel objects came from `Cavdil1200` and `CD34dil1430`,
#: so Caveolin and CD34 and NOT CD31, which never drew anything. Lymphatic from
#: `PDPNmorphClos0-5_Thr300`.
CLASSIFIER_MARKER: dict[str, tuple[str, ...]] = {
    "Nerve": ("B3TUBB",),
    "Lymphatic": ("Podoplanin",),
    "Blood vessel": ("Caveolin", "CD34"),
}

#: The `crossed_frac_*` column that measures how much of a cell each object covers.
#:
#: `build_crossed_fraction.py` (correction D) writes
#: ``crossed_frac_<col> = area(cell outline INTERSECT that structure's objects) / area(cell)``
#: from the QuPath cell polygons and the same objects the classification uses.  Its spelling is
#: the short one; `STRUCTURES` uses the display name, so the two are mapped here rather than
#: guessed at a call site.
STRUCTURE_AREA_COLUMN: dict[str, str] = {
    "Blood vessel": "crossed_frac_vessel",
    "Nerve": "crossed_frac_nerve",
    "Lymphatic": "crossed_frac_lymphatic",
}

#: The share of a cell's own area that must lie inside an object for the cell to be a CANDIDATE
#: for that object's class.
#:
#: "The majority of the cell boundaries are within the object" -- the maintainer's rule, step 2.
#: A half is the majority and nothing here is tuned: the number is what the word means.
MAJORITY_AREA = 0.5

#: Markers that can carry an identity, from `lineages.LINEAGE_MARKERS` -- "a statement of what
#: each antibody means, made before any result is seen". State markers (proliferation,
#: exhaustion, activation) are absent by that module's own design: they say what a cell is doing,
#: not what it is, so they must not be able to win a cell's top-marker argument.
IDENTITY_MARKERS: tuple[str, ...] = tuple(sorted(
    {m for markers in LIN.LINEAGE_MARKERS.values() for m in markers}))

#: WHAT `top_identity_marker` IS, AND WHY ITS NAME IS WRONG.
#:
#: It is NOT "the marker this cell expresses most". It is the argmax of the SCALED value, and
#: `S.apply_scale` maps every marker onto its OWN floor-to-p99.9 range. A marker that is broadly
#: and moderately expressed therefore sits high on its own scale for most cells and takes the
#: argmax from a marker that is sharply expressed on few. The column is the argmax of a
#: PER-MARKER-NORMALISED value, and that is the whole of what it says.
#:
#: MEASURED 2026-09-18 over the analysis-eligible cells, through `classify`'s own `values` and
#: not a standalone reading of it. Winning share of the argmax, with the winning marker's own
#: scaled p50 in brackets:
#:
#:   115  (860,863 cells)    FOXP3 21.09 % [0.0000] · Ker8_18 14.51 % [0.0000] ·
#:                           E_cadherin 11.06 % [0.0017] · EpCAM 8.80 % [0.0000] ·
#:                           CD34 6.97 % [0.0087] · CD20 4.30 % [0.0000]
#:   6442 (1,419,669 cells)  E_cadherin 21.76 % · FOXP3 21.72 % · EpCAM 9.60 % ·
#:                           Iba1 7.61 % · CD11b 7.26 % · CD34 6.99 %
#:   6534 (1,160,909 cells)  FOXP3 24.58 % · E_cadherin 10.18 % · CD34 9.71 % ·
#:                           Pan_Cytokeratin 7.75 % · CD11b 6.29 % · Ker8_18 5.46 %
#:
#: FOXP3 is in `lineages.SUPPORTING_ONLY` -- the maintainer's statement that a marker may not
#: carry a call on its own, FOXP3 because it sits above the positivity model's own boundary in
#: 90.4 % of cells -- and it is nonetheless the MODAL winner on two of those three donors. A
#: SUPPORTING_ONLY marker takes the argmax on 30.13 % (115), 31.76 % (6442) and 33.34 % (6534)
#: of cells.
#:
#: A THIRD DEFECT, MEASURED HERE AND ALSO NOT REPAIRED. The argmax does not range over these 34
#: markers at all. `classify` assembles `values` from what the SUBPANELS, the immune rules and
#: the supporting tier ask for, and the argmax then runs over whichever identity markers that
#: set happens to contain -- 28 of 34 on all three donors. CD31, Collagen_IV, Keratin_5, SMA,
#: TP63 and Vimentin are never requested, so they can never win, which removes the whole Stromal
#: lineage (SMA, Vimentin, Collagen_IV) and CD31 from a column named for identity. Read against
#: the full 34 the shares move -- on 115, Collagen_IV takes 6.30 % and Vimentin 3.66 % -- so the
#: marker set is doing part of the work the name attributes to the data. `identity_argmax_audit`
#: records `markers_available` and `markers_absent` on every donor so this stays checkable. It is
#: not repaired because widening the set would change which cells correction C admits.
#:
#: THE CONSEQUENCE, and the reason this matters now. Correction C in `_structure_calls` requires
#: the cell's entry in this column to be one of the markers that DREW the object. Read against
#: an argmax a broad marker wins, that requirement withholds structure calls on grounds that
#: have nothing to do with which structure the cell is in. Until correction C nothing read this
#: column, which is why the defect never surfaced.
#:
#: A SECOND DEFECT, REPORTED AND DELIBERATELY NOT FIXED, because fixing it would change a
#: computed value and this lane may not: `np.nanargmax(np.nan_to_num(..., nan=-1.0))` over a row
#: with NO identity marker measured at all returns index 0 -- `B3TUBB`, which is Nerve's
#: classifier marker -- instead of UNAVAILABLE. Measured 2026-09-18: 0 such cells on 115 and
#: 6442, 13 on 6534. The two alternative columns below DO return unavailable there, which is part
#: of why they disagree.

#: The three readings of the identity argmax this file computes side by side, and the ONE that
#: correction C reads. Naming them here is the point: which reading is right is a METHOD choice,
#: it belongs to the maintainer, and it should be changeable in one line rather than found inside
#: `classify`.
IDENTITY_ARGMAX_SCALED = "top_identity_marker"
IDENTITY_ARGMAX_NO_SUPPORTING = "argmax_scaled_no_supporting_only"
IDENTITY_ARGMAX_BY_RANK = "argmax_within_marker_rank"
#: The marker the cell EXPRESSES most -- the argmax of the aligned values before any scaling.
#: This one IS the rule: "if the same marker that created the object is simply higher than all
#: the rest, it gets classified", and the approved plan's "the structure's marker be the cell's
#: TOP EXPRESSED marker". The other three are normalised readings and are kept for comparison.
IDENTITY_ARGMAX_EXPRESSED = "top_expressed_marker"

#: WHICH column `_structure_calls` compares against `CLASSIFIER_MARKER`.
#:
#: THE RULE, not a choice. The maintainer's words: "IF THE SAME MARKER THAT CREATED THE OBJECT IS
#: SIMPLY HIGHER THAN ALL THE REST, IT GETS CLASSIFIED." The approved plan says the same thing --
#: "the rule asks additionally that the structure's marker be the cell's TOP EXPRESSED marker" --
#: and records that two normalised variants were tested against the calls on disk and BOTH FIT
#: WORSE (head-to-head on scaled value 97.0-98.3 %, ratio-to-level 94.8-97.8 %).
#:
#: This read `top_identity_marker` until 2026-09-19, which is the argmax AFTER `apply_scale`, so
#: it asked "highest on its own scale" -- close to the scaled head-to-head the plan had already
#: rejected. That let a broadly expressed marker take the slot from the object's own: on donor 115
#: CD56 held it on 96 % of the withheld nerve cells, and CD56 is a nerve marker in SUPPORTING_ONLY.
#: Do not move this back to a normalised column.
CORRECTION_C_COLUMN = IDENTITY_ARGMAX_EXPRESSED

#: Scaled value an immune rule's markers must clear. Recovered, not chosen.
IMMUNE_REQUIRE_MIN = 0.10

#: The supporting tier's level over `bottom_tier_immune`, cohort-wide. Recovers independently on
#: 6442 (0.8715) and 6436 (0.8716).
SUPPORTING_LEVEL = 0.8715

#: A cell inside an islet that falls through to Exocrine is not called Exocrine.
#:
#: The maintainer's reason, in substance: exocrine markers only get into an islet annotation by
#: MISTAKE, so a cell the chain can only call Exocrine while it sits inside an islet outline is
#: evidence about the annotation, not about the cell. It is ruled out of analysis and carries
#: this reason in `excluded_reason` -- a column written on every row, not only under
#: `--emit-detail`, so the count stays auditable from the table itself.
ISLET_EXOCRINE_REASON = "exocrine_inside_islet_border"

#: The shape feature the islet chain uses, and the reason it is a NUCLEUS feature.
#:
#: Cell-boundary features were measured across six donors and every class -- immune, exocrine,
#: alpha, beta, delta -- at cell circularity 0.514-0.523, a 1 % spread. Real immune and real
#: exocrine cells are not the same shape, so that is the signature of boundaries cut against
#: neighbours in packed tissue: the cell outline is a tessellation artefact and carries little
#: morphology. Nuclei are segmented from DAPI and are not tessellated. The cell features are
#: carried alongside and marked, exactly as `islet_identity.py` carries them, so the choice
#: stays checkable rather than assumed.
SHAPE_FEATURE = "nucleus_elongation"

#: Nucleus elongation is max diameter over min diameter, so its floor is 1.0, not 0.0. A level
#: fitted AT that floor would admit every cell with a nucleus -- see `S.endocrine_level`.
SHAPE_FLOOR = 1.0

#: Per-cell nucleus and cell morphology, written by `extract_qupath_fundamentals.py`.
CROSSED = S.OUT / "crossed"
FUNDAMENTALS = S.OUT / "qupath_fundamentals"

#: A cell for which not one marker was measured. It carries no evidence of any kind, so it is
#: not exocrine, not endocrine and not immune -- it is unmeasured, and inside an islet it is
#: ruled out under its OWN reason rather than swept into the exocrine-annotation-error set.
#:
#: This is not hypothetical either. Donor 6534's UniFORM-aligned table carries 205,070 rows with
#: no value for any of its 57 lineage markers, and 25,472 of its islet cells have no hormone read
#: on them at all. They must never be named Beta by an argmax over three NaNs, and they must not
#: be counted as exocrine tissue the islet annotation caught.
NO_MEASUREMENT_REASON = "no_marker_values"


def _log(): return json.loads((OUT / "subpanel_log.json").read_text())


def crossed_area(donor: str, object_ids: pd.Series) -> dict[str, np.ndarray]:
    """Per structure, the fraction of each cell's own area that object covers -- by `object_id`.

    The input is correction D's output, which measures the intersection of the QuPath cell polygon
    with the same objects the classification uses.  It is joined by `object_id` and nothing else.

    THE OVERLAP IS ASKED OF EVERY CELL IN THE TISSUE, not of the islet-associated ones.  A vessel
    runs where it runs, and whether the majority of a cell's area lies inside one is a question
    about that cell.  Correction D restricted itself to islet core and expansion until 2026-09-19;
    that was a limit of that producer's run and never a property of the rule.

    NaN IS UNAVAILABLE AND IS NOT ZERO, and the difference is load-bearing here.  A cell with a
    ROW and a 0.0 was measured: its outline was tested against that structure's objects and none
    touched it.  A cell with NO ROW was never tested.  The first is "no overlap"; the second is
    "we do not know", and a cell we do not know about cannot be made a candidate NOR silently
    passed through.  `_structure_calls` records it under its own name.
    """
    path = CROSSED / f"crossed_{donor}.parquet"
    out = {name: np.full(len(object_ids), np.nan) for name in STRUCTURES}
    if not path.exists():
        return out
    frame = pq.read_table(path).to_pandas()
    joined = object_ids.to_frame("object_id").merge(
        frame, on="object_id", how="left", validate="one_to_one")
    for name, column in STRUCTURE_AREA_COLUMN.items():
        if column in joined.columns:
            out[name] = joined[column].to_numpy(dtype=float)
    return out


def scaled_markers(donor: str, wanted: set[str], scale: pd.DataFrame,
                   object_ids: pd.Series) -> pd.DataFrame:
    """Scaled values for the markers this donor actually carries, joined to `object_ids`.

    By NAME, never by position: three panels exist and a marker absent from one donor is absent,
    not shifted. An absent marker yields a column of NaN, and `NaN >= level` is False, so a rule
    needing it simply cannot fire -- which is the correct behaviour and not zero-filling.
    """
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    have = set(pq.ParquetFile(path).schema_arrow.names)
    present = sorted(wanted & have)
    frame = pq.read_table(path, columns=["object_id"] + present).to_pandas()
    frame = object_ids.to_frame("object_id").merge(frame, on="object_id", how="left")
    values = pd.DataFrame(S.apply_scale(frame[present].to_numpy(), present, scale),
                          columns=present, index=frame.index)
    for marker in sorted(wanted - have):
        values[marker] = np.nan
    # The EXPRESSED values, before `apply_scale` put every marker on its own range. The structure
    # rule asks whether the object's own marker is "simply higher than all the rest", which is a
    # question about expression and cannot be answered from a per-marker-normalised number.
    # Carried as an attribute so no call site has to re-read the aligned table.
    raw = pd.DataFrame(frame[present].to_numpy(dtype=float), columns=present, index=frame.index)
    for marker in sorted(wanted - have):
        raw[marker] = np.nan
    values.attrs["expressed"] = raw
    return values


#: Columns the geometry half of the classification needs, whichever source supplies them.
GEOMETRY_COLUMNS = ["object_id", "X_centroid", "Y_centroid", "islet_num", "cell_region",
                    "qc_analysis_eligible"]


def donor_context(donor: str, membership: str, pixel: dict, stem: dict) -> pd.DataFrame:
    """object_id, context, cell_region for every analysis-eligible cell of one donor."""
    path, _ = S.cell_source(donor)
    run = pq.read_table(path, columns=GEOMETRY_COLUMNS).to_pandas()
    run = run[run["qc_analysis_eligible"]].reset_index(drop=True)

    if membership == "parent":
        parent = run["islet_num"].fillna("").astype(str).str.strip()
        in_islet = ((parent != "") & run["cell_region"].eq("core")).to_numpy()
        islet_num = parent.to_numpy(dtype=object)
    elif membership == "geometry":
        geo = MEMBERSHIP.membership(donor, run["X_centroid"].to_numpy(),
                                    run["Y_centroid"].to_numpy(), pixel[donor], stem[donor])
        in_islet = geo["in_islet"].to_numpy()
        # WHICH islet, from the same containment test that decided THAT it is in one. Taking the
        # number from the outlines and the membership from the outlines keeps one answer; taking
        # the number from `Parent` while the membership came from geometry would size an islet
        # from a different islet's cells.
        islet_num = geo["islet_num"].to_numpy(dtype=object)
    else:
        raise ValueError(f"unknown membership {membership!r}")
    islet_num = np.where(in_islet, islet_num, None)

    # Ask the object question of EVERY cell, islet included -- correction C. `context` keeps the
    # precedence answer; `object_context` keeps what the object actually said, which the original
    # table discarded and which is why an earlier analysis of this same lane went wrong.
    object_context = CX.object_label(run["X_centroid"].to_numpy(), run["Y_centroid"].to_numpy(),
                                     stem[donor], pixel[donor])
    context = np.where(in_islet, "islet", object_context)
    frame = pd.DataFrame({"object_id": run["object_id"], "context": context,
                          "object_context": object_context, "in_islet": in_islet,
                          "islet_num": islet_num, "cell_region": run["cell_region"]})
    # Per-structure membership, NOT the collapsed label: a cell inside both a vessel and a
    # lymphatic must be able to answer for both. See `context.memberships`.
    inside = CX.memberships(run["X_centroid"].to_numpy(), run["Y_centroid"].to_numpy(),
                            stem[donor], pixel[donor])
    for name, mask in inside.items():
        frame[f"in_{name}"] = mask
    return frame


LEVEL_CACHE = WORK / "reference_levels.json"


def fit_levels(donors: list[str], scale: pd.DataFrame, log: dict,
               membership: str = "geometry", work: pathlib.Path | None = None) -> dict:
    """Every structure's level, from the reference population pooled over the given donors.

    POOLING IS NOT THE DEFAULT. It is kept only so `--gate` can reproduce the table on disk, which
    was built this way. Measured 2026-09-15, the same reference population fitted per donor gives:

        hormone       0.0253 - 0.4936   19.5x spread   (pooled 0.1010)
        Blood vessel  0.1031 - 0.5387    5.2x          (pooled 0.3073)
        Nerve         0.0913 - 0.6631    7.3x          (pooled 0.2272)
        Lymphatic     0.0592 - 0.4175    7.0x          (pooled 0.1776)

    and the hormone level separates almost perfectly by acquisition batch -- batch 1 median
    0.0418, batch 2 median 0.3167, with only 6436 crossing. So the pooled level is a blend of two
    batch populations: 2.4x too strict for batch 1 and 3x too permissive for batch 2. Acquisition
    batch is aliased with disease status at the coarse level, so a level built that way carries
    the confound directly into the calls.

    Fitting per donor also makes the level mean what it says -- "above the ordinary tissue THIS
    structure sits in" -- since the reference population is that donor's own tissue.
    """
    work = work or WORK
    pools: dict[str, list] = {name: [] for name in ("hormone", ID.MARKER, *STRUCTURES)}
    for donor in donors:
        ctx = pd.read_parquet(work / f"{donor}__{membership}.parquet")
        tissue = ctx["context"].to_numpy() == "tissue"
        if not tissue.any():
            continue
        wanted = (set(S.HORMONES) | {m for name in STRUCTURES for m in log["subpanels"][name]}
                  | {ID.MARKER})
        values = scaled_markers(donor, wanted, scale, ctx["object_id"])
        # CD99's ORDINARY-TISSUE level is fitted here, with the others, from the same reference
        # population: the p99 of this donor's cells in no classified object and outside every
        # islet. It is NOT the level the islet chain decides on. It is carried into
        # `fit_islet_levels` as `donor_level`, recorded as `donor_reference_level` and reported
        # beside the per-stratum levels so a reader can see both numbers where a stratum could
        # not be fitted -- never substituted for them.
        for name, markers in [("hormone", list(S.HORMONES)), (ID.MARKER, [ID.MARKER])] + \
                             [(n, log["subpanels"][n]) for n in STRUCTURES]:
            present = [m for m in markers if values[m].notna().any()]
            if not present:
                continue
            pools[name].append(np.nanmax(values.loc[tissue, present].to_numpy(), axis=1))
    levels = {}
    for name, chunks in pools.items():
        if not chunks:
            raise ValueError(f"{name}: reference population is empty -- refusing to classify")
        levels[name] = S.reference_level(np.concatenate(chunks))
    return levels


def immune_rules() -> dict:
    """phenotype -> (all_of, [any_of groups], none_of), immune phenotypes only."""
    table = RULES.requirements(RULES.load_rules(), with_negatives=True)
    return {k: v for k, v in table.items() if k.startswith(("Imm-", "T-"))}


def fires(values: pd.DataFrame, rule: tuple) -> np.ndarray:
    """Does each cell satisfy one immune definition, with its declared negatives enforced.

    A marker the donor does not carry is NaN, and every comparison against NaN is False. That
    makes an all_of containing an absent marker unsatisfiable (correct -- the evidence is not
    there) and an any_of fall back to whichever members ARE present (also correct).
    """
    all_of, any_of, none_of = rule
    ok = np.ones(len(values), dtype=bool)
    for marker in all_of:
        ok &= (values[marker].to_numpy() >= IMMUNE_REQUIRE_MIN)
    for group in any_of:
        members = [m for m in group if m in values.columns]
        if not members:
            continue
        ok &= (values[members].to_numpy() >= IMMUNE_REQUIRE_MIN).any(axis=1)
    for marker in none_of:
        # A negative must be ABSENT. NaN is "not measured", which cannot demonstrate presence,
        # so it does not veto -- the alternative silently kills every rule on a partial panel.
        v = values[marker].to_numpy()
        ok &= ~(v >= IMMUNE_REQUIRE_MIN)
    return ok


def identity_argmax(values: pd.DataFrame) -> dict[str, np.ndarray]:
    """Four readings of "which identity marker speaks loudest for this cell", side by side.

    `top_expressed_marker` is THE RULE -- the marker the cell expresses most, from the aligned
    values before any scaling. The other three are per-marker-normalised readings, kept beside it
    so the disagreement stays visible and checkable.

    They are computed together and written together because THEY DISAGREE, and the disagreement
    is the finding. None of them is promoted here: `CORRECTION_C_COLUMN` says which one the chain
    reads, and it is left on the column that is on disk.

    `top_identity_marker` -- UNCHANGED, bit for bit, including its two defects.
        The argmax of the SCALED value: `S.apply_scale` put every marker on its own
        floor-to-p99.9 range, so this is "which marker is highest ON ITS OWN SCALE", not "which
        marker is highest". The full measurement is on `IDENTITY_MARKERS` above. Its second
        defect is preserved too: a cell with no identity marker measured at all gets index 0,
        `B3TUBB`, rather than nothing. Both are reported, neither is repaired, because repairing
        either would move a call and that is the maintainer's to move.

    `argmax_scaled_no_supporting_only` -- the same argmax with `lineages.SUPPORTING_ONLY`
        withheld. That frozenset is the maintainer's existing statement that M2Gal3, CD11b,
        FOXP3, CD56 and TP63 may not carry a call alone; this column is what that statement looks
        like applied to the argmax. It does NOT make the remaining argmax meaningful -- the
        scale problem is untouched and Ker8_18 and E-cadherin simply inherit the win.

    `argmax_within_marker_rank` -- the argmax over each value's RANK within its own marker,
        across this donor's analysis-eligible cells (average ranks, so the large tie at the
        scale floor shares one rank). This is the reading that asks "how far up its own
        distribution is this cell, for each marker, and which marker puts it highest", which is
        what the name `top_identity_marker` promises and does not deliver. It is a PER-DONOR
        quantity by construction -- the ranking population is that donor's cells -- so it is not
        comparable across donors and must not be pooled into one.

    UNAVAILABLE IS WRITTEN AS UNAVAILABLE in the two new columns, never as a tie-break winner:
    a cell with no member of that column's marker set measured gets "", because there is no
    answer, not because the answer is the alphabetically first marker.
    """
    present = [m for m in IDENTITY_MARKERS if m in values.columns]
    names = np.array(present, dtype=object)
    scaled = values[present].to_numpy(dtype=float)
    n = len(scaled)

    # UNCHANGED. Do not "fix" the all-NaN tie-break here -- see the docstring and the note on
    # IDENTITY_MARKERS. Anything written into this column changes what correction C withholds.
    top = names[np.nanargmax(np.nan_to_num(scaled, nan=-1.0), axis=1)]

    keep = [i for i, m in enumerate(present) if m not in LIN.SUPPORTING_ONLY]
    no_supporting = np.full(n, "", dtype=object)
    if keep:
        sub = scaled[:, keep]
        have = np.isfinite(sub).any(axis=1)
        no_supporting[have] = names[keep][
            np.nanargmax(np.nan_to_num(sub[have], nan=-1.0), axis=1)]

    ranks = np.full(scaled.shape, np.nan, dtype=np.float32)
    for j in range(scaled.shape[1]):
        ranks[:, j] = pd.Series(scaled[:, j]).rank(pct=True, na_option="keep").to_numpy()
    by_rank = np.full(n, "", dtype=object)
    have = np.isfinite(ranks).any(axis=1)
    by_rank[have] = names[np.nanargmax(np.nan_to_num(ranks[have], nan=-1.0), axis=1)]

    # The rule's reading: the marker this cell EXPRESSES most, from the values before scaling.
    # A cell with no identity marker measured is UNAVAILABLE and gets "", never a tie-break --
    # index 0 is `B3TUBB`, Nerve's own classifier marker, so a tie-break there would admit a cell
    # to Nerve on nothing at all.
    raw = values.attrs.get("expressed")
    if raw is None:
        expressed = np.full(len(values), "", dtype=object)
    else:
        cols = [m for m in present if m in raw.columns]
        M = raw[cols].to_numpy(dtype=float)
        none_measured = np.isnan(M).all(axis=1)
        expressed = np.array(cols, dtype=object)[
            np.argmax(np.where(np.isnan(M), -np.inf, M), axis=1)]
        expressed = np.where(none_measured, "", expressed)
    return {IDENTITY_ARGMAX_SCALED: top,
            IDENTITY_ARGMAX_EXPRESSED: expressed,
            IDENTITY_ARGMAX_NO_SUPPORTING: no_supporting,
            IDENTITY_ARGMAX_BY_RANK: by_rank}


def identity_argmax_audit(columns: dict[str, np.ndarray], values: pd.DataFrame) -> dict:
    """What each reading says, and how far it is from the one correction C reads.

    Counts, per donor, never a pooled rate: the three donors measured on 2026-09-18 disagree
    with each other about how big the disagreement is, which is itself worth seeing.

    `differs` counts every cell whose entry is not the current column's. `differs_both_available`
    excludes the cells where the alternative is unavailable, so the two numbers separate "this
    reading picks a different marker" from "this reading declines to pick one".

    `owns_classifier_marker` is correction C's exposure and nothing more: how many cells would
    satisfy its top-marker requirement if they were inside the matching object and cleared its
    level. It is not a count of calls and must not be read as one.

    `markers_absent` IS THE THIRD DEFECT AND IT IS RECORDED RATHER THAN REPAIRED -- see the note
    on `IDENTITY_MARKERS`. A marker that never reaches `values` can never win the argmax, and
    naming it here is what makes that checkable from the audit instead of from the source.
    """
    current = columns[CORRECTION_C_COLUMN]
    n = len(current)
    classifier = sorted({m for markers in CLASSIFIER_MARKER.values() for m in markers})
    present = [m for m in IDENTITY_MARKERS if m in values.columns]
    out = {"n_cells": n, "correction_c_reads": CORRECTION_C_COLUMN,
           "markers_declared": len(IDENTITY_MARKERS), "markers_available": present,
           "markers_absent": [m for m in IDENTITY_MARKERS if m not in values.columns],
           "columns": {}}
    for name, col in columns.items():
        unavailable = col == ""
        differs = col != current
        out["columns"][name] = {
            "unavailable": int(unavailable.sum()),
            "differs_from_correction_c_column": int(differs.sum()),
            "differs_both_available": int((differs & ~unavailable).sum()),
            "owns_classifier_marker": int(np.isin(col, classifier).sum()),
            "top_share": {str(m): int(c) for m, c in
                          pd.Series(col).value_counts().head(8).items()},
        }
    return out


def report_identity_argmax(audit: dict) -> None:
    """Print the three readings and their disagreement. Counts first, then the share."""
    A = audit.get("identity_argmax")
    if not A:
        return
    n = A["n_cells"]
    reads = A["correction_c_reads"]
    note = ("the marker the cell EXPRESSES most -- the rule" if reads == IDENTITY_ARGMAX_EXPRESSED
            else "the argmax of a PER-MARKER-NORMALISED value, NOT the top expressed marker, "
                 "which is NOT what the rule asks for")
    print(f"   identity argmax ({n:,} cells): correction C reads {reads!r} -- {note}")
    if A["markers_absent"]:
        print(f"     over {len(A['markers_available'])} of {A['markers_declared']} identity "
              f"markers; these NEVER reach the argmax and so can never win it: "
              f"{', '.join(A['markers_absent'])}")
    print(f"     {'column':<36}{'unavailable':>12}{'differs':>11}{'both avail':>12}"
          f"{'owns clf':>10}   leaders")
    for name, c in A["columns"].items():
        lead = ", ".join(f"{m} {100.0 * v / n:.1f}%" for m, v in list(c["top_share"].items())[:4])
        print(f"     {name:<36}{c['unavailable']:>12,}"
              f"{c['differs_from_correction_c_column']:>11,}"
              f"{c['differs_both_available']:>12,}{c['owns_classifier_marker']:>10,}   {lead}")


def nucleus_shape(donor: str, object_ids: pd.Series) -> pd.DataFrame:
    """Per-cell morphology for one donor, joined to `object_ids` BY object_id.

    By object_id and nothing else -- never row order, never a centroid match. The fundamentals
    table carries every detected object, the context frame only the analysis-eligible ones, so
    the two differ in length by construction and a positional join would be silently wrong.

    A cell with no nucleus measurement comes back NaN and is thereby out of every calculation
    that needs one: it cannot enter a reference fit (only finite values are given to the
    percentile) and it cannot clear a shape level, because `NaN >= level` is False. That is an
    unavailable measurement read as unavailable, and it is not zero.
    """
    columns = [SHAPE_FEATURE] + [c for c in ID.NUCLEUS_SHAPE + ID.CELL_SHAPE
                                 if c != SHAPE_FEATURE]
    path = FUNDAMENTALS / f"donor_id={donor}/data.parquet"
    if not path.exists():
        # Unavailable, and named as such: the shape criterion cannot be asked of this donor.
        return pd.DataFrame({c: np.full(len(object_ids), np.nan) for c in columns})
    frame = ID.derive_shape(pq.read_table(path).to_pandas())
    keep = ["object_id"] + [c for c in columns if c in frame.columns]
    joined = object_ids.to_frame("object_id").merge(frame[keep], on="object_id", how="left",
                                                    validate="m:1")
    for c in columns:
        if c not in joined.columns:
            joined[c] = np.nan
    return joined[columns].reset_index(drop=True)


def hormone_evidence(values: pd.DataFrame) -> dict:
    """What the three hormones say about each cell, and how confidently.

    `complete` is the point of interest. A cell one of whose hormones was not measured cannot be
    CONFIRMED endocrine -- the runner-up it would have to beat is unknown, and treating an
    unmeasured channel as a zero runner-up manufactures confidence out of an absence. Such a
    cell is simply not admitted to the reference the levels are fitted from.
    """
    present = [h for h in S.HORMONES if h in values.columns and values[h].notna().any()]
    if not present:
        raise ValueError("no hormone channel on this donor -- refusing to classify its islets")
    m = values[present].to_numpy(dtype=float)
    filled = np.nan_to_num(m, nan=-1.0)
    order = np.sort(filled, axis=1)
    top = order[:, -1]
    runner_up = order[:, -2] if m.shape[1] >= 2 else np.full(len(m), 0.0)
    which = np.array(present, dtype=object)[np.argmax(filled, axis=1)]
    complete = ~np.isnan(m).any(axis=1)
    top = np.where(top < 0, np.nan, top)
    runner_up = np.where(runner_up < 0, 0.0, runner_up)
    # islet_identity's own definition of a confirmed endocrine cell, imported rather than
    # restated: top hormone at or above REFERENCE_MIN and at least REFERENCE_MARGIN times the
    # runner-up. It is deliberately strict because this set is what "endocrine" LOOKS like, and
    # every level below is fitted from it.
    confirmed = (complete & (top >= ID.REFERENCE_MIN)
                 & (top >= ID.REFERENCE_MARGIN * np.maximum(runner_up, 1e-6)))
    # `measured` is the floor under the whole islet chain: the argmax that would name this cell
    # Beta, Alpha or Delta reads these three channels and nothing else, so with none of them
    # measured the name would come from the tie-break order rather than from the tissue.
    return {"top": top, "runner_up": runner_up, "which": which, "measured": ~np.isnan(m).all(1),
            "complete": complete, "confirmed": confirmed, "present": present}


def fit_islet_levels(donor_level: float | None, confirmed: np.ndarray, stratum: np.ndarray,
                     cd99: np.ndarray, shape: np.ndarray,
                     quantile: float = S.ENDOCRINE_QUANTILE) -> dict:
    """The endocrine levels for ONE donor, one per size stratum. Maintainer's ruling: PER DONOR
    AND PER SIZE STRATUM, restored 2026-09-18.

    Both levels are fitted from that donor's hormone-CONFIRMED cells IN THAT STRATUM -- a large
    islet's endocrine cells are the only fair yardstick for a large islet's uncertain ones, which
    is also what `islet_identity_mixture.py` does for its own calibration.

    THE FLOOR GUARD IS NOT THE STRATIFICATION AND MUST NOT BE REMOVED WITH IT. `apply_scale`
    clips at the measured background, so a marker at or below background reads exactly 0.0; on
    6534, 20.3 % of islet cells sit there, which puts a low quantile of its confirmed cells at
    0.0000 in every stratum. `CD99 >= 0.0` is true of every measured cell, so that "level" admits
    the whole islet while looking fitted. `S.endocrine_level` refuses a level at or below the
    feature's own floor and returns it as unavailable instead -- 0.0 for a scaled marker, 1.0 for
    an elongation that is a ratio of diameters.

    `donor_level` is that donor's CD99 REFERENCE LEVEL -- the p99 of its ordinary tissue, the
    estimator every other level in this classification uses. It is carried, reported beside the
    per-stratum levels, and NOT substituted for them: it is there so a reader can see both
    numbers where a stratum could not be fitted, and so the fallback is a maintainer's choice
    rather than a silent one. Measured 2026-09-18: 115 0.0518, 6442 0.0496, 6534 0.0610.
    """
    out = {"donor_reference_level": donor_level}
    for name in S.SIZE_STRATA:
        pool = confirmed & (stratum == name)
        n = int(pool.sum())
        entry = {"n_confirmed": n, "CD99": None, SHAPE_FEATURE: None}
        if n:
            entry["CD99"] = S.endocrine_level(cd99[pool], quantile)
            if name == S.LARGE_STRATUM:
                entry[SHAPE_FEATURE] = S.endocrine_level(shape[pool], quantile,
                                                         floor=SHAPE_FLOOR)
        out[name] = entry
    return out


def endocrine_decision(open_islet: np.ndarray, stratum: np.ndarray, cd99: np.ndarray,
                       elongation: np.ndarray, islet_levels: dict,
                       measured: np.ndarray) -> dict:
    """Steps 3 and 4: is this islet cell endocrine, and on which evidence.

    CD99 decides, at its own donor x stratum level. Elongation is a SECOND criterion and can only
    ever ADD: a cell the marker does not reach can still be endocrine if its nucleus is as
    elongated as that stratum's confirmed endocrine cells. It is structurally barred outside the
    large stratum -- masked here AND with no level fitted there at all.

    The mixture test's verdict on elongation (3 of 16 donors firing, median p_sep 0.119) is
    reported beside every count rather than used to remove the route: it tests ONE threshold
    against a cohort-pooled null, and a single stratification is not expected to hold across
    donors whose islets differ in shape by disease.

    `measured` -- at least one hormone channel read on this cell -- gates BOTH criteria, and it
    gates elongation for the same reason it gates CD99 even though elongation never touches a
    marker: whatever calls the cell endocrine, the argmax that then NAMES it Beta, Alpha or
    Delta has only those three channels to read. A shape is not a name. Without this gate,
    np.argmax over three NaNs filled to -1 returns index 0 -- "Beta" -- every time, which on
    6534 was 23,554 large-islet cells named from no hormone measurement at all.
    """
    n = len(open_islet)
    open_islet = np.asarray(open_islet) & np.asarray(measured)
    level = np.array([islet_levels.get(s, {}).get("CD99") for s in stratum], dtype=object)
    have_level = np.array([v is not None for v in level]) & open_islet
    cd99_ok = np.zeros(n, dtype=bool)
    if have_level.any():
        cd99_ok[have_level] = (np.asarray(cd99, dtype=float)[have_level]
                               >= np.array(list(level[have_level]), dtype=float))

    large = np.asarray(stratum) == S.LARGE_STRATUM
    shape_level = islet_levels.get(S.LARGE_STRATUM, {}).get(SHAPE_FEATURE)
    shape_ok = (np.zeros(n, dtype=bool) if shape_level is None
                else large & have_level & (np.asarray(elongation, dtype=float) >= shape_level))
    # ELONGATION IS A ROUTE: a cell the key marker does not reach can still be endocrine if its
    # nucleus is as elongated as that stratum's confirmed endocrine cells. Maintainer's ruling,
    # restored 2026-09-18 after I removed it without being asked to.
    #
    # WHAT THE MIXTURE TEST SAID, recorded because the number matters when reading any count
    # this produces: islet_identity_mixture.py over 20 donors puts nucleus elongation at 3 of 16
    # fitted donors firing, median separation p 0.119 against a Bonferroni bar of 0.0083, against
    # CD99's 5 of 10 at p_LR 0.0010 / p_sep 0.0100. That test pools each stratum's donors against
    # one null, and the maintainer's standing point is that a single stratification is not
    # expected to hold across donors -- T1D islets do not have the shape ND islets do -- so a
    # cohort-level test of one threshold is not the last word on it. `n_by_elongation_only` is
    # reported separately on every donor so the contribution is always visible.
    return {"have_level": have_level, "cd99_ok": cd99_ok, "shape_ok": shape_ok,
            "endocrine": have_level & (cd99_ok | shape_ok), "large": large}


def _structure_calls(ctx: pd.DataFrame, values: pd.DataFrame, levels: dict, log: dict,
                     top_marker: np.ndarray, require_top_marker: bool,
                     area: dict[str, np.ndarray]) -> tuple:
    """Step 1. Structures, weakest first so a stronger class overwrites.

    Asked of EVERY cell, islet included: the object evidence is real inside an islet too.

    THE MAINTAINER'S RULE, in four steps, none of them sufficient alone:

      (a) the cell's CENTROID is exactly inside the object's bounds;
      (b) THE WHOLE CELL IS THEN CHECKED -- the MAJORITY of the cell's own area must lie within
          the object (`MAJORITY_AREA`), and that is what makes it a CANDIDATE;
      (c) the marker that CREATED the object is simply higher than all the rest;
      (d) otherwise PASSTHROUGH, carrying the FRACTION of the cell's area the object covers.

    Step (b) is the rule. A centroid inside an object is necessary and nowhere near sufficient --
    an object is a region of tissue and most of what sits in it is whatever the region runs
    through. Measured 2026-09-19 on donor 115: 125,725 cells have a centroid inside a vessel
    object and only 14,661 are majority-inside. The step was ABSENT from this function until
    2026-09-19, so centroid-inside was being treated as candidacy; it is not.

    A cell whose brightest marker is Ker8_18 is not a blood vessel because it lies inside one; it
    is an exocrine cell the vessel passes through, and it keeps its own identity and carries how
    much it was crossed as a measurement.

    A SMALL NUMBER OF STRUCTURE CALLS IS THE EXPECTED RESULT AND IS NOT A DEFECT. Most vessels,
    lymphatics and nerves are thin processes running through and around cells rather than the cell
    bodies themselves, so most of what an object contains is other cells it passes. Those are
    passthroughs by design. When a CELL BODY is there we want it called, and step (b) is what
    finds it. Do not loosen `MAJORITY_AREA` because a class came out small -- the survivors being
    few is the rule working, and the pixel classifier has already done the work of saying where
    the object is.

    THE AREA MAY BE UNAVAILABLE, and that is a third outcome, not a refusal. A cell correction D
    never measured cannot be shown to be majority-inside NOR shown not to be, so it is neither
    called nor quietly passed through: it is counted under `inside_area_unavailable` and its
    passthrough fraction stays NaN. Do not read that as zero overlap.

    `top_marker` IS NOT THE CELL'S TOP EXPRESSED MARKER, whatever the column is called. It is
    the argmax of a per-marker-normalised value, so a broadly and moderately expressed marker
    takes it from a sharply expressed one -- on the three donors measured 2026-09-18 a
    `lineages.SUPPORTING_ONLY` marker wins it on 30-33 % of cells, FOXP3 alone on 21-25 %. That
    is measured in the note on `IDENTITY_MARKERS`, two alternative readings are computed and
    written beside it by `identity_argmax`, and `CORRECTION_C_COLUMN` is the one line that says
    which one this test reads. Which it SHOULD read is the maintainer's ruling and is open.

    The top-marker condition is correction C and it is a deliberate TIGHTENING, not part of the
    reconstruction. `--no-top-marker` drops it so its effect can be isolated in a diff, and the
    audit counts what it moves per structure.
    """
    n = len(ctx)
    cell_class = np.full(n, None, dtype=object)
    resolved_by = np.full(n, "", dtype=object)
    passed_through: dict[str, np.ndarray] = {}
    audit: dict[str, dict] = {}
    for name in reversed(STRUCTURES):
        markers = [m for m in log["subpanels"][name] if m in values.columns]
        if not markers:
            continue
        inside = ctx[f"in_{name}"].to_numpy()
        clears = np.nanmax(values[markers].to_numpy(), axis=1) >= levels[name]
        owns = np.isin(top_marker, CLASSIFIER_MARKER[name])
        # Step (b). NaN is unavailable, and `NaN > x` is False, so an unmeasured cell is never a
        # candidate -- but it is also never counted as "measured and not majority", which is why
        # the two are tallied apart below.
        frac = area.get(name, np.full(n, np.nan))
        known = ~np.isnan(frac)
        majority = known & (frac > MAJORITY_AREA)
        candidate = inside & majority
        hit = candidate & clears & (owns if require_top_marker else True)
        cell_class[hit] = name
        resolved_by[hit] = "structure"
        # Step (d). Everything else the object touches keeps its own identity and carries HOW MUCH
        # of it the object covers -- the fraction, not a flag. NaN where it was never measured.
        passed_through[name] = np.where(inside & ~hit, frac, np.nan)
        audit[name] = {
            "inside": int(inside.sum()),
            # Step (b), reported on its own because it is the step that was missing.
            "inside_and_majority": int(candidate.sum()),
            "inside_area_unavailable": int((inside & ~known).sum()),
            "majority_withholds": int((inside & known & ~majority).sum()),
            "inside_and_clears": int((candidate & clears).sum()),
            "called": int(hit.sum()),
            # What correction C moves: inside the object, clearing its level, but the cell's own
            # brightest identity marker is not one that drew the object.
            "top_marker_withholds": int((candidate & clears & ~owns).sum()),
            "top_marker_withholds_in_islet": int(
                (candidate & clears & ~owns & ctx["in_islet"].to_numpy()).sum()),
        }
    return cell_class, resolved_by, passed_through, audit


def _tissue_chain(values: pd.DataFrame, log: dict, rules: dict, cell_class: np.ndarray,
                  definition: np.ndarray, resolved_by: np.ndarray) -> None:
    """Steps 3-5 of the ordinary chain, UNCHANGED: immune definitions, supporting tier, Exocrine.

    Mutates the three arrays in place. Every cell the islet chain does not call endocrine comes
    through here, which is the point: the non-endocrine cells inside an islet are the analysis,
    so they are typed by the same definitions as everything else rather than discarded.
    """
    n = len(cell_class)
    unclaimed = cell_class == None  # noqa: E711
    for phenotype, rule in sorted(rules.items()):
        hit = unclaimed & fires(values, rule)
        cell_class[hit] = "Immune"
        definition[hit] = phenotype
        resolved_by[hit] = "definition"
        unclaimed &= ~hit

    # The supporting tier. A cell is eligible when its immune evidence CLEARS the level and its
    # exocrine evidence does NOT -- the same level on both sides, so the tier never overrules a
    # cell the exocrine markers speak for as loudly as the immune ones do.
    #
    # This was recovered, not chosen, and the recovery is worth recording because it is the one
    # part of the method RECONSTRUCTION.md left open. All 35,332 cells the original table claims
    # by this tier on 6436 carry exocrine evidence below 0.8715; none is above it. Applying the
    # tier to every unclaimed cell instead -- the permissive form -- wrongly claims 5,671 strongly
    # exocrine cells there (exocrine max p50 0.98). Restricting to cells with NO exocrine evidence
    # was already measured too narrow: it keeps only 19,142 of the 35,332.
    #
    # It is not exact. 1,223 cells on 6436 remain wrongly claimed, 98.4 % of them driven by PD_1
    # on cells whose exocrine evidence sits just under the level (p50 0.843 against 0.044 for the
    # cells the tier should claim). Closing that would mean picking a threshold because it fits
    # the answer, which is not recovery, so it is left open and reported.
    tier = [m for m in log["bottom_tier_immune"] if m in values.columns]
    exocrine = [m for m in S.SUBPANEL["exocrine"] if m in values.columns]
    if tier:
        clears = np.nanmax(values[tier].to_numpy(), axis=1) >= SUPPORTING_LEVEL
        speaks_exocrine = (np.nanmax(values[exocrine].to_numpy(), axis=1) >= SUPPORTING_LEVEL
                           if exocrine else np.zeros(n, dtype=bool))
        hit = unclaimed & clears & ~speaks_exocrine
        cell_class[hit] = "Immune"
        resolved_by[hit] = "supporting"
        unclaimed &= ~hit

    cell_class[unclaimed] = "Exocrine"
    resolved_by[unclaimed] = "fallthrough"


def _finish(in_islet: np.ndarray, cell_class: np.ndarray,
            exocrine_measured: np.ndarray | None = None) -> np.ndarray:
    """Step 2, applied where it can be: exocrine inside an islet outline is an annotation error.

    Exocrine markers get into an islet annotation only by mistake -- that is the maintainer's
    rule and it is a statement about the annotation, not about the cell. So a cell the chain can
    only call Exocrine while its centroid sits inside an islet outline is ruled out of analysis
    and carries the reason, rather than being counted as an exocrine cell in an islet.

    THE REASON HAS TO BE TRUE OF THE CELL IT IS WRITTEN ON. `exocrine_inside_islet_border`
    asserts that exocrine markers reached this cell, and that assertion needs those markers to
    have been MEASURED. A cell whose Ker8_18, E-cadherin and Pan-Cytokeratin are all missing
    reaches Exocrine by the same fall-through, but for the opposite cause: not exocrine
    evidence, no evidence. It is ruled out under `no_marker_values` instead.

    This is not hypothetical. Donor 6534's aligned table carries 205,070 rows with no value for
    any of the 57 lineage markers, and 25,461 of its islet cells reach Exocrine with neither a
    hormone nor an exocrine marker read on them. Counted under the exocrine reason they would
    inflate, by more than a factor of two, a number that is read as "how often the islet
    annotation caught exocrine tissue".

    Only the REASON moves here, never the class: a cell the ordinary chain called Immune on
    markers that WERE measured keeps that call, because the non-endocrine cells inside the islet
    are the analysis.
    """
    reason = np.full(len(cell_class), "", dtype=object)
    ruled_out = in_islet & (cell_class == "Exocrine")
    cell_class[ruled_out] = "Excluded"
    reason[ruled_out] = ISLET_EXOCRINE_REASON
    if exocrine_measured is not None:
        reason[ruled_out & ~np.asarray(exocrine_measured)] = NO_MEASUREMENT_REASON
    return reason


def _run_chain(values, log, rules, structure_class, structure_resolved, in_islet,
               endocrine: np.ndarray, endocrine_name: np.ndarray, endocrine_route: np.ndarray,
               exocrine_measured: np.ndarray | None = None):
    """One full pass: structures (given), the islet answer (given), then the ordinary chain."""
    cell_class = structure_class.copy()
    resolved_by = structure_resolved.copy()
    definition = np.full(len(cell_class), "", dtype=object)
    hit = endocrine & (cell_class == None)  # noqa: E711
    cell_class[hit] = endocrine_name[hit]
    resolved_by[hit] = endocrine_route[hit]
    _tissue_chain(values, log, rules, cell_class, definition, resolved_by)
    reason = _finish(in_islet, cell_class, exocrine_measured)
    return cell_class, definition, resolved_by, reason


def classify(donor: str, ctx: pd.DataFrame, levels: dict, scale: pd.DataFrame,
             log: dict, detail: bool, require_top_marker: bool = True,
             islet_chain: bool = True, quantile: float = S.ENDOCRINE_QUANTILE,
             audit: dict | None = None) -> pd.DataFrame:
    """The subpanel classification for one donor.

    THE ISLET CHAIN, for a cell whose centroid is inside an islet outline:

      1. STRUCTURE -- inside a classified object, clearing that object's level, and the object's
         own marker is the cell's entry in `CORRECTION_C_COLUMN` (see `identity_argmax`: that
         column is an argmax of per-marker-normalised values, NOT the top expressed marker).
         The classified objects separate vessels and nerves, and they do so inside islets too.
      2. EXOCRINE INSIDE AN ISLET IS AN ANNOTATION ERROR, applied at the end of the chain where
         the evidence is in: `Excluded`, with `exocrine_inside_islet_border` as the reason.
      3. ENDOCRINE OR NOT, ON CD99 -- the key differentiator -- against the level fitted from
         THIS donor's AND THIS size stratum's own hormone-confirmed cells, with the floor guard.
         The donor's ordinary-tissue CD99 level is carried and reported beside it, never
         substituted for it.
      4. ELONGATION as a second ROUTE that can only ADD, in the large stratum only.
      5. ENDOCRINE -> the hormone argmax names Beta / Alpha / Delta.
      6. NOT ENDOCRINE -> the ordinary tissue chain, unchanged.

    ISLET SIZE, which decides the stratum, is counted over the cells the analysis ADMITS: islet
    core, non-structural, not `Excluded`. That is `islet_identity.load`'s definition and it is
    reused rather than restated -- but it is a property of a COMPLETED classification, so it
    cannot be known before the chain runs. The chain therefore runs twice: once as it stood
    before this correction, purely to fix the admitted set and so the sizes, and once with the
    islet chain. The first pass is also exactly the baseline the audit reports against, so the
    "what moved" numbers cost nothing extra.

    `--no-islet-chain` returns the first pass, which is the behaviour before this change.
    """
    rules = immune_rules()
    wanted = set(S.HORMONES) | {ID.MARKER}
    for markers in log["subpanels"].values():
        wanted |= set(markers)
    for all_of, any_of, none_of in rules.values():
        wanted |= set(all_of) | set(none_of) | {m for g in any_of for m in g}
    wanted |= set(log["bottom_tier_immune"])
    values = scaled_markers(donor, wanted, scale, ctx["object_id"])

    if "islet_num" not in ctx.columns:
        raise ValueError(
            f"{donor}: the cached context predates islet sizing and carries no islet_num. "
            f"Rerun with --refit so the stratum can be counted.")

    n = len(ctx)
    in_islet = ctx["in_islet"].to_numpy()
    name_of = {v: k for k, v in log["hormones"].items()}
    A = {} if audit is None else audit
    A["donor_id"] = donor
    A["n_cells"] = n
    A["n_in_islet"] = int(in_islet.sum())

    # Three readings of the identity argmax, computed side by side and all three written. Which
    # one correction C reads is `CORRECTION_C_COLUMN`, one line at module level, and it is left
    # where it was: choosing between them moves calls and is the maintainer's to choose.
    identity_columns = identity_argmax(values)
    top_marker = identity_columns[CORRECTION_C_COLUMN]
    A["identity_argmax"] = identity_argmax_audit(identity_columns, values)

    # --- 1. structures ------------------------------------------------------------------------
    area = crossed_area(donor, ctx["object_id"])
    A["area_measured"] = {name: int((~np.isnan(v)).sum()) for name, v in area.items()}
    structure_class, structure_resolved, passed_through, A["structures"] = _structure_calls(
        ctx, values, levels, log, top_marker, require_top_marker, area)
    structural = structure_class != None  # noqa: E711

    # --- the hormone evidence, and the baseline pass -------------------------------------------
    horm = hormone_evidence(values)
    A["hormones_present"] = horm["present"]
    A["n_incomplete_hormone"] = int((~horm["complete"]).sum())
    # Two different absences, kept apart because they gate different things. `measured` is "a
    # hormone was read on this cell" and it gates the endocrine NAME -- the argmax that writes
    # Beta, Alpha or Delta reads nothing else. `exo_measured` is "an exocrine marker was read"
    # and it gates the exclusion REASON, because `exocrine_inside_islet_border` is a claim about
    # exocrine markers and cannot be written on a cell that has none.
    exo = [m for m in S.SUBPANEL["exocrine"] if m in values.columns]
    exo_measured = (values[exo].notna().any(axis=1).to_numpy() if exo
                    else np.zeros(n, dtype=bool))
    A["n_no_hormone"] = int((~horm["measured"]).sum())
    A["n_no_hormone_in_islet"] = int((in_islet & ~horm["measured"]).sum())
    A["n_no_exocrine_marker_in_islet"] = int((in_islet & ~exo_measured).sum())
    base_endocrine = in_islet & horm["measured"] & (horm["top"] >= levels["hormone"])
    endocrine_name = np.array([name_of[w] for w in horm["which"]], dtype=object)
    base_class, base_definition, base_resolved, base_reason = _run_chain(
        values, log, rules, structure_class, structure_resolved, in_islet,
        base_endocrine, endocrine_name, np.full(n, "hormone", dtype=object), exo_measured)

    # --- islet size, over the cells the analysis admits ----------------------------------------
    admitted = in_islet & ~structural & (base_class != "Excluded")
    islet_num = ctx["islet_num"].to_numpy()
    sizes = pd.Series(islet_num[admitted]).value_counts()
    islet_cells = pd.Series(islet_num).map(sizes).to_numpy(dtype=float)
    islet_cells = np.where(in_islet, np.nan_to_num(islet_cells, nan=0.0), np.nan)
    stratum = S.size_stratum(islet_cells)
    A["n_admitted"] = int(admitted.sum())
    A["stratum_cells"] = {s: int((stratum == s).sum()) for s in S.SIZE_STRATA}

    if not islet_chain:
        A["islet_chain"] = False
        frame = _emit(donor, ctx, base_class, base_definition, base_resolved, base_reason,
                      islet_cells, stratum, np.full(n, "", dtype=object),
                      identity_columns, passed_through, detail)
        return frame

    # --- 3. CD99 at the donor x stratum level, and 4. elongation in the large stratum ----------
    cd99 = (values[ID.MARKER].to_numpy(dtype=float) if ID.MARKER in values.columns
            else np.full(n, np.nan))
    shape = nucleus_shape(donor, ctx["object_id"])
    elongation = shape[SHAPE_FEATURE].to_numpy(dtype=float)
    A["n_no_nucleus"] = int((in_islet & ~np.isfinite(elongation)).sum())

    confirmed_pool = admitted & horm["confirmed"]
    A["n_confirmed_outside_admitted"] = int((in_islet & horm["confirmed"] & ~admitted).sum())
    islet_levels = fit_islet_levels(levels.get(ID.MARKER), confirmed_pool, stratum,
                                    cd99, elongation, quantile)
    A["islet_levels"] = islet_levels
    A["endocrine_quantile"] = quantile

    # DIAGNOSTIC ONLY, never consulted by the chain: the level each stratum WOULD have produced
    # with no minimum-count guard, and the confirmed/undecided medians behind it. This is what
    # answers "did stratifying change anything" -- if the four levels barely differ, the ruling
    # cost nothing, and that is worth knowing from the same run rather than from an argument.
    undecided = admitted & (horm["top"] < ID.UNDECIDED_MAX)
    unguarded: dict = {}
    for name in S.SIZE_STRATA:
        pool = confirmed_pool & (stratum == name)
        und = undecided & (stratum == name)
        unguarded[name] = {
            "n_confirmed": int(pool.sum()), "n_undecided": int(und.sum()),
            # floor=-inf on purpose: the diagnostic must SHOW a degenerate level rather than
            # refuse it, because "this stratum's level lands on the scale floor" is the finding.
            "CD99": (S.endocrine_level(cd99[pool], quantile, 1, float("-inf"))
                     if pool.any() else None),
            SHAPE_FEATURE: (S.endocrine_level(elongation[pool], quantile, 1, float("-inf"))
                            if pool.any() else None),
            "cd99_median_confirmed": float(np.nanmedian(cd99[pool])) if pool.any() else None,
            "cd99_median_undecided": float(np.nanmedian(cd99[und])) if und.any() else None,
            "elong_var_ratio": (float(np.nanvar(elongation[und]) / np.nanvar(elongation[pool]))
                                if und.any() and pool.any()
                                and np.nanvar(elongation[pool]) > 0 else None),
            "cd99_var_ratio": (float(np.nanvar(cd99[und]) / np.nanvar(cd99[pool]))
                               if und.any() and pool.any()
                               and np.nanvar(cd99[pool]) > 0 else None)}
    A["islet_levels_unguarded"] = unguarded

    open_islet = in_islet & ~structural
    D = endocrine_decision(open_islet, stratum, cd99, elongation, islet_levels,
                           horm["measured"])
    have_level, cd99_ok, shape_ok, large = (D["have_level"], D["cd99_ok"], D["shape_ok"],
                                            D["large"])

    # Where the stratum has no fittable level the question CANNOT be asked -- the level is
    # unavailable, not zero and not the cohort's. Those cells keep the hormone question they
    # were asked before, and the audit counts them so an unanswerable stratum stays visible.
    fallback = open_islet & ~have_level
    endocrine = D["endocrine"] | (fallback & base_endocrine)
    route = np.full(n, "", dtype=object)
    route[endocrine & cd99_ok] = "cd99"
    route[endocrine & ~cd99_ok & shape_ok] = "elongation"
    route[endocrine & fallback] = "hormone"

    cell_class, definition, resolved_by, reason = _run_chain(
        values, log, rules, structure_class, structure_resolved, in_islet,
        endocrine, endocrine_name, route, exo_measured)

    # --- what moved ---------------------------------------------------------------------------
    A["islet_chain"] = True
    A["n_endocrine_base"] = int((base_endocrine & open_islet).sum())
    A["n_endocrine_new"] = int(endocrine.sum())
    A["n_by_cd99"] = int((endocrine & cd99_ok).sum())
    A["n_by_elongation_only"] = int((endocrine & ~cd99_ok & shape_ok).sum())
    A["n_by_hormone_fallback"] = int((endocrine & fallback).sum())
    A["n_hormone_clearing_demoted"] = int(
        (open_islet & have_level & base_endocrine & ~endocrine).sum())
    A["n_confirmed_demoted"] = int(
        (open_islet & have_level & horm["confirmed"] & ~endocrine).sum())
    # The reading this chain did NOT take: elongation as a requirement rather than a second
    # route. Counted, not applied, so the maintainer can rule on it against a number.
    A["n_large_cd99_ok_but_round"] = int((large & have_level & cd99_ok & ~shape_ok).sum())
    changed = base_class != cell_class
    A["n_changed"] = int(changed.sum())
    # Named on its own, because the study question is the immune environment around and inside
    # islets: a cell the immune DEFINITIONS claimed before and the CD99/elongation criteria claim
    # now is an immune cell removed from that analysis, and it must never arrive as a footnote
    # inside a confusion table.
    endocrine_names = set(log["hormones"])
    A["n_immune_to_endocrine"] = int(
        ((base_class == "Immune") & np.isin(cell_class, list(endocrine_names))).sum())
    A["n_excluded_to_endocrine"] = int(
        ((base_class == "Excluded") & np.isin(cell_class, list(endocrine_names))).sum())
    A["n_endocrine_to_other"] = int(
        (np.isin(base_class, list(endocrine_names)) & ~np.isin(cell_class,
                                                               list(endocrine_names))).sum())
    A["changes"] = {f"{a} -> {b}": int(c) for (a, b), c in
                    pd.Series(list(zip(base_class[changed], cell_class[changed])))
                    .value_counts().items()}
    A["per_stratum"] = {
        s: {"cells": int((stratum == s).sum()),
            "endocrine_base": int((base_endocrine & open_islet & (stratum == s)).sum()),
            "endocrine_new": int((endocrine & (stratum == s)).sum()),
            "by_cd99": int((endocrine & cd99_ok & (stratum == s)).sum()),
            "by_elongation_only": int((endocrine & ~cd99_ok & shape_ok & (stratum == s)).sum()),
            "changed": int((changed & (stratum == s)).sum())}
        for s in S.SIZE_STRATA}

    return _emit(donor, ctx, cell_class, definition, resolved_by, reason, islet_cells, stratum,
                 route, identity_columns, passed_through, detail)


def _emit(donor, ctx, cell_class, definition, resolved_by, reason, islet_cells, stratum,
          route, identity_columns, passed_through, detail) -> pd.DataFrame:
    """The table.

    `object_context` is written on every row, not only under `--emit-detail`. It is what the
    classified objects said before the islet precedence overwrote the context, so it is the only
    place a cell that sits in a vessel INSIDE an islet can be seen as one. That gap misled an
    earlier analysis of this lane; a column that exists only under a flag would leave it open.
    """
    frame = pd.DataFrame({
        "donor_id": donor, "object_id": ctx["object_id"].to_numpy(),
        "cell_class": cell_class.astype(str), "context": ctx["context"].to_numpy(),
        "object_context": ctx["object_context"].to_numpy(),
        "cell_region": ctx["cell_region"].to_numpy(), "excluded_reason": reason.astype(str),
        "islet_num": ctx["islet_num"].to_numpy(),
        "islet_cells": islet_cells, "islet_stratum": np.asarray(stratum, dtype=object),
        "islet_route": np.asarray(route, dtype=object)})
    if detail:
        frame["immune_definition"] = definition.astype(str)
        frame["resolved_by"] = resolved_by.astype(str)
        # All three readings of the identity argmax, side by side and named for what each one
        # is. `top_identity_marker` keeps its name and its values so nothing downstream breaks;
        # the other two are there so the disagreement is on the table rather than in a message.
        # See `identity_argmax` for what separates them.
        #
        # A BARE ARRAY IS STILL ACCEPTED, and it means the one column it always meant. This
        # helper's existing callers -- `core/phenocycler/tests/current/test_build_cells_carry.py`
        # among them -- pass the top-marker array positionally, and widening a private signature
        # is not a reason to make them wrong.
        columns = (identity_columns if isinstance(identity_columns, dict)
                   else {IDENTITY_ARGMAX_SCALED: identity_columns})
        for name, column in columns.items():
            frame[name] = column
        # Step (d) of the structure rule: the passthrough carries the FRACTION of the cell's own
        # area that object covers, not a flag -- "being crossed is a measurement, not a
        # disqualification". NaN where correction D never measured the cell, which is UNAVAILABLE
        # and is not a zero overlap. It was a boolean until 2026-09-19.
        for name, fraction in passed_through.items():
            frame[f"passed_through_{name.replace(' ', '_')}"] = fraction
    return frame


def report_islet_audit(audit: dict) -> None:
    """Print what the islet chain moved, per donor and per stratum, in counts.

    Every line here is a count of cells, never a rate presented on its own: a percentage of an
    unstated denominator is how a change of 400 cells and a change of 400,000 come to read the
    same. The levels are printed whether or not they moved anything, because "the strata barely
    differ" is itself the answer to whether stratifying was needed.
    """
    if not audit or not audit.get("islet_chain"):
        return
    print(f"   islet chain: {audit['n_in_islet']:,} cells in an islet, "
          f"{audit['n_admitted']:,} admitted to the size count")
    print(f"     endocrine {audit['n_endocrine_base']:,} (hormone level) -> "
          f"{audit['n_endocrine_new']:,}   by CD99 {audit['n_by_cd99']:,}, "
          f"by elongation alone {audit['n_by_elongation_only']:,}, "
          f"hormone fallback {audit['n_by_hormone_fallback']:,}")
    print(f"     hormone-clearing cells CD99 demotes {audit['n_hormone_clearing_demoted']:,}"
          f"  (of which hormone-CONFIRMED {audit['n_confirmed_demoted']:,})")
    # KNOWN DEFECT, named here rather than fixed in a documentation pass: `islet_levels` is keyed
    # by `donor_reference_level` and the four stratum names -- there is no top-level "CD99" key --
    # so this lookup returns None on every donor and the line below always prints "unavailable".
    # The per-stratum CD99 levels and the donor reference level are in the audit JSON
    # (`--audit-out`) and are simply not reached by this print. Do not read the printed
    # "unavailable" as a statement about the data.
    cd99 = audit["islet_levels"].get("CD99")
    print(f"     CD99 level (p99 of this donor's ordinary tissue, one per donor -- the reference\n"
          f"     population lies outside every islet, so it carries no islet size): "
          + ("unavailable" if cd99 is None else f"{cd99:.4f}"))
    print(f"     {'stratum':<8}{'cells':>12}{'confirmed':>11}"
          f"{'elong level':>13}{'endo base':>11}{'endo new':>10}{'by elong':>10}")
    for name in S.SIZE_STRATA:
        lv = audit["islet_levels"][name]
        st = audit["per_stratum"][name]
        el = "-" if lv[SHAPE_FEATURE] is None else f"{lv[SHAPE_FEATURE]:.4f}"
        print(f"     {name:<8}{st['cells']:>12,}{lv['n_confirmed']:>11,}{el:>13}"
              f"{st['endocrine_base']:>11,}{st['endocrine_new']:>10,}"
              f"{st['by_elongation_only']:>10,}")
    print(f"     {'stratum':<8}{'CD99 conf':>11}{'CD99 und':>10}{'ratio':>8}"
          f"{'elong var':>11}{'CD99 var':>10}{'unguarded lvl':>15}   (diagnostic, not used)")
    for name in S.SIZE_STRATA:
        u = audit["islet_levels_unguarded"][name]
        f = lambda v, w=10, d=3: (f"{v:>{w}.{d}f}" if isinstance(v, float) else f"{'-':>{w}}")
        ratio = (u["cd99_median_confirmed"] / u["cd99_median_undecided"]
                 if u["cd99_median_confirmed"] and u["cd99_median_undecided"] else None)
        print(f"     {name:<8}{f(u['cd99_median_confirmed'], 11)}{f(u['cd99_median_undecided'])}"
              f"{f(ratio, 8, 2)}{f(u['elong_var_ratio'], 11)}{f(u['cd99_var_ratio'])}"
              f"{f(u['CD99'], 15, 4)}")
    if audit["n_changed"]:
        top = sorted(audit["changes"].items(), key=lambda kv: -kv[1])[:6]
        print(f"     {audit['n_changed']:,} cells change class: " +
              "; ".join(f"{k} {v:,}" for k, v in top))
        print(f"     immune -> endocrine {audit['n_immune_to_endocrine']:,}   "
              f"excluded -> endocrine {audit['n_excluded_to_endocrine']:,}   "
              f"endocrine -> other {audit['n_endocrine_to_other']:,}")
    for name, st in audit["structures"].items():
        if st["top_marker_withholds"]:
            print(f"     correction C withholds {name}: {st['top_marker_withholds']:,} cells "
                  f"inside and clearing ({st['top_marker_withholds_in_islet']:,} of them in an "
                  f"islet); {st['called']:,} called")


def gate(rebuilt: pd.DataFrame, donor: str) -> dict:
    """Row-for-row against the table on disk, on (donor_id, object_id, cell_class)."""
    old = pq.read_table(OUT / "subpanel_calls.parquet",
                        columns=["object_id", "cell_class", "context"],
                        filters=[("donor_id", "=", donor)]).to_pandas()
    merged = old.merge(rebuilt[["object_id", "cell_class", "context"]], on="object_id",
                       how="outer", suffixes=("_old", "_new"), indicator=True)
    both = merged["_merge"] == "both"
    agree = (merged.loc[both, "cell_class_old"] == merged.loc[both, "cell_class_new"])
    ctx_agree = (merged.loc[both, "context_old"] == merged.loc[both, "context_new"])
    confusion = (merged.loc[both & ~agree.reindex(merged.index, fill_value=False)]
                 .groupby(["cell_class_old", "cell_class_new"]).size()
                 .sort_values(ascending=False).head(6))
    return {"donor_id": donor, "rows_old": int((merged["_merge"] != "right_only").sum()),
            "rows_new": int((merged["_merge"] != "left_only").sum()),
            "only_old": int((merged["_merge"] == "left_only").sum()),
            "only_new": int((merged["_merge"] == "right_only").sum()),
            "class_match_pct": 100.0 * float(agree.mean()) if len(agree) else float("nan"),
            "context_match_pct": 100.0 * float(ctx_agree.mean()) if len(ctx_agree) else float("nan"),
            "confusion": confusion}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--membership", choices=("parent", "geometry"), default="parent",
                    help="geometry is the method -- islet membership from the islet OUTLINES "
                         "(islet_membership.membership). parent is QuPath's `Parent` column, "
                         "the SUPERSEDED behaviour, and it is still the default because stage 1 "
                         "is a reconstruction gate: a bare run has to reproduce the table on "
                         "disk, which was built with parent. THE DEFAULT IS NOT THE METHOD. "
                         "Its measured cost: the Parent labels cells `peri` that sit up to "
                         "74.6 um inside a core outline, and because a false in_islet means the "
                         "cell is never asked the islet question at all, 11,556 islet-core cells "
                         "cohort-wide (2.7 %%; 35.3 %% on 6539) fell through to Exocrine, 5,250 of "
                         "them hormone-positive -- measured 2026-09-15, recorded in "
                         "islet_membership.py's docstring and context.py's `memberships`, and "
                         "restated in RECONSTRUCTION.md under Context. A real rebuild passes "
                         "--membership geometry explicitly; --gate pins parent by itself")
    ap.add_argument("--emit-detail", action="store_true")
    ap.add_argument("--no-top-marker", action="store_true",
                    help="drop correction C's requirement, to isolate its effect in a diff")
    ap.add_argument("--no-islet-chain", action="store_true",
                    help="ask the islet only the pooled hormone question, as before the CD99 "
                         "and elongation criteria -- to isolate their effect in a diff")
    ap.add_argument("--endocrine-quantile", type=float, default=S.ENDOCRINE_QUANTILE,
                    help="quantile of the hormone-confirmed cells that sets the CD99 and "
                         "elongation levels (default: the mirror of the reference quantile)")
    ap.add_argument("--work", default=None,
                    help="directory for the cached per-donor context (default: the typing tree)")
    ap.add_argument("--audit-out", default=None,
                    help="write the per-donor islet-chain audit as JSON here")
    ap.add_argument("--levels", choices=("per-donor", "pooled"), default="per-donor",
                    help="per-donor is the method; pooled exists only to reproduce the old table")
    ap.add_argument("--refit", action="store_true", help="recompute context and levels")
    ap.add_argument("--gate", action="store_true",
                    help="compare against the table on disk instead of writing")
    ap.add_argument("--out", default=None, help="output directory (default: subpanel_rebuilt)")
    args = ap.parse_args()

    log = _log()
    scale = S.load_scale()
    pixel, stem = MEMBERSHIP._manifest()
    # Every donor that CAN be classified is one with aligned marker values. This read the
    # historical run tree's expression partitions until 2026-10-01 -- the same 20 donors while that
    # tree was the only source, and none at all in a checkout whose cells come from revised exports.
    every = sorted(p.name.split("=")[1] for p in S.ALIGNED.glob("donor_id=*"))
    donors = args.donors or every
    work = pathlib.Path(args.work) if args.work else WORK
    work.mkdir(parents=True, exist_ok=True)

    print(f"membership: {args.membership}   classifying {len(donors)} of {len(every)} donors")
    t0 = time.time()
    # Per donor, a donor needs only its own context. Pooled levels cannot be fitted without the
    # whole cohort, so that mode -- and only that mode -- pays for every donor.
    scope = every if args.levels == "pooled" else donors
    # Keyed by membership as well as donor: the two modes give different context, and a cache
    # that ignored that would silently serve one mode's answer to the other.
    key = lambda d: work / f"{d}__{args.membership}.parquet"
    stale = [d for d in scope if not key(d).exists() or args.refit]
    for donor in stale:
        ctx = donor_context(donor, args.membership, pixel, stem)
        ctx.to_parquet(key(donor), index=False)
        print(f"  context {donor}: {len(ctx):,} eligible cells, "
              f"{int(ctx['in_islet'].sum()):,} in an islet  [{S.cell_source(donor)[1]}]"
              f"   [{time.time()-t0:.0f}s]", flush=True)

    level_cache = work / LEVEL_CACHE.name
    if args.levels == "pooled":
        if level_cache.exists() and not args.refit and not stale:
            levels = json.loads(level_cache.read_text())
            print("   (pooled levels read from cache; --refit to recompute)")
        else:
            levels = fit_levels(every, scale, log, args.membership, work)
            level_cache.write_text(json.dumps(levels, indent=1))
        print("\nPOOLED reference levels -- reproduces the old table; NOT the method:")
        for name, value in levels.items():
            recorded = log["bars"].get(name) or (0.1010 if name == "hormone" else None)
            print(f"   {name:14s} {value:.4f}" +
                  (f"   recorded {recorded:.4f}" if recorded else ""))
        per_donor = {d: levels for d in donors}
    else:
        per_donor = {d: fit_levels([d], scale, log, args.membership, work) for d in donors}
        print("\nPER-DONOR reference levels (p99 of that donor's OWN ordinary tissue):")
        frame = pd.DataFrame(per_donor).T
        frame.index.name = "donor"
        print(frame.to_string(float_format=lambda v: f"{v:.4f}"))

    target = pathlib.Path(args.out) if args.out else S.OUT / "subpanel_rebuilt"
    if not args.gate:
        target.mkdir(parents=True, exist_ok=True)

    if args.gate and args.levels != "pooled":
        print("\nNOTE: --gate against per-donor levels is NOT a reconstruction gate. The table on\n"
              "disk was built with pooled levels, so what follows is the CORRECTION'S EFFECT, not\n"
              "a reproduction failure. Run `--gate --levels pooled` for the gate itself.")

    gates, audits = [], []
    for donor in donors:
        ctx = pd.read_parquet(key(donor))
        audit: dict = {}
        frame = classify(donor, ctx, per_donor[donor], scale, log, args.emit_detail,
                         not args.no_top_marker, islet_chain=not args.no_islet_chain,
                         quantile=args.endocrine_quantile, audit=audit)
        audits.append(audit)
        report_identity_argmax(audit)
        report_islet_audit(audit)
        if args.gate:
            g = gate(frame, donor)
            gates.append(g)
            print(f"\n{donor}: class {g['class_match_pct']:.4f}%  context "
                  f"{g['context_match_pct']:.4f}%   rows old {g['rows_old']:,} / "
                  f"new {g['rows_new']:,}  (only-old {g['only_old']:,}, only-new {g['only_new']:,})")
            if len(g["confusion"]):
                for (a, b), n in g["confusion"].items():
                    print(f"      {a} -> {b}: {n:,}")
        else:
            part = target / f"donor_id={donor}"
            part.mkdir(parents=True, exist_ok=True)
            frame.drop(columns=["donor_id"]).to_parquet(part / "data.parquet", index=False)
            print(f"  {donor}: {len(frame):,} cells  " +
                  "  ".join(f"{k} {v:,}" for k, v in
                            frame['cell_class'].value_counts().head(5).items()), flush=True)

    if args.audit_out:
        out = pathlib.Path(args.audit_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(audits, indent=1))
        print(f"\naudit -> {out}")

    if args.gate and gates:
        table = pd.DataFrame([{k: v for k, v in g.items() if k != "confusion"} for g in gates])
        pd.set_option("display.width", 200)
        print("\n" + table.to_string(index=False))
        worst = table["class_match_pct"].min()
        print(f"\nGATE: worst donor {worst:.4f}% on cell_class. "
              f"{'PASS' if worst >= 99.99 else 'NOT YET EXACT -- do not apply corrections'}")
    elif not args.gate:
        print(f"\n-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
