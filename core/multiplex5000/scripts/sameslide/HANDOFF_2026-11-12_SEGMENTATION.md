# HANDOFF 2026-11-12: InstanSeg vs Xenium segmentation comparison (sameslide lane)

Status: code complete and reviewed; tested only against a stub harness; **never run on real data**.
Origin: claude.ai cloud session https://claude.ai/code/session_014Hj88vwdzaGt1d4LtCGBHm (2026-11-03/04).
The four code files in this package were also attached in that chat on 2026-11-04.

## Why this exists

The question (2026-11-03): has code been developed to compare Xenium segmentation to the
InstanSeg/PhenoCycler segmentation for nuclei and cell bodies? Survey of the repo said no:

- `multiplex5000/integration/boundaries.py` loads both platforms' outlines and warps them, and
  `integration/sameslide.py::match_by_iou` pairs whole cells at IoU >= 0.25 — but as a pairing
  step for downstream transfer, not an evaluation (no thresholds sweep, no split/merge, no
  nucleus outlines, no saved output found).
- The sameslide lane (`scripts/sameslide/`) matches nuclei by **centroid only**
  (`register_points.py`, `recheck_nuclei.py`) to validate registration (median pair offsets
  0.63–1.11 µm across the four donors, table in `HANDOFF_2026-10-19.md`).
- `build_joined.py`/`lift_calls.py` compute the full overlap graph (`overlap_edges.parquet`) but
  use it only to carry protein/calls onto Xenium cells. `paired_typing.py` (2026-11-02) matches
  by centroid within 5 µm and marks the segmentation relationship "unavailable without overlap
  evidence"; `docs/evaluations/blocker_remediation_20261102/geom2.md` lists that evidence as
  open. This stage is that evidence.

## What this package adds

Two modules plus tests, all in `core/multiplex5000/scripts/sameslide/`:

- **`segcompare.py`** — standalone metric library (numpy/pandas/shapely 2.x/scipy only; no lane
  imports; reusable for restain or MACSima comparisons). Given two polygon sets in one frame:
  - overlap graph: every positive-area intersection with IoU and per-object covered fractions;
  - one-to-one matching at IoU 0.1–0.9: maximum pairs, then maximum total IoU, solved exactly
    with a sparse assignment (greedy loses pairs; verified against brute force on 3,000 random
    graphs during review);
  - panoptic quality at IoU 0.5: SQ = mean matched IoU, RQ = (matched_a+matched_b)/(n_a+n_b)
    (symmetric F1 — neither platform is ground truth, every count is reported from both sides),
    PQ = SQ×RQ;
  - splits/merges by strict-majority containment (a partner covering > half the object's area);
  - per-object unmatched reason at IoU 0.5: `no_overlap` (detection difference) > `split` >
    `merged` > `boundary`;
  - matched-pair geometry: log2 area ratio b/a, centroid offset (a mean offset away from zero is
    residual registration, not segmentation);
  - registration ceiling: each object's IoU with itself shifted by the recheck's median pair
    distance in a seeded random direction — what two *identical* segmentations would score at
    that misregistration; read observed IoU against this, not against 1. Conservative (the
    recheck distance already contains segmentation disagreement);
  - counting domain `shared_area()`: 10 µm occupancy grid per platform, gaps < ~60 µm closed,
    intersection eroded 20 µm; objects outside it stay matchable but leave every denominator.
    Area only one platform segmented is reported (`area_mm2.a_only/b_only`), never silently
    dropped;
  - per-stratum summaries, `one_side()` for masks only one side has, `hierarchy()` for whether
    matched nuclei sit in matched cells.
- **`segmentation_four.py`** — the lane stage, patterned on `join_four.py` (same receipts,
  fingerprint, `attempts/<id>/`, `current.json`, `scientifically_unavailable` propagation).
  Frame: PhenoCycler; Xenium nuclei and cells warped forward through the validated transform
  (displacement field included — never inverted). Objects: everything each platform segmented
  (preparation `geometry.parquet` `wkb`/`nucleus_wkb`; bundle `cell_boundaries.parquet` and
  `nucleus_boundaries.parquet` split on `label_id`, so a cell's second nucleus is its own
  object). Strata: tissue from the preparation's routing (`panc`→pancreas, `pln_1`→pLN; Xenium
  objects take the nearest InstanSeg cell's tissue), plus Xenium `segmentation_method` reported
  separately (an expanded nucleus is not a measured outline). Also: nucleus-in-matched-cell
  hierarchy, nucleus/cell area ratio per platform, structural self-checks that gate
  `validated_complete`, and a per-section PNG (matched fraction vs IoU with the ceiling curve,
  and matched-or-why-not bars).

Outputs per section under `<section>/segmentation/attempts/<id>/`: `segmentation.json`,
`thresholds.csv` (tidy: stratum × threshold × level), `edges_{nucleus,cell}.parquet`,
`objects_{instanseg,xenium}_{nucleus,cell}.parquet` (per-object: area, partner, pair IoU,
reason, held_by/holds_n, matched flag at every threshold — join to typing via
`object_id`/`cell_id`), `segmentation.png`, `receipt.json`.

## Verified against the real lane files (2026-11-12, staged copies from 2026-11-03)

- Imports follow `join_four.py`'s conventions: path root = `core/multiplex5000` (`HERE.parents[1]`),
  `from phenocycler.integration.transform import Transform`, `from manifest import
  load_manifest, get_section`, `from workflow import DEFAULT_OUT, section_root, identity,
  write_json, now, valid_receipt`, `from join_four import TISSUE, current`, `from build_joined
  import xenium_polygons` (defined at build_joined.py:87).
- `workflow.py` provides every helper named above; `TISSUE`/`current` at join_four.py:53/56.
- `recheck_nuclei.py:75` writes exactly `summary['pair_p50_um']`; the stage handles a missing or
  None value (ceiling columns NaN, JSON `null` — regression-tested).
- The recheck receipt is checked against the preparation sha before running, as in `join_four`.

## Check before the first real run

1. `grep -n "_read_boundary_rings" core/multiplex5000/phenocycler/integration/boundaries.py` —
   the stage imports it from the `phenocycler` package. It exists at line 191 of the
   `multiplex5000/integration/boundaries.py` copy; if the `phenocycler` package differs, point
   the import at the right module.
2. Preparation `geometry.parquet` `nucleus_wkb` holds **polygons**, not points:
   `shapely.get_type_id(shapely.from_wkb(df.nucleus_wkb.dropna().head(100)))` → all 3 or 6.
3. `prep/panc/cells.parquet` and `prep/pln_1/cells.parquet` carry `object_id` (used for tissue
   routing).
4. Bundle `cells.parquet` has `segmentation_method` (otherwise that breakdown reads
   `not_recorded`; everything else unaffected).
5. Open question for Justin: the counting domain is my 10/30/20 µm grid construction. If the
   lane already defines "area Xenium imaged" (e.g. in `reseed_four.py`, which was never read),
   switch `shared_area` to that and keep the a_only/b_only reporting.

## Test and run

```bash
cd core/multiplex5000/scripts/sameslide
python -m pytest tests/test_segcompare.py tests/test_segmentation_four.py   # 23 tests
python segmentation_four.py --manifest <manifest.json> --section <section_id>
```

The end-to-end test builds a synthetic preparation and a vendor-format bundle (long vertex
tables, a two-nuclei cell, a cell with no nucleus, an unmatched extra cell, a pure-translation
transform), so on the real tree it exercises the real `_read_boundary_rings` and
`xenium_polygons`.

Sections: the four donors with validated transforms — 112, 6347, 6580 (panc + pLN) and 6533
(panc). 6533 pLN should record `scientifically_unavailable` (its registration is).

Sizing (cloud container, synthetic packed discs, single core): 250k objects/side → 504k edges,
33 s for one level end to end. A full section (~1.0M InstanSeg × 0.7M Xenium cells, two levels)
should sit in the minutes range; memory is dominated by the polygon arrays — a standard
HiPerGator node is ample, no GPU.

## Review trail

An independent review (same session) found, and I fixed, with regression tests: a NaN ceiling
reaching the `allow_nan=False` receipt writer (now `null`); two structural checks that could
never fail (reasons now checked against the taxonomy and against partners; one-to-one checked in
both directions; empty-headline guard added); valid non-polygonal geometries (GeometryCollection
with a stray line) now reduced to their polygon parts; `hierarchy()` restricted to the counted
area like every other number; `matched_{t}` column collisions at thresholds like 0.25 vs 0.2;
ceiling share computed over the stratum's finite values; the reuse fingerprint extended to every
module that shapes the outputs (boundaries, join_four, workflow). `match()`'s optimality,
`shift_ceiling`'s vertex ordering (holes, multipolygons), `shared_area` row/col indexing, and
all side a/b bookkeeping were verified clean.
