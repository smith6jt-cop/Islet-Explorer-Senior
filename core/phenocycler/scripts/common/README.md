# `scripts/common/` — shared, model-agnostic lane utilities (Track E.3, 2026-09-05)

Moved here from the closed DeepCellTypes lane (`scripts/dct/`) so that a new typing pilot
can reuse them without importing from a closed lane, and so that the DCT review document's
"discard `scripts/dct/`" instruction is executable. Nothing under `scripts/**` is hashed into
anything: `phenocycler/pipeline.py` computes no `run_id` at all any more (the content-addressed
run framework it used to feed was removed on 2026-09-16).

| module | what | importers |
|---|---|---|
| `gates.py` | `Verdict`, `require_evidence`, `assess_self_match`, `assess_anchor_concordance` — fail-closed evidence primitives (unjudged is never a pass) | `dct/concordance_gate.py`, `dct/plane_identity.py`, `broadtype/` |
| `tiff_window.py` | `read_window`, `selftest_bit_exact` — windowed reads of tiled/stripped TIFF pages without zarr; each worker must open its own `TiffFile` | `dct/arms.py`, `dct/feasibility_tile.py`, `restain/tiles.py`, `broadtype/` |
| `window_mask.py` | `Window`, `rasterize_windows`, `verify_against_archive` — per-window label masks from the QuPath GeoJSON stream with the invariant `label = feature_index + 1`, so `object_ids[label - 1]` is the cell's UUID | `dct/arms.py`, `broadtype/frame.py` |
| `label_space.py` | `LabelSpace.restrict`, `load_label_space` — the restrict-and-renormalise idiom over a JSON crosswalk that carries a written reason per class. Its default `CROSSWALK` still points at `scripts/dct/broad_class_crosswalk.json`, the DCT lane's artifact; other lanes pass their own crosswalk path | `dct/arms.py`, `dct/concordance_gate.py`, `broadtype/arms/dct_*.py` |

Tests: `../../.venv/bin/python -m pytest -q scripts/common/tests` (synthetic TIFFs and geometry,
no cohort data). The real-data bit-identity record for `window_mask` remains
`data/development/deepcell_types_review/window_mask_verification.json` (32 windows, four donors).

Consumers add `scripts/common` to `sys.path` beside their own directory (the lane idiom), e.g.
`sys.path.insert(0, str(HERE.parent / "common"))`.
