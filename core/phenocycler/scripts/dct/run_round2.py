"""Round two: separate panel COMPOSITION from panel SIZE, on the selection donors only.

C1 screens each axis independently.  If the reduced panel moves the gate, that alone cannot
say whether 30 lineage markers helped or 30 channels helped -- the model's channel-wise
transformer has no positional encoding, so panel size is a live explanation on its own.

Five seeded random 30-channel draws answer it.  If they move the gate as far as `type_only`
does, the effect is size and the type/state split is incidental; if they do not, the
composition is doing the work and the reduced panel is a real setting to adopt.

`A1d-type-plus-exclusive` asks the question the type/state split is too crude for: a
FUNCTIONAL marker present on some compartments and absent from others discriminates them
even though it is not a lineage marker.  HLA-DR separates B cells, DCs and macrophages from
resting exocrine epithelium; Granzyme B is CD8/NK and nothing else.  Read against `all` (52)
and `type_only` (30) it answers both halves at once -- whether those 13 markers ADD over a
bare lineage panel, and whether the 9 ubiquitous ones SUBTRACT from the full one.

`A1c` / `A1e` ask the follow-on: the model's training panels essentially always carry a
nuclear channel and ours is withheld by a crosswalk decision, so if a reduced panel is what
we adopt, does it want DAPI back?

Selection donors only -- 6539, 6450, 6523.  115, 6521, 6593, 6534 and 6591 are held out for
the release verdict and are not touched here.

    TQDM_DISABLE=1 python scripts/dct/run_round2.py         # from core/phenocycler/
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms import ARMS, CHECKPOINTS, OUT, build_masks, run_arm  # noqa: E402
from arms import ROUND2_ARMS  # noqa: E402

SELECTION = ["6539", "6450", "6523"]
REGIONS = [0, 1, 4]
LOG = OUT / "factorial_round2.jsonl"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    done = set()
    if LOG.exists():
        for line in LOG.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done.add((r["donor"], r["arm"], r["checkpoint"], r["region"]))
        print(f"resuming: {len(done)} records already on disk")

    t0 = time.time()
    for donor in SELECTION:
        print(f"\n=== {donor}", flush=True)
        # every round-2 arm reads the same 2048 window as the baseline, so one pass covers it
        cache = build_masks(donor, ROUND2_ARMS + [ARMS[0]], REGIONS)
        for ckpt in CHECKPOINTS:
            for arm in ROUND2_ARMS:
                todo = [k for k in REGIONS if (donor, arm.name, ckpt, k) not in done]
                if not todo:
                    continue
                try:
                    recs = run_arm(donor, arm, ckpt, todo, cache)
                except Exception as exc:
                    print(f"  !! {donor}/{arm.name}/{ckpt}: {type(exc).__name__}: {exc}", flush=True)
                    recs = [{"donor": donor, "region": k, "arm": arm.name, "checkpoint": ckpt,
                             "error": f"{type(exc).__name__}: {exc}"} for k in todo]
                with LOG.open("a") as fh:
                    for r in recs:
                        fh.write(json.dumps(r) + "\n")
        del cache
    print(f"\nround 2 complete in {(time.time() - t0) / 60:.0f} min -> {LOG}")


if __name__ == "__main__":
    main()
