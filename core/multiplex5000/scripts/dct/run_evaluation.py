"""C2: the release verdict, on the five donors nothing has been fitted to.

The configuration was chosen on 6539, 6450 and 6523 -- checkpoint `2026-06-23-ptft`, all 52
usable channels, a 2048-px window scoring a 1024-px core, the cell mask, insulin as
`C-peptide`, and the 26-class label space.  Those three donors have seen the gate and cannot
un-see it, so the verdict is computed here instead, on 115, 6521, 6593, 6534 and 6591.

**One configuration and one checkpoint, deliberately.**  Running `2026-06-15` here as well
would cost little and answer nothing worth spending the held-out set on: the selection
evidence is 3 donors x 3 regions x 3 panel variants, every statistic, in one direction.
Evaluating the rejected alternative on held-out data invites switching to it after the fact,
which is how a held-out set stops being one.

All 8 regions per donor -- the full frozen frame, not the 3 the screen used.

    TQDM_DISABLE=1 python scripts/dct/run_evaluation.py     # from core/multiplex5000/
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms import ARMS, OUT, build_masks, run_arm  # noqa: E402

EVALUATION = ["115", "6521", "6593", "6534", "6591"]
REGIONS = list(range(8))
CHECKPOINT = "2026-06-23-ptft"
CONFIG = ARMS[0]                     # A0-baseline IS the selected configuration
LOG = OUT / "evaluation_c2.jsonl"


def main():
    assert CONFIG.name == "A0-baseline" and CONFIG.panel == "all" and CONFIG.window == 2048 \
        and CONFIG.mask == "cell" and CONFIG.ins == "C-peptide" and not CONFIG.with_dapi, \
        "the frozen configuration has drifted from what was selected"

    OUT.mkdir(parents=True, exist_ok=True)
    done = set()
    if LOG.exists():
        for line in LOG.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done.add((r["donor"], r["region"]))
        print(f"resuming: {len(done)} records already on disk")

    t0 = time.time()
    for donor in EVALUATION:
        todo = [k for k in REGIONS if (donor, k) not in done]
        if not todo:
            continue
        print(f"\n=== {donor}", flush=True)
        t = time.time()
        cache = build_masks(donor, [CONFIG], todo)
        print(f"    {len(cache)} windows in {time.time() - t:.0f}s", flush=True)
        try:
            recs = run_arm(donor, CONFIG, CHECKPOINT, todo, cache)
        except Exception as exc:
            print(f"  !! {donor}: {type(exc).__name__}: {exc}", flush=True)
            recs = [{"donor": donor, "region": k, "arm": CONFIG.name, "checkpoint": CHECKPOINT,
                     "error": f"{type(exc).__name__}: {exc}"} for k in todo]
        with LOG.open("a") as fh:
            for r in recs:
                fh.write(json.dumps(r) + "\n")
        del cache
    print(f"\nC2 complete in {(time.time() - t0) / 60:.0f} min -> {LOG}")


if __name__ == "__main__":
    main()
