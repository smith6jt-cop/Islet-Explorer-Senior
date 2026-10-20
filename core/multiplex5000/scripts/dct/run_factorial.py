"""C1 screening: which confound, if any, moves the concordance gate.

    TQDM_DISABLE=1 python scripts/dct/run_factorial.py        # from core/multiplex5000/

Three donors, one per panel variant and covering both mask provenances, so an arm that
only works on one panel is visible here rather than after the confirmation set is spent:

    6539  batch 1, cd38-59ch,        archived mask
    6450  batch 2, b-catenin1-59ch,  rasterized
    6523  batch 2, neither-58ch,     rasterized  (every index from 16 shifts)

Three regions each, stated in advance: the first two blind hash draws and the first
islet-anchored one, so the screen always sees an islet field -- the study question is the
peri-islet immune environment and the blind draw can miss islets entirely.

Both checkpoints for every arm, because the arms demonstrably interact with the checkpoint:
restricting the label space helps ptft and does not repair 2026-06-15.

Results append to `factorial_c1.jsonl` as they arrive, so a run that dies partway leaves
evidence rather than nothing.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms import ARMS, CHECKPOINTS, OUT, build_masks, run_arm  # noqa: E402

DONORS = ["6539", "6450", "6523"]
REGIONS = [0, 1, 4]                      # two blind draws, one islet-anchored
LOG = OUT / "factorial_c1.jsonl"


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
    for donor in DONORS:
        print(f"\n=== {donor} — one GeoJSON pass for every window every arm needs", flush=True)
        t = time.time()
        cache = build_masks(donor, ARMS, REGIONS)
        print(f"    {len(cache)} windows in {time.time() - t:.0f}s", flush=True)
        for ckpt in CHECKPOINTS:
            for arm in ARMS:
                todo = [k for k in REGIONS if (donor, arm.name, ckpt, k) not in done]
                if not todo:
                    continue
                try:
                    recs = run_arm(donor, arm, ckpt, todo, cache)
                except Exception as exc:                    # an arm that dies is a finding
                    print(f"  !! {donor}/{arm.name}/{ckpt}: {type(exc).__name__}: {exc}", flush=True)
                    recs = [{"donor": donor, "region": k, "arm": arm.name, "checkpoint": ckpt,
                             "error": f"{type(exc).__name__}: {exc}"} for k in todo]
                with LOG.open("a") as fh:
                    for r in recs:
                        fh.write(json.dumps(r) + "\n")
        del cache
    print(f"\nC1 complete in {(time.time() - t0) / 60:.0f} min -> {LOG}")


if __name__ == "__main__":
    main()
