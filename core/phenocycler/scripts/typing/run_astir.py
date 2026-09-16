#!/usr/bin/env python3
"""Cell typing with ASTIR -- a probabilistic marker-based typer, run in ``astir-env``.

WHY THIS AND NOT A HAND-BUILT RULE.  High background and nonspecific staining is a routine problem
in multiplex imaging, and ASTIR is built for it: it models, per marker, both a background level and
a cell-type-specific level, and infers per-cell type probabilities by variational inference with a
native ``Unknown`` outcome.  Nothing here chooses an operating point, a lift, a positive fraction
or a family-wise bar -- those were an attempt to rebuild what this model already does properly.

It takes exactly the marker -> cell-type dictionary the lane already declares
(:mod:`lineages`), so the biology stated up front is the model's input rather than something
reconstructed from a clustering afterwards.

``diagnostics_celltype`` is the part that speaks to the background problem directly: for every
cell type and every one of its markers it tests whether that marker is expressed higher in that
type than in the types that do NOT declare it, and reports the ones that fail.  That is the
measured, model-based version of "this channel is background", and it is reported with the run
rather than being decided in advance.

Input is the Nimbus per-marker confidence table by default (``--input confidence``); pass
``--input intensity`` to run the same model on the REDSEA intensities from the substrate, which
matters because the two disagree on at least one channel -- on donor 6436 the INS confidence has
no operating point that exceeds 15x islet enrichment while the intensity's top 0.05% is entirely
inside islet polygons.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"
ASTIR_PY = Path.home() / "miniconda3/envs/astir-env/bin/python"

WORKER = r'''
import json, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
job = json.load(open(sys.argv[1]))
expr = pd.read_parquet(job["expr"])
ids = expr.pop("object_id")
marker_dict = {"cell_type": job["marker_dict"]}
from astir import Astir
a = Astir(input_expr=expr, marker_dict=marker_dict, random_seed=job["seed"])
a.fit_type(max_epochs=job["max_epochs"], learning_rate=job["lr"],
           batch_size=job["batch_size"], n_init=job["n_init"])
prob = a.get_celltype_probabilities()
calls = a.get_celltypes(threshold=job["threshold"])
out = pd.DataFrame({"object_id": ids.to_numpy(),
                    "astir_type": np.asarray(calls).ravel()})
for c in prob.columns:
    out[f"p__{c}"] = prob[c].to_numpy()
out.to_parquet(job["out"], index=False)
diag = a.diagnostics_celltype()
diag.to_csv(job["diagnostics"], index=False)
json.dump({"n_cells": int(len(out)), "types": sorted(map(str, set(out.astir_type))),
           "n_diagnostics": int(len(diag))}, open(job["result"], "w"), indent=2)
print(f"astir: {len(out):,} cells, {len(diag)} marker diagnostics", flush=True)
'''


def marker_dict() -> dict:
    sys.path.insert(0, str(Path(__file__).parent))
    import lineages as L
    return {lin: list(ms) for lin, ms in L.LINEAGE_MARKERS.items()}


def run(donor: str, source: str, threshold: float, epochs: int, seed: int) -> Path:
    sys.path.insert(0, str(Path(__file__).parent))
    import pandas as pd
    import lineages as L

    md = marker_dict()
    wanted = sorted({m for ms in md.values() for m in ms})

    if source == "confidence":
        df = pd.read_parquet(DEV / f"positivity/donor_id={donor}/data.parquet")
        cols = [m for m in wanted if m in df.columns]
        expr = df[["object_id", *cols]].copy()
    else:
        df = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet")
        cols = [m for m in wanted if f"v__{m}" in df.columns]
        expr = df[["object_id"] + [f"v__{m}" for m in cols]].copy()
        expr.columns = ["object_id", *cols]
        # ASTIR expects expression on a scale where a background and a positive level are
        # separable; raw intensity is heavy-tailed, so use the model's usual arcsinh.
        import numpy as np
        expr[cols] = np.arcsinh(expr[cols].to_numpy(float) / 5.0)

    sub = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet",
                          columns=["object_id", "eligible"])
    expr = expr.merge(sub, on="object_id").query("eligible").drop(columns="eligible")
    expr = expr.dropna()
    md = {lin: [m for m in ms if m in cols] for lin, ms in md.items()}
    md = {lin: ms for lin, ms in md.items() if ms}

    dest = DEV / f"astir/donor_id={donor}"
    dest.mkdir(parents=True, exist_ok=True)
    expr_path = dest / f"expr_{source}.parquet"
    expr.to_parquet(expr_path, index=False)

    job = {"expr": str(expr_path), "marker_dict": md, "seed": seed,
           "max_epochs": epochs, "lr": 2e-3, "batch_size": 1024, "n_init": 3,
           "threshold": threshold,
           "out": str(dest / f"calls_{source}.parquet"),
           "diagnostics": str(dest / f"diagnostics_{source}.csv"),
           "result": str(dest / f"result_{source}.json")}
    (dest / f"job_{source}.json").write_text(json.dumps(job, indent=1))
    worker = dest / "_worker.py"
    worker.write_text(WORKER)

    print(f"[{donor}] ASTIR on {len(expr):,} cells x {len(cols)} markers "
          f"({source}); {len(md)} cell types", flush=True)
    r = subprocess.run([str(ASTIR_PY), str(worker), str(dest / f"job_{source}.json")],
                       capture_output=True, text=True)
    print(r.stdout[-2000:])
    if r.returncode != 0:
        raise RuntimeError(f"astir failed:\n{r.stderr[-3000:]}")
    return dest / f"calls_{source}.parquet"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--input", choices=("confidence", "intensity"), default="confidence")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="ASTIR's own probability threshold for assigning a type; below it a "
                         "cell is Unknown, which is the model's own outcome and not a marker cut")
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1234)
    a = ap.parse_args()
    for d in a.donors:
        p = run(d, a.input, a.threshold, a.epochs, a.seed)
        import pandas as pd
        t = pd.read_parquet(p)
        print(f"\n[{d}] {p}")
        for v, n in t["astir_type"].value_counts().items():
            print(f"   {v:<16} {n:>9,}  {n/len(t)*100:>5.1f}%")


if __name__ == "__main__":
    main()
