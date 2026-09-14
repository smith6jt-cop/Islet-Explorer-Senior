#!/usr/bin/env python3
"""Cell typing with TACIT -- learns one threshold per cell type from the data, in R 4.5.1.

WHY THIS TOOL FOR THIS PROBLEM.  High background is routine in multiplexed imaging and TACIT is
built to derive the cut rather than be handed one: it CLR-normalises, builds microclusters, scores
each cell against the declared TYPE x MARKER signature, then fits a segmented model over the
ordered microcluster medians and takes the cut that minimises the two misclassification rates
between the resulting high and low groups.  The threshold is a property of the score's own
bimodality in microcluster space; no expected frequency enters, and nothing in this file chooses
an operating point.

KNOWN LIMITATION OF THE TOOL, recorded rather than hidden: TACIT's microclustering calls Seurat's
``FindClusters`` with no ``algorithm`` argument, so it runs **Louvain**, which can return badly
connected and even disconnected communities -- the defect Leiden was introduced to fix.  That is
inside the published package and not a choice made here; this lane's own clustering uses Leiden.
Read a TACIT result as a second opinion, not as an arbiter.

The R half is the lane's existing, carefully documented ``arms/tacit_run.R``, reused unchanged --
it already replicates TACIT's preprocessing verbatim, calibrates the Louvain resolution to the
package's own microcluster-size guidance, and checks that the binaries it recomputes match the
package's cell for cell.  This module only supplies the signature and the expression.

INPUT SCALE MATTERS HERE.  TACIT's demotion step uses an absolute internal constant, so it expects
a normalised scale; feeding raw intensities is a known way to get it wrong.  The default input is
therefore the Nimbus per-marker confidence, which is already on a common [0, 1] scale.  Pass
``--input intensity`` to run it on arcsinh REDSEA values instead -- worth doing because the two
sources disagree on at least one channel that matters (6436's INS).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"
R_SCRIPT = REPO / "core/phenocycler/scripts/broadtype/arms/tacit_run.R"
SCRATCH = Path("/data/typing_tacit")

#: TACIT's own knobs, as the existing lane pinned them.  ``r = None`` calibrates the resolution to
#: the package's stated microcluster-size band rather than fixing it.
PARAMS = {"p": 10, "r": None, "seed": 0,
          "r_ladder": [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 80, 120],
          "band_lo": 0.001, "band_hi": 0.005,
          #: EVERY eligible cell.  The 200,000-cell subsample carried over from the earlier lane
          #: was that driver's choice, not TACIT's requirement -- the package takes whatever
          #: `data_expression` it is handed.  On 64 cores and ~490 GB free there is no reason to
          #: learn a threshold from a fraction of a slide that fits in memory whole.
          "n_subsample": None}


def signature(markers: list[str]) -> pd.DataFrame:
    sys.path.insert(0, str(Path(__file__).parent))
    import lineages as L
    rows = []
    for lin, ms in L.LINEAGE_MARKERS.items():
        present = [m for m in ms if m in markers]
        if not present:
            continue
        rows.append({"cell_type": lin, **{m: int(m in present) for m in markers}})
    return pd.DataFrame(rows)


def load_expression(donor: str, source: str) -> tuple[pd.DataFrame, list[str]]:
    sys.path.insert(0, str(Path(__file__).parent))
    import lineages as L
    wanted = L.identity_markers()
    if source == "confidence":
        df = pd.read_parquet(DEV / f"positivity/donor_id={donor}/data.parquet")
        cols = [m for m in wanted if m in df.columns]
        expr = df[["object_id", *cols]].copy()
    else:
        df = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet")
        cols = [m for m in wanted if f"v__{m}" in df.columns]
        expr = df[["object_id"] + [f"v__{m}" for m in cols]].copy()
        expr.columns = ["object_id", *cols]
        expr[cols] = np.arcsinh(expr[cols].to_numpy(float) / 5.0)
    sub = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet",
                          columns=["object_id", "eligible"])
    expr = expr.merge(sub, on="object_id").query("eligible").drop(columns="eligible")
    expr = expr.dropna().reset_index(drop=True)
    # the R half expects values already left-censored at zero
    expr[cols] = np.clip(expr[cols].to_numpy(float), 0, None)
    return expr, cols


def run(donor: str, source: str, seed: int) -> Path:
    expr, cols = load_expression(donor, source)
    sig = signature(cols)
    if PARAMS["n_subsample"] is None:
        in_sub = np.ones(len(expr), bool)          # learn on every cell
        n_sub = len(expr)
    else:
        rng = np.random.default_rng(seed)
        n_sub = min(int(PARAMS["n_subsample"]), len(expr))
        in_sub = np.zeros(len(expr), bool)
        in_sub[rng.choice(len(expr), n_sub, replace=False)] = True
    expr.insert(1, "in_subsample", in_sub)

    scratch = SCRATCH / f"{donor}__{source}"
    scratch.mkdir(parents=True, exist_ok=True)
    expr.to_parquet(scratch / "input.parquet", index=False)
    sig.to_csv(scratch / "signature.csv", index=False)
    (scratch / "params.json").write_text(json.dumps(
        {"seed": seed, "p": PARAMS["p"], "r": PARAMS["r"],
         "r_ladder": PARAMS["r_ladder"],
         "band_lo": PARAMS["band_lo"], "band_hi": PARAMS["band_hi"]}))

    print(f"[{donor}] TACIT on {len(expr):,} cells x {len(cols)} markers ({source}); "
          f"{len(sig)} cell types; learning on {n_sub:,}", flush=True)
    t0 = time.time()
    r = subprocess.run(["Rscript", str(R_SCRIPT), str(scratch)], capture_output=True, text=True)
    print(r.stdout[-2500:])
    if r.returncode != 0:
        fail = scratch / "failure.json"
        raise RuntimeError(f"tacit failed (rc={r.returncode}):\n"
                           f"{fail.read_text() if fail.exists() else ''}\n{r.stderr[-3000:]}")
    print(f"[{donor}] TACIT finished in {(time.time()-t0)/60:.1f} min", flush=True)

    scores = pd.read_parquet(scratch / "scores.parquet")
    dest = DEV / f"tacit/donor_id={donor}"
    dest.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(dest / f"calls_{source}.parquet", index=False)
    for name in ("learned.json", "tacit_stdout.log"):
        if (scratch / name).exists():
            (dest / f"{Path(name).stem}_{source}{Path(name).suffix}").write_text(
                (scratch / name).read_text())
    return dest / f"calls_{source}.parquet"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--input", choices=("confidence", "intensity"), default="confidence")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    for d in a.donors:
        p = run(d, a.input, a.seed)
        t = pd.read_parquet(p)
        col = "tacit_label" if "tacit_label" in t.columns else t.columns[-1]
        print(f"\n[{d}] {p}")
        for v, n in t[col].value_counts().items():
            print(f"   {str(v):<16} {n:>9,}  {n/len(t)*100:>5.1f}%")


if __name__ == "__main__":
    main()
