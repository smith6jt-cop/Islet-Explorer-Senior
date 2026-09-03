"""Label-free cross-check on the declared marker-exclusivity classification.

`marker_exclusivity.json` says, from biology, which non-lineage markers carry compartment
information by exclusivity.  That is a written judgement, so it wants an independent check
that does not come from the same head.

This one uses **only this run's own expression parquet**: no cell-type label, no model
output, no threshold, no expected frequency.  The question it answers is narrow and purely
distributional -- *does this marker partition the tissue at all, or is it on everything?* --
via the Gini coefficient of its per-cell intensities, which is scale-invariant and so
comparable across markers measured on different antibodies.

A marker present on 5 % of cells and absent elsewhere has a high Gini.  A housekeeping
protein present on all of them has a low one.

**It is descriptive, not decisive, and the measurement says so.**  Genuine lineage markers
span the whole range here -- Caveolin 0.753 down to TP63 0.134, FOXP3 0.202, CD4 0.203 --
because a low Gini can mean "rare or weakly stained" as easily as "on everything", and a
high one can mean "proliferating" (Ki67) rather than "a particular compartment".  So a high
Gini is neither necessary nor sufficient for compartment information.  This check can
contradict the declared list and prompt a re-read; it can never replace it, and
disagreements are reported and left standing rather than reconciled.

It has already earned its place once: it caught Caveolin declared uninformative while
carrying the highest Gini in the panel, which turned out to be a labelling error -- two
different reasons for exclusion collapsed into one flag.

Selection donors only (6539, 6450, 6523); 115, 6521, 6593, 6534 and 6591 are held out.

    python scripts/dct/marker_exclusivity_check.py          # from core/phenocycler/
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXPRESSION = ROOT / "core/phenocycler/data/runs/30f6689e4d27d295f890/expression"
OUT = ROOT / "data/development/deepcell_types_review"
SELECTION = ["6539", "6450", "6523"]


def gini(x):
    """Scale-invariant concentration of a non-negative vector.  0 = flat, 1 = all in one cell."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 100:
        return np.nan
    x = x - min(x.min(), 0.0)            # a shifted floor, not a threshold: Gini needs x >= 0
    if x.sum() <= 0:
        return np.nan
    x = np.sort(x)
    n = x.size
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def main():
    decl = json.loads((HERE / "marker_exclusivity.json").read_text())["entries"]
    reg = json.loads((ROOT / "core/phenocycler/phenocycler/marker_registry.json").read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    kind = {m["name"]: m["kind"] for m in reg}

    per = {}
    for d in SELECTION:
        df = pd.read_parquet(EXPRESSION / f"donor_id={d}" / "data.parquet")
        for c in df.columns:
            if kind.get(c) is None:
                continue
            per.setdefault(c, {})[d] = gini(df[c].to_numpy(float))

    rows = []
    for m, byd in per.items():
        vals = [v for v in byd.values() if np.isfinite(v)]
        if not vals:
            continue
        rows.append({"marker": m, "kind": kind[m], "gini_mean": float(np.mean(vals)),
                     "gini_min": float(np.min(vals)), "gini_max": float(np.max(vals)),
                     "declared_informative": decl.get(m, {}).get("lineage_informative")})
    rows.sort(key=lambda r: -r["gini_mean"])

    # The type markers are the reference class: they are lineage markers by construction, so
    # their Gini range is what "partitions the tissue" looks like in this dataset.
    ref = [r["gini_mean"] for r in rows if r["kind"] == "type"]
    lo, hi = float(np.percentile(ref, 10)), float(np.percentile(ref, 90))
    print(f"type-marker Gini reference range (10th-90th pct over {len(ref)} markers): {lo:.3f}-{hi:.3f}\n")
    print(f"{'marker':<18}{'kind':<12}{'Gini (mean of 3)':>18}{'declared':>12}   note")
    disagree = []
    for r in rows:
        d = r["declared_informative"]
        note = ""
        if d is True and r["gini_mean"] < lo:
            note = "DISAGREES: declared informative, flatter than 90% of lineage markers"
            disagree.append(r)
        elif d is False and r["gini_mean"] > hi:
            note = "DISAGREES: declared uninformative, yet partitions more than most lineage markers"
            disagree.append(r)
        print(f"{r['marker']:<18}{r['kind']:<12}{r['gini_mean']:>18.3f}"
              f"{str(d):>12}   {note}")

    out = {"method_version": "dct-marker-exclusivity-gini-2026-08-24",
           "statistic": "Gini coefficient of per-cell intensity; label-free, model-free, scale-invariant",
           "selection_donors": SELECTION,
           "type_marker_reference_range": [lo, hi],
           "descriptive_not_decisive": ("Genuine lineage markers span the whole range (Caveolin "
                                        "0.753 to TP63 0.134), so a high Gini is neither necessary "
                                        "nor sufficient for compartment information. This check can "
                                        "contradict the declared list and prompt a re-read; it can "
                                        "never replace it."),
           "rows": rows, "disagreements": disagree}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "marker_exclusivity_gini.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n{len(disagree)} disagreement(s) with the declared list; written to "
          f"{OUT / 'marker_exclusivity_gini.json'}")


if __name__ == "__main__":
    main()
