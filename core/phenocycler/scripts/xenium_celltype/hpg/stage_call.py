#!/usr/bin/env python
"""Cluster one (donor, arm) embedding across the k sweep and gate every group against its null.

The k sweep is granularity, not an answer. Annotation is by null-gated enrichment afterwards,
so a finer k costs resolution and never invents a lineage; what changes with k is whether a
small population is its own group or a minority inside a larger one. The sweep is reported
as a stability curve -- the lineage shares as a function of k, each with its binomial CI --
rather than collapsed to a chosen k, because choosing k by the composition it produces would
be tuning toward an expected frequency.

The sweep goes finer than the workstation version's [20, 30, 40] for a reason that is only
available at full n: at 1.1M cells, k=320 still leaves ~3,400 cells per group on average,
far above nullgate's 50-cell floor, so the granularity is limited by the biology rather than
by the pooling. At 50,000 cells k=320 would have averaged 156 cells per group and the floor
would have started refusing groups.

Unavailable is a result, and it is never redistributed across the lineages that did clear.
Note the reason is NOT "acinar has no axis": 02_pancreas_exocrine carries CTRL and CUZD1,
both flagged acinar in the panel file and both present, with CUZD1 detected in 65.8% of
donor 6539's cells and supplying 57% of that core's mean. What is absent is the canonical
zymogen set. See run_donor's module docstring for the measurement and the narrow claim.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import numpy as np                                    # noqa: E402
import pandas as pd                                   # noqa: E402

import run_donor                                      # noqa: E402
import xenio                                          # noqa: E402
from nullgate_fast import annotate_groups             # noqa: E402
from stage_prep import load_at_site                   # noqa: E402

#: Finer than the workstation sweep because full n supports it; see the module docstring.
DEFAULT_K = (20, 40, 80, 160, 320)
#: 50 draws puts the null's 95th percentile on 2.5 draws a side. 200 puts it on 10. The
#: linear rewrite in nullgate_fast makes the difference cost nothing, so there is no longer a
#: reason to accept a noisy threshold.
DEFAULT_DRAWS = 200


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--arm", required=True, choices=["nonspatial", "spatial"])
    ap.add_argument("--k", type=int, nargs="+", default=list(DEFAULT_K))
    ap.add_argument("--n-draws", type=int, default=DEFAULT_DRAWS)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)

    man = json.loads((a.run / "work_manifest.json").read_text())
    seed = int(man["kodama"]["seed"])
    repo = HERE.parents[2]
    dest = a.run / f"donor_id={a.donor}" / f"roi={a.roi}"
    out_calls = dest / f"group_calls_{a.arm}.csv"
    out_comp = dest / f"composition_{a.arm}.csv"
    if out_calls.is_file() and out_comp.is_file() and not a.force:
        print(f"[call] {a.donor}/{a.arm} already complete — skipping")
        return 0

    emb_path = dest / f"embedding_{a.arm}.csv"
    if not emb_path.is_file():
        raise SystemExit(f"[call] {emb_path} is absent; the merge for this arm did not finish")

    from phenocycler.config import load_config
    from phenocycler.integration.vocab import load_panel_taxonomy
    cfg = load_config(str(repo / "config.ini"))
    tax = load_panel_taxonomy(cfg, "pancreas")
    cores = {lin: tax.lineage_genes(lin)
             for lin in sorted({v for v in tax.panel_to_lineage.values() if v})}

    dd = load_at_site(a.donor, a.roi, repo / "data" / "integration" / "xenium_paths.csv")
    emb = pd.read_csv(emb_path).set_index("cell_id").loc[dd.cell_id].reset_index()
    Xn = xenio.lognorm(dd.counts)

    rows, comps = [], []
    for k in a.k:
        groups = run_donor.cluster(emb, k=k, seed=seed)
        ann = annotate_groups(Xn, dd.genes, groups, cores, n_draws=a.n_draws, seed=seed)
        ann.insert(0, "k", k); ann.insert(0, "arm", a.arm); ann.insert(0, "donor_id", a.donor)
        rows.append(ann)
        comp = run_donor.composition(ann, groups)
        comp.insert(0, "k", k); comp.insert(0, "arm", a.arm); comp.insert(0, "donor_id", a.donor)
        comps.append(comp)
        unavail = float(comp.loc[comp.lineage == "Unavailable", "frac"].sum())
        smallest = int(pd.Series(groups).value_counts().min())
        print(f"[call] {a.donor} arm={a.arm} k={k}: "
              f"{int((ann.lineage != '').sum())}/{len(ann)} groups cleared their null; "
              f"Unavailable {unavail:.1%}; smallest group {smallest:,} cells")

    pd.concat(rows, ignore_index=True).to_csv(out_calls, index=False)
    pd.concat(comps, ignore_index=True).to_csv(out_comp, index=False)
    (dest / f"call_{a.arm}.json").write_text(json.dumps(
        {"complete": True, "donor_id": a.donor, "roi": a.roi, "arm": a.arm,
         "k": list(a.k), "n_draws": a.n_draws, "seed": seed,
         "n_cells": int(len(dd.cell_id)), "n_genes": len(dd.genes),
         "lineages_scored": sorted(cores), "min_group_cells": 50}, indent=2))
    print(f"[call] -> {out_calls.name}, {out_comp.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
