#!/usr/bin/env python
"""Step 2 driver: KINTSUGI's SOM over the UniFORM-aligned block.

    python scripts/kintsugi/run_step2.py --donors 6551            # the pilot
    python scripts/kintsugi/run_step2.py --arm A --nodes 25        # cohort, one arm
    python scripts/kintsugi/run_step2.py --sweep 9 16 25 36 49 64  # the stability sweep

Must run in an environment with pyFlowSOM, which needs ``numpy < 2`` -- see the lane README.
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

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
sys.path.insert(0, str(HERE))

import features                                     # noqa: E402
import panel                                        # noqa: E402
import som_adapter                                  # noqa: E402

OUT = panel.REPO / "data" / "development" / "kintsugi"

def run_one(block: pd.DataFrame, markers: list[str], nodes: int, arm: str,
            seed: int) -> tuple[np.ndarray, pd.DataFrame, dict]:
    values = block[markers].to_numpy(dtype=np.float64)
    from sklearn.preprocessing import StandardScaler
    scaled = StandardScaler().fit_transform(values)
    check = features.scaler_is_feature_weighting(values, block["donor_id"].to_numpy(), scaled)
    if not check["passes"]:
        raise RuntimeError(
            f"the scaler moved donors relative to each other "
            f"(max relative change {check['max_relative_change']:.2e} > {check['tolerance']}). "
            "That is donor correction, which INTENSITY_SCALE.md forbids. Arm refused.")

    t0 = time.time()
    clusters, profile = som_adapter.som_cluster(
        values, markers, num_nodes=nodes,
        num_passes=som_adapter.NOTEBOOK_SOM["num_passes"],
        lr_start=som_adapter.NOTEBOOK_SOM["lr_start"],
        lr_end=som_adapter.NOTEBOOK_SOM["lr_end"], seed=seed)
    elapsed = time.time() - t0

    purity = features.donor_purity(clusters, block["donor_id"].to_numpy())
    return clusters, profile, {
        "arm": arm, "num_nodes": nodes, "seed": seed,
        "n_clusters_occupied": int(len(np.unique(clusters))),
        "seconds": round(elapsed, 1),
        "scaler_check": check,
        "donor_purity_max": float(purity["purity"].max()),
        "donor_purity_median": float(purity["purity"].median()),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--arm", default="A", choices=list(features.ARMS))
    ap.add_argument("--nodes", type=int, default=som_adapter.NOTEBOOK_SOM["num_nodes"])
    ap.add_argument("--sweep", nargs="*", type=int, default=None)
    ap.add_argument("--seed", type=int, default=som_adapter.NOTEBOOK_SOM["seed"])
    args = ap.parse_args()

    donors = args.donors or [str(d) for d in panel.cohort()["donor_id"]]
    markers = features.feature_markers(args.arm)
    print(f"arm {args.arm}: {len(markers)} markers, {len(donors)} donor(s)", flush=True)

    block, load_report = features.load_block(donors, markers, args.arm)
    print(f"  {load_report['n_cells']:,} cells "
          f"(-{load_report['dropped_qc_ineligible']:,} qc-ineligible, "
          f"-{load_report['dropped_nonfinite']:,} non-finite)", flush=True)

    out_dir = OUT / "clusters"
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = "-".join(donors) if len(donors) <= 2 else "cohort"

    reports = []
    for nodes in (args.sweep or [args.nodes]):
        clusters, profile, report = run_one(block, markers, nodes, args.arm, args.seed)
        report["donors"] = donors
        report["n_markers"] = len(markers)
        report.update(load_report)
        reports.append(report)
        print(f"  nodes={nodes:>3}  occupied={report['n_clusters_occupied']:>3}  "
              f"{report['seconds']:>6.1f}s  donor-purity max={report['donor_purity_max']:.3f} "
              f"median={report['donor_purity_median']:.3f}", flush=True)

        stem = f"{tag}_arm{args.arm}_n{nodes}"
        pd.DataFrame({"donor_id": block["donor_id"], "object_id": block["object_id"],
                      "cell_region": block["cell_region"], "islet_num": block["islet_num"],
                      "cluster": clusters}).to_parquet(
            out_dir / f"{stem}.parquet", index=False)
        profile.to_csv(out_dir / f"{stem}_profile.csv")

    manifest = {"run_id": panel.RUN_ID, "arm": args.arm,
                "markers": markers, "sweep": reports,
                **som_adapter.provenance()}
    path = out_dir / f"{tag}_arm{args.arm}_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
