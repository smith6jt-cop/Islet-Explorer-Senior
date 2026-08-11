"""Run the production RESTORE + simultaneous typing path on a few donors.

A pilot is a review artifact, not a pipeline stage: it reads an existing
content-addressed run's immutable inputs and writes only under data/pilots/.
It never writes a stage manifest and never mutates data/runs/.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from phenocycler.config import load_config  # noqa: E402
from phenocycler.marker_registry import load_registry  # noqa: E402
from phenocycler.restore_stage import (  # noqa: E402
    RestoreStageConfig,
    restore_donor_from_paths,
)
from phenocycler.simultaneous_typing import (  # noqa: E402
    load_broad_gate_policy,
    type_simultaneous_cells,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", default="3bf4fc8e391d007766d1")
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--config", default=str(REPO / "config.ini"))
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    registry = load_registry(cfg.marker_registry)
    restore_config = RestoreStageConfig.from_pipeline(cfg)
    policy = load_broad_gate_policy(registry=registry)

    run_root = REPO / "data" / "runs" / args.run
    out = REPO / "data" / "pilots" / args.label
    (out / "marker_evidence").mkdir(parents=True, exist_ok=True)
    (out / "assignments").mkdir(parents=True, exist_ok=True)
    (out / "audit").mkdir(parents=True, exist_ok=True)

    print(f"pilot            : {out}")
    print(f"source run       : {run_root}")
    print(f"divisor          : {restore_config.divisor_statistic}")
    print(f"features         : {','.join(restore_config.feature_compartments)}")
    print(f"control          : {restore_config.control_definition}")
    print(f"gates            : {policy.rules_version}")
    print(f"donors           : {' '.join(args.donors)}", flush=True)

    models = []
    for donor in args.donors:
        started = time.time()
        result = restore_donor_from_paths(
            donor_id=donor,
            cells_root=run_root / "cells",
            redsea_root=run_root / "redsea",
            geometry_root=run_root / "geometry_qc",
            registry=registry,
            config=restore_config,
            compartment_sample_root=cfg.restore_compartment_sample,
            expression_reviews_path=cfg.restore_expression_reviews,
        )
        evidence = result.evidence_wide
        assignments = type_simultaneous_cells(
            evidence, marker_registry=registry, policy=policy
        )
        evidence.to_parquet(
            out / "marker_evidence" / f"donor_id={donor}.parquet", index=False
        )
        assignments.to_parquet(
            out / "assignments" / f"donor_id={donor}.parquet", index=False
        )
        models.append(result.model_table)
        counts = assignments["broad_label"].value_counts()
        elapsed = time.time() - started
        print(
            f"[{donor}] {len(evidence):,} cells in {elapsed/60:.1f} min :: "
            + ", ".join(f"{k}={v:,}" for k, v in counts.items()),
            flush=True,
        )

    model_table = pd.concat(models, ignore_index=True)
    model_table.to_parquet(out / "audit" / "restore_pair_models.parquet", index=False)
    (out / "pilot.json").write_text(
        json.dumps(
            {
                "source_run": args.run,
                "donors": list(args.donors),
                "divisor_statistic": restore_config.divisor_statistic,
                "feature_compartments": list(restore_config.feature_compartments),
                "control_definition": restore_config.control_definition,
                "broad_rules_version": policy.rules_version,
            },
            indent=2,
        )
    )
    print("\npair states:")
    print(model_table["state"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
