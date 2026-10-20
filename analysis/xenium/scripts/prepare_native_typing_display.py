#!/usr/bin/env python3
"""Read-only, bounded-memory display summaries of production protein typing.

All counts remain in the independent typing audit. Centroid samples are uniform
within donor, for display only. Marker summaries describe the markers used by
typing and are not independent accuracy estimates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
from protein_analysis_admission import DEFAULT_EXCLUSIONS, load_protein_exclusions, protein_analysis_mask

MARKERS = ["INS", "GCG", "SST", "CD99", "CD3e", "CD4", "CD8", "CD20",
           "CD79a", "CD68", "CD11c", "CD56", "CD66", "CD31", "Podoplanin", "B3TUBB"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path("/data/io60-calibration/pancreas_release_production"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--exclusions", type=Path, default=DEFAULT_EXCLUSIONS,
                   help="mandatory exact-key protein typing exclusion registry")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summaries, samples, manifest = [], [], []
    for path in sorted(args.source.glob("donor_id=*/data.parquet")):
        native = path.parent.name.split("=", 1)[1]
        donor = str(int(native))
        names = set(pq.ParquetFile(path).schema.names)
        columns = ["cell_class", "X_centroid", "Y_centroid", "disease_status"]
        columns += [c for c in ("cell_id", "object_id", "qc_analysis_eligible", "analysis_eligible",
                                "resolved_by", "endocrine_object_class", "endocrine_subtype_status")
                    if c in names]
        columns += ["typing_scaled_" + m for m in MARKERS if "typing_scaled_" + m in names]
        cells = pd.read_parquet(path, columns=columns)
        exclusions = load_protein_exclusions(args.exclusions, donor)
        cells = cells.loc[protein_analysis_mask(cells, donor, exclusions)].copy()
        cells["analysis_eligible"] = True
        n = len(cells)
        rng = np.random.default_rng(20261005 + int(donor))
        index = np.sort(rng.choice(n, size=min(n, 100000), replace=False))
        display_columns = list(dict.fromkeys([c for c in columns if not c.startswith("typing_scaled_")]
                                            + ["analysis_eligible"]))
        sample = cells.iloc[index][display_columns].copy()
        sample["donor_id"] = donor
        samples.append(sample)
        marker_columns = ["typing_scaled_" + m for m in MARKERS if "typing_scaled_" + m in names]
        for identity, part in cells.groupby("cell_class", dropna=False, observed=True):
            for column in marker_columns:
                values = pd.to_numeric(part[column], errors="coerce")
                summaries.append({"donor_id": donor, "cell_class": identity,
                                  "marker": column.removeprefix("typing_scaled_"),
                                  "denominator_scope": "typing_analysis_eligible_cells",
                                  "n_class_cells": len(part), "n_marker_measured": int(values.notna().sum()),
                                  "mean_calibrated_marker_signal": float(values.mean())})
        manifest.append({"donor_id": donor, "source": str(path), "eligible_cells": n,
                         "mandatory_typing_exclusion_keys": len(exclusions),
                         "display_cells": len(sample), "seed": 20261005 + int(donor)})
        print(f"donor {donor}: {n:,} eligible; {len(sample):,} displayed", flush=True)
    pd.concat(samples, ignore_index=True).to_parquet(args.out / "pancreas_native_centroid_display_sample.parquet", index=False)
    pd.DataFrame(summaries).to_csv(args.out / "pancreas_native_identity_marker_profiles.csv", index=False)
    (args.out / "native_display_sampling_manifest.json").write_text(json.dumps({
        "sampling": "uniform without replacement within each donor; capped at 100000 cells; display only",
        "denominator_scope": "typing_analysis_eligible_cells",
        "mandatory_exclusion_registry": str(args.exclusions.resolve()),
        "analysis_policy": "available eligibility masks AND non-Excluded primary AND exact mandatory keys excluded; ordinary endocrine subtype uncertainty retained",
        "coordinates": "native release X_centroid/Y_centroid; no inter-donor registration or common tissue scale asserted",
        "marker_profile": "within-donor class cell-average of released typing_scaled values; measured entries only; no missing-to-zero conversion; not independent validation",
        "donors": manifest,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
