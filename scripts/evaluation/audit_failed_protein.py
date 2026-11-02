#!/usr/bin/env python3
"""Inspect all primary tables and reproduce refinement admission, without changing calls."""
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from validate_typing_release import sha, spatial, write

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"


def main():
    stage = next((RUN / "releases/protein-evaluation-v2").glob(".06_typing.incomplete-*"))
    out = RUN / "protein_primary"
    for name in ("samples", "overlays"):
        (out / name).mkdir(parents=True, exist_ok=True)
    run = json.loads((stage / "pancreas/run.json").read_text())
    assert run["failed_stage"] == "refine"
    os.environ.update(IO60_ALIGNED_ROOT=run["aligned_root"], IO60_TYPING_ROOT=run["typing_root"],
                      IO60_GEOMETRY_ROOT=run["geometry_root"], IO60_SCALE_PATH=run["scale"],
                      IO60_SUBPANEL_LOG=run["registry_criteria"], IO60_SHIFTS_PATH=run["shifts_path"])
    sys.path.insert(0, str(BASE / "core/multiplex5000/scripts/typing"))
    import refine
    markers = refine.derive_markers()
    inventory = pd.read_csv(BASE / "docs/evaluations/dynamic_portability_20261102/inventories/io60_acquisitions.csv",
                            dtype=str, keep_default_na=False)
    rows, checks = [], []
    for donor in run["donors"]:
        acquisition = inventory[(inventory.donor_id == donor) & (inventory.platform == "phenocycler") & (inventory.tissue == "pancreas")]
        assert len(acquisition) == 1
        source = Path(run["aligned_root"]) / f"donor_id={donor}/data.parquet"
        actual = sha(source)
        checks.append({"role": "input", "path": str(source), "expected": run["aligned_files"][donor],
                       "actual": actual, "status": "pass" if actual == run["aligned_files"][donor] else "fail"})
        primary = stage / f"pancreas/cells/donor_id={donor}/data.parquet"
        cells = pq.ParquetFile(primary).read(columns=["object_id", "cell_class", "qc_analysis_eligible",
                "X_centroid", "Y_centroid", "cell_region", "refine_pass", "refined_class"]).to_pandas()
        original = pq.ParquetFile(source).read(columns=["object_id", "qc_analysis_eligible", "cell_region", "X_centroid", "Y_centroid"]).to_pandas()
        assert not cells.object_id.isna().any() and not cells.object_id.duplicated().any()
        joined = cells.merge(original, on="object_id", how="outer", validate="one_to_one", indicator=True, suffixes=("", "_source"))
        aligned_match = joined._merge.eq("both")
        shared = joined.loc[aligned_match]
        expected = original.loc[original.qc_analysis_eligible & ~original.cell_region.eq("pLN"), "object_id"]
        aligned_identity = set(cells.object_id) == set(expected)
        assert np.isfinite(cells[["X_centroid", "Y_centroid"]]).all().all()
        assert cells.refine_pass.eq(0).all() and cells.refined_class.isna().all()
        values, have = refine.donor_values(donor, markers)
        frame = cells.merge(values[["object_id"] + have], on="object_id", how="left", validate="one_to_one")
        eligible_ids = values.loc[values.qc_analysis_eligible, "object_id"]
        candidate = frame.object_id.isin(eligible_ids) & frame.qc_analysis_eligible & ~frame.cell_class.isin(refine.UNTOUCHED)
        finite = np.isfinite(frame[have].to_numpy("float64"))
        measurable = candidate.to_numpy() & finite.all(axis=1)
        missing = {m: int((~finite[:, i] & candidate.to_numpy()).sum()) for i, m in enumerate(have)}
        raw = pq.ParquetFile(source).read(columns=have).to_pandas()
        raw_missing = {m: int((~np.isfinite(raw[m])).sum()) for m in have}
        geometry_source, geometry_kind = refine.S.cell_source(donor)
        geometry = pq.ParquetFile(geometry_source).read(columns=["object_id", "qc_analysis_eligible"]).to_pandas()
        row = {"donor_id": donor, "acquisition_id": acquisition.iloc[0].acquisition_id,
               "declared_segmentation_id": acquisition.iloc[0].segmentation_version, "cells": len(original),
               "eligible": int(original.qc_analysis_eligible.sum()), "primary_cells": len(cells),
               "primary_eligible": int(cells.qc_analysis_eligible.sum()),
               "expected_primary_cells": len(expected), "primary_missing_ids": int((~expected.isin(cells.object_id)).sum()),
               "primary_extra_ids": int((~cells.object_id.isin(expected)).sum()),
               "aligned_shared_ids": int(aligned_match.sum()), "aligned_only_ids": int(joined._merge.eq("right_only").sum()),
               "primary_only_ids": int(joined._merge.eq("left_only").sum()),
               "shared_coordinate_nonidentical_rows": int((shared[["X_centroid", "Y_centroid"]].to_numpy() != shared[["X_centroid_source", "Y_centroid_source"]].to_numpy()).any(axis=1).sum()),
               "geometry_source": {"path": str(geometry_source), "sha256": sha(geometry_source), "kind": geometry_kind},
               "geometry_cells": len(geometry), "primary_ids_in_geometry": bool(cells.object_id.isin(geometry.object_id).all()),
               "identity_check": "pass" if aligned_identity else "fail", "disposition": "blocked",
               "stage_status": "failed", "review_scope": "primary-only diagnostics; no completed Stage 06 receipt",
               "class_counts": cells.cell_class.value_counts().to_dict(), "tissue_counts": cells.cell_region.value_counts().to_dict(),
               "refinement_candidates_before_finite_filter": int(candidate.sum()),
               "refinement_measurable_cells": int(measurable.sum()), "markers_carried": have,
               "markers_absent": sorted(set(markers) - set(have)), "raw_nonfinite_by_marker": raw_missing,
               "candidate_nonfinite_after_fragment_mask": missing,
               "primary_table": {"path": str(primary), "sha256": sha(primary)},
               "reason": "All 25 primary tables built; first donor 0112 has zero measurable cells and stops refinement because an unavailable Ker8_18 column is treated as carried. No refined or assembled mixed-tissue output was published."}
        rows.append(row)
        frame = cells.rename(columns={"object_id": "native_cell_id", "X_centroid": "x_um", "Y_centroid": "y_um"})
        spatial(out, donor, frame[["native_cell_id", "x_um", "y_um", "cell_class", "cell_region"]], "cell_class")
        print(donor, len(cells), "measurable", int(measurable.sum()), "nonfinite", {m: n for m, n in missing.items() if n}, flush=True)
    for name in ("failure.json", "pancreas/run.json", "pancreas/subpanel_audit.json", "pancreas/refine.log",
                 "pancreas/refinement/representation_params.parquet", "pancreas/refinement/representation_provenance.json",
                 "pancreas_scale.parquet"):
        dest = out / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(stage / name, dest)
    upstream = Path(run["aligned_root"]).parents[2]
    for source_stage in ("00_inventory", "03_protein_measurement", "05_normalization"):
        for name in ("receipt.json", "gate/accepted.json", "gate/gate.json"):
            src = upstream / source_stage / name
            if src.is_file():
                dest = out / "upstream" / source_stage / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, dest)
    pd.DataFrame(checks).to_csv(out / "receipt_checks.csv", index=False)
    write(out / "per_acquisition.json", rows)
    write(out / "summary.json", {"acquisitions": len(rows), "cells": sum(r["cells"] for r in rows),
                                "eligible": sum(r["eligible"] for r in rows), "primary_cells": sum(r["primary_cells"] for r in rows),
                                "identity_failures": [r["donor_id"] for r in rows if r["identity_check"] == "fail"],
                                "hash_failures": [r for r in checks if r["status"] != "pass"],
                                "zero_measurable_acquisitions": [r["donor_id"] for r in rows if r["refinement_measurable_cells"] == 0],
                                "completed_typing_acquisitions": 0, "status": "blocked",
                                "scope": "all primary IDs/counts/coordinates plus actual refinement admission; failed stage is not a release"})


if __name__ == "__main__":
    main()
