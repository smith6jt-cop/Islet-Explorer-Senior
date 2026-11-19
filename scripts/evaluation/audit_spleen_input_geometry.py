#!/usr/bin/env python3
"""Spatial input diagnostics for eight untyped RNA and four unsupported PC33 acquisitions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from validate_typing_release import spatial


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--protein-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    for folder in ("samples", "overlays"):
        (args.out / folder).mkdir()
    rows = []
    for pair in json.loads(args.pairs.read_text())["pairs"]:
        path = Path(pair["xenium_output"]) / "cells.parquet"
        cells = pq.read_table(path, columns=["cell_id", "x_centroid", "y_centroid", "nucleus_count"]).to_pandas()
        cells = cells.rename(columns={"cell_id": "cell_uid", "x_centroid": "x_um", "y_centroid": "y_um"})
        cells["admission"] = np.where(cells.nucleus_count > 0, "nucleus present", "no nucleus")
        name = pair["xenium_acquisition_id"].replace(":", "_")
        spatial(args.out, name, cells, "admission")
        rows.append({"acquisition_id": pair["xenium_acquisition_id"], "plot": name, "cells": len(cells),
                     "finite_coordinates": bool(np.isfinite(cells[["x_um", "y_um"]]).all().all()),
                     "admission_counts": cells.admission.value_counts().to_dict(), "source": str(path),
                     "scope": "source region geometry and nucleus admission only; no RNA identities inferred"})
    receipt = json.loads((args.protein_root / "06_typing/receipt.json").read_text())
    for section in receipt["settings"]["not_typed"]["spleen_pc33"]:
        path = args.protein_root / f"03_protein_measurement/data/expression/donor_id={section}/data.parquet"
        cells = pq.ParquetFile(path).read(columns=["object_id", "X_centroid", "Y_centroid", "qc_analysis_eligible"]).to_pandas()
        cells = cells.rename(columns={"object_id": "cell_uid", "X_centroid": "x_um", "Y_centroid": "y_um"})
        cells["admission"] = np.where(cells.qc_analysis_eligible, "QC eligible", "QC excluded")
        spatial(args.out, section, cells, "admission")
        with path.open("rb") as handle:
            source_hash = hashlib.file_digest(handle, "sha256").hexdigest()
        rows.append({"section_id": section, "plot": section, "cells": len(cells),
                     "unique_native_ids": bool(cells.cell_uid.is_unique),
                     "finite_coordinates": bool(np.isfinite(cells[["x_um", "y_um"]]).all().all()),
                     "admission_counts": cells.admission.value_counts().to_dict(), "source": str(path), "sha256": source_hash,
                     "scope": "source slide geometry/QC only; PC33 lineage remains unsupported"})
    (args.out / "per_acquisition.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
