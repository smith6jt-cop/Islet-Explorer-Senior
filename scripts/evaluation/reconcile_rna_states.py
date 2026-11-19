#!/usr/bin/env python3
"""Independently count every RNA call state, including unavailable acquisitions."""
import argparse
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.resolve().is_relative_to(args.root.resolve()):
        parser.error("output must be outside the scientific release")
    rows = []
    for row in pd.read_csv(args.root / "06_typing/rna_acquisitions.csv",
                           keep_default_na=False).to_dict("records"):
        sample = row["sample_key"]
        path = args.root / f"06_typing/rna/{sample}/{sample}.calls.parquet"
        counts = dict.fromkeys(["called_cells", "broad_cells", "unresolved_cells",
                                "uncalled_cells", "unavailable_cells"], 0)
        if path.is_file():
            counts_by_label = pq.read_table(path, columns=["lineage"]).to_pandas().lineage.value_counts(dropna=False)
            for label, count in counts_by_label.items():
                if not isinstance(label, str):
                    raise TypeError(f"{sample}: non-string lineage")
                if label.startswith("Unavailable_"):
                    state = "unavailable_cells"
                elif label.startswith("Uncalled_") or label.endswith("_uncalled"):
                    state = "uncalled_cells"
                elif label.startswith("Unresolved_") or label.endswith("_unresolved"):
                    state = "unresolved_cells"
                elif label == "Non_immune":
                    state = "broad_cells"
                else:
                    state = "called_cells"
                counts[state] += int(count)
        else:
            if row["status"] != "unavailable":
                raise ValueError(f"{sample}: missing candidate calls")
            counts["unavailable_cells"] = pq.ParquetFile(
                args.root / f"04_rna_ingest/cells/{sample}.parquet").metadata.num_rows
        errors = {k: {"actual": v, "declared": row[k]} for k, v in counts.items() if v != int(row[k])}
        rows.append({"acquisition_id": row["acquisition_id"], "counts": counts,
                     "errors": errors, "status": "fail" if errors else "pass"})
        print(row["acquisition_id"], rows[-1]["status"], flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"acquisitions": rows, "failures": sum(bool(r["errors"]) for r in rows)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
