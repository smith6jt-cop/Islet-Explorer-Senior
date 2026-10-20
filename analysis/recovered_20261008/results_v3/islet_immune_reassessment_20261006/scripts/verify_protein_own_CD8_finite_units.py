#!/usr/bin/env python3
"""Reject mismatched CD8 units on finite in-block UUIDs independently of labels."""
from pathlib import Path
import argparse
import importlib.util
import json
import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("source_rows", HERE / "verify_protein_candidate_expression_sources.py")
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--n-per-donor", type=int, default=2048)
    args = ap.parse_args(); root = args.root.resolve()
    manifest = json.loads((root / "source_manifest.json").read_text())
    if manifest["arms"] != ["own_compartment"] or len(manifest["donors"]) != 25:
        raise ValueError("Finite CD8 verification requires the exact all25 own-source generation")
    shifts_path = root / "own_compartment/shifts.json"
    shifts = json.loads(shifts_path.read_text())
    rows = []
    for donor in manifest["donors"]:
        path = root / f"own_compartment/uniform/donor_id={donor}/data.parquet"
        parts = []; n = 0; outside = 0
        for batch in pq.ParquetFile(path).iter_batches(batch_size=32768, columns=["object_id", "CD8", "in_alignment_block", "qc_analysis_eligible"]):
            frame = batch.to_pandas()
            denied = ~frame.in_alignment_block.to_numpy(bool)
            if np.isfinite(frame.loc[denied, "CD8"].to_numpy(float)).any():
                raise ValueError(f"Separate-block CD8 aligned value is finite: {donor}")
            outside += int(denied.sum())
            eligible = frame.in_alignment_block.to_numpy(bool) & frame.qc_analysis_eligible.to_numpy(bool) & np.isfinite(frame.CD8.to_numpy(float))
            if n < args.n_per_donor and eligible.any():
                part = frame.loc[eligible, ["object_id", "CD8"]].iloc[:args.n_per_donor - n]
                parts.append(part); n += len(part)
        if n != args.n_per_donor:
            raise ValueError(f"Insufficient finite in-block CD8 evidence to verify units: {donor}, n={n}")
        import pandas as pd
        aligned = pd.concat(parts, ignore_index=True).set_index("object_id")
        source_path = root / f"own_compartment/expression/donor_id={donor}/data.parquet"
        source = S.exact_rows(source_path, ["CD8"], aligned.index)
        multiplier = float(shifts["shifts"]["CD8"]["scale"][donor])
        if multiplier <= 0 or not np.isfinite(multiplier):
            raise ValueError(f"Invalid matching CD8 multiplier: {donor}")
        recovered = aligned.CD8.to_numpy(float) / multiplier
        if not np.allclose(recovered, source.CD8.to_numpy(float), rtol=1e-7, atol=1e-8, equal_nan=False):
            raise ValueError(f"Inverse CD8 fails finite own-source parity: {donor}")
        rows.append(dict(donor_id=donor, sampled_finite_QC_in_block_exact_UUIDs=n,
            separate_block_rows_checked=outside, separate_block_CD8_all_unavailable=True,
            inverse_multiplier=1/multiplier, own_source_finite_unit_parity=True,
            expression_sha256=S.sha(source_path), aligned_sha256=S.sha(path)))
        print("verified finite own-source CD8 units", donor, n, "outside block", outside, flush=True)
    (root / "CD8_finite_own_source_unit_parity.json").write_text(json.dumps(dict(status="passed", rows=rows,
        shifts_sha256=S.sha(shifts_path), verifier_sha256=S.sha(Path(__file__)),
        n_exact_finite_UUIDs=sum(r["sampled_finite_QC_in_block_exact_UUIDs"] for r in rows),
        sampling="First QC-eligible in-block finite CD8 UUIDs; no label or intensity-tail selection. Every outside-block CD8 row also checked.",
        accuracy_validated=False, interpretation="Independent finite-value unit and block masking integrity; no biological validation"), indent=2)+"\n")

if __name__ == "__main__":
    main()
