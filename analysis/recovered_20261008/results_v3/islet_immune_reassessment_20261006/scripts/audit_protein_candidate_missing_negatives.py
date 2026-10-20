#!/usr/bin/env python3
"""Count exact primary calls made while a declared negative was unavailable.

This audit retains all original calls and never assigns an alternative lineage.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
sys.path[:0] = [str(ROOT / "core/phenocycler/scripts/typing"), str(ROOT / "core/phenocycler/scripts/celltype")]
import subpanel

def read(path, columns=None):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-root", type=Path, required=True)
    ap.add_argument("--assert-no-unavailable-negatives", action="store_true")
    args = ap.parse_args(); root = args.candidate_root.resolve()
    manifest = json.loads((root / "source_manifest.json").read_text())
    rules = subpanel.immune_rules()
    negative_markers = sorted({m for rule in rules.values() for m in rule[2]})
    summaries = []; sources = []; exact = []
    for method in manifest["arms"]:
        status = json.loads((root / method / "typing/run.json").read_text())["status"]
        if status not in ("primary_complete", "complete"):
            raise RuntimeError(f"{method} primary calls are incomplete: {status}")
    for method in manifest["arms"]:
        scale_path = root / method / "scale.parquet"
        scale = read(scale_path).set_index("marker")
        for donor in manifest["donors"]:
            calls_path = root / method / f"typing/cells/donor_id={donor}/data.parquet"
            values_path = root / method / f"uniform/donor_id={donor}/data.parquet"
            calls = read(calls_path, ["object_id", "cell_class", "immune_definition", "resolved_by", "qc_analysis_eligible"])
            names = pq.ParquetFile(values_path).schema.names
            have = [m for m in negative_markers if m in names]
            values = read(values_path, ["object_id", *have])
            frame = calls.merge(values, on="object_id", how="left", validate="one_to_one")
            availability = read(root / method / f"availability/donor_id={donor}/data.parquet").set_index("marker")
            for marker in negative_markers:
                if marker not in have or marker not in scale.index:
                    frame[marker] = np.nan
            definition_calls = frame.cell_class.eq("Immune") & frame.resolved_by.eq("definition")
            for phenotype, group in frame.loc[definition_calls].groupby("immune_definition"):
                if phenotype not in rules:
                    raise ValueError(f"Unknown exact producer immune definition: {phenotype}")
                negatives = rules[phenotype][2]
                missing = ~np.isfinite(group[negatives].to_numpy(float))
                any_missing = missing.any(axis=1)
                not_acquired = [m for m in negatives if not bool(availability.loc[m, "available"])]
                no_scale = [m for m in negatives if m not in scale.index]
                summaries.append(dict(method=method, donor_id=donor, immune_definition=phenotype,
                    n_primary_definition_calls=len(group), n_with_unavailable_declared_negative=int(any_missing.sum()),
                    declared_negative_markers=";".join(negatives), unacquired_negative_markers=";".join(not_acquired),
                    missing_scale_negative_markers=";".join(no_scale),
                    implication="Declared negative evidence unavailable; original call retained as audit context; no alternative identity assigned"))
                if any_missing.any():
                    flagged = group.loc[any_missing, ["object_id", "cell_class", "immune_definition", "resolved_by"]].copy()
                    flagged["donor_id"] = donor; flagged["method"] = method
                    flagged["unavailable_negative_markers"] = [";".join(np.asarray(negatives)[row]) for row in missing[any_missing]]
                    flagged["analysis_review_status"] = "declared_negative_unavailable"
                    exact.append(flagged)
            sources.append(dict(method=method, donor_id=donor, calls_path=str(calls_path), calls_sha256=sha(calls_path),
                                values_path=str(values_path), values_sha256=sha(values_path), scale_sha256=sha(scale_path)))
            print("audited declared negatives", method, donor, int(definition_calls.sum()), flush=True)
    table = root / "resume_missing_negative_primary_definition_counts.csv"
    pd.DataFrame(summaries).to_csv(table, index=False)
    keys = root / "resume_missing_negative_primary_definition_exact_keys.parquet"
    exact_table = (pd.concat(exact, ignore_index=True) if exact else pd.DataFrame(columns=[
        "object_id", "cell_class", "immune_definition", "resolved_by", "donor_id", "method",
        "unavailable_negative_markers", "analysis_review_status"]))
    exact_table.to_parquet(keys, index=False)
    if (root / "generation_config.json").exists():
        original_sources = json.loads((root / "generation_config.json").read_text())
        entry = next(row for row in original_sources["source_scripts"] if Path(row["path"]).name == "subpanel.py")
        old_subpanel = dict(original_path=entry["path"], snapshot_path=entry["snapshot"], sha256=entry["sha256"], generation="corrected_all25")
    else:
        original_sources = json.loads((root / "diagnostic_primary_source_snapshot.json").read_text())
        old_subpanel = next(entry for entry in original_sources["source_files"]
                            if Path(entry["original_path"]).name == "subpanel.py")
    if sha(old_subpanel["snapshot_path"]) != old_subpanel["sha256"]:
        raise ValueError("Preserved executed primary producer changed")
    (root / "resume_missing_negative_primary_definition_QC.json").write_text(json.dumps(dict(
        status="complete", current_rule_reader_code_sha256=sha(Path(subpanel.__file__)),
        executed_primary_subpanel=old_subpanel,
        counts_path=str(table), counts_sha256=sha(table), exact_keys_path=str(keys), exact_keys_sha256=sha(keys),
        sources=sources, n_flagged_exact_keys=sum(len(part) for part in exact),
        original_calls_modified=False, accuracy_validated=False,
        interpretation="Missing-negative integrity gap measured on exact producer calls; no relabeling, threshold fitting or biological validation"), indent=2) + "\n")
    if args.assert_no_unavailable_negatives and len(exact_table):
        raise ValueError(f"Corrected primary retains {len(exact_table)} unavailable declared negative definition calls")

if __name__ == "__main__":
    main()
