#!/usr/bin/env python3
"""Validate declared hash algorithms and add full-byte observations without rewriting receipts."""
import argparse
import hashlib
import json
import platform
import shutil
import sys
from pathlib import Path

from multiplex5000.artifacts import sha256_sampled_file


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="Receipt-bound QuPath manifest declaring full or sampled hash modes")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = json.loads(args.receipt.read_text())
    manifest_bindings = [r for r in receipt["inputs"] if r["name"] == "qupath_manifest"]
    if len(manifest_bindings) != 1 or sha(args.manifest) != manifest_bindings[0].get("sha256"):
        raise ValueError("hash-mode manifest is not uniquely bound by receipt")
    fingerprints = {}
    for image in json.loads(args.manifest.read_text())["images"]:
        for item in image.values():
            if isinstance(item, dict) and {"path", "hash_mode", "sha256"} <= item.keys():
                old = fingerprints.setdefault(item["path"], item)
                if old != item:
                    raise ValueError("conflicting fingerprints for one path")
    inputs = [row for row in receipt["inputs"] if not row.get("sha256")]
    if len({row["name"] for row in inputs}) != len(inputs):
        raise ValueError("duplicate receipt input names")
    report = {"schema": "recorded-input-byte-verification/1", "complete": False,
              "receipt": str(args.receipt), "receipt_sha256": sha(args.receipt),
              "producer_sha256": sha(Path(__file__)), "python": sys.version,
              "manifest_sha256": sha(args.manifest),
              "scope": "full current bytes plus comparison using declared historical hash mode; sampled history does not establish historical full bytes",
              "platform": platform.platform(), "inputs": []}
    for row in inputs:
        result = {"name": row["name"], "path": row["path"], "recorded_sha256": row.get("recorded_sha256")}
        try:
            fingerprint = fingerprints[row["path"]]
            if fingerprint["sha256"] != row.get("recorded_sha256"):
                raise ValueError("receipt and manifest recorded hashes differ")
            mode = fingerprint["hash_mode"]
            if mode not in {"full", "sampled"}:
                raise ValueError("undeclared hash algorithm")
            path = Path(row["path"])
            before = path.stat()
            digest = sha(path)
            comparison = (digest if mode == "full" else sha256_sampled_file(
                path, sample_size_bytes=fingerprint["sample_size_bytes"], sample_count=fingerprint["sample_count"]))
            after = path.stat()
            stable = (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
            result.update(sha256=digest, hash_mode="full", recorded_hash_mode=mode,
                          observed_recorded_mode_hash=comparison, bytes=after.st_size, mtime_ns=after.st_mtime_ns,
                          stable_during_hash=stable, size_matches_receipt=after.st_size == row["bytes"],
                          status="pass" if stable and comparison == row.get("recorded_sha256") and after.st_size == row["bytes"] else "fail")
        except (OSError, KeyError, ValueError) as exc:
            result.update(status="unavailable", reason=str(exc))
        report["inputs"].append(result)
        (args.out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
        print(result["name"], result["status"], flush=True)
    report.update(complete=True, passed=sum(x["status"] == "pass" for x in report["inputs"]),
                  declared=len(inputs), bytes_hashed=sum(x.get("bytes", 0) for x in report["inputs"]))
    report["status"] = "pass" if report["passed"] == len(inputs) else "fail"
    shutil.copyfile(__file__, args.out / "producer_verify_recorded_inputs.py")
    (args.out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "passed", "declared", "bytes_hashed")}))
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
