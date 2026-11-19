"""Run the approved native-region screen and its independent validators.

This is a dataset-analysis command. It neither accepts a typing gate nor activates
a release. Input contracts stay external and immutable; outputs use a new folder.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

from ..runtime_context import path_status, runtime_context


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mplx5k gammadelta", description=__doc__)
    p.add_argument("--inputs-root", type=Path, required=True,
                   help="folder containing approved_runtime_config.json, native_sources.json, "
                        "reviewed_criteria.json and historical membership tables")
    p.add_argument("--region", required=True,
                   help="unique inventory row key (legacy run_id is a slide ID in the 25-region study)")
    p.add_argument("--output-dir", type=Path, required=True, help="new, nonexistent output folder")
    p.add_argument("--stage", choices=("all", "selection", "coordinates"), default="all")
    p.add_argument("--selection-dir", type=Path,
                   help="previous independently validated selection, for --stage coordinates only")
    p.add_argument("--dry-run", action="store_true",
                   help="inspect contract, source availability and planned commands without writing outputs")
    p.add_argument("--chunk-cells", type=int, default=20000)
    p.add_argument("--batch-size", type=int, default=262144)
    return p


def _prepare(args: argparse.Namespace):
    if args.chunk_cells < 1 or args.batch_size < 1:
        raise ValueError("Chunk and batch sizes must be positive")
    if bool(args.selection_dir) != (args.stage == "coordinates"):
        raise ValueError("--selection-dir is required only for --stage coordinates")
    inputs = args.inputs_root.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing existing output directory: {output}")
    if inputs == output or inputs in output.parents or output in inputs.parents:
        raise ValueError("Inputs and output must be separate trees")
    # Protect the immutable historical package even if a user chooses its sibling
    # folders as a destination. The low-level producers also refuse overwrites.
    for ancestor in output.parents:
        if (ancestor / "run_manifest.json").is_file():
            raise ValueError(f"Output would change an immutable run package: {ancestor}")
    inventory = inputs / "native_sources.json"
    rows = json.loads(inventory.read_text())["runs"]
    matches = [row for row in rows if str(row["run_id"]) == args.region]
    if len(matches) != 1:
        raise ValueError("The region key must identify exactly one native bundle")
    row = matches[0]
    if Path(args.region).name != args.region or args.region in (".", ".."):
        raise ValueError("Region key cannot contain path components")
    config = inputs / "approved_runtime_config.json"
    criteria = inputs / "reviewed_criteria.json"
    historical_and = inputs / "historical_and" / (args.region + ".csv")
    historical_or = inputs / "historical_or2_cells.csv"
    files = [inventory, config, criteria, historical_and, historical_or]
    bindings = {str(path): sha(path) for path in files}
    package = Path(__file__).resolve().parent
    code_hashes = {p.name: sha(p) for p in sorted(package.iterdir())
                   if p.is_file() and p.suffix in (".py", ".json")}
    selection = args.selection_dir.resolve() if args.selection_dir else output / "selector"
    shared = ["--config", str(config), "--inventory", str(inventory)]
    selector_args = [*shared, "--criteria", str(criteria), "--region", args.region,
                     "--historical-and-cells", str(historical_and),
                     "--historical-or2-cells", str(historical_or),
                     "--chunk-cells", str(args.chunk_cells)]
    coordinates_args = [*shared, "--selection-dir", str(selection), "--region", args.region,
                        "--batch-size", str(args.batch_size)]
    steps = []
    if args.stage in ("all", "selection"):
        steps.extend([
            ("select_corrected", [*selector_args, "--out-dir", str(selection)]),
            ("validate_selection", [*selector_args, "--selection-dir", str(selection),
                                     "--out-dir", str(selection / "validation")]),
        ])
    if args.stage in ("all", "coordinates"):
        steps.extend([
            ("extract_corrected", [*coordinates_args, "--output-dir", str(output / "transcripts")]),
            ("validate_coordinates", [*coordinates_args, "--extraction-dir", str(output / "transcripts"),
                                       "--output-dir", str(output / "coordinate_validation")]),
        ])
    receipt = {
        "schema": "multiplex5000-gammadelta-pipeline/1", "status": "RUNNING",
        "started_utc": datetime.now(timezone.utc).isoformat(), "stage": args.stage,
        "inventory_row_key": args.region, "native_unit": row.get("native_unit", {}),
        "native_bundle": row["bundle"], "legacy_id_semantics": row.get("legacy_id_semantics"),
        "inputs_sha256": bindings, "producer_sha256": code_hashes,
        "python": sys.version, "executable": sys.executable, "steps": [],
        "biological_identity_accepted": False, "release_activated": False,
        "nucleus_filter_applied": False, "tissue_filter_applied": False,
    }

    return receipt, steps, output, selection, bindings, code_hashes


def _checks(args, receipt, selection):
    """Check small contracts and source metadata; full native hashes remain execution checks."""
    from .select_corrected import load_contract

    inputs = args.inputs_root.resolve()
    config, _, _ = load_contract(inputs / "approved_runtime_config.json",
                                 inputs / "reviewed_criteria.json", inputs / "native_sources.json")
    expected = {inputs / "historical_and" / (args.region + ".csv"):
                config["historical_and_sha256_by_region"][args.region],
                inputs / "historical_or2_cells.csv": config["historical_or2_sha256"]}
    for path, digest in expected.items():
        if receipt["inputs_sha256"][str(path)] != digest:
            raise ValueError(f"Historical membership hash mismatch: {path}")
    rows = json.loads((inputs / "native_sources.json").read_text())["runs"]
    row = next(row for row in rows if str(row["run_id"]) == args.region)
    sources = row["sources"]
    metadata = [name for name in sources if name in ("cells.csv", "cells.csv.gz", "cells.parquet")]
    if len(metadata) != 1:
        raise ValueError("Inventory must pin exactly one native cell metadata file")
    required = ["cell_feature_matrix.h5", metadata[0]]
    if args.stage != "selection":
        required.append("transcripts.parquet")
    checks, blockers = {}, []
    for name in required:
        if name not in sources:
            checks[name] = {"status": "missing_source_pin"}
            blockers.append(f"Missing source pin: {name}")
            continue
        pin = sources[name]
        path = Path(pin["path"])
        if not path.is_absolute():
            raise ValueError(f"Source must name an absolute pinned path: {path}")
        if name != "transcripts.parquet" and path.parent != Path(row["bundle"]):
            raise ValueError(f"H5 and cell metadata must belong to the pinned bundle: {path}")
        check = path_status(path)
        check["expected_sha256"] = pin["sha256"]
        check["full_sha256_verified"] = False
        if len(pin["sha256"]) != 64 or any(c not in "0123456789abcdef" for c in pin["sha256"]):
            raise ValueError(f"Invalid full source hash: {path}")
        if check["status"] == "available":
            check["matches_pinned_stat"] = (check["regular_file"]
                and check["bytes"] == pin["bytes"] and check["mtime_ns"] == pin["mtime_ns"])
            if not check["matches_pinned_stat"]:
                blockers.append(f"Source metadata changed: {path}")
        else:
            blockers.append(f"Source {check['status']}: {path}")
        checks[name] = check
    if args.stage == "coordinates":
        # Reuse the existing artifact-binding checks; do not reimplement scientific selection.
        from .extract_corrected import load_contract as coordinate_contract
        if (selection == args.output_dir.resolve() or selection in args.output_dir.resolve().parents
                or args.output_dir.resolve() in selection.parents):
            raise ValueError("Selection and coordinate output must be separate trees")
        coordinate_contract(argparse.Namespace(config=inputs / "approved_runtime_config.json",
                            inventory=inputs / "native_sources.json", selection_dir=selection,
                            selection_validation=None, run_id=args.region))
    return checks, blockers


def inspect(args: argparse.Namespace) -> dict:
    receipt, steps, output, selection, _, _ = _prepare(args)
    checks, blockers = _checks(args, receipt, selection)
    return {"schema": "multiplex5000-gammadelta-preflight/1",
            "status": "BLOCKED" if blockers else "READY_FOR_EXECUTION",
            "side_effects": False, "runtime": runtime_context(),
            "inventory_row_key": args.region, "native_unit": receipt["native_unit"],
            "legacy_id_semantics": receipt["legacy_id_semantics"],
            "native_bundle": receipt["native_bundle"], "stage": args.stage,
            "inputs_root": str(args.inputs_root.resolve()), "output_dir": str(output),
            "selection_dir": str(selection), "receipt": str(output / "pipeline.json"),
            "inputs_sha256": receipt["inputs_sha256"],
            "producer_sha256": receipt["producer_sha256"],
            "required_sources": checks, "blockers": blockers,
            "commands": [[sys.executable, "-m", "multiplex5000.gammadelta." + module, *arguments]
                         for module, arguments in steps],
            "native_content_verified": False, "scientific_acceptance": False,
            "scope": "Small contract hashes and source path/size/mtime only. READY means these "
                     "checks passed, not full source hash/schema validation or scientific acceptance. "
                     "The producers and independent validators check native content during execution."}


def run(args: argparse.Namespace) -> dict:
    receipt, steps, output, selection, bindings, code_hashes = _prepare(args)
    _, blockers = _checks(args, receipt, selection)
    if blockers:
        raise ValueError("; ".join(blockers))
    package = Path(__file__).resolve().parent
    output.mkdir(parents=True)
    receipt_path = output / "pipeline.json"

    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")

    save()
    try:
        for module, arguments in steps:
            command = [sys.executable, "-m", "multiplex5000.gammadelta." + module, *arguments]
            started = time.monotonic()
            result = subprocess.run(command, capture_output=True, text=True)
            (output / (module + ".stdout.txt")).write_text(result.stdout)
            (output / (module + ".stderr.txt")).write_text(result.stderr)
            receipt["steps"].append({"command": command, "returncode": result.returncode,
                                      "elapsed_seconds": time.monotonic() - started})
            save()
            if result.returncode:
                raise RuntimeError(f"{module} failed; see {output / (module + '.stderr.txt')}")
        for path, digest in bindings.items():
            if sha(Path(path)) != digest:
                raise RuntimeError(f"Input contract changed: {path}")
        for name, digest in code_hashes.items():
            if sha(package / name) != digest:
                raise RuntimeError(f"Producer changed: {name}")
        validations = [selection / "validation/selection_validation.json"]
        if args.stage != "selection":
            validations.append(output / "coordinate_validation/coordinate_validation.json")
        receipt["validation_receipts"] = {}
        for path in validations:
            check = json.loads(path.read_text())
            if check.get("status") != "PASS" or check.get("all_checks_pass") is not True:
                raise RuntimeError(f"Independent validation did not pass: {path}")
            receipt["validation_receipts"][str(path)] = sha(path)
        receipt["status"] = "COMPLETE_TECHNICAL_VALIDATION"
    except BaseException as error:
        receipt["status"] = "FAILED"
        receipt["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()
    return receipt


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.dry_run:
            plan = inspect(args)
            print(json.dumps(plan, indent=2))
            return 2 if plan["blockers"] else 0
        receipt = run(args)
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        print(f"mplx5k gammadelta: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"status": receipt["status"], "inventory_row_key": args.region,
                      "receipt": str(args.output_dir.resolve() / "pipeline.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
