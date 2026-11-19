"""Inspect, reuse and review pinned sequential Xenium–PhenoCycler analyses.

This analysis route takes an explicit recipe and a fresh output directory. It
does not use same-slide stages, default project roots, or release acceptance.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import sys

from ..runtime_context import runtime_context, path_status
from .sequential_workflow import inspect_bundle, verify_bundle, replay_assignment, digest


def _recipe(path):
    path = Path(path).expanduser().resolve()
    source = path.read_bytes()
    initial_hash = hashlib.sha256(source).hexdigest()
    recipe = json.loads(source)
    if recipe.get("schema") != "multiplex5000-sequential-recipe/1":
        raise ValueError("Unrecognized sequential recipe schema")
    for section in ("base_bundle", "current_islets"):
        item = recipe[section]
        root = Path(item["path"])
        if not root.is_absolute():
            raise ValueError(f"{section} must bind an absolute source package path")
        if digest(root / "artifact_manifest.json") != item["manifest_sha256"]:
            raise ValueError(f"Source manifest changed: {section}")
    _pins(recipe)
    if digest(path) != initial_hash:
        raise ValueError("Recipe changed during inspection")
    return path, recipe, initial_hash


def _pins(recipe):
    checked = 0
    for section in ("base_bundle", "current_islets"):
        root = Path(recipe[section]["path"]).resolve()
        for pin in recipe[section].get("pins", []):
            path = (root / pin["path"]).resolve()
            if root not in path.parents:
                raise ValueError(f"Recipe artifact escapes source package: {pin['path']}")
            if digest(path) != pin["sha256"] or path.stat().st_size != pin["bytes"]:
                raise ValueError(f"Pinned source changed: {path}")
            checked += 1
    return checked


def _fresh_output(value, recipe):
    output = Path(value).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Use a new output directory: {output}")
    roots = [Path(recipe[k]["path"]).resolve() for k in ("base_bundle", "current_islets")]
    roots += [Path(r["path"]).resolve().parent for r in recipe.get("related_reports", [])]
    specs = json.loads((roots[0] / "run_specs.json").read_text())
    roots += [Path(r["output"]).resolve() for r in specs]
    if any(output == r or r in output.parents or output in r.parents for r in roots):
        raise ValueError("Output must be separate from all source result and native construction trees")
    for ancestor in output.parents:
        if any((ancestor / name).is_file() for name in
               ("artifact_manifest.json", "run_manifest.json", "review_manifest.json", "pipeline.json")):
            raise ValueError(f"Output would alter an existing result package: {ancestor}")
    return output


def inspect_recipe(path, output=None):
    path, recipe, recipe_hash = _recipe(path)
    base = inspect_bundle(recipe["base_bundle"]["path"])
    if base["missing_artifacts"]:
        raise ValueError(f"Missing base artifacts: {base['missing_artifacts']}")
    target = _fresh_output(output, recipe) if output else None
    return {"schema": "multiplex5000-sequential-inspection/1", "status": "READY_FOR_VERIFICATION",
            "recipe": str(path), "recipe_sha256": recipe_hash, "name": recipe["name"],
            "output_dir": str(target) if target else None, "runtime": runtime_context(),
            "corrected_multicell": base, "current_islets": recipe["current_islets"],
            "pinned_recipe_files_verified": _pins(recipe), "side_effects": False,
            "native_content_verified": False, "scientific_acceptance": False,
            "planned_actions": ["Verify completed corrected multicell/native outputs and transforms",
                                "Verify current Islet-only result, native IDs and displayed evidence",
                                "Create a separate review workspace with export/resume"],
            "native_construction_planned": False, "cohort_assignment_rerun_planned": False}


def verify_recipe(path):
    from .sequential_current import verify_current_islets
    path, recipe, recipe_hash = _recipe(path)
    base = verify_bundle(recipe["base_bundle"]["path"], verify_native=True)
    current = verify_current_islets(recipe["current_islets"]["path"],
                                     expected_manifest_sha256=recipe["current_islets"]["manifest_sha256"],
                                     verify_images=True)
    native = json.loads((Path(recipe["base_bundle"]["path"]) / "native_output_reuse.json").read_text())
    native_paths = {(str(r["donor_id"]), r["role"]): (Path(r["output"]) / "objects.parquet").resolve()
                    for r in native}
    for row in current["native_sources"]:
        if Path(row["path"]).resolve() != native_paths.get((row["donor_id"], row["role"])):
            raise ValueError("Current Islet analysis uses a different native construction")
    if current["display_directory"] != recipe["current_islets"]["display_directory"]:
        raise ValueError("Recipe selects a different display generation")
    _pins(recipe)
    if digest(path) != recipe_hash:
        raise ValueError("Recipe changed during verification")
    return {"schema": "multiplex5000-sequential-verification/1", "status": "PASS",
            "recipe": str(path), "recipe_sha256": recipe_hash, "runtime": runtime_context(),
            "corrected_multicell": base, "current_islets": current, "side_effects": False,
            "raw_matrices_or_full_TIFFs_rehashed": False,
            "scientific_designations_changed": False, "production_activated": False,
            "verified_utc": datetime.now(timezone.utc).isoformat()}


def _command(*args):
    command = ["env", "PYTHONDONTWRITEBYTECODE=1", f"PYTHONPATH={Path(__file__).resolve().parents[2]}"]
    site = runtime_context().get("site_file")
    if site:
        command.append("MPLX5K_SITE_FILE=" + str(site))
    return shlex.join([*command, sys.executable, "-B", "-m", "multiplex5000.cli", "sequential", *map(str, args)])


def run_recipe(path, output):
    from .sequential_review import build_workspace
    path, recipe, recipe_hash = _recipe(path)
    target = _fresh_output(output, recipe)
    verification = verify_recipe(path)
    if verification["recipe_sha256"] != recipe_hash or digest(path) != recipe_hash:
        raise ValueError("Recipe changed before workspace creation")
    target.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "multiplex5000-sequential-workflow/1", "status": "RUNNING",
               "started_utc": datetime.now(timezone.utc).isoformat(), "recipe": str(path),
               "recipe_sha256": recipe_hash, "output_dir": str(target), "runtime": runtime_context(),
               "native_construction_executed": False, "cohort_assignments_recomputed": False,
               "scientific_designations_changed": False, "production_activated": False}
    try:
        (target / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
        (target / "recipe.json").write_text(json.dumps(recipe, indent=2) + "\n")
        receipt["review"] = build_workspace(recipe, target, verification=verification, recipe_path=path)
        _pins(recipe)
        if digest(path) != receipt["recipe_sha256"]:
            raise ValueError("Recipe changed during execution")
        receipt["commands"] = {
            "verify_sources": _command("verify", "--recipe", path),
            "verify_saved_review": _command("review-status", "--review-dir", target,
                                             "--decisions", "SAVED_PROGRESS.json"),
        }
        receipt["status"] = "COMPLETE_REUSE_AND_REVIEW_WORKSPACE"
    except BaseException as error:
        receipt["status"] = "FAILED"
        receipt["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        (target / "pipeline.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def parser():
    top = argparse.ArgumentParser(prog="mplx5k sequential", description=__doc__)
    commands = top.add_subparsers(dest="action", required=True)
    for action, help_text in (("inspect", "read-only pinned inputs and runtime inspection"),
                              ("verify", "read-only result, native reuse and current review evidence checks"),
                              ("run", "reuse completed results in a new review workspace"),
                              ("replay-multicell", "bounded donor replay of the earlier corrected multicell assignments")):
        p = commands.add_parser(action, help=help_text)
        p.add_argument("--recipe", type=Path, required=True)
        if action in ("run", "replay-multicell"):
            p.add_argument("--output-dir", type=Path, required=True)
            p.add_argument("--dry-run", action="store_true")
        if action == "replay-multicell":
            p.add_argument("--donor", required=True, help="one donor with existing PASS/WARN transform")
    review = commands.add_parser("review-status", help="validate exported reviewer notes and source freshness")
    review.add_argument("--review-dir", type=Path, required=True)
    review.add_argument("--decisions", type=Path, required=True)
    return top


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.action == "inspect":
            result = inspect_recipe(args.recipe)
        elif args.action == "verify":
            result = verify_recipe(args.recipe)
        elif args.action == "review-status":
            from .sequential_review import review_status
            result = review_status(args.review_dir, args.decisions)
        elif args.action == "replay-multicell":
            path, recipe, recipe_hash = _recipe(args.recipe)
            output = _fresh_output(args.output_dir, recipe)
            if args.dry_run:
                specs = json.loads((Path(recipe["base_bundle"]["path"]) / "run_specs.json").read_text())
                selected = [s for s in specs if str(s["donor_id"]) == args.donor]
                if sorted(s["role"] for s in selected) != ["RNA", "protein"]:
                    raise ValueError("Replay donor needs exactly one recorded RNA/protein section pair")
                result = {"status": "READY_FOR_BOUNDED_MULTICELL_REPLAY", "side_effects": False,
                          "donor": args.donor, "output_dir": str(output), "native_sections": selected,
                          "current_islet_analysis_recomputed": False, "scientific_acceptance": False}
            else:
                result = replay_assignment(recipe["base_bundle"]["path"], output, donor=args.donor)
        elif args.dry_run:
            result = inspect_recipe(args.recipe, args.output_dir)
        else:
            result = run_recipe(args.recipe, args.output_dir)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, KeyError, ImportError) as error:
        print(f"mplx5k sequential: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
