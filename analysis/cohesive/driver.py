"""One controller for the current study; never writes to source or app roots."""
from __future__ import annotations

import argparse
import fcntl
import importlib
from pathlib import Path
import sys

from .common import Context, binding, read_json, replace_json, require, run_task, utc_now


STAGES = {
    "reconcile": ("sources", "reconcile", ()),
    "qualify": ("qualification", "qualify", ("reconcile",)),
    "protocol": ("protocol", "freeze", ("qualify",)),
    "atlas": ("atlas", "build_tables", ("qualify",)),
    "h1": ("endpoints", "h1", ("protocol",)),
    "h2": ("endpoints", "h2", ("protocol",)),
    "h3": ("endpoints", "h3", ("protocol",)),
    "h4": ("endpoints", "h4", ("protocol",)),
    "h5": ("rna_endpoints", "h5", ("protocol",)),
    "h6": ("rna_endpoints", "h6", ("protocol", "h3")),
    "inference": ("inference", "infer", ("h1", "h2", "h3", "h4", "h5", "h6")),
    "crossstudy": ("crossstudy", "compare", ("protocol",)),
    "cd8": ("cd8_roc", "roc_stage", ("reconcile",)),
    "figures": ("figures", "render", ("atlas", "inference", "crossstudy", "cd8")),
    "verify": ("delivery", "verify", ("figures",)),
}


def required_stages(through: str) -> list[str]:
    needed = set()

    def visit(name):
        if name not in needed:
            for dependency in STAGES[name][2]:
                visit(dependency)
            needed.add(name)

    visit(through)
    return [name for name in STAGES if name in needed]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "status", "run"))
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--through", choices=tuple(STAGES), default="verify")
    args = parser.parse_args()
    require(args.inputs.is_file(), "An explicit existing input contract is required")
    require(args.protocol.is_file(), "An explicit existing scientific protocol is required")
    output = args.output.resolve()
    require(not args.inputs.resolve().is_relative_to(output), "Output contains the input contract")
    require(not args.protocol.resolve().is_relative_to(output), "Output contains the scientific protocol")
    require(not output.is_relative_to(Path(__file__).resolve().parents[2]),
            "Keep generated analysis products outside the source worktree")
    manifest_path = output / "analysis_manifest.json"
    stages = required_stages(args.through)

    if args.command in ("plan", "status"):
        manifest = read_json(manifest_path) if manifest_path.exists() else {"tasks": {}}
        for name in stages:
            previous = manifest["tasks"].get(name)
            state = previous["state"] if previous else "pending"
            dependencies = ", ".join(STAGES[name][2]) or "none"
            print(f"{name:12} {state:12} prerequisites={dependencies}")
        print("Status is recorded execution state, not a fresh verification of source bytes.")
        print(f"Outputs: {output}\nSame-slide processing: separately owned; not submitted by this controller.")
        return

    output.mkdir(parents=True, exist_ok=True)
    with (output / ".run.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        roots = {"inputs": binding(args.inputs), "protocol": binding(args.protocol)}
        if manifest_path.exists():
            manifest = read_json(manifest_path)
            require(manifest["roots"] == roots, "Input/protocol changed: select a new analysis output folder")
        else:
            manifest = {"schema_version": 1, "created_UTC": utc_now(), "roots": roots,
                        "production_activation": False, "tasks": {}}
        context = Context(args.inputs.resolve(), args.protocol.resolve(), output, {})
        imported = {}
        for name in stages:
            module_name, function, dependencies = STAGES[name]
            module = importlib.import_module(f".{module_name}", __package__)
            code = binding(module.__file__)
            require(imported.setdefault(module_name, code) == code,
                    f"{module_name} changed on disk during this run; start a new run")
            files = [*roots.values(), code, binding(Path(__file__).with_name("common.py"))]
            files.extend(module.input_bindings(context, function))
            files.extend(binding(context.task(dep) / "completion.json") for dep in dependencies)
            manifest["tasks"][name] = {"state": "running", "started_UTC": utc_now()}
            replace_json(manifest_path, manifest)
            completed = False
            try:
                directory, receipt = run_task(
                    output, name, files, {"operation": function},
                    lambda stage: getattr(module, function)(context, stage),
                )
                context.completed[name] = directory
                manifest["tasks"][name] = {
                    "state": "complete", "directory": str(directory),
                    "receipt": binding(directory / "completion.json"), "result": receipt["result"],
                }
                completed = True
            finally:
                if not completed:
                    error = sys.exception()
                    manifest["tasks"][name] = {"state": "failed", "failed_UTC": utc_now(),
                                               "error_type": type(error).__name__, "error": str(error)}
                replace_json(manifest_path, manifest)


if __name__ == "__main__":
    main()
