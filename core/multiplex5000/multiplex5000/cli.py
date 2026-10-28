"""``mplx5k`` -- run a stage, read and accept its gate, see where a release stands, run the graph.

    mplx5k stage run 00_inventory               # one stage, through the task contract
    mplx5k status                               # every stage: complete, gate, accepted, failed runs
    mplx5k gate show 00_inventory               # its checks, and where its figures are
    mplx5k gate accept 00_inventory --by NAME   # refused while a blocking check fails
    mplx5k graph --dry-run                      # the Snakemake route (needs the workflow extra)

Every command takes ``--project`` (an id under ``core/multiplex5000/projects`` or a path; default
``io60``) and ``--releases-root`` (default: the site's). The notebooks call the same code.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__, gates, roots, stages, tasks
from . import project as projects
from .pln_sources import PLNSourceError, load_exports

WORKFLOW = Path(__file__).resolve().parent / "workflow"


def _project(args) -> projects.Project:
    root = Path(args.releases_root) if args.releases_root else None
    return projects.load(args.project, releases_root=root)


def _stage_run(args) -> int:
    project = _project(args)
    result = tasks.run_stage(stages.get(args.stage), project, force=args.force)
    if result.status == "complete":
        print(f"{result.stage}: already complete at {result.path} (use --force to run it again)")
        return 0
    verdict = "passes" if result.gate_passed else "FAILS a blocking check"
    print(f"{result.stage}: published {result.path}; its gate {verdict}")
    if result.superseded:
        print(f"  the previous output was moved aside to {result.superseded}")
    if result.accepted_automatically:
        print("  accepted automatically (the project file marks this gate automatic)")
    else:
        print(f"  read {result.path / 'gate' / 'gate.html'}, then: mplx5k gate accept "
              f"{result.stage} --by <name>")
    return 0


def _status(args) -> int:
    project = _project(args)
    print(f"project {project.id}: release {project.release} at {project.release_dir}")
    for stage in project.stages:
        folder = project.release_dir / stage
        if not tasks.is_complete(folder):
            print(f"  {stage:<26} not run")
            continue
        gate = gates.load(folder)
        try:
            acceptance = gates.accepted(folder)
        except gates.GateError as err:
            acceptance, note = None, f" ({err})"
        else:
            note = ""
        state = ("accepted by " + acceptance["accepted_by"] + " " + acceptance["accepted_at"]
                 if acceptance else "not accepted" + note)
        print(f"  {stage:<26} complete; gate {'passes' if gate['passed'] else 'FAILS'}; {state}")
    for path in tasks.failed_attempts(project.release_dir):
        print(f"  failed or interrupted run left at {path}")
    if "06_typing" in project.stages:
        print("  research comparisons: release gates do not block separate analysis of available "
              "source data")
        typing = project.release_dir / "06_typing"
        if not tasks.is_complete(typing):
            print("  class-specific comparisons: need Stage 06 output; available source "
                  "measurements can be compared")
        else:
            try:
                accepted = gates.accepted(typing)
            except gates.GateError:
                accepted = None
            if accepted:
                print("  class-specific comparisons: Stage 06 accepted for this release")
            else:
                print("  class-specific comparisons: Stage 06 candidate calls can be analyzed; "
                      "label them provisional")
        print("  release-backed claims: require processing acceptance; accepting Stage 06 does "
              "not run comparative analysis")
    return 0


def _gate_show(args) -> int:
    project = _project(args)
    folder = project.release_dir / args.stage
    record = gates.load(folder)
    print(f"{args.stage} ({project.release}): {record['summary']}")
    for c in record["checks"]:
        mark = ("ok  " if c["passed"] else "FAIL") if c["blocking"] else "note"
        print(f"  {mark} {c['name']}: {c['value']}" + (f"  [{c['detail']}]" if c["detail"] else ""))
    print(f"  figures: {', '.join(str(folder / 'gate' / f) for f in record['figures']) or 'none'}")
    print(f"  page: {folder / 'gate' / 'gate.html'}")
    try:
        acceptance = gates.accepted(folder)
    except gates.GateError as err:
        print(f"  {err}")
    else:
        print(f"  accepted by {acceptance['accepted_by']} at {acceptance['accepted_at']}"
              if acceptance else "  not accepted")
    return 0


def _gate_accept(args) -> int:
    project = _project(args)
    path = gates.accept(project.release_dir / args.stage, by=args.by, note=args.note or "")
    print(f"{args.stage}: accepted by {args.by}; recorded in {path}")
    return 0


def _pln_sources(args) -> int:
    root = Path(args.root) if args.root else roots.source("pln_v3_export_cells")
    if root is None:
        raise PLNSourceError("name sources.pln_v3_export_cells in the site file or pass --root")
    roster = roots.INTEGRATION / "pln_node_annotations.csv"
    index = roots.REPO / "core/multiplex5000/sites/pln_qupath_exports_v3.csv"
    exports = load_exports(roster, index, source_root=root,
                           selected_images=set(args.image) if args.image else None,
                           require_sources=True)
    print(json.dumps({"images_checked": len(exports),
                      "export_cells": sum(item.export_cells for item in exports),
                      "older_roster_counts_differ": sum(
                          item.roster_cells != item.export_cells for item in exports),
                      "segmentation_ids": sorted({item.segmentation_id for item in exports}),
                      "source_root": str(root)}, indent=2))
    return 0


def graph_spec(project: projects.Project) -> dict:
    """What the Snakefile needs, so Snakemake never imports this package: each stage's declared
    inputs, the previous stage's acceptance, and its receipt."""
    out = []
    for stage in project.stages:
        previous = project.previous(stage)
        # A later stage may discover its raw inputs from an accepted earlier
        # stage's inventory. Defer that discovery until the human gate exists;
        # on the next graph invocation its full declared inputs are included.
        ready = previous is None or gates.accepted(project.release_dir / previous) is not None
        inputs = ({i.name: str(i.path) for i in stages.get(stage).inputs(project)
                   if i.required or Path(i.path).exists()} if ready else {})
        if previous is not None:
            inputs["previous_gate_accepted"] = str(project.release_dir / previous / gates.GATE_DIR
                                                   / gates.ACCEPTED_FILE)
        out.append({"id": stage, "inputs": inputs,
                    "receipt": str(project.release_dir / stage / tasks.RECEIPT)})
    # A graph invocation builds the next incomplete stage. Targeting every
    # receipt at once would make Snakemake reject a missing future human gate
    # before it could run the first stage.
    target = next((item["receipt"] for item in out
                   if not tasks.is_complete(project.release_dir / item["id"])),
                  out[-1]["receipt"])
    return {"project": str(project.path), "release": project.release,
            "releases_root": str(project.releases_root),
            "release_dir": str(project.release_dir), "target_receipt": target,
            "stages": out}


def _graph_spec(args) -> int:
    print(json.dumps(graph_spec(_project(args)), indent=2))
    return 0


def _snakemake(named) -> str | None:
    """The Snakemake to run: --snakemake, else the site file's, else the one on PATH."""
    if named:
        return named
    site = roots.SITE.get("snakemake") if roots.SITE else None
    return str(site) if site else shutil.which("snakemake")


def _graph(args, extra: list[str]) -> int:
    project = _project(args)
    snakemake = _snakemake(args.snakemake)
    if snakemake is None:
        print("the graph route needs Snakemake in its own environment (requirements-graph.txt); "
              "name it with --snakemake or the site file's \"snakemake\"", file=sys.stderr)
        return 2
    project.release_dir.mkdir(parents=True, exist_ok=True)
    command = [snakemake, "--snakefile", str(WORKFLOW / "Snakefile"),
               "--directory", str(project.release_dir),
               "--workflow-profile", str(WORKFLOW / "profiles" / args.profile),
               "--config", f"project={project.path}", f"releases_root={project.releases_root}",
               f"python={sys.executable}"]
    if args.dry_run:
        command.append("--dry-run")
    command += extra
    print(" ".join(command), flush=True)
    return subprocess.call(command)


def parser() -> argparse.ArgumentParser:
    top = argparse.ArgumentParser(prog="mplx5k", description=__doc__.split("\n\n")[0])
    top.add_argument("--version", action="version", version=f"multiplex5000 {__version__}")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--project", default="io60", help="project id or file (default io60)")
    common.add_argument("--releases-root", help="where releases are built (default: the site's)")
    sub = top.add_subparsers(dest="command", required=True)

    stage = sub.add_parser("stage", help="run one stage").add_subparsers(dest="action",
                                                                         required=True)
    run = stage.add_parser("run", parents=[common], help="run a stage through the task contract")
    run.add_argument("stage")
    run.add_argument("--force", action="store_true",
                     help="run again and replace the published output (moved aside, not deleted)")
    run.set_defaults(handler=_stage_run)

    status = sub.add_parser("status", parents=[common], help="where the release stands")
    status.set_defaults(handler=_status)

    gate = sub.add_parser("gate", help="read or accept a gate").add_subparsers(dest="action",
                                                                               required=True)
    show = gate.add_parser("show", parents=[common], help="a stage's checks and figures")
    show.add_argument("stage")
    show.set_defaults(handler=_gate_show)
    accept = gate.add_parser("accept", parents=[common], help="accept a stage's gate")
    accept.add_argument("stage")
    accept.add_argument("--by", required=True, help="who accepts it")
    accept.add_argument("--note", help="anything the acceptance should say")
    accept.set_defaults(handler=_gate_accept)

    sources = sub.add_parser("sources", help="inspect recorded source bindings")
    source_sub = sources.add_subparsers(dest="source_kind", required=True)
    pln = source_sub.add_parser("pln", help="check dedicated pLN v3 QuPath exports")
    pln.add_argument("--root", help="Parquet folder (default: site sources.pln_v3_export_cells)")
    pln.add_argument("--image", action="append", help="one roster image to check (repeatable)")
    pln.set_defaults(handler=_pln_sources)

    graph = sub.add_parser("graph", parents=[common],
                           help="run the project's stages through Snakemake (extra arguments "
                                "after -- go to Snakemake)")
    graph.add_argument("--profile", default="local", choices=["local", "hipergator"])
    graph.add_argument("--dry-run", action="store_true")
    graph.add_argument("--snakemake", help="the snakemake executable (default: the site's, else "
                                           "PATH)")
    graph.set_defaults(handler=None)

    spec = sub.add_parser("graph-spec", parents=[common],
                          help="the stages as the Snakefile reads them (JSON)")
    spec.set_defaults(handler=_graph_spec)
    return top


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    extra: list[str] = []
    if "--" in argv:
        cut = argv.index("--")
        argv, extra = argv[:cut], argv[cut + 1:]
    args = parser().parse_args(argv)
    try:
        if args.command == "graph":
            return _graph(args, extra)
        return args.handler(args)
    except (gates.GateError, tasks.TaskError, projects.ProjectError, stages.StageError,
            PLNSourceError) as err:
        print(f"mplx5k: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
