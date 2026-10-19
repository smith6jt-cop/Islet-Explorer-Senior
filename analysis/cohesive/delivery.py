"""Verify the main-cohort package and write one status view generated from the task receipts.

Checks are about the build: receipts and bytes, plotted values against the tables they came from, and
deck completeness. None of them compares a result with an expected biological finding.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import qualification
from .common import Context, binding, read_json, require, validate_outputs, write_json
from .driver import STAGES

UPSTREAM = [name for name in STAGES if name != "verify"]


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "verify", "Unknown delivery operation")
    return [binding(qualification.GENERATION)]


def verify(context: Context, stage: Path) -> dict:
    checks = []

    def check(name, passed, detail=""):
        checks.append({"check": name, "passed": bool(passed), "detail": detail})
        require(passed, f"Delivery check failed: {name} {detail}")

    receipts = {}
    for name in UPSTREAM:
        directory = context.task(name)
        receipt = read_json(directory / "completion.json")
        validate_outputs(directory, receipt)
        receipts[name] = receipt
        check(f"{name} outputs match their receipt", True, str(directory))

    figures = context.task("figures")
    slides = read_json(figures / "slides.json")["slides"]
    for slide in slides:
        check(f"{slide['name']} image exists", (figures / slide["image"]).is_file())
        for path in slide["source_data"]:
            check(f"{slide['name']} source table exists", (figures / path).is_file(), path)
    pages = sorted((figures / "qa_pages").glob("page_*.png"))
    check("PDF has a title page and one page per figure", len(pages) == len(slides) + 1, f"{len(pages)} pages")

    tests = pd.read_csv(context.task("inference") / "primary_tests.csv")
    plotted = pd.read_csv(figures / "source_data" / "06_six_hypotheses_summary__primary_tests.csv")
    check("summary slide plots the inference table exactly",
          tests[["hypothesis", "effect", "p", "q_BH6"]].equals(plotted[["hypothesis", "effect", "p", "q_BH6"]]))
    check("six hypotheses in the family", sorted(tests.hypothesis) == [f"H{i}" for i in range(1, 7)])
    check("unavailable tests display no p or q",
          tests.loc[tests.status.ne("available"), ["p", "q_BH6"]].isna().all().all())
    for h in ("H1", "H2", "H3", "H4", "H5"):
        endpoints = pd.read_csv(context.task(h.lower()) / "donor_endpoints.csv", dtype={"donor_id": str})
        primary = endpoints.loc[endpoints.configuration.eq("primary")]
        row = tests.set_index("hypothesis").loc[h]
        if row.status == "available":
            usable = primary.loc[primary.status.eq("available")]
            effect = usable.loc[usable.disease_status.eq("T1D"), "donor_effect"].mean() - \
                usable.loc[usable.disease_status.eq("ND"), "donor_effect"].mean()
            check(f"{h} effect recomputes from its donor endpoints", np.isclose(effect, row.effect), f"{effect} vs {row.effect}")
        check(f"{h} unavailable donors carry a reason",
              primary.loc[primary.status.ne("available"), "reason"].astype("string").fillna("").str.len().gt(0).all())

    manifest = read_json(context.output / "analysis_manifest.json")
    frozen = read_json(context.task("protocol") / "frozen_protocol.json")
    amendments = " and ".join(frozen[key]["ruled_UTC"] for key in ("amendment", "second_amendment") if key in frozen)
    generation = read_json(Path(qualification.GENERATION))
    lines = ["# IO60 cohesive analysis — status", "",
             f"Output: `{context.output}`  ", f"Inputs: `{manifest['roots']['inputs']['path']}`  ",
             f"Pancreas protein typing: `{generation['generation_root']}` ({generation['ruling']})  ",
             f"Protocol: `{manifest['roots']['protocol']['path']}` (frozen with the {amendments} amendments)  ",
             "Production activation: none. Retrospective/exploratory.", "", "| task | result |", "|---|---|"]
    for name in UPSTREAM:
        result = receipts[name]["result"]
        summary = {k: v for k, v in result.items() if not isinstance(v, (dict, list))} if isinstance(result, dict) else result
        lines.append(f"| {name} | {summary} |")
    lines += ["", "## Six hypotheses", "", "| | effect [95% CI] | p | BH q | reading |", "|---|---|---|---|---|"]
    for row in tests.to_dict("records"):
        if row["status"] == "available":
            lines.append(f"| {row['hypothesis']} | {row['effect']:.3g} [{row['CI95_low']:.3g}, {row['CI95_high']:.3g}] | "
                         f"{row['p']:.3g} | {row['q_BH6']:.3g} | {row['interpretation']} |")
        else:
            lines.append(f"| {row['hypothesis']} | unavailable | — | — | {row['interpretation']} |")
    agreement = pd.read_csv(context.task("crossstudy") / "agreement.csv")
    lines += ["", "## Cross-study accuracy on 6551 (audit context only)", "",
              "| comparison | TVD, % points | TCF1+TOX− | TCF1+TOX+ | TCF1−TOX+ |", "|---|---|---|---|---|"]
    for row in agreement.to_dict("records"):
        lines.append(f"| {row['comparison']} | {row['total_variation_distance_percent_points']:.1f} | "
                     + " | ".join(f"{row[f'difference_{s}']:+.1f}" for s in ("TCF1+TOX-", "TCF1+TOX+", "TCF1-TOX+")) + " |")
    lines += ["", "Same-slide four-donor integration: separately owned; not part of this package."]
    (stage / "STATUS.md").write_text("\n".join(lines) + "\n")
    pd.DataFrame(checks).to_csv(stage / "delivery_checks.csv", index=False)
    write_json(stage / "package_index.json", {
        "deck_pptx": str(figures / "deck.pptx"), "deck_pdf": str(figures / "deck.pdf"),
        "figures": str(figures / "figures"), "source_data": str(figures / "source_data"),
        "primary_tests": str(context.task("inference") / "primary_tests.csv"),
        "sensitivity_tests": str(context.task("inference") / "sensitivity_tests.csv"),
        "atlas_tables": str(context.task("atlas"))})
    return {"checks": len(checks), "all_passed": all(c["passed"] for c in checks)}
