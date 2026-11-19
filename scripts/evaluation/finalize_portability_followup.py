#!/usr/bin/env python3
"""Assemble immutable follow-up evidence and update the dated evaluation ledgers.

Run after the read-only audits finish, before recording observations and hashing
the final delivery. Existing frozen evidence may only be reused byte-for-byte;
changed results need a new dated artifact, not a rewrite of history references.
"""
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
DOC = BASE / "docs/evaluations/dynamic_portability_20261102"
FOLLOW = BASE / "outputs/handoff_completion_20261102"
MALDI = BASE / "outputs/maldi_completion_20261102"


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def freeze(path, payload):
    if path.exists() and path.read_bytes() != payload:
        raise ValueError(f"Frozen evidence differs; create a new version instead: {path}")
    path.write_bytes(payload)


def main():
    numerical = read(FOLLOW / "codex/numerical_parity.json")
    direct = read(FOLLOW / "codex/source_to_port_diagnostics.json")
    producer = read(FOLLOW / "codex/producer_environment_single_thread_tables.json")
    for result in (numerical, direct, producer):
        assert result["status"] == "complete" and result["images_completed"] == result["expected_images"] == 20
    proof = read(FOLLOW / "codex/journal_proof.json")
    assert proof["status"] == "producer_outputs_and_execution_inputs_verified"
    copies = {
        "maldi_source_evidence.json": MALDI / "source_evidence.json",
        "pln_object_availability.json": MALDI / "pln_object_availability.json",
        "codex_journal_proof.json": FOLLOW / "codex/journal_proof.json",
        "codex_completion_summary.json": FOLLOW / "codex/completion_summary.json",
        "codex_numerical_parity.json": FOLLOW / "codex/numerical_parity.json",
        "codex_source_to_port_diagnostics.json": FOLLOW / "codex/source_to_port_diagnostics.json",
        "codex_producer_single_thread_tables.json": FOLLOW / "codex/producer_environment_single_thread_tables.json",
        "parameter_source_manifest.json": FOLLOW / "parameter_source_manifest.json",
    }
    for destination, source in copies.items():
        freeze(DOC / destination, source.read_bytes())
    trace = read(FOLLOW / "parameter_traces.json")
    assert trace["summary"]["traced_declarations"] == trace["summary"]["catalog_declarations"] == 572
    for name in ("parameter_traces.json", "parameter_literal_candidates.json"):
        freeze(DOC / (name + ".gz"), gzip.compress((FOLLOW / name).read_bytes(), mtime=0))
    summary = {k: trace[k] for k in ("schema", "scope", "complete", "completion_interpretation", "summary", "method_limits", "remaining_work")}
    freeze(DOC / "parameter_trace_summary.json", (json.dumps(summary, indent=2) + "\n").encode())

    acquisitions = read(DOC / "acquisition_dispositions.json")
    by_image = {r["image"]: r for r in numerical["results"]}
    for row in acquisitions:
        if row["cohort"] == "codex_baseline":
            result = by_image[row["acquisition_id"]]
            row.update(cells=result["identity_cells"], cell_count_scope="eligible original identity table; exact native IDs reconciled",
                       disposition="completed",
                       reason="Exact original-source/port computations and identity assertions pass. Original mtime fixture remains rejected; levels/canary exact CSV assertions fail. Producer-thread/parser replay is separately diagnosed; v20 current-file bindings remain mismatched.",
                       evidence="codex_numerical_parity.json; codex_source_to_port_diagnostics.json; codex_journal_proof.json; codex_producer_single_thread_tables.json")
    write(DOC / "acquisition_dispositions.json", acquisitions)
    with (DOC / "acquisition_dispositions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for row in acquisitions for k in row)))
        writer.writeheader()
        writer.writerows(acquisitions)
    accounting = read(DOC / "cohort_accounting.json")
    # Keep original inventory partitions; the evaluation disposition is updated
    # separately below without inventing donor/section relationships.
    accounting["final_disposition_counts"] = {s: sum(r["disposition"] == s for r in acquisitions)
                                              for s in ("completed", "excluded", "unavailable", "blocked")}
    counts = Counter((r["cohort"], r["modality"], r["disposition"]) for r in acquisitions)
    accounting["by_cohort_modality_disposition"] = [{"cohort": c, "modality": m, "disposition": d, "count": n}
                                                   for (c, m, d), n in sorted(counts.items())]
    accounting["followup"] = "CODEX20 evaluated by exact original-source/port parity and identity reconciliation; original artifact assertion failures remain separately recorded."
    write(DOC / "cohort_accounting.json", accounting)

    with (DOC / "dependencies.csv").open(newline="") as handle:
        dependencies = list(csv.DictReader(handle))
    for row in dependencies:
        if row["id"] == "protein.v20":
            row["disposition"] = "Independent v20 audit: 22/24 current bindings mismatch. Recovered original source-to-port comparisons pass all60 exactly; original mtime fixture and40 exact CSV failures retained separately. See codex_completion_summary.json."
        if row["id"] == "protein.pln":
            row["disposition"] += "; receipt protein_cutoff_fraction is not passed to evaluator; actual candidate threshold reads REQUIRE_MIN. See parameter_traces.json.gz."
        if row["id"] == "rna.naming":
            row["disposition"] += "; receipt naming-policy metadata does not govern hardcoded shared_er_exact_v1 producer argument; actual report bindings retained."
        if row["id"].startswith("maldi."):
            row["disposition"] += "; follow-up source comparison and exact pLN annotation availability are recorded separately without acceptance in maldi_source_review.json and pln_object_availability.json."
    with (DOC / "dependencies-v2.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dependencies[0]))
        writer.writeheader()
        writer.writerows(dependencies)

    readiness = read(DOC / "readiness.json")
    for route in readiness["routes"]:
        if route["assay"] == "codex" and route["relationship_type"] == "independent_acquisition":
            route["axes"]["valid_inputs"] = "journal/output hashes recovered; historical endocrine inputs retain stat-only binding limitation"
            route["axes"]["computational_validation"] = "20 exact identity passes; 40 exact CSV assertion failures; all20 direct original-source/port calculations exact; separate producer-thread/parser replay diagnostic"
            route["evidence"] += [n for n in copies if n.startswith("codex_") and n not in route["evidence"]]
        if route["relationship_type"] == "serial_section":
            route["evidence"] += [n for n in ("maldi_source_review.json", "maldi_source_evidence.json", "pln_object_availability.json") if n not in route["evidence"]]
    write(DOC / "readiness.json", readiness)

    validation = read(DOC / "validation.json")
    validation["codex_journal_validated_numerical"] = {k: v for k, v in numerical.items() if k != "results"}
    validation["codex_exact_source_to_port"] = {k: v for k, v in direct.items() if k != "results"}
    validation["codex_producer_environment_replay"] = {k: v for k, v in producer.items() if k != "results"}
    validation["codex_identity_cells"] = sum(r["identity_cells"] for r in numerical["results"])
    validation["parameter_execution_trace"] = trace["summary"]
    validation["MALDI_source_visual_comparisons"] = 32
    validation["pLN_object_sources"] = {k: read(copies["pln_object_availability.json"])[k] for k in ("acquisitions", "objects", "area_objects")}
    history = BASE / "outputs/dynamic_portability_20261102/parameter-history"
    validation["parameter_recording"] = read(history / "recording_summary.json")
    validation["MCP_scientific_project"] = read(history / "scientific_mcp_check.json")
    write(DOC / "validation.json", validation)
    print(json.dumps({"frozen_artifacts": len(copies) + 3, "acquisitions": len(acquisitions),
                      "codex_identity_cells": validation["codex_identity_cells"], "trace": trace["summary"]}, indent=2))


if __name__ == "__main__":
    main()
