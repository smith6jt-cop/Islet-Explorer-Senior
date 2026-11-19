#!/usr/bin/env python3
"""Assemble the dated evaluation ledger from frozen inventories and audit outputs.

This is a review artifact builder, not a scientific producer or an approval route.
Identifiers come from explicit source rows; donor equality never creates a pair.
"""
import json
import shutil
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"
DOC = BASE / "docs/evaluations/dynamic_portability_20261102"


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main():
    sources = DOC / "inventories"
    sources.mkdir(parents=True, exist_ok=True)
    for name, path in {
        "io60_acquisitions.csv": RUN / "remote/rna/acquisitions.csv",
        "io60_relationships.csv": RUN / "remote/rna/relationships.csv",
        "maldi_sections.csv": RUN / "maldi/sections.csv",
        "spleen_cohort.json": RUN / "remote/cohort-spec-spleen-primary.json",
        "spleen_pair_sources.json": RUN / "remote/spleen_pair_sources.json",
    }.items():
        shutil.copyfile(path, sources / name)
    rna = {r["acquisition_id"]: r for r in read(RUN / "remote/rna/per_acquisition.json")}
    pln = {r["acquisition_id"].split(":", 1)[1]: r for r in read(RUN / "pln/per_acquisition.json")}
    protein_path = RUN / "protein/per_acquisition.json"
    protein = {r["acquisition_id"]: r for r in read(protein_path)} if protein_path.exists() else {}
    primary_path = RUN / "protein_primary/per_acquisition.json"
    primary = {r["acquisition_id"]: r for r in read(primary_path)} if primary_path.exists() else {}
    ledger = []
    inventory = pd.read_csv(sources / "io60_acquisitions.csv", keep_default_na=False, dtype=str)
    for row in inventory.to_dict("records"):
        item = {"cohort": "io60", "acquisition_id": row["acquisition_id"],
                "donor_id": row["donor_id"], "physical_section_id": None,
                "section_identity_evidence": "inventory slide/region; no inferred cross-assay identity",
                "slide_id": row["slide_id"] or row["xenium_slide"] or None,
                "roi": row["roi"], "modality": row["platform"], "tissue": row["tissue"],
                "panel": row["panel_id"] or None, "segmentation_id": row["segmentation_version"] or None,
                "source": "inventories/io60_acquisitions.csv", "cells": None, "spots": None,
                "disposition": "blocked", "evidence": None, "reason": None}
        if row["platform"] == "xenium":
            result = rna[row["acquisition_id"]]
            item.update(disposition="unavailable" if result["status"] == "unavailable" else "completed",
                        cells=int(result["cells"]), reason=result["reason"],
                        evidence="remote/rna/per_acquisition.json")
        elif row["platform"] == "phenocycler":
            if row["acquisition_id"] in protein:
                result = protein[row["acquisition_id"]]
                item.update(disposition="completed" if result["identity_check"] == "pass" else "blocked",
                            cells=result.get("cells"), reason=result.get("reason", "Legacy mixed-tissue candidate; not per-image port parity"),
                            evidence="protein/per_acquisition.json")
            elif row["acquisition_id"] in primary:
                result = primary[row["acquisition_id"]]
                item.update(disposition="blocked", cells=result["cells"], primary_cells=result["primary_cells"],
                            reason=result["reason"], completed_substeps=["primary-table production", "primary/source compatibility audit"],
                            evidence="protein_primary/per_acquisition.json")
            elif row["qptiff"] in pln:
                result = pln[row["qptiff"]]
                item.update(disposition="completed", cells=int(result["cells"]),
                            reason="Two-section candidate audit; original manifest recovered at producer pin; no final class acceptance",
                            evidence="pln/per_acquisition.json")
            else:
                item.update(reason="Full canonical dedicated-pLN route unavailable" if row["tissue"] != "pancreas"
                            else "Isolated existing 25-section evaluation requires final receipt",
                            evidence="protein-evaluation.log" if row["tissue"] == "pancreas" else "dependencies.csv#protein.pln")
        else:
            item.update(disposition="excluded", reason="post-Xenium IF is outside the requested typing routes",
                        evidence="inventories/io60_acquisitions.csv")
        ledger.append(item)
    for row in pd.read_csv(sources / "maldi_sections.csv", keep_default_na=False).to_dict("records"):
        ledger.append({"cohort": "io60", "acquisition_id": row["acquisition_id"], "donor_id": str(row["donor_id"]),
                       "physical_section_id": row["tissue_id"], "modality": "maldi", "tissue": None,
                       "panel": "38 N-glycans + 10 glycogen", "segmentation_id": None,
                       "source": "inventories/maldi_sections.csv", "cells": None, "spots": row["spots"],
                       "disposition": "blocked", "completed_substeps": ["source pairing", "spot accounting", "mask and denominator audit"],
                       "reason": row["registration"], "evidence": "maldi/summary.json"})
    codex = {r["image_name"]: r for r in read(sources / "codex_project_manifest.json")["entries"]}
    parity = read(RUN / "codex_parity/summary.json")
    assert set(codex) == {r["image"] for r in parity}, "parity roster must equal the source project inventory"
    for row in parity:
        ledger.append({"cohort": "codex_baseline", "acquisition_id": row["image"], "donor_id": codex[row["image"]]["case_id"],
                       "physical_section_id": None, "modality": "codex", "tissue": "pancreas", "panel": "codex-uniform32",
                       "segmentation_id": None, "cells": None, "spots": None,
                       "source": "inventories/codex_project_manifest.json (exact roster reconciled to parity tables)", "disposition": "blocked",
                       "reason": "Original parity fixture rejects table producer timestamps; numerical comparison did not execute",
                       "evidence": f"codex_parity/{row['image']}.log"})
    spleen = read(sources / "spleen_cohort.json")
    spleen_result = {r["section_id"]: r for r in read(RUN / "remote/spleen/per_acquisition.json")}
    pc33_inputs = {r["section_id"]: r for r in read(RUN / "remote/spleen_input_geometry/per_acquisition.json") if "section_id" in r}
    for row in spleen["images"]:
        result = spleen_result[row["donor_id"]]
        ledger.append({"cohort": "spleen", "acquisition_id": row["image_id"], "donor_id": row["donor_id"],
                       "physical_section_id": row["image_id"], "modality": "phenocycler", "tissue": "spleen",
                       "panel": row["panel_id"], "segmentation_id": spleen["defaults"]["segmentation_version"],
                       "cells": result.get("source_cells", pc33_inputs.get(row["donor_id"], {}).get("cells")),
                       "typed_eligible_cells": result.get("cells"), "spots": None, "source": "inventories/spleen_cohort.json",
                       "disposition": "completed" if result["identity_check"] == "pass" else "blocked",
                       "reason": result.get("reason", "Later independent PC67 producer; large raw-image input hashes unavailable"),
                       "evidence": "remote/spleen/per_acquisition.json"})
    pairs = []
    rna_inputs = {r["acquisition_id"]: r for r in read(RUN / "remote/spleen_rna_inputs.json")["acquisitions"]}
    for row in read(sources / "spleen_pair_sources.json")["pairs"]:
        measured = rna_inputs[row["xenium_acquisition_id"]]
        ledger.append({"cohort": "spleen", "acquisition_id": row["xenium_acquisition_id"], "donor_id": row["donor_id"],
                       "physical_section_id": row["protein_entry"], "modality": "xenium", "tissue": "spleen",
                       "panel": "gene_set_sha256:" + measured["gene_panel_sha256"], "segmentation_id": None,
                       "cells": measured["cells"], "spots": None,
                       "source": "inventories/spleen_pair_sources.json", "disposition": "blocked",
                       "reason": "Sources verified; full independent spleen RNA lineage route/specification unavailable (one gene-evidence pilot is not all eight)",
                       "evidence": "remote/spleen_rna_inputs.json"})
        pairs.append({**row, "pair_id": row["protein_entry"] + "__" + row["xenium_acquisition_id"],
                      "status": "blocked", "source_verification": "pass", "registration": "no accepted Stage 07 receipt",
                      "independent_rna_typing": "unavailable full-cohort specification/route", "links": None,
                      "concordance": None, "evidence": "remote/spleen_pair_sources.json"})
    relations = pd.read_csv(sources / "io60_relationships.csv", keep_default_na=False, dtype=str)
    for row in relations[(relations.kind == "same_section") & relations.evidence.str.startswith("section_modes")].to_dict("records"):
        # Match the *recorded serial* to both source inventories, never donor alone.
        serial = row["b"].rsplit(" ", 1)[-1]
        p = inventory[(inventory.platform == "phenocycler") & (inventory.xenium_slide == serial)]
        r = inventory[(inventory.platform == "xenium") & (inventory.acquisition_id == f"xenium:{serial}:panc")]
        pairs.append({"pair_id": "io60:" + serial + ":panc", "tissue": "pancreas", "relationship": "same_section",
                      "protein_acquisition_ids": p.acquisition_id.tolist(), "rna_acquisition_ids": r.acquisition_id.tolist(),
                      "source_verification": "inventory declaration; full frame/segmentation binding unavailable",
                      "evidence": row["evidence"], "status": "blocked", "registration": "no accepted Stage 07 receipt",
                      "links": None, "concordance": None})
    for row in ledger:
        row["cell_count_scope"] = "segmented input cells in this acquisition" if row.get("cells") is not None else "unavailable"
        row["assay"] = row["modality"]
        row["measurement_modality"] = {"codex": "protein", "phenocycler": "protein", "if": "protein", "xenium": "rna", "maldi": "maldi"}[row["modality"]]
    write(DOC / "acquisition_dispositions.json", ledger)
    pd.DataFrame(ledger).to_csv(DOC / "acquisition_dispositions.csv", index=False)
    write(DOC / "pair_dispositions.json", {"explicit_pairs": pairs,
          "unresolved": [{"tissue": "pancreatic_lymph_node", "status": "blocked",
                          "reason": "pLN regions occur on four recorded serials, but the inventory does not assert their physical same-section relationship; donor/slide agreement cannot create it"}],
          "not_additional_pairs": "Four donor_notes same_section rows repeat the four section_modes declarations"})
    write(DOC / "cohort_accounting.json", {"acquisitions": len(ledger),
          "by_cohort_modality_disposition": pd.DataFrame(ledger).groupby(["cohort", "modality", "disposition"]).size().reset_index(name="count").to_dict("records"),
          "explicit_same_slide_pairs": len(pairs),
          "central_exclusions": [{"donor": "6457", "platform": "phenocycler", "reason": "existing platform-specific exclusion; absent from frozen PhenoCycler inventory", "applies_to_xenium": False, "applies_to_maldi": False}],
          "unknown": "null section/panel/segmentation/count fields mean unavailable, never inferred from donor IDs"})
    readiness = []
    for (modality, tissue, panel), group in pd.DataFrame(ledger).groupby(["modality", "tissue", "panel"], dropna=False):
        if modality == "if":
            continue
        tissue = None if pd.isna(tissue) else tissue
        panel = None if pd.isna(panel) else panel
        axes = {"code_availability": "available", "valid_inputs": "partially_verified",
                "scientific_specification": "version_specific; portability review required",
                "executable_routing": "available", "computational_validation": "see per-acquisition dispositions",
                "biological_review": "not performed; no acceptance asserted"}
        remedy = []
        if modality == "codex":
            axes.update(valid_inputs="blocked: unbound/stale table producer", computational_validation="20/20 provenance guards blocked before numerical comparison")
            remedy = ["PROV-1", "BIO-1", "CORE-2"]
        elif modality == "xenium" and tissue != "spleen":
            axes.update(valid_inputs="receipt hashes and every native/qualified ID validated",
                        computational_validation="all 65 across both tissues accounted; 60 candidates and five expected unavailable; state counts independently reconciled",
                        scientific_specification="tissue-specific identity cores and fixed thresholds; actual panel coverage in reports; absent genes are not zeros")
            remedy = ["RNA-1", "RNA-2"]
        elif modality == "xenium":
            axes.update(code_availability="bounded gene-evidence pilot only", valid_inputs="all eight source bundles verified; 5,001 genes with identical gene-set hash; 4,167,720 cell IDs reconciled and nucleus admission measured",
                        executable_routing="full independent spleen RNA typing unavailable",
                        scientific_specification="full spleen RNA lineage specification missing", computational_validation="blocked for eight; pilot is not substituted")
            remedy = ["RNA-2", "GEOM-1"]
        elif modality == "maldi":
            axes.update(valid_inputs="both source hashes and paired spot keys verified; tissue/neighbour inspections blank",
                        executable_routing="ingest available; registration blocked; Xenium neighbour not built",
                        computational_validation="all 32 sections / 1,666,450 spots audited for keys, finite values, mask and denominators; no accepted transforms/objects",
                        scientific_specification="serial spot/object evidence only; physical tissue, intensity calibration and pLN object semantics unresolved")
            remedy = ["MALDI-1", "MALDI-2", "MALDI-3"]
        elif tissue == "spleen" and panel == "spleen_pc33":
            axes.update(code_availability="no PC33 registry/lineage route", scientific_specification="lineage roots absent; four sections deliberately not typed",
                        executable_routing="unsupported", computational_validation="blocked, all four retained")
            remedy = ["PANEL-2"]
        elif tissue == "spleen":
            axes.update(code_availability="later a9549cc producer; pinned df3fe1d port proposals remain separate",
                        computational_validation="all eight native eligible cell universes reconciled; 4,346,039 cells",
                        valid_inputs="77 receipt entries checked; 10 large image/geometry input hashes unavailable")
            remedy = ["PANEL-1", "BIO-3", "MEAS-1"]
        elif tissue == "pancreatic_lymph_node":
            axes.update(code_availability="legacy candidate route accepts configured rosters; per-image pLN port unsupported",
                        executable_routing="only two sections have bound upstream artifacts; remaining 17 blocked",
                        valid_inputs="two pilot sections verified; old manifest recovered exactly at producer pin",
                        computational_validation="two sections / 832,062 cells reconciled; 17 blocked; no final classes")
            remedy = ["PANEL-3", "BIO-1", "BIO-2"]
        elif modality == "phenocycler":
            axes.update(code_availability="legacy full25 mixed-tissue route; new per-image port has no PhenoCycler pancreas registry",
                        computational_validation="all 25 primary tables produced; refinement failed on 0112; primary/source compatibility audited separately; no final receipt",
                        scientific_specification="pancreas context, hormone references and donor-derived scale remain version-specific")
            remedy = ["MEAS-1", "BIO-1", "BIO-2", "CORE-2", "PROT-1"]
        readiness.append({"modality": {"codex": "protein", "phenocycler": "protein", "xenium": "rna", "maldi": "maldi"}[modality],
                          "assay": modality, "tissue": tissue, "panel": panel,
                          "relationship_type": "serial_section" if modality == "maldi" else "independent_acquisition",
                          "acquisitions": group.acquisition_id.tolist(), "axes": axes, "remediations": remedy,
                          "evidence": sorted(set(group.evidence.dropna()))})
    for tissue in ("pancreas", "pancreatic_lymph_node", "spleen"):
        readiness.append({"modality": "protein/rna", "assay": "PhenoCycler/Xenium", "tissue": tissue, "panel": "PC67/RNA" if tissue == "spleen" else None,
                          "relationship_type": "same_section", "axes": {
                              "code_availability": "one-pair spleen candidate-link pilot; shared matching functions",
                              "valid_inputs": "eight verified spleen source pairs / four pancreas declarations; pLN physical relationships unresolved",
                              "scientific_specification": "hierarchy mapping, split/merge policy and uncertainty propagation missing",
                              "executable_routing": "bounded single-pair pilot; no full-cohort route/accepted transforms",
                              "computational_validation": "contract tests and adversarial probes only; no real same-cell links validated",
                              "biological_review": "not performed; concordance cannot establish truth"},
                          "remediations": ["GEOM-1", "GEOM-2"], "evidence": ["pair_dispositions.json", "fault_probes/results.json"]})
    write(DOC / "readiness.json", {"schema": "dynamic-portability-readiness/1", "routes": readiness,
                                  "interpretation": "Axes are independent; code availability is not scientific acceptance"})


if __name__ == "__main__":
    main()
