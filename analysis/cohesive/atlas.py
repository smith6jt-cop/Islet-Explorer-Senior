"""Atlas overview: cohort coverage, qualification flow and true-Islet composition, descriptive only.

Identity is the refined hierarchical identity (broad label kept by refinement, then its subclass).
Islet density uses the actual outline area of the original true Islet; there is no expansion.
Pooled ratios and unweighted means of per-islet ratios are different estimands and both are kept.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import geometry
from .common import Context, binding, parallel_map, peak_rss_gib, read_json, require, write_json
from .endpoints import ENDOCRINE, identity
from .sources import audit_binding

POSTIF_COVERAGE = "results/v3/islet_immune_reassessment_20261006/tables/qc_postIF_islet_coverage.csv"
#: The cohort scale the typing used; T-cell state positivity is read at its 0.10 level, like the rules.
SCALE = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008/own_compartment/scale.parquet")
LEVEL = 0.10
STATE_MARKERS = ("TCF_1", "TOX", "PD_1", "CD8", "CD4", "FOXP3", "Ki67", "LAG3", "Granzyme_B")
RNA_T = {"pancreas_RNA": ("Immune_T",),
         "pLN_RNA": ("Immune_T", "T_CD4", "T_CD8", "T_FOXP3_detected_pattern", "T_CD4_CD8_both", "T_subset_unresolved")}
MEMORY_GENES = ("TCF7", "LEF1")
CYTOTOXIC_GENES = ("TOX", "GZMB", "PRF1", "NKG7")
STRUCTURES = {"Blood vessel-associated allocation": "vascular", "Nerve-associated allocation": "neural",
              "Lymphatic-associated allocation": "lymphatic"}
GIB_PER_TASK = 4
CELL_COLUMNS = ["object_id", "endocrine_object_id", "endocrine_object_class", "immune_class",
                "annotation_name", "broad_endocrine", "broad_immune", "structure_associated",
                "refinement_keeps_broad_label", "disease_status"]


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "build_tables", "Unknown atlas operation")
    return [binding(context.task("qualify") / "completion.json"), binding(Path(__file__).with_name("endpoints.py")),
            audit_binding(context, geometry.QUALIFIED_OBJECTS), binding(context.original_repo / POSTIF_COVERAGE),
            binding(SCALE)]


def category(cells: pd.DataFrame) -> np.ndarray:
    """One display category per cell from the refined identity; structure allocations stay structural."""
    label, t_cell = identity(cells)
    name = cells.annotation_name.astype("string").fillna("").to_numpy(str)
    kept = cells.refinement_keeps_broad_label.fillna(False).to_numpy(bool)
    out = np.where(label != "", label, "Identity unresolved").astype(object)
    immune = label == "Immune"
    out[immune] = cells.immune_class.astype("string").fillna("Imm-unassigned").to_numpy(str)[immune]
    for source, target in STRUCTURES.items():
        out[(name == source) & kept & cells.structure_associated.fillna(False).to_numpy(bool)] = f"structure: {target}"
    demoted = (name != "") & (name != "Identity unresolved") & ~kept
    out[demoted] = "Broad label not kept by refinement"
    return out.astype(str)


def t_states(cells: pd.DataFrame, t_cell: np.ndarray, member: np.ndarray, donor: str) -> pd.DataFrame:
    """TCF1/TOX state shares and marker-positive fractions of refined T cells, inside true Islets and outside."""
    scale = pd.read_parquet(SCALE).set_index("marker")
    positive = {m: np.clip((cells[f"redsea_aligned_{m}"].to_numpy(float) - scale.at[m, "floor"]) / scale.at[m, "span"],
                           0, 1) >= LEVEL for m in STATE_MARKERS}
    measured = {m: np.isfinite(cells[f"redsea_aligned_{m}"].to_numpy(float)) for m in STATE_MARKERS}
    rows = []
    for compartment, where in (("inside true Islets", member), ("outside true Islets", ~member)):
        t = t_cell & where
        both = t & measured["TCF_1"] & measured["TOX"]
        row = {"donor_id": donor, "compartment": compartment, "T_cells": int(t.sum()),
               "T_cells_TCF1_TOX_measured": int(both.sum())}
        for name, (a, b) in {"TCF1+TOX-": (True, False), "TCF1+TOX+": (True, True),
                             "TCF1-TOX+": (False, True), "TCF1-TOX-": (False, False)}.items():
            hit = both & (positive["TCF_1"] == a) & (positive["TOX"] == b)
            row[f"share_{name}"] = float(hit.sum() / both.sum()) if both.any() else np.nan
        for m in STATE_MARKERS:
            denominator = t & measured[m]
            row[f"fraction_{m}_positive"] = float((positive[m] & denominator).sum() / denominator.sum()) if denominator.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def donor_tables(task: dict) -> dict:
    started = time.perf_counter()
    cells = pd.read_parquet(task["cells"], columns=CELL_COLUMNS + ["immune_subtype"] +
                            [f"redsea_aligned_{m}" for m in STATE_MARKERS])
    cells["category"] = category(cells)
    tissue = cells.category.value_counts().rename_axis("category").reset_index(name="cells")
    tissue.insert(0, "donor_id", task["donor"])
    member = cells.endocrine_object_class.eq("Islet").fillna(False)
    counts = (cells.loc[member].groupby(["endocrine_object_id", "category"]).size()
              .unstack(fill_value=0))
    label, t_cell = identity(cells)
    immune = label == "Immune"
    subtypes = (cells.loc[immune].assign(compartment=np.where(member[immune], "inside true Islets", "outside true Islets"),
                                         subtype=cells.loc[immune, "immune_subtype"].astype("string").fillna("no subtype rule fired"),
                                         immune_class=cells.loc[immune, "immune_class"].astype("string").fillna("unassigned"))
                .groupby(["compartment", "immune_class", "subtype"]).size().rename("cells").reset_index())
    subtypes.insert(0, "donor_id", task["donor"])
    return {"tissue": tissue, "islet_counts": counts.assign(donor_id=task["donor"]).reset_index(),
            "immune_subtypes": subtypes, "t_states": t_states(cells, t_cell, member.to_numpy(bool), task["donor"]),
            "resource": {"donor_id": task["donor"], "peak_rss_gib": peak_rss_gib(),
                         "seconds": time.perf_counter() - started}}


def rna_t_states(task: dict) -> dict:
    """Memory (TCF7/LEF1) and cytotoxic (TOX/GZMB/PRF1/NKG7) programs of a section's released T cells."""
    started = time.perf_counter()
    genes = [f"RNA_{g}" for g in MEMORY_GENES + CYTOTOXIC_GENES]
    cells = pd.read_parquet(task["cells"], columns=["candidate_lineage", "handoff_named_lineage_use_eligible",
                                                    "raw_library_counts", *genes])
    t = cells.loc[cells.handoff_named_lineage_use_eligible.fillna(False)
                  & cells.candidate_lineage.astype(str).isin(RNA_T[task["lane"]])]
    depth = t.raw_library_counts.to_numpy(float)

    def program(names):
        return np.mean(np.column_stack([np.log1p(1e4 * t[f"RNA_{g}"].to_numpy(float) / depth) for g in names]), axis=1)

    row = {**{k: task[k] for k in ("scope_id", "lane", "donor_id", "section_id")}, "T_cells": len(t)}
    if len(t):
        memory, cytotoxic = program(MEMORY_GENES), program(CYTOTOXIC_GENES)
        row.update(mean_memory_program=float(memory.mean()), mean_cytotoxic_program=float(cytotoxic.mean()),
                   fraction_any_memory_count=float((t[[f"RNA_{g}" for g in MEMORY_GENES]].sum(axis=1) > 0).mean()),
                   fraction_any_cytotoxic_count=float((t[[f"RNA_{g}" for g in CYTOTOXIC_GENES]].sum(axis=1) > 0).mean()),
                   **{f"T_{lineage}": int((t.candidate_lineage.astype(str) == lineage).sum()) for lineage in RNA_T[task["lane"]]})
    return {"row": row, "resource": {"scope_id": task["scope_id"], "peak_rss_gib": peak_rss_gib(),
                                     "seconds": time.perf_counter() - started}}


def islet_table(objects: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    islets = objects.loc[objects.original_object_class.eq("Islet"),
                         ["donor_id", "donor_status", "endocrine_object_id", "islet_key", "outline_area_um2",
                          "outline_area_available", "n_original_cells_in_object", "n_qualified_member"]]
    table = islets.merge(counts, on=["donor_id", "endocrine_object_id"], how="left", validate="one_to_one")
    categories = [c for c in counts.columns if c not in ("donor_id", "endocrine_object_id")]
    table[categories] = table[categories].fillna(0).astype(int)
    table["admitted_nucleated_members"] = table[categories].sum(axis=1)
    endocrine = [c for c in categories if c in ENDOCRINE or c == "Endocrine (subtype unresolved)"]
    immune = [c for c in categories if c.startswith("Imm-")]
    table["refined_endocrine"] = table[endocrine].sum(axis=1)
    table["refined_subtyped_endocrine"] = table[[c for c in endocrine if c in ENDOCRINE]].sum(axis=1)
    table["refined_immune"] = table[immune].sum(axis=1)
    area_mm2 = table.outline_area_um2 / 1e6
    table["immune_per_mm2"] = table.refined_immune / area_mm2
    table["immune_fraction_of_members"] = table.refined_immune / table.admitted_nucleated_members.where(
        table.admitted_nucleated_members > 0)
    for subtype in ENDOCRINE:
        table[f"{subtype.lower()}_fraction_of_subtyped_endocrine"] = table.get(subtype, 0) / \
            table.refined_subtyped_endocrine.where(table.refined_subtyped_endocrine > 0)
    return table


def donor_summary(islets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, z in islets.groupby("donor_id", sort=True):
        area = z.outline_area_um2.sum() / 1e6
        row = {"donor_id": donor, "disease_status": z.donor_status.iloc[0], "true_islets": len(z),
               "true_islets_with_admitted_member": int(z.admitted_nucleated_members.gt(0).sum()),
               "total_outline_area_mm2": area, "admitted_nucleated_members": int(z.admitted_nucleated_members.sum()),
               "refined_endocrine": int(z.refined_endocrine.sum()), "refined_immune": int(z.refined_immune.sum()),
               "immune_per_mm2_pooled": z.refined_immune.sum() / area if area > 0 else np.nan,
               "immune_per_mm2_mean_of_islets": z.immune_per_mm2.mean(),
               "immune_fraction_pooled": z.refined_immune.sum() / z.admitted_nucleated_members.sum()
               if z.admitted_nucleated_members.sum() else np.nan,
               "immune_fraction_mean_of_islets": z.immune_fraction_of_members.mean(),
               "fraction_of_islets_with_immune": float(z.refined_immune.gt(0).mean())}
        for subtype in ENDOCRINE:
            total = z.refined_subtyped_endocrine.sum()
            row[f"{subtype.lower()}_fraction_pooled"] = z.get(subtype, pd.Series(0)).sum() / total if total else np.nan
            row[f"{subtype.lower()}_fraction_mean_of_islets"] = z[f"{subtype.lower()}_fraction_of_subtyped_endocrine"].mean()
        for name in [c for c in z.columns if c.startswith("Imm-")]:
            row[f"{name}_per_mm2_pooled"] = z[name].sum() / area if area > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def build_tables(context: Context, stage: Path) -> dict:
    qualified = context.task("qualify")
    protein = pd.read_csv(qualified / "protein_scope_qualification.csv", dtype={"donor_id": str})
    rna = pd.read_csv(qualified / "rna_scope_qualification.csv", dtype={"donor_id": str, "source_serial": str})
    tasks = [{"donor": row.donor_id, "cells": str(qualified / "protein" / f"donor_id={row.donor_id}" / "cells.parquet"),
              "rows": int(row.biological_rows)} for row in protein.itertuples()]
    results = parallel_map(donor_tables, tasks, GIB_PER_TASK, "atlas", cost=[t["rows"] for t in tasks])
    disease_of = protein.set_index("donor_id").disease_status
    pd.concat([r["immune_subtypes"] for r in results], ignore_index=True).assign(
        disease_status=lambda t: t.donor_id.map(disease_of)).to_csv(stage / "protein_immune_subtypes.csv", index=False)
    pd.concat([r["t_states"] for r in results], ignore_index=True).assign(
        disease_status=lambda t: t.donor_id.map(disease_of)).to_csv(stage / "protein_t_cell_states.csv", index=False)
    rna_tasks = [{"scope_id": r.scope_id, "lane": r.lane, "donor_id": r.donor_id, "section_id": r.section_id,
                  "cells": str(qualified / r.relative_data_path), "rows": int(r.biological_rows)} for r in rna.itertuples()]
    rna_results = parallel_map(rna_t_states, rna_tasks, GIB_PER_TASK, "atlas RNA T states", cost=[t["rows"] for t in rna_tasks])
    pd.DataFrame([r["row"] for r in rna_results]).assign(
        disease_status=lambda t: t.donor_id.map(disease_of)).to_csv(stage / "rna_t_cell_states.csv", index=False)
    tissue = pd.concat([r["tissue"] for r in results], ignore_index=True)
    tissue = tissue.merge(protein[["donor_id", "disease_status"]], on="donor_id", validate="many_to_one")
    tissue.to_csv(stage / "protein_tissue_composition.csv", index=False)
    objects = pd.read_parquet(audit_binding(context, geometry.QUALIFIED_OBJECTS)["path"])
    objects["donor_id"] = objects.donor_id.astype(str)
    counts = pd.concat([r["islet_counts"] for r in results], ignore_index=True).fillna(0)
    islets = islet_table(objects, counts)
    require(len(islets) == int(objects.original_object_class.eq("Islet").sum()), "True-Islet universe changed")
    islets.to_parquet(stage / "true_islet_composition.parquet", index=False)
    summary = donor_summary(islets)
    summary.to_csv(stage / "donor_true_islet_summary.csv", index=False)

    disease = protein.set_index("donor_id").disease_status
    rna_counts = pd.read_csv(qualified / "rna_annotation_counts.csv", dtype={"donor_id": str, "source_serial": str})
    rna_counts["disease_status"] = rna_counts.donor_id.map(disease)
    require(rna_counts.disease_status.notna().all(), "An RNA donor has no protein disease authority")
    rna_counts.to_csv(stage / "rna_section_composition.csv", index=False)
    postif = pd.read_csv(context.original_repo / POSTIF_COVERAGE, dtype={"donor_id": str, "xenium_serial": str})
    coverage = []
    for row in protein.itertuples():
        sections = rna.loc[rna.donor_id.eq(row.donor_id)]
        mask = postif.loc[postif.donor_id.eq(row.donor_id)]
        coverage.append({
            "donor_id": row.donor_id, "disease_status": row.disease_status,
            "protein_pancreas_cells": int(row.biological_rows),
            "protein_true_islets": int(summary.set_index("donor_id").true_islets.get(row.donor_id, 0)),
            "protein_immune_cells": int(row.broad_immune),
            "rna_pancreas_sections": int(sections.lane.eq("pancreas_RNA").sum()),
            "rna_pancreas_cells": int(sections.loc[sections.lane.eq("pancreas_RNA"), "biological_rows"].sum()),
            "rna_pLN_sections": int(sections.lane.eq("pLN_RNA").sum()),
            "rna_pLN_cells": int(sections.loc[sections.lane.eq("pLN_RNA"), "biological_rows"].sum()),
            "postIF_islet_masks": int(mask.n_retained_islets.fillna(0).sum()) if mask.mask_available.any() else 0,
            "postIF_state": "available" if mask.mask_available.any() else
            ("; ".join(mask.exclusion_reason.dropna().astype(str)) or "no post-IF record")})
    coverage = pd.DataFrame(coverage)
    coverage.to_csv(stage / "cohort_coverage.csv", index=False)
    flow = protein[["donor_id", "disease_status", "source_rows", "QC_rows", "admission_rows", "biological_rows"]]
    flow.to_csv(stage / "protein_qualification_flow.csv", index=False)
    rna.assign(disease_status=rna.donor_id.map(disease)).to_csv(stage / "rna_qualification_flow.csv", index=False)
    pd.DataFrame([r["resource"] for r in results]).to_csv(stage / "resource_usage.csv", index=False)
    result = {"donors": len(protein), "disease_counts": protein.disease_status.value_counts().to_dict(),
              "true_islets": len(islets), "true_islets_with_admitted_member": int(islets.admitted_nucleated_members.gt(0).sum()),
              "rna_sections": rna.lane.value_counts().to_dict(),
              "postIF_mask_donors": int(coverage.postIF_state.eq("available").sum())}
    write_json(stage / "atlas_summary.json", result)
    return result
