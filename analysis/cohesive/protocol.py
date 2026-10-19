"""Freeze the spatial_v2 protocol with the maintainer's 2026-10-17 amendment, before any endpoint."""
from __future__ import annotations

import copy
from pathlib import Path

import pandas as pd

from .common import Context, binding, read_json, require, write_json


AMENDMENT = Path(__file__).with_name("protocol_amendment_20261017.json")
CROSS_STUDY_AMENDMENT = Path(__file__).with_name("protocol_amendment_20261018.json")

# What each endpoint cannot be computed without. Not support cutoffs: a donor lacking one of these
# is UNAVAILABLE for that endpoint because the quantity is undefined, not because it is small.
COMPUTABILITY = {
    "H1": "At least one single-part true islet with a valid outline, refined immune cells in its boundary "
          "kernel, vascular and neural objects in the section, and a vascular-minus-neural proximity "
          "profile that varies around the boundary.",
    "H2": "At least one true islet with refined Beta cells and at least one other refined endocrine "
          "subtype, so the composition-preserving shuffle is not the identity. Both components are "
          "reported; the joint prediction is claimed only where their directions agree.",
    "H3": "Refined endocrine cells with a measured HLA_A, a homotypic graph with at least one edge, "
          "and residual HLA_A variance after the nuisance fit.",
    "H4": "Refined endocrine cells with measured HLA_A and PD_L1, refined T cells with measured PD_1 "
          "and Granzyme_B or CD107a within the neighbourhood, and a full-rank design with a nonzero "
          "interaction residual.",
    "H5": "Assayed CCL19 or CCL21 and LTB, released B/plasma and T identities inside a localized "
          "islet, and variable organizer, B-LTB and T states.",
    "H6": "A computable H3 endpoint for the same donor, released pLN T identities, and a conditioning "
          "design that leaves residual variance across paired donors.",
}
REQUIRED_PROTEIN_MARKERS = {"H1": (), "H2": ("INS", "GCG", "SST"), "H3": ("HLA_A",),
                            "H4": ("HLA_A", "PD_L1", "PD_1"), "H5": (), "H6": ("HLA_A",)}
# H6 is the archived script's TCF7/LEF1 versus TOX/GZMB/PRF1/NKG7 contrast (evaluate_spatial_RNA_H5_H6.py:19).
REQUIRED_RNA_GENES = {"H5": ("CCL19", "CCL21", "LTB"),
                      "H6": ("TCF7", "LEF1", "TOX", "GZMB", "PRF1", "NKG7")}


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "freeze", "Unknown protocol operation")
    return [binding(AMENDMENT), binding(CROSS_STUDY_AMENDMENT)]


def amended(original: dict, amendment: dict) -> dict:
    require(original["family"] == amendment["amends"]["family"], "Amendment names another protocol")
    frozen = copy.deepcopy(original)
    for hypothesis in frozen["hypotheses"]:
        removed = hypothesis.pop("minimum_support")
        hypothesis["removed_minimum_support"] = removed
        hypothesis["computability"] = COMPUTABILITY[hypothesis["id"]]
        hypothesis["eligibility"] = hypothesis["eligibility"].replace(
            "Native Beta/Alpha/Delta", "Refined Beta/Alpha/Delta").replace(
            "native T lineage", "refined T lineage").replace(
            "Native endocrine", "Refined endocrine")
    frozen["statistics"]["removed_minimum_disease_support"] = frozen["statistics"].pop("minimum_disease_support")
    for key in ("pancreatitis", "protein_identity"):
        frozen["sensitivity"].pop(key)
    frozen["protein_identity"] = next(row["effect"] for row in amendment["rulings"]
                                      if row["id"] == "refinement_is_phenotyping")
    frozen["unavailable"] = ("An endpoint that cannot be computed for a donor is UNAVAILABLE and shown as "
                             "such; statistics.unavailable_p applies to multiplicity bookkeeping only and is "
                             "never displayed as a measured result.")
    frozen["amendment"] = amendment
    frozen["status"] = "frozen before any new endpoint; retrospective/exploratory"
    return frozen


def cross_study(frozen: dict, amendment: dict) -> dict:
    """The 2026-10-18 rulings: 6551 stays everywhere and anchors a cross-study accuracy report."""
    require(amendment["amends"]["after"] == AMENDMENT.name, "Cross-study amendment applies after 2026-10-17")
    out = copy.deepcopy(frozen)
    effects = {row["id"]: row["effect"] for row in amendment["rulings"]}
    out["sensitivity"]["removed_source_overlap"] = out["sensitivity"].pop("source_overlap")
    out["cross_study_accuracy"] = effects["literature_overlap_is_cross_study_accuracy"]
    out["protein_identity_coverage"] = effects["four_same_slide_donors_typed_like_the_rest"]
    out["second_amendment"] = amendment
    return out


def feasibility(context: Context) -> pd.DataFrame:
    qualified = context.task("qualify")
    markers = pd.read_csv(qualified / "protein_marker_support.csv", dtype={"donor_id": str})
    measured = markers[markers.measured > 0].groupby("donor_id").marker.agg(set)
    rna = read_json(qualified / "rna_measurement_support.json")
    rows = []
    for donor, present in measured.items():
        for hypothesis, needed in REQUIRED_PROTEIN_MARKERS.items():
            missing = sorted(set(needed) - present)
            rows.append({"hypothesis": hypothesis, "modality": "protein", "scope": donor,
                         "required": ";".join(needed), "unmeasured": ";".join(missing),
                         "measurements_available": not missing})
    for scope in rna:
        for hypothesis, needed in REQUIRED_RNA_GENES.items():
            missing = sorted(set(needed) & set(scope["unavailable_genes"]))
            rows.append({"hypothesis": hypothesis, "modality": "RNA", "scope": scope["scope_id"],
                         "required": ";".join(needed), "unmeasured": ";".join(missing),
                         "measurements_available": not missing})
    return pd.DataFrame(rows)


def freeze(context: Context, stage: Path) -> dict:
    frozen = cross_study(amended(read_json(context.protocol), read_json(AMENDMENT)),
                         read_json(CROSS_STUDY_AMENDMENT))
    write_json(stage / "frozen_protocol.json", frozen)
    table = feasibility(context)
    table.to_csv(stage / "measurement_feasibility.csv", index=False)
    return {"hypotheses": [row["id"] for row in frozen["hypotheses"]],
            "support_cutoffs": 0, "scopes_lacking_a_required_measurement":
                int((~table.measurements_available).sum()),
            "feasibility_rows": len(table)}
