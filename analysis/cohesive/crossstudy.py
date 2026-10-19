"""Cross-study accuracy on donor 6551 (ruled 2026-10-18: "use it to infer cross-study accuracy").

Peters et al. (bioRxiv 10.1101/2025.05.23.655848, Figure 7) analysed the same donor's PhenoCycler image
with their own segmentation, gates and QuPath classifiers. Their analysis object
(`/data/6551_manuscript/6551_leiden_umap_tcf_tox_analysis.h5ad`) holds the 26,781 CD3e+ T cells they
classed TCF1+TOX-, TCF1+TOX+ or TCF1-TOX+, with tissue (pancreas / pLN) and signed distance to islets;
it reproduces the published shares exactly. Here the same quantities are computed from this pipeline:
T cells by the hierarchical identity in the pancreas, and by the producer's Imm-T-cell definition in
both tissues; TCF1 and TOX positive at the cohort scale's 0.10 level. Different segmentation and
islet definitions are part of what is compared.

Audit context only: agreement is reported, never used to tune, gate or release anything, and it does
not verify any biological finding.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely

from . import geometry
from .common import Context, binding, read_json, require, threaded, write_json
from .endpoints import identity
from .sources import audit_binding

PETERS = Path("/data/6551_manuscript/6551_leiden_umap_tcf_tox_analysis.h5ad")
GEN = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008")
UNIFORM = GEN / "own_compartment/uniform/donor_id=6551/data.parquet"
SCALE = GEN / "own_compartment/scale.parquet"
NUCLEUS = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/core/phenocycler/data/redsea/donor_id=6551/data_0.parquet")
LEVEL = 0.10
STATES = ("TCF1+TOX-", "TCF1+TOX+", "TCF1-TOX+")
PETERS_STATES = {"CD3e: TCF-1": "TCF1+TOX-", "CD3e: TCF-1: TOX": "TCF1+TOX+", "CD3e: TOX": "TCF1-TOX+"}
BINS = [-np.inf, -50, 0, 50, 100, np.inf]
BIN_LABELS = ["< -50", "-50 to 0", "0 to 50", "50 to 100", "> 100"]
EPITHELIAL = ("E_cadherin", "EpCAM", "Pan_Cytokeratin")


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "compare", "Unknown cross-study operation")
    return [binding(context.task("qualify") / "completion.json"), binding(Path(geometry.__file__)),
            audit_binding(context, geometry.QUALIFIED_OBJECTS),
            *threaded(binding, [PETERS, UNIFORM, SCALE, NUCLEUS])]


def scaled(frame: pd.DataFrame, scale: pd.DataFrame, marker: str, column: str) -> np.ndarray:
    return np.clip((frame[column].to_numpy(float) - scale.at[marker, "floor"]) / scale.at[marker, "span"], 0, 1)


def state(tcf1: np.ndarray, tox: np.ndarray) -> np.ndarray:
    out = np.full(len(tcf1), "TCF1-TOX-", dtype=object)
    out[(tcf1 >= LEVEL) & (tox < LEVEL)] = "TCF1+TOX-"
    out[(tcf1 >= LEVEL) & (tox >= LEVEL)] = "TCF1+TOX+"
    out[(tcf1 < LEVEL) & (tox >= LEVEL)] = "TCF1-TOX+"
    out[~(np.isfinite(tcf1) & np.isfinite(tox))] = "unmeasured"
    return out


def shares(states: np.ndarray) -> dict:
    """Shares among TCF1- and/or TOX-positive T cells: the published Figure 7 denominator."""
    kept = np.isin(states, STATES)
    n = int(kept.sum())
    return {"n": n, **{s: float((states[kept] == s).mean() * 100) if n else np.nan for s in STATES}}


def peters_reference() -> tuple[pd.DataFrame, pd.DataFrame]:
    import anndata
    obs = anndata.read_h5ad(PETERS, backed="r").obs
    states = obs.TCF_TOX_population.astype(str).map(PETERS_STATES).to_numpy()
    require(pd.notna(states).all(), "Unexpected Peters T-cell population label")
    rows = [{"study": "Peters et al.", "tissue": t, **shares(states[obs.Tissue.astype(str).eq(t).to_numpy()])}
            for t in ("Pancreas", "pLN")]
    pancreas = obs.Tissue.astype(str).eq("Pancreas").to_numpy()
    bins = pd.cut(obs.SignedDistanceIslet[pancreas], BINS, right=False, labels=BIN_LABELS)
    distance = [{"study": "Peters et al.", "distance_bin": b, **shares(states[pancreas][bins.to_numpy() == b])}
                for b in BIN_LABELS]
    return pd.DataFrame(rows), pd.DataFrame(distance)


def ours(context: Context) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    scale = pd.read_parquet(SCALE).set_index("marker")
    cells = pd.read_parquet(context.task("qualify") / "protein" / "donor_id=6551" / "cells.parquet")
    label, typed_t = identity(cells)
    s = {m: scaled(cells, scale, m, f"redsea_aligned_{m}") for m in ("CD3e", "TCF_1", "TOX", *EPITHELIAL)}
    rule_t = (s["CD3e"] >= LEVEL) & np.all([s[m] < LEVEL for m in EPITHELIAL], axis=0)
    pancreas_state = state(s["TCF_1"], s["TOX"])
    # pLN: the same-block node, aligned with the pancreas (2026-10-02); no typing producer names its cells
    # with the pancreas rules, so only the Imm-T-cell definition applies there.
    node = pq.read_table(UNIFORM, columns=["object_id", "qc_analysis_eligible", "in_alignment_block", "cell_region",
                                           "CD3e", "TCF_1", "TOX", *EPITHELIAL]).to_pandas()
    node = node.loc[node.cell_region.eq("pLN") & node.in_alignment_block.fillna(False) & node.qc_analysis_eligible.fillna(False)]
    nucleus = pq.read_table(NUCLEUS, columns=["object_id", "nucleus_area_px"]).to_pandas()
    node = node.merge(nucleus, on="object_id", how="left", validate="one_to_one")
    node = node.loc[pd.to_numeric(node.nucleus_area_px) > 0]
    n = {m: scaled(node, scale, m, m) for m in ("CD3e", "TCF_1", "TOX", *EPITHELIAL)}
    node_rule_t = (n["CD3e"] >= LEVEL) & np.all([n[m] < LEVEL for m in EPITHELIAL], axis=0)
    node_state = state(n["TCF_1"], n["TOX"])
    rows = [{"study": "this pipeline, hierarchical T identity", "tissue": "Pancreas", **shares(pancreas_state[typed_t])},
            {"study": "this pipeline, Imm-T-cell definition", "tissue": "Pancreas", **shares(pancreas_state[rule_t])},
            {"study": "this pipeline, Imm-T-cell definition", "tissue": "pLN", **shares(node_state[node_rule_t])}]
    islets = geometry.qualified_islets(Path(audit_binding(context, geometry.QUALIFIED_OBJECTS)["path"]), "6551")
    images = geometry.image_records(context)
    outlines, _ = geometry.true_islets(geometry.donor_files(context, "6551")["islets"], images["6551"]["pixel_um"],
                                       set(islets.endocrine_object_id))
    union = shapely.union_all(list(outlines.values()))
    points = shapely.points(cells.loc[typed_t, ["X_centroid", "Y_centroid"]].to_numpy(float))
    signed = np.where(shapely.covers(union, points), -shapely.distance(points, union.boundary), shapely.distance(points, union))
    bins = pd.cut(signed, BINS, right=False, labels=BIN_LABELS)
    distance = [{"study": "this pipeline, hierarchical T identity", "distance_bin": b,
                 **shares(pancreas_state[typed_t][bins.to_numpy() == b])} for b in BIN_LABELS]
    counts = {"pancreas_typed_T": int(typed_t.sum()), "pancreas_rule_T": int(rule_t.sum()),
              "pLN_rule_T": int(node_rule_t.sum()), "pLN_cells": int(len(node)),
              "pancreas_typed_T_unmeasured_TCF1_or_TOX": int((pancreas_state[typed_t] == "unmeasured").sum())}
    return pd.DataFrame(rows), pd.DataFrame(distance), counts


def compare(context: Context, stage: Path) -> dict:
    p_states, p_distance = peters_reference()
    o_states, o_distance, counts = ours(context)
    states = pd.concat([p_states, o_states], ignore_index=True)
    distance = pd.concat([p_distance, o_distance], ignore_index=True)
    states.to_csv(stage / "tcf1_tox_shares_by_tissue.csv", index=False)
    distance.to_csv(stage / "tcf1_tox_shares_by_islet_distance.csv", index=False)
    agreement = []
    for _, row in o_states.iterrows():
        ref = p_states.set_index("tissue").loc[row.tissue]
        tvd = 0.5 * sum(abs(row[s] - ref[s]) for s in STATES)
        agreement.append({"comparison": f"{row.study} vs Peters et al., {row.tissue}",
                          "total_variation_distance_percent_points": tvd,
                          **{f"difference_{s}": row[s] - ref[s] for s in STATES}})
    agreement = pd.DataFrame(agreement)
    agreement.to_csv(stage / "agreement.csv", index=False)
    write_json(stage / "crossstudy_definition.json", {
        "reference": str(PETERS), "reference_denominator": "CD3e+ T cells positive for TCF1 and/or TOX (Figure 7)",
        "this_pipeline": {"pancreas_T": "hierarchical identity: broad Immune kept by refinement, immune_class Imm-T-cell; "
                                        "and the producer's Imm-T-cell definition (CD3e >= 0.10, epithelium < 0.10)",
                          "pLN_T": "the producer's Imm-T-cell definition on 6551's same-block node (aligned with the pancreas)",
                          "positivity": "TCF_1 / TOX >= 0.10 on the cohort scale", "distance": "signed distance to the "
                          "union of the original true-Islet outlines (negative inside)"},
        "counts": counts, "use": "audit context only; never a tuning, gating or release target"})
    return {"counts": counts, "agreement": agreement.round(2).to_dict("records")}
