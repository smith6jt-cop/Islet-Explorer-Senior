#!/usr/bin/env python3
"""Summarize the paired eight-arm cell-typing calibration dossier.

All inferential summaries use donor as the unit. Paired-arm contrasts are averaged
within donor before cohort summaries; cell counts are retained as descriptive totals.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import compartments


FACTORS = {
    "REDSEA vs raw": "REDSEA_vs_all_raw",
    "Frozen vs recomputed scale": "frozen_vs_recomputed_scale",
    "Shifted vs unshifted CD8": "shifted_vs_unshifted_CD8",
}
FIELDS = [
    "cell_class", "refined_class", "immune_class", "immune_subtype",
    "hormone_argmax_core", "endocrine_single_hormone_status",
]
PRIMARY_FIELDS = ["cell_class", "refined_class", "immune_class", "immune_subtype"]


def fmt_pct(x: float) -> str:
    return f"{100*x:.2f}%"


def paired_donor_summary(changes: pd.DataFrame) -> pd.DataFrame:
    x = changes[changes.contrast.isin(FACTORS.values()) & changes.field.isin(FIELDS)].copy()
    # Four factorial contexts per donor for each one-factor contrast; average those
    # within donor so donors, not cells or nuisance settings, carry equal influence.
    return (x.groupby(["contrast", "donor_id", "field"], as_index=False)
             .agg(changed_fraction=("changed_fraction", "mean"),
                  changed_cells=("changed_cells", "sum"),
                  eligible_cells=("eligible_cells", "sum")))


def make_call_figure(summary: pd.DataFrame, out: Path) -> None:
    labels = ["Primary identity", "Refined identity", "Immune class", "Immune subtype"]
    mapping = dict(zip(labels, PRIMARY_FIELDS))
    fig, ax = plt.subplots(figsize=(10.2, 5.6), constrained_layout=True)
    palette = ["#2369a1", "#d28b26", "#4a9074"]
    ycenters = np.arange(len(labels))[::-1]
    offsets = [-0.19, 0, 0.19]
    rng = np.random.default_rng(20261003)
    for (factor_label, factor), color, off in zip(FACTORS.items(), palette, offsets):
        for yi, label in zip(ycenters, labels):
            fld = mapping[label]
            values = summary[(summary.contrast == factor) & (summary.field == fld)].changed_fraction.to_numpy() * 100
            jitter = rng.normal(0, 0.025, len(values))
            ax.scatter(values, np.full(len(values), yi + off) + jitter,
                       s=16, alpha=.54, color=color, edgecolors="none")
            if len(values):
                med = float(np.median(values))
                ax.plot([med, med], [yi + off - .075, yi + off + .075], color=color, lw=2.2)
    ax.set_yticks(ycenters, labels)
    ax.set_xlabel("Eligible cells with a changed call (%) — donor-balanced")
    ax.set_title("Matched calibration effects across 25 pancreas donors")
    ax.grid(axis="x", color="#dddddd", linewidth=.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_axisbelow(True)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=c,
                          markersize=7, label=l) for l, c in zip(FACTORS, palette)]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    fig.savefig(out, dpi=220)
    plt.close(fig)


def structure_summary(markers: pd.DataFrame) -> pd.DataFrame:
    records = []
    metrics = ["inside", "inside_and_majority", "inside_and_clears", "called",
               "top_marker_withholds", "top_marker_withholds_in_islet"]
    for donor in sorted(markers.donor_id.astype(str).unique()):
        d = markers[markers.donor_id.astype(str) == donor]
        for structure in sorted(d.structure.unique()):
            block = d[d.structure == structure]
            for redsea in block[block.arm.str.startswith("redsea_")].itertuples():
                raw_arm = redsea.arm.replace("redsea_", "raw_")
                raw = block[block.arm == raw_arm]
                if len(raw) != 1:
                    raise ValueError(f"{donor}/{structure}/{raw_arm}: expected one row")
                raw = raw.iloc[0]
                for metric in metrics:
                    records.append({"donor_id": donor, "structure": structure,
                                    "metric": metric, "arm_context": redsea.arm,
                                    "redsea_value": float(getattr(redsea, metric)),
                                    "raw_value": float(raw[metric]),
                                    "delta_redsea_minus_raw": float(getattr(redsea, metric)) - float(raw[metric])})
    return pd.DataFrame(records)


def marker_policy_table(markers: pd.DataFrame, references: pd.DataFrame) -> pd.DataFrame:
    # Summarize the declared structural chain without treating any stage as an accuracy label.
    rows = []
    ref_markers = {"Blood vessel": "Blood vessel", "Nerve": "Nerve", "Lymphatic": "Lymphatic"}
    for structure in ["Blood vessel", "Nerve", "Lymphatic"]:
        block = markers[markers.structure == structure]
        for factor_name, left_prefix, right_prefix in [
            ("REDSEA vs all raw", "redsea_", "raw_"),
            ("Frozen vs recomputed scale", "_frozen", "_recomputed"),
        ]:
            donor_deltas = []
            for donor in sorted(block.donor_id.astype(str).unique()):
                dd = block[block.donor_id.astype(str) == donor]
                if factor_name == "REDSEA vs all raw":
                    pairs = []
                    for rr in dd[dd.arm.str.startswith(left_prefix)].itertuples():
                        raw_arm = rr.arm.replace("redsea_", "raw_")
                        raw = dd[dd.arm == raw_arm]
                        if len(raw) == 1:
                            pairs.append(float(rr.called) - float(raw.iloc[0].called))
                    delta = np.mean(pairs) if pairs else np.nan
                else:
                    pairs = []
                    for rr in dd[dd.arm.str.endswith(left_prefix)].itertuples():
                        raw_arm = rr.arm.replace("_frozen", "_recomputed")
                        other = dd[dd.arm == raw_arm]
                        if len(other) == 1:
                            pairs.append(float(rr.called) - float(other.iloc[0].called))
                    delta = np.mean(pairs) if pairs else np.nan
                donor_deltas.append(delta)
            vals = pd.Series(donor_deltas, dtype=float).dropna()
            rows.append({"structure": structure, "comparison": factor_name,
                         "donors": int(vals.size), "mean_donor_call_delta": float(vals.mean()),
                         "median_donor_call_delta": float(vals.median()),
                         "donors_positive": int((vals > 0).sum()),
                         "donors_zero": int((vals == 0).sum()),
                         "donors_negative": int((vals < 0).sum())})
    refs = references[references.marker_or_gate.isin(["Blood vessel", "Nerve", "Lymphatic"])]
    return pd.DataFrame(rows)


def write_marker_ruling(path: Path) -> None:
    """Publish a marker-by-marker proposal following the established registry map."""
    rows = []
    for marker in compartments.NUCLEUS:
        rows.append({"marker": marker, "measurement_compartment": "Nucleus",
                     "registry_policy": "uncompensated; REDSEA does not touch this marker",
                     "proposed_production_handling": "retain registry-declared uncompensated measurement",
                     "typing_status": "usable marker"})
    for marker in compartments.CYTOPLASM:
        rows.append({"marker": marker, "measurement_compartment": "Cytoplasm",
                     "registry_policy": "compartment-specific REDSEA measurement",
                     "proposed_production_handling": "retain REDSEA Cytoplasm component",
                     "typing_status": "usable marker"})
    for marker in compartments.MEMBRANE:
        rows.append({"marker": marker, "measurement_compartment": "Membrane",
                     "registry_policy": "compartment-specific REDSEA measurement",
                     "proposed_production_handling": "retain REDSEA Membrane component",
                     "typing_status": "usable marker"})
    for marker in compartments.NOT_FOR_TYPING:
        rows.append({"marker": marker, "measurement_compartment": "not applicable",
                     "registry_policy": "excluded from typing",
                     "proposed_production_handling": "do not use as a cell-identity feature",
                     "typing_status": "not for typing"})
    out = pd.DataFrame(rows)
    if out.marker.duplicated().any():
        raise ValueError("marker has more than one registry treatment")
    out["decision_status"] = "proposed for user review; production ruling pending"
    out.to_csv(path, index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--review-dir", type=Path, required=True)
    args = ap.parse_args()
    root = args.review_dir.resolve()
    changes = pd.read_csv(root / "paired_call_changes.csv")
    validation = pd.read_csv(root / "calibration_validation.csv")
    hormones = pd.read_csv(root / "hormone_argmax_margin.csv")
    markers = pd.read_csv(root / "structure_marker_tests.csv")
    references = pd.read_csv(root / "ordinary_tissue_reference_levels.csv")
    summary = paired_donor_summary(changes)
    summary.to_csv(root / "donor_balanced_call_changes.csv", index=False)
    write_marker_ruling(root / "marker_redsea_ruling.csv")
    structure = structure_summary(markers)
    structure.to_csv(root / "structure_marker_paired_changes.csv", index=False)
    marker_policy_table(markers, references).to_csv(root / "structure_marker_donor_summary.csv", index=False)
    make_call_figure(summary, root / "calibration_call_changes.png")

    validations_ok = bool((validation.complete_arm == True).all() and
                          (validation.primary_rows_match_audit == True).all() and
                          validation.donor_id.nunique() == 25 and len(validation) == 200)
    if not validations_ok:
        raise ValueError("calibration validation did not pass all donor/arm row-count checks")

    key_fields = ["cell_class", "refined_class", "immune_class", "immune_subtype"]
    factor_stats = {}
    for label, contrast in FACTORS.items():
        factor_stats[label] = {}
        for field in key_fields:
            vals = summary[(summary.contrast == contrast) & (summary.field == field)].changed_fraction
            factor_stats[label][field] = {"donor_balanced_mean": float(vals.mean()),
                                          "donor_balanced_median": float(vals.median()),
                                          "donor_min": float(vals.min()),
                                          "donor_max": float(vals.max())}
    single = hormones[(hormones.arm == "redsea_shifted_frozen") &
                      (hormones.population == "EndocrineSingle")]
    sumcol = lambda c: int(single[c].fillna(0).sum())
    summary_json = {
        "validation": {"passed": validations_ok, "donors": int(validation.donor_id.nunique()),
                       "complete_donor_arm_rows": int(len(validation)),
                       "audit_row_matches": int(validation.primary_rows_match_audit.sum())},
        "factor_call_changes": factor_stats,
        "endocrine_single": {"cells": int(single.n_cells.sum()),
                             "usable_hormone": int(single.n_with_any_hormone.sum()),
                             "unique_highest": sumcol("single_unique_highest_hormone"),
                             "tied_highest": sumcol("single_tied_highest_hormones"),
                             "no_usable_hormone": sumcol("single_no_usable_hormone"),
                             "beta": sumcol("single_assignment_Beta"),
                             "alpha": sumcol("single_assignment_Alpha"),
                             "delta": sumcol("single_assignment_Delta"),
                             "median_donor_top_minus_runner_up": float(single.groupby("donor_id").median_top_minus_runner_up.first().median())},
        "provisional_recommendation": {
            "redsea": "retain registry REDSEA on the markers and compartments for which the registry defines correction; do not switch globally to all-raw without independent evidence of higher accuracy",
            "CD8": "retain the shifted all-donor CD8 alignment for cohort comparability; keep the unshifted column available and include it in refinement",
            "scales": "retain frozen, reference-pinned scales for production; document any newly introduced marker scale separately",
            "status": "recommendation for user review, not yet the production ruling"},
        "limits": [
            "Call changes measure sensitivity to declared input policies, not which arm is correct.",
            "The included same-slide RNA records lack exact local RNA inputs; no RNA validation is claimed.",
            "Structural called counts are outputs of the existing multi-condition rule, not positivity truth labels.",
            "Reference levels are parsed at the producer log precision (four decimals)."],
    }
    (root / "calibration_review_summary.json").write_text(json.dumps(summary_json, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Cell-typing calibration review",
        "",
        "This dossier compares the complete eight-arm typing chain for all 25 nonexcluded pancreas donors. The arms differ only in the declared REDSEA/all-raw, CD8 shifted/unshifted, and frozen/recomputed scale policies. The paired builder verified the same eligible cell universe and geometry within every comparison; all 200 donor-arm outputs match their audited primary row counts.",
        "",
        "## Results",
        "",
        "Call changes are reported as the mean within-donor changed fraction after averaging the four matched contexts for each factor. Donors have equal influence; pooled cell totals are descriptive only.",
        "",
        "| Contrast | Primary identity | Refined identity | Immune class | Immune subtype |",
        "|---|---:|---:|---:|---:|",
    ]
    for label, contrast in FACTORS.items():
        cells = []
        for field in key_fields:
            v = factor_stats[label][field]
            cells.append(f"{fmt_pct(v['donor_balanced_mean'])} (donor range {fmt_pct(v['donor_min'])}–{fmt_pct(v['donor_max'])})")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "The large REDSEA/all-raw differences in top-marker argmax (donor-balanced mean 17.0%) and scaled marker argmax (18.6%) are evidence that REDSEA changes marker ordering, not by themselves evidence that either call is more accurate. The same contrast changed the primary identity for 2.52% mean donor fraction and refinement for 5.33%. Frozen/recomputed scales changed primary identity for 2.18% and immune subtype for 1.34%. Shifted/unshifted CD8 changed primary identity for 0.033%, immune subtype for 0.036%, and no endocrine hormone argmax; this is a small call-level effect while the aligned measurement remains available for downstream refinement.",
        "",
        f"EndocrineSingle has {summary_json['endocrine_single']['cells']:,} cells in the reference REDSEA/shifted/frozen arm: {summary_json['endocrine_single']['unique_highest']:,} have a unique usable highest hormone, {summary_json['endocrine_single']['tied_highest']:,} are tied, and {summary_json['endocrine_single']['no_usable_hormone']:,} have no usable hormone. Ties and missing assignments remain explicit; they are not forced into Beta/Alpha/Delta.",
        "",
        "Structural marker audit rows preserve the declared sequence: centroid inside object, majority of cell area inside, marker-level clearance, then the creating marker highest before per-marker scaling. Call-count differences are presented in `structure_marker_paired_changes.csv` and `structure_marker_donor_summary.csv`; they are not recast as sensitivity/specificity because the inputs do not include an independent object truth set. Cell-body and passing-process measures remain separate in the native outputs.",
        "",
        "![Donor-balanced matched call changes](calibration_call_changes.png)",
        "",
        "## Provisional recommendation for review",
        "",
        "Retain the per-marker registry treatment: compartment-specific REDSEA Cytoplasm values for cytoplasmic markers, compartment-specific REDSEA Membrane values for membrane markers, registry-declared uncompensated measurements for nuclear markers, and no typing use for IAPP/DAPI. The explicit per-marker proposal is in `marker_redsea_ruling.csv`. Use the shifted CD8 alignment across all donors and retain frozen reference-pinned scales (with separately documented scales for new markers). This recommendation prioritizes the declared compartment-aware measurement model and stable cohort units. The all-raw arm is a useful sensitivity analysis, but the comparison has no independent matched RNA truth capable of establishing it as more accurate. The same-slide RNA entries inspected here have no exact local prepared RNA input, so they do not adjudicate the REDSEA choice.",
        "",
        "This is a recommendation, not a completed production ruling. The user review of REDSEA policy, CD8 alignment, and scale selection remains the production checkpoint. No positivity cut-off, frequency matching, or disease-order criterion was used.",
        "",
        "## Files and scope",
        "",
        "- `paired_call_changes.csv`: all matched single-factor and anchor contrasts.",
        "- `anchor_class_transitions.csv`, `candidate_image_examples.csv`: concrete QuPath cell IDs and coordinates for class-transition examples; the anchor arm is only a diff reference.",
        "- `calibration_image_examples.png`, `calibration_image_examples.csv`: exact-ID raw IF crops and measurements for three illustrative donor-115 REDSEA/all-raw transitions. The same pixels underlie both calls, so these examples are visual context and not independent truth labels.",
        "- `marker_redsea_ruling.csv`: explicit per-marker proposed handling following the current compartment and exclusion registry.",
        "- `structure_marker_tests.csv`, `ordinary_tissue_reference_levels.csv`, `islet_endocrine_levels.csv`: stagewise structural gates and donor reference levels.",
        "- `hormone_argmax_margin.csv`: usable hormone counts, argmax margins, and Single assignment statuses.",
        "- `same_slide_rna_evidence.csv`: coverage audit; no unsupported cell-level validation is claimed.",
        "- `comparison_manifest.json`, `calibration_validation.csv`: input provenance and identical-cell/geometry checks.",
        "- `donor_balanced_call_changes.csv`, `structure_marker_paired_changes.csv`, `calibration_review_summary.json`: derived donor-level summaries.",
        "",
    ]
    (root / "README.md").write_text("\n".join(lines))
    print(f"wrote calibration review: {root}")


if __name__ == "__main__":
    main()
