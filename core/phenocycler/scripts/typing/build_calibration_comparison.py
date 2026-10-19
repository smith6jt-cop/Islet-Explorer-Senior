#!/usr/bin/env python3
"""Build the consolidated, paired eight-arm typing calibration review.

The first listed arm is only a stable diff anchor. It is not a recommended production policy.
Every comparison is donor paired and uses the same output object order, eligibility flags, and
geometry; a mismatch is a hard error rather than an implicit inner join.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(HERE))
import subpanel_io as S  # noqa: E402
import subpanel as SP  # noqa: E402


EXPECTED_ARMS = [
    "redsea_shifted_frozen", "redsea_unshifted_frozen",
    "raw_shifted_frozen", "raw_unshifted_frozen",
    "redsea_shifted_recomputed", "redsea_unshifted_recomputed",
    "raw_shifted_recomputed", "raw_unshifted_recomputed",
]
CELL_FIELDS = [
    "object_id", "qc_analysis_eligible", "cell_region", "X_centroid", "Y_centroid",
    "signed_dist_islet_um", "islet_num", "endocrine_object_class", "endocrine_object_id",
    "cell_class", "refined_class", "refined_phenotype", "phenotype", "is_endocrine",
    "endocrine_route", "resolved_by", "immune_leaf", "top_identity_marker",
    "top_expressed_marker", "argmax_scaled_no_supporting_only", "argmax_within_marker_rank",
    "islet_route", "islet_stratum", "islet_cells",
]
GEOMETRY_FIELDS = [
    "qc_analysis_eligible", "cell_region", "X_centroid", "Y_centroid",
    "signed_dist_islet_um", "islet_num", "endocrine_object_class", "endocrine_object_id",
]
CALL_FIELDS = [
    "cell_class", "refined_class", "refined_phenotype", "is_endocrine",
    "endocrine_route", "resolved_by", "immune_leaf", "top_identity_marker",
    "top_expressed_marker", "argmax_scaled_no_supporting_only",
    "argmax_within_marker_rank", "islet_route",
]
IMMUNE_SUBTYPE_FIELDS = {
    "immune_class", "immune_subtype", "immune_class_undecided", "immune_subtype_undecided"
}


def read_frame(path: Path, columns: list[str]) -> pd.DataFrame:
    pf = pq.ParquetFile(path)
    have = set(pf.schema_arrow.names)
    missing = sorted(set(columns) - have)
    if missing:
        raise ValueError(f"{path}: missing expected columns {missing}")
    return pq.read_table(path, columns=columns).to_pandas()


def norm_label(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=object)
    missing = pd.isna(arr)
    if not missing.any():
        return arr
    out = arr.copy()
    out[missing] = "<NA>"
    return out


def same_values(left: np.ndarray, right: np.ndarray, name: str) -> bool:
    if left.dtype.kind in "fiu" or right.dtype.kind in "fiu":
        return bool(np.array_equal(left, right, equal_nan=True))
    return bool(np.array_equal(norm_label(left), norm_label(right)))


def changed_count(left: np.ndarray, right: np.ndarray, eligible: np.ndarray) -> int:
    if left.dtype.kind in "fiu" or right.dtype.kind in "fiu":
        equal = np.isclose(left, right, rtol=0, atol=0, equal_nan=True)
    else:
        equal = norm_label(left) == norm_label(right)
    return int((eligible & ~equal).sum())


def transition_rows(donor: str, arm: str, field: str, old: np.ndarray,
                    new: np.ndarray, eligible: np.ndarray) -> list[dict]:
    old_s, new_s = norm_label(old), norm_label(new)
    changed = eligible & (old_s != new_s)
    if not changed.any():
        return []
    pairs = pd.MultiIndex.from_arrays([old_s[changed], new_s[changed]],
                                      names=["from", "to"])
    counts = pd.Series(1, index=pairs, dtype=np.int64).groupby(
        level=["from", "to"], sort=True).sum()
    return [{"donor_id": donor, "comparison_arm": arm, "field": field,
             "from": a, "to": b, "n_cells": int(n)}
            for (a, b), n in counts.items()]


def factorial_pairs() -> list[tuple[str, str, str]]:
    """All matched, single-factor contrasts, returned as (effect, arm_a, arm_b)."""
    out = []
    for cd8 in ("shifted", "unshifted"):
        for scale in ("frozen", "recomputed"):
            out.append(("REDSEA_vs_all_raw", f"redsea_{cd8}_{scale}", f"raw_{cd8}_{scale}"))
    for source in ("redsea", "raw"):
        for scale in ("frozen", "recomputed"):
            out.append(("shifted_vs_unshifted_CD8",
                        f"{source}_shifted_{scale}", f"{source}_unshifted_{scale}"))
    for source in ("redsea", "raw"):
        for cd8 in ("shifted", "unshifted"):
            out.append(("frozen_vs_recomputed_scale",
                        f"{source}_{cd8}_frozen", f"{source}_{cd8}_recomputed"))
    return out


def json_audits(arm_root: Path) -> dict[str, dict]:
    value = json.loads((arm_root / "subpanel_audit.json").read_text())
    return {str(row["donor_id"]): row for row in value}


def parse_reference_levels(log_path: Path) -> dict[str, dict[str, float]]:
    """Read the producer's printed, donor-specific gate levels (4-decimal audit precision)."""
    lines = log_path.read_text(errors="replace").splitlines()
    try:
        start = next(i for i, line in enumerate(lines)
                     if "PER-DONOR reference levels" in line)
    except StopIteration:
        raise ValueError(f"{log_path}: per-donor reference-level table is missing")
    columns = ("hormone", "CD99", *SP.STRUCTURES)
    out = {}
    for line in lines[start + 1:]:
        fields = line.split()
        if not fields:
            if out:
                break
            continue
        if len(fields) == 6 and re.fullmatch(r"\d+", fields[0]):
            out[fields[0]] = {name: float(value) for name, value in zip(columns, fields[1:])}
    if not out:
        raise ValueError(f"{log_path}: no donor reference levels parsed")
    return out


def structure_test_rows(donor: str, arm: str, audit: dict,
                        levels: dict[str, float], log: dict) -> list[dict]:
    rows = []
    fields = ("inside", "inside_and_majority", "inside_area_unavailable",
              "majority_withholds", "inside_and_clears", "called",
              "top_marker_withholds", "top_marker_withholds_in_islet")
    for structure, values in audit.get("structures", {}).items():
        row = {"donor_id": donor, "arm": arm, "structure": structure}
        row.update({k: values.get(k) for k in fields})
        row["marker_channels"] = ";".join(log.get("subpanels", {}).get(structure, []))
        row["ordinary_tissue_reference_level_p99_printed_4dp"] = levels.get(structure)
        rows.append(row)
    return rows


def hormone_profile(donor: str, arm: str, arm_root: Path, run: dict) -> tuple[list[dict], dict]:
    """Use the exact fragment mask, scale transform, hormone evidence and single-cell rule."""
    path = arm_root / "subpanel_rebuilt" / f"donor_id={donor}" / "data.parquet"
    calls = read_frame(path, ["object_id", "cell_region", "endocrine_object_class", "cell_class"])
    core = calls.loc[calls.cell_region.eq("core")].reset_index(drop=True)
    if core.empty:
        return [], {"object_id": np.array([], dtype=object), "argmax": np.array([], dtype=object),
                    "margin": np.array([], dtype=float)}

    S.ALIGNED = Path(run["aligned_root"])
    S.SCALE_CSV = Path(run["scale"])
    scale = pd.read_parquet(S.SCALE_CSV).set_index("marker")
    values = SP.scaled_markers(donor, set(S.HORMONES), scale, core.object_id)
    evidence = SP.hormone_evidence(values)
    names = {v: k for k, v in json.loads((REPO / "data/development/typing/subpanel/subpanel_log.json").read_text())["hormones"].items()}
    single_call, single_status = SP.direct_single_identity(values, names)
    top = np.asarray(evidence["top"], dtype=float)
    runner = np.asarray(evidence["runner_up"], dtype=float)
    margin = top - runner
    population = core["endocrine_object_class"].fillna("all_core")
    rows = []
    for pop in ["all_core", "EndocrineSingle", "EndocrineCluster", "Islet"]:
        mask = (np.ones(len(core), dtype=bool) if pop == "all_core"
                else population.eq(pop).to_numpy())
        if not mask.any():
            continue
        vals = margin[mask & np.isfinite(margin)]
        usable = np.asarray(evidence["measured"])[mask]
        complete = np.asarray(evidence["complete"])[mask]
        base = {"donor_id": donor, "arm": arm, "population": pop,
                "n_cells": int(mask.sum()), "n_with_any_hormone": int(usable.sum()),
                "n_with_all_hormones": int(complete.sum()),
                "median_top_scaled": float(np.nanmedian(top[mask])) if np.isfinite(top[mask]).any() else np.nan,
                "median_runner_up_scaled": float(np.nanmedian(runner[mask])) if np.isfinite(runner[mask]).any() else np.nan,
                "median_top_minus_runner_up": float(np.nanmedian(vals)) if len(vals) else np.nan,
                "q25_top_minus_runner_up": float(np.nanquantile(vals, .25)) if len(vals) else np.nan,
                "q75_top_minus_runner_up": float(np.nanquantile(vals, .75)) if len(vals) else np.nan}
        for hormone in S.HORMONES:
            base[f"argmax_{hormone}"] = int(
                ((np.asarray(evidence["which"], dtype=object) == hormone) & mask
                 & np.asarray(evidence["measured"])).sum())
        if pop == "EndocrineSingle":
            for status in ("unique_highest_hormone", "tied_highest_hormones", "no_usable_hormone"):
                base[f"single_{status}"] = int(((single_status == status) & mask).sum())
            for label in ("Beta", "Alpha", "Delta", "Endocrine (ambiguous subtype)",
                          "Endocrine (unresolved subtype)"):
                base[f"single_assignment_{label}"] = int(((single_call == label) & mask).sum())
        rows.append(base)
    status = np.asarray(evidence["measured"], dtype=bool)
    argmax = np.asarray(evidence["which"], dtype=object)
    argmax[~status] = "<UNAVAILABLE>"
    return rows, {"object_id": core.object_id.to_numpy(dtype=object), "argmax": argmax,
                  "margin": margin, "single_mask": population.eq("EndocrineSingle").to_numpy(),
                  "single_status": np.asarray(single_status, dtype=object),
                  "single_call": np.asarray(single_call, dtype=object)}


def immune_calls(arm_root: Path, donor: str) -> pd.DataFrame:
    path = arm_root / "immune" / f"donor_id={donor}" / "data.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"immune subtype producer output is missing: {path}")
    return read_frame(path, ["object_id", "immune_class", "immune_subtype",
                             "immune_class_undecided", "immune_subtype_undecided"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms-root", type=Path, default=Path("/data/io60-calibration/arms"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--input-manifest", type=Path,
                    default=Path("results/v3/release_20261003/current_input_manifest.json"))
    ap.add_argument("--alignment-manifest", type=Path,
                    default=Path("/data/io60-calibration/alignment_arms/alignment_arms.json"))
    ap.add_argument("--anchor", default="redsea_shifted_frozen")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="development only; normal review requires all eight complete arms")
    args = ap.parse_args()
    arms_root, out = args.arms_root.resolve(), args.out.resolve()
    names = EXPECTED_ARMS
    if args.anchor not in names:
        raise ValueError(f"anchor must be one of {names}")
    runs, audits = {}, {}
    for name in names:
        arm = arms_root / name
        run = json.loads((arm / "run.json").read_text())
        if run.get("status") != "complete" and not args.allow_incomplete:
            raise RuntimeError(f"{name} status is {run.get('status')!r}; comparison awaits full typing")
        runs[name] = run
        audits[name] = json_audits(arm)
    donors = list(runs[args.anchor]["donors"])
    if "6457" in donors or len(donors) != 25:
        raise ValueError(f"unexpected donor universe {donors}")
    if any(set(run["donors"]) != set(donors) for run in runs.values()):
        raise ValueError("calibration arms do not have identical donor sets")
    out.mkdir(parents=True, exist_ok=True)

    pair_rows, transition_data, marker_rows, hormone_rows = [], [], [], []
    reference_level_rows, endocrine_level_rows = [], []
    immune_rows, distribution_rows, example_rows, validation = [], [], [], []
    all_pairs = [("anchor_vs_arm", args.anchor, other) for other in names if other != args.anchor]
    all_pairs += factorial_pairs()
    typing_log = json.loads((REPO / "data/development/typing/subpanel/subpanel_log.json").read_text())
    reference_levels = {name: parse_reference_levels(arms_root / name / "subpanel.log")
                        for name in names}

    for donor in donors:
        frames = {}
        immune = {}
        profiles = {}
        for arm in names:
            root = arms_root / arm
            frame = read_frame(root / "cells" / f"donor_id={donor}" / "data.parquet", CELL_FIELDS)
            if frame.object_id.duplicated().any():
                raise ValueError(f"{arm}/{donor}: duplicate object_id in cell release")
            frames[arm] = frame
            immune[arm] = immune_calls(root, donor)
            hormone_rows_arm, profiles[arm] = hormone_profile(donor, arm, root, runs[arm])
            hormone_rows.extend(hormone_rows_arm)
            marker_rows.extend(structure_test_rows(donor, arm, audits[arm][donor],
                                                   reference_levels[arm][donor], typing_log))
            for marker, value in reference_levels[arm][donor].items():
                reference_level_rows.append({"donor_id": donor, "arm": arm,
                    "marker_or_gate": marker, "reference_level_p99_printed_4dp": value,
                    "precision_note": "subpanel.log producer display precision is 4 decimals"})
            audit = audits[arm][donor]
            for stratum, entry in audit.get("islet_levels", {}).items():
                if not isinstance(entry, dict):
                    continue
                unguarded = audit.get("islet_levels_unguarded", {}).get(stratum, {})
                endocrine_level_rows.append({"donor_id": donor, "arm": arm, "stratum": stratum,
                    "n_hormone_confirmed_reference_cells": entry.get("n_confirmed"),
                    "cd99_level": entry.get("CD99"),
                    "nucleus_elongation_level": entry.get("nucleus_elongation"),
                    "cd99_level_unguarded": unguarded.get("CD99"),
                    "nucleus_elongation_level_unguarded": unguarded.get("nucleus_elongation")})

            validation.append({"donor_id": donor, "arm": arm,
                "cell_rows": len(frame), "qc_analysis_eligible": int(frame.qc_analysis_eligible.sum()),
                "audit_n_cells": int(audit["n_cells"]),
                "primary_rows_match_audit": len(frame) == int(audit["n_cells"]),
                "complete_arm": runs[arm].get("status") == "complete"})
            for field in ("cell_class", "refined_class", "refined_phenotype"):
                values = norm_label(frame[field])
                counts = collections.Counter(values[frame.qc_analysis_eligible.to_numpy(dtype=bool)])
                for label, count in sorted(counts.items()):
                    distribution_rows.append({"donor_id": donor, "arm": arm,
                                              "field": field, "label": label, "n_cells": count})

        anchor_frame = frames[args.anchor]
        anchor_ids = anchor_frame.object_id.to_numpy(dtype=object)
        anchor_elig = anchor_frame.qc_analysis_eligible.to_numpy(dtype=bool)
        if len(np.unique(anchor_ids)) != len(anchor_ids):
            raise ValueError(f"{donor}: duplicate IDs in anchor")
        for arm in names:
            frame = frames[arm]
            if not np.array_equal(frame.object_id.to_numpy(dtype=object), anchor_ids):
                raise ValueError(f"{donor}/{arm}: object IDs or row order differ from the anchor")
            for field in GEOMETRY_FIELDS:
                if not same_values(anchor_frame[field].to_numpy(), frame[field].to_numpy(), field):
                    raise ValueError(f"{donor}/{arm}: geometry/eligibility changed in {field}")

        # Map sparse immune output to the common object universe. Missing call rows are explicit
        # unresolved values; they are never counted as a resolved phenotype.
        immune_maps = {}
        for arm, table in immune.items():
            if table.object_id.duplicated().any():
                raise ValueError(f"{donor}/{arm}: duplicate object_id in immune subtyping")
            table = table.set_index("object_id")
            immune_maps[arm] = {
                field: table[field].reindex(anchor_ids).to_numpy(dtype=object)
                for field in ("immune_class", "immune_subtype")}
            for field in ("immune_class_undecided", "immune_subtype_undecided"):
                immune_maps[arm][field] = table[field].reindex(anchor_ids).fillna(False).to_numpy(dtype=bool)

        # Baseline comparison is merely a readable diff. The factorial rows below are the
        # controlled scientific contrasts that separate the three declared choices.
        for effect, left_name, right_name in all_pairs:
            left, right = frames[left_name], frames[right_name]
            if not np.array_equal(left.object_id.to_numpy(dtype=object),
                                  right.object_id.to_numpy(dtype=object)):
                raise ValueError(f"{donor}: {left_name} and {right_name} object universes differ")
            elig = (left.qc_analysis_eligible.to_numpy(dtype=bool)
                    & right.qc_analysis_eligible.to_numpy(dtype=bool))
            fields = list(CALL_FIELDS)
            fields += ["immune_class", "immune_subtype", "immune_class_undecided",
                       "immune_subtype_undecided"]
            for field in fields:
                a = immune_maps[left_name][field] if field in IMMUNE_SUBTYPE_FIELDS else left[field].to_numpy()
                b = immune_maps[right_name][field] if field in IMMUNE_SUBTYPE_FIELDS else right[field].to_numpy()
                n_changed = changed_count(a, b, elig)
                pair_rows.append({"donor_id": donor, "contrast": effect,
                                  "arm_a": left_name, "arm_b": right_name, "field": field,
                                  "eligible_cells": int(elig.sum()), "changed_cells": n_changed,
                                  "changed_fraction": n_changed / max(int(elig.sum()), 1)})
                if effect == "anchor_vs_arm":
                    transition_data.extend(transition_rows(donor, right_name, field, a, b, elig))

            h_a, h_b = profiles[left_name], profiles[right_name]
            if not np.array_equal(h_a["object_id"], h_b["object_id"]):
                raise ValueError(f"{donor}: hormone evidence universe differs in {left_name}/{right_name}")
            # Hormone evidence is defined over the same geometric core subset; IDs are identical
            # across arms, and this mask corresponds to the core-only profile rows.
            h_change = int((h_a["argmax"] != h_b["argmax"]).sum())
            pair_rows.append({"donor_id": donor, "contrast": effect, "arm_a": left_name,
                              "arm_b": right_name, "field": "hormone_argmax_core",
                              "eligible_cells": len(h_a["argmax"]), "changed_cells": h_change,
                              "changed_fraction": h_change / max(len(h_a["argmax"]), 1)})
            single = h_a["single_mask"] & h_b["single_mask"]
            s_change = int((norm_label(h_a["single_status"])[single]
                            != norm_label(h_b["single_status"])[single]).sum())
            pair_rows.append({"donor_id": donor, "contrast": effect, "arm_a": left_name,
                              "arm_b": right_name, "field": "endocrine_single_hormone_status",
                              "eligible_cells": int(single.sum()), "changed_cells": s_change,
                              "changed_fraction": s_change / max(int(single.sum()), 1)})

            if effect == "anchor_vs_arm":
                primary = right
                diff = (norm_label(left.cell_class) != norm_label(right.cell_class)) & elig
                # Up to three stable coordinate examples per donor/arm for each common primary
                # transition. They support later QuPath crop assembly without inventing matches.
                old_class, new_class = norm_label(left.cell_class), norm_label(right.cell_class)
                changes = np.where(diff)[0]
                groups = collections.defaultdict(list)
                for i in changes:
                    groups[(old_class[i], new_class[i])].append(i)
                for (old, new), indices in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
                    ranked = sorted(indices, key=lambda i: hashlib.sha256(
                        f"{right_name}|{donor}|{anchor_ids[i]}".encode()).hexdigest())[:3]
                    for i in ranked:
                        example_rows.append({"donor_id": donor, "arm": right_name,
                            "from_cell_class": old, "to_cell_class": new,
                            "transition_n": len(indices), "object_id": anchor_ids[i],
                            "X_centroid": right.X_centroid.iloc[i],
                            "Y_centroid": right.Y_centroid.iloc[i],
                            "cell_region": right.cell_region.iloc[i],
                            "signed_dist_islet_um": right.signed_dist_islet_um.iloc[i],
                            "top_expressed_marker": right.top_expressed_marker.iloc[i],
                            "refined_class": right.refined_class.iloc[i]})

    pd.DataFrame(pair_rows).to_csv(out / "paired_call_changes.csv", index=False)
    pd.DataFrame(transition_data).to_csv(out / "anchor_class_transitions.csv", index=False)
    pd.DataFrame(marker_rows).to_csv(out / "structure_marker_tests.csv", index=False)
    pd.DataFrame(reference_level_rows).to_csv(out / "ordinary_tissue_reference_levels.csv", index=False)
    pd.DataFrame(endocrine_level_rows).to_csv(out / "islet_endocrine_levels.csv", index=False)
    pd.DataFrame(hormone_rows).to_csv(out / "hormone_argmax_margin.csv", index=False)
    pd.DataFrame(distribution_rows).to_csv(out / "typing_distributions.csv", index=False)
    pd.DataFrame(example_rows).to_csv(out / "candidate_image_examples.csv", index=False)
    pd.DataFrame(validation).to_csv(out / "calibration_validation.csv", index=False)
    input_manifest = json.loads(args.input_manifest.read_text())
    same_slide = [row for row in input_manifest.get("records", [])
                  if row.get("section_relationship") == "same_slide"]
    same_slide_fields = ["donor_id", "tissue", "slide_id", "modality", "physical_block_id",
                         "section_relationship_evidence", "local_source_available",
                         "local_source_path", "prepared_h5ad", "prepared_h5ad_source_kind",
                         "prepared_match_is_exact", "coverage_status"]
    pd.DataFrame([{field: row.get(field) for field in same_slide_fields}
                  for row in same_slide]).to_csv(out / "same_slide_rna_evidence.csv", index=False)
    alignment_manifest = json.loads(args.alignment_manifest.read_text())
    reference_donor_rows = []
    for source_name in ("redsea", "raw"):
        source_arm = alignment_manifest["arms"][f"{source_name}_shifted"]
        shift_path = Path(source_arm["shifts"])
        shift_data = json.loads(shift_path.read_text())
        if shift_data.get("uniform_commit") != source_arm.get("uniform_commit"):
            raise ValueError(f"{source_name}: UniFORM commit differs from alignment manifest")
        for marker, item in shift_data.get("shifts", {}).items():
            reference_donor_rows.append({"input_policy": source_name, "marker": marker,
                "reference_donor": item.get("reference_donor"),
                "reference_source": item.get("reference_source"),
                "n_donors": item.get("n_donors"),
                "shift_file": str(shift_path), "shift_file_sha256": source_arm["shifts_sha256"],
                "uniform_commit": source_arm.get("uniform_commit")})
    pd.DataFrame(reference_donor_rows).to_csv(out / "uniform_reference_donors.csv", index=False)
    manifest = {
        "schema_version": 1,
        "anchor_arm": args.anchor,
        "anchor_is_diff_reference_only": True,
        "alignment_manifest": str(args.alignment_manifest.resolve()),
        "donors": donors,
        "arms": {name: {"status": runs[name].get("status"), "scale": runs[name]["scale"],
                          "scale_sha256": runs[name]["scale_sha256"],
                          "aligned_root": runs[name]["aligned_root"],
                          "criteria_sha256": runs[name]["registry_criteria_sha256"],
                          "stages": runs[name].get("stages", [])} for name in names},
        "geometry_checked_identical": True,
        "eligible_cell_universe_checked_identical": True,
        "hormone_values_fragment_masked": True,
        "primary_comparison": "paired within donor; no positivity-frequency matching or disease-order gate",
        "notes": [
            "Marker-test rows preserve the declared in-object, majority-area, fitted-level, and top-expressed-marker stages.",
            "Hormone margins use the producer scale and fragment mask; missing channels remain unavailable.",
            "same_slide_rna_evidence.csv lists every same-slide RNA record; no cell-level RNA check is run without exact local input and verified registration.",
            "The top_identity_marker field is retained as legacy output but is a per-marker-scaled argmax, not top expression.",
        ],
    }
    (out / "comparison_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote paired calibration dossier: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
