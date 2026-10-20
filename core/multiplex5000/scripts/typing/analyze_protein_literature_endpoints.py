#!/usr/bin/env python3
"""Donor-level Peters and Medina-Serpas protein endpoints from native tissue calls.

Pancreas measurements use the explicitly selected candidate typing arm. Dedicated pLN values
are read from the section's own QuPath export and use its native pLN calls. Per-cell values are
reduced to donor-level summaries before cohort inference. No cross-tissue absolute intensity
comparison is made because the panels and tissue-specific calibrations differ.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


REPO = Path(__file__).resolve().parents[4]
PANCREAS_MARKERS = {
    "CD8": "CD8", "TCF1": "TCF_1", "TOX": "TOX", "Granzyme_B": "Granzyme_B",
    "LAG3": "LAG3", "PD_1": "PD_1", "CD39": "CD39", "IFNG": "IFNG",
    "ICOS": "ICOS", "IDO1": "IDO1", "HLA_DR": "HLA_DR", "HLA_A": "HLA_A",
    "VISTA": "VISTA", "PD_L1": "PD_L1",
}
PLN_MARKERS = {**PANCREAS_MARKERS, "PD_L1": "PDL1"}
PETERS_MARKERS = ("CD8", "TCF1", "TOX", "Granzyme_B", "LAG3", "PD_1", "CD39", "IFNG", "ICOS")
INFLAMMATORY_MARKERS = ("IFNG", "IDO1", "HLA_DR", "HLA_A", "VISTA", "PD_L1")
BOOTSTRAPS = 10_000
SEED = 20261003


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_values(values: pd.Series) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(float)

def auc_cd8_vs_other_t(cd8: np.ndarray, other: np.ndarray) -> float:
    x = cd8[np.isfinite(cd8)]
    y = other[np.isfinite(other)]
    if not len(x) or not len(y):
        return np.nan
    return float(mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").statistic /
                 (len(x) * len(y)))


def auc_immune_vs_other(immune: np.ndarray, other: np.ndarray) -> float:
    x = immune[np.isfinite(immune)]
    y = other[np.isfinite(other)]
    if not len(x) or not len(y):
        return np.nan
    return float(mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").statistic /
                 (len(x) * len(y)))


def rho(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if len(x) < 2 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan
    value = spearmanr(x, y).statistic
    return float(value) if np.isfinite(value) else np.nan


def bootstrap_mean_ci(values: np.ndarray, seed: int) -> tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = rng.choice(x, size=(BOOTSTRAPS, len(x)), replace=True).mean(axis=1)
    return tuple(map(float, np.quantile(means, [0.025, 0.975])))


def one_sample_p(values: np.ndarray, reference: float) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)] - reference
    if not len(x) or np.allclose(x, 0):
        return 1.0 if len(x) else np.nan
    return float(wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)


def paired_p(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x) or np.allclose(x, 0):
        return 1.0 if len(x) else np.nan
    return float(wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)


def disease_map(input_manifest: dict) -> dict[str, str]:
    values: dict[str, set[str]] = {}
    for record in input_manifest.get("records", []):
        donor = str(record.get("donor_id", ""))
        status = str(record.get("disease_status", ""))
        if donor and status and status.lower() != "unavailable":
            values.setdefault(donor.lstrip("0") or "0", set()).add(status)
    for donor, statuses in values.items():
        if len(statuses) > 1:
            raise ValueError(f"{donor}: inconsistent donor disease labels {sorted(statuses)}")
    return {d: sorted(v)[0] for d, v in values.items()}


def disease_for(donor: str, diseases: dict[str, str]) -> str | None:
    return diseases.get(str(donor).lstrip("0") or "0")


def marker_summary(donor: str, disease: str, tissue: str, frame: pd.DataFrame,
                   t_mask: np.ndarray, cd8_mask: np.ndarray, other_t_mask: np.ndarray,
                   marker_name: str, marker_col: str, distance_col: str | None,
                   inflammatory: bool = False) -> tuple[dict, dict | None]:
    raw = finite_values(frame[marker_col]) if marker_col in frame else np.full(len(frame), np.nan)
    t_values = raw[t_mask]
    t_values = t_values[np.isfinite(t_values)]
    q = np.quantile(t_values, [0.25, 0.5, 0.75]) if len(t_values) else [np.nan] * 3
    auc = auc_cd8_vs_other_t(raw[cd8_mask], raw[other_t_mask])
    if distance_col and distance_col in frame:
        distance = finite_values(frame[distance_col])[t_mask]
        marker_t = raw[t_mask]
        distance_rho = rho(distance, marker_t)
        distance_n = int((np.isfinite(distance) & np.isfinite(marker_t)).sum())
    else:
        distance_rho, distance_n = np.nan, 0
    row = {
        "donor_id": donor, "disease_status": disease, "tissue": tissue,
        "marker": marker_name, "source_marker_column": marker_col,
        "n_qc_eligible_cells": int(frame.qc_analysis_eligible.fillna(False).sum()),
        "n_qc_eligible_cells_marker_measured": int(np.isfinite(raw[frame.qc_analysis_eligible.fillna(False).to_numpy(bool)]).sum()),
        "marker_available_in_qc_cells": bool(np.isfinite(raw[frame.qc_analysis_eligible.fillna(False).to_numpy(bool)]).any()),
        "n_t_cells": int(t_mask.sum()), "n_t_cells_marker_measured": int(np.isfinite(raw[t_mask]).sum()),
        "t_cell_marker_q25": float(q[0]), "t_cell_marker_median": float(q[1]),
        "t_cell_marker_q75": float(q[2]),
        "n_strict_cd8_calls_marker_measured": int(np.isfinite(raw[cd8_mask]).sum()),
        "n_known_noncd8_t_calls_marker_measured": int(np.isfinite(raw[other_t_mask]).sum()),
        "marker_auc_strict_cd8_vs_known_noncd8_t": auc,
        "distance_spearman_marker_vs_signed_nearest_islet_um": distance_rho,
        "n_t_cells_marker_and_distance_measured": distance_n,
    }
    inflammation = None
    if inflammatory and tissue == "pancreas" and "cell_class" in frame:
        cell_class = frame.cell_class.astype("string")
        immune = frame.qc_analysis_eligible.fillna(False).to_numpy(bool) & cell_class.eq("Immune").fillna(False).to_numpy()
        other = frame.qc_analysis_eligible.fillna(False).to_numpy(bool) & cell_class.notna().to_numpy() & cell_class.ne("Immune").fillna(False).to_numpy()
        inflammation = {
            "donor_id": donor, "disease_status": disease, "tissue": tissue,
            "marker": marker_name, "source_marker_column": marker_col,
            "n_immune_cells": int(immune.sum()), "n_immune_marker_measured": int(np.isfinite(raw[immune]).sum()),
            "n_other_native_typed_cells": int(other.sum()), "n_other_marker_measured": int(np.isfinite(raw[other]).sum()),
            "marker_auc_immune_vs_other_native_typed": auc_immune_vs_other(raw[immune], raw[other]),
            "interpretation_limit": "within-donor signal distribution; calls and markers may share evidence; not independent specificity or a disease effect",
        }
    return row, inflammation


def summarize_sample(donor: str, disease: str, tissue: str, frame: pd.DataFrame,
                     markers: dict[str, str], distance_col: str | None) -> tuple[dict, list[dict], list[dict]]:
    eligible = frame.qc_analysis_eligible.fillna(False).to_numpy(bool)
    if tissue == "pancreas":
        immune_class = frame.immune_class.astype("string")
        t_mask = eligible & immune_class.eq("Imm-T-cell").fillna(False).to_numpy()
        subtype = frame.immune_subtype.astype("string")
        strict_cd8 = t_mask & subtype.eq("T-CD8").fillna(False).to_numpy()
        known_noncd8 = t_mask & subtype.isin(["T-CD4-helper", "T-regulatory"]).fillna(False).to_numpy()
        explicit_cd8 = t_mask & subtype.isin(["T-CD8", "T-CD8-cytotoxic"]).fillna(False).to_numpy()
        cell_class = frame.cell_class.astype("string")
        immune_mask = eligible & cell_class.eq("Immune").fillna(False).to_numpy()
        distance_col = distance_col or "signed_dist_islet_um"
        intensity_prefix = "typing_scaled_"
    elif tissue == "pLN":
        call = frame.t_nk_class.astype("string")
        t_mask = eligible & call.isin(["T_CD8", "T_CD4", "T_reg", "T_NK_unresolved"]).fillna(False).to_numpy()
        strict_cd8 = t_mask & call.eq("T_CD8").fillna(False).to_numpy()
        known_noncd8 = t_mask & call.isin(["T_CD4", "T_reg"]).fillna(False).to_numpy()
        explicit_cd8 = strict_cd8.copy()
        immune_mask = eligible & frame.immune_class.astype("string").ne("").fillna(False).to_numpy()
        intensity_prefix = ""
    else:
        raise ValueError(tissue)

    distance_values = finite_values(frame[distance_col]) if distance_col and distance_col in frame else None
    if distance_values is not None:
        d = distance_values[t_mask]
        called = explicit_cd8[t_mask].astype(float)
        identity_distance_rho = rho(d, called)
        identity_distance_n = int((np.isfinite(d) & np.isfinite(called)).sum())
    else:
        identity_distance_rho, identity_distance_n = np.nan, 0
    tissue_row = {
        "donor_id": donor, "disease_status": disease, "tissue": tissue,
        "n_qc_analysis_eligible_cells": int(eligible.sum()), "n_native_t_cells": int(t_mask.sum()),
        "n_native_strict_cd8_calls": int(strict_cd8.sum()),
        "n_native_explicit_cd8_calls": int(explicit_cd8.sum()),
        "n_native_known_noncd8_t_calls": int(known_noncd8.sum()),
        "fraction_explicit_cd8_calls_among_native_t_cells": float(explicit_cd8.sum() / t_mask.sum()) if t_mask.sum() else np.nan,
        "n_native_immune_calls": int(immune_mask.sum()),
        "native_cd8_call_distance_spearman_signed_nearest_islet_um": identity_distance_rho,
        "n_t_cells_distance_measured": identity_distance_n,
        "distance_measurement": "continuous signed distance to nearest islet in micrometres; no distance bands" if distance_values is not None else "not applicable to pLN",
    }
    rows, inflammation = [], []
    for marker, source_marker in markers.items():
        column = intensity_prefix + source_marker
        row, inflammatory_row = marker_summary(
            donor, disease, tissue, frame, t_mask, strict_cd8, known_noncd8,
            marker, column, distance_col,
            inflammatory=marker in INFLAMMATORY_MARKERS)
        rows.append(row)
        if inflammatory_row:
            inflammation.append(inflammatory_row)
    return tissue_row, rows, inflammation


def donor_aggregates(marker_rows: pd.DataFrame, inflammation_rows: pd.DataFrame) -> pd.DataFrame:
    summaries = []
    for tissue, marker, family, values_col, null in [
        (r.tissue, r.marker, "Peters_T_cell_marker_AUC", "marker_auc_strict_cd8_vs_known_noncd8_t", 0.5)
        for r in marker_rows[["tissue", "marker"]].drop_duplicates().itertuples(index=False)
        if r.marker in PETERS_MARKERS
    ]:
        values = marker_rows.loc[(marker_rows.tissue == tissue) & (marker_rows.marker == marker), values_col].to_numpy(float)
        values = values[np.isfinite(values)]
        lo, hi = bootstrap_mean_ci(values, SEED + len(summaries))
        summaries.append({"tissue": tissue, "endpoint_family": family, "marker": marker,
                          "donor_balanced_mean": float(np.mean(values)) if len(values) else np.nan,
                          "donor_balanced_median": float(np.median(values)) if len(values) else np.nan,
                          "bootstrap_95ci_low": lo, "bootstrap_95ci_high": hi,
                          "n_donors": int(len(values)), "donor_level_wilcoxon_vs_null_p": one_sample_p(values, null)})
    for marker in PETERS_MARKERS:
        part = marker_rows[(marker_rows.tissue == "pancreas") & marker_rows.marker.eq(marker)]
        values = part["distance_spearman_marker_vs_signed_nearest_islet_um"].to_numpy(float)
        values = values[np.isfinite(values)]
        lo, hi = bootstrap_mean_ci(values, SEED + len(summaries))
        summaries.append({"tissue": "pancreas", "endpoint_family": "Peters_continuous_distance_association",
                          "marker": marker, "donor_balanced_mean": float(np.mean(values)) if len(values) else np.nan,
                          "donor_balanced_median": float(np.median(values)) if len(values) else np.nan,
                          "bootstrap_95ci_low": lo, "bootstrap_95ci_high": hi,
                          "n_donors": int(len(values)), "donor_level_wilcoxon_vs_null_p": one_sample_p(values, 0.0)})
    for marker in INFLAMMATORY_MARKERS:
        part = inflammation_rows[inflammation_rows.marker.eq(marker)]
        values = part.marker_auc_immune_vs_other_native_typed.to_numpy(float)
        values = values[np.isfinite(values)]
        lo, hi = bootstrap_mean_ci(values, SEED + len(summaries))
        summaries.append({"tissue": "pancreas", "endpoint_family": "Medina_pancreatic_inflammatory_marker_distribution",
                          "marker": marker, "donor_balanced_mean": float(np.mean(values)) if len(values) else np.nan,
                          "donor_balanced_median": float(np.median(values)) if len(values) else np.nan,
                          "bootstrap_95ci_low": lo, "bootstrap_95ci_high": hi,
                          "n_donors": int(len(values)), "donor_level_wilcoxon_vs_null_p": one_sample_p(values, 0.5)})
    result = pd.DataFrame(summaries)
    result["bh_q_within_endpoint_family"] = np.nan
    for _, idx in result.groupby(["endpoint_family", "tissue"], dropna=False).groups.items():
        valid = result.loc[idx, "donor_level_wilcoxon_vs_null_p"].notna()
        use = np.asarray(idx)[valid.to_numpy()]
        if len(use):
            result.loc[use, "bh_q_within_endpoint_family"] = multipletests(
                result.loc[use, "donor_level_wilcoxon_vs_null_p"].to_numpy(float), method="fdr_bh")[1]
    return result


def paired_tissue_auc(marker_rows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for marker in PETERS_MARKERS:
        part = marker_rows[marker_rows.marker.eq(marker)
                           & marker_rows.tissue.isin(["pancreas", "pLN"])]
        wide = part.pivot_table(index="donor_id", columns="tissue",
                                values="marker_auc_strict_cd8_vs_known_noncd8_t", aggfunc="first")
        if not {"pancreas", "pLN"}.issubset(wide.columns):
            continue
        wide = wide.dropna(subset=["pancreas", "pLN"])
        delta = (wide["pLN"] - wide["pancreas"]).to_numpy(float)
        lo, hi = bootstrap_mean_ci(delta, SEED + len(rows))
        rows.append({"marker": marker, "endpoint_family": "Peters_T_cell_marker_AUC",
                     "outcome": "donor_paired_difference_pLN_minus_pancreas",
                     "n_paired_donors": len(delta), "mean_paired_difference": float(np.mean(delta)) if len(delta) else np.nan,
                     "median_paired_difference": float(np.median(delta)) if len(delta) else np.nan,
                     "bootstrap_95ci_low": lo, "bootstrap_95ci_high": hi,
                     "wilcoxon_signed_rank_p": paired_p(delta), "donor_ids": ";".join(wide.index.astype(str)),
                     "interpretation_limit": "within-donor rank discrimination in native call groups; no absolute protein intensity comparison; cell calls are not independent validation truth"})
    out = pd.DataFrame(rows)
    if len(out):
        out["bh_q_within_endpoint_family"] = multipletests(out.wilcoxon_signed_rank_p, method="fdr_bh")[1]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pancreas-root", type=Path, required=True,
                    help="candidate pancreas cell release, partitioned by donor")
    ap.add_argument("--pln-cells-root", type=Path,
                    default=REPO / "core/phenocycler/data/integration/cells_pln")
    ap.add_argument("--pln-calls-root", type=Path,
                    default=REPO / "data/development/typing/pln/pheno")
    ap.add_argument("--input-manifest", type=Path,
                    default=REPO / "results/v3/release_20261003/current_input_manifest.json")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--analysis-assignments-root", type=Path, required=True,
                    help="current exact-key assignments; analysis_eligible excludes quarantined cells")
    args = ap.parse_args()
    pancreas_manifest = json.loads((args.pancreas_root / "manifest.json").read_text())
    if pancreas_manifest.get("status") == "production":
        arm_status = "production"
    else:
        arm_status = "candidate_validation_build"
    input_manifest = json.loads(args.input_manifest.read_text())
    diseases = disease_map(input_manifest)
    tissue_summaries: list[dict] = []
    marker_rows: list[dict] = []
    inflammation_rows: list[dict] = []
    source_files: list[dict] = []

    pancreas_dirs = sorted(args.pancreas_root.glob("donor_id=*/data.parquet"))
    for path in pancreas_dirs:
        donor = path.parent.name.split("=", 1)[1]
        if donor == "6457" or donor.lstrip("0") == "6457":
            raise ValueError("excluded donor 6457 is present in pancreas candidate release")
        donor_disease = disease_for(donor, diseases)
        if donor_disease is None:
            raise ValueError(f"{donor}: disease status is missing from current inventory")
        columns = ["cell_id", "qc_analysis_eligible", "cell_class", "immune_class", "immune_subtype",
                   "signed_dist_islet_um", *["typing_scaled_" + c for c in PANCREAS_MARKERS.values()]]
        available_columns = set(pq.ParquetFile(path).schema.names)
        frame = pd.read_parquet(path, columns=[c for c in columns if c in available_columns])
        assignment_path = args.analysis_assignments_root / f"donor_id={int(donor)}" / "data.parquet"
        assignments = pd.read_parquet(assignment_path, columns=["cell_id", "analysis_eligible", "cell_class"])
        if frame.cell_id.duplicated().any() or assignments.cell_id.duplicated().any():
            raise ValueError(f"{donor}: duplicated literature/assignment cell key")
        assert len(frame) == len(assignments)
        aligned = assignments.set_index("cell_id").reindex(frame.cell_id)
        assert aligned.analysis_eligible.notna().all(), f"{donor}: stale/missing cell join"
        source_qc = frame.qc_analysis_eligible.fillna(False).to_numpy(bool)
        retained = aligned.analysis_eligible.fillna(False).to_numpy(bool)
        assert not (retained & ~source_qc).any()
        assert not (retained & aligned.cell_class.eq("Excluded").to_numpy(bool)).any()
        frame["source_qc_analysis_eligible"] = source_qc
        # All summaries and comparator groups use the same retained-cell mask.
        # Original QC remains provenance, never an analytical denominator.
        frame["qc_analysis_eligible"] = retained
        frame["cell_class"] = aligned.cell_class.to_numpy()
        tissue_row, rows, inflammation = summarize_sample(
            donor, donor_disease, "pancreas", frame, PANCREAS_MARKERS, "signed_dist_islet_um")
        tissue_row["n_source_qc_eligible_cells"] = int(source_qc.sum())
        tissue_row["n_quarantined_qc_cells"] = int((source_qc & ~retained).sum())
        tissue_row["analysis_assignment_path"] = str(assignment_path)
        tissue_summaries.append(tissue_row); marker_rows.extend(rows); inflammation_rows.extend(inflammation)
        source_files.append({"donor_id": donor, "tissue": "pancreas", "path": str(path.resolve()),
                             "bytes": path.stat().st_size, "sha256": file_sha256(path),
                             "analysis_assignment_path": str(assignment_path),
                             "analysis_assignment_sha256": file_sha256(assignment_path),
                             "quarantined_qc_cells": int((source_qc & ~retained).sum())})
        print(f"pancreas {donor}: {tissue_row['n_native_t_cells']} T cells; {len(rows)} markers", flush=True)

    pln_raw_files = sorted(args.pln_cells_root.glob("donor_id=*/roi=*/data_*.parquet"))
    if not pln_raw_files:
        raise FileNotFoundError(f"no pLN source cell partitions under {args.pln_cells_root}")
    for raw_path in pln_raw_files:
        donor = raw_path.parents[1].name.split("=", 1)[1]
        roi = raw_path.parent.name.split("=", 1)[1]
        call_path = args.pln_calls_root / f"donor_id={donor}" / f"roi={roi}" / "data.parquet"
        if not call_path.is_file():
            raise FileNotFoundError(f"{donor}/{roi}: native pLN calls unavailable: {call_path}")
        names = set(PLN_MARKERS.values())
        available_columns = set(pq.ParquetFile(raw_path).schema.names)
        raw = pd.read_parquet(raw_path, columns=["object_id", "qc_analysis_eligible",
                                                 *[m for m in names if m in available_columns]])
        calls = pd.read_parquet(call_path, columns=["object_id", "immune_class", "t_nk_class"])
        if raw.object_id.duplicated().any() or calls.object_id.duplicated().any():
            raise ValueError(f"{donor}/{roi}: duplicate pLN cell key")
        frame = raw.merge(calls, on="object_id", how="left", validate="one_to_one")
        tissue_row, rows, _ = summarize_sample(donor, disease_for(donor, diseases) or "unavailable",
                                               "pLN", frame, PLN_MARKERS, None)
        tissue_row["roi"] = roi
        for row in rows:
            row["roi"] = roi
        tissue_summaries.append(tissue_row); marker_rows.extend(rows)
        source_files.append({"donor_id": donor, "tissue": "pLN", "roi": roi,
                             "raw_path": str(raw_path.resolve()), "calls_path": str(call_path.resolve()),
                             "raw_sha256": file_sha256(raw_path), "calls_sha256": file_sha256(call_path)})
        print(f"pLN {donor}/{roi}: {tissue_row['n_native_t_cells']} T cells; {len(rows)} markers", flush=True)

    tissue_df = pd.DataFrame(tissue_summaries)
    assert tissue_df.loc[tissue_df.tissue.eq("pancreas"), "n_quarantined_qc_cells"].sum() == 11484
    assert tissue_df.loc[tissue_df.tissue.eq("pancreas"), "n_qc_analysis_eligible_cells"].sum() == 31385695
    marker_df = pd.DataFrame(marker_rows)
    inflammation_df = pd.DataFrame(inflammation_rows)
    coverage = (marker_df.groupby(["tissue", "marker", "source_marker_column"], as_index=False)
                .agg(n_donors=("donor_id", "nunique"),
                     n_donors_with_qc_marker=("marker_available_in_qc_cells", "sum"),
                     n_donors_with_t_cell_marker=("n_t_cells_marker_measured", lambda x: int(x.gt(0).sum())),
                     qc_cells_measured=("n_qc_eligible_cells_marker_measured", "sum"),
                     t_cells_measured=("n_t_cells_marker_measured", "sum")))
    paired = paired_tissue_auc(marker_df)
    summaries = donor_aggregates(marker_df, inflammation_df)
    args.out.mkdir(parents=True, exist_ok=True)
    tissue_df.to_csv(args.out / "protein_literature_donor_tissue_counts.csv", index=False)
    marker_df.to_csv(args.out / "protein_literature_donor_marker_endpoints.csv", index=False)
    coverage.to_csv(args.out / "protein_literature_marker_coverage.csv", index=False)
    inflammation_df.to_csv(args.out / "protein_literature_donor_inflammation_endpoints.csv", index=False)
    paired.to_csv(args.out / "protein_literature_paired_tissue_effects.csv", index=False)
    summaries.to_csv(args.out / "protein_literature_donor_balanced_summary.csv", index=False)
    payload = {
        "title": "Donor-level Peters and Medina-Serpas protein endpoints",
        "status": arm_status, "pancreas_calibration_arm": pancreas_manifest.get("calibration_arm"),
        "pancreas_release_manifest_sha256": file_sha256(args.pancreas_root / "manifest.json"),
        "production_ruling_sha256": pancreas_manifest.get("production_ruling_sha256"),
        "current_input_manifest_sha256": file_sha256(args.input_manifest),
        "analysis_script": str(Path(__file__).resolve()), "analysis_script_sha256": file_sha256(Path(__file__).resolve()),
        "pancreas_marker_source": (
            "typing_scaled_<marker> from the selected production pancreas arm"
            if arm_status == "production" else
            "typing_scaled_<marker> from the selected validation/candidate pancreas arm"),
        "pln_marker_source": "native donor/ROI QuPath measurements joined by exact object_id to dedicated pLN calls",
        "donor_unit": True,
        "exclusion_policy": "All 11,484 pattern-reviewed protein cells are excluded from every endpoint, comparator, and analytic denominator; original QC counts are metadata only",
        "analysis_assignments_root": str(args.analysis_assignments_root.resolve()),
        "quarantined_cells": 11484,
        "cell_level_inference": False,
        "bootstrap": {"replicates": BOOTSTRAPS, "seed": SEED, "resampling_unit": "donor"},
        "multiplicity": "Benjamini-Hochberg within each tissue x predefined endpoint family",
        "pancreas_cd8_identity": "native Imm-T-cell leaves T-CD8 and T-CD8-cytotoxic, with strict T-CD8 reported separately; T-exhausted and T-progenitor-memory are not assumed CD8-only because their rules permit CD4 or CD8",
        "pLN_cd8_identity": "native t_nk_class T_CD8; known non-CD8 comparator calls are T_CD4 and T_reg; unresolved calls are not relabelled",
        "features": {"Peters": list(PETERS_MARKERS), "Medina-Serpas_pancreatic_inflammation": list(INFLAMMATORY_MARKERS)},
        "paired_tissue_estimand": "donor-paired difference in within-tissue rank AUC of marker signal in native strict CD8 versus known non-CD8 T-cell calls (pLN minus pancreas); not an absolute intensity or cell-intrinsic tissue effect",
        "distance_estimand": "within-donor Spearman association of measured marker signal or explicit native CD8 calls with continuous signed distance to nearest islet in pancreas; no distance bands",
        "guardrails": [
            "The donor, not the cell, is the inferential unit; all per-cell values first become within-donor estimates.",
            "No positivity threshold or expected disease ordering is introduced.",
            "Absolute protein intensity is not compared across pancreas and pLN panels or calibrations.",
            "AUCs involving native calls are descriptive because typing rules and marker evidence can overlap; they are not independent accuracy validation.",
            "TCF1, TOX, Granzyme B, PD-1, LAG3 and other targets remain separate marker endpoints, not a fabricated complete state signature.",
            "pLN has no pancreas-islet distance outcome; continuous distance associations are pancreas-only.",
            "No disease adjustment is fit where disease and batch are inseparable.",
        ],
        "sources": source_files,
        "outputs": sorted(str(p.resolve()) for p in args.out.glob("*.csv")),
    }
    (args.out / "protein_literature_analysis_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {len(marker_df)} donor/marker rows, {len(paired)} paired marker effects -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
