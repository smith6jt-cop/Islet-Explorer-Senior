#!/usr/bin/env python3
"""Full corrected RNA generation on frozen native QC/partitions, all25 sections.

Uses the generic producer's selectable shared_er_exact_v1 policy, with reviewed
serial-qualified uncertainty decisions. Full native-style cell calls refresh
every lineage-derived field; original values survive as native_* provenance.
Raw detection is recomputed where changed marker roles or labels require it.
The measured shared-ER features and all native membership/QC columns survive.
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

from trial_residual_pancreas_graph_refinement import module, sha

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/v3/islet_immune_reassessment_20261006/cell_typing_validation/candidate_naming_replay"
FROZEN = Path("/data/xenium_typing/panc_naming_repair_trial_20261008")
ADJ = Path("/data/xenium_typing/panc_naming_repair_adjudicated_candidate_20261008")
WORK = Path("/data/xenium_typing/panc_20261003_bundle/_work")
MF = ROOT / "results/v3/complete_cell_typing_paper_style_20261005/source_tables/xenium_immune_subpanels_20261005/manifest.json"
DERIVED = ["lineage", "common_lineage", "reason", "parent_cluster_lineage", "cell_core_argmax",
           "cell_agrees_with_cluster_call", "n_naming_detected_for_cluster_call",
           "n_naming_standouts_behind_cluster_call", "cluster_score_margin", "cluster_naming_standouts",
           "cluster_purity_cell_argmax_agrees", "cluster_cell_argmax_entropy"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("/data/xenium_typing/panc_corrected_RNA_release_20261008"))
    ap.add_argument("--review-decisions", type=Path, default=BASE / "adjudicated_candidate_20261008/all28_original_partition_uncertainty_decisions.csv")
    ap.add_argument("--review-evidence", type=Path,
                    help="bounded original-partition review JSON; verifies retained and changed dispositions against the pinned reused release")
    ap.add_argument("--report-out", type=Path, default=BASE / "corrected_RNA_release_20261008")
    ap.add_argument("--reuse-verified-naming-counts", type=Path,
                    help="immutable corrected release with identical policy/core hashes and known labels; bind hashes and reuse verified actual raw counters")
    args = ap.parse_args()
    if args.out.exists() or args.report_out.exists():
        raise ValueError("Preserve previous releases; choose fresh output directories")
    args.out.mkdir(parents=True)
    args.report_out.mkdir(parents=True)
    producer_path = ROOT / "analysis/xenium/scripts/panc_xenium_cluster_name.py"
    pc = module("corrected_full_release_producer", producer_path)
    raw = module("corrected_full_release_raw_reader", ROOT / "analysis/xenium/scripts/render_rna_candidate_roi_review.py")
    original_panels = pc.load_panels()
    panels = pc.naming_policy_panels(original_panels, "shared_er_exact_v1")
    names = list(pc.LINEAGE_CORES)
    core_sources = {str(v[0]): sha(v[0]) for v in pc.LINEAGE_CORES.values()}
    reused = None
    reused_manifest = None
    if args.reuse_verified_naming_counts:
        cache_manifest = args.reuse_verified_naming_counts / "release.json"
        reused = json.loads(cache_manifest.read_text())
        reused_manifest = reused
        if reused["producer_naming_policy"] != "shared_er_exact_v1" or reused["identity_core_sources_sha256"] != core_sources:
            raise ValueError("Verified counter cache has different naming policy or core features/roles")
        cache_snapshot = args.reuse_verified_naming_counts / "generic_worker_source_snapshot.py"
        if sha(cache_snapshot) != reused["generic_worker_sha256"]:
            raise ValueError("Verified counter cache producer snapshot hash differs")
        cache_worker = module("verified_counter_cache_producer", cache_snapshot)
        cached_roles = cache_worker.naming_policy_panels(original_panels, "shared_er_exact_v1")
        if any(cached_roles["cores"][c] != panels["cores"][c] for c in names):
            raise ValueError("Verified counter cache has different active naming/support roles")
        reused = {r["section_id"]: r for r in reused["sections"]}
    if any(original_panels["cores"][c][0] != panels["cores"][c][0] for c in names):
        raise ValueError("Measured support-core features changed")
    inventory = [r for r in json.loads(MF.read_text())["sections"] if r["tissue"] == "pancreas"]
    decisions = pd.read_parquet(ADJ / "all25_adjudicated_candidate_partition_decisions.parquet")
    frozen = pd.read_parquet(FROZEN / "frozen_partition_candidate_decisions.parquet")
    review = pd.read_csv(args.review_decisions, dtype={"section_id": str, "source_serial": str})
    serial_map = {r["section_id"]: str(r["source_serial"]).zfill(7) for r in inventory}
    expected_serials = review.section_id.map(serial_map)
    if "source_serial" in review and not review.source_serial.str.zfill(7).eq(expected_serials).all():
        raise ValueError("Supplied review decisions carry a different source acquisition serial")
    review["source_serial"] = expected_serials
    if review.source_serial.isna().any() or review.duplicated(["section_id", "source_serial", "cluster_id"]).any():
        raise ValueError("Reviewed uncertainty decisions lack unique acquisition identities")
    review_evidence, evidence_rows, additional_keys = None, None, None
    if args.review_evidence:
        if reused is None:
            raise ValueError("Bounded review evidence must bind to a verified preceding release")
        review_evidence = json.loads(args.review_evidence.read_text())
        evidence_rows = pd.DataFrame(review_evidence["decisions"])
        evidence_rows["section_id"] = evidence_rows.section_id.astype(str)
        evidence_rows["source_serial"] = evidence_rows.source_serial.astype(str).str.zfill(7)
        if len(evidence_rows) != review_evidence["n_reviewed_original_partitions"] or evidence_rows.duplicated(["section_id", "source_serial", "original_partition_id"]).any():
            raise ValueError("Bounded evidence has non-unique or incomplete reviewed partitions")
        if not evidence_rows.source_serial.eq(evidence_rows.section_id.map(serial_map)).all():
            raise ValueError("Bounded evidence serial differs from pinned inventory")
        added_path = Path(review_evidence["uncertainty_exact_keys_path"])
        if sha(added_path) != review_evidence["uncertainty_exact_keys_sha256"]:
            raise ValueError("Bounded review exact-key evidence bytes differ")
        additional_keys = pd.read_parquet(added_path)
        if len(additional_keys) != review_evidence["n_additional_exact_uncertainty_cells"] or not additional_keys.serial_qualified_cell_key.is_unique:
            raise ValueError("Bounded review exact-key count/uniqueness differs")
        previous_path = args.reuse_verified_naming_counts / "serial_qualified_review_decisions.csv"
        if sha(previous_path) != reused_manifest["serial_qualified_review_decisions_sha256"]:
            raise ValueError("Preceding reviewed decisions differ from pinned release")
        previous_review = pd.read_csv(previous_path, dtype={"section_id": str, "source_serial": str})
        review_ix = review.set_index(["section_id", "source_serial", "cluster_id"])
        previous_ix = previous_review.set_index(["section_id", "source_serial", "cluster_id"])
        if not review_ix.reindex(previous_ix.index)[previous_ix.columns].equals(previous_ix):
            raise ValueError("The preceding uncertainty decisions must remain unchanged")
        specific = evidence_rows.loc[evidence_rows.disposition.eq("specific_keyed_uncertainty")]
        added_index = pd.MultiIndex.from_frame(specific[["section_id", "source_serial", "original_partition_id"]].rename(columns={"original_partition_id": "cluster_id"}))
        if set(review_ix.index) - set(previous_ix.index) != set(added_index) or set(added_index) & set(previous_ix.index):
            raise ValueError("Additional decisions must exactly match the disjoint bounded review scope")
        if not review_ix.reindex(added_index).candidate_uncertainty_state.tolist() == specific.recommended_candidate_lineage.tolist():
            raise ValueError("Added uncertainty states differ from the explicit review")
        (args.out / "bounded_original_partition_review_snapshot.json").write_bytes(args.review_evidence.read_bytes())
    review_path = args.out / "serial_qualified_review_decisions.csv"
    review.to_csv(review_path, index=False)
    assets, summaries, checks, scope, lineage_counts, final_stats = [], [], [], [], [], []
    for report in inventory:
        sec, serial, donor = report["section_id"], serial_map[report["section_id"]], str(report["donor_id"])
        native_path = Path(report["native_calls_path"])
        matrix_path = Path(report["h5ad_path"])
        if sha(matrix_path) != report["h5ad_sha256"]:
            raise ValueError("Actual matrix bytes differ from pinned section hash")
        if sha(native_path) != report["native_calls_sha256"]:
            raise ValueError("Pinned original native calls changed")
        native = pd.read_parquet(native_path)
        calls = pd.concat([native.copy(), native[DERIVED].add_prefix("native_")], axis=1)
        if not native.cell_id.is_unique or not native.section_key.eq(sec).all():
            raise ValueError("Native section keys are not exact unique identities")
        with h5py.File(matrix_path, "r") as h:
            matrix_ids = raw.read_array(h["obs"][h["obs"].attrs["_index"]]).astype(str)
        if len(matrix_ids) != len(native) or not pd.Index(matrix_ids).is_unique or (pd.Index(matrix_ids).get_indexer(native.cell_id) < 0).any():
            raise ValueError("Hash-pinned matrix/native exact cell membership is not a bijection")
        vendor_path = Path(report["source_acquisition_path"]) / "cells.parquet"
        if "__" + serial + "__" not in str(vendor_path):
            raise ValueError("Exact acquisition source path does not bind the declared serial")
        vendor_hash = sha(vendor_path)
        vendor = pd.read_parquet(vendor_path, columns=["cell_id", "nucleus_count"])
        if not vendor.cell_id.is_unique:
            raise ValueError("Acquisition-native vendor keys are duplicated")
        vendor_ix = pd.Index(vendor.cell_id).get_indexer(native.cell_id)
        if (vendor_ix < 0).any() or not np.array_equal(vendor.nucleus_count.to_numpy()[vendor_ix], native.nucleus_count.to_numpy()):
            raise ValueError("Serial-qualified vendor/native exact membership or nucleus values differ")
        del vendor, vendor_ix, matrix_ids
        qc = (native.depth.ge(10) & native.n_genes_detected.ge(5) & native.nucleus_count.gt(0)).fillna(False).astype(bool)
        table = frozen[frozen.section_id.eq(sec)].set_index("cluster_id")
        expected = decisions[decisions.section_id.eq(sec)].set_index("cluster_id")
        gp = WORK / (sec + ".gene_stats.parquet")
        stats = pd.read_parquet(gp)
        genes = sorted(stats.gene.unique())
        order = sorted(table.index)
        gi = {g: i for i, g in enumerate(genes)}
        def mat(col):
            return stats.pivot(index="cluster_id", columns="gene", values=col).reindex(index=order, columns=genes).to_numpy()
        z, det, rest, means = [mat(c) for c in ["z_across_clusters", "det_frac", "det_frac_rest_of_section", "mean_lognorm"]]
        sizes = table.reindex(order).n_cells.to_numpy(dtype=int)
        naming = {c: np.array([gi[g] for g in panels["cores"][c][1] if g in gi], dtype=int) for c in names}
        tc, nk = [np.array([gi[g] for g in gs if g in gi], dtype=int) for gs in [pc.T_DEFINING, pc.NK_NAMING_HINT]]
        result = pc.decide(z, det, rest, pc.Z_STANDOUT, names, naming, naming, tc, nk,
                           exact_detection_fdr=.05, cluster_sizes=sizes)
        if result["calls"] != table.reindex(order).candidate_lineage.tolist():
            raise ValueError("Selectable generic repaired policy does not exactly reproduce frozen all25 trial")
        labels = dict(zip(order, result["calls"]))
        reasons = dict(zip(order, result["reasons"]))
        here = review[review.section_id.eq(sec)]
        if not here.source_serial.eq(serial).all():
            raise ValueError("Reviewed source acquisition serial differs")
        for row in here.itertuples():
            if row.cluster_id not in table.index or int(table.loc[row.cluster_id, "n_cells"]) != int(row.n_cells):
                raise ValueError("Reviewed partition absent or native membership size differs")
            labels[row.cluster_id], reasons[row.cluster_id] = row.candidate_uncertainty_state, row.uncertainty_reason
        for cid in order:
            previous = expected.loc[cid, "candidate_lineage"]
            if labels[cid] != previous:
                if previous not in names or not labels[cid].startswith("Unresolved_") or cid not in set(here.cluster_id):
                    raise ValueError("Review changes must preserve original uncertainty and only add explicit known-to-uncertainty dispositions")
        calls.loc[qc, "lineage"] = calls.loc[qc, "subcluster_id"].map(labels)
        calls.loc[qc, "reason"] = calls.loc[qc, "subcluster_id"].map(reasons)
        calls["common_lineage"] = calls.lineage.map({c: pc.LINEAGE_CORES[c][1] for c in names}).fillna("")
        # Reconstruct initial-parent statistics from exact frozen child counts and
        # means. This changes naming only; parent and final memberships are fixed.
        parents = table.reindex(order).parent_cluster_id.to_numpy()
        po = sorted(set(parents))
        pmeans, pdets, psizes = [], [], []
        for parent in po:
            which = parents == parent
            psizes.append(sizes[which].sum())
            pmeans.append(np.average(means[which], axis=0, weights=sizes[which]))
            pdets.append(np.average(det[which], axis=0, weights=sizes[which]))
        pmeans, pdets, psizes = np.array(pmeans), np.array(pdets), np.array(psizes)
        ptotaldet = np.rint(pdets * psizes[:, None]).sum(0)
        prest = (ptotaldet[None, :] - pdets * psizes[:, None]) / (psizes.sum() - psizes)[:, None]
        pz = np.divide(pmeans - pmeans.mean(0), pmeans.std(0), out=np.zeros_like(pmeans), where=pmeans.std(0) > 0)
        native_naming = {c: np.array([gi[g] for g in original_panels["cores"][c][1] if g in gi], dtype=int) for c in names}
        native_parent_replay = pc.decide(pz, pdets, prest, pc.Z_STANDOUT, names, native_naming, native_naming, tc, nk)
        original_parents_path = WORK / (sec + ".clusters.parquet")
        original_parents = pd.read_parquet(original_parents_path).set_index("cluster_id").reindex(po)
        if original_parents.n_cells.isna().any() or not np.array_equal(psizes, original_parents.n_cells.to_numpy()):
            raise ValueError("Reconstructed initial-parent counts differ from original native summaries")
        if native_parent_replay["calls"] != original_parents.lineage_before_refinement.tolist():
            raise ValueError("Reconstructed native-policy initial-parent labels differ; do not claim generic reproducibility")
        pr = pc.decide(pz, pdets, prest, pc.Z_STANDOUT, names, naming, naming, tc, nk,
                       exact_detection_fdr=.05, cluster_sizes=psizes)
        calls.loc[qc, "parent_cluster_lineage"] = calls.loc[qc, "cluster_id"].map(dict(zip(po, pr["calls"])))
        known = qc & calls.lineage.isin(names)
        # Full-core support argmax remains mathematically unchanged: no measured
        # full-core genes or their weights change. It is descriptive support, not
        # an identity label or validation truth. Agreement is refreshed below.
        calls["cell_core_argmax_policy"] = "descriptive_full_core_support_including_shared_ER_not_cell_identity"
        calls["cell_agrees_with_cluster_call"] = (known & calls.cell_core_argmax.eq(calls.lineage) & calls.cell_core_argmax.ne("")).fillna(False).astype(bool)
        calls.loc[qc & ~known, "n_naming_detected_for_cluster_call"] = np.nan
        recalc = known & (calls.lineage.ne(native.lineage) | calls.lineage.eq("Immune_B_Plasma"))
        cache_record = None
        if reused is not None:
            cache_record = reused[sec]
            cache_path = Path(cache_record["corrected_calls_path"])
            if cache_record["source_serial"] != serial or sha(cache_path) != cache_record["corrected_calls_sha256"]:
                raise ValueError("Verified counter cache acquisition or call-byte hash differs")
            if cache_record["matrix_upstream_sha256"] != report["h5ad_sha256"] or cache_record["native_calls_sha256"] != report["native_calls_sha256"] or cache_record["serial_vendor_cells_sha256"] != vendor_hash:
                raise ValueError("Verified counter cache matrix/native/vendor bytes differ from current verified inputs")
            cache = pd.read_parquet(cache_path, columns=["cell_id", "serial_qualified_cell_key", "subcluster_id", "lineage", "n_naming_detected_for_cluster_call", "naming_policy"])
            if not cache.cell_id.equals(calls.cell_id) or not cache.naming_policy.eq("shared_er_exact_v1").all():
                raise ValueError("Verified counter cache exact key order or policy differs")
            if not cache.loc[known, "lineage"].equals(calls.loc[known, "lineage"]):
                raise ValueError("Known class changed relative to cached raw counters; recompute from matrix")
            if cache.loc[known, "n_naming_detected_for_cluster_call"].isna().any():
                raise ValueError("Verified raw naming counters are unavailable for a known class")
            if evidence_rows is not None:
                for erow in evidence_rows.loc[evidence_rows.section_id.eq(sec)].itertuples():
                    if erow.source_calls_path != str(cache_path) or erow.source_calls_sha256 != cache_record["corrected_calls_sha256"]:
                        raise ValueError("Bounded evidence is not bound to the current verified cache bytes")
                    group = cache.loc[cache.subcluster_id.eq(erow.original_partition_id)]
                    if len(group) != int(erow.n_native_QC_cells) or not group.lineage.eq(erow.current_candidate_lineage).all():
                        raise ValueError("Bounded original partition membership/label differs")
                    if not calls.loc[group.index, "lineage"].eq(erow.recommended_candidate_lineage).all():
                        raise ValueError("Final calls do not implement the explicit retained/uncertainty disposition")
                    if erow.disposition == "specific_keyed_uncertainty":
                        ek = additional_keys.loc[additional_keys.subcluster_id.eq(erow.original_partition_id)]
                        if len(ek) != len(group) or set(ek.serial_qualified_cell_key) != set(group.serial_qualified_cell_key):
                            raise ValueError("Bounded reviewed exact cell keys differ from the consumed original partition")
            calls.loc[known, "n_naming_detected_for_cluster_call"] = cache.loc[known, "n_naming_detected_for_cluster_call"]
            del cache
        elif recalc.any():
            selected_ids = calls.loc[recalc, "cell_id"].tolist()
            wanted = sorted({g for c in calls.loc[recalc, "lineage"].unique() for g in panels["cores"][c][1]})
            measured, counts = raw.raw_selected(Path(report["h5ad_path"]), selected_ids, wanted)
            values = pd.DataFrame(counts, index=calls.index[recalc], columns=measured)
            dd = values.to_numpy() * 10_000 / calls.loc[recalc, "depth"].to_numpy()[:, None] >= 1
            for c in calls.loc[recalc, "lineage"].unique():
                rows = calls.loc[recalc, "lineage"].eq(c).to_numpy()
                cols = [j for j, g in enumerate(measured) if g in panels["cores"][c][1]]
                calls.loc[values.index[rows], "n_naming_detected_for_cluster_call"] = dd[rows][:, cols].sum(1).astype(np.float32)
        standout_count, model_count, margin = {}, {}, dict(zip(order, result["margin"]))
        for k, cid in enumerate(order):
            label = labels[cid]
            standout_count[cid] = int(result["naming_ok"][k, names.index(label)]) if label in names else np.nan
            model_count[cid] = int(result["naming_ok"][k, result["win"][k]])
        calls.loc[qc, "n_naming_standouts_behind_cluster_call"] = calls.loc[qc, "subcluster_id"].map(standout_count)
        calls.loc[qc, "cluster_naming_standouts"] = calls.loc[qc, "subcluster_id"].map(standout_count)
        calls.loc[known, "cluster_score_margin"] = calls.loc[known, "subcluster_id"].map(margin)
        calls.loc[qc & ~known, "cluster_score_margin"] = np.nan
        calls["naming_model_score_margin"] = calls.subcluster_id.map(margin)
        calls["naming_model_winning_core_standouts"] = calls.subcluster_id.map(model_count)
        agrees = calls.loc[qc].groupby("subcluster_id").cell_agrees_with_cluster_call.mean()
        calls.loc[known, "cluster_purity_cell_argmax_agrees"] = calls.loc[known, "subcluster_id"].map(agrees)
        calls.loc[qc & ~known, "cluster_purity_cell_argmax_agrees"] = np.nan
        freqs = calls.loc[qc].groupby(["subcluster_id", "cell_core_argmax"], observed=True).size()
        freqs = freqs / freqs.groupby(level=0).sum()
        entropy = -(freqs * np.log(freqs + 1e-12)).groupby(level=0).sum()
        calls.loc[qc, "cluster_cell_argmax_entropy"] = calls.loc[qc, "subcluster_id"].map(entropy)
        calls["source_serial"] = serial
        calls["donor_id"], calls["roi"] = donor, "panc"
        calls["cell_key"] = sec + "|" + calls.cell_id.astype(str)
        calls["serial_qualified_cell_key"] = serial + "|" + calls.cell_key
        calls["naming_policy"] = "shared_er_exact_v1"
        calls["identity_uncertainty"] = calls.lineage.str.startswith("Unresolved_")
        calls["lineage_analysis_eligible"] = known & ~calls.identity_uncertainty
        target = args.out / "calls" / ("donor_id=" + donor) / ("section=" + sec) / "data.parquet"
        target.parent.mkdir(parents=True)
        calls.to_parquet(target, index=False, compression="zstd")
        readback = pd.read_parquet(target)
        for col in native.columns:
            actual = readback["native_" + col] if col in DERIVED else readback[col]
            if not actual.equals(native[col]):
                raise ValueError("Original native value lost: " + sec + ":" + col)
        for col in DERIVED:
            if not readback.loc[~qc, col].equals(native.loc[~qc, col]):
                raise ValueError("QC-excluded derived field altered: " + sec + ":" + col)
        if readback.lineage.isna().any() or not readback.cell_key.is_unique:
            raise ValueError("Corrected full release exact-key or lineage coverage failed")
        unresolved_known = qc & ~readback.lineage.isin(names)
        if readback.loc[unresolved_known, ["n_naming_detected_for_cluster_call", "n_naming_standouts_behind_cluster_call", "cluster_naming_standouts", "cluster_score_margin", "cluster_purity_cell_argmax_agrees"]].notna().any().any():
            raise ValueError("Unresolved identities retained asserted lineage-derived numeric fields")
        if readback.loc[unresolved_known, "lineage_analysis_eligible"].any():
            raise ValueError("Unresolved identities admitted to lineage summaries")
        immune = known & calls.common_lineage.isin(["B_Plasma", "T_NK", "Myeloid", "Granulocyte"])
        unresolved_n = int(unresolved_known.sum())
        summaries.append({"donor_id": donor, "section_id": sec, "source_serial": serial,
                          "n_all_source_cells": len(calls), "n_native_QC_cells": int(qc.sum()),
                          "n_named_lineage_cells": int(known.sum()), "n_QC_identity_unresolved_cells": unresolved_n,
                          "n_keyed_program_uncertainty_cells": int((qc & calls.identity_uncertainty).sum()),
                          "known_named_immune_cells_lower_bound": int(immune.sum()),
                          "complete_immune_cell_count": None if unresolved_n else int(immune.sum()),
                          "complete_immune_count_status": "unavailable_identity_uncertainty" if unresolved_n else "measured",
                          "unresolved_identity_is_measured_zero_immune": False})
        checks.append({"section_id": sec, "source_serial": serial, "n_cells": len(calls),
                       "all_native_values_and_original_partitions_preserved": True, "original_native_labels_preserved": True,
                       "all_QC_excluded_derived_fields_unchanged": True, "corrected_policy_exactly_reproduces_frozen_trial": True,
                       "reconstructed_native_parent_labels_and_counts_exact": True,
                       "original_adjudicated_decisions_preserved_with_only_serial_qualified_review_additions": True,
                       "all_readback_contracts_pass": True,
                       "n_actual_raw_naming_count_recomputed_cells": int(recalc.sum())})
        assets.append({"donor_id": donor, "section_id": sec, "source_serial": serial,
                       "native_calls_path": str(native_path), "native_calls_sha256": report["native_calls_sha256"],
                       "corrected_calls_path": str(target), "corrected_calls_sha256": sha(target),
                       "frozen_gene_stats_path": str(gp), "frozen_gene_stats_sha256": sha(gp),
                       "matrix_path": report["h5ad_path"], "matrix_upstream_sha256": report["h5ad_sha256"], "matrix_hash_verified_this_run": True,
                       "serial_vendor_cells_path": str(vendor_path), "serial_vendor_cells_sha256": vendor_hash,
                       "matrix_native_bijection_and_vendor_nucleus_binding_verified": True,
                       "original_parent_summary_path": str(original_parents_path), "original_parent_summary_sha256": sha(original_parents_path)})
        assets[-1]["raw_naming_counter_source"] = ({"verified_corrected_calls_path": cache_record["corrected_calls_path"],
                                                    "verified_corrected_calls_sha256": cache_record["corrected_calls_sha256"],
                                                    "counters_reused_after_exact_keys_and_unchanged_known_labels_check": True}
                                                   if cache_record else {"actual_hash_verified_matrix_read_this_run": True})
        checks[-1]["n_actual_raw_naming_count_recomputed_cells"] = 0 if reused is not None else int(recalc.sum())
        checks[-1]["verified_raw_naming_counters_reused"] = reused is not None
        cs = calls.groupby(["native_lineage", "lineage"], observed=True).size().reset_index(name="n_cells")
        cs["section_id"], cs["source_serial"] = sec, serial
        scope.append(cs)
        lc = calls.groupby(["lineage", "common_lineage"], observed=True).size().reset_index(name="n_cells")
        lc["section_id"], lc["source_serial"] = sec, serial
        lineage_counts.append(lc)
        for k, cid in enumerate(order):
            rec = {"section_id": sec, "source_serial": serial, "cluster_id": cid, "n_cells": int(sizes[k]),
                   "native_lineage": table.loc[cid, "native_lineage"], "lineage": labels[cid], "reason": reasons[cid],
                   "common_lineage": pc.LINEAGE_CORES.get(labels[cid], ("", ""))[1],
                   "naming_model_winning_core": names[int(result["win"][k])], "naming_model_score_margin": float(result["margin"][k]),
                   "n_admitted_naming_genes_of_call": standout_count[cid]}
            for j, c in enumerate(names):
                rec["score__" + c] = float(result["scores"][k, j])
                rec["naming_standouts__" + c] = int(result["naming_so"][k, j])
                rec["naming_admission__" + c] = int(result["naming_ok"][k, j])
            final_stats.append(rec)
        print(sec, len(calls), "source cells", int(recalc.sum()), "verified naming counter scope", "reused" if reused is not None else "recomputed", flush=True)
        del native, calls, readback
        gc.collect()
    pd.DataFrame(summaries).to_csv(args.out / "cohort_identity_and_immune_availability.csv", index=False)
    pd.DataFrame(final_stats).to_parquet(args.out / "corrected_all6492_partition_calls_and_scores.parquet", index=False)
    pd.concat(lineage_counts).to_csv(args.out / "corrected_lineage_counts.csv", index=False)
    pd.concat(scope).to_csv(args.report_out / "all25_original_native_to_corrected_call_scope.csv", index=False)
    pd.DataFrame(checks).to_csv(args.report_out / "all25_full_native_style_readback_verification.csv", index=False)
    manifest = {"release_status": "corrected_full25_frozen_partition_generation_for_review", "producer_naming_policy": "shared_er_exact_v1",
                "n_sections": len(assets), "n_source_cells": sum(r["n_all_source_cells"] for r in summaries),
                "n_native_QC_cells": sum(r["n_native_QC_cells"] for r in summaries), "sections": assets,
                "CLI_script_path": str(Path(__file__)), "CLI_script_sha256": sha(Path(__file__)),
                "generic_worker_path": str(producer_path), "generic_worker_sha256": sha(producer_path),
                "identity_core_sources_sha256": core_sources, "native_inventory_path": str(MF), "native_inventory_sha256": sha(MF),
                "review_input_path": str(args.review_decisions), "review_input_sha256": sha(args.review_decisions),
                "bounded_review_evidence_path": str(args.review_evidence) if args.review_evidence else None,
                "bounded_review_evidence_sha256": sha(args.review_evidence) if args.review_evidence else None,
                "bounded_review_all_dispositions_verified": evidence_rows is not None,
                "n_bounded_review_original_partitions": len(evidence_rows) if evidence_rows is not None else 0,
                "n_keyed_program_uncertainty_cells": sum(r["n_keyed_program_uncertainty_cells"] for r in summaries),
                "serial_qualified_review_decisions_path": str(review_path), "serial_qualified_review_decisions_sha256": sha(review_path),
                "verified_naming_counter_cache_release_path": str(args.reuse_verified_naming_counts / "release.json") if args.reuse_verified_naming_counts else None,
                "verified_naming_counter_cache_release_sha256": sha(args.reuse_verified_naming_counts / "release.json") if args.reuse_verified_naming_counts else None,
                "original_native_derived_field_provenance": {c: "native_" + c for c in DERIVED},
                "all_original_native_values_QC_and_memberships_preserved": True, "native_generation_modified": False,
                "cell_core_argmax": "Reused exactly because measured full-core features/weights are unchanged. Descriptive support diagnostic including ER3, never an activated per-cell identity or accuracy target. Agreement and purity are recomputed against corrected lineage.",
                "statistical_inference_limit": "Expression selected the clusters. Section-wide hypergeometric/BH q is an exploratory detection contrast after selection, not a guaranteed lineage false-discovery rate or independent validation.",
                "heldout_evidence_limit": "NR5A2/HNF4A were held out of native naming cores, not excluded from clustering/features. Independent reference context is never per-cell argmax relabeling.",
                "unresolved_immune_count_behavior": "Complete immune counts unavailable when native-QC identities remain unresolved; known named immune lower bounds and all unresolved denominators/counts are explicit.",
                "rejected_graph_children_activated": False, "whole_lineage_accuracy_completed": False,
                "producer_cli_reproducibility": "Run this CLI for the frozen cohort. The generic worker supports --naming-policy shared_er_exact_v1 --uncertainty-decisions serial_qualified_review_decisions.csv --source-serial SERIAL for future section processing; native default remains provenance-compatible."}
    (args.out / "release.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (args.report_out / "corrected_full25_release_receipt.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (args.out / "source_snapshot.py").write_bytes(Path(__file__).read_bytes())
    (args.out / "generic_worker_source_snapshot.py").write_bytes(producer_path.read_bytes())


if __name__ == "__main__":
    main()
