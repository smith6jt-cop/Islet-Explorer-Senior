#!/usr/bin/env python3
"""Cohort driver and release writer for the pancreas Xenium 5K typing.

Runs ``panc_xenium_cluster_name.py`` once per pancreas section -- each in its own process, so
one section's memory is returned to the machine before the next starts -- and then assembles a
declared release from the per-section artifacts.

WHAT THE RELEASE IS
-------------------
``release.json`` (schema_version 1, its own release_id, roi ``panc``, the donor list, a sha256
per input), ``label_space.json``, and ``calls/donor_id=<d>/roi=panc/data.parquet`` keyed on
``cell_id``.  It is a DEVELOPMENT release: nothing here is production typing and no latch moves.

WHERE IT IS WRITTEN
-------------------
Under ``/data``.  ``/home`` is at 97%, and a cohort of 21 sections writes tens of GB.

CONTROLS AND COHERENCE
----------------------
* The Y-probe control is a control ONLY on a female donor, where every detection of a
  zero-template Y probe is false.  It is therefore reported per section AND aggregated over the
  FEMALE donors alone; a male donor's Y rate is a measure of that donor's Y expression and is
  reported as such, never as a background estimate.
* Depth, not gDNA, is the constraint on a 5K-panel call.  Per-donor median transcripts per cell
  is reported and a section is flagged on depth alone.
* INTERNAL COHERENCE: if the genes that MAY NAME a lineage are detected in a section and that
  lineage still receives no cells, that is a scoring failure to investigate.  It is a comparison
  of a section against itself.  No expected fraction exists for any lineage and none is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORKER = HERE / "panc_xenium_cluster_name.py"

sys.path.insert(0, str(HERE))
import importlib.util as _u
_spec = _u.spec_from_file_location("panc_worker", WORKER)
worker = _u.module_from_spec(_spec)
_spec.loader.exec_module(worker)

ROI = "panc"
#: A section whose median transcripts per cell is at or below this is flagged LOW DEPTH.  It is
#: a flag on the substrate, not a rejection: the calls are still written and the flag travels
#: with them.  Depth is what limits a 5K-panel call; it is not a gDNA statement.
LOW_DEPTH_MEDIAN = 100.0
#: The zero-template floor probes, read on FEMALE donors only.
Y_TARGET_FEMALE = 0.003


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def run_one(section: str, work: str, score_on: str) -> tuple[str, int, str]:
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, str(WORKER), "--section", section, "--work", work,
         "--score-on", score_on],
        capture_output=True, text=True,
    )
    tail = (proc.stdout + proc.stderr).strip().splitlines()
    return section, proc.returncode, " | ".join(tail[-4:]) + f" [{time.time() - t0:.0f}s]"


def stage_sections(work: Path, score_on: str, jobs: int, only: list[str] | None,
                   force: bool) -> list[tuple[str, str, Path]]:
    todo, done = [], []
    secs = worker.sections()
    if only:
        secs = [s for s in secs if s[1] in only or s[0] in only]
    for donor, section, h5ad in secs:
        if (work / f"{section}.report.json").is_file() and not force:
            done.append(section)
        else:
            todo.append(section)
    log(f"{len(secs)} sections; {len(done)} already in {work}; {len(todo)} to run on {jobs} jobs")
    if todo:
        work.mkdir(parents=True, exist_ok=True)
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(run_one, s, str(work), score_on) for s in todo]
            for f in futs:
                section, rc, tail = f.result()
                log(f"{'OK ' if rc == 0 else 'FAIL'} {section}: {tail}")
    return secs


def assemble(work: Path, out: Path, score_on: str, secs) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    reports, per_donor, frames = {}, [], []
    for donor, section, h5ad in secs:
        rp = work / f"{section}.report.json"
        if not rp.is_file():
            log(f"SKIP {section}: no report (the section did not complete)")
            continue
        reports[section] = json.loads(rp.read_text())

    log(f"hashing {len(reports)} h5ad inputs")
    with ProcessPoolExecutor(max_workers=6) as ex:
        hashes = dict(zip(
            [s for s in reports],
            ex.map(sha256_of, [Path(reports[s]["h5ad"]) for s in reports]),
        ))

    label_names = list(worker.LINEAGE_CORES)
    for donor, section, h5ad in secs:
        if section not in reports:
            continue
        r = reports[section]
        calls = pd.read_parquet(work / f"{section}.calls.parquet")
        # How DECIDED each cluster's call was, carried per cell.  The method admits the argmax
        # on one standout naming gene, so a cluster that won by 0.008 and a cluster that won by
        # 2.1 both produce a label; only one of them is resolved.  Nothing is suppressed on this
        # and no threshold is invented -- the margin and the per-cell agreement travel with the
        # call so an undecided one is visible where it is used.
        subs = pd.read_parquet(work / f"{section}.subclusters.parquet",
                               columns=["subcluster_id", "score_margin_over_runner_up",
                                        "n_naming_above_background_of_winning_core",
                                        "purity_cell_argmax_agrees",
                                        "cell_argmax_entropy"])
        subs = subs.rename(columns={
            "score_margin_over_runner_up": "cluster_score_margin",
            "n_naming_above_background_of_winning_core": "cluster_naming_standouts",
            "purity_cell_argmax_agrees": "cluster_purity_cell_argmax_agrees",
            "cell_argmax_entropy": "cluster_cell_argmax_entropy"})
        calls = calls.merge(subs, on="subcluster_id", how="left")
        # A missing measurement is UNAVAILABLE.  Older work artifacts wrote -1 where the nucleus
        # rule was unmeasurable; the release never carries a negative sentinel for it.
        measured = bool(r["nucleus_rule_measurable"])
        calls["nucleus_count"] = pd.array(
            calls["nucleus_count"].to_numpy() if measured else [pd.NA] * len(calls),
            dtype="Int32")
        thin = calls["lineage"].isin(label_names) & (calls["cluster_score_margin"] < 0.05)
        dst = out / "calls" / f"donor_id={donor}" / f"roi={ROI}"
        dst.mkdir(parents=True, exist_ok=True)
        # donor_id and roi are the hive partition keys.  Repeating them inside the file makes the
        # tree unreadable as one dataset -- the path key of a numeric donor infers as int32 while
        # the column is a string -- so they live in the path only.  section_key stays in the file
        # because it is not a partition key and a donor may one day carry more than one section.
        calls.drop(columns=["donor_id", "roi"]).to_parquet(dst / "data.parquet", index=False)
        counts = calls["lineage"].value_counts().to_dict()
        n = len(calls)
        n_called = int(sum(v for k, v in counts.items() if k in label_names))
        per_donor.append({
            "donor_id": donor, "section_key": section, "sex": r["sex"],
            "disease_status": r["disease_status"], "batch": r["batch"],
            "n_cells": n,
            "n_no_nucleus_excluded": r["n_cells_no_nucleus_excluded"],
            "n_below_cluster_qc": int(counts.get("Uncalled_below_cluster_qc", 0)),
            "n_called": n_called,
            "called_fraction_of_all_cells": round(n_called / n, 4),
            "called_fraction_of_nucleated_cells": round(
                n_called / max(n - r["n_cells_no_nucleus_excluded"], 1), 4),
            "median_transcripts_per_cell": r["median_transcripts_per_cell"],
            "median_genes_per_cell": r["median_genes_per_cell"],
            "low_depth_flag": bool(r["median_transcripts_per_cell"] <= LOW_DEPTH_MEDIAN),
            "n_clusters": r["cluster_meta"]["n_clusters"],
            "n_subclusters": r.get("n_subclusters"),
            "knn_holdout_agreement": round(r["cluster_meta"]["knn_holdout_agreement"], 4),
            "y_probe_detect_frac_any": r["y_probe_detect_frac_any"],
            "y_floor_detect_frac_any": r["y_floor_detect_frac_any"],
            "n_cd3_negative_nk_redirects": r["n_clusters_cd3_negative_nk_redirect"],
            "n_called_on_margin_below_0.05": int(thin.sum()),
            "n_called_on_one_naming_gene_in_cluster_below_depth_50": int(
                r["called_cells_on_one_naming_gene_in_low_depth_cluster"]),
            **{f"n__{lab}": int(counts.get(lab, 0)) for lab in label_names},
            **{f"n__{lab}": int(counts.get(lab, 0)) for lab in sorted(
                k for k in counts if k not in label_names)},
        })
        frames.append(pd.DataFrame({"donor_id": donor, "lineage": list(counts),
                                    "n": list(counts.values())}))
        del calls

    donors = sorted({d["donor_id"] for d in per_donor})
    females = [d for d in per_donor if str(d["sex"]).lower().startswith("f")]
    males = [d for d in per_donor if str(d["sex"]).lower().startswith("m")]

    def wmean(rows, key):
        if not rows:
            return None
        w = np.array([r["n_cells"] for r in rows], dtype=float)
        v = np.array([r[key] for r in rows], dtype=float)
        return float((v * w).sum() / w.sum())

    control = {
        "rule": "the Y probes are a background read only where the template is absent",
        "floor_group": worker.Y_FLOOR,
        "crosshyb_excluded_from_the_control": worker.Y_CROSSHYB,
        "n_female_donors": len(females), "n_male_donors": len(males),
        "female_donors": sorted({d["donor_id"] for d in females}),
        "female_y_floor_detect_frac_any_cellweighted": wmean(females, "y_floor_detect_frac_any"),
        "female_y_probe_detect_frac_any_cellweighted": wmean(females, "y_probe_detect_frac_any"),
        "female_y_floor_per_section": {d["section_key"]: d["y_floor_detect_frac_any"]
                                       for d in females},
        "target_female_floor": Y_TARGET_FEMALE,
        "male_y_floor_detect_frac_any_cellweighted": wmean(males, "y_floor_detect_frac_any"),
        "male_note": ("a male donor's Y rate measures that donor's Y expression, not background; "
                      "it is reported and never used as a control"),
    }

    coh_rows = []
    for section, r in reports.items():
        for gene, blk in r.get("coherence_genes", {}).items():
            if not blk.get("on_panel"):
                coh_rows.append({"section_key": section, "gene": gene, "on_panel": False})
                continue
            coh_rows.append({"section_key": section, "gene": gene, "on_panel": True,
                             "det_frac_section": blk["det_frac_section"],
                             "n_cells_detected": blk["n_cells_detected"],
                             "max_cluster_det_frac": blk["max_cluster_det_frac"],
                             "max_cluster_lineage": blk["max_cluster_lineage"]})
    coh = pd.DataFrame(coh_rows)

    # A lineage cannot be absent from a section while its own naming evidence is detected there.
    incoherent = []
    watch = {"Immune_T": ["CD3E", "CD3G", "CD4", "CD8A"],
             "Immune_NK_ILC": ["NKG7", "KLRD1"],
             "Immune_B_Plasma": ["MS4A1"], "Myeloid": ["CD68"]}
    for row in per_donor:
        for lineage, genes in watch.items():
            if row.get(f"n__{lineage}", 0) > 0:
                continue
            sub = coh[(coh.section_key == row["section_key"]) & (coh.gene.isin(genes))
                      & coh.on_panel]
            if not sub.empty and float(sub.det_frac_section.max()) > 0:
                incoherent.append({
                    "section_key": row["section_key"], "lineage_with_zero_cells": lineage,
                    "naming_genes_detected": {g: float(v) for g, v in
                                              zip(sub.gene, sub.det_frac_section)},
                    "verdict": "SCORING FAILURE TO INVESTIGATE, not a measurement",
                })

    panel_files = sorted({p for p, _ in worker.LINEAGE_CORES.values()}) + \
        sorted(worker.EMBEDDING_EXTRA.values())
    inputs = [{"kind": "section_h5ad", "donor_id": d, "section_key": s,
               "path": reports[s]["h5ad"], "sha256": hashes[s]}
              for d, s, _ in secs if s in reports]
    inputs += [{"kind": "identity_core", "path": str(p.relative_to(REPO)),
                "sha256": sha256_of(p)} for p in panel_files]
    inputs.append({"kind": "cohort_exclusions",
                   "path": "core/phenocycler/phenocycler/cohort.py",
                   "sha256": sha256_of(REPO / "core/phenocycler/phenocycler/cohort.py"),
                   "note": "DONOR_EXCLUSIONS is read from here; no literal donor list is used"})

    pdf = pd.DataFrame(per_donor)
    pdf.to_csv(out / "per_donor_summary.csv", index=False)
    coh.to_csv(out / "coherence_genes.csv", index=False)
    pd.concat(frames).pivot_table(index="donor_id", columns="lineage", values="n",
                                  fill_value=0).to_csv(out / "lineage_counts.csv")

    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", str(REPO), "status", "--porcelain"],
                                capture_output=True, text=True).stdout.strip())

    label_space = {
        "schema_version": 1,
        "space_id": "xenium-panc-cluster-name-v1",
        "producer_lane": "xenium-panc-cluster-name",
        "labels": [
            {"label": lab, "common_lineage": worker.LINEAGE_CORES[lab][1],
             "core": str(worker.LINEAGE_CORES[lab][0].relative_to(REPO)),
             "n_genes_on_panel": reports[next(iter(reports))]["core_sizes"][lab][
                 "n_genes_on_panel"],
             "n_naming_genes_on_panel": reports[next(iter(reports))]["core_sizes"][lab][
                 "n_naming_on_panel"],
             "reason": "named by canonical enrichment against the other clusters of the "
                       "same section, admitted on a standout naming gene of this core"}
            for lab in label_names
        ],
        "uncalled": [
            {"label": "Uncalled_no_naming_evidence",
             "reason": "the winning core has no naming gene standing out above the cluster's "
                       "own background; never reassigned to the runner-up"},
            {"label": "Uncalled_no_enrichment",
             "reason": "the winning core's score is not above the section mean at all"},
            {"label": "Uncalled_below_cluster_qc",
             "reason": f"nCount_Xenium < {worker.QC_MIN_COUNTS} or nFeature_Xenium < "
                       f"{worker.QC_MIN_FEATURES}; too little measured to partition on"},
            {"label": "Unavailable_no_nucleus",
             "reason": "nucleus_count == 0; excluded from every calculation, never zero"},
        ],
        "not_a_label": [
            {"name": "epithelial_general",
             "reason": "16_epithelial_general is an embedding extra. KRT/EPCAM are shared by "
                       "ductal, acinar and endocrine cells, so it is measured and reported and "
                       "is never a lineage."},
        ],
        "what_the_labels_do_not_support": {
            "Exocrine": "ONE label. The 5K panel carries no PRSS1/CPA1/CELA3A/CTRB1/CTRB2/"
                        "AMY2A/REG1A. No acinar and no ductal claim is made or supported.",
            "Endocrine": "INS, GCG, SST and PPY are absent from the 5,101-gene panel, measured "
                         "per section. Endocrine resolves only through the curated surrogate "
                         "panel 01_pancreas_endocrine. The call supports 'this cluster is "
                         "endocrine'; it supports no alpha / beta / delta statement.",
        },
    }
    (out / "label_space.json").write_text(json.dumps(label_space, indent=1))

    release = {
        "schema_version": 1,
        "release_id": f"xenium-panc-cluster-name-{datetime.now(timezone.utc):%Y%m%d}",
        "modality": "xenium",
        "roi": ROI,
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "producer": {"lane": "xenium-panc-cluster-name",
                     "worker": str(WORKER.relative_to(REPO)),
                     "driver": str(Path(__file__).resolve().relative_to(REPO)),
                     "commit": commit, "dirty": dirty},
        "label_space": "label_space.json",
        "calls_dir": "calls",
        "cell_key": "cell_id",
        "partition_keys": ["donor_id", "roi"],
        "partition_key_note": "donor_id and roi are in the directory path only, never repeated "
                              "as columns, so calls/ reads as one hive dataset. Read donor_id as "
                              "a STRING -- every donor in this cohort happens to be all digits, "
                              "so an inferring reader will turn it into an integer and 112 will "
                              "stop matching the donor it names.",
        "read_calls_with": "pyarrow.dataset.dataset(<release>/calls, format='parquet', "
                           "partitioning=pyarrow.dataset.partitioning(pyarrow.schema("
                           "[('donor_id', pyarrow.string()), ('roi', pyarrow.string())]), "
                           "flavor='hive'))",
        "columns": {"label": "lineage", "common_lineage": "common_lineage",
                    "reason": "reason", "cluster": "subcluster_id",
                    "parent_cluster": "cluster_id"},
        "stamp": "development",
        "latch": "no production typing is declared; this is a development release",
        "method": {
            "clustering": "within donor, per section; no joint embedding",
            "naming": "canonical enrichment against the other clusters of the SAME section; no "
                      "level, quantile or threshold crosses a slide",
            "score_on": score_on,
            "score_rule": "mean z over the genes of a core that MAY NAME its lineage "
                          "(core/panels evidence column: canonical or both)",
            "admission": "argmax core, admitted only on a standout naming gene of that core "
                         "above the cluster's own background; otherwise uncalled, never "
                         "reassigned",
            "negative_evidence": "a cluster whose argmax is Immune_T but in which neither CD3E "
                                 "nor CD3G stands out, while NK naming genes do, is "
                                 "Immune_NK_ILC",
            "detection": "lognorm >= log1p(1), lognorm = log1p(raw * 1e4 / cell depth); no raw "
                         "count floor and no depth-scaled equivalent",
            "nucleus_rule": "a cell with no nucleus is excluded from every calculation",
            "exclusions": "derived from phenocycler.cohort.DONOR_EXCLUSIONS, never a literal",
        },
        "donors": donors,
        "n_donors": len(donors),
        "n_sections": len(reports),
        "n_cells_total": int(pdf.n_cells.sum()),
        "n_cells_called": int(pdf.n_called.sum()),
        "n_cells_no_nucleus_excluded": int(pdf.n_no_nucleus_excluded.sum()),
        "inputs": inputs,
        "controls": control,
        "depth": {
            "rule": "depth, not gDNA, is the constraint on a 5K-panel call",
            "low_depth_median_threshold": LOW_DEPTH_MEDIAN,
            "sections_flagged_low_depth": [d["section_key"] for d in per_donor
                                           if d["low_depth_flag"]],
            "median_transcripts_per_cell": {d["section_key"]:
                                            d["median_transcripts_per_cell"]
                                            for d in per_donor},
        },
        "coherence": {
            "what_it_is": "a section compared against itself: a lineage cannot be absent while "
                          "its own naming evidence is detected. It is not a check against any "
                          "prior biological finding and no expected fraction is used.",
            "genes": worker.COHERENCE_GENES,
            "table": "coherence_genes.csv",
            "sections_failing": incoherent,
        },
        "decidedness": {
            "what_it_is": "the argmax is admitted on evidence, not on a margin; these counts say "
                          "how decided each admitted call was, so an undecided one is visible "
                          "where it is used. No call was suppressed on a margin and no margin "
                          "threshold operates anywhere in the method.",
            "columns_on_every_called_cell": ["cluster_score_margin",
                                             "cluster_naming_standouts",
                                             "cluster_purity_cell_argmax_agrees",
                                             "cluster_cell_argmax_entropy",
                                             "n_naming_standouts_behind_cluster_call",
                                             "cluster_median_depth"],
            "n_called_on_margin_below_0.05": int(pdf["n_called_on_margin_below_0.05"].sum()),
            "n_called_on_one_naming_gene_in_cluster_below_depth_50": int(
                pdf["n_called_on_one_naming_gene_in_cluster_below_depth_50"].sum()),
        },
        "per_donor": per_donor,
    }
    (out / "release.json").write_text(json.dumps(release, indent=1))
    log(f"release written to {out}")
    return release


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/data/xenium_typing/panc_2026-09-17/_work")
    ap.add_argument("--release", default="/data/xenium_typing/panc_2026-09-17/release")
    ap.add_argument("--score-on", choices=("naming", "all"), default=worker.SCORE_ON_DEFAULT)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--assemble-only", action="store_true")
    args = ap.parse_args()
    work, out = Path(args.work), Path(args.release)
    if args.assemble_only:
        secs = worker.sections()
        if args.sections:
            secs = [s for s in secs if s[1] in args.sections or s[0] in args.sections]
    else:
        secs = stage_sections(work, args.score_on, args.jobs, args.sections, args.force)
    assemble(work, out, args.score_on, secs)


if __name__ == "__main__":
    main()
