#!/usr/bin/env python
"""pLN arm (plan step B5): B1 metrics, B2 implied negative controls, B3 detection.

Extends the PANCREAS same-slide vs different-slide background comparison to the pancreatic
lymph node. Four pLN Xenium bundles exist for the two Pheno-Xeno donors:

    donor 6533 : 0080870 (same-slide, PhenoCycler first) | 0041393 (different slide, 2024)
    donor 112  : 0081366 (same-slide, PhenoCycler first) | 0041396 (different slide, 2024)

The two same-slide pLN bundles are NOT recorded in xenium_paths.csv (only their `panc`
siblings are), so they are named here explicitly and the provenance is written into every
output table rather than left to the reader.

CONTROL-SET RULE (the thing this script exists to get right)
------------------------------------------------------------
The control set is fixed BY BIOLOGY A PRIORI. No gene enters or leaves it because of how much
it is expressed in these data -- that is circular, and it makes `clean_neg` a blend of
background and signal rather than an estimate of background.

panel_audit.csv was written for PANCREAS. A lymph node is a different tissue, so every
category is re-derived here for lymph node, gene by gene, in LN_VERDICT below. A gene is
marked ln_control_valid = "no" when LYMPH-NODE BIOLOGY says a resident cell type transcribes
it -- macrophages, granulocytes, FRC/FDC stroma, HEV and vascular endothelium, vascular smooth
muscle, peri-nodal adipose, lymphocytes themselves. That verdict is recorded with its reason
and is independent of any number this script computes.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import xenio  # noqa: E402

OUT = pathlib.Path("/blue/maigan/smith6jt/Xenium_Analysis/data/processed/pln")
AUDIT = pathlib.Path("/blue/maigan/smith6jt/Xenium_Analysis/data/processed/panel_audit.csv")

PX = pathlib.Path("/orange/maigan/20260910__211837__6533_HDL112_pLNpanc_PhenoXeno")
NPOD = pathlib.Path("/orange/maigan/Xenium_T1D_PancPLN/MM036_Xenium_nPODBatch2")

#: donor -> section design -> bundle. `recorded` says whether xenium_paths.csv carries the row.
BUNDLES = [
    dict(donor_id="6533", donor_sex="F", serial="0080870", design="same_slide",
         phenocycler="0080870_6533_pLNpanc_Scan3.er.qptiff", recorded_in_xenium_paths=False,
         path=PX / "output-XETG00298__0080870__pLN__20260910__211918"),
    dict(donor_id="6533", donor_sex="F", serial="0041393", design="different_slide",
         phenocycler="", recorded_in_xenium_paths=True,
         path=NPOD / "output-XETG00298__0041393__pLN__20241019__205207"),
    dict(donor_id="112", donor_sex="F", serial="0081366", design="same_slide",
         phenocycler="0081366_HDL112_pLNpanc_Scan3.er.qptiff", recorded_in_xenium_paths=False,
         path=PX / "output-XETG00298__0081366__pLN__20260910__211918"),
    dict(donor_id="112", donor_sex="F", serial="0041396", design="different_slide",
         phenocycler="", recorded_in_xenium_paths=True,
         path=NPOD / "output-XETG00298__0041396__pLN__20241019__205207"),
]

# --- B2(a) zero-template Y probes -------------------------------------------------------
# Both donors are female (data/donor_metadata_panc.xlsx: 112 F, 6533 F), so every Y probe has
# NO template in any tissue. That argument is tissue-independent, so it carries from pancreas
# to lymph node unchanged.
Y_CLEAN = ("TBL1Y", "USP9Y", "EIF1AY")        # no highly expressed X paralog to cross-hyb from
Y_CROSSHYB = ("DDX3Y", "KDM5D")               # DDX3X / KDM5C paralogs -- reported SEPARATELY
Y_OTHER = ("TSPY1",)                          # Y-linked, testis-restricted; not in the established triple

#: gene -> (ln_control_valid, reason). Only genes whose LYMPH-NODE verdict differs from the
#: pancreas audit, or that carry a caveat worth recording, are listed; everything else in an
#: audited category inherits ln_control_valid = "yes".
LN_VERDICT = {
    "NLRP7":   ("no",  "NLRP7 is an inflammasome sensor of human MACROPHAGES/monocytes; a lymph node is full of them"),
    "ARG1":    ("no",  "arginase-1 is a granulocyte/myeloid-suppressor-cell gene as well as a hepatocyte gene"),
    "PCK1":    ("no",  "adipocyte glyceroneogenesis; hilar and peri-nodal adipose is part of a pLN section"),
    "GJA5":    ("no",  "connexin-40 marks ARTERIAL ENDOTHELIUM, not only cardiac conduction tissue"),
    "NPPC":    ("no",  "C-type natriuretic peptide is an ENDOTHELIAL/fibroblast gene"),
    "PLN":     ("no",  "phospholamban is expressed by SMOOTH MUSCLE (capsule, trabeculae, vessels), not only myocardium"),
    "RYR1":    ("no",  "RyR1 is reported in human B lymphocytes and dendritic cells"),
    "COL11A1": ("no",  "COL11A1 is an activated-fibroblast/CAF program; FRC and capsular fibroblasts are LN stroma"),
    "COMP":    ("no",  "COMP is expressed by fibrotic fibroblasts and vascular smooth muscle, not only cartilage"),
    "COL10A1": ("no",  "COL10A1 is reported in activated/cancer-associated fibroblasts as well as hypertrophic chondrocytes"),
    "WT1":     ("no",  "WT1 is expressed by MESOTHELIUM (retroperitoneal surfaces) and by haematopoietic progenitors"),
    "TAS2R38": ("no",  "human LEUKOCYTES express class-I/II bitter taste receptors incl. TAS2R38; also a 1-gene set"),
    # kept, but the caveat is on the record
    "KRT14":   ("yes", "basal keratinocyte / myoepithelium; a lymph node carries no epithelium"),
    "KRT1":    ("yes", "cornified epidermis; no template in lymph node (it IS expressed by pancreatic duct, so this gene is LN-only)"),
    "UMOD":    ("yes", "thick ascending limb only; no template in lymph node (pancreas expresses it, so LN-only)"),
    "AQP3":    ("yes", "ductal/centroacinar in pancreas and a live Gamma rule marker there; lymph node has no such epithelium"),
    "C8B":     ("yes", "terminal complement C8/C9 is predominantly hepatic; macrophage complement output is C1q/C3/CFB"),
    "C9":      ("yes", "terminal complement C8/C9 is predominantly hepatic; macrophage complement output is C1q/C3/CFB"),
    "OTX2":    ("yes", "retina/pineal/choroid plexus; absent from lymphoid tissue"),
    "RFX2":    ("yes", "ciliogenesis TF, broadest gene in the spermatogenesis set; lymph node has no ciliated cells"),
    "MAGEA1":  ("yes", "cancer-testis antigen; these are non-neoplastic nPOD donors"),
    "MAGEA3":  ("yes", "cancer-testis antigen; these are non-neoplastic nPOD donors"),
    "MAGEA4":  ("yes", "cancer-testis antigen; these are non-neoplastic nPOD donors"),
    "MAGEA10": ("yes", "cancer-testis antigen; these are non-neoplastic nPOD donors"),
    "MAGEA12": ("yes", "cancer-testis antigen; these are non-neoplastic nPOD donors"),
    "HCN2":    ("yes", "HCN2 is neuronal as well as cardiac; peri-nodal nerve fibres are sparse but possible"),
}


def load(b: dict):
    ids, genes, counts = xenio._read_bundle_matrix(b["path"])
    return ids, pd.Index(genes), counts.tocsc()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------- B1: every metrics_summary.csv field ----------------
    rows = []
    for b in BUNDLES:
        m = pd.read_csv(b["path"] / "metrics_summary.csv").iloc[0].to_dict()
        exp = json.loads((b["path"] / "experiment.xenium").read_text())
        rec = {"donor_id": b["donor_id"], "donor_sex": b["donor_sex"], "serial": b["serial"],
               "design": b["design"], "roi": "pLN",
               "phenocycler_scan": b["phenocycler"],
               "recorded_in_xenium_paths_csv": b["recorded_in_xenium_paths"],
               "bundle_path": str(b["path"]),
               "panel_num_targets_predesigned": exp.get("panel_num_targets_predesigned"),
               "panel_num_targets_custom": exp.get("panel_num_targets_custom"),
               "instrument_sw_version": exp.get("instrument_sw_version"),
               "analysis_sw_version": exp.get("analysis_sw_version"),
               "run_start_time": exp.get("run_start_time")}
        rec.update(m)
        rows.append(rec)
    b1 = pd.DataFrame(rows)

    lead = ["donor_id", "donor_sex", "serial", "design", "roi", "run_name", "cassette_name",
            "region_name", "phenocycler_scan", "recorded_in_xenium_paths_csv",
            # decoding quality
            "fraction_transcripts_decoded_q20", "fraction_predesigned_transcripts_decoded_q20",
            "fraction_custom_transcripts_decoded_q20", "fraction_transcripts_assigned",
            "thickness_transcripts_high_quality",
            # control rates
            "adjusted_negative_control_probe_rate", "adjusted_negative_control_codeword_rate",
            "adjusted_genomic_control_probe_rate",
            "negative_control_probe_counts_per_control_per_cell",
            "genomic_control_probe_counts_per_control_per_cell",
            "estimated_number_of_false_positive_transcripts_per_cell",
            "estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts",
            # segmentation
            "stain_definition", "segmented_cell_stain_frac", "segmented_cell_boundary_frac",
            "segmented_cell_interior_frac", "segmented_cell_nuc_expansion_frac",
            "segmented_cell_boundary_count", "segmented_cell_interior_count",
            "segmented_cell_nuc_expansion_count", "fraction_empty_cells",
            # yield / density
            "num_cells_detected", "cells_per_100um2", "region_area", "total_cell_area",
            "total_high_quality_decoded_transcripts", "decoded_transcripts_per_100um2",
            "nuclear_transcripts_per_100um2", "median_genes_per_cell",
            "median_transcripts_per_cell", "median_predesigned_genes_per_cell",
            "median_predesigned_transcripts_per_cell", "median_custom_genes_per_cell",
            "median_custom_transcripts_per_cell"]
    ordered = [c for c in lead if c in b1.columns] + [c for c in b1.columns if c not in lead]
    b1 = b1[ordered]
    b1.to_csv(OUT / "B1_pln_metrics_summary_all_fields.csv", index=False)
    print(f"[B1] {len(b1)} bundles x {b1.shape[1]} fields -> B1_pln_metrics_summary_all_fields.csv",
          flush=True)

    # ---------------- load matrices, build the 5,001-gene intersection ----------------
    mats = {}
    gene_sets = []
    for b in BUNDLES:
        ids, genes, X = load(b)
        mats[b["serial"]] = (ids, genes, X)
        gene_sets.append(set(genes))
        print(f"[load] {b['serial']} {b['design']}: {X.shape[0]:,} cells x {len(genes):,} "
              f"gene features, nnz {X.nnz:,}", flush=True)
    inter = sorted(set.intersection(*gene_sets))
    union = sorted(set.union(*gene_sets))
    print(f"[panel] intersection {len(inter):,} genes; union {len(union):,}", flush=True)
    pd.DataFrame({"gene": inter}).to_csv(OUT / "B3_gene_intersection.csv", index=False)

    # ---------------- B3: per-gene detection + per-cell summaries ----------------
    per_gene = {"gene": inter}
    percell_rows = []
    for b in BUNDLES:
        s = b["serial"]
        ids, genes, X = mats[s]
        cols = pd.Index(genes).get_indexer(inter)
        assert (cols >= 0).all()
        Xi = X[:, cols].tocsr()
        n = Xi.shape[0]
        ge1 = np.asarray((Xi >= 1).sum(axis=0)).ravel()
        ge2 = np.asarray((Xi >= 2).sum(axis=0)).ravel()
        tot = np.asarray(Xi.sum(axis=0)).ravel()
        per_gene[f"n_cells_ge1_{s}"] = ge1
        per_gene[f"det_frac_ge1_{s}"] = ge1 / n
        per_gene[f"n_cells_ge2_{s}"] = ge2
        per_gene[f"det_frac_ge2_{s}"] = ge2 / n
        per_gene[f"total_counts_{s}"] = tot
        per_gene[f"counts_per_cell_{s}"] = tot / n

        depth = np.asarray(Xi.sum(axis=1)).ravel()
        ngene = np.asarray((Xi > 0).sum(axis=1)).ravel()
        percell_rows.append(dict(
            donor_id=b["donor_id"], serial=s, design=b["design"],
            n_cells=n, n_genes_intersection=len(inter),
            total_counts_on_intersection=int(depth.sum()),
            counts_per_cell_mean=float(depth.mean()),
            counts_per_cell_median=float(np.median(depth)),
            counts_per_cell_q25=float(np.percentile(depth, 25)),
            counts_per_cell_q75=float(np.percentile(depth, 75)),
            genes_per_cell_mean=float(ngene.mean()),
            genes_per_cell_median=float(np.median(ngene)),
            genes_per_cell_q25=float(np.percentile(ngene, 25)),
            genes_per_cell_q75=float(np.percentile(ngene, 75)),
            frac_cells_zero_counts=float((depth == 0).mean()),
            mean_det_frac_ge1_all_genes=float((ge1 / n).mean()),
            median_det_frac_ge1_all_genes=float(np.median(ge1 / n)),
            mean_det_frac_ge2_all_genes=float((ge2 / n).mean()),
            median_det_frac_ge2_all_genes=float(np.median(ge2 / n)),
        ))
        mats[s] = (ids, genes, Xi)          # keep the intersected CSR, drop the wide one
        print(f"[B3] {s} done", flush=True)

    pg = pd.DataFrame(per_gene)
    pg.to_csv(OUT / "B3_per_gene_detection.csv", index=False)
    pc = pd.DataFrame(percell_rows)
    pc.to_csv(OUT / "B3_per_cell_summary.csv", index=False)
    print("[B3] per-gene + per-cell written", flush=True)

    # ---------------- B2: implied negative controls ----------------
    audit = pd.read_csv(AUDIT)
    ann = audit[audit.category.notna()][["gene", "category", "exclude_recommended"]].copy()

    ctrl = []
    for g in Y_CLEAN:
        ctrl.append(dict(gene=g, set_name="Y_zero_template_clean", category="sex_chromosome"))
    for g in Y_CROSSHYB:
        ctrl.append(dict(gene=g, set_name="Y_zero_template_crosshyb", category="sex_chromosome"))
    for g in Y_OTHER:
        ctrl.append(dict(gene=g, set_name="Y_zero_template_other", category="sex_chromosome"))
    ydone = set(Y_CLEAN) | set(Y_CROSSHYB) | set(Y_OTHER)
    for _, r in ann.iterrows():
        if r.category == "sex_chromosome" or r.gene in ydone:
            continue
        ctrl.append(dict(gene=r.gene, set_name=f"audit::{r.category}", category=r.category))
    c = pd.DataFrame(ctrl).drop_duplicates("gene")
    c = c.merge(ann[["gene", "exclude_recommended"]], on="gene", how="left")

    verdict = [LN_VERDICT.get(g, ("yes", "")) for g in c.gene]
    c["ln_control_valid"] = [v[0] for v in verdict]
    c["ln_reason"] = [v[1] for v in verdict]
    # a gene the pancreas audit itself revised out of the exclusion list is not a control here
    revised = c.exclude_recommended.eq("no")
    c.loc[revised & c.ln_control_valid.eq("yes"), "ln_reason"] = (
        "panel_audit revised exclude yes->no on evidence of real expression; not a control in any tissue")
    c.loc[revised, "ln_control_valid"] = "no"
    # the Y probes carry no exclude_recommended reasoning of that kind
    c.loc[c.set_name.str.startswith("Y_"), "ln_control_valid"] = "yes"
    c.loc[c.set_name.str.startswith("Y_") & c.ln_reason.eq(""), "ln_reason"] = (
        "Y-linked; both donors are female so there is no template in any tissue")

    c["on_intersection"] = c.gene.isin(inter)
    for b in BUNDLES:
        s = b["serial"]
        m = pg.set_index("gene")
        c[f"det_frac_ge1_{s}"] = c.gene.map(m[f"det_frac_ge1_{s}"])
        c[f"det_frac_ge2_{s}"] = c.gene.map(m[f"det_frac_ge2_{s}"])
        c[f"counts_per_cell_{s}"] = c.gene.map(m[f"counts_per_cell_{s}"])
    c = c.sort_values(["set_name", "gene"])
    c.to_csv(OUT / "B2_control_genes_per_gene.csv", index=False)

    # per-set aggregates: ALL audited genes in the set, and the LN-valid subset
    agg = []
    for (sname), grp in c.groupby("set_name"):
        for scope, sub in (("all_audited", grp),
                           ("ln_valid_only", grp[grp.ln_control_valid.eq("yes")])):
            sub = sub[sub.on_intersection]
            row = dict(set_name=sname, scope=scope, n_genes=len(sub),
                       genes=";".join(sorted(sub.gene)))
            for b in BUNDLES:
                s = b["serial"]
                row[f"mean_det_frac_ge1_{s}"] = float(sub[f"det_frac_ge1_{s}"].mean()) if len(sub) else np.nan
                row[f"mean_det_frac_ge2_{s}"] = float(sub[f"det_frac_ge2_{s}"].mean()) if len(sub) else np.nan
                row[f"mean_counts_per_cell_{s}"] = float(sub[f"counts_per_cell_{s}"].mean()) if len(sub) else np.nan
            for donor, ss, ds in (("6533", "0080870", "0041393"), ("112", "0081366", "0041396")):
                a, bq = row[f"mean_det_frac_ge1_{ss}"], row[f"mean_det_frac_ge1_{ds}"]
                row[f"fold_sameslide_over_diffslide_ge1_{donor}"] = (a / bq) if (bq and bq > 0) else np.nan
            agg.append(row)
    ag = pd.DataFrame(agg).sort_values(["set_name", "scope"])
    ag.to_csv(OUT / "B2_control_sets_summary.csv", index=False)
    print("[B2] control tables written", flush=True)

    # headline: median gene vs control floors
    head = []
    for b in BUNDLES:
        s = b["serial"]
        med = float(np.median(pg[f"det_frac_ge1_{s}"]))
        row = dict(donor_id=b["donor_id"], serial=s, design=b["design"],
                   n_cells=int(pc.set_index("serial").loc[s, "n_cells"]),
                   median_gene_det_frac_ge1=med)
        for sname in ("Y_zero_template_clean", "Y_zero_template_crosshyb", "Y_zero_template_other"):
            sub = c[(c.set_name == sname) & c.on_intersection]
            row[sname] = float(sub[f"det_frac_ge1_{s}"].mean())
        lnv = c[(c.ln_control_valid == "yes") & c.on_intersection & (~c.set_name.str.startswith("Y_"))]
        row["ln_valid_autosomal_all_categories"] = float(lnv[f"det_frac_ge1_{s}"].mean())
        row["n_ln_valid_autosomal"] = int(len(lnv))
        head.append(row)
    hd = pd.DataFrame(head)
    hd.to_csv(OUT / "B2_headline_floors.csv", index=False)
    print(hd.to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
