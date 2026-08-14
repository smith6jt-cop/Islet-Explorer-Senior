#!/usr/bin/env python3
"""Build panel role-classification tables for the Xenium subpanels.

Emits two machine-readable artifacts under data/:

  panel_roles.csv       one row per (tissue, subpanel): the 3-tier role
                        (identity / subclass / state) and intended
                        phenotyping use (broad_mask / fine_tier /
                        module_score / embedding_extra / exclude).

  panel_gene_roles.csv  gene-level curation for the *mixed* panels (those
                        that bundle identity, subclass, and state genes):
                        each gene tagged identity / subclass / state with a
                        lineage hint and short note.

Rationale and literature citations live in docs/panel_classification.md.
This is a metadata generator (stdlib only); it reads the existing subpanel
CSVs so gene symbols are exact, and applies curated role lookups. Run from
the repo root:  python3 scripts/build_panel_roles.py

(The repo's data-build scripts are otherwise R; this one is Python because
the classification is hand-curated metadata, not a numeric recompute, and
no R runtime is assumed where it was authored.)
"""
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SHARED = os.path.join(DATA, "subpanels_shared")
TISSUES = {
    "pancreas": os.path.join(DATA, "tissues", "pancreas", "subpanels"),
    "thymus":   os.path.join(DATA, "tissues", "thymus", "subpanels"),
}


def n_genes(path):
    """Count data rows (non-empty, minus header) of a subpanel CSV."""
    with open(path, newline="") as fh:
        rows = [ln for ln in fh.read().splitlines() if ln.strip()]
    return max(0, len(rows) - 1)


def gene_col(path):
    """Return the ordered list of gene symbols (first column, no header)."""
    out = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd, None)
        for row in rd:
            if row and row[0].strip():
                out.append(row[0].strip())
    return out


# ---------------------------------------------------------------------------
# Panel-level roles.  tuple = (tier, phenotyping_use, lineage_or_group, source, note)
#   source: "shared" => consumer reads subpanels_shared/<stem>.csv for the mask
# Stems not listed fall through to defaults (pathway 21-49 / residual 99*).
# ---------------------------------------------------------------------------
M = "broad_mask"; F = "fine_tier"; S = "module_score"; E = "embedding_extra"; X = "exclude"
ID = "identity"; SUB = "subclass"; ST = "state"; RES = "residual"; CUS = "custom"

PANEL_ROLES = {
"pancreas": {
 "01_pancreas_endocrine": (ID, M, "Endocrine", None, ""),
 "02_pancreas_exocrine": (ID, M, "Exocrine", None, ""),
 "03_immune_T_cell": (ID, M, "Immune_T", "shared", "consumer reads shared/ copy"),
 "04_immune_B_plasma": (ID, M, "Immune_B", "shared", "consumer reads shared/ copy"),
 "05_immune_NK_ILC": (ID, M, "Immune_T", "shared", "NK/ILC; allowlist unions into Immune_T; consumer reads shared/"),
 "06_immune_myeloid": (ID, M, "Myeloid", "shared", "consumer reads shared/ copy"),
 "07_immune_granulocyte": (ID, M, "Myeloid", "shared", "consumer reads shared/ copy"),
 "08_antigen_presentation_MHC": (ST, X, "NK_T_recognition_MHCI", None,
   "MISLABELED: NK receptors + CD4/8 + MHC-I machinery, no HLA-DR/CD74; identity genes already in 03/05/06; ~65% overlap with T mask; see panel_gene_roles.csv"),
 "09_cytokines_chemokines": (ST, S, "Cytokine_signaling", None,
   "immune-biased communication program; ~45% receptors (some lineage-informative -> fine_tier); see gene roles"),
 "10_interferon_antiviral": (ST, S, "IFN_response", None,
   "canonical inducible ISG/JAK-STAT state (Schoggins 2019); IFNA* = pDC-identity candidates"),
 "11_endothelial_vascular": (ID, M, "Vascular", "shared", "consumer reads shared/ copy"),
 "12_pericyte_smooth_muscle": (ID, M, "Stromal", "shared", "consumer reads shared/ copy"),
 "13_fibroblast_stellate": (ID, M, "Stromal", None, ""),
 "14_ECM_remodeling": (ST, S, "Matrix_remodeling", None,
   "matrisome; mesenchymal-leaning but overlaps 12/13 + basement-membrane (epi/endo) + broad integrins; secreted subset = stromal embedding-extra candidate"),
 "15_neural_glial": (ID, M, "Neural", "shared", "consumer reads shared/ copy"),
 "16_epithelial_general": (ID, E, "Epithelial", "shared",
   "overlaps endocrine/exocrine/ductal; used as embedding-extra (gene union), NOT a mask"),
 "17_cell_cycle_proliferation": (ST, S, "Proliferation", None, "cross-lineage; consumer tracks via CYCLING_MARKERS"),
 "18_apoptosis_cell_death": (ST, S, "Cell_death", None, "cross-lineage program"),
 "19_ER_stress_UPR": (ST, S, "ER_stress", None, "cross-lineage program"),
 "20_signaling_pathways_developmental": (ST, S, "Dev_signaling", None, "WNT/NOTCH/TGFb/BMP/FGF/Hedgehog programs"),
},
"thymus": {
 "01_thymic_epithelial_cortical": (ID, M, "TEC", None, ""),
 "02_thymic_epithelial_medullary": (ID, M, "TEC", None, "ensure AIRE/FEZF2 present (mTEC identity)"),
 "03_hassall_keratinized": (ID, M, "TEC", None, ""),
 "04_thymic_specialized_TEC": (ID, M, "TEC", None, ""),
 "05_thymic_myoid": (ID, M, "Stromal", None, "degraded (3 genes) per user exclusions"),
 "06_thymocyte_DN": (ID, M, "Thymocyte", None, ""),
 "07_thymocyte_DP": (ID, M, "Thymocyte", None, ""),
 "08_thymocyte_SP": (ID, M, "Thymocyte", None, ""),
 "09_thymic_treg": (ID, M, "Thymocyte", None, ""),
 "10_TCR_VDJ_machinery": (SUB, F, "Thymocyte_immature", None, "RAG1/2/DNTT/PTCRA = DN/DP subclass machinery; not a new lineage"),
 "11_positive_negative_selection": (ST, S, "Selection", None,
   "apoptosis + MHC machinery + NR4A; AIRE/FEZF2 = mTEC identity -> promote; see gene roles"),
 "12_thymic_DC_myeloid": (ID, M, "Myeloid", None, ""),
 "13_thymic_B_plasma": (ID, M, "B_lymphoid", None, ""),
 "14_thymic_NK_ILC": (ID, M, "NK", None, ""),
 "15_thymic_stroma_fibroblast": (ID, M, "Stromal", None, ""),
 "16_thymic_endothelial_pericyte": (ID, M, "Vascular", None, ""),
 "17_notch_il7_signaling": (SUB, F, "Thymocyte_DN", None,
   "IL7R/KIT/FLT3/NOTCH1 = DN-stage subclass; mixed with Notch-effector state; see gene roles"),
 "18_crosstalk_RANKL_LTBR": (ST, S, "Crosstalk_NFkB", None, "TEC<->thymocyte crosstalk; NF-kB program"),
 "19_AIRE_TRA_program": (ID, F, "TEC_mTEC", None,
   "AIRE/FEZF2 = mTEC identity (promote to mask 02/04); TRA reporters = mTEC promiscuous expression; degraded (7 genes)"),
 "20_thymic_egress_S1P": (SUB, F, "Thymocyte_SP_egress", None,
   "S1PR1/KLF2/SELL/CCR7 = SP maturation/egress; ZBTB16 = innate-T"),
 "50_custom_100_addon": (CUS, X, "Custom_T1D_GWAS", None, "standalone add-on; not a phenotyping mask"),
},
}

# Extra mask rows the consumer sources from shared/ that have no tissue-folder file.
SHARED_APPEND = {
 "thymus": [
   ("15_neural_glial", ID, M, "Neural", "sourced from subpanels_shared/ (no thymus-specific neural panel)"),
   ("16_epithelial_general", ID, E, "Epithelial", "embedding-extra; sourced from subpanels_shared/"),
 ],
}

PATHWAY_RE = re.compile(r"^(2[1-9]|3[0-9]|4[0-9])_")


def panel_default(stem):
    if PATHWAY_RE.match(stem):
        return (ST, S, "Pathway", None,
                "functional module; pathway panels overlap heavily by design; key by name not number across tissues")
    if stem.startswith("99"):
        return (RES, X, "Residual", None, "uncategorized residual; nested 99 superset of 99c superset of 99d")
    return ("UNCLASSIFIED", "review", "?", None, "no role rule matched -- review")


def build_panel_roles(n_identity=None):
    """One row per (tissue, subpanel).  `n_identity` maps (source, stem) ->
    identity_core size (from derive_identity_panels); blank where not a
    cell-type identity panel."""
    n_identity = n_identity or {}
    rows = []
    for tissue, sdir in TISSUES.items():
        stems = sorted(s[:-4] for s in os.listdir(sdir) if s.endswith(".csv"))
        for stem in stems:
            tier, use, grp, src, note = PANEL_ROLES.get(tissue, {}).get(stem) or panel_default(stem)
            resolved = "shared" if src == "shared" else tissue
            read_from = os.path.join(SHARED, stem + ".csv") if src == "shared" else os.path.join(sdir, stem + ".csv")
            srcpath = ("subpanels_shared/%s.csv" % stem) if src == "shared" \
                      else ("data/tissues/%s/subpanels/%s.csv" % (tissue, stem))
            rows.append([tissue, stem, tier, use, grp,
                         n_genes(read_from if os.path.exists(read_from) else os.path.join(sdir, stem + ".csv")),
                         n_identity.get((resolved, stem), ""),
                         srcpath, note])
        for stem, tier, use, grp, note in SHARED_APPEND.get(tissue, []):
            rows.append([tissue, stem, tier, use, grp,
                         n_genes(os.path.join(SHARED, stem + ".csv")),
                         n_identity.get(("shared", stem), ""),
                         "subpanels_shared/%s.csv" % stem, note])
    out = os.path.join(DATA, "panel_roles.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tissue", "subpanel", "tier", "phenotyping_use",
                    "lineage_or_group", "n_genes", "n_identity_genes",
                    "phenotyping_source", "note"])
        w.writerows(rows)
    return rows


# ---------------------------------------------------------------------------
# Gene-level curation for the MIXED panels.
# Each panel: a per-gene override dict + a (tier, hint, note) DEFAULT applied
# to any gene not explicitly listed.  Genes are read from the actual CSV so
# symbols are exact; an override for a gene absent from the file is ignored.
# ---------------------------------------------------------------------------
MIXED = {
 ("shared", "08_antigen_presentation_MHC", os.path.join(SHARED, "08_antigen_presentation_MHC.csv")): {
   "default": (ST, "MHC-I/stress", "MHC-I peptide-loading / stress ligand"),
   "genes": {
     "CD4": (SUB, "T_CD4", "T-helper coreceptor (also low myeloid); already in 03_immune_T_cell"),
     "CD8A": (SUB, "T_CD8/NK", "CD8 coreceptor; already in 03"),
     "CD8B": (SUB, "T_CD8", "CD8 coreceptor; already in 03"),
     "KIR2DL1": (ID, "NK", "inhibitory KIR; NK identity; already in 05_immune_NK_ILC"),
     "KIR3DL1": (ID, "NK", "inhibitory KIR; NK identity; already in 05"),
     "KLRC1": (ID, "NK/CD8", "NKG2A; NK/CD8 identity"),
     "KLRC2": (ID, "NK", "NKG2C; NK identity"),
     "KLRK1": (SUB, "NK/CD8/gd", "NKG2D; cytotoxic-lymphocyte receptor"),
     "LILRB1": (SUB, "NK/myeloid", "inhibitory receptor; NK/myeloid"),
     "LILRB2": (SUB, "myeloid", "myeloid inhibitory receptor"),
     "CIITA": (ST, "APC/MHC-II", "MHC-II master regulator: constitutive in APC, IFN-g-inducible elsewhere (Baton 2004)"),
     "NLRC5": (ST, "MHC-I", "MHC-I transactivator; IFN-inducible, broad"),
     "TAP1": (ST, "MHC-I", "MHC-I peptide transport; IFN-inducible, broad"),
     "TAP2": (ST, "MHC-I", "MHC-I peptide transport; IFN-inducible, broad"),
     "MICA": (ST, "stress-ligand", "stress-induced NKG2D ligand; inducible across cells"),
     "MICB": (ST, "stress-ligand", "stress-induced NKG2D ligand"),
     "ULBP2": (ST, "stress-ligand", "stress-induced NKG2D ligand"),
   }},
 ("shared", "09_cytokines_chemokines", os.path.join(SHARED, "09_cytokines_chemokines.csv")): {
   "default": (ST, "cytokine-signaling", "inducible cytokine/chemokine ligand or broad receptor"),
   "genes": {
     "XCR1": (ID, "cDC1", "cDC1 identity receptor"),
     "CXCR1": (ID, "neutrophil", "neutrophil identity"),
     "CXCR2": (ID, "neutrophil", "neutrophil identity"),
     "IL3RA": (ID, "pDC/basophil", "CD123; pDC/basophil"),
     "IL5RA": (ID, "eosinophil", "eosinophil identity"),
     "CCR2": (SUB, "monocyte", "inflammatory monocyte"),
     "CX3CR1": (SUB, "monocyte/NK", "patrolling monocyte / NK / Trm"),
     "CCR7": (SUB, "naive/cm-T;mDC", "naive/central-memory T, mature DC"),
     "CXCR5": (SUB, "Tfh/B", "follicular T / B homing"),
     "CXCR3": (SUB, "Th1/CD8eff", "Th1 / effector CD8"),
     "CXCR6": (SUB, "Trm/NKT", "tissue-resident T / NKT"),
     "CCR9": (SUB, "gut-T/pDC", "gut-homing T / pDC"),
     "CCR6": (SUB, "Th17", "Th17 / some DC,B"),
     "CCR4": (SUB, "Th2/Treg", "Th2 / Treg / skin T"),
     "CCR3": (SUB, "eosinophil", "eosinophil / basophil"),
     "CCR5": (SUB, "Th1/macrophage", "Th1 / macrophage"),
     "CCR8": (SUB, "Treg/Th2", "Treg / Th2"),
     "IL2RA": (SUB, "Treg/actT", "CD25; Treg / activated T"),
     "IL7R": (SUB, "T/ILC", "CD127; T / ILC"),
     "IL1RL1": (SUB, "ILC2/Th2/mast", "ST2; ILC2 / Th2 / mast"),
     "IL23R": (SUB, "Th17/ILC3", "Th17 / ILC3"),
     "IL17A": (SUB, "Th17", "Th17 effector cytokine"),
     "IL17F": (SUB, "Th17", "Th17 effector cytokine"),
     "IL4": (SUB, "Th2/ILC2", "Th2 / ILC2 cytokine"),
     "IL5": (SUB, "Th2/ILC2", "Th2 / ILC2 cytokine"),
     "IL13": (SUB, "Th2/ILC2", "Th2 / ILC2 cytokine"),
     "IL21": (SUB, "Tfh", "Tfh cytokine"),
     "IL22": (SUB, "Th17/ILC3", "Th17 / ILC3 cytokine"),
     "CXCL13": (SUB, "FDC/stroma", "follicular dendritic / lymphoid stroma"),
     "CCL19": (SUB, "FRC/stroma", "fibroblastic reticular cell / lymphoid stroma"),
     "CCL17": (SUB, "DC", "DC-derived"),
     "CCL22": (SUB, "DC", "DC-derived"),
     "TSLP": (SUB, "epithelial", "epithelial alarmin"),
     "IL33": (SUB, "epithelial/stroma", "epithelial/stromal alarmin"),
   }},
 ("shared", "10_interferon_antiviral", os.path.join(SHARED, "10_interferon_antiviral.csv")): {
   "default": (ST, "IFN-response", "ISG / sensor / JAK-STAT: inducible across nucleated cells"),
   "genes": {
     "IFNA1": (ID, "pDC(cand)", "type-I IFN ligand; pDC-biased production"),
     "IFNA2": (ID, "pDC(cand)", "type-I IFN ligand; pDC-biased"),
     "IFNA7": (ID, "pDC(cand)", "type-I IFN ligand; pDC-biased"),
     "IFNA8": (ID, "pDC(cand)", "type-I IFN ligand; pDC-biased"),
     "IFNA17": (ID, "pDC(cand)", "type-I IFN ligand; pDC-biased"),
     "IFNG": (SUB, "Th1/NK/CD8", "IFN-gamma production = Th1/NK/effector-CD8"),
     "IFNL1": (ST, "epithelial-IFN", "type-III IFN; epithelial-biased response"),
     "IFNL2": (ST, "epithelial-IFN", "type-III IFN; epithelial-biased response"),
     "IFNL3": (ST, "epithelial-IFN", "type-III IFN; epithelial-biased response"),
   }},
 ("shared", "14_ECM_remodeling", os.path.join(SHARED, "14_ECM_remodeling.csv")): {
   "default": (ST, "Stromal-matrix", "matrisome / remodeling; mesenchymal-enriched but secreted/shared"),
   "genes": {
     "ITGA2B": (ID, "megakaryocyte", "CD41; platelet/megakaryocyte identity"),
     "GP6": (ID, "platelet", "platelet collagen receptor"),
     "ITGB3": (SUB, "platelet/endo", "platelet / endothelial"),
     "ITGB7": (SUB, "lymphocyte", "gut-homing lymphocyte integrin"),
     "ITGB4": (SUB, "epithelial", "hemidesmosome; epithelial"),
     "ITGB6": (SUB, "epithelial", "epithelial integrin"),
     "ITGA6": (SUB, "epithelial", "epithelial/basal integrin"),
     "MMP8": (ID, "neutrophil", "neutrophil collagenase"),
     "MMP9": (SUB, "neutrophil/macrophage", "gelatinase; granulocyte/macrophage"),
     "MMP12": (SUB, "macrophage", "macrophage metalloelastase"),
     "CD36": (SUB, "macrophage/endo", "scavenger receptor; macrophage/endothelial/platelet"),
     "SV2C": (ST, "neural(odd)", "synaptic vesicle protein; off-target for an ECM panel"),
   }},
 ("thymus", "11_positive_negative_selection", os.path.join(TISSUES["thymus"], "11_positive_negative_selection.csv")): {
   "default": (ST, "selection-process", "selection-associated process gene"),
   "genes": {
     "AIRE": (ID, "mTEC", "medullary TEC identity TF -> promote to TEC mask"),
     "FEZF2": (ID, "mTEC", "medullary TEC identity TF -> promote to TEC mask"),
     "CIITA": (ST, "APC/MHC-II", "MHC-II master regulator"),
     "CD74": (SUB, "APC(B/DC)", "MHC-II invariant chain; B/DC/APC"),
     "B2M": (ST, "MHC-I", "MHC-I light chain; broad/IFN-inducible"),
     "ERAP1": (ST, "MHC-I", "MHC-I peptide trimming"),
     "PSMB8": (ST, "immunoproteasome", "IFN-inducible immunoproteasome"),
     "PSMB9": (ST, "immunoproteasome", "IFN-inducible immunoproteasome"),
     "PSMB10": (ST, "immunoproteasome", "IFN-inducible immunoproteasome"),
     "TAP1": (ST, "MHC-I", "MHC-I transport"),
     "TAP2": (ST, "MHC-I", "MHC-I transport"),
     "BAD": (ST, "apoptosis", "apoptosis; cross-lineage"),
     "BAX": (ST, "apoptosis", "apoptosis; cross-lineage"),
     "BCL2": (ST, "apoptosis", "apoptosis; cross-lineage"),
     "BCL2L11": (ST, "apoptosis", "BIM; negative-selection apoptosis"),
     "CASP3": (ST, "apoptosis", "apoptosis; cross-lineage"),
     "CASP9": (ST, "apoptosis", "apoptosis; cross-lineage"),
     "CYCS": (ST, "apoptosis", "cytochrome c; apoptosis"),
     "HRK": (ST, "apoptosis", "apoptosis"),
     "PMAIP1": (ST, "apoptosis", "NOXA; apoptosis"),
     "NR4A1": (ST, "TCR-activation", "immediate-early TCR signal / negative selection"),
     "NR4A3": (ST, "TCR-activation", "immediate-early TCR signal / negative selection"),
   }},
 ("thymus", "17_notch_il7_signaling", os.path.join(TISSUES["thymus"], "17_notch_il7_signaling.csv")): {
   "default": (ST, "Notch-signaling", "Notch/IL7 signaling-effector or housekeeping"),
   "genes": {
     "IL7R": (SUB, "DN/early-T", "early-thymocyte / T / ILC"),
     "KIT": (SUB, "ETP/DN1", "early thymic progenitor"),
     "FLT3": (SUB, "progenitor/DC", "progenitor / DC"),
     "NOTCH1": (SUB, "T-commitment", "T-lineage commitment (DN)"),
     "DLL1": (SUB, "TEC-niche", "Notch ligand; TEC/stroma niche"),
     "DLL4": (SUB, "TEC-niche", "Notch ligand; cTEC niche"),
     "JAG2": (SUB, "TEC-niche", "Notch ligand; TEC niche"),
     "CXCL12": (ST, "niche/stroma", "niche chemokine; stromal"),
   }},
 ("thymus", "19_AIRE_TRA_program", os.path.join(TISSUES["thymus"], "19_AIRE_TRA_program.csv")): {
   "default": (SUB, "TRA-mTEC", "tissue-restricted antigen; mTEC promiscuous expression readout"),
   "genes": {
     "AIRE": (ID, "mTEC", "medullary TEC identity TF -> promote to TEC mask"),
     "FEZF2": (ID, "mTEC", "medullary TEC identity TF -> promote to TEC mask"),
     "CIITA": (ST, "APC/MHC-II", "MHC-II master regulator"),
   }},
 ("thymus", "20_thymic_egress_S1P", os.path.join(TISSUES["thymus"], "20_thymic_egress_S1P.csv")): {
   "default": (SUB, "SP-egress", "SP-thymocyte maturation/egress state"),
   "genes": {
     "CD69": (ST, "activation", "activation/retention marker"),
     "FOXO1": (ST, "naive-TF", "naive/quiescence TF"),
     "ZBTB16": (SUB, "innate-T", "PLZF; innate-like T / NKT"),
     "EBI3": (ST, "regulatory", "IL-27/IL-35 subunit"),
   }},
}


def build_gene_roles_mixed(ct_terms):
    """Mixed/program-panel gene curation, upgraded to the unified 10-column
    schema used by derive_identity_panels.  These panels (08/09/10/14, thymus
    11/17/19/20) bundle a few lineage markers among program genes and were
    classified by hand from the literature, so evidence='canonical'; cl_terms
    is the gene's reference cell_type (passed in), kept for context."""
    rows = []
    for (source, stem, path), spec in MIXED.items():
        dtier, dhint, dnote = spec["default"]
        for g in gene_col(path):
            tier, hint, note = spec["genes"].get(g, (dtier, dhint, dnote))
            cl_terms = ";".join(ct_terms.get(g, [])[:6])
            rows.append([source, stem, g, tier, hint, "", "", cl_terms,
                         "canonical", note])
    return rows


if __name__ == "__main__":
    import sys
    from collections import Counter
    # Import the cell-type identity classifier lazily (sibling module in scripts/)
    # so importing build_panel_roles as a library has no ontology-parsing side
    # effect and there is no import-time cycle.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import derive_identity_panels as derive

    # 1. CL-based classification of every cell-type identity panel
    #    (+ identity_core/, identity_audit_*, identity_epithelial_general_*).
    res = derive.run(write_outputs=True)
    # 2. panel_roles.csv, now carrying n_identity_genes.
    pr = build_panel_roles(res["n_identity"])
    # 3. unified panel_gene_roles.csv = hand-curated mixed rows + cell-type rows.
    ct_terms, _ = derive.load_reference()
    mixed = build_gene_roles_mixed(ct_terms)
    derive.write_panel_gene_roles(res["rows"], mixed)

    # --- validation summary ---
    print("CL release:", derive.cl_data_version())
    print("panel_roles.csv: %d rows" % len(pr))
    for t in ("pancreas", "thymus"):
        c = Counter(r[3] for r in pr if r[0] == t)
        masks = sorted(r[4] for r in pr if r[0] == t and r[3] == "broad_mask")
        print("  %-8s use-counts: %s" % (t, dict(c)))
        print("           broad_mask lineages: %s" % masks)
        bad = [r[1] for r in pr if r[0] == t and r[2] == "UNCLASSIFIED"]
        if bad:
            print("  [WARN] unclassified in %s: %s" % (t, bad))
    tiers = Counter(r[3] for r in res["rows"])
    print("panel_gene_roles.csv: %d mixed + %d cell-type = %d rows"
          % (len(mixed), len(res["rows"]), len(mixed) + len(res["rows"])))
    print("  cell-type tier counts:", dict(tiers))
