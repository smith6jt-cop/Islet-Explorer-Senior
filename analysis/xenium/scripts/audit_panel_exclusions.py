"""Audit the Xenium 5K + custom panel for genes that are biologically
implausible in the target tissue — candidates to drop from cell-type
scoring panels and (optionally) from the analytical feature set entirely.

CLI:
    --tissue {pancreas, thymus}   (default: pancreas)
    --samples ID1 ID2 ...         (default: 0041323 0041326)
    --output PATH                 (default: data/processed/panel_audit.csv)
    --staging-dir DIR             (default: /tmp; reads
                                   panel_genes_union.csv +
                                   panel_detection_{sample}.csv from here)

Inputs (in --staging-dir):
    panel_genes_union.csv          union of genes across all samples
    panel_detection_{sample}.csv   per-sample n_cells / detection_pct

Output:
    panel_audit{,_thymus}.csv      gene-by-gene flagging with category,
                                     rationale, recommendation, and
                                     per-sample detection stats

Rules used in this audit:
    - Only flag genes confidently restricted to tissues not present in
      the target tissue. For pancreas this means acinar, ductal,
      endocrine, stellate, endothelial, smooth muscle, pericyte,
      fibroblast, Schwann, immune (and lymph-node residents). For
      thymus this means cTEC/mTEC (incl. Hassall's corpuscles), the
      thymocyte maturation hierarchy, plus the same stromal/vascular/
      immune scaffolding.
    - Generic stromal/epithelial markers, broadly-expressed
      transcription factors, and any gene with even a single plausible
      immune or stromal role is KEPT (better to under-flag than
      over-flag).
    - Empirical sample-asymmetry (>10× detection-rate difference
      between samples on a normalized basis) is flagged as a separate
      signal — panel cross-detection or quality-failure candidates,
      not necessarily biologically implausible.

Tissue-specific deviations from the universal set:
    pancreas: keeps the original `skin_cornification` block including
              the Hassall's-corpuscle-compatible keratins (KRT1/5/10/14,
              FLG/IVL/LOR/TGM3 etc.) because they are stratified-skin-
              specific in the pancreas context.
    thymus:   moves those Hassall's-compatible genes OUT of exclusion
              (medullary thymic epithelium genuinely expresses them);
              adds pancreatic_exocrine, pancreatic_endocrine, and
              pancreatic_islet_TF categories.

The GHRL `SINGLE_GENE_PANEL_REMOVE` note is pancreas-specific (Epsilon
scoring artefact); it is skipped for non-pancreas tissue.
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))

# =============================================================================
# Curated tissue-restricted marker sets — high-confidence exclusions
# =============================================================================

EXCLUDE_CATEGORIES = {

    # --- Photoreceptor / retinal ---
    "photoreceptor_retinal": {
        "rationale": "Restricted to retinal photoreceptors / RPE; not present "
                     "in pancreas or LN.",
        "genes": {"RHO", "OPN1SW", "OPN1MW", "OPN1LW", "GRK1", "GRK7", "RCVRN",
                   "ARR3", "GUCA1A", "GUCA1B", "NRL", "CRX", "NR2E3", "RPGR",
                   "ABCA4", "USH2A", "RP1", "RP2", "ROM1", "PDE6A", "PDE6B",
                   "PDE6G", "PDE6H", "TULP1", "AIPL1", "IMPG1", "IMPG2",
                   "FSCN2", "RD3", "GUCY2D", "GUCY2F", "KCNV2", "RBP3",
                   "PRPH2", "SAG", "GNAT1", "GNAT2", "CNGA1", "CNGB1",
                   "CNGA3", "CNGB3", "RGR", "RPE65", "LRAT", "RDH5",
                   "BEST1", "VSX2", "OTX2"},  # OTX2 also some non-retinal but mostly retina
    },

    # --- Olfactory + taste receptors ---
    "chemosensory_receptor": {
        "rationale": "Olfactory / taste receptors; nasal/oral epithelium only.",
        "genes": set(),  # populated by prefix-match below
        "prefix_match": ["OR1", "OR2", "OR3", "OR4", "OR5", "OR6", "OR7",
                          "OR8", "OR9", "OR10", "OR11", "OR12", "OR13", "OR14",
                          "OR51", "OR52", "OR56", "OR6Y1", "OR8B",
                          "TAS1R", "TAS2R"],
    },

    # --- Sperm / spermatogenesis-restricted ---
    "spermatogenesis": {
        "rationale": "Spermatogenic / testis-only expression; absent from "
                     "pancreas/LN tissues.",
        "genes": {"PRM1", "PRM2", "PRM3", "TNP1", "TNP2", "AKAP4", "ODF1",
                   "ODF2", "ODF3", "ODF4", "OAZ3", "CRISP2", "ZPBP", "ZPBP2",
                   "SPATA1", "SPATA3", "SPATA4", "SPATA8", "SPATA16", "SPATA17",
                   "SPATA22", "TEX14", "TEX15", "TEX19", "LDHC", "PGK2", "RFX2",
                   "BOLL", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "STK31", "DPY19L2",
                   "PRSS21", "TSGA10", "FATE1", "MAGEA1", "MAGEA3", "MAGEA4",
                   "MAGEA10", "MAGEA12", "PIWIL1", "PIWIL2", "PIWIL4",
                   "SYCP1", "SYCP2", "SYCP3", "MLH3", "HORMAD1", "HORMAD2",
                   "STAG3", "MEI1", "MEIOB", "MEIOC", "DMC1", "REC8",
                   "DDX25", "TDRD1", "TDRD5", "TDRD6", "TDRD9"},
    },

    # --- Oocyte / female germline ---
    "oocyte_germline": {
        "rationale": "Oocyte / female germline-specific expression.",
        "genes": {"ZP1", "ZP2", "ZP3", "ZP4", "GDF9", "BMP15", "DAZL",
                   "DDX4", "FIGLA", "NOBOX", "OOSP1", "OOSP2", "OOEP",
                   "LHX8", "PADI6", "TLE6", "KHDC3L", "NLRP5", "NLRP7",
                   "NLRP14", "JAG1", "ZAR1", "MOS", "H1FOO", "H1-7"},
    },

    # --- Embryonic pluripotency-restricted ---
    "embryonic_pluripotency": {
        "rationale": "Pluripotency factors / embryonic-restricted; not in "
                     "adult pancreas or LN.",
        "genes": {"POU5F1", "NANOG", "ZFP42", "LEFTY1", "LEFTY2", "DPPA3",
                   "DPPA5", "TDGF1", "DNMT3L", "PRDM14", "TFAP2C"},
        # SOX2 deliberately omitted (some adult progenitor expression)
    },

    # --- Skeletal-muscle-specific contractile / regulatory ---
    "skeletal_muscle": {
        "rationale": "Skeletal-muscle-specific contractile machinery; not "
                     "expressed in pancreas/LN smooth muscle or "
                     "myofibroblasts.",
        "genes": {"MYH1", "MYH2", "MYH4", "MYH8", "MYH13", "MYH15",
                   "MYL1", "MYBPC1", "MYBPC2", "NEB", "TNNI2", "TNNT3",
                   "TNNC2", "ACTA1", "ACTN2", "ACTN3", "MYOD1", "MYF5",
                   "MYF6", "MYOG", "TTN",  # titin: strong skel/cardiac
                   "RYR1",  # skeletal muscle ryanodine receptor
                   "DMD",  # dystrophin: skeletal/cardiac muscle
                   "OBSCN", "OBSL1", "MYOM1", "MYOM2", "MYOZ1", "MYOZ2",
                   "MYOZ3", "CASQ1", "PVALB",  # parvalbumin (skel/CNS interneurons)
                   "MYBPH", "TRIM63", "FBXO32",  # atrogenes
                   "MUSK", "DOK7"},  # NMJ
        # DES (desmin) deliberately kept — also myofibroblast/stellate
        # ACTA2 (smooth muscle alpha actin) kept
    },

    # --- Cardiac-muscle-specific ---
    "cardiac_muscle": {
        "rationale": "Cardiac-muscle-specific; not in pancreas/LN.",
        "genes": {"MYH6", "MYH7", "MYH7B",  # cardiac heavy chains
                   "NPPA", "NPPB", "NPPC",  # natriuretic peptides
                   "TNNT2", "TNNI3",  # cardiac troponins
                   "TNNC1",  # cardiac/slow skeletal — borderline
                   "ACTC1",  # cardiac alpha actin
                   "MYL2", "MYL7", "MYL3", "MYL4",  # cardiac myosin LC
                   "MYBPC3",  # cardiac MyBP-C
                   "HCN4", "HCN2",  # pacemaker
                   "KCNH2", "SCN5A", "KCNQ1",  # cardiac ion channels
                   "GJA5",  # connexin-40 (atrial)
                   "RYR2",  # cardiac ryanodine receptor
                   "PLN",  # phospholamban
                   "CASQ2",  # cardiac calsequestrin
                   "TBX5",  # cardiac TF
                   "NKX2-5",  # cardiac TF
                   "ANKRD1",  # cardiac ankyrin
                   "BMP10"},  # cardiac BMP
    },

    # --- Hepatocyte-restricted ---
    "hepatocyte": {
        "rationale": "Hepatocyte-specific; absent from pancreatic "
                     "exocrine/endocrine and lymphoid tissues.",
        "genes": {"ALB", "FGA", "FGB", "FGG",  # albumin, fibrinogen
                   "AFP",  # alpha-fetoprotein (fetal liver)
                   "F2", "F7", "F9", "F11", "F12", "F13A1", "F13B",  # clotting
                   "APOH",  # beta-2-glycoprotein
                   "ITIH1", "ITIH2", "ITIH3", "ITIH4",  # inter-α-trypsin
                   "SERPINC1", "SERPIND1",  # antithrombin, heparin cofactor II
                   "HPX",  # hemopexin
                   "HRG",  # histidine-rich glycoprotein
                   "AHSG",  # fetuin-A
                   "PROC", "PROS1",  # protein C, S
                   "AGT",  # angiotensinogen (hepatic)
                   "C8A", "C8B", "C9",  # complement (hepatic)
                   "CYP7A1", "CYP8B1",  # bile acid synthesis
                   "BAAT",  # bile acid CoA:amino acid N-acyltransferase
                   "G6PC1",  # glucose-6-phosphatase liver
                   "PCK1",  # PEPCK liver-restricted variant — borderline (some kidney)
                   "TAT",  # tyrosine aminotransferase (liver)
                   "HAO1",  # hydroxyacid oxidase 1 (liver)
                   "BHMT",  # betaine-homocysteine methyltransferase
                   "ARG1"},  # arginase 1 (liver, also some myeloid)
        # Deliberately KEPT: TTR (transthyretin — also choroid plexus and some
        # endocrine), HNF4A (also pancreas/intestine), APOB/APOA1/APOA2/APOA4
        # (also intestine), CYP3A4 (also intestine).
    },

    # --- Bone / hard tissue ---
    "bone_cartilage_tooth": {
        "rationale": "Bone matrix, hypertrophic chondrocyte, or tooth-specific.",
        "genes": {"BGLAP",  # osteocalcin
                   "SP7",  # osterix
                   "IBSP",  # bone sialoprotein
                   "COL10A1",  # hypertrophic chondrocyte collagen
                   "MEPE", "DMP1",  # bone/SIBLINGs
                   "PHEX",  # phosphate metab, bone
                   "ACAN",  # aggrecan (cartilage)
                   "MATN1", "MATN3",  # cartilage matrilins
                   "COMP",  # cartilage oligomeric matrix protein
                   "COL11A1", "COL11A2",  # cartilage collagens
                   "AMELX", "AMELY", "ENAM", "AMBN",  # tooth enamel
                   "MMP20",  # tooth-specific metalloprotease
                   "DSPP",  # dentin
                   "ODAM",  # odontogenic
                   "COL2A1"},  # type II collagen — mainly cartilage; some embryo
    },

    # --- Skin keratinization — pancreas variant (hard cornification +
    # stratified keratins; for pancreas only, where stratified epithelium
    # is absent). In the thymus variant below, the Hassall's-corpuscle-
    # compatible subset is removed.
    "skin_cornification": {
        "rationale": "Cornified envelope / hard skin keratinization; not in "
                     "pancreas ductal or LN.",
        "genes": {"FLG", "FLG2",  # filaggrin
                   "IVL",  # involucrin
                   "LOR",  # loricrin
                   "HRNR",  # hornerin
                   "TGM3",  # cornification transglutaminase
                   "DSG1", "DSG3", "DSG4",  # desmoglein 1/3/4 (skin layers)
                   "DSC1", "DSC3",  # desmocollin 1/3 (skin)
                   "KRT1", "KRT2", "KRT9", "KRT10",  # cornified keratins
                   "KRT14", "KRT15",  # basal skin keratins
                   "KRT5", "KRT6A", "KRT6B", "KRT6C",  # squamous epithelia
                   "KRT16", "KRT17",  # hyperproliferative skin
                   "KRT25", "KRT26", "KRT27", "KRT28",
                   "KRT31", "KRT32", "KRT33A", "KRT33B", "KRT34",
                   "KRT35", "KRT36", "KRT37", "KRT38",
                   "KRT39", "KRT40",  # hair keratins
                   "KRT71", "KRT72", "KRT73", "KRT74", "KRT75",
                   "KRT76", "KRT77", "KRT78", "KRT79",  # type II keratins (skin)
                   "KRT80", "KRT81", "KRT82", "KRT83", "KRT84",
                   "KRT85", "KRT86",  # hair follicle keratins
                   "TCHH",  # trichohyalin
                   # KRTAP* keratin-associated handled by prefix
                   },
        "prefix_match": ["KRTAP", "LCE", "SPRR"],
        # Deliberately KEPT: KRT8, KRT18 (simple epithelium / pancreas ductal),
        # KRT19 (ductal!), KRT20 (gut/Merkel), KRT7 (ductal/glandular).
    },

    # --- Late-cornified envelope / small proline-rich ---
    # Already handled via prefix in skin_cornification

    # --- Specific nephron segments ---
    "nephron_specific": {
        "rationale": "Specific nephron segment; not in pancreas/LN. Some "
                     "broadly-expressed kidney genes (e.g., RBP1, KL) "
                     "deliberately kept.",
        "genes": {"UMOD",  # uromodulin (TAL/distal nephron)
                   "SLC12A1",  # NKCC2 (loop of Henle)
                   "SLC12A3",  # NCC (DCT)
                   "NPHS1", "NPHS2",  # podocyte slit diaphragm
                   "AQP2", "AQP3",  # collecting duct (some not strictly kidney)
                   "SLC34A1",  # phosphate cotransporter (PT)
                   "SLC22A6", "SLC22A8",  # OAT1/3 (PT)
                   "KCNJ1",  # ROMK (TAL/CD)
                   "CLCNKA", "CLCNKB",  # kidney chloride channels
                   "BSND",  # barttin (TAL)
                   "WT1",  # podocyte/Wilms — but also some other tissues
                   "ATP6V1G3",  # intercalated cell
                   "PVALB",  # also skeletal but used in kidney too — already in skel
                   "SLC12A2",  # widely expressed actually — REMOVE from list
                   "FXYD2",  # kidney (and some neural)
                   "REN",  # renin — JG cells; KEEP (vascular tissue context valid)
                   },
        # Note: REN included but KEEP-flagged below; actually I'll just remove it
        # PVALB removed (handled by skeletal); SLC12A2 removed (broad)
    },

    # --- Mammary / lactation-specific ---
    "mammary": {
        "rationale": "Mammary gland / lactation-specific; not in pancreas/LN.",
        "genes": {"CSN1S1", "CSN2", "CSN3",  # caseins
                   "LALBA",  # lactalbumin
                   "WAP",  # whey acidic protein
                   "MUC15"},  # mostly mammary/lung
    },

    # --- Sex-chromosome genes (donor-sex confound between samples) ---
    "sex_chromosome": {
        "rationale": "Y-chromosome / X-inactivation genes with strong "
                     "sex-of-donor bias. Empirical signature in this dataset: "
                     "KDM5D and DDX3Y are 10.3% / 7.7% in 0041326 vs ~0% in "
                     "0041323 — donor-sex differs between samples. Exclude "
                     "from cell-type scoring, PCA/scVI input, and any "
                     "cross-sample comparison; retain in obs as a "
                     "sex-confirmation tag only.",
        "genes": {
            # Y-chromosome single-copy genes
            "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "RPS4Y1", "RPS4Y2",
            "TBL1Y", "TMSB4Y", "USP9Y", "UTY", "ZFY", "TXLNGY", "PRKY",
            # Y-cluster repeats (TSPY family etc.)
            "TSPY1", "TSPY2", "TSPY10", "TSPYL2", "BPY2", "PRY", "VCY",
            "XKRY", "CDY1", "CDY2A", "CDY2B", "DAZ1", "DAZ2", "DAZ3", "DAZ4",
            # X-inactivation marker
            "XIST", "TSIX",
        },
    },

    # --- Salivary-restricted ---
    "salivary": {
        "rationale": "Salivary-restricted secretory proteins.",
        "genes": {"STATH",  # statherin
                   "HTN1", "HTN3",  # histatins
                   "MUC7",  # salivary mucin
                   "PRH1", "PRH2",  # proline-rich
                   "CSTA",  # cystatin A — actually broad; remove
                   "AMY2A", "AMY2B"},  # pancreatic amylase isoforms — KEEP if pancreatic
        # Removing AMY2A/B from this list — pancreas-specific.
    },
}

# Manually clean up a few edge-cases we listed and want to remove:
EXCLUDE_CATEGORIES["nephron_specific"]["genes"].discard("REN")
EXCLUDE_CATEGORIES["nephron_specific"]["genes"].discard("PVALB")
EXCLUDE_CATEGORIES["nephron_specific"]["genes"].discard("SLC12A2")
EXCLUDE_CATEGORIES["salivary"]["genes"].discard("CSTA")
EXCLUDE_CATEGORIES["salivary"]["genes"].discard("AMY2A")
EXCLUDE_CATEGORIES["salivary"]["genes"].discard("AMY2B")


# =============================================================================
# Tissue-specific overrides for the universal exclusion set
# =============================================================================
# Hassall's corpuscles in the thymic medulla legitimately express
# stratified-skin keratins (KRT1/5/10/14/15/16/17) and cornified-envelope
# genes (FLG, IVL, LOR, TGM3 + the desmoglein/desmocollin pairs).
# Excluding them in a thymus audit would discard real signal.
HASSALLS_COMPATIBLE_SKIN_GENES = {
    "KRT1", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT10",
    "KRT14", "KRT15", "KRT16", "KRT17",
    "FLG", "FLG2", "IVL", "LOR", "HRNR", "TGM3",
    "DSG1", "DSG3", "DSG4", "DSC1", "DSC3",
}

# Genes that ARE on the pancreas islet/exocrine cellular program and are
# tissue-implausible in thymus. CHGA/CHGB/SCG2/SCGN deliberately KEPT —
# pan-neuroendocrine markers, expressed by some thymic medullary cells.
# NEUROD1/ISL1/INSM1/FEV deliberately KEPT — same reason.
PANCREAS_SPECIFIC_CATEGORIES = {
    "pancreatic_exocrine": {
        "rationale": "Pancreatic exocrine zymogen / digestive enzyme / "
                     "regenerating-islet-derived protein; restricted to "
                     "pancreatic acinar / centroacinar / ductal cells. "
                     "Not expressed in thymus.",
        "genes": {
            # Serine proteases (trypsinogens, chymotrypsinogens)
            "PRSS1", "PRSS2", "PRSS3", "CTRB1", "CTRB2", "CTRC",
            # Carboxypeptidases
            "CPA1", "CPA2", "CPB1",
            # Elastases
            "CELA1", "CELA2A", "CELA3A", "CELA3B",
            # Amylases (AMY1A is salivary too; AMY2A/B are pancreas-only)
            "AMY1A", "AMY2A", "AMY2B",
            # Acinar / centroacinar marker
            "CUZD1", "SYCN", "GP2",
            # Regenerating-islet-derived (acinar / ductal)
            "REG1A", "REG1B", "REG3A", "REG3G", "REG4",
            # Lipase + colipase
            "PNLIP", "PNLIPRP1", "PNLIPRP2", "CLPS",
            # Tissue kallikrein (pancreas + salivary; on this panel salivary-context
            # KLK1 already kept; we re-flag here as exocrine-secretory)
            "KLK1",
        },
    },
    "pancreatic_endocrine": {
        "rationale": "Pancreatic islet hormone; restricted to islet "
                     "endocrine lineages (β, α, δ, γ/PP, ε). Not in thymus.",
        "genes": {
            "INS", "GCG", "SST", "PPY", "IAPP", "GHRL",
        },
    },
    "pancreatic_islet_TF": {
        "rationale": "Pancreatic-islet-specific transcription factor or "
                     "β-cell-restricted transporter. Not expressed in "
                     "thymus. Pan-neuroendocrine TFs (NEUROD1, ISL1, "
                     "INSM1, FEV, CHGA, CHGB, SCG2, SCGN) deliberately "
                     "kept because some are expressed by thymic medullary "
                     "neuroendocrine-like cells.",
        "genes": {
            "PDX1", "MAFA", "MAFB", "NKX6-1",
            "SLC30A8",  # β-cell zinc transporter
            "ARX", "PAX4",
        },
    },
}


def apply_tissue_override(categories, tissue):
    """Mutate `categories` in-place for tissue-specific deviations.

    pancreas: no-op (the dict was authored for pancreas/LN).
    thymus:   remove Hassall's-compatible skin genes from
              skin_cornification; add pancreatic_{exocrine,
              endocrine, islet_TF} categories.
    """
    if tissue == "pancreas":
        return categories
    if tissue == "thymus":
        skin = categories["skin_cornification"]
        for g in list(HASSALLS_COMPATIBLE_SKIN_GENES):
            skin["genes"].discard(g)
        skin["rationale"] = (
            "Cornified envelope / hair-keratin restricted to hard "
            "skin / hair follicle. Hassall's-corpuscle-compatible "
            "stratified keratins (KRT1/5/10/14/15/16/17) and "
            "FLG/IVL/LOR/TGM3/DSG1-3/DSC1-3 are deliberately KEPT "
            "for thymus — Hassall's corpuscles in medullary thymic "
            "epithelium genuinely express them."
        )
        for cat, info in PANCREAS_SPECIFIC_CATEGORIES.items():
            categories[cat] = {
                "rationale": info["rationale"],
                "genes": set(info["genes"]),
            }
        return categories
    raise ValueError(f"unknown tissue: {tissue!r}")


# =============================================================================
# Single-gene-panel offenders — flag for SCORING-PANEL REMOVAL specifically
# =============================================================================
# These are NOT for whole-panel exclusion; they're for cell-type scoring
# panel reform. Genes that, even when valid, should not be a 1-gene
# "panel" because argmax becomes detection-floor-sensitive.
# Only applied when tissue == 'pancreas' (GHRL is already excluded outright
# in thymus under pancreatic_endocrine).
SINGLE_GENE_PANEL_REMOVE_PANCREAS = {
    "GHRL": "Sole Epsilon marker on the panel; 1-gene scoring caused "
            "150K cells in 0041326 to be argmax-mislabeled Epsilon "
            "(see HANDOFF 2026-04-28). Either drop Epsilon from the "
            "celltype dict entirely (recommended) or expand the panel "
            "with co-markers (ARX is α-also; no clean Epsilon multi-gene "
            "panel exists in the 5K).",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tissue",
        choices=("pancreas", "thymus"),
        default="pancreas",
        help="Target tissue. Controls which exclusion overrides apply.",
    )
    ap.add_argument(
        "--samples",
        nargs="+",
        default=["0041323", "0041326"],
        help="Sample IDs. Used to read detection CSVs and to name "
             "output columns.",
    )
    ap.add_argument(
        "--output",
        default=str(ROOT / "panel_audit.csv"),
        help="Output CSV path.",
    )
    ap.add_argument(
        "--staging-dir",
        default="/tmp",
        help="Directory containing panel_genes_union.csv and "
             "panel_detection_{sample}.csv (default: /tmp).",
    )
    args = ap.parse_args()

    if args.tissue != "pancreas" and len(args.samples) != 2:
        # column naming for ratio assumes exactly 2 samples; keep that contract
        if len(args.samples) < 2:
            raise SystemExit("--samples must list at least 2 sample IDs")

    staging = Path(args.staging_dir)
    categories = apply_tissue_override(EXCLUDE_CATEGORIES, args.tissue)
    single_gene_panel_remove = (
        SINGLE_GENE_PANEL_REMOVE_PANCREAS if args.tissue == "pancreas" else {}
    )

    union = pd.read_csv(staging / "panel_genes_union.csv")["gene"].tolist()
    union_set = set(union)

    # Resolve prefix-matched categories
    for cat, info in categories.items():
        if "prefix_match" in info:
            for prefix in info["prefix_match"]:
                for g in union:
                    if g.startswith(prefix):
                        info["genes"].add(g)

    # Build gene → category mapping
    gene_to_cat = {}
    gene_to_rationale = {}
    for cat, info in categories.items():
        for g in info["genes"]:
            if g in union_set:
                # If a gene already tagged, keep first match (categories
                # ordered with most-confident first); avoids overlap.
                gene_to_cat.setdefault(g, cat)
                gene_to_rationale.setdefault(g, info["rationale"])

    # Per-sample detection
    det = {}
    for s in args.samples:
        d = pd.read_csv(staging / f"panel_detection_{s}.csv", index_col=0)
        det[s] = d

    # The first two samples define the cross-sample log2 ratio direction.
    s1, s2 = args.samples[0], args.samples[1]
    ratio_col = f"log2_detection_ratio_{s2}_over_{s1}"

    # Build the audit table
    rows = []
    for g in sorted(union_set):
        cat = gene_to_cat.get(g, "")
        rationale = gene_to_rationale.get(g, "")
        single = single_gene_panel_remove.get(g, "")

        d_per_sample = {
            s: float(det[s]["detection_pct"].get(g, np.nan))
            for s in args.samples
        }
        n_per_sample = {
            s: int(det[s]["n_cells"].get(g, 0))
            for s in args.samples
        }

        # Sample-asymmetry: ratio of detection rates for the first two samples
        d1 = d_per_sample[s1]
        d2 = d_per_sample[s2]
        if d1 > 0.001 and d2 > 0.001:
            asym_log2 = float(np.log2(d2 / d1))
        else:
            asym_log2 = np.nan

        row = {
            "gene": g,
            "category": cat,
            "rationale": rationale,
            "single_gene_panel_note": single,
            "exclude_recommended": "yes" if cat else "no",
        }
        for s in args.samples:
            row[f"detection_pct_{s}"] = round(d_per_sample[s], 3)
        for s in args.samples:
            row[f"n_cells_{s}"] = n_per_sample[s]
        row[ratio_col] = (round(asym_log2, 3)
                          if np.isfinite(asym_log2) else "")
        rows.append(row)

    df = pd.DataFrame(rows)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    # Summary
    n_total = len(df)
    n_excluded = (df["exclude_recommended"] == "yes").sum()
    n_per_cat = (df[df["category"] != ""]
                 .groupby("category").size()
                 .sort_values(ascending=False))
    print(f"Panel audit ({args.tissue}): {n_total:,} unique genes across "
          f"{len(args.samples)} samples")
    print(f"  recommended exclusions: {n_excluded:,} "
          f"({100*n_excluded/n_total:.1f}%)")
    print("  by category:")
    print(n_per_cat.to_string())
    print()

    # Sample-asymmetry top-15
    asym_df = df[df[ratio_col] != ""].copy()
    asym_df[ratio_col] = pd.to_numeric(asym_df[ratio_col], errors="coerce")
    asym_df = asym_df.dropna(subset=[ratio_col])
    asym_df["abs_log2"] = asym_df[ratio_col].abs()
    print("Top 15 sample-asymmetric genes (|log2 ratio| largest):")
    cols = ["gene", "category",
            f"detection_pct_{s1}", f"detection_pct_{s2}",
            ratio_col]
    print(asym_df.nlargest(15, "abs_log2")[cols].to_string(index=False))
    print()
    if args.tissue == "pancreas":
        ghrl = df[df["gene"] == "GHRL"]
        if not ghrl.empty:
            print("GHRL specifically (the Epsilon mislabel root cause):")
            print(ghrl[["gene", f"detection_pct_{s1}",
                        f"detection_pct_{s2}",
                        ratio_col,
                        "single_gene_panel_note"]].to_string(index=False))
    print(f"\nFull table: {out}")


if __name__ == "__main__":
    main()
