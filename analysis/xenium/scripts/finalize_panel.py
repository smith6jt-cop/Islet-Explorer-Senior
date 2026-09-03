"""Stage 0 — finalize panel masks + donor_sex tag on each phenotyped/raw h5ad.

CLI:
    --paths PATH [PATH ...]   explicit h5ad paths (overrides --samples)
    --samples ID1 ID2 ...     defaults to data/processed/{ID}/{ID}_phenotyped.h5ad
    --audit-csv PATH          panel_audit.csv to apply (default pancreas)
    --dry-run                 do not save

Reads:
    --audit-csv (default data/processed/panel_audit.csv) — candidate exclusions
    each h5ad listed (default _phenotyped.h5ad for each sample)

Writes back, in-place:
    adata.var['panel_for_scoring']    bool — True for genes safe in cell-type
                                       scoring (drops non-sex tissue-restricted)
    adata.var['panel_for_embedding']  bool — additionally drops sex-chromosome
                                       genes; for HVG/PCA/scVI input
    adata.var['exclude_category']     str — audit category (or '' for kept genes)
    adata.var['subpanel_<lineage>']   bool — per-broad-lineage subpanel membership
                                       (from data/reference/subpanels; one per
                                       tissue lineage). Skipped with --no-subpanels.
    adata.var['celltype_panel']       bool — union of subpanel genes ∩ panel_for_
                                       embedding; the Tier-1 scVI gene set used by
                                       scripts/hierarchical_phenotype.py.
    adata.obs['donor_sex']            str — 'male' / 'female' from KDM5D+DDX3Y

This script is idempotent: re-running it just overwrites the flags. It writes
to whichever h5ad file is passed; the same logic applies to raw zarr ingestion
output and to phenotyped outputs. Tissue-agnostic — pass the appropriate
--audit-csv for the tissue.
"""
import argparse
import importlib
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.workflow_paths import require_existing_files
DEFAULT_SAMPLES = ("0041323", "0041326")
DEFAULT_AUDIT_CSV = ROOT / "panel_audit.csv"
DEFAULT_SUBPANEL_DIR = REPO / "data/development/xenium_panels/subpanels"

# Sex-call thresholds (lognorm scale)
# KDM5D + DDX3Y are reliably 7-10% in male, near-zero in female.
SEX_GENE_THRESH_LOGNORM = np.log1p(1.0)
SEX_DETECTION_FRAC = 0.02  # ≥ 2% of cells expressing → male


def annotate_var_flags(adata, audit_df):
    """Add panel_for_scoring, panel_for_embedding, exclude_category to var."""
    audit_lookup = audit_df.set_index("gene")["category"].to_dict()
    cat_col = adata.var_names.map(lambda g: audit_lookup.get(g, ""))
    adata.var["exclude_category"] = pd.Categorical(cat_col)
    sex_chrom_genes = set(audit_df[audit_df["category"] == "sex_chromosome"]["gene"])
    nonsex_excluded = set(audit_df[audit_df["exclude_recommended"] == "yes"]["gene"]) - sex_chrom_genes

    adata.var["panel_for_scoring"] = ~adata.var_names.isin(nonsex_excluded)
    adata.var["panel_for_embedding"] = (~adata.var_names.isin(nonsex_excluded)
                                          & ~adata.var_names.isin(sex_chrom_genes))


def annotate_subpanel_masks(adata, broad_to_genes, tissue=None, mod=None):
    """Write var['subpanel_<L>'] per broad lineage + the identity/state masks.

    `broad_to_genes` is {lineage: set(genes)} — the upstream ontology-curated
    IDENTITY_CORE per lineage (scripts/tissue_markers/_identity_core.broad_genes);
    the special key '_embedding_extra' contributes to the union ONLY (no per-lineage
    mask).

    Masks stamped (IDENTITY-CORE rewire 2026-06-08 session 2):
      subpanel_<L>      per-lineage IDENTITY mask = the clean identity_core ∩ panel
                        ∩ embedding. NO _identity_trim — the upstream cores are
                        already CL-curated + specificity-audited cross-lineage-clean.
                        Used by the broad CALL (absolute scoring + raw-count anchors)
                        AND Tier-2 within-lineage sub-clustering.
      celltype_panel    union ∩ panel_for_embedding (bookkeeping; == identity_panel
                        minus the promoted genes).
      identity_panel    the scANVI FEATURE set = (union ∪ promoted) ∩ embedding. Now
                        the TIGHT identity_core union (not the fat subpanels) — clean
                        features + clean seeds (lesson 4) are the catch-all fix.
      is_state_only     genes whose ONLY vendored membership is a state/fine_tier
                        panel (audit of the GUARD; not in any identity panel).
      identity_promoted genes promoted from a state/fine_tier panel into an identity
                        mask (thymus mTEC AIRE/FEZF2; upstream Rec 3).

    When tissue/mod are None (e.g. --no-subpanels off-path), only subpanel_<L> +
    celltype_panel are written (legacy behavior).
    """
    if "panel_for_embedding" not in adata.var.columns:
        raise ValueError("annotate_subpanel_masks requires panel_for_embedding "
                         "— call annotate_var_flags first")
    live = set(adata.var_names)
    embed = adata.var["panel_for_embedding"].to_numpy()

    # --- per-lineage IDENTITY-CORE masks (no trim — cores are pre-clean) ----------
    union = set()
    for lineage, genes in sorted(broad_to_genes.items()):
        present = genes & live
        union |= present
        if lineage.startswith("_"):
            print(f"    {lineage}: {len(present)}/{len(genes)} on panel (union-only)")
            continue
        adata.var[f"subpanel_{lineage}"] = adata.var_names.isin(present)
        print(f"    subpanel_{lineage} (identity_core): {len(present)}/{len(genes)} "
              f"on panel")
    adata.var["celltype_panel"] = adata.var_names.isin(union) & embed
    print(f"    celltype_panel (identity_core union ∩ embedding): "
          f"{int(adata.var['celltype_panel'].sum()):,}/{adata.n_vars:,}")

    if tissue is None:
        return

    # --- identity_panel (scANVI feature set) + GUARD audit masks -----------------
    from scripts.tissue_markers import _panel_roles
    promoted = _panel_roles.promoted_identity_genes(tissue)          # {gene: lineage}
    prom_present = {g for g in promoted if g in live}
    identity_genes = union | prom_present
    adata.var["identity_panel"] = adata.var_names.isin(identity_genes) & embed
    adata.var["identity_promoted"] = adata.var_names.isin(prom_present)
    # promote AIRE/FEZF2-style genes into their target subpanel mask too
    for g, lineage in promoted.items():
        col = f"subpanel_{lineage}"
        if g in live and col in adata.var.columns:
            adata.var.loc[g, col] = True

    state_side = _panel_roles.panel_genes(tissue, {"module_score", "fine_tier"})
    is_state_only = (adata.var_names.isin(state_side)
                     & ~adata.var_names.isin(identity_genes))
    adata.var["is_state_only"] = is_state_only
    n_id = int(adata.var["identity_panel"].sum())
    n_emb = int(embed.sum())
    print(f"    identity_panel (scANVI features): {n_id:,}/{n_emb:,} of "
          f"panel_for_embedding ({n_emb - n_id:,} pure-state/pathway genes excluded)")
    if prom_present:
        print(f"    identity_promoted: {sorted(prom_present)} -> "
              f"{ {g: promoted[g] for g in sorted(prom_present)} }")
    print(f"    is_state_only (audit): {int(is_state_only.sum()):,} genes")


def call_donor_sex(adata):
    """Return 'male' or 'female' from KDM5D + DDX3Y detection."""
    layer = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    sex_genes = [g for g in ("KDM5D", "DDX3Y") if g in adata.var_names]
    if not sex_genes:
        return "unknown"
    n_cells = adata.n_obs
    n_detected_max = 0
    for g in sex_genes:
        i = adata.var_names.get_loc(g)
        col = layer[:, i]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        col = np.asarray(col).ravel()
        n_det = int((col >= SEX_GENE_THRESH_LOGNORM).sum())
        n_detected_max = max(n_detected_max, n_det)
    frac = n_detected_max / max(n_cells, 1)
    return "male" if frac >= SEX_DETECTION_FRAC else "female"


def process_one(p, audit_df, subpanel_genes=None, tissue=None, mod=None,
                dry_run=False):
    print(f"\n=== {p.name} ===")
    a = sc.read_h5ad(p)
    print(f"    shape: {a.shape}")
    annotate_var_flags(a, audit_df)
    n_score = int(a.var["panel_for_scoring"].sum())
    n_embed = int(a.var["panel_for_embedding"].sum())
    print(f"    panel_for_scoring: {n_score:,} / {a.n_vars:,}")
    print(f"    panel_for_embedding: {n_embed:,} / {a.n_vars:,}")
    if subpanel_genes is not None:
        annotate_subpanel_masks(a, subpanel_genes, tissue=tissue, mod=mod)
    sex = call_donor_sex(a)
    a.obs["donor_sex"] = pd.Categorical([sex] * a.n_obs,
                                          categories=["female", "male", "unknown"])
    print(f"    donor_sex: {sex}")
    if dry_run:
        print("    DRY RUN — not saving")
        return
    a.write(p, compression="gzip")
    print("    saved")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", nargs="+", default=None,
                       help="h5ad paths to annotate (default: --samples' "
                            "_phenotyped.h5ad files under data/processed/)")
    ap.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES),
                       help=f"Sample IDs to default-resolve to "
                            f"data/processed/{{s}}/{{s}}_phenotyped.h5ad "
                            f"when --paths is not given (default: "
                            f"{' '.join(DEFAULT_SAMPLES)})")
    ap.add_argument("--audit-csv", default=str(DEFAULT_AUDIT_CSV),
                       help=f"Panel audit CSV produced by "
                            f"audit_panel_exclusions.py (default: "
                            f"{DEFAULT_AUDIT_CSV})")
    ap.add_argument("--tissue", choices=("pancreas", "thymus"),
                       default="pancreas",
                       help="tissue module supplying the subpanel allowlist "
                            "(default: pancreas)")
    ap.add_argument("--subpanel-dir", default=str(DEFAULT_SUBPANEL_DIR),
                       help=f"vendored cell-type subpanel root "
                            f"(default: {DEFAULT_SUBPANEL_DIR})")
    ap.add_argument("--no-subpanels", action="store_true",
                       help="skip writing the subpanel_* / celltype_panel masks")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    if args.paths:
        candidates = [Path(p) for p in args.paths]
    else:
        candidates = [DATA_ROOT / f"processed/{s}/{s}_phenotyped.h5ad"
                      for s in args.samples]
    try:
        paths = require_existing_files(
            candidates, context="Xenium panel finalization"
        )
    except FileNotFoundError as exc:
        ap.error(str(exc))

    audit_df = pd.read_csv(args.audit_csv)
    print(f"loaded audit: {len(audit_df):,} genes, "
          f"{(audit_df['exclude_recommended']=='yes').sum()} flagged "
          f"({args.audit_csv})")

    # Subpanel gene sets are sample-independent (intersection with each sample's
    # live var_names happens inside annotate_subpanel_masks). IDENTITY-CORE rewire
    # (2026-06-08 session 2): masks come from the upstream ontology-curated
    # identity_core (clean), NOT the fat subpanels.
    subpanel_genes = None
    mod = None
    if args.no_subpanels:
        print("[subpanels] --no-subpanels set; skipping subpanel masks")
    else:
        try:
            mod = importlib.import_module(f"scripts.tissue_markers.{args.tissue}")
            if hasattr(mod, "SUBPANEL_ALLOWLIST"):
                # Fail-loud if the vendored panel_roles.csv disagrees with the
                # hand-maintained allowlist (drift guard) before touching any sample.
                from scripts.tissue_markers import _identity_core, _panel_roles
                _panel_roles.assert_consistent(mod, args.tissue)
                broad = _identity_core.broad_genes(args.tissue, mod, args.subpanel_dir)
                extra = _identity_core.embedding_extra_genes(args.tissue, mod,
                                                             args.subpanel_dir)
                subpanel_genes = {L: set(g) for L, g in broad.items()}
                subpanel_genes["_embedding_extra"] = set(extra)
                tot = len(set().union(*subpanel_genes.values())) if subpanel_genes else 0
                print(f"[identity_core] loaded {len(broad)} broad lineages "
                      f"({tot:,} core genes pre-intersection) from "
                      f"{args.subpanel_dir}/identity_core")
            else:
                print(f"[subpanels] tissue '{args.tissue}' has no SUBPANEL_ALLOWLIST "
                      f"— skipping subpanel masks")
        except FileNotFoundError as e:
            print(f"[identity_core] vendored snapshot incomplete ({e}); "
                  f"skipping subpanel masks")

    for p in paths:
        process_one(p, audit_df, subpanel_genes=subpanel_genes,
                    tissue=args.tissue, mod=mod, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
