#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan that is also the
NAMED PRODUCER of a live artifact under ``data/app_data`` (code tier 2), which is why it is
kept rather than deleted. It cannot be re-run: the inputs named below are absent from this
checkout.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method, and ``data/app_data``
is the OTHER app domain: the two datasets share five donors and no vocabulary.

Build enriched H5AD for Islet Explorer app.

Takes the validated trajectory h5ad (data/adata_ins_root.h5ad) and enriches it with:
1. QuPath-derived data from groovy exports (markers, targets, composition)
2. Donor metadata from CODEX_Pancreas_Donors.xlsx

Output: data/islet_explorer.h5ad — single file containing all app data.

Canonical lineage:
    CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, single-cell)
      ↓ fixed_islet_aggregation.py
    islets_core_fixed.h5ad  (islet-level, core only)
      ↓ rebuild_trajectory.ipynb
    data/adata_ins_root.h5ad  (+ pseudotime + UMAP)
      ↓ THIS SCRIPT
    data/islet_explorer.h5ad  (+ groovy data + donor metadata)

Usage:
    python scripts/build_h5ad_for_app.py
    python scripts/build_h5ad_for_app.py --trajectory data/adata_ins_root.h5ad --output data/islet_explorer.h5ad
"""

import argparse
import os
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


# ── Constants ──────────────────────────────────────────────────────────────

GROOVY_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'groovy_exports')
DONOR_KEY_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'CODEX_Pancreas_Donors.xlsx')
DEFAULT_TRAJECTORY = os.path.join(os.path.dirname(__file__), '..', 'data', 'adata_ins_root.h5ad')
DEFAULT_OUTPUT = os.path.join(os.path.dirname(__file__), '..', 'data', 'islet_explorer.h5ad')


# ── Donor metadata loading ────────────────────────────────────────────────

def parse_aab_flags(aabs_value):
    """Parse autoantibody string into boolean flags."""
    flags = {
        'AAb_GADA': False, 'AAb_IA2A': False, 'AAb_ZnT8A': False,
        'AAb_IAA': False, 'AAb_mIAA': False,
    }
    if not isinstance(aabs_value, str) or not aabs_value.strip():
        return flags
    tokens = {re.sub(r'[^A-Za-z0-9]', '', p).upper()
              for p in aabs_value.split(',') if p.strip()}
    if 'GADA' in tokens:
        flags['AAb_GADA'] = True
    if 'IA2A' in tokens:
        flags['AAb_IA2A'] = True
    if 'ZNT8A' in tokens or 'ZNT8' in tokens:
        flags['AAb_ZnT8A'] = True
    if 'MIAA' in tokens:
        flags['AAb_mIAA'] = True
    if 'IAA' in tokens:
        flags['AAb_IAA'] = True
    return flags


def load_donor_key(key_path):
    """Load donor metadata from CODEX_Pancreas_Donors.xlsx."""
    key = pd.read_excel(key_path, sheet_name=0)
    # Drop trailing empty columns (Unnamed: N)
    key = key.loc[:, ~key.columns.str.startswith('Unnamed')]
    key['Case ID'] = key['Case ID'].astype(str).str.zfill(4)
    if 'Aabs' in key.columns:
        flags_df = key['Aabs'].apply(parse_aab_flags).apply(pd.Series)
        for c in flags_df.columns:
            key[c] = flags_df[c].astype(bool)
    return key


# ── Groovy TSV loading ────────────────────────────────────────────────────

def parse_case_id(filename):
    """Extract case ID from filename like islet_targets_6479.tsv."""
    m = re.search(r'_(\d+)\.tsv$', filename)
    if m:
        return m.group(1).zfill(4)
    return None


def load_groovy_targets(groovy_dir, donor_key):
    """Load all islet_targets_*.tsv files into a single DataFrame."""
    frames = []
    for f in sorted(Path(groovy_dir).glob('islet_targets_*.tsv')):
        df = pd.read_csv(f, sep='\t')
        case_id = parse_case_id(f.name)
        df.insert(0, 'Case ID', case_id)
        df.insert(1, 'source_file', f.name)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.merge(donor_key, on='Case ID', how='left')
    return combined


def load_groovy_markers(groovy_dir, donor_key):
    """Load all islet_markers_*.tsv and LGALS3_*.tsv files."""
    frames = []
    for pattern in ['islet_markers_*.tsv', 'LGALS3_*.tsv']:
        for f in sorted(Path(groovy_dir).glob(pattern)):
            df = pd.read_csv(f, sep='\t')
            case_id = parse_case_id(f.name)
            df.insert(0, 'Case ID', case_id)
            df.insert(1, 'source_file', f.name)
            frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.merge(donor_key, on='Case ID', how='left')
    return combined


def load_groovy_composition(groovy_dir, donor_key):
    """Load all islet_composition_*.tsv files."""
    frames = []
    for f in sorted(Path(groovy_dir).glob('islet_composition_*.tsv')):
        df = pd.read_csv(f, sep='\t')
        case_id = parse_case_id(f.name)
        df.insert(0, 'Case ID', case_id)
        df.insert(1, 'source_file', f.name)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.merge(donor_key, on='Case ID', how='left')
    return combined


# ── Build islet key mapping ───────────────────────────────────────────────

def extract_islet_key_from_region(region_str):
    """Extract islet key from region string like 'Islet_156_core' or 'Islet_156_band'.

    Returns (islet_name, region_type) e.g. ('Islet_156', 'islet_core').
    """
    if not isinstance(region_str, str):
        return None, None

    s = region_str.strip()
    # Match patterns like Islet_156_core, Islet_156_band
    m = re.match(r'^(Islet_\d+)_(core|band)$', s)
    if m:
        islet_name = m.group(1)
        region_type = f'islet_{m.group(2)}'
        return islet_name, region_type
    return None, None


def build_islet_key(case_id, imageid, islet_name):
    """Build a unique islet key compatible with the h5ad islet_id format.

    The h5ad uses: imageid_Islet_N (e.g., '6479_Islet_156')
    Groovy data uses: Case ID + region name
    """
    return f"{imageid}_{islet_name}"


# ── Main build function ──────────────────────────────────────────────────

def build_enriched_h5ad(trajectory_path, groovy_dir, donor_key_path, output_path):
    """Build enriched H5AD by merging trajectory data with groovy exports."""

    print("=" * 60)
    print("Building enriched H5AD for Islet Explorer")
    print("=" * 60)

    # Load trajectory h5ad
    print(f"\n1. Loading trajectory: {trajectory_path}")
    adata = ad.read_h5ad(trajectory_path)
    print(f"   Shape: {adata.shape}")
    print(f"   Obs columns: {list(adata.obs.columns)}")
    print(f"   Obsm keys: {list(adata.obsm.keys())}")

    # Build donor key from the h5ad's own obs data (imageid -> donor_status, age, gender, AAb flags)
    # The external CODEX_Pancreas_Donors.xlsx has a different donor cohort and cannot be used here.
    print(f"\n2. Building donor key from h5ad obs")
    adata.obs['case_id_str'] = adata.obs['imageid'].astype(str).str.extract(r'(\d+)')[0].str.zfill(4)

    # Build per-donor metadata from h5ad obs (one row per donor)
    donor_key_cols = ['case_id_str', 'donor_status', 'age', 'gender']
    aab_cols = ['GADA', 'ZnT8A', 'IA2A', 'mIAA']
    donor_key_cols += [c for c in aab_cols if c in adata.obs.columns]
    donor_key_df = adata.obs[donor_key_cols].drop_duplicates(subset=['case_id_str']).copy()

    # Rename to match Excel conventions
    donor_key_df = donor_key_df.rename(columns={
        'case_id_str': 'Case ID',
        'donor_status': 'Donor Status',
        'age': 'Age',
        'gender': 'Gender',
    })
    # Create AAb flag columns matching Excel format
    for aab in aab_cols:
        if aab in adata.obs.columns:
            donor_key_df[f'AAb_{aab}'] = donor_key_df[aab].apply(
                lambda x: bool(x) if not pd.isna(x) else False
            )
            donor_key_df.drop(columns=[aab], inplace=True)

    print(f"   Donors: {len(donor_key_df)}")
    print(f"   Donor Status: {dict(donor_key_df['Donor Status'].value_counts())}")

    # Load groovy data with h5ad-derived donor key
    print(f"\n3. Loading groovy exports: {groovy_dir}")

    targets_df = load_groovy_targets(groovy_dir, donor_key_df)
    print(f"   Targets: {len(targets_df)} rows")

    markers_df = load_groovy_markers(groovy_dir, donor_key_df)
    print(f"   Markers: {len(markers_df)} rows")

    comp_df = load_groovy_composition(groovy_dir, donor_key_df)
    print(f"   Composition: {len(comp_df)} rows")

    # Build mapping from h5ad islet_id to groovy data
    # h5ad islet_id format: "imageid_Islet_N" (e.g., "6479_Islet_156")
    # Groovy region format: "Islet_N_core" or "Islet_N_band"

    # We need to map Case ID (from donor key) to imageid (in h5ad)
    # Build case_id -> imageid lookup from adata
    case_to_imageid = dict(zip(
        adata.obs['case_id_str'].values,
        adata.obs['imageid'].astype(str).values
    ))
    print(f"\n   Case ID -> imageid mapping: {case_to_imageid}")

    # Store groovy DataFrames in adata.uns for the Shiny app to extract
    # This allows the R app to reconstruct the same data frames as prep_data()
    print(f"\n4. Storing groovy data in adata.uns")

    # Convert to serializable format for h5ad storage
    # h5ad supports storing DataFrames in .uns as dict-of-arrays

    def df_to_uns_dict(df, name):
        """Convert DataFrame to a dict suitable for h5ad .uns storage."""
        result = {}
        for col in df.columns:
            vals = df[col].values
            # Convert to appropriate types
            if pd.api.types.is_bool_dtype(vals):
                result[col] = vals.astype(bool)
            elif pd.api.types.is_numeric_dtype(vals):
                result[col] = vals.astype(float)
            else:
                result[col] = vals.astype(str)
        return result

    # Store targets (core region only for density, all regions for counts)
    if len(targets_df) > 0:
        # Add islet_key for matching
        targets_df['_islet_name'] = targets_df['region'].apply(
            lambda x: extract_islet_key_from_region(x)[0] if isinstance(x, str) else None
        )
        targets_df['_region_type'] = targets_df['region'].apply(
            lambda x: extract_islet_key_from_region(x)[1] if isinstance(x, str) else None
        )
        adata.uns['groovy_targets_columns'] = list(targets_df.columns)
        adata.uns['groovy_targets_n_rows'] = len(targets_df)
        # Store as compressed arrays
        for col in targets_df.columns:
            vals = targets_df[col].values
            key = f'groovy_targets_{col}'
            if pd.api.types.is_bool_dtype(vals):
                adata.uns[key] = vals.astype(bool)
            elif pd.api.types.is_numeric_dtype(vals):
                adata.uns[key] = np.array(vals, dtype=float)
            else:
                adata.uns[key] = np.array(vals, dtype=str)

    # Store markers
    if len(markers_df) > 0:
        markers_df['_islet_name'] = markers_df['region'].apply(
            lambda x: extract_islet_key_from_region(x)[0] if isinstance(x, str) else None
        )
        markers_df['_region_type'] = markers_df['region'].apply(
            lambda x: extract_islet_key_from_region(x)[1] if isinstance(x, str) else None
        )
        adata.uns['groovy_markers_columns'] = list(markers_df.columns)
        adata.uns['groovy_markers_n_rows'] = len(markers_df)
        for col in markers_df.columns:
            vals = markers_df[col].values
            key = f'groovy_markers_{col}'
            if pd.api.types.is_bool_dtype(vals):
                adata.uns[key] = vals.astype(bool)
            elif pd.api.types.is_numeric_dtype(vals):
                adata.uns[key] = np.array(vals, dtype=float)
            else:
                adata.uns[key] = np.array(vals, dtype=str)

    # Store composition
    if len(comp_df) > 0:
        adata.uns['groovy_composition_columns'] = list(comp_df.columns)
        adata.uns['groovy_composition_n_rows'] = len(comp_df)
        for col in comp_df.columns:
            vals = comp_df[col].values
            key = f'groovy_composition_{col}'
            if pd.api.types.is_bool_dtype(vals):
                adata.uns[key] = vals.astype(bool)
            elif pd.api.types.is_numeric_dtype(vals):
                adata.uns[key] = np.array(vals, dtype=float)
            else:
                adata.uns[key] = np.array(vals, dtype=str)

    # Step 4.5: Merge neighborhood metrics into adata.obs
    nbr_path = os.path.join(os.path.dirname(output_path), 'neighborhood_metrics.csv')
    if os.path.exists(nbr_path):
        print(f"\n4.5. Merging neighborhood metrics: {nbr_path}")
        nbr = pd.read_csv(nbr_path)
        print(f"     Loaded {len(nbr)} rows, {len(nbr.columns)} columns")

        # Build combined_islet_id in adata.obs for merge
        adata.obs['_combined_islet_id'] = (
            adata.obs['imageid'].astype(str) + '_' + adata.obs['base_islet_id'].astype(str)
        )

        # Merge neighborhood columns (exclude combined_islet_id itself)
        nbr_cols = [c for c in nbr.columns if c != 'combined_islet_id']
        nbr_merge = nbr.set_index('combined_islet_id')[nbr_cols]

        before_cols = len(adata.obs.columns)
        adata.obs = adata.obs.join(nbr_merge, on='_combined_islet_id', how='left')
        adata.obs.drop(columns=['_combined_islet_id'], inplace=True)
        after_cols = len(adata.obs.columns)
        print(f"     Added {after_cols - before_cols} neighborhood columns to .obs")

        n_with_data = adata.obs['total_cells_peri'].notna().sum() if 'total_cells_peri' in adata.obs.columns else 0
        print(f"     Islets with peri-islet data: {n_with_data}/{len(adata.obs)}")
    else:
        print(f"\n4.5. Neighborhood metrics not found at {nbr_path} (skipping)")

    # Step 4.7: Merge Leiden clustering from islets_core_clustered.h5ad
    clust_path = os.path.join(os.path.dirname(output_path), '..', 'islet_analysis', 'islets_core_clustered.h5ad')
    if os.path.exists(clust_path):
        print(f"\n4.7. Merging Leiden clustering: {clust_path}")
        clust = ad.read_h5ad(clust_path)
        leiden_cols = [c for c in clust.obs.columns if c.startswith("leiden_")]
        print(f"     Leiden columns: {leiden_cols}")

        # Merge Leiden cluster assignments by index alignment (both indexed by islet_id)
        for col in leiden_cols:
            adata.obs[col] = clust.obs[col].astype(str).reindex(adata.obs.index)

        # Store UMAP coords from .obsm['X_umap'] as leiden_umap_1/2 in .obs
        if 'X_umap' in clust.obsm:
            umap_df = pd.DataFrame(
                clust.obsm['X_umap'],
                index=clust.obs.index,
                columns=['leiden_umap_1', 'leiden_umap_2']
            )
            for col in umap_df.columns:
                adata.obs[col] = umap_df[col].reindex(adata.obs.index)
            print(f"     Merged {len(leiden_cols)} Leiden columns + 2 UMAP coords")
        else:
            # Fallback: check if umap_1/umap_2 are in .obs
            if all(c in clust.obs.columns for c in ['umap_1', 'umap_2']):
                adata.obs['leiden_umap_1'] = clust.obs['umap_1'].reindex(adata.obs.index)
                adata.obs['leiden_umap_2'] = clust.obs['umap_2'].reindex(adata.obs.index)
                print(f"     Merged {len(leiden_cols)} Leiden columns + 2 UMAP coords (from .obs)")
            else:
                print(f"     Merged {len(leiden_cols)} Leiden columns (no UMAP coords found)")
    else:
        print(f"\n4.7. Clustered H5AD not found at {clust_path} (skipping Leiden merge)")

    # Store provenance metadata
    adata.uns['data_provenance'] = {
        'pipeline': 'Phase 2: Data Audit, QC Validation & Pipeline Formalization',
        'trajectory_source': os.path.basename(trajectory_path),
        'groovy_source': os.path.basename(groovy_dir),
        'donor_key_source': os.path.basename(donor_key_path),
        'n_target_rows': len(targets_df),
        'n_marker_rows': len(markers_df),
        'n_composition_rows': len(comp_df),
    }

    # Store case_to_imageid mapping for the R app
    adata.uns['case_to_imageid'] = case_to_imageid

    # Ensure string/categorical columns have no NaN (h5py can't write mixed str/NaN)
    for col in adata.obs.columns:
        if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            adata.obs[col] = adata.obs[col].astype(str).fillna('')
        elif adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].fillna('').astype(str)

    # Write output
    print(f"\n5. Writing output: {output_path}")
    adata.write_h5ad(output_path)

    file_size = os.path.getsize(output_path) / (1024 * 1024)
    print(f"   Size: {file_size:.1f} MB")
    print(f"   Shape: {adata.shape}")
    print(f"   Obs columns: {list(adata.obs.columns)}")
    print(f"   Obsm keys: {list(adata.obsm.keys())}")
    print(f"   Uns keys: {[k for k in adata.uns.keys() if not k.startswith('groovy_')]}")

    print(f"\n{'='*60}")
    print(f"SUCCESS: {output_path}")
    print(f"{'='*60}")

    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build enriched H5AD for Islet Explorer app')
    parser.add_argument('--trajectory', default=DEFAULT_TRAJECTORY,
                        help='Path to validated trajectory h5ad')
    parser.add_argument('--groovy-dir', default=GROOVY_DIR,
                        help='Path to groovy_exports directory')
    parser.add_argument('--donor-key', default=DONOR_KEY_PATH,
                        help='Path to CODEX_Pancreas_Donors.xlsx')
    parser.add_argument('--output', default=DEFAULT_OUTPUT,
                        help='Output path for enriched h5ad')
    args = parser.parse_args()

    # Resolve paths relative to script location
    script_dir = os.path.dirname(os.path.abspath(__file__))

    trajectory = args.trajectory
    if not os.path.isabs(trajectory):
        trajectory = os.path.join(script_dir, trajectory)

    groovy_dir = args.groovy_dir
    if not os.path.isabs(groovy_dir):
        groovy_dir = os.path.join(script_dir, groovy_dir)

    donor_key = args.donor_key
    if not os.path.isabs(donor_key):
        donor_key = os.path.join(script_dir, donor_key)

    output = args.output
    if not os.path.isabs(output):
        output = os.path.join(script_dir, output)

    # Verify inputs exist
    for path, label in [(trajectory, 'trajectory h5ad'), (groovy_dir, 'groovy exports'), (donor_key, 'donor key')]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    build_enriched_h5ad(trajectory, groovy_dir, donor_key, output)
