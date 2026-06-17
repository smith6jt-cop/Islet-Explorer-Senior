"""
Fixed islet aggregation functions for CODEX pancreas analysis.

Key fixes:
1. Excludes peri-islet expanded regions (_exp20um) - only uses true islet cells
2. Keeps only protein expression in .X matrix
3. Stores phenotype counts and morphology separately in .obsm
4. Stores scVI latent means in .obsm for consistent use in downstream analysis
"""

import pandas as pd
import numpy as np
from scipy import sparse
import anndata as ad


def create_islet_dataset_fixed(adata_subset, region='islet_only', min_cells=0,
                                require_paired=True, _skip_diagnostic=False,
                                _valid_islet_ids=None):
    """
    Create islet-level dataset from single-cell data.

    Parameters
    ----------
    adata_subset : AnnData
        Single-cell AnnData with phenotyped cells.
    region : str, default='islet_only'
        Which cells to include:
        - 'islet_only': Only cells within islet boundaries (recommended for islet biology)
        - 'peri_only': Only cells in the 20μm expanded peri-islet region (microenvironment)
        - 'combined': Both islet and peri-islet cells combined (original behavior)
        - 'separate': Returns a tuple of (islet_adata, peri_adata) for comparative analysis
    min_cells : int, default=0
        Minimum number of cells required per islet. Islets with fewer cells are dropped.
        With default=0, all islets with at least 1 cell are kept.
    require_paired : bool, default=True
        If True, only keep islets that have BOTH core and peri-islet regions.

    Returns
    -------
    AnnData or tuple of AnnData
        Islet-level dataset(s) with:
        - .X: Mean protein expression only (n_islets × n_proteins)
        - .obsm['X_scVI_mean']: Mean scVI latent representation
        - .obsm['phenotype_counts']: Cell type counts per islet
        - .obsm['phenotype_proportions']: Cell type proportions per islet
        - .obsm['morphology']: Mean morphological features
        - .obs: Islet metadata including donor info
    """

    if region == 'separate':
        print("=== Creating SEPARATE islet and peri-islet datasets ===")
        print(f"Filters: min_cells={min_cells}, require_paired={require_paired}")

        # First, identify valid islets (have both core and peri, meet min_cells)
        parent_col = 'Parent' if 'Parent' in adata_subset.obs.columns else 'Name'
        parent_str = adata_subset.obs[parent_col].astype(str).str.strip()

        # Get base islet ID for each cell
        def get_base_id(p):
            if isinstance(p, str) and '_exp20um' in p:
                return p.replace('_exp20um', '')
            return p

        is_peri = parent_str.str.endswith('_exp20um')
        is_core = parent_str.str.startswith('Islet_') & ~is_peri

        # Add base_islet_id and imageid to create unique islet key
        base_ids = parent_str.apply(get_base_id)
        islet_keys = adata_subset.obs['imageid'].astype(str) + '_' + base_ids.astype(str)

        # Count cells per islet for core and peri
        core_counts = islet_keys[is_core].value_counts()
        peri_counts = islet_keys[is_peri].value_counts()

        print(f"\nBefore filtering:")
        print(f"  Core islets: {len(core_counts)}")
        print(f"  Peri-islet regions: {len(peri_counts)}")

        # Find islets meeting criteria
        if require_paired:
            # Must have both core and peri
            paired_islets = set(core_counts.index) & set(peri_counts.index)
            print(f"  Paired (have both core & peri): {len(paired_islets)}")
        else:
            paired_islets = set(core_counts.index) | set(peri_counts.index)

        # Apply min_cells filter
        valid_core = set(core_counts[core_counts > min_cells].index)
        valid_peri = set(peri_counts[peri_counts > min_cells].index)

        print(f"  Core with >{min_cells} cells: {len(valid_core)}")
        print(f"  Peri with >{min_cells} cells: {len(valid_peri)}")

        # Final valid set: paired AND both regions meet min_cells
        valid_islet_ids = paired_islets & valid_core & valid_peri
        print(f"\nAfter filtering: {len(valid_islet_ids)} valid islets")

        n_dropped_unpaired = len(paired_islets if not require_paired else (set(core_counts.index) | set(peri_counts.index))) - len(paired_islets if require_paired else paired_islets)
        n_dropped_small = len(paired_islets) - len(valid_islet_ids)
        print(f"  Dropped (unpaired): {len(set(core_counts.index) ^ set(peri_counts.index)) if require_paired else 0}")
        print(f"  Dropped (≤{min_cells} cells): {n_dropped_small}")

        # Now create datasets with only valid islets
        adata_islet = create_islet_dataset_fixed(
            adata_subset, region='islet_only', min_cells=min_cells,
            _skip_diagnostic=False, _valid_islet_ids=valid_islet_ids
        )
        print("\n" + "-"*60 + "\n")
        adata_peri = create_islet_dataset_fixed(
            adata_subset, region='peri_only', min_cells=min_cells,
            _skip_diagnostic=True, _valid_islet_ids=valid_islet_ids
        )
        return adata_islet, adata_peri

    morpho_features = [
        'Nucleus Area', 'Nucleus Length', 'Nucleus Circularity', 'Nucleus Solidity',
        'Nucleus: Max diameter µm', 'Nucleus: Min diameter µm',
        'Cell Area', 'Cell Length', 'Cell Circularity', 'Cell Solidity',
        'Cell: Max diameter µm', 'Cell: Min diameter µm',
        'Dist to Closest Tissue', 'Dist to Closest Capillary',
        'Dist to Closest Lymphatic', 'Dist to Closest Nerve'
    ]

    # Get parent column
    if 'Parent' in adata_subset.obs.columns:
        parent_col = 'Parent'
    elif 'Name' in adata_subset.obs.columns:
        parent_col = 'Name'
    else:
        raise KeyError("adata_subset.obs must contain a 'Parent' or 'Name' column.")

    # Convert to string and strip whitespace
    parent_str = adata_subset.obs[parent_col].astype(str).str.strip()

    # DIAGNOSTIC: Show unique Parent value patterns (only on first call)
    if not _skip_diagnostic:
        print(f"\n=== DIAGNOSTIC: Parent column analysis ===")
        print(f"Total cells in input: {adata_subset.n_obs:,}")
        unique_parents = parent_str.unique()
        print(f"Unique Parent values: {len(unique_parents)}")

        # Sample some Parent values to understand the format
        sample_parents = parent_str.value_counts().head(15)
        print(f"\nTop 15 most common Parent values:")
        for parent_val, count in sample_parents.items():
            print(f"  '{parent_val}': {count:,} cells")

        # Check if phenotype column exists and what values it has
        if 'phenotype' in adata_subset.obs.columns:
            pheno_counts = adata_subset.obs['phenotype'].value_counts()
            print(f"\nTotal phenotype distribution (ALL cells in input):")
            for pheno, count in pheno_counts.head(10).items():
                print(f"  {pheno}: {count:,} ({count/adata_subset.n_obs:.1%})")

    # Identify core vs peri-islet cells using explicit string matching
    # Core cells: Parent starts with "Islet_" and does NOT end with "_exp20um"
    # Peri cells: Parent ends with "_exp20um"
    is_peri = parent_str.str.endswith('_exp20um')
    is_core = parent_str.str.startswith('Islet_') & ~is_peri

    # Also track "other" cells that don't match either pattern
    is_other = ~is_core & ~is_peri

    n_core = is_core.sum()
    n_peri = is_peri.sum()
    n_other = is_other.sum()

    # Show filtering summary (only on first call)
    if not _skip_diagnostic:
        print(f"\n=== Cell filtering results ===")
        print(f"Core islet cells (Islet_N pattern): {n_core:,} ({n_core/adata_subset.n_obs:.1%})")
        print(f"Peri-islet cells (_exp20um suffix): {n_peri:,} ({n_peri/adata_subset.n_obs:.1%})")
        print(f"Other cells (no islet pattern): {n_other:,} ({n_other/adata_subset.n_obs:.1%})")

        # Show phenotype breakdown for each category
        if 'phenotype' in adata_subset.obs.columns:
            print(f"\n=== Phenotype distribution BY FILTER CATEGORY ===")

            if n_core > 0:
                core_phenos = adata_subset.obs.loc[is_core, 'phenotype'].value_counts()
                print(f"\nCORE ISLET cells ({n_core:,} cells):")
                for pheno, count in core_phenos.head(8).items():
                    print(f"  {pheno}: {count:,} ({count/n_core:.1%})")

            if n_peri > 0:
                peri_phenos = adata_subset.obs.loc[is_peri, 'phenotype'].value_counts()
                print(f"\nPERI-ISLET cells ({n_peri:,} cells):")
                for pheno, count in peri_phenos.head(8).items():
                    print(f"  {pheno}: {count:,} ({count/n_peri:.1%})")

            if n_other > 0:
                other_phenos = adata_subset.obs.loc[is_other, 'phenotype'].value_counts()
                print(f"\nOTHER cells ({n_other:,} cells):")
                for pheno, count in other_phenos.head(8).items():
                    print(f"  {pheno}: {count:,} ({count/n_other:.1%})")

    # Select cells based on region parameter
    if region == 'islet_only':
        cell_mask = is_core
        region_label = 'islet'
    elif region == 'peri_only':
        cell_mask = is_peri
        region_label = 'peri_islet'
    elif region == 'combined':
        cell_mask = is_core | is_peri
        region_label = 'combined'
    else:
        raise ValueError(f"region must be 'islet_only', 'peri_only', 'combined', or 'separate', got '{region}'")

    if not cell_mask.any():
        sample_values = parent_str.unique()[:20].tolist()
        raise ValueError(
            f"No cells found for region='{region}'. "
            f"Sample Parent values: {sample_values}"
        )

    # Store parent_str in the dataframe for use after filtering
    adata_subset.obs['_parent_str'] = parent_str

    # Filter to selected cells
    islet_cells = adata_subset[cell_mask].copy()
    print(f"\n=== SELECTED {islet_cells.n_obs:,} cells for region='{region}' ===")

    # Create islet IDs from Parent column
    def get_base_islet_id(parent_name):
        if isinstance(parent_name, str) and '_exp20um' in parent_name:
            return parent_name.replace('_exp20um', '')
        return parent_name

    islet_cells.obs['base_islet_id'] = islet_cells.obs['_parent_str'].apply(get_base_islet_id)
    islet_cells.obs['islet_id'] = (
        islet_cells.obs['imageid'].astype(str) + '_' + islet_cells.obs['base_islet_id'].astype(str)
    )

    # Filter to valid islets if specified (from region='separate' pre-filtering)
    if _valid_islet_ids is not None:
        n_before = islet_cells.n_obs
        islet_cells = islet_cells[islet_cells.obs['islet_id'].isin(_valid_islet_ids)].copy()
        print(f"Filtered to valid islets: {n_before:,} -> {islet_cells.n_obs:,} cells")

    unique_islets = islet_cells.obs['islet_id'].unique().tolist()
    print(f"Total unique islets: {len(unique_islets)}")

    # Get phenotype information
    phenotype_series = islet_cells.obs.get('phenotype')
    unique_phenotypes = []
    if phenotype_series is not None:
        # Convert to string to avoid categorical/string type mismatches
        phenotype_series = phenotype_series.astype(str)
        islet_cells.obs['phenotype'] = phenotype_series  # Update the column
        unique_phenotypes = sorted(phenotype_series.dropna().unique().tolist())
    print(f"Total unique phenotypes: {len(unique_phenotypes)}")

    # Print phenotype distribution for this region (VERIFICATION)
    if phenotype_series is not None:
        print(f"\n=== VERIFICATION: Phenotype distribution for region='{region}' ===")
        print(f"(Calculated from filtered islet_cells with {islet_cells.n_obs:,} cells)")

        # Double-check: recalculate directly from islet_cells.obs
        direct_dist = islet_cells.obs['phenotype'].astype(str).value_counts(normalize=True)
        print(f"Top 10 phenotypes:")
        for pheno, prop in direct_dist.head(10).items():
            count = int(prop * islet_cells.n_obs)
            print(f"  {pheno}: {prop:.1%} ({count:,} cells)")

        # Check for suspicious values - only warn for islet_only, not peri_only
        if region == 'islet_only' and 'Acinar' in direct_dist.index and direct_dist['Acinar'] > 0.30:
            print(f"\n⚠️  WARNING: Acinar is {direct_dist['Acinar']:.1%} - this is HIGH for islet core!")
            print("   Islet core should be dominated by endocrine cells (Beta, Alpha, Delta).")
            print("   This suggests the filtering is not correctly identifying islet core cells.")
            print("   Check the Parent column values above.")
        elif region == 'peri_only' and 'Acinar' in direct_dist.index and direct_dist['Acinar'] > 0.50:
            print(f"\n✓ Acinar is {direct_dist['Acinar']:.1%} - this is EXPECTED for peri-islet regions.")
            print("   (Peri-islet = surrounding exocrine tissue, which is mostly Acinar)")

    protein_markers = islet_cells.var_names.tolist()
    available_morpho = [feat for feat in morpho_features if feat in islet_cells.obs.columns]
    print(f"\nProtein markers: {len(protein_markers)}")
    print(f"Available morphological features: {len(available_morpho)}")

    # Initialize storage
    protein_data = []
    pheno_counts_data = []
    morpho_data = []
    islet_metadata = []
    scvi_means = []

    print(f"\nAggregating {len(unique_islets)} islets...")

    for i, islet_id in enumerate(unique_islets):
        if i % 1000 == 0 and i > 0:
            print(f"  Processed {i}/{len(unique_islets)} islets")

        islet_subset = islet_cells[islet_cells.obs['islet_id'] == islet_id]

        # Mean protein expression
        matrix = islet_subset.X
        if sparse.issparse(matrix):
            mean_protein = np.asarray(matrix.mean(axis=0)).ravel()
        else:
            mean_protein = np.asarray(matrix.mean(axis=0)).ravel()
        protein_data.append(mean_protein)

        # Phenotype counts
        if len(unique_phenotypes) > 0:
            pheno_series = islet_subset.obs['phenotype'].astype(str)
            pheno_counts = pheno_series.value_counts()
            pheno_totals = pheno_counts.reindex(unique_phenotypes, fill_value=0).to_numpy(dtype=float)
        else:
            pheno_totals = np.array([], dtype=float)
        pheno_counts_data.append(pheno_totals)

        # Mean morphological features
        morpho_means = []
        for feat in available_morpho:
            feat_values = pd.to_numeric(islet_subset.obs[feat], errors='coerce')
            if feat_values.notna().any():
                morpho_means.append(float(np.nanmean(feat_values)))
            else:
                morpho_means.append(np.nan)
        morpho_data.append(np.asarray(morpho_means, dtype=float))

        # Metadata
        first_cell = islet_subset.obs.iloc[0]
        islet_metadata.append({
            'islet_id': islet_id,
            'base_islet_id': first_cell['base_islet_id'],
            'region': region_label,  # Track which region this data comes from
            'donor_status': first_cell['Donor Status'],
            'age': first_cell['Age'],
            'gender': first_cell['Gender'],
            'imageid': first_cell['imageid'],
            'total_cells': int(islet_subset.n_obs),
            'GADA': first_cell.get('GADA', False) if not pd.isna(first_cell.get('GADA')) else False,
            'ZnT8A': first_cell.get('ZnT8A', False) if not pd.isna(first_cell.get('ZnT8A')) else False,
            'IA2A': first_cell.get('IA2A', False) if not pd.isna(first_cell.get('IA2A')) else False,
            'mIAA': first_cell.get('mIAA', False) if not pd.isna(first_cell.get('mIAA')) else False
        })

        # Mean scVI latent representation
        if 'X_scVI' in islet_subset.obsm:
            scvi_array = np.asarray(islet_subset.obsm['X_scVI'], dtype=float)
            scvi_means.append(scvi_array.mean(axis=0))
        else:
            scvi_means.append(None)

    # Build matrices
    protein_matrix = np.vstack(protein_data).astype(np.float32)

    morpho_matrix = (
        np.vstack(morpho_data).astype(np.float32)
        if morpho_data and available_morpho
        else np.zeros((len(protein_matrix), 0), dtype=np.float32)
    )

    pheno_counts_matrix = (
        np.vstack(pheno_counts_data).astype(np.float32)
        if pheno_counts_data and len(unique_phenotypes) > 0
        else np.zeros((len(protein_matrix), 0), dtype=np.float32)
    )

    islet_metadata_df = pd.DataFrame(islet_metadata)

    print(f"\nFinal dataset: {protein_matrix.shape[0]} islets")
    print(f"  Protein features: {protein_matrix.shape[1]}")
    print(f"  Morphology features: {morpho_matrix.shape[1]}")
    print(f"  Phenotype categories: {pheno_counts_matrix.shape[1]}")

    # KEY FIX #2: Only protein data in .X
    adata = ad.AnnData(X=protein_matrix)
    adata.var_names = protein_markers
    adata.var['feature_type'] = 'protein'

    adata.obs_names = islet_metadata_df['islet_id']
    for col in islet_metadata_df.columns:
        adata.obs[col] = islet_metadata_df[col].values

    # KEY FIX #3: Store phenotypes and morphology in .obsm
    if len(unique_phenotypes) > 0:
        adata.obsm['phenotype_counts'] = pheno_counts_matrix
        # Compute proportions
        row_sums = pheno_counts_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        adata.obsm['phenotype_proportions'] = pheno_counts_matrix / row_sums

        # Store phenotype names for reference
        adata.uns['phenotype_names'] = list(unique_phenotypes)

        # Also add proportion columns to obs for easy access/analysis by disease state
        for j, pheno in enumerate(unique_phenotypes):
            adata.obs[f'prop_{pheno}'] = adata.obsm['phenotype_proportions'][:, j]
            adata.obs[f'count_{pheno}'] = adata.obsm['phenotype_counts'][:, j]

    if available_morpho:
        adata.obsm['morphology'] = morpho_matrix
        adata.uns['morphology_names'] = available_morpho

    # KEY FIX #4: Store scVI means in .obsm for use with sc.pp.neighbors
    scvi_template = next((arr for arr in scvi_means if arr is not None), None)
    if scvi_template is not None:
        scvi_dim = scvi_template.shape[0]
        scvi_matrix = np.vstack([
            arr if arr is not None else np.full(scvi_dim, np.nan)
            for arr in scvi_means
        ]).astype(np.float32)
        adata.obsm['X_scVI_mean'] = scvi_matrix
        print(f"  scVI latent dimensions: {scvi_dim}")

    # Set categorical types
    adata.obs['donor_status'] = adata.obs['donor_status'].astype('category')
    adata.obs['gender'] = adata.obs['gender'].astype('category')
    adata.obs['region'] = adata.obs['region'].astype('category')

    # Print phenotype distribution summary
    print("\nDonor status distribution:")
    print(adata.obs['donor_status'].value_counts())

    if 'dominant_phenotype' in adata.obs.columns:
        print("\nDominant phenotype distribution:")
        print(adata.obs['dominant_phenotype'].value_counts().head(10))

    return adata


def compute_neighbors_and_umap(adata, use_scvi=True, n_neighbors=15):
    """
    Compute neighbors and UMAP using the scVI latent space.

    Parameters
    ----------
    adata : AnnData
        Islet-level AnnData from create_islet_dataset_fixed()
    use_scvi : bool
        If True, use X_scVI_mean representation (recommended)
    n_neighbors : int
        Number of neighbors for graph construction
    """
    import scanpy as sc

    if use_scvi and 'X_scVI_mean' in adata.obsm:
        print("Computing neighbors using scVI latent space...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_scVI_mean', metric='cosine')
    else:
        print("Computing neighbors using protein expression (PCA)...")
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)

    sc.tl.umap(adata, min_dist=0.1, spread=1.5)
    print("UMAP computed.")

    return adata


def compute_trajectory(adata, root_marker='INS', root_status='ND'):
    """
    Compute diffusion pseudotime trajectory.

    Parameters
    ----------
    adata : AnnData
        Islet-level AnnData with neighbors computed
    root_marker : str
        Protein marker to use for finding the root cell (highest expression)
    root_status : str
        Donor status to restrict root search to (e.g., 'ND' for healthy)
    """
    import scanpy as sc

    # Find root: highest marker expression in specified donor status
    if root_marker in adata.var_names:
        marker_idx = list(adata.var_names).index(root_marker)
        status_mask = adata.obs['donor_status'] == root_status

        if status_mask.any():
            marker_data = adata.X[status_mask, marker_idx]
            status_indices = np.where(status_mask)[0]
            best_idx = status_indices[np.argmax(marker_data)]

            adata.uns['iroot'] = best_idx
            print(f"Root set to islet {best_idx} ({root_status}, {root_marker}={adata.X[best_idx, marker_idx]:.3f})")
        else:
            print(f"Warning: No islets with donor_status='{root_status}', using global max")
            adata.uns['iroot'] = np.argmax(adata.X[:, marker_idx])
    else:
        print(f"Warning: {root_marker} not found, using first islet as root")
        adata.uns['iroot'] = 0

    # Compute diffusion map and pseudotime
    sc.tl.diffmap(adata, n_comps=10)
    sc.tl.dpt(adata)

    # Add diffusion components to obs
    adata.obs['DC1'] = adata.obsm['X_diffmap'][:, 0]
    adata.obs['DC2'] = adata.obsm['X_diffmap'][:, 1]

    # Print pseudotime statistics
    print("\nPseudotime by donor status:")
    for status in adata.obs['donor_status'].cat.categories:
        mask = adata.obs['donor_status'] == status
        mean_pt = adata.obs.loc[mask, 'dpt_pseudotime'].mean()
        std_pt = adata.obs.loc[mask, 'dpt_pseudotime'].std()
        print(f"  {status}: {mean_pt:.3f} ± {std_pt:.3f}")

    return adata


def compare_islet_vs_peri(adata_islet, adata_peri, markers=None):
    """
    Compare islet core vs peri-islet microenvironment.

    Parameters
    ----------
    adata_islet : AnnData
        Islet-level data from region='islet_only'
    adata_peri : AnnData
        Peri-islet data from region='peri_only'
    markers : list, optional
        Protein markers to compare. If None, uses all shared markers.

    Returns
    -------
    pd.DataFrame
        Comparison statistics for each marker and phenotype
    """
    from scipy import stats

    results = []

    # Compare protein expression
    shared_markers = list(set(adata_islet.var_names) & set(adata_peri.var_names))
    if markers:
        shared_markers = [m for m in markers if m in shared_markers]

    print("=== PROTEIN EXPRESSION COMPARISON ===")
    print(f"Comparing {len(shared_markers)} markers\n")

    for marker in shared_markers:
        islet_expr = adata_islet[:, marker].X.flatten()
        peri_expr = adata_peri[:, marker].X.flatten()

        # Mann-Whitney U test (non-parametric)
        stat, pval = stats.mannwhitneyu(islet_expr, peri_expr, alternative='two-sided')

        fold_change = np.mean(islet_expr) / (np.mean(peri_expr) + 1e-10)

        results.append({
            'feature': marker,
            'type': 'protein',
            'islet_mean': np.mean(islet_expr),
            'peri_mean': np.mean(peri_expr),
            'fold_change_islet_vs_peri': fold_change,
            'pvalue': pval
        })

    # Compare phenotype proportions
    if 'phenotype_names' in adata_islet.uns and 'phenotype_names' in adata_peri.uns:
        print("\n=== PHENOTYPE PROPORTION COMPARISON ===")

        islet_phenos = set(adata_islet.uns['phenotype_names'])
        peri_phenos = set(adata_peri.uns['phenotype_names'])
        shared_phenos = list(islet_phenos & peri_phenos)

        for pheno in shared_phenos:
            islet_idx = adata_islet.uns['phenotype_names'].index(pheno)
            peri_idx = adata_peri.uns['phenotype_names'].index(pheno)

            islet_prop = adata_islet.obsm['phenotype_proportions'][:, islet_idx]
            peri_prop = adata_peri.obsm['phenotype_proportions'][:, peri_idx]

            stat, pval = stats.mannwhitneyu(islet_prop, peri_prop, alternative='two-sided')

            results.append({
                'feature': pheno,
                'type': 'phenotype',
                'islet_mean': np.mean(islet_prop),
                'peri_mean': np.mean(peri_prop),
                'fold_change_islet_vs_peri': np.mean(islet_prop) / (np.mean(peri_prop) + 1e-10),
                'pvalue': pval
            })

    results_df = pd.DataFrame(results)

    # Add multiple testing correction
    from statsmodels.stats.multitest import multipletests
    _, results_df['pvalue_adj'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')

    # Sort by fold change
    results_df = results_df.sort_values('fold_change_islet_vs_peri', ascending=False)

    # Print top differences
    print("\nTop features enriched in ISLET CORE:")
    top_islet = results_df[results_df['fold_change_islet_vs_peri'] > 1].head(10)
    for _, row in top_islet.iterrows():
        sig = "***" if row['pvalue_adj'] < 0.001 else "**" if row['pvalue_adj'] < 0.01 else "*" if row['pvalue_adj'] < 0.05 else ""
        print(f"  {row['feature']}: {row['fold_change_islet_vs_peri']:.2f}x {sig}")

    print("\nTop features enriched in PERI-ISLET:")
    top_peri = results_df[results_df['fold_change_islet_vs_peri'] < 1].tail(10).iloc[::-1]
    for _, row in top_peri.iterrows():
        sig = "***" if row['pvalue_adj'] < 0.001 else "**" if row['pvalue_adj'] < 0.01 else "*" if row['pvalue_adj'] < 0.05 else ""
        print(f"  {row['feature']}: {1/row['fold_change_islet_vs_peri']:.2f}x {sig}")

    return results_df


def merge_islet_peri_datasets(adata_islet, adata_peri):
    """
    Merge islet and peri-islet datasets for joint analysis.

    The merged dataset has islets matched by base_islet_id, with separate
    columns for islet vs peri-islet features.

    Parameters
    ----------
    adata_islet : AnnData
        Islet core data
    adata_peri : AnnData
        Peri-islet data

    Returns
    -------
    AnnData
        Merged dataset with islet features in .X and peri features in .obsm
    """
    import anndata as ad

    # Find matching islets
    islet_ids = set(adata_islet.obs['base_islet_id'].astype(str) + '_' + adata_islet.obs['imageid'].astype(str))
    peri_ids = set(adata_peri.obs['base_islet_id'].astype(str) + '_' + adata_peri.obs['imageid'].astype(str))
    shared_ids = islet_ids & peri_ids

    print(f"Islets with both core and peri-islet data: {len(shared_ids)}")

    # Create matching key
    adata_islet.obs['match_key'] = adata_islet.obs['base_islet_id'].astype(str) + '_' + adata_islet.obs['imageid'].astype(str)
    adata_peri.obs['match_key'] = adata_peri.obs['base_islet_id'].astype(str) + '_' + adata_peri.obs['imageid'].astype(str)

    # Filter to shared
    adata_islet_shared = adata_islet[adata_islet.obs['match_key'].isin(shared_ids)].copy()
    adata_peri_shared = adata_peri[adata_peri.obs['match_key'].isin(shared_ids)].copy()

    # Sort both by match_key to align
    islet_order = adata_islet_shared.obs['match_key'].argsort()
    peri_order = adata_peri_shared.obs['match_key'].argsort()

    adata_islet_sorted = adata_islet_shared[islet_order].copy()
    adata_peri_sorted = adata_peri_shared[peri_order].copy()

    # Use islet as base, add peri data
    adata_merged = adata_islet_sorted.copy()

    # Add peri-islet protein expression
    adata_merged.obsm['X_peri'] = adata_peri_sorted.X

    # Add peri-islet phenotypes
    if 'phenotype_proportions' in adata_peri_sorted.obsm:
        adata_merged.obsm['phenotype_proportions_peri'] = adata_peri_sorted.obsm['phenotype_proportions']
        adata_merged.obsm['phenotype_counts_peri'] = adata_peri_sorted.obsm['phenotype_counts']

    # Add peri-islet cell counts
    adata_merged.obs['total_cells_peri'] = adata_peri_sorted.obs['total_cells'].values
    adata_merged.obs['total_cells_islet'] = adata_merged.obs['total_cells']

    # Rename for clarity
    adata_merged.uns['peri_var_names'] = list(adata_peri_sorted.var_names)
    if 'phenotype_names' in adata_peri_sorted.uns:
        adata_merged.uns['phenotype_names_peri'] = adata_peri_sorted.uns['phenotype_names']

    print(f"Merged dataset: {adata_merged.n_obs} islets with paired core/peri data")

    return adata_merged


# Example usage (illustrative demo only — the Senior pipeline drives aggregation
# via scripts/senior/aggregate_islets.py, which imports the functions above).
if __name__ == "__main__":
    import scanpy as sc
    import os

    # Load your scVI-processed single-cell data (override with ISLET_BASE_DIR)
    base_dir = os.environ.get("ISLET_BASE_DIR", os.path.expanduser("~/IO60panc2nd"))
    adata_scvi = sc.read_h5ad(os.path.join(base_dir, 'single_cell_analysis',
                                            'CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad'))

    # ==========================================================================
    # OPTION 1: Create separate islet and peri-islet datasets
    # ==========================================================================
    print("\n" + "="*60)
    print("CREATING SEPARATE ISLET AND PERI-ISLET DATASETS")
    print("="*60)

    adata_islet, adata_peri = create_islet_dataset_fixed(adata_scvi, region='separate')

    # Compare the two regions
    comparison_df = compare_islet_vs_peri(adata_islet, adata_peri)
    comparison_df.to_csv('islet_vs_peri_comparison.csv', index=False)

    # Analyze islet core (for T1D progression trajectory)
    print("\n" + "="*60)
    print("ISLET CORE ANALYSIS")
    print("="*60)
    adata_islet = compute_neighbors_and_umap(adata_islet, use_scvi=True)
    adata_islet = compute_trajectory(adata_islet, root_marker='INS', root_status='ND')
    adata_islet.write_h5ad('islets_core_only.h5ad')

    # Analyze peri-islet microenvironment separately
    print("\n" + "="*60)
    print("PERI-ISLET MICROENVIRONMENT ANALYSIS")
    print("="*60)
    adata_peri = compute_neighbors_and_umap(adata_peri, use_scvi=True)
    adata_peri.write_h5ad('islets_peri_only.h5ad')

    # ==========================================================================
    # OPTION 2: Merge for paired analysis
    # ==========================================================================
    print("\n" + "="*60)
    print("CREATING MERGED DATASET FOR PAIRED ANALYSIS")
    print("="*60)
    adata_merged = merge_islet_peri_datasets(adata_islet, adata_peri)
    adata_merged.write_h5ad('islets_merged_paired.h5ad')

    # ==========================================================================
    # OPTION 3: Combined analysis (original behavior, if needed)
    # ==========================================================================
    # adata_combined = create_islet_dataset_fixed(adata_scvi, region='combined')

    # Visualize
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    sc.pl.umap(adata_islet, color='donor_status', ax=axes[0, 0], show=False,
               title='Islet Core: Donor Status')
    sc.pl.umap(adata_islet, color='dominant_phenotype', ax=axes[0, 1], show=False,
               title='Islet Core: Dominant Phenotype')
    sc.pl.umap(adata_islet, color='dpt_pseudotime', ax=axes[0, 2], show=False,
               title='Islet Core: Pseudotime')

    sc.pl.umap(adata_peri, color='donor_status', ax=axes[1, 0], show=False,
               title='Peri-islet: Donor Status')
    sc.pl.umap(adata_peri, color='dominant_phenotype', ax=axes[1, 1], show=False,
               title='Peri-islet: Dominant Phenotype')
    sc.pl.umap(adata_peri, color='total_cells', ax=axes[1, 2], show=False,
               title='Peri-islet: Cell Count')

    plt.tight_layout()
    plt.savefig('islet_vs_peri_umap.png', dpi=150)
    plt.show()
