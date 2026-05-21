#!/usr/bin/env python
"""
Retrain scVI without Age/Gender covariates and compare diagnostics.

Purpose: The original scVI model used categorical_covariate_keys=["Age", "Gender"],
which uniquely identified each of the 15 donors (each has a unique Age value).
This script retrains with NO covariates to compare latent space quality.

Architecture kept identical:
  - n_latent=10, n_layers=2, n_hidden=128, dropout=0.1
  - gene_likelihood='nb', dispersion='gene-batch'
  - layer="counts" (recovered raw MFI from log layer)
  - Training: max_epochs=400, early_stopping, patience=45, lr=0.001, batch_size=2048
"""

import os
import sys
import time
import json
import numpy as np
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings("ignore")

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Save all output to a file
RESULTS_DIR = "scripts/scvi_comparison_results"
os.makedirs(RESULTS_DIR, exist_ok=True)
log_file = open(os.path.join(RESULTS_DIR, "diagnostics_log.txt"), "w")

def log(msg=""):
    print(msg)
    sys.stdout.flush()
    log_file.write(msg + "\n")
    log_file.flush()

# ============================================================
# STEP 1: Load single-cell data and recover raw MFI
# ============================================================
log("=" * 70)
log("STEP 1: Loading single-cell data")
log("=" * 70)

t0 = time.time()
sc_path = "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
adata = ad.read_h5ad(sc_path)
log(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars} proteins ({time.time()-t0:.1f}s)")
log(f"Layers: {list(adata.layers.keys())}")

# Recover raw MFI from log layer: log layer = log1p(raw_mfi)
log("\nRecovering raw MFI from log layer...")
raw_mfi = np.expm1(adata.layers['log'].copy())
if hasattr(raw_mfi, 'toarray'):
    raw_mfi = raw_mfi.toarray()
raw_mfi = raw_mfi.astype(np.float32)
adata.layers["counts"] = raw_mfi
log(f"  Raw MFI range: {raw_mfi.min():.2f} to {raw_mfi.max():.2f}")
log(f"  Mean: {raw_mfi.mean():.2f}")

# Store original scVI latent for comparison
original_latent = adata.obsm['X_scVI'].copy()
log(f"\nOriginal X_scVI shape: {original_latent.shape}")

# ============================================================
# STEP 2: Setup and train scVI WITHOUT covariates
# ============================================================
log("\n" + "=" * 70)
log("STEP 2: Training scVI WITHOUT Age/Gender covariates")
log("=" * 70)

import scvi
import torch

scvi.settings.seed = 0
log(f"scvi-tools version: {scvi.__version__}")
log(f"CUDA: {torch.cuda.is_available()}, GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A'}")

# Setup WITHOUT any covariates (the ONLY change from original)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    # NO batch_key
    # NO categorical_covariate_keys  <-- THIS IS THE CHANGE
    # NO continuous_covariate_keys
)

log("\nscVI setup complete (NO covariates)")
log(f"  Cells: {adata.n_obs:,}, Proteins: {adata.n_vars}")
log(f"  _scvi_batch unique: {adata.obs['_scvi_batch'].unique()}")

# Create model with IDENTICAL architecture
model = scvi.model.SCVI(
    adata,
    n_latent=10,
    n_layers=2,
    n_hidden=128,
    dropout_rate=0.1,
    dispersion='gene-batch',
    gene_likelihood='nb',
)
log(f"\nModel created")

# Train with same hyperparameters
t_train = time.time()
log("\nStarting training...")
model.train(
    max_epochs=400,
    accelerator="gpu" if torch.cuda.is_available() else "cpu",
    devices=1,
    batch_size=2048,
    train_size=0.9,
    early_stopping=True,
    early_stopping_patience=45,
    check_val_every_n_epoch=1,
    plan_kwargs={"lr": 0.001},
    datasplitter_kwargs={
        "num_workers": min(16, os.cpu_count()),
        "pin_memory": True
    }
)
train_time = time.time() - t_train
log(f"\nTraining complete in {train_time:.1f}s ({train_time/60:.1f} min)")

# ============================================================
# STEP 3: Extract new latent representations
# ============================================================
log("\n" + "=" * 70)
log("STEP 3: Extracting latent representations")
log("=" * 70)

new_latent = model.get_latent_representation()
log(f"New latent shape: {new_latent.shape}")
log(f"Original latent shape: {original_latent.shape}")

# Save latent representations to disk
np.save(os.path.join(RESULTS_DIR, "new_latent_single_cell.npy"), new_latent)
log(f"Saved new single-cell latent to {RESULTS_DIR}/new_latent_single_cell.npy")

# ============================================================
# STEP 4: Aggregate to islet level (same as fixed_islet_aggregation.py)
# ============================================================
log("\n" + "=" * 70)
log("STEP 4: Islet-level aggregation")
log("=" * 70)

# We need core-only cells (region='islet_only'), min_cells=20
parent = adata.obs['Parent'].astype(str)

# Core cells: Parent starts with Islet_ but NOT ending in _exp20um
core_mask = parent.str.startswith('Islet_') & ~parent.str.endswith('_exp20um')
log(f"Core cells: {core_mask.sum():,} / {adata.n_obs:,}")

core_obs = adata.obs[core_mask].copy()
core_obs['islet_id'] = core_obs['imageid'].astype(str) + '_' + core_obs['Parent'].astype(str)

# Get both old and new latent for core cells
core_idx = np.where(core_mask)[0]
core_new_latent = new_latent[core_idx]
core_orig_latent = original_latent[core_idx]

# Aggregate: mean latent per islet
islet_ids = core_obs['islet_id'].values
unique_islets = np.unique(islet_ids)
log(f"Unique islets (before min_cells filter): {len(unique_islets)}")

# Count cells per islet and filter
islet_counts = pd.Series(islet_ids).value_counts()
valid_islets = islet_counts[islet_counts >= 20].index
log(f"Islets with >= 20 cells: {len(valid_islets)}")

# Aggregate latent means
new_means = []
orig_means = []
islet_metadata = []

for islet in valid_islets:
    mask = islet_ids == islet
    new_means.append(core_new_latent[mask].mean(axis=0))
    orig_means.append(core_orig_latent[mask].mean(axis=0))

    first_idx = np.where(mask)[0][0]
    row = core_obs.iloc[first_idx]
    islet_metadata.append({
        'islet_id': islet,
        'imageid': str(row['imageid']),
        'donor_status': row['Donor Status'],
    })

new_islet_latent = np.array(new_means)
orig_islet_latent = np.array(orig_means)
meta_df = pd.DataFrame(islet_metadata)

log(f"\nAggregated to {len(meta_df)} islets")
log(f"New islet latent shape: {new_islet_latent.shape}")
log(f"Original islet latent shape: {orig_islet_latent.shape}")

# Save islet-level results
np.save(os.path.join(RESULTS_DIR, "new_islet_latent.npy"), new_islet_latent)
np.save(os.path.join(RESULTS_DIR, "orig_islet_latent.npy"), orig_islet_latent)
meta_df.to_csv(os.path.join(RESULTS_DIR, "islet_metadata.csv"), index=False)
log(f"Saved islet-level latents and metadata to {RESULTS_DIR}/")

# Load the app H5AD to get pseudotime for correlation
app_adata = ad.read_h5ad("data/islet_explorer.h5ad")
app_obs = app_adata.obs.copy()
if 'combined_islet_id' in app_obs.columns:
    pt_map = app_obs.set_index('combined_islet_id')['dpt_pseudotime'].to_dict()
elif 'islet_id' in app_obs.columns:
    pt_map = app_obs.set_index('islet_id')['dpt_pseudotime'].to_dict()
else:
    pt_map = app_obs['dpt_pseudotime'].to_dict()

meta_df['dpt_pseudotime'] = meta_df['islet_id'].map(pt_map)
has_pt = meta_df['dpt_pseudotime'].notna().sum()
log(f"Matched pseudotime for {has_pt}/{len(meta_df)} islets")

# ============================================================
# STEP 5: Run diagnostics on BOTH latent spaces
# ============================================================
log("\n" + "=" * 70)
log("STEP 5: COMPARATIVE DIAGNOSTICS")
log("=" * 70)

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

obs = meta_df

def run_diagnostics(latent, obs, label):
    """Run full diagnostic suite on a latent space."""
    log(f"\n{'='*60}")
    log(f"  DIAGNOSTICS: {label}")
    log(f"{'='*60}")

    results = {}

    # --- 1. Variance explained ---
    enc_donor = OneHotEncoder(sparse_output=False, drop='first')
    enc_status = OneHotEncoder(sparse_output=False, drop='first')
    X_donor = enc_donor.fit_transform(obs[['imageid']].astype(str))
    X_status = enc_status.fit_transform(obs[['donor_status']].astype(str))

    reg_d = LinearRegression().fit(X_donor, latent)
    reg_s = LinearRegression().fit(X_status, latent)
    r2_donor = reg_d.score(X_donor, latent)
    r2_status = reg_s.score(X_status, latent)

    log(f"\n  R² (overall):")
    log(f"    imageid:      {r2_donor:.4f}")
    log(f"    donor_status: {r2_status:.4f}")
    log(f"    Ratio:        {r2_donor/max(r2_status, 1e-10):.2f}x")
    results['r2_donor'] = r2_donor
    results['r2_status'] = r2_status

    # --- 2. Eta-squared per dimension ---
    def eta_squared(groups, values):
        grand_mean = np.mean(values)
        ss_total = np.sum((values - grand_mean)**2)
        ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
        return ss_between / ss_total if ss_total > 0 else 0

    eta2_d_list = []
    eta2_s_list = []
    for i in range(latent.shape[1]):
        vals = latent[:, i]
        gd = [vals[obs['imageid'].values == d] for d in obs['imageid'].unique()]
        gs = [vals[obs['donor_status'].values == s] for s in obs['donor_status'].unique()]
        eta2_d_list.append(eta_squared(gd, vals))
        eta2_s_list.append(eta_squared(gs, vals))

    log(f"\n  Mean eta² across dimensions:")
    log(f"    imageid:      {np.mean(eta2_d_list):.4f}")
    log(f"    donor_status: {np.mean(eta2_s_list):.4f}")
    log(f"    Excess donor: {np.mean(eta2_d_list) - np.mean(eta2_s_list):.4f}")
    results['eta2_donor'] = np.mean(eta2_d_list)
    results['eta2_status'] = np.mean(eta2_s_list)
    results['eta2_excess'] = np.mean(eta2_d_list) - np.mean(eta2_s_list)

    # Per-dimension detail
    log(f"\n  Per-dimension eta²:")
    log(f"  {'Dim':<6} {'eta²(donor)':<14} {'eta²(status)':<14} {'Excess':<10}")
    log(f"  {'-'*44}")
    for i in range(latent.shape[1]):
        log(f"  {i:<6} {eta2_d_list[i]:<14.4f} {eta2_s_list[i]:<14.4f} {eta2_d_list[i]-eta2_s_list[i]:<+10.4f}")

    # --- 3. Silhouette scores ---
    le_d = LabelEncoder()
    le_s = LabelEncoder()
    labels_d = le_d.fit_transform(obs['imageid'].astype(str))
    labels_s = le_s.fit_transform(obs['donor_status'].astype(str))

    sil_d = silhouette_score(latent, labels_d, metric='cosine')
    sil_s = silhouette_score(latent, labels_s, metric='cosine')

    log(f"\n  Silhouette (cosine):")
    log(f"    imageid:      {sil_d:.4f}")
    log(f"    donor_status: {sil_s:.4f}")
    log(f"    Ratio:        {sil_d/max(sil_s, 1e-10):.2f}x")
    results['sil_donor'] = sil_d
    results['sil_status'] = sil_s

    # --- 4. PCA: what drives top components ---
    pca = PCA(n_components=min(5, latent.shape[1]))
    pcs = pca.fit_transform(latent)

    log(f"\n  PCA of latent space:")
    log(f"  {'PC':<5} {'VarExpl':<10} {'eta²(donor)':<14} {'eta²(status)':<14} {'r(pseudotime)':<15}")
    log(f"  {'-'*58}")

    pt_vals = obs['dpt_pseudotime'].values
    for i in range(min(5, pcs.shape[1])):
        pc = pcs[:, i]
        gd = [pc[obs['imageid'].values == d] for d in obs['imageid'].unique()]
        gs = [pc[obs['donor_status'].values == s] for s in obs['donor_status'].unique()]
        e2d = eta_squared(gd, pc)
        e2s = eta_squared(gs, pc)

        valid = ~np.isnan(pt_vals)
        if valid.sum() > 0:
            r_pt, _ = spearmanr(pc[valid], pt_vals[valid])
            pt_str = f"{r_pt:.3f}"
        else:
            pt_str = "N/A"

        log(f"  PC{i:<3} {pca.explained_variance_ratio_[i]:<10.4f} {e2d:<14.4f} {e2s:<14.4f} {pt_str:<15}")

    results['pc0_var'] = pca.explained_variance_ratio_[0]

    # --- 5. k-NN donor mixing ---
    k = 15
    nn = NearestNeighbors(n_neighbors=k+1, metric='cosine')
    nn.fit(latent)
    _, indices = nn.kneighbors(latent)

    donors = obs['imageid'].values
    statuses = obs['donor_status'].values
    same_d = np.array([np.mean(donors[indices[i, 1:]] == donors[i]) for i in range(len(latent))])
    same_s = np.array([np.mean(statuses[indices[i, 1:]] == statuses[i]) for i in range(len(latent))])

    n_total = len(obs)
    expected_d = np.mean([(obs['imageid'] == d).sum() / n_total for d in obs['imageid'].unique()])
    expected_s = np.mean([(obs['donor_status'] == s).sum() / n_total for s in obs['donor_status'].unique()])

    log(f"\n  k-NN mixing (k={k}):")
    log(f"    Same donor:  {np.mean(same_d):.4f} (random: {expected_d:.4f}, enrichment: {np.mean(same_d)/expected_d:.2f}x)")
    log(f"    Same status: {np.mean(same_s):.4f} (random: {expected_s:.4f}, enrichment: {np.mean(same_s)/expected_s:.2f}x)")
    results['knn_donor_frac'] = np.mean(same_d)
    results['knn_status_frac'] = np.mean(same_s)
    results['knn_donor_enrichment'] = np.mean(same_d) / expected_d
    results['knn_status_enrichment'] = np.mean(same_s) / expected_s

    # --- 6. Per-disease-group donor mixing ---
    log(f"\n  Per-group k-NN mixing:")
    for status in ['ND', 'Aab+', 'T1D']:
        mask = obs['donor_status'].values == status
        log(f"    {status}: same-donor={np.mean(same_d[mask]):.4f}, same-status={np.mean(same_s[mask]):.4f}")

    return results


# Run on both latent spaces
log("\n\nRunning diagnostics on ORIGINAL latent space (with Age/Gender covariates)...")
orig_results = run_diagnostics(orig_islet_latent, obs, "ORIGINAL (Age+Gender covariates)")

log("\n\nRunning diagnostics on NEW latent space (NO covariates)...")
new_results = run_diagnostics(new_islet_latent, obs, "NEW (No covariates)")

# ============================================================
# STEP 6: Side-by-side comparison
# ============================================================
log("\n\n" + "=" * 70)
log("STEP 6: SIDE-BY-SIDE COMPARISON")
log("=" * 70)

metrics = [
    ('R² imageid (donor)', 'r2_donor', 'lower is better'),
    ('R² donor_status (biology)', 'r2_status', 'higher is better'),
    ('eta² imageid', 'eta2_donor', 'lower is better'),
    ('eta² donor_status', 'eta2_status', 'higher is better'),
    ('eta² excess (donor-status)', 'eta2_excess', 'lower is better'),
    ('Silhouette imageid (cosine)', 'sil_donor', 'lower is better'),
    ('Silhouette donor_status (cosine)', 'sil_status', 'higher is better'),
    ('k-NN same-donor frac', 'knn_donor_frac', 'lower is better'),
    ('k-NN same-status frac', 'knn_status_frac', 'higher is better'),
    ('k-NN donor enrichment', 'knn_donor_enrichment', 'lower is better'),
    ('k-NN status enrichment', 'knn_status_enrichment', 'higher is better'),
]

log(f"\n{'Metric':<35} {'Original':<12} {'New (no cov)':<14} {'Change':<12} {'Better?':<10}")
log("-" * 83)

for name, key, direction in metrics:
    o = orig_results[key]
    n = new_results[key]
    change = n - o

    if 'higher is better' in direction:
        better = "YES" if n > o else ("NO" if n < o else "SAME")
    else:
        better = "YES" if n < o else ("NO" if n > o else "SAME")

    log(f"{name:<35} {o:<12.4f} {n:<14.4f} {change:<+12.4f} {better:<10}")

# ============================================================
# STEP 7: Summary assessment
# ============================================================
log("\n\n" + "=" * 70)
log("STEP 7: SUMMARY ASSESSMENT")
log("=" * 70)

log(f"\nTraining time: {train_time:.1f}s ({train_time/60:.1f} min)")
log(f"\nKey comparisons:")
log(f"  Donor variance explained:  {orig_results['r2_donor']:.4f} -> {new_results['r2_donor']:.4f}")
log(f"  Biology variance explained: {orig_results['r2_status']:.4f} -> {new_results['r2_status']:.4f}")
log(f"  Excess donor effect (eta²): {orig_results['eta2_excess']:.4f} -> {new_results['eta2_excess']:.4f}")
log(f"  Donor silhouette:           {orig_results['sil_donor']:.4f} -> {new_results['sil_donor']:.4f}")
log(f"  k-NN donor enrichment:      {orig_results['knn_donor_enrichment']:.2f}x -> {new_results['knn_donor_enrichment']:.2f}x")

# Verdict
improvements = 0
regressions = 0
for name, key, direction in metrics:
    o = orig_results[key]
    n = new_results[key]
    if 'higher is better' in direction:
        if n > o * 1.01: improvements += 1
        elif n < o * 0.99: regressions += 1
    else:
        if n < o * 0.99: improvements += 1
        elif n > o * 1.01: regressions += 1

log(f"\nScorecard: {improvements} metrics improved, {regressions} metrics regressed, {len(metrics) - improvements - regressions} unchanged")

if new_results['eta2_excess'] < orig_results['eta2_excess']:
    log("\n  -> Removing covariates REDUCED excess donor effects")
else:
    log("\n  -> Removing covariates INCREASED excess donor effects")
    log("     (Age/Gender covariates were helping remove donor signal)")

if new_results['r2_status'] > orig_results['r2_status']:
    log("  -> Biological signal is STRONGER without covariates")
else:
    log("  -> Biological signal is WEAKER without covariates")

# Save results as JSON
all_results = {
    'original': {k: float(v) for k, v in orig_results.items()},
    'new_no_covariates': {k: float(v) for k, v in new_results.items()},
    'training_time_seconds': train_time,
}
with open(os.path.join(RESULTS_DIR, "comparison_results.json"), "w") as f:
    json.dump(all_results, f, indent=2)

log(f"\nResults saved to {RESULTS_DIR}/comparison_results.json")
log("\nDone.")
log_file.close()
