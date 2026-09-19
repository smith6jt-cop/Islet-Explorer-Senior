#!/usr/bin/env python
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan (code tier 1):
nothing in the tree imports it, and the inputs named below are absent from this checkout.
It is kept as the record of how that artifact was built, not as a description of anything
current.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method.

Supplemental Figure: scVI model comparison
  - Original (with Age+Gender categorical covariates)
  - Retrained without covariates

6-panel figure:
  A/B: PCA of latent space colored by donor / disease status
  C:   Per-dimension eta² (variance explained)
  D:   Donor silhouette distributions
  E:   k-NN neighbor mixing
  F:   Summary metrics bar chart
"""

import os
import numpy as np
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples
from sklearn.preprocessing import LabelEncoder
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore")

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = "scripts/scvi_comparison_results"

# ── Load data ──────────────────────────────────────────────
orig_latent = np.load(os.path.join(RESULTS_DIR, "orig_islet_latent.npy"))
new_latent  = np.load(os.path.join(RESULTS_DIR, "new_islet_latent.npy"))
meta        = pd.read_csv(os.path.join(RESULTS_DIR, "islet_metadata.csv"))

with open(os.path.join(RESULTS_DIR, "comparison_results.json")) as f:
    results = json.load(f)

# Pseudotime
import anndata as ad
app = ad.read_h5ad("data/islet_explorer.h5ad")
app_obs = app.obs.copy()
if 'combined_islet_id' in app_obs.columns:
    pt_map = app_obs.set_index('combined_islet_id')['dpt_pseudotime'].to_dict()
elif 'islet_id' in app_obs.columns:
    pt_map = app_obs.set_index('islet_id')['dpt_pseudotime'].to_dict()
else:
    pt_map = app_obs['dpt_pseudotime'].to_dict()
meta['dpt_pseudotime'] = meta['islet_id'].map(pt_map)

# ── Colors ─────────────────────────────────────────────────
STATUS_COLORS = {'ND': '#2ca02c', 'Aab+': '#ffcc00', 'T1D': '#9467bd'}
STATUS_ORDER = ['ND', 'Aab+', 'T1D']
DONOR_CMAP = plt.cm.tab20
donors_unique = sorted(meta['imageid'].unique())
donor_color_map = {d: DONOR_CMAP(i / len(donors_unique)) for i, d in enumerate(donors_unique)}

BLUE = '#2166ac'
RED  = '#b2182b'

# ── Precompute ─────────────────────────────────────────────
pca_orig = PCA(n_components=2).fit_transform(orig_latent)
pca_new  = PCA(n_components=2).fit_transform(new_latent)

def eta_squared(groups, values):
    gm = np.mean(values)
    ss_t = np.sum((values - gm)**2)
    ss_b = sum(len(g) * (np.mean(g) - gm)**2 for g in groups)
    return ss_b / ss_t if ss_t > 0 else 0

def per_dim_eta2(latent, col):
    out = []
    for i in range(latent.shape[1]):
        v = latent[:, i]
        groups = [v[meta[col].values == g] for g in meta[col].unique()]
        out.append(eta_squared(groups, v))
    return np.array(out)

def knn_same_frac(latent, labels, k=15):
    nn = NearestNeighbors(n_neighbors=k+1, metric='cosine').fit(latent)
    _, idx = nn.kneighbors(latent)
    return np.array([np.mean(labels[idx[i, 1:]] == labels[i]) for i in range(len(latent))])

knn_orig_donor  = knn_same_frac(orig_latent, meta['imageid'].values)
knn_new_donor   = knn_same_frac(new_latent, meta['imageid'].values)
knn_orig_status = knn_same_frac(orig_latent, meta['donor_status'].values)
knn_new_status  = knn_same_frac(new_latent, meta['donor_status'].values)

labels_donor = LabelEncoder().fit_transform(meta['imageid'].astype(str))
sil_orig = silhouette_samples(orig_latent, labels_donor, metric='cosine')
sil_new  = silhouette_samples(new_latent, labels_donor, metric='cosine')

# ================================================================
#  FIGURE
# ================================================================
fig, axes = plt.subplots(4, 2, figsize=(14, 22),
                         gridspec_kw={'height_ratios': [1, 1, 1, 1],
                                      'hspace': 0.40, 'wspace': 0.30})

# Turn off the placeholder axes for row 0 — we'll use inset axes
for ax in axes[0]:
    ax.set_visible(False)

lbl_kw = dict(fontsize=18, fontweight='bold', va='top', ha='left',
              transform=fig.transFigure)

# ── Row 0: PCA panels (4 subplots via inset axes) ─────────
positions = [
    # (left, bottom, width, height)
    (0.06, 0.77, 0.20, 0.18),   # A-left:  orig, donor
    (0.29, 0.77, 0.20, 0.18),   # A-right: new, donor
    (0.56, 0.77, 0.20, 0.18),   # B-left:  orig, status
    (0.79, 0.77, 0.20, 0.18),   # B-right: new, status
]
pca_axes = [fig.add_axes(p) for p in positions]

# A: by donor
for ax, pca_data, title in [(pca_axes[0], pca_orig, 'Original (Age+Gender)'),
                              (pca_axes[1], pca_new, 'No covariates')]:
    for donor in donors_unique:
        mask = meta['imageid'] == donor
        status = meta.loc[mask, 'donor_status'].iloc[0]
        ax.scatter(pca_data[mask, 0], pca_data[mask, 1],
                   c=[donor_color_map[donor]], s=6, alpha=0.5, rasterized=True)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=4)
    ax.set_xlabel('PC1', fontsize=8)
    ax.set_ylabel('PC2', fontsize=8)
    ax.tick_params(labelsize=7)

# Donor legend (compact)
handles = [plt.Line2D([0],[0], marker='o', color='w', markersize=5,
                       markerfacecolor=donor_color_map[d],
                       label=f'{d} ({meta.loc[meta["imageid"]==d, "donor_status"].iloc[0]})')
           for d in donors_unique]
pca_axes[1].legend(handles=handles, fontsize=5, ncol=2, loc='lower right',
                   framealpha=0.9, handletextpad=0.2, columnspacing=0.4,
                   borderpad=0.3, markerscale=1)

fig.text(0.03, 0.96, 'A', **lbl_kw)
fig.text(0.06, 0.959, 'Latent space PCA — colored by donor',
         fontsize=11, fontstyle='italic', transform=fig.transFigure)

# B: by disease status
for ax, pca_data, title in [(pca_axes[2], pca_orig, 'Original (Age+Gender)'),
                              (pca_axes[3], pca_new, 'No covariates')]:
    for status in STATUS_ORDER:
        mask = meta['donor_status'] == status
        ax.scatter(pca_data[mask, 0], pca_data[mask, 1],
                   c=STATUS_COLORS[status], s=6, alpha=0.5,
                   label=status, rasterized=True)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=4)
    ax.set_xlabel('PC1', fontsize=8)
    ax.set_ylabel('PC2', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=8, loc='upper right', framealpha=0.9, markerscale=2)

fig.text(0.53, 0.96, 'B', **lbl_kw)
fig.text(0.56, 0.959, 'Latent space PCA — colored by disease status',
         fontsize=11, fontstyle='italic', transform=fig.transFigure)

# ── Panel C: Per-dimension eta² ───────────────────────────
ax_c = axes[1, 0]
dims = np.arange(orig_latent.shape[1])
e2_od = per_dim_eta2(orig_latent, 'imageid')
e2_nd = per_dim_eta2(new_latent, 'imageid')
e2_os = per_dim_eta2(orig_latent, 'donor_status')
e2_ns = per_dim_eta2(new_latent, 'donor_status')

w = 0.2
ax_c.bar(dims - 1.5*w, e2_od, w, color=BLUE, alpha=0.85, label='Original — donor')
ax_c.bar(dims - 0.5*w, e2_os, w, color=BLUE, alpha=0.35, label='Original — disease', edgecolor=BLUE, linewidth=0.5)
ax_c.bar(dims + 0.5*w, e2_nd, w, color=RED, alpha=0.85, label='No cov — donor')
ax_c.bar(dims + 1.5*w, e2_ns, w, color=RED, alpha=0.35, label='No cov — disease', edgecolor=RED, linewidth=0.5)

ax_c.set_xlabel('Latent dimension', fontsize=11)
ax_c.set_ylabel(r'$\eta^2$ (variance explained)', fontsize=11)
ax_c.set_xticks(dims)
ax_c.legend(fontsize=7.5, ncol=2, loc='upper right', framealpha=0.9)
ax_c.set_ylim(0, 1.05)
ax_c.axhline(0.5, color='grey', ls='--', alpha=0.3, lw=0.8)
ax_c.set_title(r'C   Per-dimension $\eta^2$: donor identity vs. disease status',
               fontsize=12, fontweight='bold', loc='left', pad=8)

# ── Panel D: Silhouette distributions ─────────────────────
ax_d = axes[1, 1]
bins = np.linspace(-0.6, 0.8, 50)
ax_d.hist(sil_orig, bins=bins, alpha=0.55, color=BLUE, density=True,
          edgecolor='white', linewidth=0.3, label='Original (Age+Gender)')
ax_d.hist(sil_new, bins=bins, alpha=0.55, color=RED, density=True,
          edgecolor='white', linewidth=0.3, label='No covariates')

ax_d.axvline(np.mean(sil_orig), color=BLUE, ls='--', lw=2,
             label=f'Original mean = {np.mean(sil_orig):.3f}')
ax_d.axvline(np.mean(sil_new), color=RED, ls='--', lw=2,
             label=f'No cov mean = {np.mean(sil_new):.3f}')
ax_d.axvline(0, color='black', lw=0.5, alpha=0.4)

ax_d.set_xlabel('Donor silhouette score (cosine)', fontsize=11)
ax_d.set_ylabel('Density', fontsize=11)
ax_d.legend(fontsize=8.5, loc='upper right')
ax_d.set_title('D   Per-islet donor silhouette (lower = better mixing)',
               fontsize=12, fontweight='bold', loc='left', pad=8)

# ── Panel E: k-NN violin plots ───────────────────────────
ax_e = axes[2, 0]

positions_viol = [0.7, 1.3, 2.5, 3.1]
data_viol = [knn_orig_donor, knn_new_donor, knn_orig_status, knn_new_status]
colors_viol = [BLUE, RED, BLUE, RED]

for pos, data, col in zip(positions_viol, data_viol, colors_viol):
    parts = ax_e.violinplot([data], positions=[pos], showmeans=True,
                             showextrema=False, widths=0.45)
    for pc in parts['bodies']:
        pc.set_facecolor(col)
        pc.set_alpha(0.6)
    parts['cmeans'].set_color(col)
    parts['cmeans'].set_linewidth(2)

# Random expectations
exp_d = 1.0 / len(donors_unique)
exp_s = 1.0 / 3
ax_e.plot([0.4, 1.6], [exp_d, exp_d], 'k:', alpha=0.5, lw=1)
ax_e.plot([2.2, 3.4], [exp_s, exp_s], 'k:', alpha=0.5, lw=1)
ax_e.text(1.65, exp_d, 'random', fontsize=7, color='grey', va='center')
ax_e.text(3.45, exp_s, 'random', fontsize=7, color='grey', va='center')

# Add mean annotations
for pos, data, col in zip(positions_viol, data_viol, colors_viol):
    ax_e.text(pos, np.mean(data) + 0.03, f'{np.mean(data):.2f}',
              ha='center', va='bottom', fontsize=8, color=col, fontweight='bold')

ax_e.set_xticks([1.0, 2.8])
ax_e.set_xticklabels(['Same-donor\nneighbors', 'Same-status\nneighbors'], fontsize=10)
ax_e.set_ylabel('Fraction of k=15 nearest neighbors', fontsize=10)
ax_e.set_ylim(0, 1.08)
ax_e.legend([Patch(facecolor=BLUE, alpha=0.6), Patch(facecolor=RED, alpha=0.6)],
            ['Original (Age+Gender)', 'No covariates'],
            fontsize=9, loc='upper left')
ax_e.set_title('E   k-NN neighbor composition',
               fontsize=12, fontweight='bold', loc='left', pad=8)

# ── Panel F: Summary bar chart ────────────────────────────
ax_f = axes[2, 1]

metric_data = [
    (r'$R^2$ donor'      + '\n(lower=better)',    'r2_donor',           'lower'),
    (r'$R^2$ disease'    + '\n(higher=better)',   'r2_status',          'higher'),
    (r'$\eta^2$ excess'  + '\n(lower=better)',    'eta2_excess',        'lower'),
    ('Donor sil.\n(lower=better)',                 'sil_donor',          'lower'),
    ('k-NN donor\nenrich. /15\n(lower=better)',   'knn_donor_enrichment','lower'),
]

names  = [m[0] for m in metric_data]
o_vals = [results['original'][m[1]] / (15 if 'knn_donor' in m[1] else 1) for m in metric_data]
n_vals = [results['new_no_covariates'][m[1]] / (15 if 'knn_donor' in m[1] else 1) for m in metric_data]
dirs   = [m[2] for m in metric_data]

x = np.arange(len(names))
bw = 0.35
bars1 = ax_f.bar(x - bw/2, o_vals, bw, color=BLUE, alpha=0.8, label='Original (Age+Gender)')
bars2 = ax_f.bar(x + bw/2, n_vals, bw, color=RED,  alpha=0.8, label='No covariates')

# Star the winner
for i in range(len(names)):
    if dirs[i] == 'higher':
        win_orig = o_vals[i] >= n_vals[i]
    else:
        win_orig = o_vals[i] <= n_vals[i]

    star_x = x[i] - bw/2 if win_orig else x[i] + bw/2
    star_y = max(o_vals[i], n_vals[i]) + 0.02
    star_c = BLUE if win_orig else RED
    ax_f.text(star_x, star_y, '*', fontsize=18, fontweight='bold',
              color=star_c, ha='center', va='bottom')

ax_f.set_xticks(x)
ax_f.set_xticklabels(names, fontsize=8, ha='center')
ax_f.set_ylabel('Score', fontsize=11)
ax_f.legend(fontsize=9, loc='upper right')
ax_f.axhline(0, color='black', lw=0.5)
ax_f.text(0.02, 0.97, '* = better model for this metric',
          transform=ax_f.transAxes, fontsize=7.5, fontstyle='italic', color='grey', va='top')
ax_f.set_title('F   Summary metrics comparison',
               fontsize=12, fontweight='bold', loc='left', pad=8)

# ── Row 3: PCA colored by pseudotime (bonus panel) ───────
ax_g = axes[3, 0]
ax_h = axes[3, 1]

pt = meta['dpt_pseudotime'].values
valid = ~np.isnan(pt)

sc_g = ax_g.scatter(pca_orig[valid, 0], pca_orig[valid, 1],
                    c=pt[valid], cmap='viridis', s=8, alpha=0.6, rasterized=True)
ax_g.set_xlabel('PC1', fontsize=10)
ax_g.set_ylabel('PC2', fontsize=10)
ax_g.set_title('G   Original — PCA colored by pseudotime',
               fontsize=12, fontweight='bold', loc='left', pad=8)
cb_g = plt.colorbar(sc_g, ax=ax_g, shrink=0.8, pad=0.02)
cb_g.set_label('Pseudotime', fontsize=9)

sc_h = ax_h.scatter(pca_new[valid, 0], pca_new[valid, 1],
                    c=pt[valid], cmap='viridis', s=8, alpha=0.6, rasterized=True)
ax_h.set_xlabel('PC1', fontsize=10)
ax_h.set_ylabel('PC2', fontsize=10)
ax_h.set_title('H   No covariates — PCA colored by pseudotime',
               fontsize=12, fontweight='bold', loc='left', pad=8)
cb_h = plt.colorbar(sc_h, ax=ax_h, shrink=0.8, pad=0.02)
cb_h.set_label('Pseudotime', fontsize=9)

# ── Save ───────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    out = os.path.join(RESULTS_DIR, f"supplemental_scvi_comparison.{ext}")
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out}")

plt.close()
print("Done.")
