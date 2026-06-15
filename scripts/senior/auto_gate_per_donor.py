#!/usr/bin/env python3
"""
Phase 0 — Automated PER-DONOR gating for Islet Explorer Senior (run in the `scimap` conda env).

Replaces the manual `sm.pl.napariGater` step (59 sliders re-validated per donor) with a
scripted, per-(marker, donor) gate + a diagnostic contact-sheet for QC. Reuses the exact
scimap rules engine (`sm.tl.phenotype_cells`) and the existing rules CSV, so the downstream
`obs['phenotype']` contract is unchanged.

WHY PER-DONOR: donor/batch staining-intensity variation is the verified root cause of the
prior global-auto-thresholding failure (INS fraction-positive varies ~8x across donors). A
single global gate is wrong per donor; scimap's `rescale` already supports a markers x imageid
gate table and scales each donor independently, so per-donor gates fix the dominant failure
mode while keeping everything else identical.

Pipeline:
  raw MFI  -> per-(marker,donor) GMM gate (log1p, posterior crossover; unimodal fallback)
           -> rescale per donor (gate -> 0.5, in [0,1])  -> sm.tl.phenotype_cells(rules)
           -> obs['phenotype'] (+ obs['gate_sep_min'] coarse confidence)

Rescale backends:
  --rescale fast   (default) vectorised numpy piecewise-linear rescale, faithful to scimap's
                   per-image MinMaxScaler in the positive range (gate->0.5, max->1). ~100x
                   faster than scimap's per-column pandas sort; required at 23M scale.
  --rescale scimap call sm.pp.rescale (per-image pandas path). Use to validate `fast`.
  --validate-rescale  run BOTH on the loaded cells and report per-cell phenotype concordance.

Smoke (one donor, capped):
  python scripts/senior/auto_gate_per_donor.py --donors 6374 --max-cells-per-donor 200000 \
         --out /tmp/pheno_6374_smoke.h5ad --diagnostics-pdf /tmp/gates_6374.pdf --validate-rescale
Full:
  python scripts/senior/auto_gate_per_donor.py --out ~/IO60panc2nd/senior_data/singlecell_protein_phenotyped.h5ad \
         --gates-out ~/IO60panc2nd/senior_data/per_donor_gates.csv \
         --diagnostics-pdf ~/IO60panc2nd/senior_data/per_donor_gate_diagnostics.pdf
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm

BASE = Path.home() / "IO60panc2nd" / "senior_data"
INP = BASE / "singlecell_protein.h5ad"
RULES = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"


# ------------------------------------------------------------------ gate seeding
def gmm_gate_log_with_sep(x, min_sep=2.0):
    """2-component GMM gate on log1p(raw), at the posterior crossover; fall back to the
    99th percentile (log1p) when the marker is effectively unimodal. Mirrors
    `gmm_gate_log` in phenotype_scimap.py but also returns a separation score
    (component-mean gap in pooled-SD units) used as a coarse per-gate confidence.

    Returns (gate_log1p, sep_score). sep is np.inf for the unimodal-fallback (treated
    as 'low confidence' downstream by being flagged separately)."""
    from sklearn.mixture import GaussianMixture
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 50:
        return float("nan"), 0.0
    lx = np.log1p(np.clip(x, 0, None))
    sub = lx[np.linspace(0, lx.size - 1, min(50000, lx.size)).astype(int)].reshape(-1, 1)
    try:
        g = GaussianMixture(2, n_init=1, random_state=0, reg_covar=1e-3).fit(sub)
    except Exception:
        return float(np.quantile(lx, 0.99)), 0.0
    mu = g.means_.ravel()
    sd = np.sqrt(g.covariances_.ravel())
    hi = int(np.argmax(mu)); lo = 1 - hi
    pooled = max(0.5 * (sd[hi] + sd[lo]), 1e-6)
    sep = float((mu[hi] - mu[lo]) / pooled)
    if (mu[hi] - mu[lo]) < min_sep * pooled:
        # unimodal-ish: high percentile fallback, mark as weakly separated
        return float(np.quantile(lx, 0.99)), sep
    grid = np.linspace(mu[lo], mu[hi], 256).reshape(-1, 1)
    p = g.predict_proba(grid)[:, hi]
    cr = np.where(p >= 0.5)[0]
    gate = float(grid[cr[0], 0] if cr.size else mu[hi])
    return gate, sep


def apply_gate_overrides(gates_df, overrides_path):
    """Apply manual gate corrections after auto-seeding. CSV columns:
        marker, gate_log1p[, donor]
    A row without `donor` (or donor='all') overrides that marker for ALL donors.
    This is the Phase-0 QC loop: scan the diagnostic PDF, find the few markers whose
    auto-gate sits in the bulk (e.g. SOX2 — a nuclear TF with high CODEX background that
    univariate GMM over-calls at ~90% positive), and raise just those — not 59 sliders."""
    ov = pd.read_csv(overrides_path)
    ov.columns = [c.strip() for c in ov.columns]
    n = 0
    for _, r in ov.iterrows():
        m = str(r["marker"]).strip()
        g = float(r["gate_log1p"])
        d = str(r["donor"]).strip() if "donor" in ov.columns and not pd.isna(r.get("donor")) else "all"
        if m not in gates_df.index:
            print(f"[override] WARN marker '{m}' not in gates; skipped"); continue
        cols = list(gates_df.columns) if d in ("all", "") else [c for c in gates_df.columns if str(c) == d]
        for c in cols:
            gates_df.loc[m, c] = g; n += 1
    print(f"[override] applied {n} gate corrections from {overrides_path}")
    return gates_df


def warn_overgated(adata_rescaled, rules_path, thresh=0.5):
    """Flag rule markers whose positive fraction (>0.5 after rescale) exceeds `thresh`.
    Lineage markers (EpCAM, Ker8_18) can legitimately be high; subtype-defining markers
    (SOX2, TP63, FOXP3 ...) that are high are likely over-gated — review them in the PDF."""
    rules = pd.read_csv(rules_path)
    rule_markers = [c for c in rules.columns[2:]]
    X = adata_rescaled.X if isinstance(adata_rescaled.X, np.ndarray) else adata_rescaled.X.toarray()
    pf = pd.Series((X > 0.5).mean(axis=0), index=list(adata_rescaled.var_names))
    flagged = pf[rule_markers][pf[rule_markers] > thresh].sort_values(ascending=False)
    if len(flagged):
        print(f"\n[warn] rule markers positive in >{thresh:.0%} of cells (review in diagnostics PDF; "
              f"override the subtype-defining ones):")
        for m, v in flagged.items():
            print(f"        {m:18s} {v:.1%} positive")
    return flagged


def compute_per_donor_gates(adata, donor_col, min_sep):
    """Return (gates_df, sep_df): markers x donors tables of log1p gate values and
    separation scores. Gates feed scimap rescale (markers index + one column per imageid)."""
    markers = list(adata.var_names)
    donors = list(pd.unique(adata.obs[donor_col]))
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    didx = adata.obs[donor_col].values
    gates = pd.DataFrame(index=pd.Index(markers, name="markers"), columns=donors, dtype=float)
    seps = pd.DataFrame(index=pd.Index(markers, name="markers"), columns=donors, dtype=float)
    for d in donors:
        rows = np.where(didx == d)[0]
        Xd = X[rows]
        print(f"[gates] donor {d}: {rows.size:,} cells")
        for j, m in enumerate(markers):
            gv, sep = gmm_gate_log_with_sep(Xd[:, j], min_sep)
            gates.loc[m, d] = gv
            seps.loc[m, d] = sep
    return gates, seps


# ------------------------------------------------------------------ rescaling
def rescale_fast(adata, gates_df, donor_col, log=True):
    """Vectorised per-donor piecewise-linear rescale faithful to scimap's per-image
    MinMaxScaler: within each donor, for each marker with gate g (in the log1p scale),
        x >= g : 0.5 + 0.5 * (x - g) / (max - g)        [== scimap high MinMaxScaler]
        x <  g : 0.5 * (x - min) / (g - min)            [== scimap low  MinMaxScaler]
    Degenerate ranges collapse to 0.5/0.0 as scimap does. Returns an (n, M) float32 array
    aligned to adata row order. Positivity (>=0.5) and positive-range magnitudes match
    scimap exactly, so sm.tl.phenotype_cells assignments are preserved."""
    markers = list(adata.var_names)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    X = np.asarray(X, dtype=np.float64)
    if log:
        X = np.log1p(np.clip(X, 0, None))
    out = np.empty_like(X, dtype=np.float32)
    didx = adata.obs[donor_col].values
    for d in pd.unique(didx):
        rows = np.where(didx == d)[0]
        Xd = X[rows]                                   # (n_d, M)
        g = gates_df[d].reindex(markers).to_numpy(dtype=float)   # (M,)
        xmin = Xd.min(axis=0)
        xmax = Xd.max(axis=0)
        hi_den = np.where(xmax > g, xmax - g, np.inf)  # inf -> high frac 0
        lo_den = np.where(g > xmin, g - xmin, np.inf)  # inf -> low frac 0
        pos = Xd >= g                                  # broadcast (n_d, M)
        hi = 0.5 + 0.5 * (Xd - g) / hi_den
        lo = 0.5 * (Xd - xmin) / lo_den
        scaled = np.where(pos, hi, lo)
        scaled = np.clip(scaled, 0.0, 1.0)
        out[rows] = scaled.astype(np.float32)
    return out


def rescale_scimap(adata, gates_df):
    """scimap's native per-image rescale (validation / fallback). Mutates a copy."""
    a = adata.copy()
    if a.raw is None:
        a.raw = a
    gdf = gates_df.copy()
    gdf.index.name = "markers"
    a = sm.pp.rescale(a, gate=gdf, log=True, method="all", verbose=False)
    return np.asarray(a.X, dtype=np.float32)


# ------------------------------------------------------------------ phenotyping
def phenotype(adata_rescaled, rules_path):
    rules = pd.read_csv(rules_path)
    a = sm.tl.phenotype_cells(adata_rescaled, phenotype=rules, label="phenotype", verbose=False)
    return a


# ------------------------------------------------------------------ diagnostics
def write_diagnostics_pdf(adata, gates_df, seps_df, donor_col, pdf_path, max_donors=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    markers = list(adata.var_names)
    donors = list(gates_df.columns)[: (max_donors or len(gates_df.columns))]
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    didx = adata.obs[donor_col].values
    ncol = min(5, len(donors)); nrow = int(np.ceil(len(donors) / ncol))
    Path(pdf_path).parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        for j, m in enumerate(markers):
            fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.2 * nrow), squeeze=False)
            for k, d in enumerate(donors):
                ax = axes[k // ncol][k % ncol]
                rows = np.where(didx == d)[0]
                lx = np.log1p(np.clip(X[rows, j], 0, None))
                ax.hist(lx, bins=80, color="#4477AA", alpha=0.8)
                g = gates_df.loc[m, d]; sep = seps_df.loc[m, d]
                ax.axvline(g, color="#CC3311", lw=1.5)
                ax.set_title(f"{d}  sep={sep:.1f}", fontsize=8)
                ax.set_yticks([])
            for k in range(len(donors), nrow * ncol):
                axes[k // ncol][k % ncol].axis("off")
            fig.suptitle(f"{m} — log1p(raw) with per-donor gate", fontsize=11)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig); plt.close(fig)
    print(f"[diag] wrote {pdf_path} ({len(markers)} pages)")


# ------------------------------------------------------------------ IO
def load_subset(path, donors, max_cells_per_donor, donor_col):
    """Read the h5ad fully into memory then subset with numpy. anndata backed-mode
    fancy-indexing of scattered rows is pathologically slow on this 23M-row / 13GB file;
    a single contiguous in-memory read (X ~5.4GB float32) + numpy masking is far faster
    and the box has the RAM (scVI training loads the same array)."""
    print(f"[load] reading {path} into memory ...")
    a = ad.read_h5ad(path)
    print(f"[load] full {a.shape}")
    if donors is None and max_cells_per_donor is None:
        return a
    dvals = a.obs[donor_col].astype(str).to_numpy()
    sel = np.ones(dvals.size, dtype=bool) if donors is None else np.isin(dvals, [str(x) for x in donors])
    idx = np.where(sel)[0]
    if max_cells_per_donor is not None:
        keep = []
        for d in pd.unique(dvals[idx]):
            di = idx[dvals[idx] == d]
            if di.size > max_cells_per_donor:
                di = di[np.linspace(0, di.size - 1, max_cells_per_donor).astype(int)]
            keep.append(di)
        idx = np.sort(np.concatenate(keep))
    a = a[idx].to_memory() if a.isbacked else a[idx].copy()
    print(f"[load] subset {a.shape} | donors {pd.unique(a.obs[donor_col]).tolist()}")
    return a


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=INP)
    ap.add_argument("--rules", type=Path, default=RULES)
    ap.add_argument("--donor-col", default="imageid")
    ap.add_argument("--donors", nargs="*", default=None, help="subset of donor ids (default all)")
    ap.add_argument("--max-cells-per-donor", type=int, default=None, help="downsample per donor (smoke)")
    ap.add_argument("--min-sep", type=float, default=2.0)
    ap.add_argument("--gate-overrides", type=Path, default=None,
                    help="CSV (marker,gate_log1p[,donor]) of manual gate corrections applied after auto-seeding")
    ap.add_argument("--rescale", choices=["fast", "scimap"], default="fast")
    ap.add_argument("--validate-rescale", action="store_true", help="run both backends, report concordance")
    ap.add_argument("--out", type=Path, default=None, help="phenotyped h5ad output")
    ap.add_argument("--gates-out", type=Path, default=None, help="markers x donors gate CSV")
    ap.add_argument("--diagnostics-pdf", type=Path, default=None)
    a = ap.parse_args(argv)

    adata = load_subset(a.input, a.donors, a.max_cells_per_donor, a.donor_col)

    # validate rules markers exist
    rules = pd.read_csv(a.rules)
    rule_markers = [c for c in rules.columns[2:]]
    missing = [m for m in rule_markers if m not in set(adata.var_names)]
    if missing:
        raise SystemExit(f"[err] rule markers missing from var_names: {missing}\n var_names={list(adata.var_names)}")

    gates_df, seps_df = compute_per_donor_gates(adata, a.donor_col, a.min_sep)
    if a.gate_overrides:
        gates_df = apply_gate_overrides(gates_df, a.gate_overrides)
    if a.gates_out:
        a.gates_out.parent.mkdir(parents=True, exist_ok=True)
        gates_df.to_csv(a.gates_out)
        seps_df.to_csv(str(a.gates_out).replace(".csv", "_sepscore.csv"))
        print(f"[gates] wrote {a.gates_out}")

    if a.diagnostics_pdf:
        write_diagnostics_pdf(adata, gates_df, seps_df, a.donor_col, a.diagnostics_pdf)

    # --- rescale + phenotype
    def build_and_phenotype(Xr, tag):
        ar = ad.AnnData(X=Xr, obs=adata.obs.copy(), var=adata.var.copy())
        ar = phenotype(ar, a.rules)
        vc = ar.obs["phenotype"].value_counts()
        print(f"\n[{tag}] phenotype counts:\n{vc.to_string()}")
        return ar

    if a.validate_rescale:
        Xf = rescale_fast(adata, gates_df, a.donor_col)
        Xs = rescale_scimap(adata, gates_df)
        af = build_and_phenotype(Xf, "fast")
        as_ = build_and_phenotype(Xs, "scimap")
        agree = float((af.obs["phenotype"].values == as_.obs["phenotype"].values).mean())
        maxabs = float(np.abs(Xf.astype(np.float64) - Xs.astype(np.float64)).max())
        print(f"\n[validate] per-cell phenotype concordance fast-vs-scimap = {agree:.4%} | max|ΔX| = {maxabs:.4f}")
        chosen = af if a.rescale == "fast" else as_
    else:
        Xr = rescale_fast(adata, gates_df, a.donor_col) if a.rescale == "fast" else rescale_scimap(adata, gates_df)
        chosen = build_and_phenotype(Xr, a.rescale)

    warn_overgated(chosen, a.rules)

    # coarse per-cell confidence: the worst (min) separation among that donor's gates
    donor_min_sep = seps_df.min(axis=0)  # per donor
    chosen.obs["gate_sep_min"] = chosen.obs[a.donor_col].map(donor_min_sep).astype(float).values

    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        chosen.write_h5ad(a.out)
        print(f"[done] wrote {a.out} {chosen.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
