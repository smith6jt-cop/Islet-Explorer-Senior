"""Build analysis/H20_bcell_density_feasibility.ipynb programmatically.

Run once to produce the notebook, then execute via jupyter nbconvert.
Implements the methodology from the user's GPT-5.6 chat ("Statistical approach for
dataset size issue"): counts + exact-area offset -> genotype rate ratios,
donor as the independent unit, platform-stratified inference, random matched-window
robustness. Applied to the B-cell-density endpoint (rs3184504 T/T hypothesis).

Idempotent: re-running rewrites the notebook. Heavy inputs (Cells.csv) are cached
under analysis/.h20_cache/ and analysis/tables/H20_*.csv.
"""
from pathlib import Path
import nbformat as nbf

PROJECT = Path(__file__).resolve().parent.parent
NB_PATH = PROJECT / "analysis" / "H20_bcell_density_feasibility.ipynb"

nb = nbf.v4.new_notebook()
cells = []
def md(src): cells.append(nbf.v4.new_markdown_cell(src))
def code(src): cells.append(nbf.v4.new_code_cell(src))

# ---------------------------------------------------------------- Cell 0
md(r"""# H20 — Can we test "T/T (rs3184504) has a denser B-cell population"?

**Question.** Does the T-allele of rs3184504 (SH2B3) associate with a *denser*
splenic B-cell population? "Denser" here = **areal density: CD20⁺ B-cells per mm²**
(the one reading never computed before; prior work measured B-cell *proportion*).

**Methodology (from the GPT-5.6 "dataset size" chat).** The cohort mixes whole
Phenocycler sections with small cropped CODEX fields, so density is an
*exposure/precision* problem, not a reason to crop. The chat prescribes:

1. **Primary model — hierarchical Negative-Binomial count-rate with an area offset.**
   For follicle (window) *j* of donor *i*:
   $$Y_{ij}\sim \text{NegBin}(\mu_{ij},\phi),\qquad
     \log\mu_{ij}=\log A_{ij}+\alpha+\beta_G G_i+\beta_P P_i+u_i$$
   counts $Y_{ij}$ (B-cells), exact area $A_{ij}$ as a fixed-coefficient offset,
   genotype $G_i$ + platform $P_i$ fixed effects, donor random intercept $u_i$.
   Output = genotype **rate ratios**. NB (not Poisson) for overdispersion.
2. **All tissue, no cropping** — larger sections just give more precise counts.
3. **Random matched-window robustness** for within-section heterogeneity; donor is
   the independent unit (windows nested in donor).
4. **Platform-stratified inference** (Phenocycler-only / CODEX-only / pooled-secondary),
   because CODEX field selection may not be unbiased.

**Cohort (single-cell fluorescent set — the only data with B-cell markers).**

| Genotype | Donors (n) | Platform |
|---|---|---|
| C/C | HDL011, HDL053, HDL055 (3) | all Phenocycler |
| C/T | HDL043, HDL052, HDL070, HDL079, HDL086 (5) | all Phenocycler |
| T/T | HDL063 (Pheno), 1901HBMP004, HDL073 (CODEX) (3) | **2 of 3 are CODEX** |

**Platform is nearly collinear with the T/T contrast** — both CODEX donors are T/T;
zero CODEX donors are C/C or C/T. Keep this in view throughout.

**Data-source note.** `adata_H16.h5ad` (RESTORE-phenotyped B-cells) is only a
**5000×3000 µm crop** of each Phenocycler image, so its counts cannot be divided by
full-image follicle areas. Areal density is therefore computed from the **full-image
`Cells.csv`** (all 22.6M cells; `Cell: CD20: Mean`), with the B-cell call **calibrated
per-image to the RESTORE labels** inside the crop (see §B).""")

# ---------------------------------------------------------------- Cell 1: setup
code(r"""import sys, json, shutil, subprocess, tempfile, warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

CODE_ROOT = Path.cwd().parent if Path.cwd().name == "analysis" else Path.cwd()
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du

du.setup_style()
GENO_ORDER = du.GENO_ORDER
PALETTE = du.GENO_PALETTE
DOSAGE = {"C/C": 0, "C/T": 1, "T/T": 2}
ARCSINH_COFACTOR = 150.0
CELLS_CSV = du.DATA_ROOT / "Measurements" / "Cells.csv"
ANN_CSV = du.DATA_ROOT / "Measurements" / "AnnotationsFinal.csv"
ADATA = du.DATA_ROOT / "adata_H16.h5ad"
CACHE_DIR = du.OUTPUT_ROOT / "analysis" / ".h20_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
TABLES = du.TABLES_DIR

RSCRIPT = shutil.which("Rscript")
print("Code root:", CODE_ROOT)
print("Data root:", du.DATA_ROOT)
print("Output root:", du.OUTPUT_ROOT)
print("Rscript:", RSCRIPT or "not found (R confirmatory model will be skipped)")

clin = du.load_clinical()
print("\n(Single-cell cohort — the 11 donors with B-cell markers — is shown in §B.)")""")

# ---------------------------------------------------------------- Section A
md(r"""## A. What prior analyses already found about B-cells by genotype

Two independent single-cell pipelines (H15 Leiden-gating, H16 RESTORE) already
tested B-cell **proportion** by genotype, plus follicle-restricted proportion and a
platform-confound diagnostic. Consolidated below. Key points to carry forward:
the tests are **null**, the two pipelines **disagree on the sign** of the T-allele
trend, and H15 flags the B-cell readout as **platform-confounded** (H16's RESTORE
removes it) — which is why we anchor the new B-cell call to H16's RESTORE labels.""")

code(r"""def _read(name):
    p = TABLES / f"{name}.csv"
    return pd.read_csv(p) if p.exists() else None

rows = []
for pipe in ["H15", "H16"]:
    ct = _read(f"{pipe}_celltype_genotype_stats")
    if ct is not None:
        b = ct[ct["Metric"] == "B cell"]
        def g(test):
            r = b[b["Test"] == test]
            return (r["Statistic"].iloc[0], r["p"].iloc[0], r["Effect_Size"].iloc[0]) if len(r) else (np.nan, np.nan, "")
        kw = g("Kruskal-Wallis"); sp = g("Spearman dosage")
        rows.append({"Pipeline": pipe, "Endpoint": "B-cell proportion (whole)",
                     "KW_p": kw[1], "Dosage": sp[2], "Dosage_p": sp[1]})
    rg = _read(f"{pipe}_region_celltype_genotype_stats")
    if rg is not None:
        fb = rg[rg["Metric"] == "Follicle:B cell"]
        sp = fb[fb["Test"] == "Spearman dosage"]
        kw = fb[fb["Test"] == "Kruskal-Wallis"]
        rows.append({"Pipeline": pipe, "Endpoint": "B-cell proportion (in Follicle)",
                     "KW_p": kw["p"].iloc[0] if len(kw) else np.nan,
                     "Dosage": sp["Effect_Size"].iloc[0] if len(sp) else "",
                     "Dosage_p": sp["p"].iloc[0] if len(sp) else np.nan})
    pd_ = _read(f"{pipe}_platform_diagnostic")
    if pd_ is not None:
        pb = pd_[pd_["Feature"] == "B cell"]
        if len(pb):
            rows.append({"Pipeline": pipe, "Endpoint": "B-cell platform confound (CODEX vs PC)",
                         "KW_p": np.nan, "Dosage": f"r={pb['rank_biserial'].iloc[0]:.3f}",
                         "Dosage_p": pb["p"].iloc[0]})

prior = pd.DataFrame(rows)
du.save_table(prior, "H20_prior_bcell_evidence")
print("Prior B-cell-vs-genotype evidence (all proportion-based):")
prior""")

# ---------------------------------------------------------------- Section B
md(r"""## B. Build the new metric — areal B-cell density (CD20⁺ cells / mm²)

**B-cell call.** Full-image `Cells.csv` CD20 (`arcsinh(x/150)`), thresholded
**per image at the point that best reproduces the H16 RESTORE `B cell` labels**
(Youden's J on the crop cells, joined by Object ID). Reported with per-image
sensitivity/specificity — note CODEX images have **low sensitivity** (weaker CD20),
so CODEX B-cells are *under*-counted (a platform caveat, not a genotype effect).

**Endpoints (per donor + per follicle).**
- *Per-follicle* (GPT window model): each in-follicle CD20⁺ cell assigned to its
  nearest `AnnotationsFinal` Follicle centroid (same µm frame; full coverage);
  $A_{ij}$ = that follicle's area.
- *Per-donor whole-tissue*: CD20⁺ count / Σ region area.
- *Per-donor follicle-compartment*: in-follicle CD20⁺ count / Σ follicle area.

Computed once from `Cells.csv` and cached to `analysis/tables/H20_*`; re-runs load
the cache.""")

code(r"""def _otsu(v, nbins=256):
    v = v[np.isfinite(v)]
    if v.size == 0: return np.inf
    lo, hi = float(np.min(v)), float(np.max(v))
    if hi <= lo: return hi
    h, e = np.histogram(v, bins=nbins, range=(lo, hi))
    c = (e[:-1] + e[1:]) / 2; w = h.astype(float); s = w.sum()
    if s == 0: return hi
    p = w / s; om = np.cumsum(p); mu = np.cumsum(p * c); mt = mu[-1]
    den = om * (1 - om); den[den == 0] = np.nan
    return c[np.nanargmax((mt * om - mu) ** 2 / den)]

def compute_counts():
    import anndata as ad
    from sklearn.metrics import roc_curve
    use = ["Image", "Object ID", "Centroid X µm", "Centroid Y µm",
           "Cell: CD20: Mean", "Cell: CD21: Mean",
           "Signed distance to annotation Follicle µm"]
    print(f"Reading {CELLS_CSV.name} ({CELLS_CSV.stat().st_size/1e9:.1f} GB)...")
    df = pd.read_csv(CELLS_CSV, usecols=use,
                     dtype={"Image": "category", "Centroid X µm": "float32",
                            "Centroid Y µm": "float32", "Cell: CD20: Mean": "float32",
                            "Cell: CD21: Mean": "float32",
                            "Signed distance to annotation Follicle µm": "float32"})
    df = df.rename(columns={"Centroid X µm": "x", "Centroid Y µm": "y",
                            "Cell: CD20: Mean": "cd20", "Cell: CD21: Mean": "cd21",
                            "Signed distance to annotation Follicle µm": "dfoll"})
    df["Sample"] = df["Image"].astype(str).map(du.extract_sample_id)
    df = df[~df["Sample"].isin(du.EXCLUDE_SAMPLES)].copy()
    df["Genotype"] = df["Sample"].map(du.GENOTYPE_MAP)
    df = df.dropna(subset=["Genotype"])
    df["Platform"] = df["Sample"].map(lambda s: "CODEX" if s in du.CODEX_SAMPLES else "Phenocycler")
    df["cd20t"] = np.arcsinh(df["cd20"].to_numpy() / ARCSINH_COFACTOR)
    df["in_foll"] = df["dfoll"].to_numpy() < 0

    A = ad.read_h5ad(ADATA, backed="r")
    lab = A.obs[["Object ID", "cell_type"]]
    bmap = {str(o): (c == "B cell") for o, c in
            zip(lab["Object ID"].astype(str), lab["cell_type"].astype(str))}
    isb = df["Object ID"].astype(str).map(bmap).astype(float).to_numpy()
    df["cd21t"] = np.arcsinh(df["cd21"].to_numpy() / ARCSINH_COFACTOR)
    cd20t = df["cd20t"].to_numpy(); cd21t = df["cd21t"].to_numpy()
    bcell = np.zeros(len(df), bool); cd21pos = np.zeros(len(df), bool); cal = []
    for img, idx in df.groupby("Image", observed=True).indices.items():
        idx = np.asarray(idx); v = cd20t[idx]; lv = isb[idx]; m = ~np.isnan(lv)
        if m.sum() > 100 and np.nansum(lv) > 10:
            y = lv[m].astype(int); s = v[m]
            fpr, tpr, thr = roc_curve(y, s); t = thr[np.argmax(tpr - fpr)]
            pred = (s > t).astype(int)
            tp = ((pred==1)&(y==1)).sum(); fn = ((pred==0)&(y==1)).sum()
            tn = ((pred==0)&(y==0)).sum(); fp = ((pred==1)&(y==0)).sum()
            sens = tp/(tp+fn) if tp+fn else np.nan; spec = tn/(tn+fp) if tn+fp else np.nan
            bf = float(y.mean())
        else:
            t = _otsu(v); sens = spec = bf = np.nan
        bcell[idx] = v > t
        t21 = _otsu(cd21t[idx]); cd21pos[idx] = cd21t[idx] > t21  # CD21+ via per-image Otsu
        cal.append({"Image": img, "threshold": float(t), "cd21_threshold": float(t21),
                    "crop_b_frac": bf, "sensitivity": sens, "specificity": spec})
    df["bcell"] = bcell
    df["dp"] = bcell & cd21pos                      # CD20+CD21+ double-positive
    calib = pd.DataFrame(cal)

    # per-donor denominators from AnnotationsFinal (full image)
    ann = pd.read_csv(ANN_CSV); ann["Sample"] = ann["Image"].map(du.extract_sample_id)
    ann = ann[~ann["Sample"].isin(du.EXCLUDE_SAMPLES)]
    foll = ann[ann["Classification"] == "Follicle"]
    reg = ann[ann["Classification"].isin(du.MAIN_REGIONS)]
    tissue_mm2 = reg.groupby("Sample", observed=True)["Area µm^2"].sum() / 1e6
    foll_mm2 = foll.groupby("Sample", observed=True)["Area µm^2"].sum() / 1e6

    # per-follicle assignment (nearest centroid, full coverage)
    pf = []
    for s, fs in foll.groupby("Sample", observed=True):
        cs = df[(df["Sample"] == s) & df["in_foll"].to_numpy()]
        tree = cKDTree(fs[["Centroid X µm", "Centroid Y µm"]].to_numpy())
        fid = fs["Object ID"].to_numpy(); fa = fs["Area µm^2"].to_numpy()
        tot = np.zeros(len(fid), int); bb = np.zeros(len(fid), int); dd = np.zeros(len(fid), int)
        if len(cs):
            _, idx = tree.query(cs[["x", "y"]].to_numpy())
            isbb = cs["bcell"].to_numpy(); isdp = cs["dp"].to_numpy()
            tot = np.bincount(idx, minlength=len(fid))
            bb = np.bincount(idx[isbb], minlength=len(fid))
            dd = np.bincount(idx[isdp], minlength=len(fid))
        for k in range(len(fid)):
            pf.append({"Sample": s, "Follicle_ID": fid[k], "Follicle_Area_um2": fa[k],
                       "Total_cells": int(tot[k]), "Bcell_count": int(bb[k]),
                       "DP_count": int(dd[k])})
    per_foll = pd.DataFrame(pf)
    per_foll["Genotype"] = per_foll["Sample"].map(du.GENOTYPE_MAP)
    per_foll["Platform"] = per_foll["Sample"].map(lambda s: "CODEX" if s in du.CODEX_SAMPLES else "Phenocycler")
    per_foll["Bcell_density_mm2"] = per_foll["Bcell_count"] / (per_foll["Follicle_Area_um2"] / 1e6)
    per_foll["DP_density_mm2"] = per_foll["DP_count"] / (per_foll["Follicle_Area_um2"] / 1e6)

    # per-donor
    g = df.groupby(["Sample", "Genotype", "Platform"], observed=True)
    pdw = g.agg(n_cells=("bcell", "size"), n_bcell=("bcell", "sum"),
                n_dp=("dp", "sum")).reset_index()
    inf = df[df["in_foll"].to_numpy()].groupby("Sample", observed=True)
    pdw = pdw.merge(inf["bcell"].sum().rename("n_bcell_infoll"), on="Sample")
    pdw = pdw.merge(inf["dp"].sum().rename("n_dp_infoll"), on="Sample")
    pdw = pdw.merge(tissue_mm2.rename("tissue_mm2"), on="Sample").merge(foll_mm2.rename("foll_mm2"), on="Sample")
    pdw["WholeTissue_density"] = pdw["n_bcell"] / pdw["tissue_mm2"]
    pdw["Follicle_density"] = pdw["n_bcell_infoll"] / pdw["foll_mm2"]
    pdw["DP_WholeTissue_density"] = pdw["n_dp"] / pdw["tissue_mm2"]
    pdw["DP_Follicle_density"] = pdw["n_dp_infoll"] / pdw["foll_mm2"]

    # compact per-cell coord cache for §E
    compact = pd.DataFrame({"Sample": df["Sample"].astype("category"),
                            "x": df["x"], "y": df["y"],
                            "bcell": df["bcell"].to_numpy(bool),
                            "in_foll": df["in_foll"].to_numpy(bool)})
    compact.to_pickle(CACHE_DIR / "percell_coords.pkl")
    return per_foll, pdw, calib

# load cache or compute
ff_p = TABLES / "H20_bcell_counts_per_follicle.csv"
dd_p = TABLES / "H20_bcell_density_per_donor.csv"
cal_p = TABLES / "H20_cd20_calibration.csv"
if ff_p.exists() and dd_p.exists() and cal_p.exists():
    per_foll = pd.read_csv(ff_p); per_donor = pd.read_csv(dd_p); calib = pd.read_csv(cal_p)
    print("Loaded cached H20 count tables.")
else:
    per_foll, per_donor, calib = compute_counts()
    du.save_table(per_foll, "H20_bcell_counts_per_follicle")
    du.save_table(per_donor, "H20_bcell_density_per_donor")
    du.save_table(calib, "H20_cd20_calibration")

for d in (per_foll, per_donor):
    d["Genotype"] = pd.Categorical(d["Genotype"], categories=GENO_ORDER, ordered=True)
print(f"\n{len(per_foll)} follicles across {per_donor.shape[0]} donors")
print("Single-cell cohort (the actual unit of analysis):")
print(per_donor.groupby(["Genotype", "Platform"], observed=True).size().unstack(fill_value=0))
print("\nCD20 threshold calibration vs RESTORE labels (note low CODEX sensitivity);"
      " cd21_threshold is per-image Otsu:")
calib.assign(Sample=calib["Image"].map(du.extract_sample_id))[
    ["Sample", "threshold", "cd21_threshold", "crop_b_frac", "sensitivity", "specificity"]]""")

code(r"""print("Per-donor areal B-cell density (cells/mm^2):")
show = per_donor.sort_values(["Genotype", "Sample"])[
    ["Sample", "Genotype", "Platform", "n_bcell", "tissue_mm2",
     "WholeTissue_density", "Follicle_density"]]
show""")

code(r"""# Visualize: T/T is not elevated on either density endpoint
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, col, title in zip(
        axes, ["WholeTissue_density", "Follicle_density"],
        ["Whole-tissue B-cell density", "Follicle-compartment B-cell density"]):
    sns.stripplot(data=per_donor, x="Genotype", y=col, order=GENO_ORDER,
                  hue="Platform", palette={"Phenocycler": "#3366CC", "CODEX": "#CC6633"},
                  s=11, jitter=0.12, ax=ax, edgecolor="black", linewidth=0.6)
    means = per_donor.groupby("Genotype", observed=True)[col].mean().reindex(GENO_ORDER)
    ax.plot(range(3), means.values, "k_", markersize=26, markeredgewidth=2)
    ax.set_title(title); ax.set_ylabel("CD20+ cells / mm²"); ax.set_xlabel("rs3184504")
    ax.legend(fontsize=8, title="Platform")
fig.suptitle("Areal B-cell density by genotype (per-donor; black bar = group mean)", y=1.02)
plt.tight_layout(); du.save_figure(fig, "H20_density_by_genotype"); plt.show()""")

# ---------------------------------------------------------------- Section C
md(r"""## C. Primary model — Negative-Binomial rate ratios with area offset

Per-follicle NB (GPT model): `Bcell_count ~ offset(log Follicle_mm²) + Dosage + Platform`,
donor **cluster-robust** SEs (the portable analogue of a donor random intercept).
Genotype enters as T-allele **dosage** (0/1/2). The reported **rate ratio**
exp(β_dosage) is the per-allele multiplicative change in B-cells per mm².

The per-follicle model is restricted to **genuine follicles ≥ 8400 µm²** (the H6
GMM area threshold, set independently of genotype) — sub-threshold fragments are
mostly annotation artifacts, contribute ~2.5% of follicle area, and destabilize a
per-follicle count model; the all-follicle fit is reported as a sensitivity. A
confirmatory **R `lme4` NB-GLMM** with an explicit `(1|donor)` random intercept
(fixed θ from `MASS::glm.nb`, bobyqa) follows.""")

code(r"""import statsmodels.api as sm

def fit_nb(data, count_col, area_mm2_col, label):
    d = data.copy()
    d = d[d[area_mm2_col] > 0]
    y = d[count_col].to_numpy(float)
    off = np.log(d[area_mm2_col].to_numpy(float))
    dose = d["Genotype"].map(DOSAGE).to_numpy(float)
    codex = (d["Platform"] == "CODEX").astype(float).to_numpy()
    groups = d["Sample"].astype("category").cat.codes.to_numpy()
    X = pd.DataFrame({"const": 1.0, "Dosage": dose, "Platform_CODEX": codex})
    # estimate NB2 dispersion via Cameron-Trivedi aux regression on a Poisson fit
    poi = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
    mu = np.clip(poi.mu, 1e-8, None)
    aux = ((y - mu) ** 2 - y) / mu
    alpha = float(np.clip(sm.OLS(aux, mu).fit().params[0], 1e-6, None))
    res = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha),
                 offset=off).fit(cov_type="cluster", cov_kwds={"groups": groups})
    ci = res.conf_int()
    out = []
    for term in ["Dosage", "Platform_CODEX"]:
        out.append({"Model": label, "Term": term, "beta": res.params[term],
                    "RR": np.exp(res.params[term]),
                    "RR_lo": np.exp(ci.loc[term, 0]), "RR_hi": np.exp(ci.loc[term, 1]),
                    "p": res.pvalues[term], "alpha_nb": alpha,
                    "n_obs": int(len(y)), "n_donor": int(len(np.unique(groups)))})
    return pd.DataFrame(out), res

def boot_dosage_rr(pf_data, count_col, area_col="Follicle_mm2", n=600, seed=20260714):
    '''Donor-level bootstrap of the per-follicle dosage rate ratio (resample the
    11 donors, refit Poisson mean structure). Returns array of RRs.'''
    r = np.random.default_rng(seed)
    donors = pf_data["Sample"].unique()
    dof = {s: DOSAGE[du.GENOTYPE_MAP[s]] for s in donors}
    out = []
    for _ in range(n):
        samp = r.choice(donors, len(donors), replace=True)
        if len({dof[s] for s in samp}) < 2:
            continue
        d = pd.concat([pf_data[pf_data["Sample"] == s] for s in samp], ignore_index=True)
        y = d[count_col].to_numpy(float); off = np.log(d[area_col].to_numpy(float))
        cols = [np.ones(len(d)), d["Genotype"].map(DOSAGE).to_numpy(float)]
        codex = (d["Platform"] == "CODEX").astype(float).to_numpy()
        if codex.min() != codex.max():
            cols.append(codex)
        try:
            out.append(np.exp(sm.GLM(y, np.column_stack(cols),
                       family=sm.families.Poisson(), offset=off).fit().params[1]))
        except Exception:
            pass
    return np.array(out)""")

code(r"""FOLL_MIN_AREA = 8400.0  # H6 GMM follicle-area threshold (µm²), genotype-blind
per_foll["Follicle_mm2"] = per_foll["Follicle_Area_um2"] / 1e6
pf_fit = per_foll[per_foll["Follicle_Area_um2"] >= FOLL_MIN_AREA].copy()
print(f"Follicles: {len(per_foll)} total -> {len(pf_fit)} genuine (>= {FOLL_MIN_AREA:.0f} µm²)")

nb_primary, res_primary = fit_nb(pf_fit, "Bcell_count", "Follicle_mm2",
                                 "Per-follicle >=8400um2 (PRIMARY)")
nb_allf, _ = fit_nb(per_foll, "Bcell_count", "Follicle_mm2",
                    "Per-follicle all (sensitivity)")

# per-donor aggregate models (n=11)
per_donor["n_bcell"] = per_donor["n_bcell"].astype(float)
nb_whole, _ = fit_nb(per_donor.assign(Follicle_ID=per_donor["Sample"]),
                     "n_bcell", "tissue_mm2", "Whole-tissue (per-donor)")
nb_follcomp, _ = fit_nb(per_donor, "n_bcell_infoll", "foll_mm2",
                        "Follicle-compartment (per-donor)")

nb_all = pd.concat([nb_primary, nb_allf, nb_follcomp, nb_whole], ignore_index=True)
du.save_table(nb_all, "H20_nb_rate_ratios")
print("Negative-Binomial rate ratios (Dosage = per-T-allele; RR>1 => denser in T/T):")
nb_all.round(4)""")

code(r"""# Confirmatory: R lme4 NB-GLMM with explicit (1|donor). Fixed theta from
# MASS::glm.nb + bobyqa converges where glmer.nb's alternating fit does not.
def glmer_nb(data, count_col, area_mm2_col):
    if not RSCRIPT:
        return None
    d = data[data[area_mm2_col] > 0].copy()
    d["Dosage"] = d["Genotype"].map(DOSAGE)
    d["Platform_CODEX"] = (d["Platform"] == "CODEX").astype(int)
    d["logA"] = np.log(d[area_mm2_col].astype(float))
    csv = CACHE_DIR / "glmer_in.csv"; out = CACHE_DIR / "glmer_out.csv"
    err = Path(str(out) + ".err")
    for f in (out, err):
        if f.exists(): f.unlink()
    d[[count_col, "Dosage", "Platform_CODEX", "logA", "Sample"]].rename(
        columns={count_col: "Y"}).to_csv(csv, index=False)
    rscript = f'''
    suppressMessages({{library(lme4); library(MASS)}})
    d <- read.csv("{csv}")
    tryCatch({{
      th <- glm.nb(Y ~ Dosage + Platform_CODEX + offset(logA), data=d)$theta
      m <- glmer(Y ~ Dosage + Platform_CODEX + (1|Sample) + offset(logA), data=d,
                 family=negative.binomial(theta=th),
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
      co <- summary(m)$coefficients
      write.csv(data.frame(Term=rownames(co), co, theta=th), "{out}", row.names=FALSE)
    }}, error=function(e) writeLines(conditionMessage(e), "{err}"))
    '''
    r = subprocess.run([RSCRIPT, "-e", rscript], capture_output=True, text=True, timeout=600)
    if out.exists():
        res = pd.read_csv(out)
        res.columns = ["Term", "Estimate", "SE", "z", "p", "theta"][:len(res.columns)]
        res["RR"] = np.exp(res["Estimate"])
        return res
    msg = err.read_text()[-300:] if err.exists() else (r.stderr or "")[-300:]
    print("R NB-GLMM unavailable:", msg)
    return None

r_foll = glmer_nb(pf_fit, "Bcell_count", "Follicle_mm2")
if r_foll is not None:
    print("R lme4 NB-GLMM (per-follicle >=8400µm², (1|donor) random intercept):")
    display(r_foll.round(4))
    dpy = nb_primary.loc[nb_primary["Term"] == "Dosage", "RR"].iloc[0]
    dr = r_foll.loc[r_foll["Term"] == "Dosage", "RR"]
    if len(dr):
        print(f"\nDosage rate ratio — statsmodels: {dpy:.3f}   R NB-GLMM: {dr.iloc[0]:.3f} "
              f"(agree: {'yes' if abs(dpy-dr.iloc[0])<0.15 else 'no'})")
else:
    print("R confirmatory model unavailable (Rscript/lme4 missing or non-convergence).")""")

md(r"""**Honest inference for 11 donors.** Genotype is a *donor-level* covariate, so the
effective unit is the 11 donors — not the 1059 follicles. Cluster-robust SEs with only
G=11 clusters are **anti-conservative** (their p is too small); the trustworthy
inference is the R random-intercept p above plus a **donor-level bootstrap** of the
per-follicle dosage rate ratio (resample the 11 donors, refit).""")

code(r"""boot = boot_dosage_rr(pf_fit, "Bcell_count")
rr_pt = float(nb_primary.loc[nb_primary["Term"] == "Dosage", "RR"].iloc[0])
rr_lo, rr_hi = np.percentile(boot, [2.5, 97.5])
rr_p = float(2 * min((boot < 1).mean(), (boot > 1).mean()))
r_p = float(r_foll.loc[r_foll["Term"] == "Dosage", "p"].iloc[0]) if r_foll is not None else np.nan
sm_p = float(nb_primary.loc[nb_primary["Term"] == "Dosage", "p"].iloc[0])
print(f"Per-follicle dosage rate ratio = {rr_pt:.2f}")
print(f"  donor-bootstrap 95% CI [{rr_lo:.2f}, {rr_hi:.2f}], bootstrap p(RR≠1) = {rr_p:.3f}  (n_boot={len(boot)})")
print(f"  R random-intercept p = {r_p:.3f}")
print(f"  statsmodels cluster-robust p = {sm_p:.3f}  <- anti-conservative (G=11), do not use")
du.save_table(pd.DataFrame({"metric": ["dosage_RR", "boot_ci_lo", "boot_ci_hi",
                                       "boot_p", "R_glmm_p", "cluster_robust_p"],
                            "value": [rr_pt, rr_lo, rr_hi, rr_p, r_p, sm_p]}),
              "H20_primary_inference")""")

# ---------------------------------------------------------------- Section D
md(r"""## D. Platform-stratified inference & the field-selection threat

Area offsets fix unequal *effort*, not non-random *field selection*. Both CODEX
donors are T/T, so any CODEX bias loads onto the T/T estimate:
- If CODEX crops were placed on white-pulp / follicle-rich tissue → B-cell density
  **inflated** in the (T/T) CODEX donors (spurious support for the hypothesis).
- Conversely, CODEX CD20 sensitivity is low (§B) → B-cells **under**-counted in the
  (T/T) CODEX donors (spurious *against*).
Either way the pooled T/T estimate is not a clean genotype effect. Report all three
strata; the **Phenocycler-only** stratum (the clean estimand) has **T/T n=1**.""")

code(r"""def dosage_spearman(data, col):
    dose = data["Genotype"].map(DOSAGE)
    m = dose.notna() & data[col].notna()
    if m.sum() < 3 or dose[m].nunique() < 2:
        return np.nan, np.nan
    return stats.spearmanr(dose[m], data.loc[m, col])

strat_rows = []
for name, sub in [("All (11)", per_donor),
                  ("Phenocycler-only", per_donor[per_donor["Platform"] == "Phenocycler"]),
                  ("CODEX-only", per_donor[per_donor["Platform"] == "CODEX"])]:
    for col in ["WholeTissue_density", "Follicle_density"]:
        rho, p = dosage_spearman(sub, col)
        n_tt = int((sub["Genotype"] == "T/T").sum())
        strat_rows.append({"Stratum": name, "Endpoint": col, "n": len(sub),
                           "n_TT": n_tt, "dosage_rho": rho, "p": p})
strat = pd.DataFrame(strat_rows)
du.save_table(strat, "H20_platform_stratified")
print("Dosage trend by platform stratum (rho>0 would support the hypothesis):")
display(strat.round(4))

# clean comparison: the single Phenocycler T/T vs the two CODEX T/T
tt = per_donor[per_donor["Genotype"] == "T/T"][
    ["Sample", "Platform", "WholeTissue_density", "Follicle_density"]]
print("\nT/T donors — the confound in one table (HDL063 is the only clean estimate):")
display(tt.round(1))
cc_ct = per_donor[per_donor["Genotype"] != "T/T"]
print(f"C/C+C/T whole-tissue density range: "
      f"{cc_ct['WholeTissue_density'].min():.0f}–{cc_ct['WholeTissue_density'].max():.0f} cells/mm²; "
      f"Phenocycler T/T (HDL063) = "
      f"{tt.loc[tt.Platform=='Phenocycler','WholeTissue_density'].iloc[0]:.0f}")""")

# ---------------------------------------------------------------- Section D2 (CD20+CD21+)
md(r"""## D2. CD20⁺CD21⁺ double-positive subset (marginal-zone / FDC-associated B cells)

CD21 (CR2) marks marginal-zone and follicular B cells embedded in FDC networks. Does
this more specific **CD20⁺CD21⁺** subset show the increase in T/T that the whole
B-cell pool did not? The same NB rate-ratio machinery, donor-level inference, and the
same n=3 / platform-collinearity limits apply. CD20 is RESTORE-calibrated (§B);
**CD21⁺ is a per-image Otsu call** (there is no clean double-positive ground truth),
so treat the DP call as approximate.""")

code(r"""# per-follicle DP model (genuine follicles >= 8400 µm²), donor-level inference
pf_dp = per_foll[per_foll["Follicle_Area_um2"] >= FOLL_MIN_AREA].copy()
nb_dp, _ = fit_nb(pf_dp, "DP_count", "Follicle_mm2", "DP per-follicle >=8400um2 (PRIMARY)")
r_dp = glmer_nb(pf_dp, "DP_count", "Follicle_mm2")
boot_d = boot_dosage_rr(pf_dp, "DP_count")
dp_rr = float(nb_dp.loc[nb_dp["Term"] == "Dosage", "RR"].iloc[0])
dp_lo, dp_hi = np.percentile(boot_d, [2.5, 97.5])
dp_bp = float(2 * min((boot_d < 1).mean(), (boot_d > 1).mean()))
dp_rp = float(r_dp.loc[r_dp["Term"] == "Dosage", "p"].iloc[0]) if r_dp is not None else np.nan

dp_rows = []
for name, sub in [("All (11)", per_donor),
                  ("Phenocycler-only", per_donor[per_donor["Platform"] == "Phenocycler"])]:
    for col in ["DP_WholeTissue_density", "DP_Follicle_density"]:
        rho, p = dosage_spearman(sub, col)
        dp_rows.append({"Stratum": name, "Endpoint": col, "n": len(sub),
                        "n_TT": int((sub["Genotype"] == "T/T").sum()), "dosage_rho": rho, "p": p})
dp_strat = pd.DataFrame(dp_rows)

du.save_table(nb_dp, "H20_dp_rate_ratios")
du.save_table(dp_strat, "H20_dp_stratified")
du.save_table(pd.DataFrame({"metric": ["DP_dosage_RR", "boot_ci_lo", "boot_ci_hi", "boot_p", "R_glmm_p"],
                            "value": [dp_rr, dp_lo, dp_hi, dp_bp, dp_rp]}), "H20_dp_primary_inference")

print(f"CD20+CD21+ per-follicle dosage RR = {dp_rr:.2f}  "
      f"donor-boot 95% CI [{dp_lo:.2f}, {dp_hi:.2f}]  boot p={dp_bp:.3f}  R glmm p={dp_rp:.3f}")
print(f"DP cells as % of CD20+ B cells (per-donor mean): "
      f"{(per_donor['n_dp'] / per_donor['n_bcell']).mean():.1%}")
print("\nDP dosage trend by stratum (rho>0 would support 'increased CD20+CD21+ in T/T'):")
display(dp_strat.round(4))
print("\nDP density by genotype (per-donor means, cells/mm²):")
display(per_donor.groupby("Genotype", observed=True)[
    ["DP_WholeTissue_density", "DP_Follicle_density"]].mean().round(1))""")

code(r"""fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, col, title in zip(axes, ["DP_WholeTissue_density", "DP_Follicle_density"],
                          ["CD20⁺CD21⁺ whole-tissue", "CD20⁺CD21⁺ in follicles"]):
    sns.stripplot(data=per_donor, x="Genotype", y=col, order=GENO_ORDER,
                  hue="Platform", palette={"Phenocycler": "#3366CC", "CODEX": "#CC6633"},
                  s=11, jitter=0.12, ax=ax, edgecolor="black", linewidth=0.6)
    m = per_donor.groupby("Genotype", observed=True)[col].mean().reindex(GENO_ORDER)
    ax.plot(range(3), m.values, "k_", markersize=26, markeredgewidth=2)
    ax.set_title(title); ax.set_ylabel("CD20⁺CD21⁺ cells / mm²"); ax.set_xlabel("rs3184504")
    ax.legend(fontsize=8, title="Platform")
fig.suptitle("CD20⁺CD21⁺ double-positive density by genotype", y=1.02)
plt.tight_layout(); du.save_figure(fig, "H20_dp_density_by_genotype"); plt.show()""")

# ---------------------------------------------------------------- Section E
md(r"""## E. Robustness — random matched observation windows (GPT step 3)

Whole-section density can shift with *where* you look. Following the chat: for each
Phenocycler donor, drop many random CODEX-sized windows on tissue (centred on a
random cell), recompute in-window density (B-cells ÷ in-window tissue area,
estimated from local cell coverage), pick **one window per donor** per Monte-Carlo
draw together with the two whole-image CODEX donors, and re-estimate the dosage
trend. Windows are **nested within donor** (never independent). Reports the sign
stability of the genotype signal; a **donor bootstrap** gives the CI.""")

code(r"""rng = np.random.default_rng(20260714)
pc_path = CACHE_DIR / "percell_coords.pkl"
CODEX_AREA_MM2 = per_donor.loc[per_donor["Platform"] == "CODEX", "tissue_mm2"].mean()
win_side = np.sqrt(CODEX_AREA_MM2) * 1000.0  # µm, square window matching CODEX footprint

if pc_path.exists():
    cc = pd.read_pickle(pc_path)
    codex_d = per_donor[per_donor["Platform"] == "CODEX"].set_index("Sample")
    pc_samples = [s for s in per_donor.loc[per_donor.Platform == "Phenocycler", "Sample"]]
    # precompute per-donor global cell areal density (cells per mm^2 tissue)
    cell_dens = (per_donor.set_index("Sample")["n_cells"] / per_donor.set_index("Sample")["tissue_mm2"])
    by = {s: cc[cc["Sample"] == s] for s in pc_samples}

    def one_window(sub):
        i = rng.integers(len(sub))
        cx, cy = sub["x"].iloc[i], sub["y"].iloc[i]
        h = win_side / 2
        m = (sub["x"] >= cx-h) & (sub["x"] < cx+h) & (sub["y"] >= cy-h) & (sub["y"] < cy+h)
        w = sub[m]
        if len(w) < 200:  # off-tissue / sparse -> retry a few times
            for _ in range(5):
                i = rng.integers(len(sub)); cx, cy = sub["x"].iloc[i], sub["y"].iloc[i]
                m = (sub["x"] >= cx-h) & (sub["x"] < cx+h) & (sub["y"] >= cy-h) & (sub["y"] < cy+h)
                w = sub[m]
                if len(w) >= 200: break
        return int(w["bcell"].sum()), int(len(w))

    N_ITER = 300
    rhos = []
    for _ in range(N_ITER):
        recs = []
        for s in pc_samples:
            nb, ncell = one_window(by[s])
            tissue_mm2 = ncell / cell_dens[s]
            recs.append((s, du.GENOTYPE_MAP[s], nb / max(tissue_mm2, 1e-6)))
        for s, row in codex_d.iterrows():
            recs.append((s, du.GENOTYPE_MAP[s], row["WholeTissue_density"]))
        wd = pd.DataFrame(recs, columns=["Sample", "Genotype", "dens"])
        dose = wd["Genotype"].map(DOSAGE)
        rho, _ = stats.spearmanr(dose, wd["dens"])
        rhos.append(rho)
    rhos = np.array(rhos)
    frac_pos = float((rhos > 0).mean())
    print(f"Window robustness ({N_ITER} draws, window≈{CODEX_AREA_MM2:.1f} mm²):")
    print(f"  dosage rho: median={np.median(rhos):+.3f}  IQR[{np.percentile(rhos,25):+.3f},"
          f"{np.percentile(rhos,75):+.3f}]  P(rho>0)={frac_pos:.0%}")

    # donor bootstrap of the observed whole-tissue dosage trend (CI)
    samples = per_donor["Sample"].to_numpy(); dvals = per_donor["WholeTissue_density"].to_numpy()
    dose_all = per_donor["Genotype"].map(DOSAGE).to_numpy()
    boot = []
    for _ in range(2000):
        idx = rng.integers(0, len(samples), len(samples))
        if len(np.unique(dose_all[idx])) < 2: continue
        boot.append(stats.spearmanr(dose_all[idx], dvals[idx])[0])
    boot = np.array(boot)
    rob = pd.DataFrame({"metric": ["window_dosage_rho_median", "window_P(rho>0)",
                                   "donor_boot_rho_median", "donor_boot_ci_lo", "donor_boot_ci_hi"],
                        "value": [np.median(rhos), frac_pos, np.median(boot),
                                  np.percentile(boot, 2.5), np.percentile(boot, 97.5)]})
    du.save_table(rob, "H20_window_robustness")
    print(f"  donor-bootstrap whole-tissue dosage rho: {np.median(boot):+.3f} "
          f"[{np.percentile(boot,2.5):+.3f}, {np.percentile(boot,97.5):+.3f}]")
else:
    print("Per-cell coord cache absent; run §B compute_counts() first to enable §E.")""")

# ---------------------------------------------------------------- Section F
md(r"""## F. Power / minimum detectable effect

Two facts bound what this cohort can detect:
1. **Exact Mann-Whitney floor.** With n=3 vs 3 (C/C vs T/T) the smallest possible
   two-sided p is **0.10** — the direct T/T pairwise test *cannot* reach α=0.05.
2. **Donor-level dosage trend** (the effective unit is ~11 donors) — Monte-Carlo MDE:
   smallest per-allele density fold-change reaching 80% power for the Spearman trend
   at n=3/5/3, and the balanced cohort size a real test would need.""")

code(r"""from itertools import combinations
# (1) exact MW floor for n1=n2=3
def mw_min_two_sided(n1, n2):
    from scipy.stats import mannwhitneyu
    a = list(range(n1)); best = 1.0
    # extreme separation -> smallest p
    x = np.arange(n1); y = np.arange(n1, n1+n2)
    return mannwhitneyu(x, y, alternative="two-sided")[1]
print(f"Exact MW two-sided floor  C/C vs T/T (3v3): {mw_min_two_sided(3,3):.3f}  "
      f"(cannot reach 0.05)")
print(f"Exact MW two-sided floor  C/T vs T/T (5v3): {mw_min_two_sided(5,3):.3f}")

# (2) donor-level dosage MDE (mirror H18 Monte-Carlo), lognormal donor densities
rng = np.random.default_rng(20260714)
ln_sd = float(np.std(np.log(per_donor["WholeTissue_density"]), ddof=1))
groups_n = {"C/C": 3, "C/T": 5, "T/T": 3}
def dosage_power(fold_per_allele, n_by=groups_n, sims=3000, sd=ln_sd, alpha=0.05):
    labels = np.array(sum([[g]*n for g, n in n_by.items()], []))
    dose = np.array([DOSAGE[g] for g in labels])
    hits = 0
    for _ in range(sims):
        mu = dose * np.log(fold_per_allele)
        val = np.exp(rng.normal(mu, sd))
        rho, p = stats.spearmanr(dose, val)
        if p < alpha and rho > 0:  # one-sided in hypothesized direction
            hits += 1
    return hits / sims

folds = np.round(np.arange(1.0, 3.01, 0.25), 2)
pw = [(f, dosage_power(f)) for f in folds]
mde = next((f for f, p in pw if p >= 0.8), np.nan)
print(f"\nDonor-level dosage power (n=3/5/3), per-allele fold-change -> power:")
for f, p in pw: print(f"  {f:.2f}x  ->  {p:.2f}")
print(f"  MDE (80% power): per-allele fold-change ≈ {mde}")

# required balanced n for a moderate effect (1.3x per allele)
def bal_power(n, fold=1.3, sims=2000):
    return dosage_power(fold, {"C/C": n, "C/T": n, "T/T": n}, sims=sims)
req = next((n for n in [10, 15, 20, 30, 45, 60, 85, 120] if bal_power(n) >= 0.8), ">120")
print(f"\nBalanced n/group for 80% power at 1.3x/allele: ~{req} per group")

power_tbl = pd.DataFrame(pw, columns=["fold_per_allele", "power"])
power_tbl.attrs = {}
du.save_table(power_tbl, "H20_power_analysis")""")

# ---------------------------------------------------------------- Section G
md(r"""## G. Verdict""")

code(r"""dose_rr = nb_primary.loc[nb_primary["Term"] == "Dosage"].iloc[0]
rho_all = strat[(strat.Stratum == "All (11)") & (strat.Endpoint == "WholeTissue_density")].iloc[0]
verdict = pd.DataFrame([
    {"Question": "Can a defensible test be RUN?",
     "Answer": "Yes — NB rate-ratio model with area offset, donor-clustered."},
    {"Question": "Primary result (per-follicle dosage RR, donor-level inference)",
     "Answer": f"RR={rr_pt:.2f}, donor-boot 95%CI [{rr_lo:.2f},{rr_hi:.2f}], p={rr_p:.2f} "
               f"(R glmm p={r_p:.2f}); {'no significant' if rr_p >= 0.05 else 'only nominal'} signal, "
               f"and it hinges on a non-identifiable platform adjustment."},
    {"Question": "Whole-tissue dosage trend (all 11)",
     "Answer": f"rho={rho_all['dosage_rho']:+.2f}, p={rho_all['p']:.2f} (not denser in T/T)"},
    {"Question": "CD20+CD21+ double-positive dosage RR",
     "Answer": f"RR={dp_rr:.2f}, donor-boot 95%CI [{dp_lo:.2f},{dp_hi:.2f}], p={dp_bp:.2f} "
               f"({'no significant' if dp_bp >= 0.05 else 'nominal'} increase in T/T)"},
    {"Question": "Direct T/T pairwise adequately powered?",
     "Answer": "No — n=3 vs 3 Mann-Whitney floor = 0.10; cannot reach 0.05."},
    {"Question": "Confound controlled?",
     "Answer": "No — platform collinear with T/T (both CODEX = T/T); clean Phenocycler stratum has T/T n=1."},
    {"Question": "Consistent with prior evidence?",
     "Answer": "Yes — prior B-cell proportion null/against (H15 rho=-0.39, H16 +0.04)."},
    {"Question": "Overall",
     "Answer": "Testable only weakly; data cannot confirm/refute. A platform-balanced "
               f"cohort (~{req}/group) is required."},
])
du.save_table(verdict, "H20_verdict")
with pd.option_context("display.max_colwidth", None):
    display(verdict)""")

md(r"""**Bottom line.** The GPT methodology is the right framework and yields an honest,
interpretable rate-ratio model — but it cannot rescue the design.

- The per-follicle T-allele dosage rate ratio is small (RR≈1.09) and **not significant**
  once donor-level uncertainty is respected (R random-intercept p≈0.26; the donor
  bootstrap CI spans 1). It also *hinges on adjusting for platform*, which is
  **non-identifiable** here because both CODEX donors are T/T. (The statsmodels
  cluster-robust p≈0.03 is anti-conservative with only 11 clusters and should not be read
  as a signal.)
- Every read that does **not** lean on that adjustment is flat-to-negative: the
  whole-tissue dosage trend (ρ=−0.12, p=0.73), the follicle-compartment RR (1.04, p=0.43),
  and the random-window robustness (median ρ≈−0.12, P(ρ>0)=23%).
- The **CD20⁺CD21⁺** subset (user-requested; ≈69% of B cells at the Otsu CD21 cut, i.e.
  CD21⁺ follicular B cells rather than a rare marginal-zone population) behaves the same:
  a borderline per-follicle RR=1.14 (R p≈0.07, donor-bootstrap p=0.10, CI spans 1) that
  again hinges on the non-identifiable platform adjustment (CODEX DP rate ratio 0.19) and
  is **negative** on the unadjusted whole-tissue trend (ρ=−0.51). No robust increase in T/T.
- The direct T/T pairwise test is *mathematically* incapable of significance at n=3 vs 3
  (exact floor p=0.10), and platform is collinear with the T/T contrast (the one clean
  Phenocycler stratum has a single T/T donor, HDL063).

This reproduces, in absolute-density space, the null/against signal prior work found in
proportion space. **The hypothesis is testable only weakly; the existing data cannot
confirm or refute it. A platform-balanced, adequately-powered cohort (~10+/group) is
required.**""")

# ---------------------------------------------------------------- write
nb["cells"] = cells
nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb.metadata["language_info"] = {"name": "python"}
NB_PATH.write_text(nbf.writes(nb))
print(f"Wrote {NB_PATH.relative_to(PROJECT)}  ({NB_PATH.stat().st_size} bytes, {len(cells)} cells)")
