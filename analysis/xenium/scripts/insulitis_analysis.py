"""Per-phenotype dynamic insulitis grading with rotation-null calibration.

Inputs may be supplied as explicit ``--paths`` to ``*_phenotyped.h5ad`` files.
Each immune sidecar is resolved beside its base as
``<serial>_immune_phenotyped.h5ad``. This supports the combined repository's
``processed/<serial>/<roi>/`` layout while retaining the historical flat
``--samples`` resolution.

Each immune subtype gets:
  - per-islet per-100-endocrine rate
  - rotation-null distribution (1000 global rotations of immune positions)
  - per-(sample x size_class) 95th/99th percentile thresholds
  - per-islet grade {no_insulitis, peri_insulitis, insulitis}

The Campbell-Thompson 3/6 absolute grade is preserved as `insulitis_grade_absolute`
on `total_immune`. The single weighted-index grade from earlier drafts was
dropped per user decision (2026-04-28).

Outputs (single writer of these files):
  data/processed/islet_infiltration_per100endo.csv  per-islet per-phenotype grades
  data/processed/islet_insulitis_grades.csv         long-format prevalence summary
  data/processed/islet_insulitis_thresholds.csv     per-stratum thresholds
  data/processed/islet_insulitis_regression.csv     per-(sample,subtype) odds ratios
  data/processed/immune_proximity_summary.csv       density enrichment (UNCHANGED logic)
"""
import argparse
import os
import shutil
import sys
import time
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
SAMPLES = ("0041323", "0041326")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.workflow_paths import resolve_insulitis_jobs

ENDO_LABELS = {
    "Endocrine", "Beta", "Alpha", "Delta", "Endocrine_pan", "Gamma", "Epsilon",
    # Include cells flagged Indeterminate by the strict 30%-panel-coverage gate
    # — they're still likely endocrine, just with weaker multi-marker support
    # (especially in low-detection sample 0041326).
    "Indeterminate_Endocrine_lowcoverage",
    "Indeterminate_Beta_lowcoverage",
    "Indeterminate_Alpha_lowcoverage",
    "Indeterminate_Delta_lowcoverage",
    "Indeterminate_endocrine",
    "Indeterminate_Endocrine",
}
EPS_UM = 50.0
MIN_SAMPLES = 10

ISLET_ZONE_UM = 50.0
PERI_PROXIMAL_UM = 200.0

IMMUNE_SUBTYPES = ["B_pan", "B_plasma", "DC",
                    "Macro_M1", "Macro_M2", "Macro_resident",
                    "Monocyte", "NK",
                    "T_cytotoxic", "T_exhausted", "T_helper", "T_reg"]
NON_IMMUNE_DROP = {"Acinar", "Schwann", "Endothelial"}

CT_PERI_THRESH = 3
CT_INSULITIS_THRESH = 6

N_ROTATIONS = 1000
RNG_SEED = 0
COMPOSITION_MARGIN = 0.1
COMPOSITION_FRAC = 0.6

SIZE_BINS = [(0, 20, "tiny"),
             (21, 50, "small"),
             (51, 100, "med"),
             (101, 200, "large"),
             (201, np.inf, "xl")]
SIZE_CLASSES = [b[2] for b in SIZE_BINS]


def size_class(n):
    for lo, hi, name in SIZE_BINS:
        if lo <= n <= hi:
            return name
    return SIZE_CLASSES[-1]


def dbscan_islets(coords):
    from sklearn.cluster import DBSCAN
    return np.asarray(DBSCAN(eps=EPS_UM, min_samples=MIN_SAMPLES,
                              n_jobs=-1).fit_predict(coords))


def compose_per_islet(scores, labels, n_islets):
    """Per-islet composition from (Beta, Alpha, Delta) scores.

    Parameters
    ----------
    scores : (n_endo, 3) ndarray
    labels : (n_endo,) DBSCAN labels (-1 = noise)
    n_islets : number of clusters (max(labels)+1)

    Returns
    -------
    DataFrame with one row per islet: n_beta, n_alpha, n_delta,
    n_endo_unresolved, *_frac, composition_class.
    """
    top2 = np.partition(scores, -2, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]
    cls = np.argmax(scores, axis=1)
    cls = np.where(margin >= COMPOSITION_MARGIN, cls, -1)

    rows = []
    for k in range(n_islets):
        in_isl = labels == k
        n_total = int(in_isl.sum())
        c = cls[in_isl]
        n_beta = int((c == 0).sum())
        n_alpha = int((c == 1).sum())
        n_delta = int((c == 2).sum())
        n_unres = int((c == -1).sum())
        assert n_beta + n_alpha + n_delta + n_unres == n_total, \
            f"composition sum mismatch on islet {k}"
        bf = n_beta / n_total
        af = n_alpha / n_total
        df_frac = n_delta / n_total
        if bf >= COMPOSITION_FRAC:
            comp = "Beta_rich"
        elif af >= COMPOSITION_FRAC:
            comp = "Alpha_rich"
        elif df_frac >= COMPOSITION_FRAC:
            comp = "Delta_rich"
        else:
            comp = "Mixed"
        rows.append({"n_beta": n_beta, "n_alpha": n_alpha, "n_delta": n_delta,
                      "n_endo_unresolved": n_unres,
                      "beta_frac": bf, "alpha_frac": af, "delta_frac": df_frac,
                      "composition_class": comp})
    return pd.DataFrame(rows)


def assign_to_islet(immune_xy, seed_xy, seed_label):
    """For each immune cell, return (in_zone, nearest_islet_idx).

    in_zone = distance to nearest seed < ISLET_ZONE_UM.
    nearest_islet_idx = islet label of nearest seed (only meaningful if in_zone).
    """
    from scipy.spatial import cKDTree
    tree = cKDTree(seed_xy)
    d, idx = tree.query(immune_xy, k=1,
                          distance_upper_bound=ISLET_ZONE_UM)
    in_zone = np.isfinite(d)
    nearest = np.where(in_zone, seed_label[np.where(in_zone, idx, 0)], -1)
    return in_zone, nearest, d


def rotate_xy(r, theta, dtheta, center):
    new_theta = theta + dtheta
    return np.stack([r * np.cos(new_theta), r * np.sin(new_theta)],
                       axis=1) + center


def process_sample(sample, phenotyped_path=None, immune_path=None):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}")
    rng = np.random.default_rng(RNG_SEED)

    # ----- Load endocrine seeds + composition scores -----
    p = (Path(phenotyped_path) if phenotyped_path is not None else
         DATA_ROOT / f"processed/{sample}/{sample}_phenotyped.h5ad")
    print(f"  loading {p.name} (backed=r)...")
    t = time.time()
    a = sc.read_h5ad(p, backed="r")
    print(f"    loaded shape={a.shape} in {time.time()-t:.1f}s")
    coords_all = np.asarray(a.obsm["spatial"])
    ct_lineage = a.obs["celltype_lineage"].astype(str).values
    is_endo = np.isin(ct_lineage, list(ENDO_LABELS))
    seed_xy_all = coords_all[is_endo]
    seed_scores = np.stack([
        a.obs["Beta_score"].values[is_endo].astype(np.float32),
        a.obs["Alpha_score"].values[is_endo].astype(np.float32),
        a.obs["Delta_score"].values[is_endo].astype(np.float32),
    ], axis=1)
    print(f"    endocrine-labeled cells: {is_endo.sum():,}")
    a.file.close()

    # ----- DBSCAN islets -----
    print("  DBSCAN islet identification...")
    labels_all = dbscan_islets(seed_xy_all)
    in_cluster = labels_all >= 0
    seed_xy = seed_xy_all[in_cluster]
    seed_label = labels_all[in_cluster]
    seed_scores_in = seed_scores[in_cluster]
    n_islets = int(seed_label.max() + 1) if (seed_label >= 0).any() else 0
    print(f"    n_islets={n_islets}, clustered seeds={seed_xy.shape[0]:,}")

    # ----- Composition per islet -----
    comp_df = compose_per_islet(seed_scores_in, seed_label, n_islets)

    # ----- Per-islet metadata -----
    n_endo_per = np.bincount(seed_label, minlength=n_islets).astype(np.int64)
    cent_x = np.bincount(seed_label, weights=seed_xy[:, 0],
                            minlength=n_islets) / np.maximum(n_endo_per, 1)
    cent_y = np.bincount(seed_label, weights=seed_xy[:, 1],
                            minlength=n_islets) / np.maximum(n_endo_per, 1)
    sz_class = np.array([size_class(n) for n in n_endo_per])

    islet_meta = pd.DataFrame({
        "islet_id": [f"{sample}_islet_{k:04d}" for k in range(n_islets)],
        "sample": sample,
        "n_endocrine": n_endo_per,
        "centroid_x": cent_x,
        "centroid_y": cent_y,
        "size_class": sz_class,
    })
    islet_meta = pd.concat([islet_meta, comp_df], axis=1)

    # ----- Load immune cells -----
    p_imm = (Path(immune_path) if immune_path is not None else
             p.with_name(f"{sample}_immune_phenotyped.h5ad"))
    a_imm = sc.read_h5ad(p_imm)
    if "sample" in a_imm.obs.columns:
        a_imm = a_imm[a_imm.obs["sample"].astype(str) == sample].copy()
    immune_xy_all = np.asarray(a_imm.obsm["spatial"])
    immune_subtype_all = a_imm.obs["immune_subtype"].astype(str).values
    print(f"    immune cells (raw): {a_imm.n_obs:,}")

    sub_to_code = {s: i for i, s in enumerate(IMMUNE_SUBTYPES)}
    imm_code_all = np.array([sub_to_code.get(s, -1) for s in immune_subtype_all])
    keep_imm = imm_code_all >= 0
    immune_xy = immune_xy_all[keep_imm]
    imm_code = imm_code_all[keep_imm]
    n_imm = len(imm_code)
    print(f"    immune cells (after dropping {sorted(NON_IMMUNE_DROP)}): "
          f"{n_imm:,}")

    # Distance bin tags for distal density-enrichment table
    in_zone_all, nearest_all, d_all = assign_to_islet(immune_xy_all, seed_xy,
                                                       seed_label)
    in_proximal_all = (d_all >= ISLET_ZONE_UM) & (d_all < PERI_PROXIMAL_UM)
    in_distal_all = d_all >= PERI_PROXIMAL_UM

    # ----- Real per-phenotype counts -----
    in_zone, nearest, d = assign_to_islet(immune_xy, seed_xy, seed_label)
    real_counts = np.zeros((n_islets, len(IMMUNE_SUBTYPES)), dtype=np.int64)
    in_idx = np.where(in_zone)[0]
    np.add.at(real_counts, (nearest[in_idx], imm_code[in_idx]), 1)

    real_per100 = (100.0 * real_counts /
                    np.maximum(n_endo_per[:, None], 1)).astype(np.float32)

    # ----- Rotation null -----
    print(f"  rotation null: {N_ROTATIONS} rotations on {n_imm:,} immune cells "
          f"× {n_islets} islets × {len(IMMUNE_SUBTYPES)} subtypes...")
    t_rot = time.time()
    center = seed_xy.mean(axis=0)
    imm_centered = immune_xy - center
    r_imm = np.hypot(imm_centered[:, 0], imm_centered[:, 1])
    theta_imm = np.arctan2(imm_centered[:, 1], imm_centered[:, 0])

    from scipy.spatial import cKDTree
    seed_tree = cKDTree(seed_xy)

    null_counts = np.zeros((N_ROTATIONS, n_islets, len(IMMUNE_SUBTYPES)),
                              dtype=np.int32)
    for rot_idx in range(N_ROTATIONS):
        dtheta = rng.uniform(0.0, 2.0 * np.pi)
        rotated = rotate_xy(r_imm, theta_imm, dtheta, center)
        d_rot, idx_rot = seed_tree.query(rotated, k=1,
                                            distance_upper_bound=ISLET_ZONE_UM)
        in_zone_rot = np.isfinite(d_rot)
        if in_zone_rot.any():
            valid_idx = np.where(in_zone_rot)[0]
            isl_rot = seed_label[idx_rot[valid_idx]]
            np.add.at(null_counts[rot_idx],
                       (isl_rot, imm_code[valid_idx]), 1)
        if (rot_idx + 1) % 200 == 0:
            elapsed = time.time() - t_rot
            print(f"    rot {rot_idx+1}/{N_ROTATIONS}  "
                  f"elapsed {elapsed:.0f}s  "
                  f"eta {elapsed/(rot_idx+1)*(N_ROTATIONS-rot_idx-1):.0f}s")
    print(f"    rotation null total: {time.time()-t_rot:.0f}s")

    null_per100 = (100.0 * null_counts /
                    np.maximum(n_endo_per[None, :, None], 1)).astype(np.float32)

    # ----- Per-stratum thresholds -----
    thresh_rows = []
    thresh_lookup = {}  # (sz_class, subtype) -> (peri, insulitis)
    for sc_name in SIZE_CLASSES:
        islet_mask = sz_class == sc_name
        if islet_mask.sum() == 0:
            for s_idx, S in enumerate(IMMUNE_SUBTYPES):
                thresh_rows.append({
                    "sample": sample, "size_class": sc_name, "subtype": S,
                    "n_islets": 0, "n_null_samples": 0,
                    "thresh_peri": np.nan, "thresh_insulitis": np.nan,
                })
                thresh_lookup[(sc_name, S)] = (np.inf, np.inf)
            continue
        for s_idx, S in enumerate(IMMUNE_SUBTYPES):
            pool = null_per100[:, islet_mask, s_idx].ravel()
            thresh_peri = float(np.percentile(pool, 95))
            thresh_insul = float(np.percentile(pool, 99))
            thresh_rows.append({
                "sample": sample, "size_class": sc_name, "subtype": S,
                "n_islets": int(islet_mask.sum()),
                "n_null_samples": int(pool.size),
                "thresh_peri": thresh_peri,
                "thresh_insulitis": thresh_insul,
            })
            thresh_lookup[(sc_name, S)] = (thresh_peri, thresh_insul)

    # ----- Per-phenotype grading -----
    grade_cols = {}
    for s_idx, S in enumerate(IMMUNE_SUBTYPES):
        col = []
        for k in range(n_islets):
            n_zone = int(real_counts[k, s_idx])
            rate = real_per100[k, s_idx]
            tp, ti = thresh_lookup[(sz_class[k], S)]
            if n_zone >= 2 and rate >= ti:
                col.append("insulitis")
            elif n_zone >= 1 and rate >= tp:
                col.append("peri_insulitis")
            else:
                col.append("no_insulitis")
        grade_cols[f"{S}_grade"] = col

    # ----- Build per-islet output -----
    per = islet_meta.copy()
    for s_idx, S in enumerate(IMMUNE_SUBTYPES):
        per[f"{S}_zone"] = real_counts[:, s_idx].astype(np.int64)
    for s_idx, S in enumerate(IMMUNE_SUBTYPES):
        per[f"{S}_per100endo"] = real_per100[:, s_idx]
    for S in IMMUNE_SUBTYPES:
        per[f"{S}_grade"] = grade_cols[f"{S}_grade"]
    per["total_immune"] = real_counts.sum(axis=1).astype(np.int64)
    per["total_per100endo"] = (100.0 * per["total_immune"]
                                / np.maximum(per["n_endocrine"], 1)).astype(np.float32)
    per["insulitis_grade_absolute"] = pd.Series(per["total_immune"]).apply(
        lambda n: ("insulitis" if n >= CT_INSULITIS_THRESH
                    else "peri_insulitis" if n >= CT_PERI_THRESH
                    else "no_insulitis"))

    # ----- Density enrichment table (unchanged from previous logic) -----
    enrichment_rows = []
    n_zone_total = int(in_zone_all.sum())
    n_distal_total = int(in_distal_all.sum())
    for st in sorted(set(immune_subtype_all)):
        is_st = immune_subtype_all == st
        n_zone_st = int((is_st & in_zone_all).sum())
        n_distal_st = int((is_st & in_distal_all).sum())
        zone_frac = n_zone_st / max(n_zone_total, 1)
        distal_frac = n_distal_st / max(n_distal_total, 1)
        ratio = (zone_frac / distal_frac) if distal_frac > 0 else float("nan")
        enrichment_rows.append({
            "sample": sample, "immune_subtype": st,
            "n_in_islet_zone": n_zone_st,
            "n_in_proximal": int((is_st & in_proximal_all).sum()),
            "n_in_distal": n_distal_st,
            "zone_frac": zone_frac, "distal_frac": distal_frac,
            "density_enrichment": ratio,
            "n_islet_zone_total": n_zone_total,
            "n_distal_total": n_distal_total,
        })
    enrichment_df = pd.DataFrame(enrichment_rows)

    # ----- Per-(sample,subtype) regression -----
    reg_rows = regression_rows(per, sample)

    return per, enrichment_df, pd.DataFrame(thresh_rows), reg_rows


def regression_rows(per, sample):
    """Per-subtype descriptive logistic regression.

    is_<S>_insulitis ~ log(n_endocrine) + composition_class.
    Returns list of dicts with odds_ratio + 95% CI. NO p-values
    — comment in CSV header notes n=2 caveat.
    """
    try:
        import statsmodels.api as sm
    except ImportError:
        print("    statsmodels not available; skipping regression")
        return []
    rows = []
    log_n = np.log(np.maximum(per["n_endocrine"].astype(float), 1.0))
    comp = per["composition_class"].astype(str).values
    base = pd.DataFrame({
        "intercept": 1.0,
        "log_n_endocrine": log_n,
        "is_alpha_rich": (comp == "Alpha_rich").astype(int),
        "is_delta_rich": (comp == "Delta_rich").astype(int),
        "is_mixed": (comp == "Mixed").astype(int),
    })
    for S in IMMUNE_SUBTYPES:
        y = (per[f"{S}_grade"] == "insulitis").astype(int).values
        if y.sum() < 5 or (y == 0).sum() < 5:
            continue
        try:
            mod = sm.Logit(y, base.values).fit(disp=False, maxiter=200)
            params = mod.params
            ci = mod.conf_int()
            for i, term in enumerate(base.columns):
                rows.append({
                    "sample": sample, "subtype": S, "term": term,
                    "odds_ratio": float(np.exp(params[i])),
                    "ci_low": float(np.exp(ci[i][0])),
                    "ci_high": float(np.exp(ci[i][1])),
                })
        except Exception as e:
            print(f"    regression {S}: {type(e).__name__}: {e}")
    return rows


# =====================================================================
# Run all samples and emit outputs
# =====================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Per-phenotype insulitis grading with rotation-null calibration "
                    "(pancreas).")
    ap.add_argument("--samples", nargs="+", default=list(SAMPLES),
                    help="pancreas sample IDs to grade (default: %(default)s)")
    ap.add_argument(
        "--paths", nargs="+", default=None,
        help="explicit *_phenotyped.h5ad paths (overrides --samples; the matching "
             "immune sidecar must be in the same directory)",
    )
    args = ap.parse_args(argv)
    try:
        jobs = resolve_insulitis_jobs(
            samples=args.samples, paths=args.paths, data_root=DATA_ROOT
        )
    except (FileNotFoundError, ValueError) as exc:
        ap.error(str(exc))
    samples = [sample for sample, _, _ in jobs]

    all_per, all_enr, all_thr, all_reg = [], [], [], []
    for sample, phenotyped, immune in jobs:
        per, enr, thr, reg = process_sample(
            sample, phenotyped_path=phenotyped, immune_path=immune
        )
        all_per.append(per)
        all_enr.append(enr)
        all_thr.append(thr)
        all_reg.extend(reg)

    per_df = pd.concat(all_per, ignore_index=True)
    enr_df = pd.concat(all_enr, ignore_index=True)
    thr_df = pd.concat(all_thr, ignore_index=True)
    reg_df = pd.DataFrame(all_reg)

    # ----- Long-format grade summary -----
    grade_rows = []
    for s in samples:
        sub = per_df[per_df["sample"] == s]
        for S in IMMUNE_SUBTYPES:
            vc = sub[f"{S}_grade"].value_counts()
            n_no = int(vc.get("no_insulitis", 0))
            n_per = int(vc.get("peri_insulitis", 0))
            n_in = int(vc.get("insulitis", 0))
            tot = n_no + n_per + n_in
            grade_rows.append({
                "sample": s, "subtype": S,
                "no_insulitis": n_no, "peri_insulitis": n_per,
                "insulitis": n_in, "total_islets": tot,
                "pct_insulitis": 100 * n_in / max(tot, 1),
                "pct_peri_or_insulitis": 100 * (n_in + n_per) / max(tot, 1),
            })
        # Campbell-Thompson absolute on total_immune
        vc = sub["insulitis_grade_absolute"].value_counts()
        n_no = int(vc.get("no_insulitis", 0))
        n_per = int(vc.get("peri_insulitis", 0))
        n_in = int(vc.get("insulitis", 0))
        tot = n_no + n_per + n_in
        grade_rows.append({
            "sample": s, "subtype": "absolute_total",
            "no_insulitis": n_no, "peri_insulitis": n_per,
            "insulitis": n_in, "total_islets": tot,
            "pct_insulitis": 100 * n_in / max(tot, 1),
            "pct_peri_or_insulitis": 100 * (n_in + n_per) / max(tot, 1),
        })
    grade_df = pd.DataFrame(grade_rows)

    # ----- Archive previous CSVs -----
    out = DATA_ROOT / "processed"
    for fn in ("islet_infiltration_per100endo.csv",
                "islet_insulitis_grades.csv",
                "immune_proximity_summary.csv"):
        p = out / fn
        if p.exists() and not (out / f"{p.stem}_preweighted.csv").exists():
            shutil.copy2(p, out / f"{p.stem}_preweighted.csv")
    print("\n  archived previous CSVs as *_preweighted.csv (one-shot)")

    per_df.to_csv(out / "islet_infiltration_per100endo.csv", index=False)
    print(f"  wrote islet_infiltration_per100endo.csv  "
          f"({len(per_df):,} islets, {per_df.shape[1]} cols)")
    grade_df.to_csv(out / "islet_insulitis_grades.csv", index=False)
    print(f"  wrote islet_insulitis_grades.csv  "
          f"({len(grade_df)} rows: {len(samples)}×{len(IMMUNE_SUBTYPES)+1})")
    thr_df.to_csv(out / "islet_insulitis_thresholds.csv", index=False)
    print(f"  wrote islet_insulitis_thresholds.csv  "
          f"({len(thr_df)} rows: {len(samples)}×{len(SIZE_CLASSES)}×"
          f"{len(IMMUNE_SUBTYPES)})")
    enr_df.to_csv(out / "immune_proximity_summary.csv", index=False)
    print(f"  wrote immune_proximity_summary.csv  ({len(enr_df)} rows)")

    if len(reg_df) > 0:
        with open(out / "islet_insulitis_regression.csv", "w") as f:
            f.write("# Descriptive only — n=2 samples, do not interpret as inferential\n")
            reg_df.to_csv(f, index=False)
        print(f"  wrote islet_insulitis_regression.csv  ({len(reg_df)} rows)")

    # ----- Print summary -----
    print(f"\n{'='*72}\n=== PER-PHENOTYPE INSULITIS PREVALENCE\n{'='*72}")
    print("  Stratum: (sample × size_class), 5 size bins.")
    print("  Threshold: rotation null, 1000 rotations, 99th percentile.\n")
    pivot = grade_df.pivot(index="subtype", columns="sample",
                              values="pct_insulitis")
    print(pivot.round(1).to_string())

    print("\n=== ABSOLUTE (Campbell-Thompson 3/6 on total_immune) ===")
    abs_rows = grade_df[grade_df["subtype"] == "absolute_total"]
    print(abs_rows[["sample", "no_insulitis", "peri_insulitis",
                     "insulitis", "total_islets", "pct_insulitis"]].to_string(index=False))


if __name__ == "__main__":
    main()
