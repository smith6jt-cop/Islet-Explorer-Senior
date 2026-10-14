"""B2 follow-up: per-category and cross-tissue-comparable control-set rates on the pLN bundles.

Reads the already-computed per-gene table; no matrices are reloaded. Adds
  * CLEAN_NEG   -- the ESTABLISHED pancreas curation from eval_detection.py, verbatim, so the
                   pLN number is on the same terms as the published pancreas 0.097% -> 3.54%.
  * CLEAN_NEG_LN -- the same list re-derived for lymph node: minus NLRP7 (47 -> 46 genes; human macrophage
                   inflammasome sensor; a lymph node is full of macrophages).
  * median as well as mean over each set -- the mean is moved by a single non-specific probe.
"""
import pathlib, sys
import numpy as np, pandas as pd

OUT = pathlib.Path("/blue/maigan/smith6jt/Xenium_Analysis/data/processed/pln")
sys.path.insert(0, "/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler/scripts/xenium_celltype")

CLEAN_NEG = (
    "ABCA4","ARR3","CNGB1","CNGB3","GUCA1A","GUCY2D","IMPG1","NR2E3","OPN1LW","PDE6H",
    "PRPH2","RPE65","VSX2","DSG1","DSG3","TGM3","COL2A1","PHEX","TNNI3","TNNT2","NPPA",
    "NPPB","TBX5","AKAP4","BOLL","DDX25","DPY19L2","MAGEA1","MAGEA10","MAGEA12","MAGEA3",
    "MAGEA4","MEIOB","PIWIL1","PRM1","SPATA16","SPATA22","SPATA4","STK31","SYCP1","TEX14",
    "DAZL","GDF9","KHDC3L","NLRP5","NLRP7","TLE6",
)
CLEAN_NEG_LN = tuple(g for g in CLEAN_NEG if g != "NLRP7")

SER = [("6533","0080870","same_slide"),("6533","0041393","different_slide"),
       ("112","0081366","same_slide"),("112","0041396","different_slide")]

pg = pd.read_csv(OUT/"B3_per_gene_detection.csv").set_index("gene")
c  = pd.read_csv(OUT/"B2_control_genes_per_gene.csv")

def stats(genes, name, scope):
    g = [x for x in genes if x in pg.index]
    r = dict(set_name=name, scope=scope, n_genes=len(g), genes=";".join(sorted(g)))
    for _, s, _d in SER:
        r[f"mean_det_pct_ge1_{s}"]   = float(pg.loc[g, f"det_frac_ge1_{s}"].mean()*100)
        r[f"median_det_pct_ge1_{s}"] = float(pg.loc[g, f"det_frac_ge1_{s}"].median()*100)
        r[f"mean_det_pct_ge2_{s}"]   = float(pg.loc[g, f"det_frac_ge2_{s}"].mean()*100)
        r[f"mean_counts_per_cell_{s}"] = float(pg.loc[g, f"counts_per_cell_{s}"].mean())
    for donor, ss, ds in (("6533","0080870","0041393"),("112","0081366","0041396")):
        r[f"fold_same_over_diff_mean_{donor}"]   = r[f"mean_det_pct_ge1_{ss}"]/r[f"mean_det_pct_ge1_{ds}"]
        r[f"fold_same_over_diff_median_{donor}"] = (r[f"median_det_pct_ge1_{ss}"]/r[f"median_det_pct_ge1_{ds}"]
                                                    if r[f"median_det_pct_ge1_{ds}"]>0 else np.nan)
    return r

rows = []
rows.append(stats(("TBL1Y","USP9Y","EIF1AY"), "Y_zero_template_clean", "a_priori"))
rows.append(stats(("DDX3Y","KDM5D"), "Y_zero_template_crosshyb", "a_priori"))
rows.append(stats(("TSPY1",), "Y_zero_template_other", "a_priori"))
rows.append(stats(CLEAN_NEG, "CLEAN_NEG_pancreas_curation", "as_established"))
rows.append(stats(CLEAN_NEG_LN, "CLEAN_NEG_lymphnode_rederived", "ln_valid_only"))

for cat, grp in c[~c.set_name.str.startswith("Y_")].groupby("set_name"):
    rows.append(stats(grp.gene.tolist(), cat, "all_audited"))
    ok = grp[grp.ln_control_valid.eq("yes")].gene.tolist()
    rows.append(stats(ok, cat, "ln_valid_only"))

lnv = c[(~c.set_name.str.startswith("Y_")) & c.ln_control_valid.eq("yes")].gene.tolist()
rows.append(stats(lnv, "ALL_audit_categories_pooled", "ln_valid_only"))
rows.append(stats(c[~c.set_name.str.startswith("Y_")].gene.tolist(),
                  "ALL_audit_categories_pooled", "all_audited"))
# the whole panel, for scale
rows.append(stats(list(pg.index), "ALL_5001_PANEL_GENES", "reference_not_a_control"))

out = pd.DataFrame(rows)
out.to_csv(OUT/"B2_control_sets_summary.csv", index=False)

show = ["set_name","scope","n_genes","mean_det_pct_ge1_0080870","mean_det_pct_ge1_0041393",
        "fold_same_over_diff_mean_6533","mean_det_pct_ge1_0081366","mean_det_pct_ge1_0041396",
        "fold_same_over_diff_mean_112"]
pd.set_option("display.width", 220); pd.set_option("display.float_format", lambda v: f"{v:,.4f}")
print(out[show].to_string(index=False))

# genes whose a-priori "no template" call the data contradicts, on the PhenoCycler-FREE sections
yf = pg.loc[["TBL1Y","USP9Y","EIF1AY"], ["det_frac_ge1_0041393","det_frac_ge1_0041396"]].mean().mean()
cand = c[(~c.set_name.str.startswith("Y_")) & c.ln_control_valid.eq("yes")].copy()
cand["diff_slide_mean_pct"] = (cand[["det_frac_ge1_0041393","det_frac_ge1_0041396"]].mean(axis=1)*100)
cand["x_above_Y_floor"] = cand["diff_slide_mean_pct"]/(yf*100)
cand = cand.sort_values("x_above_Y_floor", ascending=False)
cand[["gene","set_name","diff_slide_mean_pct","x_above_Y_floor"]].to_csv(
    OUT/"B2_control_genes_vs_Y_floor.csv", index=False)
print("\nTop a-priori controls by how far they sit above the zero-template Y floor on the "
      "PhenoCycler-FREE pLN sections (Y floor = %.4f%%):" % (yf*100))
print(cand[["gene","set_name","diff_slide_mean_pct","x_above_Y_floor"]].head(15).to_string(index=False))
