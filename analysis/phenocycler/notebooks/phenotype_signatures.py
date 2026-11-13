"""Measured marker summaries for the PhenoCycler notebook series.

Pancreas structural and endocrine markers follow typing/refine.py::OWN_MARKERS.
Immune markers follow phenotype_workflow.csv, the recorded macrophage rule and
the recovered protein profiles. pLN markers follow the independent tissue rule
table in core/panels/data/tissues/pancreatic_lymph_node/phenotype_rules.csv.
"""
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


AMBIGUOUS = "Endocrine (ambiguous subtype)"
BELOW_FLOOR = "Endocrine (unresolved subtype)"
DISPLAY_NAMES = {
    AMBIGUOUS: "Endocrine, competing hormones",
    BELOW_FLOOR: "CD99 endocrine, no hormone above floor",
}
HORMONES = ["GCG", "INS", "SST"]
PAN_STRUCTURE = HORMONES + ["CD99", "Ker8_18", "E_cadherin", "Pan_Cytokeratin",
                            "CD31", "Caveolin", "CD34", "Podoplanin", "B3TUBB", "CD56"]
PLN_STRUCTURE = ["B3TUBB", "SOX2", "Vimentin", "SMA", "Collagen_IV", "Podoplanin",
                 "CD31", "CD34", "Caveolin", "E_cadherin", "Pan_Cytokeratin", "Keratin_5", "TP63"]
IMMUNE = ["CD3e", "CD20", "CD79a", "CD68", "CD163", "Iba1", "CD11c", "CD209",
          "MPO", "CD66", "CD11b", "CD56", "Granzyme_B", "CD38", "CD57", "IDO1", "HLA_DR"]
T_CELL = ["CD3e", "CD4", "CD8", "FOXP3", "ICOS", "CD39", "PD_1", "TOX", "LAG3",
          "VISTA", "TCF_1", "Bcl_2", "CD107a", "Granzyme_B", "CD44", "CD57", "Ki67", "HLA_DR", "Vimentin"]
MACROPHAGE = ["CD68", "CD163", "M2Gal3", "CD4"]
PLN_MYELOID = ["CD68", "CD163", "M2Gal3", "Iba1", "CD206", "CD11b", "CD11c",
               "HLA_DR", "iNOS", "CD209", "CD14", "CD141"]
PAN_ORDER = ["Alpha", "Beta", "Delta", BELOW_FLOOR, "Exocrine", "Blood vessel", "Lymphatic",
             "Nerve", "Immune", "Unresolved (no lineage evidence)"]
PLN_ORDER = ["Endothelial", "Stromal", "Smooth_muscle", "Neural", "Immune", "Ambiguous candidate", "No rule match"]
PLN_STROMAL = {"Neural", "Nerve_fiber", "Schwann_cell", "Stromal", "Fibroblast", "FRC",
               "Pericyte", "Basement_membrane", "Smooth_muscle", "Endothelial", "Blood_vessel",
               "Lymphatic_vessel", "HEV", "Myofibroblast_CAF", "Epithelial", "Ductal_Epithelial", "Basal_Epithelial", "Progenitor_Epithelial"}


def names_in_order(labels, tissue):
    preferred = PAN_ORDER if tissue == "Pancreas" else PLN_ORDER
    return [x for x in preferred if x in labels] + sorted(set(labels).difference(preferred))


def marker_panels(tissue, level, labels):
    """Keep every named assignment visible; group long leaf lists by the relevant panel."""
    ordered = names_in_order(labels, tissue)
    structural = PAN_STRUCTURE if tissue == "Pancreas" else PLN_STRUCTURE
    immune = list(IMMUNE)
    if tissue == "pLN":
        immune = list(dict.fromkeys(immune + PLN_MYELOID))
    if level == "primary":
        # Every named population has its own distinguishing channels on these same axes.
        anchors = (["CD3e","CD20","CD68","Iba1","CD79a","MPO","CD66","CD11c","CD209"]
                   if tissue == "Pancreas" else ["CD3e","CD20","CD68","CD11c","CD11b"])
        return [("Primary lineage markers", ordered, list(dict.fromkeys(structural + anchors)))]
    if level == "immune_class":
        return [("Immune lineages", ordered, immune)]
    if level == "t_subtype":
        return [("T-cell phenotypes", ordered, T_CELL)]
    if level == "macrophage_profile":
        return [("Macrophage protein profiles", ordered, MACROPHAGE)]
    if level != "detail":
        raise ValueError(f"Unknown hierarchy level {level}")
    groups = OrderedDict((name, []) for name in ("Tissue", "Immune", "B cell", "T cell"))
    for label in ordered:
        if (tissue == "Pancreas" and label in PAN_ORDER) or label in PLN_STROMAL:
            key = "Tissue"
        elif tissue == "pLN" and label in {"B_cell", "Memory_B_cell", "GC_B_cell", "Naive_B_cell", "Plasma_cell", "Plasmablast"}:
            key = "B cell"
        elif label.startswith("T-") or label in {"T_cell", "Tfh", "Treg", "DN_T_cell", "Progenitor_exhausted", "Terminally_exhausted"} or any(
                term in label for term in ("CD4", "CD8", "Treg")):
            key = "T cell"
        else:
            key = "Immune"
        groups[key].append(label)
    return [(k, v, {"Tissue": structural, "Immune": immune, "B cell": ["CD20", "CD79a", "CD38", "Bcl_2", "Ki67"], "T cell": T_CELL}[k]) for k, v in groups.items() if v]


def measured_profiles(counts, evidence, *, tissue="Pancreas", level="primary", disease=None,
                      donor=None, weighting="donors"):
    """Return phenotype support and exact donor/cell denominators for every displayed measurement."""
    if weighting not in ("donors", "cells"):
        raise ValueError("Choose weighting='donors' or 'cells'.")
    scope = counts.loc[counts.tissue.eq(tissue) & counts.level.eq(level)].copy()
    if disease is not None:
        scope = scope.loc[scope.disease.eq(disease)]
    if donor is not None:
        scope = scope.loc[scope.donor.astype(str).eq(str(donor))]
    if scope.empty:
        raise ValueError(f"No {tissue} section in the selected scope")
    keys = ["tissue", "level", "donor", "disease", "label"]
    frame = evidence.merge(scope[keys + ["assessable","comparable"]], on=keys, validate="many_to_one")
    eligible = "assessable" if donor is not None else "comparable"
    frame = frame.loc[frame[eligible] & frame.n_measured.gt(0)].copy()
    scope.loc[~scope[eligible], "cells"] = np.nan
    support = scope.groupby("label", observed=True).agg(
        cells=("cells", lambda s: s.sum(min_count=1)), contributing=("cells", lambda s: s.gt(0).sum()),
        assessable=("cells","count"))
    frame["signal_sum"] = frame.mean_scaled * frame.n_measured
    values = frame.groupby(["label", "marker"], observed=True).agg(
        mean=("mean_scaled", "mean"), fraction=("positive_fraction", "mean"),
        donors=("donor", "nunique"), measured=("n_measured", "sum"),
        positive=("n_positive", "sum"), signal_sum=("signal_sum","sum"))
    if weighting == "cells":
        values["mean"] = values.signal_sum / values.measured
        values["fraction"] = values.positive / values.measured
    return support, values.drop(columns="signal_sum"), scope.donor.nunique()


def draw_matrix(values, labels, markers, row_text, *, title, measure, tissue, scope_note):
    """Identical axes/scales for independently selectable dotplots and heatmaps."""
    if not labels or not markers:
        return
    idx = pd.MultiIndex.from_product([labels,markers], names=["label","marker"])
    val = values.reindex(idx)
    means = val["mean"].to_numpy().reshape(len(labels),len(markers))
    fractions = val.fraction.to_numpy().reshape(means.shape)
    if measure not in ("heatmap", "dotplot"):
        raise ValueError("Choose heatmap or dotplot")
    mean = pd.DataFrame(means, index=row_text, columns=markers)
    fraction = pd.DataFrame(fractions, index=row_text, columns=markers)
    rows = pd.DataFrame({"phenotype": pd.Categorical(row_text, categories=row_text, ordered=True)},
                        index=row_text)
    summary = ad.AnnData(means, obs=rows, var=pd.DataFrame(index=markers))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e2e2e2")
    size = (max(10.2, .48 * len(markers) + 4.3), max(4.4, .34 * len(labels) + 2.3))
    params = dict(categories_order=row_text, title=title+"\n"+scope_note, figsize=size, vmin=0, vmax=1)
    if measure == "heatmap":
        plot = sc.pl.MatrixPlot(summary, markers, "phenotype", values_df=mean, **params)
        plot.style(cmap=cmap, edge_lw=.25, edge_color="white").legend(title="Mean scaled protein")
    else:
        plot = sc.pl.DotPlot(summary, markers, "phenotype", dot_color_df=mean,
                             dot_size_df=fraction.fillna(0), **params)
        plot.style(cmap=cmap, dot_min=0, dot_max=1, size_exponent=1, largest_dot=115).legend(
            colorbar_title="Mean scaled protein", size_title="Measured cells ≥ 0.10 (%)")
    plot.make_figure()
    ax = plot.get_axes()["mainplot_ax"]
    y, x = np.where(~np.isfinite(means))
    if measure == "heatmap":
        for row, column in zip(y,x):
            ax.add_patch(Rectangle((column,row),1,1,facecolor="#c6c9cc",edgecolor="white",linewidth=.3))
    else:
        ax.scatter(x+.5, y+.5, marker="x", s=12, linewidths=.7, color="#92989a")
    ax.tick_params(axis="both", labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")
    plt.show()


def show_profiles(counts,evidence,*,tissue="Pancreas",level="primary",disease=None,donor=None,
                  weighting="donors",measures=("dotplot","heatmap")):
    support, values, total = measured_profiles(counts,evidence,tissue=tissue,level=level,disease=disease,
                                               donor=donor,weighting=weighting)
    if total == 0:
        return values
    scope = f"Donor {donor}" if donor is not None else f"{disease or 'Cohort'} · {total} donors"
    unit = "Single section" if donor is not None else ("Equal donor weight" if weighting=="donors" else "Pooled cells")
    # Include registered but absent or unassessable phenotypes for the selected section.
    labels = [x for x in support.index if x != AMBIGUOUS]
    for title, rows, markers in marker_panels(tissue,level,labels):
        captions = []
        for row in rows:
            item=support.loc[row]
            n = f"{int(item.cells):,}" if pd.notna(item.cells) else "unassessable"
            nrow = f"{n} cells" if pd.notna(item.cells) else n
            if donor is None:
                nrow += f", {int(item.contributing)}/{total} donors"
            captions.append(f"{DISPLAY_NAMES.get(row,row)}  ·  {nrow}")
        for measure in measures:
            draw_matrix(values, rows, markers, captions, title=f"{tissue} · {scope} · {title}",
                        measure=measure,tissue=tissue,scope_note=f"{unit}   |   Gray {'×' if measure=='dotplot' else 'tile'} = no reportable measurement")
    return values.loc[values.index.get_level_values("label") != AMBIGUOUS]


def donor_heatmap(counts,evidence,*,tissue,level,label,markers):
    all_sections = counts.loc[(counts.tissue==tissue)&(counts.level==level)&(counts.label==label)].copy()
    if all_sections.empty:
        raise ValueError(f"No {label} label at {level} in {tissue}")
    all_sections["group_order"] = all_sections.disease.map({"ND":0,"AAb+":1,"T1D":2})
    all_sections["donor_num"] = all_sections.donor.astype(int)
    all_sections = all_sections.sort_values(["group_order","donor_num"])
    frame = evidence.loc[(evidence.tissue==tissue)&(evidence.level==level)&(evidence.label==label)]
    # This graph compares donors, so require comparable class and channel measurements.
    frame = frame.merge(all_sections.loc[all_sections.comparable,["donor"]],on="donor",validate="many_to_one")
    frame = frame.set_index(["donor","marker"]).rename(columns={"mean_scaled":"mean","positive_fraction":"fraction"})
    support = [f"{r.disease} · {r.donor}  ·  {int(r.cells):,} cells" if pd.notna(r.cells) else f"{r.disease} · {r.donor} · unassessable"
               for r in all_sections.itertuples()]
    draw_matrix(frame, all_sections.donor.tolist(),markers,support,title=f"{tissue} · {DISPLAY_NAMES.get(label,label)}",
                measure="heatmap",tissue=tissue,scope_note="One section per row   |   Gray = no comparable measurement")
    return frame[["n_measured","n_positive","mean","fraction"]]


def conflict_summary(data,counts,*,disease=None,donor=None):
    q=pd.read_parquet(Path(data)/"endocrine_conflict_qc.parquet")
    q["donor"]=q.donor.astype(str)
    if q.object_id.duplicated().any():
        raise ValueError("Repeated native cell in endocrine conflict evidence")
    q=q.loc[q.corrected_class.eq(AMBIGUOUS)].copy()
    if "final_class" in q:
        q=q.loc[q.final_class.eq(AMBIGUOUS)]
    census=counts.loc[counts.tissue.eq("Pancreas") & counts.level.eq("primary") & counts.label.eq(AMBIGUOUS)]
    if disease is not None:
        q=q.loc[q.disease.eq(disease)];census=census.loc[census.disease.eq(disease)]
    if donor is not None:
        q=q.loc[q.donor.eq(str(donor))];census=census.loc[census.donor.eq(str(donor))]
    observed=q.groupby("donor",observed=True).size().reindex(census.donor,fill_value=0)
    np.testing.assert_array_equal(observed.to_numpy(),census.cells.fillna(0).to_numpy())
    q["competing hormones"] = q.a_marker.astype("string").fillna("unknown") + " / " + q.b_marker.astype("string").fillna("unknown")
    return q.groupby(["competing hormones","decision"],observed=True).agg(
        cells=("object_id","size"),donors=("donor","nunique"),median_area_gap_pp=("area_lead_pp","median"),
        marker_A_median=("a_interior_median_scaled","median"),
        marker_B_median=("b_interior_median_scaled","median")).reset_index()


def validate_release(data):
    """Check the declared assay release and the actual marker list saved with its UMAP."""
    import json
    import h5py
    data = Path(data)
    expected = "io60-phenocycler-20261113-no-epcam"
    receipt_file = data / "typing_release.json"
    if not receipt_file.is_file():
        raise FileNotFoundError(f"PhenoCycler release receipt missing at {receipt_file}")
    receipt = json.loads(receipt_file.read_text())
    if receipt["release_id"] != expected or "EpCAM" not in receipt["typing_excluded_markers"]:
        raise ValueError(f"Select the {expected} bundle with PHENOCYCLER_DATA")
    expected_markers = list(receipt["markers_used_for_pancreas_umap"])
    with h5py.File(data/"pancreas_umap.h5ad", "r") as panel:
        actual_id = ad.io.read_elem(panel["uns/release_id"])
        actual_markers = list(ad.io.read_elem(panel["uns/umap_markers"]))
    if (actual_id != expected or not expected_markers or len(set(expected_markers)) != len(expected_markers)
            or "EpCAM" in expected_markers or expected_markers != actual_markers):
        raise ValueError("The saved pancreas UMAP and the typing release declare different inputs")
    manifest = pd.read_csv(data/"pancreas_manifest.csv", dtype={"donor":str})
    if manifest.donor.duplicated().any() or int(receipt["pancreas_cells"]) != int(manifest.cells.sum()):
        raise ValueError("The full-section census disagrees with the typing release")
    print(f"{expected} · {len(manifest)} donors · {manifest.cells.sum():,} admitted pancreas cells")
    print("Verified pancreas UMAP markers:", ", ".join(actual_markers))
    return receipt


def _endocrine_observations(file):
    """Read three hormone columns and selected categorical annotations without loading the section."""
    import h5py
    with h5py.File(file,"r") as handle:
        obs = handle["obs"]
        primary = ad.io.read_elem(obs["primary"])
        chosen = primary.isin([AMBIGUOUS,"Alpha","Beta","Delta"])
        ids = np.flatnonzero(chosen)
        out = pd.DataFrame({"label": np.asarray(primary[chosen],dtype=str)})
        for field in ("islet_route", "endocrine_subtype_status", "refine_state"):
            out[field] = np.asarray(ad.io.read_elem(obs[field])[ids], dtype=object)
        out["object_id"] = None
        ambiguous = out.label.eq(AMBIGUOUS).to_numpy()
        if ambiguous.any():
            keys = obs["object_id"]
            if isinstance(keys,h5py.Group) and "values" in keys:
                keys = keys["values"]
            if isinstance(keys,h5py.Dataset):
                out.loc[ambiguous,"object_id"] = keys.asstr()[ids[ambiguous]]
            else:
                out.loc[ambiguous,"object_id"] = np.asarray(ad.io.read_elem(keys)[ids[ambiguous]])
        var = ad.io.read_elem(handle["var"])
        for marker in HORMONES:
            col = var.index.get_loc(marker)
            scaled = np.asarray(handle["layers/typing_scaled"][:,col])[ids]
            raw = np.asarray(handle["X"][:,col])[ids]
            span = var.at[marker,"span"]
            if not np.isfinite(span) or span <= 0:
                uncapped = np.full(len(raw),np.nan)
            else:
                uncapped = (raw.astype(float) - var.at[marker,"floor"]) / span
            if not np.allclose(np.clip(uncapped,0,1),scaled,rtol=1e-5,atol=1e-6,equal_nan=True):
                raise ValueError(f"Raw and scaled {marker} differ in {file.name}")
            out[marker] = scaled
            out[marker+"_uncapped"] = uncapped
        return out


def endocrine_evidence(data,counts):
    """Summarize measured named/ambiguous cells within donor and retain true denominators."""
    data = Path(data)
    reference = pd.read_parquet(data/"endocrine_conflict_qc.parquet")
    reference = reference.loc[reference.corrected_class.eq(AMBIGUOUS)].copy()
    if "final_class" in reference:
        reference = reference.loc[reference.final_class.eq(AMBIGUOUS)]
    reference["donor"] = reference.donor.astype(str)
    if reference.object_id.duplicated().any():
        raise ValueError("Multiple conflict records for one native cell")
    manifest = pd.read_csv(data/"pancreas_manifest.csv",dtype={"donor":str})
    scored = []
    qc = []
    for section in manifest.itertuples():
        did = section.donor
        cells = _endocrine_observations(data/section.path)
        official = counts.loc[counts.tissue.eq("Pancreas") & counts.level.eq("primary") & counts.donor.eq(did)]
        for label in ["Alpha","Beta","Delta",AMBIGUOUS]:
            if cells.label.eq(label).sum() != official.loc[official.label.eq(label),"cells"].sum():
                raise ValueError(f"Endocrine denominator differs in {did} for {label}")
        unknown = cells.loc[cells.label.eq(AMBIGUOUS)]
        ref = reference.loc[reference.donor.eq(did)]
        if set(unknown.object_id) != set(ref.object_id):
            raise ValueError(f"The final ambiguous calls and native evidence have different IDs in {did}")
        unknown = unknown.merge(ref[["object_id","a_marker","b_marker"]],on="object_id",validate="one_to_one")
        for marker,label in zip(HORMONES,["Alpha","Beta","Delta"]):
            named = cells.loc[cells.label.eq(label)]
            measured = named[marker].notna()
            low = measured & named[marker].lt(.10)
            unique = named.endocrine_subtype_status.eq("unique_highest_hormone")
            direct = named.islet_route.eq("single_hormone")
            unavailable = named.refine_state.eq("unavailable_bound")
            qc.append({"donor":did,"disease":section.disease,"named_type":label,"defining_hormone":marker,
                       "named_cells":len(named),"measured":int(measured.sum()),"below_0.10":int(low.sum()),
                       "unique_hormone_winner":int(unique.sum()),"direct_Single":int(direct.sum()),
                       "no_refinement_bound":int(unavailable.sum()),
                       "below_0.10_and_no_bound":int((low & unavailable).sum()),
                       "direct_Single_below_0.10_and_no_bound":int((direct & low & unavailable).sum()),
                       "weakest_measured_signal":named[marker].min()})
        for (a,b), ambiguous in unknown.groupby(["a_marker","b_marker"],dropna=False,observed=True):
            for marker in (a,b):
                if marker not in HORMONES:
                    continue
                label = {"GCG":"Alpha","INS":"Beta","SST":"Delta"}[marker]
                named = cells.loc[cells.label.eq(label)]
                for status,part in [("Competing hormones",ambiguous),(label,named)]:
                    measured = part[marker].dropna()
                    uncapped = part[marker+"_uncapped"]
                    scored.append(dict(donor=did,disease=section.disease,pair=f"{a} / {b}",hormone=marker,
                                       population=status,cells=len(part),measured=len(measured),
                                       positive=int(measured.ge(.1).sum()),scaled_sum=measured.astype(float).sum(),
                                       uncapped_sum=uncapped.sum()))
    pairs = pd.DataFrame(scored,columns=["donor","disease","pair","hormone","population","cells","measured",
                                         "positive","scaled_sum","uncapped_sum"])
    return pairs,pd.DataFrame(qc)


def matched_hormones(pairs,*,disease=None,donor=None):
    """For each pair and contender, use the identical measured donor set for both populations."""
    use = pairs.copy()
    if disease is not None:
        use = use.loc[use.disease.eq(disease)]
    if donor is not None:
        use = use.loc[use.donor.eq(str(donor))]
    if use.empty:
        return pd.DataFrame()
    keys = ["pair","hormone","donor"]
    use["matched"] = use.groupby(keys,observed=True).measured.transform("min").gt(0)
    rows = []
    for (pair,marker,status), all_donors in use.groupby(["pair","hormone","population"],observed=True):
        both = all_donors.loc[all_donors.matched]
        measured = int(both.measured.sum())
        rows.append({"pair":pair,"hormone":marker,"population":status,
                     "donors_with_this_pair":all_donors.donor.nunique(),
                     "donors_measured_in_both":both.donor.nunique(),"matched_population_cells":int(both.cells.sum()),
                     "matched_measured":measured,"matched_positive":int(both.positive.sum()),
                     "equal_donor_mean_scaled":(both.scaled_sum/both.measured).mean(),
                     "equal_donor_positive_pct":100*(both.positive/both.measured).mean(),
                     "equal_donor_mean_uncapped":(both.uncapped_sum/both.measured).mean(),
                     "pooled_mean_scaled":both.scaled_sum.sum()/measured if measured else np.nan,
                     "pooled_positive_pct":100*both.positive.sum()/measured if measured else np.nan,
                     "pooled_mean_uncapped":both.uncapped_sum.sum()/measured if measured else np.nan})
    return pd.DataFrame(rows)


def neighbor_support(data, counts):
    """Bind section-level concordance aggregates to the admitted donor/phenotype census."""
    frame = pd.read_parquet(Path(data)/"pancreas_refinement_qc.parquet")
    frame["donor"] = frame.donor.astype(str)
    unit = ["donor","disease","celltype"]
    actual = frame.groupby(unit,observed=True,dropna=False).cells.sum().sort_index()
    census = counts.loc[counts.tissue.eq("Pancreas") & counts.level.eq("primary")].rename(columns={"label":"celltype"})
    expected = census.groupby(unit,observed=True).cells.sum().fillna(0).sort_index()
    actual, expected = actual.align(expected,fill_value=0)
    if not actual.eq(expected).all():
        raise ValueError("Refinement quality counts differ from the admitted phenotype census")
    if not (frame.concordance_measured.between(0,frame.cells).all() and
            frame.zero_concordance.between(0,frame.concordance_measured).all()):
        raise ValueError("Neighbor support counts have inconsistent denominators")
    present = frame.concordance_measured.gt(0)
    if not frame.loc[present,"mean_concordance"].between(0,1).all():
        raise ValueError("Measured neighbor agreement must range from 0 to 1")
    if not frame.low_concordance.between(frame.zero_concordance,frame.concordance_measured).all():
        raise ValueError("Low neighbor agreement must include all measured zero-agreement cells")
    frame["agreement_sum"] = (frame.mean_concordance * frame.concordance_measured).fillna(0)
    frame["changed_since_first_pass"] = frame.cells.where(frame.celltype.ne(frame.first_pass),0)
    frame["not_graph_measured"] = frame.cells-frame.concordance_measured
    return frame


def summarize_neighbors(frame, keys):
    num = ["cells","changed_since_first_pass","concordance_measured","not_graph_measured","zero_concordance","agreement_sum"]
    num.append("low_concordance")
    summary = frame.groupby(keys,observed=True,dropna=False)[num].sum()
    summary["mean_neighbor_agreement"] = summary.pop("agreement_sum") / summary.concordance_measured.replace(0,np.nan)
    summary["low_pct_of_measured"] = summary.low_concordance.mul(100)/summary.concordance_measured.replace(0,np.nan)
    summary["zero_pct_of_measured"] = summary.zero_concordance.mul(100)/summary.concordance_measured.replace(0,np.nan)
    return summary


def plot_matched_hormones(matched, *, scope, weighting="donors"):
    """Each dot is one independently matched pair/contender/arm, with exact measured support."""
    if matched.empty:
        print(f"{scope} · no retained pair")
        return
    if weighting not in ("donors","cells"):
        raise ValueError("Choose 'donors' or 'cells'")
    mean_col, pct_col = (("equal_donor_mean_scaled","equal_donor_positive_pct") if weighting=="donors"
                         else ("pooled_mean_scaled","pooled_positive_pct"))
    rows = matched.sort_values(["pair","hormone","population"],key=lambda x:
                               x.map({"Competing hormones":0}).fillna(1) if x.name=="population" else x)
    if rows.duplicated(["pair","hormone","population"]).any():
        raise ValueError("Each pair and hormone must have one record per comparison arm")
    groups = rows[["pair","hormone"]].drop_duplicates().itertuples(index=False,name=None)
    entries = []
    for pair, hormone in groups:
        panel = rows[(rows.pair==pair) & (rows.hormone==hormone)]
        if len(panel)!=2 or panel.donors_measured_in_both.nunique()!=1:
            raise ValueError(f"Comparison arms disagree for {pair}, {hormone}")
        unknown = panel[panel.population.eq("Competing hormones")].iloc[0]
        named = panel[panel.population.ne("Competing hormones")].iloc[0]
        row_label = f"{pair}  ·  {hormone}"
        entries.append({"row":row_label, "unknown":unknown[mean_col], "named":named[mean_col],
                        "unknown_fraction":unknown[pct_col]/100, "named_fraction":named[pct_col]/100,
                        "annotation":f"{named.population}   |   {int(unknown.matched_measured):,} / {int(named.matched_measured):,}   |   "
                                     f"{int(named.donors_measured_in_both)} / {int(named.donors_with_this_pair)}"})
    frame = pd.DataFrame(entries).set_index("row")
    # The two columns use each row's same hormone and identical donor roster.
    # Quantities are precomputed once for each arm. The AnnData only supplies the Scanpy grouping axes.
    columns = ["Competing pair", "Named subtype"]
    means = frame[["unknown","named"]].set_axis(columns,axis=1)
    fractions = frame[["unknown_fraction","named_fraction"]].set_axis(columns,axis=1)
    chart_labels=[f"{label}    |    {note}" for label,note in frame.annotation.items()]
    means.index = chart_labels
    fractions.index = chart_labels
    obj = ad.AnnData(means.to_numpy(), obs=pd.DataFrame(
        {"comparison":pd.Categorical(chart_labels,categories=chart_labels,ordered=True)},index=chart_labels),
        var=pd.DataFrame(index=columns))
    mode = "Equal donor weight" if weighting=="donors" else "Pooled measured cells"
    title = f"{scope} · Hormone signal in competing and named cells\n{mode} · same cytoplasmic compartment"
    plot = sc.pl.DotPlot(obj,columns,"comparison",categories_order=chart_labels,
                         dot_color_df=means,dot_size_df=fractions.fillna(0),vmin=0,vmax=1,
                         figsize=(6.0,max(3.4,.42*len(frame)+2.2)),title=title)
    plot.style(cmap="viridis",dot_min=0,dot_max=1,size_exponent=1,largest_dot=155).legend(
        colorbar_title="Mean scaled protein",size_title="Measured cells ≥ 0.10 (%)")
    plot.make_figure()
    ax=plot.get_axes()["mainplot_ax"]
    yy,xx=np.where(~np.isfinite(means.to_numpy()))
    ax.scatter(xx+.5,yy+.5,marker="x",s=18,c="#767d80",lw=.8)
    ax.tick_params(axis="both",labelsize=8)
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    ax.text(0,-.17,"Pair · hormone  |  named subtype  |  measured cells (pair / named)  |  donors (matched / available)",
            ha="center",va="top",transform=ax.transAxes,fontsize=7.5,color="#364047")
    plt.show()
