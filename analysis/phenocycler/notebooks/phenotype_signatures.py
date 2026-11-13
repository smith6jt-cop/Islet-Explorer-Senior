"""Measured marker summaries for the PhenoCycler notebook series.

Pancreas structural and endocrine markers follow typing/refine.py::OWN_MARKERS.
Immune markers follow phenotype_workflow.csv, the recorded macrophage rule and
the recovered protein profiles. pLN markers follow the independent tissue rule
table in core/panels/data/tissues/pancreatic_lymph_node/phenotype_rules.csv.
"""
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
MACROPHAGE = ["CD68", "CD163", "M2Gal3"]
PLN_MYELOID = ["CD68", "CD163", "M2Gal3", "Iba1", "CD206", "CD11b", "CD11c",
               "HLA_DR", "iNOS", "CD209", "CD14", "CD141"]
PAN_ORDER = ["Alpha", "Beta", "Delta", AMBIGUOUS, BELOW_FLOOR, "Exocrine", "Blood vessel", "Lymphatic",
             "Nerve", "Immune", "Imm-T-cell", "Imm-B-cell", "Imm-plasma-cell", "Imm-NK-cell", "Imm-neutrophil",
             "Imm-dendritic-cell", "Imm-macrophage", "T-CD4-helper", "T-regulatory", "T-CD8", "T-CD8-cytotoxic",
             "T-exhausted", "T-progenitor-memory"]
PLN_ORDER = ["Endothelial", "Blood_vessel", "HEV", "Lymphatic_vessel", "Stromal", "FRC", "Fibroblast",
             "Basement_membrane", "Pericyte", "Smooth_muscle", "Neural", "Nerve_fiber", "Schwann_cell", "Immune"]


def names_in_order(labels, tissue):
    remaining = set(labels)
    preferred = PAN_ORDER if tissue == "Pancreas" else PLN_ORDER
    unresolved = [s for s in ["Immune, class unspecified", "Unresolved (no lineage evidence)",
                              "Ambiguous candidate", "No rule match", "Unavailable candidate"] if s in remaining]
    first = [s for s in preferred if s in remaining]
    return first + sorted(remaining.difference(first, unresolved)) + unresolved


def full_panel(evidence, tissue):
    """Return every released calibrated protein except the excluded identity channel."""
    available = set(evidence.loc[evidence.tissue.eq(tissue), "marker"].dropna()) - {"EpCAM"}
    preferred = (PAN_STRUCTURE if tissue == "Pancreas" else PLN_STRUCTURE) + IMMUNE + T_CELL + MACROPHAGE + PLN_MYELOID
    first = [name for name in dict.fromkeys(preferred) if name in available]
    return first + sorted(available.difference(first))


def hierarchy_tables(counts, evidence=None, level="subtype"):
    """Add a mutually exclusive leaf level while retaining the source hierarchy for QC."""
    if level != "subtype":
        return counts, evidence

    def leaves(frame):
        keep_detail = frame.level.eq("detail") & ~(frame.tissue.eq("Pancreas") & frame.label.eq("Imm-macrophage"))
        profiles = frame.tissue.eq("Pancreas") & frame.level.eq("macrophage_profile")
        return frame.loc[keep_detail | profiles].assign(level="subtype")

    if "subtype" not in counts.level.unique():
        deep = leaves(counts)
        deep["parent_cells"] = deep.eligible_cells
        deep["fraction_parent"] = deep.fraction_tissue
        grouped = deep.groupby(["tissue", "donor"], observed=True)
        if not grouped.cells.sum().eq(grouped.eligible_cells.first()).all():
            raise ValueError("Deepest subtype counts do not partition each complete section")
        counts = pd.concat([counts, deep], ignore_index=True)
    if evidence is not None and "subtype" not in evidence.level.unique():
        evidence = pd.concat([evidence, leaves(evidence)], ignore_index=True)
    return counts, evidence


def measured_profiles(counts, evidence, *, tissue="Pancreas", level="subtype", disease=None,
                      donor=None, weighting="donors"):
    """Return phenotype support and exact donor/cell denominators for every displayed measurement."""
    if weighting not in ("donors", "cells"):
        raise ValueError("Choose weighting='donors' or 'cells'.")
    if level == "subtype":
        counts, evidence = hierarchy_tables(counts, evidence)
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


def _scale_heatmap_columns(mean):
    """Scale finite displayed group means per marker and flag values with an undefined range."""
    measured = mean.where(np.isfinite(mean))
    minimum = measured.min(axis=0)
    span = measured.max(axis=0) - minimum
    scaled = (measured - minimum).div(span.where(span.gt(0)))
    return scaled, measured.notna() & span.eq(0)


def draw_matrix(values, labels, markers, row_text, *, title, measure, tissue, scope_note, applicable_only=False):
    """Draw Scanpy profiles from full-section means and measured positivity fractions."""
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
    size = (max(6.2, .32 * len(markers) + 3.5), max(4.5, .27 * len(labels) + 2.0))
    params = dict(categories_order=row_text, title=title+"\n"+scope_note, figsize=size, vmin=0, vmax=1)
    if measure == "heatmap":
        scaled, constant = _scale_heatmap_columns(mean)
        plot = sc.pl.MatrixPlot(summary, markers, "phenotype", values_df=scaled, **params)
        plot.style(cmap=cmap, edge_lw=.25, edge_color="white").legend(title="Mean expression\nmin–max per marker")
    else:
        plot = sc.pl.DotPlot(summary, markers, "phenotype", dot_color_df=mean,
                             dot_size_df=fraction.fillna(0), standard_scale="var", **params)
        plot.style(cmap=cmap, dot_min=0, dot_max=1, size_exponent=1, largest_dot=115).legend(
            colorbar_title="Mean expression\nscaled per marker", size_title="Measured cells\n≥ 0.10 (%)")
    plot.make_figure()
    ax = plot.get_axes()["mainplot_ax"]
    missing = ~np.isfinite(means)
    if applicable_only:
        missing &= idx.isin(values.index).reshape(means.shape)
    y, x = np.where(missing)
    if measure == "heatmap":
        for row, column in zip(y,x):
            ax.add_patch(Rectangle((column,row),1,1,facecolor="#c6c9cc",edgecolor="white",linewidth=.3))
        keys = []
        for mask, face, edge, note in (
                (constant & mean.ne(0), "#fff2da", "#987348", "Constant mean"),
                (constant & mean.eq(0), "#482173", "#ddd0ed", "Constant zero")):
            if not mask.to_numpy().any():
                continue
            style = dict(facecolor=face, edgecolor=edge, hatch="//", linewidth=.3)
            for row, column in zip(*np.where(mask)):
                ax.add_patch(Rectangle((column,row),1,1,**style))
            keys.append(Rectangle((0,0),1,1,label=note,**style))
        if keys:
            plot.get_axes()["color_legend_ax"].legend(
                handles=keys, loc="lower center", bbox_to_anchor=(.5, 2.7), borderaxespad=0,
                frameon=False, fontsize=8)
    else:
        ax.scatter(x+.5, y+.5, marker="x", s=12, linewidths=.7, color="#92989a")
    ax.tick_params(axis="both", labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=90 if len(markers) > 15 else 0, ha="center")
    plt.show()


def show_profiles(counts,evidence,*,tissue="Pancreas",level="subtype",disease=None,donor=None,
                  weighting="donors",measures=("dotplot","heatmap")):
    support, values, total = measured_profiles(counts,evidence,tissue=tissue,level=level,disease=disease,
                                               donor=donor,weighting=weighting)
    if total == 0:
        return values
    scope = f"Donor {donor}" if donor is not None else f"{disease or 'Cohort'} · {total} donors"
    unit = "Single section" if donor is not None else ("Equal donor weight" if weighting=="donors" else "Pooled cells")
    rows = names_in_order(support.index.tolist(), tissue)
    markers = MACROPHAGE if tissue == "Pancreas" and level == "macrophage_profile" else full_panel(evidence, tissue)
    captions = []
    for row in rows:
        item = support.loc[row]
        n = f"{int(item.cells):,} cells" if pd.notna(item.cells) else "unassessable"
        if donor is None:
            n += f", {int(item.contributing)}/{total} donors"
        captions.append(f"{subtype_display_name(row)}  ·  {n}")
    for measure in measures:
        draw_matrix(values, rows, markers, captions, title=f"{tissue} · {scope} · {level.replace('_', ' ')}{' candidates' if tissue == 'pLN' else ''}",
                    measure=measure, tissue=tissue,
                    scope_note=f"{unit}   |   {len(markers)} proteins   |   Gray {'×' if measure=='dotplot' else 'tile'} = no reportable measurement")
    return values.loc[values.index.get_level_values("marker").isin(markers)]



def donor_heatmap(counts,evidence,*,tissue,level,label,markers):
    counts, evidence = hierarchy_tables(counts, evidence, level)
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
    profile_labels = pd.read_parquet(data / "donor_counts.parquet", columns=["tissue", "level", "label"])
    profiles = profile_labels.loc[profile_labels.tissue.eq("Pancreas") &
                                  profile_labels.level.eq("macrophage_profile"), "label"]
    profile_types = set(profiles.dropna()).difference({"Macrophage, profile unavailable"})
    if profile_types != {"Macrophage_Gal3+", "Macrophage_Gal3-"}:
        raise ValueError("Select the release with Gal3+ and Gal3- macrophage subtypes.")
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
    keys = ["pair", "hormone", "population"]
    if matched.duplicated(keys).any():
        raise ValueError("Repeated comparison arm")
    named = {"GCG": "Alpha", "INS": "Beta", "SST": "Delta"}
    for (pair, hormone), arms in matched.groupby(["pair", "hormone"], observed=True):
        if (set(arms.population) != {"Competing hormones", named[hormone]}
                or arms.donors_measured_in_both.nunique() != 1):
            raise ValueError(f"Inconsistent matched arms for {pair}, {hormone}")
    markers = [m for m in HORMONES if m in set(matched.hormone)]
    rows, captions, records = [], [], []
    for pair, frame in matched.groupby("pair", sort=True, observed=True):
        for population in ["Competing hormones", *[named[m] for m in markers]]:
            selected = frame[frame.population.eq(population)].sort_values("hormone")
            if selected.empty:
                continue
            label = f"{pair}  ·  {population}"
            rows.append(label)
            notes = {}
            for item in selected.to_dict("records"):
                m = item["hormone"]
                notes[m] = (f'{item["matched_measured"]:,} cells, {item["donors_measured_in_both"]}'
                            f'/{item["donors_with_this_pair"]} donors')
                records.append({"label": label, "marker": m, "mean": item[mean_col], "fraction": item[pct_col]/100})
            support = next(iter(notes.values())) if len(set(notes.values())) == 1 else " · ".join(
                f"{m} {value}" for m, value in notes.items())
            captions.append(f"{label}  ·  {support}")
    values = pd.DataFrame(records).set_index(["label", "marker"])
    unit = "Equal donor weight" if weighting == "donors" else "Pooled measured cells"
    draw_matrix(values, rows, markers, captions, title=f"{scope} · Matched endocrine comparisons",
                measure="dotplot", tissue="Pancreas", applicable_only=True,
                scope_note=f"{unit}  |  Blank = outside comparison; × = no matched measurement")


def cell_subtypes(obs, tissue):
    """Return the deepest recorded assignment for each cell in the tissue hierarchy."""
    labels = obs["detail"].astype("string").fillna(obs["primary"].astype("string"))
    if tissue == "Pancreas":
        profiles = obs["macrophage_profile"].astype("string")
        labels = labels.mask(labels.eq("Imm-macrophage") & profiles.notna(), profiles)
    if labels.isna().any():
        raise ValueError(f"Missing cell assignment in the {tissue} release")
    return labels.astype("category")


def subtype_display_name(label):
    """Compact labels on axes; the notebook tables retain the release identifiers."""
    if label in {"Macrophage_Gal3+", "Macrophage_Gal3-"}:
        return "Macrophage · Gal3" + ("+" if label.endswith("+") else "−")
    return DISPLAY_NAMES.get(label, label).replace("_", " ")


def phenotype_colors(counts, tissue, extra_labels=()):
    """Stable colors from the complete tissue census, shared by all donor and group views."""
    fixed = ({"Alpha":"#cf6c30", "Beta":"#287cac", "Delta":"#4c8b58", "Immune":"#965084",
              "Exocrine":"#c1b7a6", AMBIGUOUS:"#bda559", BELOW_FLOOR:"#99853f",
              "Blood vessel":"#b34b50", "Lymphatic":"#68a69a", "Nerve":"#5b558c",
              "Unresolved (no lineage evidence)":"#868686"} if tissue == "Pancreas" else
             {"Ambiguous candidate":"#294e77", "Endothelial":"#bf6a2c", "Immune":"#348472",
              "Neural":"#9062a5", "No rule match":"#ac4e62", "Smooth_muscle":"#86734c",
              "Stromal":"#608aa5", "Unavailable candidate":"#737373"})
    present = set(counts.loc[counts.tissue.eq(tissue), "label"].dropna())
    present.update(extra_labels)
    others = sorted(present.difference(fixed))
    colors = [c for c in sc.pl.palettes.default_102 if c.lower() not in {x.lower() for x in fixed.values()}]
    if len(others) > len(colors):
        raise ValueError(f"The {tissue} census exceeds the available categorical colors")
    return {**fixed, **dict(zip(others, colors))}


def phenotype_key(ax, labels, palette, *, columns=2, title=None):
    """Show the identity key in its own figure area."""
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", linestyle="", markersize=5, color=palette[label],
                      label=subtype_display_name(label)) for label in labels]
    ax.set_axis_off()
    ax.legend(handles=handles, loc="center", ncol=columns, title=title, frameon=False,
              fontsize=8, title_fontsize=9, columnspacing=1.2, handletextpad=.35)


def pln_review_labels(group, status):
    """Keep explicit provenance on each broad review or unresolved row."""
    return group.astype("string").where(~status.eq("Cell and named cluster agree"), "Review · " + group.astype("string"))


def pln_evidence_profiles(counts, evidence, review_counts, review_markers):
    """Partition strict no-matches once using native full-section review measurements."""
    base = counts.loc[counts.tissue.eq("pLN") & counts.level.eq("subtype")].copy()
    markers = evidence.loc[evidence.tissue.eq("pLN") & evidence.level.eq("subtype")]
    keys = ["tissue", "donor", "disease"]
    no_match = base.loc[base.label.eq("No rule match")]
    added = review_counts.merge(no_match[keys+["eligible_cells"]], on=keys, validate="many_to_one")
    totals = added.groupby(keys, observed=True).cells.sum().sort_index()
    expected = no_match.set_index(keys).cells.sort_index().astype("int64")
    pd.testing.assert_series_equal(totals, expected, check_dtype=False)
    added["label"] = pln_review_labels(added.review_group, added.review_status)
    added = added.assign(level="subtype", assessable=True, comparable=True,
                         parent_cells=added.eligible_cells, fraction_tissue=added.cells / added.eligible_cells)
    added["fraction_parent"] = added.fraction_tissue
    extra = review_markers.assign(label=pln_review_labels(review_markers.review_group, review_markers.review_status),
                                  level="subtype")
    return (pd.concat([base.loc[base.label.ne("No rule match")], added], ignore_index=True),
            pd.concat([markers.loc[markers.label.ne("No rule match")], extra], ignore_index=True))
