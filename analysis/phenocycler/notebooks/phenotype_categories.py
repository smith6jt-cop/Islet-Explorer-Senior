"""Identify recorded assignment QC and review categories in the two tissue workflows."""

PANCREAS_QC = frozenset({
    "Endocrine (ambiguous subtype)", "Endocrine (unresolved subtype)",
    "Immune, class unspecified", "Unresolved (no lineage evidence)",
})
PLN_QC = frozenset({
    "Ambiguous candidate", "No rule match", "Unavailable candidate", "Unknown", "Unassessable",
    "Cell and named cluster disagree", "No named immune cluster", "Outside prior cluster analysis",
    "Immune, class unspecified", "T cell, subtype unspecified",
    "Review · B/plasma", "Review · Granulocyte", "Review · Myeloid", "Review · T cell",
})
QC_LABELS = {"Pancreas": PANCREAS_QC, "pLN": PLN_QC}
# Fixed colors keep small assignment-QC populations visible on white in every scope.
QC_COLORS = {
    "Endocrine (ambiguous subtype)": "#ad6122", "Endocrine (unresolved subtype)": "#706220",
    "Immune, class unspecified": "#8b4776", "Unresolved (no lineage evidence)": "#62676c",
    "Ambiguous candidate": "#294e77", "No rule match": "#ac4e62", "Unavailable candidate": "#64696a",
    "Unknown": "#686868", "Unassessable": "#555d65", "Cell and named cluster disagree": "#ac442b",
    "No named immune cluster": "#4f3829", "Outside prior cluster analysis": "#987211",
    "Review · B/plasma": "#3a792e", "Review · Granulocyte": "#067575", "Review · Myeloid": "#225c9e",
    "Review · T cell": "#812f80", "T cell, subtype unspecified": "#456d75",
}
PANEL_NAMES = ("Cell identities", "Assignment QC and review")


def is_qc_label(label, tissue):
    """Classify the recorded label; assay availability is a separate property."""
    return label in QC_LABELS[tissue]


def category_panels(labels, tissue):
    """Partition an ordered collection without merging any recorded label."""
    labels = list(labels)
    return [(title, [label for label in labels if is_qc_label(label, tissue) == qc])
            for qc, title in enumerate(PANEL_NAMES)]


def draw_identity_map(ax, view, *, labels, palette, title, basis="umap", key="subtype", size=6):
    """Draw one group of labels with the complete view's coordinate limits."""
    import numpy as np
    import scanpy as sc

    coordinate_key = "X_umap" if basis == "umap" else basis
    coords = np.asarray(view.obsm[coordinate_key])
    selected = view[view.obs[key].isin(labels)].copy()
    if selected.n_obs:
        selected.obs[key] = selected.obs[key].cat.remove_unused_categories()
        sc.pl.embedding(selected, basis=basis, color=key, palette=palette, size=size, ax=ax,
                        legend_loc="none", frameon=False, show=False)
    else:
        ax.text(.5, .5, "No cells in this selection", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
    if len(coords):
        limits = np.vstack([np.nanmin(coords, axis=0), np.nanmax(coords, axis=0)])
        margin = np.maximum(np.diff(limits, axis=0)[0] * .025, .001)
        ax.set_xlim(limits[0, 0] - margin[0], limits[1, 0] + margin[0])
        ax.set_ylim(limits[0, 1] - margin[1], limits[1, 1] + margin[1])
    ax.set_title(f"{title} · {selected.n_obs:,} cells")
    ax.set_aspect("equal")
