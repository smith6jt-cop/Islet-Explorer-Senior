#!/usr/bin/env python3
"""Create paper-style figures from the complete native call census, read-only.

The native producers are not changed. Plotting transforms are explicit, and
unknown/blank fields are preserved rather than interpreted as measured absence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import textwrap
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
from protein_analysis_admission import DEFAULT_EXCLUSIONS, load_protein_exclusions, protein_analysis_mask
ANALOG = REPO / "results/v3/literature_figure_reproduction_20261005"
TRAJECTORY = Path("/data/io60-calibration/nd_referenced_trajectory_production_all_exclusions_20261005")
ARIAL = Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
ARIAL_BOLD = Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf")
ARIAL_ITALIC = Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf")
PAGE = (11.69, 8.27)
EXPORT_DPI = 600
GROUP_COLORS = {"ND": "#3182bd", "AAb+": "#fed976", "T1D": "#e31a1c", "unavailable": "#999999"}
# The first five fills are actual atlas Fig3E PDF vector colors. Extension is
# explicit, because our complete identity set is larger than the paper panel.
TYPE_COLORS = ["#ff7f0e", "#ff9896", "#17becf", "#279e68", "#ffbb78",
               "#9467bd", "#1f77b4", "#8c564b", "#e377c2", "#bcbd22",
               "#aec7e8", "#c5b0d5", "#98df8a", "#c49c94", "#f7b6d2", "#dbdb8d"]
PETERS_COLORS = ["#f8766d", "#e7851e", "#d09400", "#b2a100", "#89ac00", "#45b500", "#00bc51",
                 "#00c087", "#00c0b2", "#00bcd6", "#00b3f2", "#28a3ff", "#9c8dff", "#d277ff", "#f166e8", "#ff61c7", "#ff689e"]
BLANK = {"<null>", "<empty>"}
TYPING_EXCLUSIONS = {"Excluded"}


def identity_marker_matrix(profiles: pd.DataFrame) -> pd.DataFrame:
    """Summarize admitted profiles, retaining broad endocrine subtype uncertainty.

    An old QC-only aggregate cannot be filtered by exact cell keys after aggregation. Require
    regeneration of that source table, rather than dropping every uncertainty label.
    """
    if ("denominator_scope" not in profiles
            or not profiles.denominator_scope.eq("typing_analysis_eligible_cells").all()):
        raise ValueError("Regenerate marker profiles with prepare_native_typing_display.py and the mandatory exclusion registry")
    admitted = profiles.loc[~profiles.cell_class.isin(TYPING_EXCLUSIONS)]
    return admitted.groupby(["marker", "cell_class"]).mean_calibrated_marker_signal.median().unstack()


def norm_donor(value: object) -> str:
    s = str(value)
    return str(int(s)) if s.isdigit() else s


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def configure() -> None:
    for path in (ARIAL, ARIAL_BOLD, ARIAL_ITALIC):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
                         "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "axes.linewidth": .65, "lines.linewidth": .8,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "figure.facecolor": "white"})


def clean_axes(ax: plt.Axes, grid: bool = False) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(width=.65, length=3, pad=2)
    if grid:
        ax.grid(axis="y", color="#dddddd", lw=.45)
        ax.set_axisbelow(True)


def panel(ax: plt.Axes, letter: str) -> None:
    rect = ax.get_position()
    x = (max(rect.x0 - .085 * rect.width, .035) - rect.x0) / rect.width
    ax.text(x, 1.025, letter + ")", transform=ax.transAxes,
            fontfamily="Arial", fontsize=14, weight="bold", va="bottom", ha="left")


def display_label(label: str) -> str:
    if label in BLANK:
        return "Native field blank"
    if label == "Endocrine (unresolved subtype)":
        return "Endocrine\nunresolved"
    if label == "Endocrine (ambiguous subtype)":
        return "Endocrine\nambiguous"
    if label == "Unresolved (no lineage evidence)":
        return "Unresolved lineage"
    return label.replace("_", " ")


def composition_axis_groups(manifest: dict, component: int) -> str:
    values = manifest["mfa_group_contributions_by_component"]
    groups = sorted(values, key=lambda group: (-values[group][component], group))[:3]
    return " + ".join(groups)


class Builder:
    def __init__(self, out: Path) -> None:
        self.out = out
        self.figdir = out / "figures"
        self.figdir.mkdir(parents=True, exist_ok=True)
        self.tmp = out / "tmp/pdfs"
        self.tmp.mkdir(parents=True, exist_ok=True)
        self.counts = pd.read_csv(out / "typing_audit/donor_counts.csv", dtype={"donor_id": str, "section_id": str})
        analysis_counts = out / "typing_audit/protein_identity_analysis_counts.csv"
        if analysis_counts.exists():
            self.counts = pd.concat([self.counts, pd.read_csv(analysis_counts, dtype={"donor_id": str, "section_id": str})], ignore_index=True)
        self.coverage = pd.read_csv(out / "typing_audit/coverage.csv", dtype={"donor_id": str})
        self.disease = self.coverage.dropna(subset=["disease_status"]).drop_duplicates("donor_id").set_index("donor_id").disease_status.to_dict()
        lut = pd.read_csv(out / "style_audit/bodenmiller_fig1i_native_colorbar_lut.csv")
        self.heat_cmap = colors.ListedColormap(lut[["red", "green", "blue"]].to_numpy() / 255., name="Bodenmiller_Fig1I")
        self.heat_cmap.set_bad("#d9d9d9")
        identity_fields = {"cell_class", "refined_class", "lineage", "immune_class", "immune_subtype", "t_nk_class",
                           "native_broad_identity_display", "native_tnk_identity_display"}
        identities = sorted(self.counts.loc[self.counts.label_field.isin(identity_fields), "cell_identity"].unique())
        palette = TYPE_COLORS[:5] + PETERS_COLORS + TYPE_COLORS[5:]
        self.identity_colors = {x: palette[i % len(palette)] for i, x in enumerate(identities)}
        self.identity_colors.update({"Exocrine": "#bdbdbd", "Non_immune": "#bdbdbd", "Excluded": "#555555",
                                     "Alpha": "#ff7f0e", "Beta": "#279e68", "Delta": "#17becf", "Immune": "#9467bd"})
        for identity in identities:
            if identity in BLANK or any(word in identity.lower() for word in ("uncalled", "unavailable", "unresolved", "ambiguous")):
                self.identity_colors[identity] = "#9e9e9e"
        self.figures: list[dict] = []

    def donor_sort(self, donor: str) -> tuple:
        order = {"ND": 0, "AAb+": 1, "T1D": 2}
        return order.get(self.disease.get(donor), 3), int(donor) if donor.isdigit() else donor

    def matrix(self, lane: str, field: str, scope: str, by_section: bool = False) -> tuple[pd.DataFrame, pd.Series]:
        # Native census is retained as an audit. Every current pancreas
        # analytical panel uses the exact-key retained-cell population.
        if lane == "pancreas_protein" and scope == "qc_analysis_eligible_cells":
            scope = "typing_analysis_eligible_cells"
        d = self.counts.loc[self.counts.lane.eq(lane) & self.counts.label_field.eq(field)
                            & self.counts.denominator_scope.eq(scope)].copy()
        if d.empty:
            raise ValueError(f"missing census: {lane}, {field}, {scope}")
        if by_section:
            d["unit"] = d.donor_id + " | " + d.section_id
        else:
            d["unit"] = d.donor_id
        unit_order = sorted(d.unit.unique(), key=lambda x: (self.donor_sort(x.split(" | ")[0]), x))
        count = d.groupby(["unit", "cell_identity"], observed=True).n_cells.sum().unstack(fill_value=0)
        # One source population denominator per biological section. A donor
        # with two distinct sections is counted twice, never as a mean of means.
        den = d.drop_duplicates(["donor_id", "section_id"])[["unit", "denominator"]].groupby("unit").denominator.sum()
        return count.reindex(unit_order).astype(float), den.reindex(unit_order).astype(float)

    def unit_label(self, unit: str, by_section: bool = False) -> str:
        donor = unit.split(" | ")[0]
        label = donor + " | " + str(self.disease.get(donor, "NA"))
        if by_section:
            section = unit.split(" | ", 1)[1]
            # Native full section IDs remain in source tables. Number repeated
            # donor sections for compact labels; no block identity is invented.
            sections = sorted(self.counts.loc[self.counts.donor_id.eq(donor), "section_id"].dropna().unique())
            if "/" in section:
                short = section.rsplit("/", 1)[-1]
            elif len(section) <= 22:
                short = section
            else:
                short = "section " + str(sections.index(section) + 1)
            label += " | " + short
        return label

    def compact_unit_labels(self, units: list[str]) -> tuple[list[str], list[str]]:
        number: dict[str, int] = {}
        full, short = [], []
        for unit in units:
            donor = unit.split(" | ")[0]
            number[donor] = number.get(donor, 0) + 1
            if " | " in unit:
                full.append(self.unit_label(donor) + " | section " + str(number[donor]))
                short.append(donor + ":" + str(number[donor]))
            else:
                full.append(self.unit_label(donor))
                short.append(donor)
        return full, short

    def header(self, fig: plt.Figure, title: str, subtitle: str = "") -> None:
        fig.text(.055, .975, title, fontsize=14, weight="bold", va="top")
        if subtitle:
            fig.text(.055, .942, subtitle, fontsize=8.5, va="top", color="#444444")

    def finish(self, fig: plt.Figure, name: str, title: str, caption: str, style: str, fields: list[str]) -> None:
        lines = textwrap.wrap(caption, width=175)
        fig.text(.055, .03, "\n".join(lines), fontsize=7.2, va="bottom", linespacing=1.35)
        png = self.figdir / (name + ".png")
        svg = self.figdir / (name + ".svg")
        pdf = self.tmp / (name + ".pdf")
        fig.savefig(png, dpi=EXPORT_DPI)
        # Text/axes remain vectors in PDF/SVG. Explicit DPI is also essential
        # for the rasterized million-cell scatter layers in those formats.
        fig.savefig(svg, dpi=EXPORT_DPI)
        fig.savefig(pdf, dpi=EXPORT_DPI)
        plt.close(fig)
        self.figures.append({"name": name, "title": title, "png": str(png.resolve()), "svg": str(svg.resolve()),
                             "pdf_page": str(pdf.resolve()), "caption": caption, "style_reference": style,
                             "native_fields": fields, "png_dpi": EXPORT_DPI,
                             "rasterized_point_layer_dpi_pdf_svg": EXPORT_DPI,
                             "text_and_axes_pdf_svg": "vector"})
        print(name, flush=True)

    def stack(self, ax: plt.Axes, count: pd.DataFrame, den: pd.Series, by_section: bool = False,
              exocrine_omitted: bool = False, identity_analysis: bool = False) -> None:
        classes = sorted(count.columns, key=lambda x: (x in BLANK or "Uncalled" in x or "Unavailable" in x, x))
        y = np.arange(len(count))
        left = np.zeros(len(count))
        for i, identity in enumerate(classes):
            value = count[identity].to_numpy() / den.to_numpy() * 100
            color = "#bdbdbd" if identity in BLANK or "Uncalled" in identity or "Unavailable" in identity or "Unresolved" in identity else self.identity_colors[identity]
            ax.barh(y, value, left=left, color=color, height=.76, linewidth=0, label=display_label(identity))
            left += value
        full_labels, _ = self.compact_unit_labels(list(count.index))
        ax.set_yticks(y, full_labels, fontsize=6.5 if by_section else 7)
        ax.invert_yaxis()
        ax.set_ylim(len(count) - .5, -.5)
        limit = max(5, np.ceil(float(left.max()) / 5) * 5) if exocrine_omitted else 100
        ax.set_xlim(0, limit)
        if exocrine_omitted:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4, integer=True))
        else:
            ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("Assigned cells (% of typing-analysis eligible cells)" if identity_analysis
                      else "Assigned cells (% of original source population)")
        clean_axes(ax)
        ax.legend(frameon=False, ncols=2, fontsize=6.7, loc="upper left", bbox_to_anchor=(0, -.10),
                  handlelength=1.0, handletextpad=.4, columnspacing=1.0, borderaxespad=0)

    def heat(self, fig: plt.Figure, rect: list[float], count: pd.DataFrame, den: pd.Series, title: str,
             letter: str, marker_values: bool = False, show_status_strip: bool = True) -> None:
        observed = [c for c in count.columns if c not in BLANK]
        fraction = count[observed].div(den, axis=0) if not marker_values else count
        if marker_values:
            values = fraction.T
        else:
            # Blank-only producer fields mean an unavailable call space. The
            # observed zero assigned classes in other donors remain real counts.
            available = count[observed].sum(axis=1).gt(0)
            values = fraction.T.copy()
            values.loc[:, ~available] = np.nan
        mu, sd = values.mean(axis=1), values.std(axis=1, ddof=0)
        z = values.sub(mu, axis=0).div(sd.replace(0, np.nan), axis=0)
        z.loc[sd.eq(0)] = values.loc[sd.eq(0)].where(values.loc[sd.eq(0)].isna(), 0)
        if len(z) > 1:
            # Dendrogram is a display-only clustering of the measured profile.
            # Fill centered missing display entries only for linkage; heatmap
            # retains NaN. No fitted/analytic table is imputed or changed.
            tree = linkage(z.fillna(0).to_numpy(), method="average", metric="euclidean", optimal_ordering=True)
            order = dendrogram(tree, no_plot=True)["leaves"]
            z = z.iloc[order]
        else:
            tree = None
        x, y, w, h = rect
        ax = fig.add_axes([x, y, w, h])
        image = ax.imshow(z.to_numpy(), cmap=self.heat_cmap, vmin=-2, vmax=2, aspect="auto", interpolation="nearest")
        ax.set_yticks(np.arange(len(z)), [display_label(v) for v in z.index], fontsize=6.8)
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", length=0, pad=4)
        _, short_labels = self.compact_unit_labels(list(z.columns))
        ax.set_xticks(np.arange(len(z.columns)), short_labels, rotation=90, fontsize=6.5)
        ax.tick_params(axis="x", length=0, pad=3)
        ax.set_title(title, fontsize=9, loc="left", pad=20)
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel(ax, letter)
        if show_status_strip:
            strip = fig.add_axes([x, y + h + .006, w, .012])
            strip.imshow(np.array([[colors.to_rgb(GROUP_COLORS.get(self.disease.get(v.split(" | ")[0]), "#999999")) for v in z.columns]]), aspect="auto")
            strip.set_axis_off()
        if tree is not None:
            dend = fig.add_axes([x - .043, y, .038, h])
            dendrogram(tree, orientation="left", no_labels=True, color_threshold=0, above_threshold_color="#777777", ax=dend)
            dend.invert_yaxis()
            dend.set_axis_off()
            for collection in dend.collections:
                collection.set_linewidth(.45)
        cbax = fig.add_axes([x + .07, y - .09, w * .55, .013])
        cb = fig.colorbar(image, cax=cbax, orientation="horizontal", ticks=[-2, 0, 2])
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=2, width=.5, labelsize=6.5)
        cb.set_label("Within-row standardized " + ("marker signal" if marker_values else "assigned-cell fraction"), fontsize=7)

    def composition(self, lane: str, field: str, scope: str, name: str, title: str,
                    subtitle: str, caveat: str = "", by_section: bool = False) -> None:
        count, den = self.matrix(lane, field, scope, by_section)
        identity_analysis = scope == "typing_analysis_eligible_cells"
        if identity_analysis:
            audit = self.counts.loc[self.counts.lane.eq(lane) & self.counts.label_field.eq(field)
                                   & self.counts.denominator_scope.eq("qc_analysis_eligible_cells")]
            native_count = audit.groupby(["donor_id", "cell_identity"]).n_cells.sum().unstack(fill_value=0)
            native_den = audit.drop_duplicates(["donor_id", "section_id"]).groupby("donor_id").denominator.sum()
            native_count.assign(denominator=native_den).to_csv(self.out / "source_tables" / (name + "_native_complete_counts.csv"))
        # Native complete counts are additionally exported in wide form.
        count.assign(denominator=den).to_csv(self.out / "source_tables" / (name + "_complete_counts.csv"))
        exocrine_omitted = lane.startswith("pancreas") and "Exocrine" in count.columns
        if exocrine_omitted:
            count = count.drop(columns="Exocrine")
        fig = plt.figure(figsize=PAGE)
        self.header(fig, title, subtitle)
        ax = fig.add_axes([.105 if by_section else .075, .30, .32, .58])
        self.stack(ax, count, den, by_section, exocrine_omitted, identity_analysis)
        panel(ax, "A")
        ax.set_title("Nonexocrine identities" if exocrine_omitted else "Complete source-population composition", loc="left", pad=15)
        self.heat(fig, [.49, .32, .34, .55], count, den, "Relative differences between donors", "B")
        source = "Atlas Fig3E horizontal composition; Bodenmiller Fig1I heatmap/Arial"
        caption = (f"{len(count)} donor/section populations; {int(den.sum()):,} cells in the displayed denominator. "
                   "A uses absolute assigned-cell percentages. B centers and scales each native class across the displayed populations: "
                   "red/blue indicates relatively more/fewer assigned cells, not marker positivity or independent validation. "
                   "Gray indicates a blank-only native field. Top strip: ND blue, AAb+ yellow, T1D red. Missing class calls are retained; no redistribution. "
                   + ("Typing-analysis eligibility excludes 11,459 native annotation-mismatch cells and 25 uncertain endocrine Singles. Exclusion status is not a biological type; original source counts remain in the native table/audit. " if identity_analysis else "")
                   + (("Exocrine is omitted for clarity; fractions use the full typing-analysis denominator, including Exocrine. " if identity_analysis else "Exocrine is omitted from both displays for clarity; percentages retain the original all-cell denominator. Its full counts remain in the table and census. ") if exocrine_omitted else "")
                   + ("Section suffixes enumerate the same populations left-to-right in B as top-to-bottom in A; full annotation/section IDs are in the complete-count table. " if by_section else "") + caveat)
        self.finish(fig, name, title, caption, source, [lane + "." + field, scope])

    def refinement(self) -> None:
        analysis_path = self.out / "typing_audit/pancreas_analysis_primary_refined_transitions.csv"
        d = pd.read_csv(analysis_path if analysis_path.exists() else self.out / "typing_audit/pancreas_primary_refined_transitions.csv", dtype={"donor_id": str})
        changed = d.loc[d.primary_identity.ne(d.refined_identity)]
        matrix = d.groupby(["primary_identity", "refined_identity"]).n_cells.sum().unstack(fill_value=0)
        matrix.to_csv(self.out / "source_tables/pancreas_primary_refined_transition_matrix.csv")
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Pancreas protein: refinement and immune resolution", "Refined labels remain beside the primary calls; primary identities are preserved.")
        ax = fig.add_axes([.085, .33, .43, .53])
        im = ax.imshow(np.log10(1 + matrix.to_numpy()), cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(matrix.columns)), [display_label(x) for x in matrix.columns], rotation=75, ha="right", fontsize=6.3)
        ax.set_yticks(range(len(matrix)), [display_label(x) for x in matrix.index], fontsize=7)
        ax.set_xlabel("Refined identity")
        ax.set_ylabel("Primary identity")
        panel(ax, "A")
        ax.set_title("Analysis-eligible primary-to-refined transitions", loc="left", pad=12)
        cbax = fig.add_axes([.10, .15, .35, .014])
        cb = fig.colorbar(im, cax=cbax, orientation="horizontal")
        cb.set_label("log10(1 + cells); counts retained in source table", fontsize=7)
        cb.outline.set_visible(False)
        rate = changed.groupby("donor_id").n_cells.sum().reindex(sorted(d.donor_id.unique(), key=self.donor_sort), fill_value=0)
        den = d.drop_duplicates("donor_id").set_index("donor_id").denominator.reindex(rate.index)
        bx = fig.add_axes([.64, .30, .29, .56])
        y = np.arange(len(rate))
        bx.scatter(rate / den * 100, y, s=15, c=[GROUP_COLORS.get(self.disease.get(x), "#999999") for x in rate.index], edgecolors="black", linewidths=.25)
        bx.set_yticks(y, [self.unit_label(x) for x in rate.index], fontsize=7)
        bx.invert_yaxis()
        bx.set_xlabel("Cells whose identity changed (%)")
        bx.set_title("Donor-level effect of refinement", loc="left", pad=12)
        clean_axes(bx)
        panel(bx, "B")
        self.finish(fig, "figure02_primary_refined_transitions", "Primary-to-refined transitions",
                    f"{int(d.n_cells.sum()):,} typing-analysis-eligible cells are reconciled exactly. {int(changed.n_cells.sum()):,} change identity. "
                    "11,459 native typing exclusions and 25 uncertain endocrine Singles are omitted from the biological transition display and retained in the audit. Other unresolved refined identities remain explicit. This chart describes producer behavior; it is not an accuracy score.",
                    "Bodenmiller Arial/panel-letter convention and compact heatmap; Peters donor-point style", ["cell_class", "refined_class"])

    def immune(self) -> None:
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Pancreas protein: all native immune classes and subtypes", "Primary Immune is available for all 25 donors; these finer native fields are available for 21.")
        for i, field in enumerate(["immune_class", "immune_subtype"]):
            count, den = self.matrix("pancreas_protein", field, "qc_analysis_eligible_cells")
            self.heat(fig, [.10 + .46 * i, .30, .28, .56], count, den, "Native " + field.replace("_", " "), chr(65 + i))
            count.assign(denominator=den).to_csv(self.out / "source_tables" / (field + "_all_donor_counts.csv"))
        self.finish(fig, "figure03_pancreas_immune_resolution", "Native pancreas immune resolution",
                    "Fractions use all retained pancreas cells after exclusion of the 11,484 quarantined cells, including retained cells without fine subtypes. Donors 112, 6347, 6533 and 6580 have blank fine-typing fields "
                    "and are gray, despite valid primary Immune calls. Figure 15 adds their separately qualified typing using the measured panel. Blank/uncalled values remain in census tables. Marker-defined state names are native rule outputs, "
                    "not independently established complete exhaustion/activation phenotypes.", "Bodenmiller Fig1I native LUT, Arial and compact donor strips", ["immune_class", "immune_subtype"])

    def pln_subtypes(self) -> None:
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "pLN protein: additional native subtype fields", "The producer's macrophage/CD4/CD8 field remains distinct from its T/NK refinement.")
        for i, (lane, label) in enumerate([("dedicated_pln_protein", "Dedicated pLN"), ("embedded_pln_protein", "Embedded pLN")]):
            count, den = self.matrix(lane, "immune_subtype", "qc_analysis_eligible_cells", by_section=lane.startswith("embedded"))
            self.heat(fig, [.10 + .46 * i, .31, .28, .53], count, den, label + " native immune subtype", chr(65 + i))
            count.assign(denominator=den).to_csv(self.out / "source_tables" / (lane + "_immune_subtype_counts.csv"))
        self.finish(fig, "figure06c_pln_protein_additional_subtypes", "Additional pLN native subtype fields",
                    "Every emitted nonblank immune_subtype value is shown. Fractions use all QC-eligible cells; blank fields remain in the source census and are not interpreted as absence. "
                    "The older immune_subtype field and the newer T/NK splitting field have separate definitions and are not merged into a new hierarchy. "
                    "pLN broad Myeloid/B_Plasma/Granulocyte resolution differs from pancreas macrophage/DC/B/plasma/neutrophil resolution.",
                    "Bodenmiller compact heatmap/Arial", ["immune_subtype", "dedicated and embedded pLN"])

    def markers(self) -> None:
        d = pd.read_csv(self.out / "source_tables/pancreas_native_identity_marker_profiles.csv", dtype={"donor_id": str})
        mean = identity_marker_matrix(d)
        mean.to_csv(self.out / "source_tables/pancreas_donor_balanced_identity_marker_matrix.csv")
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Pancreas protein: measured markers behind native identities", "Each entry is the donor median of the class-average released calibrated marker signal.")
        self.heat(fig, [.14, .30, .56, .53], mean.T, pd.Series(1, index=mean.columns), "Native identity-marker profiles", "A", marker_values=True, show_status_strip=False)
        # Replace generic donor tick labels and strip, since columns are classes.
        ax = fig.axes[0]
        short = {"Endocrine (ambiguous subtype)": "Endocrine\nambiguous", "Endocrine (unresolved subtype)": "Endocrine\nunresolved"}
        ax.set_xticklabels([short.get(x, display_label(x)) for x in mean.columns], rotation=45, ha="right", fontsize=7)
        self.finish(fig, "figure07_identity_marker_profiles", "Native identity-marker profiles",
                    "Marker rows are standardized across native classes solely for display. Hormone fragment masking is inherited from the released typing inputs; "
                    "missing measurements are excluded from summaries, never converted to zeros. This descriptive profile reuses classification markers and therefore "
                    "cannot independently validate typing. Ordinary endocrine uncertainty stays in the admitted broad endocrine population; it does not establish an Alpha, Beta or Delta subtype. Protein supports broad Exocrine only; acinar/ductal classes are not invented.",
                    "Bodenmiller Fig1I exact native colorbar samples and Arial", ["cell_class", "typing_scaled_* measurements"])

    def spatial(self) -> None:
        d = pd.read_parquet(self.out / "source_tables/pancreas_native_centroid_display_sample.parquet")
        sampling_path = self.out / "source_tables/native_display_sampling_manifest.json"
        sampling = json.loads(sampling_path.read_text()) if sampling_path.is_file() else {}
        exclusion_registry = Path(sampling.get("mandatory_exclusion_registry", str(DEFAULT_EXCLUSIONS)))
        d = pd.concat([
            part.loc[protein_analysis_mask(part, str(donor), load_protein_exclusions(exclusion_registry, str(donor)))]
            for donor, part in d.groupby("donor_id", sort=False)], ignore_index=True)
        donors = sorted(d.donor_id.unique(), key=self.donor_sort)
        identities = sorted(set(d.cell_class.dropna().unique()) - TYPING_EXCLUSIONS)
        palette = {x: self.identity_colors[x] for x in identities}
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Pancreas protein: native tissue cell-typing maps", "25 donor sections, uniform sampling for display; complete counts are used throughout the atlas.")
        gs = fig.add_gridspec(5, 5, left=.045, right=.96, bottom=.17, top=.88, wspace=.10, hspace=.22)
        for ax, donor in zip([fig.add_subplot(g) for g in gs], donors):
            part = d.loc[d.donor_id.eq(donor)]
            # Broad gray Exocrine first; other native classes remain visible.
            for identity in sorted(identities, key=lambda x: x != "Exocrine"):
                cells = part.loc[part.cell_class.eq(identity)]
                ax.scatter(cells.X_centroid, cells.Y_centroid, s=.08, c=palette[identity], linewidths=0, rasterized=True,
                           alpha=.10 if identity == "Exocrine" else 1.)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(float(part.X_centroid.min()), float(part.X_centroid.max()))
            ax.set_ylim(float(part.Y_centroid.min()), float(part.Y_centroid.max()))
            ax.invert_yaxis()
            ax.set_axis_off()
            ax.set_title(self.unit_label(donor), fontsize=7.5, loc="left", pad=2)
        fig.legend(handles=[Patch(facecolor=palette[x], label=display_label(x)) for x in identities], frameon=False,
                   ncols=6, loc="lower center", bbox_to_anchor=(.5, .085), fontsize=7, handlelength=1)
        self.finish(fig, "figure08_all_pancreas_native_maps", "All pancreas native maps",
                    "Up to 100,000 cells per donor sampled uniformly without replacement from the exact-key admitted population, fixed seeds. These are actual native tissue centroids, not a UMAP "
                    "or spatial registration between sections. Panels fit their own tissue extents; no common physical scale is asserted. Sampling is only for rendering and can hide rare cells. "
                    "Exocrine is drawn first as faint context. Ordinary endocrine uncertainty remains visible with subtype unavailable; actual typing exclusions are omitted from points/legends. Source coordinates and exact counts remain in the audit and native Parquet partitions.", "Atlas compact cell-identity maps; certified Arial", ["X_centroid", "Y_centroid", "cell_class"])

    def objects(self) -> None:
        d = pd.read_csv(ANALOG / "endocrine_object_class_donor_summary.csv", dtype={"donor_id": str})
        fig, axes = plt.subplots(2, 3, figsize=PAGE)
        fig.subplots_adjust(left=.075, right=.97, top=.85, bottom=.21, wspace=.32, hspace=.48)
        self.header(fig, "Endocrine objects: Single, Cluster and Islet", "Objects with retained cells are included; there is no ten-cell floor and no expansion ring.")
        rng = np.random.default_rng(20261005)
        for col, cls in enumerate(["Single", "Cluster", "Islet"]):
            part = d.loc[d.object_class.eq(cls)]
            for row, (field, label) in enumerate([("object_fraction", "Fraction of donor endocrine objects (%)"), ("median_object_area_um2", "Donor median object area (um2)")]):
                ax = axes[row, col]
                for i, status in enumerate(["ND", "AAb+", "T1D"]):
                    values = part.loc[part.donor_status.eq(status), field].dropna().to_numpy()
                    if row == 0:
                        values = values * 100
                    jitter = rng.uniform(-.18, .18, len(values))
                    ax.scatter(i + jitter, values, s=22, c=GROUP_COLORS[status], alpha=.65, edgecolors="none", zorder=3)
                    if len(values):
                        ax.errorbar(i, np.mean(values), yerr=np.std(values, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0,
                                    fmt="_", ms=15, color="black", lw=.9, capsize=3, zorder=4)
                ax.set_xticks(range(3), ["ND", "AAb+", "T1D"], rotation=45, ha="right", weight="bold")
                ax.set_ylabel(label, weight="bold")
                ax.set_title(cls, weight="bold", fontsize=9)
                clean_axes(ax)
                panel(ax, chr(65 + row * 3 + col))
        self.finish(fig, "figure09_original_object_classes", "Original endocrine object classes",
                    "Points are independent donors; black lines show mean +/- SEM across donors. Upper panels count original objects, not cells; lower panels use each donor's median "
                    "object area. Original 47,352 object records/counts/geometry are retained in the inventory; objects containing only quarantined cells are excluded from the analysis. Drotar Fig3 donor-point/mean-SEM layout is visually recreated; its figure interior is rasterized, "
                    "so its original font cannot be certified. The recreation uses exact installed Arial, not an asserted recovered Drotar font.", "Drotar Fig3 visual layout; Arial visual recreation", ["original_object_class", "object_count", "object_area_um2"])

    def paired(self) -> None:
        availability = pd.read_csv(ANALOG / "same_slide_pheno_xenium_availability.csv", dtype=str)
        paired = pd.read_csv(ANALOG / "donor_matched_section_composition.csv", dtype={"donor_id": str})
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Same-slide availability and verified sequential comparisons", "PhenoCycler + Xenium same-slide, sequential sections, source blocks and post-Xenium IF are separate relationships.")
        ax = fig.add_axes([.065, .59, .88, .28])
        ax.set_axis_off()
        rows = []
        for r in availability.itertuples():
            rows.append([r.donor_id, r.shared_slide_serial.zfill(7), r.user_confirmed_acquisition_status.replace("_", " "),
                         r.current_xenium_call_serial.zfill(7), "Not in current call release"])
        table = ax.table(cellText=rows, colLabels=["Donor", "Same-slide serial", "Acquisition status", "Available RNA serial", "Exact same-slide RNA calls"],
                         colWidths=[.08, .16, .18, .17, .33], cellLoc="left", loc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.65)
        for (row, col), cell in table.get_celld().items():
            cell.set_linewidth(.4)
            cell.set_edgecolor("#dddddd")
            if row == 0:
                cell.set_text_props(weight="bold")
                cell.set_facecolor("#eeeeee")
        panel(ax, "A")
        for i, phenotype in enumerate(["Endocrine", "Immune"]):
            bx = fig.add_axes([.16 + i * .43, .18, .28, .28])
            part = paired.loc[paired.phenotype.eq(phenotype)]
            for r in part.itertuples():
                bx.plot([0, 1], [r.protein_fraction * 100, r.rna_fraction * 100], color=GROUP_COLORS[r.disease_status], lw=.8)
                bx.scatter([0, 1], [r.protein_fraction * 100, r.rna_fraction * 100], c=GROUP_COLORS[r.disease_status], s=27, edgecolors="black", linewidths=.25)
                bx.annotate(r.donor_id, (1, r.rna_fraction * 100), xytext=(6, 0), textcoords="offset points", va="center", fontsize=7)
            bx.set_xticks([0, 1], ["Protein section", "RNA section"])
            bx.set_ylabel("Fraction of named cells (%)")
            bx.set_title(phenotype, weight="bold")
            bx.set_xlim(-.2, 1.35)
            clean_axes(bx)
            panel(bx, chr(66 + i))
        self.finish(fig, "figure10_same_slide_and_sequential", "Same-slide and sequential scope",
                    "Only donors 112 and 6533 have explicitly verified sequential pairs in the current metadata. Each line is one donor's two sections, not matched cells. "
                    "Different cell definitions/panels and section sampling can change fractions; two donor pairs support descriptive comparisons only. "
                    "Completed same-slide acquisition does not mean its exact RNA calls are available locally. Same-block/different-block nodes remain metadata, not an inferred spatial distance.",
                    "Richardson paired donor-point/Arial convention; compact native provenance table", ["confirmed sequential section fractions", "exact same-slide serial coverage"])

    def shared(self) -> None:
        d = pd.read_csv(ANALOG / "shared_marker_donor_profiles.csv", dtype={"donor_id": str})
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Shared-panel measurements across every available donor", "12 measured overlapping markers; protein and RNA are shown on separate within-marker scales.")
        for i, (field, label) in enumerate([("protein_median_raw_signal", "Protein: raw marker intensity"), ("rna_mean_lognorm_all_cells", "RNA: log-normalized gene signal")]):
            matrix = d.pivot(index="donor_id", columns="protein_marker", values=field).reindex(sorted(d.donor_id.unique(), key=self.donor_sort))
            self.heat(fig, [.10 + .46 * i, .29, .28, .57], matrix, pd.Series(1, index=matrix.index), label, chr(65 + i), marker_values=True)
        self.finish(fig, "figure11_shared_panel_marker_measurements", "Shared-panel marker measurements",
                    "All 25 currently eligible released donors are shown; 6457 is excluded. Each modality is centered/scaled within each marker across donors for display only. "
                    "Red/blue compares donors within a row, never raw protein to RNA units. These are matched-donor profiles, not same-cell concordance. Physical source block IDs remain "
                    "unresolved where metadata does not confirm them. Exact availability, source acquisition and fragment masking are in shared_marker_donor_profiles.csv.",
                    "Bodenmiller Fig1I native LUT/Arial compact heatmap", ["protein_median_raw_signal", "rna_mean_lognorm_all_cells"])

    def restain(self) -> None:
        d = pd.read_csv(ANALOG / "post_xenium_if_donor_metrics.csv", dtype={"donor_id": str})
        fig, axes = plt.subplots(1, 2, figsize=PAGE, sharey=True)
        fig.subplots_adjust(left=.11, right=.95, top=.84, bottom=.29, wspace=.17)
        self.header(fig, "Post-Xenium IF: hormone signal and broad RNA endocrine identity", "A separate restain assay with 16 verified registered sections; this is not a PhenoCycler/Xenium same-slide release.")
        donors = sorted(d.donor_id.unique(), key=self.donor_sort)
        ypos = {v: i for i, v in enumerate(donors)}
        for i, marker in enumerate(["INS", "GCG"]):
            ax = axes[i]
            part = d.loc[d.marker.eq(marker)]
            for r in part.itertuples():
                ax.scatter(r.protein_vs_rna_endocrine_identity_auc, ypos[r.donor_id], s=25,
                           c=GROUP_COLORS[r.disease_status], edgecolors="black", linewidths=.25)
            ax.axvline(.5, color="#999999", ls="--", lw=.6)
            ax.set_xlim(0, 1)
            ax.set_yticks(range(len(donors)), [self.unit_label(x) for x in donors])
            ax.set_ylim(len(donors) - .5, -.5)
            ax.set_title(marker + " protein", loc="left", weight="bold")
            ax.set_xlabel("Separation of RNA endocrine cells (AUC)")
            clean_axes(ax)
            panel(ax, chr(65 + i))
        self.finish(fig, "figure12_post_xenium_if_registered_evidence", "Post-Xenium IF endocrine evidence",
                    "Each point is one donor. AUC is the probability that a randomly chosen RNA-called endocrine cell has higher hormone signal than a named nonendocrine cell: "
                    "0.5 means no separation; below 0.5 means the direction is reversed. Panels do not measure INS/GCG RNA, so this is broad-lineage evidence, not same-hormone validation. "
                    "Donors 112/6533 failed registration; 6579 was unavailable. No failed section is silently omitted from the coverage tables.",
                    "Peters Arial quantified panel style/black spines; donor-point layout", ["protein_INS", "protein_GCG", "RNA lineage", "verified registration"])

    def reduced_panel(self) -> None:
        path = self.out / "typing_audit/panel_qualified_immune_counts.csv"
        d = pd.read_csv(path, dtype={"donor_id": str}, keep_default_na=False)
        donors = sorted(d.donor_id.unique(), key=self.donor_sort)
        den = d.drop_duplicates("donor_id").set_index("donor_id").denominator.reindex(donors)
        def field_matrix(field):
            part = d.loc[d.label_field.eq(field)]
            count = part.pivot(index="donor_id", columns="cell_identity", values="n_cells").fillna(0).reindex(donors)
            if not np.array_equal(count.sum(axis=1).to_numpy(), den.to_numpy()):
                raise ValueError("qualified tier must reconcile every primary Immune cell: " + field)
            count.assign(denominator_primary_Immune=den).to_csv(self.out / "source_tables" / (field + "_complete_counts.csv"))
            return count
        classes = field_matrix("qualified_immune_class")
        states = field_matrix("qualified_immune_subtype")
        support = field_matrix("qualified_class_confidence")
        if int(den.sum()) != 360035 or len(donors) != 4:
            raise ValueError("qualified tier population changed")
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Four 70-channel donors: completed panel-qualified immune typing",
                    "360,035 primary Immune cells evaluated | measured CD38 retained | unavailable Ker8/18 negative evidence explicitly omitted")
        ax = fig.add_axes([.095, .59, .34, .27])
        left = np.zeros(len(donors))
        for identity in sorted(classes.columns, key=lambda x: x in BLANK):
            value = classes[identity].to_numpy() / den.to_numpy() * 100
            label = "No measured rule matched" if identity in BLANK else identity.removeprefix("Imm-").replace("-", " ")
            ax.barh(np.arange(len(donors)), value, left=left, height=.65, linewidth=0,
                    color=self.identity_colors.get(identity, "#9e9e9e"), label=label)
            left += value
        ax.set_yticks(range(len(donors)), [self.unit_label(x) for x in donors])
        ax.set_ylim(len(donors) - .5, -.5)
        ax.set_xlim(0, 100)
        ax.set_xlabel("Primary Immune cells (%)")
        ax.set_title("Selected qualified class", loc="left", pad=14)
        clean_axes(ax)
        ax.legend(frameon=False, ncols=2, fontsize=7, bbox_to_anchor=(0, -.23), loc="upper left", handlelength=1)
        panel(ax, "A")
        bx = fig.add_axes([.63, .35, .28, .51])
        values = states.drop(columns=list(BLANK), errors="ignore").div(den, axis=0).T * 100
        image = bx.imshow(values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
        bx.set_xticks(range(len(donors)), donors)
        bx.set_yticks(range(len(values)), [x.replace("-", " ") for x in values.index], fontsize=7)
        bx.set_title("Selected qualified state", loc="left", pad=14)
        bx.tick_params(length=0)
        for y in range(len(values)):
            for x in range(len(donors)):
                v = values.iloc[y, x]
                number = "<0.01" if 0 < v < .01 else (f"{v:.2f}" if v < .1 and v > 0 else f"{v:.1f}")
                bx.text(x, y, number, ha="center", va="center", fontsize=7,
                        color="black" if v > 55 else "white")
        cbax = fig.add_axes([.65, .265, .24, .013])
        cb = fig.colorbar(image, cax=cbax, orientation="horizontal", ticks=[0, 50, 100])
        cb.outline.set_visible(False)
        cb.set_label("Primary Immune cells (%); shared scale", fontsize=7)
        panel(bx, "B")
        cx = fig.add_axes([.095, .19, .34, .18])
        left = np.zeros(len(donors))
        support_colors = {"single measured-rule support": "#279e68", "multiple compatible measured rules": "#ffbb78",
                          "no measured rule matched": "#9e9e9e"}
        if set(support.columns) != set(support_colors):
            raise ValueError("unexpected qualified support category")
        for identity in support.columns:
            value = support[identity].to_numpy() / den.to_numpy() * 100
            cx.barh(np.arange(len(donors)), value, left=left, height=.65, linewidth=0,
                    color=support_colors.get(identity, "#bdbdbd"), label=identity)
            left += value
        cx.set_yticks(range(len(donors)), donors)
        cx.set_ylim(len(donors) - .5, -.5)
        cx.set_xlim(0, 100)
        cx.set_xlabel("Primary Immune cells (%)")
        cx.set_title("Rule support and overlap", loc="left", pad=13)
        clean_axes(cx)
        panel(cx, "C")
        cx.legend(handles=[Patch(facecolor=support_colors[x], label=y) for x, y in [
            ("single measured-rule support", "One compatible class"),
            ("multiple compatible measured rules", "Multiple compatible classes"),
            ("no measured rule matched", "No class rule matched")]], frameon=False, ncols=1,
            fontsize=7, bbox_to_anchor=(1.03, .6), loc="upper left", handlelength=1)
        self.finish(fig, "figure15_qualified_reduced_panel_immune_evidence", "Qualified reduced-panel immune evidence",
                    "A partitions all primary Immune cells using the existing rule conflict resolver; 331,533 receive a qualified class. C retains the 69,727 cells with multiple class candidates. "
                    "All candidates remain in cell-level tables. B shows measured-rule state support within the selected qualified T/macrophage parents; no selected state is not evidence of biological absence. "
                    "These calls omit only the physically unacquired negative Ker8/18 term; measured positive CD38 and every other native marker requirement and threshold are retained. "
                    "Qualified labels do not satisfy the full-panel strict rule; state names do not independently establish function. Strict native fields, calibration and trajectories remain unchanged.",
                    "Atlas horizontal composition; Peters exact Arial and black-spine quantified panel style", ["qualified_immune_class", "qualified_immune_subtype", "qualified_class_confidence", "all rule candidates"])

    def vendor_umaps(self) -> None:
        source = self.out / "source_tables/rna_expression_umap_display_sample.parquet"
        if not source.exists():
            raise ValueError("all current RNA expression maps must be completed before the final atlas")
        d = pd.read_parquet(source)
        for lane, name, tissue_title in [("pancreas_rna", "figure16_pancreas_rna_vendor_umaps", "Pancreas RNA"),
                                        ("pln_rna", "figure17_pln_rna_vendor_umaps", "pLN RNA")]:
            part = d.loc[d.lane.eq(lane)].copy()
            populations = self.counts.loc[self.counts.lane.eq(lane), ["donor_id", "section_id"]].drop_duplicates()
            pairs = sorted(list(populations.itertuples(index=False, name=None)), key=lambda x: (self.donor_sort(x[0]), x[1]))
            fig = plt.figure(figsize=PAGE)
            self.header(fig, tissue_title + ": expression UMAPs with current native labels",
                        "All current sections | exact cell-ID joins | saved vendor maps and separately computed current-expression maps")
            nrows = int(np.ceil(len(pairs) / 5))
            gs = fig.add_gridspec(nrows, 5, left=.045, right=.97, bottom=.22, top=.88, wspace=.14, hspace=.28)
            axes = [fig.add_subplot(g) for g in gs]
            symbols = {}
            labels = sorted(part.native_lineage.unique())
            for ax, (donor, section) in zip(axes, pairs):
                cells = part.loc[part.donor_id.eq(donor) & part.section_id.eq(section)]
                ax.set_axis_off()
                section_num = symbols.get(donor, 0) + 1
                symbols[donor] = section_num
                repeated = int(populations.donor_id.eq(donor).sum()) > 1
                if cells.empty:
                    raise ValueError("missing measured RNA expression map: " + section)
                origin = cells.embedding_source.unique()
                if len(origin) != 1:
                    raise ValueError("mixed embedding source within section: " + section)
                map_label = "saved vendor map" if origin[0] == "saved_vendor" else "current-expression map"
                ax.set_title(self.unit_label(donor) + (" | section " + str(section_num) if repeated else "") + "\n" + map_label,
                             loc="left", fontsize=6.5, pad=2)
                for identity in sorted(labels, key=lambda x: x != "Exocrine" and x != "Non_immune"):
                    subset = cells.loc[cells.native_lineage.eq(identity)]
                    ax.scatter(subset.UMAP1, subset.UMAP2, c=self.identity_colors.get(identity, "#9e9e9e"),
                               s=.10, linewidths=0, rasterized=True, alpha=.10 if identity == "Exocrine" else 1.)
                ax.set_aspect("equal", adjustable="box")
            for ax in axes[len(pairs):]:
                ax.set_visible(False)
            fig.legend(handles=[Patch(facecolor=self.identity_colors.get(x, "#9e9e9e"), label=display_label(x)) for x in labels],
                       ncols=5, frameon=False, fontsize=6.6, loc="lower center", bbox_to_anchor=(.5, .095), handlelength=1)
            available = part[["donor_id", "section_id"]].drop_duplicates()
            saved = part.loc[part.embedding_source.eq("saved_vendor"), "section_id"].nunique()
            if len(available) != len(pairs):
                raise ValueError("incomplete RNA display coverage")
            self.finish(fig, name, tissue_title + " complete expression UMAPs",
                        f"All {len(pairs)} current sections are displayed: {saved} saved vendor maps and {len(pairs) - saved} newly calculated maps from exact current RNA counts. "
                        "Each panel uniformly displays 50,000 nucleated cells with fixed seeds. Saved coordinates are unchanged; new maps use normalized/log-transformed RNA, variable genes, PCA and UMAP. "
                        "Each map is independently fitted for display; positions/axes cannot be compared across panels and these maps are not the native classifier's graph. "
                        "Colors are current native calls, including unresolved labels; exocrine is faint context. They are not vendor cluster names or independent validation. Exact acquisitions, cell-ID joins, coverage and repeatability are documented. "
                        "No missing acquisition is replaced by an older section.", "Atlas/Peters compact cell-identity map style; exact installed Arial", ["saved/computed gene-expression UMAP1/2", "current native RNA lineage", "exact cell ID", "embedding_source"])

    def literature(self) -> None:
        root = Path("/data/io60-calibration/protein_literature_all_exclusions_20261005")
        donor = pd.read_csv(root / "protein_literature_donor_marker_endpoints.csv")
        inflammation = pd.read_csv(root / "protein_literature_donor_inflammation_endpoints.csv")
        summary = pd.read_csv(root / "protein_literature_donor_balanced_summary.csv")
        for path in root.glob("*.csv"):
            shutil.copy2(path, self.out / "source_tables" / path.name)
        fig, axes = plt.subplots(3, 1, figsize=PAGE)
        fig.subplots_adjust(left=.10, right=.96, top=.85, bottom=.22, hspace=.64)
        self.header(fig, "Peters and Medina-Serpas: measured protein endpoints", "Native labels and donor-balanced uncertainty; classification markers are descriptive evidence, not independent accuracy.")
        tissue_colors = {"pancreas": "#f8766d", "pLN": "#00c0b2"}
        settings = [(donor, "marker_auc_strict_cd8_vs_known_noncd8_t", "Peters_T_cell_marker_AUC", ["CD8", "TCF1", "TOX", "Granzyme_B", "LAG3", "PD_1", "CD39", "ICOS"], "CD8 vs other named T cells (AUC)", [0, 1], .5),
                    (donor.loc[donor.tissue.eq("pancreas")], "distance_spearman_marker_vs_signed_nearest_islet_um", "Peters_continuous_distance_association", ["CD8", "TCF1", "TOX", "Granzyme_B", "LAG3", "PD_1", "CD39", "ICOS"], "Association with continuous islet distance", [-1, 1], 0),
                    (inflammation.loc[inflammation.tissue.eq("pancreas")], "marker_auc_immune_vs_other_native_typed", "Medina_pancreatic_inflammatory_marker_distribution", ["IDO1", "HLA_DR", "HLA_A", "VISTA", "PD_L1"], "Immune vs other typed cells (AUC)", [0, 1], .5)]
        for ax, (frame, outcome, family, markers, ylabel, limits, null) in zip(axes, settings):
            for i, marker in enumerate(markers):
                for tissue in sorted(frame.tissue.unique()):
                    values = pd.to_numeric(frame.loc[frame.marker.eq(marker) & frame.tissue.eq(tissue), outcome], errors="coerce").dropna().to_numpy()
                    offset = (-.11 if tissue == "pancreas" else .11) if len(frame.tissue.unique()) > 1 else 0
                    x = i + offset
                    ax.scatter(x + np.linspace(-.04, .04, len(values)), values, s=9, c=tissue_colors[tissue], alpha=.6, linewidths=0)
                    s = summary.loc[summary.marker.eq(marker) & summary.tissue.eq(tissue) & summary.endpoint_family.eq(family)]
                    if len(s):
                        row = s.iloc[0]
                        ax.errorbar(x, row.donor_balanced_mean,
                                    yerr=[[max(0, row.donor_balanced_mean - row.bootstrap_95ci_low)], [max(0, row.bootstrap_95ci_high - row.donor_balanced_mean)]],
                                    fmt="_", ms=9, color="black", capsize=2, lw=.8)
            ax.axhline(null, color="#999999", lw=.5, ls="--")
            ax.set_xticks(range(len(markers)), [m.replace("_", " ") for m in markers], fontsize=7)
            ax.set_ylim(*limits)
            ax.set_ylabel(ylabel, fontsize=7.5)
            clean_axes(ax, grid=True)
            panel(ax, chr(65 + list(axes).index(ax)))
        axes[0].legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c, label=t, markersize=4) for t, c in tissue_colors.items()],
                       frameon=False, ncols=2, fontsize=7, loc="upper right", bbox_to_anchor=(1, 1.25))
        self.finish(fig, "figure18_literature_protein_endpoints", "Measured protein literature endpoints",
                    "Each point is a donor; black intervals are existing donor-bootstrap 95% CIs about the donor-balanced mean. AUC 0.5 indicates no separation. "
                    "Panel B is Spearman association with continuously measured signed micrometre distance, using no distance bins or rings. "
                    "Fine-T-cell endpoints exclude unavailable strict parents: the four reduced-panel donors remain in the typing census. IFNG protein is not measured and is not filled with zero. "
                    "Partial marker-supported states do not establish complete exhaustion, activation or TCR clonality.", "Peters quantified-panel Arial, donor points and thin gray guides", ["existing native protein endpoints", "continuous signed distance", "donor bootstrap"])
        source = REPO / "results/v3/release_20261003/literature_analysis/literature_rna_paired_tissue_effects.csv"
        effects = pd.read_csv(source)
        shutil.copy2(source, self.out / "source_tables" / source.name)
        fig, axes = plt.subplots(1, 2, figsize=PAGE, sharex=True)
        fig.subplots_adjust(left=.085, right=.96, top=.85, bottom=.23, wspace=.30)
        self.header(fig, "Peters and Medina-Serpas: donor-paired RNA evidence", "Actual measured features; tissue detection differences can reflect cell composition.")
        extent = np.nanmax(np.abs(effects[["mean_difference_bootstrap_95ci_low", "mean_difference_bootstrap_95ci_high"]].to_numpy())) * 100
        for i, paper in enumerate(["Peters", "Medina-Serpas"]):
            ax = axes[i]
            part = effects.loc[effects.paper.eq(paper)].dropna(subset=["mean_paired_difference"])
            for y, row in enumerate(part.itertuples()):
                ax.plot([row.mean_difference_bootstrap_95ci_low * 100, row.mean_difference_bootstrap_95ci_high * 100], [y, y], c="#555555", lw=.7)
                ax.scatter(row.mean_paired_difference * 100, y, s=15, c="#00bcd6", edgecolors="black", linewidths=.25)
            ax.set_yticks(range(len(part)), part.gene, fontsize=7, fontstyle="italic")
            ax.set_ylim(len(part) - .5, -.5)
            ax.axvline(0, color="#999999", lw=.6, ls="--")
            ax.set_xlim(-extent * 1.08, extent * 1.08)
            ax.set_xlabel("pLN minus pancreas detection (percentage points)", fontsize=7.5)
            ax.set_title(paper, loc="left", fontsize=10, weight="bold")
            clean_axes(ax)
            panel(ax, chr(65 + i))
        self.finish(fig, "figure19_literature_rna_paired_evidence", "Measured RNA literature endpoints",
                    "Points and intervals are the existing mean donor-paired differences and donor-bootstrap 95% CIs, restricted to exact current local acquisitions and panel coverage. "
                    "Pair counts and BH-adjusted q values for the predefined endpoint families are in literature_rna_paired_tissue_effects.csv. "
                    "LTB/CCL19/CCL21 and supported TCF7/TOX/activation features remain partial measured signatures; whole-tissue changes are not cell-intrinsic changes or complete phenotypes. "
                    "No TCR clonality is inferred. Medina-Serpas's raster figure interiors do not permit certified exact font recovery.", "Peters exact Arial quantified typography; Medina visual recreation", ["existing paired RNA endpoints", "measured panel coverage", "donor bootstrap"])

    def trajectories(self) -> None:
        d = pd.read_parquet(TRAJECTORY / "donor_endocrine_object_trajectories.parquet")
        m = json.loads((TRAJECTORY / "trajectory_manifest.json").read_text())
        d["donor_id"] = d.donor_id.map(norm_donor)
        good = d.trajectory_status.eq("ordered")
        nd = d.loc[good & d.donor_status.eq("ND")].groupby("donor_id")[["mfa_coordinate_1", "mfa_coordinate_2"]].median().median().to_numpy()
        coords = d.loc[good, ["mfa_coordinate_1", "mfa_coordinate_2"]].to_numpy()
        scale = np.maximum(np.quantile(np.abs(coords - nd), .975, axis=0), 1e-6)
        d["_x"] = (d.mfa_coordinate_1 - nd[0]) / scale[0]
        d["_y"] = (d.mfa_coordinate_2 - nd[1]) / scale[1]
        donors = sorted(d.donor_id.unique(), key=self.donor_sort)
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "Individual-donor endocrine-object ordering", "Frozen non-diabetic reference; shared display scales; the center marks typical ND composition.")
        gs = fig.add_gridspec(5, 5, left=.07, right=.92, bottom=.19, top=.87, wspace=.14, hspace=.30)
        markers = {"Islet": "o", "EndocrineCluster": "s", "EndocrineSingle": "^"}
        for ax, donor in zip([fig.add_subplot(g) for g in gs], donors):
            allpart = d.loc[d.donor_id.eq(donor)]
            part = allpart.loc[allpart.trajectory_status.eq("ordered")]
            for cls, cells in part.groupby("original_object_class"):
                ax.scatter(cells._x, cells._y, c=cells.within_donor_component_normalized_time, cmap="viridis", vmin=0, vmax=1,
                           s=3.0, marker=markers.get(cls, "o"), linewidths=0, alpha=.65, rasterized=True)
            for _, component in part.groupby("component_id"):
                roots = component.loc[component.endocrine_object_id.eq(component.component_root_object_id.iloc[0])]
                ax.scatter(roots._x, roots._y, marker="*", s=45, facecolor="white", edgecolor="black", linewidth=.65, zorder=5)
            ax.axvline(0, color="#dddddd", lw=.4, zorder=0)
            ax.axhline(0, color="#dddddd", lw=.4, zorder=0)
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
            clean_axes(ax)
            ax.set_title(self.unit_label(donor) + f" | {len(allpart):,} objects", fontsize=6.9, loc="left", pad=3)
        cax = fig.add_axes([.923, .29, .009, .45])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap="viridis"), cax=cax)
        cb.set_ticks([0, .5, 1], labels=["Start", "Middle", "Further"])
        cb.ax.tick_params(labelsize=7, length=2)
        cb.outline.set_visible(False)
        fig.text(.5, .15, composition_axis_groups(m, 0).capitalize() + " composition direction", ha="center", fontsize=8)
        fig.text(.025, .52, composition_axis_groups(m, 1).capitalize() + " composition direction", va="center", rotation=90, fontsize=8)
        fig.legend(handles=[plt.Line2D([], [], marker=marker, color="black", ls="", markersize=4, label=cls.replace("Endocrine", ""))
                            for cls, marker in markers.items()], frameon=False, ncols=3, fontsize=7,
                   loc="lower center", bbox_to_anchor=(.5, .103))
        clipped = int(((d.loc[good, ["_x", "_y"]].abs() > 1).any(axis=1)).sum())
        self.finish(fig, "figure13_all_donor_trajectory_views", "All donor endocrine-object orderings",
                    "Each point is a Single, Cluster or Islet; star = the component start closest to the frozen ND center. Color orders measured compositions within a donor and connected component, "
                    "not elapsed disease time. A robust ND-centered display and identical component scales are used for all donors; " + f"{clipped:,} extreme coordinates lie outside the display window and remain in data. "
                    f"{len(d):,} original objects retained; {int((~good).sum()):,} have no ordering, including 84 objects containing only quarantined cells (25 uncertain Singles). First two of {m['representation_rank']} reference components are displayed; ordering uses the complete frozen representation. Five measured groups enter the fit; stromal identity is unavailable. Numeric ticks and repetitive finite/stability scores are omitted.",
                    "Bodenmiller Arial/compact small-point panels; data-adapted donor paths, not the paper's developmental STREAM", ["frozen ND MFA coordinates", "within_donor_component_normalized_time"])
        contributions = pd.DataFrame({group: values[:2] for group, values in m["mfa_group_contributions_by_component"].items()}, index=["Composition direction 1", "Composition direction 2"])
        fig = plt.figure(figsize=PAGE)
        self.header(fig, "What the composition directions mean", "Measured feature contributions and cross-donor deviation are distinct from within-donor ordering.")
        ax = fig.add_axes([.10, .30, .34, .52])
        (contributions.T * 100).plot.barh(ax=ax, color=["#3182bd", "#fd8d3c"], width=.70, edgecolor="none")
        ax.set_xlabel("Share of the fitted direction (%)")
        ax.set_ylabel("")
        ax.invert_yaxis()
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.15), fontsize=7)
        clean_axes(ax)
        panel(ax, "A")
        ax.set_title("Biological groups behind the axes", loc="left", pad=12)
        bx = fig.add_axes([.62, .26, .30, .60])
        dist = d.groupby("donor_id").nd_reference_distance_mfa.agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)]).reindex(donors)
        dist.columns = ["median", "p25", "p75"]
        for y, (donor, row) in enumerate(dist.iterrows()):
            bx.plot([row.p25, row.p75], [y, y], color="#555555", lw=.7)
            bx.scatter(row["median"], y, s=18, c=GROUP_COLORS.get(self.disease.get(donor), "#999999"), edgecolors="black", linewidths=.25)
        bx.set_yticks(range(len(donors)), [self.unit_label(x) for x in donors], fontsize=7)
        bx.invert_yaxis()
        bx.set_xlabel("Distance from frozen ND composition center")
        bx.set_xticks([])
        bx.set_title("Cross-donor ND deviation", loc="left", pad=12)
        clean_axes(bx)
        panel(bx, "B")
        contributions.to_csv(self.out / "source_tables/trajectory_biological_group_contributions.csv")
        dist.to_csv(self.out / "source_tables/donor_nd_reference_distance_summary.csv")
        self.finish(fig, "figure14_trajectory_biological_interpretation", "Biological interpretation of the fitted directions",
                    "A reports contributions to the two displayed MFA directions, not manually assigned importance weights or signed biological effects. "
                    "B shows each donor's median and middle 50% of object distances in the same frozen ND representation. Larger distance means less ND-like measured composition, "
                    f"not longer disease duration or a severity score. After exclusion of all 84 objects containing only quarantined cells, ND-only validation selected rank {m['representation_rank']} and the frozen graph settings; no disease ordering was imposed.",
                    "Bodenmiller exact Arial; Peters compact donor-point/interval convention", ["mfa_group_contributions_by_component", "nd_reference_distance_mfa"])

    def run(self) -> None:
        (self.out / "source_tables").mkdir(exist_ok=True)
        self.composition("pancreas_protein", "cell_class", "typing_analysis_eligible_cells", "figure01_pancreas_protein_all_native_classes",
                         "Pancreas protein: endocrine, structural and immune typing", "25 donors | typing exclusions omitted | exocrine omitted from display; native counts retained", "Exocrine is broad: no acinar/ductal distinction is supported.")
        self.refinement()
        self.immune()
        self.composition("pancreas_rna", "lineage", "nucleated_cells", "figure04_pancreas_rna_native_lineages",
                         "Pancreas RNA: nonexocrine native lineages", "Development release | 25 donors | exocrine omitted from display; full counts retained",
                         "Xenium panels lack INS/GCG/SST/PPY and acinar enzyme targets; no RNA hormone subtype or acinar/ductal call is asserted. Donor 115 carries a B/Plasma coherence warning; zero assigned cells is not proof of absence.")
        self.composition("dedicated_pln_protein", "native_tnk_identity_display", "qc_analysis_eligible_cells", "figure05_dedicated_pln_protein_native_calls",
                         "Dedicated pLN protein: complete native immune typing", "19 slides | tissue-specific producer | native T/NK refinement assembled beside the broad calls",
                         "The display combines existing is_immune/immune_class/t_nk_class fields without new inference. Nonimmune structural subtypes remain unresolved.")
        self.composition("embedded_pln_protein", "native_tnk_identity_display", "qc_analysis_eligible_cells", "figure06_embedded_pln_protein_native_calls",
                         "Embedded pLN protein: every annotated node", "12 annotations across 10 donors | node membership and source-block relationships remain separate",
                         "All same-block/different-block assignments stay in the coverage table; no pancreas-to-node distance is inferred.", True)
        self.composition("pln_rna", "lineage", "nucleated_cells", "figure06b_pln_rna_complete_native_calls",
                         "pLN RNA: all native lineages and T/NK subsets", "Development release | 27 sections across 24 donors | nonimmune cells remain a broad compartment",
                         "Unresolved T subsets remain explicit. This producer does not name detailed nonimmune structural classes.", True)
        self.pln_subtypes()
        self.markers()
        self.spatial()
        self.objects()
        self.paired()
        self.shared()
        self.restain()
        self.trajectories()
        self.reduced_panel()
        self.vendor_umaps()
        self.literature()
        figure_manifest = {"font_files": [{"path": str(p), "sha256": sha(p)} for p in (ARIAL, ARIAL_BOLD, ARIAL_ITALIC)],
                           "identity_colors": self.identity_colors,
                           "exact_font_scope": "Arial portions of verified vector source panels; original raster-only figure fonts cannot be certified",
                           "figure_layout_scope": "data-adapted recreations of explicitly referenced panel styles, not numerical reproduction of published cohorts",
                           "stage_palette_mapping": "source Control/sAAb+/Onset colors reused for ND/AAb+/T1D visual group keys only; disease durations or single/multiple autoantibodies are not inferred",
                           "figures": self.figures}
        (self.out / "figure_manifest.json").write_text(json.dumps(figure_manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    configure()
    Builder(args.out.resolve()).run()


if __name__ == "__main__":
    main()
