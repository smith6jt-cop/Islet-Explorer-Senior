"""Shared data loading, genotype mapping, and analysis utilities.

Used by all H1–H10 hypothesis notebooks.
"""

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import kruskal, mannwhitneyu, spearmanr

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
CODE_ROOT = Path(__file__).resolve().parent.parent
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from identity import (  # noqa: E402
    SAMPLE_BY_ALIAS,
    SAMPLE_MANIFEST,
    SampleManifestError,
    load_sample_manifest as _load_sample_manifest,
    resolve_sample_id,
)
from paths import DATA_ROOT, OUTPUT_ROOT  # noqa: E402

# Backwards-compatible name used by the U-Net helpers. It means the external
# data root after absorption; tracked source lives at CODE_ROOT.
PROJECT = DATA_ROOT
DATA_CSV = PROJECT / "Measurements" / "AnnotationsFinal.csv"
GROUPS_XLSX = CODE_ROOT / "Groups.xlsx"
FIGURES_DIR = OUTPUT_ROOT / "analysis" / "figures"
TABLES_DIR = OUTPUT_ROOT / "analysis" / "tables"

GENO_ORDER = ["C/C", "C/T", "T/T"]
GENO_PALETTE = dict(zip(GENO_ORDER, sns.color_palette("Set2", 3)))
MAIN_REGIONS = ["Follicle", "PALS", "RedPulp", "Trabeculae"]
EXCLUDE_SAMPLES = {"HDL018", "HDL021", "HDL172"}

CODEX_SAMPLES = {"1901HBMP004", "HDL073", "HBMP006", "HBMP007", "HDL098"}
PHENOCYCLER_SAMPLES = {"HDL011", "HDL043", "HDL052", "HDL053", "HDL055",
                       "HDL063", "HDL070", "HDL079", "HDL086", "HDL094"}
CODEX_CSV = PROJECT / "Measurements" / "ForSH2B3.csv"

# Single-cell sources. Cells.csv covers the original 13 images; the second CSV
# adds the three annotated CODEX sections that give the single-cell cohort its
# first non-T/T CODEX donors.
CELLS_CSV = PROJECT / "Measurements" / "Cells.csv"
NEW_CODEX_CELLS_CSV = (PROJECT / "Measurements" /
                       "Spleen_Upload2_2024Sept17_v2_distances.csv")
# That CSV holds six sections; only these three carry annotations. The other
# three (20_06_CC2C, 20_07_CC2B, 20_08_CC2B) are excluded. They must be
# filtered by section name, not by sample: 20_06_CC2C is the same donor
# (HBMP006) as the annotated 20_06_CC2B.
NEW_CODEX_SECTIONS = ("20_06_CC2B.tif", "20_07_CC1A.tif", "20_08_CC3A.tif")

# CODEX region/class harmonization
#
# Two maps, deliberately different — see the Sinusoid note below.
#
# _CODEX_REGION_MAP rewrites the Parent string, i.e. it decides which region a
# child object (vessel) is counted under. _CODEX_CLASS_MAP rewrites the
# annotation's own Classification, i.e. it decides which region that
# annotation's *area* is summed into by compute_density().
#
# Sinusoid is in the Parent map but NOT the class map. In the CODEX export,
# Red_Pulp and Sinusoid are sibling annotations (both parented to the image
# root) that occupy the same territory: near-identical centroids, and Sinusoid
# carries a ~20x larger perimeter for a smaller area — the sinusoidal network
# threading through the red pulp, traced as its own overlapping annotation.
# Summing both areas double-counted the red pulp, putting total annotated
# tissue at 132-149% of the image field of view. Excluding Sinusoid's area
# brings the three CODEX images to 94-107% of FOV. Vessels parented to
# Annotation (Sinusoid) are still attributed to RedPulp, so the density
# numerator is unchanged; only the denominator is corrected.
_CODEX_REGION_MAP = {"Red_Pulp": "RedPulp", "Sinusoid": "RedPulp",
                     "Trabecula": "Trabeculae", "Peripheral_White_Pulp": "PALS"}
_CODEX_CLASS_MAP = {**{k: v for k, v in _CODEX_REGION_MAP.items()
                       if k != "Sinusoid"},
                    "Vessel": "SmallVessel",
                    "Vein": "LargeVessel", "Arteriole": "LargeVessel"}

# Clinical / SNP constants
CLINICAL_COLS = ["Age (yrs)", "Gender", "Ethnicity", "C-pep (ng/ml)", "HbA1c", "BMI"]
CONTINUOUS_CLINICAL = ["Age (yrs)", "C-pep (ng/ml)", "HbA1c", "BMI"]

def load_sample_manifest(path: Path = SAMPLE_MANIFEST) -> pd.DataFrame:
    """Load the exact alias table as a DataFrame for notebook consumers."""

    return pd.DataFrame(_load_sample_manifest(path))


_SAMPLE_MANIFEST = load_sample_manifest()
_SAMPLE_BY_ALIAS = SAMPLE_BY_ALIAS

# GeoJSON pixel → µm, per sample. Verified empirically; the OME-XML header is
# not authoritative here.
#
# - HBMP006/HBMP007/HDL098 are 0.37740007578193524. Confirmed by reproducing
#   every region area in ForSH2B3.csv to four decimal places.
# - 1901HBMP004 and HDL073 are 0.374 even though their OME-XML declares
#   0.5077663810243286 — a stale header written by convert_to_ome_tiff.py,
#   which hardcodes the Phenocycler pixel size. QuPath's own measurements and
#   signed distances for these two are at 0.374: recomputing Cells.csv signed
#   distances from GeoJSON gives r = 1.0000 (median |Δ| 0.08 µm) at 0.374,
#   versus r ≈ 0.13-0.21 at 0.5078.
_CODEX_PIXEL_SIZE = 0.37740007578193524
_PIXEL_SIZES = {
    "1901HBMP004": 0.374,
    "HDL073": 0.374,
    "HBMP006": _CODEX_PIXEL_SIZE,
    "HBMP007": _CODEX_PIXEL_SIZE,
    "HDL098": _CODEX_PIXEL_SIZE,
}
DEFAULT_PIXEL_SIZE = 0.5077663810243286


def get_pixel_size(sample: str) -> float:
    """µm per pixel for a sample's GeoJSON coordinates."""
    return _PIXEL_SIZES.get(sample, DEFAULT_PIXEL_SIZE)


# ---------------------------------------------------------------------------
# Sample ID extraction
# ---------------------------------------------------------------------------
def extract_sample_id(image_name: str) -> str:
    """Resolve an exact recorded alias to its canonical sample ID.

    Filename parsing is deliberately forbidden: a plausible prefix is not
    donor identity. Add a reviewed row to ``config/sample_manifest.csv`` when
    a new image, region stem, or metadata alias enters the workflow.
    """

    return resolve_sample_id(image_name)


# ---------------------------------------------------------------------------
# Genotype mapping
# ---------------------------------------------------------------------------
def _normalize_genotype(g):
    """Normalize a genotype string: T/C → C/T, G A → A/G."""
    if pd.isna(g):
        return np.nan
    s = str(g).strip()
    parts = re.split(r"[\s/]+", s)
    if len(parts) == 2:
        return "/".join(sorted(parts))
    return s


def _build_genotype_map() -> dict:
    """Build sample_id → genotype dictionary from Groups.xlsx."""
    g = pd.read_excel(GROUPS_XLSX)
    geno_col = "rs3184504" if "rs3184504" in g.columns else "rs3184504 (SH2B3)"
    geno_map = {}
    for _, row in g.iterrows():
        geno = row.get(geno_col)
        if pd.isna(geno):
            continue
        geno = _normalize_genotype(geno)
        # Map by HANDEL ID
        handel = row.get("HANDEL ID")
        if pd.notna(handel):
            geno_map[str(handel).strip()] = geno
        # ALT IDs are exact aliases in the manifest. No digit slicing.
        alt = row.get("ALT ID")
        if pd.notna(alt):
            alt_str = str(int(alt))
            geno_map[resolve_sample_id(alt_str)] = geno
    return geno_map


GENOTYPE_MAP = _build_genotype_map()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data() -> pd.DataFrame:
    """Load AnnotationsFinal.csv with Sample and Genotype columns.

    Drops excluded samples and rows without genotype.
    """
    df = pd.read_csv(DATA_CSV)
    df["Sample"] = df["Image"].apply(extract_sample_id)
    df["Genotype"] = df["Sample"].map(GENOTYPE_MAP)
    # Drop exclusions and unmapped
    df = df[~df["Sample"].isin(EXCLUDE_SAMPLES)]
    df = df.dropna(subset=["Genotype"])
    df["Genotype"] = pd.Categorical(df["Genotype"], categories=GENO_ORDER, ordered=True)
    return df


def load_all_data() -> pd.DataFrame:
    """Load AnnotationsFinal.csv + ForSH2B3.csv (CODEX), harmonized.

    Returns combined DataFrame with Sample, Genotype, Platform columns.
    """
    df = load_data()
    # Platform is a property of the sample, not of the source CSV. Two CODEX
    # samples (1901HBMP004, HDL073) live in AnnotationsFinal.csv, so a blanket
    # "Phenocycler" here mislabelled them. Derive it from CODEX_SAMPLES, the
    # same rule load_clinical() and build_feature_matrix() use.
    df["Platform"] = df["Sample"].apply(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")

    if CODEX_CSV.exists():
        codex = pd.read_csv(CODEX_CSV)
        codex["Classification"] = codex["Classification"].map(
            lambda c: _CODEX_CLASS_MAP.get(c, c))
        codex["Parent"] = codex["Parent"].str.replace(
            r"Annotation \((\w+)\)",
            lambda m: f"Annotation ({_CODEX_REGION_MAP.get(m.group(1), m.group(1))})",
            regex=True,
        )
        codex["Sample"] = codex["Image"].apply(extract_sample_id)
        codex["Genotype"] = codex["Sample"].map(GENOTYPE_MAP)
        codex["Platform"] = "CODEX"
        codex = codex[~codex["Sample"].isin(EXCLUDE_SAMPLES)]
        codex = codex.dropna(subset=["Genotype"])
        codex["Genotype"] = pd.Categorical(
            codex["Genotype"], categories=GENO_ORDER, ordered=True)
        df = pd.concat([df, codex], ignore_index=True)
        df["Genotype"] = pd.Categorical(
            df["Genotype"], categories=GENO_ORDER, ordered=True)

    return df


# ---------------------------------------------------------------------------
# Single-cell loading (both cohorts, harmonized)
# ---------------------------------------------------------------------------
#
# The two CSVs spell marker columns differently — "Cell: {m}: Mean" in
# Cells.csv, "{m}: Cell: Mean" in the new CODEX one — and five markers differ
# in name between the platforms. Cells.csv already merges both cohorts, so it
# can carry *both* spellings of the same marker; those are coalesced.

X_COL, Y_COL = "Centroid X µm", "Centroid Y µm"

_MARKER_ALIASES = {
    "HLA-DR": ("HLA-DR", "HLADR"),
    "FOXP3": ("FOXP3", "FoxP3"),
    "SMA": ("SMA", "SMActin"),
    "CollagenIV": ("Collagen IV", "CollagenIV"),
    "Ecadherin": ("E-cadherin", "ECAD"),
}

# Markers all-NaN in at least one of the three new CODEX sections. Verified by
# direct non-null counts:
#   20_06_CC2B  CD68, Podoplanin, CD44
#   20_07_CC1A  CD11c, PanCK, CD44
#   20_08_CC3A  CD1c, CollagenIV, CD44
MARKERS_LOST_IN_NEW_CODEX = ("CD11c", "CD1c", "CD44", "CD68",
                             "CollagenIV", "PanCK", "Podoplanin")

# Markers with data in *every* one of the 14 donors — the only set an analysis
# spanning both platforms may use. Measured per donor from both CSVs, not
# inferred from column presence: a column can exist and be entirely NaN, which
# is how CD35 and Lyve1 look Phenocycler-available while carrying no data
# there. Of the 39 canonical markers across the two cohorts, 21 fail this test.
#
# H15/H16 cluster on 22 markers, of which CD44, CD11c, CD68, Podoplanin and
# CollagenIV are unusable — leaving the 17-marker panel for the 14-donor run.
MARKERS_USABLE_ALL_DONORS = (
    "CD107a", "CD163", "CD20", "CD21", "CD31", "CD34", "CD3e", "CD4", "CD45",
    "CD45RO", "CD8", "DAPI", "Ecadherin", "FOXP3", "HLA-DR", "Ki67", "SMA",
    "Vimentin",
)


def _marker_columns(available, marker):
    """Every column spelling in `available` that carries `marker`."""
    return [pat for name in _MARKER_ALIASES.get(marker, (marker,))
            for pat in (f"Cell: {name}: Mean", f"{name}: Cell: Mean")
            if pat in available]


def load_single_cells(markers=(), images=None, extra_columns=(),
                      include_excluded=False):
    """Load both single-cell cohorts with canonical marker names.

    Parameters
    ----------
    markers : iterable of canonical marker names (e.g. "CD20", "HLA-DR")
    images : optional iterable of Image values to keep
    extra_columns : additional Cells.csv columns to carry through (e.g. the
        QuPath signed-distance columns, which exist only for that cohort)
    include_excluded : keep the quality-excluded donors (HDL018, HDL021) that
        still carry a genotype. H12 and H13 deliberately re-include them to
        enlarge the C/T group, and they need region polygons to do that on the
        region layer rather than on the broken `Parent` hierarchy. HDL172 has
        no genotype and is dropped either way.

    Returns
    -------
    DataFrame with Image, Sample, Genotype, Platform, x, y and one float
    column per requested marker. Markers absent from a cohort come back NaN.
    """
    markers = list(markers)
    frames = []

    for path, sections in ((CELLS_CSV, None),
                           (NEW_CODEX_CELLS_CSV, NEW_CODEX_SECTIONS)):
        if not path.exists():
            continue
        available = pd.read_csv(path, nrows=1).columns.tolist()
        wanted = {m: _marker_columns(available, m) for m in markers}
        extras = [c for c in extra_columns if c in available]
        use = list(dict.fromkeys(          # dedupe, preserve order
            ["Image", X_COL, Y_COL]
            + sorted({c for cols in wanted.values() for c in cols}) + extras))
        raw = pd.read_csv(path, usecols=use, engine="pyarrow")
        if sections is not None:
            raw = raw[raw["Image"].isin(sections)]
        if images is not None:
            raw = raw[raw["Image"].isin(set(images))]
        if raw.empty:
            continue

        out = pd.DataFrame({"Image": raw["Image"].astype(str),
                            "x": raw[X_COL].astype("float32"),
                            "y": raw[Y_COL].astype("float32")})
        for m, cols in wanted.items():
            if not cols:
                out[m] = np.float32("nan")
                continue
            v = raw[cols[0]].astype("float32")
            for c in cols[1:]:          # coalesce alternate spellings
                v = v.fillna(raw[c].astype("float32"))
            out[m] = v
        for c in extras:
            out[c] = raw[c].values
        frames.append(out)

    if not frames:
        return pd.DataFrame(columns=["Image", "Sample", "Genotype", "Platform",
                                     "x", "y"] + markers)

    cells = pd.concat(frames, ignore_index=True)
    cells["Sample"] = cells["Image"].apply(extract_sample_id)
    cells["Genotype"] = cells["Sample"].map(GENOTYPE_MAP)
    cells["Platform"] = cells["Sample"].apply(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")
    if not include_excluded:
        cells = cells[~cells["Sample"].isin(EXCLUDE_SAMPLES)]
    cells = cells.dropna(subset=["Genotype"])
    cells["Genotype"] = pd.Categorical(cells["Genotype"],
                                       categories=GENO_ORDER, ordered=True)
    return cells.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------
def get_regions(df: pd.DataFrame) -> pd.DataFrame:
    """Return region-level annotations (Follicle, PALS, RedPulp, Trabeculae, LargeVessel)."""
    region_classes = {"Follicle", "PALS", "RedPulp", "Trabeculae", "LargeVessel"}
    return df[df["Classification"].isin(region_classes)].copy()


def assign_regions_by_image(df: pd.DataFrame, region_layer=None,
                            x_col: str = "Centroid X µm",
                            y_col: str = "Centroid Y µm",
                            image_col: str = "Image") -> np.ndarray:
    """Region label for every row, by which polygon contains its centroid.

    The one region-assignment path. Objects outside every polygon come back
    "Unassigned" — cohort-wide that is 4-9% of SmallVessels, sitting just
    beyond the KDE tissue contour at section edges, so it is an expected
    category rather than an error. An image with no polygons at all is left
    entirely "Unassigned".

    Works for annotations and for single cells; pass x_col/y_col when the frame
    spells the centroid columns differently (load_single_cells uses x/y).
    """
    if region_layer is None:
        region_layer = load_region_layer()
    by_sample = {extract_sample_id(stem): polys
                 for stem, polys in region_layer.items()}

    out = np.full(len(df), "Unassigned", dtype=object)
    positions = pd.Series(np.arange(len(df)), index=df.index)
    for image, g in df.groupby(image_col, observed=True):
        polys = region_layer.get(_image_stem(str(image)),
                                 by_sample.get(extract_sample_id(str(image))))
        if not polys:
            continue
        coords = g[[x_col, y_col]].values.astype(float)
        out[positions.loc[g.index].values] = assign_objects_to_polygons(coords, polys)
    return out


def signed_distance_by_image(df: pd.DataFrame, region: str, region_layer=None,
                             x_col: str = "Centroid X µm",
                             y_col: str = "Centroid Y µm",
                             image_col: str = "Image",
                             boundary_step: float = 8.0) -> np.ndarray:
    """Signed distance to one region's polygons, per image, negative inside.

    The region-layer counterpart to QuPath's "Signed distance to annotation X
    µm" columns. Prefer this over those columns for Follicle and PALS: QuPath
    measured to the *unfiltered* annotations, so a cell sitting inside a
    fragment the area floor removed reads as inside a follicle while being
    labelled RedPulp. Measuring against the same polygons the labels came from
    makes the two agree by construction.

    Images with no polygons for the region come back NaN.
    """
    if region_layer is None:
        region_layer = load_region_layer()
    by_sample = {extract_sample_id(stem): polys
                 for stem, polys in region_layer.items()}

    out = np.full(len(df), np.nan, dtype="float64")
    positions = pd.Series(np.arange(len(df)), index=df.index)
    for image, g in df.groupby(image_col, observed=True):
        polys = (region_layer.get(_image_stem(str(image)),
                                  by_sample.get(extract_sample_id(str(image))))
                 or {}).get(region, [])
        coords = g[[x_col, y_col]].values.astype(float)
        out[positions.loc[g.index].values] = signed_distance_to_regions(
            coords, polys, boundary_step=boundary_step)
    return out


def get_vessels(df: pd.DataFrame, region_layer=None,
                from_parent: bool = False) -> pd.DataFrame:
    """Return SmallVessel annotations with a Region column.

    Region comes from which region-layer polygon contains the vessel centroid.
    The old route — parsing QuPath's `Parent` string — is kept behind
    ``from_parent=True`` only so the two can be compared; it must not be used
    for analysis. The hierarchy it reads is known broken: 1,014 vessels are
    orphaned at the image root, and 1901HBMP004 has 35 Follicle annotations
    with zero Follicle-parented vessels, so `Parent` silently reports no
    follicular vasculature for a section that plainly has some.

    Vessels landing outside every polygon get Region "Unassigned" rather than
    NaN, so they can be counted rather than vanishing in a groupby.
    """
    v = df[df["Classification"] == "SmallVessel"].copy()
    if from_parent:
        v["Region"] = v["Parent"].str.extract(r"Annotation \((\w+)\)")
        return v

    v["Region"] = assign_regions_by_image(v, region_layer=region_layer)
    return v


# ---------------------------------------------------------------------------
# Density computation
# ---------------------------------------------------------------------------
def compute_density(df: pd.DataFrame, region_layer=None,
                    regions=None) -> pd.DataFrame:
    """Vessel density per image per region.

    Delegates to :func:`compute_density_from_regions`. This used to sum
    annotation areas and read a vessel's region off the QuPath `Parent` string;
    both were wrong in ways that could not be patched. `Parent` loses the 1,014
    vessels orphaned at the image root, and summing annotation areas
    double-counts the CODEX Red_Pulp/Sinusoid pair. The region layer fixes both,
    so the two paths are no longer kept side by side.

    Returns DataFrame with columns:
      Image, Sample, Genotype, Region, Vessel_Count, Region_Area_um2,
      Region_Area_mm2, Density_per_mm2, RedPulp_Density, Density_Normalized
    """
    if regions is None:
        regions = MAIN_REGIONS
    return compute_density_from_regions(df, region_layer=region_layer,
                                        regions=regions)


REGIONS_DIR = OUTPUT_ROOT / "analysis" / "regions"


def load_region_layer(regions_dir=None):
    """Load the persisted region polygons.

    Built by scripts/build_region_layer_qupath.py from the QuPath pixel-
    classifier exports. Coordinates are µm, matching the annotation and cell
    centroids, so no pixel-size conversion is needed.

    Returns
    -------
    {image_stem: {region_name: [Polygon]}}
    """
    from shapely.geometry import shape

    d = Path(regions_dir) if regions_dir else REGIONS_DIR
    if not d.exists():
        raise FileNotFoundError(
            f"{d} not found — run scripts/build_region_layer_qupath.py first")

    layer = {}
    for path in sorted(d.glob("*.geojson")):
        by_region = {}
        for feat in json.load(open(path)).get("features", []):
            name = (feat.get("properties", {}) or {}).get("region")
            by_region.setdefault(name, []).append(shape(feat["geometry"]))
        layer[path.stem] = by_region
    return layer


def _image_stem(image: str) -> str:
    """Match an Image value to a region-layer key."""
    return (image.replace(".ome.tiff", "").replace(".ome.tif", "")
            .replace(".tiff", "").replace(".tif", ""))


def compute_density_from_regions(df: pd.DataFrame, region_layer=None,
                                 regions=None) -> pd.DataFrame:
    """Vessel density per image per region, using cell-derived regions.

    The polygon-based sibling of :func:`compute_density`. Same output schema,
    so downstream code is unchanged. Two differences in what it measures:

    - A vessel's region comes from which polygon contains its centroid, not
      from the QuPath `Parent` string. That is what makes the result immune to
      the broken hierarchy (1,014 vessels orphaned at the image root) and to
      the CODEX Red_Pulp/Sinusoid area double-count.
    - Region area is the cell-derived polygon area, so RedPulp is the tissue
      remainder rather than a sum of annotation areas.

    Parameters
    ----------
    df : annotation DataFrame (needs Image, Sample, Genotype, Classification,
        the centroid columns and Area µm^2) — SmallVessel rows supply the
        objects being counted
    region_layer : output of :func:`load_region_layer`; loaded if omitted
    regions : region names to report (default MAIN_REGIONS + Perifollicular)

    Returns
    -------
    DataFrame with Image, Sample, Genotype, Region, Vessel_Count,
    Region_Area_mm2, Density_per_mm2, RedPulp_Density, Density_Normalized
    """
    from shapely.ops import unary_union

    if region_layer is None:
        region_layer = load_region_layer()
    if regions is None:
        regions = MAIN_REGIONS + ["Perifollicular"]

    # Fall back to the sample when the stem does not match: ForSH2B3.csv spells
    # the CODEX sections "2006_CC2B" while the single-cell CSV (and so the
    # region layer) spells them "20_06_CC2B". extract_sample_id resolves both.
    by_sample = {extract_sample_id(stem): polys
                 for stem, polys in region_layer.items()}

    vessels = df[df["Classification"] == "SmallVessel"]
    rows = []
    for image, g in df.groupby("Image", observed=True):
        sample = g["Sample"].iloc[0]
        polys_by_region = region_layer.get(_image_stem(image),
                                           by_sample.get(sample))
        if polys_by_region is None:
            continue
        genotype = g["Genotype"].iloc[0]

        v = vessels[vessels["Image"] == image]
        coords = v[["Centroid X µm", "Centroid Y µm"]].values
        labels = (assign_objects_to_polygons(coords, polys_by_region)
                  if len(coords) else np.array([], dtype=object))
        counts = pd.Series(labels).value_counts()

        for region in regions:
            polys = polys_by_region.get(region, [])
            area_mm2 = (unary_union(polys).area / 1e6) if polys else 0.0
            n = int(counts.get(region, 0))
            rows.append({"Image": image, "Sample": sample, "Genotype": genotype,
                         "Region": region, "Vessel_Count": n,
                         "Region_Area_um2": area_mm2 * 1e6,
                         "Region_Area_mm2": area_mm2,
                         "Density_per_mm2": n / area_mm2 if area_mm2 > 0 else np.nan})

    density = pd.DataFrame(rows)
    if density.empty:
        return density

    rp = (density[density["Region"] == "RedPulp"][["Image", "Density_per_mm2"]]
          .rename(columns={"Density_per_mm2": "RedPulp_Density"}))
    density = density.merge(rp, on="Image", how="left")
    density["Density_Normalized"] = (density["Density_per_mm2"]
                                     / density["RedPulp_Density"])
    return density


# ---------------------------------------------------------------------------
# Spatial assignment (for H5)
# ---------------------------------------------------------------------------
def assign_vessels_to_follicles(df: pd.DataFrame, image: str,
                                region_layer=None) -> pd.DataFrame:
    """Assign intra-follicular SmallVessels to the nearest follicle centroid.

    Which vessels count as intra-follicular now comes from the region layer,
    not from the QuPath `Parent` string — see :func:`get_vessels`. Nearest
    follicle is still a cKDTree lookup on centroids.

    Returns vessel DataFrame with Follicle_ID and Follicle_Area columns.
    """
    img_df = df[df["Image"] == image]
    follicles = img_df[img_df["Classification"] == "Follicle"].copy()
    vessels = get_vessels(img_df, region_layer=region_layer)
    vessels = vessels[vessels["Region"] == "Follicle"].copy()

    if follicles.empty or vessels.empty:
        return pd.DataFrame()

    # Build tree from follicle centroids
    fol_coords = follicles[["Centroid X µm", "Centroid Y µm"]].values
    tree = cKDTree(fol_coords)

    # Query vessel centroids
    ves_coords = vessels[["Centroid X µm", "Centroid Y µm"]].values
    _, indices = tree.query(ves_coords)

    vessels["Follicle_ID"] = follicles["Object ID"].values[indices]
    vessels["Follicle_Area"] = follicles["Area µm^2"].values[indices]

    return vessels


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------
def rank_biserial(x, y):
    """Rank-biserial effect size: r = 1 - 2U/(n1*n2)."""
    u_stat, _ = mannwhitneyu(x, y, alternative="two-sided")
    n1, n2 = len(x), len(y)
    return 1 - 2 * u_stat / (n1 * n2)


def run_kruskal(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Run Kruskal-Wallis test across genotype groups.

    Returns (H_statistic, p_value) or (NaN, NaN) if <2 groups have data.
    """
    groups = [g[value_col].dropna().values for _, g in data.groupby(group_col, observed=True)]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return np.nan, np.nan
    h, p = kruskal(*groups)
    return h, p


def run_pairwise(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Run pairwise Mann-Whitney U tests with rank-biserial effect size.

    Returns list of dicts with keys: Comparison, U, p, r, n1, n2.
    """
    results = []
    groups = data.groupby(group_col, observed=True)[value_col]
    group_dict = {name: g.dropna().values for name, g in groups}
    pairs = [("C/C", "C/T"), ("C/C", "T/T"), ("C/T", "T/T")]
    for g1, g2 in pairs:
        x = group_dict.get(g1, np.array([]))
        y = group_dict.get(g2, np.array([]))
        if len(x) < 1 or len(y) < 1:
            results.append({"Comparison": f"{g1} vs {g2}", "U": np.nan, "p": np.nan,
                            "r": np.nan, "n1": len(x), "n2": len(y)})
            continue
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        r = 1 - 2 * u / (len(x) * len(y))
        results.append({"Comparison": f"{g1} vs {g2}", "U": u, "p": p,
                        "r": r, "n1": len(x), "n2": len(y)})
    return results


def run_dosage_trend(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Spearman correlation with ordinal genotype (C/C=0, C/T=1, T/T=2).

    Returns (rho, p_value).
    """
    ordinal = data[group_col].map({"C/C": 0, "C/T": 1, "T/T": 2})
    valid = ordinal.notna() & data[value_col].notna()
    if valid.sum() < 3:
        return np.nan, np.nan
    return spearmanr(ordinal[valid], data.loc[valid, value_col])


def full_stats_table(data: pd.DataFrame, value_col: str, label: str = ""):
    """Run all three statistical tests and return a summary DataFrame."""
    h, kw_p = run_kruskal(data, value_col)
    pw = run_pairwise(data, value_col)
    rho, sp_p = run_dosage_trend(data, value_col)

    rows = [{"Test": "Kruskal-Wallis", "Metric": label, "Statistic": h, "p": kw_p, "Effect_Size": ""}]
    for r in pw:
        rows.append({"Test": f"Mann-Whitney ({r['Comparison']})", "Metric": label,
                      "Statistic": r["U"], "p": r["p"], "Effect_Size": f"r={r['r']:.3f}"})
    rows.append({"Test": "Spearman dosage", "Metric": label,
                  "Statistic": rho, "p": sp_p, "Effect_Size": f"rho={rho:.3f}" if not np.isnan(rho) else ""})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def setup_style():
    """Configure seaborn/matplotlib defaults."""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams["figure.facecolor"] = "white"


def save_figure(fig, name: str, tight=True):
    """Save figure to analysis/figures/ as PNG."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / f"{name}.png"
    if tight:
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    else:
        fig.savefig(path, dpi=150, facecolor="white")
    print(f"Saved: {path}")


def save_table(df: pd.DataFrame, name: str):
    """Save DataFrame to analysis/tables/ as CSV."""
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLES_DIR / f"{name}.csv"
    df.to_csv(path, index=False)
    print(f"Saved: {path}")


# ---------------------------------------------------------------------------
# Clinical data loading (H8, H10)
# ---------------------------------------------------------------------------
def _resolve_sample_from_row(row) -> str | None:
    """Resolve a Groups.xlsx row to a canonical sample ID."""
    handel = row.get("HANDEL ID")
    if pd.notna(handel):
        return str(handel).strip()
    alt = row.get("ALT ID")
    if pd.notna(alt):
        return resolve_sample_id(str(int(alt)))
    return None


def load_clinical() -> pd.DataFrame:
    """Load clinical metadata from Groups.xlsx.

    Returns DataFrame with columns: Sample, Genotype, Platform,
    Age (yrs), Gender, Ethnicity, C-pep (ng/ml), HbA1c, BMI.
    """
    g = pd.read_excel(GROUPS_XLSX)
    g["Sample"] = g.apply(_resolve_sample_from_row, axis=1)
    g = g.dropna(subset=["Sample"])
    g = g[~g["Sample"].isin(EXCLUDE_SAMPLES)]

    # Genotype
    geno_col = "rs3184504" if "rs3184504" in g.columns else "rs3184504 (SH2B3)"
    g["Genotype"] = g[geno_col].apply(_normalize_genotype)
    g = g.dropna(subset=["Genotype"])
    g["Genotype"] = pd.Categorical(g["Genotype"], categories=GENO_ORDER, ordered=True)

    # Platform
    g["Platform"] = g["Sample"].apply(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")

    # Clean clinical columns
    out = g[["Sample", "Genotype", "Platform"] + CLINICAL_COLS].copy()
    for col in CONTINUOUS_CLINICAL:
        out[col] = pd.to_numeric(
            out[col].where(out[col] != "No Data"), errors="coerce")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# SNP panel loading (H9)
# ---------------------------------------------------------------------------
def load_snp_panel(project_samples=None):
    """Load SNP genotype and risk-dosage matrices from Groups.xlsx.

    Returns (geno_df, dosage_df):
      - geno_df: Sample × SNP matrix of normalized genotype strings
      - dosage_df: Sample × risk-annotated-SNP matrix of 0/1/2 dosage
    """
    g = pd.read_excel(GROUPS_XLSX)
    g["Sample"] = g.apply(_resolve_sample_from_row, axis=1)
    g = g.dropna(subset=["Sample"])
    g = g[~g["Sample"].isin(EXCLUDE_SAMPLES)]

    snp_cols = [c for c in g.columns
                if c.startswith("rs") and "_risk_het_prot" not in c]
    risk_cols = {c.replace("_risk_het_prot", ""): c
                 for c in g.columns if c.endswith("_risk_het_prot")}

    # Build genotype matrix
    geno_records = {}
    for _, row in g.iterrows():
        sid = row["Sample"]
        if project_samples and sid not in project_samples:
            continue
        rec = {}
        for snp in snp_cols:
            val = row[snp]
            rec[snp] = _normalize_genotype(val) if pd.notna(val) else np.nan
        geno_records[sid] = rec
    geno_df = pd.DataFrame.from_dict(geno_records, orient="index")
    geno_df.index.name = "Sample"

    # Build dosage matrix from risk annotations
    dosage_map = {"Protective": 0, "Het": 1, "Risk": 2}
    dosage_records = {}
    for _, row in g.iterrows():
        sid = row["Sample"]
        if project_samples and sid not in project_samples:
            continue
        rec = {}
        for snp, risk_col in risk_cols.items():
            val = row.get(risk_col)
            if pd.notna(val) and val in dosage_map:
                rec[snp] = dosage_map[val]
            else:
                rec[snp] = np.nan
        dosage_records[sid] = rec
    dosage_df = pd.DataFrame.from_dict(dosage_records, orient="index")
    dosage_df.index.name = "Sample"

    # Drop monomorphic SNPs from dosage
    polymorphic = dosage_df.columns[dosage_df.nunique() > 1]
    dosage_df = dosage_df[polymorphic]

    return geno_df, dosage_df


# ---------------------------------------------------------------------------
# Feature matrix construction (H8, H9, H10)
# ---------------------------------------------------------------------------
def build_feature_matrix(df: pd.DataFrame, region_layer=None) -> pd.DataFrame:
    """Build per-sample morphological feature matrix from annotation data.

    Returns ~23 numeric columns indexed by Sample, plus Genotype and Platform.

    The region layer is loaded once and shared: compute_density and get_vessels
    would otherwise parse ~130 MB of GeoJSON apiece.
    """
    if region_layer is None:
        region_layer = load_region_layer()
    density = compute_density(df, region_layer=region_layer)
    regions = get_regions(df)
    vessels = get_vessels(df, region_layer=region_layer)

    features = {}

    # --- Vessel density per region (raw + normalized) ---
    for region in MAIN_REGIONS:
        rd = density[density["Region"] == region]
        features[f"{region}_density"] = rd.set_index("Sample")["Density_per_mm2"]
        if region != "RedPulp":
            features[f"{region}_norm_density"] = rd.set_index("Sample")["Density_Normalized"]

    # --- Follicle metrics ---
    follicles = regions[regions["Classification"] == "Follicle"]
    fol_per_img = follicles.groupby(["Sample"], observed=True).agg(
        Follicle_count=("Object ID", "count"),
        Follicle_total_area=("Area µm^2", "sum"),
        Follicle_mean_area=("Area µm^2", "mean"),
    )
    for col in fol_per_img.columns:
        features[col] = fol_per_img[col]

    # Follicle fraction of total tissue
    total_area = regions[regions["Classification"].isin(MAIN_REGIONS)].groupby(
        "Sample", observed=True)["Area µm^2"].sum()
    fol_area = follicles.groupby("Sample", observed=True)["Area µm^2"].sum()
    features["Follicle_fraction"] = fol_area / total_area

    # --- Tissue proportions ---
    region_areas = (
        regions[regions["Classification"].isin(MAIN_REGIONS)]
        .groupby(["Sample", "Classification"], observed=True)["Area µm^2"]
        .sum()
        .unstack(fill_value=0)
    )
    region_total = region_areas.sum(axis=1)
    for r in MAIN_REGIONS:
        if r in region_areas.columns:
            features[f"{r}_fraction"] = region_areas[r] / region_total

    # --- White pulp metrics ---
    wp_area = pd.Series(0.0, index=region_total.index)
    for r in ["Follicle", "PALS"]:
        if r in region_areas.columns:
            wp_area = wp_area + region_areas[r]
    features["WP_fraction"] = wp_area / region_total
    if "Follicle" in region_areas.columns and "PALS" in region_areas.columns:
        features["Follicle_PALS_ratio"] = region_areas["Follicle"] / region_areas["PALS"].replace(0, np.nan)

    # --- Vessel morphology (median per sample) ---
    sv = df[df["Classification"] == "SmallVessel"].copy()
    if not sv.empty:
        sv["Elongation"] = sv["Max diameter µm"] / sv["Min diameter µm"].replace(0, np.nan)
        morph_agg = sv.groupby("Sample", observed=True).agg(
            Vessel_median_area=("Area µm^2", "median"),
            Vessel_median_circularity=("Circularity", "median"),
            Vessel_median_solidity=("Solidity", "median"),
            Vessel_median_elongation=("Elongation", "median"),
        )
        for col in morph_agg.columns:
            features[col] = morph_agg[col]

    # Combine all features
    feat_df = pd.DataFrame(features)
    feat_df.index.name = "Sample"

    # Add Genotype and Platform
    sample_meta = df.drop_duplicates("Sample").set_index("Sample")[["Genotype"]]
    feat_df = feat_df.join(sample_meta)
    feat_df["Platform"] = feat_df.index.map(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")

    return feat_df


# ---------------------------------------------------------------------------
# Polygenic risk score (H9)
# ---------------------------------------------------------------------------
def compute_prs(dosage_df: pd.DataFrame, min_snps: int = 10) -> pd.Series:
    """Compute unweighted polygenic risk score from risk dosage matrix.

    PRS = sum(dosage) / (2 × N_genotyped), normalized to [0,1].
    Returns NaN for samples with fewer than min_snps genotyped.
    """
    n_genotyped = dosage_df.notna().sum(axis=1)
    total_dosage = dosage_df.sum(axis=1)
    prs = total_dosage / (2 * n_genotyped)
    prs[n_genotyped < min_snps] = np.nan
    prs.name = "PRS"
    return prs


# ---------------------------------------------------------------------------
# Partial correlation (H8)
# ---------------------------------------------------------------------------
def partial_spearman(data: pd.DataFrame, x_col: str, y_col: str,
                     covariate_cols: list) -> tuple:
    """Partial Spearman correlation controlling for covariates.

    Rank-transforms all variables, regresses out covariates via OLS residuals,
    then computes Pearson correlation on residuals.
    Returns (rho, p_value).
    """
    from scipy.stats import pearsonr

    cols = [x_col, y_col] + list(covariate_cols)
    sub = data[cols].dropna()
    if len(sub) < 4:
        return np.nan, np.nan

    # Rank-transform
    ranked = sub.rank()

    # Regress out covariates
    if covariate_cols:
        C = ranked[covariate_cols].values
        C = np.column_stack([C, np.ones(len(C))])
        for col in [x_col, y_col]:
            y_vec = ranked[col].values
            beta, _, _, _ = np.linalg.lstsq(C, y_vec, rcond=None)
            ranked[col] = y_vec - C @ beta

    rho, p = pearsonr(ranked[x_col], ranked[y_col])
    return rho, p


# ---------------------------------------------------------------------------
# Platform diagnostic (H8, H10)
# ---------------------------------------------------------------------------
def platform_diagnostic(feature_df: pd.DataFrame,
                        feature_cols: list) -> pd.DataFrame:
    """Mann-Whitney test for CODEX vs Phenocycler on each feature.

    Returns DataFrame with Feature, U, p, rank_biserial, n_CODEX, n_PC.
    """
    results = []
    for feat in feature_cols:
        codex = feature_df.loc[feature_df["Platform"] == "CODEX", feat].dropna()
        pc = feature_df.loc[feature_df["Platform"] == "Phenocycler", feat].dropna()
        if len(codex) < 1 or len(pc) < 1:
            results.append({"Feature": feat, "U": np.nan, "p": np.nan,
                            "rank_biserial": np.nan, "n_CODEX": len(codex), "n_PC": len(pc)})
            continue
        u, p = mannwhitneyu(codex, pc, alternative="two-sided")
        r = 1 - 2 * u / (len(codex) * len(pc))
        results.append({"Feature": feat, "U": u, "p": p,
                        "rank_biserial": r, "n_CODEX": len(codex), "n_PC": len(pc)})
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Platform residualization (H11)
# ---------------------------------------------------------------------------
def residualize_platform(feature_df: pd.DataFrame,
                         feature_cols: list) -> pd.DataFrame:
    """Remove platform effect from features via OLS residuals.

    For each feature: fit y ~ platform_indicator, return residuals + grand mean.
    Preserves genotype variance while removing systematic platform offset.
    """
    out = feature_df.copy()
    platform_ind = (out["Platform"] == "CODEX").astype(float).values
    X = np.column_stack([np.ones(len(platform_ind)), platform_ind])
    for col in feature_cols:
        y = out[col].values.astype(float)
        mask = np.isfinite(y)
        if mask.sum() < 3:
            continue
        beta, _, _, _ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
        residuals = y.copy()
        residuals[mask] = y[mask] - X[mask] @ beta + beta[0]  # residual + intercept
        out[col] = residuals
    return out


# ---------------------------------------------------------------------------
# Cell region assignment by signed distance
# ---------------------------------------------------------------------------
def assign_region_by_distance(dist_follicle, dist_pals, dist_redpulp):
    """Assign cells to the region they are deepest inside. SUPERSEDED.

    Reads QuPath's pre-computed "Signed distance to annotation X µm" columns,
    which were measured against the pre-rebuild annotation geometry and so no
    longer agree with the region layer — most visibly for 1901HBMP004, whose
    exported geometry was stale. Use :func:`assign_objects_to_polygons` against
    :func:`load_region_layer` instead. Kept so the old labelling can be
    reproduced for comparison.

    Takes three arrays of signed distances (negative = inside).
    Returns array of region labels: Follicle, PALS, RedPulp, or Other.
    """
    regions = np.array(["Follicle", "PALS", "RedPulp"])
    dists = np.column_stack([dist_follicle, dist_pals, dist_redpulp])
    min_idx = np.argmin(dists, axis=1)
    min_val = np.min(dists, axis=1)
    result = regions[min_idx]
    result[min_val >= 0] = "Other"
    return result


def signed_distance_to_regions(coords, polys, boundary_step=8.0):
    """Distance from each point to the nearest polygon boundary, negative inside.

    The region-layer replacement for QuPath's "Signed distance to annotation X
    µm" columns, so distance-based readouts (H15/H16 cell-type distance
    profiles) can be measured against the same geometry everything else uses.

    Boundaries are sampled every ``boundary_step`` µm and queried with a
    cKDTree, which is what makes this tractable over millions of cells; the
    cost is that a reported distance can be off by up to half a step near a
    highly convoluted edge. Containment is exact.

    Returns NaN everywhere if ``polys`` is empty — a region a donor does not
    have is unmeasured, not infinitely far away.
    """
    import shapely
    from shapely.ops import unary_union

    coords = np.asarray(coords, dtype=float)
    if not polys or len(coords) == 0:
        return np.full(len(coords), np.nan)

    pts = []
    for p in polys:
        for ring in [p.exterior, *p.interiors]:
            n = max(int(ring.length / boundary_step), 8)
            t = np.linspace(0, ring.length, n, endpoint=False)
            pts.append(np.array([ring.interpolate(s).coords[0] for s in t]))
    if not pts:
        return np.full(len(coords), np.nan)

    dist, _ = cKDTree(np.vstack(pts)).query(coords, k=1)
    inside = shapely.contains_xy(unary_union(polys), coords[:, 0], coords[:, 1])
    return np.where(inside, -dist, dist)


# ---------------------------------------------------------------------------
# KDE-based region boundary utilities (H14)
# ---------------------------------------------------------------------------
def compute_marker_kde(coords, x_range, y_range, grid_spacing=10, bandwidth=75):
    """Bin cell coordinates onto a grid and smooth with Gaussian filter.

    Parameters
    ----------
    coords : (N, 2) array of (x, y) positions in µm
    x_range : (x_min, x_max)
    y_range : (y_min, y_max)
    grid_spacing : bin size in µm (default 10)
    bandwidth : Gaussian sigma in µm (default 75)

    Returns
    -------
    density : 2-D array (ny, nx) of smoothed cell density
    x_edges : 1-D array of bin edges along x
    y_edges : 1-D array of bin edges along y
    """
    from scipy.ndimage import gaussian_filter

    x_edges = np.arange(x_range[0], x_range[1] + grid_spacing, grid_spacing)
    y_edges = np.arange(y_range[0], y_range[1] + grid_spacing, grid_spacing)
    counts, _, _ = np.histogram2d(
        coords[:, 0], coords[:, 1], bins=[x_edges, y_edges])
    # histogram2d returns (nx, ny); transpose to (ny, nx) for image convention
    counts = counts.T
    sigma = bandwidth / grid_spacing
    density = gaussian_filter(counts.astype(float), sigma=sigma)
    return density, x_edges, y_edges


def polygons_at_level(density, level, grid_spacing, x_min, y_min, min_area=8400):
    """Contour a density grid at an *absolute* level and return Shapely polygons.

    The level-based sibling of :func:`extract_region_polygons`, which contours at
    a fraction of the grid maximum. Marker regions (Follicle, PALS) are naturally
    expressed relative to their own peak; tissue extent is not — its threshold
    comes from Otsu on the density itself — so both spellings are needed.

    Parameters
    ----------
    density : 2-D smoothed density array (ny, nx)
    level : absolute density value to contour at
    grid_spacing : µm per grid cell
    x_min, y_min : origin offset in µm
    min_area : discard polygons smaller than this (µm², default 8400)

    Returns
    -------
    list of shapely.geometry.Polygon in µm coordinates
    """
    from skimage.measure import find_contours
    from shapely.geometry import Polygon
    from shapely.validation import make_valid

    if level <= 0:
        return []

    contours = find_contours(density, level)
    polygons = []
    for contour in contours:
        if len(contour) < 4:
            continue
        # contour is (row, col) → convert to (x_µm, y_µm)
        coords_um = np.column_stack([
            contour[:, 1] * grid_spacing + x_min,
            contour[:, 0] * grid_spacing + y_min,
        ])
        poly = Polygon(coords_um)
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty or poly.geom_type not in ('Polygon', 'MultiPolygon'):
            continue
        if poly.geom_type == 'MultiPolygon':
            for p in poly.geoms:
                if p.area >= min_area:
                    polygons.append(p)
        elif poly.area >= min_area:
            polygons.append(poly)
    return polygons


def extract_region_polygons(density, threshold_frac, grid_spacing, x_min, y_min,
                            min_area=8400):
    """Extract contours from density grid and return Shapely polygons.

    Parameters
    ----------
    density : 2-D smoothed density array (ny, nx)
    threshold_frac : fraction of max density for contouring (e.g. 0.10)
    grid_spacing : µm per grid cell
    x_min, y_min : origin offset in µm
    min_area : discard polygons smaller than this (µm², default 8400)

    Returns
    -------
    list of shapely.geometry.Polygon in µm coordinates
    """
    return polygons_at_level(density, threshold_frac * density.max(),
                             grid_spacing, x_min, y_min, min_area=min_area)


def assign_objects_to_polygons(coords, polygon_dict):
    """Assign point coordinates to regions via spatial index.

    Parameters
    ----------
    coords : (N, 2) array of (x, y) in µm
    polygon_dict : {region_name: list_of_polygons}

    Returns
    -------
    labels : array of region name strings (len N), "Unassigned" if outside all
    """
    import shapely
    from shapely import STRtree

    n = len(coords)
    labels = np.full(n, "Unassigned", dtype=object)

    # Build flat list of (polygon, region_name) and STRtree
    all_polys = []
    poly_labels = []
    for region, polys in polygon_dict.items():
        for p in polys:
            all_polys.append(p)
            poly_labels.append(region)

    if not all_polys:
        return labels

    tree = STRtree(all_polys)

    # Vectorized query: predicate="within" means point.within(polygon)
    # (shapely 2.x evaluates query_geom.predicate(tree_geom))
    #
    # Both steps have to stay vectorised. Building the points with a Python
    # comprehension over 500k coordinates, and writing the matches back one at a
    # time, cost 27s per 500k-row chunk against 0.8s for the tree query itself —
    # so ~97% of the time was interpreter overhead, and it was paid once per
    # chunk in H11 and once per notebook everywhere else.
    coords = np.asarray(coords, dtype=float)
    points = shapely.points(coords[:, 0], coords[:, 1])
    pt_idx, tree_idx = tree.query(points, predicate="within")
    # Fancy assignment keeps last-writer-wins on duplicate point indices, which
    # is what the loop did. The layer is disjoint, so this only decides overlaps
    # that should not exist.
    labels[pt_idx] = np.asarray(poly_labels, dtype=object)[tree_idx]

    return labels


# ---------------------------------------------------------------------------
# Canonical cell-derived region layer
# ---------------------------------------------------------------------------
#
# Regions come from the single-cell data, not from QuPath annotations. The
# QuPath exports failed repeatedly (1,014 vessels orphaned at the image root, a
# CODEX Red_Pulp/Sinusoid area double-count, and a segmentation version skew
# where 1901HBMP004 carries 54 Follicle annotations in GeoJSON but 35 in
# AnnotationsFinal.csv), so one method now generates every region.
#
# Two things stay QuPath-derived because the cell data cannot supply them:
# Trabeculae (no marker in the panel identifies it) and the SmallVessel objects
# themselves (there is no vessel object in the cell tables) — only the region a
# vessel is attributed to changes.

REGION_ORDER = ["Follicle", "PALS", "Perifollicular", "RedPulp", "Trabeculae"]

# Tissue extent. Coarser grid and wider bandwidth than the marker regions: this
# is a glass-vs-tissue decision, not a boundary, and the Otsu split wants a
# smooth unimodal-per-class histogram.
TISSUE_GRID_SPACING = 20
TISSUE_BANDWIDTH = 60
TISSUE_MIN_AREA = 50_000

# Marker regions — the settings validated in H14, kept for the legacy
# percentile path (marker_region_polygons); the canonical builder no longer
# uses them.
REGION_GRID_SPACING = 10
MIN_REGION_AREA = 8400
CD20_PERCENTILE, CD20_THRESHOLD_FRAC = 90, 0.10
CD3E_PERCENTILE, CD3E_THRESHOLD_FRAC = 95, 0.25

# White pulp is where either lymphoid marker's z-scored intensity field rises
# a standard deviation above its tissue mean; it then splits into Follicle
# (CD20-dominant) and PALS (CD3e-dominant). A fixed z cut is used rather than a
# per-image Otsu: Otsu on the combined field swings between 1.22 and 2.87
# across donors and at the high end erases the follicles outright
# (20_06_CC2B → 0% Follicle). z > 1.0 is the value at which both regions come
# back non-empty for every one of the 14 donors.
WHITE_PULP_Z = 1.0

# Minimum follicle area, from a 2-component GMM on log10 area of the
# cell-derived polygons — the same procedure H6 uses, but fitted to this
# method's own output rather than to the QuPath annotations. The QuPath-fitted
# boundary (7,945 µm², essentially the old 8,400 constant) is far too permissive
# here: the cell-derived polygons have their real-follicle mode at ~170,000 µm²,
# and 8,400 admits a tail of ~200 objects holding 42-100 CD20-high cells each
# against 1,146 for a median follicle. Those are the "small areas that are not
# actually follicles".
MIN_FOLLICLE_AREA = 57_000

# Minimum area for a *QuPath* Follicle object. Their areas are bimodal — a
# 2-component GMM on log10 area splits at ~8,000 µm2 (equivalently ~101 µm
# across), which is the filter H6 already applies. Below it lie ~50% of the
# objects carrying ~2.1% of the follicle area: fragments, not follicles.
# scripts/fit_area_floors.py re-derives this; it currently reports 8,490 µm2,
# and the 7% gap moves neither the object count nor the area dropped by more
# than a few tenths of a percent, so the original constant stands.
MIN_QUPATH_FOLLICLE_AREA = 7_945

# Minimum area for a QuPath PALS object. Unlike the follicles', this is NOT a
# GMM crossover, and an earlier comment here wrongly said it was. PALS areas
# are not cleanly bimodal, so the 2-component fit splits the body of the
# distribution rather than separating fragments from sheaths: it lands at
# 4,141 µm2 and would discard 67% of the objects and 9.5% of the PALS area.
# That contradicts the validation, which requires PALS to reproduce the
# annotation CSV at 0.95-1.00x. 612 µm2 is instead a conservative geometric
# floor — roughly one arteriole's cross-section, a circle ~28 µm across — and
# removes 0.2% of the objects and 0.01% of the area. Run
# scripts/fit_area_floors.py to see both numbers.
MIN_QUPATH_PALS_AREA = 612

# Anything smaller than a single cell is a geometry artifact, not a region.
# Intersecting with the tissue mask and differencing the perifollicular band
# shatters polygons into slivers: before this floor the layer held PALS and
# Perifollicular objects of exactly 0 µm2, with 10th percentiles of 7 and 9 µm2.
SLIVER_AREA = 100

# Perifollicular zone.
PFZ_MIN_WIDTH = 25          # µm — guaranteed annulus, so the zone is never empty
PFZ_MAX_WIDTH = 100         # µm — how far marker signal may extend it
PFZ_WIDTH_SWEEP = (50, 100, 150, 200)
PFZ_PERCENTILE = 85
PFZ_THRESHOLD_FRAC = 0.20
# Both arms are restricted to MARKERS_USABLE_ALL_DONORS. CD35 and Lyve1 are the
# textbook choices here and are deliberately absent: each is entirely missing
# from all nine Phenocycler donors, so including either would define the
# perifollicular zone from different markers on different platforms — the exact
# confound this project exists to avoid. CD11c, the classic marginal-zone
# marker, is out for the same reason (all-NaN in HBMP007).
PFZ_MARKERS = ("CD21", "CD163")                  # FDC / perifollicular macrophage
PFZ_VASCULAR_MARKERS = ("CD31", "SMA", "CD8")    # sheathed capillaries, CD8+ sinusoidal endothelium

# QC. A marker that carries no regional structure produces a "region" covering
# most of the section. 20_07_CC1A's CD3e does exactly this: its top 5% of cells
# are spatially uniform, yielding a 14.2 mm² PALS in a 13.4 mm² section against
# a QuPath value of 0.119 mm². Such a region is dropped and recorded, never
# emitted as a number.
MAX_REGION_TISSUE_FRAC = 0.35


def adaptive_bandwidth(coords, k=5, scale=4, bw_min=25, bw_max=100):
    """Per-image KDE bandwidth from the k-th nearest-neighbour distance.

    Parameters
    ----------
    coords : (N, 2) array of marker-positive cell positions in µm
    k : which nearest neighbour to measure (default 5)
    scale : multiplier on the median k-NN distance (default 4)
    bw_min, bw_max : clamp range in µm

    Returns
    -------
    float bandwidth in µm
    """
    if len(coords) < k + 1:
        return (bw_min + bw_max) / 2
    dists, _ = cKDTree(coords).query(coords, k=k + 1)   # k+1: index 0 is self
    return float(np.clip(np.median(dists[:, k]) * scale, bw_min, bw_max))


def tissue_extent(coords, grid_spacing=TISSUE_GRID_SPACING,
                  bandwidth=TISSUE_BANDWIDTH, min_area=TISSUE_MIN_AREA):
    """Tissue footprint from cell density alone, via Otsu on the smoothed grid.

    Independent of any QuPath annotation. Validated against QuPath's total
    annotated area across all 16 images: ratio 0.92–1.12, median 1.00.

    Parameters
    ----------
    coords : (N, 2) array of *all* cell centroids in µm
    grid_spacing, bandwidth : KDE grid and smoothing, µm
    min_area : discard fragments below this (µm²)

    Returns
    -------
    list of shapely Polygons in µm
    """
    from skimage.filters import threshold_otsu

    if len(coords) < 10:
        return []
    pad = 2 * bandwidth
    x_range = (coords[:, 0].min() - pad, coords[:, 0].max() + pad)
    y_range = (coords[:, 1].min() - pad, coords[:, 1].max() + pad)
    density, xe, ye = compute_marker_kde(coords, x_range, y_range,
                                         grid_spacing=grid_spacing,
                                         bandwidth=bandwidth)
    occupied = density[density > 1e-9]
    if occupied.size == 0:
        return []
    level = threshold_otsu(occupied)
    return polygons_at_level(density, level, grid_spacing, xe[0], ye[0],
                             min_area=min_area)


def marker_enrichment_fields(coords, markers, grid_spacing=TISSUE_GRID_SPACING,
                             bandwidth=TISSUE_BANDWIDTH):
    """Smoothed per-marker intensity fields, z-scored within tissue.

    For each marker: local mean of log1p(intensity) on a smoothed grid, then
    standardised against its own mean and SD over tissue. Two markers are then
    directly comparable in units of their own spatial variability, which is what
    lets Follicle and PALS be split by whichever dominates locally.

    Why intensity and not a positivity gate. Both cohorts are
    background-subtracted with negatives clipped to zero — no marker anywhere
    has a negative value, and the exactly-zero fraction runs from 0.1% in the
    Phenocycler donors to 86% (CD3e, 20_06_CC2B) in the CODEX ones. A fixed
    percentile gate is platform-dependent by construction: "top 10%" means
    something different when 86% of cells sit on the floor than when 0.1% do.

    Rank/ECDF normalisation was tried and rejected: discarding magnitude makes
    spatial contrast a function of the *clipped fraction*, so the marker with
    more zeros wins the dominance test. It inverted HDL073, a visibly
    follicle-dominated section, to PALS-dominant. Working on log intensity keeps
    magnitude, and the z-score removes the scale difference without rewarding
    clipping.

    Returns
    -------
    (fields, density, x_edges, y_edges) — fields is {marker: 2-D z-score array}
    """
    from scipy.ndimage import gaussian_filter
    from skimage.filters import threshold_otsu

    coords = np.asarray(coords, dtype=float)
    pad = 2 * bandwidth
    xe = np.arange(coords[:, 0].min() - pad,
                   coords[:, 0].max() + pad + grid_spacing, grid_spacing)
    ye = np.arange(coords[:, 1].min() - pad,
                   coords[:, 1].max() + pad + grid_spacing, grid_spacing)
    sigma = bandwidth / grid_spacing

    counts, _, _ = np.histogram2d(coords[:, 0], coords[:, 1], bins=[xe, ye])
    density = gaussian_filter(counts.T.astype(float), sigma)
    occupied = density[density > 1e-9]
    mask = (density > threshold_otsu(occupied)) if occupied.size else density > 0
    if not mask.any():
        mask = density > 0

    fields = {}
    for name, values in markers.items():
        v = np.nan_to_num(np.log1p(np.clip(np.asarray(values, dtype=float),
                                           0, None)))
        num, _, _ = np.histogram2d(coords[:, 0], coords[:, 1], bins=[xe, ye],
                                   weights=v)
        local = gaussian_filter(num.T.astype(float), sigma) / np.maximum(density, 1e-9)
        mu, sd = local[mask].mean(), local[mask].std()
        fields[name] = (local - mu) / (sd if sd > 0 else 1.0)
    return fields, density, xe, ye


def marker_region_polygons(all_coords, values, percentile, threshold_frac,
                           grid_spacing=REGION_GRID_SPACING,
                           min_area=MIN_REGION_AREA):
    """Contour the density of the top-`percentile` cells for one marker.

    The percentile threshold is deliberate: thresholding on the intensity
    distribution instead (Otsu on log1p) calls 40–73% of cells positive and
    inflates every region several-fold.

    Returns
    -------
    (polygons, info) where info carries the intensity threshold, positive-cell
    count and bandwidth actually used.
    """
    values = np.asarray(values, dtype=float)
    info = {"threshold": np.nan, "n_positive": 0, "bandwidth": np.nan}
    if not np.isfinite(values).any():
        return [], info

    thresh = np.nanpercentile(values, percentile)
    hi = all_coords[values > thresh]
    info["threshold"], info["n_positive"] = float(thresh), int(len(hi))
    if len(hi) < 10:
        return [], info

    bw = adaptive_bandwidth(hi)
    info["bandwidth"] = bw
    pad = 2 * bw
    x_range = (all_coords[:, 0].min() - pad, all_coords[:, 0].max() + pad)
    y_range = (all_coords[:, 1].min() - pad, all_coords[:, 1].max() + pad)
    density, xe, ye = compute_marker_kde(hi, x_range, y_range,
                                         grid_spacing=grid_spacing,
                                         bandwidth=bw)
    polys = extract_region_polygons(density, threshold_frac, grid_spacing,
                                    xe[0], ye[0], min_area=min_area)
    return polys, info


def perifollicular_score(markers):
    """Composite perifollicular signal per cell.

    Two arms averaged after z-scoring: the marginal-zone/FDC/macrophage
    identity markers (CD21/CD35/CD163) and the vascular signature of the human
    perifollicular zone (CD31/Lyve1/SMA/CD8 — sheathed capillaries and
    CD8+ sinusoidal endothelium). Arms are averaged separately first so a
    four-marker vascular arm does not outvote a three-marker identity arm.

    Parameters
    ----------
    markers : {name: (N,) array} of per-cell mean intensities

    Returns
    -------
    (N,) array of composite score, or None if neither arm has any marker
    """
    def _arm(names):
        cols = []
        for n in names:
            v = markers.get(n)
            if v is None:
                continue
            v = np.asarray(v, dtype=float)
            if not np.isfinite(v).any():
                continue
            lv = np.log1p(np.clip(v, 0, None))
            sd = np.nanstd(lv)
            if not np.isfinite(sd) or sd == 0:
                continue
            cols.append((lv - np.nanmean(lv)) / sd)
        return np.nanmean(np.column_stack(cols), axis=1) if cols else None

    arms = [a for a in (_arm(PFZ_MARKERS), _arm(PFZ_VASCULAR_MARKERS)) if a is not None]
    if not arms:
        return None
    return np.nanmean(np.column_stack(arms), axis=1)


def perifollicular_zone(follicle_polys, exclude_polys, tissue_polys, coords,
                        score=None, min_width=PFZ_MIN_WIDTH,
                        max_width=PFZ_MAX_WIDTH, percentile=PFZ_PERCENTILE,
                        threshold_frac=PFZ_THRESHOLD_FRAC,
                        grid_spacing=REGION_GRID_SPACING,
                        min_area=MIN_REGION_AREA):
    """Perifollicular band around the cell-derived follicles.

    A `min_width` annulus is guaranteed so the zone is never empty next to a
    follicle; marker signal then extends it outward as far as `max_width`. The
    geometry is the scaffold and the markers choose the outer boundary — a
    purely marker-driven contour would collapse onto the follicle, because
    CD21/CD35 label the follicular dendritic network largely *inside* it.

    Parameters
    ----------
    follicle_polys : list of Polygons defining the follicles
    exclude_polys : list of Polygons to subtract (PALS, Trabeculae)
    tissue_polys : list of Polygons bounding the tissue
    coords : (N, 2) all cell positions in µm
    score : (N,) composite perifollicular score, or None for a fixed-width band

    Returns
    -------
    list of shapely Polygons in µm
    """
    from shapely.ops import unary_union

    if not follicle_polys:
        return []
    fol = unary_union(follicle_polys)
    tissue = unary_union(tissue_polys) if tissue_polys else None
    excl = unary_union(exclude_polys) if exclude_polys else None

    def _band(width):
        b = fol.buffer(width).difference(fol)
        if excl is not None:
            b = b.difference(excl)
        if tissue is not None:
            b = b.intersection(tissue)
        return b

    guaranteed = _band(min_width)
    scaffold = _band(max_width)
    if scaffold.is_empty:
        return _as_polygon_list(guaranteed, min_area=0.0)

    zone = guaranteed
    if score is not None and np.isfinite(score).any():
        hi = coords[score > np.nanpercentile(score, percentile)]
        if len(hi) >= 10:
            bw = adaptive_bandwidth(hi)
            pad = 2 * bw
            x_range = (coords[:, 0].min() - pad, coords[:, 0].max() + pad)
            y_range = (coords[:, 1].min() - pad, coords[:, 1].max() + pad)
            density, xe, ye = compute_marker_kde(hi, x_range, y_range,
                                                 grid_spacing=grid_spacing,
                                                 bandwidth=bw)
            signal = extract_region_polygons(density, threshold_frac,
                                             grid_spacing, xe[0], ye[0],
                                             min_area=min_area)
            if signal:
                guided = unary_union(signal).intersection(scaffold)
                zone = unary_union([zone, guided])

    return _as_polygon_list(zone, min_area=0.0)


def _as_polygon_list(geom, min_area=MIN_REGION_AREA):
    """Flatten a geometry to a list of non-trivial Polygons."""
    if geom is None or geom.is_empty:
        return []
    parts = geom.geoms if geom.geom_type.startswith("Multi") or \
        geom.geom_type == "GeometryCollection" else [geom]
    return [p for p in parts if p.geom_type == "Polygon" and p.area >= min_area]


def build_cell_derived_regions(coords, markers, trabeculae=None,
                               max_region_frac=MAX_REGION_TISSUE_FRAC,
                               pfz_max_width=PFZ_MAX_WIDTH):
    """Build the full disjoint region set for one image from its cells.

    Regions are constructed disjoint and in priority order, so
    :func:`assign_objects_to_polygons` gives each object exactly one label:
    Follicle, then PALS (minus Follicle), then Perifollicular (minus both),
    then Trabeculae, then RedPulp as the remaining tissue.

    Parameters
    ----------
    coords : (N, 2) array of all cell centroids in µm
    markers : {name: (N,) array} of per-cell mean intensities; needs at least
        CD20 and CD3e, plus the PFZ markers to guide the perifollicular border
    trabeculae : list of Polygons from QuPath, or None
    max_region_frac : drop (and record) any marker region covering more than
        this fraction of tissue — see `MAX_REGION_TISSUE_FRAC`
    pfz_max_width : outer limit of the perifollicular band, µm

    Returns
    -------
    (regions, qc) — regions is {name: [Polygon]}, qc is a dict of diagnostics
    including `tissue_mm2`, per-region areas, and `dropped` for regions that
    failed the coverage check.
    """
    from shapely.ops import unary_union

    coords = np.asarray(coords, dtype=float)
    qc = {"n_cells": len(coords), "flagged": [], "warnings": []}

    tissue = tissue_extent(coords)
    tissue_area = unary_union(tissue).area if tissue else 0.0
    qc["tissue_mm2"] = tissue_area / 1e6
    if tissue_area <= 0:
        qc["warnings"].append("empty tissue extent")
        return {r: [] for r in REGION_ORDER}, qc

    def _check(name, polys):
        """Flag a region that covers an unusual amount of tissue — never drop it.

        This used to zero the region, which is how HBMP007 and HDL098 ended up
        with no PALS at all. Every spleen has a T zone, so an empty PALS is the
        method reporting its own failure, not a property of the donor. The
        fraction and a warning are recorded instead, so an outlier stays visible
        without a fabricated zero entering the analysis.
        """
        if not polys:
            qc["warnings"].append(f"{name} is empty")
            return polys
        frac = unary_union(polys).area / tissue_area
        qc[f"{name}_tissue_frac"] = frac
        if frac > max_region_frac:
            qc["flagged"].append(name)
            qc["warnings"].append(
                f"{name} covers {frac:.0%} of tissue (> {max_region_frac:.0%}) "
                "— unusually large; kept but flagged")
        return polys

    tissue_u = unary_union(tissue)

    def _clip(polys):
        """Keep regions inside the tissue footprint — the KDE pad overshoots it."""
        return [q for p in polys for q in _as_polygon_list(p.intersection(tissue_u))]

    # White pulp, then the B/T split. Follicle and PALS are carved from one
    # mask by which marker dominates locally, so PALS cannot come back empty
    # while white pulp exists — every spleen has a T zone, and a method that
    # reports none for a donor is reporting its own failure.
    fields, density, xe, ye = marker_enrichment_fields(coords, {
        "CD20": markers.get("CD20", np.full(len(coords), np.nan)),
        "CD3e": markers.get("CD3e", np.full(len(coords), np.nan)),
    })
    e20, e3e = fields["CD20"], fields["CD3e"]
    wp_grid = np.maximum(e20, e3e) > WHITE_PULP_Z
    qc["white_pulp_frac"] = float(wp_grid.mean())

    def _from_grid(mask, min_area):
        if not mask.any():
            return []
        polys = polygons_at_level(mask.astype(float), 0.5, TISSUE_GRID_SPACING,
                                  xe[0], ye[0], min_area=min_area)
        return _clip(polys)

    fol = _from_grid(wp_grid & (e20 >= e3e), MIN_FOLLICLE_AREA)
    pals = _from_grid(wp_grid & (e3e > e20), MIN_REGION_AREA)
    fol = _check("Follicle", fol)
    pals = _check("PALS", pals)
    if fol and pals:                      # keep them disjoint
        fu = unary_union(fol)
        pals = [q for p in pals for q in _as_polygon_list(p.difference(fu))]
    qc["Follicle_info"] = {"threshold": WHITE_PULP_Z,
                           "n_positive": int((wp_grid & (e20 >= e3e)).sum()),
                           "bandwidth": TISSUE_BANDWIDTH}
    qc["PALS_info"] = {"threshold": WHITE_PULP_Z,
                       "n_positive": int((wp_grid & (e3e > e20)).sum()),
                       "bandwidth": TISSUE_BANDWIDTH}

    pfz = perifollicular_zone(fol, pals + (trabeculae or []), tissue, coords,
                              score=perifollicular_score(markers),
                              max_width=pfz_max_width)

    trab = list(trabeculae) if trabeculae else []
    claimed = fol + pals + pfz + trab
    red = _as_polygon_list(
        unary_union(tissue).difference(unary_union(claimed))) if claimed else tissue

    regions = {"Follicle": fol, "PALS": pals, "Perifollicular": pfz,
               "RedPulp": red, "Trabeculae": trab}
    for name, polys in regions.items():
        qc[f"{name}_mm2"] = (unary_union(polys).area / 1e6) if polys else 0.0
        qc[f"{name}_n"] = len(polys)
    return regions, qc
