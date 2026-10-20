"""PhenoCycler inputs for the QC figures, loaded once per donor into plain tables.

Everything here joins on ``object_id`` -- the four run tables (``cells``, ``geometry_qc``,
``expression``, ``redsea``) are NOT in the same row order (checked 2026-09-05: none of the three
pairs align positionally), so a positional join would be silently wrong.  Two derivations that
are easy to get wrong are done here once, with their reasons:

* **Nucleus area** is ``nucleus_cell_area_ratio x raster_area_ratio x cell_area`` from
  ``geometry_qc`` and ``cells`` (``cell_area`` is QuPath's um^2).  ``geometry_qc`` carries no area
  column, and the ``redsea`` table's ``nucleus_area_px``/``cell_area_px`` are ZERO-FILLED for the
  ~13 % of cells that are spillover-context ineligible (``apply_spillover_context_eligibility``
  removes them from the raster), so reading areas from there would print 0 for a missing
  measurement.  Where the nucleus is missing (``missing_nucleus_geometry`` /
  ``nucleus_raster_empty``) the area is NaN -- Unavailable, never 0 -- and counted.
* **The donor area ceiling** is ``max_cell_area_median_multiple x median(cell_area over the
  provisional set)`` where provisional follows ``geometry_qc.evaluate_geometry_metrics`` exactly:
  cell and nucleus geometry present, finite area >= min, finite solidity >= min, duplicate losers
  removed.  The recomputed count above the ceiling is printed beside the run's own
  ``cell_area_above_max`` count so the reproduction is checkable.

The REDSEA window (P05) maps mask labels to cells through the GeoJSON stream order
(``label = feature_index + 1``, the ``window_mask`` invariant).  The id list is streamed once with
the pipeline's own iterator, cached, and verified against the QuPath manifest's
``object_id_sha256`` before it is used; centroid matching is never a fallback.
"""
from __future__ import annotations

import json
import pathlib
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
import style  # noqa: E402
from tiff_window import read_window  # noqa: E402  (scripts/common)
from phenocycler.artifacts import hash_identifiers, iter_qupath_geojson_features  # noqa: E402
from phenocycler.redsea import qptiff_channels  # noqa: E402

#: Reasons under which the nucleus measurement does not exist for the cell.
NUCLEUS_MISSING_REASONS = ("missing_nucleus_geometry", "nucleus_raster_empty", "missing_cell_geometry")
GEOMETRY_COLUMNS = ["object_id", "analysis_eligible", "estimation_eligible", "spillover_context_eligible",
                    "analysis_reasons", "estimation_reasons", "spillover_context_reasons", "duplicate_group",
                    "duplicate_survivor", "raster_area_ratio", "nucleus_cell_area_ratio", "analysis_reason",
                    "estimation_reason", "spillover_context_reason"]
CELL_COLUMNS = ["object_id", "cell_region", "cell_area", "cell_solidity", "X_centroid", "Y_centroid"]


@dataclass
class PhenoTables:
    """One donor's geometry-QC view: every run row, joined on object_id."""
    donor: str
    status: str
    table: pd.DataFrame          # cells x geometry_qc x expression flags, one row per run row
    thresholds: dict             # run.json config.canonical.geometry_qc
    n_geojson_features: int
    n_ingest_rejects: int
    inputs: list = field(default_factory=list)

    @property
    def n(self) -> int:
        return int(len(self.table))


def _read(path: pathlib.Path, columns=None) -> pd.DataFrame:
    df = pq.read_table(path, columns=columns).to_pandas()
    if "object_id" in df.columns:
        df["object_id"] = df["object_id"].astype(str)
    return df


def has_reason(reasons: pd.Series, *tokens: str) -> np.ndarray:
    """Token membership in a ';'-joined reason string, never a substring match."""
    parts = reasons.fillna("").astype(str).str.split(";")
    want = set(tokens)
    return np.fromiter((bool(want & set(p)) for p in parts), dtype=bool, count=len(parts))


def nucleus_area_um2(t: pd.DataFrame) -> tuple:
    """(area, n_unavailable): NaN where the nucleus measurement does not exist."""
    missing = has_reason(t["analysis_reasons"], *NUCLEUS_MISSING_REASONS) \
        | has_reason(t["spillover_context_reasons"], *NUCLEUS_MISSING_REASONS)
    area = t["nucleus_cell_area_ratio"].to_numpy(float) * t["raster_area_ratio"].to_numpy(float) \
        * t["cell_area"].to_numpy(float)
    area = np.where(missing, np.nan, area)
    return area, int((~np.isfinite(area)).sum())


def provisional_mask(t: pd.DataFrame, thr: dict) -> np.ndarray:
    """``evaluate_geometry_metrics``'s provisional set, from the recorded decisions."""
    area = t["cell_area"].to_numpy(float)
    sol = t["cell_solidity"].to_numpy(float)
    nucleus_missing = has_reason(t["analysis_reasons"], "missing_nucleus_geometry", "missing_cell_geometry")
    dup_loser = (t["duplicate_group"].fillna("").astype(str) != "") & (~t["duplicate_survivor"].fillna(False).astype(bool))
    return (~nucleus_missing & np.isfinite(area) & (area > 0) & (area >= float(thr["min_cell_area_um2"]))
            & np.isfinite(sol) & (sol >= 0) & (sol <= 1) & (sol >= float(thr["min_cell_solidity"])) & ~dup_loser.to_numpy())


def recompute_ceiling(t: pd.DataFrame, thr: dict) -> dict:
    prov = provisional_mask(t, thr)
    med = float(np.median(t.loc[prov, "cell_area"].to_numpy(float)))
    ceiling = float(thr["max_cell_area_median_multiple"]) * med
    area = t["cell_area"].to_numpy(float)
    n_recomputed = int((np.isfinite(area) & (area > ceiling)).sum())
    n_flagged = int(has_reason(t["analysis_reasons"], "cell_area_above_max").sum())
    return {"n_provisional": int(prov.sum()), "median_cell_area_um2_provisional": med,
            "ceiling_um2": ceiling, "multiple": float(thr["max_cell_area_median_multiple"]),
            "n_above_ceiling_recomputed": n_recomputed, "n_flagged_cell_area_above_max": n_flagged,
            "reproduces_flag": n_recomputed == n_flagged}


def load_tables(donor: str) -> PhenoTables:
    cells = _read(pins.cells_parquet(donor), CELL_COLUMNS)
    geom = _read(pins.geometry_parquet(donor), GEOMETRY_COLUMNS)
    expr = _read(pins.expression_parquet(donor), ["object_id", "qc_analysis_eligible", "qc_estimation_eligible"])
    if not (len(cells) == len(geom) == len(expr)):
        raise RuntimeError(f"donor {donor}: run tables differ in length {len(cells)} / {len(geom)} / {len(expr)}")
    t = cells.merge(geom, on="object_id", how="inner", validate="one_to_one") \
             .merge(expr, on="object_id", how="inner", validate="one_to_one")
    if len(t) != len(cells):
        raise RuntimeError(f"donor {donor}: object_id join lost rows ({len(t)} of {len(cells)})")
    img = pins.manifest_image(donor)
    n_feat = int(img["cell_geometry"]["feature_count"])
    n_rej = 0
    if pins.INGEST_REJECTS.exists():
        rej = pq.read_table(pins.INGEST_REJECTS).to_pandas()
        if "donor_id" in rej.columns:
            n_rej = int((rej["donor_id"].astype(str) == str(donor)).sum())
    inputs = [pins.cells_parquet(donor), pins.geometry_parquet(donor), pins.expression_parquet(donor),
              pins.RUN_JSON, pins.QUPATH_MANIFEST] + ([pins.INGEST_REJECTS] if pins.INGEST_REJECTS.exists() else [])
    return PhenoTables(donor=str(donor), status=pins.donor_status(donor), table=t,
                       thresholds=pins.geometry_thresholds(), n_geojson_features=n_feat, n_ingest_rejects=n_rej,
                       inputs=inputs)


# --------------------------------------------------------------------------- Step-2 aligned values (P04)
@dataclass
class AlignedMarker:
    marker: str
    values: np.ndarray           # analysis-eligible aligned values (NaN dropped), or raw where not aligned
    mode_aligned: float          # Step-2 histogram mode on the aligned scale; NaN when not derivable
    scale: float                 # the donor's Step-2 multiplicative scale (1.0 = untouched)
    provenance: str              # shifts.json provenance, or a reason string
    nucleus_measured: bool
    n_eligible: int
    n_finite: int
    n_zero: int
    available: bool = True
    unavailable_reason: str = ""
    registry_kind: str = ""


def shifts() -> dict:
    return json.loads(pins.SHIFTS.read_text())


def registry_kinds() -> dict:
    reg = json.loads(pins.REGISTRY.read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    return {m["name"]: m.get("kind", "") for m in reg}


def histogram_mode_log(values_raw: np.ndarray, lo: float, hi: float, num_bins: int = pins.UNIFORM_NUM_BINS) -> float:
    """UniFORM's landmark as ``uniform_adapter.background_mode_spread`` computes it: log of the
    values >= 1, ``num_bins`` bins on the cohort's global range, tallest bin skipping bin 0."""
    v = np.asarray(values_raw, dtype=float)
    v = v[np.isfinite(v) & (v >= 1.0)]
    if v.size == 0 or not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        return float("nan")
    c = np.histogram(np.log(v), bins=num_bins, range=(lo, hi))[0]
    binw = (hi - lo) / (num_bins - 1)
    return float(lo + (int(np.argmax(c[1:])) + 1) * binw)


def load_aligned(donor: str, markers: list | None = None) -> tuple:
    """Every acquired marker on the aligned scale with its Step-2 mode; absent markers flagged.

    Returns (list[AlignedMarker], inputs).  A marker the donor never acquired (b_Catenin1 on both
    donors) is ``available=False`` with reason ``not acquired``; a registry-excluded marker (IAPP)
    is drawn from the raw expression table with scale 1 and says so -- Step 2 never aligned it.
    """
    acquired = pins.acquired_markers(donor)
    markers = list(markers) if markers is not None else acquired + [m for m in ("b_Catenin1",) if m not in acquired]
    sh = shifts()["shifts"]
    kinds = registry_kinds()
    uni_cols = set(pq.read_schema(pins.uniform_parquet(donor)).names)
    have_uni = [m for m in markers if m in uni_cols]
    uni = _read(pins.uniform_parquet(donor), ["object_id", "qc_analysis_eligible"] + have_uni)
    elig = uni["qc_analysis_eligible"].to_numpy(bool)
    need_raw = [m for m in markers if m not in uni_cols and m in acquired]
    raw = None
    if need_raw:
        expr_cols = set(pq.read_schema(pins.expression_parquet(donor)).names)
        need_raw = [m for m in need_raw if m in expr_cols]
        raw = _read(pins.expression_parquet(donor), ["object_id", "qc_analysis_eligible"] + need_raw)
    out = []
    for m in markers:
        nucleus = m in pins.NUCLEUS_MEASURED
        if m not in acquired:
            out.append(AlignedMarker(m, np.zeros(0), float("nan"), float("nan"), "not_acquired", nucleus, 0, 0, 0,
                                     available=False, unavailable_reason="not acquired on this donor",
                                     registry_kind=kinds.get(m, "")))
            continue
        if m in uni_cols:
            v = uni.loc[elig, m].to_numpy(float)
            rec = sh.get(m, {})
            scale = float(rec.get("scale", {}).get(str(donor), 1.0))
            prov = str(rec.get("provenance", "no_shift_record"))
            lo, hi = float(rec.get("global_min", np.nan)), float(rec.get("global_max", np.nan))
            fin = v[np.isfinite(v)]
            mode = histogram_mode_log(fin / scale if scale > 0 else fin, lo, hi)
            mode_aligned = float(np.exp(mode) * scale) if np.isfinite(mode) else float("nan")
            n_elig = int(elig.sum())
        else:
            v = raw.loc[raw["qc_analysis_eligible"].to_numpy(bool), m].to_numpy(float)
            fin = v[np.isfinite(v)]
            scale, prov = 1.0, f"not aligned in Step 2 (registry kind={kinds.get(m, '?')}); raw expression scale"
            lo = float(np.log(fin[fin >= 1].min())) if (fin >= 1).any() else np.nan
            hi = float(np.log(fin.max())) if fin.size else np.nan
            mode = histogram_mode_log(fin, lo, hi)
            mode_aligned = float(np.exp(mode)) if np.isfinite(mode) else float("nan")
            n_elig = int(raw["qc_analysis_eligible"].sum())
        out.append(AlignedMarker(m, fin, mode_aligned, scale, prov, nucleus, n_elig, int(fin.size),
                                 int((fin == 0).sum()), registry_kind=kinds.get(m, "")))
    inputs = [pins.uniform_parquet(donor), pins.SHIFTS, pins.REGISTRY, pins.QUPATH_MANIFEST]
    if raw is not None:
        inputs.append(pins.expression_parquet(donor))
    return out, inputs


# --------------------------------------------------------------------------- the frame and the REDSEA window (P05)
def frame_regions(donor: str) -> dict:
    doc = json.loads((pins.FRAME / "regions.json").read_text())
    return doc["donors"][str(donor)]


def first_islet_region(donor: str) -> tuple:
    """(index, region) of the first islet-anchored draw in the frozen frame."""
    regs = frame_regions(donor)["regions"]
    for i, r in enumerate(regs):
        if r.get("draw") == "islet":
            return i, r
    raise RuntimeError(f"donor {donor}: the frame has no islet-anchored region")


def frame_mask(donor: str, idx: int) -> np.ndarray:
    return np.load(pins.FRAME / "masks" / f"donor_id={donor}" / f"w{idx:02d}.npy")


def geojson_object_ids(donor: str, cache_dir: pathlib.Path | None = None) -> tuple:
    """The GeoJSON feature ids in stream order (``object_ids[label - 1]``), verified against the
    QuPath manifest's ``object_id_sha256`` and ``feature_count`` before use.  Cached because the
    stream takes about a minute per donor."""
    img = pins.manifest_image(donor)
    want_sha = str(img["cell_geometry"]["object_id_sha256"])
    want_n = int(img["cell_geometry"]["feature_count"])
    cache_dir = pathlib.Path(cache_dir) if cache_dir else pins.OUT / str(donor) / "pheno" / "_cache"
    cache = cache_dir / "geojson_object_ids.parquet"
    ids = None
    if cache.exists():
        meta = pq.read_metadata(cache).metadata or {}
        if meta.get(b"object_id_sha256", b"").decode() == want_sha:
            ids = pq.read_table(cache)["object_id"].to_numpy(zero_copy_only=False).astype(object)
    if ids is None:
        ids = np.asarray([str(ft.get("id")) for ft in iter_qupath_geojson_features(pins.geojson_path(donor))],
                         dtype=object)
        got = hash_identifiers(ids)
        if got != want_sha or len(ids) != want_n:
            raise RuntimeError(f"donor {donor}: GeoJSON id stream ({len(ids)} ids, sha {got[:12]}) does not match "
                               f"the QuPath manifest ({want_n}, {want_sha[:12]}); refusing to map mask labels")
        cache_dir.mkdir(parents=True, exist_ok=True)
        import pyarrow as pa
        tbl = pa.table({"object_id": pa.array(list(ids), type=pa.string())})
        tbl = tbl.replace_schema_metadata({"object_id_sha256": want_sha, "feature_count": str(want_n),
                                           "geojson": str(pins.geojson_path(donor)), "order": "GeoJSON stream"})
        pq.write_table(tbl, cache)
    if len(ids) != want_n:
        raise RuntimeError(f"donor {donor}: cached id list has {len(ids)} ids, manifest says {want_n}")
    return ids, cache


def plane_population(page, *, tile: int = 512, stride_tiles: int = 8) -> np.ndarray:
    """A section-wide, selection-independent pixel population: every ``stride_tiles``-th 512-px tile
    in both directions (about 1/64 of the slide), read chunk-wise.  The window's own range would
    rescale each panel to its contents; this one cannot move when the window changes."""
    H, W = page.shape[-2], page.shape[-1]
    step = tile * stride_tiles
    parts = [read_window(page, y, y + tile, x, x + tile).ravel()
             for y in range(0, max(H - tile, 1), step) for x in range(0, max(W - tile, 1), step)]
    return np.concatenate(parts) if parts else np.zeros(0)


@dataclass
class RedseaWindow:
    """One frame region for P05: raw planes, the label mask, and per-label raw / corrected values."""
    donor: str
    status: str
    region_idx: int
    region: dict
    window: dict                 # y0, x0, size (level-0 px)
    core_px: int
    px_um: float
    markers: tuple
    planes: dict                 # marker -> 2D uint16 window
    ranges: dict                 # marker -> DisplayRange over the section-wide lattice population
    mask: np.ndarray             # int32 global labels
    lut_raw: dict                # marker -> float array indexed by label (NaN = Unavailable)
    lut_corrected: dict          # marker -> float array indexed by label
    fill_ranges: dict            # marker -> DisplayRange over the donor's raw per-cell means (shared by raw/corrected)
    status_by_label: np.ndarray  # 0 background/unused, 1 eligible, 2 spillover-ineligible, 3 not in run (ingest reject)
    counts: dict
    inputs: list = field(default_factory=list)


def load_redsea_window(donor: str, markers: tuple = pins.REDSEA_PAIR) -> RedseaWindow:
    idx, region = first_islet_region(donor)
    mask = frame_mask(donor, idx)
    win = region["window"]
    frame_doc = json.loads((pins.FRAME / "regions.json").read_text())
    core_px = int(frame_doc["criteria"]["core_px"])
    ids, cache = geojson_object_ids(donor)
    labels = np.unique(mask[mask > 0])
    if labels.size and int(labels.max()) > len(ids):
        raise RuntimeError(f"donor {donor}: mask label {labels.max()} exceeds the GeoJSON id count {len(ids)}")
    oid = pd.DataFrame({"label": labels, "object_id": ids[labels - 1].astype(str)})

    raw = _read(pins.cells_parquet(donor), ["object_id", *markers])
    corr = _read(pins.expression_parquet(donor), ["object_id", *markers])
    geom = _read(pins.geometry_parquet(donor), ["object_id", "spillover_context_eligible", "spillover_context_reason"])
    j = oid.merge(raw.rename(columns={m: f"raw__{m}" for m in markers}), on="object_id", how="left") \
           .merge(corr.rename(columns={m: f"corr__{m}" for m in markers}), on="object_id", how="left") \
           .merge(geom, on="object_id", how="left")
    in_run = j["spillover_context_eligible"].notna().to_numpy()
    eligible = j["spillover_context_eligible"].fillna(False).astype(bool).to_numpy()
    n_lab = int(labels.max()) + 1 if labels.size else 1
    status = np.zeros(n_lab, dtype=np.int8)
    status[j["label"].to_numpy()] = np.where(~in_run, 3, np.where(eligible, 1, 2))
    lut_raw, lut_corr, fill_ranges = {}, {}, {}
    for m in markers:
        lr = np.full(n_lab, np.nan); lc = np.full(n_lab, np.nan)
        lr[j["label"].to_numpy()] = j[f"raw__{m}"].to_numpy(float)
        lc[j["label"].to_numpy()] = j[f"corr__{m}"].to_numpy(float)
        # corrected values of spillover-ineligible cells are NaN in expression already; keep them Unavailable
        lc[j["label"].to_numpy()[~eligible]] = np.nan
        lut_raw[m], lut_corr[m] = lr, lc
        # one white point per marker for raw AND corrected: the donor-wide raw QuPath means (selection-independent)
        fill_ranges[m] = style.unsigned_display_range(raw[m].to_numpy(float), ceiling=65535.0)

    tf, pages, names, shape, px_um = qptiff_channels(pins.qptiff_path(donor), 0)
    planes, ranges = {}, {}
    try:
        for m in markers:
            ch = pins.channel_index(donor, m)
            page = pages[ch]
            planes[m] = read_window(page, win["y0"], win["y0"] + win["size"], win["x0"], win["x0"] + win["size"])
            ranges[m] = style.unsigned_display_range(plane_population(page), ceiling=65535.0)
    finally:
        tf.close()
    counts = {"n_labels_in_window": int(labels.size), "n_in_run": int(in_run.sum()),
              "n_spillover_eligible": int(eligible.sum()), "n_spillover_ineligible": int((in_run & ~eligible).sum()),
              "n_not_in_run": int((~in_run).sum()),
              "ineligible_reasons": j.loc[in_run & ~eligible, "spillover_context_reason"].value_counts().to_dict(),
              "n_cells_donor": int(len(raw)), "n_spillover_eligible_donor": int(geom["spillover_context_eligible"].sum())}
    inputs = [pins.FRAME / "regions.json", pins.FRAME / "masks" / f"donor_id={donor}" / f"w{idx:02d}.npy",
              pins.cells_parquet(donor), pins.expression_parquet(donor), pins.geometry_parquet(donor),
              pins.qptiff_path(donor), pins.geojson_path(donor), pins.QUPATH_MANIFEST, cache]
    return RedseaWindow(donor=str(donor), status=pins.donor_status(donor), region_idx=idx, region=region, window=win,
                        core_px=core_px, px_um=float(px_um or pins.PIXEL_SIZE_UM), markers=tuple(markers), planes=planes,
                        ranges=ranges, mask=mask, lut_raw=lut_raw, lut_corrected=lut_corr, fill_ranges=fill_ranges,
                        status_by_label=status, counts=counts, inputs=inputs)
