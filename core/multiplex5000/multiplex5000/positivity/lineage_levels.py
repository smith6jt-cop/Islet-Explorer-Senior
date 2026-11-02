"""lineage_levels.py -- per-section marker levels from mutually exclusive lineages (spleen PC67, 2026-11-01).

The spleen registry is not a copy of another tissue's (maintainer 2026-11-01: "Start the registry over"; "Stop trying
to apply parameters used from a different dataset"). In lymphoid tissue nearly every marker has a population that
cannot carry it, so a marker's level is read on that population, per section, as the source method's known-negative
rule (Barlow et al. eLife 2024: positive above the 99th percentile of cells known not to express the marker) -- with
spleen populations and spleen-measured constants:

- **Lineages from the data.** The B and T lineages' required markers (CD20, CD3e) have a population above their
  background on every section. Pass 1 splits each where it begins (``lineage_split``: the density minimum above the
  background mode, Gaussian KDE with Scott's bandwidth, else a two-component mixture boundary); pass 2 replaces
  those splits with the two markers' own known-negative levels (CD20 on T cells, CD3e on B cells). B cells: the B
  marker positive and the T marker negative; T cells: the converse.
- **Known-negative levels** (registry method ``known_negative``): the LEVEL_QUANTILE-th percentile of the marker in its
  reference populations (registry ``level_reference``: ``b_cells``, ``t_cells``, ``lymphocytes``, ``non_haematopoietic``).
  A reference cell with a cell above a first-pass level within ISOLATION_UM is left out, since a neighbour's signal is
  not the marker's background; the cell itself stays, since its own tail is what the percentile measures. A reference
  with fewer than KNOWN_NEGATIVE_MIN_ABOVE cells above the percentile gives no level.
- **Non-haematopoietic cells** (CD45's reference): above their levels for a marker of the Mesenchymal, Endothelial or
  Epithelial compartment, in neither lymphoid lineage, with no lymphocyte within ISOLATION_UM.
- **Markers no lineage excludes** (registry method ``valley``): the density minimum between the two highest modes;
  ``background``: the background mode plus BACKGROUND_K half-widths (K measured on the spleen known-negative markers).
- **Readings:** each marker's positivity measurement (whole cell; nuclear antigens on the nucleus), and for a marker
  with a differentiation measurement the same method on that reading too (owner 2026-10-30).

``choose`` returns the tables ``levels.choose`` returns (levels, levels_differentiation, selection, chosen, positivity,
options, differentiation_selection, stops, immune_set_aside), so canary, coverage and the hierarchy read them unchanged.
A method that gives no level on a section is reported against the method, for every section ("flag the methods you
are using so they will work with the images and markers", owner 2026-10-28).
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde

from . import registry as RG

POPULATIONS = ("b_cells", "t_cells", "lymphocytes", "non_haematopoietic")
METHODS = ("known_negative", "valley", "background")
NON_HAEM_COMPARTMENTS = ("Mesenchymal", "Endothelial", "Epithelial")
KDE_SAMPLE = 200_000
GRID = 600
SEED = 0


class LineageLevelError(ValueError):
    """The registry or the cells cannot give a level as declared."""


def log1p(x) -> np.ndarray:
    return np.log1p(np.clip(np.nan_to_num(np.asarray(x, float), nan=0.0), 0, None))


def _kde(v: np.ndarray):
    v = v[v > 0]
    if len(v) < 1000:
        return None, None
    rng = np.random.default_rng(SEED)
    s = v if len(v) <= KDE_SAMPLE else rng.choice(v, KDE_SAMPLE, replace=False)
    grid = np.linspace(np.percentile(s, 0.1), np.percentile(s, 99.9), GRID)
    return grid, gaussian_kde(s)(grid)


def _modes(dens: np.ndarray) -> list[int]:
    return [i for i in range(1, len(dens) - 1) if dens[i] >= dens[i - 1] and dens[i] >= dens[i + 1]]


def valley(v: np.ndarray) -> tuple[float, float, float]:
    """(lower mode, upper mode, valley): the density minimum between the two highest modes; NaN if unimodal."""
    grid, dens = _kde(v)
    if grid is None:
        return np.nan, np.nan, np.nan
    peaks = _modes(dens)
    if len(peaks) < 2:
        return (float(grid[peaks[0]]) if peaks else np.nan), np.nan, np.nan
    a, b = sorted(sorted(peaks, key=lambda i: dens[i])[-2:])
    j = a + int(np.argmin(dens[a:b + 1]))
    return float(grid[a]), float(grid[b]), float(grid[j])


def lineage_split(v: np.ndarray) -> float:
    """Where a lineage marker's positives (brighter than its background) begin. The density minimum between the
    background (highest) mode and the highest mode above it; where the positives form only a shoulder, the first point
    above the mode where the density reaches twice the background mirrored about its mode (the left side, reflected),
    i.e. where a cell is as likely positive as background. Pass 1 only: pass 2 replaces it with the known-negative
    level."""
    grid, dens = _kde(v)
    if grid is None:
        return np.nan
    top = int(np.argmax(dens))
    right = [i for i in _modes(dens) if i > top]
    if right:
        b = max(right, key=lambda i: dens[i])
        return float(grid[top + int(np.argmin(dens[top:b + 1]))])
    x = grid[top + 1:]
    mirrored = np.interp(2 * grid[top] - x, grid, dens, left=0.0, right=0.0)
    over = np.flatnonzero(dens[top + 1:] >= 2 * mirrored)
    return float(x[over[0]]) if len(over) else np.nan


def background(v: np.ndarray, k: float) -> tuple[float, float, float]:
    """(mode, left half-width, mode + k half-widths) of the highest KDE mode."""
    grid, dens = _kde(v)
    if grid is None:
        return np.nan, np.nan, np.nan
    i = int(np.argmax(dens))
    j = i
    while j > 0 and dens[j] > dens[i] / 2:
        j -= 1
    if dens[j] > dens[i] / 2:
        return float(grid[i]), np.nan, np.nan
    w = float(grid[i] - grid[j]) / np.sqrt(2 * np.log(2))       # half-width at half maximum -> sigma
    return float(grid[i]), w, float(grid[i] + k * w)


def away(xy: np.ndarray, ref: np.ndarray, src: np.ndarray, cut: float) -> np.ndarray:
    """Cells of `ref` with no `src` cell other than themselves within `cut` um."""
    if cut <= 0 or not src.any() or not ref.any():
        return ref.copy()
    idx = np.flatnonzero(ref)
    near = cKDTree(xy[src]).query_ball_point(xy[idx], r=cut, return_length=True)
    out = np.zeros(len(ref), bool)
    out[idx[(np.asarray(near) - src[idx].astype(int)) == 0]] = True
    return out


def known_negative(v: np.ndarray, ref: np.ndarray, xy: np.ndarray, q: float, cut: float, min_above: float) -> tuple[float, int]:
    """(level, reference cells used): two passes, the second without reference cells next to first-pass positives."""
    def pct(mask):
        o = v[mask]
        return float(np.percentile(o, q)) if len(o) * (1 - q / 100) >= min_above else np.nan
    first = pct(ref)
    if not np.isfinite(first):
        return np.nan, int(ref.sum())
    iso = away(xy, ref, v > first, cut)
    return pct(iso), int(iso.sum())


def _lineage_marker(reg: RG.Registry, lineage: str) -> str:
    ms = [m for m in reg.broad(lineage) if reg.table.loc[m, "required"] == "yes" and m in reg.measured_markers()]
    if len(ms) != 1:
        raise LineageLevelError(f"{reg.name}: lineage {lineage} needs exactly one required marker, found {ms}")
    return ms[0]


def _reading(d: pd.DataFrame, reg: RG.Registry, m: str, meas: str) -> np.ndarray:
    return log1p(d[reg.column(m, meas)])


def image_levels(source, tag: str, reg_name: tuple) -> list[dict]:
    reg = RG.load(*reg_name)
    q, cut = reg.const("LEVEL_QUANTILE"), reg.const("ISOLATION_UM")
    k, min_above = reg.const("BACKGROUND_K"), reg.const("KNOWN_NEGATIVE_MIN_ABOVE")
    markers = reg.measured_markers()
    cols = sorted({reg.column(m, reg.measurement(m)) for m in markers}
                  | {reg.column(m, reg.differentiation(m)) for m in markers if reg.differentiation(m)})
    d = source.levels_frame(tag, cols, [])
    xy = d[["x_um", "y_um"]].to_numpy(float)
    n = len(d)
    bm, tm = _lineage_marker(reg, "B"), _lineage_marker(reg, "T")
    vb, vt = _reading(d, reg, bm, reg.measurement(bm)), _reading(d, reg, tm, reg.measurement(tm))
    cb, ct = lineage_split(vb), lineage_split(vt)
    if not (np.isfinite(cb) and np.isfinite(ct)):
        raise LineageLevelError(f"{tag}: {bm} or {tm} has no population above its background; the lineages cannot be "
                                "read from the data")
    rows: list[dict] = []
    for npass in (1, 2):
        B, T = (vb > cb) & ~(vt > ct), (vt > ct) & ~(vb > cb)
        pops = {"b_cells": B, "t_cells": T, "lymphocytes": B | T}
        levels: dict[tuple[str, str], float] = {}
        rows = []

        def one(m, meas, reading):
            method = reg.background_method(m)
            v = _reading(d, reg, m, meas)
            ref_n, center, width = 0, np.nan, np.nan
            if method == "known_negative":
                names = RG.split(reg.table.loc[m, "level_reference"].replace("+", ";"))
                if not names or set(names) - set(POPULATIONS):
                    raise LineageLevelError(f"{m}: level_reference {names} is not among {POPULATIONS}")
                if "non_haematopoietic" in names and "non_haematopoietic" not in pops:
                    return None                                     # read after the non-haematopoietic markers
                ref = np.logical_or.reduce([pops[x] for x in names])
                lvl, ref_n = known_negative(v, ref, xy, q, cut, min_above)
            elif method == "valley":
                center, _, lvl = valley(v)
            elif method == "background":
                center, width, lvl = background(v, k)
            else:
                raise LineageLevelError(f"{m}: method {method!r} is not among {METHODS}")
            levels[(m, meas)] = lvl
            pos = v > lvl if np.isfinite(lvl) else np.zeros(n, bool)
            rows.append({"image": tag, "marker": m, "measurement": meas, "method": method, "reading": reading,
                         "pass": npass, "center": center, "width": width, "level_log": lvl,
                         "level_raw": float(np.expm1(lvl)) if np.isfinite(lvl) else np.nan,
                         "positive_fraction": float(pos.mean()), "zero_fraction": float(np.mean(v <= 0)),
                         "reference_cells": ref_n, "background_center": center, "background_width": width})
            return lvl

        deferred = []
        for m in markers:
            for meas, reading in [(reg.measurement(m), "positivity")] + \
                    ([(reg.differentiation(m), "differentiation")] if reg.differentiation(m) else []):
                if one(m, meas, reading) is None:
                    deferred.append((m, meas, reading))
        nonhaem = np.zeros(n, bool)
        for m in markers:
            if reg.role(m) == "defining" and reg.table.loc[m, "compartment"] in NON_HAEM_COMPARTMENTS:
                lvl = levels.get((m, reg.measurement(m)), np.nan)
                if np.isfinite(lvl):
                    nonhaem |= _reading(d, reg, m, reg.measurement(m)) > lvl
        pops["non_haematopoietic"] = away(xy, nonhaem & ~(B | T), B | T, cut)
        for m, meas, reading in deferred:
            one(m, meas, reading)
        if npass == 1:                      # pass 2: the lineages from the two markers' own known-negative levels
            nb, nt = levels.get((bm, reg.measurement(bm)), np.nan), levels.get((tm, reg.measurement(tm)), np.nan)
            if not (np.isfinite(nb) and np.isfinite(nt)):
                break
            cb, ct = nb, nt
    for r in rows:
        r["b_cells"], r["t_cells"] = int(((vb > cb) & ~(vt > ct)).sum()), int(((vt > ct) & ~(vb > cb)).sum())
    return rows


def choose(source, reg_name: tuple, tags: list[str], workers: int = 1) -> dict:
    """Every section's levels and the per-marker decisions, in ``levels.choose``'s table shapes."""
    with ProcessPoolExecutor(max(1, min(workers, len(tags)))) as ex:
        rows = [r for rs in ex.map(image_levels, repeat(source), tags, repeat(reg_name)) for r in rs]
    t = pd.DataFrame(rows)
    cols = ["image", "marker", "measurement", "method", "center", "width", "level_log", "level_raw",
            "positive_fraction", "zero_fraction", "background_center", "background_width"]
    lv = t[t["reading"] == "positivity"]
    ld = t[t["reading"] == "differentiation"]
    decisions, diff = [], []
    for (m, reading), g in t.groupby(["marker", "reading"], sort=False):
        failing = sorted(g.loc[~np.isfinite(g["level_log"]), "image"])
        rec = {"marker": m, "decision": "stop" if failing else "use", "measurement": g["measurement"].iloc[0],
               "method": g["method"].iloc[0], "out_of_place": np.nan,
               "sensitivity": float(g["positive_fraction"].median()),
               "why": f"{g['method'].iloc[0]} gives no level on: {' '.join(failing)}" if failing else ""}
        (decisions if reading == "positivity" else diff).append(rec)
    sel, dsel = pd.DataFrame(decisions), pd.DataFrame(diff)
    reg = RG.load(*reg_name)
    # levels.choose's shapes: chosen = each marker on the measurement it is compared on (its differentiation reading
    # where it has one); positivity = the whole-cell (positivity) level of the markers compared on another reading
    chosen = {r.marker: (r.measurement, r.method) for r in sel.itertuples()}
    diff_markers = set(dsel["marker"]) if len(dsel) else set()
    positivity = {m: chosen[m] for m in diff_markers}
    chosen.update({r.marker: (r.measurement, r.method) for r in dsel.itertuples()})
    stops = list(sel.loc[sel["decision"] == "stop", "marker"])
    return {"options": t, "selection": sel, "differentiation_selection": dsel, "levels": lv[cols].reset_index(drop=True),
            "levels_differentiation": ld[cols].reset_index(drop=True), "chosen": chosen, "positivity": positivity,
            "immune_set_aside": [], "stops": stops}


def canary_levels(tables: dict, reg_name: tuple) -> pd.DataFrame:
    """canary.levels' table for a lineage-levels run. The canary null of a marker is the cells of compartments that
    cannot carry it; in spleen that is the same lymphoid population its level is read on, at the same quantile, so the
    canary level is the level and canary depth counts a cell's positive typing markers."""
    reg = RG.load(*reg_name)
    lv = tables["levels"]
    typing = set(reg.typing_markers())
    out = lv[lv["marker"].isin(typing)][["image", "marker", "level_log"]].copy()
    out["canary_log"] = out["level_log"]
    out["canary_fold_over_level"] = 1.0
    out["null_reference"] = [reg.table.loc[m, "level_reference"] or reg.background_method(m) for m in out["marker"]]
    return out.reset_index(drop=True)
