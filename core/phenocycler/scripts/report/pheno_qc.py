"""PhenoCycler QC figures P01-P05 for one donor.  Inputs are the immutable run and the Step-2 aligned
values only -- no identity release is needed, nothing is called.  Every axis starts at the measured
floor 0 (asinh when the range demands it, stated on the axis), every distribution panel prints its
denominator, a missing measurement is drawn as Unavailable (grey) and never as 0.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import pheno_io  # noqa: E402
import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord  # noqa: E402

SEC = "pheno"


def _out(donor: str) -> pathlib.Path:
    return pins.OUT / str(donor) / SEC


def _rec(fid: str, donor: str, files: list, inputs: list, facts: dict) -> FigureRecord:
    return FigureRecord(figure_id=fid, donor=str(donor), section=SEC, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _sidecar(out_dir: pathlib.Path, fid: str, facts: dict) -> pathlib.Path:
    p = out_dir / f"{fid}.json"
    p.write_text(json.dumps(facts, indent=1, default=str) + "\n")
    return p


# --------------------------------------------------------------------------- P01
def p01_counts(donor: str, prov: dict, T: pheno_io.PhenoTables | None = None) -> FigureRecord:
    T = T or pheno_io.load_tables(donor)
    t = T.table
    n = len(t)
    counts = {"objects in run": n,
              "analysis eligible": int(t["analysis_eligible"].fillna(False).astype(bool).sum()),
              "estimation eligible": int(t["estimation_eligible"].fillna(False).astype(bool).sum()),
              "spillover-context eligible": int(t["spillover_context_eligible"].fillna(False).astype(bool).sum())}
    region = t["cell_region"].astype(str).value_counts().to_dict()
    fig, axes = style.figure(1, 2, width=5.6)
    ax = axes[0, 0]
    ax.barh(list(counts)[::-1], [counts[k] for k in list(counts)[::-1]], color=[style.ROLE_COLORS.get(k.split()[0], style.ELIGIBLE) for k in list(counts)[::-1]])
    for i, k in enumerate(list(counts)[::-1]):
        ax.text(counts[k], i, f" {counts[k]:,}", va="center", fontsize=8)
    ax.set_xlabel("cells"); ax.set_title(f"P01-A  run objects and eligibility (donor {donor}, {T.status})")
    style.denominator(ax, f"denominator: {n:,} run objects; GeoJSON features {T.n_geojson_features:,}; ingest rejects {T.n_ingest_rejects:,}", loc="lower right")
    ax = axes[0, 1]
    order = [r for r in ("tissue", "peri", "core") if r in region] + [r for r in region if r not in ("tissue", "peri", "core")]
    vals = [region[r] for r in order]
    ax.bar(order, vals, color=[{"core": "#CCBB44", "peri": "#EE8866", "tissue": "#77AADD"}.get(r, "#BBBBBB") for r in order])
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}\n{v / n:.1%}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("cells"); ax.set_title("P01-B  QuPath islet proximity (audit axis: Islet_N / exp20um)")
    style.denominator(ax, f"denominator: {n:,} run objects", loc="upper right")
    style.title(fig, f"P01 cell counts and eligibility -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "P01_cell_counts"); plt.close(fig)
    facts = {"counts": counts, "cell_region": region, "n_geojson_features": T.n_geojson_features, "n_ingest_rejects": T.n_ingest_rejects}
    files.append(_sidecar(out, "P01_cell_counts", facts))
    return _rec("P01", donor, files, T.inputs, facts)


# --------------------------------------------------------------------------- P02
def p02_areas(donor: str, prov: dict, T: pheno_io.PhenoTables | None = None) -> FigureRecord:
    T = T or pheno_io.load_tables(donor)
    t, thr = T.table, T.thresholds
    elig = t["analysis_eligible"].fillna(False).astype(bool).to_numpy()
    area = t["cell_area"].to_numpy(float)
    nuc, n_nuc_unavail = pheno_io.nucleus_area_um2(t)
    ceil = pheno_io.recompute_ceiling(t, thr)
    sol = t["cell_solidity"].to_numpy(float)
    fig, axes = style.figure(1, 3, width=5.4)
    ax = axes[0, 0]
    hi = float(np.nanquantile(area, 0.999)) * 1.5
    style.hist_from_zero(ax, area[elig], hi=hi, color=style.ELIGIBLE, label=f"analysis eligible n={elig.sum():,}")
    style.hist_from_zero(ax, area[~elig], hi=hi, color=style.EXCLUDED, label=f"excluded n={(~elig).sum():,}", alpha=0.6)
    style.threshold_line(ax, float(thr["min_cell_area_um2"]), f"min {thr['min_cell_area_um2']} um2")
    style.threshold_line(ax, ceil["ceiling_um2"], f"ceiling {ceil['multiple']:.0f} x median = {ceil['ceiling_um2']:.0f} um2")
    ax.set_xlabel(style.axis_label_from_zero("cell area", "um2", False)); ax.set_ylabel("cells"); ax.legend(fontsize=7, loc="upper right")
    ax.set_title("P02-A  cell area (cells stage, QuPath)")
    style.denominator(ax, f"denominator: {len(t):,} run objects; values > {hi:.0f} folded into the last bin", loc="lower right")
    ax = axes[0, 1]
    fin = np.isfinite(nuc)
    hi2 = float(np.nanquantile(nuc, 0.999)) * 1.5
    style.hist_from_zero(ax, nuc[fin & elig], hi=hi2, color=style.ELIGIBLE, label=f"eligible, nucleus measured n={(fin & elig).sum():,}")
    style.hist_from_zero(ax, nuc[fin & ~elig], hi=hi2, color=style.EXCLUDED, label=f"excluded n={(fin & ~elig).sum():,}", alpha=0.6)
    ax.set_xlabel(style.axis_label_from_zero("nucleus area", "um2", False)); ax.set_ylabel("cells")
    ax.set_title("P02-B  nucleus area = nucleus/cell ratio x raster ratio x cell area")
    ax.legend(handles=ax.get_legend_handles_labels()[0] + [style.unavailable_patch(f"no nucleus geometry n={n_nuc_unavail:,}")], fontsize=7, loc="upper right")
    style.denominator(ax, f"denominator: {len(t):,} run objects; {n_nuc_unavail:,} Unavailable (not drawn as 0)", loc="lower right")
    ax = axes[0, 2]
    ax.hist(sol[elig & np.isfinite(sol)], bins=100, range=(0, 1), color=style.ELIGIBLE, alpha=0.85, label="eligible")
    ax.hist(sol[~elig & np.isfinite(sol)], bins=100, range=(0, 1), color=style.EXCLUDED, alpha=0.6, label="excluded")
    style.threshold_line(ax, float(thr["min_cell_solidity"]), f"min solidity {thr['min_cell_solidity']}")
    ax.set_xlabel("cell solidity (0-1)"); ax.set_ylabel("cells"); ax.legend(fontsize=7); ax.set_title("P02-C  solidity")
    style.denominator(ax, f"denominator: {len(t):,} run objects", loc="upper left")
    style.title(fig, f"P02 area and shape QC with the run's geometry thresholds -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "P02_area_distributions"); plt.close(fig)
    facts = {"thresholds": thr, "ceiling": ceil, "n_nucleus_unavailable": n_nuc_unavail,
             "cell_area_median_eligible": float(np.nanmedian(area[elig])), "nucleus_area_median_eligible": float(np.nanmedian(nuc[fin & elig]))}
    files.append(_sidecar(out, "P02_area_distributions", facts))
    return _rec("P02", donor, files, T.inputs, facts)


# --------------------------------------------------------------------------- P03
def p03_eligibility(donor: str, prov: dict, T: pheno_io.PhenoTables | None = None) -> FigureRecord:
    T = T or pheno_io.load_tables(donor)
    t = T.table
    n = len(t)
    fig, axes = style.figure(1, 2, width=6.4)
    first = t["analysis_reason"].fillna("").astype(str).replace("", "(eligible)").value_counts()
    ax = axes[0, 0]
    ax.barh(first.index[::-1], first.to_numpy()[::-1], color=[style.ELIGIBLE if k == "(eligible)" else style.EXCLUDED for k in first.index[::-1]])
    for i, (k, v) in enumerate(zip(first.index[::-1], first.to_numpy()[::-1])):
        ax.text(v, i, f" {v:,} ({v / n:.1%})", va="center", fontsize=8)
    ax.set_xlabel("cells"); ax.set_title("P03-A  first analysis reason (geometry_qc.analysis_reason)")
    style.denominator(ax, f"denominator: {n:,} run objects", loc="lower right")
    combos = t["analysis_reasons"].fillna("").astype(str)
    combos = combos[combos != ""].value_counts().head(12)
    ax = axes[0, 1]
    ax.barh(combos.index[::-1], combos.to_numpy()[::-1], color=style.EXCLUDED)
    for i, v in enumerate(combos.to_numpy()[::-1]):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    ax.set_xlabel("cells"); ax.set_title("P03-B  top reason combinations (';'-joined; overlapping)")
    style.denominator(ax, f"denominator: {int((t['analysis_reasons'].fillna('') != '').sum()):,} excluded objects", loc="lower right")
    style.title(fig, f"P03 geometry-QC eligibility by reason -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "P03_eligibility_by_reason"); plt.close(fig)
    facts = {"first_reason": first.to_dict(), "top_combinations": combos.to_dict()}
    files.append(_sidecar(out, "P03_eligibility_by_reason", facts))
    return _rec("P03", donor, files, T.inputs, facts)


# --------------------------------------------------------------------------- P04
def p04_marker_intensities(donor: str, prov: dict) -> FigureRecord:
    marks, inputs = pheno_io.load_aligned(donor)
    ncol = 6
    nrow = int(np.ceil(len(marks) / ncol))
    fig, axes = style.figure(nrow, ncol, width=3.4, height=2.6)
    facts = {}
    for k, m in enumerate(marks):
        ax = axes[k // ncol, k % ncol]
        if not m.available or m.n_finite == 0:
            style.note(ax, f"{m.marker}\nUnavailable: {m.unavailable_reason or 'no finite values'}")
            ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("#EEEEEE")
            facts[m.marker] = {"available": False, "reason": m.unavailable_reason}
            continue
        v = np.asarray(m.values, float)
        v = v[np.isfinite(v)]
        pos = v[v > 0]
        hi = float(np.quantile(pos, 0.999)) if pos.size else 1.0
        use_asinh = True
        x = np.arcsinh(np.clip(v, 0, None))
        ax.hist(x, bins=120, range=(0, np.arcsinh(hi) * 1.05), color=style.ELIGIBLE if not m.nucleus_measured else "#8c6bb1", alpha=0.9)
        if np.isfinite(m.mode_aligned):
            ax.axvline(np.arcsinh(max(m.mode_aligned, 0)), color=style.MODE, lw=1.0, ls="--")
        ax.set_yscale("log")
        ttl = m.marker + (" [nucleus-measured, raw]" if m.nucleus_measured else "") + (f" [{m.registry_kind}]" if m.registry_kind else "")
        ax.set_title(ttl, fontsize=8)
        ax.text(0.98, 0.95, f"n={m.n_finite:,}\nzero {m.n_zero / max(m.n_finite, 1):.0%}\nmode {m.mode_aligned:.3g}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
        ax.tick_params(labelsize=6)
        facts[m.marker] = {"available": True, "n_finite": m.n_finite, "n_zero": m.n_zero, "mode_aligned": m.mode_aligned,
                           "scale": m.scale, "nucleus_measured": m.nucleus_measured, "kind": m.registry_kind}
    for k in range(len(marks), nrow * ncol):
        axes[k // ncol, k % ncol].axis("off")
    fig.text(0.5, 0.955, "x: asinh(aligned intensity) from the measured floor 0; y: log count; red dashed = Step-2 background mode; "
             "nucleus-measured markers (purple) are raw QuPath nucleus means, not REDSEA-corrected", ha="center", fontsize=8)
    style.title(fig, f"P04 per-marker intensity distributions on the aligned scale -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig, top=0.93)
    out = _out(donor); files = style.savefig(fig, out, "P04_marker_intensities"); plt.close(fig)
    files.append(_sidecar(out, "P04_marker_intensities", facts))
    return _rec("P04", donor, files, inputs, {"n_markers": len(marks)})


# --------------------------------------------------------------------------- P05
def p05_redsea_window(donor: str, prov: dict) -> FigureRecord:
    W = pheno_io.load_redsea_window(donor)
    markers = list(W.markers)
    fig, axes = style.figure(len(markers), 3, width=4.6, height=4.4)
    lut_status = W.status_by_label
    for i, m in enumerate(markers):
        ax = axes[i, 0]
        rng_ = W.ranges[m]
        img = W.planes[m].astype(float)
        ax.imshow(np.clip(img, rng_.vmin, rng_.vmax), cmap="gray", vmin=rng_.vmin, vmax=rng_.vmax, interpolation="nearest")
        ax.set_title(f"{m} raw plane  [window {rng_.vmin:.0f}-{rng_.vmax:.0f}; {getattr(rng_, 'frac_at_zero', float('nan')):.1%} at 0]", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        for j, (kind, lut) in enumerate((("raw QuPath mean", W.lut_raw[m]), ("REDSEA-corrected mean", W.lut_corrected[m]))):
            ax = axes[i, 1 + j]
            fr = W.fill_ranges[m]
            filled = np.full(W.mask.shape, np.nan)
            lab = W.mask
            ok = lab > 0
            vals = lut[lab[ok]]
            filled[ok] = vals
            cmap = matplotlib.colormaps["viridis"].copy(); cmap.set_bad(style.UNAVAILABLE)
            ax.imshow(np.clip(filled, fr.vmin, fr.vmax), cmap=cmap, vmin=fr.vmin, vmax=fr.vmax, interpolation="nearest")
            # spillover-ineligible cells: grey overlay
            st = lut_status[lab] if lut_status is not None else None
            if st is not None:
                grey = np.where(ok & np.isin(st, [2, 3]), 1.0, np.nan)
                ax.imshow(grey, cmap=matplotlib.colors.ListedColormap([style.UNAVAILABLE]), interpolation="nearest", alpha=0.9)
            ax.set_title(f"{m} cells filled by {kind}  [0-{fr.vmax:.3g}]", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
    fig.legend(handles=[style.unavailable_patch("spillover-ineligible or not in run (grey)")], loc="lower right", fontsize=8)
    style.title(fig, f"P05 REDSEA window: raw planes vs per-cell raw / corrected means -- donor {donor}, frame region {W.region_idx} ({W.region.get('draw', '')}) at {W.px_um:.4f} um/px")
    style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "P05_redsea_window"); plt.close(fig)
    facts = {"region_idx": W.region_idx, "window": W.window, "markers": markers, "counts": W.counts}
    files.append(_sidecar(out, "P05_redsea_window", facts))
    return _rec("P05", donor, files, W.inputs, facts)


FIGURES = {"P01": p01_counts, "P02": p02_areas, "P03": p03_eligibility, "P04": p04_marker_intensities, "P05": p05_redsea_window}


def run(donor: str, prov: dict, only: list | None = None) -> list:
    recs = []
    T = None
    for fid, fn in FIGURES.items():
        if only and fid not in only:
            continue
        if fid in ("P01", "P02", "P03"):
            T = T or pheno_io.load_tables(donor)
            recs.append(fn(donor, prov, T))
        else:
            recs.append(fn(donor, prov))
        print(f"[report] {fid} donor {donor}: {recs[-1].files[0]}", flush=True)
    return recs
