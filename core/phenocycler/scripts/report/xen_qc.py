"""Xenium QC figures X01-X03: every bundle cell, the base builder's QC recomputed, the negative controls."""
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

import pins  # noqa: E402
import style  # noqa: E402
import xen_io  # noqa: E402
from manifest import FigureRecord  # noqa: E402

SEC = "xen"


def _rec(fid, donor, files, inputs, facts, section=SEC):
    return FigureRecord(figure_id=fid, donor=str(donor), section=section, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _sidecar(out_dir: pathlib.Path, fid: str, facts: dict) -> pathlib.Path:
    p = out_dir / f"{fid}.json"
    p.write_text(json.dumps(facts, indent=1, default=str) + "\n")
    return p


def x01_cell_qc(donor: str, prov: dict, X: xen_io.XenTables | None = None) -> FigureRecord:
    X = X or xen_io.load_tables(donor)
    c = X.cells
    keep = c["retained_base"].to_numpy(bool)
    qc = X.qc_recomputed
    panels = [("n_counts_matrix", "transcripts per cell (gene features)", xen_io.MIN_COUNTS, "min_counts"),
              ("n_genes", "genes detected per cell", xen_io.MIN_GENES, "min_genes"),
              ("cell_area", "cell area (um2)", None, None),
              ("density", "transcripts per um2", qc.get("density_cap_counts_per_um2"), "density cap (q98)")]
    fig, axes = style.figure(1, 5, width=4.2)
    facts = {}
    for k, (col, label, thr, thr_name) in enumerate(panels):
        ax = axes[0, k]
        v = c[col].to_numpy(float)
        hi = float(np.nanquantile(v[np.isfinite(v)], 0.999)) * 1.3 if np.isfinite(v).any() else 1.0
        asinh = style.range_demands_asinh(v)
        style.hist_from_zero(ax, v[keep], hi=hi, asinh=asinh, color=style.ELIGIBLE, label=f"retained n={keep.sum():,}")
        style.hist_from_zero(ax, v[~keep], hi=hi, asinh=asinh, color=style.EXCLUDED, label=f"dropped n={(~keep).sum():,}", alpha=0.6)
        if thr is not None:
            style.threshold_line(ax, np.arcsinh(thr) if asinh else thr, f"{thr_name} = {thr:.3g}")
        ax.set_xlabel(style.axis_label_from_zero(label, "", asinh)); ax.set_ylabel("cells"); ax.legend(fontsize=7)
        ax.set_title(f"X01-{'ABCD'[k]}  {label}", fontsize=9)
        style.denominator(ax, f"denominator: {len(c):,} bundle cells", loc="lower right")
        facts[col] = {"median_retained": float(np.nanmedian(v[keep])), "median_dropped": float(np.nanmedian(v[~keep])) if (~keep).any() else None, "threshold": thr}
    ax = axes[0, 4]
    ratio = c["nucleus_area"].to_numpy(float) / np.where(c["cell_area"].to_numpy(float) > 0, c["cell_area"].to_numpy(float), np.nan)
    ax.hist(ratio[keep & np.isfinite(ratio)], bins=100, range=(0, 1), color=style.ELIGIBLE, alpha=0.85, label="retained")
    ax.hist(ratio[~keep & np.isfinite(ratio)], bins=100, range=(0, 1), color=style.EXCLUDED, alpha=0.6, label="dropped")
    ax.set_xlabel("nucleus area / cell area"); ax.set_ylabel("cells"); ax.legend(fontsize=7); ax.set_title("X01-E  nucleus / cell area", fontsize=9)
    style.denominator(ax, f"denominator: {len(c):,} bundle cells; {int((~np.isfinite(ratio)).sum()):,} Unavailable", loc="upper left")
    style.title(fig, f"X01 Xenium per-cell QC -- donor {donor} (serial {X.serial}, {X.status}); retained {keep.sum():,} of {len(c):,}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X01_cell_qc_distributions"); plt.close(fig)
    facts["qc_recomputed"] = qc
    files.append(_sidecar(out, "X01_cell_qc_distributions", facts))
    return _rec("X01", donor, files, X.inputs, facts)


def x02_drop_reasons(donor: str, prov: dict, X: xen_io.XenTables | None = None) -> FigureRecord:
    X = X or xen_io.load_tables(donor)
    c = X.cells
    f = c[["fail_min_genes", "fail_min_counts", "fail_density"]].astype(bool)
    dropped = ~c["retained_base"].to_numpy(bool)
    combos = f.apply(lambda r: "+".join([k.replace("fail_", "") for k, v in r.items() if v]) or "(none)", axis=1)
    vc = combos[dropped].value_counts()
    singles = f.sum()
    fig, axes = style.figure(1, 2, width=6.0)
    ax = axes[0, 0]
    ax.barh(vc.index[::-1], vc.to_numpy()[::-1], color=style.EXCLUDED)
    for i, v in enumerate(vc.to_numpy()[::-1]):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    ax.set_xlabel("dropped cells"); ax.set_title("X02-A  overlapping QC failures among dropped cells (exact set per cell)")
    style.denominator(ax, f"denominator: {int(dropped.sum()):,} dropped of {len(c):,} bundle cells", loc="lower right")
    ax = axes[0, 1]
    ax.bar(list(singles.index.str.replace("fail_", "")), singles.to_numpy(), color=style.EXCLUDED)
    for i, v in enumerate(singles.to_numpy()):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("cells failing the rule"); ax.set_title(f"X02-B  per-rule counts (sum {int(singles.sum()):,} > actual drop {int(dropped.sum()):,}: the masks overlap)")
    style.denominator(ax, f"denominator: {len(c):,} bundle cells", loc="upper right")
    style.title(fig, f"X02 QC-drop reasons -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X02_qc_drop_reasons"); plt.close(fig)
    facts = {"combinations": vc.to_dict(), "per_rule": singles.to_dict(), "n_dropped": int(dropped.sum()), "n_cells": int(len(c)),
             "retained_sets_identical": X.qc_recomputed.get("retained_sets_identical")}
    files.append(_sidecar(out, "X02_qc_drop_reasons", facts))
    return _rec("X02", donor, files, X.inputs, facts)


def x03_negative_controls(donors: list, prov: dict, tables: dict | None = None) -> FigureRecord:
    tables = tables or {d: xen_io.load_tables(d) for d in donors}
    types = list(pins.XENIUM_FEATURE_TYPES) if hasattr(pins, "XENIUM_FEATURE_TYPES") else sorted(next(iter(tables.values())).feature_types)
    fig, axes = style.figure(1, 2, width=6.2)
    ax = axes[0, 0]
    w = 0.8 / len(donors)
    facts = {}
    for i, d in enumerate(donors):
        X = tables[d]
        vals = [X.per_type_mean_per_cell.get(t, np.nan) for t in types]
        ax.bar(np.arange(len(types)) + i * w, vals, width=w, label=f"{d} ({X.status}, {pins.XENIUM_SERIAL.get(str(d), '')})", color=style.status_color(X.status))
        facts[d] = {"per_type_mean_per_cell": X.per_type_mean_per_cell, "feature_types": X.feature_types, "per_type_total": X.per_type_total}
    ax.set_xticks(np.arange(len(types)) + w * (len(donors) - 1) / 2)
    ax.set_xticklabels([f"{t}\n(n={tables[donors[0]].feature_types.get(t, '?')})" for t in types], fontsize=7)
    ax.set_yscale("log"); ax.set_ylabel("mean decoded counts per cell per feature (log)"); ax.legend(fontsize=7)
    ax.set_title("X03-A  per-feature-type mean counts per cell (all bundle cells)")
    style.denominator(ax, "denominator: every bundle cell; per feature, so feature-set size does not inflate a type", loc="lower left")
    ax = axes[0, 1]
    rows = []
    for d in donors:
        ms = tables[d].metrics_summary
        if ms is None or len(ms) == 0:
            rows.append((d, "metrics_summary.csv not transferred"))
            continue
        cols = [c for c in ms.columns if any(k in c.lower() for k in ("negative", "unassigned", "q20", "median_transcripts", "decod", "cells_detected"))]
        rows.append((d, "\n".join(f"{c}: {ms.iloc[0][c]}" for c in cols[:10])))
    ax.axis("off")
    ax.text(0.0, 1.0, "\n\n".join(f"donor {d}:\n{txt}" for d, txt in rows), va="top", ha="left", fontsize=7, family="monospace", transform=ax.transAxes)
    ax.set_title("X03-B  instrument metrics_summary (decoding QC) where transferred")
    style.title(fig, "X03 negative controls and decoding QC -- donors " + ", ".join(donors))
    style.footer(fig, prov); style.finish(fig)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, "X03_negative_controls"); plt.close(fig)
    files.append(_sidecar(out, "X03_negative_controls", facts))
    inputs = [p for d in donors for p in tables[d].inputs]
    return _rec("X03", "_pair", files, inputs, facts, section="_pair")


FIGURES = {"X01": x01_cell_qc, "X02": x02_drop_reasons}


def run(donor: str, prov: dict, only: list | None = None) -> list:
    recs = []
    X = xen_io.load_tables(donor)
    for fid, fn in FIGURES.items():
        if only and fid not in only:
            continue
        recs.append(fn(donor, prov, X))
        print(f"[report] {fid} donor {donor}: {recs[-1].files[0]}", flush=True)
    return recs
