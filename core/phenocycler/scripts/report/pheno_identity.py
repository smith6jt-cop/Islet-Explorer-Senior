"""PhenoCycler identity figures P06-P11 from the development release (data/development/broadtype/release).

Nothing here re-types a cell: labels are read from the release, `evidence_positive` is read as written
(True only for a verified label), Ambiguous is hatched, Unavailable is grey, and every panel prints its
denominator.  A missing release renders a stated 'not available' panel.
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

import identity_io as io  # noqa: E402
import pheno_io  # noqa: E402
import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord  # noqa: E402

SEC = "pheno"
MOD = "phenocycler"


def _rec(fid, donor, files, inputs, facts, section=SEC):
    return FigureRecord(figure_id=fid, donor=str(donor), section=section, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _side(out, fid, facts):
    p = out / f"{fid}.json"; p.write_text(json.dumps(facts, indent=1, default=str) + "\n"); return p


def _unavailable(fid: str, donor: str, prov: dict, why: str, section=SEC) -> FigureRecord:
    fig, axes = style.figure(1, 1, width=7, height=3)
    axes[0, 0].axis("off"); style.note(axes[0, 0], f"{fid}: not available\n{why}")
    style.title(fig, f"{fid} -- donor {donor}: not available"); style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / (str(donor) if section != "_pair" else "_pair") / (SEC if section != "_pair" else "pheno")
    files = style.savefig(fig, out, f"{fid}_not_available"); plt.close(fig)
    return _rec(fid, donor, files, [], {"not_available": why}, section=section)


def _calls(donor: str) -> pd.DataFrame | None:
    return io.release_calls(MOD, donor)


def label_order(common: dict) -> list:
    return [L for L in common] + [io.AMBIGUOUS, io.UNAVAILABLE]


# ------------------------------------------------------------------ P06 marker x label heatmap
def p06_marker_by_label(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("P06", donor, prov, "no PhenoCycler release yet (data/development/broadtype/release)")
    tab = pd.read_parquet(pins.uniform_parquet(donor))
    tab["object_id"] = tab["object_id"].astype(str)
    dec = json.loads((pins.PKG / "scripts" / "broadtype" / "declarations.json").read_text())
    ident = [m for m in dec["identity_markers"] if m in tab.columns]
    state = [m for m in dec["state_markers"] if m in tab.columns]
    markers = ident + state
    df = tab[["object_id", *markers]].merge(calls[["object_id", "lineage_native"]], on="object_id", how="inner")
    labs = [L for L in io.common_of(MOD)] + [io.AMBIGUOUS]
    df = df[df["lineage_native"].isin(labs)]
    X = np.arcsinh(np.clip(df[markers].to_numpy(float), 0, None))
    med = np.nanmedian(X, axis=0); mad = np.nanmedian(np.abs(X - med), axis=0); mad = np.where(mad > 0, mad, 1.0)
    Z = (X - med) / mad
    lab = df["lineage_native"].to_numpy()
    counts = {L: int((lab == L).sum()) for L in labs}
    H = np.full((len(markers), len(labs)), np.nan)
    for j, L in enumerate(labs):
        sel = lab == L
        if sel.sum() >= 50:
            H[:, j] = np.nanmedian(Z[sel], axis=0)
    # label-permutation null on a hash subsample: the max |median z| over labels x markers
    rng = np.random.default_rng(0)
    sub = np.sort(rng.choice(len(lab), min(200_000, len(lab)), replace=False))
    null_max = []
    for _ in range(50):
        perm = lab[sub][rng.permutation(len(sub))]
        mx = 0.0
        for L in labs:
            sel = perm == L
            if sel.sum() >= 50:
                mx = max(mx, float(np.nanmax(np.abs(np.nanmedian(Z[sub][sel], axis=0)))))
        null_max.append(mx)
    null95 = float(np.quantile(null_max, 0.95))
    fig, axes = style.figure(1, 1, width=9, height=0.28 * len(markers) + 2)
    ax = axes[0, 0]
    vmax = float(np.nanmax(np.abs(H))) if np.isfinite(H).any() else 1.0
    im = ax.imshow(H, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(labs))); ax.set_xticklabels([f"{L}\nn={counts[L]:,}" for L in labs], fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels([m + ("  [state]" if m in state else "") for m in markers], fontsize=7)
    ax.axhline(len(ident) - 0.5, color="black", lw=1.0)
    ax.text(len(labs) - 0.5, len(ident) - 0.5, " identity | state ", va="center", ha="left", fontsize=7)
    masked = np.abs(H) < null95
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            if np.isfinite(H[i, j]) and masked[i, j]:
                ax.text(j, i, "·", ha="center", va="center", fontsize=6, color="#555555")
    fig.colorbar(im, ax=ax, fraction=0.02, label="median robust z of asinh(aligned) within label")
    ax.set_title(style.wrap(f"P06  marker x label medians; '·' = below the label-permutation null (max |z| p95 = {null95:.2f}, 50 draws on 200k cells); state rows are reported, never inputs"), fontsize=8)
    style.title(fig, f"P06 marker-by-label heatmap -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig, top=0.92)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "P06_marker_by_label_heatmap"); plt.close(fig)
    facts = {"labels": counts, "null_max_p95": null95, "n_identity": len(ident), "n_state": len(state)}
    files.append(_side(out, "P06_marker_by_label_heatmap", facts))
    return _rec("P06", donor, files, [pins.uniform_parquet(donor), *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ P07 spatial map
def p07_spatial(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("P07", donor, prov, "no PhenoCycler release yet")
    cells = pd.read_parquet(pins.cells_parquet(donor), columns=["object_id", "X_centroid", "Y_centroid"])
    cells["object_id"] = cells["object_id"].astype(str)
    df = cells.merge(calls[["object_id", "lineage_native"]], on="object_id", how="inner")
    common, colors = io.common_of(MOD), io.common_colors()
    labs = label_order(common)
    px = 20.0
    gx = (df["X_centroid"].to_numpy() / px).astype(int); gy = (df["Y_centroid"].to_numpy() / px).astype(int)
    W, Hh = gx.max() + 1, gy.max() + 1
    code = pd.Categorical(df["lineage_native"], categories=labs).codes
    raster = np.full((Hh, W), -1, dtype=int)
    counts = np.zeros((Hh, W, len(labs)), dtype=np.int32)
    np.add.at(counts, (gy, gx, code), 1)
    filled = counts.sum(axis=2) > 0
    raster[filled] = counts.argmax(axis=2)[filled]
    fig, axes = style.figure(1, 2, width=7, height=6)
    ax = axes[0, 0]
    rgb = np.ones((Hh, W, 3))
    for j, L in enumerate(labs):
        c = matplotlib.colors.to_rgb(io.label_color(L, common, colors))
        rgb[raster == j] = c
    ax.imshow(rgb, origin="upper", interpolation="nearest")
    ax.set_title(style.wrap(f"P07-A  20 um majority raster of lineage_native ({len(df):,} cells)"), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=f"{L} n={(df['lineage_native'] == L).sum():,}", hatch="///" if L == io.AMBIGUOUS else None) for L in labs], fontsize=6, loc="lower left")
    ax = axes[0, 1]
    reg_idx, reg = pheno_io.first_islet_region(donor)
    w = reg["window"]
    x0, y0, s = w["x0"] * pins.PIXEL_SIZE_UM, w["y0"] * pins.PIXEL_SIZE_UM, w["size"] * pins.PIXEL_SIZE_UM
    sel = (df["X_centroid"] >= x0) & (df["X_centroid"] < x0 + s) & (df["Y_centroid"] >= y0) & (df["Y_centroid"] < y0 + s)
    sub = df[sel]
    ax.scatter(sub["X_centroid"], sub["Y_centroid"], s=3, c=[io.label_color(L, common, colors) for L in sub["lineage_native"]], linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(style.wrap(f"P07-B  frame region {reg_idx} ({reg.get('draw')}), {len(sub):,} cells, {s:.0f} um window"), fontsize=9)
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    style.title(fig, f"P07 spatial map of the released broad labels -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "P07_spatial_map"); plt.close(fig)
    facts = {"labels": df["lineage_native"].value_counts().to_dict(), "region_idx": reg_idx}
    files.append(_side(out, "P07_spatial_map", facts))
    return _rec("P07", donor, files, [pins.cells_parquet(donor), *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ P08 UMAP diagnostic
def p08_umap(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("P08", donor, prov, "no PhenoCycler release yet")
    import umap
    tab = pd.read_parquet(pins.uniform_parquet(donor))
    tab["object_id"] = tab["object_id"].astype(str)
    dec = json.loads((pins.PKG / "scripts" / "broadtype" / "declarations.json").read_text())
    ident = [m for m in dec["identity_markers"] if m in tab.columns]
    df = tab[["object_id", *ident]].merge(calls[["object_id", "lineage_native"]], on="object_id", how="inner")
    df = df[~df["lineage_native"].eq(io.UNAVAILABLE)]
    gate = io.broadtype_lane()["gate"]
    idx = gate.hash_subsample(df["object_id"].to_numpy(), 200_000, salt="p08")
    sub = df.iloc[idx]
    X = np.arcsinh(np.clip(sub[ident].to_numpy(float), 0, None))
    X = (X - np.nanmean(X, axis=0)) / np.where(np.nanstd(X, axis=0) > 0, np.nanstd(X, axis=0), 1)
    X = np.nan_to_num(X)
    emb = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0, n_jobs=1).fit_transform(X)
    common, colors = io.common_of(MOD), io.common_colors()
    fig, axes = style.figure(1, 1, width=7, height=6.5)
    ax = axes[0, 0]
    for L in label_order(common):
        sel = sub["lineage_native"].to_numpy() == L
        if sel.any():
            ax.scatter(emb[sel, 0], emb[sel, 1], s=1, c=io.label_color(L, common, colors), label=f"{L} n={sel.sum():,}", linewidths=0, alpha=0.7)
    ax.legend(fontsize=6, markerscale=6, loc="best"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(style.wrap(f"P08  diagnostic UMAP of {len(sub):,} hash-sampled cells on {len(ident)} identity markers (asinh, standardised); coloured by released label -- a picture, not evidence"), fontsize=8)
    style.title(fig, f"P08 diagnostic UMAP -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "P08_umap_diagnostic"); plt.close(fig)
    facts = {"n_cells": int(len(sub)), "markers": ident, "umap": {"n_neighbors": 30, "min_dist": 0.3, "random_state": 0}}
    files.append(_side(out, "P08_umap_diagnostic", facts))
    return _rec("P08", donor, files, [pins.uniform_parquet(donor), *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ P09 proportions by donor (pair)
def p09_proportions(donors: list, prov: dict) -> FigureRecord:
    common, colors = io.common_of(MOD), io.common_colors()
    if not common:
        return _unavailable("P09", "_pair", prov, "no PhenoCycler release yet", section="_pair")
    labs = label_order(common)
    rows_all, rows_ev, ns = {}, {}, {}
    for d in donors:
        c = _calls(d)
        if c is None:
            continue
        name = f"{d} ({style.status_display(pins.donor_status(d))})"
        vc = c["lineage_native"].value_counts()
        rows_all[name] = {L: vc.get(L, 0) / len(c) for L in labs}
        ev = c[c["evidence_positive"].astype(bool)]
        rows_ev[name] = {L: (ev["lineage_native"] == L).sum() / max(len(ev), 1) for L in common}
        ns[name] = (len(c), len(ev))
    A = pd.DataFrame(rows_all).T; B = pd.DataFrame(rows_ev).T
    fig, axes = style.figure(1, 2, width=6.2, height=5)
    ax = axes[0, 0]; io.stacked_bars(ax, A, labs, common, colors)
    ax.set_ylabel("fraction of ALL run objects"); ax.set_title(style.wrap("P09-A  every object: labels + Ambiguous (hatched) + Unavailable (grey)"), fontsize=9)
    style.denominator(ax, "; ".join(f"{k}: n={v[0]:,}" for k, v in ns.items()), loc="upper right")
    ax = axes[0, 1]; io.stacked_bars(ax, B, list(common), common, colors)
    ax.set_ylabel("fraction of evidence-positive cells"); ax.set_title(style.wrap("P09-B  evidence-positive only (a VERIFIED label on that donor)"), fontsize=9)
    style.denominator(ax, "; ".join(f"{k}: n_ev={v[1]:,}" for k, v in ns.items()), loc="upper right")
    fig.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in labs], loc="lower center", bbox_to_anchor=(0.5, 0.05), ncol=5, fontsize=7)
    style.title(fig, f"P09 broad-label proportions by donor (ND / Aab+ n=0 / T1D) -- producer {io.producer(MOD)}")
    style.footer(fig, prov); style.finish(fig, bottom=0.21)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, "P09_proportions_by_donor"); plt.close(fig)
    facts = {"all": A.to_dict(), "evidence_positive": B.to_dict(), "n": ns, "note": "no Aab+ donor in this pair; the axis slot is empty"}
    files.append(_side(out, "P09_proportions_by_donor", facts))
    return _rec("P09", "_pair", files, io.release_inputs(MOD), facts, section="_pair")


# ------------------------------------------------------------------ P10 proportions by region
def p10_regions(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("P10", donor, prov, "no PhenoCycler release yet")
    common, colors = io.common_of(MOD), io.common_colors()
    labs = label_order(common)
    reg = calls["cell_region"].astype(str) if "cell_region" in calls.columns else pd.Series("tissue", index=calls.index)
    order = [r for r in ("tissue", "peri", "core") if r in set(reg)]
    T = pd.DataFrame({r: calls.loc[reg == r, "lineage_native"].value_counts(normalize=True).reindex(labs).fillna(0) for r in order}).T
    T.index = [f"{r}\nn={int((reg == r).sum()):,}" for r in order]
    fig, axes = style.figure(1, 2, width=6.2, height=5)
    ax = axes[0, 0]; io.stacked_bars(ax, T, labs, common, colors)
    ax.set_title(style.wrap("P10-A  by QuPath islet proximity (audit axis: Islet_N core / 20 um ring / tissue)"), fontsize=9); ax.set_ylabel("fraction of objects in region")
    ax = axes[0, 1]
    ct = io.cells_table(MOD, donor)
    zone_ok = ct is not None and "structure_id" in ct.columns and ct["structure_id"].astype(str).ne("").any()
    if zone_ok:
        ct = ct[["cell_id", "structure_id", "x_um", "y_um"]]
        ct["cell_id"] = ct["cell_id"].astype(str)
        m = calls[["object_id", "lineage_native"]].merge(ct, left_on="object_id", right_on="cell_id", how="inner")
        core = m["structure_id"].astype(str).ne("")
        # peri = within 20 um of any islet member cell (identity-derived), tissue = rest
        from scipy.spatial import cKDTree
        tree = cKDTree(m.loc[core, ["x_um", "y_um"]].to_numpy())
        d, _ = tree.query(m.loc[~core, ["x_um", "y_um"]].to_numpy(), distance_upper_bound=20.0)
        zone = np.where(core, "islet core", "tissue"); zone[np.flatnonzero(~core)[np.isfinite(d)]] = "peri (<=20 um)"
        Z = pd.DataFrame({z: m.loc[zone == z, "lineage_native"].value_counts(normalize=True).reindex(labs).fillna(0) for z in ("tissue", "peri (<=20 um)", "islet core")}).T
        Z.index = [f"{z}\nn={int((zone == z).sum()):,}" for z in ("tissue", "peri (<=20 um)", "islet core")]
        io.stacked_bars(ax, Z, labs, common, colors)
        ax.set_title(style.wrap("P10-B  by identity-derived islet zone (integration S3 structures)"), fontsize=9)
    else:
        ax.axis("off"); style.note(ax, "P10-B not available: identity-derived structures (integration S3) not yet on disk for this donor")
    fig.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in labs], loc="lower center", bbox_to_anchor=(0.5, 0.05), ncol=5, fontsize=7)
    style.title(fig, f"P10 broad-label proportions by islet proximity -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig, bottom=0.21)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "P10_proportions_by_region"); plt.close(fig)
    facts = {"by_cell_region": T.to_dict(), "identity_zones_available": bool(zone_ok)}
    files.append(_side(out, "P10_proportions_by_region", facts))
    return _rec("P10", donor, files, io.release_inputs(MOD), facts)


# ------------------------------------------------------------------ P11 Ki67 rank per label (no cut)
def p11_ki67(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None or "ki67_donor_rank" not in calls.columns:
        return _unavailable("P11", donor, prov, "no PhenoCycler release with ki67_donor_rank yet")
    common, colors = io.common_of(MOD), io.common_colors()
    labs = [L for L in label_order(common) if (calls["lineage_native"] == L).sum() >= 50]
    fig, axes = style.figure(1, 1, width=9, height=4.5)
    ax = axes[0, 0]
    data = [calls.loc[calls["lineage_native"] == L, "ki67_donor_rank"].dropna().to_numpy() for L in labs]
    parts = ax.violinplot(data, showmedians=True, widths=0.8)
    for body, L in zip(parts["bodies"], labs):
        body.set_facecolor(io.label_color(L, common, colors)); body.set_alpha(0.8)
    ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels([f"{L}\nn={len(v):,}" for L, v in zip(labs, data)], fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("Ki67 within-donor rank (0-1)"); ax.set_ylim(0, 1)
    ax.axhline(0.5, color="#999999", lw=0.8, ls=":")
    ax.set_title(style.wrap("P11  proliferation as a STATE: the continuous within-donor rank of nucleus Ki67 per released label; no cut, no positive call (plan: PCNA/Ki67 declaration open)"), fontsize=8)
    style.denominator(ax, "denominator per violin: cells of that label with a Ki67 measurement", loc="upper left")
    style.title(fig, f"P11 Ki67 rank by broad label -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "P11_ki67_rank"); plt.close(fig)
    facts = {L: {"n": int(len(v)), "median_rank": float(np.median(v)) if len(v) else None} for L, v in zip(labs, data)}
    files.append(_side(out, "P11_ki67_rank", facts))
    return _rec("P11", donor, files, io.release_inputs(MOD), facts)


PER_DONOR = {"P06": p06_marker_by_label, "P07": p07_spatial, "P08": p08_umap, "P10": p10_regions, "P11": p11_ki67}
PAIR = {"P09": p09_proportions}
