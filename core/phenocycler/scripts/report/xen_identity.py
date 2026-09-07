"""Xenium identity figures X04-X09 from the development release (data/development/xenium_typing/release).

Labels are read from the release; nothing here re-types a cell.  Identity cores come from the panel
snapshot (identity subpanels only -- state subpanels are never drawn as identity evidence).  QC-dropped
cells are drawn Unavailable.  A missing release renders a stated 'not available' panel.
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
import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord  # noqa: E402

sys.path.insert(0, str(pins.PKG / "scripts" / "xenium_typing"))
import pins as xpins  # noqa: E402  (scripts/xenium_typing/pins -- shadows nothing: imported under another name)

SEC = "xen"
MOD = "xenium"
CYCLING_RULE = ">= 2 of MKI67/TOP2A/PCNA at raw count >= 1 (analysis/xenium hierarchical_phenotype rule)"


def _rec(fid, donor, files, inputs, facts, section=SEC):
    return FigureRecord(figure_id=fid, donor=str(donor), section=section, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _side(out, stem, facts):
    p = out / f"{stem}.json"; p.write_text(json.dumps(facts, indent=1, default=str) + "\n"); return p


def _unavailable(fid: str, donor: str, prov: dict, why: str, section=SEC) -> FigureRecord:
    fig, axes = style.figure(1, 1, width=7, height=3)
    axes[0, 0].axis("off"); style.note(axes[0, 0], f"{fid}: not available\n{why}")
    style.title(fig, f"{fid} -- donor {donor}: not available"); style.footer(fig, prov, None if section == "_pair" else donor); style.finish(fig)
    out = pins.OUT / ("_pair" if section == "_pair" else str(donor)) / SEC
    files = style.savefig(fig, out, f"{fid}_not_available"); plt.close(fig)
    return _rec(fid, donor, files, [], {"not_available": why}, section=section)


def _calls(donor):
    return io.release_calls(MOD, donor)


def label_order(common):
    return list(common) + [io.AMBIGUOUS, io.UNAVAILABLE]


def _cells(donor: str) -> pd.DataFrame:
    c = pd.read_parquet(xpins.cells_parquet(donor), columns=["cell_id", "x_centroid", "y_centroid", "cell_area"])
    c["cell_id"] = c["cell_id"].astype(str)
    return c


# ------------------------------------------------------------------ X04 identity-core dotplot + nullgate
def x04_cores(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("X04", donor, prov, "no Xenium release yet (data/development/xenium_typing/release)")
    sys.path.insert(0, str(pins.PKG / "scripts" / "xenium_typing")); import gate as xgate  # noqa: E402
    base = xgate.load_base_for_gate(donor)
    cores = xpins.identity_cores()
    ids = base["obs_names"]; X = base["lognorm"]; genes = base["var_names"]
    lab = calls.set_index("cell_id")["lineage_native"].reindex(ids).fillna(io.UNAVAILABLE).to_numpy()
    common = io.common_of(MOD); labs = [L for L in label_order(common) if (lab == L).sum() >= 50]
    gi = {g: i for i, g in enumerate(genes)}
    core_names = list(cores)
    M = np.full((len(core_names), len(labs)), np.nan); F = np.full_like(M, np.nan)
    for i, cn in enumerate(core_names):
        cols = [gi[g] for g in cores[cn] if g in gi]
        if not cols:
            continue
        sub = X[:, cols]
        score = np.asarray(sub.mean(axis=1)).ravel()
        det = np.asarray((sub > 0).mean(axis=1)).ravel()
        for j, L in enumerate(labs):
            sel = lab == L
            M[i, j] = float(score[sel].mean()); F[i, j] = float(det[sel].mean())
    fig, axes = style.figure(1, 2, width=6.5, height=5.2)
    ax = axes[0, 0]
    mx = np.nanmax(M) if np.isfinite(M).any() else 1.0
    for i in range(len(core_names)):
        for j in range(len(labs)):
            if np.isfinite(M[i, j]):
                ax.scatter(j, i, s=20 + 380 * F[i, j], c=[[M[i, j] / mx, 0.2, 1 - M[i, j] / mx]], linewidths=0)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels([f"{L}\nn={(lab == L).sum():,}" for L in labs], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(core_names))); ax.set_yticklabels([f"{c} ({len(cores[c])} genes)" for c in core_names], fontsize=7)
    ax.invert_yaxis(); ax.set_xlim(-0.6, len(labs) - 0.4)
    ax.set_title(style.wrap("X04-A  identity-core mean lognorm (colour, 0 -> max) and detected fraction (size) per released label; state subpanels are not drawn"), fontsize=8)
    ax = axes[0, 1]
    arm = io.producer(MOD); v = io.verdict(MOD, arm, donor) if arm else None
    if v and v.get("nullgate", {}).get("evaluable"):
        per = v.get("lineages", {})
        names = [L for L in xpins.LINEAGES]
        real = [per.get(L, {}).get("award", {}).get("real", np.nan) for L in names]
        p95 = [per.get(L, {}).get("award", {}).get("null_p95", np.nan) for L in names]
        x = np.arange(len(names))
        ax.bar(x - 0.2, real, 0.4, color=[io.label_color(L, common, io.common_colors()) for L in names], label="real (own core on own cells)")
        ax.bar(x + 0.2, p95, 0.4, color="#BBBBBB", label="matched null p95")
        for k, L in enumerate(names):
            a = per.get(L, {}).get("award", {})
            ax.text(k, max(np.nan_to_num(real[k]), np.nan_to_num(p95[k])), "awarded" if a.get("awarded") else ("unjudgeable" if not a.get("judgeable", True) else "not awarded"), ha="center", va="bottom", fontsize=6, rotation=90)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7); ax.set_ylabel("nullgate contrast (standardised core mean)")
        ax.axhline(0, color="black", lw=0.6); ax.legend(fontsize=7, loc="upper left")
        ax.set_title(style.wrap(f"X04-B  null-gated award per lineage label ({v.get('null_draws')} draws; {arm})"), fontsize=8)
        facts_b = {L: per.get(L, {}).get("award", {}) for L in names}
    else:
        ax.axis("off"); style.note(ax, f"X04-B not available: no gate record for producer {arm!r} on donor {donor}")
        facts_b = {"not_available": "no gate record"}
    style.title(fig, f"X04 identity cores by released label -- donor {donor}, producer {arm}")
    style.footer(fig, prov, donor); style.finish(fig, top=0.9)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X04_identity_cores"); plt.close(fig)
    facts = {"labels": {L: int((lab == L).sum()) for L in labs}, "cores": {c: len(g) for c, g in cores.items()}, "award": facts_b,
             "identity_state_overlap_on_base": base["identity_state_overlap"]}
    files.append(_side(out, "X04_identity_cores", facts))
    return _rec("X04", donor, files, [pathlib.Path(base["path"]), *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ X05 diagnostic UMAP from the producer latent
def x05_umap(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("X05", donor, prov, "no Xenium release yet")
    arm = io.producer(MOD)
    jm = io.XENIUM_TYPING / "joint_model"
    sys.path.insert(0, str(pins.PKG / "scripts" / "xenium_typing")); import gate as xgate  # noqa: E402
    if arm == xpins.ARM_JOINT and (jm / "train_latent.npy").exists():
        Z = np.load(jm / "train_latent.npy"); obs = pd.read_parquet(jm / "train_obs.parquet")
        sel = obs["donor_id"].astype(str).to_numpy() == str(donor)
        ids = obs.index.astype(str).to_numpy()[sel]
        ids = np.array([i.rsplit("-", 1)[0] if i.count("-") > 1 else i for i in ids])
        Zd = Z[sel]; latent_name = "joint scVI latent (batch = run)"; inputs = [jm / "train_latent.npy", jm / "train_obs.parquet"]
    elif arm == xpins.ARM_SCANVI and xpins.sidecar_path(donor).exists():
        import h5py
        from anndata.io import read_elem
        with h5py.File(xpins.sidecar_path(donor), "r") as f:
            ids_all = np.asarray(read_elem(f["obs"]).index).astype(str)
            Z_all = np.asarray(read_elem(f["obsm/X_scvi_broad"]), dtype=np.float32)
        idx = xgate.hash_subsample(ids_all, 200_000, salt="x05")
        ids, Zd = ids_all[idx], Z_all[idx]; latent_name = "per-donor scANVI broad latent (X_scvi_broad, hash-sampled 200k)"; inputs = [xpins.sidecar_path(donor)]
    else:
        return _unavailable("X05", donor, prov, f"producer {arm!r} keeps no per-cell latent on disk for donor {donor}")
    lab = calls.set_index("cell_id")["lineage_native"].reindex(ids).fillna(io.UNAVAILABLE).to_numpy()
    import umap
    emb = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0, n_jobs=1).fit_transform(Zd)
    common, colors = io.common_of(MOD), io.common_colors()
    fig, axes = style.figure(1, 1, width=7, height=6.5); ax = axes[0, 0]
    for L in label_order(common):
        m = lab == L
        if m.any():
            ax.scatter(emb[m, 0], emb[m, 1], s=1.5, c=io.label_color(L, common, colors), label=f"{L} n={m.sum():,}", linewidths=0, alpha=0.7)
    ax.legend(fontsize=6, markerscale=6); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(style.wrap(f"X05  diagnostic UMAP of {len(Zd):,} cells of this donor in the producer's latent: {latent_name}; coloured by released label -- a picture, not evidence"), fontsize=8)
    style.title(fig, f"X05 diagnostic UMAP -- donor {donor}, producer {arm}"); style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X05_umap_diagnostic"); plt.close(fig)
    facts = {"n_cells": int(len(Zd)), "latent_dim": int(Zd.shape[1]), "latent": latent_name, "umap": {"n_neighbors": 30, "min_dist": 0.3, "random_state": 0}}
    files.append(_side(out, "X05_umap_diagnostic", facts))
    return _rec("X05", donor, files, [*inputs, *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ X06 spatial map
def x06_spatial(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("X06", donor, prov, "no Xenium release yet")
    cells = _cells(donor)
    df = cells.merge(calls[["cell_id", "lineage_native", "reason"]], on="cell_id", how="left")
    df["lineage_native"] = df["lineage_native"].fillna(io.UNAVAILABLE)
    common, colors = io.common_of(MOD), io.common_colors(); labs = label_order(common)
    px = 20.0
    gx = (df["x_centroid"].to_numpy() / px).astype(int); gy = (df["y_centroid"].to_numpy() / px).astype(int)
    W, H = gx.max() + 1, gy.max() + 1
    code = pd.Categorical(df["lineage_native"], categories=labs).codes
    counts = np.zeros((H, W, len(labs)), dtype=np.int32); np.add.at(counts, (gy, gx, code), 1)
    filled = counts.sum(axis=2) > 0; raster = np.full((H, W), -1); raster[filled] = counts.argmax(axis=2)[filled]
    rgb = np.ones((H, W, 3))
    for j, L in enumerate(labs):
        rgb[raster == j] = matplotlib.colors.to_rgb(io.label_color(L, common, colors))
    n_qc = int(df["reason"].fillna("").str.startswith("qc_dropped").sum())
    fig, axes = style.figure(1, 2, width=7, height=6)
    ax = axes[0, 0]; ax.imshow(rgb, origin="upper", interpolation="nearest"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(style.wrap(f"X06-A  20 um majority raster of lineage_native ({len(df):,} bundle cells; {n_qc:,} QC-dropped drawn Unavailable)"), fontsize=8)
    ax.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=f"{L} n={(df['lineage_native'] == L).sum():,}", hatch="///" if L == io.AMBIGUOUS else None) for L in labs], fontsize=6, loc="lower left")
    ax = axes[0, 1]
    # a 1 mm window centred on the densest 20 um bin of Endocrine calls, else the tissue centre
    endo = [L for L in labs if common.get(L) == "Endocrine"]
    if endo and (df["lineage_native"].isin(endo)).sum() > 0:
        e = df[df["lineage_native"].isin(endo)]
        cx, cy = e["x_centroid"].median(), e["y_centroid"].median()
    else:
        cx, cy = df["x_centroid"].median(), df["y_centroid"].median()
    half = 500.0
    s = df[(df["x_centroid"].between(cx - half, cx + half)) & (df["y_centroid"].between(cy - half, cy + half))]
    ax.scatter(s["x_centroid"], s["y_centroid"], s=2, c=[io.label_color(L, common, colors) for L in s["lineage_native"]], linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.set_title(style.wrap(f"X06-B  1 mm window at the median Endocrine position, {len(s):,} cells"), fontsize=8)
    style.title(fig, f"X06 spatial map of the released Xenium labels -- donor {donor}, producer {io.producer(MOD)}")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X06_spatial_map"); plt.close(fig)
    facts = {"labels": df["lineage_native"].value_counts().to_dict(), "n_qc_dropped": n_qc, "window_center_um": [float(cx), float(cy)]}
    files.append(_side(out, "X06_spatial_map", facts))
    return _rec("X06", donor, files, [xpins.cells_parquet(donor), *io.release_inputs(MOD)], facts)


# ------------------------------------------------------------------ X07 proportions by donor (pair)
def x07_proportions(donors: list, prov: dict) -> FigureRecord:
    common, colors = io.common_of(MOD), io.common_colors()
    if not common:
        return _unavailable("X07", "_pair", prov, "no Xenium release yet", section="_pair")
    labs = label_order(common)
    A, B, LS, ns = {}, {}, {}, {}
    for d in donors:
        c = _calls(d)
        if c is None:
            continue
        name = f"{d} ({style.status_display(pins.donor_status(d))})"
        vc = c["lineage_native"].value_counts(); A[name] = {L: vc.get(L, 0) / len(c) for L in labs}
        ev = c[c["evidence_positive"].astype(bool)]; B[name] = {L: (ev["lineage_native"] == L).sum() / max(len(ev), 1) for L in common}
        typed = c[c["lineage_native"].isin(list(common))]
        low = typed["reason"].astype(str).str.contains("low_signal")
        LS[name] = {L: float(low[typed["lineage_native"] == L].mean()) if (typed["lineage_native"] == L).sum() else np.nan for L in common}
        ns[name] = (len(c), len(ev), int(low.sum()))
    A, B, LS = pd.DataFrame(A).T, pd.DataFrame(B).T, pd.DataFrame(LS).T
    fig, axes = style.figure(1, 3, width=5, height=5)
    ax = axes[0, 0]; io.stacked_bars(ax, A, labs, common, colors); ax.set_ylabel("fraction of ALL bundle cells")
    ax.set_title(style.wrap("X07-A  every bundle cell (QC-dropped = Unavailable)"), fontsize=8); style.denominator(ax, "; ".join(f"{k}: n={v[0]:,}" for k, v in ns.items()))
    ax = axes[0, 1]; io.stacked_bars(ax, B, list(common), common, colors); ax.set_ylabel("fraction of evidence-positive cells")
    ax.set_title(style.wrap("X07-B  evidence-positive only (VERIFIED lineage, not low-signal)"), fontsize=8); style.denominator(ax, "; ".join(f"{k}: n_ev={v[1]:,}" for k, v in ns.items()))
    ax = axes[0, 2]
    x = np.arange(len(common)); w = 0.8 / max(len(LS), 1)
    for i, (name, row) in enumerate(LS.iterrows()):
        ax.bar(x + i * w, row.reindex(list(common)).to_numpy(float), w, color=style.status_color(pins.donor_status(name.split(" ")[0])), label=name)
    ax.set_xticks(x + w * (len(LS) - 1) / 2); ax.set_xticklabels(list(common), rotation=45, ha="right", fontsize=7); ax.set_ylim(0, 1)
    ax.set_ylabel("share of typed cells flagged low_signal"); ax.legend(fontsize=7); ax.set_title(style.wrap("X07-C  low-signal share per lineage"), fontsize=8)
    fig.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in labs], loc="lower center", bbox_to_anchor=(0.5, 0.05), ncol=5, fontsize=7)
    style.title(fig, f"X07 Xenium lineage proportions by donor (ND / Aab+ n=0 / T1D) -- producer {io.producer(MOD)}")
    style.footer(fig, prov); style.finish(fig, bottom=0.21)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, "X07_proportions_by_donor"); plt.close(fig)
    facts = {"all": A.to_dict(), "evidence_positive": B.to_dict(), "low_signal_share": LS.to_dict(), "n": ns}
    files.append(_side(out, "X07_proportions_by_donor", facts))
    return _rec("X07", "_pair", files, io.release_inputs(MOD), facts, section="_pair")


# ------------------------------------------------------------------ X08 proportions by identity-derived islet zone
def x08_zones(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None:
        return _unavailable("X08", donor, prov, "no Xenium release yet")
    ct = io.cells_table(MOD, donor)
    if ct is None or "structure_id" not in ct.columns or not ct["structure_id"].astype(str).ne("").any():
        return _unavailable("X08", donor, prov, "identity-derived Xenium structures (integration S3) not yet on disk for this donor; no QuPath annotation is used on Xenium")
    common, colors = io.common_of(MOD), io.common_colors(); labs = label_order(common)
    ct = ct[["cell_id", "structure_id", "x_um", "y_um"]].copy(); ct["cell_id"] = ct["cell_id"].astype(str)
    m = calls[["cell_id", "lineage_native"]].merge(ct, on="cell_id", how="inner")
    core = m["structure_id"].astype(str).ne("").to_numpy()
    from scipy.spatial import cKDTree
    tree = cKDTree(m.loc[core, ["x_um", "y_um"]].to_numpy())
    d, _ = tree.query(m.loc[~core, ["x_um", "y_um"]].to_numpy(), distance_upper_bound=20.0)
    zone = np.where(core, "islet core", "tissue"); zone[np.flatnonzero(~core)[np.isfinite(d)]] = "peri (<=20 um)"
    order = ("tissue", "peri (<=20 um)", "islet core")
    Z = pd.DataFrame({z: m.loc[zone == z, "lineage_native"].value_counts(normalize=True).reindex(labs).fillna(0) for z in order}).T
    Z.index = [f"{z}\nn={int((zone == z).sum()):,}" for z in order]
    fig, axes = style.figure(1, 1, width=6, height=5); ax = axes[0, 0]
    io.stacked_bars(ax, Z, labs, common, colors); ax.set_ylabel("fraction of cells in zone")
    ax.set_title(style.wrap("X08  by identity-derived islet zone (S3 Xenium structures; peri = within 20 um of a structure member)"), fontsize=8)
    fig.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in labs], loc="lower center", bbox_to_anchor=(0.5, 0.05), ncol=5, fontsize=7)
    style.title(fig, f"X08 Xenium proportions by islet zone -- donor {donor}, producer {io.producer(MOD)}"); style.footer(fig, prov, donor); style.finish(fig, bottom=0.21)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X08_proportions_by_zone"); plt.close(fig)
    facts = {"by_zone": Z.to_dict(), "n_structures": int(m.loc[core, "structure_id"].nunique())}
    files.append(_side(out, "X08_proportions_by_zone", facts))
    return _rec("X08", donor, files, io.release_inputs(MOD), facts)


# ------------------------------------------------------------------ X09 cycling per label
def x09_cycling(donor: str, prov: dict) -> FigureRecord:
    calls = _calls(donor)
    if calls is None or "cycling" not in calls.columns:
        return _unavailable("X09", donor, prov, "no Xenium release with a cycling column yet")
    common, colors = io.common_of(MOD), io.common_colors()
    labs = [L for L in label_order(common) if (calls["lineage_native"] == L).sum() >= 50 and L != io.UNAVAILABLE]
    typed = calls[calls["lineage_native"].isin(labs)]
    frac = typed.groupby("lineage_native")["cycling"].mean().reindex(labs)
    n = typed.groupby("lineage_native").size().reindex(labs)
    fig, axes = style.figure(1, 1, width=8, height=4.5); ax = axes[0, 0]
    ax.bar(range(len(labs)), frac.to_numpy(float), color=[io.label_color(L, common, colors) for L in labs])
    for i, (f, k) in enumerate(zip(frac, n)):
        ax.text(i, f, f"{f:.3f}\n(n={int(k):,})", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("fraction cycling"); ax.set_ylim(0, max(0.05, float(np.nanmax(frac)) * 1.3))
    ax.set_title(style.wrap(f"X09  cycling as a detection STATE per released label; rule: {CYCLING_RULE}"), fontsize=8)
    style.title(fig, f"X09 cycling per lineage -- donor {donor}, producer {io.producer(MOD)}"); style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "X09_cycling"); plt.close(fig)
    facts = {"rule": CYCLING_RULE, "fraction": frac.to_dict(), "n": n.astype(int).to_dict()}
    files.append(_side(out, "X09_cycling", facts))
    return _rec("X09", donor, files, io.release_inputs(MOD), facts)


PER_DONOR = {"X04": x04_cores, "X05": x05_umap, "X06": x06_spatial, "X08": x08_zones, "X09": x09_cycling}
PAIR = {"X07": x07_proportions}
