"""K01-K03: the documented CORAL/KRONOS2 run, laid out for reading.

Reads only ``data/development/kronos2_documented/`` -- the output of ``scripts/kronos2_doc/run.py``,
which is CORAL Tutorial 3 (per-cell KRONOS2 + the mean-marker baseline) followed by Tutorial 4
(neighbours 15 on the raw embeddings, Leiden 0.5, UMAP, seed 0).  No label is assigned anywhere here:
Tutorial 4 assigns none, and CORAL states that naming clusters from their marker signature is the
reader's step.  These figures are that signature.

    K01  UMAP per feature space, coloured by Leiden cluster, with cluster sizes
    K02  cluster x marker profile (mean-marker intensity, z across clusters) -- the readout
    K03  the clusters in space, one frame window per donor
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

import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord  # noqa: E402

SEC = "kronos2_doc"
ROOT = pins.REPO / "data" / "development" / "kronos2_documented"
SPACES = ("kronos2", "mean_marker")
NICE = {"kronos2": "KRONOS2 768-d cell embedding", "mean_marker": "mean-marker baseline (per-cell intensity)"}


def _rec(fid, donor, files, inputs, facts):
    return FigureRecord(figure_id=fid, donor=str(donor), section=SEC, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _side(out, stem, facts):
    p = out / f"{stem}.json"; p.write_text(json.dumps(facts, indent=1, default=str) + "\n"); return p


def _unavailable(fid, donor, prov, why):
    fig, axes = style.figure(1, 1, width=7, height=3)
    axes[0, 0].axis("off"); style.note(axes[0, 0], f"{fid}: not available\n{why}")
    style.title(fig, f"{fid} -- donor {donor}: not available"); style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, f"{fid}_not_available"); plt.close(fig)
    return _rec(fid, donor, files, [], {"not_available": why})


def load(donor: str) -> dict | None:
    d = ROOT / f"donor_id={donor}"
    rec = d / "coral" / "features_record.json"
    crec = d / "clusters" / "cluster_record.json"
    if not (rec.exists() and crec.exists()):
        return None
    out = {"features": json.loads(rec.read_text()), "clusters": json.loads(crec.read_text()), "dir": d}
    for s in SPACES:
        out[s] = pd.read_parquet(d / "clusters" / f"clusters__{s}.parquet")
        out[f"prof_{s}"] = pd.read_csv(d / "clusters" / f"profiles__{s}.csv")
    return out


def _caption(D: dict) -> str:
    f, c = D["features"], D["clusters"]
    return (f"CORAL Tutorial 3 + 4 as documented: {f['n_cells']:,} cells, {f['patch_size']} px cell-centred tiles "
            f"({f['patch_um']} um at {f['mpp']:.4f} um/px; the tutorial's own slide is 24 um), batch {f['batch_size_used']} "
            f"(CORAL's gold standard), then neighbours 15 on the raw embedding, Leiden {c['resolution']}, seed {c['seed']}")


# ------------------------------------------------------------------ K01 UMAPs
def k01_umap(donor: str, prov: dict) -> FigureRecord:
    D = load(donor)
    if D is None:
        return _unavailable("K01", donor, prov, "the documented run has not finished for this donor")
    fig, axes = style.figure(1, 2, width=6.2, height=5.6)
    facts = {}
    for ax, s in zip(axes[0], SPACES):
        df = D[s]; k = int(df["cluster"].max()) + 1
        cmap = matplotlib.colormaps["tab20"]
        for cl in range(k):
            m = df["cluster"].to_numpy() == cl
            ax.scatter(df.loc[m, "umap_1"], df.loc[m, "umap_2"], s=1.5, color=cmap(cl % 20), linewidths=0, alpha=0.75)
            cx, cy = df.loc[m, "umap_1"].median(), df.loc[m, "umap_2"].median()
            ax.text(cx, cy, str(cl), fontsize=8, weight="bold", ha="center", va="center")
        ax.set_xticks([]); ax.set_yticks([])
        sizes = df["cluster"].value_counts().sort_index()
        ax.set_title(style.wrap(f"K01  {NICE[s]}: {k} Leiden clusters on {len(df):,} cells. "
                                f"Sizes {', '.join(f'{i}:{n:,}' for i, n in sizes.items())}", 78), fontsize=7)
        facts[s] = {"n_clusters": k, "sizes": sizes.to_dict()}
    facts["kronos2_vs_mean_marker"] = D["clusters"]["kronos2_vs_mean_marker"]
    style.title(fig, f"K01 documented Leiden clusters -- donor {donor} ({style.status_display(pins.donor_status(donor))}); "
                     f"KRONOS2 vs the mean-marker baseline agree at ARI {D['clusters']['kronos2_vs_mean_marker']['ari']:.2f}")
    fig.text(0.5, 0.925, style.wrap(_caption(D), 150), ha="center", va="top", fontsize=6.5)
    style.footer(fig, prov, donor); style.finish(fig, top=0.87)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "K01_umap_clusters"); plt.close(fig)
    files.append(_side(out, "K01_umap_clusters", facts))
    return _rec("K01", donor, files, [D["dir"] / "clusters" / "cluster_record.json"], facts)


# ------------------------------------------------------------------ K02 cluster x marker profile
def k02_profiles(donor: str, prov: dict) -> FigureRecord:
    D = load(donor)
    if D is None:
        return _unavailable("K02", donor, prov, "the documented run has not finished for this donor")
    fig, axes = style.figure(1, 2, width=7, height=6.4)
    facts = {}
    for ax, s in zip(axes[0], SPACES):
        prof = D[f"prof_{s}"]
        zc = [c for c in prof.columns if c.startswith("z_")]
        markers = [c[2:] for c in zc]
        Z = prof[zc].to_numpy(float).T                      # markers x clusters
        n = prof["n_cells"].to_numpy(int)
        vmax = float(np.nanmax(np.abs(Z))) if np.isfinite(Z).any() else 1.0
        im = ax.imshow(Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_yticks(range(len(markers))); ax.set_yticklabels(markers, fontsize=6)
        ax.set_xticks(range(len(n))); ax.set_xticklabels([f"{i}\nn={v:,}" for i, v in enumerate(n)], fontsize=6)
        for j in range(Z.shape[1]):
            i = int(np.nanargmax(Z[:, j])) if np.isfinite(Z[:, j]).any() else 0
            ax.text(j, i, "*", ha="center", va="center", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="mean marker per cluster, z across clusters")
        ax.set_title(style.wrap(f"K02  clusters from the {NICE[s]}; '*' marks each cluster's top marker. "
                                f"Values are the mean-marker intensity CORAL calls 'the standard input to cell "
                                f"phenotyping'; no label is assigned.", 80), fontsize=7)
        facts[s] = {"markers": markers, "n_cells": n.tolist(),
                    "top_marker_per_cluster": {int(j): markers[int(np.nanargmax(Z[:, j]))] for j in range(Z.shape[1])}}
    style.title(fig, f"K02 what each documented cluster is made of -- donor {donor} ({style.status_display(pins.donor_status(donor))})")
    fig.text(0.5, 0.935, style.wrap(_caption(D), 150), ha="center", va="top", fontsize=6.5)
    style.footer(fig, prov, donor); style.finish(fig, top=0.89)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "K02_cluster_marker_profile"); plt.close(fig)
    files.append(_side(out, "K02_cluster_marker_profile", facts))
    for s in SPACES:
        (out / f"K02_profiles__{s}.csv").write_text((D["dir"] / "clusters" / f"profiles__{s}.csv").read_text())
        files.append(out / f"K02_profiles__{s}.csv")
    return _rec("K02", donor, files, [D["dir"] / "clusters" / f"profiles__{s}.csv" for s in SPACES], facts)


# ------------------------------------------------------------------ K03 clusters in space
def k03_spatial(donor: str, prov: dict) -> FigureRecord:
    D = load(donor)
    if D is None:
        return _unavailable("K03", donor, prov, "the documented run has not finished for this donor")
    import pheno_io
    cells = pd.read_parquet(pins.cells_parquet(donor), columns=["object_id", "X_centroid", "Y_centroid"])
    cells["object_id"] = cells["object_id"].astype(str)
    reg_idx, reg = pheno_io.first_islet_region(donor)
    fig, axes = style.figure(1, 2, width=6, height=5.6)
    facts = {"region_idx": reg_idx}
    for ax, s in zip(axes[0], SPACES):
        by = pd.read_parquet(D["dir"] / "clusters" / f"clusters__{s}__by_object_id.parquet")
        m = by[by["region_idx"] == reg_idx].merge(cells, on="object_id", how="inner")
        cmap = matplotlib.colormaps["tab20"]
        ax.scatter(m["X_centroid"], m["Y_centroid"], s=6, c=[cmap(int(c) % 20) for c in m["cluster"]], linewidths=0)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
        ax.set_title(style.wrap(f"K03  {NICE[s]}: frame region {reg_idx} ({reg.get('draw')}), {len(m):,} cells, "
                                f"coloured by Leiden cluster", 72), fontsize=7.5)
        facts[s] = m["cluster"].value_counts().sort_index().to_dict()
    style.title(fig, f"K03 documented clusters in space -- donor {donor} ({style.status_display(pins.donor_status(donor))})")
    style.footer(fig, prov, donor); style.finish(fig)
    out = pins.OUT / str(donor) / SEC; files = style.savefig(fig, out, "K03_spatial"); plt.close(fig)
    files.append(_side(out, "K03_spatial", facts))
    return _rec("K03", donor, files, [pins.cells_parquet(donor)], facts)


PER_DONOR = {"K01": k01_umap, "K02": k02_profiles, "K03": k03_spatial}
PAIR: dict = {}
