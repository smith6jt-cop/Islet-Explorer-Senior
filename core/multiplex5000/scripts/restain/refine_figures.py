#!/usr/bin/env python
"""What the candidate registrations did, as pictures rather than as a claim.

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/refine_figures.py --donors 6436 6442

Reads ``refine_<roi>.json`` and each donor's point sets. Three figures per donor:

``refine_candidates.png``   every candidate on one axis: the per-cell pair distance and the
                            fraction of Xenium nuclei within 2 um, each against its own null.
                            Losers are drawn, not dropped.
``refine_pairs.png``        the winning candidate's per-cell distance distribution against the
                            placement null, with the quantiles marked. This is the number the
                            lane owes: what an individual cell suffers, not a window median.
``refine_residual.png``     where what is left points, per 500 um bin. A field of arrows all
                            pointing one way is a translation still to remove; a field with no
                            direction is the detector's own precision.

Distances are distances, so every axis starts at 0 -- the same rule ``display.py`` enforces on
intensity, for the same reason.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                       # headless: no display on this box

import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import align as AL                                                # noqa: E402
import refine as RF                                               # noqa: E402
from run_register import OUT                                      # noqa: E402

from phenocycler.integration.transform import Transform           # noqa: E402

WIN = "#1B7837"          # the candidate that was kept
LOSE = "#999999"         # the candidates that were not
NULL = "#B2182B"         # the placement null
BIN_UM = 500.0


def candidates_figure(rec: dict, path: Path, bars: dict | None = None) -> Path:
    cands = rec["candidates"]
    names = [c["name"] for c in cands]
    p50 = [c["score"]["pair_distance_um"].get("p50", np.nan) for c in cands]
    p90 = [c["score"]["pair_distance_um"].get("p90", np.nan) for c in cands]
    hit = [100 * c["score"]["matched_within"]["2um"] for c in cands]
    nul = [100 * c["score"]["matched_within_null"]["2um"] for c in cands]
    colors = [WIN if n == rec["winner"] else LOSE for n in names]
    y = np.arange(len(names))

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 0.55 * len(names) + 2.6), sharey=True)
    axes[0].barh(y, p50, color=colors)
    axes[0].plot(p90, y, "|", color="black", markersize=11, label="p90")
    axes[0].set_xlim(left=0)
    axes[0].set_xlabel("per-cell distance to the nearest IF nucleus (um)")
    axes[0].set_title("pair distance: median bar, p90 tick", fontsize=10)
    import matplotlib.transforms as mtransforms
    blend = mtransforms.blended_transform_factory(axes[0].transData, axes[0].transAxes)
    for label, v, ls in (bars or {}).values():
        axes[0].axvline(v, color="#2166AC", linestyle=ls, linewidth=1.1)
        axes[0].text(v, 0.02, f" {label} {v:.2f}", color="#2166AC", fontsize=7,
                     rotation=90, va="bottom", ha="left", transform=blend)
    axes[0].legend(loc="lower right", fontsize=8, frameon=False)

    axes[1].barh(y, hit, color=colors)
    axes[1].plot(nul, y, "o", color=NULL, markersize=4, label="placement null")
    ceiling = 100 * rec.get("detection_ceiling", float("nan"))
    if np.isfinite(ceiling):
        # No registration can match a Xenium nucleus the IF detector never found, so this is
        # where the bars run out of room. Without it a 53% bar reads as a failure rather than as
        # 94% of what was available.
        axes[1].axvline(ceiling, color="#2166AC", linestyle="--", linewidth=1.1)
        axes[1].text(ceiling, 0.02, f" IF detection ceiling {ceiling:.0f}%", color="#2166AC",
                     fontsize=7, rotation=90, va="bottom", ha="left",
                     transform=mtransforms.blended_transform_factory(
                         axes[1].transData, axes[1].transAxes))
    axes[1].set_xlim(0, 100)
    axes[1].set_xlabel("% of Xenium nuclei with an IF nucleus within 2 um")
    axes[1].set_title("matched fraction, against chance", fontsize=10)
    axes[1].legend(loc="lower right", fontsize=8, frameon=False)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(names, fontsize=9)
    axes[0].invert_yaxis()
    fig.suptitle(f"donor {rec['donor_id']} - every candidate registration, scored alike "
                 f"(kept: {rec['winner']})", fontsize=11)
    fig.text(0.5, 0.005,
             f"{rec['n_xenium_nuclei']:,} Xenium nuclei vs {rec['n_if_nuclei']:,} detected IF "
             f"nuclei; pairs are mutual nearest neighbours within "
             f"{rec['constants']['PAIR_MAX_UM']:g} um; the null displaces the whole field by "
             f"{rec['constants']['NULL_OFFSET_UM']:g} um.",
             ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def pairs_figure(rec: dict, d_all: np.ndarray, d_null: np.ndarray, path: Path) -> Path:
    # 0.5% of Xenium nuclei sit hundreds of microns from any IF detection -- they are outside
    # the drawn IF tissue, so there is nothing there to be registered to and the distance is not
    # a registration error. Letting them set the axis would compress the part that is. The
    # window covers the null's own bulk, and the overflow is COUNTED on the panel rather than
    # cropped away.
    hi = float(max(10.0, np.percentile(d_null, 90)))
    over = int((d_all > hi).sum())
    bins = np.linspace(0.0, hi, 160)                 # from 0: a distance has a real floor
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.hist(np.clip(d_all, 0, hi), bins=bins, color=WIN, alpha=0.85,
            label=f"registered ({rec['winner']})")
    ax.hist(np.clip(d_null, 0, hi), bins=bins, histtype="step", color=NULL, linewidth=1.4,
            label=f"placement null (+{rec['constants']['NULL_OFFSET_UM']:g} um)")
    top = ax.get_ylim()[1]
    for k, (q, style) in enumerate(((50, "-"), (90, "--"), (99, ":"))):
        v = float(np.percentile(d_all, q))
        if v <= hi:
            ax.axvline(v, color="black", linestyle=style, linewidth=1.0)
            ax.text(v, top * (0.96 - 0.09 * k), f" p{q} {v:.2f} um", va="top", ha="left",
                    fontsize=8)
        else:
            ax.text(0.99, 0.96 - 0.09 * k, f"p{q} {v:.1f} um (off scale) ",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.set_xlim(0, hi)
    if over:
        ax.text(0.99, 0.62, f"{over:,} nuclei ({100*over/len(d_all):.2f}%) beyond {hi:.0f} um, "
                            f"piled at the right edge -- mostly nuclei outside the drawn IF "
                            f"tissue, where no IF nucleus was detected at all",
                transform=ax.transAxes, ha="right", fontsize=7, color="#444444")
    ax.set_xlabel("distance from a Xenium nucleus to the nearest IF nucleus (um)")
    ax.set_ylabel("Xenium nuclei")
    ax.set_title(f"donor {rec['donor_id']} - per-cell registration error, all "
                 f"{len(d_all):,} nuclei", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def residual_figure(rec: dict, pts: np.ndarray, res: np.ndarray, path: Path) -> Path:
    from scipy.stats import binned_statistic_2d

    x0, y0 = pts.min(axis=0)
    x1, y1 = pts.max(axis=0)
    xe = np.arange(x0, x1 + BIN_UM, BIN_UM)
    ye = np.arange(y0, y1 + BIN_UM, BIN_UM)
    cnt, *_ = binned_statistic_2d(pts[:, 1], pts[:, 0], None, "count", bins=[ye, xe])
    dx, *_ = binned_statistic_2d(pts[:, 1], pts[:, 0], res[:, 0], "median", bins=[ye, xe])
    dy, *_ = binned_statistic_2d(pts[:, 1], pts[:, 0], res[:, 1], "median", bins=[ye, xe])
    ok = (cnt >= 20) & np.isfinite(dx) & np.isfinite(dy)
    cy, cx = np.meshgrid(0.5 * (ye[:-1] + ye[1:]), 0.5 * (xe[:-1] + xe[1:]), indexing="ij")
    mag = np.hypot(dx, dy)

    # Arrow length is scaled so the 99th-percentile bin reaches about one bin width. Fixing the
    # scale instead would make a good registration draw an empty page and a bad one draw a mess,
    # which reads as "no data" rather than "no residual" -- so the FACTOR is printed instead.
    hi = float(np.percentile(mag[ok], 99)) if ok.any() else 1.0
    arrow_x = BIN_UM / max(hi, 1e-9)
    fig, ax = plt.subplots(figsize=(7.6, 9.6))
    q = ax.quiver(cx[ok], cy[ok], dx[ok], dy[ok], mag[ok], cmap="magma",
                  angles="xy", scale_units="xy", scale=1.0 / arrow_x, width=0.004)
    q.set_clim(0.0, hi)
    cb = fig.colorbar(q, ax=ax, shrink=0.55)
    cb.set_label(f"median residual in the bin (um), arrows x{arrow_x:.0f}")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um, IF frame)")
    ax.set_ylabel("y (um, IF frame)")
    med = float(np.median(np.linalg.norm(res, axis=1)))
    const = float(np.median(np.linalg.norm(res - np.median(res, axis=0), axis=1)))
    ax.set_title(f"donor {rec['donor_id']} - what {rec['winner']} leaves, per {BIN_UM:g} um bin\n"
                 f"median {med:.3f} um; {const:.3f} um after removing a constant "
                 f"({np.median(res[:, 0]):+.3f}, {np.median(res[:, 1]):+.3f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def figures_for(rec: dict, out_dir: Path) -> list[Path]:
    from scipy.spatial import cKDTree

    donor_dir = out_dir / "valis" / f"donor_id={rec['donor_id']}"
    xen = np.load(donor_dir / "xenium_nuclei_um.npy")
    ifn = np.load(donor_dir / "if_nuclei_um.npy")
    tf = Transform.load(donor_dir, stem="refined")

    warped = tf.apply(xen)
    tree = cKDTree(ifn)
    d_all, _ = tree.query(warped, k=1)
    d_null, _ = tree.query(warped + np.array([rec["constants"]["NULL_OFFSET_UM"], 0.0]), k=1)
    i, j, _d = RF.mutual_pairs(warped, ifn, max_dist_um=rec["constants"]["PAIR_MAX_UM"])

    # The three displacements at which a disc still covers 75% of its own footprint, for the
    # three objects a reader might mean by "a cell": the gate's own cell diameter, the NUCLEUS
    # this measurement actually matches, and the hormone sampler's 8 um aperture. Which one is
    # the bar is a decision, so all three are drawn rather than one being chosen here.
    bars = {k: (k, AL.self_overlap_offset_um(d_um, 0.75), ls)
            for k, d_um, ls in (("cell", rec.get("median_cell_diameter_um", 9.3), "--"),
                                ("nucleus", rec.get("median_nucleus_diameter_um", 6.15), ":"),
                                ("8um aperture", 8.0, "-."))}

    d = out_dir / "figures" / f"donor_id={rec['donor_id']}" / "refine"
    return [candidates_figure(rec, d / "refine_candidates.png", bars),
            pairs_figure(rec, d_all, d_null, d / "refine_pairs.png"),
            residual_figure(rec, warped[i], ifn[j] - warped[i], d / "refine_residual.png")]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    p = args.out / f"refine_{args.roi}.json"
    if not p.exists():
        print(f"{p} is absent -- run run_refine.py first", file=sys.stderr)
        return 1
    for rec in json.loads(p.read_text()):
        if args.donors and rec["donor_id"] not in set(args.donors):
            continue
        for f in figures_for(rec, args.out):
            print(f"-> {f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
