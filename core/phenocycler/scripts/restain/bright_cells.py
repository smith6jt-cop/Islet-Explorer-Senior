#!/usr/bin/env python
"""Start from the clearly bright IF cells. Do they carry the endocrine genes you would expect?

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/bright_cells.py --donors 6436 6442

**The measurement runs protein -> RNA, not RNA -> protein.** The maintainer can see the bright
insulin and glucagon cells directly in the image, so the protein is the unambiguous side and it
is what selects the cells. The question is then descriptive and needs no truth set: of the cells
that are clearly bright for a hormone, how many carry the endocrine genes associated with it,
and how does that compare with every other cell in the same section?

**The cut is on absolute counts, and the number of bright cells is an OUTPUT.** Selecting "the
top 1%" would impose how many hormone-positive cells a donor has, which is exactly the
assumption this repository forbids -- and it would impose it hardest on 6436, a T1D donor that
should have fewer insulin-bright cells than the ND donor. So the ladder is in the channel's own
8-bit counts and the resulting count is reported, never chosen. Nothing is tuned to it.

**No cell is labelled and no classifier is scored.** The gene side stays a raw detection rate:
the fraction of bright cells with at least one transcript, against the fraction among all other
cells. ``subtype_hint`` comes from ``core/panels``' endocrine identity core and is a hint in its
own file's words -- it groups the rows of a table, and asserts nothing about any cell.

The whole curve is reported rather than one threshold, because where "clearly bright" begins is
a judgement, and a reader should be able to see whether the answer depends on it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import hormone as HZ                                             # noqa: E402
import validate as VAL                                           # noqa: E402
import xenio                                                     # noqa: E402
from run_hormone import SECTIONS                                 # noqa: E402
from run_register import (GEOJSON_DIR, OUT, SUBPANELS,           # noqa: E402
                          XENIUM_PATHS, discover, donor_of)

#: Mean 8-bit counts inside a cell's own polygon, above which that cell is called clearly bright
#: for a channel. A ladder rather than a value: where "clearly bright" begins is a judgement, and
#: reporting the curve lets a reader see whether the answer turns on it. 8-bit data, so 100 is
#: 40% of full scale and nothing at that level is ambiguous.
BRIGHT_LADDER = (10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 140.0)


class BrightCellsError(RuntimeError):
    """The bright-cell comparison cannot be made from what is on disk."""


def detection_rates(counts, cols: list[int], mask: np.ndarray) -> np.ndarray:
    """Fraction of the cells in ``mask`` with at least one transcript, per gene in ``cols``."""
    import scipy.sparse as sp

    if not mask.any():
        return np.full(len(cols), np.nan)
    X = sp.csc_matrix(counts)[mask][:, cols]
    X.data = np.ones_like(X.data)
    return np.asarray(X.sum(axis=0)).ravel() / float(mask.sum())


def compare_one(donor: str, serial: str, *, roi: str = "panc", out_dir: Path = OUT,
                compartment: str = "cell", ladder=BRIGHT_LADDER) -> dict:
    dd = xenio.load_donor(donor, XENIUM_PATHS, roi=roi, prefer="bundle")
    _core, by_hint = VAL.endocrine_core(SUBPANELS)
    gi = {g: i for i, g in enumerate(dd.genes)}
    wanted = [(g, hint or "pan") for hint, genes in by_hint.items()
              for g in genes if g in gi]
    cols = [gi[g] for g, _ in wanted]

    rec = {"donor_id": donor, "xenium_serial": serial, "roi": roi,
           "compartment": compartment, "n_cells": int(len(dd.cell_id)),
           "direction": "protein selects the cells; RNA is then described",
           "selection": "mean 8-bit counts inside the cell's own Xenium polygon",
           "no_ground_truth": ("no cell is labelled and no classifier is scored; the gene side "
                               "is a raw detection rate"),
           "genes": [{"gene": g, "subtype_hint": h} for g, h in wanted],
           "channels": {}}

    for marker in ("INS", "GCG"):
        p = out_dir / "protein_rna" / f"{donor}_{serial}_{compartment}_{marker}.npy"
        if not p.exists():
            raise BrightCellsError(f"{p} is absent -- run protein_rna.py first")
        v = np.load(p)
        measured = np.isfinite(v)
        steps = []
        for thr in ladder:
            bright = measured & (v >= thr)
            rest = measured & ~bright
            r_b = detection_rates(dd.counts, cols, bright)
            r_r = detection_rates(dd.counts, cols, rest)
            per_gene = {}
            for (g, hint), b, r in zip(wanted, r_b, r_r):
                # None, never NaN: an empty bright set has no detection rate, and writing a
                # number there would make "nothing was measured" indistinguishable from zero.
                ok = np.isfinite(b) and np.isfinite(r)
                per_gene[g] = {"subtype_hint": hint,
                               "detected_in_bright": float(b) if np.isfinite(b) else None,
                               "detected_in_rest": float(r) if np.isfinite(r) else None,
                               "ratio": float(b / r) if (ok and r > 0) else None}
            def mean_ratio(h):
                vals = [x["ratio"] for x in per_gene.values()
                        if x["subtype_hint"] == h and x["ratio"] is not None
                        and np.isfinite(x["ratio"])]
                return float(np.mean(vals)) if vals else None
            steps.append({
                "threshold_counts": float(thr),
                "n_bright": int(bright.sum()),
                "pct_of_measured": float(100.0 * bright.sum() / max(measured.sum(), 1)),
                "genes": per_gene,
                "mean_ratio_beta": mean_ratio("beta"),
                "mean_ratio_alpha": mean_ratio("alpha"),
                "mean_ratio_pan": mean_ratio("pan"),
            })
        rec["channels"][marker] = {"n_measured": int(measured.sum()), "ladder": steps}

        print(f"\n[bright] {donor} ({'T1D' if donor == '6436' else ''}) {marker}: "
              f"{int(measured.sum()):,} cells measured", flush=True)
        print(f"[bright]   {'thr':>5} {'n bright':>9} {'% cells':>8}   "
              f"{'beta x':>7} {'alpha x':>8} {'pan x':>7}", flush=True)
        for st in steps:
            def f(x):
                return "  n/a" if x is None else f"{x:5.2f}"
            print(f"[bright]   {st['threshold_counts']:5.0f} {st['n_bright']:9,} "
                  f"{st['pct_of_measured']:7.3f}%   {f(st['mean_ratio_beta']):>7} "
                  f"{f(st['mean_ratio_alpha']):>8} {f(st['mean_ratio_pan']):>7}", flush=True)
    return rec


def figure(recs: list, out: Path) -> Path:
    """The whole ladder, so the answer can be seen not to depend on where "bright" begins.

    The bands are means over the subpanel's ``beta``- and ``alpha``-hinted genes. Both are large
    on every panel because the comparison set is the whole section, 97% of which is exocrine
    where an endocrine gene is near-absent -- so the height of a line is mostly "endocrine or
    not". **The subtype-informative part is the GAP between the two lines within one panel**,
    and that is why the panels are drawn per channel rather than pooled.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(recs), 2, figsize=(13.5, 5.4 * len(recs)), squeeze=False)
    for r, rec in enumerate(recs):
        for c, marker in enumerate(("INS", "GCG")):
            ax = axes[r][c]
            steps = rec["channels"][marker]["ladder"]
            thr = [st["threshold_counts"] for st in steps]
            for hint, colour, lab in (("beta", "#1c7ed6", "beta-hinted genes"),
                                      ("alpha", "#f76707", "alpha-hinted genes"),
                                      ("pan", "#adb5bd", "pan-endocrine")):
                y = [st[f"mean_ratio_{hint}"] for st in steps]
                xs = [t for t, v in zip(thr, y) if v is not None]
                ys = [v for v in y if v is not None]
                ax.plot(xs, ys, "o-", color=colour, lw=2.2 if hint != "pan" else 1.2,
                        label=lab, zorder=3 if hint != "pan" else 1)
            ax.set_yscale("log")
            ax.set_xlabel("cell is 'clearly bright' at ≥ this many 8-bit counts")
            if c == 0:
                ax.set_ylabel("detection rate in bright cells\n÷ rate in every other cell")
            for st in steps:
                if st["n_bright"]:
                    ax.annotate(f"{st['n_bright']:,}", (st["threshold_counts"], 1.05),
                                fontsize=7, rotation=90, ha="center", va="bottom",
                                color="#495057")
            ax.axhline(1.0, color="#212529", lw=0.8)
            ax.set_title(f"donor {rec['donor_id']} — cells bright for "
                         f"{'INS (Opal 690)' if marker == 'INS' else 'GCG (Opal 520)'}\n"
                         f"grey numbers = how many cells that is", fontsize=10)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 0 and c == 0:
                ax.legend(fontsize=8, frameon=False)
    fig.suptitle("Selected on the protein, described by the RNA — no cell is labelled\n"
                 "Threshold is on absolute counts, so the number of bright cells is an "
                 "output, never a target", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[bright] {out}", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--compartment", default="cell", choices=["cell", "nucleus"])
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    recs = []
    for serial, _gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        if sections[sections.xenium_serial == serial].empty:
            continue
        recs.append(compare_one(donor, serial, roi=args.roi, out_dir=args.out,
                                compartment=args.compartment))
    if not recs:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / f"bright_cells_{args.roi}_{args.compartment}.json"
    p.write_text(json.dumps(recs, indent=2, allow_nan=False))
    figure(recs, args.out / "protein_rna" / f"bright_cells_{args.compartment}.png")
    print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
