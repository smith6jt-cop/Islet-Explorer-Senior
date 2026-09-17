#!/usr/bin/env python
"""Score every detection mode, including the context-conditional one, against ground truth.

The comparison is the one `detect.py` established and is kept identical so the new mode is
measured on the same terms rather than on terms chosen for it:

* **Y floor** -- both donors are female, so `TBL1Y`, `USP9Y` and `EIF1AY` have no template of
  any kind and every count they register is background. `DDX3Y` and `KDM5D` cross-hybridise
  with their highly expressed X paralogs, so they measure something else and are reported
  separately rather than averaged in.
* **clean autosomal negatives** -- an EXPLICIT, per-gene curated list (``CLEAN_NEG``) of panel
  genes with no expression in pancreas in ANY compartment. Real gene probes, so what they pick
  up IS the gene-level background. See "Why ``CLEAN_NEG`` is a gene list" below -- this set is
  deliberately NOT derived from ``panel_audit.csv`` categories.
* **median gene** -- the sensitivity that background is being traded against.
* **selectivity** -- (median gene rate / clean negative rate), and the gain over `raw >= 1`.
  A mode that halves background and halves sensitivity has gained nothing, and this is the
  number that says so.

    python scripts/xenium_celltype/eval_detection.py --bundle <dir> --calls <parquet>

Why ``CLEAN_NEG`` is a gene list and NOT a ``panel_audit.csv`` category filter
-----------------------------------------------------------------------------
This set was previously built by substring-matching the ``category`` column of
``panel_audit.csv`` against ``CLEAN_NEG_CATEGORIES = ("skin", "nephron")``.  That selection
was WRONG: it admitted genes pancreas expresses, so their contribution to ``clean_neg`` was
real signal booked as background, and every clean_neg figure this script has ever produced is
contaminated by it.

A ``panel_audit.csv`` category records where a gene's CANONICAL biology sits -- it is a
panel-design annotation, not a statement that the gene is silent in pancreas.  The
``nephron_specific`` category alone contributed three genes that pancreas genuinely
expresses:

* ``AQP3`` -- ductal / centroacinar epithelium.  16.6% and 7.3% detection on PhenoCycler-free
  pancreas sections (7.2% / 11.9% in the audit's own ``detection_pct_0041323`` /
  ``detection_pct_0041326`` columns).  It is also a live Gamma-cell rule marker in the Xenium
  phenotyping pipeline (``tissue_markers/pancreas.py``), where an attempt to drop it on
  exactly this "the panel annotation says kidney" reasoning was made and reverted.
* ``UMOD`` -- 1.6% / 3.2% detection in the same audit columns, well clear of the Y_CLEAN
  optical floor measured on these sections.
* ``WT1``  -- mesothelium, which lines the pancreatic capsule, plus stellate cells.

``skin_cornification`` was the same failure: it carries ``KRT1`` and ``KRT14``, keratins that
pancreatic ductal epithelium expresses.

The consequence is not a small bias, it is a category error: once genes the tissue expresses
are in the set, ``clean_neg`` is no longer an estimate of background at all, it is a blend of
background and real signal.  ``clean_neg`` is the DENOMINATOR of ``selectivity``, so every
selectivity number this script printed was computed against a quantity that does not mean what
its name says.  The old selection was also thin -- on the base hAtlas_v1.1 5K panel it
resolved to just 7 genes, of which 2 (``UMOD``, ``WT1``) were contaminants -- so a single
expressed probe moved the whole denominator.  A category label cannot answer "is this gene
silent in pancreas"; only per-gene curation can.

``CLEAN_NEG`` below is that curation: photoreceptor, cornified epidermis, bone/cartilage,
cardiac sarcomere, spermatogenesis and oocyte genes with no expression in pancreas in any
compartment -- islet, acinar, ductal, stellate, endothelium, nerve, or infiltrating immune.
Do NOT reintroduce a category- or substring-based selection here.  Add or remove genes one at
a time, with a reason.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import context as CTX                                  # noqa: E402
import detect as DET                                   # noqa: E402
import gdna as GD                                      # noqa: E402

#: Zero-template on a female donor, and without a highly expressed X paralog to cross-hybridise
#: from. These are the honest floor.
Y_CLEAN = ("TBL1Y", "USP9Y", "EIF1AY")
#: Zero-template but cross-hybridising (DDX3X / KDM5C), so a different measurement.
Y_CROSSHYB = ("DDX3Y", "KDM5D")
#: Clean autosomal negatives: an EXPLICIT curated gene list, one gene at a time, NOT a
#: ``panel_audit.csv`` category filter -- see the module docstring for why the category
#: selection was wrong (it admitted AQP3, UMOD and WT1, which pancreas expresses).
#: Every gene here has no expression in pancreas in any compartment: islet, acinar, ductal,
#: stellate, endothelium, nerve, or infiltrating immune.  Real gene probes, so whatever they
#: register IS the gene-level background.
CLEAN_NEG = (
    # photoreceptor / retinal pigment epithelium
    "ABCA4", "ARR3", "CNGB1", "CNGB3", "GUCA1A", "GUCY2D", "IMPG1", "NR2E3", "OPN1LW",
    "PDE6H", "PRPH2", "RPE65", "VSX2",
    # cornified epidermis (NOT the KRT1/KRT14 keratins the old "skin" category swept in --
    # ductal epithelium expresses those)
    "DSG1", "DSG3", "TGM3",
    # bone / cartilage
    "COL2A1", "PHEX",
    # cardiac sarcomere / natriuretic peptides
    "TNNI3", "TNNT2", "NPPA", "NPPB", "TBX5",
    # spermatogenesis / cancer-testis
    "AKAP4", "BOLL", "DDX25", "DPY19L2", "MAGEA1", "MAGEA10", "MAGEA12", "MAGEA3", "MAGEA4",
    "MEIOB", "PIWIL1", "PRM1", "SPATA16", "SPATA22", "SPATA4", "STK31", "SYCP1", "TEX14",
    # oocyte / subcortical maternal complex
    "DAZL", "GDF9", "KHDC3L", "NLRP5", "NLRP7", "TLE6",
)


def positive_sets(spearman: pathlib.Path, genes: pd.Index) -> dict:
    """Per reference, the RNA counterparts of the markers its rule REQUIRES.

    Mandatory markers only. The OR-groups are the broad arms -- ``Imm-macrophage``'s second
    group is "any one of nine pan-immune markers", which would drag MS4A1 and MPO into a
    macrophage positive control and measure nothing.

    This is the sensitivity side of the ground truth: a cell the PROTEIN calls Endo-vascular
    should carry PECAM1 in RNA, and a rule that suppresses it has not gained selectivity, it
    has lost the signal.
    """
    sys.path.insert(0, str(HERE.parent / "celltype"))
    import rules as R                                  # noqa: PLC0415

    sp_ = pd.read_csv(spearman)
    m2g = {r.marker: [g for g in str(r.gene).split("+") if g in genes]
           for r in sp_.itertuples() if isinstance(r.gene, str)}
    out = {}
    for ref, (all_of, _groups) in R.requirements(R.load_rules()).items():
        g = sorted({x for m in all_of for x in m2g.get(m, [])})
        if g:
            out[ref] = g
    return out


def clean_negatives(genes: pd.Index) -> list:
    """``CLEAN_NEG`` restricted to the genes this bundle's panel actually carries.

    Intersection only -- no category matching, no audit CSV.  A panel that is missing some of
    these (the base hAtlas_v1.1 5K panel does not carry every one) simply measures background
    on fewer probes.
    """
    return sorted(set(CLEAN_NEG) & set(genes))


def rate(det: sp.csr_matrix, cols: np.ndarray) -> float:
    if len(cols) == 0:
        return float("nan")
    return float(np.asarray(det[:, cols].mean(axis=0)).ravel().mean())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True, type=pathlib.Path)
    ap.add_argument("--calls", required=True, type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args(argv)

    t0 = time.time()
    ad = sc.read_10x_h5(args.bundle / "cell_feature_matrix.h5")
    ad.var_names_make_unique()
    ad.obs_names = ad.obs_names.astype(str)
    X = sp.csr_matrix(ad.X)
    print(f"[eval] counts {X.shape[0]:,} x {X.shape[1]:,}, nnz {X.nnz:,} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # Stage 1's own reader: genomic_control_counts comes from cells.parquet, shrunk toward an
    # area-conditioned prior, and verify_against_vendor refuses a >5% disagreement with the
    # vendor's stated false-positive rate. Re-deriving it here would be a second opinion on a
    # number that is already checked.
    g = GD.summarise(args.bundle)
    lam = pd.Series(g["_lambda"], index=pd.Index(g["_cell_id"]).astype(str)) \
            .reindex(ad.obs_names).to_numpy()
    if not np.isfinite(lam).all():
        raise SystemExit("some cells in the matrix have no lambda from cells.parquet")
    print(f"[eval] lambda mean {lam.mean():.6f} over {g['n_probes']} probes; "
          f"vendor check {g['vendor_check']}", flush=True)

    calls = pd.read_parquet(args.calls)
    calls["cell_id"] = calls["cell_id"].astype(str)
    calls = calls.set_index("cell_id").reindex(ad.obs_names)
    calls["call"] = calls["call"].fillna(CTX.POOLED)
    labels = CTX.class_labels(calls)
    labels.index = pd.RangeIndex(len(labels))
    vc = labels.value_counts()
    print(f"[eval] {len(vc)} context classes; pooled {int(vc.get(CTX.POOLED, 0)):,} cells "
          f"({vc.get(CTX.POOLED, 0)/len(labels):.1%})", flush=True)

    genes = ad.var_names
    gidx = {g: i for i, g in enumerate(genes)}
    sets = {
        "y_clean": np.array([gidx[g] for g in Y_CLEAN if g in gidx]),
        "y_crosshyb": np.array([gidx[g] for g in Y_CROSSHYB if g in gidx]),
        "clean_neg": np.array([gidx[g] for g in clean_negatives(genes)]),
    }
    real = np.array([i for g, i in gidx.items()
                     if not g.startswith(("Genomic", "Negative", "Unassigned", "Deprecated"))])
    print(f"[eval] negative sets: " +
          ", ".join(f"{k} {len(v)}" for k, v in sets.items()), flush=True)

    # A spatial hold-out, the same device run_refine.py uses: pi and mu are FITTED on one
    # side of a median x split and the context rule is EVALUATED on the other. Without it the
    # prior is estimated from the very cells it then scores, and the positive control would
    # measure that circularity rather than the rule.
    cells = pd.read_parquet(args.bundle / "cells.parquet")
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id").reindex(ad.obs_names)
    xcol = next((c for c in ("x_centroid", "x_um", "x") if c in cells.columns), None)
    if xcol is None:
        raise SystemExit(f"cells.parquet has no x centroid column: {list(cells.columns)[:8]}")
    xs = cells[xcol].to_numpy()
    split = float(np.median(xs))
    left = xs < split
    pi_l, mu_l, meta = CTX.class_tables(X[left], lam[left], labels[left].reset_index(drop=True))
    pi_r, mu_r, _ = CTX.class_tables(X[~left], lam[~left], labels[~left].reset_index(drop=True))
    pi, mu, _ = CTX.class_tables(X, lam, labels)
    print(f"[eval] class tables; spatial split at x={split:.0f} um "
          f"({int(left.sum()):,} left / {int((~left).sum()):,} right) ({time.time()-t0:.0f}s)",
          flush=True)

    def context_heldout() -> sp.csr_matrix:
        """Each half scored with the prior fitted on the OTHER half."""
        out = sp.lil_matrix(X.shape, dtype=bool)
        for mask, pi_o, mu_o in ((left, pi_r, mu_r), (~left, pi_l, mu_l)):
            lab = labels[mask].reset_index(drop=True)
            keep = [c for c in lab.unique() if c in pi_o.index]
            lab = lab.where(lab.isin(keep), CTX.POOLED)
            if CTX.POOLED not in pi_o.index:        # pooled must exist to fall back to
                continue
            d = CTX.detected(X[mask], lam[mask], lab, pi_o, mu_o)
            out[np.flatnonzero(mask)] = d
        return sp.csr_matrix(out)

    pos = positive_sets(args.calls.parent / "protein_gene_spearman.csv", genes)
    pos_cols = {ref: np.array([gidx[g] for g in gs]) for ref, gs in pos.items()}
    print(f"[eval] positive control: {len(pos_cols)} references with required-marker genes",
          flush=True)

    def positive_rate(det: sp.csr_matrix) -> float:
        """Detection of a class's OWN required-marker genes, in cells of that class.

        Pooled across references weighted by class size, so a large class is not outvoted by
        a small one.
        """
        num = den = 0.0
        for ref, cols in pos_cols.items():
            rows = np.flatnonzero((labels == ref).to_numpy())
            if len(rows) < CTX.MIN_CLASS_CELLS:
                continue
            num += float(det[rows][:, cols].mean()) * len(rows)
            den += len(rows)
        return num / den if den else float("nan")

    modes = {
        "raw>=1": lambda: DET.detected(X, mode="raw", floor=1),
        "count>=2": lambda: DET.detected(X, mode="count", floor=2),
        "poisson a=.01": lambda: DET.detected(X, mode="poisson", gdna_lambda=lam),
        "lognorm log1p(1)": None,     # needs the normalised layer; filled below
        "context": lambda: CTX.detected(X, lam, labels, pi, mu),
        "context heldout": context_heldout,
    }
    ln = ad.copy()
    sc.pp.normalize_total(ln)
    sc.pp.log1p(ln)
    modes["lognorm log1p(1)"] = lambda: DET.detected(sp.csr_matrix(ln.X), mode="lognorm",
                                                     thr=float(np.log1p(1.0)))
    rows = []
    for name, fn in modes.items():
        det = fn()
        med = float(np.median(np.asarray(det[:, real].mean(axis=0)).ravel()))
        r = {"mode": name, "median_gene": med, "positive": positive_rate(det)}
        for k, cols in sets.items():
            r[k] = rate(det, cols)
        # Selectivity on the POSITIVE control, not the median gene: suppressing a gene a class
        # does not express is the intended behaviour, so counting it as lost sensitivity
        # measures the design as though it were a defect.
        r["selectivity"] = r["positive"] / r["clean_neg"] if r["clean_neg"] else float("nan")
        rows.append(r)
        print(f"  {name:<18} y_clean {r['y_clean']:.5%}  clean_neg {r['clean_neg']:.4%}  "
              f"median_gene {med:.4%}  POSITIVE {r['positive']:.4%}  "
              f"selectivity {r['selectivity']:.2f}", flush=True)

    df = pd.DataFrame(rows)
    base = df.loc[df["mode"] == "raw>=1", "selectivity"].iloc[0]
    df["selectivity_gain"] = df["selectivity"] / base
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    (args.out.with_suffix(".classes.json")).write_text(json.dumps(meta, indent=2, default=float))
    print("\n" + df.to_string(index=False))
    print(f"\n[eval] wrote {args.out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
