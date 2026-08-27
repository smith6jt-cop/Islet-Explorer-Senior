#!/usr/bin/env python
"""Per-donor Xenium cell typing: PCA -> KODAMA -> cluster -> null-gated annotation.

Replaces the per-cell identity-core argmax, which was measured and refuted: on donor 6539,
against gene sets matched on core size and per-gene detection decile, five of eight lineages
sat inside their null and two (Endocrine, Immune_B) did worse than random. At ~96% per-cell
zeros an argmax over absolute means reports which core happens to hold a high expresser, and
a matched random set has that property too. Pooling cells into groups first is what buys the
signal a 316-transcript cell does not have; :mod:`nullgate` then refuses any group whose best
core cannot beat its own matched null.

Two arms are always run. ``spatial`` lets a cell's neighbourhood inform its representation;
``nonspatial`` does not. For a study asking whether immune cells sit near islets, a type call
that already knows where the cell is makes the answer partly circular -- so ``nonspatial`` is
the one that may carry a spatial claim, and the difference between the arms is reported as
the measurement of how much location leaks into identity.

Known limit, stated narrowly because an earlier draft overstated it. The claim was "acinar
cells have no axis", on the reading that ``02_pancreas_exocrine`` is ductal. **That is wrong,
and the file says so**: of its eight genes, ``CTRL`` and ``CUZD1`` carry ``subtype_hint=acinar``
and both are on the 5,101-gene panel. Measured on donor 6539, all 317,610 cells, CUZD1 is
detected in 65.8% at mean count 10.69 and CTRL in 40.7%; ``CUZD1|CTRL`` covers 66.5% of cells
and the two supply 71% of the core's mean. So the Exocrine core is already an acinar call,
driven by CUZD1 -- not a ductal one.

What is true is narrower: the canonical *zymogen* set is absent (PRSS1, PRSS2, CPA1, CPB1,
CELA2A/3A, CTRB1/2, PNLIP, CLPS, AMY2A, REG1A, GP2, PTF1A, RBPJL -- checked against the panel),
as are INS/GCG/SST/PPY. So acinar identity here rests on two genes whose acinar specificity in
THIS data has not been measured. Measure it before leaning on the Exocrine share; do not
report the share as an abundance estimate, and do not restate the absent-zymogen fact as a
claim that the compartment cannot be resolved.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))          # core/phenocycler, for `phenocycler`

import xenio                                          # noqa: E402
from nullgate import annotate_groups                   # noqa: E402

DEV_ROOT = HERE.parents[3] / "data" / "development" / "xenium_celltype"
ARMS = ("nonspatial", "spatial")

#: Cells embedded per donor. KODAMA's cost grows as ~n^1.44 in this configuration -- measured
#: on 50 PCs at landmarks=500, M=5: 5k = 44.7 s, 20k = 316.6 s, 50k = 1222.2 s -- so the full
#: 317k-1.4M cells per donor is days per arm and the cohort is weeks.
#:
#: The deliverable is donor-level COMPOSITION, not a label on every cell, and a subsample
#: estimates that directly: at 50,000 cells a population at 1% frequency is measured from
#: ~500 cells (binomial SE 0.14 percentage points). So this is a sampling decision with a
#: computable error, reported as a CI on every share -- not a tuning knob, and not a
#: statement about what any frequency ought to be. Raise it if a rarer population matters
#: and accept the runtime; the error bars will show what it bought.
DEFAULT_MAX_CELLS = 50_000


def subsample(dd: "xenio.DonorData", max_cells: int, seed: int) -> tuple["xenio.DonorData", dict]:
    """Take at most `max_cells` cells at random, recording exactly what was taken."""
    n = len(dd.cell_id)
    if max_cells <= 0 or n <= max_cells:
        return dd, {"subsampled": False, "n_available": n, "n_used": n, "subsample_seed": seed}
    idx = np.sort(np.random.default_rng(seed).choice(n, max_cells, replace=False))
    import dataclasses
    sub = dataclasses.replace(dd, cell_id=dd.cell_id[idx],
                              counts=dd.counts[idx].tocsr(), xy=dd.xy[idx])
    return sub, {"subsampled": True, "n_available": n, "n_used": int(max_cells),
                 "subsample_seed": seed}


def composition(ann: pd.DataFrame, groups: np.ndarray) -> pd.DataFrame:
    """Cells per awarded lineage, with a binomial CI -- the donor-level deliverable.

    Groups that cleared no null contribute to ``Unavailable``; they are never redistributed
    across the lineages that did clear, which would manufacture the very shares the gate
    refused to award.
    """
    sizes = pd.Series(groups).value_counts()
    lab = ann.set_index("group").lineage.reindex(sizes.index).fillna("")
    lab = lab.replace("", "Unavailable")
    n = int(sizes.sum())
    tot = sizes.groupby(lab).sum()
    frac = tot / n
    se = np.sqrt(frac * (1 - frac) / n)
    return pd.DataFrame({"n_cells": tot, "frac": frac,
                         "ci_lo": (frac - 1.96 * se).clip(lower=0),
                         "ci_hi": (frac + 1.96 * se).clip(upper=1)}).reset_index(
                             names="lineage").sort_values("frac", ascending=False)


def prepare(dd: "xenio.DonorData", out: Path, *, n_comps: int = 50, seed: int = 0) -> dict:
    """lognorm -> PCA -> the three files run_kodama.R reads.

    PCA is what makes KODAMA tractable here: measured at identical n/landmarks/M, 5,101
    genes took 1583.2 s and 50 components took 44.7 s. It is unsupervised and
    identity-agnostic, so it cannot encode the answer -- but ``n_comps`` is a choice and is
    recorded here so a later run can vary it and show the calls do not depend on it.
    """
    from sklearn.decomposition import TruncatedSVD
    out.mkdir(parents=True, exist_ok=True)
    Xn = xenio.lognorm(dd.counts)
    svd = TruncatedSVD(n_components=n_comps, random_state=seed)
    pcs = svd.fit_transform(Xn)
    sio.mmwrite(str(out / "pcs.mtx"), sp.csr_matrix(pcs.astype("float64")))
    np.savetxt(out / "coords.csv", dd.xy, delimiter=",", header="x,y", comments="")
    pd.DataFrame({"cell_id": dd.cell_id}).to_csv(out / "cell_id.csv", index=False)
    meta = {"donor_id": dd.donor_id, "source": dd.source, "n_cells": int(len(dd.cell_id)),
            "n_genes": len(dd.genes), "n_comps": int(n_comps), "seed": seed,
            "explained_variance_ratio_sum": float(svd.explained_variance_ratio_.sum()),
            "n_dropped_no_spatial": dd.n_dropped_no_spatial}
    (out / "prep.json").write_text(json.dumps(meta, indent=2))
    return meta


def run_kodama(inp: Path, out: Path, arm: str, *, M: int, landmarks: int,
               n_cores: int, seed: int = 0) -> Path:
    cmd = ["Rscript", str(HERE / "run_kodama.R"), "--input", str(inp), "--output", str(out),
           "--arm", arm, "--M", str(M), "--landmarks", str(landmarks),
           "--n-cores", str(n_cores), "--seed", str(seed)]
    subprocess.run(cmd, check=True)
    emb = out / f"embedding_{arm}.csv"
    if not emb.exists():
        raise RuntimeError(f"run_kodama.R exited 0 but wrote no {emb.name}")
    return emb


def cluster(emb: pd.DataFrame, *, k: int, seed: int = 0) -> np.ndarray:
    """KMeans on the KODAMA embedding.

    `k` is granularity, not an answer: annotation is by null-gated enrichment afterwards, so
    over-clustering costs resolution and never invents a lineage. Several k are run and the
    stability of the resulting lineage shares is what gets reported.
    """
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, n_init=10, random_state=seed)
    return km.fit_predict(emb[["dim1", "dim2"]].to_numpy())


def annotate(dd, groups, tax, *, n_draws: int, seed: int) -> pd.DataFrame:
    cores = {lin: tax.lineage_genes(lin)
             for lin in sorted({v for v in tax.panel_to_lineage.values() if v})}
    # the gate compares the real core against random cores on the caller's scale; lognorm
    # keeps a per-cell depth difference from masquerading as a lineage contrast
    return annotate_groups(xenio.lognorm(dd.counts), dd.genes, groups, cores,
                           n_draws=n_draws, seed=seed)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--xenium-paths", type=Path,
                    default=HERE.parents[1] / "data" / "integration" / "xenium_paths.csv")
    ap.add_argument("--out", type=Path, default=DEV_ROOT)
    ap.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
    ap.add_argument("--k", type=int, nargs="+", default=[20, 30, 40])
    ap.add_argument("--n-comps", type=int, default=50)
    ap.add_argument("--M", type=int, default=20)
    ap.add_argument("--landmarks", type=int, default=2000)
    ap.add_argument("--n-cores", type=int, default=16)
    ap.add_argument("--n-draws", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=DEFAULT_MAX_CELLS,
                    help="cells embedded per donor; 0 = all (days per arm, see module docs)")
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)

    from phenocycler.config import load_config
    from phenocycler.integration.vocab import load_panel_taxonomy
    cfg = load_config(str(HERE.parents[1] / "config.ini"))
    tax = load_panel_taxonomy(cfg, "pancreas")

    dest = a.out / f"donor_id={a.donor}" / f"roi={a.roi}"
    dd = xenio.load_donor(a.donor, a.xenium_paths, roi=a.roi)
    dd, sub_meta = subsample(dd, a.max_cells, a.seed)
    meta = prepare(dd, dest / "input", n_comps=a.n_comps, seed=a.seed)
    meta.update(sub_meta)
    (dest / "input" / "prep.json").write_text(json.dumps(meta, indent=2))
    print(f"[prep] {meta['donor_id']} {meta['n_used']:,} of {meta['n_available']:,} cells "
          f"from {meta['source']}; {a.n_comps} PCs explain "
          f"{meta['explained_variance_ratio_sum']:.1%}")
    if a.prepare_only:
        return 0

    rows, comps = [], []
    for arm in a.arms:
        emb = pd.read_csv(run_kodama(dest / "input", dest, arm, M=a.M,
                                     landmarks=a.landmarks, n_cores=a.n_cores, seed=a.seed))
        emb = emb.set_index("cell_id").loc[dd.cell_id].reset_index()
        for k in a.k:
            groups = cluster(emb, k=k, seed=a.seed)
            ann = annotate(dd, groups, tax, n_draws=a.n_draws, seed=a.seed)
            ann.insert(0, "k", k); ann.insert(0, "arm", arm)
            rows.append(ann)
            comp = composition(ann, groups)
            comp.insert(0, "k", k); comp.insert(0, "arm", arm); comp.insert(0, "donor_id", a.donor)
            comps.append(comp)
            n_called = int((ann.lineage != "").sum())
            unavail = float(comp.loc[comp.lineage == "Unavailable", "frac"].sum())
            print(f"[annotate] arm={arm} k={k}: {n_called}/{len(ann)} groups cleared their "
                  f"null; Unavailable {unavail:.1%} of cells")
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(dest / "group_calls.csv", index=False)
    pd.concat(comps, ignore_index=True).to_csv(dest / "composition.csv", index=False)
    print(f"[done] -> {dest/'group_calls.csv'} and composition.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
