#!/usr/bin/env python3
"""Rebuild the Trajectory tab's H5AD on the WITHIN-DONOR ordering, for every donor it covers.

THE SINGLE PRODUCER of `<app>/adata_ins_root.h5ad`, as of 2026-09-18. Nothing else writes that
path. Four scripts in this tree name a file by that name and this is the only one that owns it:

* `scripts/reaggregate_islets.py` and `scripts/recompute_viz_umap.py` write
  `<repo>/data/adata_ins_root.h5ad` from `islet_analysis/` and
  `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` -- the legacy CODEX lane, and
  none of those three paths exists in this checkout, so neither can run.
* `core/multiplex5000/scripts/typing/trajectory_islet.py` used to target it and no longer does. It
  fits the POOLED islet-composition ordering, and on 2026-09-17 it REFUSED to write the deployed
  file: that ordering failed its own `pseudotime_is_not_mostly_donor` check (donor eta-squared
  0.444 within T1D, 0.314 within Aab+, against a 0.15 ceiling). It is the behaviour that script
  was built to have, so the split follows it -- trajectory_islet.py keeps `islet_leiden.parquet`
  and writes its own `islet_composition.h5ad` / `.rejected.h5ad`, and THIS script reads that
  rejected frame and writes the deployed file.

Before that split the deployed file was stuck at ten donors, two disease statuses and a superseded
typing generation, while the ordering the maintainer decided to ship -- the WITHIN-DONOR
single-cell trajectory, one destiny/slingshot fit per donor, no pooled embedding -- had been
rebuilt for all twenty. This script joins the two: it takes the twenty-donor FRAME out of the rejected file
(composition, ALR, the islet-composition UMAP) and replaces the one quantity that was rejected --
the pooled ordering -- with the within-donor one, summarised per islet.

WHAT IS AND IS NOT POOLED IN THE FILE THIS WRITES. Two constructions live in it and they must not
be read as one:

* `dpt_pseudotime` / `dpt_pseudotime_combined` are WITHIN-DONOR. Each donor was ordered on its own
  diffusion map; the stored value is that donor's own quantile, rank/(n+1). Two donors' values are
  not on one axis and a cohort-pooled comparison of them measures nothing.
* `obsm['X_umap']` is the POOLED islet-composition embedding, carried over unchanged from the
  islet-level lane, which measured and decided it separately. It is a description of composition,
  not an ordering and not a disease axis.

The legacy column NAMES are kept because the app reads them by name and its server module belongs
to another lane. What is in them is recorded in `uns['trajectory_construction']`, and the honest
names -- `sc_pseudotime_core`, `sc_pseudotime_peri`, `sc_pseudotime_peri_minus_core` -- are carried
beside them in `obs`.

NOTHING IS IMPUTED. An islet with no ordered cell in a region has no value there and carries the
reason it has none, in `pseudotime_unavailable_reason`. Zero is not a pseudotime.

    .venv/bin/python scripts/rebuild_adata_ins_root.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "core" / "multiplex5000"))

from phenocycler.cohort import DONOR_EXCLUSIONS, REINSTATED  # noqa: E402

TYPING = REPO / "data" / "development" / "typing"
TRAJ = TYPING / "trajectory"
APP = TYPING / "app"

#: The same floor the islet summary applies: below it an islet's "median" is one or two cells'
#: position, which is not a summary of anything. Restated here only so the combined core+peri
#: value this script computes uses the identical rule; it is not a second policy.
MIN_SAMPLED_CELLS = 5

#: Donor status order. ND < Aab+ < T1D, everywhere, including the categorical written to disk.
STATUS_ORDER = ("ND", "Aab+", "T1D")

#: Where the rejected pooled frame is looked for, in THIS order. An explicit list, not a glob and
#: not the newest mtime: `touch` must not be able to choose which frame the deployed file is built
#: from. The first name is what trajectory_islet.py writes now; the second is what it wrote before
#: 2026-09-18 and is the file currently on disk.
FRAME_CANDIDATES = (
    TYPING / "umap_rebuild_2026-09-17" / "islet_composition.rejected.h5ad",
    TYPING / "umap_rebuild_2026-09-17" / "adata_ins_root.rejected.h5ad",
)

#: The one path this script writes. Named here so the refusal below can check it.
OUT_NAME = "adata_ins_root.h5ad"


def resolve_frame() -> Path:
    """The first FRAME_CANDIDATES entry that exists, or a refusal naming all of them."""
    for candidate in FRAME_CANDIDATES:
        if candidate.exists():
            return candidate
    raise SystemExit(
        "no rejected pooled frame found. Looked, in order, at:\n  "
        + "\n  ".join(str(c) for c in FRAME_CANDIDATES)
        + "\nRun core/multiplex5000/scripts/typing/trajectory_islet.py --out-dir "
        + str(FRAME_CANDIDATES[0].parent) + " first, or pass --frame."
    )


def _norm_donor(value) -> str:
    text = str(value).strip()
    text = text.lstrip("0") or "0"
    return text


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float | None:
    """Between-group sum of squares over total. None when it is not defined."""
    ok = np.isfinite(values)
    values, groups = values[ok], groups[ok]
    if values.size < 2:
        return None
    levels = pd.unique(groups)
    if len(levels) < 2:
        return None
    grand = values.mean()
    total = float(((values - grand) ** 2).sum())
    if total <= 0:
        return None
    between = 0.0
    for level in levels:
        here = values[groups == level]
        between += here.size * (here.mean() - grand) ** 2
    return float(between / total)


def combined_medians(pseudotime_path: Path) -> pd.DataFrame:
    """Median within-donor quantile over an islet's ordered core AND peri cells together.

    The app's "Core + Peri" mode needs one value per islet over both regions. The islet summary
    carries the regions apart, so it is computed here from the same ordered cells, under the same
    sampled-cell floor. Cells outside the subsample were never ordered and are simply absent --
    they are not counted as anything.

    TWO conditions, not one, and the second is the point. A median over the pooled cells needs at
    least MIN_SAMPLED_CELLS of them, as anywhere else. It ALSO needs at least one ordered cell in
    EACH region, because a value computed from twelve peri cells and no core cell is a peri value,
    and labelling it "Core + Peri" would tell the reader it contains a compartment it does not.
    Measured on this build: the cell-count rule alone admits 15,265 islets, the two rules together
    admit 13,423, and of the 1,842 that the second rule turns away, 1,840 have no ordered core cell
    at all and 2 have no ordered peri cell. Those 1,840 are the ones that would have been labelled
    "Core + Peri" while containing no core.
    """
    cells = pd.read_parquet(pseudotime_path,
                            columns=["donor_id", "islet_key", "region", "pseudotime"])
    cells = cells[np.isfinite(cells["pseudotime"])]
    cells = cells[cells["region"].isin(["core", "peri"])]
    cells["case_id"] = cells["donor_id"].map(_norm_donor)
    grouped = cells.groupby(["case_id", "islet_key"], observed=True)
    out = grouped["pseudotime"].agg(["median", "size"]).rename(
        columns={"median": "sc_pseudotime_combined", "size": "sc_n_combined"})
    regions = grouped["region"].nunique().rename("sc_n_regions_combined")
    out = out.join(regions)
    unusable = ((out["sc_n_combined"] < MIN_SAMPLED_CELLS)
                | (out["sc_n_regions_combined"] < 2))
    out.loc[unusable, "sc_pseudotime_combined"] = np.nan
    return out.reset_index()


def unavailable_reason(row) -> str:
    """The named reason an islet has no core ordering. Never blank when the value is missing."""
    if not row["_in_summary"]:
        return "no_ordered_cell_in_either_region_for_this_islet"
    n_core = row["sc_n_core"]
    if not np.isfinite(n_core) or n_core <= 0:
        return "no_core_cell_in_this_donor_ordered_subsample"
    if n_core < MIN_SAMPLED_CELLS:
        return f"below_{MIN_SAMPLED_CELLS}_sampled_core_cell_floor"
    return "core_pseudotime_not_finite"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frame", type=Path, default=None,
                        help="the 20-donor islet frame whose pooled ordering is being replaced "
                             "(default: the first of FRAME_CANDIDATES that exists)")
    parser.add_argument("--summary", type=Path, default=TRAJ / "singlecell_islet_summary.parquet")
    parser.add_argument("--cells", type=Path, default=TRAJ / "singlecell_pseudotime.parquet")
    parser.add_argument("--fit", type=Path, default=TRAJ / "singlecell_trajectory.json")
    parser.add_argument("--relation", type=Path, default=TRAJ / "crossdonor_relation.json")
    parser.add_argument("--marker-means", type=Path, default=APP / "islet_marker_means.parquet")
    parser.add_argument("--out", type=Path, default=APP / OUT_NAME)
    parser.add_argument("--no-backup", action="store_true",
                        help="do not copy the existing output aside first")
    args = parser.parse_args()
    if args.frame is None:
        args.frame = resolve_frame()

    fit = json.loads(args.fit.read_text())
    if fit.get("fit_scope") != "within_donor":
        raise SystemExit(f"refusing: {args.fit} records fit_scope={fit.get('fit_scope')!r}, "
                         "and a pooled ordering is not what this file is for")

    frame = ad.read_h5ad(args.frame)
    print(f"frame   {args.frame.name}: {frame.n_obs:,} islets, "
          f"{frame.obs['imageid'].nunique()} donors")

    before = None
    if args.out.exists():
        old = ad.read_h5ad(args.out)
        before = {
            "n_obs": int(old.n_obs),
            "imageid": [str(x) for x in pd.Index(old.obs["imageid"]).unique().sort_values()],
            "donor_status": [str(x) for x in pd.Index(old.obs["donor_status"]).unique().sort_values()],
            "reasons": {str(k): int(v) for k, v in
                        old.obs["pseudotime_unavailable_reason"].value_counts().items()},
        }
        print(f"\nBEFORE  {args.out}")
        print(f"  islets        {before['n_obs']:,}")
        print(f"  imageid       {len(before['imageid'])}: {', '.join(before['imageid'])}")
        print(f"  donor_status  {', '.join(before['donor_status'])}")
        del old
        if not args.no_backup:
            stamp = datetime.fromtimestamp(args.out.stat().st_mtime).strftime("%Y%m%d-%H%M")
            keep = args.out.with_name(f"{args.out.stem}.superseded-{stamp}.h5ad")
            if not keep.exists():
                shutil.copy2(args.out, keep)
            print(f"  kept aside    {keep.name}")

    obs = frame.obs.copy()
    obs["case_id"] = obs["imageid"].map(_norm_donor).astype(str)
    obs["islet_key"] = obs["base_islet_id"].astype(str)

    # ---- the within-donor ordering, per islet -------------------------------------------------
    summary = pd.read_parquet(args.summary)
    summary["case_id"] = summary["case_id"].map(_norm_donor).astype(str)
    summary["_in_summary"] = True
    joined = obs.merge(summary, on=["case_id", "islet_key"], how="left", validate="one_to_one")
    joined["_in_summary"] = joined["_in_summary"].fillna(False).astype(bool)
    joined.index = obs.index

    combined = combined_medians(args.cells)
    joined = joined.merge(combined, on=["case_id", "islet_key"], how="left", validate="one_to_one")
    joined.index = obs.index

    for column in ("sc_pseudotime_core", "sc_pseudotime_peri", "sc_pseudotime_peri_minus_core",
                   "sc_pseudotime_mean_core", "sc_pseudotime_mean_peri",
                   "sc_n_core", "sc_n_peri", "sc_pseudotime_combined", "sc_n_combined",
                   "sc_n_regions_combined"):
        obs[column] = joined[column].to_numpy(dtype=float)

    obs["dpt_pseudotime"] = obs["sc_pseudotime_core"].to_numpy(dtype=float)
    obs["dpt_pseudotime_combined"] = obs["sc_pseudotime_combined"].to_numpy(dtype=float)

    missing = ~np.isfinite(obs["dpt_pseudotime"].to_numpy(dtype=float))
    reasons = np.array([""] * len(obs), dtype=object)
    if missing.any():
        reasons[missing] = joined.loc[missing].apply(unavailable_reason, axis=1).to_numpy()
    obs["pseudotime_unavailable_reason"] = pd.Categorical(reasons)

    # ---- donor status, in the cohort's order --------------------------------------------------
    present = [s for s in STATUS_ORDER if s in set(obs["donor_status"].astype(str))]
    unexpected = sorted(set(obs["donor_status"].astype(str)) - set(STATUS_ORDER))
    if unexpected:
        raise SystemExit(f"unrecognised donor_status values: {unexpected}")
    obs["donor_status"] = pd.Categorical(obs["donor_status"].astype(str),
                                         categories=present, ordered=True)

    excluded = sorted(set(obs["case_id"]) & set(DONOR_EXCLUSIONS))
    if excluded:
        raise SystemExit(f"centrally excluded donor present: {excluded}")

    out = ad.AnnData(X=frame.X.copy(), obs=obs.drop(columns=["case_id", "islet_key"]),
                     var=frame.var.copy())
    out.obs_names = frame.obs_names
    for key, value in frame.obsm.items():
        out.obsm[key] = value

    # ---- per-islet marker means, which the rejected frame does not carry ----------------------
    means = pd.read_parquet(args.marker_means)
    means["case_id"] = means["case_id"].map(_norm_donor).astype(str)
    means = means.set_index(means["case_id"] + "_" + means["islet_key"].astype(str))
    marker_cols = [c for c in means.columns if c.startswith("marker_mean_")]
    aligned = means.reindex(out.obs_names)
    matched = int(aligned[marker_cols[0]].notna().sum()) if marker_cols else 0
    out.obsm["X_marker_mean"] = aligned[marker_cols].to_numpy(dtype=float)
    out.uns["marker_mean_names"] = np.array(marker_cols, dtype=object)
    print(f"\nmarker means: {matched:,} of {out.n_obs:,} islets matched, "
          f"{len(marker_cols)} columns; unmatched islets carry NaN, not zero")

    out.uns["phenotype_names"] = frame.uns["phenotype_names"]
    if "embedding_batch" in frame.uns:
        out.uns["embedding_batch"] = frame.uns["embedding_batch"]

    # ---- what this file is, stated in it ------------------------------------------------------
    relation = json.loads(args.relation.read_text()) if args.relation.exists() else {}
    out.uns["trajectory_construction"] = {
        "fit_scope": "within_donor",
        "algorithm": fit["algorithm"],
        "source": fit["source"],
        "sampling": fit["sampling"],
        "pseudotime_column": fit["pseudotime_column"],
        "orientation": fit["orientation"],
        "n_population_cells": int(fit["n_population"]),
        "n_ordered_cells": int(fit["n_ordered"]),
        "n_donors_ordered": int(fit["n_donors_ordered"]),
        "min_sampled_cells_per_islet_region": MIN_SAMPLED_CELLS,
        "umap_is_not_the_ordering": (
            "obsm['X_umap'] is the POOLED islet-composition embedding from the islet-level lane, "
            "carried over unchanged. The ordering in dpt_pseudotime is within-donor and shares no "
            "construction with it. The UMAP describes composition and is not a disease axis."),
        "cross_donor_comparability": (
            "Each donor's pseudotime is its own quantile. Values from two donors are not on one "
            "axis; peri-minus-core is a within-islet difference and is the quantity that carries "
            "across donors."),
        "typing_generation": "data/development/typing/cells/",
        "donor_exclusions": dict(DONOR_EXCLUSIONS),
        "reinstated": dict(REINSTATED),
    }
    if relation:
        out.uns["crossdonor_relation"] = {
            "what_this_is": relation["what_this_is"],
            "unit_permuted": relation["unit_permuted"],
            "n_permutations": int(relation["n_permutations"]),
            "relation": relation["relation"],
        }

    # ---- the method record the app's caption reads --------------------------------------------
    pt = out.obs["dpt_pseudotime"].to_numpy(dtype=float)
    donors = out.obs["imageid"].astype(str).to_numpy()
    statuses = out.obs["donor_status"].astype(str).to_numpy()
    per_status_eta = {}
    per_status_donors = {}
    for status in present:
        here = statuses == status
        value = eta_squared(pt[here], donors[here])
        per_status_eta[status] = value
        per_status_donors[status] = int(len(set(donors[here])))
    ordered_mask = np.isfinite(pt)

    beta = None
    if "Beta_cell" in list(out.var_names):
        counts = np.asarray(out.X, dtype=float)
        total = counts.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            proportion = counts[:, list(out.var_names).index("Beta_cell")] / total
        both = ordered_mask & np.isfinite(proportion)
        if both.sum() > 2:
            from scipy import stats as _stats
            beta = float(_stats.spearmanr(pt[both], proportion[both]).statistic)

    out.uns["trajectory_validation"] = {
        "fit_scope": "within_donor",
        "gate_every_donor_has_its_own_fit":
            bool(int(fit["n_donors_ordered"]) == int(fit["n_donors_population"])),
        "gate_no_pooled_embedding": True,
        "gate_ordered_islets_all_finite": bool(np.isfinite(pt[ordered_mask]).all()),
        "gate_every_missing_value_carries_a_reason":
            bool((out.obs["pseudotime_unavailable_reason"].astype(str).str.len() > 0)
                 .to_numpy()[~ordered_mask].all()),
        "hard_pass": None,  # filled below
        "islets_ordered": int(ordered_mask.sum()),
        "islets_unavailable": int((~ordered_mask).sum()),
        "islets_ordered_but_absent_from_this_frame": None,  # filled below
        "islets_with_paired_peri_minus_core":
            int(np.isfinite(out.obs["sc_pseudotime_peri_minus_core"].to_numpy(dtype=float)).sum()),
        "n_population_cells": int(fit["n_population"]),
        "n_ordered_cells": int(fit["n_ordered"]),
        "sampled_fraction_cells": float(fit["n_ordered"]) / float(fit["n_population"]),
        "spearman_pseudotime_beta_proportion": beta,
        "within_status_donor_count": per_status_donors,
        # DELIBERATELY NOT WRITTEN UNDER `within_status_donor_eta2`. The app's caption prints that
        # key under the words "within-status donor eta-squared", which on a POOLED ordering is the
        # batch check with a 0.15 ceiling. On a within-donor fit it is not that check and must not
        # be printed as though it were: every donor is rank-uniform on its OWN axis, so what this
        # number measures is where each donor's core cells sit inside that donor's own spread --
        # a fact about the summarising, not about batch. Printed beside a "gates pass" badge it
        # would read as a ceiling that was cleared. It is carried in full under its own name, and
        # the cross-donor question is asked separately and properly, on the relation between the
        # twenty orderings with a permutation null over donors -- uns['crossdonor_relation'].
        "within_status_donor_eta2_of_islet_medians":
            {k: v for k, v in per_status_eta.items() if v is not None},
        "within_status_donor_eta2": None,
        "within_status_donor_eta2_reason": (
            "not_applicable_to_a_within_donor_fit(each_donor_is_ordered_on_its_own_axis); "
            "see within_status_donor_eta2_of_islet_medians and crossdonor_relation"),
        "batch_correction": "none",
        "batch_correction_reason": (
            "No correction is applied and none is needed at this level: the ordering is fit "
            "within donor, so no cross-donor axis was constructed for a batch effect to live on. "
            "Acquisition batch is aliased with disease status at the coarse level in this cohort "
            "and is reported, never corrected on."),
        "cohort_limit": (
            f"{', '.join(f'{per_status_donors[s]} {s}' for s in present)} donors. Two of the "
            "three Aab+ donors are in batch 2, which contains no ND donor, and 6521 is the only "
            "Aab+ donor in batch 1. The cohort supports SHOWING Aab+ islets; it does not support "
            "an ND < Aab+ < T1D ordering claim, and this figure is not one."),
    }
    # Islets the within-donor ordering covers that this islet-composition frame does not contain.
    # They are measured and they exist; what they are not is drawable here. Saying so beats letting
    # them vanish in a join.
    frame_keys = set(obs["case_id"].astype(str) + "_" + obs["islet_key"].astype(str))
    summary_keys = set(summary["case_id"].astype(str) + "_" + summary["islet_key"].astype(str))
    out.uns["trajectory_validation"]["islets_ordered_but_absent_from_this_frame"] = \
        int(len(summary_keys - frame_keys))

    gates = [v for k, v in out.uns["trajectory_validation"].items()
             if k.startswith("gate_") and v is not None]
    out.uns["trajectory_validation"]["hard_pass"] = bool(all(gates))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(args.out)

    print(f"\nAFTER   {args.out}")
    print(f"  islets        {out.n_obs:,}")
    after_ids = sorted(set(out.obs['imageid'].astype(str)), key=lambda s: (len(s), s))
    print(f"  imageid       {len(after_ids)}: {', '.join(after_ids)}")
    print(f"  donor_status  {', '.join(list(out.obs['donor_status'].cat.categories))} "
          f"(ordered={out.obs['donor_status'].cat.ordered})")
    print(f"                {dict(out.obs['donor_status'].value_counts().reindex(present))}")
    print(f"  ordered       {int(ordered_mask.sum()):,}  "
          f"unavailable {int((~ordered_mask).sum()):,}")
    print("  reasons:")
    for reason, count in out.obs["pseudotime_unavailable_reason"].value_counts().items():
        label = reason if str(reason) else "(ordered)"
        print(f"    {label:<52s} {count:>7,}")
    print(f"  hard_pass     {out.uns['trajectory_validation']['hard_pass']}")
    medians = out.uns["trajectory_validation"]["within_status_donor_eta2_of_islet_medians"]
    print("  donor eta2 of islet-median pseudotime, within status "
          "(reported, not a batch check on a within-donor fit): "
          + ", ".join(f"{k}={v:.4f}" for k, v in medians.items()))
    print(f"  ordered islets absent from this frame: "
          f"{out.uns['trajectory_validation']['islets_ordered_but_absent_from_this_frame']:,}")
    if before is not None:
        gained = sorted(set(after_ids) - set(before["imageid"]), key=lambda s: (len(s), s))
        print(f"\n  donors gained {len(gained)}: {', '.join(gained)}")
        print(f"  statuses      {before['donor_status']} -> {list(out.obs['donor_status'].cat.categories)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
