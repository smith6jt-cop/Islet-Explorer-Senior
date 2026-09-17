"""Arm B ``som-kintsugi``: KINTSUGI's pyFlowSOM clustering of the type-marker block, annotated
behind the label-free null gate, with the node count chosen by seed stability alone.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/arms/som_kintsugi.py \
        --donor 6436 --scope frame [--seed N] [--ablate MARKER] [--keep-scratch]

Method (frozen in ``scripts/broadtype/README.md``; a change here is a new arm name):

1. features = ``arm_base.within_donor_rank(ctx.X)`` -- the declared identity (type) markers
   only, each ranked in [0, 1] within the loaded cells.  Rank is the statistic
   ``celltype/README.md`` records as safe on the shape-mismatched markers, and it is what the
   kintsugi lane's arm B uses, so this arm is comparable to that lane.
2. The SOM runs in ``external/.venv-som`` (pyFlowSOM needs numpy < 2) through
   ``arms/som_kintsugi_worker.py`` -> ``scripts/kintsugi/som_adapter.som_cluster``, which
   asserts the pinned KINTSUGI checkout.  Grid: nodes in {9, 16, 25, 36, 49} x seeds
   {ctx.seed, ..., ctx.seed + 4}; num_passes / lr from ``som_adapter.NOTEBOOK_SOM``.
3. Every fit is annotated with ``gate.null_gated_cluster_annotation`` on
   ``Z = gate.standardise_within_donor(features)``: a cluster is awarded a label only where its
   declared-marker contrast beats the p95 of 200 size-matched random cell sets AND ranks first
   by a 25 % margin; otherwise Ambiguous(multi_lineage) / Unavailable(no_lineage_beats_null) /
   Unavailable(cluster_below_min_cells).
4. The node count is the argmax over the grid of the mean pairwise ARI of the ANNOTATED label
   vectors across the 5 seeds -- stability of the answer this arm gives, never agreement with
   another arm, a QuPath class, or an expected composition.  Ties go to the fewest nodes.
5. Primary calls = the fit at ``ctx.seed`` with the chosen node count.  ``score`` = the awarded
   label's contrast minus its null p95 (cluster-level; NaN for Ambiguous/Unavailable).

False-award rate, measured beside every table as the protocol requires.  One per-column
shuffled copy of the feature block is made (marginals kept, every cell-level association
destroyed) and two rates are recorded from it:

* ``annotation_only`` -- the annotation re-run on the shuffled ``Z`` with the PRIMARY clusters
  (the task's literal specification).  Cluster membership is then independent of the features,
  so this is the gate's Type-I rate on random cell sets of the clusters' sizes.
* ``arm`` -- the SOM re-fitted on the shuffled block over the same nodes x seeds grid (a second
  worker call) and annotated with the same gate; awarded / eligible clusters pooled over the 5
  seeds per node count.  This is what the whole arm does on structureless input, and the pooled
  rate at the chosen node count is the headline ``false_award_rate``.  Pooling is what makes it
  usable as a reference: one fit of 16 nodes estimates a rate to +/- 0.12.

DEVIATION, recorded here and in the report: the plan named only the first measurement.  On
synthetic iid noise it reads ~0.25 while the arm itself awards ~0.7 of its clusters, because a
SOM cluster is *selected* to be co-extreme in the features the gate then tests, and the gate's
random-cell-set null does not model that selection.  The literal measurement therefore cannot
be the ceiling the synthetic noise test asks for; the arm-level one can, so both are recorded
and the test reads the arm-level rate.  Nothing in ``gate.py`` is changed by this arm.

Missing measurements: a cell with a non-finite value in any identity feature cannot enter a
Euclidean SOM and is ``Unavailable(non_finite_identity_input)`` -- never imputed, never zero.
(On donor 6436's frame those are the 5,732 cells the geometry stage refused; observed after
the fact, not used as an input.)  A label whose declared markers are all absent from
``ctx.markers`` (ablation or panel variant) is excluded from ``declared`` and can never be
awarded; the exclusion is recorded.

Recording: ``arm_base.cli`` passes ``DEFAULT_PARAMS`` to the record by reference AFTER
``run`` returns, so ``run`` fills ``DEFAULT_PARAMS["result"]`` (chosen nodes, ARI table,
false-award rates, per-fit summaries, the primary annotation table) and the same dict is placed
in ``ctx.extra["som_kintsugi"]`` for a caller that writes its own record.  Scratch parquet
lives under ``data/development/broadtype/scratch/som-kintsugi/`` and is deleted when the run
ends unless ``--keep-scratch`` is given.
"""
from __future__ import annotations

import itertools
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))               # scripts/broadtype
import arm_base  # noqa: E402
import declarations  # noqa: E402
import gate  # noqa: E402
import pins  # noqa: E402
from criteria import CRITERIA  # noqa: E402
from records import AMBIGUOUS, CALL_COLUMNS, UNAVAILABLE  # noqa: E402

sys.path.insert(0, str(pins.SCRIPTS / "kintsugi"))
import som_adapter  # noqa: E402  (pin constants and NOTEBOOK_SOM only; pyFlowSOM never imports in .venv)

ARM = "som-kintsugi"
STOCHASTIC = True   # honours ctx.seed; the screen runs seeds 1..5 as replicates
WORKER = HERE / "som_kintsugi_worker.py"
SOM_PYTHON = pins.REPO / "external" / ".venv-som" / "bin" / "python"
SCRATCH_ROOT = pins.OUT / "scratch" / ARM
KEEP_SCRATCH = False

NODES_GRID = (9, 16, 25, 36, 49)
N_SEEDS = 5
NON_FINITE_REASON = "non_finite_identity_input"

DEFAULT_PARAMS = {
    "features": "arm_base.within_donor_rank(ctx.X): declared type markers only, per-marker rank in [0, 1] "
                "within the loaded cells; a cell with any non-finite feature is Unavailable, never imputed",
    "som": {
        "adapter": "scripts/kintsugi/som_adapter.som_cluster (KINTSUGI Kanalysis.visualize_som core)",
        "nodes_grid": list(NODES_GRID),
        "n_seeds": N_SEEDS,
        "seeds": "ctx.seed .. ctx.seed + 4 (pyFlowSOM node-initialisation seed)",
        "num_passes": som_adapter.NOTEBOOK_SOM["num_passes"],
        "lr_start": som_adapter.NOTEBOOK_SOM["lr_start"],
        "lr_end": som_adapter.NOTEBOOK_SOM["lr_end"],
        "scaling": "StandardScaler inside som_cluster -- per-marker feature weighting for Euclidean distance",
    },
    "annotation": {
        "fn": "gate.null_gated_cluster_annotation",
        "Z": "gate.standardise_within_donor(features)",
        "declared": "label -> column indices of its declared markers present in ctx.markers; a label with none is excluded",
        "null_draws": CRITERIA["null_draws"],
        "null_quantile": CRITERIA["null_quantile"],
        "margin_frac": 0.25,
        "min_cells": CRITERIA["min_calls"],
        "seed": "ctx.seed (the same null seed for every fit, so seed stability isolates the SOM)",
    },
    "node_choice": "argmax over nodes_grid of the mean pairwise adjusted Rand index (sklearn) of the annotated "
                   "per-cell label vectors across the 5 seeds, over clustered cells; ties -> fewest nodes; "
                   "never by agreement with another arm, class or composition",
    "primary": "the fit at ctx.seed with the chosen node count",
    "score": "awarded label's cluster contrast minus its null p95; NaN for Ambiguous/Unavailable",
    "false_award_measurement": {
        "shuffled_copy": "every feature column independently permuted (rng seed ctx.seed), then standardised",
        "annotation_only": "annotation on the shuffled Z with the PRIMARY clusters; awarded / clusters >= min_cells "
                           "(the gate's Type-I rate on random cell sets)",
        "arm": "SOM re-fitted on the shuffled block over the same nodes x seeds grid, annotated with the same gate; "
               "awarded / clusters >= min_cells pooled over the 5 seeds per node count "
               "(what the arm does on structureless input)",
        "headline": "false_award_rate = pooled arm rate at the chosen node count",
    },
    "pins": {
        "kintsugi_commit": som_adapter.KINTSUGI_COMMIT,
        "kintsugi_version": som_adapter.KINTSUGI_VERSION,
        "som_python": str(SOM_PYTHON),
        "worker": str(WORKER.relative_to(pins.REPO)),
    },
    "result": None,        # filled by run(): chosen nodes, ARI table, false-award rates, fit summaries
}


# --------------------------------------------------------------------------- helpers
def declared_columns(markers: list) -> tuple[dict, dict]:
    """label -> indices into ``markers`` of its declared markers; labels with none are excluded.

    The exclusion is what lets an ablation run transparently: withholding a label's only
    remaining marker leaves the arm nothing to test that label with, so it can never be awarded
    and the record says so, instead of the label silently dropping out of the table.
    """
    dec = declarations.load()
    pos = {m: j for j, m in enumerate(markers)}
    declared, excluded = {}, {}
    for L, spec in dec["labels"].items():
        cols = [pos[m] for m in spec["markers"] if m in pos]
        if cols:
            declared[L] = cols
        else:
            excluded[L] = "no declared marker present in the identity input"
    return declared, excluded


def mean_pairwise_ari(vectors: list) -> tuple[float, list]:
    """Mean ARI over all pairs of integer-coded label vectors (descriptive of stability only)."""
    from sklearn.metrics import adjusted_rand_score
    pairs = []
    for (i, a), (j, b) in itertools.combinations(enumerate(vectors), 2):
        pairs.append({"i": i, "j": j, "ari": float(adjusted_rand_score(a, b))})
    return (float(np.mean([p["ari"] for p in pairs])) if pairs else float("nan")), pairs


def choose_nodes(ari_by_nodes: dict) -> int:
    """argmax of mean pairwise ARI; ties -> fewest nodes (a stated, content-blind rule)."""
    best = max(ari_by_nodes.values())
    return min(n for n, v in ari_by_nodes.items() if v == best)


def shuffle_columns(F: np.ndarray, seed: int) -> np.ndarray:
    """Independently permute each column: marginals kept, every cell-level association destroyed."""
    rng = np.random.default_rng(seed)
    out = np.empty_like(F)
    for j in range(F.shape[1]):
        out[:, j] = F[rng.permutation(F.shape[0]), j]
    return out


def award_rate(ann: pd.DataFrame) -> dict:
    """Fraction of annotatable clusters (>= min_cells) the gate awarded a label; cells beside it."""
    eligible = ann["reason"].fillna("") != "cluster_below_min_cells"
    awarded = ~ann["label"].isin([AMBIGUOUS, UNAVAILABLE])
    n_cells = int(ann["n"].sum())
    return {"n_clusters": int(len(ann)), "n_eligible": int(eligible.sum()), "n_awarded": int(awarded.sum()),
            "n_ambiguous": int((ann["label"] == AMBIGUOUS).sum()),
            "rate": float(awarded.sum() / eligible.sum()) if eligible.any() else float("nan"),
            "frac_cells_awarded": float(ann.loc[awarded, "n"].sum() / n_cells) if n_cells else float("nan")}


def run_worker(features: pd.DataFrame, ids: np.ndarray, fits: list, scratch: pathlib.Path) -> dict:
    """Write the block, fit every (nodes, seed) in .venv-som, return {(nodes, seed): clusters}."""
    if not SOM_PYTHON.exists():
        raise RuntimeError(
            f"{SOM_PYTHON} is missing. The SOM needs numpy < 2 (pyFlowSOM); build it as "
            "scripts/kintsugi/README.md says:\n  python3 -m venv external/.venv-som\n"
            "  external/.venv-som/bin/pip install -c external/KINTSUGI/constraints.txt "
            "pyFlowSOM pandas scipy scikit-learn pyarrow shapely")
    scratch.mkdir(parents=True, exist_ok=True)
    feat_path = scratch / "features.parquet"
    block = features.copy()
    block.insert(0, "object_id", pd.Series(ids).astype(str).to_numpy())
    block.to_parquet(feat_path, index=False, row_group_size=100_000)
    fit_dir = scratch / "fits"
    cmd = [str(SOM_PYTHON), str(WORKER), "--features", str(feat_path), "--out-dir", str(fit_dir),
           "--num-passes", str(som_adapter.NOTEBOOK_SOM["num_passes"]),
           "--lr-start", str(som_adapter.NOTEBOOK_SOM["lr_start"]),
           "--lr-end", str(som_adapter.NOTEBOOK_SOM["lr_end"]),
           "--fits", *[f"{n}:{s}" for n, s in fits]]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(pins.PKG))
    if proc.returncode != 0:
        raise RuntimeError(f"SOM worker failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
    manifest = json.loads((fit_dir / "worker_manifest.json").read_text())
    want = pd.Index(pd.Series(ids).astype(str))
    clusters = {}
    for rec in manifest["fits"]:
        tab = pd.read_parquet(rec["path"])
        got = pd.Index(tab["object_id"].astype(str))
        if not want.equals(got):
            raise RuntimeError(f"fit nodes={rec['nodes']} seed={rec['seed']}: object_id order does not match the input")
        clusters[(int(rec["nodes"]), int(rec["seed"]))] = tab["cluster"].to_numpy(np.int32)
    return {"clusters": clusters, "manifest": manifest}


def annotate_fit(Z: np.ndarray, clusters: np.ndarray, declared: dict, seed: int) -> tuple[pd.DataFrame, dict]:
    """Annotate one fit; return (the per-cluster table, cluster -> its row as a dict)."""
    ann = gate.null_gated_cluster_annotation(Z, clusters, declared, seed=seed)
    by_cluster = {int(r["cluster"]): r for r in ann.to_dict("records")}
    return ann, by_cluster


def label_vector(clusters: np.ndarray, by_cluster: dict, vocab: dict) -> np.ndarray:
    """Per-cell integer label codes for one fit (Ambiguous/Unavailable are labels here too)."""
    lut = np.empty(int(clusters.max()) + 1, dtype=np.int32)
    for c, r in by_cluster.items():
        lut[c] = vocab[r["label"]]
    return lut[clusters]


# --------------------------------------------------------------------------- the arm
def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    t0 = time.time()
    dec = declarations.load()
    ids = np.asarray(ctx.ids).astype(str)
    n = len(ids)
    result = {"donor_id": ctx.donor_id, "scope": ctx.scope, "config_id": ctx.config_id, "seed": int(ctx.seed),
              "n_cells": int(n), "identity_markers": list(ctx.markers)}
    DEFAULT_PARAMS["result"] = result
    ctx.extra["som_kintsugi"] = result

    declared, excluded_labels = declared_columns(list(ctx.markers))
    result["declared_columns"] = {L: [ctx.markers[j] for j in cols] for L, cols in declared.items()}
    result["excluded_labels"] = excluded_labels
    if not declared:
        result["outcome"] = "no declared label has a marker in the identity input"
        return pd.DataFrame({"object_id": ids, "label": UNAVAILABLE, "score": np.nan,
                             "reason": "ablation_removes_required_input"})[list(CALL_COLUMNS)]

    # 1. features: within-donor rank of the type markers; non-finite cells cannot enter a Euclidean SOM
    feats = arm_base.within_donor_rank(ctx.X)
    F = feats.to_numpy(np.float64)
    finite = np.isfinite(F).all(axis=1)
    result["n_non_finite"] = int((~finite).sum())
    result["non_finite_by_marker"] = {m: int(k) for m, k in zip(ctx.markers, (~np.isfinite(F)).sum(axis=0)) if k}
    calls = pd.DataFrame({"object_id": ids, "label": UNAVAILABLE, "score": np.nan, "reason": NON_FINITE_REASON})
    n_fit = int(finite.sum())
    if n_fit < CRITERIA["min_calls"]:
        result["outcome"] = f"only {n_fit} cells with finite identity features"
        calls.loc[finite, "reason"] = "too_few_finite_cells_to_cluster"
        return calls[list(CALL_COLUMNS)]
    Ff, idf = F[finite], ids[finite]
    Fs = shuffle_columns(Ff, int(ctx.seed))           # the structureless twin, for the false-award rates

    # 2. the SOM grid (and the shuffled twin's grid at ctx.seed), in .venv-som
    seeds = [int(ctx.seed) + k for k in range(N_SEEDS)]
    fits = [(nn, s) for nn in NODES_GRID for s in seeds]
    scratch = SCRATCH_ROOT / f"{ctx.donor_id}__{ctx.scope}__{ctx.config_id.replace(':', '_')}__{os.getpid()}"
    t_som = time.time()
    try:
        cols = list(ctx.markers)
        worker = run_worker(pd.DataFrame(Ff, columns=cols), idf, fits, scratch / "real")
        shuffled = run_worker(pd.DataFrame(Fs, columns=cols), idf, fits, scratch / "shuffled")
    finally:
        if not KEEP_SCRATCH:
            shutil.rmtree(scratch, ignore_errors=True)
    result["som_seconds"] = round(time.time() - t_som, 1)
    result["worker"] = {k: v for k, v in worker["manifest"].items() if k != "fits"}
    result["fits"] = worker["manifest"]["fits"]

    # 3. annotate every fit against the same standardised block and null seed
    Z = gate.standardise_within_donor(Ff)
    vocab = {L: i for i, L in enumerate([*dec["labels"], AMBIGUOUS, UNAVAILABLE])}
    inv_vocab = {v: k for k, v in vocab.items()}
    t_ann = time.time()
    annotated, vectors = {}, {}
    for key in fits:
        ann, by_cluster = annotate_fit(Z, worker["clusters"][key], declared, int(ctx.seed))
        annotated[key] = (ann, by_cluster)
        vectors[key] = label_vector(worker["clusters"][key], by_cluster, vocab)
    for rec in result["fits"]:
        key = (rec["nodes"], rec["seed"])
        ann = annotated[key][0]
        rec["award"] = award_rate(ann)
        rec["cluster_labels"] = ann["label"].value_counts().to_dict()
        rec["cell_labels"] = pd.Series(vectors[key]).map(inv_vocab).value_counts().to_dict()
    result["annotation_seconds"] = round(time.time() - t_ann, 1)

    # 4. node count by seed stability of the annotated labels (and nothing else)
    ari_table = {}
    for nn in NODES_GRID:
        mean_ari, pairs = mean_pairwise_ari([vectors[(nn, s)] for s in seeds])
        ari_table[nn] = {"mean_pairwise_ari": mean_ari, "pairs": pairs}
    chosen = choose_nodes({nn: v["mean_pairwise_ari"] for nn, v in ari_table.items()})
    result["ari_by_nodes"] = {str(nn): v for nn, v in ari_table.items()}
    result["chosen_nodes"] = int(chosen)

    # 5. primary calls from the fit at ctx.seed with the chosen node count
    key = (chosen, int(ctx.seed))
    ann, by_cluster = annotated[key]
    clusters = worker["clusters"][key]
    lab = np.empty(n_fit, dtype=object)
    score = np.full(n_fit, np.nan)
    reason = np.empty(n_fit, dtype=object)
    for c, r in by_cluster.items():
        sel = clusters == c
        lab[sel] = r["label"]
        reason[sel] = r.get("reason", "") or ""
        if r["label"] not in (AMBIGUOUS, UNAVAILABLE):
            score[sel] = float(r[f"contrast_{r['label']}"]) - float(r[f"p95_{r['label']}"])
    calls.loc[finite, "label"] = lab
    calls.loc[finite, "score"] = score
    calls.loc[finite, "reason"] = reason
    result["primary_annotation"] = json.loads(ann.to_json(orient="records"))
    result["primary_award"] = award_rate(ann)

    # 6. the false-award rates on the structureless twin, beside the table
    Zs = gate.standardise_within_donor(Fs)
    ann_only = gate.null_gated_cluster_annotation(Zs, clusters, declared, seed=int(ctx.seed))
    arm_rates = {}
    for nn in NODES_GRID:
        per_seed = [award_rate(gate.null_gated_cluster_annotation(Zs, shuffled["clusters"][(nn, s)], declared, seed=int(ctx.seed)))
                    for s in seeds]
        n_el, n_aw = sum(r["n_eligible"] for r in per_seed), sum(r["n_awarded"] for r in per_seed)
        arm_rates[str(nn)] = {"n_eligible": n_el, "n_awarded": n_aw, "rate": float(n_aw / n_el) if n_el else float("nan"),
                              "frac_cells_awarded": float(np.mean([r["frac_cells_awarded"] for r in per_seed])),
                              "per_seed": [{"seed": s, **r} for s, r in zip(seeds, per_seed)]}
    result["false_award"] = {
        "annotation_only": award_rate(ann_only),
        "arm_by_nodes": arm_rates,
        "false_award_rate": arm_rates[str(chosen)]["rate"],
        "definition": DEFAULT_PARAMS["false_award_measurement"],
    }
    result["label_counts"] = calls["label"].value_counts().to_dict()
    result["seconds"] = round(time.time() - t0, 1)
    return calls[list(CALL_COLUMNS)]


def _pop_flag(flag: str) -> bool:
    """Consume a flag arm_base.cli does not know, so the shared CLI stays the shared CLI."""
    if flag in sys.argv:
        sys.argv = [a for a in sys.argv if a != flag]
        return True
    return False


if __name__ == "__main__":
    KEEP_SCRATCH = _pop_flag("--keep-scratch")
    raise SystemExit(arm_base.cli(ARM, run, env=str(SOM_PYTHON), default_params=DEFAULT_PARAMS))
