#!/home/smith6jt/miniconda3/envs/scvi-env/bin/python
"""Arm X2 ``joint-run-aware``: one scVI over a hash subsample of every donor, batch = the Xenium
instrument run, Leiden clusters on the joint latent, lineages awarded to clusters by the null gate,
labels carried to every cell by a kNN vote in the same latent.

    cd core/multiplex5000
    /home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/xenium_typing/joint_model.py --dry-run
    /home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/xenium_typing/joint_model.py --n-per-donor 50000

Why joint and why the run is the batch.  A per-donor model gives twenty latent spaces whose clusters
correspond only by name; a joint model gives one space in which a cluster IS the same population
across donors, which is what a cross-donor comparison needs.  ``batch_key = xenium_run`` (the 11
instrument runs, the technical unit) and never the donor: six of the eleven runs mix ND with non-ND
donors, so correcting on run does not alias disease, while correcting on donor would remove the very
biology the study is about (decision #13's reasoning, applied at cohort scale).

Identity vs state (decision #17): the embedding uses ``var['identity_panel']`` (478 curated identity
genes; zero overlap with ``is_state_only``, asserted per base); state panels are scored AFTER typing
(``state_panels.score_state_panels``) and written as ``state_*_score`` columns; ``cycling`` is the
existing detection rule (>= 2 of MKI67/TOP2A/PCNA at count >= 1).  No gene is selected by variance
(decision #2), no library-size covariate is passed (decision #1), ``random_state`` is fixed
(decision #6), ``dl_num_workers = 0`` (decision #16).

Annotation of a cluster: ``gate.nullgate_table`` on the subsample's lognorm -- a lineage is awarded
only when its identity core beats a size- and detection-matched random core and it is the ONLY
lineage that does (two -> Ambiguous(multi_lineage), none -> Unavailable).  Resolution is chosen by seed
stability (mean pairwise ARI of the ANNOTATED labels across 5 seeds), never by agreement with another
arm.  Per-cell label = majority of the 15 nearest training cells' cluster labels in the joint latent;
vote fraction < 0.5 -> Ambiguous(low_margin).  Replicates = the same at the other four seeds.  Fine
tier = the identity-core ``subtype_hint`` panels of the cell's broad lineage (detection-gated argmax,
the typer's tier-2 rule in its simplest form; recorded as such).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gate  # noqa: E402
import pins  # noqa: E402
sys.path.insert(0, str(pins.XENIUM_ANALYSIS / "scripts"))
sys.path.insert(0, str(pins.XENIUM_ANALYSIS))

ARM = pins.ARM_JOINT
DEFAULT = {
    "method_version": "xenium-typing-joint-run-aware-2026-09-05",
    "n_per_donor": 50_000, "batch_key": "xenium_run", "genes": "var['identity_panel'] (478)",
    "n_latent": 20, "n_layers": 2, "gene_likelihood": "nb", "max_epochs": 400, "early_stopping_patience": 40,
    "batch_size": 4096, "dl_num_workers": 0, "knn_k": 15, "resolutions": [0.3, 0.5, 1.0], "seeds": [0, 1, 2, 3, 4],
    "vote_margin": 0.5, "min_cluster_cells": 50, "fine_tier": "identity-core subtype_hint panels, detection-gated argmax within the broad lineage",
    "cycling_rule": ">= 2 of MKI67/TOP2A/PCNA at raw count >= 1",
    "leiden": "igraph Graph.community_leiden(objective_function='modularity', resolution=r, n_iterations=2, beta=0.01), seed via igraph RNG; "
              "replaces leidenalg n_iterations=-1, which had not completed one fit on 1,000,000 nodes after 3.5 h (stopped 2026-09-06 before any result)",
}
OUT = pins.OUT / "joint_model"


def log(msg: str) -> None:
    print(f"[joint {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_base(donor: str, *, genes_mask_col: str = "identity_panel"):
    import anndata as ad
    a = ad.read_h5ad(pins.base_path(donor))
    if "identity_panel" in a.var.columns and "is_state_only" in a.var.columns:
        overlap = int((a.var["identity_panel"].to_numpy(bool) & a.var["is_state_only"].to_numpy(bool)).sum())
        if overlap:
            raise RuntimeError(f"donor {donor}: {overlap} state-only genes inside identity_panel -- tiering violated")
    a.obs["donor_id"] = str(donor)
    a.obs["xenium_run"] = pins.xenium_run(donor)
    return a


def subsample(a, n: int, salt: str):
    idx = gate.hash_subsample(np.asarray(a.obs_names).astype(str), n, salt=salt)
    return a[np.sort(idx)].copy()


def train_joint(train, *, seed: int, p: dict, use_gpu: bool):
    import scvi
    import torch
    scvi.settings.seed = seed
    scvi.settings.dl_num_workers = 0
    if use_gpu:
        torch.set_float32_matmul_precision("high")
    scvi.model.SCVI.setup_anndata(train, layer="counts", batch_key=p["batch_key"])
    model = scvi.model.SCVI(train, n_latent=p["n_latent"], n_layers=p["n_layers"], gene_likelihood=p["gene_likelihood"])
    t0 = time.time()
    model.train(max_epochs=p["max_epochs"], early_stopping=True, early_stopping_patience=p["early_stopping_patience"],
                batch_size=p["batch_size"], accelerator="gpu" if use_gpu else "cpu", devices=1,
                load_sparse_tensor=bool(use_gpu))
    log(f"scVI trained in {(time.time() - t0) / 60:.1f} min; epochs run {len(model.history['elbo_train'])}")
    return model


LEIDEN_ITERATIONS = 2


def leiden_fits(Z: np.ndarray, *, k: int, resolutions: list, seeds: list) -> tuple:
    """kNN graph (pynndescent) + Leiden through igraph's C implementation (the Xenium lane's decision #5
    flavour).  The first full-size run used ``leidenalg.find_partition(..., n_iterations=-1)`` and had not
    finished ONE fit after 3.5 h on 1,000,000 nodes; it was stopped before any result existed and this
    implementation recorded in its place (modularity objective with a resolution parameter is the
    RBConfiguration partition; ``n_iterations=2`` is scanpy's igraph-flavour default)."""
    import random
    import igraph as ig
    import pynndescent
    t0 = time.time()
    index = pynndescent.NNDescent(Z.astype(np.float32), n_neighbors=k + 1, metric="euclidean", random_state=0, n_jobs=8)
    nbr, _ = index.neighbor_graph
    n = Z.shape[0]
    src = np.repeat(np.arange(n), k)
    dst = nbr[:, 1:k + 1].ravel()
    g = ig.Graph(n=n, edges=np.column_stack([src, dst]).astype(np.int64).tolist(), directed=False).simplify()
    log(f"kNN graph: {n:,} nodes, {g.ecount():,} edges ({time.time() - t0:.0f}s)")
    fits = {}
    for r in resolutions:
        for s in seeds:
            t1 = time.time()
            ig.set_random_number_generator(random.Random(int(s)))
            part = g.community_leiden(objective_function="modularity", resolution=float(r), n_iterations=LEIDEN_ITERATIONS, beta=0.01)
            fits[(r, s)] = np.asarray(part.membership, dtype=np.int64)
            log(f"leiden r={r} seed={s}: {len(set(part.membership))} clusters ({time.time() - t1:.0f}s)")
    return fits, index


def annotate(X_lognorm, genes, clusters: np.ndarray, cores: dict, *, seed: int) -> tuple:
    long, _ = gate.nullgate_table(X_lognorm, genes, clusters.astype(str), cores, seed=seed)
    lab = {}
    for g_, sub in long.groupby("group"):
        if not sub["judgeable"].any():
            lab[g_] = (gate.UNAVAILABLE, "cluster_below_min_cells")
            continue
        clear = sub.loc[sub["clears"], "lineage"].tolist()
        lab[g_] = (clear[0], "") if len(clear) == 1 else ((gate.UNAVAILABLE, "no_lineage_beats_null") if not clear else (gate.AMBIGUOUS, "multi_lineage"))
    labels = np.array([lab[str(c)][0] for c in clusters], dtype=object)
    reasons = np.array([lab[str(c)][1] for c in clusters], dtype=object)
    return labels, reasons, long


def mean_pairwise_ari(vectors: list) -> float:
    from sklearn.metrics import adjusted_rand_score
    vals = [adjusted_rand_score(vectors[i], vectors[j]) for i in range(len(vectors)) for j in range(i + 1, len(vectors))]
    return float(np.mean(vals)) if vals else float("nan")


def fine_tier(a, broad: np.ndarray, cores_fine: dict, fine_to_broad: dict, *, floor: int = 1) -> np.ndarray:
    """Detection-gated argmax over the fine panels of the cell's broad lineage."""
    genes = {g: i for i, g in enumerate(np.asarray(a.var_names).astype(str))}
    L = a.layers["lognorm"]
    C = a.layers["counts"]
    fine = broad.astype(object).copy()
    for B in set(fine_to_broad.values()):
        fts = [f for f, b in fine_to_broad.items() if b == B and f in cores_fine and any(g in genes for g in cores_fine[f])]
        rows = np.flatnonzero(broad == B)
        if not rows.size or not fts:
            continue
        S = np.full((rows.size, len(fts)), -np.inf)
        for j, f in enumerate(fts):
            cols = [genes[g] for g in cores_fine[f] if g in genes]
            score = np.asarray(L[rows][:, cols].mean(axis=1)).ravel()
            det = np.asarray((C[rows][:, cols] >= floor).sum(axis=1)).ravel()
            S[:, j] = np.where(det >= 1, score, -np.inf)
        best = S.argmax(axis=1)
        ok = np.isfinite(S[np.arange(rows.size), best])
        fine[rows[ok]] = np.asarray(fts, dtype=object)[best[ok]]
    return fine


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=None)
    ap.add_argument("--n-per-donor", type=int, default=DEFAULT["n_per_donor"])
    ap.add_argument("--max-epochs", type=int, default=DEFAULT["max_epochs"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out-suffix", default="")
    ap.add_argument("--reuse-model", action="store_true",
                    help="load <out>/scvi_model (trained earlier on the identical hash-drawn training set) instead of training again")
    a = ap.parse_args(argv)
    p = dict(DEFAULT, n_per_donor=a.n_per_donor, max_epochs=a.max_epochs)
    donors = list(a.donors) if a.donors else list(pins.donors())
    snap = pins.panel_snapshot()
    out = pathlib.Path(str(OUT) + a.out_suffix)
    out.mkdir(parents=True, exist_ok=True)
    import anndata as ad
    import torch
    use_gpu = torch.cuda.is_available()
    log(f"{len(donors)} donors; n_per_donor {p['n_per_donor']}; gpu={use_gpu}; out={out}")

    # A. subsample every donor
    subs, runs = [], {}
    for d in donors:
        base = read_base(d)
        runs[d] = pins.xenium_run(d)
        sub = subsample(base, p["n_per_donor"], salt=f"joint|{a.seed}")
        subs.append(sub)
        log(f"donor {d} ({runs[d]}): {base.n_obs:,} cells -> subsample {sub.n_obs:,}")
        del base
    if a.dry_run:
        log("dry run: inputs load; stopping")
        return 0
    train_full = ad.concat(subs, join="inner", merge="same")
    train_full.obs_names_make_unique()
    del subs
    ident = train_full.var["identity_panel"].to_numpy(bool)
    ident_genes = [str(g_) for g_ in np.asarray(train_full.var_names)[ident]]
    train = train_full[:, ident].copy()
    log(f"training set {train.n_obs:,} cells x {train.n_vars} identity genes; runs {sorted(set(train.obs['xenium_run']))}")

    # B/C. scVI + latent
    if a.reuse_model and (out / "scvi_model").exists():
        import scvi
        scvi.settings.seed = a.seed
        if use_gpu:
            torch.set_float32_matmul_precision("high")      # same matmul precision as the training run
        scvi.model.SCVI.setup_anndata(train, layer="counts", batch_key=p["batch_key"])
        model = scvi.model.SCVI.load(str(out / "scvi_model"), adata=train, accelerator="gpu" if use_gpu else "cpu")
        log(f"reused {out / 'scvi_model'} (trained earlier on the same hash-drawn training set; no training this run)")
    else:
        model = train_joint(train, seed=a.seed, p=p, use_gpu=use_gpu)
        model.save(str(out / "scvi_model"), overwrite=True)
    Z = model.get_latent_representation()
    np.save(out / "train_latent.npy", Z)
    train_full.obs[["donor_id", "xenium_run"]].to_parquet(out / "train_obs.parquet")

    # D/E. clusters on the joint latent, annotated by the null gate, resolution by seed stability
    cores = pins.identity_cores()
    fits, index = leiden_fits(Z, k=p["knn_k"], resolutions=p["resolutions"], seeds=p["seeds"])
    ann = {}
    for key, cl in fits.items():
        labels, reasons, long = annotate(train_full.layers["lognorm"], np.asarray(train_full.var_names).astype(str), cl, cores, seed=a.seed)
        ann[key] = (labels, reasons)
        long.to_parquet(out / f"nullgate__r{key[0]}_s{key[1]}.parquet", index=False)
        log(f"fit r={key[0]} s={key[1]}: {pd.Series(labels).value_counts().to_dict()}")
    stability = {r: mean_pairwise_ari([ann[(r, s)][0] for s in p["seeds"]]) for r in p["resolutions"]}
    chosen = max(stability, key=stability.get)
    log(f"seed stability (mean pairwise ARI of annotated labels): {stability}; chosen resolution {chosen}")
    (out / "fits.json").write_text(json.dumps({"stability": stability, "chosen_resolution": chosen, "params": p,
                                              "panel": snap, "n_train": int(train.n_obs), "runs": runs}, indent=1, default=str))

    # F-J. every cell of every donor
    from tissue_markers import pancreas as mod
    from tissue_markers._identity_core import hint_fine_panels
    import state_panels
    fine_to_broad = {f: B for B, fs in mod.BROAD_TO_FINE.items() for f in fs}
    for d in donors:
        t0 = time.time()
        base = read_base(d)
        cores_fine = hint_fine_panels("pancreas", mod, reference_dir=pins.SUBPANEL_DIR, live_genes=list(base.var_names))
        # subset by NAME: donor bases differ by a few genes (5,098 vs 5,101), so a boolean mask from the
        # inner-joined training var index does not line up with every base
        missing = [g_ for g_ in ident_genes if g_ not in set(base.var_names)]
        if missing:
            raise RuntimeError(f"donor {d}: {len(missing)} training identity genes absent from its base: {missing[:5]}")
        sub_d = base[:, ident_genes].copy()
        import scvi
        scvi.model.SCVI.setup_anndata(sub_d, layer="counts", batch_key=p["batch_key"])
        Zd = model.get_latent_representation(sub_d)
        nbr, _ = index.query(Zd.astype(np.float32), k=p["knn_k"])
        primary_labels, prim_reasons = ann[(chosen, a.seed)]
        def vote(lab_train):
            votes = lab_train[nbr]                                    # (n, k) object
            out_l = np.empty(len(votes), dtype=object); conf = np.empty(len(votes)); rs = np.empty(len(votes), dtype=object)
            for i in range(len(votes)):
                u, c = np.unique(votes[i], return_counts=True)
                j = c.argmax(); frac = c[j] / c.sum()
                if frac < p["vote_margin"]:
                    out_l[i], conf[i], rs[i] = gate.AMBIGUOUS, frac, "low_margin"
                else:
                    out_l[i], conf[i], rs[i] = u[j], frac, ""
            return out_l, conf, rs
        lab, conf, rs = vote(primary_labels)
        # a cluster the gate left Ambiguous/Unavailable carries its reason to the cells it labels
        reason_of = dict(zip(primary_labels, prim_reasons))
        rs = np.where((rs == "") & np.isin(lab, [gate.AMBIGUOUS, gate.UNAVAILABLE]), [reason_of.get(x, "") for x in lab], rs)
        fine = fine_tier(base, lab, cores_fine, fine_to_broad)
        genes = {g: i for i, g in enumerate(np.asarray(base.var_names).astype(str))}
        cyc_cols = [genes[g] for g in mod.CYCLING_MARKERS if g in genes]
        cycling = np.asarray((base.layers["counts"][:, cyc_cols] >= 1).sum(axis=1)).ravel() >= 2
        state_panels.score_state_panels(base, "pancreas", reference_dir=pins.SUBPANEL_DIR)
        calls = pd.DataFrame({"cell_id": np.asarray(base.obs_names).astype(str), "lineage_native": lab, "celltype_fine": fine,
                              "confidence": conf, "reason": rs, "cycling": cycling})
        for c in base.obs.columns:
            if c.startswith("state_") and c.endswith("_score"):
                calls[c] = base.obs[c].astype(float).to_numpy()
        reps = {}
        for s in p["seeds"]:
            if s == a.seed:
                continue
            rl, _c, _r = vote(ann[(chosen, s)][0])
            reps[f"seed:{s}"] = rl
        rep_df = pd.DataFrame({"cell_id": calls["cell_id"], **reps})
        rec = gate.ArmRecord(arm=ARM, donor_id=str(d), config_id="primary", identity_inputs={L: len(g) for L, g in cores.items()},
                             n_cells=int(len(calls)), n_called=int(calls["lineage_native"].isin(pins.LINEAGES).sum()),
                             n_ambiguous=int((calls["lineage_native"] == gate.AMBIGUOUS).sum()),
                             n_unavailable=int((calls["lineage_native"] == gate.UNAVAILABLE).sum()),
                             seconds=time.time() - t0, env=sys.executable, stochastic=True,
                             pinned={"base": str(pins.base_path(d)), "panel": snap, "model": str(out / "scvi_model"), "fits": str(out / "fits.json")},
                             params={**p, "chosen_resolution": chosen, "stability": stability, "batch_of_donor": runs[d]},
                             notes="joint scVI (batch = xenium_run) + Leiden + null-gated cluster annotation; labels by kNN vote in the joint latent")
        gate.write_calls(calls, rec, expected_ids=calls["cell_id"].to_numpy(), replicates=rep_df)
        log(f"donor {d}: {calls['lineage_native'].value_counts().to_dict()} ({(time.time() - t0) / 60:.1f} min)")
        del base, sub_d
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
