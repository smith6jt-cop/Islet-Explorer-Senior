"""The label-free gate for Xenium typing arms, and the one record shape they all write.

Every statistic here is a measurement against a null drawn from the same cells; none consults
an expected frequency, a disease label, another arm's answer, or any QuPath object.

Blocking, per arm and donor, per lineage label:
  (a) matched-null lineage award on the LABELS -- the arm's label vector is the grouping and
      ``nullgate`` (``scripts/xenium_celltype``) asks whether the label's OWN identity core
      reads higher inside the label than outside, against 200 size- and detection-matched
      random gene sets.  A label is *awarded* iff its own core beats ``max(p95, 0)``.  The
      statistic is exactly ``hpg/nullgate_fast.annotate_groups`` (proved equal to the
      reference to 1e-18), re-run here in a form that keeps every (group, lineage) row.
  (b) spatial coherence -- Moran's I of the one-hot label on a k = 12 kNN graph over the
      section's micron centroids, against 200 within-section label permutations.
  (c) stability -- the spatial verdict is invariant across the arm's replicates (seeds x
      resolutions, or the Analyzer's k ladder).  An arm that declares itself stochastic but
      supplies no replicates is UNJUDGED on stability, and unjudgeable is never a pass.
Descriptive, never blocking:
  (d) cross-arm agreement (ARI / AMI) on cells both arms labelled;
  (e) rank association of the same-section IF INS and GCG per-cell values (restain lane) with
      Endocrine membership -- BOTH hormones, AUROC and Spearman, no cell called from protein.
  (f) ``judge`` combines them: verified = awarded AND coherent AND stable, per lineage, and
      prints the family-wise false-award caveat beside every table (``FAMILY_WISE_CAVEAT``).

The bar (``BAR``) names evidence floors on a statistic's n and null draws -- never a share of
cells, a prior, or an ordering.
"""
from __future__ import annotations

import json
import pathlib
import sys
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
import gates  # noqa: E402  (scripts/common/gates.py -- fail-closed evidence floors)
import usability  # noqa: E402  (scripts/celltype/usability.py -- kNN graph + Moran's I)
from nullgate import (MIN_GROUP_CELLS, GroupCall, NullGateError, build_null_pools,  # noqa: E402
                      decide_group, detection_deciles, match_strata, resolve_cores,
                      standardise_genes)
from nullgate_fast import group_gene_deltas  # noqa: E402

AMBIGUOUS = "Ambiguous"
UNAVAILABLE = "Unavailable"
LINEAGES = pins.LINEAGES
CALL_COLUMNS = ("cell_id", "lineage_native", "celltype_fine", "confidence", "reason", "cycling")
STATE_SUFFIX = "_score"
STATE_PREFIX = "state_"

#: The bar.  Floors are on a statistic's n; nothing here is a share of cells.
BAR = {
    "min_calls": 50,                 # fewest cells a label needs before any statistic is judged
    "min_group_cells": MIN_GROUP_CELLS,
    "null_draws": 200,
    "null_quantile": 0.95,
    "morans_k": 12,
    "spatial_max_cells": 200_000,    # hash subsample for the kNN graph; the same cells for every label
    "knn_vote_margin": 0.5,          # a kNN majority below this fraction is Ambiguous(low_margin)
}

#: Printed beside every table, verbatim from the nullgate module's own measurement.
FAMILY_WISE_CAVEAT = (
    "Family-wise false-award caveat: 'real > p95' is a 5% false-positive rate per core, applied "
    "once per core offered. On pure synthetic sparse noise (random groups, no structure; 50 draws, "
    "12 groups, 8 seeds) the gate awarded a lineage to 6.2% of groups with 1 core offered and 30.2% "
    "with the 8 the pancreas taxonomy offers (nullgate.py docstring). An award on one group in one "
    "run is at that rate; the evidence is a lineage awarded consistently across resolutions, seeds, "
    "arms and donors. Measured on noise with no structure -- a floor on the false-positive rate, and "
    "no statement about power on tissue."
)


# =============================================================================== the record
@dataclass
class ArmRecord:
    """Provenance beside every ``calls.parquet``.  ``identity_inputs`` is what the tiering test
    reads: the identity-core gene sets (by lineage, with counts) that entered the label."""
    arm: str
    donor_id: str
    config_id: str                       # 'primary' | 'dry-run' | a replicate name
    identity_inputs: dict                # lineage -> n identity-core genes used (type genes only)
    n_cells: int
    n_called: int
    n_ambiguous: int
    n_unavailable: int
    seconds: float
    env: str                             # interpreter path
    stochastic: bool                     # does the method depend on a seed?
    replicates: list = field(default_factory=list)   # replicate names present in replicate_labels.parquet
    pinned: dict = field(default_factory=dict)
    params: dict = field(default_factory=dict)
    notes: str = ""
    code: dict = field(default_factory=pins.git_head)
    written_utc: str = field(default_factory=pins.utc_now)


def state_columns(df: pd.DataFrame) -> list:
    return [c for c in df.columns if c.startswith(STATE_PREFIX) and c.endswith(STATE_SUFFIX)]


def validate_calls(df: pd.DataFrame, expected_ids: np.ndarray | None = None) -> pd.DataFrame:
    missing = [c for c in CALL_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"calls missing columns {missing}")
    if df["cell_id"].duplicated().any():
        raise ValueError("duplicate cell_id in calls")
    lab = df["lineage_native"].astype(str)
    unknown = set(lab) - set(LINEAGES) - {AMBIGUOUS, UNAVAILABLE}
    if unknown:
        raise ValueError(f"undeclared labels in calls: {sorted(unknown)}")
    bad = df.loc[lab.isin([AMBIGUOUS, UNAVAILABLE]) & (df["reason"].fillna("").astype(str) == "")]
    if len(bad):
        raise ValueError(f"{len(bad)} Ambiguous/Unavailable calls carry no reason")
    if expected_ids is not None:
        want = pd.Index(pd.Series(expected_ids).astype(str))
        got = pd.Index(df["cell_id"].astype(str))
        if not want.equals(got):
            lost, extra = want.difference(got), got.difference(want)
            raise ValueError(f"calls do not cover the base cell set in order: {len(lost)} missing, "
                             f"{len(extra)} extra, order_equal={set(want) == set(got)}")
    return df


def write_calls(df: pd.DataFrame, record: ArmRecord, *, expected_ids: np.ndarray | None = None,
                replicates: pd.DataFrame | None = None) -> pathlib.Path:
    """``calls.parquet`` + a JSONL record; optional ``replicate_labels.parquet`` (cell_id + one
    label column per replicate) that ``judge`` reads for the stability criterion."""
    d = pins.arm_dir(record.arm, record.donor_id)
    d.mkdir(parents=True, exist_ok=True)
    cols = list(CALL_COLUMNS) + state_columns(df)
    df = validate_calls(df[cols].copy(), expected_ids)
    df["cell_id"] = df["cell_id"].astype(str)
    df["confidence"] = df["confidence"].astype("float32")
    df["reason"] = df["reason"].fillna("").astype(str)
    df["cycling"] = df["cycling"].astype("boolean")
    path = d / ("calls.parquet" if record.config_id == "primary" else f"calls__{record.config_id}.parquet")
    df.to_parquet(path, index=False)
    if replicates is not None:
        if not replicates["cell_id"].astype(str).equals(df["cell_id"]):
            raise ValueError("replicate_labels must list the same cells in the same order as calls")
        replicates.to_parquet(d / "replicate_labels.parquet", index=False)
        record.replicates = [c for c in replicates.columns if c != "cell_id"]
    with (d / "records.jsonl").open("a") as fh:
        fh.write(json.dumps({**asdict(record), "calls": str(path)}, default=str) + "\n")
    return path


def read_calls(arm: str, donor: str) -> pd.DataFrame:
    p = pins.arm_dir(arm, donor) / "calls.parquet"
    if not p.exists():
        raise FileNotFoundError(f"{arm} has no primary calls for donor {donor}: {p}")
    return pd.read_parquet(p)


def read_record(arm: str, donor: str, config_id: str = "primary") -> dict:
    p = pins.arm_dir(arm, donor) / "records.jsonl"
    if not p.exists():
        raise FileNotFoundError(f"{arm} donor {donor}: no records.jsonl")
    recs = [json.loads(line) for line in p.read_text().splitlines() if line.strip()]
    recs = [r for r in recs if r.get("config_id") == config_id]
    if not recs:
        raise FileNotFoundError(f"{arm} donor {donor}: no {config_id!r} record")
    return recs[-1]


def read_replicates(arm: str, donor: str) -> pd.DataFrame | None:
    p = pins.arm_dir(arm, donor) / "replicate_labels.parquet"
    return pd.read_parquet(p) if p.exists() else None


# =============================================================================== helpers
def hash_subsample(ids: np.ndarray, n: int, salt: str = "xenium_typing") -> np.ndarray:
    """Deterministic, content-keyed subsample (sorted indices), independent of row order."""
    import hashlib
    ids = np.asarray(ids).astype(str)
    if len(ids) <= n:
        return np.arange(len(ids))
    keys = np.fromiter((int(hashlib.sha256(f"{i}|{salt}".encode()).hexdigest()[:16], 16) for i in ids),
                       dtype=np.uint64, count=len(ids))
    return np.sort(np.argsort(keys, kind="stable")[:n])


def groupcall_labels(calls: pd.DataFrame) -> dict:
    """GroupCall rows -> {group: (label, reason)} for a cluster-level arm.

    An awarded lineage is the label (``multiple_cores_cleared`` stays as a note on an award);
    cores that cleared but are not separable -> Ambiguous; no core beating its null, or too few
    cells -> Unavailable.  Never the nearest lineage.
    """
    out = {}
    for r in calls.itertuples(index=False):
        if r.lineage:
            out[r.group] = (str(r.lineage), str(r.reason))
        elif r.reason == "cores_cleared_but_are_not_separable":
            out[r.group] = (AMBIGUOUS, str(r.reason))
        else:
            out[r.group] = (UNAVAILABLE, str(r.reason) or "no_core_beats_its_matched_null")
    return out


# =============================================================================== (a) nullgate
def nullgate_table(X, genes, groups: np.ndarray, cores: dict, *, n_draws: int = BAR["null_draws"],
                   seed: int = 0, min_cells: int = BAR["min_group_cells"]) -> tuple:
    """The nullgate statistic for every (group, lineage), plus the reference's per-group decision.

    This is ``hpg/nullgate_fast.annotate_groups`` with the inner loop kept open: same
    ``standardise_genes``, ``match_strata``, ``build_null_pools`` (a gene never stands in for
    itself), same RNG consumption order, same ``decide_group``.  ``tests/test_gate.py`` asserts
    the second return value equals ``annotate_groups`` frame for frame.  Returns
    ``(long, calls)``: ``long`` has one row per (group, lineage) with ``real``, ``null_mean``,
    ``null_p95``, ``z``, ``clears`` (= real > max(p95, 0)); ``calls`` is the GroupCall frame.
    """
    groups = np.asarray(groups).astype(str)
    Xs = standardise_genes(X)
    dec = match_strata(Xs)
    coarse = detection_deciles(Xs)
    rng = np.random.default_rng(seed)
    resolved = resolve_cores(genes, cores)
    pools = build_null_pools(genes, resolved, dec, coarse)
    labels_u, sizes, delta = group_gene_deltas(Xs, groups)
    order = {lab: i for i, lab in enumerate(labels_u)}
    long_rows, calls = [], []
    for g in pd.unique(groups):
        i = order[g]
        n = int(sizes[i])
        if n < min_cells:
            calls.append(GroupCall(g, n, "", "too_few_cells", np.nan, np.nan, np.nan, np.nan, ""))
            for lin in resolved:
                long_rows.append({"group": g, "n": n, "lineage": lin, "real": np.nan, "null_mean": np.nan,
                                  "null_p95": np.nan, "z": np.nan, "clears": False, "judgeable": False})
            continue
        d = delta[i]
        stats, nulls = {}, {}
        for lin, cols in resolved.items():
            cols = np.asarray(cols)
            real = float(d[cols].mean())
            picks = np.stack([[rng.choice(pools[lin][c]) for c in cols] for _ in range(n_draws)])
            null = d[picks].mean(axis=1)
            sd = float(null.std(ddof=1))
            p95 = float(np.percentile(null, 95))
            z = (real - float(null.mean())) / sd if sd > 0 else np.nan
            stats[lin] = (real, float(null.mean()), p95, z)
            nulls[lin] = null
            long_rows.append({"group": g, "n": n, "lineage": lin, "real": real, "null_mean": float(null.mean()),
                              "null_p95": p95, "z": z, "clears": bool(real > max(p95, 0.0)), "judgeable": True})
        calls.append(decide_group(g, n, stats, nulls))
    return pd.DataFrame(long_rows), pd.DataFrame([vars(c) for c in calls])


def award_table(long: pd.DataFrame, calls: pd.DataFrame, labels=LINEAGES) -> pd.DataFrame:
    """Per lineage LABEL: is its own core awarded on the cells the arm gave that label?"""
    rows = []
    present = set(long["group"].astype(str)) if len(long) else set()
    for L in labels:
        if L not in present:
            rows.append({"label": L, "n_cells": 0, "judgeable": False, "awarded": False,
                         "reason": "label_not_called"})
            continue
        own = long[(long["group"].astype(str) == L) & (long["lineage"] == L)]
        call = calls[calls["group"].astype(str) == L].iloc[0]
        if own.empty or not bool(own["judgeable"].iloc[0]):
            rows.append({"label": L, "n_cells": int(call["n_cells"]), "judgeable": False, "awarded": False,
                         "reason": str(call["reason"]) or "own core not judgeable"})
            continue
        o = own.iloc[0]
        rows.append({"label": L, "n_cells": int(o["n"]), "judgeable": True,
                     "real_own": float(o["real"]), "null_p95_own": float(o["null_p95"]), "z_own": float(o["z"]),
                     "awarded": bool(o["clears"]),
                     "decision_lineage": str(call["lineage"]), "decision_reason": str(call["reason"]),
                     "runners_up": str(call["runners_up"]),
                     "reason": "" if bool(o["clears"]) else "own_core_inside_matched_null"})
    return pd.DataFrame(rows)


# =============================================================================== (b) spatial
def spatial_coherence(labels: np.ndarray, label: str, xy: np.ndarray, *, sections: np.ndarray | None = None,
                      k: int = BAR["morans_k"], draws: int = BAR["null_draws"], q: float = BAR["null_quantile"],
                      seed: int = 0, max_cells: int = BAR["spatial_max_cells"], ids: np.ndarray | None = None,
                      nbr_cache: dict | None = None) -> dict:
    """Moran's I of the one-hot label against within-section label permutations.

    Judged cells are those the arm did not mark Unavailable.  When the section holds more than
    ``max_cells`` judged cells the graph is built on a hash subsample keyed by ``ids`` -- the
    SAME cells for every label (``nbr_cache`` shares the graph across labels).
    """
    labels = np.asarray(labels).astype(str)
    judged = labels != UNAVAILABLE
    idx = np.flatnonzero(judged)
    if len(idx) > max_cells:
        key_ids = np.asarray(ids)[idx] if ids is not None else idx.astype(str)
        idx = idx[hash_subsample(key_ids, max_cells, salt="moran")]
    lab = labels[idx]
    n_L = int((lab == label).sum())
    reason = gates.require_evidence(n_L, BAR["min_calls"], f"{label}: calls")
    if reason:
        return {"label": label, "I": np.nan, "null_p": np.nan, "n_calls": n_L, "judgeable": False,
                "clear": False, "reason": reason}
    pts = np.asarray(xy, dtype=float)[idx]
    if nbr_cache is not None and "nbr" in nbr_cache and nbr_cache.get("n") == len(idx):
        nbr = nbr_cache["nbr"]
    else:
        nbr = usability.neighbour_graph(pts, k=k)
        if nbr_cache is not None:
            nbr_cache.update(nbr=nbr, n=len(idx))
    one = (lab == label).astype(float)
    I_obs = float(usability.morans_i(one, nbr))
    g = np.asarray(sections)[idx] if sections is not None else np.zeros(len(idx), dtype=int)
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for d in range(draws):
        perm = one.copy()
        for gv in np.unique(g):
            sel = np.flatnonzero(g == gv)
            perm[sel] = one[rng.permutation(sel)]
        null[d] = usability.morans_i(perm, nbr)
    p = float(np.nanquantile(null, q))
    return {"label": label, "I": I_obs, "null_p": p, "null_mean": float(np.nanmean(null)), "n_calls": n_L,
            "n_cells_used": int(len(idx)), "judgeable": True, "clear": bool(I_obs > p), "reason": ""}


# =============================================================================== (c) stability
def stability(primary: dict, replicates: dict) -> dict:
    """``primary``: label -> spatial verdict dict; ``replicates``: name -> {label -> verdict}.
    Invariant iff every replicate's ``clear`` equals the primary's.  No replicates -> None."""
    out = {}
    for L, prim in primary.items():
        per = []
        for name, reps in replicates.items():
            r = reps.get(L, {})
            per.append({"replicate": name, "spatial_clear": bool(r.get("clear", False)),
                        "judgeable": bool(r.get("judgeable", False)), "n_calls": int(r.get("n_calls", 0))})
        out[L] = {"per_replicate": per,
                  "invariant": (all(p["spatial_clear"] == bool(prim.get("clear", False)) for p in per)
                                if per else None)}
    return out


# =============================================================================== (d) agreement
def agreement(a: np.ndarray, b: np.ndarray) -> dict:
    """ARI / AMI on cells neither side marked Unavailable (Ambiguous is a label)."""
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    a = np.asarray(a).astype(str)
    b = np.asarray(b).astype(str)
    both = (a != UNAVAILABLE) & (b != UNAVAILABLE)
    n = int(both.sum())
    if n < BAR["min_calls"]:
        return {"ari": np.nan, "ami": np.nan, "n": n}
    return {"ari": float(adjusted_rand_score(a[both], b[both])),
            "ami": float(adjusted_mutual_info_score(a[both], b[both])), "n": n}


def cross_arm_agreement(donor: str, arms=pins.ARMS) -> dict:
    """Descriptive only: never chooses an arm, a k, or a threshold."""
    have = {}
    for arm in arms:
        try:
            c = read_calls(arm, donor)
        except FileNotFoundError:
            continue
        have[arm] = c.set_index(c["cell_id"].astype(str))["lineage_native"].astype(str)
    out = {}
    names = list(have)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            common = have[a].index.intersection(have[b].index)
            out[f"{a}|{b}"] = agreement(have[a].loc[common].to_numpy(), have[b].loc[common].to_numpy())
    return out


# =============================================================================== (e) IF association
def if_association(labels: np.ndarray, values: np.ndarray, label: str = "Endocrine") -> dict:
    """Rank association of a same-section IF per-cell value with membership of ``label``.

    AUROC (Mann-Whitney U on ranks; 0.5 = no association) and Spearman rho of the value against
    the one-hot label, over cells with a finite value.  DESCRIPTIVE ONLY: a protein intensity is
    evidence about the cells the arm called, not a call, and no threshold is applied to it.
    """
    from scipy.stats import rankdata
    labels = np.asarray(labels).astype(str)
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v) & (labels != UNAVAILABLE)
    is_L = labels == label
    n1, n0 = int((finite & is_L).sum()), int((finite & ~is_L).sum())
    out = {"label": label, "n_finite": int(finite.sum()), "n_in": n1, "n_out": n0, "descriptive_only": True}
    reason = gates.require_evidence(n1, BAR["min_calls"], f"{label}: cells with IF value") \
        or gates.require_evidence(n0, BAR["min_calls"], f"not-{label}: cells with IF value")
    if reason:
        out.update(judgeable=False, reason=reason, auroc=np.nan, spearman_rho=np.nan)
        return out
    r = rankdata(v[finite])
    pos = is_L[finite]
    U = float(r[pos].sum() - n1 * (n1 + 1) / 2.0)
    rho = float(np.corrcoef(r, pos.astype(float))[0, 1])
    per_label = {}
    for L in np.unique(labels[finite]):
        per_label[str(L)] = {"n": int((labels[finite] == L).sum()),
                             "median_value": float(np.median(v[finite][labels[finite] == L]))}
    out.update(judgeable=True, reason="", auroc=U / (n1 * n0), spearman_rho=rho,
               median_in=float(np.median(v[finite & is_L])), median_out=float(np.median(v[finite & ~is_L])),
               per_label_median=per_label)
    return out


# =============================================================================== data for judge
def load_base_for_gate(donor: str) -> dict:
    """obs_names, var_names, lognorm (CSR), spatial (microns) from the donor's base -- by h5py so
    the counts layer is never read here."""
    import h5py
    from anndata.io import read_elem
    p = pins.base_path(donor)
    if not p.exists():
        raise FileNotFoundError(f"base for donor {donor} is absent: {p}")
    with h5py.File(p, "r") as f:
        obs_names = np.asarray(read_elem(f["obs"]).index).astype(str)
        var = read_elem(f["var"])
        var_names = np.asarray(var.index).astype(str)
        X = read_elem(f["layers/lognorm"])
        xy = np.asarray(read_elem(f["obsm/spatial"]), dtype=float)
        tier = int((var["identity_panel"].to_numpy(bool) & var["is_state_only"].to_numpy(bool)).sum()) \
            if {"identity_panel", "is_state_only"} <= set(var.columns) else None
    return {"obs_names": obs_names, "var_names": var_names, "lognorm": X, "xy": xy,
            "identity_state_overlap": tier, "path": str(p)}


def load_if_values(donor: str, base_ids: np.ndarray, markers=("INS", "GCG")) -> dict:
    """Restain per-cell IF values aligned to ``base_ids`` through the bundle's cells.parquet
    order (that is the order the .npy arrays were written in).  Absent -> {}."""
    out = {}
    cells_p = pins.cells_parquet(donor)
    if not cells_p.exists():
        return out
    order = None
    for m in markers:
        p = pins.if_values_path(donor, m)
        if not p.exists():
            continue
        vals = np.load(p)
        if order is None:
            cid = pd.read_parquet(cells_p, columns=["cell_id"])["cell_id"].astype(str).to_numpy()
            if len(cid) != len(vals):
                raise ValueError(f"{p} has {len(vals)} values for {len(cid)} bundle cells")
            order = pd.Series(np.arange(len(cid)), index=cid).reindex(np.asarray(base_ids).astype(str))
            if order.isna().any():
                raise ValueError(f"{int(order.isna().sum())} base cells of donor {donor} are not in cells.parquet")
            order = order.to_numpy().astype(int)
        out[m] = {"values": vals[order], "path": str(p)}
    return out


# =============================================================================== (f) judge
def judge_arrays(*, ids: np.ndarray, labels: np.ndarray, X, genes, xy: np.ndarray, cores: dict,
                 replicates: dict | None = None, if_values: dict | None = None, sections: np.ndarray | None = None,
                 stochastic: bool = True, seed: int = 0, draws: int = BAR["null_draws"], arm: str = "",
                 donor: str = "") -> tuple:
    """Blocking + descriptive statistics for every lineage label of one arm on one donor.

    ``labels`` aligned to ``ids`` (the base cells); ``X`` the base lognorm over ``genes``;
    ``replicates`` name -> labels aligned to ``ids``.  Returns ``(verdict dict, nullgate long)``.
    """
    labels = np.asarray(labels).astype(str)
    ids = np.asarray(ids).astype(str)
    n = len(labels)
    out = {"arm": arm, "donor_id": donor, "n_cells": n, "seed": seed, "null_draws": draws,
           "label_counts": {k: int(v) for k, v in pd.Series(labels).value_counts().items()},
           "bar": dict(BAR), "caveat": FAMILY_WISE_CAVEAT}
    # (a) nullgate on the labels
    long = pd.DataFrame()
    n_groups = len(np.unique(labels))
    if n_groups < 2:
        out["nullgate"] = {"evaluable": False, "reason": f"{n_groups} label group(s): the contrast has no outside"}
        award = pd.DataFrame([{"label": L, "n_cells": int((labels == L).sum()), "judgeable": False,
                               "awarded": False, "reason": "nullgate_not_evaluable"} for L in LINEAGES])
    else:
        try:
            long, calls = nullgate_table(X, genes, labels, cores, n_draws=draws, seed=seed)
            award = award_table(long, calls)
            out["nullgate"] = {"evaluable": True, "groups": calls.to_dict(orient="records")}
        except NullGateError as exc:
            out["nullgate"] = {"evaluable": False, "reason": str(exc)}
            award = pd.DataFrame([{"label": L, "n_cells": int((labels == L).sum()), "judgeable": False,
                                   "awarded": False, "reason": "nullgate_not_evaluable"} for L in LINEAGES])
    award = award.set_index("label")
    # (b) spatial coherence on the primary labels, one shared graph
    cache: dict = {}
    spatial = {L: spatial_coherence(labels, L, xy, sections=sections, seed=seed, draws=draws, ids=ids,
                                    nbr_cache=cache) for L in LINEAGES}
    # (c) stability over replicates (each replicate on its own labels, same graph cells)
    reps_spatial = {}
    for name, rl in (replicates or {}).items():
        rl = np.asarray(rl).astype(str)
        rc: dict = {}
        reps_spatial[name] = {L: spatial_coherence(rl, L, xy, sections=sections, seed=seed, draws=draws, ids=ids,
                                                   nbr_cache=rc) for L in LINEAGES}
    stab = stability(spatial, reps_spatial)
    # (d) agreement with each replicate, descriptive
    out["replicate_agreement"] = {name: agreement(labels, np.asarray(rl)) for name, rl in (replicates or {}).items()}
    # (e) IF association, descriptive, both hormones
    out["if_association"] = {m: if_association(labels, v["values"]) | {"path": v["path"]}
                             for m, v in (if_values or {}).items()}
    # verdict per lineage
    per = {}
    for L in LINEAGES:
        a = award.loc[L].to_dict() if L in award.index else {"judgeable": False, "awarded": False}
        a = {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in a.items()}
        s = spatial[L]
        st = stab[L]
        if replicates:
            stable = st["invariant"]
            stable_reason = ""
        elif stochastic:
            stable, stable_reason = None, "stochastic arm with no replicates: stability unjudged"
        else:
            stable, stable_reason = True, "deterministic arm: stability holds by construction"
        blocking = [bool(a.get("awarded", False)), bool(s["clear"]), bool(stable)]
        per[L] = {"award": a, "spatial": s, "stability": {**st, "reason": stable_reason},
                  "verified": all(blocking),
                  "verified_without_stability": bool(a.get("awarded", False)) and bool(s["clear"])}
    out["lineages"] = per
    out["n_verified"] = int(sum(v["verified"] for v in per.values()))
    out["verified_lineages"] = sorted(L for L, v in per.items() if v["verified"])
    return out, long


def judge(arm: str, donor: str, *, seed: int = 0, draws: int = BAR["null_draws"], write: bool = True) -> dict:
    """Judge one arm on one donor from its written calls; writes ``gate/<arm>/donor_id=<d>.json``
    and the nullgate long table beside it.  Prints the family-wise caveat."""
    calls = read_calls(arm, donor)
    rec = read_record(arm, donor)
    base = load_base_for_gate(donor)
    ids = base["obs_names"]
    c = calls.set_index(calls["cell_id"].astype(str))
    if set(c.index) != set(ids):
        raise ValueError(f"{arm} donor {donor}: calls cover {len(c)} cells, base has {len(ids)}")
    labels = c["lineage_native"].reindex(ids).astype(str).to_numpy()
    reps = None
    rl = read_replicates(arm, donor)
    if rl is not None:
        rl = rl.set_index(rl["cell_id"].astype(str)).reindex(ids)
        reps = {col: rl[col].astype(str).to_numpy() for col in rl.columns if col != "cell_id"}
    cores = pins.identity_cores()
    ifv = load_if_values(donor, ids)
    verdict, long = judge_arrays(ids=ids, labels=labels, X=base["lognorm"], genes=base["var_names"], xy=base["xy"],
                                 cores=cores, replicates=reps, if_values=ifv, stochastic=bool(rec.get("stochastic", True)),
                                 seed=seed, draws=draws, arm=arm, donor=donor)
    verdict["inputs"] = {"calls": str(pins.arm_dir(arm, donor) / "calls.parquet"), "base": base["path"],
                         "identity_state_overlap": base["identity_state_overlap"],
                         "identity_cores": {L: len(g) for L, g in cores.items()},
                         "record_written_utc": rec.get("written_utc"), "panel_snapshot": pins.panel_snapshot()}
    verdict["cross_arm_agreement"] = cross_arm_agreement(donor)
    verdict["code"] = pins.git_head()
    verdict["written_utc"] = pins.utc_now()
    if write:
        d = pins.gate_dir(arm)
        d.mkdir(parents=True, exist_ok=True)
        pins.gate_json(arm, donor).write_text(json.dumps(verdict, indent=1, default=_json_default))
        if len(long):
            long.to_parquet(d / f"donor_id={donor}__nullgate.parquet", index=False)
    print(f"[gate] {arm} donor {donor}: verified {verdict['n_verified']}/{len(LINEAGES)} "
          f"{verdict['verified_lineages']}")
    for L, v in verdict["lineages"].items():
        a, s = v["award"], v["spatial"]
        print(f"    {L:<10} n={a.get('n_cells', 0):>8,}  awarded={a.get('awarded')!s:<5} "
              f"real={a.get('real_own')}  p95={a.get('null_p95_own')}  "
              f"I={s['I']:.4f} p95={s['null_p']:.4f} coherent={s['clear']!s:<5} "
              f"stable={v['stability']['invariant']}  verified={v['verified']}")
    print("  " + FAMILY_WISE_CAVEAT)
    return verdict


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="judge a Xenium typing arm on donors")
    ap.add_argument("--arm", required=True, choices=pins.ARMS)
    ap.add_argument("--donors", nargs="+", default=["6436", "6442"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--draws", type=int, default=BAR["null_draws"])
    a = ap.parse_args(argv)
    pins.assert_pins()
    for d in a.donors:
        judge(a.arm, d, seed=a.seed, draws=a.draws)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
