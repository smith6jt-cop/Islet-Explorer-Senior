"""Shared helpers for every arm, so each arm is only its method.

An arm module exposes

    ARM = "<name>"                                      # the arm name (a new configuration = a new name)
    IDENTITY_MARKERS(donor) -> list[str]                # the type markers it uses (declarations.identity_markers ∩ donor panel)
    run(ctx: ArmContext) -> pd.DataFrame                # calls with records.CALL_COLUMNS for every ctx.ids

and a CLI built from ``cli(main_fn)``: ``--donor D --scope frame|donor [--ablate MARKER] [--seed N]
[--subsample-frac F] [--replicate NAME]``.  ``ArmContext`` carries the aligned Step-2 values for the
frozen cell set (frame) or the whole donor, ALREADY restricted to identity (type) markers and with
the ablated marker removed; the arm never sees a state marker.  ``write`` appends the record.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402
from gate import hash_subsample  # noqa: E402

CONTEXT_COLUMNS = ("object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid", "cell_region", "islet_num")


def frame_cells(donor: str) -> pd.DataFrame:
    p = pins.OUT / "frame" / "cells" / f"donor_id={donor}.parquet"
    if not p.exists():
        raise FileNotFoundError(f"frame not built for donor {donor}: run scripts/broadtype/frame.py first")
    return pd.read_parquet(p)


def identity_markers_for(donor: str) -> list:
    """Declared identity (type) markers this donor actually carries, in declaration order."""
    dec = declarations.load()
    import panel
    have = set(panel.donor_markers(donor))
    return [m for m in dec["identity_markers"] if m in have]


def state_markers_for(donor: str) -> list:
    dec = declarations.load()
    import panel
    have = set(panel.donor_markers(donor))
    return [m for m in dec["state_markers"] if m in have]


@dataclass
class ArmContext:
    donor_id: str
    scope: str                       # 'frame' | 'donor'
    ids: np.ndarray                  # object_id, the frozen cell set (every one gets a row)
    X: pd.DataFrame                  # aligned values, columns = identity markers (ablated marker removed)
    context: pd.DataFrame            # CONTEXT_COLUMNS aligned to ids
    markers: list                    # X.columns as a list
    ablated: str | None = None
    seed: int = CRITERIA["primary_seed"]
    subsample_frac: float | None = None
    replicate: str | None = None
    region_idx: np.ndarray | None = None
    extra: dict = field(default_factory=dict)

    @property
    def config_id(self) -> str:
        if self.ablated:
            return f"ablate:{self.ablated}"
        if self.replicate:
            return self.replicate
        return "primary"


def load_context(donor: str, scope: str, *, ablate: str | None = None, seed: int = CRITERIA["primary_seed"],
                 subsample_frac: float | None = None, replicate: str | None = None) -> ArmContext:
    tab = pd.read_parquet(pins.uniform_parquet(donor))
    tab["object_id"] = tab["object_id"].astype(str)
    markers = identity_markers_for(donor)
    region_idx = None
    if scope == "frame":
        fc = frame_cells(donor)
        fc = fc.loc[fc["in_expression"]]
        tab = tab.merge(fc[["object_id", "region_idx"]], on="object_id", how="inner")
        region_idx = tab["region_idx"].to_numpy()
    elif scope != "donor":
        raise ValueError(f"scope must be frame or donor, got {scope!r}")
    if subsample_frac is not None and 0 < subsample_frac < 1:
        keep = hash_subsample(tab["object_id"].to_numpy(), int(round(subsample_frac * len(tab))),
                              salt=f"{replicate or 'subsample'}|{seed}")
        tab = tab.iloc[keep].reset_index(drop=True)
        region_idx = region_idx[keep] if region_idx is not None else None
    if ablate:
        if ablate not in markers:
            raise KeyError(f"{ablate!r} is not an identity marker on donor {donor}")
        markers = [m for m in markers if m != ablate]
    X = tab[markers].astype(float).reset_index(drop=True)
    ctx = tab[[c for c in CONTEXT_COLUMNS if c in tab.columns]].reset_index(drop=True)
    return ArmContext(donor_id=str(donor), scope=scope, ids=tab["object_id"].to_numpy(), X=X, context=ctx,
                      markers=markers, ablated=ablate, seed=seed, subsample_frac=subsample_frac,
                      replicate=replicate, region_idx=region_idx)


def within_donor_rank(X: pd.DataFrame) -> pd.DataFrame:
    """Rank in [0, 1] per marker within the loaded cells; NaN stays NaN."""
    out = X.copy()
    for c in out.columns:
        v = out[c].to_numpy(float)
        ok = np.isfinite(v)
        r = np.full(v.shape, np.nan)
        if ok.sum() > 1:
            r[ok] = (np.argsort(np.argsort(v[ok], kind="mergesort")) + 0.5) / ok.sum()
        out[c] = r
    return out


def state_columns(donor: str, ids: np.ndarray) -> pd.DataFrame:
    """The state tier, computed AFTER identity: within-donor rank of each state marker for ``ids``."""
    tab = pd.read_parquet(pins.uniform_parquet(donor))
    tab["object_id"] = tab["object_id"].astype(str)
    have = set(state_markers_for(donor))
    sm = [c for c in tab.columns if c in have]          # the table's own column order
    ranks = within_donor_rank(tab[sm])
    ranks.columns = [f"state_{m}_donor_rank" for m in sm]
    ranks.insert(0, "object_id", tab["object_id"].to_numpy())
    out = ranks.set_index("object_id").reindex(pd.Index(ids.astype(str), name="object_id"))
    return out.reset_index()


def write(ctx: ArmContext, arm: str, calls: pd.DataFrame, *, seconds: float, env: str, params: dict | None = None,
          notes: str = "", extra_pins: dict | None = None) -> pathlib.Path:
    dec = declarations.load()
    calls = records.fill_unavailable(calls, ctx.ids)
    rec = records.ArmRecord(
        arm=arm, donor_id=ctx.donor_id, config_id=ctx.config_id, identity_markers=list(ctx.markers),
        scope=ctx.scope, n_cells=int(len(calls)), n_called=int((~calls["label"].isin([records.AMBIGUOUS, records.UNAVAILABLE])).sum()),
        n_ambiguous=int((calls["label"] == records.AMBIGUOUS).sum()), n_unavailable=int((calls["label"] == records.UNAVAILABLE).sum()),
        seconds=float(seconds), env=env, params=dict(params or {}), notes=notes,
        pinned={"run_id": pins.RUN_ID, "expression_content_id": pins.EXPRESSION_CONTENT_ID,
                "declarations": dec["method_version"], **(extra_pins or {})})
    return records.write_calls(calls, rec, labels=set(dec["labels"]), expected_ids=ctx.ids)


def cli(arm: str, run_fn, *, env: str, default_params: dict | None = None) -> int:
    ap = argparse.ArgumentParser(description=f"broadtype arm {arm}")
    ap.add_argument("--donor", required=True)
    ap.add_argument("--scope", choices=("frame", "donor"), default="frame")
    ap.add_argument("--ablate", default=None)
    ap.add_argument("--seed", type=int, default=CRITERIA["primary_seed"])
    ap.add_argument("--subsample-frac", type=float, default=None)
    ap.add_argument("--replicate", default=None, help="record name for a seed/subsample replicate")
    a = ap.parse_args()
    pins.assert_pins()
    ctx = load_context(a.donor, a.scope, ablate=a.ablate, seed=a.seed, subsample_frac=a.subsample_frac,
                       replicate=a.replicate)
    t0 = time.time()
    calls = run_fn(ctx)
    path = write(ctx, arm, calls, seconds=time.time() - t0, env=env, params=default_params)
    print(f"{arm} {a.donor} {ctx.config_id}: wrote {path} ({len(ctx.ids):,} cells, {time.time() - t0:.1f}s)")
    return 0


def geojson_object_ids(donor: str) -> np.ndarray:
    """The donor's GeoJSON feature ids in stream order, so ``oids[label - 1]`` is a mask label's UUID.

    One full stream of the 3-4 GB GeoJSON (~45 s); cached under ``data/development/broadtype/scratch``
    keyed by the manifest's recorded GeoJSON sha256, so a regenerated GeoJSON invalidates the cache
    instead of silently mislabelling every cell.
    """
    from phenocycler.artifacts import iter_qupath_geojson_features
    im = pins.manifest_image(donor)
    sha = str(im["cell_geojson"].get("sha256", ""))[:16] or "nosha"
    cache = pins.OUT / "scratch" / "oids" / f"{donor}__{sha}.npy"
    if cache.exists():
        return np.load(cache, allow_pickle=True)
    ids = np.array([str(ft.get("id")) for ft in iter_qupath_geojson_features(pins.geojson_path(donor))], dtype=object)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache, ids, allow_pickle=True)
    return ids


def frame_regions(donor: str) -> list:
    """The donor's frozen regions with their 2048-px windows, from ``frame/regions.json``."""
    import json
    reg = json.loads((pins.OUT / "frame" / "regions.json").read_text())
    return reg["donors"][str(donor)]["regions"]


def frame_mask(donor: str, region_idx: int) -> np.ndarray:
    return np.load(pins.OUT / "frame" / "masks" / f"donor_id={donor}" / f"w{region_idx:02d}.npy")
