"""Arm D ``celesta``: CELESTA (Zhang et al. 2022) at its package DEFAULTS, one frame region at a time.

CELESTA turns each marker into an activation probability (a two-component Gaussian mixture on
``asinh(x / 10)``, a sigmoid through the components' crossing, min-max scaled to [0, 1]), scores
every cell against a 1/0/NA signature per type (``1 - mean squared deviation`` over the
signature's non-NA entries, normalised across types), anchors the cells whose argmax type
passes the anchor thresholds (all its ``1`` markers >= 0.7, all its ``0`` markers <= 0.9), and
then assigns the rest iteratively with a spatial term -- the type probabilities of the five
nearest neighbours, weighted by distance to the nearest anchored cell of that type within a
100-unit bandwidth -- under the looser index thresholds (>= 0.5 / <= 1.0), refitting the prior
matrix from the assigned cells each pass.  Cells that never pass stay ``Unknown``.

What this arm fixes, and why:

* **Signature from the declarations, nothing else.**  One lineage round; all seven declared
  labels are top level (lineage codes ``1_0_1 .. 1_0_7``, the flat round-1 form the package's
  own ``prior_marker_info`` uses for its five top-level types).  Code ``1`` for a label's
  declared identity markers, ``0`` for the markers in its ``negatives`` list, ``NA`` otherwise.
  Every column is a registry ``type`` marker because ``ctx.markers`` already is.
* **Package defaults, verbatim, no FilterCells.**  The R half writes the defaults it actually
  ran (``formals()``) to ``env.json``; this side refuses to record ``CELESTA_DEFAULTS`` if the
  package disagrees, so the recorded parameters can never drift from the executed ones.
* **In-method prior, flagged:** ``FitGmmModel`` fits Rmixmod ``Gaussian_p_Lk_Ck`` with
  ``equal.proportions = TRUE`` -- a 50/50 mixture per marker.  That is a frequency-shaped
  assumption inside the package, not one this lane makes; the protocol flags it and this arm
  records it (``params["in_method_prior"]``) rather than tuning around it.
* **Per frame region.**  CELESTA's neighbourhood machinery and its per-cell R loops are built
  for one field of view; the frame's regions are the pre-registered contiguous cell sets.  The
  GMM is therefore fit per region, so a marker's activation scale is region-local -- a property
  of the method as run, recorded here.  A region with <= 10 finite cells is
  ``Unavailable(region_too_small)``: 10 is CELESTA's own floor (``CreateCelestaObject`` returns
  NULL there), not a choice of ours.
* **Coordinates in microns.**  ``X_centroid`` / ``Y_centroid`` are microns (the dct lane
  verified this against the mask), so the default bandwidth 100 is 100 um and ``k = 5``
  neighbours are the five nearest by Euclidean distance in the plane.
* **Missing is unavailable.**  A cell with any non-finite identity value is
  ``Unavailable(missing_measurement)`` and is left out of its region's CELESTA input -- and so
  out of that region's neighbour graph, which is recorded per region.
* **Negatives left-censored at 0, counted.**  ``asinh(x / 10)`` is defined on negatives, but a
  negative aligned value is an alignment artefact, not a fluorescence.  The censor count per
  marker is written to the diagnostics; on the cohort's aligned values it is zero.
* **Row order fixed.**  ``FitGmmModel`` thins exact zeros by dropping the FIRST n zero indices,
  so the fit depends on row order; each region is sorted by ``object_id`` before it is written,
  making the result invariant to upstream row order.
* **Seeded, because the package is not.**  ``FitGmmModel`` calls ``set.seed(1)``, but Rmixmod's
  ``mixmodCluster(seed = NULL)`` seeds its C++ EM from its own source, so two identical calls
  return different fits (verified 2026-09-05; the scores moved in the third decimal, the labels
  did not, on the synthetic plant).  The R half supplies ``seed = ctx.seed`` to
  ``mixmodCluster`` through a global wrapper -- CELESTA's NAMESPACE imports nothing, so its
  unqualified call resolves there first -- and changes nothing else.  The arm is therefore
  stochastic in the contract's sense: the primary runs at seed 0 and the screen's seed
  replicates measure the EM-initialisation sensitivity.
* **Score** = the assigned type's entry of ``@current_cell_prob`` (the normalised score for
  anchored cells, the neighbour-updated probability for the rest); NaN when not assigned.
  ``Unknown`` -> ``Ambiguous(celesta_unassigned)``.
* **Wall-clock cap** 60 min for the frame.  The R half writes each region's result as it
  finishes, so at the cap the finished regions stand and the rest are
  ``Unavailable(celesta_wall_clock_cap)``.

Recorded deviations / extensions (the task froze the frame path only):

* ``scope == "donor"`` has no ``region_idx``; the arm tiles the donor into ``core_px``-sized
  squares (1024 px = ~520 um, the frame's core size) on the centroid grid and treats each tile
  as a region, so the per-region method is the same object at both scopes.  The cap scales
  with the tile count (60 min per 16 regions).
* The pinned CELESTA commit's NAMESPACE does not import ``%<-%`` from zeallot although every
  constructor line uses it; the R half attaches zeallot explicitly.  A packaging defect of the
  pin, recorded, not patched.

Pins: CELESTA 0.0.0.9000 at ``smith6jt-cop/CELESTA@e5e92ffcef5ccada8053c2e042f34c446b6b3565``,
spdep 1.4.2, Rmixmod 2.1.10, R 4.5.1 (``/usr/local/bin/Rscript``).  The R half is
``arms/celesta_run.R`` beside this file.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/arms/celesta.py --donor 6436 --scope frame
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402

ARM = "celesta"
RSCRIPT = pathlib.Path("/usr/local/bin/Rscript")
R_SCRIPT = HERE / "celesta_run.R"

#: The CELESTA install this arm was written against.  ``r_env()`` refuses any other commit.
CELESTA_REMOTE = "smith6jt-cop/CELESTA"
CELESTA_COMMIT = "e5e92ffcef5ccada8053c2e042f34c446b6b3565"
CELESTA_VERSION = "0.0.0.9000"
SPDEP_VERSION = "1.4.2"
RMIXMOD_VERSION = "2.1.10"

#: The package defaults this arm runs at.  Read back from ``formals()`` by the R half and
#: compared to these before anything is recorded; a mismatch stops the arm.
CELESTA_DEFAULTS = {
    "max_iteration": 10,
    "cell_change_threshold": 0.01,
    "min_diff": 0,
    "min_probability": 0,
    "high_expression_threshold_anchor": 0.7,
    "low_expression_threshold_anchor": 0.9,
    "high_expression_threshold_index": 0.5,
    "low_expression_threshold_index": 1,
    "scale_factor": 5,
    "bandwidth": 100,
    "cofactor_argument": 10,
    "cofactor_effective": 10,
    "transform_type": 1,
    "number_of_neighbors": 5,
    "neighbor_bandwidth": 100,
}

#: ``CreateCelestaObject`` returns NULL at ``length(cell_ID) <= 10``; a region at or below it is
#: Unavailable(region_too_small).  The package's floor, not ours.
MIN_REGION_CELLS = 10
#: The frame's cap; donor scope scales it by tile count (see ``wall_clock_cap_s``).
WALL_CLOCK_CAP_S = 60 * 60
CAP_REGIONS = CRITERIA["n_blind"] + CRITERIA["n_islet"]
#: Donor-scope tiling: the frame's core size, in microns.
DONOR_TILE_UM = CRITERIA["core_px"] * pins.PIXEL_SIZE_UM

SCRATCH_ROOT = pins.OUT / "scratch" / ARM

DEFAULT_PARAMS = {
    "method": "CELESTA per region at package defaults; one lineage round; 7 top-level labels 1_0_1..1_0_7",
    "signature": "1 = declared identity markers, 0 = declarations.json negatives, NA otherwise",
    "celesta_defaults": CELESTA_DEFAULTS,
    "filter_cells": False,
    "in_method_prior": "Rmixmod equal.proportions=TRUE (a 50/50 mixture per marker); flagged, not tuned",
    "transform": "asinh(x / 10) inside CELESTA (cofactor literal 10)",
    "values": "ctx.X aligned values, negatives left-censored at 0 (count recorded in diagnostics)",
    "coordinates": "X_centroid / Y_centroid in microns; bandwidth 100 = 100 um; k = 5 nearest",
    "regions": "frame: ctx.region_idx; donor: %.1f um grid tiles (core_px * pixel size)" % DONOR_TILE_UM,
    "min_region_cells": MIN_REGION_CELLS,
    "row_order": "sorted by object_id within region (FitGmmModel's zero thinning is order-dependent)",
    "missing": "any non-finite identity value -> Unavailable(missing_measurement); excluded from the region input",
    "unknown": "CELESTA 'Unknown' -> Ambiguous(celesta_unassigned)",
    "score": "assigned type's entry of @current_cell_prob; NaN when not assigned",
    "wall_clock_cap_s_per_%d_regions" % CAP_REGIONS: WALL_CLOCK_CAP_S,
    "seed": "ctx.seed -> Rmixmod mixmodCluster(seed=) via a global wrapper; the package's set.seed(1) does not reach Rmixmod",
    "pins": {"celesta": f"{CELESTA_REMOTE}@{CELESTA_COMMIT}", "celesta_version": CELESTA_VERSION,
             "spdep": SPDEP_VERSION, "rmixmod": RMIXMOD_VERSION, "rscript": str(RSCRIPT)},
    "zeallot_attached_explicitly": True,
}

REASON_MISSING = "missing_measurement"
REASON_SMALL = "region_too_small"
REASON_UNASSIGNED = "celesta_unassigned"
REASON_CAP = "celesta_wall_clock_cap"
REASON_ERROR = "celesta_error"
REASON_ACTIVATION = "celesta_error:activation_nan"
REASON_ABLATION = "ablation_removes_required_input"


# --------------------------------------------------------------------------- the signature
def build_prior(markers: list, dec: dict | None = None) -> pd.DataFrame:
    """The ``prior_marker_info`` frame: col 1 cell type, col 2 lineage code, cols 3+ markers.

    Every declared label with at least one declared marker among ``markers`` gets a row; a label
    whose declared markers are all absent (an ablation could do this only to a one-marker label,
    and none is declared) is dropped, because a row with no ``1`` would be assignable on its
    ``0`` markers alone.  Codes: 1 declared, 0 negative, NA otherwise; columns in ``markers``
    order so the R side and this side agree on the marker set without renaming.
    """
    dec = dec or declarations.load()
    rows = []
    k = 0
    for label, spec in dec["labels"].items():
        pos = [m for m in spec["markers"] if m in markers]
        if not pos:
            continue
        k += 1
        row = {"cell_type": label, "lineage_level": f"1_0_{k}"}
        neg = set(spec["negatives"]) - set(pos)
        for m in markers:
            row[m] = 1 if m in pos else (0 if m in neg else pd.NA)
        rows.append(row)
    out = pd.DataFrame(rows, columns=["cell_type", "lineage_level", *markers])
    for m in markers:
        out[m] = out[m].astype("Int64")
    return out


def dropped_labels(markers: list, dec: dict | None = None) -> list:
    dec = dec or declarations.load()
    return [L for L, s in dec["labels"].items() if not any(m in markers for m in s["markers"])]


# --------------------------------------------------------------------------- regions
def region_ids(ctx: arm_base.ArmContext) -> np.ndarray:
    """Frame: the frozen ``region_idx``.  Donor: ``core_px``-sized grid tiles on the centroids."""
    if ctx.region_idx is not None:
        return np.array([f"{int(r):02d}" for r in ctx.region_idx], dtype=object)
    x = ctx.context["X_centroid"].to_numpy(float)
    y = ctx.context["Y_centroid"].to_numpy(float)
    ix = np.floor(x / DONOR_TILE_UM).astype(int)
    iy = np.floor(y / DONOR_TILE_UM).astype(int)
    return np.array([f"t{a:03d}_{b:03d}" for a, b in zip(iy, ix)], dtype=object)


def wall_clock_cap_s(n_regions: int) -> int:
    return WALL_CLOCK_CAP_S * max(1, math.ceil(n_regions / CAP_REGIONS))


# --------------------------------------------------------------------------- R side
_R_OK: bool | None = None


def r_available() -> bool:
    """Rscript present and CELESTA (with zeallot) loadable; tests skip cleanly otherwise."""
    global _R_OK
    if _R_OK is None:
        if not RSCRIPT.exists():
            _R_OK = False
        else:
            try:
                p = subprocess.run([str(RSCRIPT), "-e", "suppressMessages({library(zeallot); library(CELESTA); library(spdep)})"],
                                   capture_output=True, text=True, timeout=300)
                _R_OK = p.returncode == 0
            except (OSError, subprocess.TimeoutExpired):
                _R_OK = False
    return _R_OK


def r_env() -> str:
    """The ``env`` string for the record, from the installed packages -- and the commit check."""
    code = ("suppressMessages({library(zeallot); library(CELESTA); library(spdep); library(Rmixmod)});"
            "d <- read.dcf(system.file('DESCRIPTION', package='CELESTA'));"
            "cat(R.version.string, as.character(packageVersion('CELESTA')), d[1,'RemoteSha'],"
            "    as.character(packageVersion('spdep')), as.character(packageVersion('Rmixmod')), sep='|')")
    p = subprocess.run([str(RSCRIPT), "-e", code], capture_output=True, text=True, check=True, timeout=300)
    rv, cv, sha, sv, mv = p.stdout.strip().splitlines()[-1].split("|")
    if sha != CELESTA_COMMIT:
        raise RuntimeError(f"installed CELESTA is {sha}, this arm is pinned to {CELESTA_COMMIT}")
    if (cv, sv, mv) != (CELESTA_VERSION, SPDEP_VERSION, RMIXMOD_VERSION):
        raise RuntimeError(f"installed versions CELESTA {cv} / spdep {sv} / Rmixmod {mv} differ from the pins "
                           f"{CELESTA_VERSION} / {SPDEP_VERSION} / {RMIXMOD_VERSION}")
    return f"{rv} ({RSCRIPT}); CELESTA {cv} {CELESTA_REMOTE}@{sha[:8]}; spdep {sv}; Rmixmod {mv}"


def _check_defaults(env: dict) -> None:
    got = env.get("defaults", {})
    bad = {k: (v, got.get(k)) for k, v in CELESTA_DEFAULTS.items()
           if got.get(k) is None or not math.isclose(float(got[k]), float(v))}
    if bad:
        raise RuntimeError(f"CELESTA package defaults differ from the recorded ones: {bad}")


def _run_r(scratch: pathlib.Path, cap_s: int, seed: int) -> dict:
    """Run the R half under the cap; a timeout is a recorded outcome, not an exception."""
    t0 = time.time()
    cmd = [str(RSCRIPT), str(R_SCRIPT), "--scratch", str(scratch), "--seed", str(int(seed))]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=cap_s)
        return {"returncode": p.returncode, "timed_out": False, "seconds": time.time() - t0,
                "stderr_tail": p.stderr[-2000:], "stdout_tail": p.stdout[-500:]}
    except subprocess.TimeoutExpired as exc:
        err = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        return {"returncode": None, "timed_out": True, "seconds": time.time() - t0, "stderr_tail": err[-2000:]}


# --------------------------------------------------------------------------- the arm
def run(ctx: arm_base.ArmContext, *, scratch_root: pathlib.Path | None = None,
        diagnostics_dir: pathlib.Path | None = None) -> pd.DataFrame:
    ids = np.asarray(ctx.ids).astype(str)
    n = len(ids)
    label = np.full(n, records.UNAVAILABLE, dtype=object)
    score = np.full(n, np.nan)
    reason = np.full(n, "", dtype=object)
    diag: dict = {"arm": ARM, "donor_id": ctx.donor_id, "scope": ctx.scope, "config_id": ctx.config_id,
                  "identity_markers": list(ctx.markers), "n_cells": int(n)}

    markers = list(ctx.markers)
    dec = declarations.load()
    prior = build_prior(markers, dec)
    diag["labels_dropped_no_marker"] = dropped_labels(markers, dec)
    if not markers or prior.empty:
        reason[:] = REASON_ABLATION
        _write_diag(ctx, {**diag, "outcome": "no signature row survives the marker set"}, diagnostics_dir)
        return pd.DataFrame({"object_id": ids, "label": label, "score": score, "reason": reason})

    # values: finite rows only; negatives censored at 0, counted
    X = ctx.X[markers].to_numpy(float)
    finite = np.isfinite(X).all(axis=1)
    reason[~finite] = REASON_MISSING
    diag["n_missing_measurement"] = int((~finite).sum())
    neg = (X < 0) & finite[:, None]
    diag["n_censored_negative_by_marker"] = {m: int(neg[:, j].sum()) for j, m in enumerate(markers) if neg[:, j].any()}
    Xc = np.where(neg, 0.0, X)

    regions = region_ids(ctx)
    xy = ctx.context[["X_centroid", "Y_centroid"]].to_numpy(float)
    uniq = sorted(set(regions))
    scratch = (scratch_root or SCRATCH_ROOT) / f"{ctx.donor_id}__{ctx.scope}__{ctx.config_id.replace(':', '_')}__{os.getpid()}"
    (scratch / "in").mkdir(parents=True, exist_ok=True)
    prior.to_csv(scratch / "prior.csv", index=False)
    sent = []
    region_diag: dict = {}
    for r in uniq:
        sel = np.flatnonzero((regions == r) & finite)
        n_all = int((regions == r).sum())
        if len(sel) <= MIN_REGION_CELLS:
            reason[sel] = REASON_SMALL
            region_diag[r] = {"n_cells": n_all, "n_finite": int(len(sel)), "status": "too_small"}
            continue
        order = sel[np.argsort(ids[sel], kind="mergesort")]           # order-invariant input
        df = pd.DataFrame({"object_id": ids[order], "X": xy[order, 0], "Y": xy[order, 1]})
        for j, m in enumerate(markers):
            df[m] = Xc[order, j]
        df.to_csv(scratch / "in" / f"region_{r}.csv", index=False)
        sent.append(r)
        region_diag[r] = {"n_cells": n_all, "n_finite": int(len(sel)), "status": "sent"}

    if sent:
        cap = wall_clock_cap_s(len(sent))
        rres = _run_r(scratch, cap, ctx.seed)
        diag["r"] = {**rres, "wall_clock_cap_s": cap, "seed": int(ctx.seed)}
        env_path = scratch / "env.json"
        if env_path.exists():
            env = json.loads(env_path.read_text())
            _check_defaults(env)
            if env.get("celesta_sha") != CELESTA_COMMIT:
                raise RuntimeError(f"R half ran CELESTA {env.get('celesta_sha')}, pinned {CELESTA_COMMIT}")
            if int(env.get("rmixmod_seed", -1)) != int(ctx.seed):
                raise RuntimeError(f"R half seeded Rmixmod with {env.get('rmixmod_seed')}, ctx.seed is {ctx.seed}")
            diag["env"] = env
        elif not rres["timed_out"]:
            raise RuntimeError(f"celesta_run.R wrote no env.json (rc={rres['returncode']}):\n{rres['stderr_tail']}")
        pos = {oid: i for i, oid in enumerate(ids)}
        labels_ok = set(dec["labels"])
        for r in sent:
            sel = np.flatnonzero((regions == r) & finite)
            st_path = scratch / "out" / f"region_{r}.status.json"
            out_path = scratch / "out" / f"region_{r}.csv"
            if not st_path.exists():
                reason[sel] = REASON_CAP if rres["timed_out"] else REASON_ERROR
                region_diag[r]["status"] = "cap" if rres["timed_out"] else "no_status"
                continue
            st = json.loads(st_path.read_text())
            region_diag[r].update({k: v for k, v in st.items() if k != "region"})
            if st.get("status") == "too_small":
                reason[sel] = REASON_SMALL
                continue
            if st.get("status") == "activation_nan":
                reason[sel] = REASON_ACTIVATION
                continue
            if st.get("status") != "ok" or not out_path.exists():
                reason[sel] = REASON_ERROR
                continue
            res = pd.read_csv(out_path, dtype={"object_id": str})
            idx = np.array([pos[o] for o in res["object_id"]])
            lab = res["celesta_label"].astype(str).to_numpy()
            prob = res["celesta_prob"].to_numpy(float)
            unknown = lab == "Unknown"
            bad = ~unknown & ~np.isin(lab, list(labels_ok))
            if bad.any():
                raise RuntimeError(f"region {r}: CELESTA returned undeclared labels {sorted(set(lab[bad]))}")
            label[idx[~unknown]] = lab[~unknown]
            score[idx[~unknown]] = prob[~unknown]
            label[idx[unknown]] = records.AMBIGUOUS
            reason[idx[unknown]] = REASON_UNASSIGNED
            left = np.setdiff1d(sel, idx)                                # cells R did not return (should be none)
            reason[left] = REASON_ERROR
            region_diag[r]["n_not_returned"] = int(len(left))
        diag["celesta_stdout_log_lines"] = _count_lines(scratch / "celesta_stdout.log")
    unreasoned = (label == records.UNAVAILABLE) & (reason == "")
    if unreasoned.any():
        raise RuntimeError(f"{int(unreasoned.sum())} cells left Unavailable without a reason -- a driver bug")
    diag["regions"] = region_diag
    diag["label_counts"] = pd.Series(label).value_counts().to_dict()
    diag["reason_counts"] = pd.Series(reason[reason != ""]).value_counts().to_dict()
    _write_diag(ctx, diag, diagnostics_dir)
    if os.environ.get("BROADTYPE_KEEP_SCRATCH"):
        print(f"{ARM}: scratch kept at {scratch}")
    else:
        shutil.rmtree(scratch, ignore_errors=True)
    return pd.DataFrame({"object_id": ids, "label": label, "score": score, "reason": reason})


def _count_lines(p: pathlib.Path) -> int | None:
    if not p.exists():
        return None
    with p.open("rb") as fh:
        return sum(1 for _ in fh)


def _write_diag(ctx: arm_base.ArmContext, diag: dict, diagnostics_dir: pathlib.Path | None) -> pathlib.Path:
    """Beside the calls: what each region did, the censor counts, the R env.  Not a call table."""
    d = diagnostics_dir or records.out_dir(ARM, ctx.donor_id, ctx.scope)
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"diagnostics__{ctx.config_id.replace(':', '_')}.json"
    p.write_text(json.dumps(diag, indent=1, default=str) + "\n")
    return p


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=r_env(), default_params=DEFAULT_PARAMS))
