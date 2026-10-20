"""Arm E ``nimbus-rules``: Nimbus-Inference per-marker positivity on the frame windows, then
the ``rules_broad8`` logic on the resulting binary table.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/arms/nimbus_rules.py --donor 6436 --scope frame

Method (frozen 2026-09-05; any change after a result is seen is a NEW arm name)
--------------------------------------------------------------------------------
1. For each frame region ``i`` of the donor the FOV is the 2048-px window
   ``frame/regions.json -> donors[d].regions[i].window`` (level-0 pixels).  Every plane of
   ``ctx.markers`` -- the declared identity markers minus the ablated one, and nothing else --
   is read from the donor qptiff with ``common/tiff_window.read_window`` on the level-0 page
   whose OME ``<Biomarker>`` equals the manifest ``channel_name`` for that marker (asserted
   plane by plane, the way ``scripts/dct/arms.py`` does) and written as a ``uint16``
   ``<registry marker>.tiff`` in a folder per FOV.  The mask is the frame's own
   ``frame/masks/donor_id=<d>/w<i>.npy`` (int32, GLOBAL labels), written as a label TIFF.
2. ``nimbus_worker.py`` (subprocess in ``external/.venv-nimbus``) runs Nimbus-Inference
   ``@241281ef`` on those FOVs: normalization = Nimbus's own foreground 0.999 quantile averaged
   over all of the donor's FOVs; ``magnification=20`` against ``model_magnification=10`` (the
   0.5078 um/px window enters the network at 1024 px, ~1.016 um/px); test-time augmentation on;
   checkpoint ``V1.pt`` pinned by sha256.  It returns one confidence per cell per marker.
3. Positivity = confidence >= 0.5 -- Nimbus's own decision boundary for a two-class output,
   recorded here and applied nowhere else.  No quantile, no per-marker cut, no frequency.
4. Label rule = the ``rules_broad8`` row logic on the binary table: a cell is a CANDIDATE for
   label L iff any of L's declared markers present in ``ctx.markers`` is positive and none of
   L's declared negatives (``declarations.json``) that is present is positive.  Exactly one
   candidate -> that label, score = the max confidence over L's declared markers.  Several ->
   ``Ambiguous(multi_lineage)``.  None -> ``Ambiguous(no_positive_marker)`` when no marker is
   positive at all, ``Ambiguous(all_candidates_negated)`` when positives exist but every label
   they point at is blocked by one of its negatives (the task's single ``no_positive_marker``
   reason, split in two so the audit can tell "Nimbus saw nothing" from "Nimbus saw a
   contradiction"; the label is Ambiguous either way).
5. Cells are joined by MASK LABEL: ``object_ids[label - 1]`` with the donor's GeoJSON stream
   order (the same list ``frame.py`` used), matched on ``(object_id, region_idx)`` so a cell that
   also appears in a neighbouring window's halo is read from the window whose core holds it.
   A cell in ``ctx.ids`` with no Nimbus row -> ``Unavailable(not_in_fov)``; a row with a missing
   confidence -> ``Unavailable(marker_confidence_missing)``.  Never 0, never a default label.

Ablation is transparent: the withheld marker is simply not written as a channel, so it can
neither vote for its label nor act as a negative.  The arm is deterministic (``ctx.seed`` is
passed to the worker and recorded, but no step draws from it).

Scope: image arm, frame only.  The protocol runs image arms whole-slide only as survivors, so
``--scope donor`` is refused here rather than emitting 1.9 M ``Unavailable`` rows.

Scratch: ``data/development/broadtype/scratch/nimbus/job__<donor>__<config>__<pid>/`` (FOV
TIFFs, masks, Nimbus output), deleted when the run ends unless ``BROADTYPE_KEEP_SCRATCH=1``
(the shared ``arm_base.cli`` has no ``--keep-scratch`` flag).  ``scratch/nimbus/cache/`` keeps
the per-donor object-id list (keyed by the manifest's GeoJSON sha256, so a regenerated GeoJSON
cannot serve stale ids) because streaming the 3.5 GB GeoJSON costs ~1-2 min per run and the
screen runs this arm ~34 times per donor.

Beside every calls parquet the arm writes ``nimbus__<config>.json``: the worker's facts
(commit, checkpoint sha, device, normalization dict, OOM ladder, seconds), the plane map, the
per-marker positivity rates and the label counts -- descriptive, never a target.
"""
from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent          # scripts/broadtype/arms
BROADTYPE = HERE.parent
if str(BROADTYPE) not in sys.path:
    sys.path.insert(0, str(BROADTYPE))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import pins  # noqa: E402  (also puts scripts/common and core/multiplex5000 on sys.path)
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402
from tiff_window import read_window  # noqa: E402  (scripts/common)

ARM = "nimbus-rules"
WORKER = HERE / "nimbus_worker.py"
NIMBUS_PY = pins.REPO / "external" / ".venv-nimbus" / "bin" / "python"
SCRATCH = pins.OUT / "scratch" / "nimbus"
CACHE = SCRATCH / "cache"

#: External pins.  ``nimbus_worker.check`` refuses to run when the environment disagrees.
PIN = {
    "nimbus_repo": "smith6jt-cop/Nimbus-Inference",
    "nimbus_commit": "241281ef85ba34d9857cc03886cc45aedd0f0b34",
    "nimbus_version": "0.0.5",
    "torch": "2.5.1+cu121",
    "numpy": "1.26.4",
    "checkpoint": "V1.pt",
    "checkpoint_sha256": "0651a3852fbafc03be415a1ec80ef947cc1913a081cf0bc649301ade887a82a9",
    "checkpoint_source": "hf://JLrumberger/Nimbus-Inference/V1.pt",
}

#: Every parameter, recorded with every record.  None is a share, a prior or an ordering.
PARAMS = {
    "window_px": CRITERIA["window_px"],
    "pixel_size_um": pins.PIXEL_SIZE_UM,
    "dataset_magnification": 20,
    "model_magnification": 10,
    "effective_model_um_per_px": pins.PIXEL_SIZE_UM * 20 / 10,
    "normalization": "nimbus prepare_normalization_dict: foreground (>0) quantile per FOV and channel, mean over FOVs",
    "normalization_quantile": 0.999,
    "normalization_n_subset": "all FOVs of the donor",
    "normalization_clip": [0.0, 1.0],
    "test_time_aug": True,
    "batch_size": 4,
    "input_shape": [1024, 1024],
    "device": "auto",
    "mixed_precision": False,
    "compile_model": False,
    "oom_ladder": "halve batch_size to 1, then halve input_shape to 256 (recorded in the sidecar)",
    "positivity_threshold": 0.5,
    "positivity_rule": "confidence >= threshold (Nimbus's own two-class decision boundary)",
    "label_rule": ("rules_broad8 on the binary table: candidate iff any declared marker positive and no "
                   "declared negative positive; one candidate -> label (score = max declared confidence); "
                   "several -> Ambiguous(multi_lineage); none -> Ambiguous(no_positive_marker | all_candidates_negated)"),
    "join": "mask label -> object_ids[label-1] (donor GeoJSON stream order), matched on (object_id, region_idx)",
    "channel_file_names": "registry marker names",
    "mask_dtype": "int32 global labels (frame/masks)",
    "plane_dtype": "uint16 (qptiff level 0, unmodified)",
    "scope": "frame only (image arm)",
    "pins": PIN,
}
ENV = f"{NIMBUS_PY} (Nimbus-Inference {PIN['nimbus_version']} @ {PIN['nimbus_commit'][:8]}, torch {PIN['torch']}, numpy {PIN['numpy']})"


# --------------------------------------------------------------------------- the label rule (pure)
def calls_from_confidence(ids: np.ndarray, conf: pd.DataFrame, markers: list, labels: dict,
                          threshold: float = PARAMS["positivity_threshold"]) -> tuple:
    """The ``rules_broad8`` logic on a per-cell confidence table.  Returns ``(calls, info)``.

    ``conf`` is aligned row-for-row to ``ids`` (a reindex result: all-NaN rows are cells Nimbus
    never scored) with one column per marker in ``markers``.  ``labels`` is
    ``declarations.json["labels"]``.  Only ``markers`` are consulted: a column present in
    ``conf`` but absent from ``markers`` (the ablated marker, say) is invisible here, which is
    what makes the ablation transparent.
    """
    ids = np.asarray(ids).astype(str)
    n = len(ids)
    C = conf[list(markers)].to_numpy(dtype=float) if len(markers) else np.zeros((n, 0))
    finite = np.isfinite(C)
    in_table = finite.any(axis=1)
    present = finite.all(axis=1)
    pos = np.where(finite, C >= threshold, False)
    col = {m: j for j, m in enumerate(markers)}

    names = list(labels)
    cand = np.zeros((n, len(names)), dtype=bool)
    best = np.full((n, len(names)), np.nan)
    unreachable = {}
    for k, L in enumerate(names):
        decl = [col[m] for m in labels[L]["markers"] if m in col]
        negs = [col[m] for m in labels[L]["negatives"] if m in col]
        if not decl:
            unreachable[L] = "none of its declared markers is in the identity input"
            continue
        c = pos[:, decl].any(axis=1)
        if negs:
            c &= ~pos[:, negs].any(axis=1)
        cand[:, k] = c
        sub = np.where(finite[:, decl], C[:, decl], -np.inf).max(axis=1)     # -inf, not NaN: no all-NaN warning
        best[:, k] = np.where(np.isfinite(sub), sub, np.nan)

    n_cand = cand.sum(axis=1)
    label = np.full(n, records.AMBIGUOUS, dtype=object)
    reason = np.full(n, "", dtype=object)
    score = np.full(n, np.nan)
    any_pos = pos.any(axis=1)
    row_max = np.where(in_table, np.nanmax(np.where(finite, C, -np.inf), axis=1), np.nan) if C.shape[1] else np.full(n, np.nan)

    one = present & (n_cand == 1)
    k1 = cand[one].argmax(axis=1)
    label[one] = np.asarray(names, dtype=object)[k1]
    score[one] = best[one, k1]

    multi = present & (n_cand > 1)
    reason[multi] = "multi_lineage"
    if multi.any():
        score[multi] = np.where(cand[multi], np.nan_to_num(best[multi], nan=-np.inf), -np.inf).max(axis=1)

    none = present & (n_cand == 0)
    reason[none & any_pos] = "all_candidates_negated"
    reason[none & ~any_pos] = "no_positive_marker"
    score[none] = row_max[none]

    missing = in_table & ~present
    label[missing] = records.UNAVAILABLE
    reason[missing] = "marker_confidence_missing"

    absent = ~in_table
    label[absent] = records.UNAVAILABLE
    reason[absent] = "not_in_fov"

    calls = pd.DataFrame({"object_id": ids, "label": label.astype(str), "score": score, "reason": reason.astype(str)})
    info = {
        "threshold": threshold,
        "n_cells": int(n),
        "n_scored": int(present.sum()),
        "n_not_in_fov": int(absent.sum()),
        "n_confidence_missing": int(missing.sum()),
        "positivity_rate": {m: float(pos[present, j].mean()) if present.any() else float("nan") for m, j in col.items()},
        "n_positive_markers_per_cell": {str(k): int(v) for k, v in
                                        zip(*np.unique(pos[present].sum(axis=1), return_counts=True))},
        "label_counts": {str(k): int(v) for k, v in calls["label"].value_counts().items()},
        "reason_counts": {str(k): int(v) for k, v in calls.loc[calls["reason"] != "", "reason"].value_counts().items()},
        "unreachable_labels": unreachable,
    }
    return calls[list(records.CALL_COLUMNS)], info


# --------------------------------------------------------------------------- planes and FOVs
def plane_index(channel_map: list, ome_names: list, markers: list) -> dict:
    """marker -> level-0 plane index, asserting OME <Biomarker> == manifest channel_name per plane.

    The manifest ``channel_map`` (``channel_index`` / ``channel_name`` / ``marker``) is the
    only place the display form (``E-cadherin``) meets the registry name (``E_cadherin``); the
    qptiff carries the display form.  A plane whose OME name disagrees with the manifest is an
    error, not a warning: reading the wrong plane under the right name is the failure this
    guards.
    """
    by_marker = {c["marker"]: c for c in channel_map}
    out = {}
    for m in markers:
        if m not in by_marker:
            raise KeyError(f"marker {m!r} has no channel_map row on this image")
        c = by_marker[m]
        ch = int(c["channel_index"])
        if ch >= len(ome_names) or ome_names[ch] != c["channel_name"]:
            got = ome_names[ch] if ch < len(ome_names) else None
            raise RuntimeError(f"OME name mismatch at plane {ch}: qptiff {got!r} != manifest {c['channel_name']!r} ({m})")
        out[m] = ch
    return out


def write_fov(pages, planes: dict, window: dict, mask: np.ndarray, fov_dir: pathlib.Path,
              mask_path: pathlib.Path) -> dict:
    """One Nimbus FOV: ``<marker>.tiff`` per plane (uint16, unmodified values) + the label TIFF."""
    import tifffile
    y0, x0, size = int(window["y0"]), int(window["x0"]), int(window["size"])
    if mask.shape != (size, size):
        raise ValueError(f"mask {mask.shape} does not match window size {size}")
    fov_dir.mkdir(parents=True, exist_ok=True)
    mask_path.parent.mkdir(parents=True, exist_ok=True)
    for m, ch in planes.items():
        arr = read_window(pages[ch], y0, y0 + size, x0, x0 + size)
        if arr.shape != (size, size):
            raise ValueError(f"window {window} is clipped by the image edge: plane {ch} read {arr.shape}")
        tifffile.imwrite(fov_dir / f"{m}.tiff", np.ascontiguousarray(arr), photometric="minisblack", compression="zlib")
    tifffile.imwrite(mask_path, np.ascontiguousarray(mask.astype(np.int32)), photometric="minisblack", compression="zlib")
    return {"name": fov_dir.name, "path": str(fov_dir), "mask": str(mask_path),
            "window": {"y0": y0, "x0": x0, "size": size}, "n_labels": int(np.unique(mask[mask > 0]).size)}


def regions_for(donor: str) -> list:
    summary = json.loads((pins.OUT / "frame" / "regions.json").read_text())
    regs = summary["donors"][str(donor)]["regions"]
    for r in regs:
        if int(r["window"]["size"]) != CRITERIA["window_px"]:
            raise RuntimeError(f"frame window size {r['window']['size']} != criteria window_px {CRITERIA['window_px']}")
    return regs


def object_ids(donor: str) -> np.ndarray:
    """The donor's GeoJSON stream order: ``object_ids[label - 1]`` is the UUID behind a mask label.

    Cached under ``scratch/nimbus/cache/`` keyed by the manifest's GeoJSON sha256; the feature
    count is checked against the manifest on every load so a short or stale list fails loud.
    """
    im = pins.manifest_image(donor)
    sha = str(im["cell_geojson"]["sha256"])[:16]
    n_expected = int(im["cell_geometry"]["feature_count"])
    cache = CACHE / f"object_ids__donor_id={donor}__geojson_{sha}.parquet"
    if cache.exists():
        oids = pd.read_parquet(cache)["object_id"].astype(str).to_numpy(dtype=object)
    else:
        from phenocycler.artifacts import iter_qupath_geojson_features
        t0 = time.time()
        oids = np.asarray([str(ft.get("id")) for ft in iter_qupath_geojson_features(pins.geojson_path(donor))], dtype=object)
        CACHE.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"object_id": pd.array(oids, dtype="string")}).to_parquet(cache, index=False)
        print(f"{ARM} {donor}: streamed {len(oids):,} GeoJSON ids in {time.time() - t0:.0f}s -> {cache}", flush=True)
    if len(oids) != n_expected:
        raise RuntimeError(f"donor {donor}: {len(oids):,} GeoJSON ids != manifest feature_count {n_expected:,}; "
                           f"delete {cache} if the GeoJSON was regenerated")
    return oids


# --------------------------------------------------------------------------- the arm
def _job_dir(ctx: arm_base.ArmContext) -> pathlib.Path:
    stem = ctx.config_id.replace(":", "_")
    return SCRATCH / f"job__{ctx.donor_id}__{ctx.scope}__{stem}__{os.getpid()}"


def run_worker(job: dict, job_dir: pathlib.Path) -> tuple:
    """Write ``job.json``, run the worker in ``external/.venv-nimbus``, read the table and facts."""
    if not NIMBUS_PY.exists():
        raise FileNotFoundError(f"{NIMBUS_PY} is missing -- the Nimbus environment is not installed")
    job_path = job_dir / "job.json"
    job_path.write_text(json.dumps(job, indent=1))
    log = job_dir / "worker.log"
    with log.open("w") as fh:
        proc = subprocess.run([str(NIMBUS_PY), str(WORKER), str(job_path)], cwd=str(job_dir),
                              stdout=fh, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        tail = "".join(log.read_text().splitlines(keepends=True)[-40:])
        raise RuntimeError(f"nimbus_worker failed (exit {proc.returncode}); log tail:\n{tail}")
    out = pathlib.Path(job["output_dir"])
    table = pd.read_parquet(out / "cell_table.parquet")
    result = json.loads((out / "result.json").read_text())
    return table, result


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    if ctx.scope != "frame":
        raise ValueError(f"{ARM} is an image arm and runs on the frame; whole-slide only as a survivor "
                         f"(README 'Frame and cell set'); got --scope {ctx.scope}")
    if ctx.region_idx is None:
        raise RuntimeError("frame context carries no region_idx")
    from phenocycler.redsea import qptiff_channels

    dec = declarations.load()
    markers = list(ctx.markers)
    regs = regions_for(ctx.donor_id)
    im = pins.manifest_image(ctx.donor_id)
    job_dir = _job_dir(ctx)
    keep_scratch = os.environ.get("BROADTYPE_KEEP_SCRATCH", "") not in ("", "0")
    t0 = time.time()
    try:
        job_dir.mkdir(parents=True, exist_ok=False)
        # 1. planes -> FOV folders, masks -> label TIFFs (only the regions holding frame cells)
        tf, pages, ome_names, (H, W), px_um = qptiff_channels(pins.qptiff_path(ctx.donor_id), 0)
        if px_um is not None and abs(px_um - pins.PIXEL_SIZE_UM) > 1e-6:
            raise RuntimeError(f"qptiff pixel size {px_um} != pinned {pins.PIXEL_SIZE_UM}")
        planes = plane_index(im["channel_map"], ome_names, markers)
        used_regions = sorted(int(i) for i in np.unique(ctx.region_idx))
        fovs = []
        mdir = pins.OUT / "frame" / "masks" / f"donor_id={ctx.donor_id}"
        for i in used_regions:
            mask = np.load(mdir / f"w{i:02d}.npy")
            fovs.append(write_fov(pages, planes, regs[i]["window"], mask,
                                  job_dir / "fovs" / f"w{i:02d}", job_dir / "masks" / f"w{i:02d}_whole_cell.tiff"))
        tf.close()
        t_fov = time.time() - t0
        print(f"{ARM} {ctx.donor_id} {ctx.config_id}: wrote {len(fovs)} FOVs x {len(markers)} planes in {t_fov:.0f}s", flush=True)

        # 2. Nimbus in its own environment
        job = {
            "arm": ARM, "donor_id": ctx.donor_id, "config_id": ctx.config_id, "seed": int(ctx.seed),
            "fovs": fovs, "channels": markers, "output_dir": str(job_dir / "nimbus_out"),
            "dataset_magnification": PARAMS["dataset_magnification"], "model_magnification": PARAMS["model_magnification"],
            "normalization_quantile": PARAMS["normalization_quantile"], "test_time_aug": PARAMS["test_time_aug"],
            "batch_size": PARAMS["batch_size"], "input_shape": PARAMS["input_shape"], "device": PARAMS["device"],
            "pin": PIN,
        }
        table, result = run_worker(job, job_dir)
        t_nimbus = time.time() - t0 - t_fov
        print(f"{ARM} {ctx.donor_id} {ctx.config_id}: nimbus {len(table):,} rows on {result['device_used']} in {t_nimbus:.0f}s", flush=True)

        # 3. label -> object_id (GeoJSON stream order), matched on (object_id, region_idx)
        oids = object_ids(ctx.donor_id)
        lab = table["label"].to_numpy(dtype=np.int64)
        if lab.min() < 1 or lab.max() > len(oids):
            raise RuntimeError(f"mask label out of range [1, {len(oids)}]: {lab.min()}..{lab.max()}")
        table["object_id"] = oids[lab - 1].astype(str)
        table["region_idx"] = table["fov"].str.slice(1).astype(int)
        if table.duplicated(["object_id", "region_idx"]).any():
            raise RuntimeError("a (object_id, region_idx) pair appears twice in the Nimbus table")
        key = pd.MultiIndex.from_arrays([pd.Index(ctx.ids.astype(str)), pd.Index(ctx.region_idx.astype(int))])
        conf = table.set_index(["object_id", "region_idx"])[markers].reindex(key)

        # 4. positivity + the rule
        calls, info = calls_from_confidence(ctx.ids, conf, markers, dec["labels"], PARAMS["positivity_threshold"])
        info.update({
            "arm": ARM, "donor_id": ctx.donor_id, "config_id": ctx.config_id, "scope": ctx.scope,
            "identity_markers": markers, "ablated": ctx.ablated, "planes": planes,
            "regions_used": used_regions, "fovs": fovs, "image_shape": [int(H), int(W)], "pixel_size_um": px_um,
            "worker": result, "seconds_fov_write": round(t_fov, 1), "seconds_nimbus": round(t_nimbus, 1),
            "seconds_total": round(time.time() - t0, 1), "params": PARAMS,
        })
        side = records.out_dir(ARM, ctx.donor_id, ctx.scope)
        side.mkdir(parents=True, exist_ok=True)
        (side / f"nimbus__{ctx.config_id.replace(':', '_')}.json").write_text(json.dumps(info, indent=1, default=str))
        rates = "  ".join(f"{m}={info['positivity_rate'][m]:.3f}" for m in markers)
        print(f"{ARM} {ctx.donor_id} {ctx.config_id}: positivity  {rates}", flush=True)
        print(f"{ARM} {ctx.donor_id} {ctx.config_id}: labels {info['label_counts']}  reasons {info['reason_counts']}", flush=True)
        return calls
    finally:
        if job_dir.exists() and not keep_scratch:
            shutil.rmtree(job_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=ENV, default_params=PARAMS))
