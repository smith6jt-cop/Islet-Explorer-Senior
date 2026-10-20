"""Nimbus-Inference worker for arm E ``nimbus-rules`` -- runs INSIDE ``external/.venv-nimbus``.

    external/.venv-nimbus/bin/python scripts/broadtype/arms/nimbus_worker.py <job.json>
    external/.venv-nimbus/bin/python scripts/broadtype/arms/nimbus_worker.py --check

This file is the only code that imports ``nimbus_inference`` or ``torch``.  The driver
(``nimbus_rules.py``, ``.venv``) writes the FOV folders and a ``job.json``, calls this script as
a subprocess, and reads back ``cell_table.parquet`` + ``result.json``.  Keeping the two
interpreters apart is not cosmetic: the Nimbus environment pins ``numpy 1.26`` and
``torch 2.5.1+cu121`` against the repo ``.venv``'s stack, and the arm contract (README
"The arm contract", item 5) says such a method is driven as a subprocess exchanging files under
``data/development/broadtype/scratch/<arm>/``.

What Nimbus is given, and why each choice is what it is:

* **FOV layout** -- ``utils.MultiplexDataset`` with a *folder per FOV* holding one single-channel
  TIFF per marker named ``<marker>.tiff`` (``get_channels`` lists ``fov_paths[0]``;
  ``get_channel_single`` opens ``<fov>/<channel><suffix>``).  The driver names the files by the
  REGISTRY marker name, so the cell table comes back keyed the way ``declarations.json`` is.
  The segmentation is supplied through ``segmentation_naming_convention``, a callable from the
  FOV path to the label TIFF; ``get_segmentation`` reads it with ``imageio`` and casts to
  ``uint32``, so the frame's GLOBAL int32 labels (up to ~1.94 M on donor 6436) survive intact.
* **Pins verified before anything runs** -- the installed commit (from the wheel's
  ``direct_url.json``), the torch / numpy versions, and the sha256 of the checkpoint file.  The
  checkpoint is loaded by NAME from the package ``assets/`` directory (``checkpoint="V1.pt"``),
  never as ``"latest"``: ``load_latest_checkpoint`` consults the Hugging Face Hub on every call
  and would silently upgrade the model the day a ``V2.pt`` appears.  A missing or mismatched file
  is an error with the download command in the message, never a silent fetch.
* **Normalization** -- Nimbus's own ``prepare_normalization_dict``: per FOV and channel the
  ``quantile`` of the FOREGROUND (> 0) pixels, averaged over the FOVs (``n_subset`` = all of the
  donor's FOVs, so no random subset is drawn), then ``get_channel_normalized`` divides and clips
  to ``[0, 1]``.  That is a per-channel intensity scale learned from the images themselves; it
  says nothing about how many cells should be positive.
* **Magnification** -- ``MultiplexDataset(magnification=20)`` and
  ``Nimbus(model_magnification=10)``: ``predict_fovs`` resizes every input by 10/20, so the
  2048-px, 0.5078 um/px window enters the network as 1024 px at ~1.016 um/px, which is the
  half-resolution the released weights were trained at (the checkpoint's own name says
  ``halfres``).  Both numbers are recorded in ``result.json``.
* **Test-time augmentation on** (Nimbus's default: 4 rotations + 2 flips, averaged) --
  ``predict_fovs`` runs each augmentation as its own forward pass at ``input_shape``; the
  ``batch_size`` only matters once the input exceeds ``input_shape`` and Nimbus tiles it.
* **OOM ladder** -- on ``torch.cuda.OutOfMemoryError`` the batch size halves to 1, then
  ``input_shape`` halves (Nimbus then tiles with ``valid`` padding and stitches).  Every attempt
  is recorded; the values that finally ran are in ``result.json``.

The worker returns CONFIDENCES (Nimbus's per-cell mean of its pixel-wise output); the 0.5
decision boundary and the label rule live in the driver, so the threshold is one recorded
number in one place.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import random
import sys
import time

#: The Nimbus-Inference fork this arm was built against, installed from git into
#: ``external/.venv-nimbus``.  The driver carries the same three constants; ``check`` fails if
#: the environment disagrees with them.
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


def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def installed_commit() -> str | None:
    """The git commit pip recorded for the installed Nimbus-Inference wheel, or None."""
    from importlib.metadata import PackageNotFoundError, distribution
    try:
        txt = distribution("Nimbus-Inference").read_text("direct_url.json")
    except PackageNotFoundError:
        return None
    if not txt:
        return None
    return json.loads(txt).get("vcs_info", {}).get("commit_id")


def assets_dir() -> pathlib.Path:
    import nimbus_inference
    return pathlib.Path(nimbus_inference.__file__).resolve().parent / "assets"


def check(pin: dict = PIN) -> dict:
    """Verify every pin; raise with the fix in the message when one is off.  Returns the facts."""
    import numpy
    import torch
    import nimbus_inference
    facts = {
        "python": sys.executable,
        "nimbus_version": getattr(nimbus_inference, "__version__", None),
        "nimbus_commit": installed_commit(),
        "torch": torch.__version__,
        "numpy": numpy.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
    }
    bad = []
    if facts["nimbus_commit"] != pin["nimbus_commit"]:
        bad.append(f"Nimbus-Inference commit {facts['nimbus_commit']!r} != pinned {pin['nimbus_commit']!r}")
    if facts["torch"] != pin["torch"]:
        bad.append(f"torch {facts['torch']!r} != pinned {pin['torch']!r}")
    if facts["numpy"] != pin["numpy"]:
        bad.append(f"numpy {facts['numpy']!r} != pinned {pin['numpy']!r}")
    ckpt = assets_dir() / pin["checkpoint"]
    if not ckpt.exists():
        bad.append(f"checkpoint {ckpt} is missing; fetch it ONCE with "
                   f"hf_hub_download(repo_id='JLrumberger/Nimbus-Inference', filename={pin['checkpoint']!r}, "
                   f"local_dir={str(assets_dir())!r}) and re-check its sha256 against the pin")
    else:
        sha = _sha256(ckpt)
        facts["checkpoint_path"] = str(ckpt)
        facts["checkpoint_sha256"] = sha
        if sha != pin["checkpoint_sha256"]:
            bad.append(f"checkpoint sha256 {sha} != pinned {pin['checkpoint_sha256']}")
    if bad:
        raise RuntimeError("nimbus_worker pins do not match the environment:\n  " + "\n  ".join(bad))
    return facts


def run_job(job: dict) -> dict:
    """Predict every FOV in ``job``; write ``cell_table.parquet`` and ``result.json`` beside it."""
    import numpy as np
    import pandas as pd
    import torch
    from nimbus_inference.nimbus import Nimbus
    from nimbus_inference.utils import MultiplexDataset

    facts = check(job.get("pin", PIN))
    out = pathlib.Path(job["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    fovs = job["fovs"]
    fov_paths = [f["path"] for f in fovs]
    mask_of = {os.path.basename(f["path"]): f["mask"] for f in fovs}

    def segmentation_naming_convention(fov_path: str) -> str:
        # Nimbus calls this with the FOV path and reads whatever path comes back.
        return mask_of[os.path.basename(fov_path)]

    seed = int(job.get("seed", 0))
    random.seed(seed)           # prepare_normalization_dict shuffles fov_paths before the (full) subset
    torch.manual_seed(seed)     # inference has no stochastic op; recorded for completeness

    t0 = time.time()
    ds = MultiplexDataset(
        fov_paths, segmentation_naming_convention=segmentation_naming_convention,
        include_channels=list(job["channels"]), suffix=".tiff", silent=True,
        magnification=int(job["dataset_magnification"]), output_dir=str(out))
    missing = sorted(set(job["channels"]) - set(ds.channels))
    if missing:
        raise RuntimeError(f"FOV folders lack channel files for {missing}")
    ds.prepare_normalization_dict(
        quantile=float(job["normalization_quantile"]), n_subset=len(fov_paths),
        multiprocessing=False, overwrite=True)
    t_norm = time.time() - t0

    attempts = []
    batch_size = int(job["batch_size"])
    input_shape = [int(v) for v in job["input_shape"]]
    while True:
        t1 = time.time()
        try:
            nimbus = Nimbus(
                ds, str(out), save_predictions=False, batch_size=batch_size,
                test_time_aug=bool(job["test_time_aug"]),
                model_magnification=int(job["model_magnification"]), input_shape=input_shape,
                device=job.get("device", "auto"), checkpoint=job.get("pin", PIN)["checkpoint"],
                compile_model=False, mixed_precision=False)
            table = nimbus.predict_fovs()
            attempts.append({"batch_size": batch_size, "input_shape": list(input_shape), "ok": True,
                             "seconds": round(time.time() - t1, 1)})
            break
        except torch.cuda.OutOfMemoryError as e:  # noqa: PERF203
            attempts.append({"batch_size": batch_size, "input_shape": list(input_shape), "ok": False,
                             "error": str(e)[:200]})
            torch.cuda.empty_cache()
            if batch_size > 1:
                batch_size //= 2
            elif input_shape[0] > 256:
                input_shape = [s // 2 for s in input_shape]
            else:
                raise
            print(f"nimbus_worker: CUDA OOM; retrying with batch_size={batch_size} input_shape={input_shape}",
                  flush=True)

    # Nimbus's table: label, fov, centroid_x, centroid_y, area, <channel>...  Keep label/fov/area
    # and the confidences; centroids are dropped on purpose -- object_id is the only join key and
    # the driver joins on the mask label, never on a centroid.
    keep = ["label", "fov", "area", *job["channels"]]
    table = table[keep].copy()
    table["label"] = table["label"].astype(np.int64)
    table.to_parquet(out / "cell_table.parquet", index=False)

    result = {
        **facts,
        "device_used": str(nimbus.device),
        "checkpoint_loaded": str(nimbus.checkpoint_path),
        "dataset_magnification": int(job["dataset_magnification"]),
        "model_magnification": int(job["model_magnification"]),
        "normalization_quantile": float(job["normalization_quantile"]),
        "normalization_n_fovs": len(fov_paths),
        "normalization_dict": {k: float(v) for k, v in ds.normalization_dict.items()},
        "normalization_clip": [0.0, 1.0],      # get_channel_normalized clips to [0, 1]
        "test_time_aug": bool(job["test_time_aug"]),
        "attempts": attempts,
        "batch_size_used": batch_size,
        "input_shape_used": list(input_shape),
        "n_fovs": len(fov_paths),
        "n_channels": len(job["channels"]),
        "n_rows": int(len(table)),
        "seconds_normalization": round(t_norm, 1),
        "seconds_total": round(time.time() - t0, 1),
        "seed": seed,
    }
    (out / "result.json").write_text(json.dumps(result, indent=1))
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Nimbus-Inference worker (runs in external/.venv-nimbus)")
    ap.add_argument("job", nargs="?", help="path to job.json written by nimbus_rules.py")
    ap.add_argument("--check", action="store_true", help="verify the pins and print the facts as JSON")
    a = ap.parse_args(argv)
    if a.check or not a.job:
        print(json.dumps(check(), indent=1))
        return 0
    job = json.loads(pathlib.Path(a.job).read_text())
    res = run_job(job)
    print(f"nimbus_worker: {res['n_rows']:,} rows over {res['n_fovs']} FOVs x {res['n_channels']} channels "
          f"on {res['device_used']} in {res['seconds_total']:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
