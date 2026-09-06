"""Arm G ``dct-ptft-broad8-typeonly``: deepcell-types on the declared type markers only.

The closed DeepCellTypes lane (``scripts/dct``, FAIL at its own bar 2026-08-24) is not reopened.
This is a NEW arm under the broad-typing protocol: the same pretrained checkpoint
(``2026-06-23-ptft``), the closed factorial's ``A1-type-only`` channel shape -- exactly
``ctx.markers``, i.e. the declared identity (type) markers minus any ablated one -- and a NEW
label-space crosswalk (``dct_broad8_crosswalk.json``, authored before any output) that collapses
the model's 51 classes onto the 7 declared labels with a written reason per class.

Per frame region the 2048-px window's channel planes are read straight from the qptiff
(``common/tiff_window.read_window``; OME <Biomarker> name asserted against the manifest's
channel map, as the closed lane did), the frame's int32 label mask is handed to ``predict`` with
``mpp = 0.5`` declared (native 0.5078 um/px; the recorded 1.553 % misstatement keeps the model's
resample a value-bit-exact identity), and only cells whose label lies in the 1024-px core -- the
frozen frame cell set -- are scored.  The 51-way softmax is restricted with
``label_space.LabelSpace.restrict`` (renormalised over kept classes; a cell whose removed mass
outweighs its best kept class is Ambiguous) and summed per declared label; the label with the
largest summed probability is the call and that sum is the score.

Deterministic: ``STOCHASTIC = False`` (the screen uses 90 % subsample replicates).  Frame scope
only for now: at ``--scope donor`` every cell is ``Unavailable(whole_slide_not_run)`` -- a
whole-slide run is a survivor's step, recorded when it happens.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/arms/dct_ptft_broad8_typeonly.py --donor 6436 --scope frame
"""
from __future__ import annotations

import json
import pathlib
import sys
import time
import warnings

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
import label_space  # noqa: E402  (scripts/common)
from tiff_window import read_window  # noqa: E402  (scripts/common)

ARM = "dct-ptft-broad8-typeonly"
STOCHASTIC = False
CHECKPOINT = pathlib.Path.home() / ".deepcell" / "models" / "deepcell-types_2026-06-23_resmlp_ptft.pt"
CROSSWALK = HERE / "dct_broad8_crosswalk.json"
DCT_NAMES = pins.SCRIPTS / "dct" / "deepcell_types_crosswalk.json"   # OME name -> DCT vocabulary name
MPP_DECLARED = 0.5
R_NOT_PREDICTED = "not_predicted_in_window"
R_WHOLE_SLIDE = "whole_slide_not_run"

DEFAULT_PARAMS = {
    "method_version": "broadtype-dct-ptft-broad8-typeonly-2026-09-05",
    "checkpoint": str(CHECKPOINT),
    "deepcell_types_commit": "59d74a97711848dabc30fc01e242b6b1e4d3d3dc",
    "channels": "ctx.markers only (declared type markers minus the ablated one); the closed factorial's A1-type-only shape",
    "window_px": 2048, "core_px": 1024,
    "mpp_declared": MPP_DECLARED, "native_um_px": pins.PIXEL_SIZE_UM,
    "label_space": str(CROSSWALK.name),
    "call": "argmax of the restricted softmax summed per declared label; score = that sum",
    "ct_abstention_k": 0,
}


def dct_channel_names(donor: str, markers: list) -> list:
    """(plane index, registry marker, DCT name) for the declared markers this donor carries."""
    from deepcell_types.utils import resolve_supported_marker
    im = pins.manifest_image(donor)
    xw = json.loads(DCT_NAMES.read_text())["entries"]
    sel, dropped = [], []
    want = set(markers)
    for c in im["channel_map"]:
        nm, mk = c["channel_name"], c["marker"]
        if mk not in want:
            continue
        if nm in xw:
            dct = xw[nm]["target"]
        else:
            dct = resolve_supported_marker(nm)
        if dct is None:
            dropped.append(mk)
            continue
        sel.append((int(c["channel_index"]), mk, dct))
    names = [d for _, _, d in sel]
    if len(names) != len(set(names)):
        dup = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"duplicate DCT marker names {dup}; the dataset would drop one silently")
    return sel, dropped


def collapse(rr, ls, dec_labels: dict) -> tuple:
    """Restricted per-class probabilities -> per-label sums, call, score, reason."""
    groups = ls.broad_groups(restricted=True)
    labels = [L for L in dec_labels if L in groups]
    kept_idx = {c: i for i, c in enumerate(ls.kept)}
    S = np.zeros((len(rr.probabilities), len(labels)))
    for j, L in enumerate(labels):
        cols = [kept_idx[c] for c in groups[L]]
        S[:, j] = np.nansum(rr.probabilities[:, cols], axis=1)
    amb = rr.ambiguous_reason != ""
    top = S.argmax(axis=1)
    call = np.asarray(labels, dtype=object)[top]
    score = S[np.arange(len(S)), top]
    call[amb] = records.AMBIGUOUS
    reason = np.where(amb, rr.ambiguous_reason, "")
    return call, score, reason


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    n = len(ctx.ids)
    if ctx.scope != "frame":
        return pd.DataFrame({"object_id": ctx.ids.astype(str), "label": records.UNAVAILABLE,
                             "score": np.nan, "reason": R_WHOLE_SLIDE})
    import torch
    from deepcell_types import predict
    from deepcell_types.config import DCTConfig
    from deepcell_types.predict import ABSTENTION_LABEL
    from phenocycler.redsea import qptiff_channels

    cfg = DCTConfig()
    if cfg.zarr_path is not None:
        raise RuntimeError(f"DCTConfig().zarr_path is {cfg.zarr_path!r}; the vocabulary would differ from the crosswalk's")
    if not CHECKPOINT.exists():
        raise FileNotFoundError(CHECKPOINT)
    classes = tuple(k for k, _ in sorted(cfg.ct2idx.items(), key=lambda kv: int(kv[1])))
    ls = label_space.load_label_space(classes, path=CROSSWALK)
    dec = declarations.load()
    sel, dropped = dct_channel_names(ctx.donor_id, list(ctx.markers))
    names = [d for _, _, d in sel]
    if not names:
        return pd.DataFrame({"object_id": ctx.ids.astype(str), "label": records.UNAVAILABLE,
                             "score": np.nan, "reason": "ablation_removes_required_input"})
    im = pins.manifest_image(ctx.donor_id)
    tf, pages, ome, (H, W), px_um = qptiff_channels(pins.qptiff_path(ctx.donor_id), 0)
    for ch, _, _ in sel:
        want = next(c["channel_name"] for c in im["channel_map"] if c["channel_index"] == ch)
        assert ome[ch] == want, f"OME name mismatch at plane {ch}: {ome[ch]!r} != {want!r}"
    oids = arm_base.geojson_object_ids(ctx.donor_id)
    regions = arm_base.frame_regions(ctx.donor_id)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    pad = (DEFAULT_PARAMS["window_px"] - DEFAULT_PARAMS["core_px"]) // 2

    out_label = np.full(n, records.UNAVAILABLE, dtype=object)
    out_score = np.full(n, np.nan)
    out_reason = np.full(n, R_NOT_PREDICTED, dtype=object)
    pos = {o: i for i, o in enumerate(ctx.ids.astype(str))}
    region_of = ctx.region_idx
    seconds = []
    for k, r in enumerate(regions):
        if region_of is not None and not (region_of == k).any():
            continue
        w = r["window"]
        mask = arm_base.frame_mask(ctx.donor_id, k)
        raw = np.empty((len(sel), w["size"], w["size"]), dtype=np.float32)
        for i, (ch, _, _) in enumerate(sel):
            raw[i] = read_window(pages[ch], w["y0"], w["y0"] + w["size"], w["x0"], w["x0"] + w["size"])
        t0 = time.time()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            res = predict(raw, mask.astype(np.int32), names, MPP_DECLARED, model_name=str(CHECKPOINT),
                          device=dev, return_probabilities=True)
        seconds.append(time.time() - t0)
        dropped_ch = [str(x.message) for x in caught if "channel" in str(x.message).lower()]
        if dropped_ch:
            raise ValueError(f"{ARM} {ctx.donor_id} r{k}: the dataset dropped channels: {dropped_ch}")
        if ABSTENTION_LABEL in res.cell_types:
            raise ValueError(f"{ARM} {ctx.donor_id} r{k}: cells vanished on resample")
        core = mask[pad:pad + DEFAULT_PARAMS["core_px"], pad:pad + DEFAULT_PARAMS["core_px"]]
        keep = set(np.unique(core[core > 0]).tolist())
        take = np.array([int(i) in keep for i in res.cell_indices])
        if not take.any():
            continue
        Pm = np.asarray(res.probabilities, float)[take]
        Pm = Pm / Pm.sum(1, keepdims=True)
        rr = ls.restrict(Pm, classes)
        call, score, reason = collapse(rr, ls, dec["labels"])
        oid = [str(oids[int(i) - 1]) for i in np.asarray(res.cell_indices)[take]]
        for o, c, s_, rs in zip(oid, call, score, reason):
            j = pos.get(o)
            if j is None:
                continue                      # a GeoJSON-only feature with no expression row
            out_label[j], out_score[j], out_reason[j] = c, float(s_), rs
    ctx.extra["dct"] = {"n_windows": len(regions), "channels": names, "dropped_markers": dropped,
                        "predict_seconds": [round(s_, 1) for s_ in seconds], "device": dev,
                        "n_predicted": int((out_reason != R_NOT_PREDICTED).sum())}
    return pd.DataFrame({"object_id": ctx.ids.astype(str), "label": out_label, "score": out_score, "reason": out_reason})


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=sys.executable, default_params=DEFAULT_PARAMS))
