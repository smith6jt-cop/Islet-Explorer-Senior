"""The confound factorial: one variable per arm, everything else held.

The 2026-08-24 feasibility check found the model calling this tissue badly, but it could
not say WHAT was wrong, because four things were confounded at once: we hand the model a
52-channel panel against a training median of 34, on the densest 2048-px tile of the
cohort's least representative donor, with an uncorrected whole-cell mask, and with insulin
fed under the name `C-peptide`.  A verdict recorded on that configuration is a statement
about the configuration.

So: one variable per arm, over the frozen 8-donor / 8-region frame in ``regions.json``.

Two design points that are not cosmetic:

* **Concentric windows around a fixed anchor.**  Every arm scores the same 1024-px core.
  A window-size arm therefore changes the model's per-FOV normalization context and the
  neighbourhood its crops see, and NOT the cell population being scored.  The superseded
  ``--median-tile`` arm changed the tile and the tissue at once and so could not attribute
  its own result.
* **The halo arm folded into the window arm.**  With a 1024 core inside a 2048 window the
  baseline already carries a 512-px halo, far above the ``CROP_SIZE//2 = 16`` the model
  needs.  The interesting contrast is therefore "no halo at all" (window == core), which is
  the 1024 arm, not a separate 128-px one.

Scoring is always in the restricted label space with the unrestricted numbers beside it;
the restriction is a fixed property of the gate, never an arm, so it can never be switched
on where it happens to help.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time
import warnings
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from concordance_gate import GateInputs, run_gate               # noqa: E402
from label_space import load_label_space                        # noqa: E402
from regions import CORE_PX                                     # noqa: E402
from tiff_window import read_window                             # noqa: E402
from window_mask import Window, rasterize_windows               # noqa: E402

__all__ = ["Arm", "ARMS", "CHECKPOINTS", "run_arm", "plan_windows"]

MANIFEST = ROOT / "core/phenocycler/data/qupath_manifest.json"
REGIONS = HERE / "regions.json"
RUN_ID = "3bf4fc8e391d007766d1"
EXPRESSION = ROOT / f"core/phenocycler/data/runs/{RUN_ID}/expression"
OUT = ROOT / "data/development/deepcell_types_review"
ARCHIVE_MASKS = pathlib.Path("/data/codex_phenocycler/archive_2026-07-25/masks")

CHECKPOINTS = {
    "2026-06-15": pathlib.Path.home() / ".deepcell/models/deepcell-types_2026-06-15_resmlp.pt",
    "2026-06-23-ptft": pathlib.Path.home() / ".deepcell/models/deepcell-types_2026-06-23_resmlp_ptft.pt",
}

#: Declared mpp.  Native is 0.5077663810243286 um; declaring 0.5 makes the model's
#: resample a value-bit-exact identity and guarantees no cell can vanish to
#: ABSTENTION_LABEL on a mask rescale.  Cost: a recorded 1.553 % scale misstatement.
MPP = 0.5

#: A vocabulary marker this panel does not carry and this tissue has no use for -- the
#: trained-but-wrong-name control for the INS arm.  Not SMA: that name is already in the
#: panel, and dataset.py drops the duplicate rather than mislabelling it.
WRONG_NAME = "AQP2"


@dataclass(frozen=True)
class Arm:
    name: str
    hypothesis: str
    window: int = 2048          # concentric read around the anchor; core is always CORE_PX
    panel: str = "all"          # "all" | "type_only"
    mask: str = "cell"          # "cell" | "nucleus"
    ins: str | None = "C-peptide"   # None drops insulin entirely
    with_dapi: bool = False


ARMS = [
    Arm("A0-baseline", "reference: what the 2026-08-24 check ran, on the fixed frame"),
    Arm("A1-type-only", "panel composition is out of distribution: 52 channels, 20 of them "
        "state markers, against a training median of 34", panel="type_only"),
    Arm("A2-window-1024", "no halo at all, and the FOV the model normalizes over is the core "
        "itself", window=CORE_PX),
    Arm("A2b-window-1500", "match the training FOV in CELL COUNT rather than in pixels",
        window=1500),
    Arm("A3-window-4096", "bracket the window arms in the other direction", window=4096),
    Arm("A4-nucleus-mask", "DCT is segmentation-gated (dataset.py:197 zeroes every pixel "
        "outside the footprint), so uncorrected boundary spillover is read as signal",
        mask="nucleus"),
    Arm("A6a-ins-dropped", "the frozen text embedding is indexed by marker NAME; measure the "
        "cost of having no insulin channel at all", ins=None),
    Arm("A6b-ins-wrong-name", f"trained-but-wrong-name control: insulin pixels under "
        f"{WRONG_NAME!r}", ins=WRONG_NAME),
]


def _load(donor):
    man = json.loads(MANIFEST.read_text())
    im = next(i for i in man["images"] if i["donor_id"] == donor)
    regs = json.loads(REGIONS.read_text())["regions"][donor]
    return im, regs


def _window_for(region, size):
    """Concentric around the region's core centre, so the scored cells never move."""
    cy = region["y0"] + region["core_px"] // 2
    cx = region["x0"] + region["core_px"] // 2
    return Window(cy - size // 2, cx - size // 2, size)


def plan_windows(donor, arms, region_idx):
    """Every distinct window one GeoJSON pass must draw for this donor."""
    _, regs = _load(donor)
    want = {}
    for a in arms:
        for k in region_idx:
            w = _window_for(regs[k], a.window)
            want[(w.y0, w.x0, w.size, a.mask == "nucleus")] = w
    return want


def _channels(im, arm, registry_kind):
    """(plane index, registry marker, DCT name) for this arm's panel."""
    from deepcell_types.utils import resolve_supported_marker
    xw = json.loads((HERE / "deepcell_types_crosswalk.json").read_text())["entries"]
    sel = []
    for c in im["channel_map"]:
        nm, mk = c["channel_name"], c["marker"]
        if nm == "DAPI":
            if not arm.with_dapi:
                continue
            dct = "dsDNA"
        elif nm == "INS":
            if arm.ins is None:
                continue
            dct = arm.ins
        elif nm in xw:
            dct = xw[nm]["target"]
            if dct is None:
                continue
        else:
            dct = resolve_supported_marker(nm)
            if dct is None:
                continue
        if arm.panel == "type_only" and registry_kind.get(mk) != "type":
            continue
        sel.append((c["channel_index"], mk, dct))
    return sel


def run_arm(donor, arm, checkpoint, region_idx, masks_cache, *, device=None):
    """One arm on one donor's chosen regions.  Returns a list of per-region records."""
    import torch
    from deepcell_types import predict
    from deepcell_types.config import DCTConfig
    from deepcell_types.predict import ABSTENTION_LABEL
    from phenocycler.redsea import qptiff_channels

    reg_json = json.loads(REGIONS.read_text())["regions"][donor]
    im, _ = _load(donor)
    registry = json.loads((ROOT / "core/phenocycler/phenocycler/marker_registry.json").read_text())["markers"]
    registry = registry if isinstance(registry, list) else [dict(name=k, **v) for k, v in registry.items()]
    kind = {m["name"]: m["kind"] for m in registry}

    cfg = DCTConfig()
    if cfg.zarr_path is not None:                     # fail closed: DEEPCELL_TYPES_ZARR_PATH
        raise RuntimeError(f"DCTConfig().zarr_path is {cfg.zarr_path!r}; the registry would "
                           "differ from the one the crosswalk was gated against")
    ckpt = CHECKPOINTS[checkpoint]
    if ckpt.suffix != ".pt" or not ckpt.exists():
        raise FileNotFoundError(ckpt)

    sel = _channels(im, arm, kind)
    names = [d for _, _, d in sel]
    if len(names) != len(set(names)):
        dup = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"{arm.name}: duplicate DCT marker names {dup}; dataset.py would drop one")

    tf, pages, ome, (H, W), px_um = qptiff_channels(im["qptiff"]["path"], 0)
    for ch, _, _ in sel:
        assert ome[ch] == next(c["channel_name"] for c in im["channel_map"]
                               if c["channel_index"] == ch), f"OME name mismatch at plane {ch}"

    E = pd.read_parquet(EXPRESSION / f"donor_id={donor}" / "data.parquet")
    E["object_id"] = E["object_id"].astype(str)
    classes = tuple(k for k, _ in sorted(cfg.ct2idx.items(), key=lambda kv: int(kv[1])))
    ls = load_label_space(classes)
    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")

    records = []
    for k in region_idx:
        r = reg_json[k]
        w = _window_for(r, arm.window)
        mask, oids = masks_cache[(w.y0, w.x0, w.size, arm.mask == "nucleus")]

        raw = np.empty((len(sel), w.size, w.size), dtype=np.float32)
        for i, (ch, _, _) in enumerate(sel):
            raw[i] = read_window(pages[ch], w.y0, w.y1, w.x0, w.x1)

        t0 = time.time()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            res = predict(raw, mask.astype(np.int32), names, MPP, model_name=str(ckpt),
                          device=dev, return_probabilities=True)
        elapsed = time.time() - t0
        # predict() emits at least two unrelated warnings; classify rather than lump them.
        dropped = [str(x.message) for x in caught if "channel" in str(x.message).lower()]
        vanished = [str(x.message) for x in caught if "resampl" in str(x.message).lower()]
        if dropped:
            raise ValueError(f"{arm.name} {donor} r{k}: dataset dropped channels: {dropped}")
        if ABSTENTION_LABEL in res.cell_types:
            raise ValueError(f"{arm.name} {donor} r{k}: {ABSTENTION_LABEL} present with "
                             f"ct_abstention_k=0 -- cells vanished on resample: {vanished}")

        # score ONLY the common core, so window size never changes the cell population
        cy0, cx0 = r["y0"] - w.y0, r["x0"] - w.x0
        core = mask[cy0:cy0 + r["core_px"], cx0:cx0 + r["core_px"]]
        keep = set(np.unique(core[core > 0]).tolist())
        take = np.array([int(i) in keep for i in res.cell_indices])

        P = np.asarray(res.probabilities, float)[take]
        P = P / P.sum(1, keepdims=True)
        oid = np.array([str(oids[int(i) - 1]) for i in res.cell_indices])[take]

        J = pd.DataFrame(P, columns=list(classes)).assign(object_id=oid).merge(
            E, on="object_id", how="inner")
        Pj = J[list(classes)].to_numpy(float)
        Pj = Pj / Pj.sum(1, keepdims=True)
        rr = ls.restrict(Pj, classes)
        meas = {c: J[c].to_numpy(float) for c in J.columns
                if c not in classes and pd.api.types.is_numeric_dtype(J[c])}

        lab = f"{donor}/r{k}/{arm.name}/{checkpoint}"
        gi = GateInputs(Pj, classes, meas, label=lab)
        rec = {
            "donor": donor, "region": k, "arm": arm.name, "checkpoint": checkpoint,
            "core_origin": [r["y0"], r["x0"]], "window": [w.y0, w.x0, w.size],
            "mask": arm.mask, "panel": arm.panel, "n_channels": len(sel),
            "channels": names, "ins_as": arm.ins, "with_dapi": arm.with_dapi,
            "n_predicted": int(len(res.cell_types)), "n_scored": int(len(J)),
            "seconds": round(elapsed, 1), "mpp_declared": MPP, "native_um_px": px_um,
            "p_impossible_mass_mean": float(rr.impossible_mass.mean()),
            "p_impossible_outweighed": float((rr.ambiguous_reason == "mass_on_impossible_classes").mean()),
            "gate_restricted": {k2: v[0] for k2, v in run_gate(gi, restricted=True).criteria.items()},
            "gate_full": {k2: v[0] for k2, v in run_gate(gi, restricted=False).criteria.items()},
        }
        out = OUT / donor
        out.mkdir(parents=True, exist_ok=True)
        J.assign(restricted_call=rr.restricted_call, unrestricted_call=rr.unrestricted_call,
                 p_impossible_mass=rr.impossible_mass, ambiguous_reason=rr.ambiguous_reason
                 ).to_parquet(out / f"arm_{arm.name}_{checkpoint}_r{k}.parquet", index=False)
        records.append(rec)
        print(f"  {lab:<52} {len(J):>6,} cells  {elapsed:>6.1f}s  "
              f"impossible {rec['p_impossible_mass_mean']:.3f}", flush=True)
    return records


def build_masks(donor, arms, region_idx):
    """One GeoJSON pass per donor for every window every arm needs."""
    im, _ = _load(donor)
    want = plan_windows(donor, arms, region_idx)
    cell = sorted({k for k in want if not k[3]})
    nuc = sorted({k for k in want if k[3]})
    cache = {}
    if cell:
        wins = [want[k] for k in cell]
        masks, oids = rasterize_windows(im["cell_geojson"]["path"], wins)
        for k, m in zip(cell, masks):
            cache[k] = (m, oids)
    if nuc:
        wins = [want[k] for k in nuc]
        masks, oids, nucs = rasterize_windows(im["cell_geojson"]["path"], wins, with_nucleus=True)
        for k, n in zip(nuc, nucs):
            cache[k] = (n, oids)
    return cache
