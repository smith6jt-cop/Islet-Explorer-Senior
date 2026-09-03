#!/usr/bin/env python
"""Register the Xenium DAPI to the post-Xenium IF DAPI with VALIS.

Runs in the ISOLATED VALIS environment, not the pipeline one:

    ~/.venvs/valis/bin/python scripts/restain/valis_register.py --donors 6442

``valis-wsi`` pins scikit-image, pandas, pyvips, opencv and pillow to versions the pipeline
cannot use, so it lives in its own venv and talks to the rest of the lane through files:
``scripts/restain/export_dapi.py`` writes the two single-channel DAPI OME-TIFFs and the Xenium
cell centroids, and this writes back where those centroids land in the IF frame. It therefore
imports NOTHING from ``phenocycler`` and needs only numpy beyond VALIS itself.

The IF slide is the reference, because the islet and tissue annotations are drawn in its frame
and the lane's transform is defined Xenium -> IF.

Points are carried across with ``Slide.warp_xy_from_to``, which maps one slide's pixels
directly into another's and so never depends on VALIS's internal registered-canvas size. The
micron<->pixel conversions use the ``origin_um`` each side recorded in ``export.json``, because
the Xenium pyramid is decimated and the IF pyramid is box-averaged and neither is
``micron = pixel * um_per_px``.

**This is the SEED, not the answer.** ``run_refine.py`` matches the two nucleus fields
directly, and all it needs from here is a transform inside half the spacing between nuclei
(~4.5 um) so that nearest-neighbour matching is unambiguous. VALIS clears that comfortably.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# VALIS 1.2.0 matches features with LightGlue on the GPU and then calls `.numpy()` on a CUDA
# tensor without moving it to host, so `rigid_register` dies with "can't convert cuda:0 device
# type tensor to numpy". Hiding the GPU makes torch build CPU tensors and the same code path
# works. The cost is small: VALIS matches on its own downsampled images (a few hundred pixels
# on a side), not on the full section. Set VALIS_USE_GPU=1 to opt back in.
if os.environ.get("VALIS_USE_GPU", "0") != "1":
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

import numpy as np

# Keep this isolated environment, but load the repository's one authoritative donor gate.
# Registration directories are mutable working data, so scanning them without this check could
# revive a stale directory for a centrally excluded donor.
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from phenocycler.cohort import ensure_eligible_donors, filter_eligible_donors  # noqa: E402

XENIUM_SLIDE = "01_xenium_dapi"
IF_SLIDE = "02_if_dapi"


def _slide(registrar, stem: str):
    """VALIS keys slides by filename stem, but be forgiving about the exact key."""
    try:
        s = registrar.get_slide(stem)
        if s is not None:
            return s
    except Exception:
        pass
    for key, s in registrar.slide_dict.items():
        if Path(str(key)).stem.split(".")[0] == stem or stem in str(key):
            return s
    raise KeyError(f"no slide matching {stem!r}; have {sorted(registrar.slide_dict)}")


def register_one(donor_dir: Path, *, micro: bool = True, micro_rigid: bool = True) -> dict:
    from valis import registration

    donor_dir = Path(donor_dir)
    src = donor_dir / "slides"
    dst = donor_dir / "valis_out"
    meta = json.loads((donor_dir / "export.json").read_text())
    xen_um, if_um = float(meta["xenium_um_per_px"]), float(meta["if_um_per_px"])
    for key in ("xenium_origin_um", "if_origin_um"):
        if key not in meta:
            raise KeyError(f"{donor_dir/'export.json'} has no {key} -- re-run export_dapi.py")
    xen_o, if_o = float(meta["xenium_origin_um"]), float(meta["if_origin_um"])

    kwargs = dict(
        image_type="fluorescence",          # both planes are DAPI, not brightfield
        reference_img_f=f"{IF_SLIDE}.ome.tiff",
        align_to_reference=True,            # every slide to the reference, not slide-to-slide
    )
    if micro_rigid:
        # THE RIGID FIT IS OTHERWISE NEVER REFINED. `Valis.register` estimates M on the
        # processed thumbnails -- 512 px on a 20,880 px slide, about 40 um/px -- and the only
        # thing in 1.2.0 that re-fits M afterwards is `MicroRigidRegistrar`, which
        # `Valis.__init__` defaults to None (registration.py:1790) and which `register()` runs
        # between the rigid and non-rigid stages (registration.py:4743-4745). `register_micro`
        # is NOT that: every path in it ends at `bk_dxdy = current + new` and it never writes M,
        # so it cannot remove a translation however high its resolution.
        from valis.micro_rigid_registrar import MicroRigidRegistrar
        kwargs["micro_rigid_registrar_cls"] = MicroRigidRegistrar

    registrar = registration.Valis(str(src), str(dst), **kwargs)
    rigid, non_rigid, err_df = registrar.register()
    out = {"donor_dir": str(donor_dir), "xenium_um_per_px": xen_um, "if_um_per_px": if_um,
           "xenium_origin_um": xen_o, "if_origin_um": if_o,
           "reference": IF_SLIDE, "micro": bool(micro), "micro_rigid": bool(micro_rigid)}
    try:
        out["error_summary"] = json.loads(err_df.to_json())
    except Exception:
        out["error_summary"] = None

    if micro:
        # `reference_img_f` IS REQUIRED HERE. Without it `register_micro` forwards None to
        # `serial_non_rigid.register_images`, where `warp_tools.get_ref_img_idx` returns index 0
        # for two images (warp_tools.py:57-59) -- the XENIUM. That makes the Xenium the fixed
        # image, so it is handed a field of zeros (serial_non_rigid.py:778-781) and the field
        # computed for the IF is then skipped by the update loop, which skips the VALIS object's
        # own reference -- the IF (registration.py:4941-4943). The two references disagree, the
        # computed field is discarded, and the only surviving effect is that the OLD field is
        # resampled: measured, non_rigid_D went 20.497 -> 21.331 on donor 6436, i.e. slightly
        # worse. Passing it makes the pass do what its name says.
        # A sparse/blank peripheral tile can make VALIS return no deformation field.  That is
        # a failed optional micro pass, not a reason to discard the rigid seed: the downstream
        # nucleus fitter will score and reject an inadequate seed on its own measurements.
        try:
            moving = _slide(registrar, XENIUM_SLIDE)
            if getattr(moving, "bk_dxdy", None) is None:
                raise RuntimeError("VALIS produced no non-rigid field for the moving slide")
            registrar.register_micro(max_non_rigid_registration_dim_px=4096,
                                     reference_img_f=f"{IF_SLIDE}.ome.tiff",
                                     align_to_reference=True)
            out["micro_status"] = "completed"
        except Exception as exc:  # noqa: BLE001 - retain the independently scored rigid seed
            out["micro_status"] = f"skipped_after_error: {type(exc).__name__}: {exc}"
            print(f"[{donor_dir.name}] optional micro pass skipped: {exc}", file=sys.stderr,
                  flush=True)
    else:
        out["micro_status"] = "not_requested"

    xen, ifs = _slide(registrar, XENIUM_SLIDE), _slide(registrar, IF_SLIDE)
    xy_um = np.load(donor_dir / "xenium_centroids_um.npy")
    warped_if_um = np.empty_like(xy_um)
    CHUNK = 200_000
    for i in range(0, len(xy_um), CHUNK):
        # Pixel coordinate 0 is at `origin_um`, not at micron 0: the Xenium pyramid KEEPS
        # level-0 pixel f*k and the IF pyramid AVERAGES f*k..f*k+f-1, so their pixel centres sit
        # at different fractions of a level-0 pixel. Ignoring it leaves +0.39 um per axis in the
        # warped result -- a constant, which is the one kind of error this lane cannot absorb.
        px = (xy_um[i:i + CHUNK] - xen_o) / xen_um
        warped_if_um[i:i + CHUNK] = xen.warp_xy_from_to(px, ifs) * if_um + if_o
    np.save(donor_dir / "warped_if_um.npy", warped_if_um)

    d = warped_if_um - xy_um
    out["n_points"] = int(len(xy_um))
    out["displacement_um"] = {
        "median": float(np.median(np.hypot(d[:, 0], d[:, 1]))),
        "dx_mean": float(d[:, 0].mean()), "dy_mean": float(d[:, 1].mean()),
        "dx_sd": float(d[:, 0].std()), "dy_sd": float(d[:, 1].std()),
    }
    (donor_dir / "valis_result.json").write_text(json.dumps(out, indent=2))
    print(f"[{donor_dir.name}] warped {out['n_points']:,} centroids -> {donor_dir/'warped_if_um.npy'}",
          flush=True)
    print(f"           median |displacement| {out['displacement_um']['median']:.1f} um "
          f"(dx {out['displacement_um']['dx_mean']:+.1f}+/-{out['displacement_um']['dx_sd']:.1f}, "
          f"dy {out['displacement_um']['dy_mean']:+.1f}+/-{out['displacement_um']['dy_sd']:.1f})",
          flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path,
                    default=Path(__file__).resolve().parents[4]
                    / "data" / "development" / "restain" / "valis")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--no-micro", action="store_true")
    ap.add_argument("--no-micro-rigid", action="store_true",
                    help="skip MicroRigidRegistrar, leaving M at the thumbnail fit")
    args = ap.parse_args(argv)

    dirs = sorted(d for d in args.root.glob("donor_id=*") if d.is_dir())
    by_donor = {d.name.split("=", 1)[1]: d for d in dirs}
    if args.donors:
        want = set(ensure_eligible_donors(args.donors, context="VALIS registration"))
        absent = sorted(want - set(by_donor))
        if absent:
            raise FileNotFoundError(f"VALIS donor directories are absent: {absent}")
        dirs = [by_donor[donor] for donor in sorted(want, key=lambda value: int(value))]
    else:
        allowed = set(filter_eligible_donors(by_donor))
        dirs = [by_donor[donor] for donor in sorted(allowed, key=lambda value: int(value))]
    if not dirs:
        print(f"no donor directories under {args.root}", file=sys.stderr)
        return 1
    for d in dirs:
        register_one(d, micro=not args.no_micro, micro_rigid=not args.no_micro_rigid)
    try:
        from valis import registration
        registration.kill_jvm()
    except Exception:
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
