#!/usr/bin/env python
"""Joint rotation + translation search for the PhenoCycler->Xenium registration, per donor.

``pheno_rotation_figures.py`` asks whether a donor's registration is merely mis-oriented, and it
asks by sweeping a residual rotation about the warped Xenium centroid. That sweep cannot answer
the question it poses. Rotating about the moving centroid holds that centroid fixed, so the
sweep traces a one-parameter curve through a three-parameter space and no point on it carries a
translation. Both metrics it plots also saturate: islet Dice dilates each side by ~75 um and so
reads zero past ~100 um of residual, and islet RMSE is capped at 300 um on mutual pairs, so
under a wrong transform it collects chance pairs and converges to the null. A registration that
is wrong in angle *and* offset therefore produces exactly the flat curves that were read as
"these two sections do not share an islet field".

This script makes translation free at every angle -- see
:func:`phenocycler.integration.register.residual_vote` for the method -- and reports the peak
against the same donor's own landscape. It also scores the winning hypothesis on **ducts and
vessels**, which no part of the search consumes and which are called without reference to any
hormone, so they are not confounded by the donor's disease state the way a hormone-seeded islet
call is.

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/restain/pheno_residual_search.py --donors 6436 6442
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from phenocycler.config import load_config  # noqa: E402
from phenocycler.integration.register import _load_structures, residual_vote  # noqa: E402
from phenocycler.integration.transform import Transform, from_rigid  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"
REG_ROOT = CORE / "data" / "integration" / "registration"

#: What each structure kind is scored at. Islets and ducts are compact; a vessel centroid is a
#: looser thing, so it gets the same bound as a duct rather than the islet one.
MATCH_TOL_UM = {"islet": 100.0, "duct": 150.0, "vessel": 150.0}


def _xy(cfg, modality: str, donor: str, kind: str) -> np.ndarray:
    df = _load_structures(cfg, modality, donor, "panc", kind=kind)
    if df is None or len(df) == 0 or "x_um" not in df.columns:
        return np.zeros((0, 2))
    return df.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)


def _matched(tf: Transform, mv_xy: np.ndarray, fx_xy: np.ndarray, tol_um: float) -> tuple:
    """How many warped moving structures have a fixed structure within ``tol_um``, and the RMSE.

    One-way, not mutual: the question here is "does this hypothesis put anatomy on anatomy",
    and a mutual test would silently drop the many-to-one cases that a section thickness
    difference legitimately produces.
    """
    if len(mv_xy) == 0 or len(fx_xy) == 0:
        return 0, 0, float("nan")
    from scipy.spatial import cKDTree

    d, _i = cKDTree(fx_xy).query(tf.apply(mv_xy))
    hit = d <= tol_um
    rmse = float(np.sqrt(np.mean(d[hit] ** 2))) if hit.any() else float("nan")
    return int(hit.sum()), int(len(mv_xy)), rmse


def _chance(mv_xy: np.ndarray, fx_xy: np.ndarray, tol_um: float) -> float:
    """Expected matches if the fixed structures were scattered at random over their own extent.

    A Poisson bound, deliberately crude: it exists so a match count can be read against
    something, not to be a calibrated null. The landscape contrast is the calibrated statistic.
    """
    if len(mv_xy) == 0 or len(fx_xy) == 0:
        return float("nan")
    area = float(np.ptp(fx_xy[:, 0]) * np.ptp(fx_xy[:, 1]))
    if area <= 0:
        return float("nan")
    return len(mv_xy) * (1.0 - np.exp(-len(fx_xy) * np.pi * tol_um ** 2 / area))


def search_donor(cfg, donor: str, *, step_deg: float, span_um: float, bin_um: float,
                 kind: str = "islet") -> dict:
    base = replace(Transform.load(REG_ROOT / f"donor_id={donor}/roi=panc", stem="transform"),
                   displacement=None)
    mv, fx = _xy(cfg, "xenium", donor, kind), _xy(cfg, "phenocycler", donor, kind)
    warped = base.apply(mv)

    vote = residual_vote(warped, fx, angles_deg=np.arange(-180.0, 180.0, step_deg),
                         bin_um=bin_um, span_um=span_um)

    # The vote rotates about the warped moving centroid and then offsets, so the corrected
    # transform is exactly that residual composed onto the stored one.
    centre = tuple(float(v) for v in warped.mean(axis=0))
    residual = from_rigid(vote["angle_deg"], (vote["dx_um"], vote["dy_um"]), center_um=centre)
    corrected = residual.compose(base)

    evidence = {}
    for k, tol in MATCH_TOL_UM.items():
        m, f = _xy(cfg, "xenium", donor, k), _xy(cfg, "phenocycler", donor, k)
        evidence[k] = {
            "tol_um": tol,
            "chance": _chance(m, f, tol),
            "stored": _matched(base, m, f, tol),
            "corrected": _matched(corrected, m, f, tol),
        }

    return {"donor": donor, "vote": vote, "base": base, "corrected": corrected,
            "centre_um": centre, "evidence": evidence}


def _report(res: dict) -> None:
    v, donor = res["vote"], res["donor"]
    b, c = res["base"], res["corrected"]
    print(f"\ndonor {donor} — {v['n_moving']} Xenium vs {v['n_fixed']} PhenoCycler islets")
    print(f"  peak contrast {v['contrast']:.2f}x  (peak {v['peak']:.2f}, p95 {v['p95']:.2f})"
          f"  at residual {v['angle_deg']:+.0f} deg"
          f"{' MIRRORED' if v['mirror'] else ''}, ({v['dx_um']:+.0f}, {v['dy_um']:+.0f}) um")
    print(f"  stored     rot {b.rotation_deg:+7.2f}  scale {b.scale:.4f}  "
          f"t ({b.affine[0, 2]:+8.1f}, {b.affine[1, 2]:+8.1f})")
    print(f"  corrected  rot {c.rotation_deg:+7.2f}  scale {c.scale:.4f}  "
          f"t ({c.affine[0, 2]:+8.1f}, {c.affine[1, 2]:+8.1f})")
    print(f"  {'kind':8s} {'tol':>5s} {'chance':>7s} {'stored':>16s} {'corrected':>16s}")
    for k, e in res["evidence"].items():
        sh, sn, sr = e["stored"]
        ch, cn, cr = e["corrected"]
        if sn == 0:
            print(f"  {k:8s} {e['tol_um']:5.0f} {'--':>7s}   none on one side")
            continue
        print(f"  {k:8s} {e['tol_um']:5.0f} {e['chance']:7.1f} "
              f"{sh:6d}/{sn:<4d} {sr:6.1f}um {ch:6d}/{cn:<4d} {cr:6.1f}um")


def figure(results: list[dict], out: Path) -> None:
    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(6.4 * n, 9.6), constrained_layout=True,
                             squeeze=False)
    for col, res in enumerate(results):
        v = res["vote"]
        ax = axes[0][col]
        for mi, (label, style) in enumerate((("no mirror", "-"), ("mirrored", "--"))):
            if v["landscape"].shape[1] > mi:
                ax.plot(v["angles_deg"], v["landscape"][:, mi], style, linewidth=1.3,
                        label=label, color="#1f77b4" if mi == 0 else "#7f7f7f")
        ax.axhline(v["p95"], color="#2ca02c", linestyle=":", linewidth=1.4,
                   label=f"p95 of the landscape = {v['p95']:.2f}")
        ax.axvline(0.0, color="#d62728", linestyle="--", linewidth=1.2,
                   label="the stored transform")
        ax.plot([v["angle_deg"]], [v["peak"]], marker="o", markersize=7, color="#ff7f0e",
                linestyle="none", label=f"peak {v['contrast']:.1f}x p95")
        ax.set_xlabel("residual rotation applied on top of the stored transform (deg)",
                      fontsize=9)
        ax.set_ylabel("best vote at this angle\n(translation free, "
                      "normalised to the landscape median)", fontsize=9)
        ax.set_title(f"donor {res['donor']} — {v['n_moving']} Xenium vs {v['n_fixed']} "
                     f"PhenoCycler islets\ntranslation is re-optimised at every angle",
                     fontsize=10)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)

        axm = axes[1][col]
        cen = v["bin_centres_um"]
        # Linear from zero: the floor of this map is a real measured vote count, so a
        # percentile black point would hide how flat the background actually is.
        im = axm.imshow(v["best_map"].T, origin="lower", cmap="magma", vmin=0.0,
                        vmax=float(v["best_map"].max()),
                        extent=(cen[0], cen[-1], cen[0], cen[-1]), interpolation="nearest")
        axm.plot([0], [0], marker="+", markersize=14, markeredgewidth=2, color="#00e5ff",
                 linestyle="none", label="zero residual (the stored transform)")
        axm.plot([v["dx_um"]], [v["dy_um"]], marker="o", markersize=9, markerfacecolor="none",
                 markeredgewidth=2, color="#7CFC00", linestyle="none", label="peak")
        axm.set_xlabel("residual x offset (µm)", fontsize=9)
        axm.set_ylabel("residual y offset (µm)", fontsize=9)
        axm.set_title(f"offset votes at the winning angle "
                      f"({v['angle_deg']:+.0f}°{', mirrored' if v['mirror'] else ''})",
                      fontsize=10)
        axm.legend(fontsize=8, loc="upper right", framealpha=0.85)
        axm.tick_params(labelsize=8)
        fig.colorbar(im, ax=axm, fraction=0.046, label="votes / landscape median")

    fig.suptitle(
        "Is the registration in the right place at all?  Rotation AND translation, searched "
        "jointly.\nA rotation-only sweep holds the warped moving centroid fixed, so it cannot "
        "see an offset — which is why it read flat.",
        fontsize=12)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\n[figure] {out}", flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Joint rotation+translation residual search")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--step-deg", type=float, default=2.0)
    ap.add_argument("--span-um", type=float, default=2500.0)
    ap.add_argument("--bin-um", type=float, default=50.0)
    ap.add_argument("--kind", default="islet", help="structure kind driving the vote")
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args(argv)

    cfg = load_config(CORE / "config.ini")
    results = [search_donor(cfg, d, step_deg=args.step_deg, span_um=args.span_um,
                            bin_um=args.bin_um, kind=args.kind) for d in args.donors]
    for res in results:
        _report(res)

    figure(results, OUT_DIR / f"pheno_residual_search_{'_'.join(args.donors)}.png")

    if args.json_out is not None:
        payload = {}
        for res in results:
            v = res["vote"]
            payload[res["donor"]] = {
                "contrast": v["contrast"], "peak": v["peak"], "p95": v["p95"],
                "angle_deg": v["angle_deg"], "mirror": v["mirror"],
                "dx_um": v["dx_um"], "dy_um": v["dy_um"],
                "bin_um": v["bin_um"], "span_um": v["span_um"], "smooth_um": v["smooth_um"],
                "stored_transform": res["base"].to_dict(),
                "corrected_transform": res["corrected"].to_dict(),
                "evidence": {k: {"tol_um": e["tol_um"], "chance": e["chance"],
                                 "stored": list(e["stored"]),
                                 "corrected": list(e["corrected"])}
                             for k, e in res["evidence"].items()},
            }
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2))
        print(f"[json] {args.json_out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
