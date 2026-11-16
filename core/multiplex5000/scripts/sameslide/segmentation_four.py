#!/usr/bin/env python3
"""InstanSeg against the Xenium Analyzer, nucleus and cell, on one registered same-slide section.

The two platforms segmented the same physical section independently: InstanSeg
(``fluorescence_nuclei_and_cells``, QuPath) on the PhenoCycler qptiff, and the Xenium Analyzer on
its own morphology image. The join (``join_four.py``) carries protein onto Xenium cells by area
overlap and needs no agreement between the two; this stage measures the agreement itself, with
``segcompare.compare`` -- geometry only, both directions, neither side taken as truth.

* **Frame.** The PhenoCycler frame, as everywhere in this lane: Xenium nuclei and cells are
  warped forward through the validated registration transform, displacement field included.
* **Objects.** Every object each platform segmented, not the analysis-eligible subset -- QC is a
  downstream choice and would make the comparison depend on it. InstanSeg nuclei and cells are
  the preparation's ``nucleus_wkb`` and ``wkb``; Xenium's are ``nucleus_boundaries.parquet``
  (split on ``label_id``, so a cell's second nucleus is its own object) and
  ``cell_boundaries.parquet``.
* **Area.** Counted only where both platforms segmented (``segcompare.shared_area``); the area
  only one of them segmented is reported, not dropped silently.
* **Strata.** Tissue for every object: an InstanSeg object's from the preparation's routing
  (``panc`` -> pancreas, ``pln_1`` -> pLN), a Xenium object's from the nearest InstanSeg cell. For
  Xenium cells, also the Analyzer's own ``segmentation_method`` where the bundle records it --
  boundary stain, interior stain, or nuclear expansion -- because an expanded nucleus is not a
  measured cell outline and should not be read as one.
* **Registration ceiling.** The recheck's median nucleus-pair distance is used as the shift in
  ``shift_ceiling``: what two identical segmentations would score at that misregistration. That
  distance includes segmentation disagreement, so the ceiling is conservative (lower than a
  registration-only ceiling).
* **Hierarchy.** Of nuclei matched across platforms, how many sit in cells matched to each other.

Receipts as ``join_four.py``: registration, its recheck against THIS preparation, and the
preparation itself; an upstream that is scientifically unavailable makes this stage unavailable
for the section, recorded as current.

Outputs (``<section>/segmentation/attempts/<id>/``): ``segmentation.json`` (summary),
``thresholds.csv``, per level ``edges_<level>.parquet``, ``objects_instanseg_<level>.parquet``,
``objects_xenium_<level>.parquet``, and ``segmentation.png``.

    python segmentation_four.py --manifest $M --section <section>
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback
from uuid import uuid4

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from phenocycler.integration.boundaries import _read_boundary_rings  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402
from build_joined import xenium_polygons  # noqa: E402
from join_four import TISSUE, current  # noqa: E402
from manifest import load_manifest, get_section  # noqa: E402
from workflow import DEFAULT_OUT, section_root, identity, write_json, now, valid_receipt  # noqa: E402
import segcompare as sc  # noqa: E402

LEVELS = ("nucleus", "cell")
PLURAL = {"nucleus": "nuclei", "cell": "cells"}
REASONS = {"matched", "boundary", "split", "merged", "no_overlap"}


def warp(polys, tf):
    """Xenium geometry forward into the PhenoCycler frame, vertex by vertex."""
    return shapely.transform(polys, lambda c: tf.apply(c))


def instanseg_objects(prep: Path):
    """InstanSeg cells and nuclei (microns, PhenoCycler frame), with tissue routing.

    Returns cells (object_id, geometry, tissue) and nuclei (object_id, geometry, parent position
    into cells). A cell without ``nucleus_wkb`` has no nucleus object; it stays a cell.
    """
    geo = pd.read_parquet(prep / "geometry.parquet", columns=["object_id", "wkb", "nucleus_wkb"])
    geo = geo.assign(object_id=geo.object_id.astype(str))
    if geo.object_id.duplicated().any():
        raise ValueError("duplicate object_id in preparation geometry")
    geo = geo[geo.wkb.notna()].reset_index(drop=True)
    tissue = pd.Series("unrouted", index=geo.object_id, dtype=object)
    for roi, name in TISSUE.items():
        path = prep / roi / "cells.parquet"
        if not path.exists():
            continue
        ids = pq.read_table(path, columns=["object_id"]).column("object_id").to_numpy(
            zero_copy_only=False).astype(str)
        hit = tissue.index.isin(ids)
        if (tissue[hit] != "unrouted").any():
            raise ValueError(f"an InstanSeg object is routed to two tissues ({roi})")
        tissue[hit] = name
    cells = pd.DataFrame({"object_id": geo.object_id.to_numpy(),
                          "geometry": shapely.from_wkb(geo.wkb.to_numpy()),
                          "tissue": tissue.to_numpy()})
    has_nuc = geo.nucleus_wkb.notna().to_numpy()
    nuclei = pd.DataFrame({"object_id": geo.object_id.to_numpy()[has_nuc],
                           "geometry": shapely.from_wkb(geo.nucleus_wkb.to_numpy()[has_nuc]),
                           "parent": np.flatnonzero(has_nuc)})
    return cells, nuclei


def xenium_objects(bundle: Path, tf):
    """Xenium cells and nuclei warped into the PhenoCycler frame, with each nucleus's cell."""
    cids, cpolys = xenium_polygons(bundle)
    cids = np.asarray(cids).astype(str)
    if len(set(cids)) != len(cids):
        raise ValueError("duplicate Xenium cell_id in cell_boundaries")
    cells = pd.DataFrame({"cell_id": cids, "geometry": warp(np.asarray(cpolys, dtype=object), tf)})
    meta_cols = [c for c in ("cell_id", "segmentation_method", "nucleus_count")
                 if c in pq.ParquetFile(bundle / "cells.parquet").schema_arrow.names]
    meta = pd.read_parquet(bundle / "cells.parquet", columns=meta_cols)
    meta = meta.assign(cell_id=meta.cell_id.astype(str)).set_index("cell_id")
    cells["segmentation_method"] = (
        meta["segmentation_method"].reindex(cells.cell_id).fillna("not_recorded").astype(str).to_numpy()
        if "segmentation_method" in meta else "not_recorded")

    path = bundle / "nucleus_boundaries.parquet"
    nid, nlab, npolys = _read_boundary_rings(path)
    nid = np.asarray(nid).astype(str)
    parent = pd.Index(cids).get_indexer(nid)
    nuclei = pd.DataFrame({"cell_id": nid, "label_id": np.asarray(nlab),
                           "geometry": warp(np.asarray(npolys, dtype=object), tf), "parent": parent})
    return cells, nuclei


def nearest_tissue(xy_from, xy_to, tissue_to):
    ok = np.isfinite(xy_from).all(1)
    out = np.full(len(xy_from), "unrouted", dtype=object)
    good = np.isfinite(xy_to).all(1)
    _, j = cKDTree(xy_to[good]).query(xy_from[ok])
    out[ok] = np.asarray(tissue_to, dtype=object)[good][j]
    return out


def centroids(geoms):
    return np.column_stack(sc.centroid_xy(geoms))


def nc_ratio(nuc_geoms, nuc_parent, cell_geoms):
    """Per cell: summed nucleus area over cell area; NaN for a cell with no nucleus."""
    n = len(cell_geoms)
    ok = nuc_parent >= 0
    s = np.bincount(nuc_parent[ok], weights=shapely.area(nuc_geoms[ok]), minlength=n)
    cnt = np.bincount(nuc_parent[ok], minlength=n)
    area = shapely.area(cell_geoms)
    return np.where((cnt > 0) & (area > 0), s / np.where(area > 0, area, 1), np.nan)


def figure(summary, thresholds, path, title):
    """Agreement by IoU threshold, both directions, beside the registration ceiling."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Categorical slots 1-2 identify the platforms; the ceiling is a reference, not a series.
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    colors = {"instanseg": "#2a78d6", "xenium": "#eb6834", "ceiling": "#8a8984"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    for ax, level in zip(axes[:2], LEVELS):
        t = thresholds[(thresholds.level == level) & (thresholds.stratum == "all")]
        ax.plot(t.iou_threshold, t.matched_frac_a, color=colors["instanseg"], lw=2, marker="o", ms=4,
                label="InstanSeg objects matched")
        ax.plot(t.iou_threshold, t.matched_frac_b, color=colors["xenium"], lw=2, marker="s", ms=4,
                label="Xenium objects matched")
        if t.ceiling_frac_a.notna().any():
            ax.plot(t.iou_threshold, t.ceiling_frac_a, color=colors["ceiling"], lw=1.5, ls="--",
                    label="ceiling at recheck residual")
        h = summary["levels"][level]["strata"]["all"]["headline"]
        f2 = {k: ("n/a" if h[k] is None else f"{h[k]:.2f}") for k in ("PQ", "SQ", "RQ")}
        ax.set_title(f"{PLURAL[level]}   PQ {f2['PQ']} (SQ {f2['SQ']} x RQ {f2['RQ']})",
                     fontsize=10, color=ink, loc="left")
        ax.set_xlabel("IoU threshold", color=muted)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(0.05, 0.95)
        ax.grid(axis="y", color=grid, lw=0.6)
        ax.tick_params(colors=muted, labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(grid)
    axes[0].set_ylabel("share of objects with a one-to-one partner", color=muted)
    axes[0].legend(frameon=False, fontsize=8, loc="lower left", labelcolor=ink)

    # Reasons take categorical slots 3-7 so they never reuse the platforms' blue and orange.
    # Three of them sit under 3:1 on white, so every segment of 6 % or more is labelled.
    ax = axes[2]
    order = ["matched", "boundary", "split", "merged", "no_overlap"]
    shades = ["#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
    bars = [(f"{PLURAL[lv]}\nInstanSeg", summary["levels"][lv]["strata"]["all"]["reason_a"]) for lv in LEVELS] + \
           [(f"{PLURAL[lv]}\nXenium", summary["levels"][lv]["strata"]["all"]["reason_b"]) for lv in LEVELS]
    for i, (_, r) in enumerate(bars):
        total = sum(r.values()) or 1
        left = 0.0
        for k, c in zip(order, shades):
            v = r.get(k, 0) / total
            if v <= 0:
                continue
            ax.barh(i, v, left=left, color=c, height=0.62, edgecolor="white", linewidth=2,
                    label=k.replace("_", " "))
            if v >= 0.06:
                ax.text(left + v / 2, i, f"{v:.0%}", ha="center", va="center", fontsize=7.5, color=ink)
            left += v
    handles, labels = ax.get_legend_handles_labels()
    seen = dict(zip(labels, handles))
    ax.legend([seen[k.replace("_", " ")] for k in order if k.replace("_", " ") in seen],
              [k.replace("_", " ") for k in order if k.replace("_", " ") in seen],
              frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              labelcolor=ink)
    ax.set_yticks(range(len(bars)), [b[0] for b in bars], fontsize=8, color=ink)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.tick_params(colors=muted, labelsize=8)
    ax.set_xlabel("share of counted objects at IoU 0.5", color=muted)
    ax.set_title("matched, or why not", fontsize=10, color=ink, loc="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(grid)
    fig.suptitle(title, fontsize=11, color=ink, x=0.01, ha="left")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def compare_section(section, prep: Path, tf, ceiling_um: float, out: Path):
    """Both levels for one section. Writes tables into ``out``; returns (summary, checks, thresholds)."""
    t0 = time.time()
    ins_cells, ins_nuc = instanseg_objects(prep)
    xen_cells, xen_nuc = xenium_objects(Path(section["bundle"]), tf)
    print(f"[seg] InstanSeg {len(ins_cells):,} cells / {len(ins_nuc):,} nuclei; Xenium "
          f"{len(xen_cells):,} cells / {len(xen_nuc):,} nuclei ({time.time() - t0:.0f}s)", flush=True)

    ins_xy = centroids(ins_cells.geometry.to_numpy())
    xen_xy = centroids(xen_cells.geometry.to_numpy())
    # One shared area for both levels, from cell centroids: a nucleus counts where its cell's
    # neighbourhood is shared.
    area = sc.shared_area(ins_xy, xen_xy)
    xen_cells["tissue"] = nearest_tissue(xen_xy, ins_xy, ins_cells.tissue.to_numpy())

    sides = {
        "cell": (ins_cells.geometry.to_numpy(), xen_cells.geometry.to_numpy(),
                 area["a"], area["b"], ins_cells.tissue.to_numpy(), xen_cells.tissue.to_numpy()),
    }
    nuc_par_ok = xen_nuc.parent.to_numpy() >= 0
    if not nuc_par_ok.all():
        print(f"[seg] {int((~nuc_par_ok).sum()):,} Xenium nuclei name no cell in cell_boundaries", flush=True)
    xn_parent = xen_nuc.parent.to_numpy()
    ins_nuc_xy = centroids(ins_nuc.geometry.to_numpy())
    xen_nuc_xy = centroids(xen_nuc.geometry.to_numpy())
    xen_nuc_tissue = np.where(nuc_par_ok, xen_cells.tissue.to_numpy()[np.maximum(xn_parent, 0)],
                              nearest_tissue(xen_nuc_xy, ins_xy, ins_cells.tissue.to_numpy()))
    # Same counted area for nuclei: a nucleus counts if its own centroid falls in it.
    sides["nucleus"] = (ins_nuc.geometry.to_numpy(), xen_nuc.geometry.to_numpy(),
                        sc.inside(area, ins_nuc_xy), sc.inside(area, xen_nuc_xy),
                        ins_cells.tissue.to_numpy()[ins_nuc.parent.to_numpy()], xen_nuc_tissue)

    levels, tables, checks, results = {}, [], {}, {}
    for level in LEVELS:
        a, b, a_dom, b_dom, a_str, b_str = sides[level]
        t1 = time.time()
        r = sc.compare(a, b, a_in_domain=a_dom, b_in_domain=b_dom, a_strata=a_str, b_strata=b_str,
                       ceiling_shift_um=ceiling_um, seed=0)
        results[level] = r
        oa, ob, e = r["objects_a"], r["objects_b"], r["edges"]
        if level == "cell":
            oa.insert(0, "object_id", ins_cells.object_id.to_numpy())
            ob.insert(0, "cell_id", xen_cells.cell_id.to_numpy())
            ob["segmentation_method"] = xen_cells.segmentation_method.to_numpy()
            e = e.assign(object_id=ins_cells.object_id.to_numpy()[e.a], cell_id=xen_cells.cell_id.to_numpy()[e.b])
        else:
            oa.insert(0, "object_id", ins_nuc.object_id.to_numpy())
            ob.insert(0, "cell_id", xen_nuc.cell_id.to_numpy())
            ob.insert(1, "label_id", xen_nuc.label_id.to_numpy())
            e = e.assign(object_id=ins_nuc.object_id.to_numpy()[e.a], label_id=xen_nuc.label_id.to_numpy()[e.b])
        oa.to_parquet(out / f"objects_instanseg_{level}.parquet", index=False)
        ob.to_parquet(out / f"objects_xenium_{level}.parquet", index=False)
        e.to_parquet(out / f"edges_{level}.parquet", index=False)
        tables.append(r["thresholds"].assign(level=level))
        levels[level] = r["summary"]

        # Structural checks, each derivable without the module's own bookkeeping.
        pa, pb = oa.partner.to_numpy(), ob.partner.to_numpy()
        i = np.flatnonzero(pa >= 0)
        j = np.flatnonzero(pb >= 0)
        rng = np.random.default_rng(1)
        sample = rng.choice(len(e), size=min(2000, len(e)), replace=False) if len(e) else np.empty(0, int)
        ga, gb = r["a"], r["b"]
        direct =(shapely.area(shapely.intersection(ga[e.a.to_numpy()[sample]], gb[e.b.to_numpy()[sample]])) /
                  shapely.area(shapely.union(ga[e.a.to_numpy()[sample]], gb[e.b.to_numpy()[sample]])))
        checks[level] = {
            "one_to_one": bool((pb[pa[i]] == i).all() and len(np.unique(pa[i])) == len(i)
                               and (pa[pb[j]] == j).all() and len(j) == len(i)),
            "headline_pairs_meet_threshold": bool((oa.pair_iou.to_numpy()[i] >= sc.HEADLINE_IOU).all()),
            "iou_in_unit_interval": bool(((e.iou > 0) & (e.iou <= 1 + 1e-9)).all()),
            "shares_in_unit_interval": bool(((e.frac_of_a <= 1 + 1e-6) & (e.frac_of_b <= 1 + 1e-6)).all()),
            "edge_iou_matches_direct_sample": bool(np.allclose(e.iou.to_numpy()[sample], direct, rtol=1e-6, atol=1e-9)),
            "reasons_match_partners": bool(
                set(oa.reason) <= REASONS and set(ob.reason) <= REASONS
                and np.array_equal(oa.reason.to_numpy() == "matched", pa >= 0)
                and np.array_equal(ob.reason.to_numpy() == "matched", pb >= 0)),
            "counted_objects_nonempty": bool(a_dom.any() and b_dom.any()),
            "headline_pairs_exist": bool(len(i) > 0 and len(e) > 0),
        }
        print(f"[seg] {level}: {len(e):,} edges, PQ {levels[level]['strata']['all']['headline']['PQ']} "
              f"({time.time() - t1:.0f}s)", flush=True)

    # Xenium cells by the Analyzer's own segmentation method.
    ob = results["cell"]["objects_b"]
    by_method = {m: sc.one_side(ob, (ob.segmentation_method == m) & ob.in_domain)
                 for m in sorted(ob.segmentation_method.unique())}
    # Nuclei matched across platforms: in cells matched to each other?
    hier = sc.hierarchy(results["nucleus"]["objects_a"].partner.to_numpy(),
                        ins_nuc.parent.to_numpy(), xn_parent,
                        results["cell"]["objects_a"].partner.to_numpy(),
                        counted_a=sides["nucleus"][2])
    ins_nc = nc_ratio(ins_nuc.geometry.to_numpy(), ins_nuc.parent.to_numpy(), ins_cells.geometry.to_numpy())
    xen_nc = nc_ratio(xen_nuc.geometry.to_numpy(), xn_parent, xen_cells.geometry.to_numpy())
    summary = {
        "section_id": section["section_id"], "donor_id": section["donor_id"], "roi": section["roi"],
        "a": "InstanSeg (PhenoCycler, QuPath)", "b": "Xenium Analyzer (warped into the PhenoCycler frame)",
        "objects": {"instanseg_cells": len(ins_cells), "instanseg_nuclei": len(ins_nuc),
                    "instanseg_cells_without_nucleus": int(len(ins_cells) - len(ins_nuc)),
                    "xenium_cells": len(xen_cells), "xenium_nuclei": len(xen_nuc),
                    "xenium_nuclei_without_cell": int((~nuc_par_ok).sum()),
                    "xenium_cells_by_nucleus_count": {str(k): int(v) for k, v in pd.Series(
                        np.bincount(xn_parent[nuc_par_ok], minlength=len(xen_cells))).value_counts().sort_index().items()}},
        "shared_area": {k: area[k] for k in ("area_mm2", "grid")},
        "counted": {"instanseg_cells": int(area["a"].sum()), "xenium_cells": int(area["b"].sum()),
                    "instanseg_nuclei": int(sides["nucleus"][2].sum()), "xenium_nuclei": int(sides["nucleus"][3].sum())},
        "ceiling_shift_um": float(ceiling_um) if np.isfinite(ceiling_um) else None,
        "levels": levels,
        "xenium_cells_by_segmentation_method": by_method,
        "nucleus_in_matched_cells": hier,
        "nucleus_to_cell_area_ratio": {"instanseg": sc.quantiles(ins_nc[area["a"]]),
                                       "xenium": sc.quantiles(xen_nc[area["b"]])},
    }
    thresholds = pd.concat(tables, ignore_index=True)
    thresholds.to_csv(out / "thresholds.csv", index=False)
    return summary, checks, thresholds


def run(manifest_path, sid, root):
    section = get_section(load_manifest(manifest_path), sid)
    donor = section["donor_id"]
    base = section_root(root, section)
    root = Path(root)
    receipts = {"registration": base / "registration/current.json",
                "registration_recheck": base / "registration_recheck/current.json",
                "preparation": root / "protein" / f"donor_id={donor}" / "preparation/current.json"}
    rec = {k: current(p, k) for k, p in receipts.items()}
    if rec["registration_recheck"]["inputs"]["preparation"]["sha256"] != identity(receipts["preparation"])["sha256"]:
        raise ValueError("registration_recheck was made against another preparation; rerun recheck_nuclei.py")
    inputs = {"section": section, "receipts": {k: identity(p) for k, p in receipts.items()},
              "code": [identity(HERE / n) for n in ("segmentation_four.py", "segcompare.py", "build_joined.py",
                                                    "join_four.py", "workflow.py")]
                      + [identity(sys.modules[_read_boundary_rings.__module__].__file__)]}
    fingerprint = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    stage = base / "segmentation"
    stage.mkdir(parents=True, exist_ok=True)
    if valid_receipt(stage / "current.json", fingerprint):
        print(f"{sid}: segmentation comparison reused", flush=True)
        return
    out = stage / "attempts" / (fingerprint[:16] + "-" + uuid4().hex[:8])
    out.mkdir(parents=True)
    report = {"status": "running", "section_id": sid, "fingerprint": fingerprint, "inputs": inputs,
              "started_utc": now(), "attempt": str(out), "job_id": os.environ.get("SLURM_JOB_ID")}
    write_json(out / "receipt.json", report)
    unavailable = {k: r.get("reason", "") for k, r in rec.items() if r["status"] == "scientifically_unavailable"}
    if unavailable:
        write_json(out / "unavailable.json", unavailable)
        report.update(status="scientifically_unavailable",
                      reason="; ".join(f"{k}: {v}" for k, v in unavailable.items()), finished_utc=now(),
                      validation={"upstream_unavailable": sorted(unavailable)},
                      outputs=[identity(out / "unavailable.json")])
        write_json(out / "receipt.json", report)
        write_json(stage / "current.json", report)
        print(f"{sid}: segmentation comparison unavailable: {report['reason']}", flush=True)
        return
    try:
        reg = rec["registration"]
        tf = Transform.load(Path(reg["transform_directory"]), stem=reg["transform_stem"])
        ceiling = rec["registration_recheck"].get("summary", {}).get("pair_p50_um")
        ceiling = float(ceiling) if ceiling is not None else float("nan")
        summary, checks, thresholds = compare_section(section, Path(rec["preparation"]["attempt"]), tf, ceiling, out)
        summary["transform"] = {"directory": reg["transform_directory"], "stem": reg["transform_stem"]}
        write_json(out / "segmentation.json", summary)
        figure(summary, thresholds, out / "segmentation.png", f"{sid}: InstanSeg vs Xenium segmentation")
        ok = all(all(c.values()) for c in checks.values())
        head = {lv: summary["levels"][lv]["strata"]["all"]["headline"] for lv in LEVELS}
        report.update(status="validated_complete" if ok else "failed", validation=checks, finished_utc=now(),
                      summary={"headline": head, "counted": summary["counted"],
                               "shared_area_mm2": summary["shared_area"]["area_mm2"],
                               "nucleus_in_matched_cells": summary["nucleus_in_matched_cells"]},
                      outputs=[identity(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "receipt.json"])
        write_json(out / "receipt.json", report)
        if not ok:
            write_json(stage / "last_failure.json", report)
            raise ValueError(f"validation failed: {checks}")
        write_json(stage / "current.json", report)
        print(f"{sid}: segmentation " + "; ".join(
            f"{lv} PQ {h['PQ']} matched InstanSeg {h['matched_frac_a']} Xenium {h['matched_frac_b']}"
            for lv, h in head.items()), flush=True)
    except Exception as exc:
        if report["status"] != "failed":
            report.update(status="failed", error=str(exc), traceback=traceback.format_exc(), finished_utc=now())
            write_json(out / "receipt.json", report)
            write_json(stage / "last_failure.json", report)
        raise


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--section", required=True)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    run(a.manifest, a.section, a.out)
