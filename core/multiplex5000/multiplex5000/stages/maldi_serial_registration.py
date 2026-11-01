"""Serial Registration, MALDI route: each inspected MALDI section onto the slide it neighbours
(2026-10-31).

A MALDI section is a serial section from the same block series as the PhenoCycler and Xenium
slides (maintainer, 2026-10-31): a registration places its spots in the neighbour's frame, never
on the neighbour's cells. The registration images are the ion images themselves.

For each MALDI section the accepted Source Inventory lists with an inspection row:

* **moving** -- its on-tissue spots (MALDI Ingest's ``images/<Tissue_id>.npz``), in the MALDI
  stage frame (um);
* **fixed** -- the inspected neighbour, in its own um frame (PhenoCycler: level-0 pixels x the
  platform profile's pixel size):

  - a PhenoCycler pancreas scan: the outline of its annotation of the inspected tissue (the
    ``Tissue`` annotation for a pancreas section, the ``pLN`` annotation for a pLN section), from
    the site's ``pancreas_annotations`` export;
  - a dedicated PhenoCycler pLN slide: its DAPI tissue mask (the slide's own Otsu level on the
    log of DAPI at ~50 um);
  - a Xenium slide: not built here; the section is reported unavailable with that reason.

* **method** (``multiplex5000.integration.register``) -- the moving outline is moved onto the
  fixed outline's centroid; ``coarse_orientation`` searches rotation and mirror (serial sections
  can be flipped); the best candidates are refined by ``refine_ecc`` on the smoothed outlines,
  rigid first, then affine; the affine is kept only where ``scale_is_physical`` holds, else the
  rigid one. Among candidates the one whose outline overlaps best (Dice) is kept.

Each transform maps MALDI stage um to the neighbour's um and is written with its QC: tissue
Dice, the outline's mutual-nearest-neighbour RMSE (um), scale, anisotropy, shear, rotation and
mirror. The QC is reported, never tuned to; the gate is a review stop -- a registration cannot be
verified by its own metrics, and the maintainer accepts the overlays.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

from .. import gates, maldi, profiles
from ..integration.rasterize import GridSpec, image_to_grid, read_qptiff_channel
from ..integration.register import (coarse_orientation, dice, nearest_neighbour_rmse,
                                    refine_ecc, scale_is_physical, scale_report)
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..integration.transform import Transform
from ..tasks import Input, Task

STAGE = "07_registration"
#: The working grid of the registration: the MALDI spot pitch.
PIXEL_UM = maldi.SPOT_PITCH_UM
#: Coarse candidates refined.
CANDIDATES = 6
#: The annotation class that outlines each tissue on a PhenoCycler pancreas scan.
OUTLINE_CLASS = {PANCREAS: "Tissue", LYMPH_NODE: "pLN"}
NOT_BUILT_XENIUM = "registration to a Xenium neighbour is not built"


class MaldiRegistrationError(ValueError):
    """A MALDI section cannot be registered as recorded."""


def _release(project) -> Path:
    return Path(project.release_dir)


def _rows(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _pixel_um() -> float:
    return float(profiles.load_all()[("platform", "phenocycler")].pixel_size_um)


def _export_tag(image_name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image_name)


def _sections(project) -> list[dict]:
    """The MALDI sections to register: those MALDI Ingest wrote, each with its inspection."""
    inventory = _release(project) / "00_inventory"
    listed = {r["tissue_id"]: r for r in _rows(inventory / "maldi_acquisitions.csv")}
    ingested = [r["tissue_id"] for r in _rows(_release(project) / "04b_maldi_ingest/sections.csv")]
    acquisitions = {r["acquisition_id"]: r for r in _rows(inventory / "acquisitions.csv")}
    out = []
    for tissue_id in ingested:
        row = dict(listed[tissue_id])
        if not row.get("tissue") or not row.get("neighbour_acquisition"):
            raise MaldiRegistrationError(f"{tissue_id}: no inspection in the accepted inventory")
        neighbour = acquisitions.get(row["neighbour_acquisition"])
        if neighbour is None:
            raise MaldiRegistrationError(f"{tissue_id}: neighbour {row['neighbour_acquisition']} "
                                         "is not an inventoried acquisition")
        row["neighbour"] = neighbour
        out.append(row)
    return out


def _neighbour_inputs(project, row: dict) -> list[Input]:
    n = row["neighbour"]
    tid = row["tissue_id"]
    if n["platform"] != "phenocycler":
        return []
    if n["tissue"] == PANCREAS:
        qptiff = Path(n["path"])
        annotations = project.source("pancreas_annotations")
        if annotations is None:
            raise MaldiRegistrationError("the pancreas_annotations source is not named")
        image_id = n["acquisition_id"].split(":", 1)[1]
        return [Input(f"qptiff:{tid}", qptiff, hash=False),
                Input(f"annotations:{tid}", Path(annotations) / f"{_export_tag(image_id)}.geojson")]
    folder = project.source("pln_qptiffs")
    if folder is None:
        raise MaldiRegistrationError("the pln_qptiffs source is not named")
    return [Input(f"qptiff:{tid}", Path(folder) / n["qptiff"], hash=False)]


def settings(project) -> dict:
    return {"route": "maldi_serial", "pixel_um": PIXEL_UM, "candidates": CANDIDATES,
            "outline_class": OUTLINE_CLASS, "scope": "registration_review_candidate"}


def inputs(project) -> list[Input]:
    release = _release(project)
    found = []
    for stage in ("00_inventory", "04b_maldi_ingest"):
        found += [Input(f"{stage}_receipt", release / stage / "receipt.json"),
                  Input(f"{stage}_gate", release / stage / "gate" / "gate.json"),
                  Input(f"{stage}_acceptance", release / stage / "gate" / "accepted.json")]
    found += [Input("acquisitions", release / "00_inventory/acquisitions.csv"),
              Input("maldi_acquisitions", release / "00_inventory/maldi_acquisitions.csv"),
              Input("maldi_sections", release / "04b_maldi_ingest/sections.csv")]
    for row in _sections(project):
        found.append(Input(f"maldi_image:{row['tissue_id']}",
                           release / f"04b_maldi_ingest/images/{row['tissue_id']}.npz"))
        found += _neighbour_inputs(project, row)
    return found


# ------------------------------------------------------------------------------- the images ----
def _polygons(path: Path, wanted: str) -> list[list[np.ndarray]]:
    """Each polygon of the wanted annotation class as [exterior, holes...], in level-0 pixels."""
    data = json.loads(Path(path).read_text())
    features = data.get("features", []) if isinstance(data, dict) else data
    out = []
    for feature in features:
        name = ((feature.get("properties") or {}).get("classification") or {}).get("name", "")
        if name != wanted:
            continue
        geometry = feature.get("geometry") or {}
        polygons = {"Polygon": [geometry.get("coordinates", [])],
                    "MultiPolygon": geometry.get("coordinates", [])}.get(geometry.get("type"), [])
        out += [[np.asarray(ring, float) for ring in polygon] for polygon in polygons if polygon]
    return out


def _fill(polygons_um: list[list[np.ndarray]], grid: GridSpec) -> np.ndarray:
    from skimage.draw import polygon as fill

    mask = np.zeros(grid.shape, bool)
    for rings in polygons_um:
        for k, ring in enumerate(rings):
            px = grid.to_pixels(ring)
            rr, cc = fill(px[:, 1], px[:, 0], grid.shape)
            mask[rr, cc] = k == 0       # the exterior fills, each hole clears
    return mask


def _centred_grid(centre: np.ndarray, half_um: float) -> GridSpec:
    side = int(np.ceil(2 * half_um / PIXEL_UM))
    return GridSpec(float(centre[0] - side * PIXEL_UM / 2), float(centre[1] - side * PIXEL_UM / 2),
                    PIXEL_UM, side, side)


def _dapi_mask(image: np.ndarray) -> np.ndarray:
    """A slide's tissue from its own DAPI: log, smoothed, Otsu's level, holes filled, pieces
    under 0.05 mm^2 dropped."""
    from scipy import ndimage
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects

    values = ndimage.gaussian_filter(np.log1p(np.maximum(image, 0)), 1.0)
    mask = values > threshold_otsu(values)
    mask = ndimage.binary_fill_holes(mask)
    return remove_small_objects(mask, min_size=20)


def _fixed(task: Task, row: dict, pixel_um: float) -> dict:
    """The neighbour's outline (and its DAPI, for the overlay) on a grid in its own um frame."""
    tid, n = row["tissue_id"], row["neighbour"]
    dapi, um_per_px = read_qptiff_channel(task.input(f"qptiff:{tid}"), 0, target_um=PIXEL_UM,
                                          base_um_per_px=pixel_um)
    if n["tissue"] == PANCREAS:
        polygons = [[ring * pixel_um for ring in rings]
                    for rings in _polygons(task.input(f"annotations:{tid}"),
                                           OUTLINE_CLASS[row["tissue"]])]
        if not polygons:
            raise MaldiRegistrationError(f"{tid}: {n['acquisition_id']} has no "
                                         f"{OUTLINE_CLASS[row['tissue']]} annotation")
        points = np.vstack([rings[0] for rings in polygons])
        lo, hi = points.min(axis=0), points.max(axis=0)
        outline_grid = GridSpec(float(lo[0]), float(lo[1]), PIXEL_UM,
                                int(np.ceil((hi[1] - lo[1]) / PIXEL_UM)) + 1,
                                int(np.ceil((hi[0] - lo[0]) / PIXEL_UM)) + 1)
        mask = _fill(polygons, outline_grid)
        origin = (outline_grid.origin_x_um, outline_grid.origin_y_um)
        source = f"{OUTLINE_CLASS[row['tissue']]} annotation"
    else:
        mask = _dapi_mask(dapi)
        polygons = None
        # a level pixel i covers [i, i + 1) level pixels of the slide: its centre is half a pixel in
        origin = (um_per_px / 2, um_per_px / 2)
        outline_grid = GridSpec(*origin, um_per_px, *mask.shape)
        source = "DAPI tissue mask (Otsu)"
    rows_, cols_ = np.nonzero(mask)
    centre = outline_grid.to_microns(np.column_stack([cols_.mean(), rows_.mean()]))[0]
    extent = np.hypot(*(np.array(mask.shape[::-1]) * outline_grid.pixel_um)) / 2
    return dict(mask=mask, mask_grid=outline_grid, polygons=polygons, dapi=dapi,
                dapi_um_per_px=um_per_px, centre=centre, extent=extent, origin=origin,
                source=source)


def _moving(task: Task, tid: str) -> dict:
    with np.load(task.input(f"maldi_image:{tid}")) as data:
        summed, on = data["summed"], data["on_tissue"]
        grid = GridSpec.from_dict(json.loads(str(data["grid"])))
    rows_, cols_ = np.nonzero(on)
    xy = grid.to_microns(np.column_stack([cols_, rows_]))
    return dict(summed=summed, on=on, grid=grid, xy=xy, centre=xy.mean(axis=0),
                extent=float(np.max(np.hypot(*(xy - xy.mean(axis=0)).T))))


def _onto(values_xy: np.ndarray, grid: GridSpec, values=None) -> np.ndarray:
    return maldi.spot_image(values_xy[:, 0], values_xy[:, 1],
                            np.ones(len(values_xy)) if values is None else values, grid, fill=0.0)


def _outline_points(mask: np.ndarray, grid: GridSpec) -> np.ndarray:
    from skimage.segmentation import find_boundaries

    rr, cc = np.nonzero(find_boundaries(mask, mode="inner"))
    return grid.to_microns(np.column_stack([cc, rr]))


def register_section(moving: dict, fixed: dict) -> tuple[Transform, dict]:
    """The MALDI-stage-um -> neighbour-um transform of one section, and its QC."""
    from scipy import ndimage

    half = max(fixed["extent"], moving["extent"]) + 1000.0
    grid = _centred_grid(fixed["centre"], half)
    fixed_mask = image_to_grid(fixed["mask"].astype(np.float32), fixed["mask_grid"].pixel_um, grid,
                               origin_um=fixed["origin"]) > 0.5
    t0 = Transform(affine=np.array([[1.0, 0, fixed["centre"][0] - moving["centre"][0]],
                                    [0, 1.0, fixed["centre"][1] - moving["centre"][1]]]))
    moved = t0.apply(moving["xy"])
    moving_mask = _onto(moved, grid) > 0
    smooth = {k: ndimage.gaussian_filter(m.astype(np.float32), 1.5)
              for k, m in (("fixed", fixed_mask), ("moving", moving_mask))}
    best = None
    for candidate in coarse_orientation(smooth["moving"], smooth["fixed"], grid)[:CANDIDATES]:
        rigid, ecc_rigid = refine_ecc(smooth["moving"], smooth["fixed"], grid, candidate,
                                      motion="euclidean")
        affine, ecc_affine = refine_ecc(smooth["moving"], smooth["fixed"], grid, rigid,
                                        motion="affine")
        physical, why = scale_is_physical(affine)
        chosen, motion, value = ((affine, "affine", ecc_affine)
                                 if physical and np.isfinite(ecc_affine)
                                 else (rigid, "rigid", ecc_rigid))
        warped = _onto(chosen.apply(moved), grid) > 0
        warped = ndimage.binary_closing(warped, iterations=1)
        score = dice(warped, fixed_mask)
        if best is None or score > best[0]:
            best = (score, chosen, motion, value, why, candidate.meta, warped)
    score, chosen, motion, value, why, coarse, warped = best
    total = chosen.compose(t0)
    rmse, pairs = nearest_neighbour_rmse(_outline_points(warped, grid),
                                         _outline_points(fixed_mask, grid))
    report = scale_report(total)
    total.meta.update({"moving_frame": "MALDI stage um", "fixed_frame": "neighbour um",
                       "motion": motion})
    qc = dict(tissue_dice=score, outline_rmse_um=rmse, outline_pairs=pairs, ecc=value,
              motion=motion, affine_scale_check=why, rotation_deg=total.rotation_deg,
              mirrored=total.is_mirrored, coarse_angle_deg=coarse.get("coarse_angle_deg"),
              coarse_mirror=coarse.get("coarse_mirror"), **report)
    return total, dict(qc, grid=grid, fixed_mask=fixed_mask, warped=warped)


def run(task: Task) -> dict:
    pixel_um = _pixel_um()
    rows, results = [], []
    for row in _sections(task.project):
        tid, n = row["tissue_id"], row["neighbour"]
        record = dict(tissue_id=tid, donor_id=row["donor_id"], tissue=row["tissue"],
                      neighbour_acquisition=n["acquisition_id"], neighbour_platform=n["platform"],
                      neighbour_tissue=n["tissue"])
        if n["platform"] != "phenocycler":
            record.update(status="unavailable", reason=NOT_BUILT_XENIUM)
            rows.append(record)
            continue
        fixed = _fixed(task, row, pixel_um)
        moving = _moving(task, tid)
        transform, qc = register_section(moving, fixed)
        transform.save(task.output(f"transforms/{maldi.safe_key(tid)}"))
        physical, _ = scale_is_physical(transform)
        record.update(status="registered", fixed_outline=fixed["source"],
                      scale_physical=physical,
                      **{k: v for k, v in qc.items() if k not in ("grid", "fixed_mask", "warped")})
        rows.append(record)
        results.append((record, fixed, moving, transform, qc))
    _write(task.output("registration.csv"), rows)
    _figure(task, results)
    registered = [r for r in rows if r["status"] == "registered"]
    return {"sections": len(rows), "registered": len(registered),
            "unavailable": sorted(f"{r['tissue_id']}: {r['reason']}" for r in rows
                                  if r["status"] != "registered"),
            "not_physical": sorted(r["tissue_id"] for r in registered if not r["scale_physical"]),
            "qc": {r["tissue_id"]: dict(dice=round(r["tissue_dice"], 3),
                                        rmse_um=round(r["outline_rmse_um"], 1),
                                        motion=r["motion"], mirrored=r["mirrored"],
                                        rotation_deg=round(r["rotation_deg"], 1))
                   for r in registered}}


def _write(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        columns += [k for k in row if k not in columns]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or ["empty"])
        writer.writeheader()
        writer.writerows(rows)


def _figure(task: Task, results: list) -> None:
    """Each registered section: the neighbour's DAPI with its outline (green), and the MALDI
    summed-ion image carried through the transform, with its outline (orange)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = max(len(results), 1)
    per_row = min(n, 4)
    rows = -(-n // per_row)
    fig, axes = plt.subplots(rows, per_row, figsize=(4.2 * per_row, 4.4 * rows), squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for ax, (record, fixed, moving, transform, qc) in zip(axes.ravel(), results):
        grid = qc["grid"]
        dapi = image_to_grid(fixed["dapi"], fixed["dapi_um_per_px"], grid)
        ax.imshow(np.log1p(np.maximum(dapi, 0)), cmap="gray")
        g = moving["grid"]
        rows_, cols_ = np.nonzero(np.isfinite(moving["summed"]))
        xy = g.to_microns(np.column_stack([cols_, rows_]))
        summed = np.log1p(np.maximum(moving["summed"][rows_, cols_], 0))
        carried = maldi.spot_image(*transform.apply(xy).T, summed, grid)
        ax.imshow(np.ma.masked_invalid(carried), cmap="magma", alpha=0.45)
        ax.contour(qc["fixed_mask"], levels=[0.5], colors="#2ca02c", linewidths=0.8)
        ax.contour(qc["warped"], levels=[0.5], colors="#ff7f0e", linewidths=0.8)
        ax.set_title(f"{record['tissue_id']} -> {record['neighbour_acquisition'].replace(' - resolution #1', '')}\n"
                     f"Dice {record['tissue_dice']:.3f}, outline RMSE {record['outline_rmse_um']:.0f} um, "
                     f"{record['motion']}{', mirrored' if record['mirrored'] else ''}", fontsize=7)
    fig.suptitle("MALDI sections on the slide each neighbours: neighbour outline green, MALDI "
                 "outline orange, MALDI summed ions over the neighbour's DAPI", fontsize=10)
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "maldi_registration.png"), dpi=110)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    C = gates.Check
    checks = [
        C("every inspected MALDI section is registered or has its reason", summary["sections"],
          "each registered, or unavailable with its reason",
          summary["sections"] == summary["registered"] + len(summary["unavailable"]),
          detail="; ".join(summary["unavailable"]) or "all registered"),
        C("every kept transform's scale is physical", summary["not_physical"], "none",
          not summary["not_physical"],
          detail="an affine whose scale is not physical is replaced by the rigid fit; a rigid fit "
                 "is physical by construction (register.scale_is_physical)"),
        C("registration QC per section", summary["qc"], "reported", True, blocking=False,
          detail="tissue Dice and outline RMSE are evidence for review, never a target; the "
                 "RMSE also sets which objects Serial Integration can evaluate"),
        C("MALDI sections not registered", summary["unavailable"], "reported", True,
          blocking=False),
    ]
    return gates.Gate(STAGE, "Each MALDI section on the slide it neighbours: a review stop. The "
                      "transforms place spots in the neighbour's frame; serial sections are "
                      "different cells.", checks, ["maldi_registration.png"])
