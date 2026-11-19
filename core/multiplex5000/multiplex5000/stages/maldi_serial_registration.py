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
  - an exact Xenium region: its instrument-declared DAPI focus image, in the region's
    native micron frame, masked by the same DAPI/Otsu method. Missing or damaged
    regions remain unavailable; another region is never substituted.

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
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .. import gates, maldi, profiles
from ..integration.rasterize import GridSpec, image_to_grid, read_qptiff_channel
from ..integration.register import (
    coarse_orientation,
    dice,
    nearest_neighbour_rmse,
    refine_ecc,
    scale_is_physical,
    scale_report,
)
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..integration.transform import Transform
from ..resolved_parameters import write_resolved_parameters
from ..tasks import Input, Task, sha256_file

STAGE = "07_registration"
#: The working grid of the registration: the MALDI spot pitch.
PIXEL_UM = maldi.SPOT_PITCH_UM
#: Coarse candidates refined.
CANDIDATES = 6
DAPI_SMOOTH_SIGMA_PX = 1.0
DAPI_MIN_COMPONENT_PIXELS = 20
OUTLINE_SMOOTH_SIGMA_PX = 1.5
REGISTRATION_PADDING_UM = 1000.0
#: The annotation class that outlines each tissue on a PhenoCycler pancreas scan.
OUTLINE_CLASS = {PANCREAS: "Tissue", LYMPH_NODE: "pLN"}


class MaldiRegistrationError(ValueError):
    """A MALDI section cannot be registered as recorded."""


class MaldiEvidenceUnavailable(MaldiRegistrationError):
    """A declared section lacks usable image/outline evidence; retain its disposition."""


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
    listed_rows = _rows(inventory / "maldi_acquisitions.csv")
    listed = {r["tissue_id"]: r for r in listed_rows}
    ingested = [r["tissue_id"] for r in _rows(_release(project) / "04b_maldi_ingest/sections.csv")]
    acquisition_rows = _rows(inventory / "acquisitions.csv")
    acquisitions = {r["acquisition_id"]: r for r in acquisition_rows}
    if len(listed) != len(listed_rows) or len(acquisitions) != len(acquisition_rows) or len(set(ingested)) != len(ingested):
        raise MaldiRegistrationError("duplicate section or acquisition identifiers")
    if set(ingested) != set(listed):
        raise MaldiRegistrationError("ingested sections must cover the declared MALDI roster exactly")
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
    if n["platform"] == "xenium":
        source = _xenium_source(project, n)
        found = [Input(f"xenium_experiment:{tid}", source["experiment"], required=False)]
        if source.get("image") is not None:
            # The new route binds the actual morphology bytes, not only file metadata.
            found.append(Input(f"xenium_dapi:{tid}", source["image"], required=False))
        return found
    if n["platform"] != "phenocycler":
        return []
    if n["tissue"] == PANCREAS:
        qptiff = Path(n["path"])
        annotations = project.source("pancreas_annotations")
        if annotations is None:
            raise MaldiRegistrationError("the pancreas_annotations source is not named")
        image_id = n["acquisition_id"].split(":", 1)[1]
        return [Input(f"qptiff:{tid}", qptiff, hash=False, required=False),
                Input(f"annotations:{tid}", Path(annotations) / f"{_export_tag(image_id)}.geojson", required=False)]
    folder = project.source("pln_qptiffs")
    if folder is None:
        raise MaldiRegistrationError("the pln_qptiffs source is not named")
    return [Input(f"qptiff:{tid}", Path(folder) / n["qptiff"], hash=False, required=False)]


def settings(project) -> dict:
    return {"route": "maldi_serial", "pixel_um": PIXEL_UM, "candidates": CANDIDATES,
            "outline_class": OUTLINE_CLASS, "scope": "registration_review_candidate",
            "xenium_outline": "instrument-declared DAPI focus; existing Otsu mask on 50 um grid",
            "dapi_mask": {"smooth_sigma_pixels": DAPI_SMOOTH_SIGMA_PX, "min_component_pixels": DAPI_MIN_COMPONENT_PIXELS,
                          "threshold": "section-fitted Otsu on smoothed log1p DAPI", "fill_holes": True},
            "outline_smooth_sigma_pixels": OUTLINE_SMOOTH_SIGMA_PX,
            "registration_padding_um": REGISTRATION_PADDING_UM,
            "xenium_frame": "exact region-native um; original image origin/orientation; no cell typing"}


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
    index = project.source("xenium_source_index")
    if index is not None:
        found.append(Input("xenium_source_index", Path(index)))
    for row in _sections(project):
        found.append(Input(f"maldi_image:{row['tissue_id']}",
                           release / f"04b_maldi_ingest/images/{row['tissue_id']}.npz"))
        found += _neighbour_inputs(project, row)
    return found


# ------------------------------------------------------------------------------- the images ----
def _xenium_source(project, neighbour: dict) -> dict:
    """Resolve only the exact selected output, and bind its instrument identity to Stage 00."""
    acquisition = neighbour["acquisition_id"]
    output = neighbour.get("output", "")
    bundle = Path(neighbour.get("path", ""))
    index = project.source("xenium_source_index")
    if index is not None:
        rows = [r for r in _rows(Path(index)) if r["acquisition_id"] == acquisition]
        if len(rows) != 1 or rows[0]["output"] != output:
            raise MaldiRegistrationError(f"{acquisition}: source index must name one exact output")
        bundle = Path(rows[0]["source_path"])
    # A declared relocation may use a different directory name. Its exact
    # acquisition/output mapping and inventoried experiment hash establish
    # identity below; directory spelling is not scientific identity.
    if not output or not bundle.is_absolute() or (index is None and bundle.name != output):
        raise MaldiRegistrationError(f"{acquisition}: missing or mismatched exact Xenium output")
    result = {"experiment": bundle / "experiment.xenium", "image": None,
              "status": "unavailable", "reason": "exact Xenium experiment record is unavailable"}
    if neighbour.get("condition") == "damaged":
        result["reason"] = "recorded damaged Xenium tissue"
        return result
    if not result["experiment"].is_file():
        return result
    experiment = json.loads(result["experiment"].read_text())
    if acquisition != f"xenium:{experiment.get('slide_id')}:{experiment.get('region_name')}":
        raise MaldiRegistrationError(f"{acquisition}: experiment names another slide or region")
    receipt = json.loads((_release(project) / "00_inventory/receipt.json").read_text())
    bindings = [item.get("sha256") for item in receipt["inputs"]
                if item["name"] == f"xenium:{output}:experiment.xenium"]
    if len(bindings) != 1:
        raise MaldiRegistrationError(f"{acquisition}: experiment needs one unique inventoried binding")
    expected = bindings[0]
    if not expected or sha256_file(result["experiment"]) != expected:
        raise MaldiRegistrationError(f"{acquisition}: experiment differs from inventoried input")
    pixel = experiment.get("pixel_size")
    if type(pixel) not in (int, float) or not np.isfinite(pixel) or pixel <= 0:
        raise MaldiRegistrationError(f"{acquisition}: missing physical pixel size")
    relative = (experiment.get("images") or {}).get("morphology_focus_filepath")
    if not isinstance(relative, str) or not relative:
        result["reason"] = "experiment does not declare a morphology focus image"
        return result
    image = (bundle / relative).resolve()
    if Path(relative).is_absolute() or not image.is_relative_to(bundle.resolve()):
        raise MaldiRegistrationError(f"{acquisition}: morphology focus escapes the exact bundle")
    result.update(image=image, pixel_um=float(pixel), experiment_sha256=expected,
                  source_binding={"inventory_output": output, "bundle_path": str(bundle),
                                  "explicit_source_index": str(index) if index is not None else None,
                                  "source_index_sha256": sha256_file(Path(index)) if index is not None else None})
    result.update(status="ready" if image.is_file() else "unavailable",
                  reason="" if image.is_file() else "instrument-declared Xenium DAPI focus is unavailable")
    return result


def _read_xenium_dapi(path: Path, pixel_um: float) -> tuple[np.ndarray, float, dict]:
    """Read the declared DAPI plane only, never implicit OME companion files or a z-stack."""
    import tifffile

    with tifffile.TiffFile(path) as tif:
        if not tif.ome_metadata:
            raise MaldiRegistrationError("Xenium DAPI needs OME physical/channel metadata")
        ome = ET.fromstring(tif.ome_metadata)
        pixel_nodes = ome.findall(".//{*}Pixels")
        pixels = pixel_nodes[0] if len(pixel_nodes) == 1 else None
        if pixels is None or pixels.get("SizeZ") != "1" or pixels.get("SizeT") != "1":
            raise MaldiRegistrationError("Xenium focus must be one z/time plane")
        for axis in ("X", "Y"):
            unit = pixels.get(f"PhysicalSize{axis}Unit")
            size = float(pixels.get(f"PhysicalSize{axis}", "nan"))
            if unit not in {"µm", "um"} or size != pixel_um:
                raise MaldiRegistrationError("Xenium OME physical pixel size differs from experiment")
        channels = pixels.findall("{*}Channel")
        dapi = [i for i, c in enumerate(channels) if c.get("Name", "").upper() == "DAPI"]
        if len(dapi) != 1:
            raise MaldiRegistrationError("Xenium focus must explicitly identify one DAPI channel")
        planes = [p for p in pixels.findall("{*}TiffData")
                  if int(p.get("FirstC", "0")) == dapi[0]
                  and int(p.get("FirstZ", "0")) == 0 and int(p.get("FirstT", "0")) == 0]
        if len(planes) != 1:
            raise MaldiRegistrationError("Xenium focus has no unique DAPI plane mapping")
        uuid = planes[0].find("{*}UUID")
        if uuid is not None and uuid.get("FileName", path.name) != path.name:
            raise MaldiRegistrationError("declared focus file does not contain the DAPI plane")
        if int(planes[0].get("PlaneCount", "1")) != 1:
            raise MaldiRegistrationError("DAPI mapping must identify one exact TIFF plane")
        page = tif.pages[int(planes[0].get("IFD", "0"))]
        levels = [page, *list(page.pages or [])]
        shape = tuple(page.shape)
        if len(shape) != 2:
            raise MaldiRegistrationError("Xenium DAPI plane must be a two-dimensional raster")
        # Xenium's native focus pyramid uses factor-two reductions, with edge truncation.
        # Check that layout explicitly instead of deriving a slightly different x/y scale
        # from rounded image dimensions or reading undeclared companion channels.
        for i, level in enumerate(levels):
            if tuple(level.shape) != tuple(n // (2 ** i) for n in shape):
                raise MaldiRegistrationError("unsupported Xenium focus pyramid geometry")
        selected = min(range(len(levels)), key=lambda i: abs(pixel_um * (2 ** i) - PIXEL_UM))
        image = np.asarray(levels[selected].asarray(), dtype=np.float32)
    if not np.isfinite(image).all() or np.any(image < 0):
        raise MaldiRegistrationError("Xenium DAPI contains invalid measurements")
    return image, pixel_um * (2 ** selected), {"dapi_channel": dapi[0], "pyramid_level": selected,
                                              "base_pixel_um": pixel_um, "base_shape": list(shape)}


def _fixed_xenium(task: Task, row: dict) -> dict:
    source = _xenium_source(task.project, row["neighbour"])
    if source["status"] != "ready":
        raise MaldiRegistrationError(source["reason"])
    tid = row["tissue_id"]
    if (task.input(f"xenium_experiment:{tid}") != source["experiment"] or
            task.input(f"xenium_dapi:{tid}").resolve() != source["image"]):
        raise MaldiRegistrationError(f"{tid}: Xenium source changed after input declaration")
    image, um_per_px, metadata = _read_xenium_dapi(source["image"], source["pixel_um"])
    # Original focus images use the region-local image corner as origin, with no
    # externally applied crop/rotation. The grid and overlay share pixel centres.
    grid = GridSpec(PIXEL_UM / 2, PIXEL_UM / 2, PIXEL_UM,
                    int(np.ceil(image.shape[0] * um_per_px / PIXEL_UM)),
                    int(np.ceil(image.shape[1] * um_per_px / PIXEL_UM)))
    dapi = image_to_grid(image, um_per_px, grid, origin_um=(um_per_px / 2, um_per_px / 2))
    mask, mask_parameters = _dapi_mask(dapi, with_parameters=True)
    rr, cc = np.nonzero(mask)
    if not len(rr):
        raise MaldiEvidenceUnavailable(f"{tid}: Xenium DAPI tissue mask is empty")
    centre = grid.to_microns(np.array([[cc.mean(), rr.mean()]]))[0]
    return dict(mask=mask, mask_grid=grid, polygons=None, dapi=dapi, dapi_um_per_px=PIXEL_UM,
                dapi_origin=(grid.origin_x_um, grid.origin_y_um), centre=centre,
                extent=float(np.hypot(*np.array(mask.shape[::-1]) * PIXEL_UM) / 2),
                origin=(grid.origin_x_um, grid.origin_y_um), source="Xenium DAPI tissue mask (Otsu)",
                frame_metadata={"source": "instrument-declared original morphology focus", "units": "um",
                                "experiment_sha256": source["experiment_sha256"],
                                "source_binding": source["source_binding"],
                                "mask_parameters": mask_parameters,
                                "region": row["neighbour"]["acquisition_id"], "crop_offset_um": [0.0, 0.0],
                                "orientation": "native instrument orientation", **metadata})


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


def _dapi_mask(image: np.ndarray, *, with_parameters: bool = False):
    """A slide's tissue from its own DAPI: log, smoothed, Otsu's level, holes filled, pieces
    under 0.05 mm^2 dropped."""
    from scipy import ndimage
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects

    values = ndimage.gaussian_filter(np.log1p(np.maximum(image, 0)), DAPI_SMOOTH_SIGMA_PX)
    threshold = float(threshold_otsu(values))
    mask = values > threshold
    mask = ndimage.binary_fill_holes(mask)
    mask = remove_small_objects(mask, min_size=DAPI_MIN_COMPONENT_PIXELS)
    parameters = {"method": "Otsu on smoothed log1p DAPI", "threshold_log1p": threshold,
                  "smooth_sigma_pixels": DAPI_SMOOTH_SIGMA_PX,
                  "min_component_pixels": DAPI_MIN_COMPONENT_PIXELS, "fill_holes": True}
    return (mask, parameters) if with_parameters else mask


def _fixed(task: Task, row: dict, pixel_um: float) -> dict:
    """The neighbour's outline (and its DAPI, for the overlay) on a grid in its own um frame."""
    tid, n = row["tissue_id"], row["neighbour"]
    if not task.input(f"qptiff:{tid}").is_file():
        raise MaldiEvidenceUnavailable(f"{tid}: exact neighbour image is unavailable")
    dapi, um_per_px = read_qptiff_channel(task.input(f"qptiff:{tid}"), 0, target_um=PIXEL_UM,
                                          base_um_per_px=pixel_um)
    if n["tissue"] == PANCREAS:
        if not task.input(f"annotations:{tid}").is_file():
            raise MaldiEvidenceUnavailable(f"{tid}: exact neighbour annotation export is unavailable")
        polygons = [[ring * pixel_um for ring in rings]
                    for rings in _polygons(task.input(f"annotations:{tid}"),
                                           OUTLINE_CLASS[row["tissue"]])]
        if not polygons:
            raise MaldiEvidenceUnavailable(f"{tid}: {n['acquisition_id']} has no "
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
        mask, mask_parameters = _dapi_mask(dapi, with_parameters=True)
        polygons = None
        # a level pixel i covers [i, i + 1) level pixels of the slide: its centre is half a pixel in
        origin = (um_per_px / 2, um_per_px / 2)
        outline_grid = GridSpec(*origin, um_per_px, *mask.shape)
        source = "DAPI tissue mask (Otsu)"
    rows_, cols_ = np.nonzero(mask)
    if not len(rows_):
        raise MaldiEvidenceUnavailable(f"{tid}: fixed tissue outline is empty")
    centre = outline_grid.to_microns(np.column_stack([cols_.mean(), rows_.mean()]))[0]
    extent = np.hypot(*(np.array(mask.shape[::-1]) * outline_grid.pixel_um)) / 2
    return dict(mask=mask, mask_grid=outline_grid, polygons=polygons, dapi=dapi,
                dapi_um_per_px=um_per_px, dapi_origin=(um_per_px / 2, um_per_px / 2),
                centre=centre, extent=extent, origin=origin,
                source=source, frame_metadata={"source": source, "pixel_um": pixel_um, "units": "um",
                                               "orientation": "native scan orientation",
                                               "mask_pixel_um": outline_grid.pixel_um,
                                               "mask_parameters": mask_parameters if polygons is None else None})


def _moving(task: Task, tid: str) -> dict:
    with np.load(task.input(f"maldi_image:{tid}")) as data:
        summed, on = data["summed"], data["on_tissue"]
        grid = GridSpec.from_dict(json.loads(str(data["grid"])))
    rows_, cols_ = np.nonzero(on)
    if not len(rows_):
        raise MaldiEvidenceUnavailable(f"{tid}: MALDI signal mask contains no spots")
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

    if (not np.asarray(fixed["mask"]).any() or not len(moving["xy"])
            or not np.isfinite(moving["xy"]).all()
            or not np.isfinite(fixed["centre"]).all()):
        raise MaldiRegistrationError("registration needs finite spots and a nonempty fixed outline")

    half = max(fixed["extent"], moving["extent"]) + REGISTRATION_PADDING_UM
    grid = _centred_grid(fixed["centre"], half)
    fixed_mask = image_to_grid(fixed["mask"].astype(np.float32), fixed["mask_grid"].pixel_um, grid,
                               origin_um=fixed["origin"]) > 0.5
    t0 = Transform(affine=np.array([[1.0, 0, fixed["centre"][0] - moving["centre"][0]],
                                    [0, 1.0, fixed["centre"][1] - moving["centre"][1]]]))
    moved = t0.apply(moving["xy"])
    moving_mask = _onto(moved, grid) > 0
    smooth = {k: ndimage.gaussian_filter(m.astype(np.float32), OUTLINE_SMOOTH_SIGMA_PX)
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
        try:
            if n["platform"] == "xenium":
                source = _xenium_source(task.project, n)
                if source["status"] != "ready":
                    raise MaldiEvidenceUnavailable(source["reason"])
                fixed = _fixed_xenium(task, row)
            elif n["platform"] == "phenocycler":
                fixed = _fixed(task, row, pixel_um)
            else:
                record.update(status="unavailable", reason=f"unsupported neighbour platform: {n['platform']}")
                rows.append(record)
                continue
            moving = _moving(task, tid)
            transform, qc = register_section(moving, fixed)
        except MaldiEvidenceUnavailable as exc:
            record.update(status="unavailable", reason=str(exc))
            rows.append(record)
            continue
        frame = f"{n['acquisition_id']}:um"
        transform.meta.update(fixed_acquisition_id=n["acquisition_id"], fixed_frame=frame,
                              moving_acquisition_id=f"maldi:{tid}", moving_frame=f"maldi:{tid}:stage_um",
                              coordinate_units="um", fixed_geometry=fixed.get("frame_metadata", {
                                  "source": fixed["source"], "pixel_um": pixel_um,
                                  "units": "um", "orientation": "native scan orientation"}))
        transform.save(task.output(f"transforms/{maldi.safe_key(tid)}"))
        physical, _ = scale_is_physical(transform)
        record.update(status="registered", fixed_outline=fixed["source"],
                      mask_parameters=fixed.get("frame_metadata", {}).get("mask_parameters"),
                      fixed_frame=frame, moving_frame=f"maldi:{tid}:stage_um", coordinate_units="um",
                      scale_physical=physical,
                      **{k: v for k, v in qc.items() if k not in ("grid", "fixed_mask", "warped")})
        rows.append(record)
        results.append((record, fixed, moving, transform, qc))
    _write(task.output("registration.csv"), rows)
    _figure(task, results)
    registered = [r for r in rows if r["status"] == "registered"]
    summary = {"sections": len(rows), "registered": len(registered),
            "unavailable": sorted(f"{r['tissue_id']}: {r['reason']}" for r in rows
                                  if r["status"] != "registered"),
            "not_physical": sorted(r["tissue_id"] for r in registered if not r["scale_physical"]),
            "qc": {r["tissue_id"]: dict(dice=round(r["tissue_dice"], 3),
                                        rmse_um=round(r["outline_rmse_um"], 1),
                                        motion=r["motion"], mirrored=r["mirrored"],
                                        rotation_deg=round(r["rotation_deg"], 1))
                   for r in registered}}
    summary["parameter_recording"] = write_resolved_parameters(
        task, operation="serial_registration",
        context={"modality": "maldi", "relationship_type": "serial_section", "units": "um",
                 "approval": "registration_review_candidate"},
        parameters=[{"name": "spot_pitch_um", "role": "fixed", "value": PIXEL_UM, "units": "um"},
                    {"name": "coarse_candidates", "role": "fixed", "value": CANDIDATES, "units": "candidates"},
                    {"name": "mask_and_registration_settings", "role": "fixed", "value": settings(task.project), "units": "per-setting pixels/um"},
                    {"name": "method", "role": "selected_method", "value": "orientation/mirror, rigid/affine ECC, physical-scale constraint", "units": None}],
        acquisitions=[{k: None if isinstance(v, (float, np.floating)) and not np.isfinite(v)
                        else v.item() if isinstance(v, np.generic) else v for k, v in row.items()} for row in rows],
        input_names=tuple(name for name in ("00_inventory_receipt", "04b_maldi_ingest_receipt", "acquisitions", "maldi_acquisitions") if name in task.inputs),
        producer_paths=(Path(__file__),))
    return summary


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
    for ax, (record, fixed, moving, transform, qc) in zip(axes.ravel(), results, strict=False):
        grid = qc["grid"]
        dapi = image_to_grid(fixed["dapi"], fixed["dapi_um_per_px"], grid,
                             origin_um=fixed.get("dapi_origin", (0.0, 0.0)))
        ax.imshow(np.log1p(np.maximum(dapi, 0)), cmap="gray")
        g = moving["grid"]
        rows_, cols_ = np.nonzero(np.isfinite(moving["summed"]))
        xy = g.to_microns(np.column_stack([cols_, rows_]))
        summed = maldi.display(moving["summed"],
                               *maldi.signal_range(moving["summed"], moving["on"]))[rows_, cols_]
        carried = maldi.spot_image(*transform.apply(xy).T, summed, grid)
        ax.imshow(np.ma.masked_invalid(carried), cmap="magma", vmin=0, vmax=1, alpha=0.45)
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
