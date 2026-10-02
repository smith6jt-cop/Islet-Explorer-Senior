"""The tissue block a cell belongs to -- the unit UniFORM aligns.

THE RULING (maintainer, 2026-10-02): "UniFORM should be done on a same-block basis."

Ten pancreas-project images carry a lymph node inside a `pLN` annotation.  A SAME-BLOCK node sits
in the pancreas block itself -- one piece of tissue, one fixation, one section -- so its cells
belong to that block's sample.  A DIFFERENT-BLOCK node is another block mounted on the same slide:
a different sample.  Its cells take no part in estimating the pancreas block's shift, in the
usability statistics, or in the cohort scale, and they receive no aligned value from the pancreas
alignment: in `celltype/uniform` their marker values are unavailable (NaN), never the raw value
passed off as aligned.  They are the pLN lane's, which names cells within a section and imports
no level from anywhere.

Which node is which is the maintainer's record (`ProcessingNotes.txt`), tabulated in
`pln_blocks.csv` and confirmed by them "Correct as shown" on 2026-10-01.  Nothing here infers it:
a `pLN` annotation with no row in that table stops the build, and so does a row whose annotation
is not on the image.

WHY A SIDECAR.  The source workflow's `cell_region` says only `tissue` or `pLN`.  Two images
(6347, 6580) carry one node of each kind, so WHICH annotation holds a cell has to be measured:

    the `pLN` annotation whose outline contains the cell's centroid;
    else, for a cell detected in a node whose centroid falls just outside its outline, the
    nearest `pLN` annotation within EDGE_TOLERANCE_UM;
    else `unassigned` -- out of the alignment, and counted in the manifest.

    python scripts/celltype/alignment_block.py --spec data/cohort-spec-v3.json \\
        --blocks ~/IO60panc2nd/results/v3/pln_blocks.csv \\
        --annotations ~/IO60panc2nd/results/v3/export/annotations_2026-10-01/manual_annotations
    -> data/development/celltype/alignment_block/donor_id=<d>.parquet   (one row per pLN-region cell)
       data/development/celltype/alignment_block/manifest.json

A donor with no `pLN` cells needs no sidecar and gets none.  A donor WITH `pLN` cells and no
sidecar stops every step that asks for its population -- there is no default.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                      # noqa: E402

OUT = panel.OUT / "alignment_block"
PLN_REGION = "pLN"
SAME, DIFFERENT, UNASSIGNED = "same-block", "different-block", "unassigned"

#: The calibration the QuPath export wrote its um centroids in (`setPixelSizeMicrons(0.5078, 0.5078)`
#: in every project script).  The annotation outlines are in pixels, so this is the one conversion.
EXPORT_PIXEL_UM = 0.5078

#: A cell is a child of the annotation it was DETECTED in; its centroid can fall just outside that
#: outline.  Within this distance it takes the nearest node; beyond it nothing is assumed.
EDGE_TOLERANCE_UM = 10.0


def image_tag(image_name: str) -> str:
    """The file stem the v3 export uses for an image."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image_name)


def read_blocks(path: pathlib.Path) -> dict[str, dict[str, str]]:
    """image name -> {pLN annotation id: 'same-block' | 'different-block'}, exactly as recorded."""
    out: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            block = (row.get("block_from_notes") or "").strip()
            if block not in (SAME, DIFFERENT):
                raise SystemExit(f"{path}: {row.get('image')!r} annotation {row.get('object_id')!r} "
                                 f"records block {block!r}; it must be {SAME!r} or {DIFFERENT!r}")
            image = out.setdefault(row["image"], {})
            if row["object_id"] in image:
                raise SystemExit(f"{path}: annotation {row['object_id']} is listed twice")
            image[row["object_id"]] = block
    return out


def pln_annotations(geojson: pathlib.Path) -> list[tuple[str, object]]:
    """(annotation id, outline in pixels) for every `pLN` annotation on one image."""
    import shapely                                  # noqa: PLC0415
    from shapely.geometry import shape              # noqa: PLC0415
    raw = json.loads(geojson.read_text())
    feats = raw["features"] if isinstance(raw, dict) else raw
    out = []
    for f in feats:
        if ((f.get("properties") or {}).get("classification") or {}).get("name") != PLN_REGION:
            continue
        g = shape(f["geometry"])
        if not g.is_valid:
            g = shapely.make_valid(g)
        out.append((str(f.get("id")), g))
    return out


def assign(x_um: np.ndarray, y_um: np.ndarray, annotations: list[tuple[str, object]],
           tolerance_um: float = EDGE_TOLERANCE_UM) -> tuple[np.ndarray, np.ndarray]:
    """Per cell: the id of its `pLN` annotation ('' = none) and how it was found."""
    import shapely                                  # noqa: PLC0415
    from shapely import STRtree                     # noqa: PLC0415
    n = len(x_um)
    ids = np.full(n, "", dtype=object)
    how = np.full(n, "no annotation within tolerance", dtype=object)
    if not annotations or n == 0:
        return ids, how
    geoms = [g for _, g in annotations]
    names = np.array([a for a, _ in annotations], dtype=object)
    points = shapely.points(x_um / EXPORT_PIXEL_UM, y_um / EXPORT_PIXEL_UM)
    tree = STRtree(geoms)
    hit, which = tree.query(points, predicate="within")
    ids[hit] = names[which]
    how[hit] = "centroid inside the outline"
    rest = np.flatnonzero(ids == "")
    if len(rest):
        near = tree.nearest(points[rest])
        dist = shapely.distance(points[rest], np.array(geoms, dtype=object)[near]) * EXPORT_PIXEL_UM
        ok = dist <= tolerance_um
        ids[rest[ok]] = names[near[ok]]
        how[rest[ok]] = f"nearest outline, within {tolerance_um:g} um"
    return ids, how


def build_donor(donor: str, image_name: str, expression: pathlib.Path, geojson: pathlib.Path,
                blocks: dict[str, str]) -> pd.DataFrame | None:
    """The sidecar for one donor, or None when it has no `pLN` cells."""
    t = pq.ParquetFile(expression).read(
        columns=["object_id", "cell_region", "X_centroid", "Y_centroid"]).to_pandas()
    t = t[t["cell_region"] == PLN_REGION]
    annotations = pln_annotations(geojson) if geojson.is_file() else []
    on_image = {a for a, _ in annotations}
    if on_image != set(blocks):
        raise SystemExit(
            f"{image_name}: the pLN annotations on the image and the recorded blocks disagree -- "
            f"on the image only: {sorted(on_image - set(blocks))}; recorded only: "
            f"{sorted(set(blocks) - on_image)}. The block of a node is recorded, never inferred.")
    if t.empty:
        if on_image:
            raise SystemExit(f"{image_name}: {len(on_image)} pLN annotation(s) and no pLN cell in "
                             f"{expression} -- the stage and the export disagree")
        return None
    if not annotations:
        raise SystemExit(f"{image_name}: {len(t):,} pLN cells and no pLN annotation in {geojson}")
    ids, how = assign(t["X_centroid"].to_numpy("float64"), t["Y_centroid"].to_numpy("float64"),
                      annotations)
    block = np.array([blocks.get(a, UNASSIGNED) if a else UNASSIGNED for a in ids], dtype=object)
    return pd.DataFrame({"donor_id": donor, "object_id": t["object_id"].to_numpy(),
                         "pln_annotation_id": ids, "block": block, "found_by": how})


def sidecar_path(donor: str, root: pathlib.Path | None = None) -> pathlib.Path | None:
    """The donor's sidecar under either spelling of its id (0112 == 112), or None."""
    root = OUT if root is None else root
    for spelling in panel._donor_spellings(str(donor)):
        p = root / f"donor_id={spelling}.parquet"
        if p.is_file():
            return p
    return None


def population(donor: str, object_id, cell_region, root: pathlib.Path | None = None) -> np.ndarray:
    """True for the cells of the pancreas block: tissue, and any SAME-BLOCK node's cells."""
    region = np.asarray(cell_region, dtype=object)
    is_pln = region == PLN_REGION
    keep = ~is_pln
    if not is_pln.any():
        return keep
    path = sidecar_path(donor, root)
    if path is None:
        raise SystemExit(
            f"donor {donor}: {int(is_pln.sum()):,} cells lie in a pLN region and there is no "
            f"alignment_block sidecar under {OUT if root is None else root}. UniFORM is done on a "
            "same-block basis (maintainer, 2026-10-02): run scripts/celltype/alignment_block.py "
            "first. There is no default for which block a node belongs to.")
    side = pq.read_table(path, columns=["object_id", "block"]).to_pandas()
    block = pd.Series(side["block"].to_numpy(), index=side["object_id"].to_numpy())
    got = block.reindex(np.asarray(object_id, dtype=object)[is_pln])
    if got.isna().any():
        raise SystemExit(f"donor {donor}: {int(got.isna().sum()):,} pLN cells are not in {path} -- "
                         "the sidecar is older than the expression stage; rebuild it")
    keep[np.flatnonzero(is_pln)[got.to_numpy() == SAME]] = True
    return keep


def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", type=pathlib.Path, required=True,
                    help="cohort specification: which image is which donor")
    ap.add_argument("--blocks", type=pathlib.Path, required=True,
                    help="pln_blocks.csv: the maintainer's record of each node's block")
    ap.add_argument("--annotations", type=pathlib.Path, required=True,
                    help="directory of <image tag>.geojson annotation exports")
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    args = ap.parse_args(argv)

    blocks = read_blocks(args.blocks)
    images = json.loads(args.spec.read_text())["images"]
    declared = {row["image_id"] for row in images}
    stray = sorted(set(blocks) - declared)
    if stray:
        raise SystemExit(f"{args.blocks}: images not in the cohort specification: {stray}")
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in images:
        donor, image = str(row["donor_id"]), row["image_id"]
        geojson = args.annotations / f"{image_tag(image)}.geojson"
        side = build_donor(donor, image, panel.expression_path(donor), geojson, blocks.get(image, {}))
        target = args.out / f"donor_id={donor}.parquet"
        if side is None:
            if target.exists():
                target.unlink()                       # the donor has no pLN cells any more
            continue
        side.to_parquet(target, index=False)
        counts = side["block"].value_counts().to_dict()
        found = side["found_by"].value_counts().to_dict()
        rows.append({"donor_id": donor, "image": image, "pln_cells": int(len(side)),
                     "annotations": {a: blocks[image][a] for a in sorted(blocks[image])},
                     SAME: int(counts.get(SAME, 0)), DIFFERENT: int(counts.get(DIFFERENT, 0)),
                     UNASSIGNED: int(counts.get(UNASSIGNED, 0)), "found_by": found})
        print(f"  {donor}: {len(side):,} pLN cells -> same-block {counts.get(SAME, 0):,}, "
              f"different-block {counts.get(DIFFERENT, 0):,}, unassigned {counts.get(UNASSIGNED, 0):,}",
              flush=True)
    manifest = {
        "ruling": "UniFORM should be done on a same-block basis. (maintainer, 2026-10-02)",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "blocks_file": str(args.blocks), "blocks_sha256": _sha256(args.blocks),
        "spec": str(args.spec), "annotations": str(args.annotations),
        "export_pixel_um": EXPORT_PIXEL_UM, "edge_tolerance_um": EDGE_TOLERANCE_UM,
        "donors_without_pln_cells": sorted(str(r["donor_id"]) for r in images
                                           if str(r["donor_id"]) not in {x["donor_id"] for x in rows}),
        "donors": rows,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"wrote {args.out}  ({len(rows)} donors with pLN cells)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
