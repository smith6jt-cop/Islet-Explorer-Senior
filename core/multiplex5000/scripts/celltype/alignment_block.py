"""Historical CLI for the packaged same-block assignment method.

The scientific decision is implemented once in
``multiplex5000.normalization_block``. This command keeps the old lane's
output layout and reads its expression partitions through ``panel``.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from datetime import datetime, timezone

from multiplex5000 import normalization_block as NB
from multiplex5000.artifacts import sha256_file

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel  # noqa: E402

OUT = panel.OUT / "alignment_block"
PLN_REGION = NB.PLN_REGION
SAME, DIFFERENT, UNASSIGNED = NB.SAME, NB.DIFFERENT, NB.UNASSIGNED
EXPORT_PIXEL_UM = NB.EXPORT_PIXEL_UM
EDGE_TOLERANCE_UM = NB.EDGE_TOLERANCE_UM
image_tag = NB.image_tag
read_blocks = NB.read_blocks
pln_annotations = NB.pln_annotations
assign = NB.assign
build_donor = NB.build_donor


def sidecar_path(donor: str, root: pathlib.Path | None = None) -> pathlib.Path | None:
    return NB.sidecar_path(donor, OUT if root is None else root)


def population(donor: str, object_id, cell_region,
               root: pathlib.Path | None = None):
    return NB.population(donor, object_id, cell_region, OUT if root is None else root)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=NB.__doc__,
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
                target.unlink()
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
        "blocks_file": str(args.blocks), "blocks_sha256": sha256_file(args.blocks),
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
