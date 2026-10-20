#!/usr/bin/env python
"""Publish PhenoCycler QuPath annotations onto the Xenium grid, for Xenium Explorer.

The serial-section counterpart of ``run_islet_annotation_export.py``.  That driver carries
same-slide post-Xenium IF annotations through a nucleus-matched transform; this one carries the
55-plex PhenoCycler section's drawn anatomy through the S4 serial-section registration in
``data/integration/registration/donor_id=<donor>/roi=panc``.

It is annotation-only and additive.  It never rewrites the registered OME-TIFF, the alignment
files, keypoints, cell groups, or the IF islet sidecar already in the bundle; it stages new
files, validates every one of them against the Explorer 4.1 contract, and only then moves them
into the bundle, rolling back the whole set if any step fails.

What it will not do is claim same-slide accuracy.  The two sections are serial, so these
outlines locate anatomy approximately -- they are not a cell-level correspondence, and the
per-donor QC file records the registration residual alongside every count.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/restain/run_pheno_annotation_export.py --donors 6436 6442
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import pheno_annotations as PA  # noqa: E402

from phenocycler.cohort import ensure_eligible_donors  # noqa: E402
from phenocycler.integration.contract import read_cell_table  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
BUNDLE_ROOT = REPO / "data" / "development" / "restain" / "explorer_bundle"
REGISTRATION_ROOT = CORE / "data" / "integration" / "registration"
CELLS_XEN = CORE / "data" / "integration" / "cells_xen"
COHORT_MANIFEST = CORE / "data" / "qupath_manifest.json"
DEFAULT_ANNOTATION_ROOT = Path.home() / "IO60panc2nd" / "results" / "Geojson"

QC_FILENAME = "pheno_annotations_qc.json"


class PublishError(RuntimeError):
    """A bundle cannot safely receive the PhenoCycler annotation sidecars."""


def _find_bundle(donor: str) -> Path:
    matches = sorted(BUNDLE_ROOT.glob(f"donor_{donor}_*"))
    if len(matches) != 1:
        raise PublishError(
            f"donor {donor}: expected exactly 1 Explorer bundle under {BUNDLE_ROOT}, "
            f"found {len(matches)}"
        )
    return matches[0]


def _load_morphology(bundle: Path) -> dict[str, Any]:
    manifest_path = bundle / "manifest.json"
    if not manifest_path.exists():
        raise PublishError(f"{bundle} has no manifest.json; it is not a rendered bundle")
    manifest = json.loads(manifest_path.read_text())
    morphology = manifest.get("morphology")
    if not isinstance(morphology, dict):
        raise PublishError(f"{manifest_path} carries no morphology block")
    for key in ("width", "height", "um_per_px"):
        if key not in morphology:
            raise PublishError(f"{manifest_path} morphology block has no {key!r}")
    return morphology


def _load_transform(donor: str, *, use_displacement: bool) -> tuple[Transform, dict[str, Any]]:
    directory = REGISTRATION_ROOT / f"donor_id={donor}" / "roi=panc"
    if not (directory / "transform.json").exists():
        raise PublishError(
            f"donor {donor}: no S4 registration at {directory}. Run "
            f"`python -m phenocycler.integration.pipeline --only register --donor {donor}` first."
        )
    full = Transform.load(directory, stem="transform")
    qc = json.loads((directory / "qc.json").read_text())
    scale = qc.get("scale_report", {}) or {}
    evidence = {
        "registration_dir": str(directory),
        "track": qc.get("track"),
        # WHICH transform placed these annotations, not just which directory it came from.
        # Donor 6436 shipped for a month under a transform that was 12.5 deg and 2.6 mm out,
        # and nothing in the published bundle recorded enough to notice. A directory path is
        # not provenance -- its contents change.
        "transform": {
            "rotation_deg": full.rotation_deg,
            "scale": full.scale,
            "scale_x": scale.get("scale_x"),
            "scale_y": scale.get("scale_y"),
            "anisotropy": scale.get("anisotropy"),
            "translation_um": list(full.translation_um),
            "affine": full.affine.tolist(),
        },
        # The out-of-sample evidence the verdict rests on. tissue_dice is carried for
        # continuity and is deliberately not part of that verdict: it moves against
        # correctness when one section is a sub-field of the other.
        "footprint_ncc": qc.get("footprint_ncc"),
        "footprint_ncc_z": qc.get("footprint_ncc_z"),
        "footprint_containment": qc.get("footprint_containment"),
        "islet_peak_ratio": qc.get("islet_peak_ratio"),
        "islet_peak_offset_um": qc.get("islet_peak_offset_um"),
        "tissue_dice": qc.get("tissue_dice"),
        "islet_dice": qc.get("islet_dice"),
        "rmse_um": qc.get("rmse_um"),
        "score": qc.get("score"),
        "stages": qc.get("stages"),
        "notes": qc.get("notes"),
        "mirrored": qc.get("mirrored"),
        "transform_has_displacement": full.displacement is not None,
        "transform_max_displacement_um": float(full.max_displacement_um),
        "displacement_report": qc.get("displacement_report"),
    }
    if use_displacement:
        return full, {**evidence, "displacement_used": True}
    # The non-rigid stage is off for this lane (reg_nonrigid = false) because its field was not
    # invertible where it mattered: measured on the real annotation vertices it failed the
    # round trip on 15.0% (6436) and 6.1% (6442) of them, costing 90 and 48 publishable islets.
    # If a field is present anyway it is dropped here and the fact recorded, rather than being
    # quietly applied.
    return replace(full, displacement=None), {**evidence, "displacement_used": False}


def _export_donor(
    donor: str,
    *,
    annotation_root: Path,
    classes: tuple[PA.PhenoClass, ...],
    min_diameter_um: float,
    use_displacement: bool,
    max_features: int,
    dry_run: bool,
    outline_pixel_um: float = PA.OUTLINE_PIXEL_UM,
) -> dict[str, Any]:
    bundle = _find_bundle(donor)
    serial = bundle.name.rsplit("_", 1)[-1]
    morphology = _load_morphology(bundle)
    width = int(morphology["width"])
    height = int(morphology["height"])
    output_um_per_px = float(morphology["um_per_px"])

    manifest = json.loads(COHORT_MANIFEST.read_text())
    stem = PA.image_stem_for(manifest, donor)
    input_um_per_px = PA.pixel_size_for(manifest, donor)
    transform, evidence = _load_transform(donor, use_displacement=use_displacement)

    # Where the Xenium section actually has data. Annotations are selected against this rather
    # than against the morphology image's bounding rectangle, which is 38% larger and mostly
    # blank -- see PA.xenium_tissue_outline.
    xen_xy = read_cell_table(CELLS_XEN, donor, "panc", modality="xenium").xy()
    tissue_outline, outline_provenance = PA.xenium_tissue_outline(
        xen_xy, pixel_um=outline_pixel_um)
    print(f"  [{donor}] tissue outline: {outline_provenance['area_um2'] / 1e6:.1f} mm2 from "
          f"{outline_provenance['n_cells']:,} cells", flush=True)

    staging = bundle / ".pheno_annotations_staging"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    report: dict[str, Any] = {
        "status": "generated",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "donor_id": donor,
        "xenium_serial": serial,
        "bundle": str(bundle),
        "section_relationship": "serial",
        "same_slide": False,
        "source_root": str(annotation_root),
        "source_image_stem": stem,
        "source_um_per_px": input_um_per_px,
        "morphology": {
            "width_px": width,
            "height_px": height,
            "um_per_px": output_um_per_px,
            "field_mm": [round(width * output_um_per_px / 1000, 3),
                         round(height * output_um_per_px / 1000, 3)],
        },
        "registration": evidence,
        "tissue_outline": outline_provenance,
        "min_diameter_um": float(min_diameter_um),
        "max_features_per_file": int(max_features),
        "classes": [],
    }

    try:
        for cls in classes:
            source = PA.source_path(annotation_root, cls, stem)
            print(f"  [{donor}] {cls.name:<14} <- {source.name}", flush=True)
            result = PA.export_class(
                source,
                staging,
                transform,
                cls,
                serial=serial,
                input_um_per_px=input_um_per_px,
                output_um_per_px=output_um_per_px,
                width=width,
                height=height,
                min_diameter_um=min_diameter_um,
                max_features=max_features,
                used_displacement=use_displacement,
                tissue_outline=tissue_outline,
            )
            report["classes"].append(result)
            audit = result["placement_audit"]
            print(
                f"      {result['n_features_total']:>7,} polygons in {result['n_files']} file(s)"
                f"  | outside the Xenium field: "
                f"{audit['excluded_outside_xenium_field']:,}"
                f"  | off the tissue outline: "
                f"{audit['excluded_outside_tissue_outline']:,}"
                f"  | below {min_diameter_um:g} um: {audit['excluded_below_minimum_size']:,}",
                flush=True,
            )

        report["n_classes"] = len(report["classes"])
        report["n_files_total"] = sum(c["n_files"] for c in report["classes"])
        report["n_features_total"] = sum(c["n_features_total"] for c in report["classes"])

        if dry_run:
            print(f"  [{donor}] dry run: staged {report['n_files_total']} file(s), publishing "
                  f"nothing", flush=True)
            report["status"] = "dry_run"
            return report

        destination = bundle / "annotations"
        destination.mkdir(parents=True, exist_ok=True)

        # Retire this lane's previous output for the classes being republished. A class that
        # splits into a different number of parts than last time would otherwise leave stale
        # part files behind -- still valid GeoJSON, still importable, but built from a
        # superseded transform. They are moved aside, not deleted, until the swap succeeds.
        retired = staging / "retired"
        retired.mkdir(parents=True, exist_ok=True)
        stale: list[tuple[Path, Path]] = []
        for cls in classes:
            pattern = f"{serial}_pheno_{cls.name}_on_xenium_grid*.geojson"
            for old_file in sorted(destination.glob(pattern)):
                aside = retired / old_file.name
                shutil.move(str(old_file), str(aside))
                stale.append((old_file, aside))

        moved: list[Path] = []
        try:
            for staged in sorted(staging.glob("*.geojson")):
                target = destination / staged.name
                shutil.move(str(staged), str(target))
                moved.append(target)
            report["replaced_files"] = [original.name for original, _ in stale]
            (bundle / QC_FILENAME).write_text(json.dumps(report, indent=2) + "\n")
        except Exception:
            for target in moved:
                target.unlink(missing_ok=True)
            for original, aside in stale:                 # put the previous output back
                if aside.exists():
                    shutil.move(str(aside), str(original))
            raise
        report["published_files"] = [str(p.relative_to(bundle)) for p in moved]
        return report
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Publish PhenoCycler QuPath annotations on the Xenium grid for Explorer"
    )
    ap.add_argument("--donors", nargs="+", required=True, help="donor ids to export")
    ap.add_argument("--annotation-root", type=Path, default=DEFAULT_ANNOTATION_ROOT,
                    help="QuPath export tree, one directory per drawn class")
    ap.add_argument("--classes", nargs="+", default=[c.name for c in PA.PHENO_CLASSES],
                    help="QuPath classification names to publish")
    ap.add_argument("--min-diameter-um", type=float, default=0.0,
                    help="drop objects below this equivalent diameter (0 keeps every object)")
    ap.add_argument("--max-features-per-file", type=int, default=PA.MAX_FEATURES_PER_FILE,
                    help="Explorer's per-file polygon cap; a class over it is split into parts")
    ap.add_argument("--use-displacement", action="store_true",
                    help="apply the S4 B-spline field. Its inverse does not converge on the "
                         "registered donors; affine-only is the default for that reason")
    ap.add_argument("--dry-run", action="store_true",
                    help="build and validate everything, publish nothing")
    args = ap.parse_args(argv)

    unknown = [c for c in args.classes if c not in PA.CLASSES_BY_NAME]
    if unknown:
        ap.error(f"unknown class(es) {unknown}; known: {sorted(PA.CLASSES_BY_NAME)}")
    classes = tuple(PA.CLASSES_BY_NAME[c] for c in args.classes)

    # An explicit request for a centrally excluded donor fails rather than being filtered
    # away silently, matching the rest of the lane.
    try:
        donors = ensure_eligible_donors(args.donors, context="PhenoCycler annotation export")
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    failures = 0
    for donor in donors:
        print(f"[pheno-annotations] donor {donor}", flush=True)
        try:
            report = _export_donor(
                donor,
                annotation_root=args.annotation_root,
                classes=classes,
                min_diameter_um=args.min_diameter_um,
                use_displacement=args.use_displacement,
                max_features=args.max_features_per_file,
                dry_run=args.dry_run,
            )
        except (PublishError, PA.PhenoAnnotationError) as exc:
            print(f"  [{donor}] FAILED: {exc}", file=sys.stderr, flush=True)
            failures += 1
            continue
        print(
            f"  [{donor}] {report['status']}: {report['n_features_total']:,} polygons, "
            f"{report['n_files_total']} file(s), {report['n_classes']} class(es)",
            flush=True,
        )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
