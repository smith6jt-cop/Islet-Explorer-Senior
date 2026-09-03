#!/usr/bin/env python
"""Fit one low-dimensional DAPI correction and reserve half the tissue for its gate.

The nucleus matcher can converge to neighbouring nuclei while reporting a deceptively small
nearest-neighbour error.  ``run_align.py`` exposes that failure as a coherent field of DAPI
offsets.  This stage may use one origin-anchored checkerboard fold of those offsets to fit a
single affine residual field; the other fold is never read by the fitter and is reserved for
an image-space gate::

    python scripts/restain/run_align.py --donors 115 --transform refined
    python scripts/restain/run_intensity_refine.py --donors 115
    python scripts/restain/run_align.py --donors 115 \
        --transform intensity_refined --window-fold odd

The output is a new, explicitly named transform.  It never overwrites ``refined`` and it is
not exportable until the matching held-out gate says PASS.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import HuberRegressor

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import align as AL  # noqa: E402
import registration_inputs as RI  # noqa: E402
import run_align as RA  # noqa: E402
from run_register import OUT  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

FIT_COORD_SCALE_UM = 5_000.0
FIELD_PIXEL_UM = 250.0
HUBER_EPSILON = 1.35
MIN_TRAIN_WINDOWS = 30
MAX_CORRECTION_GRADIENT = 0.02
MIN_CORRECTION_DETERMINANT = 0.95
MAX_CORRECTION_DETERMINANT = 1.05
MAX_CORRECTION_UM = AL.SEARCH_RADIUS_UM


class IntensityRefineError(RuntimeError):
    """The correction cannot be fitted under the declared hold-out contract."""


def _safe_stem(stem: str, *, role: str) -> str:
    """Require one filename token before using a caller-supplied transform stem."""
    value = str(stem)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value) is None:
        raise IntensityRefineError(f"unsafe {role} transform stem {value!r}")
    return value


def _same_path(left: Path | str, right: Path | str) -> bool:
    """Compare declared paths without requiring their targets to exist."""
    return Path(left).expanduser().resolve(strict=False) == Path(right).expanduser().resolve(
        strict=False)


def _verify_registration_input(gate: dict, gate_path: Path,
                               record: RI.RegistrationInput) -> None:
    """Bind a source gate to the exact inventory row selected by the batch CLI."""
    expected_identity = {
        "donor_id": record.donor_id,
        "serial": record.serial,
        "roi": record.roi,
    }
    for key, expected in expected_identity.items():
        if str(gate.get(key, "")) != str(expected):
            raise IntensityRefineError(
                f"{gate_path}: {key}={gate.get(key)!r}, expected {expected!r}")

    registration_input = gate.get("registration_input")
    if not isinstance(registration_input, dict):
        raise IntensityRefineError(
            f"{gate_path}: missing object-valued registration_input")
    for key, expected in (("qptiff", record.qptiff), ("bundle", record.bundle)):
        declared = registration_input.get(key)
        if not declared or not _same_path(declared, expected):
            raise IntensityRefineError(
                f"{gate_path}: registration_input.{key}={declared!r}, expected {expected}")

    expected_mapping = {
        "pixel_um": float(record.pixel_um),
        "pixel_size_source": "qptiff XML PixelSizeMicrons",
        "dapi_page": int(record.dapi_page),
        "dapi_channel_name": record.channel_names[record.dapi_page],
    }
    for key, expected in expected_mapping.items():
        if registration_input.get(key) != expected:
            raise IntensityRefineError(
                f"{gate_path}: registration_input.{key}="
                f"{registration_input.get(key)!r}, expected {expected!r}")

    certified = gate.get("alignment", {}).get("certifies_mapping")
    if not isinstance(certified, dict):
        raise IntensityRefineError(
            f"{gate_path}: missing alignment.certifies_mapping")
    certified_expected = {
        "pixel_um_full_res": float(record.pixel_um),
        "source": "qptiff XML PixelSizeMicrons",
        "dapi_page": int(record.dapi_page),
        "dapi_channel_name": record.channel_names[record.dapi_page],
    }
    for key, expected in certified_expected.items():
        if certified.get(key) != expected:
            raise IntensityRefineError(
                f"{gate_path}: alignment.certifies_mapping.{key}="
                f"{certified.get(key)!r}, expected {expected!r}")


def _atomic_json(path: Path, value: dict) -> None:
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False))
    tmp.replace(path)


def _fold(cx_um: float, cy_um: float) -> str:
    ix = int(np.floor(float(cx_um) / RA.WINDOW_STEP_UM))
    iy = int(np.floor(float(cy_um) / RA.WINDOW_STEP_UM))
    return "even" if (ix + iy) % 2 == 0 else "odd"


def _verified_source(donor_dir: Path, stem: str, *,
                     record: RI.RegistrationInput) -> tuple[Transform, dict, Path]:
    gate_path = donor_dir / f"alignment_{stem}.json"
    if not gate_path.exists():
        raise IntensityRefineError(
            f"{gate_path} is absent -- first gate {stem!r} on all occupied windows")
    gate = json.loads(gate_path.read_text())
    a = gate.get("alignment", {})
    if a.get("transform_source") != stem:
        raise IntensityRefineError(
            f"{gate_path} measured {a.get('transform_source')!r}, not {stem!r}")
    if a.get("window_set") != "occupied" or a.get("window_fold", "all") != "all":
        raise IntensityRefineError(
            f"{gate_path} must cover all occupied windows before a fold can be reserved")
    _verify_registration_input(gate, gate_path, record)

    transform_path = donor_dir / f"{stem}.json"
    if not transform_path.exists():
        raise IntensityRefineError(f"source transform is absent: {transform_path}")
    actual = RA._sha256(transform_path)
    if a.get("transform_sha256") != actual:
        raise IntensityRefineError(
            f"{gate_path} does not hash the current {transform_path.name}; rerun run_align.py")
    payload = json.loads(transform_path.read_text())
    displacement_path = donor_dir / f"{stem}_displacement.npz"
    if payload.get("has_displacement"):
        if not displacement_path.exists():
            raise IntensityRefineError(
                f"{transform_path.name} declares a field but {displacement_path.name} is absent")
        if a.get("displacement_sha256") != RA._sha256(displacement_path):
            raise IntensityRefineError(
                f"{gate_path} does not hash the current displacement field; rerun run_align.py")
    elif a.get("displacement_sha256") is not None:
        raise IntensityRefineError(f"{gate_path} hashes a field the transform does not use")
    return Transform.load(donor_dir, stem=stem), gate, gate_path


def _training_windows(gate: dict, fold: str) -> list[dict]:
    out = []
    for w in gate["alignment"].get("windows", []):
        if _fold(w["cx_um"], w["cy_um"]) != fold:
            continue
        # The peak must be demonstrably useful: moving by it increases DAPI contrast.  This
        # filter sees only the training fold, never the fold used by the release gate.
        if (w.get("reason") != "measured" or w.get("truncated")
                or w.get("contrast") is None or w.get("contrast_shifted") is None
                or not w["contrast_shifted"] > w["contrast"]):
            continue
        out.append(w)
    if len(out) < MIN_TRAIN_WINDOWS:
        raise IntensityRefineError(
            f"only {len(out)} reliable {fold}-fold windows; need {MIN_TRAIN_WINDOWS}")
    return out


def _fit_correction(windows: list[dict]) -> dict:
    q = np.asarray([[w["cx_um"], w["cy_um"]] for w in windows], dtype=np.float64)
    d = np.asarray([[w["dx_um"], w["dy_um"]] for w in windows], dtype=np.float64)
    center = q.mean(axis=0)
    x = (q - center) / FIT_COORD_SCALE_UM
    models = []
    for j in range(2):
        model = HuberRegressor(epsilon=HUBER_EPSILON, alpha=0.0, max_iter=1_000)
        model.fit(x, d[:, j])
        models.append(model)
    pred = np.column_stack([m.predict(x) for m in models])
    residual = np.linalg.norm(d - pred, axis=1)
    baseline = np.linalg.norm(d, axis=1)
    gradient = np.vstack([m.coef_ for m in models]) / FIT_COORD_SCALE_UM
    correction_linear = np.eye(2) + gradient
    determinant = float(np.linalg.det(correction_linear))
    if float(np.max(np.abs(gradient))) > MAX_CORRECTION_GRADIENT:
        raise IntensityRefineError(
            f"correction gradient {gradient.tolist()} exceeds {MAX_CORRECTION_GRADIENT}")
    if not MIN_CORRECTION_DETERMINANT <= determinant <= MAX_CORRECTION_DETERMINANT:
        raise IntensityRefineError(
            f"correction determinant {determinant:.6f} is outside "
            f"[{MIN_CORRECTION_DETERMINANT}, {MAX_CORRECTION_DETERMINANT}]")
    if float(np.median(residual)) >= float(np.median(baseline)):
        raise IntensityRefineError("the robust affine correction does not improve its training fold")
    return {
        "center_um": center,
        "models": models,
        "gradient": gradient,
        "correction_determinant": determinant,
        "baseline_median_um": float(np.median(baseline)),
        "residual_median_um": float(np.median(residual)),
        "residual_p90_um": float(np.percentile(residual, 90)),
    }


def _predict(fit: dict, xy_um: np.ndarray) -> np.ndarray:
    x = (np.asarray(xy_um, dtype=np.float64) - fit["center_um"]) / FIT_COORD_SCALE_UM
    return np.column_stack([m.predict(x) for m in fit["models"]])


def _field_grid(source: Transform, windows: list[dict]) -> tuple[GridSpec, np.ndarray]:
    if source.displacement is not None:
        return source.grid, np.asarray(source.displacement, dtype=np.float32)
    q = np.asarray([[w["cx_um"], w["cy_um"]] for w in windows], dtype=np.float64)
    margin = RA.WINDOW_UM
    x0 = np.floor((q[:, 0].min() - margin) / FIELD_PIXEL_UM) * FIELD_PIXEL_UM
    y0 = np.floor((q[:, 1].min() - margin) / FIELD_PIXEL_UM) * FIELD_PIXEL_UM
    x1 = np.ceil((q[:, 0].max() + margin) / FIELD_PIXEL_UM) * FIELD_PIXEL_UM
    y1 = np.ceil((q[:, 1].max() + margin) / FIELD_PIXEL_UM) * FIELD_PIXEL_UM
    grid = GridSpec(float(x0), float(y0), FIELD_PIXEL_UM,
                    int(round((y1 - y0) / FIELD_PIXEL_UM)) + 1,
                    int(round((x1 - x0) / FIELD_PIXEL_UM)) + 1)
    return grid, np.zeros((2, grid.height, grid.width), dtype=np.float32)


def _corrected_transform(source: Transform, fit: dict, all_windows: list[dict]) -> Transform:
    grid, old_field = _field_grid(source, all_windows)
    rows, cols = np.mgrid[:grid.height, :grid.width]
    u = grid.to_microns(np.column_stack([cols.ravel(), rows.ravel()]))
    old = np.column_stack([old_field[0].ravel(), old_field[1].ravel()])
    # Transform.apply computes q = u + D(u).  The new field is therefore
    # D_new(u) = D(u) + correction(q), preserving the old transform exactly before adding the
    # low-dimensional intensity correction.
    new = old + _predict(fit, u + old)
    max_correction = float(np.max(np.linalg.norm(new - old, axis=1)))
    if max_correction > MAX_CORRECTION_UM:
        raise IntensityRefineError(
            f"predicted correction reaches {max_correction:.2f} um, above the "
            f"{MAX_CORRECTION_UM:.2f} um search radius")
    field = np.stack([new[:, 0].reshape(grid.shape), new[:, 1].reshape(grid.shape)]).astype(
        np.float32)
    meta = dict(source.meta)
    meta.update({
        "source": "one-fold DAPI window correction",
        "correction_family": "robust affine residual represented on a 250-um field",
        "training_fold": fit["training_fold"],
        "withheld_gate_fold": fit["withheld_fold"],
        "max_added_correction_um": max_correction,
    })
    return Transform(affine=np.array(source.affine), displacement=field, grid=grid, meta=meta)


def _save_transform_atomic(tf: Transform, donor_dir: Path, stem: str) -> tuple[Path, Path]:
    tmp_stem = f".{stem}.tmp-{os.getpid()}"
    tmp_json = tf.save(donor_dir, stem=tmp_stem)
    tmp_npz = donor_dir / f"{tmp_stem}_displacement.npz"
    dst_json = donor_dir / f"{stem}.json"
    dst_npz = donor_dir / f"{stem}_displacement.npz"
    tmp_npz.replace(dst_npz)
    tmp_json.replace(dst_json)  # JSON last: readers never see a declaration without its field.
    return dst_json, dst_npz


def refine_one(record: RI.RegistrationInput, donor_dir: Path, *,
               source_stem: str, output_stem: str,
               training_fold: str, overwrite: bool = False) -> dict:
    donor_dir = Path(donor_dir)
    source_stem = _safe_stem(source_stem, role="source")
    output_stem = _safe_stem(output_stem, role="output")
    if source_stem == output_stem:
        raise IntensityRefineError("source and output stems must differ")
    output_path = donor_dir / f"{output_stem}.json"
    if output_path.exists() and not overwrite:
        raise IntensityRefineError(f"{output_path} exists; pass --overwrite to replace it")

    source, gate, gate_path = _verified_source(
        donor_dir, source_stem, record=record)
    train = _training_windows(gate, training_fold)
    fit = _fit_correction(train)
    fit["training_fold"] = training_fold
    fit["withheld_fold"] = "odd" if training_fold == "even" else "even"
    corrected = _corrected_transform(source, fit, gate["alignment"]["windows"])
    transform_path, displacement_path = _save_transform_atomic(
        corrected, donor_dir, output_stem)

    models = fit.pop("models")
    report = {
        "donor_id": gate["donor_id"], "serial": gate["serial"], "roi": gate["roi"],
        "source": {
            "stem": source_stem, "gate_file": str(gate_path),
            "gate_sha256": RA._sha256(gate_path),
            "transform_sha256": gate["alignment"]["transform_sha256"],
            "displacement_sha256": gate["alignment"].get("displacement_sha256"),
            "verdict": gate["alignment"]["verdict"],
            "median_offset_um": gate["alignment"]["median_offset_um"],
        },
        "fit": {
            "family": "Huber affine residual", "epsilon": HUBER_EPSILON,
            "coordinate_scale_um": FIT_COORD_SCALE_UM,
            "training_fold": training_fold, "withheld_gate_fold": fit["withheld_fold"],
            "n_training_windows": len(train),
            "selection": "measured, non-truncated, and DAPI contrast improves after shift",
            "center_um": fit["center_um"].tolist(),
            "intercept_um": [float(m.intercept_) for m in models],
            "coefficient_per_scaled_coordinate": [m.coef_.tolist() for m in models],
            "gradient_um_per_um": fit["gradient"].tolist(),
            "correction_determinant": fit["correction_determinant"],
            "training_baseline_median_um": fit["baseline_median_um"],
            "training_residual_median_um": fit["residual_median_um"],
            "training_residual_p90_um": fit["residual_p90_um"],
        },
        "output": {
            "stem": output_stem, "transform_file": str(transform_path),
            "transform_sha256": RA._sha256(transform_path),
            "displacement_file": str(displacement_path),
            "displacement_sha256": RA._sha256(displacement_path),
            "max_displacement_um": corrected.max_displacement_um,
            "required_gate_command": (
                f"run_align.py --donors {gate['donor_id']} --transform {output_stem} "
                f"--window-fold {fit['withheld_fold']}")
        },
        "constants": {
            "MIN_TRAIN_WINDOWS": MIN_TRAIN_WINDOWS,
            "MAX_CORRECTION_GRADIENT": MAX_CORRECTION_GRADIENT,
            "MIN_CORRECTION_DETERMINANT": MIN_CORRECTION_DETERMINANT,
            "MAX_CORRECTION_DETERMINANT": MAX_CORRECTION_DETERMINANT,
            "MAX_CORRECTION_UM": MAX_CORRECTION_UM,
            "FIELD_PIXEL_UM": FIELD_PIXEL_UM,
        },
    }
    report_path = donor_dir / f"intensity_refinement_{output_stem}.json"
    _atomic_json(report_path, report)
    report["report_file"] = str(report_path)
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--source-stem", default="refined")
    ap.add_argument("--output-stem", default="intensity_refined")
    ap.add_argument("--training-fold", default="even", choices=["even", "odd"])
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--inventory", type=Path, default=RI.REGISTRATION_INPUTS)
    ap.add_argument("--images-dir", type=Path, default=RI.IMAGES)
    ap.add_argument("--geojson-dir", type=Path, default=RI.GEOJSON_DIR)
    ap.add_argument("--xenium-paths", type=Path, default=RI.XENIUM_PATHS)
    ap.add_argument("--xenium-root", type=Path, default=RI.XENIUM_ROOT)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    records = RI.load_registration_inputs(
        args.inventory, args.xenium_paths, args.images_dir, args.geojson_dir,
        args.xenium_root, roi=args.roi)
    selected = RI.select_registration_inputs(records, args.donors,
                                              context="run_intensity_refine")
    for record in selected:
        RI.require_registration_files(record)
    rows = []
    for record in selected:
        row = refine_one(record, args.out / "valis" / f"donor_id={record.donor_id}",
                         source_stem=args.source_stem, output_stem=args.output_stem,
                         training_fold=args.training_fold, overwrite=args.overwrite)
        rows.append(row)
        f = row["fit"]
        print(f"[{record.donor_id}] fitted {f['n_training_windows']} "
              f"{f['training_fold']}-fold windows: median "
              f"{f['training_baseline_median_um']:.2f} -> "
              f"{f['training_residual_median_um']:.2f} um; "
              f"{f['withheld_gate_fold']} fold remains held out", flush=True)
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    summary = args.out / f"intensity_refinement_{args.roi}_{args.output_stem}.json"
    _atomic_json(summary, {"rows": rows})
    print(f"-> {summary}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
