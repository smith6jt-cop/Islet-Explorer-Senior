"""Spatial hold-out and provenance contracts for DAPI intensity refinement."""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import registration_inputs as RI  # noqa: E402
import run_align as RA  # noqa: E402
import run_intensity_refine as IR  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402


def _record(donor_dir: Path) -> RI.RegistrationInput:
    return RI.RegistrationInput(
        donor_id="115", serial="0059902", roi="panc",
        qptiff=donor_dir / "source.qptiff", bundle=donor_dir / "xenium_bundle",
        geojson=None, pixel_um=0.25, dapi_page=0,
        channel_names=("DAPI", "Opal570"), pyramid_shapes=((64, 48), (32, 24)),
        provenance_source="test fixture",
    )


def _source(donor_dir: Path, *, odd_offset: float = 0.0,
            with_displacement: bool = False) -> RI.RegistrationInput:
    donor_dir.mkdir(parents=True)
    record = _record(donor_dir)
    grid = GridSpec(-500.0, -500.0, 250.0, 22, 22) if with_displacement else None
    displacement = None
    if grid is not None:
        displacement = np.zeros((2, grid.height, grid.width), dtype=np.float32)
        displacement[0] = 3.0
        displacement[1] = -4.0
    transform_path = Transform(
        affine=np.array([[1.0, 0.0, 10.0], [0.0, 1.0, -20.0]]),
        displacement=displacement, grid=grid).save(donor_dir, stem="refined")
    displacement_path = donor_dir / "refined_displacement.npz"
    windows = []
    for iy in range(8):
        for ix in range(8):
            x, y = (ix + 0.5) * 500.0, (iy + 0.5) * 500.0
            # A small affine residual that the production family can represent exactly.
            dx = 2.0 + 0.001 * x - 0.0005 * y
            dy = -1.0 + 0.0002 * x + 0.0008 * y
            if (ix + iy) % 2:
                dx += odd_offset
                dy -= odd_offset
            windows.append({
                "cx_um": x, "cy_um": y, "dx_um": dx, "dy_um": dy,
                "reason": "measured", "truncated": False,
                "contrast": 1.2, "contrast_shifted": 2.0,
            })
    gate = {
        "donor_id": "115", "serial": "0059902", "roi": "panc",
        "registration_input": {
            "qptiff": str(record.qptiff), "bundle": str(record.bundle),
            "pixel_um": 0.25, "pixel_size_source": "qptiff XML PixelSizeMicrons",
            "dapi_page": 0, "dapi_channel_name": "DAPI",
        },
        "alignment": {
            "transform_source": "refined", "window_set": "occupied",
            "window_fold": "all", "transform_sha256": RA._sha256(transform_path),
            "displacement_sha256": (
                RA._sha256(displacement_path) if with_displacement else None),
            "certifies_mapping": {
                "pixel_um_full_res": 0.25,
                "source": "qptiff XML PixelSizeMicrons",
                "dapi_page": 0, "dapi_channel_name": "DAPI",
            },
            "verdict": "FAIL",
            "median_offset_um": 5.0, "windows": windows,
        },
    }
    (donor_dir / "alignment_refined.json").write_text(json.dumps(gate))
    return record


def test_fit_uses_only_declared_training_fold(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    record_a = _source(a, odd_offset=0.0)
    record_b = _source(b, odd_offset=1_000.0)

    ra = IR.refine_one(record_a, a, source_stem="refined", output_stem="intensity_refined",
                       training_fold="even")
    IR.refine_one(record_b, b, source_stem="refined", output_stem="intensity_refined",
                  training_fold="even")

    assert ra["fit"]["n_training_windows"] == 32
    assert ra["fit"]["withheld_gate_fold"] == "odd"
    np.testing.assert_allclose(ra["fit"]["gradient_um_per_um"],
                               [[0.001, -0.0005], [0.0002, 0.0008]], atol=1e-8)
    # Changing every held-out offset by a kilometre-equivalent does not change one output
    # field byte: the fitter never looked at that fold.
    za = np.load(a / "intensity_refined_displacement.npz")["displacement"]
    zb = np.load(b / "intensity_refined_displacement.npz")["displacement"]
    np.testing.assert_array_equal(za, zb)
    assert Transform.load(a, stem="intensity_refined").displacement is not None


def test_source_gate_must_hash_exact_transform(tmp_path):
    donor_dir = tmp_path / "donor"
    record = _source(donor_dir)
    gate_path = donor_dir / "alignment_refined.json"
    gate = json.loads(gate_path.read_text())
    gate["alignment"]["transform_sha256"] = "0" * 64
    gate_path.write_text(json.dumps(gate))
    with pytest.raises(IR.IntensityRefineError, match="does not hash the current"):
        IR.refine_one(record, donor_dir, source_stem="refined",
                      output_stem="intensity_refined",
                      training_fold="even")


def test_correction_is_added_in_fixed_space_after_existing_displacement(tmp_path):
    donor_dir = tmp_path / "donor"
    record = _source(donor_dir, with_displacement=True)
    IR.refine_one(record, donor_dir, source_stem="refined",
                  output_stem="intensity_refined", training_fold="even")
    source = Transform.load(donor_dir, stem="refined")
    corrected = Transform.load(donor_dir, stem="intensity_refined")
    moving = np.array([[200.0, 300.0], [1_200.0, 1_700.0], [3_500.0, 3_000.0]])
    fixed = source.apply(moving)
    expected_delta = np.column_stack([
        2.0 + 0.001 * fixed[:, 0] - 0.0005 * fixed[:, 1],
        -1.0 + 0.0002 * fixed[:, 0] + 0.0008 * fixed[:, 1],
    ])
    np.testing.assert_allclose(corrected.apply(moving), fixed + expected_delta, atol=1e-5)


@pytest.mark.parametrize(
    ("changes", "match"),
    [
        ({"donor_id": "999"}, "donor_id"),
        ({"serial": "0099999"}, "serial"),
        ({"roi": "other"}, "roi"),
        ({"qptiff": Path("/wrong/source.qptiff")}, "registration_input.qptiff"),
        ({"bundle": Path("/wrong/xenium_bundle")}, "registration_input.bundle"),
        ({"pixel_um": 0.5}, "registration_input.pixel_um"),
        ({"dapi_page": 1}, "registration_input.dapi_page"),
        ({"channel_names": ("NOT_DAPI", "Opal570")},
         "registration_input.dapi_channel_name"),
    ],
)
def test_selected_inventory_record_must_match_source_gate(tmp_path, changes, match):
    donor_dir = tmp_path / "donor"
    record = _source(donor_dir)
    wrong = replace(record, **changes)
    with pytest.raises(IR.IntensityRefineError, match=match):
        IR.refine_one(wrong, donor_dir, source_stem="refined",
                      output_stem="intensity_refined", training_fold="even")
    assert not (donor_dir / "intensity_refined.json").exists()


@pytest.mark.parametrize(
    ("source_stem", "output_stem", "match"),
    [
        ("../refined", "intensity_refined", "unsafe source"),
        ("refined", "../intensity_refined", "unsafe output"),
        ("refined", ".hidden", "unsafe output"),
        ("refined", "nested/name", "unsafe output"),
    ],
)
def test_transform_stems_are_single_safe_tokens(
        tmp_path, source_stem, output_stem, match):
    donor_dir = tmp_path / "donor"
    record = _source(donor_dir)
    with pytest.raises(IR.IntensityRefineError, match=match):
        IR.refine_one(record, donor_dir, source_stem=source_stem,
                      output_stem=output_stem, training_fold="even")


def test_training_rejects_unreliable_or_too_few_windows():
    gate = {"alignment": {"windows": [
        {"cx_um": 250.0, "cy_um": 250.0, "reason": "measured", "truncated": True,
         "contrast": 1.0, "contrast_shifted": 2.0}
    ]}}
    with pytest.raises(IR.IntensityRefineError, match="only 0 reliable"):
        IR._training_windows(gate, "even")
