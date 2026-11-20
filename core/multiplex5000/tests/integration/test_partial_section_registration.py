"""Registration must follow interior morphology and preserve reviewed section products."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from multiplex5000.integration import register
from multiplex5000.integration.rasterize import GridSpec
from multiplex5000.integration.transform import Transform, affine_from_pixels


def test_partial_field_uses_its_internal_structure_over_a_brighter_duplicate_outline():
    import cv2
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(182)
    h, w = 180, 280
    yy, xx = np.mgrid[:h, :w]
    mask = ((xx - 45) / 34) ** 2 + ((yy - 63) / 47) ** 2 < 1
    texture = gaussian_filter(rng.normal(size=(h, w)).astype("float32"), 1.65)
    texture /= texture.std()
    moving = ((0.60 + 0.085 * texture) * mask).astype("float32")
    truth = cv2.getRotationMatrix2D((45.0, 63.0), -5.0, 1)
    truth[:, 2] += [180, 36]
    present = cv2.warpAffine(moving, truth, (w, h))
    region = cv2.warpAffine(mask.astype("float32"), truth, (w, h))
    fixed = (0.31 + 0.10 * gaussian_filter(rng.normal(size=(h, w)).astype("float32"), 1.2)).astype(
        "float32"
    )
    fixed += 0.55 * present - 0.31 * region
    decoy = np.roll(
        (0.58 + 0.08 * gaussian_filter(rng.normal(size=(h, w)), 2)) * mask, (18, 12), axis=(0, 1)
    )
    fixed += decoy.astype("float32")
    grid = GridSpec(134, -79, 16.0, h, w)
    probes = grid.to_microns(np.array([[30, 75], [41, 90], [67, 49]]))
    result = register.register_image(
        moving, fixed, grid, moving_mask=mask, nonrigid=False, max_rotation_deg=30
    )
    expected = affine_from_pixels(truth, grid).apply(probes)
    assert np.max(np.linalg.norm(result.transform.apply(probes) - expected, axis=1)) < grid.pixel_um
    assert result.transform.meta["coarse_method"] == "regional_texture"


def test_stage_preserves_reviewed_registration_and_checks_its_binding(tmp_path, monkeypatch):
    cfg = SimpleNamespace(registration_dir=tmp_path, fixed_modality="phenocycler")
    row = {
        "donor_id": "example",
        "roi": "panc",
        "tissue": "pancreas",
        "pair_status": "paired",
        "pheno_image": str(tmp_path / "scan.tif"),
        "xenium_bundle": str(tmp_path / "bundle"),
    }
    monkeypatch.setattr(register, "load_manifest", lambda cfg: pd.DataFrame([row]))
    out = register.registration_dir(cfg, "example", "panc")
    transform = Transform(affine=np.array([[1.0, 0.0, 2.0], [0.0, 1.0, 4.0]]))
    result = register.RegistrationResult(transform, "image", stages=["landmark_review"])
    review = result.to_dict()
    review["section_landmark_validation"] = {
        "source": {
            "pc": row["pheno_image"],
            "xenium": str(Path(row["xenium_bundle"]) / "focus" / "image.tif"),
        },
        "landmarks": [
            {"xenium_native_um": [3.0, 8.0], "pc_native_um": [5.0, 12.0], "updated_error_um": 0.0}
        ],
    }

    def publish_review():
        transform.save(out)
        (out / "qc.json").write_text(json.dumps(review))

    # Simulate a review arriving between the initial read and persistence.
    def concurrent_review(*args, **kwargs):
        publish_review()
        return register.RegistrationResult(Transform(), "image")

    monkeypatch.setattr(register, "register_donor", concurrent_review)
    assert register.run_register(cfg).iloc[0]["stages"] == "landmark_review"
    saved = {p.name: p.read_bytes() for p in out.iterdir()}

    def forbidden(*args, **kwargs):
        pytest.fail("Existing reviewed section was recomputed")

    monkeypatch.setattr(register, "register_donor", forbidden)
    assert register.run_register(cfg).iloc[0]["stages"] == "landmark_review"
    row["xenium_bundle"] = str(tmp_path / "different-section")
    with pytest.raises(ValueError, match="different inputs or direction"):
        register.run_register(cfg)
    assert saved == {p.name: p.read_bytes() for p in out.iterdir()}
    row["xenium_bundle"] = str(tmp_path / "bundle")
    review["section_landmark_validation"]["landmarks"][0]["updated_error_um"] = 50.0
    (out / "qc.json").write_text(json.dumps(review))
    with pytest.raises(ValueError, match="recorded section landmarks"):
        register.run_register(cfg)
