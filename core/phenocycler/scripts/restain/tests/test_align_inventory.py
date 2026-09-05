"""Annotation-free registration-gate contracts."""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import align as AL  # noqa: E402
import registration_inputs as RI  # noqa: E402
import run_align as RA  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402


def _record(tmp_path: Path, *, geojson=None) -> RI.RegistrationInput:
    return RI.RegistrationInput(
        donor_id="6539", serial="0059865", roi="panc",
        qptiff=tmp_path / "0059865_Scan1.qptiff", bundle=tmp_path / "bundle",
        geojson=geojson, pixel_um=0.5, dapi_page=0, channel_names=("DAPI",),
        pyramid_shapes=((32, 32),), provenance_source="if_transfer_provenance.csv")


def _window_signature(windows):
    return [(w["cx_um"], w["cy_um"], w["bounds"], len(w["member_indices"]))
            for w in windows]


def test_occupied_windows_are_origin_anchored_deterministic_and_partition_points():
    xy = np.array([
        [-0.1, -0.1], [0.0, 0.0], [499.999, 499.999],
        [500.0, 0.0], [100.0, 500.0], [np.nan, 20.0],
    ])
    windows = RA.occupied_windows(xy)
    assert _window_signature(windows) == [
        (-250.0, -250.0, (-500.0, -500.0, 0.0, 0.0), 1),
        (250.0, 250.0, (0.0, 0.0, 500.0, 500.0), 2),
        (750.0, 250.0, (500.0, 0.0, 1000.0, 500.0), 1),
        (250.0, 750.0, (0.0, 500.0, 500.0, 1000.0), 1),
    ]
    assigned = np.concatenate([w["member_indices"] for w in windows])
    assert sorted(assigned.tolist()) == [0, 1, 2, 3, 4]

    shuffled = xy[[3, 0, 4, 2, 5, 1]]
    assert _window_signature(RA.occupied_windows(shuffled)) == _window_signature(windows)


def test_occupied_windows_refuse_a_non_partitioning_grid():
    with pytest.raises(AL.AlignError, match="size_um == step_um"):
        RA.occupied_windows(np.zeros((2, 2)), size_um=500.0, step_um=250.0)
    with pytest.raises(AL.AlignError, match=r"shape \(n, 2\)"):
        RA.occupied_windows(np.zeros((2, 3)))


def test_annotation_window_request_without_geojson_is_refused_before_io(tmp_path):
    with pytest.raises(AL.AlignError, match="no GeoJSON"):
        RA.align_one(_record(tmp_path), windows="tissue", transform="refined")


def test_tissue_window_request_with_islet_only_geojson_is_refused_before_io(
        tmp_path, monkeypatch):
    path = tmp_path / "islets.geojson"
    monkeypatch.setattr(RA.geo, "load_annotations",
                        lambda _path: {RA.geo.ISLET: [object()]})
    with pytest.raises(AL.AlignError, match="no 'Tissue' annotation"):
        RA.align_one(_record(tmp_path, geojson=path), windows="tissue", transform="refined")


def test_annotation_derived_ecc_refuses_islet_only_geojson_before_registration(
        tmp_path, monkeypatch):
    path = tmp_path / "islets.geojson"
    monkeypatch.setattr(RA.geo, "load_annotations",
                        lambda _path: {RA.geo.ISLET: [object()]})
    monkeypatch.setattr(
        RA, "register_one",
        lambda *_args, **_kwargs: pytest.fail("ECC registration started without Tissue"))
    with pytest.raises(AL.AlignError, match="ECC transform without a 'Tissue' annotation"):
        RA.align_one(_record(tmp_path, geojson=path), windows="auto", transform="ecc")


@pytest.mark.parametrize("islet_only", [False, True], ids=["no-geojson", "islet-only"])
def test_refined_alignment_uses_inventory_and_never_calls_geojson_registration(
        tmp_path, monkeypatch, islet_only):
    import shapely
    import tifffile

    geojson = tmp_path / "islets.geojson" if islet_only else None
    record = _record(tmp_path, geojson=geojson)
    if islet_only:
        monkeypatch.setattr(RA.geo, "load_annotations",
                            lambda _path: {RA.geo.ISLET: [SimpleNamespace(area=0.0)]})
    # Eighty nuclei in one occupied 500-um bin, enough for the split-half precision arm.
    geoms = np.asarray([shapely.box(10 + i, 10, 11 + i, 11) for i in range(80)], dtype=object)
    donor_dir = tmp_path / "valis" / "donor_id=6539"
    Transform(affine=np.eye(2, 3)).save(donor_dir, stem="refined")
    monkeypatch.setattr(RA, "register_one",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError(
                            "refined mode must not call the GeoJSON/ECC registration")))
    monkeypatch.setattr(RA, "_median_cell_diameter_um", lambda bundle: 10.0)
    monkeypatch.setattr(RA.B, "load_xenium_nucleus_polygons",
                        lambda bundle: (SimpleNamespace(values=geoms), np.arange(len(geoms))))
    monkeypatch.setattr(RA.TL, "window_grid",
                        lambda bounds, pixel_um: (GridSpec(0.0, 0.0, pixel_um, 8, 8),
                                                  (0, 8, 0, 8)))
    monkeypatch.setattr(RA.TL, "read_window_on_grid",
                        lambda page, grid, pixel_um: (np.ones(grid.shape),
                                                      np.ones(grid.shape, dtype=bool)))
    monkeypatch.setattr(RA, "rasterize_binary",
                        lambda geoms, grid: np.ones(grid.shape, dtype=bool))

    def measured(_dapi, _valid, _template, *, pixel_um, islet_index=-1, n_nuclei=0,
                 **_kwargs):
        return AL.Offset(islet_index=islet_index, n_nuclei=n_nuclei, dx_um=0.0, dy_um=0.0,
                         magnitude_um=0.0, ncc_peak=1.0, peak_z=10.0, reason="measured")

    monkeypatch.setattr(RA.AL, "offset_by_ncc", measured)
    monkeypatch.setattr(RA.AL, "contrast_inside_vs_gaps", lambda *a, **k: 2.0)

    page = object()
    level = SimpleNamespace(shape=(32, 32), pages=[page])
    fake_tif = SimpleNamespace(series=[SimpleNamespace(levels=[level])])

    class FakeTiffFile:
        def __init__(self, _path):
            pass

        def __enter__(self):
            return fake_tif

        def __exit__(self, *_exc):
            return False

    monkeypatch.setattr(tifffile, "TiffFile", FakeTiffFile)

    rec = RA.align_one(record, transform="refined", windows="auto", out_dir=tmp_path)
    alignment = rec["alignment"]
    assert rec["geojson"] == (str(geojson) if islet_only else None)
    assert alignment["window_set"] == "occupied"
    assert alignment["verdict"] == "PASS"
    assert alignment["certifies_mapping"] == {
        "pixel_um_full_res": 0.5, "source": "qptiff XML PixelSizeMicrons",
        "dapi_page": 0, "dapi_channel_name": "DAPI",
    }
    assert alignment["transform_file"].endswith("donor_id=6539/refined.json")
    assert alignment["transform_sha256"] == RA._sha256(donor_dir / "refined.json")
    assert alignment["displacement_sha256"] is None


def test_checkerboard_window_folds_are_disjoint_and_exhaustive():
    xy = np.array([[10.0, 10.0], [510.0, 10.0], [10.0, 510.0], [510.0, 510.0]])
    probes = RA.occupied_windows(xy)
    even = RA._select_window_fold(probes, "even")
    odd = RA._select_window_fold(probes, "odd")
    assert {w["index"] for w in even}.isdisjoint({w["index"] for w in odd})
    assert {w["index"] for w in even + odd} == {w["index"] for w in probes}
    assert len(even) == len(odd) == 2


def test_main_writes_atomic_per_donor_refined_gate_and_cohort_summary(tmp_path, monkeypatch):
    record = _record(tmp_path)
    monkeypatch.setattr(RA.RI, "load_registration_inputs", lambda *a, **k: {"6539": record})
    selected = {}

    def select(records, donors, *, context):
        selected.update(records=records, donors=donors, context=context)
        return [record]

    monkeypatch.setattr(RA.RI, "select_registration_inputs", select)
    monkeypatch.setattr(RA.RI, "require_registration_files", lambda _record: None)
    result = {
        "donor_id": "6539", "serial": "0059865", "roi": "panc",
        "alignment": {
            "verdict": "PASS", "verdict_why": "synthetic pass",
            "median_offset_um": 0.0, "threshold_pass_um": 2.0,
            "threshold_fail_um": 4.0, "transform_source": "refined",
            "n_measured": 1, "n_windows": 1, "window_set": "occupied",
            "n_truncated": 0,
        },
    }
    monkeypatch.setattr(RA, "align_one", lambda record, **kwargs: result)

    code = RA.main(["--out", str(tmp_path), "--donors", "6539"])
    assert code == 0
    assert selected["donors"] == ["6539"] and selected["context"] == "run_align"
    gate = tmp_path / "valis" / "donor_id=6539" / "alignment_refined.json"
    cohort = tmp_path / "alignment_panc_refined.json"
    assert gate.exists() and cohort.exists()
    import json
    gate_rec = json.loads(gate.read_text())
    cohort_rec = json.loads(cohort.read_text())
    assert gate_rec["alignment"]["gate_file"] == str(gate)
    assert cohort_rec == [gate_rec]
    assert not list(gate.parent.glob(".*.tmp"))
