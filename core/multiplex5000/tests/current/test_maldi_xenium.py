"""Exact-region morphology, declared planes and physical-frame contracts."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import tifffile

from multiplex5000.stages import maldi_serial_registration as R
from multiplex5000.tasks import Input, Task, sha256_file


def _focus(path, *, channel="DAPI", pixel=10.0, unit="µm"):
    image = np.zeros((128, 128), dtype=np.uint16)
    image[25:105, 20:80] = 1000
    image[25:60, 80:115] = 1000
    with tifffile.TiffWriter(path, ome=True) as writer:
        writer.write(image, subifds=2, metadata={"axes": "YX", "Channel": {"Name": channel},
                     "PhysicalSizeX": pixel, "PhysicalSizeY": pixel,
                     "PhysicalSizeXUnit": unit, "PhysicalSizeYUnit": unit})
        writer.write(image[::2, ::2], subfiletype=1)
        writer.write(image[::4, ::4], subfiletype=1)
    return image


def _source(tmp_path):
    bundle = tmp_path / "output-XETG00001__0000001__Region_1__20260101__000000"
    bundle.mkdir()
    experiment = bundle / "experiment.xenium"
    experiment.write_text(json.dumps({"slide_id": "0000001", "region_name": "Region_1",
                          "pixel_size": 10.0, "images": {"morphology_focus_filepath": "dapi.ome.tif"}}))
    image = bundle / "dapi.ome.tif"
    _focus(image)
    release = tmp_path / "release"
    (release / "00_inventory").mkdir(parents=True)
    receipt = release / "00_inventory/receipt.json"
    receipt.write_text(json.dumps({"inputs": [{"name": f"xenium:{bundle.name}:experiment.xenium",
                                             "sha256": sha256_file(experiment)}]}))
    project = SimpleNamespace(release_dir=release, source=lambda _: None)
    neighbour = {"acquisition_id": "xenium:0000001:Region_1", "output": bundle.name,
                 "path": str(bundle), "platform": "xenium", "tissue": "pancreas"}
    return project, neighbour, experiment, image


def test_reader_uses_declared_dapi_and_known_pyramid_scale(tmp_path):
    path = tmp_path / "dapi.ome.tif"
    raw = _focus(path)
    image, pixel, metadata = R._read_xenium_dapi(path, 10.0)
    np.testing.assert_array_equal(image, raw[::4, ::4])
    assert pixel == 40.0 and metadata["pyramid_level"] == 2
    assert metadata["base_shape"] == [128, 128]


@pytest.mark.parametrize("channel,pixel,unit", [("boundary", 10.0, "µm"),
                                              ("DAPI", 10.0, "nm"), ("DAPI", 11.0, "µm")])
def test_wrong_channel_or_physical_frame_rejected(tmp_path, channel, pixel, unit):
    path = tmp_path / "dapi.ome.tif"
    _focus(path, channel=channel, pixel=pixel, unit=unit)
    with pytest.raises(R.MaldiRegistrationError):
        R._read_xenium_dapi(path, 10.0)


def test_region_and_receipt_are_checked_before_using_image(tmp_path):
    project, neighbour, experiment, _ = _source(tmp_path)
    assert R._xenium_source(project, neighbour)["status"] == "ready"
    with pytest.raises(R.MaldiRegistrationError, match="another slide or region"):
        R._xenium_source(project, dict(neighbour, acquisition_id="xenium:0000001:Region_2"))
    data = json.loads(experiment.read_text())
    data["pixel_size"] = 11.0
    experiment.write_text(json.dumps(data))
    with pytest.raises(R.MaldiRegistrationError, match="differs from inventoried input"):
        R._xenium_source(project, neighbour)


@pytest.mark.parametrize("conflicting", [False, True])
def test_experiment_receipt_binding_must_be_unique(tmp_path, conflicting):
    project, neighbour, _, _ = _source(tmp_path)
    path = Path(project.release_dir) / "00_inventory/receipt.json"
    receipt = json.loads(path.read_text())
    second = dict(receipt["inputs"][0])
    if conflicting:
        second["sha256"] = "0" * 64
    receipt["inputs"].append(second)
    path.write_text(json.dumps(receipt))
    with pytest.raises(R.MaldiRegistrationError, match="unique inventoried binding"):
        R._xenium_source(project, neighbour)


def test_unavailable_and_damaged_exact_region_do_not_substitute(tmp_path):
    project, neighbour, _, image = _source(tmp_path)
    assert "damaged" in R._xenium_source(project, dict(neighbour, condition="damaged"))["reason"]
    image.unlink()
    result = R._xenium_source(project, neighbour)
    assert result["status"] == "unavailable" and result["image"] == image


def test_explicit_relocation_uses_bound_identity_not_directory_spelling(tmp_path):
    project, neighbour, experiment, _ = _source(tmp_path)
    original = experiment.parent
    relocated = original.with_name("renamed-bundle")
    original.rename(relocated)
    index = tmp_path / "source_index.csv"
    index.write_text("acquisition_id,output,source_path\n" +
                     f"{neighbour['acquisition_id']},{neighbour['output']},{relocated}\n")
    project.source = lambda name: index if name == "xenium_source_index" else None
    source = R._xenium_source(project, neighbour)
    assert source["status"] == "ready"
    assert source["source_binding"]["source_index_sha256"] == sha256_file(index)
    changed = relocated / "experiment.xenium"
    changed.write_text(changed.read_text().replace('"Region_1"', '"Region_2"'))
    with pytest.raises(R.MaldiRegistrationError, match="another slide or region"):
        R._xenium_source(project, neighbour)


def test_undeclared_rename_and_wrong_index_output_rejected(tmp_path):
    project, neighbour, experiment, _ = _source(tmp_path)
    renamed = experiment.parent.with_name("renamed-bundle")
    experiment.parent.rename(renamed)
    with pytest.raises(R.MaldiRegistrationError, match="mismatched exact"):
        R._xenium_source(project, dict(neighbour, path=str(renamed)))
    index = tmp_path / "wrong_index.csv"
    index.write_text("acquisition_id,output,source_path\n" +
                     f"{neighbour['acquisition_id']},wrong-output,{renamed}\n")
    project.source = lambda _: index
    with pytest.raises(R.MaldiRegistrationError, match="one exact output"):
        R._xenium_source(project, neighbour)


def test_declared_focus_cannot_escape_region_bundle(tmp_path):
    project, neighbour, experiment, _ = _source(tmp_path)
    data = json.loads(experiment.read_text())
    data["images"]["morphology_focus_filepath"] = "../other.ome.tif"
    experiment.write_text(json.dumps(data))
    receipt = Path(project.release_dir) / "00_inventory/receipt.json"
    bound = json.loads(receipt.read_text())
    bound["inputs"][0]["sha256"] = sha256_file(experiment)
    receipt.write_text(json.dumps(bound))
    with pytest.raises(R.MaldiRegistrationError, match="escapes"):
        R._xenium_source(project, neighbour)


def test_fixed_grid_uses_native_region_origin_and_records_source(tmp_path):
    project, neighbour, experiment, image = _source(tmp_path)
    row = {"tissue_id": "sample", "neighbour": neighbour}
    declared = R._neighbour_inputs(project, row)
    assert next(i for i in declared if i.name == "xenium_dapi:sample").hash
    task = Task("07_registration", "candidate", tmp_path / "out",
                {i.name: i for i in declared}, project=project)
    fixed = R._fixed_xenium(task, row)
    assert fixed["origin"] == fixed["dapi_origin"] == (25.0, 25.0)
    assert fixed["mask_grid"].pixel_um == 50.0 and fixed["mask"].any()
    assert fixed["frame_metadata"]["region"] == neighbour["acquisition_id"]
    assert fixed["frame_metadata"]["experiment_sha256"] == sha256_file(experiment)
    other = tmp_path / "other.ome.tif"
    other.write_bytes(image.read_bytes())
    task.inputs["xenium_dapi:sample"] = Input("xenium_dapi:sample", other)
    with pytest.raises(R.MaldiRegistrationError, match="changed after input declaration"):
        R._fixed_xenium(task, row)


def test_empty_moving_mask_is_explicitly_unavailable(tmp_path):
    from multiplex5000.integration.rasterize import GridSpec

    path = tmp_path / "empty.npz"
    np.savez(path, summed=np.zeros((3, 3)), on_tissue=np.zeros((3, 3), bool),
             grid=json.dumps(GridSpec(0, 0, 50, 3, 3).as_dict()))
    task = Task("07_registration", "candidate", tmp_path / "out",
                {"maldi_image:sample": Input("maldi_image:sample", path)})
    with pytest.raises(R.MaldiRegistrationError, match="contains no spots"):
        R._moving(task, "sample")
