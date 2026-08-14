from __future__ import annotations

import json

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts.review.qupath_fundamental_inputs import (
    SOURCE_MEASUREMENTS,
    load_qupath_projection_donor,
    load_qupath_sources,
    refresh_qupath_projection,
    validate_qupath_projection,
)


def _export_frame() -> pd.DataFrame:
    rows = []
    for image, prefix in (("image-a", "a"), ("image-b", "b"), ("not-selected", "x")):
        for index in range(3):
            rows.append(
                {
                    "Image": image,
                    "Object ID": f"{prefix}-{index}",
                    "Centroid X µm": 10 + index,
                    "Centroid Y µm": 20 + index,
                    "Cell: Area µm^2": 30 + index,
                    "Cell: Solidity": 0.8,
                    "Cell: DAPI: Mean": 100 + index,
                    "Cell: CD3e: Mean": 200 + index,
                    "Cell: CD79a: Mean": 300 + index,
                    "Cell: Iba1: Mean (norm)": [-0.5, 0.25, 1.5][index],
                    "Cell: Ker8-18: Mean": 400 + index,
                    # These tempting alternatives must never enter the projection.
                    "Cell: Iba1: Mean": 9999,
                    "Iba1__restore_normalized": 8888,
                }
            )
    return pd.DataFrame.from_records(rows)


def _write_sources(tmp_path):
    source = tmp_path / "cells.csv"
    _export_frame().to_csv(source, index=False)
    spec = tmp_path / "cohort-spec.json"
    spec.write_text(
        json.dumps(
            {
                "defaults": {"measurement_csv": [str(source)]},
                "images": [
                    {"donor_id": "D1", "image_id": "image-a", "measurement_csv_index": 0},
                    {"donor_id": "D2", "image_id": "image-b", "measurement_csv_index": 0},
                ],
            }
        ),
        encoding="utf-8",
    )
    return source, spec


def test_live_spec_resolves_exact_selected_source(tmp_path):
    source, spec = _write_sources(tmp_path)

    observed = load_qupath_sources(spec, ("D2", "D1"))

    assert [item.donor_id for item in observed] == ["D2", "D1"]
    assert {item.measurement_csv for item in observed} == {source}
    assert [item.image_id for item in observed] == ["image-b", "image-a"]


def test_live_spec_rejects_one_image_assigned_to_two_donors(tmp_path):
    _source, spec = _write_sources(tmp_path)
    payload = json.loads(spec.read_text(encoding="utf-8"))
    payload["images"][1]["image_id"] = payload["images"][0]["image_id"]
    spec.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="share image_id"):
        load_qupath_sources(spec, ("D1", "D2"))


def test_projection_reads_exact_qupath_norm_and_filters_exact_images(tmp_path):
    _source, spec = _write_sources(tmp_path)
    output = tmp_path / "development" / "qupath-inputs"

    refresh_qupath_projection(spec, output, ("D1", "D2"), block_size=1024)

    first = pd.read_parquet(output / "donor_id=D1" / "data.parquet")
    schema = pq.ParquetFile(output / "donor_id=D1" / "data.parquet").schema_arrow
    assert first["object_id"].tolist() == ["a-0", "a-1", "a-2"]
    assert first["Iba1_norm"].tolist() == [-0.5, 0.25, 1.5]
    assert tuple(first.columns) == ("donor_id", *SOURCE_MEASUREMENTS)
    assert all(
        schema.field(column).type == pa.float64()
        for column in SOURCE_MEASUREMENTS
        if column not in {"image", "object_id"}
    )
    assert "Iba1" not in first
    assert "Iba1__restore_normalized" not in first
    assert "CD79a" not in first
    summary = json.loads((output / "projection_summary.json").read_text(encoding="utf-8"))
    assert summary["restore_outputs_read"] is False
    assert summary["source_measurements"]["Iba1_norm"] == "Cell: Iba1: Mean (norm)"
    assert summary["donors"] == {"D1": 3, "D2": 3}
    assert summary["donor_images"] == {"D1": "image-a", "D2": "image-b"}


@pytest.mark.parametrize("fault", ["count", "donor", "image"])
def test_projection_consumer_checks_declared_partition_identity(tmp_path, fault):
    _source, spec = _write_sources(tmp_path)
    output = tmp_path / "development" / "qupath-inputs"
    refresh_qupath_projection(spec, output, ("D1",), block_size=1024)
    path = output / "donor_id=D1" / "data.parquet"
    frame = pd.read_parquet(path)
    if fault == "count":
        frame = frame.iloc[:-1]
    elif fault == "donor":
        frame.loc[0, "donor_id"] = "WRONG"
    else:
        frame.loc[0, "image"] = "wrong-image"
    frame.to_parquet(path, index=False)

    with pytest.raises(ValueError, match="row count|internal donor_id|declared image"):
        load_qupath_projection_donor(output, "D1", ("Iba1_norm", "Ker8_18"))


def test_projection_contract_rejects_duplicate_donor_image_inventory(tmp_path):
    _source, spec = _write_sources(tmp_path)
    output = tmp_path / "development" / "qupath-inputs"
    refresh_qupath_projection(spec, output, ("D1", "D2"), block_size=1024)
    summary_path = output / "projection_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["donor_images"]["D2"] = summary["donor_images"]["D1"]
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(ValueError, match="one image to multiple donors"):
        validate_qupath_projection(output, ("D1", "D2"))


def test_projection_refresh_is_mutable_and_transactional(tmp_path):
    source, spec = _write_sources(tmp_path)
    output = tmp_path / "development" / "qupath-inputs"
    sibling = output.parent / "keep"
    sibling.mkdir(parents=True)
    (sibling / "user.txt").write_text("keep", encoding="utf-8")
    refresh_qupath_projection(spec, output, ("D1",), block_size=1024)

    changed = _export_frame()
    changed.loc[changed["Image"].eq("image-a"), "Cell: Iba1: Mean (norm)"] = 7.0
    changed.to_csv(source, index=False)
    refresh_qupath_projection(spec, output, ("D1",), block_size=1024)

    observed = pd.read_parquet(output / "donor_id=D1" / "data.parquet")
    assert observed["Iba1_norm"].tolist() == [7.0, 7.0, 7.0]
    assert (sibling / "user.txt").read_text(encoding="utf-8") == "keep"

    broken = changed.drop(columns="Cell: Iba1: Mean (norm)")
    broken.to_csv(source, index=False)
    with pytest.raises(ValueError, match=r"Cell: Iba1: Mean \(norm\)"):
        refresh_qupath_projection(spec, output, ("D1",), block_size=1024)
    pd.testing.assert_frame_equal(
        observed,
        pd.read_parquet(output / "donor_id=D1" / "data.parquet"),
    )


def test_projection_refuses_unowned_output(tmp_path):
    _source, spec = _write_sources(tmp_path)
    output = tmp_path / "development" / "user-data"
    output.mkdir(parents=True)
    (output / "keep.txt").write_text("mine", encoding="utf-8")

    with pytest.raises(ValueError, match="not owned"):
        refresh_qupath_projection(spec, output, ("D1",), block_size=1024)

    assert (output / "keep.txt").read_text(encoding="utf-8") == "mine"
