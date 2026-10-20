"""UniFORM is done on a same-block basis (maintainer, 2026-10-02).

A same-block lymph node's cells belong to the pancreas block's sample; a different-block node on
the same slide is another sample.  These tests pin that the block is READ from the maintainer's
record and never inferred, and that all three places the lane picks its cells -- the usability
table, the shift, the cohort scale -- pick the same ones.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from shapely.geometry import box, mapping

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import alignment_block as AB
import panel
import run_step1 as S1
import run_step2 as S2
import run_step3 as S3

PX = AB.EXPORT_PIXEL_UM
IMAGE = "0080502_6580_pLNpanc_Scan1.er.qptiff - resolution #1"
SAME_ID, DIFF_ID = "node-in-the-pancreas-block", "node-from-another-block"


def _um(px: float) -> float:
    return px * PX


@pytest.fixture
def slide(tmp_path, monkeypatch):
    """One slide: pancreas tissue, a same-block node at x 1000-1100 px, a different-block node at 5000-5100."""
    flat = tmp_path / "expression"
    monkeypatch.setattr(panel, "EXPRESSION_FLAT", flat)
    monkeypatch.setattr(panel, "EXPRESSION", tmp_path / "retired")
    monkeypatch.setattr(panel, "donor_markers", lambda d=None: ["CD3e", "INS"])
    monkeypatch.setattr(panel, "cohort", lambda: pd.DataFrame({"donor_id": ["6580", "6374"]}))
    monkeypatch.setattr(AB, "OUT", tmp_path / "alignment_block")
    S2._sample_rows.cache_clear()

    rows = []

    def add(oid, region, x_px, y_px, cd3e, eligible=True):
        rows.append({"donor_id": "6580", "object_id": oid, "qc_analysis_eligible": eligible, "cell_region": region,
                     "X_centroid": _um(x_px), "Y_centroid": _um(y_px), "islet_num": None, "CD3e": cd3e, "INS": 50.0})
    for i in range(6):
        add(f"t{i}", "tissue", 100 + i, 100, 10.0 + i)
    for i in range(4):
        add(f"s{i}", "pLN", 1010 + i, 1050, 200.0 + i)                 # inside the same-block node
    for i in range(5):
        add(f"d{i}", "pLN", 5010 + i, 5050, 9000.0 + i)                # inside the different-block node
    add("edge", "pLN", 1104, 1050, 300.0)                              # 4 px = 2 um outside the same-block node
    add("lost", "pLN", 3000, 3000, 400.0)                              # in no node, near none
    add("ineligible", "pLN", 1020, 1020, 500.0, eligible=False)
    df = pd.DataFrame(rows)
    part = flat / "donor_id=6580"
    part.mkdir(parents=True)
    df.to_parquet(part / "data.parquet", index=False)

    other = pd.DataFrame({"donor_id": "6374", "object_id": [f"o{i}" for i in range(5)], "qc_analysis_eligible": True,
                          "cell_region": "tissue", "X_centroid": 1.0, "Y_centroid": 1.0, "islet_num": None,
                          "CD3e": [11.0, 12.0, 13.0, 14.0, 15.0], "INS": 40.0})
    (flat / "donor_id=6374").mkdir()
    other.to_parquet(flat / "donor_id=6374" / "data.parquet", index=False)

    annotations = tmp_path / "annotations"
    annotations.mkdir()
    feats = [{"type": "Feature", "id": "tissue-outline", "geometry": mapping(box(0, 0, 4000, 4000)),
              "properties": {"objectType": "annotation", "classification": {"name": "Tissue"}}},
             {"type": "Feature", "id": SAME_ID, "geometry": mapping(box(1000, 1000, 1100, 1100)),
              "properties": {"objectType": "annotation", "classification": {"name": "pLN"}}},
             {"type": "Feature", "id": DIFF_ID, "geometry": mapping(box(5000, 5000, 5100, 5100)),
              "properties": {"objectType": "annotation", "classification": {"name": "pLN"}}}]
    (annotations / f"{AB.image_tag(IMAGE)}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    (annotations / "6374.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats[:1]}))

    blocks = tmp_path / "pln_blocks.csv"
    blocks.write_text("image,object_id,block_from_notes\n"
                      f"{IMAGE},{SAME_ID},same-block\n{IMAGE},{DIFF_ID},different-block\n")
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"images": [{"donor_id": "6580", "image_id": IMAGE}, {"donor_id": "6374", "image_id": "6374"}]}))
    return {"spec": spec, "blocks": blocks, "annotations": annotations, "frame": df, "tmp": tmp_path}


def _build(slide):
    assert AB.main(["--spec", str(slide["spec"]), "--blocks", str(slide["blocks"]),
                    "--annotations", str(slide["annotations"])]) == 0


def test_a_cell_takes_the_block_of_the_node_that_holds_it(slide):
    _build(slide)
    side = pd.read_parquet(AB.OUT / "donor_id=6580.parquet").set_index("object_id")
    assert set(side.index) == {r for r in slide["frame"]["object_id"] if not r.startswith("t")}, "pLN-region cells only"
    assert (side.loc[[f"s{i}" for i in range(4)], "block"] == AB.SAME).all()
    assert (side.loc[[f"d{i}" for i in range(5)], "block"] == AB.DIFFERENT).all()
    assert side.loc["s0", "pln_annotation_id"] == SAME_ID and side.loc["d0", "pln_annotation_id"] == DIFF_ID
    # detected in the node, centroid 2 um outside its outline: the nearest node, and it says so
    assert side.loc["edge", "block"] == AB.SAME and "nearest" in side.loc["edge", "found_by"]
    # in no node and near none: nothing is assumed
    assert side.loc["lost", "block"] == AB.UNASSIGNED and side.loc["lost", "pln_annotation_id"] == ""
    assert not (AB.OUT / "donor_id=6374.parquet").exists(), "a donor with no pLN cells gets no sidecar"
    manifest = json.loads((AB.OUT / "manifest.json").read_text())
    row, = manifest["donors"]
    assert (row[AB.SAME], row[AB.DIFFERENT], row[AB.UNASSIGNED]) == (6, 5, 1)       # s0-3, edge, ineligible | d0-4 | lost
    assert manifest["donors_without_pln_cells"] == ["6374"]


def test_the_pancreas_block_is_tissue_plus_the_same_block_node(slide):
    _build(slide)
    f = slide["frame"]
    keep = AB.population("6580", f["object_id"], f["cell_region"])
    got = set(f.loc[keep, "object_id"])
    assert got == {f"t{i}" for i in range(6)} | {f"s{i}" for i in range(4)} | {"edge", "ineligible"}
    # a donor without pLN cells needs no sidecar: everything is its pancreas block
    assert AB.population("6374", ["a", "b"], ["tissue", "tissue"]).all()


def test_no_sidecar_is_a_refusal_not_a_default(slide):
    f = slide["frame"]
    with pytest.raises(SystemExit, match="same-block basis"):
        AB.population("6580", f["object_id"], f["cell_region"])
    _build(slide)
    with pytest.raises(SystemExit, match="older than the expression stage"):
        AB.population("6580", ["a-new-cell"], ["pLN"])


def test_a_node_without_a_recorded_block_stops_the_build(slide):
    slide["blocks"].write_text(f"image,object_id,block_from_notes\n{IMAGE},{SAME_ID},same-block\n")
    with pytest.raises(SystemExit, match="recorded, never inferred"):
        _build(slide)
    slide["blocks"].write_text("image,object_id,block_from_notes\n"
                               f"{IMAGE},{SAME_ID},same-block\n{IMAGE},{DIFF_ID},attached\n")
    with pytest.raises(SystemExit, match="must be 'same-block' or 'different-block'"):
        _build(slide)


def test_all_three_steps_use_the_block(slide):
    _build(slide)
    block_ids = {f"t{i}" for i in range(6)} | {f"s{i}" for i in range(4)} | {"edge"}      # eligible, in the block

    # step 1: the usability table's cells
    frame = S1.donor_frame("6580", ["CD3e", "INS"])
    assert set(frame["object_id"]) == block_ids
    assert frame.attrs["counts"] == {"rows": 18, "eligible": 17, "in_block": 11, "complete": 11}

    # step 2: the values a shift is estimated from -- no different-block cell (they read 9000+)
    vals = S2.marker_values("CD3e", ["6580", "6374"])
    assert len(vals["6580"]) == 11 and vals["6580"].max() < 1000 and len(vals["6374"]) == 5

    # step 2: the aligned table -- the different-block cells have no aligned value
    shifts = {m: {"scale": {"6580": 2.0, "6374": 1.0}} for m in ("CD3e", "INS")}
    out = slide["tmp"] / "uniform"
    S2.apply_shift(shifts, ["CD3e", "INS"], ["6580", "6374"], out=out)
    u = pd.read_parquet(out / "donor_id=6580" / "data.parquet").set_index("object_id")
    assert len(u) == 18, "one row per cell of the expression stage"
    assert u.loc["t0", "CD3e"] == 20.0 and u.loc["s0", "CD3e"] == 400.0 and u.loc["edge", "CD3e"] == 600.0
    outside = [f"d{i}" for i in range(5)] + ["lost"]
    assert u.loc[outside, ["CD3e", "INS"]].isna().all().all() and not u.loc[outside, "in_alignment_block"].any()
    assert u.loc[sorted(block_ids), "in_alignment_block"].all() and u["CD3e"].notna().sum() == 12

    # step 3: the cohort scale's sample
    S3_ALIGNED = S3.ALIGNED
    try:
        S3.ALIGNED = out
        seen = {}

        def marker_scale(values, floors):
            seen.update(values)
            return pd.DataFrame()
        real, S3.SC.marker_scale = S3.SC.marker_scale, marker_scale
        S3.build_scale(["6580", "6374"], ["CD3e"], ratio=False)
    finally:
        S3.ALIGNED, S3.SC.marker_scale = S3_ALIGNED, real
    assert len(seen["CD3e"]) == 11 + 5 and np.isfinite(seen["CD3e"]).all()
