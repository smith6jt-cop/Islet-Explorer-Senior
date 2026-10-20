"""The call schema every arm writes, so the gate reads one shape.

``validate_calls`` must refuse what the contract forbids -- a missing reason on an
Ambiguous/Unavailable row, a duplicated ``object_id``, a label no declaration carries, and a
call set that is not the frozen cell set -- and ``fill_unavailable`` must give EVERY frozen
cell a row.  A cell an arm forgets is ``Unavailable(arm_no_output)``, never absent.
"""
from __future__ import annotations

import json
import uuid

import numpy as np
import pandas as pd
import pytest

import pins
import records
from records import AMBIGUOUS, CALL_COLUMNS, UNAVAILABLE, fill_unavailable, validate_calls


def _ids(n, seed=0):
    rng = np.random.default_rng(seed)
    return np.array([str(uuid.UUID(bytes=rng.bytes(16), version=4)) for _ in range(n)], dtype=object)


def _calls(ids, labels, reasons=None, scores=None):
    return pd.DataFrame({
        "object_id": ids,
        "label": labels,
        "score": scores if scores is not None else np.linspace(0, 1, len(ids)),
        "reason": reasons if reasons is not None else [""] * len(ids),
    })


@pytest.fixture
def labels(dec):
    return set(dec["labels"])


# --------------------------------------------------------------------------- validate_calls
def test_a_well_formed_frame_passes(dec, labels):
    ids = _ids(6)
    df = _calls(ids, ["Endocrine", "T_NK", AMBIGUOUS, UNAVAILABLE, "Myeloid", "Stromal"],
                reasons=["", "", "multi_type", "ablation_removes_required_input", "", ""])
    out = validate_calls(df, expected_ids=ids, labels=labels)
    assert list(out.columns) == list(CALL_COLUMNS)
    assert len(out) == 6


def test_missing_columns_are_rejected():
    with pytest.raises(ValueError, match="missing columns"):
        validate_calls(pd.DataFrame({"object_id": _ids(2), "label": ["T_NK", "T_NK"]}))


@pytest.mark.parametrize("bucket", [AMBIGUOUS, UNAVAILABLE])
def test_an_ambiguous_or_unavailable_call_without_a_reason_is_rejected(bucket):
    ids = _ids(3)
    df = _calls(ids, ["T_NK", bucket, "Myeloid"], reasons=["", "", ""])
    with pytest.raises(ValueError, match="carry no reason"):
        validate_calls(df)
    df["reason"] = ["", None, ""]           # a NaN reason is an empty reason
    with pytest.raises(ValueError, match="carry no reason"):
        validate_calls(df)


def test_a_called_cell_needs_no_reason_and_a_bucket_with_one_passes():
    ids = _ids(2)
    validate_calls(_calls(ids, ["T_NK", AMBIGUOUS], reasons=["", "low_margin"]))


def test_duplicate_object_ids_are_rejected():
    ids = _ids(3)
    ids[2] = ids[0]
    with pytest.raises(ValueError, match="duplicate object_id"):
        validate_calls(_calls(ids, ["T_NK"] * 3))


def test_an_undeclared_label_is_rejected(labels):
    ids = _ids(3)
    with pytest.raises(ValueError, match=r"undeclared labels in calls: \['Neural'\]"):
        validate_calls(_calls(ids, ["T_NK", "Neural", "Myeloid"]), labels=labels)
    with pytest.raises(ValueError, match="undeclared labels"):
        validate_calls(_calls(ids, ["T_NK", "Immune", "Myeloid"]), labels=labels)


def test_the_buckets_are_never_undeclared(labels):
    """Ambiguous and Unavailable are part of the schema, not labels an arm must declare."""
    ids = _ids(2)
    validate_calls(_calls(ids, [AMBIGUOUS, UNAVAILABLE], reasons=["a", "b"]), labels=labels)


def test_a_cell_set_mismatch_is_rejected_missing_and_extra():
    ids = _ids(5)
    other = _ids(2, seed=1)
    df = _calls(np.concatenate([ids[:3], other]), ["T_NK"] * 5)
    with pytest.raises(ValueError, match="2 missing, 2 extra"):
        validate_calls(df, expected_ids=ids)
    with pytest.raises(ValueError, match="1 missing, 0 extra"):
        validate_calls(_calls(ids[:4], ["T_NK"] * 4), expected_ids=ids)
    with pytest.raises(ValueError, match="0 missing, 1 extra"):
        validate_calls(_calls(np.concatenate([ids, other[:1]]), ["T_NK"] * 6), expected_ids=ids)


@pytest.mark.xfail(strict=True, reason=(
    "records.validate_calls compares the call set with pandas Index.equals, which is order-"
    "sensitive: a COMPLETE set in another order is rejected as '0 missing, 0 extra'. Harmless "
    "on the driver path (arm_base.write reorders through fill_unavailable first) but the message "
    "misdescribes the failure; recorded here, not corrected in the spine."))
def test_a_complete_set_in_another_order_is_accepted():
    ids = _ids(5)
    validate_calls(_calls(ids[::-1], ["T_NK"] * 5), expected_ids=ids)


# --------------------------------------------------------------------------- fill_unavailable
def test_fill_unavailable_covers_the_frame_in_frame_order():
    ids = _ids(8)
    sub = ids[[6, 2, 4]]
    df = _calls(sub, ["T_NK", AMBIGUOUS, "Stromal"], reasons=["", "multi_type", ""], scores=[0.9, 0.5, 0.7])
    out = fill_unavailable(df, ids)
    assert list(out.columns) == list(CALL_COLUMNS)
    assert out["object_id"].tolist() == list(ids), "every frozen cell, in the frozen order"
    by = out.set_index("object_id")
    assert by.loc[ids[6], "label"] == "T_NK" and by.loc[ids[6], "reason"] == ""
    assert by.loc[ids[2], "label"] == AMBIGUOUS and by.loc[ids[2], "reason"] == "multi_type"
    assert by.loc[ids[4], "score"] == pytest.approx(0.7)
    missing = by.drop(index=list(sub))
    assert (missing["label"] == UNAVAILABLE).all()
    assert (missing["reason"] == "arm_no_output").all()
    assert missing["score"].isna().all()
    assert out["score"].dtype == float
    # and the result passes the validator against the same frame
    validate_calls(out, expected_ids=ids)


def test_fill_unavailable_with_no_calls_at_all_is_all_unavailable():
    ids = _ids(4)
    empty = pd.DataFrame({c: pd.Series(dtype=t) for c, t in
                          zip(CALL_COLUMNS, (object, object, float, object))})
    out = fill_unavailable(empty, ids, reason="ablation_removes_required_input")
    assert len(out) == 4
    assert (out["label"] == UNAVAILABLE).all()
    assert (out["reason"] == "ablation_removes_required_input").all()
    validate_calls(out, expected_ids=ids)


def test_fill_unavailable_never_invents_a_cell_outside_the_frame():
    """A call for a cell the frame does not hold is dropped, not appended: the frame is the set."""
    ids = _ids(4)
    stranger = _ids(1, seed=7)
    out = fill_unavailable(_calls(np.concatenate([ids[:2], stranger]), ["T_NK"] * 3), ids)
    assert out["object_id"].tolist() == list(ids)
    assert stranger[0] not in set(out["object_id"])


def test_fill_unavailable_coerces_object_id_to_str():
    ids = np.arange(5).astype(object)
    df = _calls(np.array([3, 1]), ["T_NK", "Myeloid"])
    out = fill_unavailable(df, ids)
    assert out["object_id"].tolist() == ["0", "1", "2", "3", "4"]
    assert out.set_index("object_id").loc["3", "label"] == "T_NK"


# --------------------------------------------------------------------------- write_calls
def test_write_calls_writes_the_parquet_and_appends_the_record(tmp_path, monkeypatch, dec, labels):
    """The JSONL beside the calls carries ``identity_markers`` -- what the tiering test reads."""
    monkeypatch.setattr(pins, "OUT", tmp_path / "broadtype")
    ids = _ids(5)
    calls = fill_unavailable(_calls(ids[:3], ["T_NK", "Myeloid", "Endocrine"]), ids)
    rec = records.ArmRecord(arm="toy", donor_id="0000", config_id="ablate:CD3e",
                            identity_markers=[m for m in dec["identity_markers"] if m != "CD3e"],
                            scope="frame", n_cells=5, n_called=3, n_ambiguous=0, n_unavailable=2,
                            seconds=0.1, env="test", params={"k": 1})
    path = records.write_calls(calls, rec, labels=labels, expected_ids=ids)
    assert path == records.out_dir("toy", "0000", "frame") / "calls__ablate_CD3e.parquet"
    assert path.exists()
    back = pd.read_parquet(path)
    assert back["object_id"].tolist() == list(ids)
    lines = (path.parent / "records.jsonl").read_text().splitlines()
    assert len(lines) == 1
    row = json.loads(lines[0])
    assert row["arm"] == "toy" and row["config_id"] == "ablate:CD3e" and row["calls"] == str(path)
    assert "CD3e" not in row["identity_markers"]
    assert set(row["identity_markers"]) <= set(dec["identity_markers"])
    assert not set(row["identity_markers"]) & set(dec["state_markers"])
    assert "commit" in row["code"] and "written_utc" in row
    # a second write appends, never truncates
    records.write_calls(calls, rec, labels=labels, expected_ids=ids)
    assert len((path.parent / "records.jsonl").read_text().splitlines()) == 2


def test_write_calls_refuses_a_bad_frame_before_touching_the_record(tmp_path, monkeypatch, labels):
    monkeypatch.setattr(pins, "OUT", tmp_path / "broadtype")
    ids = _ids(3)
    rec = records.ArmRecord(arm="toy", donor_id="0000", config_id="primary", identity_markers=[],
                            scope="frame", n_cells=3, n_called=3, n_ambiguous=0, n_unavailable=0,
                            seconds=0.0, env="test")
    with pytest.raises(ValueError, match="undeclared labels"):
        records.write_calls(_calls(ids, ["T_NK", "Neural", "T_NK"]), rec, labels=labels, expected_ids=ids)
    d = records.out_dir("toy", "0000", "frame")
    assert not (d / "records.jsonl").exists()
