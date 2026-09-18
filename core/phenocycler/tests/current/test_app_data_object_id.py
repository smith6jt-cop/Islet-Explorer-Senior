"""The app's per-cell artifacts must carry the QuPath UUID they were built from.

``cells/`` and ``donors/`` had no ``object_id`` until 2026-09-18.  A drill-down cell therefore
traced back to nothing: the app could draw it, colour it and report it, and no reader could say
which QuPath detection it was.  The repo rule is that the UUID is preserved throughout and is
never stood in for by a centroid, a row number or a filename, so these tests hold the column in
place -- present, first, non-null and unique within the donor -- and hold the refusal in place
for the two ways it can stop being a key.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

TYPING = Path(__file__).resolve().parents[2] / "scripts/typing"
sys.path.insert(0, str(TYPING))

import build_app_data as B  # noqa: E402
import subpanel_io as S  # noqa: E402

MARKERS = ["CD3e", "INS"]


def _cell_rows(n: int, *, ids: list[str] | None = None) -> pd.DataFrame:
    ids = ids if ids is not None else [f"uuid-{i:04d}" for i in range(n)]
    return pd.DataFrame({
        "object_id": ids,
        "X_centroid": [float(i) for i in range(n)],
        "Y_centroid": [float(i) * 2 for i in range(n)],
        "phenotype": ["Beta"] * n,
        "cell_region_20um": ["core", "expansion", "tissue"] * (n // 3) + ["core"] * (n % 3),
        "islet_key": ["Islet_1", "Islet_1", None] * (n // 3) + ["Islet_1"] * (n % 3),
        "cell_class": ["Endocrine"] * n,
    })


@pytest.fixture()
def lane(tmp_path, monkeypatch):
    """A two-donor typing tree and an aligned tree, wired into the module's own globals."""
    out = tmp_path / "typing"
    aligned = tmp_path / "uniform"
    for donor, n in (("6374", 9), ("6461", 6)):
        frame = _cell_rows(n)
        target = out / f"cells/donor_id={donor}"
        target.mkdir(parents=True)
        frame.to_parquet(target / "data.parquet", index=False)

        values = pd.DataFrame({"object_id": frame["object_id"]})
        for marker in MARKERS:
            values[marker] = 1.0
        aligned_target = aligned / f"donor_id={donor}"
        aligned_target.mkdir(parents=True)
        values.to_parquet(aligned_target / "data.parquet", index=False)

    monkeypatch.setattr(S, "OUT", out)
    monkeypatch.setattr(S, "ALIGNED", aligned)
    monkeypatch.setattr(S, "load_scale",
                        lambda: pd.DataFrame(index=pd.Index(MARKERS, name="marker")))
    monkeypatch.setattr(B, "APP", out / "app")
    return out


def test_donor_csv_leads_with_object_id(lane):
    B.write_donor_csv(["6374"])
    written = pd.read_csv(lane / "app/donors/6374.csv", dtype={"object_id": str})
    assert list(written.columns) == ["object_id", "X_centroid", "Y_centroid", "phenotype",
                                     "cell_region", "islet_name"]
    assert written["object_id"].is_unique
    assert written["object_id"].notna().all()


def test_drilldown_csv_leads_with_object_id(lane):
    B.write_drilldown_cells(["6374"])
    written = pd.read_csv(lane / "app/cells/6374_Islet_1.csv", dtype={"object_id": str})
    assert list(written.columns)[:5] == ["object_id", "X_centroid", "Y_centroid", "phenotype",
                                         "cell_region"]
    assert set(MARKERS) <= set(written.columns)
    assert written["object_id"].is_unique


def test_the_uuid_survives_the_join_it_was_made_on(lane):
    """The written UUIDs are the ones in the cell table, not a re-derived identifier."""
    B.write_drilldown_cells(["6374"])
    source = pd.read_parquet(lane / "cells/donor_id=6374/data.parquet")
    written = pd.read_csv(lane / "app/cells/6374_Islet_1.csv", dtype={"object_id": str})
    assert set(written["object_id"]) <= set(source["object_id"])


def test_a_duplicated_uuid_stops_the_writer(lane):
    """Two rows claiming one detection would multiply cells in any join. Refuse instead."""
    target = lane / "cells/donor_id=6374/data.parquet"
    frame = pd.read_parquet(target)
    frame.loc[frame.index[1], "object_id"] = frame.loc[frame.index[0], "object_id"]
    frame.to_parquet(target, index=False)
    with pytest.raises(ValueError, match="duplicated object_id"):
        B.write_donor_csv(["6374"])


def test_a_missing_uuid_stops_the_writer(lane):
    """A null UUID is a lost key, not a cell to write without one."""
    target = lane / "cells/donor_id=6374/data.parquet"
    frame = pd.read_parquet(target)
    frame.loc[frame.index[2], "object_id"] = None
    frame.to_parquet(target, index=False)
    with pytest.raises(ValueError, match="carry no object_id"):
        B.write_donor_csv(["6374"])


def test_an_empty_uuid_string_counts_as_missing(lane):
    """QuPath writes an empty field for an absent value; it is not an identifier."""
    target = lane / "cells/donor_id=6374/data.parquet"
    frame = pd.read_parquet(target)
    frame.loc[frame.index[3], "object_id"] = ""
    frame.to_parquet(target, index=False)
    with pytest.raises(ValueError, match="carry no object_id"):
        B.write_donor_csv(["6374"])


def test_the_check_runs_after_the_excluded_filter(lane):
    """A duplicate both of whose rows are Excluded never reaches a written file.

    The check is on the frame that is written, so it reports what the artifact contains rather
    than what the source table happened to hold.
    """
    target = lane / "cells/donor_id=6461/data.parquet"
    frame = pd.read_parquet(target)
    frame.loc[frame.index[0], "cell_class"] = "Excluded"
    frame.loc[frame.index[1], "cell_class"] = "Excluded"
    frame.loc[frame.index[1], "object_id"] = frame.loc[frame.index[0], "object_id"]
    frame.to_parquet(target, index=False)
    B.write_donor_csv(["6461"])
    written = pd.read_csv(lane / "app/donors/6461.csv", dtype={"object_id": str})
    assert len(written) == 4
    assert written["object_id"].is_unique
