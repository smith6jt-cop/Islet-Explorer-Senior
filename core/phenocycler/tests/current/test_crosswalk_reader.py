"""The vocabulary crosswalk's lineage rows are read one way, strictly.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_crosswalk_reader.py -q

Until 2026-09-24 three builders each read ``vocab_crosswalk.csv`` their own way; one matched
``is_immune`` case-insensitively and wrote a blank scope as the string ``"nan"``.  They share
``scripts/senior/crosswalk.py`` now.  The crosswalk is a tracked file, so the checks against the
real one run in the default lane.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior import build_compartment_gradient as compartments  # noqa: E402
from scripts.senior import build_grid_composition as grid  # noqa: E402
from scripts.senior import crosswalk as XW  # noqa: E402
from scripts.senior import modality_concordance as concordance  # noqa: E402

HEADER = "tissue,kind,common,phenocycler,xenium,is_immune,immune_tissue_scope\n"


def _write(tmp_path, body: str) -> Path:
    path = tmp_path / "vocab_crosswalk.csv"
    path.write_text(HEADER + body)
    return path


def test_a_blank_is_an_empty_string_and_a_value_is_what_it_says(tmp_path):
    path = _write(tmp_path, "pancreas,lineage,T_NK,Immune,T_NK,yes,\n"
                            "pancreas,lineage,B_Plasma,Immune,B_Plasma,Yes,pancreas\n"
                            "pancreas,marker,CD3,CD3,CD3E,,\n")
    rows = XW.lineage_rows(path)
    assert list(rows["common"]) == ["T_NK", "B_Plasma"]
    assert rows.loc[rows["common"].eq("T_NK"), "immune_tissue_scope"].item() == ""
    assert not rows.isna().any().any()
    assert list(rows.loc[rows["is_immune"].eq("yes"), "common"]) == ["T_NK"], (
        "'Yes' is not 'yes'; the crosswalk is generated and never spells it that way")


def test_a_crosswalk_with_no_lineage_rows_is_refused_by_every_reader(tmp_path, monkeypatch):
    path = _write(tmp_path, "pancreas,marker,CD3,CD3,CD3E,,\n")
    with pytest.raises(XW.CrosswalkError, match="no lineage rows"):
        XW.lineage_rows(path)
    monkeypatch.setattr(grid, "CROSSWALK", path)
    with pytest.raises(grid.GridError, match="no lineage rows"):
        grid.load_crosswalk()
    monkeypatch.setattr(concordance, "CROSSWALK", path)
    with pytest.raises(XW.CrosswalkError):
        concordance.load_crosswalk()
    monkeypatch.setattr(compartments, "CROSSWALK", path)
    with pytest.raises(XW.CrosswalkError):
        compartments.load_crosswalk_immune()


def test_the_three_builders_read_the_real_crosswalk_through_the_one_reader():
    for module in (grid, concordance, compartments):
        assert Path(module.CROSSWALK) == XW.CROSSWALK, module.__name__
    rows = XW.lineage_rows()
    pd.testing.assert_frame_equal(grid.load_crosswalk(), rows)
    pd.testing.assert_frame_equal(concordance.load_crosswalk(), rows)
    immune, scope = compartments.load_crosswalk_immune()
    assert immune == set(rows.loc[rows["is_immune"].eq("yes"), "common"])
    assert set(scope) == immune
    assert "nan" not in scope.values()


def test_no_builder_reads_the_crosswalk_file_itself():
    """A fourth private reader is how three came to disagree."""
    for module in (grid, concordance, compartments):
        tree = ast.parse(Path(module.__file__).read_text())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "read_csv" and node.args
                    and isinstance(node.args[0], ast.Name) and node.args[0].id == "CROSSWALK"):
                pytest.fail(f"{module.__name__} reads CROSSWALK itself at line {node.lineno}")
