"""Which cells the downstream immune-marker rules should evaluate."""
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/typing"))
import immune_subtype


def test_final_immune_subtyping_uses_promotions_and_excludes_demotions(tmp_path):
    path = tmp_path / "cells.parquet"
    pd.DataFrame({"object_id": ["gained", "lost", "kept"], "cell_class": ["Exocrine", "Immune", "Immune"],
                  "refined_class": ["Immune", "Exocrine", "Immune"]}).to_parquet(path, index=False)
    cells, selected, stage = immune_subtype.immune_parent_cells(path)
    assert stage == "refined" and cells.loc[selected, "object_id"].tolist() == ["gained", "kept"]
    cells, selected, stage = immune_subtype.immune_parent_cells(path, primary=True)
    assert stage == "primary" and cells.loc[selected, "object_id"].tolist() == ["lost", "kept"]


def test_pre_refinement_inputs_are_primary_and_partial_results_fail(tmp_path):
    path = tmp_path / "cells.parquet"
    frame = pd.DataFrame({"object_id": ["a", "b"], "cell_class": ["Immune", "Exocrine"]})
    for old_format in (True, False):
        if not old_format:
            frame["refined_class"] = None
            frame["refine_state"] = "not_run"
        frame.to_parquet(path, index=False)
        _, mask, stage = immune_subtype.immune_parent_cells(path)
        assert stage == "primary" and mask.tolist() == [True, False]
    frame.loc[0, "refined_class"] = "Immune"
    frame.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="incomplete final"):
        immune_subtype.immune_parent_cells(path)
