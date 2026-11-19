"""Refit the CODEX cohort and verify the released ordering by islet identifier."""
from pathlib import Path
import numpy as np
import pandas as pd
import paga_pseudotime as engine

root = Path(__file__).resolve().parent
engine._FEATURES_CSV = str(root / "barlow_islet_features.csv")
result = engine.compute_pseudotime(min_cells=10)
features = result["features"]
assert len(features) == 56
context = engine._context(10, None)
actual = pd.Series(result["pseudotime"], index=context["sub"].islet_id)
expected = pd.read_csv(root / "diffusion_coordinates.csv").set_index("islet_id").pseudotime_oriented
assert actual.index.is_unique and set(actual.index) == set(expected.index)
np.testing.assert_allclose(actual.loc[expected.index], expected, atol=1e-6, rtol=1e-6)
print(f"Verified {len(actual)} CODEX islets and {len(features)} graph inputs.")
