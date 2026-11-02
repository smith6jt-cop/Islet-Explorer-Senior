"""positivity.lineage_levels: per-section levels from mutually exclusive lineages (spleen, maintainer 2026-11-01)."""

from __future__ import annotations

import csv

import numpy as np
import pandas as pd
import pytest

from multiplex5000.positivity import lineage_levels as LL
from multiplex5000.positivity import registry as RG

COLUMNS = ["marker", "channel", "on_every_image", "role", "compartment", "lineage", "shared_by", "required",
           "measurement", "measurement_candidates", "differentiation", "positivity_class", "background_method",
           "level_reference", "carrier_lineages", "source", "proposal"]
ROWS = [  # marker, role, compartment, lineage, required, measurement, differentiation, method, reference
    ("CD20", "defining", "Hematopoietic", "B", "yes", "whole_cell", "", "known_negative", "t_cells"),
    ("CD3e", "defining", "Hematopoietic", "T", "yes", "whole_cell", "membrane", "known_negative", "b_cells"),
    ("CD8", "defining", "Hematopoietic", "T:CD8", "no", "whole_cell", "", "known_negative", "b_cells"),
    ("CD45", "defining", "Hematopoietic", "Leukocyte", "yes", "whole_cell", "", "known_negative", "non_haematopoietic"),
    ("SMA", "defining", "Mesenchymal", "Mural", "no", "whole_cell", "", "known_negative", "lymphocytes"),
    ("Ki67", "state", "", "", "no", "nucleus", "", "valley", ""),
    ("CD44", "state", "", "", "no", "whole_cell", "", "background", ""),
]
CONSTANTS = {"ISOLATION_UM": 10, "LEVEL_QUANTILE": 99, "KNOWN_NEGATIVE_MIN_ABOVE": 10, "BACKGROUND_K": 3,
             "POSITIVE_LINEAGE_MIN": 0.5}


@pytest.fixture()
def reg_name(tmp_path):
    with (tmp_path / "test__spleen.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for m, role, comp, lin, req, meas, diff, method, ref in ROWS:
            w.writerow(dict(dict.fromkeys(COLUMNS, ""), marker=m, channel=m, on_every_image="yes", role=role,
                            compartment=comp, lineage=lin, required=req, measurement=meas, measurement_candidates=meas,
                            differentiation=diff, positivity_class=method, background_method=method,
                            level_reference=ref))
    with (tmp_path / "test__spleen__constants.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "value", "source", "reason", "proposal"])
        for k, v in CONSTANTS.items():
            w.writerow([k, v, "test", "", ""])
    return ("test", "spleen", str(tmp_path))


def _section(seed=0, n_b=4000, n_t=4000, n_s=1500):
    """B cells on the left, T cells on the right, stromal (SMA+, CD45-) cells in a strip at x 1000-1100 um."""
    rng = np.random.default_rng(seed)
    n = n_b + n_t + n_s
    lin = np.array(["B"] * n_b + ["T"] * n_t + ["S"] * n_s)
    x = np.where(lin == "B", rng.uniform(0, 900, n), np.where(lin == "T", rng.uniform(1200, 2100, n), rng.uniform(1000, 1100, n)))
    y = rng.uniform(0, 1000, n)
    bg = lambda: np.exp(rng.normal(5.0, 0.3, n))
    hi = lambda: np.exp(rng.normal(8.0, 0.3, n))
    d = pd.DataFrame({"x_um": x, "y_um": y})
    d["Cell__CD20"] = np.where(lin == "B", hi(), bg())
    d["Cell__CD3e"] = np.where(lin == "T", hi(), bg())
    d["Membrane__CD3e"] = np.where(lin == "T", hi(), bg())
    cd8 = (lin == "T") & (rng.random(n) < 0.4)
    d["Cell__CD8"] = np.where(cd8, hi(), bg())
    d["Cell__CD45"] = np.where(lin == "S", bg(), hi())
    d["Cell__SMA"] = np.where(lin == "S", hi(), bg())
    ki = rng.random(n) < 0.15
    d["Nucleus__Ki67"] = np.where(ki, hi(), bg())
    d["Cell__CD44"] = bg()
    return d, lin, cd8


class Source:
    def __init__(self, frames):
        self.frames = frames

    def levels_frame(self, tag, cols, ctx):
        return self.frames[tag][["x_um", "y_um"] + cols].copy()


def test_levels_come_from_the_lineage_that_cannot_carry_the_marker(reg_name):
    d, lin, cd8 = _section()
    t = LL.choose(Source({"s1": d}), reg_name, ["s1"])
    lv = t["levels"].set_index("marker")
    assert set(t["selection"]["decision"]) == {"use"} and not t["stops"]
    # CD3e's level is the 99th percentile of CD3e on B cells: ~1 % of B cells above it, every T cell above it
    pos3 = np.log1p(d["Cell__CD3e"]) > lv.loc["CD3e", "level_log"]
    assert pos3[lin == "T"].mean() > 0.99 and pos3[lin == "B"].mean() < 0.03
    pos8 = np.log1p(d["Cell__CD8"]) > lv.loc["CD8", "level_log"]
    assert pos8[cd8].mean() > 0.99 and pos8[lin == "B"].mean() < 0.03
    # CD45 on non-haematopoietic cells (the SMA+ strip, away from lymphocytes): every lymphocyte above it
    pos45 = np.log1p(d["Cell__CD45"]) > lv.loc["CD45", "level_log"]
    assert pos45[lin != "S"].mean() > 0.99 and pos45[lin == "S"].mean() < 0.05
    # Ki67 by its valley, CD44 (no positives) by its background mode + K half-widths
    pki = np.log1p(d["Nucleus__Ki67"]) > lv.loc["Ki67", "level_log"]
    assert abs(pki.mean() - 0.15) < 0.02
    assert lv.loc["CD44", "positive_fraction"] < 0.01
    # the differentiation reading of CD3e has its own level, and both passes found the lineages
    assert list(t["levels_differentiation"]["marker"]) == ["CD3e"]
    assert t["chosen"]["CD3e"] == ("membrane", "known_negative") and t["positivity"]["CD3e"][0] == "whole_cell"
    opts = t["options"]
    assert opts["b_cells"].iloc[0] > 3900 and opts["t_cells"].iloc[0] > 3900


def test_a_neighbour_spillover_is_isolated_away(reg_name):
    """B cells touching T cells carry CD3e spillover; the level ignores B cells within ISOLATION_UM of a positive."""
    d, lin, _ = _section(seed=1)
    near = (lin == "B") & (d["x_um"] > 880)                   # the B cells next to the T zone's edge strip
    rng = np.random.default_rng(2)
    t_edge = (lin == "B") & (d["x_um"] > 895)
    d.loc[t_edge, "Cell__CD3e"] = np.exp(rng.normal(8.0, 0.3, t_edge.sum()))     # a few bright 'T-like' cells
    d.loc[near & ~t_edge, "Cell__CD3e"] = np.exp(rng.normal(7.0, 0.3, (near & ~t_edge).sum()))   # their spillover
    lvl = LL.choose(Source({"s1": d}), reg_name, ["s1"])["levels"].set_index("marker").loc["CD3e", "level_log"]
    far_b = (lin == "B") & (d["x_um"] < 800)
    assert (np.log1p(d["Cell__CD3e"]) > lvl)[far_b].mean() < 0.02
    assert lvl < 7.0                                          # the spillover band did not set the level


def test_a_method_without_a_level_is_reported_for_every_section(reg_name):
    d1, _, _ = _section(seed=3)
    d2, _, _ = _section(seed=4)
    d2["Nucleus__Ki67"] = np.exp(np.random.default_rng(5).normal(5.0, 0.3, len(d2)))      # unimodal: no valley
    t = LL.choose(Source({"a": d1, "b": d2}), reg_name, ["a", "b"])
    sel = t["selection"].set_index("marker")
    assert sel.loc["Ki67", "decision"] == "stop" and "valley gives no level on: b" in sel.loc["Ki67", "why"]
    assert t["stops"] == ["Ki67"]
    can = LL.canary_levels(t, reg_name)
    assert set(can["marker"]) == {"CD20", "CD3e", "CD8", "CD45", "SMA"}        # typing markers only
    assert (can["canary_log"] == can["level_log"]).all()


def test_lineages_must_be_bimodal(reg_name):
    d, _, _ = _section(seed=6)
    d["Cell__CD20"] = np.exp(np.random.default_rng(7).normal(5.0, 0.3, len(d)))
    with pytest.raises(LL.LineageLevelError, match="no population above its background"):
        LL.image_levels(Source({"s": d}), "s", reg_name)


def test_a_shoulder_population_still_splits_the_lineage():
    """HDL080's CD20: the B cells are a shoulder of the background, not a second peak."""
    rng = np.random.default_rng(8)
    v = np.log1p(np.concatenate([np.exp(rng.normal(5.8, 0.45, 80000)), np.exp(rng.normal(7.0, 0.4, 15000)),
                                 np.exp(rng.normal(1.0, 0.5, 3000))]))     # plus a low tail of REDSEA near-zeros
    cut = LL.lineage_split(v)
    assert 6.0 < cut < 7.0
