"""The UniFORM audit (Step 8e) on synthetic donors whose answer is planted.

Nothing here reads the cohort: every donor, shift and metadata row is built in the test.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import uniform_audit as UA                          # noqa: E402

from phenocycler.cohort import DONOR_EXCLUSIONS     # noqa: E402

HORMONES = ["INS", "GCG", "SST"]
BATCH_GAIN = 3.0


def _donor_cells(donor: str, gain: float, rng: np.random.Generator, n: int = 3000) -> pd.DataFrame:
    """Beta, alpha and delta cells (one hormone each) plus hormone-negative reference tissue.

    Every value is multiplied by ``gain`` -- the planted batch effect.
    """
    groups = {"INS": n, "GCG": n, "SST": n, None: 2 * n}
    frames = []
    for hot, k in groups.items():
        t = {m: rng.lognormal(0.0, 0.3, k) for m in HORMONES}
        if hot is not None:
            t[hot] = rng.lognormal(math.log(8.0), 0.3, k)
        t["reference_population"] = np.full(k, hot is None)
        frames.append(pd.DataFrame(t))
    df = pd.concat(frames, ignore_index=True)
    for m in HORMONES:
        df[m] *= gain
    df["donor_id"] = donor
    df["object_id"] = [f"{donor}-{i}" for i in range(len(df))]
    df["qc_analysis_eligible"] = True
    return df


@pytest.fixture
def cohort(tmp_path):
    """Eight donors in two stain batches; batch S2 is acquired 3x brighter on every marker."""
    rng = np.random.default_rng(7)
    donors = {f"90{i:02d}": ("S1" if i < 4 else "S2") for i in range(8)}
    disease = {d: ("ND" if i % 2 == 0 else "T1D") for i, d in enumerate(donors)}
    root = tmp_path / "expression"
    for d, batch in donors.items():
        part = root / f"donor_id={d}"
        part.mkdir(parents=True)
        _donor_cells(d, BATCH_GAIN if batch == "S2" else 1.0, rng).drop(
            columns="donor_id").to_parquet(part / "data.parquet", index=False)
    shift = {d: (math.log(BATCH_GAIN) if b == "S2" else 0.0) for d, b in donors.items()}
    shifts = {"uniform_commit": "synthetic", "shifts": {
        m: {"provenance": "automatic_registration", "shift_log": shift,
            "scale": {d: math.exp(-s) for d, s in shift.items()}} for m in HORMONES}}
    shifts["shifts"]["CD38"] = {"provenance": "refused_reagent_change", "shift_log": {},
                                "scale": {d: 1.0 for d in donors}}
    (tmp_path / "shifts.json").write_text(json.dumps(shifts))
    pd.DataFrame({"donor_id": list(donors), "donor_status": [disease[d] for d in donors]}
                 ).to_csv(tmp_path / "meta.csv", index=False)
    pd.DataFrame({"donor_id": list(donors), "acq_run": [f"2026-01-0{i}" for i in range(8)],
                  "stain_group": list(donors.values())}).to_csv(tmp_path / "runs.csv",
                                                                 index=False)
    pd.DataFrame({"marker_a": ["INS", "INS", "GCG"], "marker_b": ["GCG", "SST", "SST"]}
                 ).to_csv(tmp_path / "pairs.csv", index=False)
    pd.DataFrame({"marker": HORMONES, "floor": 2.0, "top": 12.0, "span": 10.0}
                 ).to_parquet(tmp_path / "scale.parquet", index=False)
    return tmp_path, root


def _run(tmp_path, root, **kw):
    return UA.run([root], tmp_path / "shifts.json", tmp_path / "meta.csv", tmp_path / "out",
                  acquisition_runs=tmp_path / "runs.csv", pairs_path=tmp_path / "pairs.csv",
                  scale_path=tmp_path / "scale.parquet",
                  reference_column="reference_population", n_perm=199, **kw)


# ------------------------------------------------------------------- (a) residual batch
def test_planted_batch_shift_is_detected_and_gone_after_alignment(cohort):
    tmp_path, root = cohort
    res = _run(tmp_path, root, apply_shifts=True)
    s = res["restore_summary"].set_index(["values", "marker", "partner"])
    before = s.xs("pre_alignment", level="values")
    after = s.xs("aligned", level="values")
    assert before["residual_batch"].all(), before
    assert (before["batch_spread_log_ratio"] > 0.8).all(), "a 3x gain is ~1.1 in log units"
    assert not after["residual_batch"].any(), after
    assert (after["batch_spread_log_ratio"] < before["batch_spread_log_ratio"] / 4).all()
    th = res["thresholds"]
    pre = th.loc[th["values"] == "pre_alignment"]
    assert pre["frac_delta"].abs().max() > 0.1, "pooled threshold misses a batch's positives"
    assert (res["partner_levels"]["level"] > 0).all()
    assert (tmp_path / "out" / "SUMMARY.md").read_text().startswith("# UniFORM audit")


def test_copositivity_uses_a_different_partner_and_stays_low_on_exclusive_markers(cohort):
    tmp_path, root = cohort
    cp = _run(tmp_path, root)["copositivity"]
    assert len(cp) and (cp["copos_partner"] != cp["partner"]).all()
    assert (cp["copos_partner"] != cp["marker"]).all()
    assert cp["copos_donor"].median() < 0.05


def test_circularity_guard_refuses_j2_equal_to_j1_or_the_marker():
    with pytest.raises(UA.CircularPartnerError):
        UA.check_partners("INS", "GCG", "GCG")
    with pytest.raises(UA.CircularPartnerError):
        UA.check_partners("INS", "GCG", "INS")
    with pytest.raises(UA.CircularPartnerError):
        UA.check_partners("INS", "INS", None)
    triples = UA.partner_triples(UA.load_pairs())
    assert (triples["copos_partner"].dropna() != triples.loc[
        triples["copos_partner"].notna(), "partner"]).all()
    bad = pd.DataFrame({"marker": ["INS"], "partner": ["GCG"], "copos_partner": ["GCG"]})
    frame = pd.DataFrame({"donor_id": ["1"], "INS": [1.0], "GCG": [1.0]})
    with pytest.raises(UA.CircularPartnerError):
        UA.restore_test(frame, bad, pd.DataFrame({"GCG": [True]}))


def test_a_marker_with_one_partner_gets_no_copositivity_fallback():
    triples = UA.partner_triples(pd.DataFrame({"marker_a": ["CD3e"], "marker_b": ["CD20"]}))
    assert triples["copos_partner"].isna().all()


def test_partner_positivity_needs_a_reference(cohort):
    tmp_path, root = cohort
    with pytest.raises(ValueError, match="reference"):
        UA.run([root], tmp_path / "shifts.json", tmp_path / "meta.csv", tmp_path / "out")


# ------------------------------------------------------------------- exclusivity ratio
def test_exclusivity_ratio_near_one_for_exclusive_and_near_zero_for_correlated():
    rng = np.random.default_rng(0)
    n = 2000
    a = np.r_[rng.uniform(0.5, 1.0, n), np.zeros(n)]
    b = np.r_[np.zeros(n), rng.uniform(0.5, 1.0, n)]
    assert UA.exclusivity_ratio(a, b) > 0.9
    t = rng.uniform(0.0, 1.0, 2 * n)
    assert UA.exclusivity_ratio(t, t + rng.normal(0, 0.01, 2 * n)) < 0.05
    assert math.isnan(UA.exclusivity_ratio(a[:10], b[:10])), "too few cells is unavailable"


def test_exclusivity_on_cohort_scale_for_synthetic_hormones(cohort):
    tmp_path, root = cohort
    ex = _run(tmp_path, root, apply_shifts=True)["exclusivity"]
    assert len(ex) == 3 * 8
    assert (ex["r"] > 0.8).all(), ex


# ------------------------------------------------------------------- (b) invariance
def _shift_frame(effect: dict[str, float], donors: pd.DataFrame, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    batch_off = {b: 0.4 * k for k, b in enumerate(sorted(donors["batch"].unique()))}
    rows = []
    for m, eff in effect.items():
        for d, dis, b in donors[["donor_id", "disease", "batch"]].itertuples(index=False):
            rows.append({"marker": m, "donor_id": d, "scale": 1.0, "provenance": "x",
                         "shift_log": batch_off[b] + (eff if dis == "T1D" else 0.0)
                         + rng.normal(0, 0.03)})
    return pd.DataFrame(rows)


def _mixed_design() -> pd.DataFrame:
    ids = [f"91{i:02d}" for i in range(12)]
    return pd.DataFrame({"donor_id": ids,
                         "disease": ["ND", "ND", "T1D", "T1D"] * 3,
                         "batch": ["S1"] * 4 + ["S2"] * 4 + ["S3"] * 4})


def test_planted_disease_shift_is_flagged_and_a_batch_only_shift_is_not():
    donors = _mixed_design()
    inv = UA.invariance_test(_shift_frame({"HLA_A": 0.8, "CD3e": 0.0}, donors), donors,
                             n_perm=999).set_index("marker")
    assert inv.loc["HLA_A", "flag"] == "disease_associated"
    assert inv.loc["HLA_A", "effect_T1D"] == pytest.approx(0.8, abs=0.1)
    assert bool(inv.loc["HLA_A", "watchlist"])
    assert inv.loc["CD3e", "flag"] == "not_flagged"
    assert inv.loc["CD3e", "status"] == "tested"


def test_confounded_design_is_reported_not_passed():
    donors = pd.DataFrame({"donor_id": [f"92{i:02d}" for i in range(8)],
                           "disease": ["ND"] * 4 + ["T1D"] * 4,
                           "batch": ["S1"] * 4 + ["S2"] * 4})
    inv = UA.invariance_test(_shift_frame({"HLA_A": 0.8}, donors), donors, n_perm=99)
    assert (inv["status"] == "confounded").all()
    assert (inv["flag"] == "confounded").all()
    assert inv["p_perm"].isna().all()
    design = UA.design_table(donors).set_index("batch")
    assert not design.loc["S2", "has_reference_disease"]
    assert not design["within_batch_contrast"].any()
    md = UA.summary_markdown(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), inv,
                             design.reset_index(), pd.DataFrame(), {})
    assert "not estimable" in md and "Batches with no ND donor" in md


def test_partially_confounded_batch_is_named_but_others_still_test():
    donors = _mixed_design()
    donors.loc[donors["batch"] == "S3", "disease"] = "T1D"
    design = UA.design_table(donors).set_index("batch")
    assert not design.loc["S3", "has_reference_disease"]
    inv = UA.invariance_test(_shift_frame({"HLA_A": 0.8}, donors), donors, n_perm=199)
    assert inv.loc[0, "status"] == "tested"
    assert inv.loc[0, "batches_with_contrast"] == "S1,S2"


def test_shifts_json_refused_markers_are_listed_not_tested(cohort):
    tmp_path, _ = cohort
    tidy, skipped = UA.load_shifts(tmp_path / "shifts.json")
    assert set(tidy["marker"]) == set(HORMONES)
    assert skipped["marker"].tolist() == ["CD38"]


# ------------------------------------------------------------------- input contracts
def test_unqualified_cells_and_excluded_donors_never_enter(tmp_path):
    excluded = sorted(DONOR_EXCLUSIONS)[0]
    rng = np.random.default_rng(3)
    for d in ("9300", excluded):
        part = tmp_path / f"donor_id={d}"
        part.mkdir()
        df = _donor_cells(d, 1.0, rng, n=10).drop(columns="donor_id")
        df.loc[:4, "qc_analysis_eligible"] = False
        df.loc[5, "INS"] = np.nan
        df.to_parquet(part / "data.parquet", index=False)
    f = UA.load_expression([tmp_path], HORMONES + ["CD45"])
    assert set(f["donor_id"]) == {"9300"}
    assert f.attrs["dropped"] == [excluded]
    assert len(f) == 50 - 5
    assert f["CD45"].isna().all(), "an unmeasured marker is unavailable, not zero"
    assert f["INS"].isna().sum() == 1


def test_missing_qc_column_is_refused(tmp_path):
    part = tmp_path / "donor_id=9301"
    part.mkdir()
    pd.DataFrame({"object_id": ["a"], "INS": [1.0]}).to_parquet(part / "data.parquet")
    with pytest.raises(KeyError, match="qc_analysis_eligible"):
        UA.load_expression([tmp_path], ["INS"])


def test_fundamental_pair_json_is_an_accepted_pair_table():
    path = pathlib.Path(UA.__file__).resolve().parents[2] / "phenocycler" / \
        "fundamental_pair_diagnostics.json"
    pairs = UA.load_pairs(path)
    assert {tuple(p) for p in pairs.to_numpy()} >= {("CD3e", "Ker8_18")}
