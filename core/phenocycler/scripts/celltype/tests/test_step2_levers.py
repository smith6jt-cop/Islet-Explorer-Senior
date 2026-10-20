"""The two settings of the UniFORM step adjusted on 2026-10-02.

The maintainer: "UniFORM should be able to handle the different pixel intensity distributions, check if there are
levers that need to be adjusted."  Two were:

* a marker whose level steps with the acquisition group is left unaligned only when UniFORM's own shift ALSO leaves
  its background modes over the bar -- before, the step alone refused it, which with 25 donors in one cohort left
  CD8 unaligned although UniFORM brings it to 0.065;
* each marker can keep the reference donor of an earlier run, so the aligned values stay in that run's unit and a
  scale measured on it still applies.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import panel
import run_step2 as S2
import uniform_adapter as UA
import usability as U


@pytest.fixture(scope="module")
def pinned_uniform():
    if not UA.UNIFORM_PATH.exists():
        pytest.skip(f"UniFORM checkout not on disk at {UA.UNIFORM_PATH}")
    return UA._import_uniform()


@pytest.fixture
def cohort(monkeypatch):
    """Six donors in two acquisition groups. `level` is the same distribution on every donor, 6x lower on group B:
    a step with the group that is a pure scale. `shape` cannot be brought together by any shift: every donor has the
    same broad bulk, which is what a cross-correlation of whole histograms aligns, and a narrow peak -- the tallest
    bin, so the "background mode" the bar is stated on -- that sits at a different place on group B's donors."""
    rng = np.random.default_rng(0)
    donors = ["a1", "a2", "a3", "b1", "b2", "b3"]
    base = rng.lognormal(4.0, 0.5, 60_000)
    level = {d: base * (1.0 if d[0] == "a" else 1 / 6.0) * f for d, f in zip(donors, (1.0, 1.2, 0.9, 1.0, 1.1, 0.95))}
    peak_at = {"a1": 4.0, "a2": 4.0, "a3": 4.0, "b1": 3.2, "b2": 4.0, "b3": 4.8}
    shape = {d: np.concatenate([rng.lognormal(4.0, 0.3, 57_000), rng.lognormal(peak_at[d], 0.012, 4_500)]) for d in donors}
    values = {"level": level, "shape": shape}
    monkeypatch.setattr(S2, "marker_values", lambda m, ds, unavailable=None: {d: values[m][d] for d in ds})
    monkeypatch.setattr(S2.UA, "gmm_cutoff", lambda *a, **k: float("nan"))
    usab = pd.DataFrame({"marker": ["level", "shape"], "acq_ratio": [5.0, 5.0]})        # both over the flag of 2.0
    assert (usab["acq_ratio"] >= U.BARS["acq_ratio_flag"]).all()
    return donors, usab


@pytest.mark.data
def test_a_stepped_marker_uniform_aligns_within_the_bar_is_aligned(pinned_uniform, cohort, tmp_path):
    donors, usab = cohort
    r = S2.derive(["level"], donors, usab, save_figs=False, figs=tmp_path / "f")["level"]
    assert r["acq_ratio_flagged"] is True and r["mode_sd_after"] <= S2.BAR_MODE_SD
    assert r["provenance"] == "automatic_registration"
    ref = r["reference_donor"]
    for d in ("b1", "b2", "b3"):                   # group B comes up ~6x onto group A's level (or A down, if B holds the reference)
        assert r["scale"][d] / r["scale"]["a1"] == pytest.approx(6.0, rel=0.25), (d, r["scale"])
    assert r["scale"][ref] == 1.0


@pytest.mark.data
def test_a_stepped_marker_uniform_cannot_align_is_left_as_it_is(pinned_uniform, cohort, tmp_path):
    donors, usab = cohort
    r = S2.derive(["shape"], donors, usab, save_figs=False, figs=tmp_path / "f")["shape"]
    assert r["mode_sd_after"] > S2.BAR_MODE_SD
    assert r["provenance"] == S2.REFUSED == "refused_reagent_change"        # the label refine.py and kintsugi read
    assert set(r["scale"].values()) == {1.0}, "a refused marker is left unshifted, not shifted and labelled"
    # the same distributions WITHOUT the step with the group are carried over the bar, not refused
    calm = usab.assign(acq_ratio=[0.5, 0.5])
    r2 = S2.derive(["shape"], donors, calm, save_figs=False, figs=tmp_path / "f")["shape"]
    assert r2["provenance"] == "automatic_registration_over_bar" and r2["acq_ratio_flagged"] is False


@pytest.mark.data
def test_a_pinned_reference_donor_sets_the_unit(pinned_uniform, cohort, tmp_path):
    donors, usab = cohort
    free = S2.derive(["level"], donors, usab, save_figs=False, figs=tmp_path / "f")["level"]
    assert free["reference_source"] == "uniform_choice"
    other = next(d for d in donors if d != free["reference_donor"] and d[0] != free["reference_donor"][0])
    pinned = S2.derive(["level"], donors, usab, save_figs=False, figs=tmp_path / "f", references={"level": other})["level"]
    assert pinned["reference_donor"] == other and pinned["reference_source"] == "pinned" and pinned["scale"][other] == 1.0
    # the same alignment in another unit: every donor's scale moves by one common factor
    ratio = [pinned["scale"][d] / free["scale"][d] for d in donors]
    assert max(ratio) / min(ratio) == pytest.approx(1.0, abs=0.03) and ratio[0] != pytest.approx(1.0, abs=0.5)
    # a pinned donor that has no histogram for the marker falls back to UniFORM's choice
    gone = S2.derive(["level"], donors, usab, save_figs=False, figs=tmp_path / "f", references={"level": "zz"})["level"]
    assert gone["reference_source"] == "uniform_choice" and gone["reference_donor"] == free["reference_donor"]


def test_references_are_read_from_an_earlier_run_and_donor_names_are_matched(tmp_path, monkeypatch):
    monkeypatch.setattr(panel, "cohort", lambda: pd.DataFrame({"donor_id": ["0112", "115", "6591"]}))
    monkeypatch.setattr(panel, "donor_markers", lambda d: ["CD3e", "CD8", "NEW"])
    monkeypatch.setattr(panel, "expression_path", lambda d: tmp_path / f"{d}.parquet")
    monkeypatch.setattr(panel, "OUT", tmp_path)
    pd.DataFrame({"marker": ["CD3e", "CD8", "NEW"], "acq_ratio": [0.5, 3.0, 0.4]}).to_parquet(tmp_path / "marker_donor_usability.parquet")
    earlier = tmp_path / "earlier_shifts.json"
    earlier.write_text(json.dumps({"shifts": {"CD3e": {"reference_donor": "112"}, "CD8": {"reference_donor": "6523"}}}))
    seen = {}

    def derive(markers, donors, usab, *, figs, references=None):
        seen.update(references=references, markers=markers)
        return {m: {"provenance": "automatic_registration", "mode_sd_after": 0.01, "mode_sd_before": 0.5, "n_donors": 3,
                    "scale": {d: 1.0 for d in donors}} for m in markers}
    monkeypatch.setattr(S2, "derive", derive)
    out = tmp_path / "uniform_x"
    assert S2.main(["--out", str(out), "--dry-run", "--references-from", str(earlier)]) == 0
    # "112" in the earlier run is this cohort's "0112"; 6523 is not in this cohort, so CD8 is left to UniFORM
    assert seen["references"] == {"CD3e": "0112"}
    assert json.loads((out / "shifts.json").read_text())["references_from"] == str(earlier)
