"""The gate must verify a planted truth AND refuse noise -- or the bars measure the gate.

The load-bearing test is the positive control: a synthetic donor whose labels are patchy in
space and whose declared markers carry real tails has to come out VERIFIED for every planted
label, through the same ``judge_label`` the screen driver calls.  A gate that only ever fails
is not a gate ("or the bars are measuring the gate rather than the arm").

The operational no-expected-frequency test is the scramble: rebuild that donor ten times with
label shares drawn from a Dirichlet and demand the verdicts be identical every time.  Every
statistic here compares a label to a null drawn from the same cells; a share enters only
through the evidence floors (counts on a statistic's n), so scrambling it must change nothing.

Everything else is the fail-closed half: no structure and no tails -> not verified; too few
calls or too few tail cells -> unjudgeable, which never clears; the null-gated cluster
annotation recovers planted blocks and its false-award rate on pure noise is MEASURED and
PRINTED, not tuned away.  Synthetic data only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import gate
from conftest import signature_arm, subsample_replicate
from criteria import CRITERIA
from records import AMBIGUOUS, UNAVAILABLE


# --------------------------------------------------------------------------- helpers
def judge_all(dec, donor, arm, *, n_replicates: int, seed: int = 0, with_cell_region: bool = True):
    """Run ``arm`` (a callable on the aligned values) as the screen driver would: primary, one
    ablation per identity marker, ``n_replicates`` 90 % hash-subsample replicates; then
    ``judge_label`` for every declared label.  Returns (primary labels, {label: result})."""
    primary = np.asarray(arm(donor.X)).astype(str)
    ablations = {m: np.asarray(arm(donor.X.drop(columns=[m]))).astype(str) for m in donor.X.columns}
    replicates = {f"subsample:{i}": subsample_replicate(primary, donor.ids, CRITERIA["subsample_frac"],
                                                        salt=f"subsample:{i}|{i}")
                  for i in range(1, n_replicates + 1)}
    tails = {m: gate.tail_mask(donor.X[m].to_numpy()) for m in donor.X.columns}
    out = {}
    for L, spec in dec["labels"].items():
        other = {m: t for m, t in tails.items() if m not in spec["markers"]}
        out[L] = gate.judge_label(
            L, primary, ablations, replicates, values=donor.X, xy=donor.xy, groups=donor.groups,
            ids=donor.ids, declared_markers=spec["markers"], held_out_gated=spec["held_out_gated"],
            other_tails=other, cell_region=donor.cell_region if with_cell_region else None, seed=seed)
    return primary, out


def _verdict_tuple(r: dict) -> tuple:
    return (r["heldout"]["verified"], r["spatial"]["clear"], r["stability"]["invariant"], r["verified"])


def _summary_rows(res: dict) -> list:
    rows = []
    for L, r in res.items():
        ho = r["heldout"]
        margins = [x["h"] - x["null_p"] for x in ho["per_marker"] if x.get("judgeable")]
        rows.append({
            "label": L, "verified": r["verified"],
            "heldout": (f"{ho['n_clear']}/{ho['n_judged']} clear, min margin {min(margins):+.2f}"
                        if ho["verified"] is not None and margins else "not gated"),
            "spatial": f"I={r['spatial']['I']:.3f} p95={r['spatial']['null_p']:.3f} clear={r['spatial']['clear']}",
            "stability": r["stability"]["invariant"],
        })
    return rows


def _draw_shares(rng, labels, floor: float, cap: float) -> dict:
    """Dirichlet shares with a floor and a cap.  The floor keeps every label above the
    EVIDENCE floors (min_calls, tail_min_k) so the only thing varying is prevalence; the cap
    keeps the enrichment ratio's own ceiling 1/s from meeting the null's -- a property of the
    statistic, not an expected prevalence."""
    k = len(labels)
    while True:
        d = rng.dirichlet(np.ones(k))
        s = floor + (1 - k * floor) * d
        if s.max() <= cap:
            return dict(zip(labels, s))


# --------------------------------------------------------------------------- (a) positive control
def test_positive_control_planted_labels_are_verified(dec, planted_donor, capsys):
    """~20k cells, 4 regions, patchy planted labels, planted marker tails; the toy arm reads them
    back from expression alone and every planted label must be VERIFIED."""
    donor = planted_donor(n=20_000, n_regions=4, seed=0)
    primary, res = judge_all(dec, donor, lambda X: signature_arm(X, dec),
                             n_replicates=CRITERIA["n_replicates"])
    acc = float((primary == donor.labels).mean())
    assert acc > 0.97, f"the toy arm did not recover the planted truth ({acc:.3f}); the fixture is broken"
    with capsys.disabled():
        print(f"\n[positive control] n={donor.n} regions={len(np.unique(donor.groups))} "
              f"toy-arm accuracy={acc:.4f}")
        print(pd.DataFrame(_summary_rows(res)).to_string(index=False))
    for L, spec in dec["labels"].items():
        r = res[L]
        assert r["verified"], f"{L} not verified: {r['heldout'].get('reason')!r} / {r['spatial']['reason']!r}"
        assert r["spatial"]["judgeable"] and r["spatial"]["clear"]
        assert r["spatial"]["I"] > r["spatial"]["null_p"] > 0
        if spec["held_out_gated"]:
            ho = r["heldout"]
            assert ho["verified"] is True
            assert ho["n_judged"] == len(spec["markers"]), "every declared marker was judgeable"
            assert ho["n_clear"] >= ho["n_required"] >= 2
            for m in ho["per_marker"]:
                assert m["judgeable"], f"{L}/{m['marker']}: {m['reason']}"
                assert m["n_tail"] >= CRITERIA["min_tail_in_pool"] and m["n_calls"] >= CRITERIA["min_calls"]
        else:
            assert r["heldout"]["verified"] is None       # Endocrine: spatial + stability only
        st = r["stability"]
        assert len(st["per_replicate"]) == CRITERIA["n_replicates"]
        assert st["invariant"] is True
        assert all(x["judgeable"] and x["spatial_clear"] for x in st["per_replicate"])
        assert r["contamination"]["judgeable"] and "clear" not in r["contamination"], "descriptive only"
    assert "islet_core_audit" in res["Endocrine"] and res["Endocrine"]["islet_core_audit"]["audit_only"]


def test_the_oracle_arm_is_verified_too(dec, planted_donor):
    """The planted labels handed straight to the gate (an arm that is right by construction)
    are verified as well -- so a verification is not an artefact of the toy arm's z-scores."""
    donor = planted_donor(n=8_000, seed=11)
    _, res = judge_all(dec, donor, lambda X: donor.labels, n_replicates=1)
    assert all(res[L]["verified"] for L in dec["labels"]), {L: _verdict_tuple(res[L]) for L in res}


def test_islet_core_audit_is_audit_only(dec, planted_donor):
    """The QuPath parentage column changes a reported number and NEVER a verdict."""
    donor = planted_donor(n=6_000, seed=5)
    arm = lambda X: signature_arm(X, dec)  # noqa: E731
    _, with_cr = judge_all(dec, donor, arm, n_replicates=1, with_cell_region=True)
    _, without = judge_all(dec, donor, arm, n_replicates=1, with_cell_region=False)
    for L in dec["labels"]:
        assert _verdict_tuple(with_cr[L]) == _verdict_tuple(without[L])
    assert "islet_core_audit" in with_cr["Endocrine"] and "islet_core_audit" not in without["Endocrine"]
    assert "islet_core_audit" not in with_cr["T_NK"]
    audit = with_cr["Endocrine"]["islet_core_audit"]
    assert audit["audit_only"] is True and np.isfinite(audit["enrichment"])


# --------------------------------------------------------------------------- (b) scrambled prevalence
def test_verdicts_are_invariant_to_scrambled_label_shares(dec, planted_donor, capsys):
    """Ten rebuilds of the positive control with Dirichlet-scrambled shares: identical verdicts."""
    rng = np.random.default_rng(2026)
    labels = list(dec["labels"])
    verdicts, achieved = [], []
    for i in range(10):
        shares = _draw_shares(rng, labels, floor=0.03, cap=0.5)
        donor = planted_donor(n=8_000, seed=100 + i, shares=shares)
        _, res = judge_all(dec, donor, lambda X: signature_arm(X, dec), n_replicates=1)
        verdicts.append({L: _verdict_tuple(res[L]) for L in labels})
        achieved.append({L: float((donor.labels == L).mean()) for L in labels})
    ach = pd.DataFrame(achieved)
    with capsys.disabled():
        print("\n[scrambled prevalence] achieved shares over 10 rebuilds (min / max per label):")
        print(pd.DataFrame({"min": ach.min(), "max": ach.max()}).T.round(3).to_string())
        print(f"verdict identical across rebuilds: {all(v == verdicts[0] for v in verdicts)}")
    # the shares genuinely moved -- otherwise invariance is vacuous
    assert (ach.max() - ach.min()).max() > 0.2, "the scramble did not move the shares"
    assert ach.min().min() >= 0.02 and ach.max().max() <= 0.55
    # every rebuild, every label, the same four-part verdict ...
    for i, v in enumerate(verdicts):
        assert v == verdicts[0], f"rebuild {i} changed a verdict: {v} != {verdicts[0]}"
    # ... and that verdict is VERIFIED, so the invariance is not seven unjudgeables agreeing
    assert all(verdicts[0][L][3] for L in labels), verdicts[0]


# --------------------------------------------------------------------------- (c) no structure, no tails
def test_a_label_with_no_spatial_structure_and_no_tail_enrichment_is_not_verified(dec, planted_donor):
    """Labels drawn cell by cell (no patches) on pure-noise expression (no tails), emitted by an
    'arm' that simply repeats them: enough calls to judge, and NOT verified.  Note the spatial
    bar's nominal false-positive rate is exactly the null's 5 % per label; the seeds pin one
    draw of that, and the held-out bar has to fail independently on every gated label."""
    donor = planted_donor(n=6_000, seed=3, spatial=False, signal=0.0)
    _, res = judge_all(dec, donor, lambda X: donor.labels, n_replicates=1)
    for L, spec in dec["labels"].items():
        r = res[L]
        assert r["spatial"]["judgeable"], "this must be a FAIL on evidence, not an unjudgeable"
        assert not r["verified"], f"{L} verified on noise: {_verdict_tuple(r)}"
        if spec["held_out_gated"]:
            ho = r["heldout"]
            assert ho["n_judged"] == len(spec["markers"])
            assert ho["verified"] is False
            for m in ho["per_marker"]:
                assert m["judgeable"] and np.isfinite(m["h"])
    # the descriptive statistic still reports (it never blocks)
    assert all(np.isfinite(res[L]["contamination"]["share"]) for L in dec["labels"])


def test_no_structure_means_moran_i_sits_in_the_null_band(dec, planted_donor):
    donor = planted_donor(n=6_000, seed=8, spatial=False, signal=0.0)
    r = gate.spatial_coherence(donor.labels, "Myeloid", donor.xy, donor.groups)
    assert r["judgeable"] and abs(r["I"]) < 5 * r["null_p"]


# --------------------------------------------------------------------------- (d) unjudgeable
def test_too_few_calls_is_unjudgeable_and_never_clears(dec, planted_donor):
    donor = planted_donor(n=4_000, seed=4)
    labels = donor.labels.copy()
    keep = np.flatnonzero(labels == "Stromal")[:30]
    labels[(labels == "Stromal")] = "Vascular"
    labels[keep] = "Stromal"                       # 30 Stromal calls: below min_calls = 50
    sp = gate.spatial_coherence(labels, "Stromal", donor.xy, donor.groups)
    assert sp == {**sp, "judgeable": False, "clear": False, "n_calls": 30}
    assert "only 30 usable observations, need >= 50" in sp["reason"]
    assert np.isnan(sp["I"]) and np.isnan(sp["null_p"])
    ho = gate.heldout_enrichment(labels, "Stromal", gate.tail_mask(donor.X["SMA"].to_numpy()), donor.groups)
    assert ho["judgeable"] is False and ho["clear"] is False and ho["n_calls"] == 30
    assert np.isnan(ho["h"])
    spec = dec["labels"]["Stromal"]
    r = gate.judge_label("Stromal", labels, {m: labels for m in spec["markers"]}, {"subsample:1": labels},
                         values=donor.X, xy=donor.xy, groups=donor.groups, ids=donor.ids,
                         declared_markers=spec["markers"], held_out_gated=True, other_tails={},
                         cell_region=None)
    assert r["verified"] is False
    assert r["heldout"]["n_judged"] == 0 and r["heldout"]["verified"] is False
    assert r["heldout"]["reason"] == "fewer than 2 judgeable held-out markers"
    assert r["spatial"]["judgeable"] is False and r["spatial"]["clear"] is False
    assert r["stability"]["invariant"] is True     # False == False: invariance alone never verifies


def test_too_few_tail_cells_is_unjudgeable_and_never_clears(dec, planted_donor):
    donor = planted_donor(n=4_000, seed=6)
    small = np.zeros(donor.n, dtype=bool)
    small[np.flatnonzero(donor.labels == "T_NK")[:10]] = True          # 10 < min_tail_in_pool
    r = gate.heldout_enrichment(donor.labels, "T_NK", small, donor.groups)
    assert r["judgeable"] is False and r["clear"] is False and r["n_tail"] == 10
    assert "tail cells in pool: only 10" in r["reason"]
    # tail cells the arm marked Unavailable do not count as evidence
    tail = np.zeros(donor.n, dtype=bool)
    idx = np.flatnonzero(donor.labels == "T_NK")[:25]
    tail[idx] = True
    labels = donor.labels.copy()
    labels[idx[:10]] = UNAVAILABLE
    r = gate.heldout_enrichment(labels, "T_NK", tail, donor.groups)
    assert r["judgeable"] is False and r["clear"] is False and r["n_tail"] == 15
    # the same 25 all judged is judgeable (and, being the label's own cells, enriched)
    r = gate.heldout_enrichment(donor.labels, "T_NK", tail, donor.groups)
    assert r["judgeable"] is True and r["n_tail"] == 25 and r["clear"] is True


def test_an_absent_label_is_unjudgeable(dec, planted_donor):
    donor = planted_donor(n=2_000, seed=9)
    labels = np.where(donor.labels == "B_Plasma", "Myeloid", donor.labels).astype(str)
    r = gate.heldout_enrichment(labels, "B_Plasma", gate.tail_mask(donor.X["CD20"].to_numpy()), donor.groups)
    assert r["judgeable"] is False and r["clear"] is False and r["n_calls"] == 0
    s = gate.spatial_coherence(labels, "B_Plasma", donor.xy, donor.groups)
    assert s["judgeable"] is False and s["clear"] is False and s["n_calls"] == 0


def test_a_missing_ablation_or_unacquired_marker_never_clears(dec, planted_donor):
    donor = planted_donor(n=4_000, seed=7)
    spec = dec["labels"]["Vascular"]                          # CD31, Caveolin, Podoplanin
    arm = lambda X: signature_arm(X, dec)  # noqa: E731
    primary = arm(donor.X)
    # only one ablation run, and one marker not in the values table
    ablations = {"CD31": arm(donor.X.drop(columns=["CD31"]))}
    values = donor.X.drop(columns=["Podoplanin"])
    r = gate.judge_label("Vascular", primary, ablations, {}, values=values, xy=donor.xy, groups=donor.groups,
                         ids=donor.ids, declared_markers=spec["markers"], held_out_gated=True,
                         other_tails={}, cell_region=None)
    per = {m["marker"]: m for m in r["heldout"]["per_marker"]}
    assert per["Caveolin"]["reason"] == "ablation_not_run" and per["Caveolin"]["clear"] is False
    assert per["Podoplanin"]["reason"] in ("ablation_not_run", "marker_not_acquired")
    assert per["Podoplanin"]["judgeable"] is False and per["Podoplanin"]["clear"] is False
    assert per["CD31"]["judgeable"] is True
    assert r["heldout"]["n_judged"] == 1 and r["heldout"]["verified"] is False
    assert r["heldout"]["reason"] == "fewer than 2 judgeable held-out markers"
    assert r["spatial"]["clear"] is True, "spatial alone clears here ..."
    assert r["verified"] is False, "... and is not enough for a held-out-gated label"
    assert r["stability"]["invariant"] is None and r["stability"]["per_replicate"] == []


def test_ablation_of_a_marker_the_arm_needs_is_unavailable_and_unjudgeable(dec, planted_donor):
    """An arm that cannot run without the withheld marker writes Unavailable for all cells; the
    held-out statistic then has no judged cells and must not clear."""
    donor = planted_donor(n=4_000, seed=10)
    all_unavailable = np.full(donor.n, UNAVAILABLE, dtype=object).astype(str)
    r = gate.heldout_enrichment(all_unavailable, "T_NK", gate.tail_mask(donor.X["CD3e"].to_numpy()), donor.groups)
    assert r["judgeable"] is False and r["clear"] is False and r["n_tail"] == 0 and r["n_calls"] == 0


# --------------------------------------------------------------------------- (e) null-gated cluster annotation
def _declared_columns(dec) -> dict:
    markers = list(dec["identity_markers"])
    return {L: [markers.index(m) for m in spec["markers"]] for L, spec in dec["labels"].items()}


def test_null_gated_cluster_annotation_recovers_planted_blocks(dec):
    rng = np.random.default_rng(0)
    labels = list(dec["labels"])
    declared = _declared_columns(dec)
    per = 300
    X = rng.normal(size=(per * len(labels) + 20, len(dec["identity_markers"])))
    clusters = np.concatenate([np.repeat(np.arange(len(labels)), per), np.full(20, 99)])
    for i, L in enumerate(labels):
        X[np.ix_(clusters == i, declared[L])] += 3.0
    Z = gate.standardise_within_donor(X)
    df = gate.null_gated_cluster_annotation(Z, clusters, declared).set_index("cluster")
    for i, L in enumerate(labels):
        row = df.loc[i]
        assert row["label"] == L, f"block {i} planted {L}, awarded {row['label']!r} ({row['reason']!r})"
        assert row["reason"] == "" and row["n_above_null"] >= 1
        assert row[f"contrast_{L}"] > row[f"p95_{L}"] > 0
    small = df.loc[99]
    assert small["label"] == UNAVAILABLE and small["reason"] == "cluster_below_min_cells" and small["n"] == 20


def test_a_cluster_mixing_two_lineages_is_ambiguous_not_awarded(dec):
    rng = np.random.default_rng(1)
    declared = _declared_columns(dec)
    # sizes: a small pure T block, a 300-cell half/half mixture, and a large structureless remainder.
    # The mixture is only detectable by an in-vs-out contrast when neither lineage is half the population
    # (with 600 of 1,200 cells T-high the mixed cluster sits at the population mean for T_NK markers).
    n = 1200
    X = rng.normal(size=(n, len(dec["identity_markers"])))
    clusters = np.concatenate([np.zeros(100, int), np.ones(300, int), np.full(800, 2)])
    X[np.ix_(clusters == 0, declared["T_NK"])] += 3.0          # a clean block
    mixed = np.flatnonzero(clusters == 1)
    X[np.ix_(mixed[: len(mixed) // 2], declared["T_NK"])] += 3.0
    X[np.ix_(mixed[len(mixed) // 2:], declared["Myeloid"])] += 3.0
    df = gate.null_gated_cluster_annotation(gate.standardise_within_donor(X), clusters, declared).set_index("cluster")
    assert df.loc[0, "label"] == "T_NK"
    assert df.loc[1, "label"] == AMBIGUOUS and df.loc[1, "reason"] == "multi_lineage"
    # the mixture is caught either by two labels above their nulls or by the co-lineage rule (the second
    # lineage's null is contaminated by the first's markers, so it may sit below its null yet be enriched)
    assert df.loc[1, "n_above_null"] >= 2 or len(df.loc[1, "co_lineages"]) >= 1
    assert df.loc[2, "label"] == UNAVAILABLE and df.loc[2, "reason"] == "no_lineage_beats_null"


def test_noise_false_award_rate_is_measured_and_printed(dec, capsys):
    """200 random clusterings of pure noise: the share of clusters awarded ANY label.  Recorded
    beside every table that uses the annotation; asserted below 0.35, never corrected."""
    rng = np.random.default_rng(42)
    declared = _declared_columns(dec)
    n, k = 2_000, 4
    Z = gate.standardise_within_donor(rng.normal(size=(n, len(dec["identity_markers"]))))
    awarded = ambiguous = total = 0
    for i in range(200):
        clusters = rng.integers(0, k, n)
        df = gate.null_gated_cluster_annotation(Z, clusters, declared, seed=i)
        total += len(df)
        awarded += int((~df["label"].isin([AMBIGUOUS, UNAVAILABLE])).sum())
        ambiguous += int((df["label"] == AMBIGUOUS).sum())
    rate = awarded / total
    with capsys.disabled():
        print(f"\n[null_gated_cluster_annotation] noise false-award rate = {rate:.3f} "
              f"({awarded}/{total} clusters awarded a label over 200 random clusterings of "
              f"{n} pure-noise cells into {k} clusters; ambiguous rate {ambiguous / total:.3f}) "
              f"-- recorded, not corrected")
    assert rate < 0.35, f"noise false-award rate {rate:.3f} is not below 0.35"


def test_standardise_within_donor_handles_nan_and_constant_columns():
    X = np.array([[1.0, 5.0, np.nan], [2.0, 5.0, 1.0], [3.0, 5.0, 2.0]])
    Z = gate.standardise_within_donor(X)
    assert np.isfinite(Z).all()
    assert np.allclose(Z[:, 1], 0.0)                     # constant column: sd -> 1, centred to 0
    assert Z[0, 2] == 0.0                                # NaN -> 0 after standardising
    assert np.isclose(Z[:, 0].mean(), 0) and np.isclose(Z[:, 0].std(), 1)


# --------------------------------------------------------------------------- (f) tail_mask
@pytest.mark.parametrize("n,expected_k", [(5_000, 20), (19_999, 20), (20_000, 20), (50_000, 50), (123_456, 123)])
def test_tail_mask_k_is_max_20_or_a_thousandth(n, expected_k):
    rng = np.random.default_rng(n)
    v = rng.normal(size=n)
    t = gate.tail_mask(v)
    assert t.dtype == bool and t.shape == (n,)
    assert int(t.sum()) == expected_k == max(CRITERIA["tail_min_k"], int(CRITERIA["tail_q"] * n))
    assert v[t].min() > v[~t].max(), "the tail is exactly the top-k values"


def test_tail_mask_never_contains_nan():
    rng = np.random.default_rng(1)
    v = rng.normal(size=3_000)
    v[rng.choice(3_000, 500, replace=False)] = np.nan
    v[:5] = 1e6                                          # the loudest values ...
    v[5:10] = np.nan                                     # ... next to NaNs
    t = gate.tail_mask(v)
    assert int(t.sum()) == 20
    assert not np.isnan(v[t]).any()
    assert t[:5].all()
    fin = np.flatnonzero(np.isfinite(v))
    top = fin[np.argsort(-v[fin])[:20]]
    assert set(np.flatnonzero(t)) == set(top)


def test_tail_mask_with_fewer_finite_than_k_takes_every_finite_value():
    v = np.array([np.nan, 3.0, np.nan, 1.0, 2.0])
    t = gate.tail_mask(v)
    assert t.tolist() == [False, True, False, True, True]
    assert gate.tail_mask(np.full(7, np.nan)).sum() == 0
    assert gate.tail_mask(np.zeros(0)).shape == (0,)


# --------------------------------------------------------------------------- (g) hash_subsample
def test_hash_subsample_is_deterministic_and_order_independent():
    rng = np.random.default_rng(0)
    ids = np.array([f"cell-{i:05d}" for i in rng.permutation(5_000)], dtype=object)
    a = gate.hash_subsample(ids, 1_234)
    b = gate.hash_subsample(ids, 1_234)
    assert np.array_equal(a, b)
    assert len(a) == 1_234 and np.all(np.diff(a) > 0), "sorted, unique indices"
    perm = rng.permutation(len(ids))
    c = gate.hash_subsample(ids[perm], 1_234)
    assert set(ids[a]) == set(ids[perm][c]), "the SAME cells, whatever the row order"
    assert set(ids[gate.hash_subsample(ids, 1_234, salt="other")]) != set(ids[a]), "a different salt draws a different set"
    assert np.array_equal(gate.hash_subsample(ids, len(ids)), np.arange(len(ids)))
    assert np.array_equal(gate.hash_subsample(ids[:10], 50), np.arange(10))


def test_hash_subsample_keys_on_content_not_position():
    ids = np.array([f"c{i}" for i in range(2_000)], dtype=object)
    pick = gate.hash_subsample(ids, 400)
    chosen = set(ids[pick])
    # append strangers: every previously chosen cell is still chosen or displaced only by a
    # stranger with a smaller key -- never by a re-numbering of the originals
    more = np.concatenate([ids, np.array([f"x{i}" for i in range(500)], dtype=object)])
    pick2 = set(more[gate.hash_subsample(more, 400)])
    assert (pick2 & chosen) == (pick2 & set(ids))
    assert len(chosen - pick2) == len(pick2 - chosen)
