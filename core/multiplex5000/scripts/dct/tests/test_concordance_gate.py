"""The concordance gate must reject a wrong model AND accept a right one.

A gate that only ever fails is not a gate.  The load-bearing test here is
``test_a_correct_model_passes_every_blocking_criterion``: a synthetic model that really
does track the measurements has to clear every bar, or the bars are measuring the gate
rather than the model.  Everything else checks that absent or under-powered evidence
fails closed, and that the anchor set is derived rather than typed.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from concordance_gate import (  # noqa: E402
    CRITERIA, GateInputs, derive_anchors, load_anchors, run_gate)
from label_space import load_label_space  # noqa: E402
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

CROSSWALK = json.loads((HERE / "broad_class_crosswalk.json").read_text())
CLASSES = tuple(CROSSWALK["entries"])
ANCHORS = load_anchors()
#: The cohort manifest is gitignored; the two tests that read it (one through
#: ``derive_anchors``) carry ``pytest.mark.data`` and skip, naming it, when it is absent.
MANIFEST = HERE.parents[1] / "data" / "qupath_manifest.json"
BLOCKING = {"sign", "specificity", "detectability", "argmax_coherence",
            "impossible_mass", "null_separation"}


# ----------------------------------------------------------------- the anchors
@pytest.mark.data
def test_anchors_are_derived_not_typed():
    """anchors.json must be exactly what the rule produces, or it has been hand-edited."""
    if not MANIFEST.exists():
        pytest.skip(f"{MANIFEST} not present (data is gitignored)")
    frozen = json.loads((HERE / "anchors.json").read_text())
    fresh = derive_anchors()
    assert fresh == frozen
    assert frozen["method_version"] == "dct-anchors-2026-08-24"


def test_every_anchor_is_a_registry_type_marker_with_a_compartment():
    reg = json.loads((HERE.parents[1] / "multiplex5000" / "marker_registry.json").read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    by = {m["name"]: m for m in reg}
    ls = load_label_space(CLASSES)
    compartments = set(ls.broad_groups(restricted=False))
    for marker, comp in ANCHORS:
        assert by[marker]["kind"] == "type", marker
        assert by[marker]["broad_compartment"] == comp, marker
        assert comp in compartments, comp


def test_anchor_set_size_is_pinned():
    """The 22-donor intersection loses Ker8_18: neither new image has that channel.

    A silent change here moves every count in the pre-registered bar. The rule and
    fractions stay fixed; its measured panel intersection determines the count.
    """
    assert len(ANCHORS) == 27


@pytest.mark.data
def test_every_anchor_is_acquired_on_every_cohort_image():
    if not MANIFEST.exists():
        pytest.skip(f"{MANIFEST} not present (data is gitignored)")
    manifest = json.loads(MANIFEST.read_text())
    missing = {
        im["donor_id"]: sorted({m for m, _ in ANCHORS}
                               - {c["marker"] for c in im["channel_map"]})
        for im in manifest["images"]
    }
    assert not {d: markers for d, markers in missing.items() if markers}


@pytest.mark.parametrize("donor", sorted(DONOR_EXCLUSIONS))
def test_anchor_derivation_rejects_centrally_excluded_donors(tmp_path, donor):
    (tmp_path / "multiplex5000").mkdir()
    (tmp_path / "multiplex5000" / "marker_registry.json").write_text('{"markers": []}')
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "qupath_manifest.json").write_text(json.dumps({
        "images": [{"donor_id": "0" + donor, "channel_map": []}],
    }))
    with pytest.raises(ValueError, match="excluded donor"):
        derive_anchors(repo_root=tmp_path)


# ------------------------------------------------------ the positive control
def _synthetic(n=4000, seed=0, informative=True):
    """A model whose probabilities really do follow the measurements.

    Each cell gets a true compartment; its markers are bright for that compartment and
    dim otherwise; the model puts its mass on a class in the true compartment.  Nothing
    here is tuned -- if the gate rejects this, the gate is wrong.
    """
    rng = np.random.default_rng(seed)
    ls = load_label_space(CLASSES)
    groups = ls.broad_groups(restricted=True)
    comps = sorted({c for _, c in ANCHORS})
    truth = rng.choice(comps, size=n)

    P = rng.random((n, len(CLASSES))) * 0.01
    for i, comp in enumerate(truth):
        pool = groups.get(comp) or groups[comps[0]]
        pick = CLASSES.index(pool[rng.integers(len(pool))])
        P[i, pick] += 20.0 if informative else 0.0
    if not informative:                      # same marginal shape, no cell-level signal
        P += rng.random((n, len(CLASSES))) * 0.01
    P /= P.sum(1, keepdims=True)

    measured = {}
    for marker, comp in ANCHORS:
        base = rng.lognormal(2.0, 0.6, n)
        measured[marker] = base * np.where(truth == comp, 6.0, 1.0)
    return GateInputs(P, CLASSES, measured, label="synthetic"), truth


def test_a_correct_model_passes_every_blocking_criterion():
    inp, _ = _synthetic()
    r = run_gate(inp, restricted=True)
    failed = {k: v[1] for k, v in r.criteria.items() if k in BLOCKING and not v[0]}
    assert not failed, f"the gate rejects a model that tracks the measurements: {failed}"
    assert r.passed


def test_an_uninformative_model_is_rejected():
    inp, _ = _synthetic(informative=False)
    r = run_gate(inp, restricted=True)
    assert not r.passed
    assert not r.criteria["null_separation"][0] or not r.criteria["detectability"][0]


def test_permutation_null_is_centred_on_no_information():
    inp, _ = _synthetic()
    r = run_gate(inp, restricted=True)
    assert r.null["draws"] > 0
    assert r.null["rho_abs"] < 0.10, r.null
    assert r.null["auroc_dev"] < 0.10, r.null


# ------------------------------------------------------------- fail closed
def test_absent_measurements_fail_rather_than_pass():
    """Every anchor unmeasured must be a FAIL. `nan > bar` is False, which is the trap."""
    inp, _ = _synthetic()
    blank = GateInputs(inp.probabilities, CLASSES,
                       {m: np.full(len(inp.probabilities), np.nan) for m, _ in ANCHORS})
    r = run_gate(blank, restricted=True)
    assert not r.passed
    assert not r.criteria["sign"][0] and not r.criteria["detectability"][0]
    assert not r.anchor_verdict.passed


def test_under_powered_evidence_fails_rather_than_passes():
    inp, _ = _synthetic(n=int(CRITERIA["min_obs"]) - 5)
    r = run_gate(inp, restricted=True)
    assert not r.passed


def test_confidence_criterion_never_blocks_release():
    """It licenses stratification; it must not be able to fail a lane on its own."""
    inp, _ = _synthetic()
    r = run_gate(inp, restricted=True)
    r.criteria["confidence_informativeness"] = (False, "forced")
    assert r.passed, "a confidence failure must not block release"


def test_no_criterion_consults_an_expected_frequency():
    """Pinned after the first draft of criterion 4 did exactly that.

    An earlier `no_collapse` capped each compartment's share of calls at 45 % and floored
    it at 0.5 %.  Roughly 85 % of exocrine pancreas is acinar, so a correct model would
    have failed the cap, and the floor was an expectation about how many nerve cells a
    pancreas holds.  Both are the tuning CLAUDE.md's first non-negotiable forbids.
    """
    for key in ("collapse_max_share", "collapse_min_share"):
        assert key not in CRITERIA, f"{key} encodes an expected cell-type frequency"
    inp, _ = _synthetic()
    r = run_gate(inp, restricted=True)
    assert "argmax_coherence" in r.criteria and "no_collapse" not in r.criteria


def test_restricted_and_full_label_spaces_both_run():
    inp, _ = _synthetic()
    a = run_gate(inp, restricted=True)
    b = run_gate(inp, restricted=False)
    assert set(a.criteria) == set(b.criteria) == BLOCKING | {"confidence_informativeness"}
