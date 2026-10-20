"""The bar is frozen: exactly these keys, exactly these values, and no key that names a share.

``criteria.CRITERIA`` was written before any arm ran (protocol Section 3).  A new key, a
changed value, or a key whose name smells of an expected frequency all fail here, so the bar
cannot drift toward a result.  The value pins repeat the protocol text verbatim on purpose:
the document and the code must say the same thing.
"""
from __future__ import annotations

import pytest

import criteria
from criteria import CRITERIA, FORBIDDEN_KEY_TOKENS, assert_frequency_free

#: The protocol's Section 3 list, verbatim.
FROZEN = {
    "tail_q": 0.001,
    "tail_min_k": 20,
    "min_calls": 50,
    "min_tail_in_pool": 20,
    "null_draws": 200,
    "null_quantile": 0.95,
    "morans_k": 12,
    "n_replicates": 5,
    "primary_seed": 0,
    "subsample_frac": 0.9,
    "n_blind": 8,
    "n_islet": 8,
    "core_px": 1024,
    "window_px": 2048,
    "heldout_frac": 2 / 3,
    "survive_frac": 2 / 3,
}


def test_the_key_set_is_exactly_the_frozen_set():
    assert set(CRITERIA) == set(FROZEN)


def test_the_values_are_the_protocol_values():
    for k, v in FROZEN.items():
        assert CRITERIA[k] == v, f"{k}: {CRITERIA[k]!r} != protocol {v!r}"


def test_assert_frequency_free_passes_on_the_bar():
    assert_frequency_free(CRITERIA)
    assert_frequency_free()                     # the import-time call, repeated


@pytest.mark.parametrize("bad_key", ["min_share", "expected_prevalence", "label_frequency",
                                     "mixture_prior", "disease_order", "status_order", "collapse_rule"])
def test_assert_frequency_free_raises_on_a_share_like_key(bad_key):
    with pytest.raises(AssertionError, match="expected-frequency term"):
        assert_frequency_free({**CRITERIA, bad_key: 0.1})


def test_assert_frequency_free_raises_on_min_share_alone():
    """The case named in the task: a key called ``min_share`` is refused outright."""
    with pytest.raises(AssertionError, match="min_share"):
        assert_frequency_free({"min_share": 0.02, **{k: v for k, v in CRITERIA.items()}})


@pytest.mark.parametrize("frac_key", ["heldout_frac", "survive_frac", "subsample_frac"])
@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_fractions_must_be_fractions_of_judged_items(frac_key, bad):
    with pytest.raises(AssertionError, match="fraction of judged items"):
        assert_frequency_free({**CRITERIA, frac_key: bad})


def test_forbidden_tokens_include_every_hidden_form_named_in_the_rules():
    """share / prevalence / frequency / prior / expected / disease / status / order -- all present."""
    for tok in ("share", "prevalence", "frequency", "prior", "expected", "disease", "status", "order"):
        assert tok in FORBIDDEN_KEY_TOKENS


def test_no_criterion_encodes_a_share_of_cells():
    """Evidence floors are counts on a statistic's n; the only fractions are of JUDGED items."""
    for k, v in CRITERIA.items():
        if k.endswith("_frac"):
            assert 0 < v <= 1
        elif k in ("tail_q", "null_quantile"):
            assert 0 < v < 1                    # a quantile, not a prevalence
        else:
            assert isinstance(v, int) and not isinstance(v, bool), f"{k} should be an integer count/size"


def test_module_import_runs_the_check():
    """``assert_frequency_free()`` is called at import, so a bad key never reaches an arm."""
    src = (criteria.__file__ and open(criteria.__file__).read())
    assert src.rstrip().endswith("assert_frequency_free()")
