"""Regression tests for the fail-closed rule.

The first plane-identity run reported PASS on 52 markers whose correlations were
all NaN, because `nan <= nan` is False. These tests exist so that cannot recur.
"""
import sys, pathlib
import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from gates import assess_self_match, require_evidence

L3 = ["CD3e", "CD68", "CD20"]
N3 = np.full(3, 5000)


def _good():
    return np.array([[0.9, 0.3, 0.2],
                     [0.3, 0.9, 0.1],
                     [0.2, 0.1, 0.9]])


def test_clean_diagonal_passes():
    v = assess_self_match(L3, _good(), N3)
    assert v.passed and bool(v) is True


def test_all_nan_fails():                       # the bug that shipped
    v = assess_self_match(L3, np.full((3, 3), np.nan), N3)
    assert not v.passed
    assert "not finite" in v.reason


def test_one_nan_fails_only_that_marker():
    S = _good(); S[1, 1] = np.nan
    v = assess_self_match(L3, S, N3)
    assert not v.passed
    assert [r[0] for r in v.rows if not r[5]] == ["CD68"]


def test_swapped_planes_fail():
    S = _good()[[1, 0, 2]]                      # CD3e and CD68 planes transposed
    assert not assess_self_match(L3, S, N3).passed


def test_empty_input_fails():
    assert not assess_self_match([], np.zeros((0, 0)), np.zeros(0)).passed


def test_too_few_observations_fails():
    v = assess_self_match(L3, _good(), np.array([5000, 10, 5000]), min_obs=50)
    assert not v.passed and "usable observations" in v.reason


def test_tie_is_not_a_pass():
    S = _good(); S[0, 1] = 0.9                  # rival ties the diagonal
    assert not assess_self_match(L3, S, N3).passed


def test_margin_is_enforced():
    S = _good(); S[0, 1] = 0.89
    assert assess_self_match(L3, S, N3, margin=0.0).passed
    assert not assess_self_match(L3, S, N3, margin=0.05).passed


def test_self_score_floor():
    S = np.array([[0.2, 0.01, 0.0], [0.01, 0.9, 0.0], [0.0, 0.0, 0.9]])
    assert assess_self_match(L3, S, N3).passed          # wins, but weakly
    assert not assess_self_match(L3, S, N3, min_self=0.5).passed


def test_shape_mismatch_fails():
    assert not assess_self_match(L3, np.zeros((2, 2)), np.zeros(2)).passed


def test_infinity_is_not_finite():
    S = _good(); S[2, 2] = np.inf
    assert not assess_self_match(L3, S, N3).passed


@pytest.mark.parametrize("n,lo,expect_fail", [(49, 50, True), (50, 50, False), (0, 1, True)])
def test_require_evidence(n, lo, expect_fail):
    assert (require_evidence(n, lo, "x") is not None) is expect_fail


# --- the reader self-test must also refuse to pass on no evidence -----------

def test_reader_selftest_refuses_vacuous_cases():
    """A window entirely outside the page compares empty-to-empty; that is not a pass."""
    import numpy as _np
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
    from tiff_window import selftest_bit_exact

    class _Page:                      # 100x100 single-chunk page, chunk stored sparse
        shape, dtype, is_tiled = (100, 100), _np.dtype("uint16"), True
        tilelength = tilewidth = 100
        chunked = (1, 1)
        dataoffsets, databytecounts = (0,), (0,)      # 0 bytes -> read_window skips, no I/O
        parent = None
        def asarray(self): return _np.zeros((100, 100), _np.uint16)

    # size larger than the page forces out-of-bounds / empty reference slices
    assert selftest_bit_exact(_Page(), size=400, verbose=False) is False
