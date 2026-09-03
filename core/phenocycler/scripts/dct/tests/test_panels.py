"""The panel arms must differ by exactly what they claim to differ by.

The three-point comparison only answers the question if the arms are clean:

    all (52)  ->  type + exclusive (43)  ->  type only (30)

`all` vs `type_plus_exclusive` must differ by exactly the markers declared present-on-
everything or orthogonal-to-identity, and `type_plus_exclusive` vs `type_only` by exactly
the markers declared lineage-informative. If either boundary leaks, the arm measures two
things at once and cannot say whether functional markers help.
"""
from __future__ import annotations

import json
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

MANIFEST = HERE.parents[1] / "data" / "qupath_manifest.json"
REGISTRY = HERE.parents[1] / "phenocycler" / "marker_registry.json"
EXCL = json.loads((HERE / "marker_exclusivity.json").read_text())["entries"]

pytest.importorskip("deepcell_types", reason="the DCT extra is not installed")
if not MANIFEST.exists():
    pytest.skip("cohort manifest not present in this checkout", allow_module_level=True)

import arms as A  # noqa: E402


def _kind():
    reg = json.loads(REGISTRY.read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    return {m["name"]: m["kind"] for m in reg}


def _panel(donor, arm):
    man = json.loads(MANIFEST.read_text())
    im = next(i for i in man["images"] if i["donor_id"] == donor)
    return {m for _, m, _ in A._channels(im, arm, _kind())}


def _by_name(name):
    return next(a for a in A.ARMS + A.ROUND2_ARMS if a.name == name)


DONORS = ["6539", "6523"]          # a 59-channel and a 58-channel panel variant


@pytest.mark.parametrize("donor", DONORS)
def test_type_only_is_exactly_the_registry_type_markers(donor):
    kind = _kind()
    got = _panel(donor, _by_name("A1-type-only"))
    assert got and all(kind[m] == "type" for m in got)


@pytest.mark.parametrize("donor", DONORS)
def test_exclusive_panel_adds_only_declared_informative_markers(donor):
    kind = _kind()
    base = _panel(donor, _by_name("A1-type-only"))
    got = _panel(donor, _by_name("A1d-type-plus-exclusive"))
    added = got - base
    assert added, "the exclusivity arm added nothing -- it would duplicate type_only"
    for m in added:
        assert kind[m] != "type"
        assert EXCL[m]["lineage_informative"] is True, m
        assert EXCL[m]["in_dct_vocabulary"] is True, m
        assert kind[m] != "excluded", f"{m} is registry-excluded and must not enter a panel arm"


@pytest.mark.parametrize("donor", DONORS)
def test_full_panel_drops_only_declared_uninformative_markers(donor):
    """all -> type_plus_exclusive must remove nothing that was declared informative."""
    full = _panel(donor, _by_name("A0-baseline"))
    got = _panel(donor, _by_name("A1d-type-plus-exclusive"))
    kind = _kind()
    for m in full - got:
        assert kind[m] != "type", f"{m} is a lineage marker and must not be dropped"
        informative = EXCL.get(m, {}).get("lineage_informative")
        excluded = kind[m] == "excluded"
        assert not informative or excluded, f"{m} was declared informative but got dropped"


@pytest.mark.parametrize("donor", DONORS)
def test_dapi_is_admitted_by_its_own_flag_not_by_the_panel(donor):
    """Pinned: the panel filter used to override with_dapi, so the arm measured nothing."""
    for base, plus in [("A1-type-only", "A1c-type-plus-dna"),
                       ("A1d-type-plus-exclusive", "A1e-type-plus-exclusive-dna")]:
        b, p = _panel(donor, _by_name(base)), _panel(donor, _by_name(plus))
        assert p - b == {"DAPI"}, f"{plus} should add DAPI and nothing else, got {p - b}"


@pytest.mark.parametrize("donor", DONORS)
def test_random30_is_size_matched_to_the_type_panel(donor):
    n_type = len(_panel(donor, _by_name("A1-type-only")))
    for a in A.ROUND2_ARMS:
        if a.panel == "random30":
            assert len(_panel(donor, a)) == n_type


def test_random30_draws_differ_between_seeds():
    """A control that drew the same channels every time would control for nothing."""
    draws = [frozenset(_panel("6539", a)) for a in A.ROUND2_ARMS if a.panel == "random30"]
    assert len(set(draws)) == len(draws) > 1


@pytest.mark.parametrize("donor", DONORS)
def test_every_arm_resolves_to_distinct_dct_names(donor):
    """Two channels resolving to one DCT name means dataset.py silently drops one."""
    man = json.loads(MANIFEST.read_text())
    im = next(i for i in man["images"] if i["donor_id"] == donor)
    for a in A.ARMS + A.ROUND2_ARMS:
        names = [d for _, _, d in A._channels(im, a, _kind())]
        assert len(names) == len(set(names)), f"{a.name} on {donor}: duplicate {names}"


def test_exclusivity_artifact_separates_its_two_reasons():
    """Pinned after the first draft collapsed them and mislabelled Caveolin."""
    for m, v in EXCL.items():
        assert isinstance(v["lineage_informative"], bool), m
        assert isinstance(v["in_dct_vocabulary"], bool), m
        assert v["reason"].strip(), m
    assert EXCL["Caveolin"]["lineage_informative"] is True
    assert EXCL["Caveolin"]["in_dct_vocabulary"] is False
