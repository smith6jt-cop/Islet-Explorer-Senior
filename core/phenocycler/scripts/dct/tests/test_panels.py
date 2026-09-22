"""The panel arms must differ by exactly what they claim to differ by.

The three-point comparison only answers the question if the arms are clean:

    all (52)  ->  type + exclusive (43)  ->  type only (30)

`all` vs `type_plus_exclusive` must differ by exactly the markers declared present-on-
everything or orthogonal-to-identity, and `type_plus_exclusive` vs `type_only` by exactly
the markers declared lineage-informative. If either boundary leaks, the arm measures two
things at once and cannot say whether functional markers help.

Every test that resolves a donor's channels reads the cohort manifest
(``core/phenocycler/data/qupath_manifest.json``, gitignored) and carries ``pytest.mark.data``;
the manifest is opened inside a module-scoped fixture so the file collects the same ids
whether or not the manifest is present.  ``test_exclusivity_artifact_separates_its_two_reasons``
reads only the tracked ``marker_exclusivity.json`` and stays in the synthetic lane.
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

import arms as A  # noqa: E402


@pytest.fixture(scope="module")
def manifest() -> dict:
    """The cohort manifest, or a skip that names it -- never a module-level skip.

    ``pytest.skip(..., allow_module_level=True)`` removed these ids from collection when the
    manifest was absent, so a run's collected count depended on the data tree and could not
    be reconciled against a deselection.
    """
    if not MANIFEST.exists():
        pytest.skip(f"cohort manifest not present in this checkout ({MANIFEST})")
    return json.loads(MANIFEST.read_text())


def _kind():
    reg = json.loads(REGISTRY.read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    return {m["name"]: m["kind"] for m in reg}


def _panel(man, donor, arm):
    im = next(i for i in man["images"] if i["donor_id"] == donor)
    return {m for _, m, _ in A._channels(im, arm, _kind())}


def _by_name(name):
    return next(a for a in A.ARMS + A.ROUND2_ARMS if a.name == name)


DONORS = ["6539", "6523"]          # a 59-channel and a 58-channel panel variant


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_type_only_is_exactly_the_registry_type_markers(manifest, donor):
    kind = _kind()
    got = _panel(manifest, donor, _by_name("A1-type-only"))
    assert got and all(kind[m] == "type" for m in got)


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_exclusive_panel_adds_only_declared_informative_markers(manifest, donor):
    kind = _kind()
    base = _panel(manifest, donor, _by_name("A1-type-only"))
    got = _panel(manifest, donor, _by_name("A1d-type-plus-exclusive"))
    added = got - base
    assert added, "the exclusivity arm added nothing -- it would duplicate type_only"
    for m in added:
        assert kind[m] != "type"
        assert EXCL[m]["lineage_informative"] is True, m
        assert EXCL[m]["in_dct_vocabulary"] is True, m
        assert kind[m] != "excluded", f"{m} is registry-excluded and must not enter a panel arm"


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_full_panel_drops_only_declared_uninformative_markers(manifest, donor):
    """all -> type_plus_exclusive must remove nothing that was declared informative."""
    full = _panel(manifest, donor, _by_name("A0-baseline"))
    got = _panel(manifest, donor, _by_name("A1d-type-plus-exclusive"))
    kind = _kind()
    for m in full - got:
        assert kind[m] != "type", f"{m} is a lineage marker and must not be dropped"
        informative = EXCL.get(m, {}).get("lineage_informative")
        excluded = kind[m] == "excluded"
        assert not informative or excluded, f"{m} was declared informative but got dropped"


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_dapi_is_admitted_by_its_own_flag_not_by_the_panel(manifest, donor):
    """Pinned: the panel filter used to override with_dapi, so the arm measured nothing."""
    for base, plus in [("A1-type-only", "A1c-type-plus-dna"),
                       ("A1d-type-plus-exclusive", "A1e-type-plus-exclusive-dna")]:
        b, p = _panel(manifest, donor, _by_name(base)), _panel(manifest, donor, _by_name(plus))
        assert p - b == {"DAPI"}, f"{plus} should add DAPI and nothing else, got {p - b}"


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_random30_is_size_matched_to_the_type_panel(manifest, donor):
    n_type = len(_panel(manifest, donor, _by_name("A1-type-only")))
    for a in A.ROUND2_ARMS:
        if a.panel == "random30":
            assert len(_panel(manifest, donor, a)) == n_type


@pytest.mark.data
def test_random30_draws_differ_between_seeds(manifest):
    """A control that drew the same channels every time would control for nothing."""
    draws = [frozenset(_panel(manifest, "6539", a)) for a in A.ROUND2_ARMS if a.panel == "random30"]
    assert len(set(draws)) == len(draws) > 1


@pytest.mark.data
@pytest.mark.parametrize("donor", DONORS)
def test_every_arm_resolves_to_distinct_dct_names(manifest, donor):
    """Two channels resolving to one DCT name means dataset.py silently drops one."""
    im = next(i for i in manifest["images"] if i["donor_id"] == donor)
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
