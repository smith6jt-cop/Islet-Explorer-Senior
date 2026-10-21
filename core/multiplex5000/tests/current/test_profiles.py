"""Profiles: the schema refuses what it cannot read, and the tissue profiles are the authority.

Until D2 (2026-10-20) each tissue profile mirrored a constant in ``integration/tissues.py``. Since
D2 that module reads the profiles, so these tests hold the other direction: it keeps no copy of
its own, and the values it serves are pinned to the ones it had before D2, so a change to a
tissue is a deliberate edit to its profile and to the pin below, never a drift.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_profiles.py -q
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from multiplex5000 import profiles
from multiplex5000.integration import tissues

REPO = Path(__file__).resolve().parents[4]
SHIPPED = profiles.load_all()


def _tissue(tissue_id: str) -> profiles.TissueProfile:
    return SHIPPED[("tissue", tissue_id)]


def test_every_shipped_profile_loads():
    kinds = {kind for kind, _ in SHIPPED}
    assert {"tissue", "platform", "panel", "cohort"} <= kinds


def test_every_tissue_the_package_knows_has_a_profile_and_no_other():
    assert {i for k, i in SHIPPED if k == "tissue"} == set(tissues.TISSUES)


@pytest.mark.parametrize("tissue_id", tissues.TISSUES)
def test_tissues_py_serves_the_profile_and_keeps_no_copy(tissue_id):
    t = _tissue(tissue_id)
    assert tissues.PROFILES[tissue_id] == t
    assert t.display_name == tissues.DISPLAY_NAME[tissue_id]
    assert all(tissues.normalize(alias) == tissue_id for alias in t.aliases)
    assert {p for p, tid in tissues.ROI_PREFIX_TO_TISSUE.items() if tid == tissue_id} == set(
        t.roi_prefixes)
    assert t.panel_tissue_dir == tissues.PANEL_TISSUE_DIR[tissue_id]
    assert set(t.tissue_only_panels) == tissues.TISSUE_ONLY_PANELS[tissue_id]
    assert t.structures == tissues.STRUCTURE_KINDS[tissue_id]
    assert t.comparable_lineages == tissues.COMPARABLE_LINEAGES[tissue_id]


#: What ``integration/tissues.py`` held as literals before D2, copied from 4f802d4. Serving the
#: profiles changed none of it.
BEFORE_D2 = {
    "TISSUES": ("pancreas", "pancreatic_lymph_node"),
    "DISPLAY_NAME": {"pancreas": "pancreas", "pancreatic_lymph_node": "pancreatic lymph node"},
    "ROI_PREFIX_TO_TISSUE": {"pan": "pancreas", "pln": "pancreatic_lymph_node",
                             "ln": "pancreatic_lymph_node"},
    "PANEL_TISSUE_DIR": {"pancreas": "pancreas", "pancreatic_lymph_node": None},
    "PANCREAS_ONLY_PANELS": frozenset({"01_pancreas_endocrine", "02_pancreas_exocrine",
                                       "13_fibroblast_stellate"}),
    "STRUCTURE_KINDS": {"pancreas": ("islet", "duct", "vessel"), "pancreatic_lymph_node": ()},
    "COMPARABLE_LINEAGES": {
        "pancreas": ("Endocrine", "Exocrine", "Stromal", "Vascular",
                     "T_NK", "B_Plasma", "Myeloid", "Neural"),
        "pancreatic_lymph_node": ("Stromal", "Vascular", "T_NK", "B_Plasma", "Myeloid")},
}


@pytest.mark.parametrize("name", sorted(BEFORE_D2))
def test_the_values_are_the_ones_tissues_py_held_before_d2(name):
    value = getattr(tissues, name)
    assert value == BEFORE_D2[name]
    if isinstance(value, dict):
        assert list(value) == list(BEFORE_D2[name]), "order matters: for_roi tries prefixes in turn"


def test_only_the_pancreas_carries_tissue_only_panels():
    assert set(_tissue(tissues.PANCREAS).tissue_only_panels) == set(tissues.PANCREAS_ONLY_PANELS)
    assert _tissue(tissues.LYMPH_NODE).tissue_only_panels == ()
    assert tissues.applies(tissues.PANCREAS, "01_pancreas_endocrine")
    assert not tissues.applies(tissues.LYMPH_NODE, "01_pancreas_endocrine")


@pytest.mark.parametrize("ini", ["config.ini", "config_metabolic.ini"])
def test_the_configured_tissues_are_the_cohort_profiles(ini):
    """config.ini's ``tissues`` names the cohort's tissues; the cohort profile is their authority."""
    import configparser
    parser = configparser.ConfigParser()
    parser.read(REPO / "core/multiplex5000" / ini)
    configured = [t.strip() for t in parser.get("integration", "tissues").split(",") if t.strip()]
    assert configured == list(SHIPPED[("cohort", "io60_phenocycler")].tissues)


def test_every_spelling_normalize_accepts_is_an_alias_of_its_tissue():
    """The profile's alias list is complete: nothing normalize() accepts is missing from it."""
    accepted = {"pln", "lymph_node", "ln", "pancreatic_ln", "panc_lymph_node", "panc", "pan",
                "pancreatic", *tissues.TISSUES}
    listed = {a for tid in tissues.TISSUES for a in _tissue(tid).aliases}
    assert accepted == listed


def test_the_phenocycler_pixel_size_mirrors_the_app_config():
    config = yaml.safe_load((REPO / "apps/islet/config/run20_default.yml").read_text())
    assert SHIPPED[("platform", "phenocycler")].pixel_size_um == config["pixel_size_um"]


def test_the_post_xenium_if_pixel_size_mirrors_restain_sections():
    import csv

    path = REPO / "core/phenocycler/data/integration/restain_sections.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    declared = SHIPPED[("platform", "post_xenium_if")].pixel_size_um
    assert rows and {float(r["pixel_size_um"]) for r in rows} == {declared}


def test_panels_and_cohorts_reference_their_authority_and_a_known_platform():
    platforms = {i for k, i in SHIPPED if k == "platform"}
    for (kind, _), profile in SHIPPED.items():
        if kind in ("panel", "cohort"):
            assert profile.platform in platforms
    cohort = SHIPPED[("cohort", "io60_phenocycler")]
    assert set(cohort.tissues) <= set(tissues.TISSUES)
    assert cohort.exclusions == "multiplex5000.cohort.DONOR_EXCLUSIONS"


def _write(tmp_path, kind, name, body) -> Path:
    folder = tmp_path / f"{kind}s"
    folder.mkdir(exist_ok=True)
    path = folder / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"schema": profiles.PROFILE_SCHEMA, "kind": kind, **body}))
    return path


def test_a_cohort_profile_that_lists_donors_is_refused(tmp_path):
    path = _write(tmp_path, "cohort", "c", {
        "id": "c", "platform": "phenocycler", "tissues": ["pancreas"],
        "membership": "core/multiplex5000/cohort-spec.json",
        "exclusions": "multiplex5000.cohort.DONOR_EXCLUSIONS", "donors": ["6457"]})
    with pytest.raises(profiles.ProfileError, match="unknown key"):
        profiles.load_profile(path)


@pytest.mark.parametrize("mutation, message", [
    (lambda b: b.pop("structures"), "missing key"),
    (lambda b: b.update(structures="islet"), "list of strings"),
    (lambda b: b.update(id="other"), "file name"),
    (lambda b: b.update(plugin=3), "string or null"),
    (lambda b: b.update(typo_key=1), "unknown key"),
])
def test_the_schema_refuses_a_malformed_tissue_profile(tmp_path, mutation, message):
    body = {"id": "t", "display_name": "t", "aliases": ["t"], "roi_prefixes": [],
            "panel_tissue_dir": None, "tissue_only_panels": [], "structures": [],
            "comparable_lineages": [], "plugin": None}
    mutation(body)
    with pytest.raises(profiles.ProfileError, match=message):
        profiles.load_profile(_write(tmp_path, "tissue", "t", body))


def test_the_schema_refuses_a_wrong_schema_kind_folder_or_reference(tmp_path):
    path = tmp_path / "tissues" / "t.yaml"
    path.parent.mkdir()
    path.write_text(yaml.safe_dump({"schema": "other/1", "kind": "tissue"}))
    with pytest.raises(profiles.ProfileError, match="schema"):
        profiles.load_profile(path)
    path.write_text(yaml.safe_dump({"schema": profiles.PROFILE_SCHEMA, "kind": "organ"}))
    with pytest.raises(profiles.ProfileError, match="kind"):
        profiles.load_profile(path)
    misplaced = _write(tmp_path, "platform", "p", {"id": "p", "display_name": "p",
                                                   "pixel_size_um": 0.5})
    moved = tmp_path / "tissues" / "p.yaml"
    misplaced.rename(moved)
    with pytest.raises(profiles.ProfileError, match="folder"):
        profiles.load_profile(moved)
    panel = _write(tmp_path, "panel", "x", {"id": "x", "platform": "phenocycler",
                                            "authority": "no/such/file.json"})
    with pytest.raises(profiles.ProfileError, match="repository-relative file"):
        profiles.load_profile(panel)
    platform = _write(tmp_path, "platform", "q", {"id": "q", "display_name": "q",
                                                  "pixel_size_um": -1})
    with pytest.raises(profiles.ProfileError, match="positive"):
        profiles.load_profile(platform)


@pytest.mark.parametrize("aliases, message", [
    ({"HDL112": 112}, "quoted donor ids"),          # an unquoted 0112 is the number 112
    ({"0112": "0112"}, "cannot also be an alias"),
    (["HDL112"], "quoted donor ids"),
])
def test_donor_aliases_are_quoted_text_and_never_an_id(tmp_path, aliases, message):
    path = _write(tmp_path, "cohort", "c", {
        "id": "c", "platform": "phenocycler", "tissues": ["pancreas"],
        "membership": "core/multiplex5000/cohort-spec.json",
        "exclusions": "multiplex5000.cohort.DONOR_EXCLUSIONS", "donor_aliases": aliases})
    with pytest.raises(profiles.ProfileError, match=message):
        profiles.load_profile(path)


def test_the_cohort_declares_its_donor_aliases():
    assert SHIPPED[("cohort", "io60_phenocycler")].donor_aliases == {
        "112": "0112", "HDL112": "0112", "HDL115": "115"}
