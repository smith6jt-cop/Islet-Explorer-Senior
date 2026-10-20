"""Profiles: the schema refuses what it cannot read, and every shipped value mirrors the constant
it names, so a profile can never become a second, silently diverging authority before phase D2
makes profiles the authority.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_profiles.py -q
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from phenocycler import profiles
from phenocycler.integration import tissues

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
def test_a_tissue_profile_mirrors_tissues_py(tissue_id):
    t = _tissue(tissue_id)
    assert t.display_name == tissues.DISPLAY_NAME[tissue_id]
    assert all(tissues.normalize(alias) == tissue_id for alias in t.aliases)
    assert {p for p, tid in tissues.ROI_PREFIX_TO_TISSUE.items() if tid == tissue_id} == set(
        t.roi_prefixes)
    assert t.panel_tissue_dir == tissues.PANEL_TISSUE_DIR[tissue_id]
    assert t.structures == tissues.STRUCTURE_KINDS[tissue_id]
    assert t.comparable_lineages == tissues.COMPARABLE_LINEAGES[tissue_id]


def test_only_the_pancreas_carries_tissue_only_panels():
    assert set(_tissue(tissues.PANCREAS).tissue_only_panels) == set(tissues.PANCREAS_ONLY_PANELS)
    assert _tissue(tissues.LYMPH_NODE).tissue_only_panels == ()


def test_every_spelling_normalize_accepts_is_an_alias_of_its_tissue():
    """The profile's alias list is complete: nothing normalize() accepts is missing from it."""
    accepted = {"pln", "lymph_node", "ln", "pancreatic_ln", "panc_lymph_node", "panc", "pan",
                "pancreatic", *tissues.TISSUES}
    listed = {a for tid in tissues.TISSUES for a in _tissue(tid).aliases}
    assert accepted == listed


def test_the_phenocycler_pixel_size_mirrors_the_app_config():
    config = yaml.safe_load((REPO / "apps/islet/config/run20_default.yml").read_text())
    assert SHIPPED[("platform", "phenocycler")].pixel_size_um == config["pixel_size_um"]


def test_panels_and_cohorts_reference_their_authority_and_a_known_platform():
    platforms = {i for k, i in SHIPPED if k == "platform"}
    for (kind, _), profile in SHIPPED.items():
        if kind in ("panel", "cohort"):
            assert profile.platform in platforms
    cohort = SHIPPED[("cohort", "io60_phenocycler")]
    assert set(cohort.tissues) <= set(tissues.TISSUES)
    assert cohort.exclusions == "phenocycler.cohort.DONOR_EXCLUSIONS"


def _write(tmp_path, kind, name, body) -> Path:
    folder = tmp_path / f"{kind}s"
    folder.mkdir(exist_ok=True)
    path = folder / f"{name}.yaml"
    path.write_text(yaml.safe_dump({"schema": profiles.PROFILE_SCHEMA, "kind": kind, **body}))
    return path


def test_a_cohort_profile_that_lists_donors_is_refused(tmp_path):
    path = _write(tmp_path, "cohort", "c", {
        "id": "c", "platform": "phenocycler", "tissues": ["pancreas"],
        "membership": "core/phenocycler/cohort-spec.json",
        "exclusions": "phenocycler.cohort.DONOR_EXCLUSIONS", "donors": ["6457"]})
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
