"""The Xenium readers ask the exclusion for the Xenium platform; the PhenoCycler ones do not.

The maintainer, 2026-10-22: "the xenium for 6457 is not excluded". 6457's recorded reason --
poor DAPI that defeated the PhenoCycler channel registration -- is a PhenoCycler one, so
``cohort.EXCLUDED_PLATFORMS`` scopes it to ``phenocycler`` and every reader of Xenium data asks
``cohort.excluded_on(donor, "xenium")``. Each test below drives one such reader with 6457 and
checks it is admitted there while the PhenoCycler side still refuses it; the fail-closed rule
(an exclusion with no declared platform covers every platform) is checked beside it.

The builders with their own fixtures carry their tests there: availability and sections in
``test_availability_layer.py``, the re-stain overlay in ``test_restain_layer.py``, the grid in
``test_grid_composition.py``, the distance profiles in ``test_distance_profiles.py``, phenospace
in ``tests/integration``, and the two Xenium typing lanes under their own ``scripts/*/tests``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from phenocycler import cohort  # noqa: E402
from phenocycler.cohort import ensure_eligible_donors, excluded_on  # noqa: E402
from multiplex5000.tissues import pancreas_xenium_worker as panc_worker  # noqa: E402


def _load(path: Path, name: str):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


pln_v3 = _load(REPO / "analysis/xenium/scripts/pln_immune_typing_v3.py", "pln_v3_xenium_scope")
cluster = _load(REPO / "core/multiplex5000/scripts/pln/cluster_pln.py", "cluster_pln_xenium_scope")
linkage = _load(REPO / "scripts/senior/same_section_linkage.py", "linkage_xenium_scope")
from scripts.senior import modality_concordance as concordance  # noqa: E402


@pytest.fixture
def fully_excluded(monkeypatch):
    """6457 as it was before 2026-10-22: excluded with no declared platform."""
    monkeypatch.delitem(cohort.EXCLUDED_PLATFORMS, "6457")
    assert excluded_on("6457", "xenium")


def test_the_ruling_this_file_rests_on():
    assert excluded_on("6457", "phenocycler")
    assert not excluded_on("6457", "xenium") and not excluded_on("6457", "if")
    with pytest.raises(ValueError, match="excluded donor"):
        ensure_eligible_donors(["6457"], context="PhenoCycler")


def _h5ads(root: Path, names: list[str]) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for name in names:
        (root / name).touch()
    return root


def test_the_pancreas_xenium_typer_types_6457(tmp_path, monkeypatch):
    monkeypatch.setattr(panc_worker, "H5AD_DIR", _h5ads(
        tmp_path, ["6457_PancXenium.h5ad", "6533_PancXenium.h5ad", "HDL112_PancXenium.h5ad"]))
    assert [(d, s) for d, s, _ in panc_worker.sections()] == [
        ("6457", "6457"), ("6533", "6533"), ("112", "HDL112")]


def test_the_pancreas_xenium_typer_still_fails_closed(tmp_path, monkeypatch, fully_excluded):
    monkeypatch.setattr(panc_worker, "H5AD_DIR", _h5ads(
        tmp_path, ["6457_PancXenium.h5ad", "6533_PancXenium.h5ad"]))
    assert [d for d, _, _ in panc_worker.sections()] == ["6533"]


def test_the_pln_xenium_typer_types_both_6457_sections(tmp_path, monkeypatch):
    monkeypatch.setattr(pln_v3, "H5AD_DIR", _h5ads(
        tmp_path, ["6457a_pLNXenium.h5ad", "6457b_pLNXenium.h5ad", "6533_pLNXenium.h5ad"]))
    assert [(d, s) for d, s, _ in pln_v3.sections()] == [
        ("6457", "6457a"), ("6457", "6457b"), ("6533", "6533")]
    assert pln_v3.refused_donors() == []


def test_the_pln_xenium_typer_still_fails_closed(tmp_path, monkeypatch, fully_excluded):
    monkeypatch.setattr(pln_v3, "H5AD_DIR", _h5ads(
        tmp_path, ["6457a_pLNXenium.h5ad", "6533_pLNXenium.h5ad"]))
    assert [d for d, _, _ in pln_v3.sections()] == ["6533"]
    assert pln_v3.refused_donors() == ["6457"]


def test_the_pln_xenium_clustering_admits_6457():
    assert cluster.xen_eligible({"6457", "6533"}) == ({"6457", "6533"}, [])


def test_the_pln_xenium_clustering_still_fails_closed(fully_excluded):
    assert cluster.xen_eligible({"6457", "6533"}) == ({"6533"}, ["6457"])


def test_a_restain_linkage_reads_6457_and_a_same_slide_one_would_not():
    """The re-stain is Xenium and IF only; the same-slide lane lifts a PhenoCycler call."""
    assert linkage.excluded_platforms("6457", same_slide=False) == []
    assert linkage.excluded_platforms("6457", same_slide=True) == ["phenocycler"]
    row = pd.Series({"centrally_excluded": True, "excluded_platforms": "phenocycler"})
    assert linkage.verdict(row) == ("nothing", "centrally excluded on phenocycler")


def test_the_linkage_still_fails_closed(fully_excluded):
    assert linkage.excluded_platforms("6457", same_slide=False) == ["xenium", "if"]


def test_concordance_reads_6457s_rna_and_refuses_its_protein():
    by_side = concordance.excluded_by_modality()
    assert "6457" in by_side["protein"] and "6457" not in by_side["rna"]


def test_concordance_still_fails_closed(fully_excluded):
    by_side = concordance.excluded_by_modality()
    assert "6457" in by_side["protein"] and "6457" in by_side["rna"]


class _Release:
    """Just what ``pln_rna_counts`` asks of a release before it opens a file."""

    available = True
    release_id = "fixture-pln"

    def require_column(self, role):
        return role

    def column(self, role):
        return role if role in {"cell_key", "common_lineage"} else None

    def calls_schema(self):
        return type("Schema", (), {"names": ["cell_id", "common_lineage"]})()

    def donor_partitions(self):
        return {"6457": []}


def test_concordance_admits_6457_to_the_rna_side_only(tmp_path, monkeypatch):
    """The protein pLN side names 6457 refused without reading its partition; the RNA side
    admits it and goes on to the next question (here, that the release holds no partition)."""
    by_side = concordance.excluded_by_modality()
    pheno = tmp_path / "pheno"
    pheno.mkdir()
    pd.DataFrame({"donor_id": ["6457"], "state": ["typed"], "reason": [""]}).to_csv(
        pheno / "per_donor_status.csv", index=False)
    monkeypatch.setattr(concordance, "PLN_PHENO", pheno)
    _, protein_states = concordance.pln_protein_counts(by_side["protein"])
    assert protein_states == {"6457": "donor centrally excluded"}

    monkeypatch.setattr(concordance, "xenium_source", lambda tissue: _Release())
    _, rna_states = concordance.pln_rna_counts(by_side["rna"])
    assert rna_states == {"6457": "no partition under fixture-pln"}


roster = _load(REPO / "scripts/senior/build_cohort_roster.py", "cohort_roster_xenium_scope")


def test_the_app_roster_records_the_platforms_each_exclusion_covers():
    built = roster.build_roster()
    assert built["platforms"] == list(cohort.PLATFORMS)
    entry = next(e for e in built["excluded_donors"] if e["donor_id"] == "6457")
    assert entry["platforms"] == ["phenocycler"]
    # The candidates are the PhenoCycler pancreas roster, so eligibility stays donor-wide.
    assert "6457" in built["all_donors"] and "6457" not in built["eligible_donors"]


def test_the_app_roster_records_an_unscoped_exclusion_as_every_platform(fully_excluded):
    entry = next(e for e in roster.build_roster()["excluded_donors"] if e["donor_id"] == "6457")
    assert entry["platforms"] == list(cohort.PLATFORMS)
