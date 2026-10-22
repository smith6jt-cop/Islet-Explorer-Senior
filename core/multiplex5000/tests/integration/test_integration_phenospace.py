"""Declared phenotype spaces/releases: strict mappings and provenance refusal tests."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from phenocycler.config import PipelineConfig
from phenocycler.integration.contract import coerce_cell_table, contract_profile, has_identity
from phenocycler.integration.export_pheno import ExportError, attach_identity
from phenocycler.integration.phenospace import (
    PhenotypeSpaceError, attach_release_identity, load_phenotype_space, load_release,
)


def _space(path, *, beta_common="Endocrine"):
    payload = {
        "schema_version": 1,
        "space_id": f"space-{beta_common}",
        "producer_lane": "synthetic-current-lane",
        "labels": [
            {"label": "Endo-vascular", "common_lineage": "Vascular"},
            {"label": "Endo-beta", "common_lineage": beta_common},
            {"label": "Endo-beta-cycling", "parent": "Endo-beta", "inherit": True},
            {"label": "Ambiguous", "common_resolution": "ambiguous", "reason": "tie"},
            {"label": "Unavailable", "common_resolution": "unavailable", "reason": "absent"},
            {"label": "AMBIGUOUS-CD11b-pos-epithelial",
             "common_resolution": "not_mappable", "reason": "winning hybrid call"},
        ],
    }
    path.write_text(json.dumps(payload))
    return path


def _release(
    tmp_path,
    *,
    donors=("6436",),
    content_id="content-ok",
    input_path="expression",
    modality="phenocycler",
):
    release_dir = tmp_path / "release"
    release_dir.mkdir(parents=True, exist_ok=True)
    _space(release_dir / "label_space.json")
    payload = {
        "schema_version": 1,
        "release_id": "synthetic-release-v1",
        "modality": modality,
        "producer": {"lane": "synthetic-current-lane", "commit": "abc123"},
        "label_space": "label_space.json",
        "calls_dir": "calls",
        "cell_key": "object_id" if modality == "phenocycler" else "cell_id",
        "donors": list(donors),
        "inputs": [{
            "kind": "expression",
            "run_id": "run1",
            "content_id": content_id,
            "path": input_path,
        }],
        "columns": {"label": "phenotype", "common_lineage": "common_lineage"},
    }
    (release_dir / "release.json").write_text(json.dumps(payload))
    return release_dir


def _cfg(tmp_path, release_dir, *, modality="phenocycler"):
    data = tmp_path / "data"
    data.mkdir(parents=True, exist_ok=True)
    kwargs = {f"{modality}_identity_release": release_dir}
    return PipelineConfig(data_dir=data, **kwargs)


def test_real_prefix_trap_is_resolved_by_parent_declaration(tmp_path):
    space = load_phenotype_space(_space(tmp_path / "space.json"))
    assert space.common_for("Endo-vascular") == ("Vascular", "mapped")
    assert space.common_for("Endo-beta") == ("Endocrine", "mapped")
    assert space.common_for("Endo-beta-cycling") == ("Endocrine", "mapped")


def test_ambiguous_unavailable_and_not_mappable_stay_distinct(tmp_path):
    space = load_phenotype_space(_space(tmp_path / "space.json"))
    assert space.common_for("Ambiguous") == ("", "ambiguous")
    assert space.common_for("Unavailable") == ("", "unavailable")
    assert space.common_for("AMBIGUOUS-CD11b-pos-epithelial") == ("", "not_mappable")


def test_undeclared_label_raises_and_spaces_do_not_mutate_module_state(tmp_path):
    a = load_phenotype_space(_space(tmp_path / "a.json", beta_common="Endocrine"))
    b = load_phenotype_space(_space(tmp_path / "b.json", beta_common="Other"))
    assert a.common_for("Endo-beta")[0] == "Endocrine"
    assert b.common_for("Endo-beta")[0] == "Other"
    assert a.common_for("Endo-beta")[0] == "Endocrine"
    with pytest.raises(PhenotypeSpaceError, match="undeclared"):
        a.common_for("Endo-by-prefix")


def test_release_round_trip_records_whether_provenance_could_be_checked(tmp_path):
    """content_id is no longer compared to anything, and that outcome is now VISIBLE.

    It used to be checked against data/runs/<run_id>/manifests/expression.json. Stage
    manifests were removed on 2026-09-16, so a release naming a wrong content id loads with
    no error. That it was possible was tolerable; that it was SILENT was not. The release now
    carries `provenance`, and a synthetic release whose declared input does not resolve says
    so rather than looking checked.
    """

    release_dir = _release(tmp_path)
    release = load_release(_cfg(tmp_path, release_dir), "phenocycler")
    assert release.release_id == "synthetic-release-v1"
    assert release.run_id == "run1"
    assert release.expression_content_id == "content-ok"
    assert release.provenance.startswith("unverifiable:"), release.provenance
    assert "does not resolve" in release.provenance

    # A wrong content id still loads -- nothing compares it -- but it is not passed off as
    # verified either.
    bad_dir = _release(tmp_path / "bad", content_id="wrong")
    loaded = load_release(_cfg(tmp_path / "bad", bad_dir), "phenocycler")
    assert loaded.expression_content_id == "wrong"
    assert loaded.provenance.startswith("unverifiable:")


def test_a_resolvable_input_is_checked_for_donor_coverage(tmp_path):
    """When the declared input exists, it must contain every donor the release claims."""

    expression = tmp_path / "expr"
    for donor in ("6436", "6442"):
        (expression / f"donor_id={donor}").mkdir(parents=True)
    release_dir = _release(tmp_path, donors=("6436", "6442"),
                           input_path=str(expression))
    release = load_release(_cfg(tmp_path, release_dir), "phenocycler")
    assert release.provenance == f"verified:{expression}"


def test_a_release_claiming_a_donor_its_input_lacks_is_refused(tmp_path):
    """The realistic wrong-id failure, and the one case worth failing on.

    No configuration explains a release claiming calls for a donor its own declared input
    does not contain. This needs no content hash and no manifest -- and unlike the check it
    replaces, it cannot fire because package source changed.
    """

    expression = tmp_path / "expr"
    (expression / "donor_id=6436").mkdir(parents=True)
    release_dir = _release(tmp_path, donors=("6436", "6442"),
                           input_path=str(expression))
    with pytest.raises(PhenotypeSpaceError, match="does not contain"):
        load_release(_cfg(tmp_path, release_dir), "phenocycler")


def test_donor_coverage_tolerates_the_zero_padding_split(tmp_path):
    """0112 and 112 are one donor; the pLNpanc partitions pad and the cohort does not."""

    expression = tmp_path / "expr"
    (expression / "donor_id=0112").mkdir(parents=True)
    release_dir = _release(tmp_path, donors=("112",), input_path=str(expression))
    release = load_release(_cfg(tmp_path, release_dir), "phenocycler")
    assert release.provenance.startswith("verified:")


def test_attach_identity_is_exact_uuid_join_and_stamps_release(tmp_path):
    release_dir = _release(tmp_path)
    calls_dir = release_dir / "calls" / "donor_id=6436"
    calls_dir.mkdir(parents=True)
    pd.DataFrame({
        "object_id": ["u2", "u1", "u3"],  # deliberately not geometry order
        "phenotype": ["Endo-vascular", "Endo-beta", "Unavailable"],
        "common_lineage": ["Vascular", "Endocrine", ""],
        "evidence_positive": [True, True, pd.NA],
    }).to_parquet(calls_dir / "data.parquet", index=False)
    release = load_release(_cfg(tmp_path, release_dir), "phenocycler")
    geometry = coerce_cell_table(
        pd.DataFrame({
            "cell_id": ["u1", "u2", "u3"],
            "x_um": [1.0, 2.0, 3.0],
            "y_um": [4.0, 5.0, 6.0],
        }),
        modality="phenocycler", donor_id="6436", roi="panc", profile="geometry",
    )

    attached = attach_identity(geometry, release, "6436")
    assert list(attached["lineage_native"]) == ["Endo-beta", "Endo-vascular", "Unavailable"]
    assert contract_profile(attached) == "identity"
    assert has_identity(attached)
    assert attached.attrs["identity_release_id"] == release.release_id

    pd.DataFrame({
        "object_id": ["u1", "u2"],
        "phenotype": ["Endo-beta", "Endo-vascular"],
    }).to_parquet(calls_dir / "data.parquet", index=False)
    with pytest.raises(ExportError, match="UUID mismatch"):
        attach_identity(geometry, release, "6436")


def test_xenium_release_uses_the_same_exact_cell_set_attachment(tmp_path):
    release_dir = _release(tmp_path, modality="xenium")
    calls_dir = release_dir / "calls" / "donor_id=6436" / "roi=panc"
    calls_dir.mkdir(parents=True)
    pd.DataFrame({
        "cell_id": ["x2", "x1"],
        "phenotype": ["Endo-vascular", "Endo-beta"],
        "common_lineage": ["Vascular", "Endocrine"],
        "evidence_positive": [True, False],
    }).to_parquet(calls_dir / "data.parquet", index=False)
    release = load_release(_cfg(tmp_path, release_dir, modality="xenium"), "xenium")
    geometry = coerce_cell_table(
        pd.DataFrame({
            "cell_id": ["x1", "x2"],
            "x_um": [1.0, 2.0],
            "y_um": [3.0, 4.0],
        }),
        modality="xenium", donor_id="6436", roi="panc", profile="geometry",
    )

    attached = attach_release_identity(geometry, release, "6436")
    assert list(attached["lineage_common"]) == ["Endocrine", "Vascular"]
    assert list(attached["evidence_positive"]) == [False, True]
    assert has_identity(attached)


def test_roi_partition_never_falls_through_to_sibling_calls(tmp_path):
    release_dir = _release(tmp_path, modality="xenium")
    sibling = release_dir / "calls" / "donor_id=6436" / "roi=pln_1"
    sibling.mkdir(parents=True)
    pd.DataFrame({
        "cell_id": ["x1"],
        "phenotype": ["Endo-beta"],
    }).to_parquet(sibling / "data.parquet", index=False)
    release = load_release(_cfg(tmp_path, release_dir, modality="xenium"), "xenium")
    geometry = coerce_cell_table(
        pd.DataFrame({"cell_id": ["x1"], "x_um": [1.0], "y_um": [2.0]}),
        modality="xenium", donor_id="6436", roi="panc", profile="geometry",
    )

    with pytest.raises(PhenotypeSpaceError, match="no calls for 6436/panc"):
        attach_release_identity(geometry, release, "6436")


@pytest.mark.parametrize(
    ("donors", "input_path", "message"),
    [
        (("6436", "6457"), "expression", "excluded donor"),
        (("6436",), "runs/run1/assignments/data.parquet", "assignments"),
        (("6436",), "restore_gated_redsea/calls.parquet", "retired"),
    ],
)
def test_release_refuses_exclusions_and_retired_inputs(tmp_path, donors, input_path, message):
    release_dir = _release(tmp_path, donors=donors, input_path=input_path)
    with pytest.raises(PhenotypeSpaceError, match=message):
        load_release(_cfg(tmp_path, release_dir), "phenocycler")


def test_a_phenocycler_only_exclusion_admits_the_donor_to_a_xenium_release(tmp_path):
    """6457's central exclusion covers its PhenoCycler data only (maintainer, 2026-10-22).

    The release's own modality is what is asked: a Xenium release may carry 6457, a
    PhenoCycler release still may not.
    """
    xenium = _release(tmp_path / "x", donors=("6436", "6457"), modality="xenium")
    release = load_release(_cfg(tmp_path / "x", xenium, modality="xenium"), "xenium")
    assert release.donors == ("6436", "6457")

    pheno = _release(tmp_path / "p", donors=("6436", "6457"))
    with pytest.raises(PhenotypeSpaceError, match="excluded donor.*on phenocycler.*6457"):
        load_release(_cfg(tmp_path / "p", pheno), "phenocycler")


def test_an_exclusion_without_a_platform_refuses_a_xenium_release(tmp_path, monkeypatch):
    from multiplex5000 import cohort

    monkeypatch.setitem(cohort.DONOR_EXCLUSIONS, "9999", "a reason")
    xenium = _release(tmp_path, donors=("6436", "9999"), modality="xenium")
    with pytest.raises(PhenotypeSpaceError, match="excluded donor.*on xenium: 9999: a reason"):
        load_release(_cfg(tmp_path, xenium, modality="xenium"), "xenium")
