"""A donor's pancreas and lymph nodes are separate sections in one cohort manifest."""

from __future__ import annotations

import pytest

from phenocycler.artifacts import ContractError, QuPathCohortManifest

from test_artifacts import _cohort_image


def _cohort(tmp_path, keys):
    images = [
        _cohort_image(
            tmp_path,
            donor_id=key,
            image_id=f"{key}_Scan1.er.qptiff",
        )
        for key in keys
    ]
    return QuPathCohortManifest.create(
        cohort_name="islet-phenocycler",
        expected_donors=list(keys),
        images=images,
    )


def test_pancreas_and_four_lymph_nodes_are_representable(tmp_path):
    keys = ["6539", "6539pln", "6539pln2", "6539pln3", "6539pln4"]
    cohort = _cohort(tmp_path, keys)
    cohort.verify_identity()
    assert list(cohort.section_keys) == keys
    assert cohort.donors() == ("6539",)


def test_sections_resolve_to_their_donor_and_roi(tmp_path):
    """Stored keys are canonical, not specification-input, order."""
    from phenocycler.sections import parse

    cohort = _cohort(tmp_path, ["6539", "6539pln", "6539pln2", "6414"])
    assert cohort.donors() == ("6414", "6539")
    assert list(cohort.section_keys) == ["6414", "6539", "6539pln", "6539pln2"]
    assert [parse(key).roi for key in cohort.section_keys] == [
        "panc",
        "panc",
        "pln_1",
        "pln_2",
    ]


def test_repeated_bare_donor_id_teaches_the_section_key_fix(tmp_path):
    images = [
        _cohort_image(
            tmp_path,
            donor_id="6539",
            image_id=f"6539{suffix}_Scan1.er.qptiff",
        )
        for suffix in ("", "pLN")
    ]
    with pytest.raises(ContractError) as excinfo:
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=["6539"],
            images=images,
        )
    message = str(excinfo.value)
    assert "SECTION key" in message
    assert "6539pln" in message
    assert "lymph nodes" in message


def test_section_keys_is_an_alias_not_a_schema_rename(tmp_path):
    cohort = _cohort(tmp_path, ["6539", "6539pln"])
    assert cohort.section_keys == cohort.donor_ids


def test_content_id_distinguishes_a_second_section(tmp_path):
    pancreas = _cohort(tmp_path / "a", ["6539"])
    with_node = _cohort(tmp_path / "b", ["6539", "6539pln"])
    assert pancreas.content_id != with_node.content_id


def test_unknown_section_is_refused_at_manifest_creation(tmp_path):
    from phenocycler.qupath_manifest import _validate_section_keys

    good = [
        _cohort_image(tmp_path / "good", donor_id=key, image_id=f"{key}.qptiff")
        for key in ("6539", "6539pln")
    ]
    _validate_section_keys(good)

    bad = [
        _cohort_image(
            tmp_path / "bad",
            donor_id="6539spleen",
            image_id="x.qptiff",
        )
    ]
    with pytest.raises(ValueError, match="do not name a section"):
        _validate_section_keys(bad)


def test_template_is_consistent_and_shows_indexed_and_unindexed_nodes():
    from phenocycler.qupath_manifest import example_specification

    spec = example_specification()
    expected = set(spec["expected_donors"])
    assert expected == {row["donor_id"] for row in spec["images"]}
    assert {"6539pln", "6539pln2"} <= expected


def test_multi_section_cohort_flows_through_ingest_and_redsea(tmp_path):
    from phenocycler import load_config, redsea
    from phenocycler.artifacts import image_id_to_donor_id
    from phenocycler.cells_parquet import donor_case_expression
    from phenocycler.cohort import ensure_eligible_donors

    keys = ["6539", "6539pln", "6539pln2", "6414"]
    cohort = _cohort(tmp_path, keys)

    mapping = image_id_to_donor_id(cohort)
    assert sorted(mapping.values()) == sorted(keys)
    assert len(set(mapping.values())) == len(keys)

    sql = donor_case_expression(mapping)
    for key in keys:
        assert f"THEN '{key}'" in sql

    cfg = load_config()
    cfg.data_dir = tmp_path / "data"
    for key in keys:
        (cfg.cells_dir / f"donor_id={key}").mkdir(parents=True)
    assert cfg.discover_sections() == sorted(keys)
    assert cfg.discover_donors() == ["6414", "6539"]

    cfg.qupath_manifest = tmp_path / "qupath_manifest.json"
    cohort.write_json(cfg.qupath_manifest)
    for key in keys:
        assert redsea.donor_image(cfg, key) == f"{key}_Scan1.er.qptiff"

    with pytest.raises(ValueError, match="excluded donor"):
        ensure_eligible_donors([*keys, "6579pln"], context="test")
