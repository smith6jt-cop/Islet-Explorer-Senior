from __future__ import annotations

import copy
import json

import pytest

from phenocycler import marker_taxonomy
from phenocycler.marker_registry import (
    DEFAULT_REGISTRY_PATH,
    NO_VERIFIED_TYPING_STATUS,
    NoVerifiedTypingError,
    load_registry,
    registry_from_dict,
)
from phenocycler.restore_stage import production_pairs


def _payload() -> dict:
    return json.loads(DEFAULT_REGISTRY_PATH.read_text(encoding="utf-8"))


def test_registry_covers_measurements_without_claiming_verified_pairs():
    registry = load_registry()

    assert registry.production_status == NO_VERIFIED_TYPING_STATUS
    assert {marker.name for marker in registry.type_markers} == set(marker_taxonomy.TYPE)
    assert {marker.name for marker in registry.state_markers} == set(marker_taxonomy.PROCESS)
    assert registry.active_markers == ()
    assert registry.broad_markers == ()
    assert registry.reference_graph == {}
    assert registry.families == {}
    assert all(marker.reference is None for marker in registry.markers)
    assert all(not marker.mutually_exclusive for marker in registry.markers)

    with pytest.raises(
        NoVerifiedTypingError, match="no verified production marker pairs"
    ):
        production_pairs(registry)


def test_no_marker_selects_the_retired_qupath_normalized_source():
    """Iba1, CD68 and M2Gal3 were contracted to QuPath's signed ``Mean (norm)`` for RESTORE.
    RESTORE is retired; the exception outlived it and put three markers on a different scale
    from the other 57 inside the same comparison, while demanding ``QuPathNorm__`` columns no
    export we hold actually carries.  Every marker now reads one source.
    """
    registry = load_registry()

    observed = {
        marker.name
        for marker in registry.markers
        if marker.measurement_source == "qupath_norm"
    }
    assert observed == set()
    for name in ("Iba1", "CD68", "M2Gal3"):
        marker = registry.marker(name)
        assert marker.measurement_source == "redsea"
        assert marker.measurement_column == f"Cell__{name}"
        # Uniform source, unchanged calibration status: still no verified reference.
        assert marker.reference is None
        assert not marker.active
    with pytest.raises(ValueError, match="no active calibration"):
        registry.marker_alpha("Iba1")


def test_registry_rejects_reintroducing_an_active_or_operative_pair():
    active = _payload()
    cd3e = next(row for row in active["markers"] if row["name"] == "CD3e")
    cd3e.update(
        calibration_status="active",
        reference="Ker8_18",
        families=["broad_immune"],
        mutually_exclusive=True,
    )
    active["families"] = {
        "broad_immune": {"familywise_alpha": 0.01, "members": ["CD3e"]}
    }
    with pytest.raises(ValueError, match="cannot contain active calibrations"):
        registry_from_dict(active)

    reference_only = _payload()
    cd79a = next(
        row for row in reference_only["markers"] if row["name"] == "CD79a"
    )
    cd79a["reference"] = "Ker8_18"
    with pytest.raises(ValueError, match="inactive marker cannot carry"):
        registry_from_dict(reference_only)


def test_registry_rejects_wrong_workflow_status_or_normalized_source_set():
    wrong_status = _payload()
    wrong_status["production_status"] = "production"
    with pytest.raises(ValueError, match="development_no_verified_typing"):
        registry_from_dict(wrong_status)

    reintroduced_source = copy.deepcopy(_payload())
    iba1 = next(row for row in reintroduced_source["markers"] if row["name"] == "Iba1")
    iba1["measurement_source"] = "qupath_norm"
    with pytest.raises(ValueError, match="no marker may select qupath_norm"):
        registry_from_dict(reintroduced_source)


def test_cd11b_remains_assay_qc_only():
    registry = load_registry()
    cd11b = registry.marker("CD11b")

    assert cd11b.kind == "excluded"
    assert cd11b.uses == ("assay_qc",)
    assert cd11b.calibration_status == "standalone"
    assert cd11b.reference is None
    assert cd11b.families == ()
