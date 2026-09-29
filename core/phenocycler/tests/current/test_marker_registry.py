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
        # redsea from a30a517 to 2026-09-29; the raw Cell mean, like every marker, since.
        assert marker.measurement_source == "qupath_raw"
        assert marker.compartment == "Cell"
        # Uniform source, unchanged calibration status: still no verified reference.
        assert marker.reference is None
        assert not marker.active
    with pytest.raises(ValueError, match="no active calibration"):
        registry.marker_alpha("Iba1")


def test_every_marker_reads_its_raw_qupath_mean():
    """The 2026-09-29 policy: no marker earned REDSEA, so none is compensated or read from it.

    Before it, 63 markers were ``redsea_corrected``, Collagen_IV was ``redsea_passthrough``
    and the eight nuclear markers were already ``uncompensated``.
    """
    registry = load_registry()

    assert len(registry.markers) == 72
    assert {marker.spillover_policy for marker in registry.markers} == {"uncompensated"}
    assert {marker.measurement_source for marker in registry.markers} == {"qupath_raw"}
    assert registry.redsea_markers == ()
    compartments = [marker.compartment for marker in registry.markers]
    assert compartments.count("Cell") == 63
    assert sorted(m.name for m in registry.markers if m.compartment == "Nucleus") == [
        "DAPI", "FOXP3", "Histone_H3_(p_Ser28)", "Ki67", "PCNA", "SOX2", "TCF_1", "TOX", "TP63",
    ]
    # REDSEA therefore compensates nothing, and the vocabulary still offers it.
    assert not marker_taxonomy.spillover_corrected_mask(
        [marker.name for marker in registry.markers]
    ).any()
    assert {"redsea_corrected", "redsea_passthrough"} <= set(
        registry.allowed_spillover_policies
    )


@pytest.mark.parametrize("compartment", ["Cell", "Nucleus"])
def test_uncompensated_is_accepted_where_ingest_carries_the_raw_mean(compartment):
    payload = _payload()
    cd3e = next(row for row in payload["markers"] if row["name"] == "CD3e")
    cd3e["compartment"] = compartment

    marker = registry_from_dict(payload).marker("CD3e")
    assert (marker.spillover_policy, marker.measurement_source) == (
        "uncompensated",
        "qupath_raw",
    )


@pytest.mark.parametrize("compartment", ["Cytoplasm", "Membrane"])
def test_uncompensated_is_refused_where_there_is_no_raw_mean(compartment):
    """Relaxed from Nucleus-only on 2026-09-29, but a raw source must still exist: ingest
    keeps only the Cell and Nucleus means, so Cytoplasm and Membrane are REDSEA outputs."""

    payload = _payload()
    cd3e = next(row for row in payload["markers"] if row["name"] == "CD3e")
    cd3e["compartment"] = compartment
    with pytest.raises(ValueError, match="ingest carries only for \\['Cell', 'Nucleus'\\]"):
        registry_from_dict(payload)

    # Routed to REDSEA, the same compartment is valid: REDSEA writes it.
    cd3e["spillover_policy"] = "redsea_corrected"
    assert registry_from_dict(payload).marker("CD3e").measurement_column == (
        f"{compartment}__CD3e"
    )


def test_a_marker_can_be_routed_back_to_redsea_by_the_registry_alone():
    payload = _payload()
    for row in payload["markers"]:
        if row["name"] in {"CD3e", "Collagen_IV"}:
            row["spillover_policy"] = (
                "redsea_corrected" if row["name"] == "CD3e" else "redsea_passthrough"
            )
    registry = registry_from_dict(payload)

    assert [marker.name for marker in registry.redsea_markers] == ["CD3e", "Collagen_IV"]
    assert {marker.measurement_source for marker in registry.redsea_markers} == {"redsea"}

    # A REDSEA policy cannot be paired with the raw source, nor the raw policy with REDSEA.
    mixed = _payload()
    cd3e = next(row for row in mixed["markers"] if row["name"] == "CD3e")
    cd3e["measurement_source"] = "redsea"
    with pytest.raises(ValueError, match="redsea measurement source requires a REDSEA"):
        registry_from_dict(mixed)
    cd3e.update(spillover_policy="redsea_corrected", measurement_source="qupath_raw")
    with pytest.raises(ValueError, match="qupath_raw measurement source requires uncompensated"):
        registry_from_dict(mixed)


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
