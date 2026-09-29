"""Expression selects one registry-declared value per marker.

Since 2026-09-29 every marker is ``uncompensated``: the selected value is QuPath's raw mean in
the marker's declared compartment, read from ``cells`` under the name ingest gave it, and the
stage reads no REDSEA partition.  The REDSEA route is still exercised here through a registry
that sends a marker back to it, because the vocabulary keeps it.
"""

from __future__ import annotations

import csv
import json

import duckdb
import pandas as pd
import pytest

from phenocycler import cells_parquet
from phenocycler.artifacts import ContractError
from phenocycler.expression import (
    build_selected_expression,
    raw_measurement_column,
    read_single_partition,
)
from phenocycler.marker_registry import (
    DEFAULT_REGISTRY_PATH,
    RAW_QUPATH_COMPARTMENTS,
    load_registry,
    registry_from_dict,
)


def _frames():
    # Column names are ingest's: the Cell mean is the bare marker name, the Nucleus mean is
    # ``Nucleus__<marker>``. Each marker also carries its OTHER compartment, with values that
    # would be caught if the wrong one were read.
    cells = pd.DataFrame(
        {
            "object_id": ["c2", "c1", "c3"],
            "image": ["img", "img", "img"],
            "X_centroid": [2.0, 1.0, 3.0],
            "Y_centroid": [12.0, 11.0, 13.0],
            "CD3e": [20.0, 10.0, 30.0],
            "Nucleus__CD3e": [200.0, 100.0, 300.0],
            "E_cadherin": [2.0, 1.0, 3.0],
            "Nucleus__E_cadherin": [0.2, 0.1, 0.3],
            "Ki67": [40.0, 50.0, 60.0],
            "Nucleus__Ki67": [4.0, 5.0, 6.0],
        }
    )
    qc = pd.DataFrame(
        {
            "object_id": ["c3", "c1", "c2"],
            "analysis_eligible": [True, False, True],
            "estimation_eligible": [True, False, False],
        }
    )
    return cells, qc


def _registry_with(**policies):
    """The shipped registry with some markers' spillover policy replaced."""

    payload = json.loads(DEFAULT_REGISTRY_PATH.read_text(encoding="utf-8"))
    for row in payload["markers"]:
        if row["name"] in policies:
            row["spillover_policy"] = policies.pop(row["name"])
    assert not policies, f"not registered: {sorted(policies)}"
    return registry_from_dict(payload)


def test_every_marker_is_read_raw_from_its_declared_compartment():
    cells, qc = _frames()
    registry = load_registry()
    expression, availability = build_selected_expression(
        cells,
        None,
        qc,
        donor_id="D1",
        acquired_markers={"CD3e", "E_cadherin", "Ki67"},
        registry=registry,
    )

    # Cell order is authoritative and every joined source is aligned to it.
    assert expression["object_id"].tolist() == ["c2", "c1", "c3"]
    assert expression["qc_analysis_eligible"].tolist() == [True, False, True]
    # Cell markers read the bare Cell mean, the nuclear marker its Nucleus mean.
    assert expression["CD3e"].tolist() == [20.0, 10.0, 30.0]
    assert expression["E_cadherin"].tolist() == [2.0, 1.0, 3.0]
    assert expression["Ki67"].tolist() == [4.0, 5.0, 6.0]

    audit = availability.set_index("marker")
    assert audit.loc["CD3e", "source"] == "qupath:CD3e"
    assert audit.loc["E_cadherin", "source"] == "qupath:E_cadherin"
    assert audit.loc["Ki67", "source"] == "qupath:Nucleus__Ki67"
    assert set(audit["measurement_source"]) == {"qupath_raw"}
    assert set(audit["spillover_policy"]) == {"uncompensated"}
    assert not audit["source"].str.startswith("redsea:").any()

    # A panel-absent marker is explicit null, not zero or a dropped column.
    assert expression["SST"].isna().all()
    assert not audit.loc["SST", "available"]
    assert audit.loc["SST", "reason"] == "marker_not_acquired_in_panel"
    assert len(expression.columns) == 4 + len(registry.markers) + 3  # keys/QC + markers + context


def test_histone_h3_reads_the_column_ingest_wrote_not_its_registry_spelling():
    """The registry keeps the qptiff spelling ``Histone_H3_(p_Ser28)``; ingest's sanitiser drops
    the parentheses, so the raw nucleus mean is ``Nucleus__Histone_H3_p_Ser28``.  The old lookup
    asked for ``Nucleus__Histone_H3_(p_Ser28)``, which no cells table carries."""

    registry = load_registry()
    histone = registry.marker("Histone_H3_(p_Ser28)")
    assert histone.compartment == "Nucleus"
    assert raw_measurement_column(histone) == "Nucleus__Histone_H3_p_Ser28"

    cells, qc = _frames()
    cells["Nucleus__Histone_H3_p_Ser28"] = [7.0, 8.0, 9.0]
    cells["Histone_H3_p_Ser28"] = [70.0, 80.0, 90.0]
    expression, availability = build_selected_expression(
        cells,
        None,
        qc,
        donor_id="D1",
        acquired_markers={"Histone_H3_(p_Ser28)"},
        registry=registry,
    )
    assert expression["Histone_H3_(p_Ser28)"].tolist() == [7.0, 8.0, 9.0]
    audit = availability.set_index("marker").loc["Histone_H3_(p_Ser28)"]
    assert audit["source"] == "qupath:Nucleus__Histone_H3_p_Ser28"

    # Only the unsanitised spelling present: that is a different column, not a fallback.
    unsanitised = cells.drop(
        columns=["Nucleus__Histone_H3_p_Ser28", "Histone_H3_p_Ser28"]
    ).assign(**{"Nucleus__Histone_H3_(p_Ser28)": [7.0, 8.0, 9.0]})
    with pytest.raises(ContractError, match="Nucleus__Histone_H3_p_Ser28"):
        build_selected_expression(
            unsanitised,
            None,
            qc,
            donor_id="D1",
            acquired_markers={"Histone_H3_(p_Ser28)"},
            registry=registry,
        )


def test_the_raw_lookup_finds_every_marker_where_ingest_put_it(tmp_path):
    """Through the real ingest SQL: QuPath names with punctuation, both compartments."""

    header = [
        "Image", "Object ID", "Centroid X µm", "Centroid Y µm", "Parent",
        "Cell: Area µm^2", "Cell: Solidity",
        "Cell: Histone H3 (p Ser28): Mean", "Nucleus: Histone H3 (p Ser28): Mean",
        "Cell: Ki67: Mean", "Nucleus: Ki67: Mean",
        "Cell: E-cadherin: Mean", "Nucleus: E-cadherin: Mean",
        "Cell: Ker8/18: Mean", "Nucleus: Ker8/18: Mean",
        "Cell: Pan-Cytokeratin: Mean", "Nucleus: Pan-Cytokeratin: Mean",
    ]
    rows = []
    for index, object_id in enumerate(("u1", "u2", "u3")):
        row = {
            "Image": "img.qptiff", "Object ID": object_id, "Centroid X µm": "1",
            "Centroid Y µm": "2", "Parent": "Annotation (Tissue)",
            "Cell: Area µm^2": "60", "Cell: Solidity": "0.9",
        }
        for offset, column in enumerate(header[7:]):
            row[column] = str(100 * (offset + 1) + index)
        rows.append(row)
    export = tmp_path / "cells.csv"
    with open(export, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    out = tmp_path / "data"
    out.mkdir()  # PARTITION_BY needs the parent, as build_cells_parquet ensures
    sql, _ = cells_parquet.build_sql(
        [export], out, None, 0.0, image_to_donor={"img.qptiff": "D1"}
    )
    duckdb.connect().execute(sql)
    cells = read_single_partition(out / "cells", "D1")

    qc = pd.DataFrame(
        {
            "object_id": ["u1", "u2", "u3"],
            "analysis_eligible": True,
            "estimation_eligible": True,
        }
    )
    expected = {  # registry name -> the QuPath column it must equal
        "Histone_H3_(p_Ser28)": "Nucleus: Histone H3 (p Ser28): Mean",
        "Ki67": "Nucleus: Ki67: Mean",
        "E_cadherin": "Cell: E-cadherin: Mean",
        "Ker8_18": "Cell: Ker8/18: Mean",
        "Pan_Cytokeratin": "Cell: Pan-Cytokeratin: Mean",
    }
    expression, _ = build_selected_expression(
        cells,
        None,
        qc,
        donor_id="D1",
        acquired_markers=set(expected),
        registry=load_registry(),
    )
    by_object = expression.set_index("object_id").loc[["u1", "u2", "u3"]]
    for marker, qupath_column in expected.items():
        assert by_object[marker].tolist() == [float(row[qupath_column]) for row in rows], marker


def test_every_registered_marker_has_one_distinct_raw_column():
    registry = load_registry()
    columns = [raw_measurement_column(marker) for marker in registry.markers]
    assert len(set(columns)) == len(columns) == len(registry.markers)
    assert {marker.compartment for marker in registry.markers} <= set(RAW_QUPATH_COMPARTMENTS)
    # Only Histone H3's registry spelling differs from ingest's.
    renamed = {
        marker.name
        for marker in registry.markers
        if cells_parquet.sanitize(marker.name) != marker.name
    }
    assert renamed == {"Histone_H3_(p_Ser28)"}


def test_declared_acquired_marker_missing_from_its_source_is_fatal():
    cells, qc = _frames()
    with pytest.raises(ContractError, match="raw QuPath Cell mean, cells column 'SST'"):
        build_selected_expression(
            cells,
            None,
            qc,
            donor_id="D1",
            acquired_markers={"CD3e", "E_cadherin", "Ki67", "SST"},
        )


@pytest.mark.parametrize("marker", ["Iba1", "CD68", "M2Gal3"])
def test_former_qupath_norm_markers_read_the_raw_mean_like_every_other_marker(marker):
    """These three took QuPath's signed ``Mean (norm)`` until the RESTORE-era exception was
    retired, then REDSEA's ``Cell__<marker>``.  They now read the bare raw Cell mean like every
    other marker, and neither a stray ``QuPathNorm__`` column nor a REDSEA-shaped ``Cell__``
    column in ``cells`` may be selected instead.
    """
    cells, qc = _frames()
    cells[f"QuPathNorm__Cell__{marker}"] = [-0.25, 0.4, 1.3]
    cells[f"Cell__{marker}"] = [900.0, 800.0, 700.0]
    cells[marker] = [5.0, 6.0, 7.0]
    expression, availability = build_selected_expression(
        cells,
        None,
        qc,
        donor_id="D1",
        acquired_markers={"CD3e", "E_cadherin", "Ki67", marker},
    )

    assert expression[marker].tolist() == [5.0, 6.0, 7.0]
    audit = availability.set_index("marker").loc[marker]
    assert audit["source"] == f"qupath:{marker}"
    assert audit["measurement_source"] == "qupath_raw"
    assert audit["spillover_policy"] == "uncompensated"


def test_a_supplied_redsea_partition_is_not_read_when_no_marker_uses_it():
    cells, qc = _frames()
    redsea = pd.DataFrame(
        {"object_id": ["c1", "c2", "c3"], "Cell__CD3e": [-1.0, -1.0, -1.0]}
    )
    with_redsea, _ = build_selected_expression(
        cells, redsea, qc, donor_id="D1", acquired_markers={"CD3e"}
    )
    without, _ = build_selected_expression(
        cells, None, qc, donor_id="D1", acquired_markers={"CD3e"}
    )
    pd.testing.assert_frame_equal(with_redsea, without)
    assert without["CD3e"].tolist() == [20.0, 10.0, 30.0]


def test_a_marker_routed_back_to_redsea_reads_it_and_needs_the_partition():
    """The REDSEA policies stay in the vocabulary; a marker sent back to one reads REDSEA."""

    registry = _registry_with(CD3e="redsea_corrected")
    assert [marker.name for marker in registry.redsea_markers] == ["CD3e"]
    cells, qc = _frames()
    redsea = pd.DataFrame(
        {"object_id": ["c1", "c2", "c3"], "Cell__CD3e": [11.0, 22.0, 33.0]}
    )
    expression, availability = build_selected_expression(
        cells,
        redsea,
        qc,
        donor_id="D1",
        acquired_markers={"CD3e", "Ki67"},
        registry=registry,
    )
    # redsea rows are c1,c2,c3 and the authoritative cell order is c2,c1,c3.
    assert expression["CD3e"].tolist() == [22.0, 11.0, 33.0]
    assert expression["Ki67"].tolist() == [4.0, 5.0, 6.0]
    audit = availability.set_index("marker")
    assert audit.loc["CD3e", "source"] == "redsea:Cell__CD3e"
    assert audit.loc["CD3e", "measurement_source"] == "redsea"

    with pytest.raises(ContractError, match="no REDSEA partition was supplied"):
        build_selected_expression(
            cells, None, qc, donor_id="D1", acquired_markers={"CD3e"}, registry=registry
        )
    with pytest.raises(ContractError, match="requires REDSEA column 'Cell__CD3e'"):
        build_selected_expression(
            cells,
            redsea.drop(columns="Cell__CD3e"),
            qc,
            donor_id="D1",
            acquired_markers={"CD3e"},
            registry=registry,
        )
    with pytest.raises(ContractError, match="object universe differs"):
        build_selected_expression(
            cells,
            redsea.iloc[:-1],
            qc,
            donor_id="D1",
            acquired_markers={"CD3e"},
            registry=registry,
        )


def test_expression_rejects_partial_object_universes():
    cells, qc = _frames()
    with pytest.raises(ContractError, match="object universe differs"):
        build_selected_expression(
            cells,
            None,
            qc.iloc[:-1],
            donor_id="D1",
            acquired_markers={"CD3e", "E_cadherin", "Ki67"},
        )
