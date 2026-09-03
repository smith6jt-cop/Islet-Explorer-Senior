"""Cell-table contract: schema, dtypes, guards, partition round-trip."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.contract import (
    CELL_TABLE_COLUMNS, CONTRACT_METADATA_IDENTITY_KEY, CONTRACT_METADATA_PROFILE_KEY,
    CONTRACT_METADATA_SCHEMA_KEY,
    CONTRACT_PROFILE_GEOMETRY, CONTRACT_PROFILE_IDENTITY, CONTRACT_SCHEMA_VERSION,
    CellTable, ContractError, coerce_cell_table, contract_profile, discover_partitions,
    empty_cell_table, feature_columns, has_identity, normalize_donor_id, partition_path,
    read_cell_table, validate_cell_table, write_cell_table,
)


def _minimal(n=10, modality="phenocycler"):
    return pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(n)],
        "x_um": np.linspace(0, 500, n),
        "y_um": np.linspace(0, 400, n),
        "lineage_native": ["Endocrine"] * n,
    })


def _identity(df, *, modality="phenocycler", donor_id="1", roi="panc", **kwargs):
    return coerce_cell_table(
        df,
        modality=modality,
        donor_id=donor_id,
        roi=roi,
        profile="identity",
        identity_release_id="synthetic-test-release",
        **kwargs,
    )


def test_coerce_fills_schema_and_dtypes():
    df = _identity(_minimal(), donor_id="6539")
    assert (df["tissue"] == "pancreas").all(), "tissue must be derived from the roi"
    for col, dtype in CELL_TABLE_COLUMNS.items():
        assert col in df.columns, col
        if dtype.startswith("float"):
            assert df[col].dtype == np.float32
        elif dtype == "boolean":
            assert str(df[col].dtype) == "boolean"
    assert (df["donor_id"] == "6539").all()
    assert (df["modality"] == "phenocycler").all()
    assert (df["roi"] == "panc").all()
    assert df["evidence_positive"].isna().all(), "missing evidence stays unavailable"


def test_coerce_preserves_extra_columns_and_sorts_features():
    raw = _minimal()
    raw["feat__CD31"] = 1.0
    raw["feat__ACTA2"] = 2.0
    raw["immune_subclass"] = "T"
    df = _identity(raw)
    # Features sorted and typed, provenance columns kept rather than dropped.
    assert feature_columns(df) == ["feat__ACTA2", "feat__CD31"]
    assert "immune_subclass" in df.columns
    assert df["feat__CD31"].dtype == np.float32


def test_missing_required_column_raises():
    bad = _minimal().drop(columns=["x_um"])
    with pytest.raises(ContractError, match="missing"):
        _identity(bad)


def test_geometry_only_table_is_valid_without_identity():
    """A table with no ``lineage_native`` is a legitimate product, not a broken one.

    This is what ``runs/<id>/expression`` can supply while no pair roster is verified, and the
    geometry stages must be able to consume it.
    """
    geom = _minimal().drop(columns=["lineage_native"])
    out = coerce_cell_table(
        geom, modality="phenocycler", donor_id="1", roi="panc", profile="geometry"
    )
    assert "lineage_native" in out.columns          # filled, not required
    assert (out["lineage_native"] == "").all()
    assert not has_identity(out)
    assert out["evidence_positive"].isna().all()
    assert contract_profile(out) == CONTRACT_PROFILE_GEOMETRY


def test_has_identity_is_about_declared_release_not_label_shaped_columns():
    labelled = _identity(_minimal())
    assert has_identity(labelled)

    blanked = labelled.copy()
    blanked["lineage_native"] = ""
    assert has_identity(blanked), "an all-Unavailable release is still an identity release"

    unproven = labelled.copy()
    unproven.attrs.clear()
    assert not has_identity(unproven)
    with pytest.raises(ContractError, match="without identity provenance"):
        validate_cell_table(unproven)


def test_profile_and_identity_release_are_explicit():
    with pytest.raises(TypeError, match="profile"):
        coerce_cell_table(_minimal(), modality="phenocycler", donor_id="1", roi="panc")
    with pytest.raises(ContractError, match="identity_release_id"):
        coerce_cell_table(
            _minimal(), modality="phenocycler", donor_id="1", roi="panc", profile="identity"
        )


def test_unknown_modality_raises():
    with pytest.raises(ContractError, match="unknown modality"):
        _identity(_minimal(), modality="visium")


def test_validate_rejects_mixed_modalities_and_duplicates():
    a = _identity(_minimal())
    b = _identity(_minimal(), modality="xenium")
    with pytest.raises(ContractError, match="mixes modalities"):
        validate_cell_table(pd.concat([a, b], ignore_index=True))

    dup = pd.concat([a, a], ignore_index=True)
    dup["modality"] = "phenocycler"
    with pytest.raises(ContractError, match="duplicate cell_id"):
        validate_cell_table(dup)


def test_validate_catches_pixels_masquerading_as_microns():
    """A tissue section is at most a few cm. A 10 cm span means someone passed pixels.

    This is the exact failure mode of the prototype being replaced, which treated Xenium's
    micron coordinates as 40x pixels — so it is worth failing loudly rather than producing a
    confidently misaligned overlay.
    """
    raw = _minimal(n=50)
    raw["x_um"] = np.linspace(0, 500_000, 50)
    with pytest.raises(ContractError, match="look like pixels"):
        _identity(raw, modality="xenium")


@pytest.mark.parametrize("value,expected", [
    (6539, "6539"), (6539.0, "6539"), ("6539", "6539"), ("6539.0", "6539"),
    (" 6539 ", "6539"), (np.int64(6539), "6539"), (None, ""), (np.nan, ""),
])
def test_donor_id_normalisation(value, expected):
    """The workbook stores donor ids numerically and PhenoCycler derives them by regex.

    Without normalisation `6539` and `6539.0` are different join keys and the whole
    cross-modality merge silently yields zero paired donors.
    """
    assert normalize_donor_id(value) == expected


def test_partition_roundtrip(tmp_path):
    df = _identity(_minimal(20), modality="xenium", donor_id="6539")
    df["feat__CD31"] = np.arange(20, dtype=float)
    df = _identity(df, modality="xenium", donor_id="6539")

    path = partition_path(tmp_path, "6539", "panc")
    write_cell_table(df, path)
    assert path.exists()

    table = read_cell_table(tmp_path, "6539", "panc", modality="xenium")
    assert table.n_cells == 20
    assert table.features == ["feat__CD31"]
    assert table.xy().shape == (20, 2)
    assert discover_partitions(tmp_path) == [("6539", "panc")]


def test_partition_stamps_schema_version_and_profile(tmp_path):
    import pyarrow.parquet as pq

    identity = _identity(_minimal(3), modality="xenium", donor_id="6539")
    identity_path = write_cell_table(identity, partition_path(tmp_path, "6539", "panc"))
    identity_meta = pq.read_schema(identity_path).metadata or {}
    assert identity_meta[CONTRACT_METADATA_SCHEMA_KEY] == str(CONTRACT_SCHEMA_VERSION).encode()
    assert identity_meta[CONTRACT_METADATA_PROFILE_KEY] == CONTRACT_PROFILE_IDENTITY.encode()
    assert identity_meta[CONTRACT_METADATA_IDENTITY_KEY] == b"synthetic-test-release"

    geometry = coerce_cell_table(
        _minimal(3).drop(columns="lineage_native"),
        modality="phenocycler", donor_id="6436", roi="panc", profile="geometry",
    )
    geometry_path = write_cell_table(geometry, partition_path(tmp_path, "6436", "panc"))
    geometry_meta = pq.read_schema(geometry_path).metadata or {}
    assert geometry_meta[CONTRACT_METADATA_PROFILE_KEY] == CONTRACT_PROFILE_GEOMETRY.encode()


def test_reader_refuses_unversioned_legacy_cell_table(tmp_path):
    df = _identity(_minimal(3), modality="xenium", donor_id="6539")
    path = partition_path(tmp_path, "6539", "panc")
    path.parent.mkdir(parents=True)
    df.to_parquet(path, index=False)
    with pytest.raises(ContractError, match="legacy integration caches"):
        read_cell_table(tmp_path, "6539", "panc", modality="xenium")


def test_registered_coordinates_guard(tmp_path):
    """Asking for registered coordinates before registration must raise, not fall back.

    Silently returning unregistered coordinates would produce a plausible-looking overlay of
    two unaligned sections, which is far worse than an error.
    """
    df = _identity(_minimal(), modality="xenium")
    t = CellTable(df=df, modality="xenium", donor_id="1", roi="panc")
    with pytest.raises(ContractError, match="register"):
        t.xy(registered=True)


def test_registered_coordinate_error_names_real_stages():
    """The guidance must name stages that exist, so it cannot rot as the pipeline changes.

    It previously said "run the register + transform stages", but `transform` is a library
    module, not a stage — and `x_um_reg` is actually written by `grid`.
    """
    from phenocycler.integration.pipeline import STAGES

    df = _identity(_minimal(), modality="xenium")
    t = CellTable(df=df, modality="xenium", donor_id="1", roi="panc")
    with pytest.raises(ContractError) as exc:
        t.xy(registered=True)

    message = str(exc.value)
    stage_names = {name for name, _pattern, _label in STAGES}
    quoted = set(re.findall(r"'([a-z_]+)'", message))
    assert quoted, "the error should name the stages to run, in quotes"
    assert quoted <= stage_names, f"error names non-existent stage(s): {quoted - stage_names}"
    assert {"register", "grid"} <= quoted


def test_empty_table_has_full_schema():
    df = empty_cell_table()
    assert list(df.columns) == list(CELL_TABLE_COLUMNS)
    validate_cell_table(df)      # zero rows is valid


def test_write_validates_before_persisting(tmp_path):
    bad = _identity(_minimal())
    bad.loc[0, "cell_id"] = bad.loc[1, "cell_id"]      # duplicate
    with pytest.raises(ContractError):
        write_cell_table(bad, tmp_path / "x.parquet")
    assert not (tmp_path / "x.parquet").exists()
