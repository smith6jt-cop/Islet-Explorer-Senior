from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pyarrow as pa
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    DonorDatasetSnapshot,
    PartialArtifactError,
    SchemaField,
    StageManifest,
    StaleArtifactError,
    ContractError,
    hash_identifiers,
)
from phenocycler.expression import write_donor_partition
from phenocycler.recovery import (
    _link_manifest_dataset,
    validate_recovered_redsea_partition,
)


REPOSITORY = Path(__file__).resolve().parents[2]


def _expected(donor: str, object_ids: list[str]) -> DonorDatasetSnapshot:
    return DonorDatasetSnapshot(
        donor_id=donor,
        files=(),
        row_count=len(object_ids),
        schema=(SchemaField("object_id", "string", True),),
        schema_sha256="unused",
        object_id_sha256=hash_identifiers(object_ids),
        content_sha256="unused",
    )


def test_recovered_redsea_partition_requires_complete_canonical_content(
    tmp_path,
):
    path = tmp_path / "data_0.parquet"
    frame = pd.DataFrame(
        {
            "object_id": ["a", "b"],
            "donor_id": ["1", "1"],
            "cell_area_px": [10.0, 11.0],
        }
    )
    frame.to_parquet(path, index=False)

    check, schema = validate_recovered_redsea_partition(
        path,
        donor_id="1",
        expected=_expected("1", ["a", "b"]),
        expected_columns=tuple(frame.columns),
    )
    assert check.donor_id == "1"
    assert check.row_count == 2
    assert schema.field("cell_area_px").type == pa.float64()

    changed = frame.copy()
    changed.loc[1, "object_id"] = "other"
    changed.to_parquet(path, index=False)
    with pytest.raises(StaleArtifactError, match="object universe differs"):
        validate_recovered_redsea_partition(
            path,
            donor_id="1",
            expected=_expected("1", ["a", "b"]),
            expected_columns=tuple(frame.columns),
        )


def test_recovered_redsea_partition_rejects_truncated_parquet(tmp_path):
    path = tmp_path / "data_0.parquet"
    path.write_bytes(b"PAR1incomplete")
    with pytest.raises(PartialArtifactError, match="cannot open"):
        validate_recovered_redsea_partition(
            path,
            donor_id="1",
            expected=_expected("1", ["a"]),
            expected_columns=("object_id", "donor_id", "cell_area_px"),
        )


def test_manifested_geometry_is_hardlinked_idempotently(tmp_path):
    source_root = tmp_path / "source" / "geometry_qc"
    for donor in ("1", "2"):
        write_donor_partition(
            pd.DataFrame(
                {
                    "donor_id": [donor],
                    "object_id": [f"cell-{donor}"],
                    "value": [1.0],
                }
            ),
            source_root,
            donor_id=donor,
        )
    manifest = StageManifest.create(
        stage="geometry",
        method_version="test",
        expected_donors=("1", "2"),
        inputs=(),
        output_root=source_root,
        config={},
        code=CodeMetadata.capture(
            REPOSITORY,
            source_paths=("phenocycler/geometry_qc.py",),
        ),
    )
    target_root = tmp_path / "target" / "geometry_qc"
    context = SimpleNamespace(
        donors=("1", "2"),
        output_root=lambda _stage: target_root,
    )

    links = _link_manifest_dataset(
        context,
        stage="geometry",
        source=manifest,
    )
    assert len(links) == 2
    assert _link_manifest_dataset(
        context,
        stage="geometry",
        source=manifest,
    ) == links
    for donor in context.donors:
        source = source_root / f"donor_id={donor}" / "data.parquet"
        target = target_root / f"donor_id={donor}" / "data.parquet"
        assert source.stat().st_ino == target.stat().st_ino


def test_recover_cli_is_dispatched_before_the_run_parser(monkeypatch):
    from phenocycler import pipeline, recovery

    called = {}

    def fake_main(arguments):
        called["arguments"] = arguments
        return 7

    monkeypatch.setattr(recovery, "main", fake_main)
    assert (
        pipeline.main(
            [
                "recover",
                "--from-run",
                "a" * 20,
                "--prepare-only",
            ]
        )
        == 7
    )
    assert called["arguments"] == [
        "--from-run",
        "a" * 20,
        "--prepare-only",
    ]


def test_recovery_blocks_pair_prefix_before_preparing_outputs(monkeypatch):
    from phenocycler import recovery

    context = SimpleNamespace(
        registry=SimpleNamespace(
            production_status="development_no_verified_typing"
        )
    )
    monkeypatch.setattr(recovery, "load_config", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(recovery, "resolve_run_context", lambda _cfg: context)
    monkeypatch.setattr(
        recovery,
        "prepare_interrupted_run_recovery",
        lambda *_args, **_kwargs: pytest.fail("recovery mutation must not begin"),
    )

    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        recovery.recover(
            config_path=None,
            source_run_id="a" * 20,
            through="type",
        )
