"""Exact sample identity for the spleen-follicle workflow.

Identity is curated in ``config/sample_manifest.csv``.  Filenames may be
aliases, but their prefixes are never treated as donor identifiers.
"""

from __future__ import annotations

import csv
from pathlib import Path


CODE_ROOT = Path(__file__).resolve().parent
SAMPLE_MANIFEST = CODE_ROOT / "config" / "sample_manifest.csv"
REQUIRED_COLUMNS = {"alias", "sample_id", "alias_type", "evidence"}


class SampleManifestError(ValueError):
    """An image or metadata alias is absent or ambiguous in the manifest."""


def load_sample_manifest(path: Path = SAMPLE_MANIFEST) -> list[dict[str, str]]:
    """Load and validate the exact alias table."""

    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        missing = sorted(REQUIRED_COLUMNS - fields)
        if missing:
            raise SampleManifestError(f"{path} lacks columns: {missing}")
        rows = [
            {
                column: str(row.get(column, "")).strip()
                for column in reader.fieldnames or ()
            }
            for row in reader
        ]

    seen: dict[str, str] = {}
    duplicates: set[str] = set()
    conflicts: set[str] = set()
    for row in rows:
        alias = row["alias"]
        sample_id = row["sample_id"]
        if not alias or not sample_id:
            raise SampleManifestError(f"{path} contains empty alias/sample rows")
        if alias in seen:
            duplicates.add(alias)
            if seen[alias] != sample_id:
                conflicts.add(alias)
        seen[alias] = sample_id
    if conflicts:
        raise SampleManifestError(
            f"{path} maps aliases to multiple samples: {sorted(conflicts)}"
        )
    if duplicates:
        raise SampleManifestError(f"{path} repeats aliases: {sorted(duplicates)}")
    return rows


SAMPLE_MANIFEST_ROWS = load_sample_manifest()
SAMPLE_BY_ALIAS = {row["alias"]: row["sample_id"] for row in SAMPLE_MANIFEST_ROWS}


def resolve_sample_id(alias: str, mapping: dict[str, str] = SAMPLE_BY_ALIAS) -> str:
    """Resolve one exact recorded alias, failing closed when it is unknown."""

    key = str(alias).strip()
    try:
        return mapping[key]
    except KeyError as exc:
        raise SampleManifestError(
            f"unlisted sample alias {key!r}; add it to {SAMPLE_MANIFEST}"
        ) from exc
