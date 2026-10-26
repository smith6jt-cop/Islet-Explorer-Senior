"""Resolve Xenium region outputs from an accepted Stage 00 inventory.

The inventory is the roster, including Xenium-only and same-slide sections. The
records manifest names exact HiPerGator output directories for bundles that are
only partly present on the workstation. Neither a donor's older section nor
another region of the same Xenium slide is a substitute for a missing bundle.

This is a source readiness check for Stage 04, not an RNA typing producer.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from . import gates
from .integration.manifest import parse_bundle_path, sample_key
from .tasks import sha256_file

INSTRUMENT_FILES = ("experiment.xenium", "gene_panel.json", "metrics_summary.csv")
CELL_FILES = ("cells.parquet", "cell_feature_matrix.h5", "cell_boundaries.parquet",
              "nucleus_boundaries.parquet")


class XeniumSourceError(ValueError):
    """The accepted roster or a recorded source is inconsistent."""


@dataclass(frozen=True)
class XeniumSource:
    acquisition_id: str
    sample_key: str
    donor_id: str
    roi: str
    tissue: str
    output: str
    inventory_path: str
    source_path: str
    status: str                 # ready, remote_only, incomplete, missing
    missing_files: tuple[str, ...]
    condition: str             # acquisition-specific damage, if recorded


def _rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")),
                                   delimiter=delimiter))


def _accepted_inventory(stage_dir: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Bind the acquisition table and instrument hashes to the accepted receipt."""
    stage_dir = Path(stage_dir)
    acceptance = gates.require_accepted(stage_dir)
    receipt_path = stage_dir / "receipt.json"
    if acceptance.get("receipt_sha256") != sha256_file(receipt_path):
        raise XeniumSourceError("Stage 00 receipt changed after gate acceptance")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("stage") != "00_inventory" or receipt.get("status") != "complete":
        raise XeniumSourceError("source roster requires a complete Stage 00 receipt")
    output = next((item for item in receipt["outputs"]
                   if item["path"] == "acquisitions.csv"), None)
    table = stage_dir / "acquisitions.csv"
    if output is None or output["sha256"] != sha256_file(table):
        raise XeniumSourceError("Stage 00 acquisitions.csv differs from its receipt")
    inputs = {item["name"]: item.get("sha256") for item in receipt["inputs"]}
    return _rows(table), inputs


def _remote_sources(path: Path | None, input_hashes: dict[str, str]) -> dict[str, str]:
    if path is None:
        return {}
    path = Path(path)
    if input_hashes.get("xenium_records_manifest") != sha256_file(path):
        raise XeniumSourceError("Xenium records manifest differs from accepted Stage 00 input")
    sources: dict[str, str] = {}
    files: dict[str, set[str]] = {}
    for row in _rows(path, "\t"):
        if row["file"] not in INSTRUMENT_FILES:
            continue
        output = row["output"]
        seen = files.setdefault(output, set())
        if row["file"] in seen:
            raise XeniumSourceError(f"{output}/{row['file']}: repeated manifest record")
        seen.add(row["file"])
        expected = input_hashes.get(f"xenium:{output}:{row['file']}")
        if expected != row["sha256"]:
            raise XeniumSourceError(f"{output}/{row['file']} differs from Stage 00 input")
        source = row["source"]
        if not source.startswith("hpg:/") or Path(source[4:]).name != row["file"]:
            raise XeniumSourceError(f"{output}/{row['file']}: invalid recorded HPG source")
        directory = str(Path(source[4:]).parent)
        if Path(directory).name != output:
            raise XeniumSourceError(f"{output}/{row['file']}: source names another output")
        if output in sources and sources[output] != directory:
            raise XeniumSourceError(f"{output}: instrument records name different directories")
        sources[output] = directory
    for output, seen in files.items():
        if seen != set(INSTRUMENT_FILES):
            raise XeniumSourceError(f"{output}: incomplete instrument-record source set")
    return sources


def _missing_files(directory: Path) -> tuple[str, ...]:
    return tuple(name for name in (*INSTRUMENT_FILES, *CELL_FILES)
                 if not (directory / name).is_file() or (directory / name).stat().st_size == 0)


def build_roster(stage_dir: Path, records_manifest: Path | None = None) -> list[XeniumSource]:
    """One row per Xenium region, with an exact source and its present-host readiness.

    A ``remote_only`` row has a recorded HPG directory but no complete bundle on this
    host. It is deliberately not ``ready`` until its cell files are checked there.
    """
    acquisitions, hashes = _accepted_inventory(stage_dir)
    remote = _remote_sources(records_manifest, hashes)
    roster: list[XeniumSource] = []
    keys: set[str] = set()
    outputs: set[str] = set()
    for row in acquisitions:
        if row["platform"] != "xenium":
            continue
        output = row["output"]
        serial, _ = parse_bundle_path(output)
        key = sample_key(serial, row["roi"])
        if not serial or not output or row["path"].rstrip("/").split("/")[-1] != output:
            raise XeniumSourceError(f"{row['acquisition_id']}: invalid output identity")
        if key in keys or output in outputs:
            raise XeniumSourceError(f"{row['acquisition_id']}: repeated region or output")
        keys.add(key)
        outputs.add(output)

        local = Path(row["path"])
        exact_remote = Path(remote[output]) if output in remote else None
        candidates = [local] + ([exact_remote] if exact_remote and exact_remote != local else [])
        ready = next((path for path in candidates if not _missing_files(path)), None)
        if ready is not None:
            for name in INSTRUMENT_FILES:
                expected = hashes.get(f"xenium:{output}:{name}")
                if expected is None or sha256_file(ready / name) != expected:
                    raise XeniumSourceError(f"{output}/{name} differs from accepted Stage 00")
            path, status, missing = ready, "ready", ()
        elif exact_remote is not None and not exact_remote.is_dir():
            path, status, missing = exact_remote, "remote_only", _missing_files(exact_remote)
        else:
            path = exact_remote if exact_remote and exact_remote.is_dir() else local
            missing = _missing_files(path)
            status = "incomplete" if path.is_dir() else "missing"
        roster.append(XeniumSource(
            acquisition_id=row["acquisition_id"], sample_key=key,
            donor_id=row["donor_id"], roi=row["roi"], tissue=row["tissue"],
            output=output, inventory_path=str(local), source_path=str(path),
            status=status, missing_files=missing, condition=row.get("condition", "")))
    return roster


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-dir", type=Path, required=True)
    parser.add_argument("--records-manifest", type=Path)
    args = parser.parse_args()
    print(json.dumps([asdict(row) for row in build_roster(args.inventory_dir,
                                                        args.records_manifest)], indent=2))


if __name__ == "__main__":
    main()
