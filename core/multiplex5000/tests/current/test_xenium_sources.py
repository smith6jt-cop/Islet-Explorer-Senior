"""Stage 04 must start from every accepted Xenium region, including remote ones."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from multiplex5000 import gates
from multiplex5000.artifacts import sha256_file
from multiplex5000.xenium_sources import (
    CELL_FILES,
    INSTRUMENT_FILES,
    XeniumSourceError,
    build_roster,
)


def _bundle(root: Path, serial: str, label: str) -> Path:
    path = root / f"output-XETG00298__{serial}__{label}__20261026__120000"
    path.mkdir(parents=True)
    for name in (*INSTRUMENT_FILES, *CELL_FILES):
        (path / name).write_bytes(f"{serial} {label} {name}".encode())
    return path


def _accepted(tmp_path: Path, records: list[dict], remote: dict[str, Path] | None = None):
    stage = tmp_path / "00_inventory"
    stage.mkdir()
    table = stage / "acquisitions.csv"
    columns = ("platform", "acquisition_id", "donor_id", "roi", "tissue", "output", "path",
               "condition")
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    inputs = []
    manifest = tmp_path / "MANIFEST.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("output", "file", "source", "sha256", "bytes",
                                                     "copied_utc"), delimiter="\t")
        writer.writeheader()
        for record in records:
            output = record["output"]
            for name in INSTRUMENT_FILES:
                digest = sha256_file(Path(record["record_bundle"]) / name)
                inputs.append({"name": f"xenium:{output}:{name}", "sha256": digest})
                if remote and output in remote:
                    writer.writerow({"output": output, "file": name,
                                     "source": f"hpg:{remote[output] / name}",
                                     "sha256": digest, "bytes": 10, "copied_utc": ""})
    inputs.append({"name": "xenium_records_manifest", "sha256": sha256_file(manifest)})
    receipt = {"stage": "00_inventory", "status": "complete", "inputs": inputs,
               "outputs": [{"path": "acquisitions.csv", "sha256": sha256_file(table)}]}
    (stage / "receipt.json").write_text(json.dumps(receipt))
    gate = stage / "gate"
    gate.mkdir()
    (gate / "gate.json").write_text(json.dumps({"stage": "00_inventory", "release": "test",
                                                "checks": [], "passed": True}))
    gates.accept(stage, by="test")
    return stage, manifest


def _record(bundle: Path, donor: str, roi: str, *, record_bundle: Path | None = None,
            condition: str = "") -> dict:
    serial = bundle.name.split("__")[1]
    return {"platform": "xenium", "acquisition_id": f"xenium:{serial}:{roi}",
            "donor_id": donor, "roi": roi,
            "tissue": "pancreas" if roi == "panc" else "pancreatic_lymph_node",
            "output": bundle.name, "path": str(bundle), "condition": condition,
            "record_bundle": str(record_bundle or bundle)}


def _site_index(path: Path, records: list[dict], sources: list[Path]) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("acquisition_id", "output", "source_path"))
        writer.writeheader()
        for record, source in zip(records, sources, strict=True):
            writer.writerow({"acquisition_id": record["acquisition_id"],
                             "output": record["output"], "source_path": str(source)})
    return path


def test_same_slide_regions_and_xenium_only_are_separate(tmp_path):
    local = _bundle(tmp_path / "local", "0030001", "panc")
    new_panc = _bundle(tmp_path / "copies", "0040001", "panc")
    new_pln = _bundle(tmp_path / "copies", "0040001", "pLN")
    xenium_only = _bundle(tmp_path / "copies", "0040002", "Region_1")
    rows = [_record(local, "6533", "panc"),
            _record(tmp_path / "absent" / new_panc.name, "6533", "panc",
                    record_bundle=new_panc),
            _record(tmp_path / "absent" / new_pln.name, "6533", "pln_2",
                    record_bundle=new_pln, condition="damaged"),
            _record(tmp_path / "absent" / xenium_only.name, "6512", "panc",
                    record_bundle=xenium_only)]
    remote = {r["output"]: Path("/orange/maigan") / r["output"] for r in rows[1:]}
    stage, manifest = _accepted(tmp_path, rows, remote)
    roster = build_roster(stage, manifest)

    assert [(r.sample_key, r.status) for r in roster] == [
        ("0030001__panc", "ready"), ("0040001__panc", "remote_only"),
        ("0040001__pln_2", "remote_only"), ("0040002__panc", "remote_only")]
    assert roster[1].source_path == str(remote[new_panc.name])
    assert roster[2].condition == "damaged"
    assert roster[3].donor_id == "6512"


def test_modified_instrument_record_is_rejected(tmp_path):
    local = _bundle(tmp_path / "local", "0040001", "panc")
    stage, manifest = _accepted(tmp_path, [_record(local, "6512", "panc")])
    (local / "experiment.xenium").write_text("changed")
    with pytest.raises(XeniumSourceError, match="differs from accepted Stage 00"):
        build_roster(stage, manifest)


def test_exact_remote_bundle_is_ready_on_the_host_that_holds_it(tmp_path):
    local = tmp_path / "old_location" / "output-XETG00298__0040001__panc__20261026__120000"
    remote = _bundle(tmp_path / "orange" / "maigan", "0040001", "panc")
    stage, manifest = _accepted(tmp_path, [_record(local, "6512", "panc",
                                                   record_bundle=remote)],
                                {remote.name: remote})
    roster = build_roster(stage, manifest)
    assert roster[0].status == "ready"
    assert roster[0].source_path == str(remote)


def test_site_index_selects_a_renamed_region_and_never_falls_back(tmp_path):
    local = _bundle(tmp_path / "local", "0040001", "panc")
    source = _bundle(tmp_path / "site", "0040001", "PanT")
    recorded = tmp_path / "old_location" / source.name.replace("PanT", "pLN")
    rows = [_record(local, "6533", "panc"),
            _record(recorded, "6533", "pln_2", record_bundle=source)]
    stage, _ = _accepted(tmp_path, rows)
    index = _site_index(tmp_path / "site.csv", rows, [local, source])
    roster = build_roster(stage, site_index=index)
    assert [r.status for r in roster] == ["ready", "ready"]
    assert roster[1].sample_key == "0040001__pln_2"
    assert roster[1].source_path == str(source)

    index = _site_index(index, rows, [local, local])
    with pytest.raises(XeniumSourceError, match="assigned to another region"):
        build_roster(stage, site_index=index)

    # The index is the run's source authority; the complete local copy cannot hide its loss.
    index = _site_index(index, rows, [tmp_path / "missing" / local.name, source])
    roster = build_roster(stage, site_index=index)
    assert roster[0].status == "missing"


def test_site_index_requires_every_accepted_region_and_correct_serial(tmp_path):
    first = _bundle(tmp_path / "local", "0040001", "panc")
    second = _bundle(tmp_path / "local", "0040002", "panc")
    rows = [_record(first, "6512", "panc"), _record(second, "6618", "panc")]
    stage, _ = _accepted(tmp_path, rows)
    index = _site_index(tmp_path / "site.csv", rows[:1], [first])
    with pytest.raises(XeniumSourceError, match="site index is missing"):
        build_roster(stage, site_index=index)
    wrong_slide = _bundle(tmp_path / "other", "0040001", "Region_2")
    index = _site_index(index, rows, [first, wrong_slide])
    with pytest.raises(XeniumSourceError, match="another Xenium slide"):
        build_roster(stage, site_index=index)


def test_tracked_hpg_index_covers_the_tracked_xenium_outputs():
    repo = Path(__file__).resolve().parents[4]
    catalog = repo / "core/phenocycler/data/integration/xenium_paths.csv"
    index = repo / "core/multiplex5000/sites/hipergator/xenium_sources_20261026.csv"
    catalog_outputs = {Path(row["source_path"].strip().rstrip("/")).name
                       for row in csv.DictReader(catalog.open())}
    indexed = list(csv.DictReader(index.open()))
    assert len(indexed) == 65
    assert len({row["acquisition_id"] for row in indexed}) == 65
    assert {row["output"] for row in indexed} == catalog_outputs
    assert all(row["source_path"].startswith("/orange/maigan/") for row in indexed)
    assert all(row["output"].split("__")[1] == Path(row["source_path"]).name.split("__")[1]
               for row in indexed)


def test_accepted_inventory_cannot_be_rewritten(tmp_path):
    local = _bundle(tmp_path / "local", "0040001", "panc")
    stage, manifest = _accepted(tmp_path, [_record(local, "6512", "panc")])
    with (stage / "acquisitions.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(XeniumSourceError, match="acquisitions.csv differs"):
        build_roster(stage, manifest)


def test_incomplete_bundle_does_not_become_ready(tmp_path):
    local = _bundle(tmp_path / "local", "0040001", "panc")
    stage, manifest = _accepted(tmp_path, [_record(local, "6512", "panc")])
    (local / "cell_feature_matrix.h5").unlink()
    roster = build_roster(stage, manifest)
    assert roster[0].status == "incomplete"
    assert roster[0].missing_files == ("cell_feature_matrix.h5",)
