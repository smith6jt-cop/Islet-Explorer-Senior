"""Bind an RNA typing run to one completed Stage 04 Xenium acquisition.

The same donor can have several sections, including same-slide and serial pLN.
Typing must use the exact instrument bundle admitted by Stage 04, rather than
searching by donor for a convenient ``cells.parquet`` or h5ad file.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from . import gates
from .artifacts import sha256_file
from .xenium_region import XeniumRegion, load_region
from .xenium_sources import XeniumSource


class XeniumTypingInputError(ValueError):
    """The recorded Stage 04 output and its raw acquisition do not agree."""


def source_from_stage04(row: dict[str, str]) -> XeniumSource:
    """Use the recorded source path and identity, without donor-level discovery."""
    required = ("acquisition_id", "sample_key", "donor_id", "roi", "tissue",
                "source_path", "condition")
    missing = [name for name in required if name not in row]
    if missing:
        raise XeniumTypingInputError(f"Stage 04 region row lacks {missing}")
    folder = Path(row["source_path"])
    if not folder.is_absolute() or not all(row[name] for name in required[:6]):
        raise XeniumTypingInputError("Stage 04 region needs complete identity and absolute source")
    return XeniumSource(row["acquisition_id"], row["sample_key"], row["donor_id"],
                        row["roi"], row["tissue"], folder.name, str(folder), str(folder),
                        "ready", (), row["condition"])


def _recorded_digest(records: list[dict], name: str, path: Path) -> str:
    matches = [item for item in records if item.get("name", item.get("path")) == name]
    if len(matches) != 1:
        raise XeniumTypingInputError(f"Stage 04 receipt needs one {name!r} record")
    record = matches[0]
    if "name" in record and Path(record.get("path", "")) != path:
        raise XeniumTypingInputError(f"Stage 04 {name} points to another acquisition")
    digest = record.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64:
        raise XeniumTypingInputError(f"Stage 04 {name} has no whole-file SHA256")
    if sha256_file(path) != digest:
        raise XeniumTypingInputError(f"Stage 04 {name} changed after ingest")
    return digest


def load_accepted_region(row: dict[str, str], receipt_path: Path, regions_path: Path,
                         stage04_cells: Path, raw_matrix: Path, raw_cells: Path) -> XeniumRegion:
    """Replay the exact Stage 04 cell join and verify its receipt before typing."""
    source = source_from_stage04(row)
    folder = Path(source.source_path)
    if raw_matrix != folder / "cell_feature_matrix.h5" or raw_cells != folder / "cells.parquet":
        raise XeniumTypingInputError(f"{source.acquisition_id}: raw inputs point elsewhere")
    stage04_dir = receipt_path.parent
    if (regions_path != stage04_dir / "regions.csv"
            or stage04_cells != stage04_dir / "cells" / f"{source.sample_key}.parquet"):
        raise XeniumTypingInputError(f"{source.acquisition_id}: Stage 04 outputs point elsewhere")
    try:
        receipt = json.loads(receipt_path.read_text())
    except (OSError, ValueError) as err:
        raise XeniumTypingInputError("Stage 04 receipt is unreadable") from err
    if receipt.get("stage") != "04_rna_ingest" or receipt.get("status") != "complete":
        raise XeniumTypingInputError("Stage 04 has no completed RNA ingest receipt")
    acceptance = gates.require_accepted(receipt_path.parent)
    if acceptance.get("receipt_sha256") != sha256_file(receipt_path):
        raise XeniumTypingInputError("Stage 04 acceptance does not bind this receipt")
    inputs = receipt.get("inputs", [])
    outputs = receipt.get("outputs", [])
    for filename, path in (("cell_feature_matrix.h5", raw_matrix), ("cells.parquet", raw_cells)):
        _recorded_digest(inputs, f"raw:{source.acquisition_id}:{filename}", path)
    _recorded_digest(outputs, "regions.csv", regions_path)
    with regions_path.open(newline="") as handle:
        recorded = [item for item in csv.DictReader(handle)
                    if item.get("acquisition_id") == source.acquisition_id]
    if len(recorded) != 1 or recorded[0] != row:
        raise XeniumTypingInputError(f"{source.acquisition_id}: Stage 04 row changed")
    _recorded_digest(outputs, f"cells/{source.sample_key}.parquet", stage04_cells)
    region = load_region(source)
    table = pq.read_table(stage04_cells, columns=["cell_uid", "cell_id", "acquisition_id",
                                                   "donor_id", "tissue", "roi", "typing_eligible",
                                                   "gene_counts", "detected_genes"])
    data = table.to_pydict()
    depth = np.asarray(region.counts.sum(axis=1)).ravel()
    detected = np.diff(region.counts.indptr)
    expected = {"cell_uid": region.cell_uid.tolist(), "cell_id": region.cell_id.tolist(),
                "acquisition_id": [source.acquisition_id] * region.n_cells,
                "donor_id": [source.donor_id] * region.n_cells,
                "tissue": [source.tissue] * region.n_cells,
                "roi": [source.roi] * region.n_cells,
                "typing_eligible": region.eligible.tolist(),
                "gene_counts": depth.tolist(), "detected_genes": detected.tolist()}
    if data != expected:
        raise XeniumTypingInputError(f"{source.acquisition_id}: raw cells differ from Stage 04")
    return region


def write_worker_h5ad(region: XeniumRegion, destination: Path, *, sex: str,
                      disease_status: str, batch: str) -> None:
    """Give the current RNA workers their raw-count h5ad, preserving exact cell IDs.

    This is a disposable worker input, never a second published data generation.
    All metadata are required from recorded acquisition and donor sources.
    """
    if not all(isinstance(value, str) and value.strip()
               for value in (sex, disease_status, batch)):
        raise XeniumTypingInputError("sex, disease status and batch need recorded values")
    import anndata as ad
    import pandas as pd

    depth = np.asarray(region.counts.sum(axis=1)).ravel()
    detected = np.diff(region.counts.indptr)
    obs = pd.DataFrame({"nCount_Xenium": depth, "nFeature_Xenium": detected,
                        "sex": sex, "disease.status": disease_status,
                        "batch": batch}, index=pd.Index(region.cell_id, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(region.genes, name="gene"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=region.counts.copy(), obs=obs, var=var).write_h5ad(destination)
