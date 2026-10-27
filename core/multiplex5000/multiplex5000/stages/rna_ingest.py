"""Stage 04: raw Xenium RNA ingest and cell-level QC for accepted regions.

One recorded source per selected Stage 00 acquisition is read. Raw matrices remain in
their instrument bundles and are bound by SHA256 in this stage's receipt; the
stage publishes region-qualified cell metadata, feature lists and a QC gate.
It does not normalize RNA or call cell types. A damaged tissue designation is
carried forward and does not suppress the other acquisitions.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import gates
from ..tasks import Input, Task
from ..xenium_region import load_region
from ..xenium_sources import INSTRUMENT_FILES, XeniumSource, build_roster, inventory_table

STAGE = "04_rna_ingest"
_KEY = re.compile(r"^[A-Za-z0-9_]+$")
_RAW = (*INSTRUMENT_FILES, "cells.parquet", "cell_feature_matrix.h5")
_COLUMNS = ("acquisition_id", "sample_key", "donor_id", "tissue", "roi", "source_path",
            "condition", "cells", "genes", "typing_eligible", "no_nucleus",
            "nucleus_missing", "median_gene_counts", "median_detected_genes", "source_status")


class RNAIngestError(ValueError):
    """An accepted Xenium acquisition cannot be ingested as recorded."""


def _inventory(project) -> Path:
    return Path(project.release_dir) / "00_inventory"


def _locations(project) -> tuple[Path | None, Path | None]:
    index = project.source("xenium_source_index")
    records = project.source("xenium_records")
    return (Path(index) if index else None,
            Path(records) / "MANIFEST.tsv" if records else None)


def _selection(project) -> set[str] | None:
    chosen = getattr(project, "settings", {}).get("rna_acquisitions")
    if chosen is None:
        return None
    if (not isinstance(chosen, list) or not chosen
            or not all(isinstance(item, str) and item for item in chosen)
            or len(set(chosen)) != len(chosen)):
        raise RNAIngestError("settings.rna_acquisitions needs distinct acquisition IDs")
    return set(chosen)


def settings(project) -> dict:
    chosen = _selection(project)
    return {"rna_acquisitions": sorted(chosen) if chosen is not None else "all_accepted"}


def _ready(project) -> list[XeniumSource]:
    index, manifest = _locations(project)
    roster = build_roster(_inventory(project), manifest, site_index=index,
                          acquisition_ids=_selection(project))
    not_ready = [(r.acquisition_id, r.status, r.missing_files)
                 for r in roster if r.status != "ready" and r.condition != "damaged"]
    if not_ready:
        raise RNAIngestError(f"Stage 04 needs every selected source ready: {not_ready[:5]}")
    for source in roster:
        if not _KEY.fullmatch(source.sample_key):
            raise RNAIngestError(f"unsafe sample key: {source.sample_key}")
    return roster


def _raw_name(source: XeniumSource, filename: str) -> str:
    return f"raw:{source.acquisition_id}:{filename}"


def inputs(project) -> list[Input]:
    """Declare accepted evidence, location records and every selected raw input."""
    inventory = _inventory(project)
    found = [Input("inventory_sources", inventory_table(inventory)),
             Input("inventory_receipt", inventory / "receipt.json"),
             Input("inventory_gate", inventory / "gate" / "gate.json"),
             Input("inventory_acceptance", inventory / "gate" / "accepted.json")]
    index, manifest = _locations(project)
    if index is not None:
        found.append(Input("xenium_source_index", index))
    if manifest is not None:
        found.append(Input("xenium_records_manifest", manifest))
    for source in _ready(project):
        if source.status != "ready":
            continue  # Recorded damaged tissue has no raw files to declare.
        for filename in _RAW:
            found.append(Input(_raw_name(source, filename),
                               Path(source.source_path) / filename))
    return found


def run(task: Task) -> dict:
    index = task.input("xenium_source_index") if "xenium_source_index" in task.inputs else None
    manifest = (task.input("xenium_records_manifest")
                if "xenium_records_manifest" in task.inputs else None)
    selected = task.settings["rna_acquisitions"]
    roster = build_roster(_inventory(task.project), manifest, site_index=index,
                          acquisition_ids=set(selected) if isinstance(selected, list) else None)
    rows = []
    total_cells = total_eligible = 0
    for source in roster:
        if not _KEY.fullmatch(source.sample_key):
            raise RNAIngestError(f"{source.acquisition_id}: unsafe sample key")
        if source.status != "ready":
            if source.condition != "damaged":
                raise RNAIngestError(f"{source.acquisition_id}: indexed source is not ready")
            rows.append({"acquisition_id": source.acquisition_id,
                         "sample_key": source.sample_key, "donor_id": source.donor_id,
                         "tissue": source.tissue, "roi": source.roi,
                         "source_path": source.source_path, "condition": source.condition,
                         "source_status": source.status})
            continue
        for filename in _RAW:
            declared = task.input(_raw_name(source, filename))
            if declared != Path(source.source_path) / filename:
                raise RNAIngestError(f"{source.acquisition_id}: source changed after declaration")
        region = load_region(source)
        counts = np.asarray(region.counts.sum(axis=1)).ravel()
        detected = np.diff(region.counts.indptr)
        no_nucleus = (region.nucleus_count == 0) & ~region.nucleus_missing
        cell_table = pa.table({
            "cell_uid": pa.array(region.cell_uid),
            "cell_id": pa.array(region.cell_id),
            "acquisition_id": pa.array([source.acquisition_id] * region.n_cells),
            "donor_id": pa.array([source.donor_id] * region.n_cells),
            "tissue": pa.array([source.tissue] * region.n_cells),
            "roi": pa.array([source.roi] * region.n_cells),
            "x_um": pa.array(region.x_um), "y_um": pa.array(region.y_um),
            "nucleus_count": pa.array(region.nucleus_count, mask=region.nucleus_missing),
            "typing_eligible": pa.array(region.eligible),
            "gene_counts": pa.array(counts), "detected_genes": pa.array(detected),
        })
        pq.write_table(cell_table, task.output(f"cells/{source.sample_key}.parquet"),
                       compression="zstd")
        task.output(f"features/{source.sample_key}.txt").write_text(
            "\n".join(region.genes) + "\n")
        rows.append({
            "acquisition_id": source.acquisition_id, "sample_key": source.sample_key,
            "donor_id": source.donor_id, "tissue": source.tissue, "roi": source.roi,
            "source_path": source.source_path, "condition": source.condition,
            "cells": region.n_cells, "genes": len(region.genes),
            "typing_eligible": int(region.eligible.sum()),
            "no_nucleus": int(no_nucleus.sum()),
            "nucleus_missing": int(region.nucleus_missing.sum()),
            "median_gene_counts": float(np.median(counts)) if len(counts) else 0.0,
            "median_detected_genes": float(np.median(detected)) if len(detected) else 0.0,
            "source_status": "ready",
        })
        total_cells += region.n_cells
        total_eligible += int(region.eligible.sum())
        del region, cell_table
    with task.output("regions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    with task.input("inventory_sources").open(newline="") as handle:
        accepted_regions = sum(row["platform"] == "xenium" for row in csv.DictReader(handle))
    _figure(task, rows, accepted_regions)
    return {"regions": len(rows), "accepted_regions": accepted_regions,
            "scope": "selected" if selected != "all_accepted" else "all_accepted",
            "cells": total_cells, "typing_eligible": total_eligible,
            "by_tissue": {tissue: sum(r["cells"] for r in rows if r["tissue"] == tissue
                                     and "cells" in r)
                          for tissue in sorted({r["tissue"] for r in rows if "cells" in r})},
            "no_nucleus": sum(r.get("no_nucleus", 0) for r in rows),
            "nucleus_missing": sum(r.get("nucleus_missing", 0) for r in rows),
            "damaged_regions": [r["acquisition_id"] for r in rows if r["condition"]],
            "unavailable_regions": [r["acquisition_id"] for r in rows
                                    if r["source_status"] != "ready"]}


def _figure(task: Task, rows: list[dict], accepted_regions: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ordered = sorted(rows, key=lambda row: (row["donor_id"], row["tissue"], row["roi"]))
    n_regions = len(ordered)
    height = max(9, 5 + 0.34 * n_regions)
    fig = plt.figure(figsize=(18, height), layout="constrained")
    grid = fig.add_gridspec(2, 6, height_ratios=[max(3, 0.34 * n_regions), 3])
    admission = fig.add_subplot(grid[0, :2])
    genes = fig.add_subplot(grid[0, 2:4], sharey=admission)
    transcripts = fig.add_subplot(grid[0, 4:], sharey=admission)
    gene_cdf = fig.add_subplot(grid[1, :3])
    transcript_cdf = fig.add_subplot(grid[1, 3:])
    tissue_colors = {"pancreas": "#287ba8", "pancreatic_lymph_node": "#ce6b27"}
    labels = [f"{row['donor_id']} {row['roi']} "
              f"({'pLN' if row['tissue'] == 'pancreatic_lymph_node' else 'pancreas'})"
              for row in ordered]
    y = np.arange(n_regions)
    eligible = np.array([row.get("typing_eligible", 0) for row in ordered])
    no_nucleus = np.array([row.get("no_nucleus", 0) for row in ordered])
    missing = np.array([row.get("nucleus_missing", 0) for row in ordered])
    admission.barh(y, eligible, color="#21867a", label="eligible, nucleus measured")
    admission.barh(y, no_nucleus, left=eligible, color="#e0a047", label="no nucleus")
    admission.barh(y, missing, left=eligible + no_nucleus, color="#9b72ad",
                   label="nucleus count missing")
    admission.set(yticks=y, yticklabels=labels, xlabel="Segmented cells (n)",
                  title="Cell admission by acquisition")
    admission.invert_yaxis()
    admission.legend(loc="lower right", fontsize=8)
    for axis, title, xlabel in (
            (genes, "Detected genes in eligible cells", "Detected genes per cell (n)"),
            (transcripts, "Raw RNA in eligible cells", "Transcripts per cell (n)")):
        axis.set(title=title, xlabel=xlabel)
        axis.tick_params(axis="y", left=False, labelleft=False)
        axis.grid(axis="x", color="#dddddd", linewidth=0.6)
    for pos, row in enumerate(ordered):
        if "cells" not in row:
            admission.text(0, pos, "source unavailable; cell count unknown", va="center",
                           fontsize=8)
            continue
        table = pq.read_table(task.staging / "cells" / f"{row['sample_key']}.parquet",
                              columns=["typing_eligible", "gene_counts", "detected_genes"])
        valid = table["typing_eligible"].to_numpy(zero_copy_only=False).astype(bool)
        color = tissue_colors.get(row["tissue"], "#555555")
        for column, interval_axis, cdf_axis in (
                ("detected_genes", genes, gene_cdf),
                ("gene_counts", transcripts, transcript_cdf)):
            values = table[column].to_numpy(zero_copy_only=False)[valid]
            if not len(values):
                continue
            p10, p25, median, p75, p90 = np.quantile(values, [.1, .25, .5, .75, .9])
            interval_axis.hlines(pos, p10, p90, color=color, linewidth=2)
            interval_axis.hlines(pos, p25, p75, color=color, linewidth=7)
            interval_axis.plot(median, pos, marker="|", markersize=13, color="black")
            quantiles = np.quantile(values, np.linspace(0, 1, 301))
            cdf_axis.plot(quantiles, np.linspace(0, 1, 301), color=color,
                          alpha=0.9 if n_regions <= 8 else 0.35,
                          label=f"{row['donor_id']} {row['roi']}: n={len(values):,}, "
                                f"median={median:,.0f}" if n_regions <= 8 else None)
    for axis, xlabel in ((gene_cdf, "Detected genes per eligible cell (n)"),
                         (transcript_cdf, "Raw transcripts per eligible cell (n)")):
        axis.set(xlabel=xlabel, ylabel="Cumulative fraction of eligible cells",
                 ylim=(0, 1), xlim=(0, None))
        axis.grid(color="#dddddd", linewidth=0.6)
        if n_regions <= 8:
            axis.legend(loc="lower right", fontsize=8)
        else:
            axis.legend(handles=[Line2D([0], [0], color=color, label=tissue)
                                 for tissue, color in tissue_colors.items()
                                 if tissue in {row["tissue"] for row in ordered}],
                        loc="lower right", fontsize=8)
    fig.suptitle(f"Xenium RNA ingest QC: {sum('cells' in row for row in rows)}/"
                 f"{n_regions} selected regions measured ({accepted_regions} accepted); "
                 f"{sum(row.get('cells', 0) for row in rows):,} measured cells\n"
                 "Intervals: 10th–90th percentile, thick bar: 25th–75th, "
                 "black tick: median; curves use all eligible cells", fontsize=14)
    fig.savefig(gates.figure_path(task, "rna_ingest_qc.png"), dpi=160)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    expected = len(_ready(task.project))
    checks = [
        gates.Check("every selected Xenium region accounted", summary["regions"],
                    f"{expected}", summary["regions"] == expected),
        gates.Check("Xenium region scope", f"{summary['regions']}/{summary['accepted_regions']}",
                    "reported; a selected pilot does not claim full-cohort coverage", True,
                    blocking=False, detail=summary["scope"]),
        gates.Check("cell accounting", summary["cells"],
                    "typing-eligible plus no-nucleus plus missing-nucleus",
                    summary["cells"] == (summary["typing_eligible"]
                                         + summary["no_nucleus"]
                                         + summary["nucleus_missing"])),
        gates.Check("damage carried forward", summary["damaged_regions"],
                    "reported; a damaged region does not stop other regions", True,
                    blocking=False),
        gates.Check("damaged sources unavailable", summary["unavailable_regions"],
                    "reported with unknown cell counts, never zero-filled", True,
                    blocking=False),
    ]
    return gates.Gate(STAGE, "Raw Xenium counts and cell QC by accepted acquisition; "
                      "no RNA typing or gene-positive cutoff.", checks,
                      ["rna_ingest_qc.png"])
