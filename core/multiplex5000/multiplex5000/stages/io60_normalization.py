"""IO60 Stage 05: same-block UniFORM shifts from the accepted Stage 03 cells.

The historical October reference-donor and marker-usability records are explicit
migration inputs. They determine the reference unit and the old conservative
refusal, while every histogram is recomputed from this release's measurements.
The output is a candidate for Stage 06 review, not a biological validation.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from .. import gates, normalization_block, normalize, profiles, roots
from ..tasks import Input, Task

STAGE = "05_normalization"
MODE_SD_BAR = 0.10
ACQ_RATIO_FLAG = 2.0
_META = ("donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
         "cell_region", "islet_num")
_SCALE_COLUMNS = ("donor_id", "marker", "reference_donor", "scale", "shift_log",
                  "n_fit_cells", "status", "acq_ratio", "mode_sd_before", "mode_sd_after")


class IO60NormalizationError(ValueError):
    """The declared IO60 normalization inputs cannot be combined safely."""


def _spec(project) -> Path:
    cohort = profiles.load_all()[("cohort", project.cohort)]
    return project.source("cohort_spec") or profiles.REPO / cohort.membership


def _images(project) -> list[dict]:
    images = json.loads(_spec(project).read_text())["images"]
    donors = [str(row["donor_id"]) for row in images]
    if len(donors) != len(set(donors)):
        raise IO60NormalizationError("IO60 cohort has more than one image per donor")
    return images


def _source(project, name: str) -> Path:
    path = project.source(name)
    if path is None:
        raise IO60NormalizationError(
            f"IO60 Stage 05 needs source {name!r} in the project or site file")
    return Path(path)


def _uniform(project) -> Path:
    return project.source("uniform_checkout") or roots.REPO / "external/UniFORM"


def _expression(project, donor: str) -> Path:
    return (project.release_dir / "03_protein_measurement/data/expression"
            / f"donor_id={donor}" / "data.parquet")


def _annotation(project, image: str) -> Path:
    return (_source(project, "pancreas_annotations")
            / f"{normalization_block.image_tag(image)}.geojson")


def settings(project) -> dict:
    return {"donors": [str(row["donor_id"]) for row in _images(project)],
            "uniform_commit": normalize.UNIFORM_COMMIT,
            "uniform_checkout": str(_uniform(project)),
            "mode_sd_bar": MODE_SD_BAR, "acq_ratio_flag": ACQ_RATIO_FLAG,
            "num_bins": normalize.NUM_BINS}


def inputs(project) -> list[Input]:
    stage = project.release_dir / "03_protein_measurement"
    found = [Input("protein_receipt", stage / "receipt.json"),
             Input("protein_acceptance", stage / "gate/accepted.json"),
             Input("cohort_spec", _spec(project)),
             Input("pln_blocks", _source(project, "pln_blocks")),
             Input("reference_record", _source(project, "normalization_references")),
             Input("usability_record", _source(project, "normalization_usability"))]
    for row in _images(project):
        donor, image = str(row["donor_id"]), row["image_id"]
        found.extend((Input(f"expression:{donor}", _expression(project, donor)),
                      Input(f"annotations:{donor}", _annotation(project, image),
                            required=False)))
    package = _uniform(project) / "UniFORM"
    if not package.is_dir():
        raise IO60NormalizationError(f"pinned UniFORM package missing at {package}")
    found.extend(Input(f"uniform_code:{path.name}", path)
                 for path in sorted(package.glob("*.py")))
    return found


def _references(path: Path, donors: list[str]) -> tuple[list[str], dict[str, str]]:
    record = json.loads(path.read_text())
    markers = record.get("markers")
    shifts = record.get("shifts")
    if (not isinstance(markers, list) or len(markers) != len(set(markers))
            or not isinstance(shifts, dict) or set(markers) != set(shifts)):
        raise IO60NormalizationError("reference record has no unique marker roster")
    if record.get("donors") != donors:
        raise IO60NormalizationError("reference record donor roster differs from Stage 03")
    references = {marker: shifts[marker]["reference_donor"] for marker in markers}
    unknown = {marker: ref for marker, ref in references.items() if ref not in donors}
    if unknown:
        raise IO60NormalizationError(f"reference donors outside the cohort: {unknown}")
    return markers, references


def _ratios(path: Path, markers: list[str], donors: list[str]) -> dict[str, float | None]:
    table = pq.read_table(path, columns=["donor_id", "marker", "acq_ratio"]).to_pandas()
    if not set(table["donor_id"].astype(str)) <= set(donors):
        raise IO60NormalizationError("usability record names a donor outside the cohort")
    out: dict[str, float | None] = {}
    for marker, group in table.groupby("marker"):
        ratios = group["acq_ratio"].dropna().unique()
        if len(ratios) > 1:
            raise IO60NormalizationError(f"{marker}: conflicting acquisition ratios")
        out[str(marker)] = float(ratios[0]) if len(ratios) else None
    if set(markers) != set(out):
        raise IO60NormalizationError("usability markers differ from the reference record")
    return out


def _population(task: Task, rows: list[dict], blocks: dict) -> dict[str, int]:
    declared = {row["image_id"] for row in rows}
    if set(blocks) - declared:
        raise IO60NormalizationError("pLN blocks name an image outside the cohort")
    counts = {"same-block": 0, "different-block": 0, "unassigned": 0}
    root = task.staging / "alignment_block"
    root.mkdir(parents=True, exist_ok=True)
    for row in rows:
        donor, image = str(row["donor_id"]), row["image_id"]
        side = normalization_block.build_donor(
            donor, image, task.input(f"expression:{donor}"),
            task.input(f"annotations:{donor}"), blocks.get(image, {}))
        if side is None:
            continue
        side.to_parquet(root / f"donor_id={donor}.parquet", index=False)
        for kind in counts:
            counts[kind] += int((side["block"] == kind).sum())
    return counts


def _mask(path: Path, donor: str, sidecars: Path) -> np.ndarray:
    frame = pq.ParquetFile(path).read(
        columns=["object_id", "qc_analysis_eligible", "cell_region"]).to_pandas()
    block = normalization_block.population(
        donor, frame["object_id"], frame["cell_region"], root=sidecars)
    return frame["qc_analysis_eligible"].to_numpy(bool) & block


def _values(paths: dict[str, Path], masks: dict[str, np.ndarray], marker: str) -> dict[str, np.ndarray]:
    values = {}
    for donor, path in paths.items():
        column = pq.ParquetFile(path).read(columns=[marker]).to_pandas()[marker].to_numpy("float64")
        sample = column[masks[donor]]
        if np.any(np.isfinite(sample) & (sample >= 1.0)):
            values[donor] = sample
    return values


def _aligned(task: Task, paths: dict[str, Path], markers: list[str],
             scales: dict[str, dict[str, float | None]], sidecars: Path) -> dict[str, int]:
    counts = {}
    for donor, path in paths.items():
        frame = pq.ParquetFile(path).read(columns=[*_META, *markers]).to_pandas()
        block = normalization_block.population(
            donor, frame["object_id"], frame["cell_region"], root=sidecars)
        for marker in markers:
            scale = scales[marker].get(donor)
            frame[marker] = (frame[marker].to_numpy("float64") * scale
                             if scale is not None else np.nan)
            frame.loc[~block, marker] = np.nan
        frame["in_alignment_block"] = block
        counts[donor] = int(block.sum())
        frame.to_parquet(task.output(f"data/aligned/donor_id={donor}/data.parquet"),
                         index=False, compression="zstd")
    return counts


def run(task: Task) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    receipt = json.loads(task.input("protein_receipt").read_text())
    donors = task.settings["donors"]
    if (receipt.get("stage") != "03_protein_measurement"
            or receipt.get("status") != "complete"
            or receipt.get("summary", {}).get("donors") != donors):
        raise IO60NormalizationError("accepted Stage 03 donor roster differs from the cohort")
    images = json.loads(task.input("cohort_spec").read_text())["images"]
    if [str(row["donor_id"]) for row in images] != donors:
        raise IO60NormalizationError("cohort changed after input declaration")
    markers, references = _references(task.input("reference_record"), donors)
    ratios = _ratios(task.input("usability_record"), markers, donors)
    blocks = normalization_block.read_blocks(task.input("pln_blocks"))
    block_counts = _population(task, images, blocks)
    sidecars = task.staging / "alignment_block"
    paths = {donor: task.input(f"expression:{donor}") for donor in donors}
    for donor, path in paths.items():
        missing = set([*_META, *markers]) - set(pq.ParquetFile(path).schema_arrow.names)
        if missing:
            raise IO60NormalizationError(f"{donor}: expression is missing {sorted(missing)}")
    masks = {donor: _mask(path, donor, sidecars) for donor, path in paths.items()}
    registration, normalization = normalize.load_uniform(Path(task.settings["uniform_checkout"]))
    records, scales, assessed, refused = [], {}, 0, []
    for marker in markers:
        sample = _values(paths, masks, marker)
        ratio = ratios[marker]
        reference = references[marker]
        result = None
        if len(sample) >= 2:
            with contextlib.redirect_stdout(io.StringIO()):
                result = normalize.fit_marker(
                    marker, sample, registration, normalization, reference=reference)
            plt.gcf().savefig(task.output(f"gate/markers/{marker}.png"),
                              dpi=110, bbox_inches="tight")
            plt.close("all")
        elif len(sample) == 1:
            # One sample has a measured marker, but no between-sample shift exists.
            reference = next(iter(sample))
        flagged = bool(ratio is not None and ratio >= ACQ_RATIO_FLAG)
        rejected = bool(result and flagged and result.mode_sd_after > MODE_SD_BAR)
        if rejected:
            refused.append(marker)
        scales[marker] = {}
        for donor in donors:
            scale = (result.scale.get(donor) if result and not rejected else
                     1.0 if donor in sample and result is None else None)
            status = ("refused" if rejected and donor in sample else
                      "aligned" if result and donor in sample else
                      "single_sample" if donor in sample else "unavailable")
            scales[marker][donor] = scale
            records.append({"donor_id": donor, "marker": marker,
                            "reference_donor": reference, "scale": scale,
                            "shift_log": result.shift_log.get(donor) if result else None,
                            "n_fit_cells": result.n_cells.get(donor, 0) if result else
                            int(np.count_nonzero(np.isfinite(sample.get(donor, [])))),
                            "status": status, "acq_ratio": ratio,
                            "mode_sd_before": result.mode_sd_before if result else None,
                            "mode_sd_after": result.mode_sd_after if result else None})
        assessed += 1
    with task.output("scales.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_SCALE_COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    in_block = _aligned(task, paths, markers, scales, sidecars)
    _figure(task, records, markers, len(donors))
    return {"donors": len(donors), "markers": len(markers), "markers_assessed": assessed,
            "different_block_cells": block_counts["different-block"],
            "unassigned_pln_cells": block_counts["unassigned"],
            "in_alignment_block": in_block, "refused_markers": refused,
            "CD8_shifted_donors": sum(r["marker"] == "CD8" and r["status"] == "aligned"
                                     for r in records)}


def _figure(task: Task, records: list[dict], markers: list[str], donors: int) -> None:
    import matplotlib.pyplot as plt

    rows = {marker: next(row for row in records if row["marker"] == marker)
            for marker in markers}
    before = [rows[marker]["mode_sd_before"]
              if rows[marker]["mode_sd_before"] is not None else np.nan for marker in markers]
    after = [rows[marker]["mode_sd_after"]
             if rows[marker]["mode_sd_after"] is not None else np.nan for marker in markers]
    fig, ax = plt.subplots(figsize=(11, max(5, len(markers) * 0.22 + 1)),
                           layout="constrained")
    y = np.arange(len(markers))
    ax.scatter(before, y, s=13, label="Before")
    ax.scatter(after, y, s=13, label="After")
    ax.axvline(MODE_SD_BAR, color="crimson", linestyle="--", label="review bar 0.10")
    ax.set(yticks=y, yticklabels=markers,
           xlabel="Background-mode spread across donors (SD, log intensity)",
           title=f"UniFORM: {len(markers)} markers across {donors} donor blocks")
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    fig.savefig(gates.figure_path(task, "mode_spread.png"), dpi=140)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    checks = [
        gates.Check("all markers assessed", summary["markers_assessed"],
                    str(summary["markers"]),
                    summary["markers_assessed"] == summary["markers"]),
        gates.Check("CD8 shifted on every donor block", summary["CD8_shifted_donors"],
                    str(summary["donors"]),
                    summary["CD8_shifted_donors"] == summary["donors"]),
        gates.Check("refused marker shifts", summary["refused_markers"],
                    "unavailable in aligned cells; review marker figures", True, blocking=False),
    ]
    figures = ["mode_spread.png"] + [f"markers/{marker}.png" for marker in
                                    sorted({row["marker"] for row in _read_scales(task)
                                            if row["mode_sd_after"]})]
    return gates.Gate(STAGE, "Same-block IO60 UniFORM shifts from Stage 03 expression; "
                      "historical references and usability are declared migration inputs.",
                      checks, figures)


def _read_scales(task: Task) -> list[dict]:
    with (task.staging / "scales.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))
