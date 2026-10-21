"""What each acquisition's instrument recorded: its panel, its run and the run's date.

* **PhenoCycler** -- the first page of a qptiff carries the instrument's description: the slide id,
  the experiment the slide was scanned under (``ModIO60_Pancreas_Slide1`` is experiment
  ``ModIO60_Pancreas``, slide position 1), the experiment's date, the acquisition software and,
  in its ScanProfile, every marker. The TIFF ``DateTime`` tag is when the scan was written. Only
  that first page is read. Slides scanned under one experiment on one day were run together --
  the maintainer, 2026-10-21: "The phenocycler was always run two slides per run, the only way a
  slide was run alone is if there happened to be an odd number in the staining batch."
* **Xenium** -- every region output carries ``experiment.xenium`` (run name and start time, slide,
  region, panel), ``gene_panel.json`` (the panel's design) and ``metrics_summary.csv`` (the run and
  panel again). Slides with one run name were run together. The three are read together and any
  difference between them is returned, never settled.

MultiPlex5000 D3 (2026-10-21).
"""

from __future__ import annotations

import csv
import html
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


class AcquisitionRecordError(ValueError):
    """An instrument record that cannot be read as the instrument writes it."""


@dataclass(frozen=True)
class QptiffHeader:
    path: str
    slide_id: str
    experiment: str          # the experiment name without its slide position
    position: str            # "Slide1", or "" when the name carries none
    experiment_date: str     # the description's Date, as written
    scan_time: str           # the TIFF DateTime tag, "YYYY-MM-DD HH:MM:SS"
    software: str
    markers: tuple[str, ...]

    @property
    def scan_date(self) -> str:
        return self.scan_time[:10]

    @property
    def run_key(self) -> tuple[str, str, str]:
        """Scanned under one experiment on one day: run together."""
        return (self.experiment, self.experiment_date, self.scan_date)


def _field(description: str, tag: str) -> str:
    found = re.search(rf"<{tag}>([^<]*)</{tag}>", description)
    return found.group(1).strip() if found else ""


def read_qptiff_header(path: Path) -> QptiffHeader:
    import tifffile

    try:
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]
            description = page.description or ""
            stamp = page.tags.get("DateTime")
            stamp = str(stamp.value) if stamp is not None else ""
    except (tifffile.TiffFileError, ValueError, OSError) as err:
        raise AcquisitionRecordError(f"{path}: not a readable TIFF ({err})") from err
    found = re.search(r"<ScanProfile>(.*?)</ScanProfile>", description, re.S)
    if not found:
        raise AcquisitionRecordError(f"{path}: the first page carries no ScanProfile")
    try:
        experiment = json.loads(html.unescape(found.group(1)))["experimentDescription"]
    except (ValueError, KeyError) as err:
        raise AcquisitionRecordError(f"{path}: the ScanProfile cannot be read ({err})") from err
    named = re.fullmatch(r"(.+)_(Slide\d+)", experiment.get("name", ""))
    markers: list[str] = []
    for well in experiment.get("wells", []):
        for item in well.get("items", []):
            marker = item.get("markerName")
            if marker and marker not in ("--", "Blank", "None") and marker not in markers:
                markers.append(marker)
    scan_time = (stamp[:10].replace(":", "-") + stamp[10:]) if stamp else ""
    header = QptiffHeader(
        path=str(path), slide_id=_field(description, "SlideID"),
        experiment=named.group(1) if named else experiment.get("name", ""),
        position=named.group(2) if named else "", experiment_date=_field(description, "Date"),
        scan_time=scan_time, software=_field(description, "AcquisitionSoftware"),
        markers=tuple(markers))
    missing = [name for name in ("experiment", "experiment_date", "scan_time")
               if not getattr(header, name)]
    if missing or not markers:
        raise AcquisitionRecordError(f"{path}: the first page does not record {missing or 'markers'}")
    return header


def phenocycler_runs(headers: list[QptiffHeader]) -> dict[tuple[str, str, str], list[QptiffHeader]]:
    """Slides by run: one experiment, one scan day."""
    runs: dict[tuple[str, str, str], list[QptiffHeader]] = defaultdict(list)
    for header in headers:
        runs[header.run_key].append(header)
    return {key: sorted(group, key=lambda h: h.scan_time) for key, group in sorted(runs.items())}


@dataclass(frozen=True)
class XeniumRecord:
    run_name: str
    run_start: str           # run_start_time as recorded, UTC
    slide_id: str
    region: str
    experiment: dict         # experiment.xenium as read; panel_names reads it
    differences: tuple[str, ...]   # where gene_panel.json or metrics_summary.csv says otherwise

    @property
    def run_date(self) -> str:
        return self.run_start[:10]


def read_xenium_records(experiment: Path, gene_panel: Path | None,
                        metrics: Path | None) -> XeniumRecord:
    try:
        record = json.loads(Path(experiment).read_text())
        run_name, run_start = record["run_name"], record["run_start_time"]
        slide_id, region = str(record["slide_id"]), record["region_name"]
    except (ValueError, KeyError) as err:
        raise AcquisitionRecordError(f"{experiment}: not an experiment record ({err})") from err
    differences = []
    design = record.get("panel_design_id")
    if gene_panel is not None:
        panel = json.loads(Path(gene_panel).read_text())["payload"]["panel"]
        if panel["identity"]["design_id"] != design:
            differences.append(f"gene_panel.json design {panel['identity']['design_id']} "
                               f"!= experiment.xenium {design}")
    if metrics is not None:
        with Path(metrics).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        for key, mine in (("run_name", run_name), ("panel_design_id", design),
                          ("region_name", region)):
            if key in row and row[key] != mine:
                differences.append(f"metrics_summary.csv {key} {row[key]} != {mine}")
    return XeniumRecord(run_name, run_start, slide_id, region, record, tuple(differences))
