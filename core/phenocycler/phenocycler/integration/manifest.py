"""
S0 — the donor <-> Xenium-run pairing manifest.

This is the module that resolves the crosswalk problem, and it is worth being explicit
about what that problem is, because it is not obvious from either repo in isolation:

1. ``donor_id`` is a perfect join key *on paper*. Both repos already share
   ``donor_metadata_panc.xlsx`` verbatim, with the same ``donor_id`` and ``disease.status``
   columns and the same nPOD donors. PhenoCycler derives ``donor_id`` by regex from the
   qptiff image name and hive-partitions everything by it. (``xenium_paths.csv`` lists 26
   donors; the run cohort is a 20-donor subset, which nothing in this repo currently pins —
   the extra Xenium donors surface as ``xenium_only`` rather than being filtered out.)

2. But the join key never reaches the Xenium data. No Xenium h5ad carries ``obs['donor_id']``
   or ``obs['roi']`` — Xenium artifacts are keyed purely by the XETG slide serial
   (``0041323``), which encodes no donor, no region, and no disease status.

3. And ``xenium_paths.csv``, the one file that *does* map donors to Xenium runs, is read by
   nothing in the Xenium repo. It is hand-maintained and unenforced. Worse, the samples that
   actually have processed h5ads (``0041323``, ``0041326``) do not appear in it at all — the
   manifest cohort and the processed cohort are disjoint.

4. Finally, the slide serial is not a section identifier. One serial can carry both a
   pancreas and a pancreatic-lymph-node region (donor 6539 -> serial ``0059865`` with both a
   ``__Panc__`` and a ``__pLN__`` bundle), so keying anything by serial alone silently
   collides two different tissues.

So this module makes the manifest load-bearing rather than decorative: it normalises donor
ids to strings on both sides, keys Xenium sections by ``{serial}__{roi}`` instead of the
bare serial, classifies every row with a ``pair_status`` so unresolved cases are *visible*
rather than silently dropped, and provides a ``--stamp`` path that writes the join key into
the Xenium cell tables where it belongs.

Rows whose donor cannot be resolved (the processed-sample case) are not an error and do not
halt anything. They are recorded as ``donor_unknown`` and excluded from paired analyses
until someone supplies the mapping in ``donor_overrides.csv`` — a two-column file that is
the right place for knowledge that lives in a lab notebook rather than in the data.

CLI::

    python -m phenocycler.integration.manifest --report
    python -m phenocycler.integration.manifest --build
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..config import PipelineConfig, load_config
from ..sections import SectionError
from ..sections import load_overrides as load_section_overrides
from ..sections import resolve as resolve_section
from .contract import normalize_donor_id
from .tissues import for_roi as tissue_for_roi

#: Columns of ``data/integration/manifest.csv``, in order.
MANIFEST_COLUMNS = [
    "donor_id",
    "roi",
    "disease_status",
    "pair_status",
    "pheno_image",
    "pheno_geojson",
    "xenium_sample_key",
    "xenium_serial",
    "xenium_batch",
    "xenium_bundle",
    "xenium_zarr",
    "xenium_h5ad",
    "section_gap_um",
    "block_id",
    "tissue",
]

#: ``pair_status`` vocabulary.
PAIRED = "paired"                  # both modalities present for this donor
PHENO_ONLY = "pheno_only"          # PhenoCycler has it, Xenium does not
XENIUM_ONLY = "xenium_only"        # Xenium has it, PhenoCycler does not
DONOR_UNKNOWN = "donor_unknown"    # a Xenium run whose donor is not recorded anywhere

#: Default ROI for callers that intentionally request a single section. Pairability is
#: established row by row from observed manifests, never from the configured tissue list.
#: The current PhenoCycler run is pancreas-only even though the Xenium catalog includes pLN.
PAIRABLE_ROI = "panc"

#: ``output-XETG00298__0059865__Panc__20250328__202338`` -> serial ``0059865``, label ``Panc``.
_BUNDLE_RE = re.compile(r"output-[A-Za-z0-9]+__(?P<serial>\d+)__(?P<label>.*?)__\d{8}__\d+/?$")


class ManifestError(RuntimeError):
    """The manifest cannot be built."""


@dataclass
class ManifestSummary:
    """Counts and diagnostics from a manifest build — what ``--report`` prints."""

    n_rows: int = 0
    by_status: dict[str, int] = field(default_factory=dict)
    paired_donors: list[str] = field(default_factory=list)
    pheno_only_donors: list[str] = field(default_factory=list)
    xenium_only_donors: list[str] = field(default_factory=list)
    unknown_keys: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def render(self) -> str:
        lines = [
            "PhenoCycler <-> Xenium pairing manifest",
            "=" * 60,
            f"rows: {self.n_rows}",
        ]
        for status in (PAIRED, PHENO_ONLY, XENIUM_ONLY, DONOR_UNKNOWN):
            lines.append(f"  {status:<14} {self.by_status.get(status, 0)}")
        lines.append("")
        lines.append(f"paired donors ({len(self.paired_donors)}): "
                     f"{', '.join(self.paired_donors) or '-'}")
        if self.pheno_only_donors:
            lines.append(f"PhenoCycler only ({len(self.pheno_only_donors)}): "
                         f"{', '.join(self.pheno_only_donors)}")
        if self.xenium_only_donors:
            lines.append(f"Xenium only ({len(self.xenium_only_donors)}): "
                         f"{', '.join(self.xenium_only_donors)}")
        if self.unknown_keys:
            lines.append("")
            lines.append(f"UNRESOLVED Xenium runs ({len(self.unknown_keys)}) — no donor_id "
                         "recorded; add rows to donor_overrides.csv to pair them:")
            for key in self.unknown_keys:
                lines.append(f"    {key}")
        if self.warnings:
            lines.append("")
            lines.append("warnings:")
            for w in self.warnings:
                lines.append(f"    {w}")
        return "\n".join(lines)


def parse_bundle_path(path: str) -> tuple[str, str]:
    """``.../output-XETG00298__0059865__Panc__20250328__202338`` -> ``('0059865', 'Panc')``.

    Returns ``('', '')`` when the path does not look like a Xenium bundle. The serial is the
    only durable identifier in the whole Xenium layout — everything downstream of ingest
    keeps it and discards the rest — so it is worth recovering carefully.
    """
    text = str(path).strip().rstrip("/")
    if not text or text.upper() == "NA":
        return "", ""
    m = _BUNDLE_RE.search(text + "/")
    if not m:
        return "", ""
    return m.group("serial"), m.group("label")


def sample_key(serial: str, roi: str) -> str:
    """The Xenium section key: ``{serial}__{roi}``.

    Never key by serial alone. A single slide serial can carry a pancreas region and a
    lymph-node region, so ``data/raw/{serial}.zarr`` collides between two different tissues.
    """
    serial = str(serial).strip()
    roi = str(roi).strip()
    if not serial:
        return ""
    return f"{serial}__{roi}" if roi else serial


def normalize_roi(label: str) -> str:
    """Map a free-text Xenium bundle label to the ``roi`` vocabulary.

    The label field is entirely unnormalised upstream — ``Panc``, ``Pan``, ``Pancreas``,
    ``PanT``, ``6380_Pan``, ``pLN``, ``pLN1``, ``6457_pLN1``, ``Region_1`` all occur. Only
    the pancreas rows are pairable against PhenoCycler, so getting this mapping right is
    what decides whether a donor is integrable.
    """
    text = str(label).strip().lower()
    if not text:
        return ""
    if "pln" in text or "_ln" in text:
        m = re.search(r"pln[_]?(\d+)", text)
        return f"pln_{m.group(1)}" if m else "pln_1"
    if "pan" in text:
        return PAIRABLE_ROI
    return text


def read_xenium_paths(path: Path) -> pd.DataFrame:
    """Read ``xenium_paths.csv`` -> normalised rows.

    The file is hand-maintained and exported from a spreadsheet, so it can carry blank
    spacer rows and an ``NA`` source path for a run nobody has located yet. Both are
    dropped here rather than propagated — but a row dropped for a *missing path* is a real
    gap in the cohort, not noise, so it is recorded in ``df.attrs['dropped_no_path']`` and
    surfaced by ``build_manifest`` as a warning. Silently discarding it is how a donor
    section goes missing and reappears downstream as "no registration produced".
    """
    cols = ["donor_id", "roi", "batch", "source_path", "xenium_serial", "xenium_sample_key"]
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=cols)

    try:
        raw = pd.read_csv(path, dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        # A zero-byte or header-only file is a legitimate state (nothing vendored yet), not
        # an error — build_manifest reports it as a warning and carries on.
        return pd.DataFrame(columns=cols)

    raw.columns = [c.strip() for c in raw.columns]
    if "donor_id" not in raw.columns:
        raise ManifestError(
            f"{path} has no 'donor_id' column (found {list(raw.columns)}); expected "
            f"donor_id,roi,batch,source_path")
    raw = raw[raw["donor_id"].astype(str).str.strip() != ""].copy()
    if raw.empty:
        return pd.DataFrame(columns=cols)

    raw["donor_id"] = raw["donor_id"].map(normalize_donor_id)
    raw["roi"] = raw["roi"].astype(str).str.strip()
    raw["batch"] = raw.get("batch", "").astype(str).str.strip()
    raw["source_path"] = raw["source_path"].astype(str).str.strip()
    no_path = raw["source_path"].str.upper().isin({"", "NA"})
    dropped = [f"{d} {r or '?'}" for d, r in
               zip(raw.loc[no_path, "donor_id"], raw.loc[no_path, "roi"])]
    raw = raw[~no_path].copy()

    parsed = raw["source_path"].map(parse_bundle_path)
    raw["xenium_serial"] = [p[0] for p in parsed]
    # Trust the manifest's own roi column; fall back to the bundle label when it is blank.
    labels = [p[1] for p in parsed]
    raw["roi"] = [r if r else normalize_roi(lab) for r, lab in zip(raw["roi"], labels)]
    raw["xenium_sample_key"] = [sample_key(s, r)
                                for s, r in zip(raw["xenium_serial"], raw["roi"])]
    out = raw.reset_index(drop=True)
    out.attrs["dropped_no_path"] = dropped
    return out


def read_donor_overrides(path: Path) -> pd.DataFrame:
    """Read ``donor_overrides.csv``: ``xenium_sample_key, donor_id, roi[, section_gap_um, block_id]``.

    This is where knowledge that exists only outside the data goes — most immediately, which
    nPOD donors ``0041323`` and ``0041326`` are. Overrides win over ``xenium_paths.csv``.
    """
    path = Path(path)
    cols = ["xenium_sample_key", "donor_id", "roi", "section_gap_um", "block_id"]
    if not path.exists():
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=cols)
    df.columns = [c.strip() for c in df.columns]
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df["donor_id"] = df["donor_id"].map(normalize_donor_id)
    df["xenium_sample_key"] = df["xenium_sample_key"].astype(str).str.strip()
    df["roi"] = df["roi"].astype(str).str.strip()
    return df.loc[:, cols].reset_index(drop=True)


def read_donor_metadata(cfg: PipelineConfig) -> pd.DataFrame:
    """``donor_id -> disease_status`` from the shared workbook.

    Both repos point at the same file with the same column names, so this is the one piece
    of donor-level metadata that needs no reconciliation — only a dtype fix, because the
    workbook stores donor ids numerically.
    """
    path = Path(cfg.donor_metadata)
    if not path.exists():
        return pd.DataFrame(columns=["donor_id", "disease_status"])
    try:
        df = pd.read_excel(path)
    except Exception as exc:  # noqa: BLE001 - missing openpyxl, corrupt file, ...
        print(f"[manifest] WARNING: could not read {path}: {exc}", flush=True)
        return pd.DataFrame(columns=["donor_id", "disease_status"])

    donor_col = cfg.metadata_donor_col
    status_col = cfg.metadata_status_col
    if donor_col not in df.columns:
        raise ManifestError(f"{path} has no {donor_col!r} column (got {list(df.columns)})")
    out = pd.DataFrame({
        "donor_id": df[donor_col].map(normalize_donor_id),
        "disease_status": (df[status_col].astype(str).str.strip()
                           if status_col in df.columns else ""),
    })
    return out.drop_duplicates(subset="donor_id").reset_index(drop=True)


def _pheno_paths(cfg: PipelineConfig, donor_id: str, roi: str = PAIRABLE_ROI) -> tuple[str, str]:
    """Best-effort ``(qptiff, geojson)`` from the authoritative QuPath manifest.

    Filename globbing left 8/20 cohort images unresolved and could substitute another
    section. The manifest already owns donor, source paths and image identity; exact tracked
    overrides handle its two irregular names. Missing mounted paths remain empty because the
    tabular stages do not need imaging data.
    """

    try:
        payload = json.loads(Path(cfg.qupath_manifest).read_text())
        overrides = load_section_overrides(cfg.section_overrides_csv)
    except (OSError, json.JSONDecodeError, SectionError):
        return "", ""

    matches = []
    for row in payload.get("images", []):
        image_id = str(row.get("image_id", "")).strip()
        try:
            section = resolve_section(
                image_id,
                donor_id=normalize_donor_id(row.get("donor_id")),
                overrides=overrides,
            )
        except SectionError:
            continue
        if (section.donor_id, section.roi) == (str(donor_id), str(roi)):
            matches.append(row)
    if len(matches) != 1:
        return "", ""

    def _path(field: str) -> str:
        value = matches[0].get(field, "")
        raw = value.get("path", "") if isinstance(value, dict) else value
        path = Path(str(raw)) if raw else None
        return str(path) if path is not None and path.exists() else ""

    return _path("qptiff"), _path("cell_geojson")


def _discover_pheno_sections(cfg: PipelineConfig,
                             summary: "ManifestSummary") -> Optional[set[tuple[str, str]]]:
    """``{(donor_id, roi)}`` from the current immutable run and QuPath manifest.

    ``None`` means "nothing exported yet" and is distinct from an empty set: the manifest
    still builds a Xenium-shaped picture rather than claiming every donor is unpaired.

    The old implementation globbed ``data/cells/donor_id=*``. Current core outputs live at
    ``data/runs/<run_id>/expression`` and a donor directory alone cannot state which physical
    section it contains. The image id in that parquet is therefore resolved against
    ``qupath_manifest.json`` plus exact, tracked section overrides.
    """

    latest = cfg.runs_dir / "LATEST"
    if not latest.exists():
        return None
    run_id = latest.read_text().strip()
    expression = cfg.runs_dir / run_id / "expression"
    if not expression.is_dir():
        summary.warnings.append(
            f"{latest} names run {run_id!r}, but it has no expression directory"
        )
        return None

    try:
        manifest_payload = json.loads(Path(cfg.qupath_manifest).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        summary.warnings.append(f"cannot read QuPath manifest {cfg.qupath_manifest}: {exc}")
        return None
    by_image = {
        str(row.get("image_id", "")).strip(): row
        for row in manifest_payload.get("images", [])
        if str(row.get("image_id", "")).strip()
    }
    try:
        overrides = load_section_overrides(cfg.section_overrides_csv)
    except SectionError as exc:
        summary.warnings.append(str(exc))
        return None

    import pyarrow.parquet as pq

    out: set[tuple[str, str]] = set()
    for part in sorted(expression.glob("donor_id=*")):
        files = sorted(part.rglob("*.parquet"))
        if not files:
            continue
        try:
            images: set[str] = set()
            for file in files:
                parquet = pq.ParquetFile(file)
                image_index = parquet.schema_arrow.names.index("image")
                needs_scan = False
                for group_index in range(parquet.metadata.num_row_groups):
                    stats = parquet.metadata.row_group(group_index).column(image_index).statistics
                    if stats is None or not stats.has_min_max:
                        needs_scan = True
                        break
                    images.add(str(stats.min).strip())
                    images.add(str(stats.max).strip())
                if needs_scan:
                    # A full image column can contain millions of repeated strings. Only
                    # take this slower path for parquet writers that omitted statistics.
                    column = parquet.read(columns=["image"])["image"]
                    images.update(str(value.as_py()).strip() for value in column.unique())
            images.discard("")
            if len(images) != 1:
                raise SectionError(
                    f"run {run_id} partition {part.name!r} carries {len(images)} image ids "
                    f"({sorted(images)}); one partition must be one physical section"
                )
            image_id = next(iter(images))
            record = by_image.get(image_id)
            if record is None:
                raise SectionError(
                    f"run {run_id} image {image_id!r} is absent from {cfg.qupath_manifest}"
                )
            sec = resolve_section(
                image_id,
                donor_id=normalize_donor_id(record.get("donor_id")),
                overrides=overrides,
            )
        except SectionError as exc:
            summary.warnings.append(f"PhenoCycler partition {part.name!r} skipped: {exc}")
            continue
        out.add((sec.donor_id, sec.roi))
    return out


def _rewrite_root(path: str, xenium_root: str) -> str:
    """Resolve a recorded bundle path against ``xenium_root``.

    ``xenium_paths.csv`` records paths **relative** to the Xenium tree::

        HuPrime5kPanTissue/<batch>/<run>/output-XETG00298__0059865__Panc__...

    which is all that identity needs — :func:`parse_bundle_path` reads the serial off the
    bundle directory name and the prefix contributes nothing. Keeping the vendored manifest
    site-neutral means no lab's filesystem layout is baked into a tracked file, and the tree
    can live wherever it was transferred to.

    Absolute paths carrying a ``/xenium/`` segment are still re-rooted, so a manifest that
    predates the migration to relative paths keeps working. An absolute path without that
    segment is returned unchanged and reported as ``unrewritten`` by
    :func:`verify_bundle_paths` — silently resolving it against the recorded root is how a
    wrong ``xenium_root`` turns into a "missing file" for the wrong reason.
    """
    if not path:
        return path
    if not path.startswith("/"):
        # Relative — the normal case. Without a root there is nothing to resolve against;
        # returning it unchanged would silently resolve against the process's cwd.
        return str(Path(xenium_root) / path) if xenium_root else path
    if not xenium_root:
        return path
    marker = "/xenium/"
    idx = path.find(marker)
    if idx == -1:
        return path
    return str(Path(xenium_root) / path[idx + len(marker):])


def _find_bundle(resolved: str, xenium_root: str) -> str:
    """Locate a bundle by its directory name when the recorded path does not resolve.

    The recorded relative path assumes the ``HuPrime5kPanTissue/<batch>/<run>/`` levels
    survived whatever copy or transfer put the tree under ``xenium_root``. They may not
    have. The bundle directory name is the part that cannot change without changing the
    sample's identity, so searching for it recovers a tree that was reorganised, while
    still refusing to invent one that is absent.

    Returns ``""`` when there is no unambiguous match — an ambiguous hit is a refusal, not
    a coin flip.
    """
    name = Path(resolved).name
    if not name or not xenium_root:
        return ""
    root = Path(xenium_root)
    if not root.is_dir():
        return ""
    hits = [p for p in root.rglob(name) if p.is_dir()]
    return str(hits[0]) if len(hits) == 1 else ""


#: Files every Xenium output bundle carries. Used to tell "the directory exists" from
#: "the directory is the bundle we think it is".
BUNDLE_MARKERS = ("cells.parquet", "cell_boundaries.parquet", "cell_feature_matrix.h5")


def verify_bundle_paths(cfg: PipelineConfig,
                        donors: Optional[set[str]] = None) -> list[dict]:
    """Resolve every recorded bundle path onto this filesystem and classify what happens.

    ``xenium_paths.csv`` is hand-maintained and, upstream, read by nothing — so nothing has
    ever checked that a path in it points at a bundle that exists. This is that check.

    Two statuses carry more information than "it is not there":

    ``no_root``
        the manifest records relative paths (the normal case) and ``xenium_root`` is unset,
        so there is nothing to resolve against. That is a configuration gap, not a missing
        bundle, and calling it ``missing`` would send the reader to the wrong problem.

    ``unrewritten``
        an *absolute* recorded path carries no ``/xenium/`` segment, so ``_rewrite_root``
        left it pointing at the filesystem it was recorded on. Setting ``xenium_root``
        against such a manifest silently resolves against the original root and reports
        "missing" for the wrong reason.

    ``relocated`` is not a failure: the recorded relative path did not resolve, but the
    bundle directory was found elsewhere under ``xenium_root``. Reported rather than
    silently accepted, because it means the recorded path is stale.
    """
    xen = read_xenium_paths(cfg.xenium_paths_csv)
    out: list[dict] = []
    for _, r in xen.iterrows():
        donor = r["donor_id"]
        if donors and donor not in donors:
            continue
        recorded = r["source_path"]
        resolved = _rewrite_root(recorded, cfg.xenium_root)
        relative = not recorded.startswith("/")
        relocated = False
        if not Path(resolved).is_dir():
            found = _find_bundle(resolved, cfg.xenium_root)
            if found:
                resolved, relocated = found, True

        if relative and not cfg.xenium_root:
            status = "no_root"
        elif cfg.xenium_root and not relative and resolved == recorded:
            status = "unrewritten"
        elif not Path(resolved).is_dir():
            status = "missing"
        elif not any((Path(resolved) / m).exists() for m in BUNDLE_MARKERS):
            status = "not_a_bundle"
        else:
            status = "relocated" if relocated else "ok"
        out.append({"donor_id": donor, "roi": r["roi"], "batch": r["batch"],
                    "recorded": recorded, "resolved": resolved, "status": status})
    for gap in xen.attrs.get("dropped_no_path", []):
        d, _, roi = gap.partition(" ")
        if not donors or d in donors:
            out.append({"donor_id": d, "roi": roi, "batch": "", "recorded": "",
                        "resolved": "", "status": "no_path"})
    return out


def render_path_verification(rows: list[dict], cfg: PipelineConfig) -> str:
    """Human-readable verification report. Every non-``ok`` row is named, not counted."""
    order = ("ok", "relocated", "missing", "no_root", "unrewritten", "not_a_bundle",
             "no_path")
    counts = {s: sum(1 for r in rows if r["status"] == s) for s in order}
    lines = [
        "Xenium bundle path verification",
        "=" * 60,
        f"manifest    : {cfg.xenium_paths_csv}",
        f"xenium_root : {cfg.xenium_root or '(unset — paths used as recorded)'}",
        f"rows        : {len(rows)}",
        "",
    ]
    lines += [f"  {s:<13} {counts[s]}" for s in order if counts[s]]
    for status in order[1:]:
        bad = [r for r in rows if r["status"] == status]
        if not bad:
            continue
        lines += ["", f"{status.upper()} ({len(bad)}):"]
        lines += [f"  {r['donor_id']:<6} {r['roi']:<7} {r['resolved'] or r['recorded'] or '-'}"
                  for r in bad]
    if counts["no_root"]:
        lines += ["", "'no_root' means the manifest records relative paths — which is the "
                      "normal case — but [integration] xenium_root is unset, so there is "
                      "nothing to resolve them against. Point it at the local Xenium tree."]
    if counts["relocated"]:
        lines += ["", "'relocated' means the bundle was found under xenium_root but NOT at "
                      "its recorded relative path. The data is fine; the manifest is stale. "
                      "Update source_path to the path shown above."]
    if counts["unrewritten"]:
        lines += ["", "'unrewritten' means xenium_root is set but these ABSOLUTE paths carry "
                      "no '/xenium/' segment, so they were left pointing at the root they "
                      "were recorded on. Make them relative to the Xenium tree — do not "
                      "ignore this."]
    return "\n".join(lines)


def build_manifest(
    cfg: PipelineConfig,
    *,
    pheno_donors: Optional[list[str]] = None,
    pheno_sections: Optional[list[tuple[str, str]]] = None,
    xenium_processed_root: Optional[Path] = None,
) -> tuple[pd.DataFrame, ManifestSummary]:
    """Build the pairing manifest.

    Parameters
    ----------
    cfg
        Resolved pipeline config (supplies both CSV paths, the donor workbook and
        ``xenium_root``).
    pheno_donors
        Override the PhenoCycler donor list. Pairing then falls back to "this donor has
        PhenoCycler data, so every ROI of a known tissue is pairable" — the best that can be
        said without knowing which sections exist.
    pheno_sections
        Override the PhenoCycler ``(donor_id, roi)`` list. This is the precise form and the
        default: sections are discovered from ``data/cells/donor_id=*``, whose keys carry the
        region (``6539`` pancreas, ``6539pln`` lymph node). With it, a donor whose lymph node
        was never scanned is ``xenium_only`` for that ROI instead of being reported as a
        pairing that does not exist.
    xenium_processed_root
        Where ``{sample}/{sample}_phenotyped.h5ad`` and ``{sample}.zarr`` live, if mounted.
        Used only to fill the ``xenium_h5ad`` / ``xenium_zarr`` columns.

    Returns
    -------
    (DataFrame, ManifestSummary)
    """
    summary = ManifestSummary()

    # `pheno_pairs` is the precise statement (which sections exist); `pheno_set` is the
    # fallback (which donors exist). `pheno_pairs is None` means "sections unknown — pair on
    # the donor and the ROI's tissue", which is what an explicit `pheno_donors` list implies.
    pheno_pairs: Optional[set[tuple[str, str]]] = None
    if pheno_sections is not None:
        pheno_pairs = {(normalize_donor_id(d), str(r)) for d, r in pheno_sections}
        pheno_set = {d for d, _ in pheno_pairs if d}
    elif pheno_donors is not None:
        pheno_set = {normalize_donor_id(d) for d in pheno_donors if normalize_donor_id(d)}
    else:
        discovered = _discover_pheno_sections(cfg, summary)
        if discovered is None:
            pheno_set = set()
        else:
            pheno_pairs = discovered
            pheno_set = {d for d, _ in pheno_pairs if d}

    xen = read_xenium_paths(cfg.xenium_paths_csv)
    overrides = read_donor_overrides(cfg.donor_overrides_csv)
    meta = read_donor_metadata(cfg)
    status_by_donor = dict(zip(meta["donor_id"], meta["disease_status"]))

    for gap in xen.attrs.get("dropped_no_path", []):
        summary.warnings.append(
            f"{cfg.xenium_paths_csv.name}: {gap} has no source_path — that donor section is "
            f"absent from the manifest and cannot be paired. Fill the path or record why."
        )

    if xen.empty and overrides.empty:
        summary.warnings.append(
            f"no Xenium runs found — {cfg.xenium_paths_csv} is missing or empty, and "
            f"{cfg.donor_overrides_csv} has no rows"
        )

    rows: list[dict] = []
    seen_keys: set[str] = set()

    override_by_key = {r["xenium_sample_key"]: r for _, r in overrides.iterrows()
                       if r["xenium_sample_key"]}

    for _, r in xen.iterrows():
        key = r["xenium_sample_key"]
        ov = override_by_key.get(key, {})
        donor = normalize_donor_id(ov.get("donor_id") or r["donor_id"])
        roi = str(ov.get("roi") or r["roi"])
        bundle = _rewrite_root(r["source_path"], cfg.xenium_root)

        rows.append({
            "donor_id": donor,
            "roi": roi,
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": "",                       # filled below
            "pheno_image": "",
            "pheno_geojson": "",
            "xenium_sample_key": key,
            "xenium_serial": r["xenium_serial"],
            "xenium_batch": r["batch"],
            "xenium_bundle": bundle,
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": str(ov.get("section_gap_um", "") or ""),
            "block_id": str(ov.get("block_id", "") or ""),
            "tissue": tissue_for_roi(roi),
        })
        seen_keys.add(key)

    # Override rows that name a Xenium sample absent from xenium_paths.csv. This is exactly
    # the 0041323 / 0041326 case: they have processed h5ads and no manifest row at all.
    for _, r in overrides.iterrows():
        key = r["xenium_sample_key"]
        if not key or key in seen_keys:
            continue
        serial = key.split("__", 1)[0]
        donor = normalize_donor_id(r["donor_id"])
        rows.append({
            "donor_id": donor,
            "roi": r["roi"] or (key.split("__", 1)[1] if "__" in key else ""),
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": "",
            "pheno_image": "",
            "pheno_geojson": "",
            "xenium_sample_key": key,
            "xenium_serial": serial,
            "xenium_batch": "",
            "xenium_bundle": "",
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": str(r.get("section_gap_um", "") or ""),
            "block_id": str(r.get("block_id", "") or ""),
            "tissue": tissue_for_roi(r["roi"] or (key.split("__", 1)[1] if "__" in key else "")),
        })
        seen_keys.add(key)

    df = pd.DataFrame(rows, columns=MANIFEST_COLUMNS) if rows else pd.DataFrame(
        columns=MANIFEST_COLUMNS)

    # Fill processed-artifact paths where the Xenium outputs are mounted.
    if xenium_processed_root is not None and len(df):
        root = Path(xenium_processed_root)
        h5ads, zarrs = [], []
        # Prefer a per-(serial, roi) file and only fall back to the serial-level one when
        # the serial carries a single ROI. One slide serial can hold a pancreas region and
        # several lymph-node regions — donor 6579's serial carries panc + three pLNs — so a
        # serial-keyed lookup gave all four rows the same h5ad. Because `_resolve_source`
        # prefers h5ad over the (correct, per-ROI) bundle, that silently imported one
        # section's labels for four different pieces of tissue.
        rois_per_serial = df.groupby("xenium_serial")["roi"].nunique().to_dict()
        for serial, roi in zip(df["xenium_serial"], df["roi"]):
            candidates = [
                root / serial / f"{serial}__{roi}_phenotyped.h5ad",
                root / serial / roi / f"{serial}_phenotyped.h5ad",
            ]
            if rois_per_serial.get(serial, 0) <= 1:
                candidates.append(root / serial / f"{serial}_phenotyped.h5ad")
            h = next((c for c in candidates if c.exists()), None)

            z_candidates = [root.parent / "raw" / f"{serial}__{roi}.zarr"]
            if rois_per_serial.get(serial, 0) <= 1:
                z_candidates.append(root.parent / "raw" / f"{serial}.zarr")
            z = next((c for c in z_candidates if c.exists()), None)

            h5ads.append(str(h) if h is not None else "")
            zarrs.append(str(z) if z is not None else "")
        df["xenium_h5ad"] = h5ads
        df["xenium_zarr"] = zarrs
        ambiguous = sorted({s for s, n in rois_per_serial.items() if n > 1})
        if ambiguous:
            print(f"[manifest] serial(s) {ambiguous} carry several ROIs; used per-ROI "
                  f"processed files only. Rows left blank fall back to their own "
                  f"xenium_bundle, which is already per-ROI.", flush=True)

    # PhenoCycler-side paths + pairing status.
    xen_donors = {d for d in df["donor_id"] if d}
    if len(df):
        statuses, images, geojsons = [], [], []
        for _, r in df.iterrows():
            donor, roi = r["donor_id"], r["roi"]
            if not donor:
                statuses.append(DONOR_UNKNOWN)
                images.append("")
                geojsons.append("")
                continue
            paired = ((donor, roi) in pheno_pairs if pheno_pairs is not None
                      else (donor in pheno_set and bool(tissue_for_roi(roi))))
            if paired:
                q, g = _pheno_paths(cfg, donor, roi)
                statuses.append(PAIRED)
                images.append(q)
                geojsons.append(g)
            else:
                statuses.append(XENIUM_ONLY)
                images.append("")
                geojsons.append("")
        df["pair_status"] = statuses
        df["pheno_image"] = images
        df["pheno_geojson"] = geojsons

    # PhenoCycler sections with no Xenium counterpart get their own rows, so the manifest is a
    # complete picture of the cohort rather than a Xenium-shaped view of it. With sections
    # known, this is exact: a donor whose pancreas paired but whose lymph node has no Xenium
    # run gets one `pheno_only` lymph-node row, rather than the whole donor being counted
    # either way on the strength of its pancreas.
    have = {(str(r["donor_id"]), str(r["roi"])) for _, r in df.iterrows()} if len(df) else set()
    if pheno_pairs is not None:
        orphans = sorted(p for p in pheno_pairs if p not in have and p[0])
    else:
        orphans = [(d, PAIRABLE_ROI) for d in sorted(pheno_set - xen_donors)]

    for donor, roi in orphans:
        q, g = _pheno_paths(cfg, donor, roi)
        df = pd.concat([df, pd.DataFrame([{
            "donor_id": donor,
            "roi": roi,
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": PHENO_ONLY,
            "pheno_image": q,
            "pheno_geojson": g,
            "xenium_sample_key": "",
            "xenium_serial": "",
            "xenium_batch": "",
            "xenium_bundle": "",
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": "",
            "block_id": "",
            "tissue": tissue_for_roi(roi),
        }])], ignore_index=True)

    if len(df):
        df = df.sort_values(["donor_id", "roi", "xenium_sample_key"],
                            kind="stable").reset_index(drop=True)

    # Collision guard. If this ever fires, two different sections are sharing an identity and
    # every downstream artifact for them would overwrite the other.
    if len(df):
        keys = df.loc[df["xenium_sample_key"] != "", "xenium_sample_key"]
        dupes = sorted(keys[keys.duplicated()].unique())
        if dupes:
            raise ManifestError(
                f"duplicate xenium_sample_key values {dupes} — two Xenium runs share a "
                "(serial, roi) identity, which would silently overwrite one another"
            )

    summary.n_rows = len(df)
    summary.by_status = df["pair_status"].value_counts().to_dict() if len(df) else {}
    summary.paired_donors = sorted(set(df.loc[df["pair_status"] == PAIRED, "donor_id"]))
    summary.pheno_only_donors = sorted(set(df.loc[df["pair_status"] == PHENO_ONLY, "donor_id"]))
    summary.xenium_only_donors = sorted(
        set(df.loc[df["pair_status"] == XENIUM_ONLY, "donor_id"]) - set(summary.paired_donors)
    )
    summary.unknown_keys = sorted(
        df.loc[df["pair_status"] == DONOR_UNKNOWN, "xenium_sample_key"])

    if summary.unknown_keys:
        summary.warnings.append(
            f"{len(summary.unknown_keys)} Xenium run(s) have no donor_id. They are excluded "
            f"from paired analyses. Add them to {cfg.donor_overrides_csv} to include them."
        )
    if not summary.paired_donors:
        summary.warnings.append(
            "no paired donors — donor-level, niche and PhenoCycler-only outputs still build, "
            "but registration, islet matching and pseudo-cell linking have nothing to run on"
        )

    return df, summary


def write_manifest(cfg: PipelineConfig, df: pd.DataFrame) -> Path:
    out = cfg.manifest_csv
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    return out


def load_manifest(cfg: PipelineConfig) -> pd.DataFrame:
    """Read the built manifest, rebuilding it on demand if absent."""
    if cfg.manifest_csv.exists():
        df = pd.read_csv(cfg.manifest_csv, dtype=str).fillna("")
        df["donor_id"] = df["donor_id"].map(normalize_donor_id)
        return df
    df, _ = build_manifest(cfg)
    return df


def rows_for_scope(df: pd.DataFrame, roi: Optional[str] = None,
                   tissue: Optional[str] = None) -> pd.DataFrame:
    """Apply the ROI/tissue scope rule. The single definition of what a run covers.

    Three shapes, matching the CLI: one ROI, one tissue's ROIs, or — with neither — every
    ROI in the table. Every stage narrows through here rather than writing its own
    ``df["roi"] == roi``; that hand-rolled form is what made the lymph node invisible to
    half the pipeline, and it is the kind of thing a new stage copies without noticing.
    """
    sub = df
    if roi:
        sub = sub[sub["roi"] == roi]
    elif tissue:
        sub = sub[sub["tissue"] == tissue]
    return sub.reset_index(drop=True)


def paired_rows(df: pd.DataFrame, roi: Optional[str] = None,
                tissue: Optional[str] = None) -> pd.DataFrame:
    """The subset that can actually be integrated, narrowed to the run's scope.

    Defaults to **every** paired section, not to the pancreas. A tissue-shaped default in the
    one function whose job is scoping is precisely how the lymph node went missing from half
    the pipeline: every caller that forgot to pass an argument quietly analysed one tissue and
    reported it as the cohort. Pass ``tissue=`` or ``roi=`` to narrow.
    """
    return rows_for_scope(df[df["pair_status"] == PAIRED], roi=roi, tissue=tissue)


def run_manifest(cfg: PipelineConfig, *, write: bool = True,
                 xenium_processed_root: Optional[Path] = None) -> pd.DataFrame:
    """Stage entry point used by ``pipeline.py``."""
    df, summary = build_manifest(cfg, xenium_processed_root=xenium_processed_root)
    print(summary.render(), flush=True)
    if write:
        path = write_manifest(cfg, df)
        print(f"\n[manifest] wrote {path} ({len(df)} rows)", flush=True)
    return df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Build/report the PhenoCycler<->Xenium manifest")
    ap.add_argument("--config", default=None, help="path to config.ini")
    ap.add_argument("--report", action="store_true",
                    help="print the pairing table without writing it")
    ap.add_argument("--build", action="store_true", help="write data/integration/manifest.csv")
    ap.add_argument("--xenium-processed", default=None,
                    help="Xenium data/processed root, to resolve h5ad/zarr paths")
    ap.add_argument("--verify-paths", action="store_true",
                    help="resolve every recorded bundle path on this filesystem and exit "
                         "non-zero if any is missing, unrewritten or not a bundle")
    ap.add_argument("--donor", action="append", dest="donors",
                    help="restrict --verify-paths to these donors (repeatable)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)

    if args.verify_paths:
        rows = verify_bundle_paths(cfg, set(args.donors) if args.donors else None)
        print(render_path_verification(rows, cfg))
        return 0 if all(r["status"] == "ok" for r in rows) else 1

    root = Path(args.xenium_processed) if args.xenium_processed else None
    df, summary = build_manifest(cfg, xenium_processed_root=root)
    print(summary.render())
    if args.build or not args.report:
        path = write_manifest(cfg, df)
        print(f"\nwrote {path} ({len(df)} rows)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
