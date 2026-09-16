"""PhenoCycler integration boundary: geometry always, identity only when declared.

Two exports behind one stage, split along what the caller actually needs.

:func:`export_geometry` is always available. It reads ``runs/<id>/expression`` — the
terminal stage of the current workflow — and emits a contract cell table carrying
coordinates, area, ``cell_region``/``qupath_islet_num``, the ``qc_analysis_eligible`` flag
as an audit column, and identity-safe intensities, with **no** ``lineage_native``. Hormone
intensities are withheld because their integration consumer applies an unverified
positivity floor. None of this claims what any cell is.

The **stage** (:func:`run_export_geometry`) exports the COMPLETE cell set of every
partition, QC-ineligible cells included, and carries the flag beside them. Two reasons.
A declared identity release covers every cell the arm saw and marks the rest
``Unavailable(reason)``; attaching it by exact ``object_id`` needs the frame to be the whole
frame, or the attach refuses on a UUID mismatch. And a downstream consumer that wants only
eligible cells filters on the flag it can see, rather than inheriting a filter it cannot.
The library function keeps ``analysis_eligible_only=True`` as its default so a direct caller
gets the analysis set it asked for.

After export the stage attaches the configured PhenoCycler release for every donor that
release declares (:func:`~.phenospace.load_release_or_none` then :func:`attach_identity`,
the same idiom ``import_xenium`` uses for Xenium). Donors the release does not declare, and
every donor while no release exists, are written geometry only and announced as such.

Identity export by any other route stays fail-closed. There is no verified marker-pair or
phenotype workflow, and historical ``restore_gated_redsea`` outputs are not accepted even
through an opt-in. :func:`require_current_pheno` is the refusal, and the identity-dependent
stages call it directly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import pandas as pd

from ..config import PipelineConfig, load_config
from ..sections import SectionError
from ..sections import load_overrides as load_section_overrides
from ..sections import parse as parse_section
from ..sections import resolve as resolve_section
from .contract import (
    as_feature, coerce_cell_table, partition_path, write_cell_table,
)


class ExportError(RuntimeError):
    """A phenotype export was requested without a verified workflow."""


NO_VERIFIED_WORKFLOW_MESSAGE = (
    "PhenoCycler phenotype export is unavailable because no marker-pair or "
    "typing workflow is currently verified. Historical phenotype/broad and "
    "restore_gated_redsea outputs are retired and have no compatibility opt-in. "
    "Use the upstream QuPath input diagnostics while pair design is in development."
)


def require_current_pheno(
    operation: str = "phenotype export",
    *,
    cfg: Optional[PipelineConfig] = None,
):
    """Return a validated declared release, or fail before any cache can be read."""

    if cfg is not None:
        from .phenospace import PhenotypeSpaceError, require_declared_phenotype

        try:
            return require_declared_phenotype(cfg, "phenocycler")
        except PhenotypeSpaceError as exc:
            detail = str(exc)
    else:
        detail = "no pipeline configuration/release was supplied"
    raise ExportError(
        f"{NO_VERIFIED_WORKFLOW_MESSAGE} Requested operation: {operation}. "
        f"Cached files are ignored and never reusable as current inputs. {detail}"
    )


def require_current_identity(
    operation: str,
    *,
    cfg: PipelineConfig,
    modalities: tuple[str, ...] = ("phenocycler", "xenium"),
):
    """Validate every identity release an operation consumes before reading caches."""

    from .phenospace import PhenotypeSpaceError, require_declared_phenotype

    releases = {}
    for modality in modalities:
        if modality == "phenocycler":
            releases[modality] = require_current_pheno(operation, cfg=cfg)
            continue
        try:
            releases[modality] = require_declared_phenotype(cfg, modality)
        except PhenotypeSpaceError as exc:
            raise ExportError(
                f"{modality.capitalize()} identity is unavailable for {operation}: {exc}. "
                "Cached files are ignored; configure a validated release."
            ) from exc
    return releases


def _refuse() -> None:
    require_current_pheno()


def export_donor(
    donor_id: str,
    cfg: PipelineConfig,
) -> pd.DataFrame:
    """Refuse donor phenotype export before reading any artifact."""

    del donor_id, cfg
    _refuse()


def run_export_pheno(*_args, **_kwargs) -> pd.DataFrame:
    """Refuse cohort phenotype export before section discovery."""

    _refuse()



#: Non-marker columns in ``expression``. Everything else in that table is a marker intensity.
_EXPRESSION_META = (
    "donor_id", "object_id", "image", "X_centroid", "Y_centroid",
    "cell_region", "islet_num", "qc_analysis_eligible", "qc_estimation_eligible",
)

# Measured, but unsafe in a geometry profile: structures.py treats their presence as
# permission to apply the retired >=5 hormone-positivity threshold.
_IDENTITY_ONLY_MARKERS = frozenset({"INS", "GCG", "SST"})


def _run_root(cfg: PipelineConfig, run_id: Optional[str] = None) -> Path:
    """Resolve the stage root to export from.

    Defaults to the flat ``data/`` tree. An explicit ``run_id`` still reads one of the
    retired content-addressed runs under ``data/runs/``, which is how the 20-donor
    cohort output from before 2026-09-16 stays reachable.
    """
    runs = cfg.runs_dir
    if not run_id:
        return cfg.data_dir
    if run_id:
        root = runs / run_id
        if not root.is_dir():
            raise ExportError(f"run {run_id!r} not found under {runs}")
        return root
    latest = runs / "LATEST"
    if not latest.exists():
        raise ExportError(f"no run_id given and {latest} does not exist")
    named = latest.read_text().strip()
    root = runs / named
    if not root.is_dir():
        raise ExportError(f"{latest} names {named!r}, which is not a directory under {runs}")
    return root


def _read_donor_stage(root: Path, stage: str, donor_id: str,
                      columns: Optional[list[str]] = None) -> pd.DataFrame:
    part = root / stage / f"donor_id={donor_id}"
    files = sorted(part.rglob("*.parquet"))
    if not files:
        raise ExportError(f"no {stage} parquet for donor {donor_id} under {part}")
    return pd.concat([pd.read_parquet(f, columns=columns) for f in files], ignore_index=True)


def _resolved_section(cfg: PipelineConfig, expr: pd.DataFrame, donor_id: str):
    """Resolve the expression partition's one image through the tracked source manifest."""

    if "image" not in expr.columns:
        raise ExportError(f"donor {donor_id}: expression carries no image provenance")
    images = sorted(expr["image"].dropna().astype(str).str.strip().unique())
    if len(images) != 1:
        raise ExportError(
            f"donor {donor_id}: expression carries {len(images)} image ids {images}; one "
            "partition must be one physical section"
        )
    try:
        payload = json.loads(Path(cfg.qupath_manifest).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ExportError(f"cannot read QuPath manifest {cfg.qupath_manifest}: {exc}") from exc
    image = images[0]
    matches = [
        row for row in payload.get("images", [])
        if str(row.get("image_id", "")).strip() == image
    ]
    if len(matches) != 1:
        raise ExportError(
            f"donor {donor_id}: image {image!r} has {len(matches)} exact rows in "
            f"{cfg.qupath_manifest}; expected one"
        )
    try:
        section = resolve_section(
            image,
            donor_id=str(matches[0].get("donor_id", "")).strip(),
            overrides=load_section_overrides(cfg.section_overrides_csv),
        )
    except SectionError as exc:
        raise ExportError(f"donor {donor_id}: {exc}") from exc
    if section.donor_id != str(donor_id):
        raise ExportError(
            f"run partition donor {donor_id} resolves through the QuPath manifest to "
            f"donor {section.donor_id}; refusing to re-key cells"
        )
    return section


def export_geometry(
    cfg: PipelineConfig,
    donor_id: str,
    *,
    run_id: Optional[str] = None,
    roi: Optional[str] = None,
    analysis_eligible_only: bool = True,
) -> pd.DataFrame:
    """Build a contract cell table from a run's ``expression`` stage. No identity.

    ``expression`` carries coordinates, ``cell_region``, ``islet_num``, both QC flags and one
    authoritative intensity per marker — but **no area**, which the contract wants as
    ``area_um2``. That comes from the same run's ``cells`` stage (``cell_area``), joined on
    ``object_id``. The join is on the QuPath UUID and nothing else: no centroid matching, no
    donor-name parsing. A cell with no area keeps a null rather than a zero, because a missing
    measurement is unavailable, not zero.

    ``lineage_native`` is left absent. ``coerce_cell_table`` fills it with the empty string,
    which is what :func:`~.contract.has_identity` reads as "no labels", so every
    identity-dependent consumer refuses this table on its own evidence rather than on a list.

    ``qc_analysis_eligible`` is carried as a nullable boolean audit column whether or not
    the table was filtered on it, so the same schema comes out either way and a reader can
    always tell which cells the run's geometry QC admitted. A null flag stays null: an
    unrecorded QC decision is unavailable, not a pass and not a fail.
    """
    root = _run_root(cfg, run_id)
    donor = str(donor_id)

    expr = _read_donor_stage(root, "expression", donor)
    section = _resolved_section(cfg, expr, donor)
    if roi is not None and roi != section.roi:
        raise ExportError(
            f"donor {donor}: requested roi {roi!r}, but source image resolves to "
            f"{section.roi!r}"
        )
    resolved_roi = section.roi
    markers = [
        c for c in expr.columns
        if c not in _EXPRESSION_META and c not in _IDENTITY_ONLY_MARKERS
    ]
    if not markers:
        raise ExportError(f"expression for donor {donor} carries no marker columns")

    if analysis_eligible_only and "qc_analysis_eligible" in expr.columns:
        keep = expr["qc_analysis_eligible"].fillna(False).astype(bool)
        expr = expr.loc[keep].reset_index(drop=True)
        if expr.empty:
            raise ExportError(f"donor {donor}: no analysis-eligible cells in expression")

    cells = _read_donor_stage(root, "cells", donor, columns=["object_id", "cell_area"])
    area = (cells.drop_duplicates("object_id")
                 .set_index("object_id")["cell_area"]
                 .reindex(expr["object_id"].astype(str)))

    out = pd.DataFrame({
        "cell_id": expr["object_id"].astype(str),
        "donor_id": donor,
        "roi": resolved_roi,
        "x_um": pd.to_numeric(expr["X_centroid"], errors="coerce"),
        "y_um": pd.to_numeric(expr["Y_centroid"], errors="coerce"),
        "area_um2": pd.to_numeric(area.to_numpy(), errors="coerce"),
        "region": expr["cell_region"].fillna("").astype(str) if "cell_region" in expr else "",
        # Audit/registration context only. S3 is the sole producer of structure_id.
        "qupath_islet_num": (expr["islet_num"].fillna("").astype(str)
                             if "islet_num" in expr else ""),
    })
    for m in markers:
        out[as_feature(m)] = pd.to_numeric(expr[m], errors="coerce").to_numpy()

    # Audit context, never a filter applied on the reader's behalf: which cells the run's
    # geometry QC admitted to analysis. Nullable so an unrecorded decision stays unrecorded.
    if "qc_analysis_eligible" in expr.columns:
        out["qc_analysis_eligible"] = pd.array(
            expr["qc_analysis_eligible"].to_numpy(dtype=object), dtype="boolean"
        )

    # Provenance: which immutable run these coordinates came from. Cheap to carry, and the
    # alternative is a table nobody can trace back to a content-addressed source.
    if "image" in expr.columns:
        out["image"] = expr["image"].astype(str).to_numpy()

    n_missing_area = int(out["area_um2"].isna().sum())
    if n_missing_area:
        print(f"[export_geometry] donor {donor}: {n_missing_area} of {len(out)} cells have no "
              f"cell_area in the run's cells stage; left null, not zeroed.", flush=True)

    return coerce_cell_table(
        out,
        modality="phenocycler",
        donor_id=donor,
        profile="geometry",
        roi=resolved_roi,
    )


def run_export_geometry(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    run_id: Optional[str] = None,
    roi: Optional[str] = None,
    analysis_eligible_only: bool = False,
) -> list[Path]:
    """Export every requested donor into ``cells_pheno_dir``, attaching a declared release.

    The stage default is the COMPLETE cell set (``analysis_eligible_only=False``), with the
    QC flag carried as a column -- see the module docstring for why. A configured PhenoCycler
    release is attached by exact ``object_id`` for each donor it declares; the attach refuses
    (``ExportError``) rather than drop or invent a cell, so a release that covers only part
    of a partition stops the stage instead of silently narrowing the frame. Donors the
    release does not declare are written geometry only and said so.
    """
    from .phenospace import load_release_or_none

    root = _run_root(cfg, run_id)
    available = sorted(p.name.split("=", 1)[1]
                       for p in (root / "expression").glob("donor_id=*") if p.is_dir())
    wanted = [str(d) for d in donors] if donors else available
    unknown = [d for d in wanted if d not in available]
    if unknown:
        raise ExportError(f"run {root.name} has no expression for donor(s) {unknown}")

    release = load_release_or_none(cfg, "phenocycler")
    if release is None:
        print("[export_geometry] no PhenoCycler identity release available; every donor is "
              "exported geometry only", flush=True)

    written: list[Path] = []
    for donor in wanted:
        df = export_geometry(cfg, donor, run_id=run_id, roi=roi,
                             analysis_eligible_only=analysis_eligible_only)
        if "qc_analysis_eligible" in df.columns:
            eligible = df["qc_analysis_eligible"]
            n_ineligible = int((eligible == False).sum())  # noqa: E712 - nullable boolean
            n_unknown = int(eligible.isna().sum())
        else:
            n_ineligible = n_unknown = 0

        if release is not None and donor in set(release.donors):
            df = attach_identity(df, release, donor)
            identity_note = f"identity attached from release {release.release_id}"
        elif release is not None:
            identity_note = (f"geometry only: release {release.release_id} does not declare "
                             f"donor {donor}")
        else:
            identity_note = "geometry only"

        resolved_roi = str(df["roi"].iloc[0])
        path = write_cell_table(df, partition_path(cfg.cells_pheno_dir, donor, resolved_roi))
        written.append(path)
        print(f"[export_geometry] {donor}: {len(df):,} cells ({n_ineligible:,} QC-ineligible "
              f"kept with the flag, {n_unknown:,} unrecorded) -> {path}; {identity_note}",
              flush=True)
    return written


def attach_identity(
    geometry: pd.DataFrame,
    release,
    donor_id: str,
) -> pd.DataFrame:
    """Attach one declared PhenoCycler release by exact QuPath UUID."""

    from .phenospace import PhenotypeSpaceError, attach_release_identity

    if release.modality != "phenocycler":
        raise ExportError(f"cannot attach {release.modality!r} release to PhenoCycler geometry")
    try:
        return attach_release_identity(geometry, release, donor_id)
    except PhenotypeSpaceError as exc:
        # Preserve the exporter-facing exception type used by the stage/CLI.
        raise ExportError(str(exc)) from exc


def _sections_in(
    cfg: PipelineConfig,
    wanted_rois: set[str],
    keep_donors: set[str] | None,
) -> list[str]:
    """Retain non-phenotyping section discovery for integration manifests."""

    output: list[str] = []
    for key in cfg.discover_sections(cfg.broad_dir):
        try:
            section = parse_section(key)
        except SectionError as exc:
            print(f"[export_pheno] SKIPPING partition {key!r}: {exc}", flush=True)
            continue
        if wanted_rois and section.roi not in wanted_rois:
            # Loudly, not silently. A section on disk whose ROI is outside the configured
            # enumeration is almost always a lymph node past `[integration] max_pln`, and
            # dropping it here made it reappear downstream as "no registration produced" —
            # an enumeration failure wearing a registration failure's clothes.
            print(f"[export_pheno] SKIPPING section {section.key!r}: roi {section.roi!r} is outside "
                  f"the configured scope {sorted(wanted_rois)}. If this is a real section, "
                  f"raise [integration] max_pln (currently {cfg.max_pln}).", flush=True)
            continue
        if keep_donors and section.donor_id not in keep_donors:
            continue
        output.append(section.key)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--donor", action="append", dest="donors")
    args = parser.parse_args(argv)
    run_export_pheno(load_config(args.config), args.donors)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "ExportError",
    "NO_VERIFIED_WORKFLOW_MESSAGE",
    "attach_identity",
    "export_donor",
    "export_geometry",
    "require_current_pheno",
    "require_current_identity",
    "run_export_geometry",
    "run_export_pheno",
]
