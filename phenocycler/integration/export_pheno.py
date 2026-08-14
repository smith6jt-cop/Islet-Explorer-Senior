"""Fail-closed PhenoCycler integration boundary.

There is currently no verified marker-pair or phenotype workflow.  Historical
``restore_gated_redsea`` outputs are not accepted, even through an opt-in.
Upstream QuPath/REDSEA development may continue, but phenotype export remains
unavailable until a new pair roster and its downstream typing are validated.
"""

from __future__ import annotations

import argparse

import pandas as pd

from ..config import PipelineConfig, load_config
from ..sections import SectionError
from ..sections import parse as parse_section


class ExportError(RuntimeError):
    """A phenotype export was requested without a verified workflow."""


NO_VERIFIED_WORKFLOW_MESSAGE = (
    "PhenoCycler phenotype export is unavailable because no marker-pair or "
    "typing workflow is currently verified. Historical phenotype/broad and "
    "restore_gated_redsea outputs are retired and have no compatibility opt-in. "
    "Use the upstream QuPath input diagnostics while pair design is in development."
)


def require_current_pheno(operation: str = "phenotype export") -> None:
    """Fail before a public stage can read cached retired phenotype tables."""

    raise ExportError(
        f"{NO_VERIFIED_WORKFLOW_MESSAGE} Requested operation: {operation}. "
        "Cached files are ignored and never reusable as current inputs."
    )


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
    "export_donor",
    "require_current_pheno",
    "run_export_pheno",
]
