"""Path contracts shared by the local Xenium workflow entry points.

The absorbed Xenium lane historically assumed a HiPerGator-only flat layout::

    processed/<serial>/<serial>_phenotyped.h5ad

The combined repository's builder has to distinguish ROIs, so its authoritative
layout is nested instead::

    processed/<serial>/<roi>/<serial>_phenotyped.h5ad

Entry points therefore accept explicit paths and validate the *entire* input set
before doing any expensive work or modifying an h5ad.  The old ``--samples``
resolution remains available for the original lane.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable


PHENOTYPED_SUFFIX = "_phenotyped.h5ad"


def require_existing_files(paths: Iterable[Path | str], *, context: str) -> list[Path]:
    """Return normalized paths, or fail once with every missing required file."""
    normalized = [Path(path).expanduser() for path in paths]
    missing = [path for path in normalized if not path.is_file()]
    if missing:
        rendered = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"{context}: required input file(s) not found:\n{rendered}")
    return normalized


def sample_from_phenotyped_path(path: Path | str) -> str:
    """Extract a serial from the workflow's exact phenotyped-base filename."""
    path = Path(path)
    if not path.name.endswith(PHENOTYPED_SUFFIX):
        raise ValueError(
            f"expected a '*{PHENOTYPED_SUFFIX}' input, got {path}"
        )
    sample = path.name[: -len(PHENOTYPED_SUFFIX)]
    if not sample:
        raise ValueError(f"phenotyped input has no sample serial: {path}")
    return sample


def resolve_phenotyping_jobs(
    *,
    samples: Iterable[str],
    paths: Iterable[Path | str] | None,
    data_root: Path,
    context: str = "Xenium hierarchical phenotyping",
) -> list[tuple[str, Path]]:
    """Resolve and preflight hierarchical-typing inputs.

    Explicit paths support the combined repository's ROI nesting.  Without them,
    the historical flat layout is retained for compatibility.
    """
    if paths is not None:
        candidates = [Path(path).expanduser() for path in paths]
        jobs = [(sample_from_phenotyped_path(path), path) for path in candidates]
    else:
        jobs = [
            (str(sample), data_root / "processed" / str(sample) / f"{sample}{PHENOTYPED_SUFFIX}")
            for sample in samples
        ]

    if not jobs:
        raise ValueError("at least one Xenium phenotyped input is required")

    duplicates = sorted(
        sample for sample, count in Counter(sample for sample, _ in jobs).items() if count > 1
    )
    if duplicates:
        raise ValueError(
            "one input per sample is required; duplicate serial(s): "
            + ", ".join(duplicates)
        )

    checked = require_existing_files(
        [path for _, path in jobs], context=context
    )
    return [(sample, path) for (sample, _), path in zip(jobs, checked, strict=True)]


def resolve_insulitis_jobs(
    *,
    samples: Iterable[str],
    paths: Iterable[Path | str] | None,
    data_root: Path,
) -> list[tuple[str, Path, Path]]:
    """Resolve each phenotyped base and its same-directory immune sidecar."""
    base_jobs = resolve_phenotyping_jobs(
        samples=samples,
        paths=paths,
        data_root=data_root,
        context="Xenium insulitis analysis",
    )
    jobs = [
        (
            sample,
            phenotyped,
            phenotyped.with_name(f"{sample}_immune_phenotyped.h5ad"),
        )
        for sample, phenotyped in base_jobs
    ]
    require_existing_files(
        [immune for _, _, immune in jobs], context="Xenium insulitis analysis"
    )
    return jobs
