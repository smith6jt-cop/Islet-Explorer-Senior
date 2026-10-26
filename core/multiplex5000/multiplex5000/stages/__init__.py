"""The stages of a release, in order, and the phase that builds each.

A stage is a module with ``STAGE`` (its id), ``inputs(project)`` (every file it reads, as
``tasks.Input``), ``run(task)`` (writes its outputs into ``task.staging`` and returns a summary)
and ``gate(task, summary)`` (its checks and figures). The Snakemake graph and the notebooks both
call these, through ``tasks.run_stage``; neither has a second implementation.

MultiPlex5000 D3 (2026-10-20). Stage 00 is the first one built.
"""

from __future__ import annotations

import importlib
from types import ModuleType

#: Every stage of the plan, in order, with the phase that builds it.
PLAN: tuple[tuple[str, str], ...] = (
    ("00_inventory", "D3"),
    ("01_image_preparation", "D9"),
    ("02_segmentation", "D9"),
    ("03_protein_measurement", "D4"),
    # RNA ingest was D8 in the original code-migration order. The maintainer's
    # 2026-10-26 typing scope requires it before D6 can type Xenium-only and
    # same-slide sections, so it is now a D6 prerequisite.
    ("04_rna_ingest", "D6 prerequisite"),
    ("05_normalization", "D5"),
    ("06_typing", "D6"),
    ("07_registration", "D10"),
    ("08_integration", "D10"),
    ("09_spatial_analysis", "D11"),
    ("10_cohort_analysis", "D11"),
    ("11_release_build", "D7"),
    ("12_release_verification", "D7"),
)

#: The stages that exist as code, by module.
BUILT: dict[str, str] = {
    "00_inventory": "multiplex5000.stages.inventory",
    "03_protein_measurement": "multiplex5000.stages.measurement",
    "04_rna_ingest": "multiplex5000.stages.rna_ingest",
    "05_normalization": "multiplex5000.stages.normalization",
    "06_typing": "multiplex5000.stages.typing",
}


class StageError(LookupError):
    """A stage that is not in the plan, or not built yet."""


def get(stage_id: str) -> ModuleType:
    phases = dict(PLAN)
    if stage_id not in phases:
        raise StageError(f"no stage {stage_id!r}; the plan's stages are {[s for s, _ in PLAN]}")
    if stage_id not in BUILT:
        raise StageError(f"stage {stage_id} is built in phase {phases[stage_id]}, not yet")
    module = importlib.import_module(BUILT[stage_id])
    if getattr(module, "STAGE", None) != stage_id:
        raise StageError(f"{BUILT[stage_id]} does not declare STAGE = {stage_id!r}")
    return module


def built() -> list[str]:
    return [s for s, _ in PLAN if s in BUILT]
