"""Route Stage 06 to the tissue pilot named by the project, never by file location."""

from __future__ import annotations

from ..tasks import Task

STAGE = "06_typing"
REQUIRE_SETTINGS_MATCH_ON_REUSE = True


class TypingRouteError(ValueError):
    """A project has no unique Stage 06 assay/tissue route."""


def _route(project):
    explicit = project.settings.get("typing_route")
    if explicit == "io60_xenium_pilot":
        from . import io60_rna_typing_pilot

        return io60_rna_typing_pilot
    if explicit == "io60_mixed_protein_pilot":
        from . import io60_mixed_protein_pilot

        return io60_mixed_protein_pilot
    if explicit == "io60_joint_typing":
        from . import io60_joint_typing

        return io60_joint_typing
    if explicit is not None:
        raise TypingRouteError(f"unknown Stage 06 typing_route: {explicit}")
    same_slide = bool(project.settings.get("same_slide_source_manifest"))
    dedicated_pln = bool(project.settings.get("pln_export_images"))
    if same_slide == dedicated_pln:
        raise TypingRouteError("Stage 06 needs exactly one recorded tissue pilot route")
    if same_slide:
        if project.cohort != "spleen_same_slide_pilot":
            raise TypingRouteError("same-slide Stage 06 is built only for the spleen pilot")
        from . import spleen_typing_pilot

        return spleen_typing_pilot
    from . import pln_typing_pilot

    return pln_typing_pilot


def settings(project) -> dict:
    return _route(project).settings(project)


def inputs(project) -> list:
    return _route(project).inputs(project)


def run(task: Task) -> dict:
    return _route(task.project).run(task)


def gate(task: Task, summary: dict):
    return _route(task.project).gate(task, summary)
