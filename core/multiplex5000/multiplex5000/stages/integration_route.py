"""Route the Integration stage to the method the project names, never by file location.

``settings.integration_route: maldi_serial`` places each registered MALDI section's spots on the
neighbouring slide's objects (``maldi_serial_integration``); a project that names no route keeps the
same-slide pilot (``same_slide_integration``), as before this router existed.
"""

from __future__ import annotations

from ..tasks import Task

STAGE = "08_integration"
ROUTES = {"maldi_serial": "maldi_serial_integration", None: "same_slide_integration"}


class IntegrationRouteError(ValueError):
    """A project names an Integration route that is not built."""


def _route(project):
    import importlib

    name = project.settings.get("integration_route")
    if name not in ROUTES:
        raise IntegrationRouteError(f"unknown integration_route: {name}")
    return importlib.import_module(f"{__package__}.{ROUTES[name]}")


def settings(project) -> dict:
    return _route(project).settings(project)


def inputs(project) -> list:
    return _route(project).inputs(project)


def run(task: Task) -> dict:
    return _route(task.project).run(task)


def gate(task: Task, summary: dict):
    return _route(task.project).gate(task, summary)
