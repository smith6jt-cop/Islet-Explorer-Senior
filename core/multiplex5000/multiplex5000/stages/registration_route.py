"""Route the Registration stage to the method the project names, never by file location.

``settings.registration_route: maldi_serial`` registers each inspected MALDI section onto the
slide it neighbours (``maldi_serial_registration``); a project that names no route keeps the
same-slide pilot (``same_slide_registration``), as before this router existed.
"""

from __future__ import annotations

from ..tasks import Task

STAGE = "07_registration"
ROUTES = {"maldi_serial": "maldi_serial_registration", None: "same_slide_registration"}


class RegistrationRouteError(ValueError):
    """A project names a Registration route that is not built."""


def _route(project):
    import importlib

    name = project.settings.get("registration_route")
    if name not in ROUTES:
        raise RegistrationRouteError(f"unknown registration_route: {name}")
    return importlib.import_module(f"{__package__}.{ROUTES[name]}")


def settings(project) -> dict:
    return _route(project).settings(project)


def inputs(project) -> list:
    return _route(project).inputs(project)


def run(task: Task) -> dict:
    return _route(task.project).run(task)


def gate(task: Task, summary: dict):
    return _route(task.project).gate(task, summary)
