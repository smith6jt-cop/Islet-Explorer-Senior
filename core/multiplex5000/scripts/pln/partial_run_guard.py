"""Compatibility imports for the package-owned partial-run guard."""

from multiplex5000.partial_run_guard import (
    PARTIAL_DIR,
    RunScope,
    refuse_shrink,
    scope_slug,
    write_scope_note,
)

__all__ = ("PARTIAL_DIR", "RunScope", "refuse_shrink", "scope_slug", "write_scope_note")
