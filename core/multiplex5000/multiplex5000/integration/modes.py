"""
Per-section-pair integration mode.

``[integration] mode`` in ``config.ini`` is a single global switch. That was correct while
every donor in the cohort had the same acquisition design, and it stops being correct the
moment one cohort contains both serial pairs and same-slide pairs. This module makes the
mode a property of the *section pair*, declared as data, so the two can coexist in one run.

**Why this is a manifest and not a flag.** ``--mode`` asserts an acquisition fact. Whether
two images are of one physical section is not a runtime preference — it is a property of how
the slide was handled, and the person typing the command is not a reliable source for it.
The dangerous direction is specific: running *serial* sections as ``same_slide`` claims
whole-cell correspondence for every cell, which serial sections cannot support, and
``sameslide.require_same_slide`` cannot catch it because the flag it checks already says
``same_slide``. Declaring the fact in a manifest, with evidence, removes the flag's ability
to invent it.

**Fail-safe default.** A pair absent from the manifest falls back to the configured global
mode (normally ``sequential``). That default is the safe direction: sequential analysis of a
same-slide pair is merely weaker than it could be, while the reverse fabricates a
measurement. So an omission costs resolution, never correctness.

The manifest is ``data/integration/section_modes.csv``; ``#`` comment lines are ignored.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from .contract import normalize_donor_id

MODES = ("sequential", "same_slide")

#: Columns the manifest must carry. ``evidence`` is required and not decorative: a
#: same_slide declaration is an experimental claim, so it has to say what backs it.
REQUIRED_COLUMNS = ("donor_id", "roi", "mode", "evidence")

#: Optional key column. One donor+roi can hold BOTH designs — donor 6533 has a serial
#: pair on Xenium 0041393 and a same-slide pair on 0080870, both roi 'panc' — so the
#: Xenium serial is part of the key. A blank serial is a wildcard for that donor+roi.
SERIAL_COLUMN = "xenium_serial"

DEFAULT_MANIFEST = Path("data/integration/section_modes.csv")


class ModeManifestError(RuntimeError):
    """The mode manifest is missing, malformed, or contradicted by a flag."""


@dataclass(frozen=True)
class SectionMode:
    """One declared (donor, roi) acquisition design."""

    donor_id: str
    roi: str
    mode: str
    evidence: str
    note: str = ""
    xenium_serial: str = ""        # "" = applies to every section of this donor+roi

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.donor_id, self.roi, self.xenium_serial)


def manifest_path(cfg) -> Path:
    """Resolve the manifest path, honouring an explicit config override."""
    raw = getattr(cfg, "section_modes_csv", None) or DEFAULT_MANIFEST
    return Path(raw).expanduser()


def load_manifest(cfg) -> dict[tuple[str, str], SectionMode]:
    """Read the manifest. A missing file is not an error — it means 'no exceptions'."""
    path = manifest_path(cfg)
    if not path.exists():
        return {}

    raw = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in raw.columns]
    if missing:
        raise ModeManifestError(
            f"{path} is missing required column(s) {missing}; expected at least "
            f"{list(REQUIRED_COLUMNS)} (optional: note)."
        )

    out: dict[tuple[str, str], SectionMode] = {}
    for i, row in raw.iterrows():
        donor = normalize_donor_id(str(row["donor_id"]).strip())
        roi = str(row["roi"]).strip()
        mode = str(row["mode"]).strip()
        evidence = str(row["evidence"]).strip()
        if not donor or not roi:
            raise ModeManifestError(f"{path} row {i + 2}: donor_id and roi are both required.")
        if mode not in MODES:
            raise ModeManifestError(
                f"{path} row {i + 2}: mode must be one of {list(MODES)}, got {mode!r}."
            )
        if mode == "same_slide" and not evidence:
            raise ModeManifestError(
                f"{path} row {i + 2}: a same_slide declaration requires non-empty 'evidence'. "
                f"It asserts that one physical section was imaged twice — record what "
                f"establishes that (run log, slide ID, acquisition order)."
            )
        entry = SectionMode(donor, roi, mode, evidence, str(row.get("note", "")).strip(),
                            str(row.get(SERIAL_COLUMN, "")).strip())
        if entry.key in out and out[entry.key].mode != mode:
            raise ModeManifestError(
                f"{path}: ({donor}, {roi}, {entry.xenium_serial or '*'}) is declared twice "
                f"with conflicting modes "
                f"({out[entry.key].mode!r} and {mode!r})."
            )
        out[entry.key] = entry
    return out


def mode_for(cfg, donor: str, roi: str = "panc",
             serial: Optional[str] = None) -> str:
    """The acquisition mode for one pair: manifest first, configured default otherwise."""
    manifest = load_manifest(cfg)
    donor_n, roi_n = normalize_donor_id(str(donor).strip()), str(roi).strip()
    serial_n = str(serial).strip() if serial else ""
    # Exact section first, then the donor+roi wildcard. Resolving a same-slide serial row
    # requires the caller to say WHICH section; without a serial the wildcard (normally
    # sequential) applies, which is the fail-safe direction this module is built on --
    # it costs resolution, not correctness.
    entry = manifest.get((donor_n, roi_n, serial_n)) if serial_n else None
    if entry is None:
        entry = manifest.get((donor_n, roi_n, ""))
    if entry:
        return entry.mode
    # A pair may only reach same_slide by being DECLARED. Keyed on the manifest FILE
    # existing, not on it having rows: an empty manifest is still an active manifest,
    # and gating on row count would leave the hole open until someone added one. Letting the global default supply
    # it would reopen exactly the hole this module closes: `--mode same_slide` would then
    # apply whole-cell correspondence to every undeclared pair in the run, silently.
    if cfg.integration_mode == "same_slide" and manifest_path(cfg).exists():
        raise ModeManifestError(
            f"donor {donor}/{roi} is not declared in {manifest_path(cfg)}, so it cannot run "
            f"same_slide even though [integration] mode is 'same_slide'. Same-slide status "
            f"is per section pair: add a row with evidence, or run this pair sequential."
        )
    return cfg.integration_mode


def declared_same_slide(cfg) -> list[SectionMode]:
    """Every pair the manifest declares same_slide — the exceptions, in one list."""
    return [e for e in load_manifest(cfg).values() if e.mode == "same_slide"]


def resolve_mode(cfg, donor: str, roi: str = "panc",
                 requested: Optional[str] = None,
                 serial: Optional[str] = None) -> str:
    """Return the pair's mode, refusing a flag that contradicts the manifest.

    ``requested`` is the ``--mode`` value (``None`` when not passed). The flag may *select*
    among declared truth; it may not *override* it, in either direction.
    """
    donor_n = normalize_donor_id(str(donor).strip())
    roi_n = str(roi).strip()
    serial_n = str(serial).strip() if serial else ""
    manifest = load_manifest(cfg)
    entry = manifest.get((donor_n, roi_n, serial_n)) if serial_n else None
    if entry is None:
        entry = manifest.get((donor_n, roi_n, ""))
    declared = entry.mode if entry else None
    effective = declared if declared is not None else cfg.integration_mode

    if requested is None:
        return effective
    if requested not in MODES:
        raise ModeManifestError(f"--mode must be one of {list(MODES)}, got {requested!r}.")
    if requested == effective:
        return effective

    if requested == "same_slide":
        raise ModeManifestError(
            f"--mode same_slide was requested for donor {donor_n} / {roi_n}, but "
            + (f"{manifest_path(cfg)} declares it {declared!r}."
               if declared else
               f"no row in {manifest_path(cfg)} declares it same_slide.")
            + "\nWhether two images are of one physical section is an acquisition fact, not "
              "a runtime option: same-slide mode claims whole-cell correspondence for every "
              "cell, and serial sections cannot support that. Declare the pair in the "
              "manifest with evidence, or run it sequential."
        )
    raise ModeManifestError(
        f"--mode sequential was requested for donor {donor_n} / {roi_n}, but "
        f"{manifest_path(cfg)} declares it same_slide (evidence: "
        f"{entry.evidence if entry else 'n/a'}).\nRunning a same-slide pair sequentially "
        f"discards a real cell-to-cell measurement. If that is genuinely intended, change "
        f"the manifest row rather than the command line, so the choice is recorded."
    )
