"""Profiles: tissue, platform, panel and cohort facts as data (MultiPlex5000 plan, phase D1).

Today these facts are Python constants scattered across the package (tissue names alone are spelled
seven ways). A profile is one small YAML file under ``core/multiplex5000/profiles/<kind>s/<id>.yaml``
holding the facts for one thing, so a new tissue or platform is a new file, not an edit in seven
places. This module defines the schema and refuses any file that does not meet it.

**Profiles do not yet replace the constants.** Until phase D2 makes them the authority, every value
a profile carries is a *mirror* of the constant it names under ``sources``, and
``tests/current/test_profiles.py`` fails if the two ever disagree -- so there is never a second,
silently diverging copy.

Two kinds of fact are deliberately *referenced*, never restated, because they already have one
authority: a cohort names its membership file and the exclusion constant
(``phenocycler.cohort.DONOR_EXCLUSIONS``); a panel names its marker registry or the panel-role table
in ``core/panels``. A cohort profile that lists donors, or a panel profile that lists markers, is
refused.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

#: The one profile format this module reads.
PROFILE_SCHEMA = "multiplex5000-profile/1"

#: Where the shipped profiles live: ``core/multiplex5000/profiles``.
PROFILES_DIR = Path(__file__).resolve().parents[1] / "profiles"

#: The repository root, for resolving the paths profiles reference.
REPO = Path(__file__).resolve().parents[3]


class ProfileError(ValueError):
    """A profile file that does not meet the schema."""


@dataclass(frozen=True)
class TissueProfile:
    id: str
    display_name: str
    aliases: tuple[str, ...]
    roi_prefixes: tuple[str, ...]
    panel_tissue_dir: Optional[str]
    tissue_only_panels: tuple[str, ...]
    structures: tuple[str, ...]
    comparable_lineages: tuple[str, ...]
    plugin: Optional[str]
    sources: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PlatformProfile:
    id: str
    display_name: str
    pixel_size_um: float
    sources: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PanelProfile:
    id: str
    platform: str
    authority: str          # repository-relative path of the one file that defines the markers
    sources: dict = field(default_factory=dict)


@dataclass(frozen=True)
class CohortProfile:
    id: str
    platform: str
    tissues: tuple[str, ...]
    membership: str         # repository-relative path of the membership authority
    exclusions: str         # dotted name of the exclusion constant; never a list
    sources: dict = field(default_factory=dict)


#: kind -> (class, required keys, optional keys). Anything else in a file is refused.
_KINDS: dict[str, tuple[type, tuple[str, ...], tuple[str, ...]]] = {
    "tissue": (TissueProfile,
               ("id", "display_name", "aliases", "roi_prefixes", "panel_tissue_dir",
                "tissue_only_panels", "structures", "comparable_lineages", "plugin"),
               ("sources",)),
    "platform": (PlatformProfile, ("id", "display_name", "pixel_size_um"), ("sources",)),
    "panel": (PanelProfile, ("id", "platform", "authority"), ("sources",)),
    "cohort": (CohortProfile, ("id", "platform", "tissues", "membership", "exclusions"),
               ("sources",)),
}

_LISTS = {"aliases", "roi_prefixes", "tissue_only_panels", "structures", "comparable_lineages",
          "tissues"}
_NULLABLE = {"panel_tissue_dir", "plugin"}


def load_profile(path: Path) -> Any:
    """Read and check one profile file; return its dataclass."""
    path = Path(path)
    try:
        data = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError) as err:
        raise ProfileError(f"{path}: cannot be read as YAML ({err})") from err
    if not isinstance(data, dict):
        raise ProfileError(f"{path}: a profile is a mapping")
    if data.get("schema") != PROFILE_SCHEMA:
        raise ProfileError(f"{path}: expected schema: {PROFILE_SCHEMA}")
    kind = data.get("kind")
    if kind not in _KINDS:
        raise ProfileError(f"{path}: kind must be one of {sorted(_KINDS)}, not {kind!r}")
    cls, required, optional = _KINDS[kind]
    body = {k: v for k, v in data.items() if k not in ("schema", "kind")}
    unknown = sorted(set(body) - set(required) - set(optional))
    if unknown:
        raise ProfileError(f"{path}: unknown key(s) {unknown} for a {kind} profile")
    missing = [k for k in required if k not in body]
    if missing:
        raise ProfileError(f"{path}: missing key(s) {missing} for a {kind} profile")
    for key, value in body.items():
        if key in _LISTS:
            if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
                raise ProfileError(f"{path}: {key} must be a list of strings")
            body[key] = tuple(value)
        elif key in _NULLABLE:
            if value is not None and not isinstance(value, str):
                raise ProfileError(f"{path}: {key} must be a string or null")
        elif key == "pixel_size_um":
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
                raise ProfileError(f"{path}: pixel_size_um must be a positive number")
            body[key] = float(value)
        elif key == "sources":
            if not isinstance(value, dict):
                raise ProfileError(f"{path}: sources must be a mapping of key -> where it is defined")
        elif not isinstance(value, str) or not value.strip():
            raise ProfileError(f"{path}: {key} must be a non-empty string")
    if body["id"] != path.stem:
        raise ProfileError(f"{path}: id {body['id']!r} must equal the file name {path.stem!r}")
    if path.parent.name != f"{kind}s":
        raise ProfileError(f"{path}: a {kind} profile lives in a '{kind}s' folder")
    if kind in ("panel", "cohort"):
        ref = body["authority" if kind == "panel" else "membership"]
        if Path(ref).is_absolute() or not (REPO / ref).is_file():
            raise ProfileError(f"{path}: {ref} must be a repository-relative file that exists")
    return cls(**body)


def load_all(root: Path = PROFILES_DIR) -> dict[tuple[str, str], Any]:
    """Every profile under ``root``, keyed by (kind, id)."""
    found: dict[tuple[str, str], Any] = {}
    for path in sorted(Path(root).glob("*s/*.yaml")):
        profile = load_profile(path)
        kind = path.parent.name[:-1]
        found[(kind, profile.id)] = profile
    return found
