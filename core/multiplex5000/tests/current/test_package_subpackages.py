"""Tagged wheel installs must contain every maintained Python subpackage."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_every_importable_subpackage_is_in_the_wheel_manifest():
    pyproject = (ROOT / "pyproject.toml").read_text()
    section = pyproject.split("[tool.setuptools]", 1)[1].split(
        "[tool.setuptools.package-data]", 1)[0]
    manifest = re.search(r"(?m)^packages\s*=\s*\[(.*?)\]", section, re.S)
    assert manifest is not None
    declared = set(re.findall(r'"([\w.]+)"', manifest.group(1)))
    actual = {".".join(path.parent.relative_to(ROOT).parts)
              for path in ROOT.rglob("__init__.py")
              if path.parent.relative_to(ROOT).parts[0] in {"multiplex5000", "phenocycler"}}
    assert actual <= declared, f"wheel omits importable subpackages: {sorted(actual - declared)}"
