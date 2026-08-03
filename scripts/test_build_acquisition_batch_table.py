"""Pin the qptiff filename contract in build_acquisition_batch_table.

The regex is the only thing standing between a mis-named export and a fabricated donor id in the
acquisition-batch table, so the accept / skip / reject split is asserted here rather than left to
review. Every name below is a real one, taken from the authoritative donor -> filename array in
`senior/convert_viewer_ometiff.sh`.

The module is imported when its dependencies are present; when pandas or the Phenocycler submodule
are missing the constants are read straight out of the source instead, so the contract is still
enforced in a bare checkout.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent
SCRIPT = SCRIPTS / "build_acquisition_batch_table.py"


def _constants() -> tuple[str, str]:
    sys.path.insert(0, str(SCRIPTS))
    try:
        import build_acquisition_batch_table as mod
        return mod.QPTIFF_NAME, mod.NON_COHORT_NAME
    except Exception:                      # pandas / phenocycler submodule unavailable
        wanted = {"QPTIFF_NAME", "NON_COHORT_NAME"}
        found = {}
        for node in ast.parse(SCRIPT.read_text()).body:
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                if node.targets[0].id in wanted:
                    found[node.targets[0].id] = ast.literal_eval(node.value)
        missing = wanted - set(found)
        if missing:
            raise AssertionError(f"constants not found in {SCRIPT.name}: {sorted(missing)}")
        return found["QPTIFF_NAME"], found["NON_COHORT_NAME"]


QPTIFF_NAME, NON_COHORT_NAME = _constants()

# donor id -> exact export name, mirroring senior/convert_viewer_ometiff.sh.
COHORT = [
    ("6374", "6374_Scan1.er.qptiff"), ("6380", "6380_Scan1.er.qptiff"),
    ("6436", "6436_Scan1.er.qptiff"), ("6442", "6442_Scan1.er.qptiff"),
    ("6457", "6457_Scan1.er.qptiff"), ("6461", "6461_Scan1.er.qptiff"),
    ("6482", "6482_Scan1.er.qptiff"), ("6493", "6493_Scan1.er.qptiff"),
    ("6521", "6521_Scan1.er.qptiff"), ("6539", "6539_Scan1.er.qptiff"),
    ("6551", "6551_Scan1.er.qptiff"), ("6593", "6593_Scan1.er.qptiff"),
    # rescans
    ("6476", "6476_Scan2.er.qptiff"), ("6579", "6579_Scan2.er.qptiff"),
    # batch 2 carries a tissue suffix, with one typo on 6534
    ("6450", "6450panc_Scan1.er.qptiff"), ("6505", "6505panc_Scan1.er.qptiff"),
    ("6523", "6523panc_Scan1.er.qptiff"), ("6566", "6566panc_Scan1.er.qptiff"),
    ("6591", "6591panc_Scan1.er.qptiff"), ("6623", "6623panc_Scan1.er.qptiff"),
    ("6534", "6534pacn_Scan1.er.qptiff"),
]

# Names the old loose pattern silently mis-parsed into a donor id.
MISPARSED = [
    ("Scan1.qptiff", "1"),
    ("Scan1.er.qptiff", "1"),
    ("HDL115pancLN_Scan1.er.qptiff", "115"),
]


@pytest.mark.parametrize("donor,name", COHORT, ids=[n for _, n in COHORT])
def test_cohort_names_parse_to_their_donor(donor, name):
    m = re.fullmatch(QPTIFF_NAME, name)
    assert m is not None, f"{name} must be accepted as a cohort scan"
    assert m.group(1) == donor


def test_every_batch2_donor_is_covered():
    """The panc/pacn suffix is what distinguishes batch 2; none of them may fall through."""
    batch2 = {"6450", "6505", "6523", "6534", "6566", "6591", "6623"}
    parsed = {re.fullmatch(QPTIFF_NAME, n).group(1) for d, n in COHORT if d in batch2}
    assert parsed == batch2


def test_lymph_node_is_skipped_not_parsed():
    name = "HDL115pancLN_Scan1.er.qptiff"
    assert re.fullmatch(NON_COHORT_NAME, name), "the HDL115 lymph node must be recognised and skipped"
    assert not re.fullmatch(QPTIFF_NAME, name), "it must never parse as a cohort donor"


@pytest.mark.parametrize("name", [
    "control.qptiff", "6450.qptiff", "Scan1.qptiff", "Scan1.er.qptiff",
    "6450_Scan1.tiff", "_6450_Scan1.er.qptiff", "6450panc_Scan1.er.qptiff.bak",
])
def test_unrecognised_names_are_rejected(name):
    assert not re.fullmatch(QPTIFF_NAME, name)
    assert not re.fullmatch(NON_COHORT_NAME, name)


@pytest.mark.parametrize("name,old_donor", MISPARSED, ids=[n for n, _ in MISPARSED])
def test_previously_misparsed_names_no_longer_yield_a_donor(name, old_donor):
    """The old pattern matched, so the fail-closed guard never fired on these."""
    assert re.match(r"(?:[A-Za-z]+)?(\d+)", name).group(1) == old_donor   # what it used to do
    assert not re.fullmatch(QPTIFF_NAME, name)                            # and no longer does


def test_no_cohort_name_is_treated_as_non_cohort():
    for _, name in COHORT:
        assert not re.fullmatch(NON_COHORT_NAME, name), f"{name} must not be skipped"
