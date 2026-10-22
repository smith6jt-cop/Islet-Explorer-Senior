"""Cohort exclusions must hold for a donor's SECTIONS, not just its bare id.

The cohort is one pancreas plus 1-4 lymph nodes per donor, and every one of those is its own
section key (`6579`, `6579pln`, `6579pln2`, ...). Those keys are what the manifest, the
partitions and the stage loops actually carry, so an exclusion that only recognises `6579`
blocks the identifier nobody passes and admits the ones that carry the cells.
"""

from __future__ import annotations

import pytest

from phenocycler.cohort import (DONOR_EXCLUSIONS, donor_of, ensure_eligible_donors,
                                filter_eligible_donors)


def test_donor_of_resolves_sections_to_their_donor():
    assert donor_of("6579") == "6579"
    assert donor_of("6579pln") == "6579"
    assert donor_of("6579pln2") == "6579"
    assert donor_of("6579xenpln") == "6579"      # post-Xenium re-stain of a lymph node
    assert donor_of("6539pln3") == "6539"


def test_donor_of_passes_through_what_it_cannot_parse():
    """A caller handing in plain ids, or something odd, must not be silently rewritten."""
    assert donor_of("not-a-section-key") == "not-a-section-key"
    assert donor_of("") == ""


#: Built from DONOR_EXCLUSIONS rather than written out, so reinstating or excluding a donor
#: needs no edit here. 6579 was excluded "until further notice" and reinstated 2026-09-17; the
#: rule this file defends -- every SECTION of an excluded donor is refused, not just the bare id
#: -- is unchanged by that, and hardcoding the donor made a cohort decision look like a
#: regression in twelve tests.
EXCLUDED_KEYS = [k for d in DONOR_EXCLUSIONS
                 for k in (d, f"{d}pln", f"{d}pln2", f"{d}pln4", f"{d}xenpln")]


@pytest.mark.parametrize("key", EXCLUDED_KEYS)
def test_every_section_of_an_excluded_donor_is_refused(key):
    """The regression this file exists for.

    Exclusion matched the raw string, so `6579pln2` — one of donor 6579's lymph nodes —
    passed both gates. Under section-keyed partitions that is every section of an excluded
    donor being readmitted, silently, while the bare id nobody passes stays blocked.
    """
    with pytest.raises(ValueError, match="excluded donor"):
        ensure_eligible_donors([key], context="test")
    assert filter_eligible_donors([key]) == []


def test_the_error_names_both_the_key_and_the_donor():
    """`<donor>pln2` does not look like `<donor>` at a glance; the message has to bridge that."""
    donor = sorted(DONOR_EXCLUSIONS)[0]
    key = f"{donor}pln2"
    with pytest.raises(ValueError) as excinfo:
        ensure_eligible_donors([key], context="test")
    message = str(excinfo.value)
    assert key in message and donor in message
    assert DONOR_EXCLUSIONS[donor].split(":")[0] in message


def test_eligible_sections_survive():
    keys = ["6539", "6539pln", "6539pln2", "6539pln3", "6539pln4"]
    assert filter_eligible_donors(keys) == keys
    assert ensure_eligible_donors(keys, context="test") == tuple(keys)


def test_a_mixed_cohort_reports_only_the_excluded_sections():
    bad = f"{sorted(DONOR_EXCLUSIONS)[0]}pln"
    with pytest.raises(ValueError) as excinfo:
        ensure_eligible_donors(["6539", "6539pln", bad], context="test")
    assert bad in str(excinfo.value)
    assert "6539pln" not in str(excinfo.value)


def test_6457s_exclusion_covers_its_phenocycler_data_only():
    """The maintainer, 2026-10-22: "the xenium for 6457 is not excluded" -- its reason is the
    PhenoCycler channel registration its DAPI defeated."""
    from multiplex5000.cohort import excluded_on

    assert excluded_on("6457", "phenocycler") and excluded_on("6457pln", "phenocycler")
    assert not excluded_on("6457", "xenium") and not excluded_on("6457", "if")
    assert not excluded_on("6374", "phenocycler")
    with pytest.raises(ValueError, match="excluded donor"):     # donor-wide gates still refuse it
        ensure_eligible_donors(["6457"])


def test_an_exclusion_with_no_declared_scope_covers_every_platform(monkeypatch):
    from multiplex5000 import cohort

    monkeypatch.setitem(cohort.DONOR_EXCLUSIONS, "9999", "a reason")
    assert all(cohort.excluded_on("9999", p) for p in ("phenocycler", "xenium", "if"))


def test_every_declared_scope_belongs_to_an_excluded_donor():
    from multiplex5000.cohort import EXCLUDED_PLATFORMS

    assert set(EXCLUDED_PLATFORMS) <= set(DONOR_EXCLUSIONS)

