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


@pytest.mark.parametrize("key", ["6579", "6579pln", "6579pln2", "6579pln4",
                                 "6579xenpln", "6457", "6457pln"])
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
    """`6579pln2` does not look like `6579` at a glance; the message has to bridge that."""
    with pytest.raises(ValueError) as excinfo:
        ensure_eligible_donors(["6579pln2"], context="test")
    message = str(excinfo.value)
    assert "6579pln2" in message and "6579" in message
    assert DONOR_EXCLUSIONS["6579"].split(":")[0] in message


def test_eligible_sections_survive():
    keys = ["6539", "6539pln", "6539pln2", "6539pln3", "6539pln4"]
    assert filter_eligible_donors(keys) == keys
    assert ensure_eligible_donors(keys, context="test") == tuple(keys)


def test_a_mixed_cohort_reports_only_the_excluded_sections():
    with pytest.raises(ValueError) as excinfo:
        ensure_eligible_donors(["6539", "6539pln", "6579pln"], context="test")
    assert "6579pln" in str(excinfo.value)
    assert "6539pln" not in str(excinfo.value)
