"""Cohort eligibility rules shared by every phenocycler analysis stage."""

from __future__ import annotations

from collections.abc import Iterable


DONOR_EXCLUSIONS: dict[str, str] = {
    "6457": (
        "Excluded until further notice: poor DAPI staining prevented accurate "
        "Phenocycler channel registration."
    ),
    "6579": (
        "Excluded until further notice by maintainer decision: pancreatitis outlier "
        "must not be included in analysis."
    ),
}


def donor_of(key: str) -> str:
    """The donor a key belongs to, whether it is a donor id or a section key.

    Exclusions are decisions about a *donor*, but the identifiers that reach these gates
    are often *section* keys — ``6579pln2`` is one of donor 6579's lymph nodes, and a
    cohort where each donor contributes a pancreas and 1-4 lymph nodes is made entirely of
    such keys. Matching the raw string would let every section of an excluded donor
    through while blocking only the bare id, which is the wrong way round for a
    fail-closed rule: the sections are the things that actually carry cells.

    Anything that does not parse as a section key is returned unchanged, so a caller
    passing plain donor ids is unaffected.
    """
    from .sections import SectionError, parse

    try:
        return parse(key).donor_id
    except SectionError:
        # Only "this is not a section key" is tolerable: the caller passed a plain donor
        # id, which is returned unchanged. Any other exception is a defect, and swallowing
        # it here would silently readmit a donor the exclusion list is meant to refuse.
        return str(key)


def ensure_eligible_donors(
    donors: Iterable[str],
    *,
    context: str = "analysis",
) -> tuple[str, ...]:
    """Return normalized donor IDs or reject any explicitly excluded donor.

    Rejects a section of an excluded donor as well as the donor itself — see
    :func:`donor_of`.
    """
    normalized = tuple(str(donor) for donor in donors)
    blocked_keys = sorted({key for key in normalized if donor_of(key) in DONOR_EXCLUSIONS})
    if blocked_keys:
        details = "; ".join(
            # Name the key that was passed *and* the donor it resolves to, so a
            # section key does not look like an unrelated identifier in the error.
            f"{key} (donor {donor_of(key)}): {DONOR_EXCLUSIONS[donor_of(key)]}"
            for key in blocked_keys
        )
        raise ValueError(f"{context} includes excluded donor(s): {details}")
    return normalized


def filter_eligible_donors(donors: Iterable[str]) -> list[str]:
    """Remove centrally excluded donors from automatic cohort discovery.

    Drops sections of an excluded donor too — see :func:`donor_of`.
    """
    return [
        str(donor)
        for donor in donors
        if donor_of(donor) not in DONOR_EXCLUSIONS
    ]
