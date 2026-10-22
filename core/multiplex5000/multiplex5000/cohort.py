"""Cohort eligibility rules shared by every phenocycler analysis stage.

THIS MODULE IS THE ONLY AUTHORITY FOR THE EXCLUSION LIST.  Derive from
:data:`DONOR_EXCLUSIONS`; no lane may restate it from a literal, and configuration
must not bypass it.  As of 2026-09-18 exactly ONE donor is centrally excluded --
6457, for a technical reason -- and :data:`REINSTATED` records 6579's lifting.
"""

from __future__ import annotations

from collections.abc import Iterable


DONOR_EXCLUSIONS: dict[str, str] = {
    "6457": (
        "Excluded until further notice: poor DAPI staining prevented accurate "
        "Phenocycler channel registration."
    ),
}

#: What each exclusion covers, by platform. 6457's reason is a PhenoCycler one -- its DAPI defeated
#: the PhenoCycler channel registration -- and the maintainer confirmed on 2026-10-22 that "the
#: xenium for 6457 is not excluded": its Xenium acquisitions, and the IF scanned on its Xenium
#: slide, stand. An exclusion not listed here covers every platform (fail closed), and every
#: function below still refuses the donor outright: a lane that reads a platform's data asks
#: :func:`excluded_on` instead.
EXCLUDED_PLATFORMS: dict[str, tuple[str, ...]] = {
    "6457": ("phenocycler",),
}


def excluded_on(donor: str, platform: str) -> bool:
    """Whether a donor's acquisitions on one platform (``phenocycler``, ``xenium``, ``if``) are
    centrally excluded."""
    donor = donor_of(donor)
    if donor not in DONOR_EXCLUSIONS:
        return False
    return platform in EXCLUDED_PLATFORMS.get(donor, (platform,))


#: Reinstated 2026-09-17 by maintainer decision. 6579 was excluded as "a pancreatitis outlier
#: [that] must not be included in analysis", "until further notice" -- and this is that notice.
#:
#: The pancreatitis is not retracted, only the exclusion. Anything that pools 6579 with the rest
#: is pooling a donor with a different disease process, so a result that moves when 6579 joins
#: has told you about 6579 rather than about the cohort. Stratify or report it separately; do not
#: treat the larger n as free.
#:
#: 6457 is NOT reinstated -- its exclusion is technical (poor DAPI prevented accurate channel
#: registration), not a judgement call, and nothing has changed about that.
REINSTATED: dict[str, str] = {
    "6579": "pancreatitis outlier; exclusion lifted 2026-09-17, the disease context is not",
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
