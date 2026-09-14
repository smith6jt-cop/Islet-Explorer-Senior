"""Which donors this build covers, and why.

The batch-2 export carries four defects (duplicated islet detections, missed islets, inflated
annotation counts, absent galectin regions), and three batch-1 donors have a small islet deficit of
their own.  Rather than assert a cohort, this records the measurement that produced it.

Missed islets were counted as clusters of hormone-high cells lying outside every islet annotation,
swept over six settings -- linking distance 15/25/40 um x within-donor hormone percentile 99/99.5.
The *counts* are parameter-sensitive; membership of the zero set is not.

    donor   settings non-zero   max missed
    115,6374,6436,6442,6461,6476,6482,6493,6539,6593        0 of 6      0
    6521                                                    2 of 6      2   (of 395 islets)
    6551                                                    2 of 6      3   (of 424)
    6380                                                    3 of 6      6   (of  67)
    every batch-2 donor                                     non-zero    10-42

So `CLEAN` is the set that survives every setting AND carries no duplicated detections.
"""
from __future__ import annotations

#: Zero missed islets under all six settings, zero duplicated islet-core detections.
CLEAN: tuple[str, ...] = ("115", "6374", "6436", "6442", "6461", "6476",
                          "6482", "6493", "6539", "6593")

#: Batch-1 donors with a small, real islet deficit -- includable at a stated cost, not clean.
#: 6521 is the only Aab+ donor among them, so adding it is the only way to get an Aab+ donor
#: into a build that excludes batch 2.
MARGINAL: dict[str, str] = {
    "6521": "1-2 missed islets of 395 (0.3-0.5%); the only Aab+ donor outside batch 2",
    "6551": "2-3 missed islets of 424 (0.5-0.7%)",
    "6380": "6 missed islets of 67 (9.0%), 13% of its islet-core cells",
}

#: Every defect in the batch-2 export is donor-wide; none of these is usable until re-exported.
BATCH2_BLOCKED: tuple[str, ...] = ("6450", "6505", "6523", "6534", "6566", "6591", "6623")


def resolve(donors: list[str] | None = None, include_marginal: bool = False) -> list[str]:
    """The donor list to build over, refusing anything knowingly broken."""
    if donors:
        blocked = sorted(set(donors) & set(BATCH2_BLOCKED))
        if blocked:
            raise ValueError(
                f"batch-2 donors are blocked until the QuPath re-export: {blocked}. "
                "Their islet detections are duplicated and their islet annotations incomplete.")
        return list(donors)
    return list(CLEAN) + (sorted(MARGINAL) if include_marginal else [])


#: What this cohort cannot answer, stated once so no figure has to discover it.
COHORT_LIMITS = """
CLEAN is 8 ND + 2 T1D and contains NO Aab+ donor, because all three Aab+ donors are either in
batch 2 (6450, 6505) or carry a small islet deficit (6521).  So no disease-progression comparison
is available on this cohort -- only ND against T1D, on two T1D donors.  Adding 6521 buys one Aab+
donor at a 0.3-0.5% islet deficit; that is a maintainer decision, not a default.
"""
