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
#:
#: The seven re-exported donors joined on 2026-09-15, when the sweep was finally implemented
#: (`missed_islets.py` -- its results were recorded here in 2026-09 and its code never was) and
#: re-run against the new annotations. All seven return zero at every setting, so the re-export
#: fixed the missed islets as well as the duplicated detections, which nothing had checked.
#:
#: 6521 stays MARGINAL despite returning zero in that run. The reimplementation is measurably
#: conservative -- it reproduces the recorded zero set exactly but under-reports the marginal
#: donors, 6521 among them -- so its zero is not strong enough to overturn the original's 1-2.
CLEAN: tuple[str, ...] = ("115", "6374", "6436", "6442", "6450", "6461", "6476", "6482",
                          "6493", "6505", "6523", "6534", "6539", "6566", "6591", "6593",
                          "6623")

#: Batch-1 donors with a small, real islet deficit -- includable at a stated cost, not clean.
#: 6521 is the only Aab+ donor among them, so adding it is the only way to get an Aab+ donor
#: into a build that excludes batch 2.
MARGINAL: dict[str, str] = {
    "6521": "1-2 missed islets of 395 (0.3-0.5%); the only Aab+ donor outside batch 2",
    "6551": "2-3 missed islets of 424 (0.5-0.7%)",
    "6380": "6 missed islets of 67 (9.0%), 13% of its islet-core cells",
}

#: Donors whose ORIGINAL export carries the duplicated islet detections. Being on this list is
#: not what blocks a donor -- the absence of a rebuild is. All seven were re-exported and rebuilt
#: on 2026-09-15, so the guard below now passes them.
BATCH2_REEXPORTED: tuple[str, ...] = ("6450", "6505", "6523", "6534", "6566", "6591", "6623")


def _rebuilt(donor: str) -> bool:
    """Does an authoritative rebuild exist for this donor?

    Checked against the filesystem rather than a hand-maintained list, because a list of what is
    blocked has to be remembered and a check does not. The previous form was a literal tuple with
    the note "until the QuPath re-export"; the re-export landed and the tuple would have kept
    refusing the donors it named.
    """
    import subpanel_io as S
    return (S.REVISED / f"donor_id={donor}/data.parquet").exists()


def resolve(donors: list[str] | None = None, include_marginal: bool = False) -> list[str]:
    """The donor list to build over, refusing anything knowingly broken."""
    if donors:
        blocked = sorted(d for d in donors if d in BATCH2_REEXPORTED and not _rebuilt(d))
        if blocked:
            raise ValueError(
                f"no rebuild from the revised export for {blocked}. Their run partitions still "
                "carry the duplicated islet detections, so they cannot be built from the run.")
        return list(donors)
    return list(CLEAN) + (sorted(MARGINAL) if include_marginal else [])


#: What this cohort cannot answer, stated once so no figure has to discover it.
COHORT_LIMITS = """
CLEAN carried no Aab+ donor until 2026-09-15, because all three were either in batch 2 or carried
a small islet deficit, so no disease-progression comparison was available -- only ND against T1D
on two T1D donors.  The re-export resolves that: 6450 and 6505 are Aab+ and both now qualify.

What remains: every batch-2 donor's marker values are REUSED from the original REDSEA run, joined
on object_id.  That is sound -- the deduplicated cells carry the same measurements, and every
analysis-eligible revised cell has an aligned row, checked at 0 missing on all seven -- but it
means their intensities were not recomputed against the corrected segmentation.  And acquisition
batch is still aliased with disease status at the coarse level; correct on stain_group or
acq_run, never on batch.
"""
