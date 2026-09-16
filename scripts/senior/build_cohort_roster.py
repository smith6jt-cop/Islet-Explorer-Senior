#!/usr/bin/env python3
"""Emit `apps/islet/config/cohort.json` — one cohort roster, read by R and Python alike.

The authority for *who is excluded* is `core/phenocycler/phenocycler/cohort.py`
(`DONOR_EXCLUSIONS`), pinned by `core/phenocycler/tests/current/test_cohort.py`. This builder
never retypes an exclusion or its reason; it imports the authority module and asks it.

The authority for *who is a candidate* is two git-tracked, non-generated sources:

  - `core/phenocycler/cohort-spec.json` `expected_donors` — the vetted panc-cohort roster
    (already excludes 6457/6579 by omission, since the spec was written after they were pulled).
  - `DONOR_EXCLUSIONS.keys()` from the authority itself — so an excluded donor still appears in
    the roster (as an excluded entry, with its measured reason) rather than silently vanishing.
    A card must be able to say "6457 was pulled, here is why", which requires the donor to be
    *in* the roster.

`donor_id` is a string and is not normalised across trees (0112 vs 112) -- see
docs/ISLET_LAYERS_SPEC.md ss 0.3. Every value here is normalised on read with `norm_donor`
before being written, and nothing downstream should have to normalise again.

Usage:
    python scripts/senior/build_cohort_roster.py [--out apps/islet/config/cohort.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PHENOCYCLER_ROOT = REPO_ROOT / "core" / "phenocycler"
COHORT_SPEC_PATH = PHENOCYCLER_ROOT / "cohort-spec.json"
DEFAULT_OUT = REPO_ROOT / "apps" / "islet" / "config" / "cohort.json"

sys.path.insert(0, str(PHENOCYCLER_ROOT))
from phenocycler.cohort import DONOR_EXCLUSIONS, filter_eligible_donors  # noqa: E402


def norm_donor(value: str) -> str:
    """'0112' -> '112'. Never int(): 6533 must not become 6533.0 through a float column.

    The one shared normaliser docs/ISLET_LAYERS_SPEC.md ss 0.3 requires every builder to use.
    """
    return str(value).strip().lstrip("0") or "0"


def load_expected_donors(spec_path: Path) -> list[str]:
    spec = json.loads(spec_path.read_text())
    return [norm_donor(d) for d in spec.get("expected_donors", [])]


def build_roster(spec_path: Path = COHORT_SPEC_PATH) -> dict:
    candidate_donors = load_expected_donors(spec_path)
    excluded_ids = sorted({norm_donor(d) for d in DONOR_EXCLUSIONS})

    # The full roster is the vetted candidates plus the known-excluded donors: an excluded donor
    # must be *in* the roster (as an excluded entry) so a card can name it and say why, rather
    # than a donor that silently never appears anywhere.
    all_donors = sorted(set(candidate_donors) | set(excluded_ids))

    eligible_donors = filter_eligible_donors(all_donors)

    excluded_donors = [
        {"donor_id": donor_id, "reason": DONOR_EXCLUSIONS[donor_id]}
        for donor_id in excluded_ids
    ]

    return {
        "spec_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": {
            "candidate_donors": str(spec_path.relative_to(REPO_ROOT)),
            "exclusion_authority": str(
                (PHENOCYCLER_ROOT / "phenocycler" / "cohort.py").relative_to(REPO_ROOT)
            ),
        },
        "all_donors": all_donors,
        "excluded_donors": excluded_donors,
        "eligible_donors": eligible_donors,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", type=Path, default=COHORT_SPEC_PATH,
                    help="cohort-spec.json to read expected_donors from")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="where to write cohort.json")
    args = ap.parse_args(argv)

    roster = build_roster(args.spec)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(roster, indent=2) + "\n")

    print(
        f"[done] {args.out}: {len(roster['all_donors'])} donors, "
        f"{len(roster['excluded_donors'])} excluded, "
        f"{len(roster['eligible_donors'])} eligible"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
