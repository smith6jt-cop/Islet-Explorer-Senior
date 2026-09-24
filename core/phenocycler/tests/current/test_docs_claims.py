"""What the docs claim about the tree, checked against the tree.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_docs_claims.py -q

On 2026-09-21 thirty-one sentences in this repository's docs were contradicted by the tree -- led
by "no donor has been re-run with the islet chain", three days after all twenty had been.  A
sentence nobody re-measures goes stale silently.  These are the claims a new session reads first,
held to what they describe: every script a CLAUDE.md names, the exclusions, the cohort the spec
declares, and (in the data lane) the typing numbers the root CLAUDE.md quotes.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

from phenocycler.cohort import DONOR_EXCLUSIONS, REINSTATED

REPO = Path(__file__).resolve().parents[4]

#: Names a CLAUDE.md cites that are not in this repository, each with the doc's own account of
#: why.  A name here must still be absent and still be named by that doc.
ABSENT_BY_DESIGN = {
    "analysis/xenium/CLAUDE.md": {
        # "moved to hpc/" and "are gone", stated in the doc's own out-of-date note
        "scripts/run_pipeline.sh", "scripts/env_rapids.sh", "scripts/R/xenium_analysis.R",
        # "have never been in this repository", the same note
        "select_embedding_features.py", "scripts/select_embedding_features.py",
        "lineage_phenotype.py", "scripts/lineage_phenotype.py",
        "refine_rare_celltypes.py", "scripts/refine_rare_celltypes.py",
        "resolve_indeterminate.py", "scripts/resolve_indeterminate.py",
        "fix_acinar_argmax.py", "validate_hierarchical_pilot.py",
        "scripts/validate_hierarchical_pilot.py", "program_panels.py", "scripts/program_panels.py",
    },
    "core/panels/CLAUDE.md": {"test-panel_validate.R"},  # "planned, never written"
}


def _tracked(pattern: str = "") -> list[str]:
    args = ["git", "ls-files"] + ([pattern] if pattern else [])
    return subprocess.run(args, cwd=REPO, capture_output=True, text=True, check=True).stdout.split()


def _named_scripts(text: str) -> set[str]:
    """Backticked script paths; a glob or a placeholder (`*`, `<`, `{`) names no one file."""
    return {t for t in re.findall(r"`([^`\s]+\.(?:py|R|sh))`", text)
            if not any(c in t for c in "*<{")}


def test_every_script_a_claude_md_names_exists():
    tracked = set(_tracked())
    basenames = {Path(t).name for t in tracked}
    missing = {}
    for doc in _tracked("*CLAUDE.md"):
        here = Path(doc).parent
        absent = ABSENT_BY_DESIGN.get(doc, set())
        for name in sorted(_named_scripts((REPO / doc).read_text())):
            found = (name in tracked or str(here / name) in tracked
                     or ("/" not in name and name in basenames)
                     or ("/" in name and any(t.endswith("/" + name) for t in tracked)))
            if not found and name not in absent:
                missing.setdefault(doc, []).append(name)
            if found and name in absent:
                pytest.fail(f"{doc}: {name} is listed as absent and exists; take it off the list")
    assert not missing, f"CLAUDE.md files name scripts that do not exist: {missing}"


def test_every_absence_the_list_excuses_is_still_named_by_its_doc():
    for doc, names in ABSENT_BY_DESIGN.items():
        text = (REPO / doc).read_text()
        stale = sorted(n for n in names if f"`{n}`" not in text)
        assert not stale, f"{doc} no longer names {stale}; take them off ABSENT_BY_DESIGN"


def test_the_exclusions_the_docs_state_are_cohort_pys():
    """6457 is the one central exclusion; 6579's was lifted on 2026-09-17."""
    assert set(DONOR_EXCLUSIONS) == {"6457"}
    assert "6579" in REINSTATED
    # The forms the stale sentences took, backticked or bold or neither: "6457 and 6579 remain
    # centrally excluded", "... are centrally excluded", "... 6579, which cohort.py centrally
    # excludes", and "6579 is/remains excluded".
    stale = re.compile(r"6457`? and `?6579`?(\*\*)?,? (remains?|are|which)\b.{0,70}?exclud"
                       r"|6579`?(\*\*)? (remains|is) (centrally )?excluded", re.IGNORECASE)
    undated = [d for d in _tracked("*.md")
               if not re.search(r"(^|/)(PLAN_)?[A-Z0-9_]*20\d\d-\d\d-\d\d", d)]
    wrong = []
    for doc in undated:
        text = " ".join((REPO / doc).read_text().split())
        wrong += [f"{doc}: ...{m.group(0)}..." for m in stale.finditer(text)]
    assert not wrong, "a live doc still calls 6579 centrally excluded:\n  " + "\n  ".join(wrong)


def test_the_cohort_spec_declares_the_acquisition_cohort():
    """The same 23 the live manifest declares: 22 sequential-slide pancreas donors less the one
    central exclusion, plus the two same-slide PhenoXeno sections UniFORM needs."""
    spec = json.loads((REPO / "core/phenocycler/cohort-spec.json").read_text())
    expected = [str(d) for d in spec["expected_donors"]]
    assert len(expected) == len(set(expected)) == 23
    assert set(expected) == {str(row["donor_id"]) for row in spec["images"]}
    assert {"0112", "6533", "6579"} <= set(expected)
    assert not set(expected) & set(DONOR_EXCLUSIONS)


def _stated(pattern: str, text: str) -> tuple[int, ...]:
    match = re.search(pattern, text)
    assert match, f"the root CLAUDE.md no longer states /{pattern}/"
    return tuple(int(g.replace(",", "")) for g in match.groups())


@pytest.mark.data
def test_the_typing_numbers_the_root_claude_md_states_are_on_disk():
    cells = REPO / "data/development/typing/cells"
    if not cells.is_dir():
        pytest.skip("typing/cells is not mounted")
    duckdb = pytest.importorskip("duckdb")
    text = " ".join((REPO / "CLAUDE.md").read_text().split())
    rows, = _stated(r"\*\*([\d,]+) cells\*\*", text)
    triple = _stated(r"the endocrine triple is Beta ([\d,]+) / Alpha ([\d,]+) / Delta ([\d,]+)", text)
    routes = _stated(r"cd99 ([\d,]+) · elongation ([\d,]+) · hormone ([\d,]+)", text)
    source = (f"read_parquet('{cells}/donor_id=*/data.parquet', hive_partitioning=true)")
    con = duckdb.connect()
    n, donors = con.execute(f"select count(*), count(distinct donor_id) from {source}").fetchone()
    assert (n, donors) == (rows, 20)
    measured = dict(con.execute(f"select cell_class, count(*) from {source} "
                                f"where cell_class in ('Beta', 'Alpha', 'Delta') group by 1").fetchall())
    assert (measured["Beta"], measured["Alpha"], measured["Delta"]) == triple
    by_route = dict(con.execute(f"select islet_route, count(*) from {source} where islet_route "
                                f"in ('cd99', 'elongation', 'hormone') group by 1").fetchall())
    assert (by_route["cd99"], by_route["elongation"], by_route["hormone"]) == routes
    routed = con.execute(f"select count(distinct donor_id) from {source} "
                         f"where islet_route in ('cd99', 'elongation', 'hormone')").fetchone()[0]
    assert routed == donors, "the doc says the chain is applied to every donor"
