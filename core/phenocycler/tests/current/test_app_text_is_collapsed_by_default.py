"""Large blocks of explanation are collapsed by default, behind a descriptive toggle.

The maintainer's rule, 2026-09-18: "if large blocks of explanatory text are needed then they should
be hidden by default with a descriptive icon to hide/unhide."

Measured when the rule was set: 4,665 words of standing explanation across 53 blocks in 12 files,
printed above the figures on every render. Statistics alone carried 2,329 -- one block of 1,061
words. None of it was wrong. What was wrong is that a caveat which is always true was given the
same weight, and the same vertical space, as a number that changes.

`ieCore::ie_explain()` is the pattern: a native `<details>`/`<summary>` disclosure, collapsed,
keyboard reachable, no server round-trip and no observer per card.

WHAT THIS TEST DOES NOT DO. It does not cap MEASURED state -- a count, a donor list, a withheld
reason, an availability note. Those are facts about what the reader is looking at and they stay
visible however long they get. It counts only prose literals, which is a static approximation: it
cannot see a string built at runtime, and it cannot tell a caveat from a caption. So the budgets
below are a RATCHET, not a specification. They record what each file carries today; lowering one is
always allowed, raising one needs a reason in the commit message.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(subprocess.run(["git", "rev-parse", "--show-toplevel"],
                           cwd=Path(__file__).resolve().parent,
                           capture_output=True, text=True, check=True).stdout.strip())
APP_R = REPO / "apps/islet/R"

#: Words of always-visible prose each file may carry. Recorded 2026-09-18 after the first pass.
#: A file absent from here may carry none. Lower these; do not raise them without saying why.
BUDGET: dict[str, int] = {
    "mod_statistics_server.R": 800,
    "mod_statistics_ui.R": 560,
    "mod_inflammation_server.R": 540,
    "card_s3_compartments.R": 430,
    "card_s2_trajectory.R": 430,
    "00_globals.R": 270,
    "mod_trajectory_ui.R": 230,
    "mod_inflammation_ui.R": 210,
    "mod_plot_server.R": 190,
    "card_s1_size_composition.R": 150,
    "card_g1_grid.R": 150,
    "mod_spatial_server.R": 140,
}

#: Style values, class names and ids are not prose, and a literal with a colon early is a CSS rule.
_NOT_PROSE = ("ie-", "font-", "margin", "padding", "background", "border", "color:", "#", "http",
              "width", "height", "display", "cursor", "list-style")

_PROSE_CALL = re.compile(
    r'(?:shiny::)?(?:tags\$p|tags\$small|helpText|p)\((.*?)\n\s*\)', re.S)


def _strip_explained(source: str) -> str:
    """Everything outside an `ie_explain(...)` call -- i.e. what renders before anything is opened."""
    out, i = [], 0
    while True:
        start = source.find("ie_explain(", i)
        if start == -1:
            out.append(source[i:])
            return "".join(out)
        out.append(source[i:start])
        j, depth = source.find("(", start), 0
        while j < len(source):
            if source[j] == "(":
                depth += 1
            elif source[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        i = j + 1


def _visible_words(source: str) -> int:
    total = 0
    for call in _PROSE_CALL.finditer(_strip_explained(source)):
        literals = [s for s in re.findall(r'"([^"]*)"', call.group(1))
                    if len(s) > 25 and not s.startswith(_NOT_PROSE)]
        words = " ".join(literals).split()
        if len(words) >= 12:          # a short caption is not a block of explanation
            total += len(words)
    return total


def test_ie_explain_exists_and_is_collapsed_by_default():
    helper = (REPO / "apps/R/ieCore/R/ie_explain.R").read_text()
    assert "tags$details" in helper, "the disclosure must be a native <details>, not a modal"
    assert "tags$summary" in helper
    # `open` is opt-in: the default argument must be FALSE, or nothing is collapsed.
    assert re.search(r"open\s*=\s*FALSE", helper), "ie_explain must default to collapsed"
    namespace = (REPO / "apps/R/ieCore/NAMESPACE").read_text()
    assert "export(ie_explain)" in namespace


@pytest.mark.parametrize("name", sorted(BUDGET))
def test_no_file_exceeds_its_visible_prose_budget(name):
    path = APP_R / name
    if not path.exists():
        pytest.skip(f"{name} no longer exists")
    words = _visible_words(path.read_text())
    assert words <= BUDGET[name], (
        f"{name} carries {words} words of always-visible explanation, over its budget of "
        f"{BUDGET[name]}. Wrap the standing text in ieCore::ie_explain() -- measured state "
        f"(counts, donor lists, withheld reasons) stays visible and is not what this counts.")


def test_files_outside_the_budget_carry_no_block_of_explanation():
    offenders = {
        p.name: _visible_words(p.read_text())
        for p in sorted(APP_R.glob("*.R")) if p.name not in BUDGET
    }
    over = {k: v for k, v in offenders.items() if v > 120}
    assert not over, (
        "these files grew a block of standing explanation and are not in BUDGET: "
        + ", ".join(f"{k} ({v} words)" for k, v in sorted(over.items(), key=lambda kv: -kv[1])))
