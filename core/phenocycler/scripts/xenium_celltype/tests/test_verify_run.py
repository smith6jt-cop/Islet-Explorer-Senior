"""Check 8 has to survive the vocabulary the composition writer actually uses.

`composition()` never writes a bare ``Unavailable``: a refusal keeps its reason, so the row
is ``Unavailable:too_few_cells`` or ``Unavailable:no_core_beats_its_matched_null`` (see
`test_nullgate.py::test_composition_never_renders_an_absence_as_a_zero`). Check 8 compared
those rows for equality with ``"Unavailable"``, which is never true -- so every refusal was
read as a lineage the composition *claims*, `claimed - awarded` was never empty, and the
check failed on a run that was in fact complete. `50_assemble.sh` exits on that status, so a
finished cohort reported as "run NOT complete".

The verifier is the thing that decides whether 20 donors are a result, so this builds a small
but complete run tree and puts it through `verify()` rather than calling the comparison
directly.
"""
import json
import pathlib
import sys

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hpg"))

import verify_run as V                                 # noqa: E402

DONOR, ROI, ARM, K = "6539", "panc", "nonspatial", 20
N_CELLS, M, N_SHARDS = 4, 4, 2


def _complete_run(tmp_path: pathlib.Path, comp: pd.DataFrame, calls: pd.DataFrame):
    """Everything checks 1-7 and 9 need, so a failure here can only be check 8."""
    run, repo = tmp_path / "run", tmp_path / "repo"
    (repo / "data").mkdir(parents=True)
    (repo / "data" / "qupath_manifest.json").write_text(
        json.dumps({"expected_donors": [DONOR]}))
    (run).mkdir(parents=True)
    (run / "work_manifest.json").write_text(json.dumps({
        "donors": [{"donor_id": DONOR, "roi": ROI}],
        "kodama": {"M": M, "n_shards": N_SHARDS, "seed": 0},
        "arms": [ARM],
    }))

    dest = run / f"donor_id={DONOR}" / f"roi={ROI}"
    (dest / "input").mkdir(parents=True)
    (dest / "input" / "prep.json").write_text(
        json.dumps({"complete": True, "n_cells": N_CELLS}))
    ids = [f"cell-{i}-1" for i in range(N_CELLS)]
    pd.DataFrame({"cell_id": ids}).to_csv(dest / "input" / "cell_id.csv", index=False)
    for f in ("pcs.mtx", "coords.csv"):
        (dest / "input" / f).touch()

    # two shards tiling 1..M against one input digest
    for i, (first, last) in enumerate([(1, M // 2), (M // 2 + 1, M)], start=1):
        (dest / f"res_{ARM}_shard{i}.json").write_text(json.dumps(
            {"complete": True, "k_global_first": first, "k_global_last": last,
             "input_digest": "deadbeef", "kodama_version": "x"}))
    pd.DataFrame({"cell_id": ids,
                  "dim1": [0.0, 1.0, 2.0, 3.0],
                  "dim2": [3.0, 2.0, 1.0, 0.0]}).to_csv(dest / f"embedding_{ARM}.csv",
                                                        index=False)
    (dest / f"kodama_{ARM}.json").write_text(json.dumps(
        {"complete": True, "n_shards": N_SHARDS}))
    calls.to_csv(dest / f"group_calls_{ARM}.csv", index=False)
    comp.to_csv(dest / f"composition_{ARM}.csv", index=False)
    return run, repo


def _frames():
    """One awarded lineage and one refusal, written the way `composition()` writes them."""
    calls = pd.DataFrame({
        "donor_id": DONOR, "arm": ARM, "k": K,
        "group": [0, 1],
        "lineage": ["Exocrine", ""],          # group 1 cleared nothing
        "reason": ["", "too_few_cells"],
        "runners_up": ["", ""],
    })
    comp = pd.DataFrame({
        "donor_id": DONOR, "arm": ARM, "k": K,
        "lineage": ["Exocrine", "Unavailable:too_few_cells"],
        "status": ["awarded", "unavailable"],
        "n_cells": [800, 200],
        "frac": [0.8, 0.2],
    })
    return comp, calls


def test_a_reasoned_unavailable_row_is_not_read_as_a_claimed_lineage(tmp_path):
    run, repo = _complete_run(tmp_path, *_frames())
    r = V.Report()
    V.verify(run, repo, r)

    eight = [p for p in r.problems if p.startswith("(8)")]
    assert not eight, eight
    assert r.ok(), r.problems


def test_the_unavailable_share_is_still_measured_per_k(tmp_path):
    """The same comparison feeds the reported share, which was silently 0 for every run."""
    run, repo = _complete_run(tmp_path, *_frames())
    r = V.Report()
    summary = V.verify(run, repo, r)

    assert summary[DONOR][ARM]["unavailable_by_k"] == {K: 0.2}
    assert any("Unavailable" in n for n in r.notes)
