"""The negative-evidence redirect must be counted on the partition whose labels ship.

The defect this pins: ``panc_xenium_cluster_name.py`` computed its redirect count from ``main``,
the PARENT partition, while every label in the release comes from ``fdec``, the final subcluster
partition, and ``panc_xenium_release.py`` published that parent number under a name that did not
say so.  The parent count is 0 in every section of the pancreas cohort, so ``6482_c18_s07`` --
215 cells, Immune_NK_ILC, named by the redirect -- was reported as 0 redirects for every donor.

Two events are separated here and they are not interchangeable:

* APPLIED -- the argmax was Immune_T, neither CD3E nor CD3G stood out, NK genes did, so the
  winner moved.
* SHIPPED -- that moved winner then cleared the enrichment gate, so the cluster went out as
  Immune_NK_ILC carrying ``REDIRECT_REASON``.  A moved winner whose NK score is not above the
  partition mean is uncalled instead, and counting only the move reports a label that never
  existed.

Nothing here asserts a biological proportion or reproduces a previous finding.  The synthetic
tests are statements about the counting rule.  The artifact tests are internal coherence: a
number a report publishes about a frame must equal what reading that frame says, which is a
statement about the build.

Run from the repository root with ``.venv/bin/python -m pytest`` on this file.
"""

from __future__ import annotations

import glob
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.tissues import pancreas_rna_naming

REPO = Path(__file__).resolve().parents[4]
WORKER_PATH = REPO / "analysis/xenium/scripts/panc_xenium_cluster_name.py"
IMPLEMENTATION_PATH = (REPO / "core/multiplex5000/multiplex5000/tissues"
                       / "pancreas_xenium_worker.py")
RELEASE_PATH = REPO / "analysis/xenium/scripts/panc_xenium_release.py"


def _load(name: str, path: Path):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


worker = _load("panc_worker_under_test", WORKER_PATH)
release = _load("panc_release_under_test", RELEASE_PATH)


def test_release_worker_uses_package_naming_decisions():
    for name in ("naming_policy_panels", "exact_detection_q", "decide", "redirect_counts"):
        assert getattr(worker, name) is getattr(pancreas_rna_naming, name)

#: The work directory of the pancreas run.  Overridable so the artifact tests can be pointed at
#: another run rather than restating a path in a second place.
WORK = Path(os.environ.get("PANC_XENIUM_WORK", "/data/xenium_typing/panc_2026-09-17/_work"))


def _sections() -> list[tuple[str, Path, Path]]:
    out = []
    for f in sorted(glob.glob(str(WORK / "*.subclusters.parquet"))):
        sub = Path(f)
        calls = sub.with_name(sub.name.replace(".subclusters.", ".calls."))
        if calls.is_file():
            out.append((sub.name.split(".")[0], calls, sub))
    return out


@pytest.fixture(scope="module")
def sections() -> list[tuple[str, Path, Path]]:
    """The per-section artifacts, read only when a test asks for them.

    Collecting this module must never touch ``WORK``: the glob runs inside the fixture, so a
    checkout without the run collects the same two ids and skips them with the reason below.
    """
    found = _sections()
    if not found:
        pytest.skip(f"no per-section artifacts under {WORK}")
    return found


# --------------------------------------------------------------- the counting rule itself

def _dec(reasons, redirected):
    return {"reasons": list(reasons), "redirected": list(redirected)}


def test_counts_the_partition_it_is_given_and_separates_applied_from_shipped():
    """The exact shape of the defect: a parent with no redirect, a child that shipped one."""
    parent = _dec(["", "", ""], [False, False, False])
    child = _dec(
        [worker.REDIRECT_REASON,
         "winning_core_Immune_NK_ILC_score_-0.069_not_above_partition_mean", "", ""],
        [True, True, False, False])
    sizes_child = np.array([215.0, 211.0, 1000.0, 5.0])

    got_child = worker.redirect_counts(child, sizes_child)
    assert got_child == {"clusters_applied": 2, "clusters_shipped": 1, "cells_shipped": 215}

    got_parent = worker.redirect_counts(parent, np.array([1200.0, 900.0, 31.0]))
    assert got_parent == {"clusters_applied": 0, "clusters_shipped": 0, "cells_shipped": 0}

    # the two partitions disagree, which is why one number cannot stand for both
    assert got_parent["clusters_applied"] != got_child["clusters_shipped"]


def test_cells_shipped_is_zero_when_the_move_was_never_admitted():
    """A redirect that was applied and then failed the gate contributes no cells and no label."""
    applied_only = _dec(
        ["winning_core_Immune_NK_ILC_score_-0.03_not_above_partition_mean"], [True])
    assert worker.redirect_counts(applied_only, np.array([763.0])) == {
        "clusters_applied": 1, "clusters_shipped": 0, "cells_shipped": 0}


# ------------------------------------------------------- decide() writes the reason it counts

def _decide_on(genes, cluster_means, cluster_dets, d_other):
    """Run ``decide`` over a hand-built partition of three clusters and two cores."""
    idx = {g: i for i, g in enumerate(genes)}
    means = np.array(cluster_means, dtype=float)
    dets = np.array(cluster_dets, dtype=float)
    mu, sd = means.mean(axis=0), means.std(axis=0, ddof=0)
    z = np.nan_to_num((means - mu) / np.where(sd > 0, sd, np.nan), nan=0.0)
    core_names = ["Immune_T", "Immune_NK_ILC"]
    cols = {"Immune_T": np.array([idx["CD3E"], idx["CD3G"]]),
            "Immune_NK_ILC": np.array([idx["NKG7"], idx["KLRD1"]])}
    t_cols = np.array([idx[g] for g in worker.T_DEFINING if g in idx])
    nk_cols = np.array([idx[g] for g in worker.NK_NAMING_HINT if g in idx])
    return worker.decide(z, dets, np.array(d_other, dtype=float), worker.Z_STANDOUT,
                         core_names, cols, cols, t_cols, nk_cols)


#: One partition, three clusters, four genes: CD3E, CD3G, NKG7, KLRD1.  In cluster 0 the T core
#: still WINS the argmax -- its means are the partition's high side -- but neither CD3 gene
#: stands out, because its detection fraction there is below the rest of the partition's.  The
#: NK genes do stand out and clear the cluster's own background.  That is the redirect.
GENES = ["CD3E", "CD3G", "NKG7", "KLRD1"]
MEANS = [[4.0, 4.0, 2.0, 2.0], [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0]]
DETS = [[0.1, 0.1, 0.6, 0.5], [0.0, 0.0, 0.2, 0.2], [0.0, 0.0, 0.0, 0.0]]
D_OTHER = [[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]


def test_decide_writes_the_shared_reason_constant_when_the_redirect_is_admitted():
    dec = _decide_on(GENES, MEANS, DETS, D_OTHER)
    # the T core won the argmax and was overruled by the negative evidence
    assert dec["win"][0] == 0
    assert dec["calls"][0] == "Immune_NK_ILC"
    assert bool(dec["redirected"][0]) is True
    assert dec["reasons"][0] == worker.REDIRECT_REASON
    # and the counter sees exactly that one cell block
    assert worker.redirect_counts(dec, np.array([215.0, 1.0, 1.0])) == {
        "clusters_applied": 1, "clusters_shipped": 1, "cells_shipped": 215}


def test_an_applied_redirect_that_fails_the_enrichment_gate_ships_no_label():
    """Same partition, but only ONE of the two NK genes is on cluster 0's high side.

    NKG7 still stands out, so the negative evidence is unchanged and the winner still moves;
    the NK core's score -- the mean over both its naming genes -- is then not above the
    partition mean, and the gate refuses the move.  The cluster goes out uncalled, counting as
    applied and not as labelled.  This is the shape of 6347_c20_s03, 6380_c11_s11 and
    6523_c19_s04 in the run on disk.
    """
    means = [[4.0, 4.0, 3.0, 0.0], [0.0, 0.0, 0.0, 3.0], [0.0, 0.0, 1.0, 3.0]]
    dec = _decide_on(GENES, means, DETS, D_OTHER)
    assert bool(dec["redirected"][0]) is True
    assert dec["calls"][0].startswith("Uncalled_")
    assert dec["reasons"][0] != worker.REDIRECT_REASON
    assert worker.redirect_counts(dec, np.array([763.0, 1.0, 1.0])) == {
        "clusters_applied": 1, "clusters_shipped": 0, "cells_shipped": 0}


# ------------------------------------------------------------------- the fields that ship

def test_the_worker_report_names_each_partition_and_drops_the_ambiguous_name():
    src = IMPLEMENTATION_PATH.read_text()
    assert '"n_clusters_cd3_negative_nk_redirect"' not in src, (
        "the report key that did not say which partition it counted is back")
    assert "redirect_counts(fdec, fsizes)" in src, (
        "the shipped count must come from the final subcluster partition")
    assert "redirect_counts(main, sizes)" in src
    for key in ('"n_subclusters_labelled_by_cd3_negative_nk_redirect"',
                '"n_cells_labelled_by_cd3_negative_nk_redirect"',
                '"n_parent_clusters_cd3_negative_nk_redirect_applied"'):
        assert key in src, key


def test_the_release_does_not_publish_the_parent_count_under_the_old_name():
    src = RELEASE_PATH.read_text()
    assert '"n_cd3_negative_nk_redirects"' not in src, (
        "the published field that carried the parent count is back")
    assert 'r["n_clusters_cd3_negative_nk_redirect"]' not in src


# ------------------------------------------------- the same rule, on the artifacts on disk

@pytest.mark.data
def test_release_counter_equals_the_cells_carrying_the_redirect_reason(sections):
    """What the release publishes must equal what reading the shipped frames says.

    This is the assertion that would have caught it: with the parent count in that field the
    published number was 0 while the calls frame carried 215 cells whose reason is the redirect.
    """
    shipped_anywhere = 0
    for section, calls_path, sub_path in sections:
        calls = pd.read_parquet(calls_path, columns=["reason", "lineage", "subcluster_id"])
        subs = pd.read_parquet(
            sub_path, columns=["subcluster_id", "reason", "cd3_negative_nk_redirect", "n_cells",
                               "lineage"])
        got = release.redirect_summary(calls, subs)

        by_cell = int(calls["reason"].eq(worker.REDIRECT_REASON).sum())
        shipped = subs["reason"].eq(worker.REDIRECT_REASON)
        by_subcluster = int(subs.loc[shipped, "n_cells"].sum())

        assert got["n_cells_labelled_by_cd3_negative_nk_redirect"] == by_cell, section
        assert by_cell == by_subcluster, (
            f"{section}: the cells carrying the redirect reason and the sizes of the "
            f"subclusters carrying it disagree")
        assert got["n_subclusters_labelled_by_cd3_negative_nk_redirect"] == int(shipped.sum())
        # every subcluster that shipped the reason went out as NK; nothing else may carry it
        assert set(subs.loc[shipped, "lineage"]) <= {"Immune_NK_ILC"}, section
        assert set(calls.loc[calls["reason"].eq(worker.REDIRECT_REASON), "lineage"]) \
            <= {"Immune_NK_ILC"}, section
        # a label can only ship where the move was applied
        assert got["n_subclusters_labelled_by_cd3_negative_nk_redirect"] \
            <= got["n_subclusters_cd3_negative_nk_redirect_applied"], section
        shipped_anywhere += by_cell

    assert shipped_anywhere > 0, (
        "no section in this work directory shipped a redirected label, so the check above "
        "passed without exercising anything")


@pytest.mark.data
def test_the_parent_partition_is_not_a_substitute_for_the_subcluster_partition(sections):
    """On these artifacts the two partitions give different answers, which is the whole point."""
    parent_applied = subcluster_shipped = 0
    for section, calls_path, sub_path in sections:
        cl_path = sub_path.with_name(sub_path.name.replace(".subclusters.", ".clusters."))
        if cl_path.is_file():
            parents = pd.read_parquet(cl_path, columns=["cd3_negative_nk_redirect"])
            parent_applied += int(parents["cd3_negative_nk_redirect"].fillna(False).sum())
        calls = pd.read_parquet(calls_path, columns=["reason"])
        subcluster_shipped += int(calls["reason"].eq(worker.REDIRECT_REASON).sum())
    assert subcluster_shipped > 0
    assert parent_applied != subcluster_shipped
