"""The spleen PC67 typing registry: acquired, registered, internally consistent, and reviewed.

Since the maintainer's review of 2026-11-01 the registry holds only what positivity needs; cells are named by the
spleen hierarchy (positivity.hierarchy), so there is no spleen labels file."""

from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from multiplex5000 import roots
from multiplex5000.positivity import levels as LV
from multiplex5000.positivity import registry as RG
from multiplex5000.positivity import tissue as TS

SPLEEN = ("phenocycler-pc67", "spleen")
HEADER = roots.REPO / "core/panels/data/tissues/spleen/pc67_export_header.tsv"
MARKER_REGISTRY = roots.REPO / "core/multiplex5000/multiplex5000/marker_registry.json"
NOT_ON_PC67 = {"B3TUBB", "CD99", "GCG", "IAPP", "INS", "Iba1", "Ker8_18", "Keratin_5", "SST", "b_Catenin1"}
ROLES = {"defining", "shared", "state", "none", "never_used"}


def _ies_name(channel: str) -> str:
    """IES marker name of a QuPath channel (redsea.SANITIZE, and the manifest's PDL1 -> PD_L1 alias)."""
    name = re.sub(r"[ /\-]", "_", channel)
    return {"PDL1": "PD_L1"}.get(name, name)


@pytest.fixture(scope="module")
def reg():
    return RG.load(*SPLEEN)


@pytest.fixture(scope="module")
def channels():
    names = HEADER.read_text().splitlines()[0].split("\t")
    return [m.group(1) for c in names if (m := re.match(r"^Cell: (.+): Mean$", c))]


def test_every_pc67_channel_has_one_row_and_nothing_else(reg, channels):
    assert len(channels) == 67
    assert set(reg.table["channel"]) == set(channels)
    assert set(reg.table.index) == {_ies_name(c) for c in channels}
    assert (reg.table["on_every_image"] == "yes").all()
    assert not set(reg.table.index) & NOT_ON_PC67


def test_every_marker_is_an_ies_registry_marker(reg):
    known = {m["name"] for m in json.loads(MARKER_REGISTRY.read_text())["markers"]}
    assert set(reg.table.index) <= known | {"DAPI"}


def test_roles_are_valid_and_compartments_belong_to_defining_markers(reg):
    assert set(reg.table["role"]) <= ROLES
    for m in reg.markers("defining"):
        assert reg.table.loc[m, "compartment"] and reg.lineages(m), m
    assert set(reg.compartments()) == {"Hematopoietic", "Endothelial", "Mesenchymal", "Epithelial"}


#: the columns that name markers (free-text provenance such as "the sheet's Iba1 row" is not a reference)
MARKER_COLUMNS = ("channel", "requires", "coexpressed_with")


def test_no_unacquired_marker_is_referenced(reg):
    named = {*reg.table.index}
    for col in MARKER_COLUMNS:
        named |= {m for v in reg.table[col] for m in RG.split(v)}
    assert not named & NOT_ON_PC67
    assert named - {"single_lineage"} <= set(reg.table.index) | {c for c in reg.table["channel"]}


def test_names_come_from_the_spleen_hierarchy_not_a_labels_file(reg):
    """Every sheet marker PC67 acquires has a registry row, so the hierarchy reads a call for it."""
    from multiplex5000.positivity import hierarchy as H
    from multiplex5000.tissues import spleen

    assert not (reg.directory / f"{reg.name}__labels.csv").exists()
    source = spleen.load_rule_source(roots.REPO / "core/panels/data/tissues/spleen/manifest.yml")
    names = H.sheet_names(source, reg.table.index)
    sheet = {m for r in source.rules for m in r.markers}
    assert sheet - set(names.values()) == {"B3TUBB", "Keratin 5"}
    assert all(reg.role(m) != "never_used" for m in names)


def test_check_bars_are_constants_and_the_method_reads_its_constants(reg):
    for bar in filter(None, reg.table["check_bar"]):
        assert reg.const(bar) == reg.const(bar)        # present and a number
    for name in ("HORMONE_SCALE_TOP_QUANTILE", "CD45_AMONG_HORMONE_MAX", "EMPTY_NULL_QUANTILE", "ISOLATION_UM"):
        reg.const(name)                                  # read unconditionally by levels/identity


def test_the_spleen_hook_accepts_it_and_reads_epithelial_and_stromal_references(reg):
    """P7 (maintainer 2026-11-01): the leukocyte checks' known-negative cells are epithelial and stromal cells; a
    marker with no chosen level is read at its own background peak, a chosen one at its level."""
    import numpy as np

    hook = TS.hooks("spleen")
    TS.check_registry(reg, hook)
    roles = LV.roles(reg)
    assert roles["core"] == ["CD3e", "CD45", "CD68", "CD163"] and not roles["hormones"]
    assert set(reg.markers("state")) >= {"Ki67", "PD_1", "PD_L1", "HLA_DR"}
    assert hook.CHECK_POPULATION == "reference"
    read = []
    background = {"Pan_Cytokeratin": np.array([1, 0, 0, 0], bool), "E_cadherin": np.zeros(4, bool),
                  "EpCAM": np.zeros(4, bool)}
    cells = hook.check_population(np.zeros(4), reg, roles, {"SMA": np.array([0, 1, 0, 0], bool)},
                                  lambda m: read.append(m) or background[m])
    assert cells.tolist() == [True, True, False, False]
    assert sorted(read) == ["E_cadherin", "EpCAM", "Pan_Cytokeratin"]    # SMA had its level; no CD45, no Vimentin


def test_every_proposal_named_exists_and_all_are_reviewed(reg):
    table = RG.proposals(reg)
    constants = pd.read_csv(reg.directory / f"{reg.name}__constants.csv", keep_default_na=False, dtype=str)
    named = {p for col in (reg.table["proposal"], constants["proposal"]) for v in col for p in RG.split(v)}
    assert named <= set(table["id"])
    assert dict(zip(table["id"], table["status"])) == {
        "P1": "accepted", "P2": "accepted", "P3": "accepted", "P4": "accepted", "P5": "accepted",
        "P6": "rejected", "P7": "accepted", "P8": "rejected", "P9": "rejected", "P10": "accepted",
        "P11": "accepted"}
    assert table["decision"].str.startswith("maintainer 2026-11-01").all()
    RG.require_reviewed(reg)
    RG.require_reviewed(RG.load("codex-uniform-32", "pancreas"))   # no proposals file: nothing to wait for


def test_cd45_and_cd163_take_the_majority_marker_levels(reg):
    """P10/P11 (maintainer 2026-11-01): nearly every spleen cell carries CD45, so its densest peak is not background."""
    assert LV.methods_for(reg, "CD45") == ["edge", "gate"] == LV.methods_for(reg, "CD163")
    assert LV.check_isolation(reg, "CD45") == "none" and LV.check_isolation(reg, "CD3e") == ""
    assert LV.check_isolation(RG.load("codex-uniform-32", "pancreas"), "CD45") == ""     # no column: the source's rule
