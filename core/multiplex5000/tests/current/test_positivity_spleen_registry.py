"""The spleen PC67 typing registry: acquired, registered, measured on spleen, and reviewed.

Maintainer 2026-11-01: "Start the registry over"; "Stop trying to apply parameters used from a different dataset."
The registry (written by SpleenFollicleCounterQP scripts/build_spleen_typing_registry.py) holds each marker's role,
reading and lineage-levels method; cells are named by the spleen hierarchy (positivity.hierarchy)."""

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


def test_every_typing_and_state_marker_has_a_lineage_levels_method(reg):
    from multiplex5000.positivity import lineage_levels as LL

    for m in reg.table.index:
        method = reg.background_method(m)
        if reg.role(m) in ("defining", "shared", "state"):
            assert method in LL.METHODS, m
            if method == "known_negative":
                assert set(RG.split(reg.table.loc[m, "level_reference"].replace("+", ";"))) <= set(LL.POPULATIONS), m
        else:
            assert reg.role(m) == "none" and method == "", m
    for lineage, marker in (("B", "CD20"), ("T", "CD3e")):
        assert [m for m in reg.broad(lineage) if reg.table.loc[m, "required"] == "yes"] == [marker]


def test_nuclear_antigens_are_read_on_the_nucleus_and_everything_else_on_the_whole_cell(reg):
    nuclear = {"Ki67", "PCNA", "Histone_H3_(p_Ser28)", "FOXP3", "TCF_1", "TOX", "TP63", "SOX2"}
    for m in reg.table.index:
        assert reg.measurement(m) == ("nucleus" if m in nuclear else "whole_cell"), m
        assert reg.differentiation(m) in ("", "membrane"), m


def test_nothing_is_carried_from_another_dataset(reg):
    consts = pd.read_csv(reg.directory / f"{reg.name}__constants.csv", keep_default_na=False, dtype=str)
    assert not consts["source"].str.contains("carried from").any()
    assert not reg.table["source"].str.contains("carried from").any()
    for name in ("ISOLATION_UM", "LEVEL_QUANTILE", "BACKGROUND_K", "KNOWN_NEGATIVE_MIN_ABOVE", "CANARY_MARGIN",
                 "CANARY_MIN_DEPTH", "CANARY_TIE_INTERIOR_LEAD"):
        reg.const(name)
    assert consts.set_index("name").loc["ISOLATION_UM", "source"].startswith("measured on the 8 PC67 sections")


def test_every_proposal_named_exists_and_all_are_reviewed(reg):
    table = RG.proposals(reg)
    constants = pd.read_csv(reg.directory / f"{reg.name}__constants.csv", keep_default_na=False, dtype=str)
    named = {p for col in (reg.table["proposal"], constants["proposal"]) for v in col for p in RG.split(v)}
    assert named <= set(table["id"])
    assert dict(zip(table["id"], table["status"])) == {
        "P1": "accepted", "P2": "accepted", "P3": "accepted", "P4": "accepted", "P5": "accepted",
        "P6": "rejected", "P7": "accepted", "P8": "rejected", "P9": "rejected", "P10": "accepted",
        "P11": "accepted", "P12": "accepted", "R1": "accepted"}
    assert table["decision"].str.startswith("maintainer 2026-11-01").all()
    RG.require_reviewed(reg)
    RG.require_reviewed(RG.load("codex-uniform-32", "pancreas"))   # no proposals file: nothing to wait for

