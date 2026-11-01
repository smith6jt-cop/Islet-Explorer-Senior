"""The spleen PC67 typing registry: acquired, registered, internally consistent, and gated on its proposals."""

from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from multiplex5000 import roots
from multiplex5000.positivity import identity as ID
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
LABEL_MARKER_COLUMNS = ("named_by", "profile")


def test_no_unacquired_marker_is_referenced(reg):
    labels = ID.labels_table(reg)
    named = {*reg.table.index}
    for col in MARKER_COLUMNS:
        named |= {m for v in reg.table[col] for m in RG.split(v)}
    for col in LABEL_MARKER_COLUMNS:
        named |= {m for v in labels[col] for m in RG.split(v.replace("|", ";"))}
    assert not named & NOT_ON_PC67
    assert named - {"single_lineage", "broad", "broad_or_subtype"} <= set(reg.table.index) | {
        c for c in reg.table["channel"]}


def test_labels_name_only_registry_lineages_and_markers(reg):
    labels = ID.labels_table(reg)
    lineages = {lin for m in reg.table.index if reg.role(m) == "defining" for lin, _ in reg.lineage_entries(m)}
    entries = {f"{lin}:{sub}" if sub else lin for m in reg.markers("defining") for lin, sub in reg.lineage_entries(m)}
    for r in labels.itertuples():
        if r.kind == "structure":
            assert r.named_by in reg.table.index and r.compartment in reg.compartments()
            continue
        assert r.lineage in entries or r.lineage.split(":")[0] in lineages, r.lineage
        if r.named_by not in ("", "broad", "broad_or_subtype"):
            assert r.named_by in reg.table.index
        for m in RG.split(r.profile.replace("|", ";")):
            assert m in reg.table.index, m


def test_check_bars_are_constants_and_the_method_reads_its_constants(reg):
    for bar in filter(None, reg.table["check_bar"]):
        assert reg.const(bar) == reg.const(bar)        # present and a number
    for name in ("HORMONE_SCALE_TOP_QUANTILE", "CD45_AMONG_HORMONE_MAX", "EMPTY_NULL_QUANTILE", "ISOLATION_UM"):
        reg.const(name)                                  # read unconditionally by levels/identity


def test_the_spleen_hook_accepts_it_and_the_method_reads_its_plan(reg):
    TS.check_registry(reg, TS.hooks("spleen"))
    roles = LV.roles(reg)
    assert roles["core"] == ["CD3e", "CD45", "CD68", "CD163"] and not roles["hormones"]
    plan = ID.build_plan(reg, reg.typing_markers())
    assert plan["comp"]["Endothelial"]["structure_only"]          # CD31 places no cell (P5)
    assert "CD31" in plan["comp"]["Hematopoietic"]["carried"]       # but adds depth (shared evidence)
    assert set(reg.markers("state")) >= {"Ki67", "PD_1", "PD_L1", "HLA_DR"}


def test_every_proposal_named_exists_and_typing_waits_for_their_review(reg):
    table = RG.proposals(reg)
    labels = ID.labels_table(reg)
    constants = pd.read_csv(reg.directory / f"{reg.name}__constants.csv", keep_default_na=False, dtype=str)
    named = {p for col in (reg.table["proposal"], labels["proposal"], constants["proposal"]) for v in col
             for p in RG.split(v)}
    assert named <= set(table["id"])
    assert set(table["status"]) == {"proposed"}
    with pytest.raises(RG.RegistryError, match="awaiting review"):
        RG.require_reviewed(reg)
    RG.require_reviewed(RG.load("codex-uniform-32", "pancreas"))   # no proposals file: nothing to wait for
