"""Cell keys: acquisition + segmentation + native id, encoded so no two cells can share one.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_identity.py -q

The cases are the ones the 2026-10-08 identity contract names: native ids that collide across
acquisitions, separators inside ids, missing parts, and uniqueness checked rather than assumed.
"""

from __future__ import annotations

import pandas as pd
import pytest

from multiplex5000 import identity as I

AWKWARD = ["plain-id_1.2~", "has:colon", "has%percent", "slash/and space", "hash #1",
           "Ünïcode", "8f14e45f-ceea-467f-a0e6-0b3c1a2b3c4d"]


@pytest.mark.parametrize("native", AWKWARD)
def test_a_key_survives_the_round_trip(native):
    key = I.CellKey("HDL115pancLN_Scan1.er.qptiff - resolution #1", "qupath-0.7.0:instanseg", native)
    assert I.parse_cell_uid(key.uid) == key
    assert key.uid.count(I.SEPARATOR) == 2


def test_a_separator_inside_a_part_cannot_make_two_cells_one():
    assert I.cell_uid("a:b", "c", "d") != I.cell_uid("a", "b:c", "d")
    assert I.cell_uid("a", "b", "c:d") != I.cell_uid("a", "b:c", "d")


def test_native_ids_that_collide_across_acquisitions_are_two_cells():
    """Xenium numbers its cells per bundle: the same cell_id in two bundles is two cells."""
    table = pd.DataFrame({"acquisition_id": ["xenium.0041117", "xenium.0041116"],
                          "segmentation_id": ["xoa-2.0", "xoa-2.0"],
                          "native_cell_id": ["aaaabbbb-1", "aaaabbbb-1"]})
    table["cell_uid"] = I.cell_uids(table)
    assert table["cell_uid"].nunique() == 2
    I.check_cell_keys(table)


def test_a_duplicated_key_is_refused_with_its_count():
    table = pd.DataFrame({"acquisition_id": ["a", "a", "a"], "segmentation_id": ["s"] * 3,
                          "native_cell_id": ["1", "1", "2"]})
    table["cell_uid"] = I.cell_uids(table)
    with pytest.raises(I.IdentityError, match=r"1 cell_uid\(s\) appear more than once \(2 rows\)"):
        I.check_cell_keys(table)


@pytest.mark.parametrize("bad", [None, "", float("nan"), 7])
def test_a_missing_or_non_text_part_is_refused(bad):
    with pytest.raises(I.IdentityError):
        I.CellKey("a", "s", bad)
    table = pd.DataFrame({"acquisition_id": ["a"], "segmentation_id": ["s"],
                          "native_cell_id": [bad]})
    with pytest.raises(I.IdentityError):
        I.cell_uids(table)


@pytest.mark.parametrize("uid", ["a%3ab:s:n", "%61:s:n", "a:s", "a:s:n:extra", "a::n"])
def test_a_uid_that_is_not_canonical_or_not_three_parts_is_refused(uid):
    with pytest.raises(I.IdentityError):
        I.parse_cell_uid(uid)


def test_the_table_builder_matches_the_one_key_builder():
    table = pd.DataFrame({"acquisition_id": ["acq one"] * len(AWKWARD),
                          "segmentation_id": ["seg"] * len(AWKWARD), "native_cell_id": AWKWARD})
    assert list(I.cell_uids(table)) == [I.cell_uid("acq one", "seg", n) for n in AWKWARD]


def test_a_uid_that_does_not_encode_its_own_row_is_refused():
    table = pd.DataFrame({"acquisition_id": ["a", "b"], "segmentation_id": ["s", "s"],
                          "native_cell_id": ["1", "2"]})
    table["cell_uid"] = [I.cell_uid("a", "s", "1"), I.cell_uid("a", "s", "2")]
    with pytest.raises(I.IdentityError, match="do not encode their row"):
        I.check_cell_keys(table)


def _release_rows() -> pd.DataFrame:
    table = pd.DataFrame({"acquisition_id": ["a", "a"], "segmentation_id": ["s", "s"],
                          "native_cell_id": ["1", "2"]})
    table["cell_uid"] = I.cell_uids(table)
    for field in I.IDENTITY_FIELDS:
        if field not in table:
            table[field] = pd.NA
    table["tissue"] = ["pancreas", pd.NA]
    return table


def test_identity_columns_allow_unknown_but_not_absent():
    table = _release_rows()
    I.check_identity_columns(table)                     # NA = recorded as unknown: allowed
    with pytest.raises(I.IdentityError, match="not carried: block_id"):
        I.check_identity_columns(table.drop(columns=["block_id"]))


def test_a_tissue_that_is_not_a_tissue_id_is_refused():
    table = _release_rows()
    table.loc[1, "tissue"] = "pLN"                      # a spelling, not the id
    with pytest.raises(I.IdentityError, match="not tissue ids"):
        I.check_identity_columns(table)


def test_every_declared_field_says_what_it_means():
    assert set(I.REQUIRED_FIELDS) <= set(I.IDENTITY_FIELDS)
    assert all(text.strip() for text in I.IDENTITY_FIELDS.values())
