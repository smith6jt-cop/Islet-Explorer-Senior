"""Cell identity: one key per cell, made from what was recorded, never from a filename.

A cell is the cell its segmentation drew in one acquisition. So its identity is three things, and
nothing less is unique: the **acquisition** it was measured in (one PhenoCycler scan, one Xenium
bundle, one IF image), the **segmentation** that drew it, and its **native id** inside that
segmentation (QuPath's UUID, kept as ``object_id``; Xenium's ``cell_id``). Two acquisitions may
reuse native ids -- Xenium numbers its cells per bundle -- and a native id outside its
acquisition names no cell at all.

``cell_uid`` writes the three as one string, ``acquisition:segmentation:native``, each part
percent-encoded so that a ``:`` (or anything else) inside an id can never make two different
cells share a key, and parsed back exactly. Uniqueness is checked, not assumed.

The other identity fields a cell carries (``IDENTITY_FIELDS``) are declared here once. A field
that was not recorded is missing -- unknown stays unknown -- and is never filled with a guess.

MultiPlex5000 D2 (2026-10-20). The acquisition and segmentation ids are assigned by the source
inventory (stage 00, D3) from recorded identifiers; the stages adopt the key from D4 on.
"""

from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import quote, unquote

import pandas as pd

#: Joins the three parts. Inside a part it is always encoded, so splitting on it is exact.
SEPARATOR = ":"
#: Characters a part keeps as they are (RFC 3986 unreserved); everything else is encoded.
_UNRESERVED = "-._~"
#: The three parts, in order, under the column names a table carries them by.
KEY_PARTS = ("acquisition_id", "segmentation_id", "native_cell_id")

#: The identity a release carries for every cell, and what each field means.
IDENTITY_FIELDS: dict[str, str] = {
    "cell_uid": "acquisition:segmentation:native, each part percent-encoded (cell_uid())",
    "acquisition_id": "the acquisition the cell was measured in: one scan, bundle or image",
    "segmentation_id": "the segmentation that drew the cell in that acquisition",
    "native_cell_id": "the id that segmentation gave it (QuPath's UUID as object_id; Xenium's "
                      "cell_id)",
    "donor_id": "the donor, as recorded; aliases resolved by an explicit mapping, never by "
                "reading a filename",
    "block_id": "the tissue block the section was cut from",
    "section_id": "the physical section on the slide",
    "tissue_roi_id": "the annotated tissue region the cell lies in (its annotation UUID)",
    "tissue": "a tissue id from the tissue profiles",
    "assay_order": "the acquisition's place among the assays run on that section, 1 first",
    "panel_id": "the marker or gene panel of the acquisition",
    "frame_id": "the coordinate frame the cell's geometry is in",
    "source_id": "the recorded source the row was read from (path and digest)",
}

#: The fields that must be present on every cell; the rest may be unknown (missing).
REQUIRED_FIELDS = ("cell_uid", *KEY_PARTS)


class IdentityError(ValueError):
    """A cell key that is missing, malformed, non-canonical or not unique."""


def _encode(part: str) -> str:
    return quote(part, safe=_UNRESERVED)


def _require_part(name: str, value) -> str:
    if value is None or (not isinstance(value, str)) or value == "":
        raise IdentityError(f"{name} must be a non-empty string, got {value!r}")
    return value


@dataclass(frozen=True)
class CellKey:
    """The three recorded parts of a cell's identity."""

    acquisition_id: str
    segmentation_id: str
    native_cell_id: str

    def __post_init__(self) -> None:
        for name in KEY_PARTS:
            _require_part(name, getattr(self, name))

    @property
    def uid(self) -> str:
        return SEPARATOR.join(_encode(getattr(self, name)) for name in KEY_PARTS)

    @classmethod
    def parse(cls, uid: str) -> "CellKey":
        """The key a ``cell_uid`` encodes. Refuses a uid that is not in canonical form, so one
        cell has one spelling."""
        if not isinstance(uid, str):
            raise IdentityError(f"a cell_uid is a string, got {uid!r}")
        parts = uid.split(SEPARATOR)
        if len(parts) != len(KEY_PARTS):
            raise IdentityError(f"{uid!r}: a cell_uid has {len(KEY_PARTS)} parts, "
                                f"acquisition:segmentation:native")
        decoded = [unquote(p) for p in parts]
        key = cls(*decoded)
        if key.uid != uid:
            raise IdentityError(f"{uid!r} is not in canonical form; it is {key.uid!r}")
        return key


def cell_uid(acquisition_id: str, segmentation_id: str, native_cell_id: str) -> str:
    """The cell's key as one string: ``acquisition:segmentation:native``."""
    return CellKey(acquisition_id, segmentation_id, native_cell_id).uid


def parse_cell_uid(uid: str) -> CellKey:
    return CellKey.parse(uid)


def cell_uids(table: pd.DataFrame) -> pd.Series:
    """``cell_uid`` for every row of a table that carries the three key parts as columns.

    A missing or empty part is an error naming how many rows lack it: a cell without its
    acquisition, segmentation or native id is not a cell this pipeline can name.
    """
    missing = [c for c in KEY_PARTS if c not in table.columns]
    if missing:
        raise IdentityError(f"the table has no {', '.join(missing)} column")
    encoded = []
    for column in KEY_PARTS:
        values = table[column]
        bad = values.isna() | (values.astype(str) == "")
        if bad.any():
            raise IdentityError(f"{int(bad.sum())} row(s) have no {column}")
        if not all(isinstance(v, str) for v in values):
            raise IdentityError(f"{column} must hold strings (an id is text, not a number)")
        encoded.append(values.map(_encode))
    uids = encoded[0].str.cat(encoded[1:], sep=SEPARATOR)
    uids.name = "cell_uid"
    return uids


def check_cell_keys(table: pd.DataFrame, *, examples: int = 5) -> None:
    """Validate a table's cell keys: present, matching their parts, and unique.

    Raises on a missing ``cell_uid``, a uid that does not encode the row's own parts (when the
    parts are carried), and any uid that appears twice -- the last names how many and shows a
    few, because a duplicated key means two rows claim to be one cell.
    """
    if "cell_uid" not in table.columns:
        raise IdentityError("the table has no cell_uid column")
    uids = table["cell_uid"]
    if uids.isna().any():
        raise IdentityError(f"{int(uids.isna().sum())} row(s) have no cell_uid")
    if all(c in table.columns for c in KEY_PARTS):
        expected = cell_uids(table)
        wrong = expected.values != uids.values
        if wrong.any():
            first = uids[wrong].iloc[0]
            raise IdentityError(f"{int(wrong.sum())} cell_uid(s) do not encode their row's "
                                f"parts; the first is {first!r}")
    dup = uids[uids.duplicated(keep=False)]
    if len(dup):
        shown = ", ".join(repr(u) for u in dup.drop_duplicates().head(examples))
        raise IdentityError(f"{dup.nunique()} cell_uid(s) appear more than once "
                            f"({len(dup)} rows), e.g. {shown}")


def check_identity_columns(table: pd.DataFrame, known_tissues=None) -> None:
    """A release's identity columns: all declared, the key parts filled, the tissue a known id.

    A field other than the key may be missing on a row (unknown); it may not be absent as a
    column, so a consumer can always tell "not recorded" from "not carried".
    """
    absent = [f for f in IDENTITY_FIELDS if f not in table.columns]
    if absent:
        raise IdentityError(f"identity column(s) not carried: {', '.join(absent)}")
    check_cell_keys(table)
    if known_tissues is None:
        from .integration import tissues as _tissues

        known_tissues = _tissues.TISSUES
    tissue = table["tissue"].dropna()
    unknown = sorted(set(tissue) - set(known_tissues))
    if unknown:
        raise IdentityError(f"tissue value(s) that are not tissue ids: {unknown}")
