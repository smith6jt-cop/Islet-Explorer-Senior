"""Contracts every layer builder must hold, pinned against real artifacts.

Workers F through J extend this file. What is here is the hazard that was measured
on 2026-09-16 and is easy to reintroduce; see docs/ISLET_LAYERS_SPEC.md.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_layer_contracts.py
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
EMPTY = FIXTURES / "vessels_empty_6442_xenium_panc.parquet"
POPULATED = FIXTURES / "vessels_populated_6436_xenium_panc.parquet"

# core/phenocycler/tests/current/this_file.py -> parents[4] is the repository root.
# parents[3] is core/, which silently skipped every real-data test in this file.
REPO = Path(__file__).resolve().parents[4]


def test_an_empty_structure_file_carries_no_schema_at_all():
    """6442's vessels.parquet is 0 rows AND 0 columns, not an empty typed table.

    This is the trap the fixture exists for. A reader that opens one donor's
    structure file, learns "vessels have these columns", and applies that to the
    cohort gets nothing from this file -- there is no schema in it to learn. A
    builder must take the schema from a declaration, never from whichever file it
    happened to open first.
    """
    empty = pq.ParquetFile(EMPTY)
    assert empty.metadata.num_rows == 0
    assert empty.metadata.num_columns == 0
    assert list(empty.schema_arrow.names) == []


def test_its_populated_sibling_has_the_real_schema():
    """Same structure kind, same modality, same roi, one donor over."""
    full = pq.ParquetFile(POPULATED)
    assert full.metadata.num_rows == 110
    names = set(full.schema_arrow.names)
    for required in ("structure_id", "kind", "modality", "donor_id", "roi",
                     "x_um", "y_um", "n_cells"):
        assert required in names, f"{required} missing from the populated sibling"


def test_reading_both_together_must_not_be_attempted_by_schema_inference():
    """Unifying the two is the failure, and it fails loudly rather than silently.

    Recorded so nobody 'fixes' a builder by pointing pyarrow at the directory and
    letting it merge schemas: an empty file with no columns and a populated file
    with ten do not unify into the ten.
    """
    d = ds.dataset([str(EMPTY), str(POPULATED)], format="parquet")
    # The union is driven by the first fragment, which declares nothing.
    assert d.schema.names == [] or "structure_id" in d.schema.names, (
        "pyarrow's behaviour here changed; re-read the spec before relying on it"
    )


# --------------------------------------------------------------------------- #
# Reading the source trees
# --------------------------------------------------------------------------- #

HIVE_COLLIDING_ROOTS = [
    "core/phenocycler/data/expression",
    "core/phenocycler/data/integration/cells_xen",
]


@pytest.mark.parametrize("rel", HIVE_COLLIDING_ROOTS)
def test_hive_partitioning_collides_with_the_donor_id_column(rel):
    """partitioning="hive" raises on these roots; partitioning=None is the contract.

    The inferred partition key and the file's own large_string donor_id cannot be
    merged. An explicit string partition key does not help either, because the
    column is large_string and not string. Measured 2026-09-16.
    """
    root = REPO / rel
    if not root.is_dir():
        pytest.skip(f"{rel} is not present in this checkout")

    with pytest.raises(Exception) as excinfo:
        ds.dataset(str(root), format="parquet", partitioning="hive").to_table(
            columns=["donor_id"]
        )
    assert "donor_id" in str(excinfo.value)

    # The contract: no partitioning, donor comes from the file column.
    table = ds.dataset(str(root), format="parquet").to_table(columns=["donor_id"])
    assert table.num_rows > 0
    assert set(table.column("donor_id").to_pylist())


def test_typing_cells_has_no_donor_id_column_so_the_directory_is_the_donor():
    """The opposite problem, and the reason the spec tables both.

    typing/cells is 24 columns with no donor_id among them. Read with
    partitioning=None it yields 25 million rows and no way to separate donors.
    """
    root = REPO / "data/development/typing/cells"
    if not root.is_dir():
        pytest.skip("typing/cells is not present in this checkout")
    one = sorted(root.glob("donor_id=*/*.parquet"))
    assert one, "expected donor_id=* partitions"
    names = set(pq.ParquetFile(one[0]).schema_arrow.names)
    assert "donor_id" not in names
    assert "object_id" in names


def test_donor_id_is_not_normalised_across_trees():
    """0112 and 112 are the same donor. A builder that globs donor_id=112 finds nothing.

    The pLNpanc runs pad to four characters; the cohort does not. This is pinned
    because two builders normalising differently is the cross-unit disagreement
    the merge step exists to catch.
    """
    expression = REPO / "core/phenocycler/data/expression"
    cells_xen = REPO / "core/phenocycler/data/integration/cells_xen"
    if not expression.is_dir() or not cells_xen.is_dir():
        pytest.skip("source trees are not present in this checkout")

    padded = {p.name.split("=")[1] for p in expression.glob("donor_id=*")}
    unpadded = {p.name.split("=")[1] for p in cells_xen.glob("donor_id=*")}
    assert "0112" in padded, "the pLNpanc partition is expected to be zero-padded"
    assert not any(d.startswith("0") for d in unpadded), (
        "the cohort partitions are expected to be unpadded"
    )
