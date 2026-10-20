"""
phenocycler.integration — PhenoCycler (protein) <-> Xenium (RNA) integration.

Both modalities first enter a versioned cell-table contract as geometry. PhenoCycler reads
an immutable core run through ``expression`` and a same-run UUID area join; Xenium ignores
embedded labels. Identity is attached only from a declared release with an exact native-key
call set and a strict label space. Retired phenotype layouts have no compatibility opt-in.

Two modes, and the distinction is the whole point:

``sequential``
    Serial sections cut from one FFPE block. The two modalities do **not** share cells, so
    cell-to-cell pairing is invalid no matter how good the registration is. Integration
    happens at four honest levels — donor, islet/structure, niche, and (explicitly labelled
    as inference) pseudo-cell.

``same_slide``
    One physical section imaged twice. Cells *are* the same cells, so a genuine paired
    protein+RNA matrix is recoverable.

``sameslide.py`` and the sequential-only steps refuse to run in the wrong mode. That guard
is a hard error rather than a warning because the validity boundary of every downstream
claim rests on it.

Stage order (``pipeline.py``)::

    manifest      S0   donor <-> Xenium-run pairing, donor_id/roi stamping
    export_pheno  S1a  immutable-run geometry; optional declared identity -> cells_pheno/
    import_xenium S1b  h5ad | zarr | bundle geometry; optional identity   -> cells_xen/
    postxen       S1c  re-stain pairing                         [Xenium identity]
    vocab         S2   release-declared vocabulary + protein<->gene crosswalk
    structures    S3   tissue structures                        [identity, pancreas]
    register      S4   coarse rigid -> affine -> optional non-rigid       [geometry]
    match         S6   islet <-> islet assignment               [identity, sequential]
    cellmatch     S6b  endocrine cell matching                  [identity, sequential]
    sameslide     S11  cell <-> cell assignment                 [identity, same_slide]
    grid          S7   composition niches + registered grid     [identity]
    donor         S8   donor-level concordance                  [identity]
    crossmodal    S9   pseudo-cell linking                      [identity, inference]
    qc/figures    S10  geometry always; identity sections only when attached

Everything after S1 speaks the modality-agnostic cell-table contract in ``contract.py``.
Parquet metadata distinguishes ``geometry`` from ``identity`` and stamps the release id, so
nullable identity columns cannot masquerade as verified calls.

The heavy readers (anndata / spatialdata / SimpleITK / opencv) are optional: install with
``pip install -e '.[integration]'`` or use ``environment-integration.yml``. Importing this
package does not pull them in — each module raises a clear error naming the extra only when
the dependency is actually needed, so the core pipeline keeps working in the lean env.
"""

from __future__ import annotations

from .contract import (
    CELL_TABLE_COLUMNS,
    MODALITIES,
    CellTable,
    feature_columns,
    read_cell_table,
    validate_cell_table,
    write_cell_table,
)

__all__ = [
    "CELL_TABLE_COLUMNS",
    "MODALITIES",
    "CellTable",
    "feature_columns",
    "read_cell_table",
    "validate_cell_table",
    "write_cell_table",
]
