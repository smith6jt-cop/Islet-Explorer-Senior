"""Portable per-gene Xenium detection evidence, without a cell-type claim.

This extracts the detection predicate recorded by the current pLN RNA lane:
``log1p(raw_gene_count * 1e4 / cell_gene_depth) >= log1p(1)``. The input is
the accepted region's raw gene-expression matrix. An unmeasured gene, a cell
without a nucleus, and a zero-depth cell are unavailable, never negative.
Tissue naming and cluster evidence remain separate Stage 06 decisions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

TARGET_SUM = 10_000.0
UNAVAILABLE = np.int8(-1)
UNDETECTED = np.int8(0)
DETECTED = np.int8(1)


class RNAEvidenceError(ValueError):
    """The raw matrix or cell admission cannot support detection evidence."""


@dataclass(frozen=True)
class RNAEvidence:
    """One tri-state array per requested gene, ordered like the input cells."""

    states: dict[str, np.ndarray]
    missing_genes: tuple[str, ...]
    eligible: np.ndarray
    depth: np.ndarray


def detect_genes(counts: sparse.csr_matrix, genes: tuple[str, ...],
                 requested: tuple[str, ...], nucleus_eligible: np.ndarray) -> RNAEvidence:
    """Detect selected genes at fixed 1e4 depth, preserving unavailable cells.

    The comparison is algebraically equal to the recorded log1p predicate and
    avoids materializing an all-gene normalized matrix. Gene availability is
    based on the assayed feature list, before looking at any cell's counts.
    """
    if not sparse.isspmatrix_csr(counts):
        raise RNAEvidenceError("counts must be a raw cells-by-genes CSR matrix")
    if counts.shape[1] != len(genes) or len(set(genes)) != len(genes):
        raise RNAEvidenceError("gene feature names do not match unique matrix columns")
    if not requested or len(set(requested)) != len(requested):
        raise RNAEvidenceError("requested genes must be nonempty and unique")
    nuclei = np.asarray(nucleus_eligible)
    if nuclei.shape != (counts.shape[0],) or nuclei.dtype != np.bool_:
        raise RNAEvidenceError("nucleus eligibility must be one Boolean per cell")
    if not np.issubdtype(counts.dtype, np.integer) or np.any(counts.data < 0):
        raise RNAEvidenceError("raw gene counts must be nonnegative integers")
    depth = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()
    if not np.isfinite(depth).all():
        raise RNAEvidenceError("cell gene depth is nonfinite")
    eligible = nuclei & (depth > 0)
    positions = {gene: index for index, gene in enumerate(genes)}
    measured = [gene for gene in requested if gene in positions]
    # A Stage 06 panel can request hundreds of genes from a large raw matrix.
    # Slice once, then use CSC's column index instead of rescanning the full
    # cells-by-genes CSR matrix separately for every requested gene.
    selected = counts[:, [positions[gene] for gene in measured]].tocsc()
    selected_position = {gene: index for index, gene in enumerate(measured)}
    states: dict[str, np.ndarray] = {}
    missing = []
    for gene in requested:
        state = np.full(counts.shape[0], UNAVAILABLE, dtype=np.int8)
        if gene not in selected_position:
            missing.append(gene)
        else:
            raw = selected.getcol(selected_position[gene]).toarray().ravel().astype(np.float64)
            state[eligible] = np.where(raw[eligible] * TARGET_SUM >= depth[eligible],
                                       DETECTED, UNDETECTED)
        states[gene] = state
    return RNAEvidence(states, tuple(missing), eligible, depth)
