#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan that is also the
NAMED PRODUCER of a live artifact under ``data/app_data`` (code tier 2), which is why it is
kept rather than deleted. It cannot be re-run: the inputs named below are absent from this
checkout.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method, and ``data/app_data``
is the OTHER app domain: the two datasets share five donors and no vocabulary.

3B(+) Senior neighborhood metrics — reuses scripts/compute_neighborhood_metrics.py
(Junior's, computes per-islet peri-islet composition / immune infiltration /
enrichment-z / min-distance) with Senior's phenotype labels and data paths. No
groovy. Run in the `scimap` conda env.

Output: data/neighborhood_metrics.csv
"""
import sys
from pathlib import Path

SEN = Path(__file__).resolve().parents[2]  # repo root: scripts/senior/<file> -> ../../
sys.path.insert(0, str(SEN / "scripts"))
import compute_neighborhood_metrics as cnm  # noqa: E402

SD = Path(__file__).resolve().parents[2] / "data"

# Senior's immune phenotype labels (drives the aggregate immune_* metrics).
cnm.IMMUNE_TYPES = [
    "CD8 T cell", "CD4 T cell", "T cell", "Regulatory T cell", "B cell",
    "Macrophage", "M2 macrophage", "Dendritic cell", "NK cell", "Granulocyte", "Immune",
]
# Per-type metrics for all phenotypes except the ambiguous scimap bucket.
cnm.PER_TYPE_EXCLUDE = {"likely-Regulatory T cell"}
# Senior centroids are already in µm (QuPath "Centroid X µm"); script computes
# distances in native coordinate units, so no pixel conversion needed.

if __name__ == "__main__":
    cnm.compute_metrics(
        str(SD / "singlecell_protein_phenotyped.h5ad"),
        str(SD / "islets_core_senior.h5ad"),
        str(SD / "neighborhood_metrics.csv"),
    )
