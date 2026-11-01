"""Per-image positivity and cell identity, ported from the CODEX Islet-Explorer typing method.

Source: Islet-Explorer ``scripts/uniform/typing/`` at commit :data:`SOURCE_COMMIT` (branch
``spatial-review-20261020``, owner sign-off 2026-11-01). The method moved intact: each marker's level
per section from its own background or known-negative reference (``levels``), a stricter canary level
(``canary``), the share of a cell's pixels above a level (``coverage``), and the canary-ranked
compartment, lineage and label of every cell (``identity``). What the code may use comes from a
registry per panel and tissue (``registry``); the code names compartments and lineages, never markers.

What changed in the move is input and output only:
- inputs come from a :class:`~multiplex5000.positivity.sources.Source`. IES runs read Protein Measurement
  (``redsea`` and ``expression`` partitions, the manifest's geometry and qptiff); a CODEX run directory
  can be read for the computation-parity check;
- the tissue-specific parts (start class, the reference populations a registry may name, the sections
  methods are chosen on, the structure objects) come from a tissue hook in ``multiplex5000.tissues``;
- functions return tables instead of writing a run directory; the stage that calls them records them.
"""

SOURCE_REPOSITORY = "smith6jt-cop/Islet-Explorer"
SOURCE_COMMIT = "578767eeaf0f95bb1df7272238f361888a0ae75d"
SOURCE_PATH = "scripts/uniform/typing"
