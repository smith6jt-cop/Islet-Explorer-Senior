# Vendored ontology + annotation provenance

Inputs vendored here so the gene-level identity classification
(`scripts/derive_identity_panels.py` → `data/panel_gene_roles.csv`,
`data/**/identity_core/`) is reproducible without network access and is
pinned to a known ontology release.

## Cell Ontology — `cl-basic.obo`

| Field | Value |
|---|---|
| File | `cl-basic.obo` |
| `data-version` | `cl/releases/2026-06-08/cl-basic.owl` |
| Release tag | `v2026-06-08` |
| Terms (CL) | 3,540 |
| Format | OBO 1.2 (the `-basic` variant: is_a-navigable, no cross-ontology imports) |
| Source | OBO Foundry Cell Ontology, <https://obophenotype.github.io/cell-ontology/> |
| Retrieved from | `https://github.com/obophenotype/cell-ontology/releases/download/v2026-06-08/cl-basic.obo` |

The canonical PURL is
`http://purl.obolibrary.org/obo/cl/cl-basic.obo`; the byte-identical release
asset above was used because the PURL host was unreachable from the build
environment. To refresh, download the `cl-basic.obo` asset of a newer
`obophenotype/cell-ontology` release, drop it here, update the table, and
re-run `python3 scripts/build_panel_roles.py` (anchor resolution fails loudly
if a term was renamed/obsoleted, so a regression surfaces immediately).

Only `[Term]` blocks, their `is_a:` edges, primary `name:`, and
`EXACT`/`NARROW` synonyms are parsed; obsolete terms are dropped. The classifier
maps each gene's `cell_type` terms to CL ids by name/synonym and tests
descendancy under per-lineage CL anchors (see `docs/panel_classification.md`).

## CellxGene cell-type annotation

The `cell_type` column the classifier reads from
`data/reference_5k/xenium5k_genes.csv` is **not** computed here — it is carried
verbatim from the 10x Genomics panel metadata.

| Field | Value |
|---|---|
| Annotation column | `cell_type` (`;`-separated CL term names) |
| Origin file | `XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv` (repo root) |
| Panel | 10x Genomics **Xenium Prime 5K Human Pan Tissue and Pathways** |
| Upstream source | CZ CELLxGENE Discover gene→cell-type associations (CL-termed), as shipped by 10x in the panel metadata |
| Coverage | 2,859 / 4,992 genes annotated (≈57%); ≈43% carry no `cell_type` |
| Term→CL match rate | 526 / 531 distinct terms map to a CL id |

Because ≈43% of the panel — including textbook markers (IAPP, MAFA, NKX6-1,
KRT19, CHGA, TRAC, …) — has no `cell_type`, the classifier merges a curated,
literature-backed canonical marker list (`CANONICAL` in
`scripts/derive_identity_panels.py`) that forces `identity` regardless of CL.
Each gene's `evidence` column records `cellxgene`, `canonical`, or `both`.
