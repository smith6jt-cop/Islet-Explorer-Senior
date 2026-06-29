# Archived legacy phenotyping files (2026-06-29)

These files were the **predecessors** of the current unified Senior phenotyping pipeline. They are
kept here for provenance only — **none are part of the live pipeline**, and nothing in
`scripts/senior/` imports them. Git history is preserved (moved with `git mv`).

## Why they were retired

The Senior 59-plex data is phenotyped by a single unified pipeline (Jun 2026):

```
REDSEA (pixel spillover)  →  RESTORE (per-image norm)  →  broad-lineage (hierarchical)  →  QuPath UUID
scripts/senior/redsea_full.py   restore_normalize.py        assign_broad_lineage.py        export_broad_class_for_qupath.py
```

This replaced three earlier approaches that all produced the fake **"Unassigned"** fraction the new
pipeline was built to eliminate (the user mandate is *fix the upstream cause, never invent a junk
category*): per-donor **auto-gating**, global **Leiden clustering**, and an **Astir/consensus**
multi-typer framework.

## What moved here, and what replaced it

| Archived file | Was | Superseded by |
|---|---|---|
| `auto_gate_per_donor.py` | Phase 0 per-donor GMM auto-gating | RESTORE per-image normalization (`restore_normalize.py`) |
| `phenotype_scimap.py` | Phase 0 scimap seeding + napariGater | RESTORE + hierarchical `assign_broad_lineage.py` |
| `cluster_annotate.py` | Phase 1 Leiden clustering + rule annotation | direct marker-rule typing in `assign_broad_lineage.py` |
| `finalize_phenotype.py` | Phase 1 classifier propagation + multi-method fusion | single hierarchical pass (no consensus needed) |
| `run_astir.py` | Phase 1 Astir probabilistic typer | — (consensus framework dropped) |
| `astir_markers.yaml` | Astir marker spec | — |
| `rules_to_astir_yaml.py` | rules CSV → Astir YAML | — |
| `export_phenotype_for_qupath.py` | centroid/measurement CSV export for QuPath | `export_broad_class_for_qupath.py` (UUID-keyed, exact) |
| `set_broad_lineage_class.groovy` | measurement-based PathClass setter | `import_broad_lineage.groovy` (UUID match, sets class + color) |
| `import clusters and cell neighborhoods.groovy` | CytoMAP-style centroid importer | `import_broad_lineage.groovy` |
| `make_phenotype_rules.py` | generated `phenotype_rules_senior.csv` | broad-lineage uses hardcoded hierarchical rules |
| `marker_reliability_audit.py` | 59-plex gating-failure audit | panel finalized; one-off diagnostic |
| `redsea_approx.py` | table-level spillover **proof** (centroid adjacency) | pixel-level `redsea_full.py` (actual correction) |
| `validate_phenotypes.py` | Phase 2 post-clustering validation | `redsea_full_validate.py` / `redsea_cohort_validate.py` (live) |
| `phenotype_rules_senior.csv` | legacy hierarchical rules (scimap format) | hardcoded rules in `assign_broad_lineage.py` |
| `CODEX_Panc_scimap_Analysis.ipynb` | Junior-era scimap phenotyping notebook | this pipeline / `notebooks/senior_phenotyping_redsea_restore.ipynb` |
| `CODEX_Panc_Analysis.ipynb` | Junior-era CODEX analysis notebook | — |
| `napari_gating_senior.ipynb` | manual napariGater gate-refinement notebook | automated RESTORE gating |
| `Comb_Xen_Phen.ipynb` | experimental Xenium-combine phenotyping notebook | — |
| `phenotype_rules1.csv`, `phenotype_rules2.csv` | older rules iterations | — |

## Still live (NOT archived)

`redsea_full.py`, `restore_normalize.py`, `assign_broad_lineage.py`,
`export_broad_class_for_qupath.py`, `redsea_full_validate.py`, `redsea_cohort_validate.py`,
`scripts/groovy/import_broad_lineage.groovy`, `scripts/groovy/view phenotypes.groovy`, and the
`build_*/compute_*/aggregate_*/extract_*` app-data scripts. The scVI/trajectory notebooks
(`single_cell_analysis/SCVI_CODEX_*`, `notebooks/rebuild_trajectory.ipynb`) were left in place — they
are not phenotyping.

See `scripts/CLAUDE.md` for the live pipeline reference.
