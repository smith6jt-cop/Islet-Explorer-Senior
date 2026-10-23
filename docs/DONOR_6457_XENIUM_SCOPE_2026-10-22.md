# 6457's Xenium on the next run of each producer (2026-10-22)

The maintainer, 2026-10-22: "the xenium for 6457 is not excluded". 6457's recorded reason, poor
DAPI that defeated the PhenoCycler channel registration, is a PhenoCycler one, so
`cohort.EXCLUDED_PLATFORMS` scopes it to `phenocycler` (PR #38). This change moves every reader
of Xenium data to `cohort.excluded_on(donor, "xenium")`, and post-Xenium IF to `"if"`.
PhenoCycler readers still refuse 6457. **Nothing has been run.** No producer, release or layer
was rebuilt, so every output below is unchanged on disk until its producer runs again.

## What 6457 has, measured 2026-10-22

- `/data/xenium/h5ad_out/`: `6457_PancXenium.h5ad` (692,855 cells), `6457a_pLNXenium.h5ad`
  (23,530) and `6457b_pLNXenium.h5ad` (202,011), all on the 5,101-gene panel.
- `xenium_paths.csv`: three rows (`panc`, `pln_1`, `pln_2`), all on slide `0037637`, run
  `MM038_Xenium_nPodBatch4`. **The bundle directories are not on this machine.**
  `/data/xenium/HuPrime5kPanTissue/6KTRAF/MM038_Xenium_nPodBatch4/` holds only 6476's two
  bundles. The complete outputs are on HiPerGator `/orange`.
- The integration manifest records the three 6457 rows as `xenium_only`. `cells_xen` has no
  6457 partition, because the ingest imports only paired rows unless it is given `--all`.
- No `restain_sections.csv` or `restain_registration_inputs.csv` row, and no VALIS directory.
- No partition in either current Xenium release (pancreas `panc_2026-09-17`, 21 donors;
  pLN `release_pln_immune_v3`, 20 donors).

## Per producer: what changes on its next run

| producer | outputs that gain 6457 |
|---|---|
| `analysis/xenium/scripts/panc_xenium_cluster_name.py`, run by `panc_xenium_release.py` | The pancreas release gains `calls/donor_id=6457/roi=panc` (692,855 cells), a `per_donor_summary.csv` row, and `lineage_counts.csv` / `coherence_genes.csv` rows; `donors` and `n_sections` go from 21 to 22. **With no local bundle, 6457's nucleus rule is UNMEASURABLE**, so the release's cohort `n_cells_no_nucleus_excluded` (775,234 today, every section measured) becomes null, and the partial sum is published on its own. Transfer the bundle first if that total should stand. |
| `analysis/xenium/scripts/pln_immune_typing_v3.py` | `partition/` and `calls/` gain sections 6457a and 6457b (225,541 cells). Both get rows in `per_donor_summary.csv` and `per_donor_tnk.csv`, and `release.json`'s `donors_refused_by_cohort` goes from `["6457"]` (as it reads today) to `[]`. The nucleus rule is UNMEASURABLE for both sections, for the same reason. In `cross`, 6457 has no PhenoCycler pLN to compare against, and that absence is named. |
| `core/multiplex5000/scripts/pln/cluster_pln.py --modality xenium` | `typing/pln/clusters/xenium/` gains profiles and assignments for 6457a and 6457b, and `cluster_report.json`'s `donors_refused_by_cohort` (`["6457"]` today) becomes `[]`. The PhenoCycler side is unchanged, because the pLN ingest never wrote 6457. |
| `core/multiplex5000/scripts/xenium_typing/pins.py` (quarantined lane) | **Nothing by default.** Its cohort is the PhenoCycler QuPath manifest's `expected_donors`, which does not list 6457. A manifest that did list it would now be accepted. |
| `core/multiplex5000/scripts/xenium_celltype/build_phenotyped_base.py` | **Nothing by default.** Its cohort is `cohort-spec.json`, which does not list 6457. `--donors 6457` is now accepted, but it needs the bundle, which is not local. |
| `multiplex5000.integration.phenospace` | Writes nothing. A Xenium identity release that declares 6457 now loads; a PhenoCycler one still raises. |
| `scripts/senior/build_sections_layer.py` | `sections/data.parquet` gains three Xenium rows for 6457 (serial `0037637`: panc, pln_1, pln_2), and no `restain_if` row. |
| `scripts/senior/build_availability.py` | `availability/data.parquet` gains donor 6457 with four rows. `phenocycler_panc`, `islet_outline` and `trajectory` are `withheld`, with reason `centrally_excluded_on_phenocycler: <cohort reason>`. `restain_if` is `unavailable / no_restain_pairing`. |
| `scripts/senior/build_restain_layer.py` | `restain/availability.parquet` gains one 6457 status row, unpaired: there is no `cells_xen` partition and no warp or section. It gains no overlay cells. |
| `scripts/senior/build_distance_profiles.py` | Once the pLN release carries 6457, the pLN Xenium arm gains compositional rows for 6457/pln_1 and pln_2. Once the pancreas release carries it, the pancreas Xenium arm gains one 6457 row, which stays unavailable at link 2 (`XEN_PANC_NO_COORDINATES`) while `cells_xen` has no 6457 partition. The PhenoCycler arms are unchanged. |
| `scripts/senior/build_grid_composition.py` | The universe gains 6457 once a Xenium release carries it. In `presence.csv`, its protein pancreas and protein pLN cells are `withheld` with the cohort reason, and their files are never read. Its RNA cells are measured once the releases are rerun (`blank` until then). `composition.csv` and `fine.csv` gain its RNA rows. `crossmodal_tnk` gains nothing, since that needs both modalities. Its RNA batch comes from the pancreas release's `per_donor_summary.csv`, as every Xenium donor's does. The universe gains 6457 only through a release: the domain's `donor_metadata.csv` drops it. |
| `scripts/senior/modality_concordance.py` | Once the releases carry it, `pln_per_donor_counts_rna.csv`, `pancreas_per_donor_counts_rna.csv` and `pln_rna_controls.csv` gain 6457 rows. Every comparison is unchanged, because a comparison joins the protein and RNA sides on `donor_id` and 6457 has no protein side. `manifest.json` gains `donor_exclusions_by_modality`. |
| `scripts/senior/same_section_linkage.py` | `donor_linkage.csv` gains an `excluded_platforms` column. 6457 still gets no row: it has no VALIS directory and no same-slide acquisition. |

## Left as they were

These read PhenoCycler data only, or gate on a PhenoCycler cohort:

- `ingest_pln_cells.py` and `build_pln_layer.py` read the PhenoCycler pLN export and roster.
- `build_islet_summary.py` and `build_trajectory_layer.py` read PhenoCycler typing.
- `build_all_tissuumaps_projects.py` reads qptiffs.
- `build_compartment_gradient.py` gates on `typing/cells` partitions. Its Xenium pLN rows are kept
  only for that cohort, so 6457's will be listed under `pln_xenium_donors_outside_cohort`, like
  the other Xenium-only donors.
- `build_cohort_roster.py`'s eligibility: the candidates come from `cohort-spec.json`, the
  PhenoCycler pancreas roster, so 6457 stays an excluded entry there. Since 2026-10-23 the entry
  also records which platforms its exclusion covers; see below.

## Moved on 2026-10-23, at the maintainer's word ("Yes stop refusing obviously. Yes record")

- **Re-stain registration.** `core/multiplex5000/scripts/restain/registration_inputs.py`
  (`select_registration_inputs`, which every registration step goes through) and
  `valis_register.py` now ask for the Xenium and IF they read. 6457 is selected the same way
  as any donor, but neither re-stain inventory (`restain_registration_inputs.csv`,
  `restain_sections.csv`) has a 6457 row and there is no VALIS directory for it. So nothing is
  registered until 6457's post-Xenium IF qptiff is recorded there and its Xenium bundle is
  local.
- **The app roster.** `build_cohort_roster.py` writes `platforms` (the vocabulary,
  `cohort.PLATFORMS`) and, on each excluded entry, the platforms its exclusion covers. 6457's
  entry says `["phenocycler"]`. `apps/islet/config/cohort.json` was rebuilt, and the rebuild
  also brings donor 112 into `all_donors` and `eligible_donors`: `cohort-spec.json` has listed
  `0112` since after the last build, on 2026-09-17. `ieCore`'s `ie_cohort_filter` /
  `ie_cohort_is_excluded` / `ie_cohort_excluded_ids` take an optional `platform`. With none
  given they remove every excluded donor, as before. No app code calls them outside the tests.
- **`cohort.excluded_on` refuses an unknown platform name.** A misspelling such as
  `"PhenoCycler"` matched no scoped exclusion, and so readmitted 6457.

## Still refusing 6457

- `run_pheno_annotation_export.py` carries PhenoCycler annotations onto Xenium, so it is right
  to refuse 6457.
- `card_x1_gradient.R` selects the pancreas arm by donor alone (the xfail in
  `test_distance_profiles.py`). A donor with a Xenium row and no PhenoCycler row, as 6457 would
  be, shows its Xenium absence there.
