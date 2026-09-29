# Declared producers: one per tissue x modality

Step 4 of the 2026-09-21 plan, scope **"Ledger only"** (maintainer, 2026-09-21). This page
declares the single producer of a class column for each of the four tissue x modality cells,
turns every other runnable typer in the tree into a recorded disposition with its measured
live uses, and records the one place a second method reached the app and how it was closed.

Measurements: the producers census of **2026-09-21** and the bare-name reference walk of
**2026-09-22** (116 ledger paths, classified LIVE / TEST-ONLY / RETIRED-SET-ONLY / DOC-ONLY /
NONE). Every `path:line` below was re-opened on **2026-09-22** at commit `f62dccd`; where the
census's line had drifted, the line printed here is the current one. Nothing on this page
moves a file. Quarantine runs through `scripts/retire/` in the main checkout after the branches
merge, and that execution is the integrator's.

The guard for this page is `core/phenocycler/tests/current/test_declared_producers.py`.

## The four producers

| tissue x modality | producer | command and interpreter | inputs | output root and class column | artifact-graph consumers |
|---|---|---|---|---|---|
| **pancreas, PhenoCycler** | `core/phenocycler/scripts/typing/subpanel.py` then `core/phenocycler/scripts/typing/build_cells.py`, with `core/phenocycler/scripts/typing/refine.py`'s pass joined back by `build_cells.py --refine` (ruled in 2026-09-23) | `.venv/bin/python core/phenocycler/scripts/typing/subpanel.py --membership geometry --emit-detail` (subpanel.py:125; the bare default is the reconstruction gate, not the rebuild -- subpanel.py:103-109), then `.venv/bin/python core/phenocycler/scripts/typing/build_cells.py [--donors ...] [--include-marginal]` (build_cells.py:67; typing/README.md:78), then `.venv/bin/python core/phenocycler/scripts/typing/refine.py --out data/development/typing/refine` (`--out` required, no default) and `.venv/bin/python core/phenocycler/scripts/typing/build_cells.py --include-marginal --refine data/development/typing/refine` (build_cells.py:507; the join at :212-277). The plain `build_cells.py` run is needed only after new `subpanel.py` calls, because refine.py reads `typing/cells` (refine.py:374-393). Interpreter: the repo `.venv` (Makefile:28; build_layers.py:116). | subpanel_io.py:68 `RUN` = the historical run `30f6689e4d27d295f890` (the only 20-donor expression output); :69 `ALIGNED` = data/development/celltype/uniform; :73 `SCALE_CSV` = data/development/broadtype/scratch/cosine_scale.parquet (read, never recomputed; subpanel_io.py:186-188); :82 `GEOJSON`. build_cells.py:517 reads data/donor_acquisition_runs.csv. | subpanel.py:1412 writes `S.OUT/subpanel_rebuilt`; build_cells.py:88 `OUT` = **data/development/typing/cells**, `donor_id=<d>/data.parquet`. Class column **`cell_class`** (subpanel.py:1214). Second pass, same lane: immune_subtype.py:98 `OUT` = data/development/typing/immune, column `immune_class` (:577); it never modifies `cell_class` (root CLAUDE.md:167-170). The refinement's verdict -- `refined_class`, `demoted`, `refine_note`, `refine_pass`, `concordance`, `concordance_before`, `immune_leaf` (build_cells.py:148) -- is joined BESIDE `cell_class` under `--refine`; `cell_class`, `is_endocrine` and `endocrine_route` stay the islet chain's, untouched (build_cells.py:212-277). The pass's own output is `data/development/typing/refine/donor_id=<d>/data.parquet` plus `refine_log.json`. **[ON DISK]** joined 2026-09-24 on all 20 partitions, `cell_class` unchanged; rerun 2026-09-26 after 6374 was retyped on its re-exported structure objects, and every layer and app artifact rebuilt from it the same day (Step 7) -- reading `cell_class`. | build_layers.py declares `typing/cells` as an input of fourteen nodes (:287, :304, :331, :339, :366, :397, :430, :437, :446, :472, :521, :540, :622, :648) and names build_cells.py as its producer (:151-153); `typing/immune` is an input of grid_composition (:649, producer at :156-157). Column readers: build_grid_composition.py:444 (`cell_class`), :589 (`immune_class`, `immune_subtype`); build_compartment_gradient.py:133 imports build_cells as the declared `cell_class` table, :150 `CELLS`. No R card reads `cell_class`: card_s3_compartments.R:21, card_g1_grid.R:33, card_x1_gradient.R:121 and card_s1_size_composition.R:19 read layers. |
| **pancreas, Xenium** | `analysis/xenium/scripts/panc_xenium_cluster_name.py`, run per section by `analysis/xenium/scripts/panc_xenium_release.py` | docs/PANC_XENIUM_TYPING_2026-09-17.md:650-655: `.venv/bin/python analysis/xenium/scripts/panc_xenium_release.py --work /data/xenium_typing/panc_2026-09-17/_work --release /data/xenium_typing/panc_2026-09-17/release --score-on naming --jobs 4`; one section alone at :661-662. Interpreter: the repo `.venv` in the write-up's own block, not scvi-env. The worker is loaded by path (panc_xenium_release.py:64-70), never imported by name. | panc_xenium_cluster_name.py:107 `H5AD_DIR` = /data/xenium/h5ad_out; :108-109 the identity cores under `core/panels/data`; :110 `XENIUM_ROOT` = /data/xenium; :111 xenium_paths.csv. | panc_xenium_release.py:683-684 default `--release` = **/data/xenium_typing/panc_2026-09-17/release** (written under /data because /home is at 97%, :32). `release.json` read 2026-09-22: `release_id` `xenium-panc-cluster-name-20260917`, `roi` `panc`, `written_utc` 2026-09-17T20:36:54+00:00. Class column **`lineage`** with **`common_lineage`** beside it (:552). It declares that it supersedes `data/development/xenium_typing/release`, the scANVI two-donor release (:86-96). | build_layers.py:651 declares the release as an input of grid_composition, the only absolute pin outside the repo (:670); :182-183 names panc_xenium_cluster_name.py. Resolved, never spelled, by `xenium_release.resolve` (xenium_release.py:690) in build_distance_profiles.py:488 and :1516, build_compartment_gradient.py:148, build_grid_composition.py:118. Column reader: build_grid_composition.py:515. Card: card_g1_grid.R:33 through grid_composition. |
| **pLN, Xenium** | `analysis/xenium/scripts/pln_immune_typing_v3.py` | Docstring :65-68 and docs/PLN_IMMUNE_XENIUM_2026-09-17.md:471-475: `.venv/bin/python analysis/xenium/scripts/pln_immune_typing_v3.py --stage partition`, then `--stage call`, `--stage tnk`, `--stage cross`, `--stage compare`. Every cohort-shaped write goes through `partial_run_guard.RunScope` (:47-48, :100). | :89 `H5AD_DIR` = /data/xenium/h5ad_out; :92 `SHARED_CLUSTERS` = data/development/typing/pln/clusters/xenium (read, deliberately not overwritten, :45); :93 `IDENTITY_CORE`; :94 `GDNA` triage parquet; :95 `XENIUM_ROOT`; :96 xenium_paths.csv; :97 `PHENO_CALLS` = data/development/typing/pln/pheno (it reads producer 4 for the compositional cross-modality comparison). Loads `core/phenocycler/scripts/pln/cluster_pln.py` by path (:561-563). | :90 `OUT_ROOT` = **data/development/xenium_typing/release_pln_immune_v3**; calls tree at :1287-1288. Class column **`lineage`** (:1128) with **`common_lineage`** (:1129-1131), `stage1`, `cluster_lineage` and `tnk_subcluster_call` beside it. Declares `supersedes: release_pln_immune_v2` (:1484-1486). `release.json` read 2026-09-22: `release_id` `xenium-pln-immune-tnk-split-dev-20260917` (RELEASE_ID at :103), `tissue` `pln`, `generated_at_utc` 2026-09-17T20:52:22Z. | build_layers.py:652 declares it for grid_composition; :548 declares the whole `data/development/xenium_typing` tree for distance_profiles; :615 declares it for compartments (until 2026-09-23 that line declared the **superseded** `release_pln_immune_v2`, which no builder reads); :179-181 names pln_immune_typing_v3.py. xenium_release.py:37-40 honours the supersedes declaration. Column readers: build_grid_composition.py:515; build_distance_profiles.py:491-500. Cards: card_n2_pln_immune.R:13, card_g1_grid.R:33, card_s3_compartments.R:21. |
| **pLN, PhenoCycler** | `core/phenocycler/scripts/typing/build_pln_immune_pheno.py` | Docstring :131-133: `--stage standout`, `--donors 6442`, or plain; arguments at :1519-1520. Interpreter pinned at docs/PLN_IMMUNE_PHENOCYCLER_2026-09-17.md:15: `.venv/bin/python core/phenocycler/scripts/typing/build_pln_immune_pheno.py`. | :163 `PLN_CELLS` = core/phenocycler/data/integration/cells_pln; :164 the ingest report; :165 `IDENTITY_CORE`; :166 `CLUSTERS` = data/development/typing/pln/clusters/phenocycler, written by `core/phenocycler/scripts/pln/cluster_pln.py` (cited at :22 as `scripts/pln/cluster_pln.py`). Imports `lineages` and `subpanel` (:156-157). | :168 `OUT` = **data/development/typing/pln/pheno**, `donor_id=<d>/roi=pln_1/data.parquet` (:1756-1758). Class columns **`immune_class`** (:1733), `immune_subtype` (:1734), `t_nk_class` (:1746, empty outside the T_NK compartment) and `is_immune` (:1732, nullable). The only file it deletes is its own stale `UNAVAILABLE.json` (:1759-1763). | build_layers.py declares `typing/pln/pheno` for distance_profiles (:547), compartments (:624) and grid_composition (:650); :174-175 names the producer. Column readers: build_compartment_gradient.py:655 (`is_immune`, `immune_class`), :739; build_grid_composition.py:554 (`is_immune`, `immune_class`, `t_nk_class`), :604 (`per_donor_status.csv`). Also read by producer 3 (pln_immune_typing_v3.py:97). Cards through layers: card_s3_compartments.R:21, card_g1_grid.R:33, card_n2_pln_immune.R:13. |

Two things the table settles. First, the interpreter for all four is the repo `.venv`; root
CLAUDE.md:11 says the `analysis/xenium` lane runs in `scvi-env`, and both Xenium write-ups
spell `.venv/bin/python` in their reproduce blocks (docs/PANC_XENIUM_TYPING_2026-09-17.md:650,
docs/PLN_IMMUNE_XENIUM_2026-09-17.md:471). The one `scvi-env` entry point in this inventory is
`xenium_typing/joint_model.py` (:7-8), which is not a producer. Second, until 2026-09-23 three
nodes reached the Xenium releases without declaring them -- `availability` and
`distance_profiles` named only `cells_xen`, and `compartments` named the superseded `_v2` -- so
the staleness check could not see a new release. Step 5 gave the first two the release roots
(`availability` has since dropped its two Xenium capabilities, which no card read, and with
them every release read and declaration -- Step 5c);
`compartments` now declares the v3 root it resolves (build_layers.py:615), on the maintainer's
go of the same day, and `test_no_node_declares_the_superseded_pln_release` keeps it so.

## Every other runnable typer

Ledger buckets are `scripts/retire/dispositions.py`'s: `CODE_TIER1-4`, `KEEP_BESIDE_CODE`
(:286-), `KEEP_BESIDE_CODE_TREES` (:281-284). `retire.py` refuses everything on
`KEEP_BESIDE_CODE`, `CODE_TIER3` and `CODE_TIER4` outright (dispositions.py:270-277).
"Not in the ledger" means the path is in none of those lists, measured 2026-09-21.

| path | class column and output | live uses, named | ledger disposition after this step |
|---|---|---|---|
| `core/phenocycler/scripts/celltype/run_step3.py` | `call` at data/development/celltype/calls (run_step3.py:25; panel.py:54). Coverage: donor 6450 only, older than its Step 2 inputs (docs/WORKFLOW_2026-09-02_MAP.md:77). | No reader of its `calls/` data anywhere in the tree. The **module** is imported by bare name for `build_scale` / `REQUIRE_MIN` / `NUCLEAR_RATIO` / `ALIGNED`: run_islet.py:45, run_sameslide_protein.py:73, figures_ambiguous.py:36, broadtype/arms/cosine_broad8.py:67, typing/refine.py:130. Its lane's tests run from Makefile:33 (`LANES`). | **LIBRARY** with those five importers named. Its `calls/` output root is a retirement candidate as data (one donor, older than its inputs, zero readers); the module stays. `KEEP_BESIDE_CODE` (dispositions.py:298, "changed on merge/hpg-sameslide"), so no move until that branch lands. |
| `core/phenocycler/scripts/celltype/run_islet.py` | `call` at data/development/celltype/islet_residents (:47) and `islet_residents_all.parquet` (:227). | Data readers: kintsugi/concordance.py:40 (`COSINE_ISLET`), kintsugi/run_step4.py:36 (runnable, :61-62), celltype/figures_ambiguous.py:171. Module importer: run_sameslide_protein.py:72. | **AUDIT** instrument with its reader named: the kintsugi method comparison and figures_ambiguous.py. `islet_residents_all.parquet` is a two-lane comparison artifact and is never a source of `cell_class` for any node. `KEEP_BESIDE_CODE` (dispositions.py:308-309). |
| `core/phenocycler/scripts/celltype/run_sameslide_protein.py` | `call` at data/development/celltype/sameslide_protein/donor_id=<d> (:82), carried with `reason`, `cell_region`, `structure`, `islet_endocrine`, `islet_num`. Donors: 6533 and 112/0112, normalised to one donor (:23-24). | sameslide/lift_calls.py:27-28 names the path in its usage block and :51 `CARRY` is exactly that column set, lifted onto Xenium cells. Imports run_step3 (:73), run_islet (:72). | **BRIDGE for donors 6533 and 0112**: it gives the same-slide Xenium cells a protein prior with its own uncertainty and reaches no build_layers node and no R card. **RULED KEPT 2026-09-23** with the whole `sameslide/` tree and `restain/run_field_refine.py`: "the 'one typing' goal does not cover sequential slide and same slide typing" (`CODE_LIBRARY`, dispositions.py). |
| `core/phenocycler/scripts/typing/refine.py` | `refined_class` (:688, :698) in a keyed sidecar at `--out`, required, no default (:781). Label `UNRESOLVED` = "Unresolved (no lineage evidence)" (:150). | Joined into `typing/cells` by `build_cells.py --refine` (the producer-1 row above); no builder reads `refined_class` -- they read `cell_class`, and whether they should read the refined class instead is open (see "Not decided here"). The module is also live as the **vocabulary authority**: vocab.py:206-209 (`PHENO_UNRESOLVED`), build_grid_composition.py:255-272 (`PANC_CLASS_NONCALL`), modality_concordance.py:186, run20_default.yml:41-42, 00_globals.R:652-653; and as the import-time guard on the producerless `refinement_log.json` (:19-25). | **Part of producer 1** since the maintainer's ruling of 2026-09-23 ("needs to be included in the cell typing for phenocycler data"), through the `--refine` join only -- not a fifth producer, and it touches no route of the islet chain. Still the vocabulary authority, and a declared library in the ledger (`scripts/retire/dispositions.py`). Its 5th-percentile bound, which demotes about 5% of every class by construction (REFINEMENT.md section 6), is unchanged. |
| `core/phenocycler/scripts/typing/call.py` | The named-cluster label at data/development/typing/calls (:53, :497-501). Header :4-10: "RETIRED -- NOT A FALLBACK ... It did not produce the calls on disk." | One importer: typing/compare_clustering.py:105 (`from call import TOP_N, _statistic, name_clusters`). Its output tree is already retired as data (dispositions.py:191-193). | **QUARANTINED 2026-09-23** with compare_clustering.py, on the ruling "Retire, but if data was ever produced with this, I want to see it regardless." What it produced is kept and was shown: `typing/calls` for donors 6436 and 6442 (2026-09-11; data tier 2 since 2026-09-16) and `typing/clustering_comparison.json` (2026-09-08, donor 6442, recommendation leiden). The 20 lane tests of its naming statistic left with it. |
| `core/phenocycler/scripts/broadtype/` (whole tree, incl. `arms/cosine_broad8.py`, `release.py`) | Methods write `label` under data/development/broadtype/arms (records.py:29, :53; pins.py:40). `release.py` turns them into data/development/broadtype/release with `lineage_native` + `common_lineage` (release.py:93-95). Run as cosine_broad8.py:44. | (a) core/phenocycler/config.ini:152 bound `broadtype/release` as the integration lane's protein identity release until 2026-09-23, attached through export_pheno.py:330 and :381-386 onto `cells_pheno`. (b) The report lane: identity_io.py:17-19, pheno_identity.py:62, :111, :176, rendered by run_report.py:77-78. (c) kronos2_doc/run.py:35-40, rendered by run_report.py:83-84. (d) subpanel_io.py:73 and :186-188 read `broadtype/scratch/cosine_scale.parquet` as producer 1's cohort floor/top; superseded.py:52-57 and :390-393 exempt that directory; build_layers.py:448 declares the parquet for app_drilldown_cells and :184-186 names it. | Retired as a producer (root CLAUDE.md:198-200, the comparison withdrawn 2026-09-08); declared **LIBRARY + AUDIT** with (b), (c) and (d) named. `CODE_TIER3`, overturned (dispositions.py:552-564); its release, frame and scratch directories are data `TIER3` (:210-215). The config.ini:152 binding was blanked 2026-09-23 (Step 5) -- see the closed-leak section. |
| `core/phenocycler/phenocycler/simultaneous_typing.py` | No output root: `type_simultaneous_cells` returns a frame and the `type` stage that would write it is blocked by `NO_VERIFIED_TYPING` (:3-7). | Two importers inside the installed package: pipeline.py:80-83, production_typing.py:21-26. | **LIBRARY behind the latch**, with both importers named and `NO_VERIFIED_TYPING` named as what keeps it from producing. Retiring it would break the package's own `pipeline` import. Not in the ledger. |
| `core/phenocycler/scripts/pln/cluster_pln.py` | Writes `cluster_id` (:595, :610, :784) and `meta_cluster` (:1016, :1041) under data/development/typing/pln/clusters (:173). **No class column** -- a partition, not a typer. Run as :112-113. | Both pLN producers depend on it: build_pln_immune_pheno.py:166 reads its output tree (:22 "consumed here"); pln_immune_typing_v3.py:561-563 loads the file by path. A by-path `importlib` load is invisible to a grep for import statements. | **LIBRARY** upstream of producers 3 and 4, both named. Not in the ledger. |
| `analysis/xenium/scripts/hierarchical_phenotype.py` | `celltype`, `celltype_broad`, `celltype_fine` and confidences into `data/processed/<s>/<s>_hierarchical.h5ad` or in place (:47-52). Its docstring :38-40 describes a split of the exocrine label that root CLAUDE.md:201-202 forbids; :82-84 records that the current pancreas release carries a single `Exocrine` lineage. | analysis/xenium/hpc/run_pipeline.sh:163-164 runs it as the HiPerGator `phenotype` stage. Cross-lane: xenium_typing/nullgate_typer.py:43-44 reads its sidecar's `celltype_broad` / `celltype` as the scANVI method's `lineage_native` / `celltype_fine`. | **RETIRED-CANDIDATE with a prerequisite**: producer 2 replaces it for the cohort, but deleting it makes nullgate_typer.py's input producerless. Close the `xenium_typing` chain below first (blank config.ini:153). The HiPerGator lane's own use stays named. Not in the ledger. |
| `core/phenocycler/scripts/xenium_typing/{release,nullgate_typer,joint_model,baseline_analyzer}.py` (one entry) | Methods write `lineage_native` + `celltype_fine` under data/development/xenium_typing/arms (gate.py:51; pins.py:60, :208-209). release.py:175 writes **data/development/xenium_typing/release** with `lineage_native` + `common_lineage` (:149); `release.json` read 2026-09-22: `release_id` `xenium-typing-scanvi-per-donor-dev-20260907`, `written_utc` 2026-09-07. joint_model.py:7-8 is the one `scvi-env` entry point (:62 output, :289 columns); baseline_analyzer.py:68 writes `Unavailable`. | **Was the load-bearing one, until 2026-09-23**: core/phenocycler/config.ini:153 bound this release as `xenium_identity_release`; import_xenium.py:454-456 loaded it and :427-429 attached it; phenospace.py:532-533 wrote `lineage_native` and `lineage_common` onto `cells_xen`; then the availability and distance_profiles builders read those columns as a fallback -- see the closed-leak section. The binding is blank and both readers are gone. Report lane: xen_identity.py:157-160 reads `joint_model/train_latent.npy`. panc_xenium_release.py:86-96 declares this release superseded. | **RETIRED-CANDIDATE, prerequisite met 2026-09-23**, and **QUARANTINED the same day** on the maintainer's ruling: release.py, nullgate_typer.py, scanvi_replicates.py and baseline_analyzer.py are under `data/ARCHIVE_2026-09-15/quarantine/` (the two lane tests that imported them left with them). joint_model.py stays while xen_identity.py:157-160 reads its latent. hierarchical_phenotype.py is still held by its HPC lane use (`analysis/xenium/hpc/run_pipeline.sh:163`). |

`typing/compartments.py`, `common/window_mask.py`, `common/label_space.py` and the rest of the
"not on any tier" set in the reference walk write no class column and are not in this table.

## Where a second method reached the app, and how it was closed

Measured 2026-09-21, re-read 2026-09-22, closed 2026-09-23 (Step 5). The record is kept so the
signature is recognisable if it comes back.

1. **The bindings.** core/phenocycler/config.ini:151 says "Optional declared identity
   releases. Blank is the fail-closed state." Until 2026-09-23 the next two lines were not
   blank: :152 `phenocycler_identity_release` bound `../../data/development/broadtype/release`
   (the withdrawn broad-typing comparison, retired as a producer at root CLAUDE.md:198-200)
   and :153 `xenium_identity_release` bound `../../data/development/xenium_typing/release`
   (the scANVI two-donor release that producer 2 declares superseded,
   panc_xenium_release.py:86-96). **Both are blank now**, in
   `config.ini` and in `config_metabolic.ini`, and `test_both_identity_release_bindings_are_blank`
   in `core/phenocycler/tests/current/test_declared_producers.py` keeps them so.
2. **How :153 reached the cells.** import_xenium.py:454-456 calls
   `phenospace.load_release_or_none(cfg, "xenium")`, which returns None on a blank key, and
   :427-429 `attach_release_identity`; phenospace.py:532-533 wrote `lineage_native` and
   `lineage_common` onto `core/phenocycler/data/integration/cells_xen`. The tree on disk still
   carries those columns from its last ingest; nothing reads them, and the next ingest writes
   none.
3. **Which nodes read it, and what they read now.** `availability` read `lineage_native` off
   cells_xen through `xenium_ingest_availability`, called only when a tissue's release did not
   resolve. That function is gone: a tissue with no resolvable release is `unavailable` with
   the resolver's own reason and the roots it searched, on every donor. `distance_profiles`
   read `lineage_common` off cells_xen through `_pln_arm`, the pLN fallback when no typed pLN
   release resolved. That function is gone: an unresolved pLN release records
   `xenium_pln_release_absent` in the artifact's provenance and emits no pLN Xenium row. Both
   nodes now declare the two declared release roots as inputs in build_layers.py, so a new
   release makes them stale; both still declare cells_xen, for the IF pairing and for a
   pancreas cell's coordinates, which only the ingest carries. Later the same day (Step 5c)
   `availability` dropped `xenium_panc` and `xenium_pln` altogether -- no card read either --
   so it now reads no release and declares none; `distance_profiles` still does both.
4. **Which card rendered it.** apps/islet/R/card_x1_gradient.R:143 requires `immune_subtype`
   and `resolved_by` from `distance_profiles/data.parquet` and :271 pastes them into the pair
   label. `resolved_by` now takes only the release arms' values (`all`, `identity_core`,
   `stage1_unresolved`, `common_lineage`); `lineage_common` is not spelled by any builder.
5. **The one reader outside the graph.** scripts/senior/same_section_linkage.py's
   `agreement()` read `lineage_native` and `celltype_fine` off cells_xen for the per-cell
   protein-vs-RNA figure (docs/SAME_SECTION_LINKAGE_2026-09-17.md section 4). It now reads the
   declared pancreas release joined on `cell_id`; the `endocrine_fine` strata are no longer
   produced because the declared release carries no fine endocrine type. Its recorded numbers
   are from the old source until the script is re-run, and the document says so.
6. **What did not leak.** `cells_pheno` (the :152 twin) is a declared input of no
   build_layers node; it reaches only the report lane (identity_io.py:19; pheno_identity.py:62,
   :111, :176). `restain` declares `cells_xen` but build_restain_layer.py names no lineage
   column. `app_drilldown_cells` declares `cosine_scale.parquet` as a scale, not a class.

The guard `test_class_columns_reach_the_graph_only_from_declared_producer_roots` is an
ordinary assertion now (its strict xfail came off with the fix). Its signature is unchanged: a
node declares a root **and** its builder names one of that root's class columns in code
(docstrings excluded). The ledger's four `xenium_typing` entries held on this prerequisite now
pass the binding half of their check (`scripts/retire/dispositions.py`, `XTYPING_PREREQ`); the
files have not been moved.

## Not decided here

**Status 2026-09-24.** Steps 0, 5, 5a, 5b and 5c of the resumption plan are landed: one writer
per app artifact, both identity-release bindings blank, both ingest fallbacks removed, the release
roots declared as inputs of `distance_profiles` (and the v3 pLN root of `compartments`), the
strict xfail off, the same-section agreement re-pointed, the ledger
executed (25 files quarantined -- six on the day's rulings; the other nineteen, the
`xenium_celltype` detection trio and its `hpg/` apparatus, were moved in 41c71ac BEFORE any
ruling and ruled "Yes to all 19" afterwards, once listed; the same-slide lane declared kept; no
candidate remains), the refinement wired into producer 1 through `build_cells.py --refine`,
and the graph remainder (5c): `--check` names every stale input and an unfinished run, reports
`typing/immune` against what it is built from, `islet_summary` stops writing its unread outline
column, and `availability` its four unread capabilities. Separately, 6579 joined both cohort
manifests with a real image entry (`c259ca8`). Nothing under `data/development` or
`/data/islet_layers` was written; every served artifact is still the pre-fix build until the
DATA step runs on a go. Step 7's DATA step then wrote `typing/refine`, re-joined `typing/cells`
and rewrote `typing/immune` (2026-09-24); its layer rebuild stopped at `islet_geometry` on
6374's missing Nerve export, so `/data/islet_layers` and `typing/app` are unchanged and still
pre-join. **Step 7 finished 2026-09-26**, once 6374's Nerve / Blood vessel / Lymphatic re-exports
were placed: 6374 was retyped (crossed fraction, `subpanel.py`, build, `refine.py` over all 20, the
join, `immune_subtype.py`; 2,477 of its cells changed `cell_class`, the other 19 donors identical
row for row), then `GO=1 make rebuild` ran all 23 stale nodes with exit 0 and `make check` reports
24 fresh. Every builder still reads `cell_class`. Step 6 follow-up debt: seven docs still say the join has not run --
`apps/islet/CLAUDE.md`, `apps/islet/README.md`, `apps/islet/docs/capability_status.md`,
`core/phenocycler/README.md`, `core/phenocycler/scripts/typing/RECONSTRUCTION.md`,
`core/phenocycler/scripts/typing/REFINEMENT.md` and `scripts/CLAUDE.md`.

**Ruled 2026-09-23, not open any more:** `refine.py` IS to be included in the PhenoCycler
typing -- scoped the same day to the `build_cells.py --refine` join, beside `cell_class`, with
the islet chain untouched (the join ran on disk 2026-09-24; builders still read `cell_class`);
`call.py` and `compare_clustering.py` retired; the four `xenium_typing` entries quarantined;
the nineteen other candidates quarantined ("Yes to all 19"); the same-slide lane kept;
`compartments` declares the v3 pLN release.

Open for the maintainer:

- **Which class the builders serve.** After the join, `typing/cells` carries `refined_class`
  beside `cell_class`, and every builder still reads `cell_class`. Serving the refined class
  instead needs decisions that touch the islet chain -- what a demoted endocrine cell's
  `is_endocrine` becomes, and what a demoted cell inside an islet outline that the pass routes
  to Exocrine becomes -- so it is not part of the wiring and nothing here decides it. Two
  answers were given on 2026-09-23 while a different design (the pass inside `build_cells.py`,
  serving the refined class) was on the table, and are recorded, not applied: "is_endocrine can
  be Alpha/Beta/Delta or CD99. CD99 is an endocrine cell that we don't have the specific marker
  for. If the CD99 tests indicate another cell type it does not keep the is_endocrine flag";
  and, for the in-islet case, `Excluded` with its own reason.
- **Two counts of islet size.** `refine.py` fits its bounds per size group counted over the
  signed-distance core (`cell_region_20um == "core"`, as `islet_identity.load` counts);
  the islet chain counts over outline membership (`in_islet`). Measured on a scratch join of
  donor 6539 on 2026-09-23: they disagree about the size group for 643 of 364,201 cells (0.18 %)
  and about the islet size for 8,723. The join keeps the pass's own as `refine_size_stratum` /
  `refine_islet_cells` beside the chain's, by design. Recorded, not decided.

- `analysis/xenium/scripts/hierarchical_phenotype.py`: tier 2 with a live HPC-lane use
  (`analysis/xenium/hpc/run_pipeline.sh:163`); the ledger refuses it until that lane's owner
  decides.
