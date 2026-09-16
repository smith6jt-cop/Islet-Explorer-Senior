"""Declared dispositions for the 2026-09-15 data retirement.

Every entry carries a CHECK: a named predicate re-verified at manifest time. Nothing is retired on
the strength of this table alone -- `build_manifest.py` runs the check and records the result, and
`retire.py` refuses to act on any entry whose check did not pass.

Tiers
-----
1  confirmed dead -- duplicate or strict subset, nothing reads it. Deleted.
2  superseded but the replacement is NARROWER. Carve-out first, then quarantined (moved), not deleted.
3  proposed for retirement and OVERTURNED by the still-read / narrower adversarial passes. Untouched.
4  maintainer's call. Untouched.
"""
from __future__ import annotations

REPO = "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior"

# --- tier 1 -----------------------------------------------------------------
# check kinds:
#   orphan_temp(evidence)      : path exists; `evidence` also exists (the completed product)
#   subset_dir(other)          : every relative file path under `path` exists under `other`
#   identical_to(other)        : md5(path) == md5(other)
#   names_exist_under(other)   : every basename directly under `path` exists under `other`
#   donor_local_cover(pattern) : every donor_id in the json has {pattern} on disk
#   plain                      : existence only (logs, backups named as such)
TIER1 = [
    dict(path="/data/tmp_b2r", check="orphan_temp",
         arg="/data/islet_ome_tiff/6593.ome.tiff",
         superseded_by="/data/islet_ome_tiff/6593.ome.tiff",
         evidence="Single orphaned mktemp dir 6593.0zsFPo; the conversion completed -- "
                  "6593.ome.tiff is 89.2 GB with its offsets sidecar, both 2026-09-14 11:06."),
    dict(path=f"{REPO}/data/development/restain/valis_baseline", check="subset_dir",
         arg=f"{REPO}/data/development/restain/valis",
         superseded_by=f"{REPO}/data/development/restain/valis",
         evidence="Same 2 donors; valis/ adds if_nuclei_um.npy, if_nuclei_z.npy, "
                  "refined_displacement.npz, refined.json, xenium_nuclei_cell_ids.npy, "
                  "xenium_nuclei_um.npy. Strict subset."),
    dict(path="/data/tissuumaps/channels/6539/.tissuumaps", check="names_exist_under",
         arg="/data/tissuumaps/channels/6539",
         superseded_by="/data/tissuumaps/channels/6539/*.tif",
         evidence="Hidden cache of 4 tifs (00_DAPI, 13_INS, 35_SST, 46_GCG); all four names exist "
                  "among the 59 deliberate marker tifs beside it."),
    dict(path=f"{REPO}/data/Islets/pre-rerun-2026-07-09", check="names_exist_under",
         arg=f"{REPO}/data/Islets",
         superseded_by=f"{REPO}/data/Islets/",
         evidence="All five re-exported into the parent 2026-09-14/15 at 3.2-4.0x the size."),
] + [
    dict(path=f"{REPO}/core/phenocycler/data/integration/{d}", check="subset_dir",
         arg=f"{REPO}/core/phenocycler/data/integration/registration",
         superseded_by=f"{REPO}/core/phenocycler/data/integration/registration",
         evidence="Every donor present in live registration/ (20 donors); backup content predates it.")
    for d in ("registration.backup_before_rotation_20260909",
              "registration.backup_before_coarse_widen_20260910",
              "registration.backup_pre_final_20260910",
              "registration.backup_pre_prewarp_20260910",
              "registration.backup_pre_rotation_bound_20260911")
] + [
    dict(path=f"{REPO}/data/development/restain/figures/donor_id={d}/valis", check="names_exist_under",
         arg=f"{REPO}/data/development/restain/figures/donor_id={d}/refined",
         superseded_by=f"{REPO}/data/development/restain/figures/donor_id={d}/refined",
         evidence="alignment_field.png and alignment_windows.png both exist under refined/, one "
                  "generation later on the transform that shipped.")
    for d in ("6436", "6442")
] + [
    dict(path=f"{REPO}/data/development/restain/alignment_panc_valis_baseline.json",
         check="identical_to",
         arg=f"{REPO}/data/development/restain/alignment_panc_valis.json",
         superseded_by=f"{REPO}/data/development/restain/alignment_panc_valis.json",
         evidence="Byte-identical twin (md5 cd4d9638fe882a3dee73ae126c4af82a). Keep the non-_baseline copy."),
    dict(path=f"{REPO}/data/development/restain/alignment_panc_refined.json",
         check="donor_local_cover", arg="alignment_refined.json",
         superseded_by="restain/valis/donor_id=*/alignment_refined.json",
         evidence="Single record, donor 115; valis/donor_id=115/alignment_refined.json exists."),
    dict(path=f"{REPO}/data/development/restain/alignment_panc_image_seed.json",
         check="donor_local_cover", arg="alignment_image_seed.json",
         superseded_by="restain/valis/donor_id=*/alignment_image_seed.json",
         evidence="Single record, donor 112; valis/donor_id=112/alignment_image_seed.json exists."),
    dict(path=f"{REPO}/data/development/restain/refine_panc.json",
         check="donor_local_cover", arg="refinement.json",
         superseded_by="restain/valis/donor_id=*/refinement.json",
         evidence="Single record, donor 115; valis/donor_id=115/refinement.json exists. NOT to be "
                  "confused with refine_panc_0.5um.json, which covers 6436/6442 and is kept."),
    dict(path=f"{REPO}/data/development/restain/nuclei_panc.json",
         check="donor_local_cover", arg="nuclei.json",
         superseded_by="restain/valis/donor_id=*/nuclei.json",
         evidence="4 records (6461 6476 6533 6591), all four donor dirs hold nuclei.json. "
                  "nuclei_panc_0.5um.json covers 6436/6442, which have none -- kept."),
    dict(path=f"{REPO}/data/development/restain/hormone_panc.json", check="identical_keys_subset",
         arg=f"{REPO}/data/development/restain/hormone_panc_refined.json",
         superseded_by=f"{REPO}/data/development/restain/hormone_panc_refined.json",
         evidence="Same 2 donors; the refined file's hormone block is a strict key superset "
                  "(adds pixel_size_um_recorded, preregistered_auroc_bar, transform_source)."),
    dict(path=f"{REPO}/data/development/restain/explorer_bundle/bundles_panc.json", check="plain",
         superseded_by="restain/explorer_bundle/donor_*/manifest.json",
         evidence="2 records; both per-bundle manifest.json exist and are supersets (add "
                  "islet_annotations). Stale 2026-09-02 index."),
    dict(path=f"{REPO}/data/development/restain/protein_rna/figures.json", check="plain",
         superseded_by=f"{REPO}/data/development/restain/protein_rna/figures_ins.json",
         evidence="figures_ins.json covers the same 2 donors and is a strict key superset."),
    dict(path=f"{REPO}/data/development/restain/protein_rna/protein_rna_6436_0059907.png", check="plain",
         superseded_by="protein_rna_6436_0059907_ins.png", evidence="Superseded figure generation."),
    dict(path=f"{REPO}/data/development/restain/protein_rna/protein_rna_6442_0059947.png", check="plain",
         superseded_by="protein_rna_6442_0059947_ins.png", evidence="Superseded figure generation."),
    # The ECC generation at the donor root (2026-08-31 16:31/16:32); the same two filenames exist
    # under refined/ on the transform that shipped.
    dict(path=f"{REPO}/data/development/restain/figures/donor_id=6436/alignment_field.png",
         check="sibling_name_under", arg=f"{REPO}/data/development/restain/figures/donor_id=6436/refined",
         superseded_by="figures/donor_id=6436/refined/alignment_field.png",
         evidence="ECC generation; the same filename exists under refined/."),
    dict(path=f"{REPO}/data/development/restain/figures/donor_id=6436/alignment_windows.png",
         check="sibling_name_under", arg=f"{REPO}/data/development/restain/figures/donor_id=6436/refined",
         superseded_by="figures/donor_id=6436/refined/alignment_windows.png",
         evidence="ECC generation; the same filename exists under refined/."),
    dict(path=f"{REPO}/data/development/restain/figures/donor_id=6442/alignment_field.png",
         check="sibling_name_under", arg=f"{REPO}/data/development/restain/figures/donor_id=6442/refined",
         superseded_by="figures/donor_id=6442/refined/alignment_field.png",
         evidence="ECC generation; the same filename exists under refined/."),
    dict(path=f"{REPO}/data/development/restain/figures/donor_id=6442/alignment_windows.png",
         check="sibling_name_under", arg=f"{REPO}/data/development/restain/figures/donor_id=6442/refined",
         superseded_by="figures/donor_id=6442/refined/alignment_windows.png",
         evidence="ECC generation; the same filename exists under refined/."),
    dict(path="/data/IO60pLN/project.qpproj.backup", check="plain",
         superseded_by="/data/IO60pLN/project.qpproj", evidence="Named backup of the live project."),
    dict(path=f"{REPO}/data/gate4_apply.log", check="plain", superseded_by="-",
         evidence="Log of the retired Gate-4 run."),
    dict(path=f"{REPO}/data/gate4_lineage.log", check="plain", superseded_by="-",
         evidence="Log of the retired Gate-4 run."),
    dict(path=f"{REPO}/data/redsea_v10_run.log", check="plain", superseded_by="-",
         evidence="Log of the run that produced the superseded data/cells_redsea."),
]

# --- tier 1 traps: verified present, never retired ---------------------------
KEEP_BESIDE_TIER1 = [
    (f"{REPO}/data/development/restain/alignment_panc_valis.json",
     "the only on-disk VALIS-arm certification for 6436/6442; neither donor dir has alignment_refined.json"),
    (f"{REPO}/data/development/restain/nuclei_panc_0.5um.json",
     "covers 6436/6442, which have NO donor-local nuclei.json"),
    (f"{REPO}/data/development/restain/refine_panc_0.5um.json",
     "covers 6436/6442, which have NO donor-local refinement.json"),
    (f"{REPO}/data/development/restain/alignment_panc.json",
     "ECC-arm record for 6436/6442; no donor-local successor exists"),
    (f"{REPO}/data/development/restain/registration_panc.json",
     "the only record of the endocrine validation; validate.py has not been re-run"),
]

# --- tier 2: carve out, then quarantine -------------------------------------
TIER2 = [
    dict(path=f"{REPO}/data/cells", carve="donors_and_columns",
         superseded_by=f"{REPO}/core/phenocycler/data/runs/30f6689e4d27d295f890/{{cells,geometry_qc}}",
         evidence="donor_id=6457 (728,698 rows) and 6579 (2,041,786) exist here and in no replacement; "
                  "14 columns absent from runs/.../cells."),
    dict(path=f"{REPO}/data/cells_redsea", carve="columns",
         superseded_by=f"{REPO}/core/phenocycler/data/runs/30f6689e4d27d295f890/{{expression,redsea}}",
         evidence="cell_area_px survives only in runs/.../redsea; v10 values are not byte-reproducible "
                  "from the run (INS corr 0.9979, max abs diff 1661.6)."),
    dict(path=f"{REPO}/data/cells_preqc", carve="columns",
         superseded_by=f"{REPO}/core/phenocycler/data/runs/30f6689e4d27d295f890/cells",
         evidence="Rows and donors match exactly; the replacement drops the 14-column QuPath "
                  "morphometry block and the two signed distances."),
    dict(path=f"{REPO}/data/development/typing/positivity", carve="none",
         superseded_by=f"{REPO}/data/development/typing/cells",
         evidence="Withdrawn two-donor lane. 39 per-marker pixel-positivity columns have no "
                  "counterpart in cells/; retired by ruling, not replaced."),
    dict(path=f"{REPO}/data/development/typing/qupath_export", carve="none",
         superseded_by=f"{REPO}/data/development/typing/cells",
         evidence="typed-orbit retires export_qupath.py. 11 of 13 columns (confidence, "
                  "ranked_*_probabilities) have no counterpart; retired by ruling."),
    dict(path=f"{REPO}/data/development/typing/substrate", carve="none",
         superseded_by="regenerable from celltype/uniform/ + broadtype/scratch/cosine_scale.parquet",
         evidence="104 v__/r__ value+rank columns for 52 markers. Regenerable but NOT byte-recoverable "
                  "(measured corr 0.945 against uniform/, max abs diff 22,511)."),
    dict(path=f"{REPO}/data/cell_qc", carve="donors",
         superseded_by=f"{REPO}/core/phenocycler/data/runs/30f6689e4d27d295f890/geometry_qc",
         evidence="22 -> 20 donors; 9 of 11 columns have no like-named successor. qc_pass / "
                  "qc_exclude_reason are the retired gate's verdict, not re-derived."),
    dict(path=f"{REPO}/data/development/typing/scratch/subpanel_context", carve="columns",
         superseded_by=f"{REPO}/data/development/typing/{{subpanel_rebuilt,cells}}",
         evidence="in_Nerve and in_Lymphatic are not reproducible: reconstruction from cells/ loses "
                  "18.4% of in_Nerve and 325 cells of in_Lymphatic."),
    dict(path=f"{REPO}/data/development/typing/astir", carve="none",
         superseded_by=f"{REPO}/data/development/typing/cells",
         evidence="Withdrawn two-donor bake-off arm; per-class posteriors have no counterpart."),
    dict(path=f"{REPO}/data/development/typing/tacit", carve="none",
         superseded_by=f"{REPO}/data/development/typing/cells",
         evidence="Withdrawn two-donor bake-off arm; per-class posteriors have no counterpart."),
    dict(path=f"{REPO}/data/development/typing/calls", carve="none",
         superseded_by=f"{REPO}/data/development/typing/cells",
         evidence="Withdrawn two-donor bake-off arm."),
] + [
    dict(path=f"{REPO}/data/development/celltype/{n}", carve="none",
         superseded_by=f"{REPO}/data/development/celltype/islet_residents_all.parquet",
         evidence="Three same-day variants of one experiment; identical 10-col / 1,351,845-row schema. "
                  "The survivor is islet_residents_all.parquet, NOT typing/cells/.")
    for n in ("islet_residents_v2", "islet_residents_v3", "islet_residents_leaf",
              "islet_residents_all_v2.parquet", "islet_residents_all_v3.parquet",
              "islet_residents_all_leaf.parquet")
]

# --- tier 3: overturned, untouched -------------------------------------------
TIER3 = [
    (f"{REPO}/data/development/typing/subpanel/subpanel_calls.parquet",
     "typed-orbit Phase 0 names it as the row-for-row reconstruction gate; Phase 0 cannot run without it"),
    (f"{REPO}/data/development/typing/app/json",
     "region source of the ten built TissUUmaps projects; only in-repo source of IsletNonEndocrine outlines"),
    (f"{REPO}/data/development/broadtype/release",
     "bound in config.ini and read as data by the report lane -- 13 rendered figures"),
    (f"{REPO}/data/development/broadtype/frame",
     "two files are declared inputs in the --verify manifest; arm_base.py fails closed without them"),
    (f"{REPO}/data/development/broadtype/scratch",
     "five live modules name it as their cache root; arm_base.py reads/writes an object-id cache inside it"),
    (f"{REPO}/data/development/report/_pair/bakeoff",
     "registered OUTPUT records in the --verify manifest; removal makes run_report.py --verify fail"),
    (f"{REPO}/data/development/typing/manifests/qupath_export.json",
     "the QuPath importer resolves <export_dir>/../manifests/qupath_export.json"),
    (f"{REPO}/data/restore_gated_redsea",
     "app does not read it, but the CLAUDE.md carve-out is live standing text -- amend the text first"),
    (f"{REPO}/data/phenotype/broad",
     "same live CLAUDE.md carve-out; also the only artifact carrying 25 <marker>_state columns"),
    (f"{REPO}/data/phenotype/qupath_class",
     "a live doc names it the surviving 20-donor adjudication instrument; replacement covers 2 of 20"),
    ("/data/xenium/h5ad_out",
     "xenio.load_donor defaults to prefer='h5ad'; also the only 24 pLN h5ads"),
    (f"{REPO}/data/redsea_scratch/intermediates",
     "read by boundary_contact_audit.py --intermediates; only per-cell polygons in this checkout"),
]

# --- tier 4: maintainer's call ------------------------------------------------
TIER4 = [
    ("/data/PancLN_Xenium", "zero repo references, no coordinates, no cell types, includes excluded "
                            "6457/6579 and a fourth disease status. AA_female*.h5ad is a separate call."),
    ("/data/IO60pLN", "RAW acquisition plus its live QuPath project -- increment 9's input, not dead"),
    ("/data/QP_Xenium", "two repo statements conflict and neither settles it"),
    ("/data/islet_ome_tiff/6457.ome.tiff", "centrally excluded donor"),
    ("/data/islet_ome_tiff/6579.ome.tiff", "centrally excluded donor"),
]
