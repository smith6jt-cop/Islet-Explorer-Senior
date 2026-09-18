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

Two KINDS of entry live here. Everything down to TIER4 is DATA. Everything below the CODE TIER
banner is CODE (and two artifacts that have no entrypoint on either side), added 2026-09-18; its
own tier meanings are written at that banner. `build_manifest.py` writes a `kind` column so the two
never get counted together, and `retire.py --kind` selects one.
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


# ============================================================================
# CODE TIER -- added 2026-09-18 (Lane U)
# ============================================================================
# Root CLAUDE.md: "When replacing a workflow, remove the obsolete entrypoint rather than
# leaving two runnable ones." The tiers above cover DATA. These cover CODE, with the same
# fail-closed contract: `build_manifest.py` runs each entry's CHECK now and `retire.py`
# refuses anything that did not PASS.
#
# Code tiers
# ----------
# c1  verified orphan. No import, source() or subprocess call reaches it from outside the
#     retiring set, AND it names no output that still exists under a live domain root.
#     Quarantined (moved). Deletion is a second, reviewable commit.
# c2  verified orphan, but a NAMED PRODUCER of an artifact that still exists under a live
#     domain root with no current replacement producer. Quarantine ONLY -- never deleted by
#     this lane, because deleting it makes a live artifact producerless, which is the exact
#     failure `typing/umap/refinement_log.json` is already in.
# c3  proposed and OVERTURNED: a live importer outside the retiring set still reaches it.
#     Untouched. The importer is named in `evidence`.
# c4  maintainer's call. Untouched.
#
# check kinds (code)
#   code_orphan          : no import / source() / subprocess reference from outside
#                          CODE_RETIRING (the union of c1 and c2) and outside the file itself
#   no_reference_anywhere: the basename appears nowhere in tracked code, config or docs
#   protected            : the path exists and is on KEEP_BESIDE_CODE (asserted, never acted on)

# --- the HiPerGator exclusion, encoded so a future run cannot retire it by accident -------
#
# `merge/hpg-sameslide` is the branch the VS Code tunnel session is developing on RIGHT NOW
# (approved plan, Context section 5). Measured 2026-09-18 with
#   git diff --name-only stage0/groundwork...merge/hpg-sameslide
# which named 21 files, plus the modules those files import (read out of the branch's own
# copies), plus the three trees the plan names wholesale. `retire.py` refuses every path on
# this list outright, the same way it refuses TIER3 and TIER4.
HPG_BRANCH = "merge/hpg-sameslide"

#: Whole trees the plan takes off the list. A candidate under any of these is refused.
KEEP_BESIDE_CODE_TREES = [
    f"{REPO}/core/phenocycler/scripts/sameslide",
    f"{REPO}/core/phenocycler/scripts/xenium_celltype",
]

KEEP_BESIDE_CODE = [
    # -- named wholesale by the plan --------------------------------------------------
    (f"{REPO}/core/phenocycler/scripts/sameslide",
     f"{HPG_BRANCH} develops it (lift_calls.py changed on the branch); whole tree off the list"),
    (f"{REPO}/core/phenocycler/scripts/xenium_celltype",
     f"{HPG_BRANCH} develops context.py, eval_detection.py, hpg/verify_run.py and "
     "tests/test_verify_run.py; whole tree including hpg/ off the list"),
    (f"{REPO}/core/phenocycler/scripts/restain/run_field_refine.py",
     f"{HPG_BRANCH} adds it (398 lines)"),
    # -- the rest of the branch's own diff ---------------------------------------------
    (f"{REPO}/core/phenocycler/scripts/celltype/panel.py", f"changed on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/run_step1.py", f"changed on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/run_step3.py", f"changed on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/run_sameslide_protein.py",
     f"added by {HPG_BRANCH} (601 lines)"),
    (f"{REPO}/analysis/xenium/scripts/tissue_markers/pancreas.py", f"changed on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/tests/current/test_no_duplicate_tracked_data.py",
     f"added by {HPG_BRANCH}"),
    (f"{REPO}/docs/XENIUM_PIPELINE.md", f"changed on {HPG_BRANCH}"),
    # -- what those files import, read out of the branch's own copies -------------------
    (f"{REPO}/core/phenocycler/scripts/celltype/rules.py",
     f"imported by run_sameslide_protein.py and eval_detection.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/run_islet.py",
     f"imported by run_sameslide_protein.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/score.py",
     f"imported by run_sameslide_protein.py and run_step3.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/celltype/usability.py",
     f"imported by run_sameslide_protein.py and run_step1.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/restain/align.py",
     f"imported by run_field_refine.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/restain/registration_inputs.py",
     f"imported by run_field_refine.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/restain/run_align.py",
     f"imported by run_field_refine.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/restain/run_intensity_refine.py",
     f"imported by run_field_refine.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/scripts/restain/run_register.py",
     f"imported by run_field_refine.py on {HPG_BRANCH}"),
    (f"{REPO}/core/phenocycler/phenocycler/integration/transform.py",
     f"imported by lift_calls.py and run_field_refine.py on {HPG_BRANCH}"),
]


def _c(path, tier, produces=(), superseded_by="-", evidence="", check="code_orphan", prereq="",
       blocked_by=()):
    """One code disposition.

    `blocked_by` names entries that must move in the SAME invocation. The direction is
    importer-first: if A imports B, moving B alone breaks A, so B is blocked_by A. `retire.py`
    refuses an entry whose blocker is not also being acted on, which is what keeps a
    within-set dependency from becoming a broken import halfway through.
    """
    return dict(path=path, tier=tier, produces=tuple(produces), superseded_by=superseded_by,
                evidence=evidence, check=check, prereq=prereq, blocked_by=tuple(blocked_by))


TYPING = f"{REPO}/core/phenocycler/scripts/typing"
XSCRIPTS = f"{REPO}/analysis/xenium/scripts"
SENIOR = f"{REPO}/scripts/senior"
ROOTS = f"{REPO}/scripts"

#: The two live domain roots (apps/islet/config/{islet_default,run20_default}.yml). An entry that
#: names an output which still EXISTS under one of these is c2, not c1.
LIVE_DOMAIN_ROOTS = (f"{REPO}/data/app_data", f"{REPO}/data/development/typing/app")

# --- code tier 1: verified orphan, no surviving output ----------------------------------
CODE_TIER1 = [
    # -- the retired typing chain (root CLAUDE.md: "Retired, and not to be revived as fallbacks")
    _c(f"{TYPING}/substrate.py", 1,
       superseded_by=f"{TYPING}/subpanel.py",
       evidence="Zero importers. Its only mention in the tree is a historical line in "
                "docs/HANDOFF_2026-09-08_TYPING_REBUILD.md:115. Its output tree "
                "data/development/typing/substrate is already TIER2 above."),
    _c(f"{TYPING}/positivity.py", 1,
       superseded_by=f"{TYPING}/subpanel.py",
       evidence="Zero references of any kind in the tree. Output tree "
                "data/development/typing/positivity is already TIER2 above."),
    _c(f"{TYPING}/run_astir.py", 1,
       superseded_by=f"{TYPING}/subpanel.py",
       evidence="Zero references. Withdrawn producer comparison; data/development/typing/astir "
                "is already TIER2 above."),
    _c(f"{TYPING}/run_tacit.py", 1,
       superseded_by=f"{TYPING}/subpanel.py",
       evidence="Zero references. Line 44 drives broadtype/arms/tacit_run.R, which STAYS (broadtype "
                "is overturned below), so retiring this file takes nothing else with it. "
                "data/development/typing/tacit is already TIER2 above."),
    _c(f"{TYPING}/tacit_run.R", 1,
       superseded_by=f"{REPO}/core/phenocycler/scripts/broadtype/arms/tacit_run.R",
       evidence="An unused second copy. run_tacit.py:44 points at broadtype/arms/tacit_run.R, not "
                "at this one, and the two differ (md5 b97f0254... vs 512d8147...). Nothing in the "
                "tree names this path."),
    # -- superseded pLN Xenium typers -------------------------------------------------------
    _c(f"{XSCRIPTS}/pln_typing.py", 1,
       superseded_by=f"{XSCRIPTS}/pln_immune_typing_v3.py",
       evidence="Zero references outside its own usage block. docs/PLN_IMMUNE_XENIUM_2026-09-17.md "
                "names _v3 as the release producer."),
    _c(f"{XSCRIPTS}/pln_immune_typing.py", 1,
       superseded_by=f"{XSCRIPTS}/pln_immune_typing_v3.py",
       evidence="Zero references outside its own usage block. _v3 does not import it (its only "
                "cross-module import is core/phenocycler/scripts/pln/partial_run_guard.py)."),
    _c(f"{XSCRIPTS}/pln_immune_typing_v2.py", 1,
       superseded_by=f"{XSCRIPTS}/pln_immune_typing_v3.py",
       evidence="Zero references outside its own usage block. build_compartment_gradient.py:101 "
                "pinned its RELEASE TREE (release_pln_immune_v2), which is Lane X's repin, not an "
                "import of this file."),
    # -- legacy CODEX app chain, root scripts/, no surviving output ------------------------
    _c(f"{ROOTS}/recompute_viz_umap.py", 1,
       superseded_by=f"{ROOTS}/rebuild_adata_ins_root.py",
       evidence="Zero importers. Reads and rewrites data/adata_ins_root.h5ad IN PLACE; it produces "
                "nothing of its own. Both that path and its input data/islet_explorer.h5ad are "
                "absent from this checkout."),
    _c(f"{ROOTS}/compile_donors.py", 1,
       evidence="Zero references of any kind in the tree, including README.md and "
                "data/app_data/DATA_PROVENANCE.md. Names no output that exists."),
    _c(f"{ROOTS}/fix_tcell_subtypes.py", 1,
       evidence="Zero references. Its only input and output is "
                "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad, which is absent "
                "(single_cell_analysis/ holds one unexecuted notebook)."),
    _c(f"{ROOTS}/retrain_scvi_no_covariates.py", 1,
       evidence="Zero references. Reads the absent CODEX h5ad and data/islet_explorer.h5ad; a "
                "one-off retraining experiment that names no surviving output."),
    _c(f"{ROOTS}/plot_scvi_comparison.py", 1,
       evidence="Zero references. Reads data/islet_explorer.h5ad (absent) and writes figures."),
    _c(f"{ROOTS}/compute_pseudotime.R", 1,
       superseded_by=f"{TYPING}/trajectory_singlecell.R",
       evidence="Zero references. Writes data/pseudotime_results.rds, which does not exist here "
                "and is not read by the app under either domain."),
    # -- orphaned scripts/senior/, no surviving output ------------------------------------
    _c(f"{SENIOR}/build_app_h5ad.py", 1,
       evidence="Zero importers. Writes data/islet_explorer_senior.h5ad from "
                "islets_{core,peri}_senior.h5ad, islet_targets.parquet, islet_composition.parquet "
                "and neighborhood_metrics.csv -- none of those five paths exists in this checkout."),
    _c(f"{SENIOR}/build_phenotyped_h5ad.py", 1,
       evidence="Zero importers. In data/singlecell_protein.h5ad -> out "
                "data/singlecell_protein_phenotyped.h5ad; neither exists."),
    _c(f"{SENIOR}/build_singlecell_anndata.py", 1,
       evidence="Zero importers. Writes data/singlecell_protein.h5ad, which does not exist."),
    _c(f"{SENIOR}/compute_islet_targets.py", 1,
       evidence="Zero importers. Writes data/islet_targets.parquet, absent here and absent under "
                "both live domain roots."),
    _c(f"{SENIOR}/aggregate_islets.py", 1,
       prereq="line 22 imports create_islet_dataset_fixed from fixed_islet_aggregation.py, which "
              "is on this list; both move together or neither does",
       evidence="Zero importers. Writes data/islets_{core,peri}_senior.h5ad, neither of which "
                "exists. NOT the same file as "
                "core/phenocycler/scripts/typing/aggregate_islets.py, which build_layers.py:269 "
                "runs and which stays."),
    _c(f"{SENIOR}/compute_islet_markers.py", 1,
       evidence="Zero importers; the only mentions are docstring cross-references in two other "
                "entries on this list. Reads islets_{core,peri}_senior.h5ad, both absent."),
    _c(f"{SENIOR}/fixed_islet_aggregation.py", 1,
       blocked_by=(f"{SENIOR}/aggregate_islets.py", f"{ROOTS}/reaggregate_islets.py"),
       prereq="moving this alone breaks its two importers, both on this list -- "
              "scripts/senior/aggregate_islets.py (c1) and scripts/reaggregate_islets.py (c2), so "
              "the c2 pass runs first",
       evidence="Importers are scripts/reaggregate_islets.py:46 and "
                "scripts/senior/aggregate_islets.py:22 -- both on this list, so the reference is "
                "internal to the retiring set. Reads the absent CODEX h5ad."),
    _c(f"{SENIOR}/compute_trajectory_senior.py", 1,
       superseded_by=f"{TYPING}/trajectory_islet.py",
       evidence="Zero importers; mod_trajectory_server.R:587 and the trajectory audits name it in "
                "prose only. Reads islets_{core,peri}_senior.h5ad, both absent. Root CLAUDE.md "
                "forbids its gate shape: it rejects a build unless pseudotime anti-correlates "
                "with INS and orders ND<T1D, which is a result used as a verification target."),
    _c(f"{SENIOR}/train_scvi_senior.py", 1,
       evidence="Zero importers. Reads data/singlecell_protein_phenotyped.h5ad, which is absent; "
                "docs/EXTERNAL_PACKAGES.md:44 already records this lane as input-starved."),
    # -- the hold framework with no instance -------------------------------------------------
    _c(f"{REPO}/apps/R/ieCore/R/ie_holds.R", 1,
       prereq="TWO EDITS FIRST, neither in this lane's files: drop "
              "\"apps/R/ieCore/R/ie_holds.R\" from ALLOWLISTED_FILES in "
              "apps/islet/tests/testthat/test-guards.R:20 and "
              "apps/islet/tests/testthat/test-no-domain-pixel.R:15 (both read the file by path), "
              "and regenerate apps/R/ieCore/NAMESPACE with roxygen so ie_assert_hold_lifted and "
              "ie_hold_is_lifted stop being exported -- Lane G2 owns NAMESPACE",
       evidence="ie_assert_hold_lifted() and ie_hold_is_lifted() have ZERO callers: grep over "
                "every .R, .py and .md in the tree returns only their own definition. There is "
                "also no INSTANCE of the contract -- no file in the repo reads 'HOLD LIFTED' and "
                "no path is registered as a retired_pin. Not wired to the restain gate on "
                "purpose: that gate refuses by MEASUREMENT (build_restain_layer.py's five states, "
                "e.g. auroc_below_preregistered_bar, and donor 115's 7.198 um residual against "
                "3.67), and putting a hand-edited text file in front of it would let an edit "
                "override a measurement -- the opposite of the NO_VERIFIED_TYPING latch it says "
                "it borrows from, which is code no file can lift."),
]

# --- code tier 2: verified orphan, but a named producer of a LIVE artifact ----------------
# Quarantine only. Never deleted here: data/app_data is the CODEX domain's live root
# (apps/islet/config/islet_default.yml) and these are the only scripts in the tree that name its
# contents as output. Retiring them is honest -- every one is already unrunnable, because the
# chain's root input single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad and every
# intermediate it names under data/ are absent -- but DELETING them would leave a live domain with
# no recorded producer at all.
CODE_TIER2 = [
    _c(f"{ROOTS}/reaggregate_islets.py", 2,
       produces=("data/app_data/adata_ins_root.h5ad",),
       superseded_by=f"{ROOTS}/rebuild_adata_ins_root.py (run20 domain only)",
       evidence="Zero importers. Named producer of adata_ins_root.h5ad (line 40); the CODEX "
                "domain's copy exists at data/app_data/adata_ins_root.h5ad (2026-04-25). Its own "
                "declared output path data/adata_ins_root.h5ad does not exist, and its input "
                "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad is absent, so it "
                "cannot be re-run."),
    _c(f"{ROOTS}/build_h5ad_for_app.py", 2,
       produces=("data/app_data/islet_explorer.h5ad",),
       evidence="Zero importers. Sole named producer of islet_explorer.h5ad, which exists at "
                "data/app_data/islet_explorer.h5ad (2026-04-25, 65.6 MB) and is loaded by "
                "apps/islet/R/data_loading.R:406. Unrunnable: data/adata_ins_root.h5ad absent."),
    _c(f"{ROOTS}/convert_h5ad_to_parquet.py", 2,
       prereq="apps/islet/R/00_globals.R:167 tells the operator to run this script by name when "
              "islets.parquet is missing; that message must name a producer that exists",
       produces=("data/app_data/parquet/islets.parquet",
                 "data/app_data/parquet/uns_composition.parquet",
                 "data/app_data/parquet/uns_markers.parquet",
                 "data/app_data/parquet/uns_targets.parquet",
                 "data/app_data/parquet/tissue"),
       evidence="Zero importers, but apps/islet/R/00_globals.R:167 tells the operator to run it by "
                "name when islets.parquet is missing, and all five outputs exist under "
                "data/app_data/parquet (2026-06-02). Unrunnable: its input "
                "data/islet_explorer.h5ad is absent."),
    _c(f"{ROOTS}/compute_neighborhood_metrics.py", 2,
       produces=("data/app_data/neighborhood_metrics.csv",),
       blocked_by=(f"{SENIOR}/compute_neighborhood_senior.py",),
       prereq="moving this alone breaks compute_neighborhood_senior.py:15, which is on this list",
       evidence="Imported only by scripts/senior/compute_neighborhood_senior.py:15, which is on "
                "this same list. apps/islet/R/mod_inflammation_server.R:62 names it by name in an "
                "operator message and apps/islet/docs/spatial.md:18,26 documents its two constants "
                "as the definition of 'immune' for those columns. Output exists at "
                "data/app_data/neighborhood_metrics.csv (2026-07-02). Unrunnable: reads the absent "
                "CODEX h5ad."),
    _c(f"{ROOTS}/extract_per_islet_cells.py", 2,
       produces=("data/app_data/cells",),
       evidence="Zero importers. One of two named producers of the CODEX drill-down cells/ "
                "directory, which exists with 5,214 files (2026-03-05). Unrunnable: reads the "
                "absent CODEX h5ad."),
    _c(f"{ROOTS}/extract_per_donor_tissue.py", 2,
       produces=("data/app_data/donors",),
       evidence="Zero importers. Sole named producer of the CODEX donors/ directory, which exists "
                "with 15 files (2026-02-20) and is the input convert_h5ad_to_parquet.py's tissue "
                "partitions are built from. Unrunnable: reads the absent CODEX h5ad."),
    _c(f"{SENIOR}/compute_neighborhood_senior.py", 2,
       produces=("data/app_data/neighborhood_metrics.csv",),
       prereq="line 15 imports scripts/compute_neighborhood_metrics.py as cnm; that file is on "
              "this list and the two move together",
       evidence="Zero importers. A second named producer of neighborhood_metrics.csv, whose CODEX "
                "copy exists. Reads data/singlecell_protein_phenotyped.h5ad, which is absent."),
    _c(f"{SENIOR}/extract_islet_cells.py", 2,
       produces=("data/app_data/cells",),
       evidence="Zero importers. The second named producer of data/app_data/cells (line 24 writes "
                "there directly). Reads data/singlecell_protein_phenotyped.h5ad, which is absent."),
    _c(f"{SENIOR}/build_composition_and_tissue.py", 2,
       produces=("data/app_data/parquet/tissue",),
       evidence="Zero importers. A second named producer of parquet/tissue/, which the Spatial tab "
                "reads through apps/islet/R/00_globals.R:165 and spatial_helpers.R:6, and which "
                "exists (2026-06-02). Reads data/singlecell_protein_phenotyped.h5ad, absent."),
    _c(f"{SENIOR}/build_viewer_assets.py", 2,
       produces=("data/app_data/islet_spatial_lookup.csv", "data/app_data/json"),
       evidence="Zero importers, but it is the SOLE named producer of islet_spatial_lookup.csv, "
                "which exists (2025-09-30) and which apps/islet/R/segmentation_helpers.R:229-230 "
                "loads through domain_data_root() to resolve a clicked islet to a centroid. It is "
                "also the declared source of <case>_structures.geojson, the DEFAULT --structures "
                "argument of the live scripts/senior/build_all_tissuumaps_projects.py:154 (that "
                "directory does not currently exist). Unrunnable: its input "
                "data/islet_targets.parquet is absent."),
]

# --- code tier 3: proposed and OVERTURNED. Untouched. -------------------------------------
CODE_TIER3 = [
    (f"{REPO}/core/phenocycler/scripts/broadtype",
     "STILL IMPORTED BY LIVE CODE, measured 2026-09-18. "
     "core/phenocycler/scripts/kronos2_doc/run.py:35-40 puts scripts/broadtype and "
     "scripts/broadtype/arms on sys.path and imports declarations, pins, kronos2_coral and "
     "kronos_clusters. kronos2_doc is not a closed lane: its output "
     "data/development/kronos2_documented/ is read by core/phenocycler/scripts/report/"
     "kronos2_doc_fig.py:33, which run_report.py:83-84 renders as a report section. "
     "The approved plan lists broadtype as 6 dead entrypoints; measured, it is 47 files, 20 with "
     "a __main__ block, and one of them is a live dependency. NOTE the near miss: "
     "scripts/report/{bakeoff,identity_io}.py also `import pins`, but they insert their own "
     "directory first and mean scripts/report/pins.py -- three lanes carry a pins.py, and a "
     "basename match would have read those as broadtype importers too."),
    (f"{REPO}/core/phenocycler/scripts/dct",
     "OVERTURNED BY THE ABOVE. dct/regions.py is imported by "
     "core/phenocycler/scripts/broadtype/frame.py:34 (`import regions as dct_regions`) and "
     "broadtype/pins.py:20 puts scripts/dct on sys.path; with broadtype staying, dct is imported "
     "by code that stays. Everything else about it is clear: 22 files, 7 with a __main__ block, "
     "no importer outside dct/ and broadtype/, and the prerequisite in scripts/CLAUDE.md:42-43 is "
     "ALREADY SATISFIED -- gates.py and tiff_window.py live only in "
     "core/phenocycler/scripts/common/ (dct/ holds neither), and dct/{arms,concordance_gate,"
     "feasibility_tile,plane_identity}.py already resolve them by putting scripts/common on "
     "sys.path. So the move scripts/CLAUDE.md asks for is done; what blocks dct now is broadtype."),
    (f"{TYPING}/call.py",
     "STILL IMPORTED BY LIVE CODE: core/phenocycler/scripts/typing/compare_clustering.py:98 does "
     "`from call import TOP_N, _statistic, name_clusters`, and compare_clustering.py is on no "
     "retirement list -- it is the lane's Leiden/LDA/KODAMA method comparison. Retiring call.py "
     "would break it. Off the list until compare_clustering.py's own disposition is settled."),
    (f"{TYPING}/export_qupath.py",
     "Imported at four sites in core/phenocycler/scripts/typing/tests/test_typing.py (203, 214, "
     "257, 617). That file is a single suite that also covers compartments.py, lineages.py, "
     "metacluster.py and call.py, which are NOT retiring, so it cannot move with the module. "
     "Prerequisite before this can be promoted to c1: split those four tests out of test_typing.py."),
    (f"{TYPING}/validate.py",
     "Imported at two sites in core/phenocycler/scripts/typing/tests/test_typing.py (457, 467); "
     "same prerequisite as export_qupath.py. NOTE: the six `import validate as VAL` lines under "
     "core/phenocycler/scripts/restain/ resolve to restain/validate.py, a DIFFERENT module, and "
     "are not evidence about this file."),
]

# --- code tier 4: maintainer's call. Untouched. -------------------------------------------
CODE_TIER4 = [
    (f"{TYPING}/compare_clustering.py",
     "57 KB, the only importer of call.py, and itself referenced by nothing but "
     "cluster_leiden.py's docstring. Retiring it would unblock call.py; keeping it keeps call.py. "
     "Either way it is a ruling about the clustering comparison, not an orphan finding."),
    (f"{REPO}/scripts/build_master_excel.py",
     "Named producer of data/app_data/master_results.xlsx (2025-09-17), which exists. Not on the "
     "plan's list; recorded here only so the next pass does not have to rediscover it."),
]

#: Every code path this lane proposes to move. The `code_orphan` check treats a reference from
#: inside this set as internal, so a within-set import (broadtype -> dct/regions.py,
#: aggregate_islets -> fixed_islet_aggregation, compute_neighborhood_senior ->
#: compute_neighborhood_metrics) does not block the group.
CODE_RETIRING = [e["path"] for e in CODE_TIER1 + CODE_TIER2]

# --- producerless AND consumerless artifacts ----------------------------------------------
# Not code; listed here because they are the same finding: a file with no entrypoint on either
# side. Quarantine only -- nothing in this lane deletes them.
ARTIFACT_TIER2 = [
    _c(f"{REPO}/data/development/typing/app/master_results_run20.xlsx", 2,
       check="no_reference_anywhere",
       evidence="No producer and no consumer. build_workbook.py:192 writes "
                "master_results_run10.xlsx, not this name, and grep over every .py/.R/.md/.yml/.sh "
                "in the tree finds the string 'master_results_run20' nowhere. Written 2026-09-11 "
                "17:15, the same minute as the file below, by something no longer in the tree."),
    _c(f"{REPO}/data/development/typing/app/islet_celltype_proportions_run20.csv", 2,
       check="no_reference_anywhere",
       superseded_by=f"{REPO}/data/development/typing/app/islet_celltype_proportions.csv",
       evidence="No producer and no consumer. build_workbook.py:196 writes "
                "islet_celltype_proportions.csv (which exists, 2026-09-17 18:00, and is what "
                "trajectory_islet.py:112 reads); nothing names the _run20 spelling. Written "
                "2026-09-11 17:15."),
]
