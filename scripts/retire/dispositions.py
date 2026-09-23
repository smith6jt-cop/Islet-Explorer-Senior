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

Two lists below the code tiers are NOT tiers. `CODE_LIBRARY` is every path the 2026-09-22 bare-name
reference walk found REACHED BY LIVE CODE, with the importers named; `retire.py` refuses to move one
outright, and a directory entry protects everything under it. `CODE_CANDIDATE` is every path that
walk found with NO live reference and that is not already on a tier: recorded, with the walk's
evidence, for the maintainer's decision -- neither refused nor acted on. Both replaced the
`merge/hpg-sameslide` branch protection (KEEP_BESIDE_CODE), which lapsed when that branch merged.
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
#   code_orphan_and_config_blank(key)
#                        : code_orphan, AND core/phenocycler/config.ini has `key` blank. For an
#                          apparatus whose class column reaches the graph through a bound
#                          RELEASE DIRECTORY, which a walk over .py names cannot see.
#   candidate_orphan     : code_orphan with CODE_CANDIDATE treated as internal too; run on every
#                          CODE_CANDIDATE row, so the manifest carries the measurement the
#                          maintainer's decision rests on
#   no_reference_anywhere: the basename appears nowhere in tracked code, config or docs
#   library              : the path exists and is on CODE_LIBRARY (asserted, never acted on).
#                          Absence FAILS: every library was measured present on main.
#
# The `merge/hpg-sameslide` protection (KEEP_BESIDE_CODE, KEEP_BESIDE_CODE_TREES, HPG_BRANCH) was
# removed 2026-09-22. Its stated reason -- the HiPerGator session "is developing on it RIGHT NOW"
# -- has passed: `git branch --merged main` lists the branch and it has no commit main lacks. Every
# path it protected was walked for bare-name references on 2026-09-22 and became either a
# CODE_LIBRARY row (reached by live code; importers named) or a CODE_CANDIDATE row (no live
# reference; recorded, not moved). Both lists follow CODE_TIER4 below.


def _c(path, tier, produces=(), superseded_by="-", evidence="", check="code_orphan", prereq="",
       blocked_by=(), arg=None):
    """One code disposition.

    `blocked_by` names entries that must move in the SAME invocation. The direction is
    importer-first: if A imports B, moving B alone breaks A, so B is blocked_by A. `retire.py`
    refuses an entry whose blocker is not also being acted on, which is what keeps a
    within-set dependency from becoming a broken import halfway through.
    """
    return dict(path=path, tier=tier, produces=tuple(produces), superseded_by=superseded_by,
                evidence=evidence, check=check, arg=arg, prereq=prereq,
                blocked_by=tuple(blocked_by))


TYPING = f"{REPO}/core/phenocycler/scripts/typing"
XSCRIPTS = f"{REPO}/analysis/xenium/scripts"
SENIOR = f"{REPO}/scripts/senior"
ROOTS = f"{REPO}/scripts"
CELLTYPE = f"{REPO}/core/phenocycler/scripts/celltype"
RESTAIN = f"{REPO}/core/phenocycler/scripts/restain"
XCT = f"{REPO}/core/phenocycler/scripts/xenium_celltype"
HPG = f"{REPO}/core/phenocycler/scripts/xenium_celltype/hpg"
XTYPING = f"{REPO}/core/phenocycler/scripts/xenium_typing"
SAMESLIDE = f"{REPO}/core/phenocycler/scripts/sameslide"
PKG = f"{REPO}/core/phenocycler/phenocycler"

#: The one prerequisite the Xenium pancreas comparison apparatus (CODE_TIER2, xenium_typing/*)
#: shares. Met 2026-09-23 (Step 5); the check still reads the binding on every run and
#: refuses again if it is ever re-bound.
XTYPING_PREREQ = (
    "Prerequisite MET 2026-09-23 (Step 5): `xenium_identity_release` at core/phenocycler/"
    "config.ini:153 is blank, so nothing attaches this apparatus's class column to cells_xen "
    "(import_xenium.py:454-456 -> phenospace.load_release_or_none returns None), and the "
    "availability and distance_profiles nodes read Xenium calls from the declared releases "
    "only (build_availability.py xenium_availability; build_distance_profiles.py "
    "_pln_xenium_release_arm; the ingest fallbacks were removed the same day). The check "
    "code_orphan_and_config_blank now passes on the binding; what remains is the reference "
    "walk, and moving the files is `retire.py --kind code --tier 2 --execute`, which is the "
    "maintainer's to run.")

#: The two live domain roots (apps/islet/config/{islet_default,run20_default}.yml). An entry that
#: names an output which still EXISTS under one of these is c2, not c1.
LIVE_DOMAIN_ROOTS = (f"{REPO}/data/app_data", f"{REPO}/data/development/typing/app")

# --- code tier 1: verified orphan, no surviving output ----------------------------------
CODE_TIER1 = [
    # -- the retired typing chain (root CLAUDE.md: "Retired, and not to be revived as fallbacks")
    _c(f"{TYPING}/substrate.py", 1,
       superseded_by=f"{TYPING}/subpanel.py",
       evidence="Zero importers (walk 2026-09-22: DOC-ONLY). Named in prose only: the chain "
                "docstring shared by call.py:6, export_qupath.py:6, validate.py:6, "
                "metacluster.py:6 and positivity.py:6, typing/README.md:157, and "
                "docs/HANDOFF_2026-09-08_TYPING_REBUILD.md:139 (also :1, :15, :137). Its output "
                "tree data/development/typing/substrate is already TIER2 above."),
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
       evidence="Zero references outside the retiring set (tacit_run.R:2 names it in a comment). "
                "Line 50 (R_SCRIPT) drives broadtype/arms/tacit_run.R, which STAYS (broadtype "
                "is overturned below), so retiring this file takes nothing else with it. "
                "data/development/typing/tacit is already TIER2 above."),
    _c(f"{TYPING}/tacit_run.R", 1,
       superseded_by=f"{REPO}/core/phenocycler/scripts/broadtype/arms/tacit_run.R",
       evidence="An unused second copy. run_tacit.py:50 points at broadtype/arms/tacit_run.R, not "
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
       evidence="Zero references outside its own usage block. "
                "scripts/senior/build_compartment_gradient.py:144 records that it USED TO pin its "
                "RELEASE TREE (release_pln_immune_v2), which is Lane X's repin, not an import of "
                "this file."),
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
       prereq="line 34 imports create_islet_dataset_fixed from fixed_islet_aggregation.py, which "
              "is on this list; both move together or neither does",
       evidence="Zero importers. Writes data/islets_{core,peri}_senior.h5ad, neither of which "
                "exists. NOT the same file as "
                "core/phenocycler/scripts/typing/aggregate_islets.py, which "
                "scripts/senior/build_layers.py:302-304 runs (node aggregate_islets) and which "
                "stays."),
    _c(f"{SENIOR}/compute_islet_markers.py", 1,
       evidence="Zero importers; the only mentions are docstring cross-references in two other "
                "entries on this list. Reads islets_{core,peri}_senior.h5ad, both absent."),
    _c(f"{SENIOR}/fixed_islet_aggregation.py", 1,
       blocked_by=(f"{SENIOR}/aggregate_islets.py", f"{ROOTS}/reaggregate_islets.py"),
       prereq="moving this alone breaks its two importers, both on this list -- "
              "scripts/senior/aggregate_islets.py (c1) and scripts/reaggregate_islets.py (c2), so "
              "the c2 pass runs first",
       evidence="Importers are scripts/reaggregate_islets.py:59-60 and "
                "scripts/senior/aggregate_islets.py:34 -- both on this list, so the reference is "
                "internal to the retiring set. Reads the absent CODEX h5ad."),
    _c(f"{SENIOR}/compute_trajectory_senior.py", 1,
       superseded_by=f"{TYPING}/trajectory_islet.py",
       evidence="Zero importers; mod_trajectory_server.R:594 and the trajectory audits name it in "
                "prose only. Reads islets_{core,peri}_senior.h5ad, both absent. Root CLAUDE.md "
                "forbids its gate shape: it rejects a build unless pseudotime anti-correlates "
                "with INS and orders ND<T1D, which is a result used as a verification target."),
    _c(f"{SENIOR}/train_scvi_senior.py", 1,
       evidence="Zero importers. Reads data/singlecell_protein_phenotyped.h5ad, which is absent; "
                "docs/EXTERNAL_PACKAGES.md:45 already records this lane as input-starved."),
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
       evidence="Zero importers. Named producer of adata_ins_root.h5ad (line 53, "
                "DEFAULT_TRAJECTORY_OUT); the CODEX "
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
       prereq="moving this alone breaks compute_neighborhood_senior.py:28, which is on this list",
       evidence="Imported only by scripts/senior/compute_neighborhood_senior.py:28, which is on "
                "this same list. apps/islet/R/mod_inflammation_server.R:73 names it by name in an "
                "operator message, README.md:294 gives its run command, and "
                "apps/islet/docs/spatial.md:36,44 documents its two constants as the definition of "
                "'immune' for those columns. Output exists at "
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
       prereq="line 28 imports scripts/compute_neighborhood_metrics.py as cnm; that file is on "
              "this list and the two move together",
       evidence="Zero importers. A second named producer of neighborhood_metrics.csv, whose CODEX "
                "copy exists. Reads data/singlecell_protein_phenotyped.h5ad, which is absent."),
    _c(f"{SENIOR}/extract_islet_cells.py", 2,
       produces=("data/app_data/cells",),
       evidence="Zero importers. The second named producer of data/app_data/cells (line 37, "
                "CELLS_OUT, writes there directly). Reads data/singlecell_protein_phenotyped.h5ad, "
                "which is absent."),
    _c(f"{SENIOR}/build_composition_and_tissue.py", 2,
       produces=("data/app_data/parquet/tissue",),
       evidence="Zero importers. A second named producer of parquet/tissue/, which the Spatial tab "
                "reads through apps/islet/R/00_globals.R:165 and spatial_helpers.R:6, and which "
                "exists (2026-06-02). Reads data/singlecell_protein_phenotyped.h5ad, absent."),
    _c(f"{SENIOR}/build_viewer_assets.py", 2,
       produces=("data/app_data/islet_spatial_lookup.csv", "data/app_data/json"),
       evidence="Zero importers, but it is the SOLE named producer of islet_spatial_lookup.csv, "
                "which exists (2025-09-30) and which apps/islet/R/segmentation_helpers.R:234-235 "
                "loads through domain_data_root() to resolve a clicked islet to a centroid. It is "
                "also the declared source of <case>_structures.geojson, the DEFAULT --structures "
                "argument of the live scripts/senior/build_all_tissuumaps_projects.py:167-170 "
                "(that directory does not currently exist). Unrunnable: its input "
                "data/islet_targets.parquet is absent."),
    # -- the Xenium pancreas comparison apparatus (arms X1 / X3 and its release), declared
    #    2026-09-22 from the producers census. Superseded IN WRITING: analysis/xenium/scripts/
    #    panc_xenium_release.py:86-91 declares that it SUPERSEDES data/development/xenium_typing/
    #    release (the scANVI-per-donor arm over 6436 and 6442, column lineage_native). Until
    #    2026-09-23 the release was BOUND, which is why the group is c2 with a prerequisite and
    #    not c1: core/phenocycler/config.ini:153 `xenium_identity_release = ../../data/development/
    #    xenium_typing/release` -> integration/import_xenium.py:454-456 load_release_or_none /
    #    attach_release_identity -> integration/phenospace.py:532-533 wrote lineage_native and
    #    lineage_common onto cells_xen -> the availability and distance_profiles nodes read them
    #    (build_availability.py's xenium_ingest_availability, build_distance_profiles.py's
    #    _pln_arm, both since removed) -> apps/islet/R/card_x1_gradient.R:143,271 rendered
    #    resolved_by. Step 5 blanked :153 and removed both readers on 2026-09-23. The check
    #    `code_orphan_and_config_blank` reads the binding itself and now passes on it; the
    #    entries stay c2 and move only when the maintainer executes the tier.
    #    joint_model.py (X2) is CODE_TIER3: the report lane reads its latent.
    _c(f"{XTYPING}/release.py", 2, check="code_orphan_and_config_blank",
       arg="xenium_identity_release",
       produces=("data/development/xenium_typing/release",),
       superseded_by=f"{XSCRIPTS}/panc_xenium_release.py",
       blocked_by=(f"{XSCRIPTS}/pln_typing.py",),
       prereq=XTYPING_PREREQ,
       evidence="Writes the bound release: release.py:175 --out defaults to pins.OUT/'release' "
                "(pins.py:60 OUT = data/development/xenium_typing), :149 declares the columns "
                "lineage_native / common_lineage / evidence_positive. Live readers of that "
                "release are the config.ini:153 binding chain above and the report lane "
                "(identity_io.py:19 RELEASE['xenium']; xen_identity.py:75,175 read "
                "lineage_native). Live importers of the MODULE: the default-lane suite "
                "xenium_typing/tests/test_xenium_typing.py:129 (`import release`; Makefile:32 "
                "LANES includes xenium_typing) and analysis/xenium/scripts/pln_typing.py:216 "
                "(CODE_TIER1, internal to the retiring set). It imports the declared libraries "
                "xenio (:46) and build_phenotyped_base (:47)."),
    _c(f"{XTYPING}/nullgate_typer.py", 2, check="code_orphan_and_config_blank",
       arg="xenium_identity_release",
       produces=("data/development/xenium_typing/arms/scanvi-per-donor",),
       superseded_by=f"{XSCRIPTS}/panc_xenium_cluster_name.py",
       blocked_by=(f"{XSCRIPTS}/pln_typing.py", f"{XTYPING}/scanvi_replicates.py"),
       prereq=XTYPING_PREREQ + " Moves with hierarchical_phenotype.py, whose sidecar it reads.",
       evidence="Arm X1 (nullgate_typer.py:31 ARM = pins.ARM_SCANVI, 'scanvi-per-donor' at "
                "pins.py:51). Nothing is re-run: :3-5 and :43 read lineage_native straight off "
                "analysis/xenium/scripts/hierarchical_phenotype.py's <serial>_hierarchical.h5ad "
                "sidecar (obs['celltype_broad']); :68 records that source. Importers: "
                "scanvi_replicates.py:102 (on this list) and analysis/xenium/scripts/"
                "pln_typing.py:204 (CODE_TIER1). It is what makes the X1 class column the head "
                "of the chain that reaches cells_xen."),
    _c(f"{XTYPING}/scanvi_replicates.py", 2, check="code_orphan_and_config_blank",
       arg="xenium_identity_release",
       superseded_by=f"{XSCRIPTS}/panc_xenium_cluster_name.py",
       prereq=XTYPING_PREREQ,
       evidence="The X1 replicate-stability runner: :102-109 imports nullgate_typer, runs the "
                "seeds and writes the replicate table through nullgate_typer.run. Walk "
                "2026-09-22: DOC-ONLY -- its only inbound reference is the docstring "
                "nullgate_typer.py:57; no importer, no runner, no test. Outbound it imports the "
                "declared library tissue_markers/pancreas.py by importlib at :46. Listed with the "
                "group because it is runnable only against arm X1."),
    _c(f"{XTYPING}/baseline_analyzer.py", 2, check="code_orphan_and_config_blank",
       arg="xenium_identity_release",
       produces=("data/development/xenium_typing/arms/xenium-analyzer-baseline",),
       superseded_by=f"{XSCRIPTS}/panc_xenium_cluster_name.py",
       prereq=XTYPING_PREREQ + " Its one importer is the lane test; that test leaves with it.",
       evidence="Arm X3 (baseline_analyzer.py:32 ARM = pins.ARM_BASELINE; xenium_typing/"
                "README.md:14). :68 writes lineage_native = gate.UNAVAILABLE for every cell until "
                "the Globus copy lands, so no transferred output exists. Sole importer: the "
                "default-lane suite xenium_typing/tests/test_xenium_typing.py:140 "
                "(`import baseline_analyzer as ba`). The census read this one as c1; it is "
                "listed c2 so the three arms and the release carry one prerequisite and move "
                "together."),
    # -- the Xenium notebook/HPC typer, declared 2026-09-22 from the producers census -----------
    _c(f"{XSCRIPTS}/hierarchical_phenotype.py", 2,
       produces=("data/processed/<serial>/<serial>_hierarchical.h5ad",),
       superseded_by=f"{XSCRIPTS}/panc_xenium_cluster_name.py",
       blocked_by=(f"{XTYPING}/nullgate_typer.py",),
       prereq="Close the xenium_typing X1 chain first (nullgate_typer.py:3-5,43 reads this "
              "file's sidecar; moving this alone makes that read producerless), and keep the "
              "notebook/HPC lane's use named: analysis/xenium/hpc/run_pipeline.sh:163-164 runs it "
              "as the `phenotype` stage, analysis/xenium/scripts/refresh_broad_umap.py:49 imports "
              "it by importlib, render_stage02_phenotyping.py:213,362 reads its celltype_broad, "
              "and analysis/xenium/CLAUDE.md:196 documents that lane as separate from the cohort "
              "lane. Those are the maintainer's decisions, not this lane's.",
       evidence="Retired as a Xenium pancreas typer because the declared producer "
                "(panc_xenium_cluster_name.py under panc_xenium_release.py) replaces it and says "
                "so at panc_xenium_release.py:86-91. Writes data/processed/<s>/<s>_hierarchical."
                "h5ad (:48, :1105) with celltype_broad (:720) and celltype (:944). It is the only "
                "lane in the tree that splits Exocrine into the two sub-labels the root CLAUDE.md "
                "forbids (:38-40; tissue_markers/pancreas.py:318), and its own :83-84 records that "
                "the current pancreas release carries a single Exocrine lineage. Live uses today "
                "are the HPC lane and X1 named in prereq; tissue_markers/pancreas.py:427,490,609 "
                "declare their tables consumed by it. c2, not c1: it is the named producer of "
                "the sidecar arm X1 still reads."),
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
     "no importer outside dct/ and broadtype/, and the prerequisite in scripts/CLAUDE.md:111-112 "
     "is ALREADY SATISFIED -- gates.py and tiff_window.py live only in "
     "core/phenocycler/scripts/common/ (dct/ holds neither), and dct/{arms,concordance_gate,"
     "feasibility_tile,plane_identity}.py already resolve them by putting scripts/common on "
     "sys.path. So the move scripts/CLAUDE.md asks for is done; what blocks dct now is broadtype."),
    (f"{TYPING}/call.py",
     "STILL IMPORTED BY LIVE CODE, and a library for exactly one consumer: "
     "core/phenocycler/scripts/typing/compare_clustering.py:105 does "
     "`from call import TOP_N, _statistic, name_clusters` -- the naming statistic the whole "
     "Leiden/LDA/KODAMA comparison is built on -- and cluster_kodama.py:45,67 pins its own "
     "graph_k / resolution to call.py's defaults. The default-lane suite "
     "core/phenocycler/scripts/typing/tests/test_typing.py imports it at 20 sites (e.g. :183). "
     "Its OUTPUT tree data/development/typing/calls is already TIER2 (withdrawn). Retiring the "
     "module would break the comparison; compare_clustering.py is a CODE_CANDIDATE below, so the "
     "disposition that closes call.py is that one, not another reference walk."),
    (f"{XTYPING}/joint_model.py",
     "STAYS while the report lane reads what it wrote. It is arm X2 of the retired Xenium "
     "pancreas comparison (joint_model.py:51 ARM = pins.ARM_JOINT; :62 OUT = "
     "pins.OUT/'joint_model'; :289 writes lineage_native + celltype_fine; the one scvi-env "
     "entrypoint in the set, :7-8). "
     "No live module imports it, but core/phenocycler/scripts/report/xen_identity.py:159-164 "
     "reads joint_model/train_latent.npy and train_obs.parquet to render the X05 latent figure "
     "that run_report.py drives. Its sibling arms X1/X3 and release.py are CODE_TIER2 with a "
     "named prerequisite; this file is held back separately because moving it makes a rendered "
     "report input producerless. It also imports the declared library "
     "analysis/xenium/scripts/tissue_markers/pancreas.py at :250."),
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
    # typing/compare_clustering.py moved to CODE_CANDIDATE on 2026-09-22 (walk: DOC-ONLY).
    (f"{REPO}/scripts/build_master_excel.py",
     "Named producer of data/app_data/master_results.xlsx (2025-09-17), which exists. Not on the "
     "plan's list; recorded here only so the next pass does not have to rediscover it."),
]

# --- declared libraries: reached by live code, measured 2026-09-22. Refused outright. ---------
# One row per path the bare-name reference walk classified LIVE (or TEST-ONLY through a test the
# default lane runs), with the importers named. A directory entry protects everything under it.
# `retire.py` refuses these the way it refuses TIER3 and TIER4; `build_manifest.py` writes them
# with tier "LIBRARY" and the `library` check, so a library that disappears from disk FAILS the
# manifest instead of silently dropping out of the guard. Nothing here is a producer ruling: the
# four producers (one per tissue x modality) are named in the test suite and in the census, and
# every alternative typer below is declared a library WITH the use that keeps it, or a candidate.
CODE_LIBRARY = [
    # -- celltype: the panel authority and the rule/score/usability modules every lane imports --
    (f"{CELLTYPE}/panel.py",
     "33 `import panel` sites resolve here through each lane's sys.path insert of "
     "scripts/celltype: broadtype/arm_base.py:46, broadtype/arms/cosine_broad8.py:65, "
     "report/pins.py:40, typing/refine.py:127, celltype/run_sameslide_protein.py:70, and the "
     "default-lane test celltype/tests/test_usability.py:25 (Makefile:32 LANES includes "
     "celltype)."),
    (f"{CELLTYPE}/rules.py",
     "Imported by both PhenoCycler producers' lanes and three comparison lanes: "
     "typing/subpanel.py:145 and typing/immune_subtype.py:95 (`import rules as RULES`), "
     "typing/refine.py:128, broadtype/arms/cosine_broad8.py:66, xenium_celltype/"
     "eval_detection.py:123, celltype/run_sameslide_protein.py:71; "
     "celltype/marker_name_crosswalk.json:4 names it as the regenerate command."),
    (f"{CELLTYPE}/score.py",
     "Imported by typing/refine.py:129, broadtype/arms/cosine_broad8.py:68 and "
     "celltype/run_sameslide_protein.py:74; celltype/README.md:362 documents it as the scale, "
     "cosine and requirement gate."),
    (f"{CELLTYPE}/usability.py",
     "Imported by xenium_typing/gate.py:42 and broadtype/gate.py:28 (both `import usability` "
     "with scripts/celltype on sys.path), celltype/run_step1.py:25, "
     "celltype/run_sameslide_protein.py:146 and the default-lane test "
     "celltype/tests/test_usability.py:26; broadtype/criteria.py:10 quotes its tail definition."),
    (f"{CELLTYPE}/run_step3.py",
     "A LIBRARY, not a producer. Its live use is build_scale / REQUIRE_MIN / NUCLEAR_RATIO / "
     "ALIGNED, imported by typing/refine.py:130, celltype/figures_ambiguous.py:36, "
     "celltype/run_islet.py:45, celltype/run_sameslide_protein.py:73, "
     "broadtype/arms/cosine_broad8.py:67 and the default-lane test "
     "tests/current/test_refine.py:87; typing/RECONSTRUCTION.md:500 and REFINEMENT.md:60 cite "
     "REQUIRE_MIN as 'imported and not restated'. Its OUTPUT root, run_step3.py:25 OUT = "
     "panel.OUT/'calls' (celltype/README.md:363), is a pancreas typing arm with zero readers "
     "and, per docs/WORKFLOW_2026-09-02_MAP.md:77, donor 6450 only and older than its own "
     "inputs -- a DATA question for the maintainer, recorded here and not acted on."),
    (f"{CELLTYPE}/run_islet.py",
     "An AUDIT instrument with its use named, not a competing producer. Module importers: "
     "celltype/figures_ambiguous.py:37, celltype/run_sameslide_protein.py:72 and the "
     "default-lane test celltype/tests/test_usability.py:649. DATA readers of what it writes "
     "(run_islet.py:47 OUT = panel.OUT/'islet_residents'; :227 islet_residents_all.parquet): "
     "kintsugi/concordance.py:40 COSINE_ISLET, read by kintsugi/run_step4.py:36 (runnable, "
     "__main__ at :61 -- the kintsugi method comparison) and celltype/figures_ambiguous.py:171. "
     "scripts/senior/build_layers.py declares NO node input under data/development/celltype "
     "other than uniform/ (:160,:340,:447), so islet_residents_all.parquet is a two-lane "
     "comparison artifact and never a source of cell_class for any app artifact. Its v2/v3/leaf "
     "variants are TIER2 above."),
    (f"{CELLTYPE}/run_step1.py",
     "The celltype lane's Step 1 entrypoint (celltype/README.md:20; docs/WORKFLOW_2026-09-02_"
     "MAP.md:75). Imported by the default-lane test celltype/tests/test_usability.py:685 "
     "(`import run_step1 as S1`); its command is printed by run_sameslide_protein.py:144."),
    # -- restain: the registration modules the live restain producers and tests import --------
    (f"{RESTAIN}/align.py",
     "Imported by restain/align_figures.py:41, run_align.py:39, run_field_refine.py:53 and the "
     "default-lane tests restain/tests/test_align_inventory.py:15 and test_restain.py:406 "
     "(Makefile:32 LANES includes restain)."),
    (f"{RESTAIN}/registration_inputs.py",
     "Imported by ten restain producers and six tests, e.g. restain/explorer.py:443, "
     "run_align.py:41, run_field_refine.py:54, tests/test_align_inventory.py:16, "
     "tests/test_intensity_refine.py:16, tests/test_restain.py:701; "
     "docs/WORKFLOW_2026-09-02_MAP.md:116 documents it."),
    (f"{RESTAIN}/run_align.py",
     "Its symbols are imported by restain/run_hormone.py:88 (_refined_transform, "
     "_valis_transform), run_field_refine.py:55, tests/test_align_inventory.py:17 and "
     "tests/test_intensity_refine.py:17; run_field_refine.py:16 gives its command."),
    (f"{RESTAIN}/run_intensity_refine.py",
     "Imported by restain/run_field_refine.py:61 (`import run_intensity_refine as IR`) and the "
     "default-lane test restain/tests/test_intensity_refine.py:18; restain/README.md:100 gives "
     "its run command."),
    (f"{RESTAIN}/run_register.py",
     "17 import sites pull OUT / GEOJSON_DIR / XENIUM_PATHS / discover / donor_of / "
     "register_one: restain/channel_check.py:39, align_figures.py:49, run_align.py:44, "
     "run_field_refine.py:62 among them; restain/README.md:871 documents it. (The package "
     "function `run_register` in phenocycler/integration is a different symbol.)"),
    # -- the integration package's transform, imported package-wide -----------------------------
    (f"{PKG}/integration/transform.py",
     "Package-internal imports (integration/register.py:51, integration/sameslide.py:61), the "
     "restain lane (run_align.py:167, run_field_refine.py:64), the sameslide tree "
     "(build_joined.py:34, lift_calls.py:48, register_points.py:42), the integration tests "
     "(tests/integration/test_integration_register.py:23) and the live "
     "scripts/senior/same_section_linkage.py:191,251."),
    # -- xenium_celltype: the four modules imported from outside the tree, their own imports,
    #    and the two test directories the default lane runs ------------------------------------
    (f"{XCT}/xenio.py",
     "The Xenium loader every lane reaches: restain/run_register.py:47, protein_rna.py:78, "
     "align_figures.py:46 (each with scripts/xenium_celltype on sys.path), "
     "xenium_typing/release.py:46, and in-tree run_donor.py:62, build_phenotyped_base.py:63, "
     "hpg/stage_call.py:41, hpg/stage_prep.py:36."),
    (f"{XCT}/nullgate.py",
     "Imported from outside the tree by xenium_typing/gate.py:43 and kintsugi/annotate.py:53 "
     "(sys.path insert at :50; kintsugi/annotate.py:21 says it is 'imported and called, never "
     "copied'); in-tree by run_donor.py:63, tests/test_nullgate.py:7 and "
     "hpg/tests/test_nullgate_fast.py:16."),
    (f"{XCT}/build_phenotyped_base.py",
     "Stage 2 of the Xenium method of record (docs/XENIUM_PIPELINE.md:151). Imported from "
     "outside the tree by report/xen_io.py:29 (the report lane) and "
     "xenium_typing/release.py:47; in-tree by tests/test_build_phenotyped_base.py:17 (its :234 "
     "says the test pins it). analysis/xenium/scripts/pln_typing.py:157 also imports it, but "
     "that file is CODE_TIER1."),
    (f"{HPG}/nullgate_fast.py",
     "The production gate: imported from outside the tree by xenium_typing/gate.py:46 "
     "(`from nullgate_fast import group_gene_deltas`), in-tree by hpg/stage_call.py:42 and "
     "hpg/tests/test_nullgate_fast.py:17; nullgate.py:244,249 and docs/XENIUM_PIPELINE.md:407 "
     "record that the two gates are asserted to agree."),
    (f"{XCT}/gdna.py",
     "Imported by the declared library build_phenotyped_base.py:200 (`import gdna as gdna_mod`) "
     "and by eval_detection.py:79; docs/XENIUM_PIPELINE.md:138 gives its operator command "
     "(`python3 scripts/xenium_celltype/gdna.py`, stage 1 of the method of record)."),
    (f"{XCT}/run_donor.py",
     "Imported by the default-lane test tests/test_nullgate.py:164 (`from run_donor import "
     "composition`) and in-tree by hpg/stage_call.py:40 and hpg/stage_prep.py:37; "
     "docs/XENIUM_PIPELINE.md:346 documents run_donor.composition's status vocabulary."),
    (f"{XCT}/run_kodama.R",
     "Called by the declared library run_donor.py:173 (`cmd = ['Rscript', str(HERE / "
     "'run_kodama.R'), ...]`); docs/WORKFLOW_2026-09-02_MAP.md:103 and "
     "docs/HANDOFF_2026-08-26_XENIUM_CELLTYPE.md:182 document it. DUP BASENAME: "
     "analysis/xenium/scripts/benchmark/methods/run_kodama.R is a different file."),
    (f"{XCT}/tests",
     "The lane's own suite, run by the default lane: Makefile:32 LANES includes xenium_celltype "
     "and Makefile:93-96 runs scripts/<lane>/tests. It imports build_phenotyped_base "
     "(test_build_phenotyped_base.py:17), nullgate and run_donor (test_nullgate.py:7,164) and "
     "hpg/verify_run.py (test_verify_run.py:24). A directory entry: everything under it is "
     "protected."),
    (f"{HPG}/tests",
     "Run by the default lane by name: Makefile:99 `pytest scripts/xenium_celltype/hpg/tests`. "
     "test_nullgate_fast.py:16-17 imports nullgate and nullgate_fast and is the proof "
     "nullgate.py:249 cites that the two gates agree. A directory entry."),
    (f"{HPG}/verify_run.py",
     "Imported by the default-lane test xenium_celltype/tests/test_verify_run.py:24 "
     "(`import verify_run as V`); invoked in-tree by hpg/slurm/50_assemble.sh:20 and "
     "hpg/submit_all.sh:141; docs/XENIUM_PIPELINE.md:349 names its check (8)."),
    # -- the rest of the former branch protection -------------------------------------------
    (f"{XSCRIPTS}/tissue_markers/pancreas.py",
     "Imported by xenium_typing/joint_model.py:250 (`from tissue_markers import pancreas as "
     "mod`) and by importlib at xenium_typing/scanvi_replicates.py:46 "
     "('scripts.tissue_markers.pancreas'); its tables are declared consumed by "
     "hierarchical_phenotype.py at pancreas.py:427,490,609."),
    (f"{REPO}/core/phenocycler/tests/current/test_no_duplicate_tracked_data.py",
     "A live default-lane test: collected by Makefile:59 (`cd core/phenocycler && $(PYTEST) "
     "-q -rfEs tests/current`). Tests are collected, not imported; nothing needs to name it."),
    (f"{REPO}/docs/XENIUM_PIPELINE.md",
     "The Xenium method of record, cited by code and records: "
     "analysis/xenium/scripts/pln_immune_typing.py:26, docs/PANC_XENIUM_TYPING_2026-09-17.md:422; "
     "it carries the operator commands for gdna.py (:138) and run_sameslide_protein.py "
     "(:245-246) and the stage-2 heading for build_phenotyped_base.py (:151)."),
    # -- alternative typers the producers census names, declared libraries with the use that
    #    keeps them (2026-09-22). None is a fifth producer. -----------------------------------
    (f"{PKG}/simultaneous_typing.py",
     "A library behind the NO_VERIFIED_TYPING latch, not a competing PhenoCycler pancreas "
     "typer. Its header (:3-6) reads RETIRED -- NOT A FALLBACK, and what stops it producing is "
     "code: the `type` stage is blocked by NO_VERIFIED_TYPING (root CLAUDE.md, 'that latch has "
     "not moved'). Live importers, both inside the installed package by bare module name: "
     "phenocycler/pipeline.py:75-78 (`from .simultaneous_typing import BroadGatePolicy, "
     "load_broad_gate_policy`; :182 records it as a source string) and "
     "phenocycler/production_typing.py:21-25 (METHOD_VERSION, BroadGatePolicy, "
     "type_simultaneous_cells). Retiring it breaks the package's own pipeline import -- a "
     "failure the code_orphan check only sees when the walk searches the bare module name."),
    (f"{REPO}/core/phenocycler/scripts/pln/cluster_pln.py",
     "The clearest library in the set: it writes a PARTITION (cluster_pln.py:173 OUT_ROOT = "
     "data/development/typing/pln/clusters; cluster_id at :595,:784) and no class column, and "
     "it is upstream of BOTH declared pLN producers. build_pln_immune_pheno.py:166 CLUSTERS "
     "reads its phenocycler tree (:22 'consumed here'; docs/PLN_IMMUNE_PHENOCYCLER_"
     "2026-09-17.md:32), and pln_immune_typing_v3.py:556-562 loads the FILE BY PATH "
     "(importlib spec_from_file_location 'cluster_pln_for_v3') for leiden_on_subsample and "
     "stable_seed (:43, :1491). Tests read it by path: tests/current/test_pln_layer.py:572,628 "
     "and test_pln_immune_pheno.py:1087. Run as docs/PLN_REGIONAL_CLUSTERING_2026-09-17.md:"
     "569-570. A by-path importlib load is invisible to an import-statement walk."),
    (f"{TYPING}/refine.py",
     "A library and the VOCABULARY AUTHORITY for the ninth class, not a fifth producer. "
     "`refined_class` has ZERO consumers (refine.py:5; modality_concordance.py:29,212-213); "
     "build_cells.py:507 `--refine` is the join that would wire it and has not been run "
     "(build_cells.py:148 REFINE_CARRY). What is live is the label string refine.py:150 "
     "UNRESOLVED = 'Unresolved (no lineage evidence)', read as the committed spelling by "
     "integration/vocab.py:206-209 (PHENO_UNRESOLVED), scripts/senior/build_grid_composition.py:"
     "272 (PANC_CLASS_NONCALL), modality_concordance.py:186 ('AUTHORITY. The class vocabulary "
     "is ... refine.py'), apps/islet/config/run20_default.yml:42,57, apps/islet/R/00_globals.R:"
     "653 and apps/R/ieCore/R/ie_availability.R:18; and the process-wide guard it installs at "
     "import (:19) against the producerless refinement_log.json. It must NOT be wired as a "
     "second PhenoCycler pancreas producer while refined_class has no consumer and its "
     "5th-percentile bound demotes about 5% of every class by construction -- root CLAUDE.md "
     "calls that an open question, not a decision. Its own imports are panel/rules/score/"
     "run_step3 (:127-130), all declared libraries above."),
]

# --- retirement candidates: NO live reference measured 2026-09-22. Recorded, not moved. -----
# One row per path the walk classified NONE, DOC-ONLY, or reached only from other rows on this
# list (in-tree only), and that is on no tier. `retire.py` neither refuses nor acts on a
# candidate -- it reads tiers 1 and 2 only -- and `build_manifest.py` writes each with tier
# "CANDIDATE" and the `code_orphan` check, so the manifest shows today's measurement beside the
# evidence. Promotion to a tier is the maintainer's decision at the end of Step 4.
CODE_CANDIDATE = [
    (f"{SAMESLIDE}",
     "Whole tree (5 files), walk 2026-09-22: no importer outside it. lift_calls.py: DOC-ONLY -- "
     "docstring of the live scripts/senior/same_section_linkage.py:16 ('read not run here') and "
     "docs/SAME_SECTION_LINKAGE_2026-09-17.md:26; it imports build_joined:47 and "
     "integration/transform:48. register_points.py: NONE -- zero inbound references; its :219 "
     "and :224 help text name read_pheno_tsv and extract_pheno_geometry.py, the only references "
     "to those two files. build_joined.py: in-tree only, imported by lift_calls.py:47. The "
     "tree has no tests and no lane entry in Makefile:32. The former protection reason "
     "(lift_calls.py changed on merge/hpg-sameslide) lapsed with the merge. Outbound it "
     "imports the declared library integration/transform.py (build_joined.py:34, "
     "register_points.py:42)."),
    (f"{RESTAIN}/run_field_refine.py",
     "TEST-ONLY via the retire guard itself: before 2026-09-22 its only references were "
     "tests/current/test_retire_code_tier.py asserting it on the keep-beside list and "
     "build_manifest.py tolerating its absence. No importer, no invocation, no doc. It imports "
     "six declared libraries (align:53, registration_inputs:54, run_align:55, "
     "run_intensity_refine:61, run_register:62, integration/transform:64) and its :4 docstring "
     "contrasts itself with run_intensity_refine.py. The former reason ('the branch adds it, "
     "398 lines') lapsed with the merge."),
    (f"{CELLTYPE}/run_sameslide_protein.py",
     "DOC-ONLY: docs/XENIUM_PIPELINE.md:239,245-246 describe it and give its command; "
     "docs/XENIUM_PIPELINE.md:367 cites run_sameslide_protein.islet_regions. No importer, no "
     "test, no invocation. Its output path (run_sameslide_protein.py:82 OUT = panel.OUT/"
     "'sameslide_protein') is named literally by sameslide/lift_calls.py:28, itself a "
     "candidate above, and lift_calls.py:52 CARRY lifts its `call` column onto Xenium cells. "
     "It is the two-donor SAME-SLIDE BRIDGE (6533 and 112), not a cohort producer: its call "
     "column reaches no build_layers node and no R card. It imports the declared libraries "
     "panel/rules/run_islet/run_step3/score/usability (:70-74, :146)."),
    (f"{TYPING}/compare_clustering.py",
     "DOC-ONLY (moved here from CODE_TIER4, where it was 'the maintainer's call'): no importer, "
     "no invocation, no test; docstrings of cluster_kodama.py, cluster_lda.py and "
     "cluster_leiden.py name it. It is the load-bearing referrer of call.py "
     "(compare_clustering.py:105 `from call import TOP_N, _statistic, name_clusters`), so "
     "retiring it would unblock call.py (CODE_TIER3) and keeping it keeps call.py. Either way "
     "it is a ruling about the clustering comparison, not an orphan finding."),
    # -- xenium_celltype: the detection-evaluation trio, in-tree only ------------------------
    (f"{XCT}/eval_detection.py",
     "DOC-ONLY: no importer, no runner, no test. Cited by name in the tracked "
     "analysis/xenium/panel_audit.csv (rows 261, 2459, 2460, 4855, 4959 give 'eval_detection.py's "
     "CLEAN_NEG curation' as the reason for an exclude flip). It imports context:77, "
     "detect:78, gdna:79 (a library) and rules:123 (a library) with its own dir on sys.path "
     "(:76). The former reason ('the branch develops it') lapsed with the merge."),
    (f"{XCT}/context.py",
     "In-tree only: the sole reference is eval_detection.py:77 (`import context as CTX`), a "
     "candidate. DUP BASENAME with core/phenocycler/scripts/typing/context.py, which owns every "
     "typing-lane hit."),
    (f"{XCT}/detect.py",
     "In-tree only: the sole reference is eval_detection.py:78 (`import detect as DET`), a "
     "candidate; eval_detection.py:4 and context.py:4 cite it in prose."),
    # -- xenium_celltype/hpg: the HiPerGator submission apparatus, in-tree only ---------------
    #    Every reference between these files stays inside this list; the only outside mentions
    #    are docs/XENIUM_PIPELINE.md:234,349,351 (prose) and the pinned-installer notes
    #    docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md:213 / docs/EXTERNAL_PACKAGES.md:16. Whether
    #    the submission path is still operated on HiPerGator is a fact the tree cannot measure,
    #    which is why these are candidates and not a tier.
    (f"{HPG}/fetch_results.sh",
     "NONE: zero inbound references (its only hits are its own usage lines); it names "
     "submit_all.sh at :15."),
    (f"{HPG}/sync_code.sh",
     "NONE: zero inbound references (its own usage lines only). It sources env.sh:22 and names "
     "setup_envs.sh:53, submit_all.sh:11,58 and verify_run.py:54."),
    (f"{HPG}/slurm",
     "Whole tree (7 scripts), in-tree only: each is submitted by hpg/submit_all.sh -- "
     "00_gate.sh at :104, 10_prep.sh :108, 05_calibrate.sh :112, 20_shard.sh :117, "
     "30_merge.sh :123, 40_call.sh :127, 50_assemble.sh :131. They source env.sh "
     "(00_gate.sh:17) and run stage_prep.py (10_prep.sh:37), stage_call.py (40_call.sh:34), "
     "verify_run.py (50_assemble.sh:20, a library), check_shard_equivalence.R (00_gate.sh:29), "
     "precompute_knn.R (10_prep.sh:43), kodama_shard.R (20_shard.sh:44, 05_calibrate.sh:58) and "
     "kodama_merge.R (30_merge.sh:38). Outside prose: docs/XENIUM_PIPELINE.md:351 "
     "(50_assemble.sh) and tests/test_verify_run.py:8."),
    (f"{HPG}/submit_all.sh",
     "In-tree pointers only: setup_envs.sh:147, sync_code.sh:11,58, fetch_results.sh:15, "
     "check_shard_equivalence.R:12, slurm/20_shard.sh:15. It is the submitter of the seven "
     "slurm scripts (:104-131) and runs preflight.py (:86) and verify_run.py (:141)."),
    (f"{HPG}/env.sh",
     "Sourced in-tree only: slurm/00_gate.sh:17 (and the six other slurm scripts), "
     "sync_code.sh:22; :20 exports XENIUM_SITES=sites.json."),
    (f"{HPG}/setup_envs.sh",
     "In-tree: sync_code.sh:53, its own :147 next-step echo, and comments in "
     "kodama_merge.R (:24 of setup_envs.sh names it). Outside prose: "
     "docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md:213 (the pinned KODAMA installer) and "
     "docs/EXTERNAL_PACKAGES.md:16. DUP BASENAME with analysis/xenium/hpc/setup_envs.sh."),
    (f"{HPG}/sites.py",
     "In-tree imports only: preflight.py:31 and stage_prep.py:35 (`import sites`); "
     "sites.json:6 names it in prose. :25 reads sites.json."),
    (f"{HPG}/sites.json",
     "Read in-tree by sites.py:25 and exported by env.sh:20; named in preflight.py:50 help "
     "text."),
    (f"{HPG}/preflight.py",
     "In-tree invocation only: submit_all.sh:86. DUP BASENAME with "
     "analysis/follicle/scripts/preflight.py, which owns the follicle docs' hits."),
    (f"{HPG}/fit_calibration.py",
     "One in-tree reference, an echoed operator instruction: slurm/05_calibrate.sh:66. No "
     "import, no automated call."),
    (f"{HPG}/kodama_merge.R",
     "In-tree Rscript call slurm/30_merge.sh:38; comments in check_shard_equivalence.R, "
     "precompute_knn.R:8, setup_envs.sh:24 and verify_run.py:8."),
    (f"{HPG}/kodama_shard.R",
     "In-tree Rscript calls slurm/20_shard.sh:44 and slurm/05_calibrate.sh:58; cited from the "
     "live typing lane as the source of a gate in a path-qualified comment, "
     "typing/cluster_kodama.py:76 (prose, not a call)."),
    (f"{HPG}/precompute_knn.R",
     "One in-tree Rscript call: slurm/10_prep.sh:43."),
    (f"{HPG}/check_shard_equivalence.R",
     "In-tree Rscript call slurm/00_gate.sh:29; comments kodama_merge.R:13, kodama_shard.R:14, "
     "setup_envs.sh:27. Outside prose: docs/XENIUM_PIPELINE.md:234."),
    (f"{HPG}/stage_call.py",
     "One in-tree invocation: slurm/40_call.sh:34. It imports the declared libraries "
     "run_donor:40, xenio:41 and nullgate_fast:42, and stage_prep:43."),
    (f"{HPG}/stage_prep.py",
     "In-tree: slurm/10_prep.sh:37 invocation and stage_call.py:43 (`from stage_prep import "
     "load_at_site`). It imports sites:35, xenio:36 and run_donor:37."),
]

#: The candidate paths, for `retire.py`'s tests: on no tier, in no refusal set.
CODE_CANDIDATE_PATHS = [p for p, _ in CODE_CANDIDATE]


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
       evidence="No producer and no consumer. build_workbook.py:234 writes "
                "master_results_run10.xlsx, not this name, and grep over every .py/.R/.md/.yml/.sh "
                "in the tree finds the string 'master_results_run20' nowhere. Written 2026-09-11 "
                "17:15, the same minute as the file below, by something no longer in the tree."),
    _c(f"{REPO}/data/development/typing/app/islet_celltype_proportions_run20.csv", 2,
       check="no_reference_anywhere",
       superseded_by=f"{REPO}/data/development/typing/app/islet_celltype_proportions.csv",
       evidence="No producer and no consumer. build_workbook.py:238 writes "
                "islet_celltype_proportions.csv (which exists, 2026-09-17 18:00, and is what "
                "trajectory_islet.py:139 reads); nothing names the _run20 spelling. Written "
                "2026-09-11 17:15."),
]
