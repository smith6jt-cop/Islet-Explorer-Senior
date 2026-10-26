# Stage 06 typing scope: both assays and same-slide acquisitions

Maintainer clarification, 2026-10-26: tissue-specific typing covers PhenoCycler-only,
Xenium-only and PhenoCycler–Xenium same-slide acquisitions. A protein-only Stage 06
is incomplete. This contract applies to the first pancreas and pancreatic lymph
node release.

## Inputs that Stage 06 must consume

| Input | Producer / authority | Role |
|---|---|---|
| Protein measurements and availability | Stage 03 | Measured protein evidence; absence is not zero |
| Protein normalization | Stage 05 | Donor-consistent marker units and the CD8 shift |
| Xenium cells and RNA QC | Stage 04 (must be built before Stage 06) | RNA evidence on Xenium-only, serial and same-slide sections |
| Tissue, ROI, panel, condition and relationships | Accepted Stage 00 inventory | Route each acquisition without inferring a pair from co-imaging |
| pLN protein phenotype rules | `core/panels/data/tissues/pancreatic_lymph_node/manifest.yml` | Selected full sheet and reviewed pLN adaptation |

The older migration schedule assigned RNA ingest to D8, after D6 typing.
Stage 04 is now a D6 prerequisite: its accepted RNA inputs are needed before
Stage 06 can claim the clarified all-modality typing scope.

Stage 03's measured roster is the 25 pancreas
PhenoCycler images in `cohort-spec-v3.json`. Its accepted receipt does not
cover the 19 dedicated PhenoCycler pLN images. The next measurement increment
must admit those images under the same compartment and exclusion policy before
Stage 06 can call their cells; it must not silently read an older pLN output.

The 2026-10-25 accepted inventory lists **28 donors**: 22 sequential, four
same-slide and two Xenium-only (6512 and 6618). It has 25 PhenoCycler
pancreas images, 19 dedicated PhenoCycler pLN images, 32 Xenium pancreas
regions and 33 Xenium pLN regions. These are acquisitions, not cell counts.
Embedded pLN/lymphoid on a pancreas image must be routed by its recorded
same-block or different-block ROI designation; it is not another acquisition.
The PhenoCycler exclusion for 6457 does not exclude its Xenium data.

The current `io60_v3` cohort spec lists only the 25 PhenoCycler image donors;
it cannot be Stage 04 or Stage 06's sole roster. Stage 00 also records the
Xenium-only 6512/6618 outputs and the four same-slide runs at their
HiPerGator `/orange/maigan` source locations. Those complete outputs are not
on this workstation. Stage 04 must use the accepted acquisition inventory
and resolve those sources on HiPerGator or from a verified transfer. It must
not substitute an older serial Xenium section for a same-slide output.
The current pancreas RNA release worker fixes its roster at 25 sections, and
the pLN v3 RNA worker discovers `/data/xenium/h5ad_out`; that local folder has
neither the two Xenium-only donors nor the four newer same-slide serials.
Reusing either worker's old roster would leave the clarified scope untyped.
The project currently fixes `settings.xenium_root` at local `/data/xenium`;
the relative Xenium bundle paths need a site-specific root on HiPerGator.
The absolute same-slide paths already recorded under `/orange/maigan` must
remain absolute rather than being joined to either root.

`multiplex5000.xenium_sources` reads the accepted Stage 00 acquisition table,
its receipt, and the checksum-bound Xenium records manifest. It uses
`serial__roi` keys, checks complete bundles and their instrument-record hashes,
and reports exact HiPerGator paths for outputs not complete on this host. The
2026-10-25 inventory resolves to 30 ready local region bundles and 35
`remote_only` bundles. `remote_only` means the HiPerGator path is recorded,
not that its large cell files have been verified there; Stage 04 must check
those files on that host before ingest. The tool is a readiness check, not
an RNA producer or a typing result.

`multiplex5000.xenium_region.load_region` is the first Stage 04 reader. It
loads one ready raw bundle's Gene Expression counts and joins its matrix
barcodes to `cells.parquet` by native ID, failing on missing or repeated IDs.
It emits acquisition-qualified `cell_uid` values and keeps cells without a
nucleus in the accounting while marking them ineligible for typing. It makes
no cell-type call or gene-positive cutoff. On the local 6539 pilot, it read
319,005 pancreas cells (41,814 no nucleus) and 222,443 pLN cells (6,344 no
nucleus), both with 5,101 genes. This pilot did not publish Stage 04 or run on
the 35 HPG-only bundles.

The Stage 04 runner declares the location record, accepted Stage 00 evidence,
and each selected region's raw matrix, cell table and instrument records as inputs. It
publishes one region-qualified cell-QC Parquet and feature list per acquisition,
plus a gate figure and receipt; the raw matrices stay at their indexed sources.
Its synthetic same-slide test passes with colliding native IDs, a zero-nucleus
cell, a missing nucleus count, and a damaged pLN region. A project can name
`settings.rna_acquisitions` as an explicit accepted subset for a pilot; the
gate reports its fraction of the full inventory and never claims full coverage.
With no subset it requires all accepted Xenium regions. On the workstation it
can use the accepted local paths and records manifest; on HiPerGator it uses
the complete exact site index. A 65-region HPG run
still needs the accepted Stage 00 evidence staged there and a site file naming
`xenium_source_index`; no full-cohort Stage 04 gate has been accepted.

A read-only HiPerGator check on 2026-10-26 found all seven required files and
matching Stage 00 instrument hashes for all 65 outputs. The 30 outputs that
are complete locally also exist under `/orange/maigan/Xenium_T1D_PancPLN`,
but their exact HiPerGator paths are not all in `MANIFEST.tsv`. The exact
65-region locator is now
`core/multiplex5000/sites/hipergator/xenium_sources_20261026.csv`; Stage 04
must bind it in its own receipt. On HiPerGator all 65 indexed paths passed
the source-file and instrument-hash check. In that tree,
serials 0059865 and 0059907 have pLN outputs named `PanT` rather than the
local `pLN`; the accepted instrument hashes match. Serials 0041116 and
0041117 each have two directory copies whose seven required files compare
byte-identically. The site index should choose one path explicitly, never
pick the first filename match found by a search.

Stage 06 must produce tissue-specific calls, evidence and unavailable reasons
for both modalities. A same-slide section needs both assay results; the
cell-to-cell relationship is established only from the recorded section mode
and registration, then used in Stage 08 integration. Serial sections never
receive a same-cell claim. The damaged 6533 pLN Xenium region may be marked
unavailable with its damage reason and must not hold up other acquisitions.
The current `section_modes.csv` records the four same-slide pancreas serials,
but has no pLN rows for their `pln_2` regions. That gap does not prevent typing
each assay's pLN cells; it prevents a same-cell pLN link until its relationship
is explicitly recorded and checked.

The selected 69-row pLN protein sheet was written for a wider pancreas-plus-pLN
panel. The maintainer approved an explicit pLN adaptation of pancreas-only
negative parent constraints. On the dedicated 55-marker pLN panel, 14 rows
still require unmeasured markers; CD38-dependent B-cell subtype rows remain
unavailable. These protein rules do not type Xenium RNA cells. Stage 06 must
use the tissue-specific RNA producer and map its calls into the release's
class vocabulary. Where the hierarchical protein method writes both
`cell_class` and `refined_class`, only the final step is served.

Stage 06 is the first manual processing review stop. It must show per-tissue,
per-assay and per-acquisition counts, marker or gene support, controls and
unavailable calls. The `NO_VERIFIED_TYPING` latch stays down until the
maintainer accepts processing. No H1–H6 inference runs before that decision.

## Spleen in the first portable pilot (maintainer, 2026-10-26)

The first **workflow pilot** also includes executable spleen same-slide typing;
this does not add spleen to the first active IO60 release. The spleen protein
reference is `core/panels/data/tissues/spleen/manifest.yml`. It derives from
the selected full 6551 hierarchy, removes pancreas-only branches and object
labels, retains applicable lymphoid rules, and adds descriptive CD21-positive
B, CD45RO-positive T, CD141-positive cDC and SMA-positive stromal leaves. The
PC67 inventory is the exact, SHA256-checked header from all eight export files.
Of the 58 retained rules, 55 use measured PC67 markers; Neural, Schwann cell
and Macrophage require B3TUBB or Iba1 and stay unavailable. These counts are
rule availability, **not typed-cell counts or validated positivity cutoffs**.
The HPG UniFORM scale table, SHA256 `7927e33b…`, covers only 33 classifier
channels per PC67 section, not the full 67 measured channels. Only five rules
have all inherited markers in that scaled set (Proliferating, PD-L1 expressing,
MHC I high, MHC II high and CD44 high); none establishes a T-cell, B-cell or
FRC lineage. Stage 05 must derive and check the remaining required marker
scales, or mark those lineages unavailable. Raw means in `cells_full` cannot
be substituted for normalized marker evidence. The existing scale table is
comparison evidence, not the completed Stage 05 producer.

The paired protein images and crop-local cell exports are under
`/blue/maigan/smith6jt/SpleenFollicleCounterQP/` on HiPerGator; the eight
Xenium regions are under
`/orange/maigan/20260826__224549__SpleenSH2B3postPheno8sample/`.
The project is an active read-only input for this pilot and is not modified.
Each Xenium region and its corresponding PC67 section must enter the accepted
inventory with a recorded same-slide relationship and coordinate-frame offset.
The two assays receive their own Stage 06 evidence before Stage 08 links cells;
the matching donor code alone never establishes a cell pair. The older spleen
integration plan said RNA cell typing would not run. The maintainer's newer
all-modality typing direction supersedes that scope statement.

The spleen pilot project now routes its selected PC67 full-cell QuPath export
through Stage 03 before the selected Xenium region enters Stage 04. This Stage
03 adapter preserves `object_id`, uses an acquisition- and export-qualified
`cell_uid`, records all cell means as raw measurements, and excludes zero or
missing nucleus area from protein typing eligibility. It verifies the export
against the accepted Stage 00 digest and labels its gate as raw measurement
evidence. It does not claim the pancreas REDSEA path or normalized protein
values for spleen. The first bounded pilot uses HDL090's one physical section
on both assays; the other seven recorded sections remain in the inventory.

`multiplex5000.tissues.evaluate.evaluate_rules` is the shared Boolean
hierarchy evaluator for the selected pLN sheet and its spleen adaptation. It
takes explicit per-cell positive, negative or unavailable marker calls; it
does not derive positivity thresholds from intensities or publish a final
cell class. Missing markers in any ancestor make a rule unavailable. Stage 05
and the tissue-specific Stage 06 producers still have to supply and validate
those marker calls before this can type a real cell.
