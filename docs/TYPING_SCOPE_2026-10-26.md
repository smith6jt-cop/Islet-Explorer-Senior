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
