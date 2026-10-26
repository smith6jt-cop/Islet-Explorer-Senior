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

Stage 06 must produce tissue-specific calls, evidence and unavailable reasons
for both modalities. A same-slide section needs both assay results; the
cell-to-cell relationship is established only from the recorded section mode
and registration, then used in Stage 08 integration. Serial sections never
receive a same-cell claim. The damaged 6533 pLN Xenium region may be marked
unavailable with its damage reason and must not hold up other acquisitions.

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
