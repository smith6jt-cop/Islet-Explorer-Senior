# HiPerGator Xenium source locator

`xenium_sources_20261026.csv` gives one exact HiPerGator directory for each of the
65 Xenium region acquisitions in the accepted 2026-10-25 Stage 00 inventory.
It is a **location index**, keyed by that inventory's `acquisition_id` and
output name. It does not define donors, tissues, region relationships or typing.
Stage 04 must bind this file as an input in its receipt and refuse a changed or
missing source; it must not search the filesystem for a similar output.

The 35 outputs incomplete on the workstation use the exact `/orange/maigan`
paths recorded in `/data/xenium_records/MANIFEST.tsv` (SHA256
`bcad804d4d04d7caf60d063b54a22dc67d4546f290be9706465cfe5ad952aee5`).
The other 30 use the top-level output directories under
`/orange/maigan/Xenium_T1D_PancPLN/<run>/`. Two exceptions are explicit:

- Local pLN outputs for serials `0059865` and `0059907` are named `PanT` on
  HiPerGator. Each selected directory's three instrument records match the
  accepted pLN output hashes.
- Serial `0041116` (`Region_1`) and serial `0041117` (`PanT`) each have a second
  MM035 directory inside the run's timestamped subdirectory. This index uses
  the top-level directory. The seven files Stage 04 requires in each pair
  compared byte-identically on 2026-10-26.

Read-only HiPerGator checks on 2026-10-26 found seven nonempty files and the
accepted Stage 00 hashes for the three instrument records in **all 65** chosen
directories. This establishes source availability, not RNA QC or typing. The
recorded damage to the 6533 pLN region still applies and need not hold up
other regions.

`multiplex5000.xenium_sources.build_roster(..., site_index=...)` checks the
index against an accepted Stage 00 table. With an index, it never falls back
to another copy when the selected path is absent. On HiPerGator every row
must report `ready` before Stage 04 ingests it. Stage 04 remains responsible
for schema, counts, cell keys, QC, scope and its own receipt and gate.

## Spleen same-slide pilot sources

`spleen_same_slide_20261026.yaml` names eight explicit pairs of a PC67
PhenoCycler section and the Xenium region run on that same physical section.
The protein cell export, 99/87 GB source image, and Xenium output have exact
HiPerGator paths. Crop origins are read from the one SHA256-bound
`sections_manifest.csv`, rather than copied into another table. The eight
protein export headers share the recorded SHA256. A read-only check on
2026-10-26 found all seven required Xenium files and both protein sources for
every pair, and matched all eight export headers and section-table identities.

`multiplex5000.same_slide_sources.load_pairs(path, require_sources=True)`
performs that check on the host holding the data. Without `require_sources`,
it exposes only the recorded paths; its `crop_origin_px` remains unset. The
manifest records **physical section pairing**, not a matched cell result.
Stage 00 must admit these acquisitions and their relationship in a spleen
pilot project before Stage 04/06 can run on them. Registration and measured
cell overlap are separate Stage 07/08 work. The SpleenFollicleCounterQP repo
and its data remain unchanged.

The tracked `projects/spleen_pilot.yaml` takes this source manifest through the
same Stage 00 function, receipt and gate as other projects. Its Stage 04 selection
is one representative region (`HDL090`, 1/8), stated in the project file. On
HiPerGator, set `MPLX5K_RELEASES_ROOT` to a writable `/blue/maigan` project
folder, then run `mplx5k stage run 00_inventory --project spleen_pilot`.
Review its `gate.html` and accept the gate before Stage 04. Neither stage
declares a cell-type result; Stage 05 and 06 enter this project when built.
