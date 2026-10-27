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
Stage 00 records two PhenoCycler whole-slide images, eight section exports,
two Xenium slide IDs in one run, eight Xenium region output folders, and the
eight physical-section relationships in a spleen pilot project before Stage
04/06 can run on them. Registration and measured
cell overlap are separate Stage 07/08 work. The SpleenFollicleCounterQP repo
and its data remain unchanged.

The tracked `projects/spleen_pilot.yaml` takes this source manifest through the
same Stage 00 function, receipt and gate as other projects. Stage 03 reads raw
QuPath means from all eight PC67 sections, Stage 04 selects one representative
Xenium region (`HDL090`, 1/8), and Stage 05 fits full-panel UniFORM shifts
across the eight protein sections. These stages have review notebooks under
`core/multiplex5000/notebooks/`. Stage 03's pilot adapter does not yet perform
REDSEA, and neither assay has a cell-type result. On HiPerGator, set
`MPLX5K_RELEASES_ROOT` to a writable `/blue/maigan` project folder, and set
the site `tmp_dir` to a sufficiently large writable `/blue/maigan` folder.
The Stage 00 output made before commit `b642305` is invalid: it counted eight
protein section exports and eight Xenium region folders as 16 acquisitions.
After installing main at or after `93e995c`, rerun that stage with
`mplx5k stage run 00_inventory --project spleen_pilot --force`. The runner
moves the old output and receipt aside. Review the new `gate.html`, which must
show two whole-slide images, eight section exports, eight Xenium regions and
one Xenium run. Do not accept the old gate or advance the pilot on it. Stage 06
produces separate review-only protein rule matches and Xenium gene evidence
for HDL090. It
publishes no final cell class or cross-assay cell link. Stage 06's gate waits
for manual review; its notebook is `notebooks/06_typing.ipynb`.

Stages 07 and 08 now use the packaged
`multiplex5000.integration.paired_typing.link_candidate_cells` function.
Stage 07 admits one supplied registration, checks the exact section and assay
identities, and draws a coordinate overlay for manual review. Name its root
as `sources.same_slide_registrations` in the site file. For the selected
section, put `transform.json` and, when the transform declares one,
`transform_displacement.npz` under `<root>/<protein_entry>/`. The transform
metadata must give the exact `protein_entry`, `xenium_acquisition_id`,
`moving_frame: region_local_um`, `fixed_frame: crop_local_um`, and
`reviewed: true`. That flag declares prior review; Stage 07's gate is still
manual. Stage 08 requires that gate's acceptance, verifies its receipt and
transform hash, then writes mutual-nearest centroid **candidates**, distances,
and unmatched counts. It makes no verified same-cell or final typing claim.
No spleen transform has been supplied to this project yet, so neither stage
has produced a real HPG result.

The HPG login shell exposes Python 3.9 by default. Load its Python 3.12
module (`source /etc/profile.d/modules.sh; module load python/3.12`) and make
a **new** virtual environment under a writable `/blue/maigan` pilot folder.
Install this repository's merged checkout from `core/multiplex5000`; keep
Snakemake in a separate environment using `requirements-graph.txt`. The
UniFORM dependency for Stage 05 must be at commit
`7c33396d12e5a19dc0b36e46a94c3293de0cd814` in that checkout's
`external/UniFORM`. Put the site file, release root and task scratch under
the same `/blue/maigan` pilot folder; set `MPLX5K_SITE_FILE` to the site file.
The Slurm command is available at `/opt/slurm/bin/sbatch` if it is absent
from a noninteractive shell's `PATH`. None of this setup needs to write in
`SpleenFollicleCounterQP` or the crowded home directory.
