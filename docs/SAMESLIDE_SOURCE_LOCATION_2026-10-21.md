# Four-slide same-slide sources on HiPerGator

Verified by read-only SSH on 2026-10-21. The tracked
[`sameslide_four.json`](../core/multiplex5000/config/sameslide_four.json) is the
section-by-section locator, not a template for putting all files under one
checkout. Its SHA-256 is `58f4b937305b62b974395cee31bf9632b1a93dbbd0e635e710c223c12ec30d28`
both here and at
`/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler/config/sameslide_four.json`
on HiPerGator. The cluster checkout is still at `37f9311` and uses the old
`core/phenocycler` name; its results have not moved with the local package rename.

| Payload | HiPerGator location |
| --- | --- |
| 112 and 6533 same-slide Xenium outputs | `/orange/maigan/20260910__211837__6533_HDL112_pLNpanc_PhenoXeno/` |
| 6347 and 6580 same-slide Xenium outputs | `/orange/maigan/20261002__185939__PhenoXeno_pLNpanc_6580_6347/` |
| Xenium-only comparison sections | The per-section bundles named in the manifest, principally under `/orange/maigan/Xenium_T1D_PancPLN/` |
| Four-slide PhenoCycler qptiffs | `/blue/maigan/smith6jt/Xenium_Analysis/data/phenocycler/` |
| Protein measurements and geometry | Per-section paths in the manifest: `/blue/maigan/smith6jt/IO60panc2nd_v3/results/v3/export/` or `/blue/maigan/smith6jt/io60_sameslide_inputs_20261019/segmentation/` |
| Delivered annotations, islets and cell calls | `/blue/maigan/smith6jt/io60_sameslide_inputs_20261019/` |
| Lane outputs, figures, and `current.json` receipts | `/blue/maigan/smith6jt/Islet-Explorer-Senior-main/data/development/sameslide_four/` |

Each of the eight same-slide region bundles (pancreas and pLN for each of four
donors) exists at the manifest's exact path and contains `experiment.xenium`,
`gene_panel.json`, `metrics_summary.csv`, and `cells.parquet`. In particular,
the four region bundles for 6347 (`0080579`) and 6580 (`0080502`) are present
under the **October 2 run directory above**. Searching only
`/orange/maigan/Xenium_T1D_PancPLN/` and
`/blue/maigan/smith6jt/Xenium_Analysis/data/` misses them. The two October 2
pancreas bundle names are:

- `output-XETG00298__0080579__panc__20261002__190021`
- `output-XETG00298__0080502__panc__20261002__190021`

Each has a sibling `__pLN__` bundle. These are source acquisitions, not
unavailable data. By contrast, the *analysis* for 6533's pLN remains
scientifically unavailable under the registration gate. The lane's
`completion/current.json` reports `validated_complete` for the completion
task, with 115 validated items, 6 scientifically unavailable items, and 1 not
run. That receipt does not turn the 6533 pLN result into a passing result.

HiPerGator home is limited to 40 GB, and [UF Research Computing's storage
guidance](https://docs.rc.ufl.edu/quickstart/practical_storage/) directs job
I/O away from `/home`; the split across `/blue` and `/orange` is deliberate.
The workstation's `/data/xenium` is not a complete mirror, and a local path
check cannot establish that a HiPerGator acquisition is absent.

## Required Stage 00 correction

`xenium_paths.csv` currently lists only 112 and 6533's same-slide **pancreas**
outputs (via `/blue` copies). Their same-slide pLN outputs are absent from that
table, as are both regions for 6347 and 6580: **six region rows are missing**.
The table does contain separate 2024 pLN sections for all four donors; those
must not be mistaken for the 2026 same-slide pLN regions. PR #38 describes
6347/6580 as having no HiPerGator output; that statement is false. Before
accepting Stage 00, reconcile all eight same-slide region bundles into its
acquisition inventory with exact serial, tissue and source path. Check whether
other `xenium_paths.csv` readers assume one row per donor and ROI before adding
the six rows there; a dedicated acquisition inventory may be necessary. If
Stage 00 runs on the workstation, copy only the small instrument-record files
it reads, record their HiPerGator source paths and SHA-256 values, and check the local copies
against those values. Keep `remote source exists` distinct from `full bundle
present locally`; regenerate run, slide and coverage counts and figures after
the six rows are admitted.

The lane's old `registration_candidates` paths are historical optional
candidates; four no longer exist. They are not the current registration
receipts and are not evidence that the Xenium bundles are missing. Existing
HiPerGator results, their inputs and their receipts must remain in place
while consolidation accounts for them; no bulk copy, relocation or rerun is
needed to correct the Stage 00 record.
