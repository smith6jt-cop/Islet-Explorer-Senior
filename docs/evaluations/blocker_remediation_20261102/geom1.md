# GEOM-1: independent same-slide batch coordination

The new [batch coordinator](../../../core/multiplex5000/multiplex5000/same_slide_batch.py)
executes the existing Stage 07/08 route independently for every explicitly
configured pair. It preserves each pair's output root, receipts and review gates,
and records the entire declared roster in one ledger. One missing transform or
failed input contract does not prevent another ready pair from progressing.

The [evaluation configuration](../../../core/multiplex5000/projects/same_slide_batch_evaluation.json)
binds the frozen 12-pair roster by SHA256. Its
[eight spleen project files](../../../core/multiplex5000/projects/same_slide_evaluation)
name exact source-manifest entries and RNA acquisitions. These are pending
evaluation configurations: they supply no assay output, transform, hierarchy map
or acceptance. The four pancreas declarations have no implemented section-keyed
adapter, so their project entries remain null. The roster's unresolved pLN
relationships are retained separately; they are not invented pairs.

## Run and inspect

From the repository root:

```bash
PYTHONPATH=core/multiplex5000 .venv/bin/python -m multiplex5000.same_slide_batch \
  preflight \
  --config core/multiplex5000/projects/same_slide_batch_evaluation.json \
  --out outputs/same_slide_preflight_new
```

Every `--out` must be a new directory. It contains `ledger.json` and an artifact
manifest. The preflight reads configured sources and receipts, writes only this
new ledger, and never creates or modifies project stages. Exit code **1** means
the declared cohort is incomplete, with per-pair requirements in the ledger;
**0** means every declared pair completed its requested candidate computation or
has an existing exclusion, with no unresolved relationship declarations. Exit
code **2** indicates an invalid configuration or command. None grants scientific
approval.

Once the exact upstream contracts exist in the configured isolated pair projects,
use the same command with action `registration` to prepare each eligible pair's
Stage 07 figures and receipt. It stops at that pair's existing manual gate. Review
the actual pair before accepting through the ordinary `mplx5k gate` workflow.
The coordinator has no approval action and refuses automatic Stage 06/07/08 gates.

After a pair's actual Stage 07 acceptance, action `integration` runs Stage 08 for
that pair when its required assay evidence is available. Other pairs remain
blocked or unavailable in the ledger. Stage 08 candidate computation and its
final candidate review are reported separately. Existing outputs are reused only
after rechecking all declared input hashes, settings and output hashes; stale
outputs fail closed. The coordinator never invokes `--force` or replaces a stage.

## Exact configuration and evidence modes

The JSON schema is `multiplex5000-same-slide-batch/1` and has four fields:

- `roster`: a file `path` and its `sha256`. The file contains `explicit_pairs`
  and any `unresolved` relationship declarations.
- `workspace_root`: the explicitly selected isolated candidate workspace.
- `pairs`: exactly one entry for every roster `pair_id`, with an explicit
  `project` path or null, and mandatory `required_evidence`.
- `schema`: the schema identifier above.

Config and roster paths may be relative to the config file. Project release
roots must be distinct, nonoverlapping descendants of `workspace_root`, and
ledger directories must be outside them. Repeated physical pairs, duplicate or
missing configurations, changed roster hashes, incompatible source identities,
mixed coordinate frames and redirected stage output paths are refused.

`required_evidence: independent_calls` requires actual independent protein and
RNA identity call columns with at least some nonempty calls, plus an explicit
tissue/hierarchy mapping. Per-cell unresolved states still remain unavailable
and are preserved in the Stage 08 outputs. Gene-only RNA evidence cannot satisfy
this mode. Geometry may be prepared independently while those biological inputs
remain unavailable, but identity integration is not substituted.

`required_evidence: spatial_gene_evidence` explicitly requests the narrower
spatial/gene-evidence route already supported by the pilot. Its ledger reports
identity comparison as unavailable when no mapping exists. The supplied
evaluation config requests `independent_calls` for every pair; it does not
silently downgrade any acquisition to the narrower route.

The coordinator verifies accepted Stage 06 source/cell bindings and call output
hashes before geometry work, and preserves the existing Stage 07/08 receipt
checks. It does not create or copy upstream artifacts or approvals. Building or
reusing upstream inputs still requires their original measurement, cell-key,
segmentation and provenance contracts. No scientific thresholds or Stage 07/08
implementation changed for this addition.

## Full declared-roster preflight

The recorded run is
[`full_cohort_preflight/ledger.json`](../../../outputs/blocker_remediation_20261102/geom1/full_cohort_preflight/ledger.json).
All **12 declared pairs are accounted for**:

| Scope | Pairs | Disposition | Concrete requirement |
| --- | ---: | --- | --- |
| Spleen PC67/Xenium | 8 | Blocked | The exact per-pair accepted Stage 06 receipt and supplied transform paths are absent in the new evaluation projects. Actual RNA identity calls and a declared hierarchy map are also required by the configured target. |
| Pancreas protein/Xenium | 4 | Unavailable | An adapter must establish exact section cell keys, coordinate frames and segmentation provenance for the donor-partitioned IO60 outputs. |
| pLN relationships | Not asserted as pairs | Unresolved | Inventory evidence does not establish their physical same-section relationships. |

These results describe the explicit pending evaluation projects; they are not a
claim that the same artifacts were searched for at every possible site or path.
The coordinator is executable for supported contracts, as demonstrated by two
ready synthetic pairs completing while a third pair remains blocked. It does
not make full-cohort independent spleen RNA biology or an IO60 section adapter
available.

Verification and source hashes are recorded in
[`geom1`](../../../outputs/blocker_remediation_20261102/geom1), including focused
batch tests, the unchanged same-slide regression suite, lint results and the
full-roster ledger. **98 tests passed**: 20 coordinator tests and 78 existing
same-slide/crossmodal tests, with clean lint and whitespace checks. Both the
repository `.venv` CLI and the scientific `scvi-env` CLI produced the same
12-pair preflight dispositions. The explicitly labelled
[`synthetic_batch/ledger`](../../../outputs/blocker_remediation_20261102/geom1/synthetic_batch/ledger)
retains the two-completed/one-blocked demonstration; its test gate acceptances are
not cohort approvals. The earlier
[read-only assessment](geom1_batch_assessment.md) records why simply removing
the singleton selectors would have coupled unrelated pair gates.
