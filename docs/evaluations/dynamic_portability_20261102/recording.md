# Project-only KINTSUGI observations

The only production-code addition in this evaluation is
`multiplex5000.parameter_observations`. It loads KINTSUGI's actual
`ParameterLearningEngine`, `ParameterRecord`, SQLite schema and insertion method
at commit `716d1ebd87db916ab04910c6374e04d2d9ea508e`. Both Git HEAD and the
engine file SHA256 are checked before opening storage. The engine hash is
`9c29be884f5dd9d405033a9224b998778fe3905f99093036849a0b9436be6db5`.

The upstream MCP success endpoint defaults to global history and user approval,
and its quality defaults cannot represent missing measurements. The evaluation
MCP adapter therefore wraps the pinned engine with three restricted operations:

- `record_parameter_observation`: proposed, evaluated, failed or blocked only.
- `get_parameter_observation_history`: inspect all revisions or one observation.
- `export_parameter_observations`: export JSON inside the selected project.

An actual stdio MCP client/server handshake exercised all three operations,
including duplicate retries, revision history and rejection of approval. See
`outputs/dynamic_portability_20261102/recording-protocol/protocol_check.json`.
This is a callable project-scoped server; no global MCP configuration was changed.
The populated scientific project was also read, retried idempotently and exported
through MCP; its export matches the CLI export byte-for-byte. See
`parameter-history/scientific_mcp_check.json` in the evaluation root.

## Setup and connection

Use a separate environment for the optional MCP dependency. The scientific
interpreter does not require MCP to import or use the recorder through its CLI.
The `parameter-recording` package extra pins `mcp==1.30.0`.

```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git external/KINTSUGI
git -C external/KINTSUGI checkout 716d1ebd87db916ab04910c6374e04d2d9ea508e
python3 -m venv --system-site-packages external/.venv-parameter-recording
external/.venv-parameter-recording/bin/pip install 'mcp==1.30.0'
```

The checkout already exists in this workspace. Use the actual repository URL
reported by `git -C external/KINTSUGI remote get-url origin` if provisioning
elsewhere. Launch from the repository root, replacing `/absolute/repo` below:

```bash
PYTHONPATH=core/multiplex5000 external/.venv-parameter-recording/bin/python \
  -m multiplex5000.parameter_observations \
  --project /absolute/repo/outputs/dynamic_portability_20261102/parameter-history \
  --kintsugi /absolute/repo/external/KINTSUGI serve
```

An MCP client uses that executable, argument list and `PYTHONPATH` environment.
The project is fixed at process launch. Clients cannot redirect storage. Replace
`serve` with `history`, `export`, `record --input observation.json`, or `retry`
for CLI inspection, recording or durable-outbox replay. CLI failures that retain
a retryable observation exit with status 2. Invalid input fails explicitly.
Queued revisions replay in observation/revision order; later revisions remain
queued if their predecessor still fails. Immutable-content conflicts require a
corrected revision and are not reported as transient storage failures.

## Record contract

The schema is `multiplex5000-parameter-observation/1`. Each record includes an
observation ID, monotonically appended revision, preceding revision, status,
operation, compatibility context, typed subject, settings, quality, outcome and
hashed evidence references. Context explicitly includes modality, assay, tissue,
panel, relationship, acquisition/section/segmentation IDs, measurement units,
producer version, registry/panel hashes and input receipts. Unknown values are
`null` or empty identifier/hash collections with an explanation in the outcome.

Subjects can be markers, genes, feature sets, acquisition pairs, acquisitions or
geometry. Marker-independent work stores SQL NULL in the upstream marker fields.
Settings distinguish `fixed`, `selected_method`, `fitted` and `review_criterion`.
The original value units are retained; diagnostic counts are observations, not
quality scores. Missing quality and image characteristics remain NULL. Detailed
tables and images stay in versioned artifacts referenced by path, ID and SHA256.

The populated evaluation database is
`outputs/dynamic_portability_20261102/parameter-history/.claude_workflow/parameter_learning.db`.
The tracked [JSON export](parameter_history.json) is immediately inspectable
without Python, MCP or KINTSUGI. Synthetic protocol records live in a separate
project and are not included in the scientific observation export.
The final export contains 351 observation IDs and 499 immutable revisions:
41 proposed, 150 evaluated, 46 failed and 114 blocked latest observations.
All remain unapproved. Exact source-to-port CODEX comparisons and failed original
artifact assertions are separate operations, so one does not erase the other.

For this dated delivery, source audits run before
`scripts/evaluation/finalize_portability_followup.py`, which freezes their reviewed
copies and compressed traces. `record_portability_observations.py` then appends
only changed observations and verifies duplicate retries. Run the scientific MCP
check, rerun the finalizer to refresh summary counts, then generate the artifact
manifest and run the delivery verifier. Do not rerun the preliminary
`assemble_portability_evidence.py` over this finalized bundle: it predates the
follow-up. A changed audit needs a new versioned artifact; previously referenced
evidence must retain its bytes and hash.

## Isolation, durability and acceptance

Storage requires an explicit absolute project root and `project_only` mode.
The original engine receives `global_db=False`. Storage and export reject path
escape, including preexisting symlinks. The supplied project directory must be
trusted; this is not a sandbox against a concurrent hostile filesystem actor.

An atomic JSON outbox write precedes database work. If recording fails, science
outputs remain untouched and the response identifies the retryable observation.
If the outbox itself is unwritable, `pending_path` is null and the caller must
retain its submitted input. Canonical payload hashing and a unique SQLite index
make identical ID/revision retries idempotent. A changed payload requires a new
revision naming its predecessor. A POSIX file lock serializes local writers;
SQLite and file-lock behavior on other filesystems must be verified before use.
History conflicts fail explicitly and are not queued as temporary failures.

SQL triggers prevent alteration, deletion or approval of observation rows.
All observations are unapproved; the upstream approved-only export excludes
them. The MCP surface offers no recommendation, application or approval tool.
Scientific routes do not import the recorder: a recording failure cannot change
scientific outputs. Existing successful, user-approved records remain a separate
upstream category, never inferred from `evaluated` status.

Future reuse requires a separately reviewed compatibility policy covering at
least modality, tissue, panel, units, method version and relationship type.
Unknown compatibility fields cannot establish a reusable match.
This evaluation implements no automatic reuse, recommendation-driven tuning,
registry approval, production deployment or release activation.

Tests exercise the actual pinned SQLite engine: restart/export round trips,
concurrent idempotent retries, revision preservation, null measurements,
project/global isolation, path restrictions, approval denial and durable retry.
