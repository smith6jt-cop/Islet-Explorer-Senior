# Spleen large-input hash verification

All ten previously unverified large inputs were read in full: **207,864,852,220
bytes**. Eight cell-geometry exports match their historical full SHA-256 hashes.
The two QPTIFF images match the historical **sampled** fingerprints declared in
the exact receipt-bound QuPath manifest. New full SHA-256 hashes are recorded
separately for both images.

| Input group | Historical algorithm | Current verification |
|---|---|---|
| Eight cell GeoJSON exports | Full-file SHA-256 | Eight exact matches |
| `0110979_SpleenPCXen_Scan1.er.qptiff` | Three deterministic 1 MiB windows | Exact sampled match; new full hash `49e0eacb15d8468df8f59586bd2029b3cc88877b2894eec745c85e195acf87db` |
| `0110618_SpleenPCXen_Scan1.er.qptiff` | Three deterministic 1 MiB windows | Exact sampled match; new full hash `a815b2d31b9c5a0dcc6e4cc0f18684df37a49de4fc2b7fcb7e2f3a6fd1bd4889` |

The manifest SHA-256 is
`f19332192f1ff43aba2685d1507af341767a983c5a99465a03b3ed66572db324`,
matching the existing Stage 06 receipt. Each source file retained its size and
modification time during the full-byte read. The
[final report](../../../outputs/blocker_remediation_20261102/hpg_maldi_xenium/spleen_algorithm_verification/verification.json)
records the expected algorithm, observed comparison digest, full-file digest,
sizes, source paths and producer identity. No original receipt was rewritten.

The first attempt compared full hashes with all `recorded_sha256` values and
reported two mismatches. That comparison was inapplicable to the two sampled
fingerprints; it was not evidence of changed image content. Its log/report/source
remain preserved under `spleen_byte_verification`. The corrected checker requires
the hash-mode manifest to match its receipt, tests sampled and full algorithms
separately, rejects a changed manifest or a digest inconsistent with its declared
full/sampled binding, and fails on missing algorithm declarations.

These checks close the missing current-byte verification. Historical sampled
fingerprints still cannot establish historical full-byte equality. They also do
not establish biological applicability of PC67 proposals, a PC33 registry, or an
independent spleen RNA hierarchy. Those scientific requirements remain separate.

Reproduction:

```bash
PYTHONPATH=core/multiplex5000 .venv/bin/python scripts/evaluation/verify_recorded_inputs.py \
  --receipt /absolute/path/to/receipt.json \
  --manifest /absolute/path/to/receipt-bound-qupath-manifest.json \
  --out /absolute/path/to/new-verification-directory
```

Run this where the receipt's exact source paths are mounted. The output directory
must be new; sources are read only. No path or donor-name substitution is used.
