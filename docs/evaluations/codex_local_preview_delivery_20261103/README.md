# Final local preview delivery and retention

Received and read the final [runbook](/data/panc_CODEX/Islet-Explorer/docs/local_duckdb_review_20261103.md) and [runtime receipt](/data/panc_CODEX/Islet-Explorer/docs/evidence/local_duckdb_review_20261103.json), delivered through operational commit `43eb115` and documentation commit `e3f4fcc23bd5e504c88e244cead48342d4d8a2a7`. The receipt SHA256 matches `736f955c62009f3a3df2482c4dc24c29fb79baaa22b7f044e31fffde5e0cadc4`. This is a delivery acknowledgement and retention note, not a new review.

The supplied preview URL is **http://127.0.0.1:3845/islet-explorer/**. Application source `7ac421e` and accepted data producer `5e6580f` retain separate identities. The sender reports completed local security, decoder, census, font and UI checks, with public dev and production unchanged. Those runtime results were not re-executed or independently audited for this acknowledgement. The earlier [static source receipt](../codex_local_display_fixes_receipt_20261103/README.md) remains unchanged.

Preserve these active backing resources while the preview is used or may be relaunched; **they are not generic cleanup candidates**:

- Data: `/data/panc_CODEX/Islet-Explorer/tmp/duckdb-dev-ui-release-20261103/final-5e6580f/artifacts/release-5e6580f/data`
- Runtime, launchers and evidence: `/data/panc_CODEX/Islet-Explorer/tmp/duckdb-local-review-20261103/r2`
- Image evidence: `/data/panc_CODEX/Islet-Explorer/tmp/duckdb-local-review-20261103/r2-image`
- Immutable image: `sha256:46b81d665cacb0f8abcd4b9fe34d7f8c786cb7916d6f9cd30e63a1223358b006`

The launcher is `r2/start-local-service.sh`; the recorded service is `islet-explorer-local-review-20261103-r2.service`. The runbook says it survives shell/chat exit but is not enabled across reboot. No stop, relaunch, cleanup, port, build or public-service action was taken by the recipient. Local-only direction remains in force.
