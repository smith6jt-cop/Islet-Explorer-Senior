# R-lane import provenance

## LymphVizkit

- Source commit: `820bb382280feb9f4af0231379943afd3ca211d0`
- Verified backup bundle:
  `/home/smith6jt/repo-backups/2026-08-13/LymphVizkit.bundle`
- Bundle SHA-256:
  `dbb2d0053c0eff3fd3a91a0a37e6222eb8d66f9644bb8ba70d230199d84975ef`
- Original `00_domain_config.R` SHA-256:
  `237418b2d67025974b2d2a66bca6a1726fb8a5ed0ec815597f892b9b81321a96`
- Original and imported `lymphoid_default.yml` SHA-256:
  `6cc842f8558fd6292295119a7e42af36cfcc7e39561fe4f0fbb00dcb91c61f7e`
- Imported synthetic-data manifest SHA-256 (sorted `sha256sum` records):
  `fa5d57388287a5c0bda905c5efb6e21b7bc5cc58fcbb410bcaaaa0b70cff78d2`

The generic loader/accessors in `R/domain_config.R` are an adaptation of
LymphVizkit's proven domain layer. Its lymph YAML and synthetic fixture are
byte-preserved. The LymphVizkit app fork itself was not copied.

## Islet Explorer (Junior)

- Publication branch: `main` at
  `f2a5c540ab95584687b19683f2c23b8bd952beaa`
- Clean local review-hardening branch observed 2026-09-04:
  `codex/ufrc-risk-review-20260902` at
  `e58996168198187840636901acdaad723317fdb5`
- Modernist source tip:
  `f4e5b2225c5b50efd1ebba97477f1e4480ab7475`

The four planned Junior harvests were already present in Senior history at
`259a1d5`, `621e0a4`, `60b470a`, and `661fa27`. The Modernist shell was merged
into Junior at `f827106` and explicitly reverted at `9cf3d50`; the later
publication release stayed on the old layout. D.2 therefore preserves the
current Senior shell and records the Modernist row as superseded. Junior's
later deployment-hardening commit remains explicit D.4 source-reconciliation
work; D.2 does not mutate or claim to replace its publication deployment.

## Recovery point

Before D.2 layout work, the complete monorepo history at `53aa300` was bundled
to `/tmp/Islet-Explorer-Senior-pre-d2-53aa300.bundle` (SHA-256
`ce7cd7e3960c7162774d84bc83f06106f1fccd534a946366b000287926da80b4`).
