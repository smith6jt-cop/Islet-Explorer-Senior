# Docs index — what to read when

The repo's documentation, in reading order:

| Doc | Location | Read it when |
|---|---|---|
| **README.md** | repo root | First contact — install, quick start, the raw→final workflow, Xenium Explorer export. |
| **HANDOFF.md** | repo root | **Resuming work.** Per-session pickup note; the top section is the current state + how to run. |
| **CLAUDE.md** | repo root | The durable architectural guide + the "do-not-regress" decisions (#1–#17). Pancreas-centric. |
| **DATA_README.md** | repo root | Data layout, h5ad layer/obs/obsm contracts, the headline CSV outputs. |
| **THYMUS_PIPELINE.md** | `docs/` | Thymus-arm specifics (samples 0041435/0041534): scVI params, zonation. *Authoritative for thymus parameters; its driver/post-proc commands are superseded by `run_pipeline.sh` — see its top banner.* |
| **HANDOFF_archive.md** | `docs/` | Historical/forensic only — pre-2026-05-27 session notes. Not loaded on session start. |
| **anchor_markers_for_review.md** | `docs/` | The curated broad/fine anchor markers behind phenotyping. |
| **upstream_panel_repair_spec.md** | `docs/` | Historical — the spec handed to XeniumPanelExplorer for the identity/subclass/state panel classification (now vendored). |
| **installed_freeze.yml** | `docs/` | A point-in-time `conda env export` snapshot (reference only; the live spec is `environment.yml`). |
| **scripts/benchmark/README.md** | `scripts/benchmark/` | The self-contained spatial-method benchmark sub-project (BANKSY/UTAG/KODAMA/MISTy). |

**The one command to run the whole thing:**
`STAGES="01 02 phenotype 03 final" TISSUE=<pancreas|thymus> bash scripts/run_pipeline.sh`
