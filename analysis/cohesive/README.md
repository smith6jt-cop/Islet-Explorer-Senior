# Cohesive analysis lane

One controller that finishes the study from the October 14 analytical handoff: reconcile inputs,
qualify them, freeze the spatial protocol, compute the six spatial hypotheses per donor, test them at
the donor level, and deliver figures, source tables and a deck. It writes only under the output folder
it is given; it never writes to source, typing or app roots, and it never activates production typing.

## The typing it reads

```bash
python -m analysis.cohesive.typing_generation      # from the repository root; writes /data/io60-analysis/typing-20261019/run
```

builds the pancreas protein typing generation and pins every donor's release, admission and candidate
file in `<root>/contract.json`, which `qualification.GENERATION` reads. It is the October 8 chain, run with
the producers recovered into this branch (`recovered_v3_producers_20261018.json`), with two rulings:

- **2026-10-18** — Ker8_18 is a veto in no immune definition, so the four same-slide donors are typed
  like every donor (`subpanel_io.UNSHARED_NEGATIVES`).
- **2026-10-19** — CD8 carries UniFORM's shift on all 25 donors, in the refinement as in the primary
  calls. The October 8 chain divided it back out for the refinement only.

## Run

```bash
make analysis-run \
  ANALYSIS_PY=$PWD/.venv/bin/python \
  ANALYSIS_INPUTS=/data/typing-validation/current_defensible_handoff_20261014_183339/LATEST_DEFENSIBLE_INPUTS.json \
  ANALYSIS_PROTOCOL=<archived spatial_hypothesis_plan_v2.json, sha256 cf648db0…> \
  ANALYSIS_OUTPUT=<a new folder under /data/io60-analysis>
```

`analysis-plan` / `analysis-status` print the stages and their recorded state; `ANALYSIS_THROUGH=<stage>`
stops early. **No hypothesis stage (h1–h6, inference) runs until the data are properly processed**
(maintainer, 2026-10-19: "the technical issues mask results"). The `.venv` must resolve `phenocycler`
to the v3 checkout, because the protein and RNA readers are path-bound and run where they live. The
deck is assembled by `deck.py` in a separate environment (`/data/io60-analysis/env/presentation`, from
`presentation-requirements.txt`) so the shared `.venv` is not altered.

## Stages

| stage | module | writes |
|---|---|---|
| reconcile | `sources.py` | selected inputs, reader contract checks, source snapshot |
| qualify | `qualification.py` | per-donor protein cells (refined hierarchical identity), per-section RNA cells with exact native expression, marker/immune availability |
| protocol | `protocol.py` | frozen protocol with the 2026-10-17 and 2026-10-18 amendments, measurement feasibility |
| atlas | `atlas.py` | coverage, qualification flow, true-Islet composition (outline area, no expansion), protein and RNA T-cell states |
| h1–h4 | `endpoints.py` | protein islet- and donor-level endpoints |
| h5, h6 | `rna_endpoints.py` | RNA organizer motifs (post-IF masks for localization only); pLN T-state architecture and H3 pairing |
| inference | `inference.py` | the six-test family (T1D − ND by donor permutation, H6 partial rank association), BH over six |
| crossstudy | `crossstudy.py` | 6551 against Peters et al.: TCF1/TOX T-cell shares by tissue and islet distance (audit context only) |
| cd8 | `cd8_roc.py` | per-donor CD8 ROC, the evidence the CD8 ruling was made on |
| figures | `figures.py`, `deck.py` | SVG/PDF/PNG figures, the exact plotted tables, deck.pptx and deck.pdf |
| verify | `delivery.py` | receipt and reconciliation checks, `STATUS.md` |

Each stage is a transaction: its fingerprint is the SHA-256 of every input file, its own code and its
settings; outputs are staged and renamed into place; a failure leaves `failure.json`. A changed input
or code file makes a new task folder, and the manifest points at the current one.

## Rules this lane holds

- **Normalization has two jobs.** Between donors, UniFORM aligns each marker. Within a donor, typing that
  involves more than one marker relies on the magnitude differences between them (maintainer,
  2026-10-19), so every marker in a comparison must be read the same way — which is why CD8 is no longer
  read unshifted in the refinement while the markers it is weighed against are shifted.
- **No support cutoffs** (ruled 2026-10-17). An endpoint is unavailable only when it cannot be computed
  — zero denominator, no variance, rank-deficient design — and it is reported with its reason, never as
  a measured zero. Bookkeeping p = 1 for multiplicity is never displayed.
- **Refined identity is phenotyping, not a sensitivity** (ruled 2026-10-17): the accepted broad label,
  kept only where `refined_class` agrees, then its subclass (Alpha/Beta/Delta, or `immune_class`).
- **Every donor is typed the same way** (ruled 2026-10-18); 6579 and 6551 are in every analysis. 6551
  also anchors the cross-study comparison, which is never used to tune, gate or release anything.
- Cells without a positive measured nucleus are excluded everywhere; cells are joined by exact UUID.
- The donor is the unit of inference; AAb+ is shown, not tested. Results are retrospective/exploratory.

The four-donor same-slide integration is owned by another session and is not produced here.
