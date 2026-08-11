# Review lane

Producers for the evidence under `data/pilots/`. They call the **production
APIs** on a few donors and write only review artifacts — they never write a stage
manifest and never mutate `data/runs/`.

These exist because the alternative kept happening: a session measures something,
puts a figure on disk, and the script that made it disappears. Every number in
`docs/METHOD_HISTORY.md` and `docs/INTENSITY_SCALE.md` should be reproducible
from here.

Run from the `Phenocycler_Analysis` directory. Paths resolve relative to the
package, so no absolute paths are baked in.

## Pilots

```bash
python scripts/review/run_pilot.py --donors 6436 6450 6566 --label baseline_p999
```

Runs `restore_donor_from_paths` + `type_simultaneous_cells` per donor, reading an
existing run's `cells/`, `redsea/` and `geometry_qc/` — the only three stages a
pilot needs. **~1.2 min per million cells**, so three donors take about five
minutes and all twenty well under an hour. If a session is about to spend hours
recomputing REDSEA to get labels, it is doing the wrong thing.

The CLI has no donor selector, which is why the pilot goes through the API
directly. That is the supported path, not a workaround.

## Method evidence

| script | what it answers |
|---|---|
| `divisor_survey.py` | what every candidate divisor would choose on each pair's own control, and how each compares to that pair's NNMF anchor |
| `pixel_anchors.py` | per donor-marker pixel saturation and the zero-clamp fraction, from tiled level-0 reads (never a whole channel — that stalls the disk) |
| `scan_acquisition.py` | exposure, responsivity and calibration epoch per donor-channel, read from the cohort manifest's declared qptiffs |

## Figures

| script | figure |
|---|---|
| `make_figures.py` | `fig1` divisor landscape per donor, `fig2` rule-vs-anchor, `fig3` composition and ambiguity attribution |
| `fig_scale.py` | `fig4` exposure spread and REDSEA zero-clamp per donor-marker |
| `fig_floor_vs_successor.py` | `fig5` Otsu floor vs the order-statistic successor at three prevalences |

Figures follow the standing rule: window on the workflow's own measured
quantities, never a percentile black point or gamma < 1, and state what falls
outside the axis rather than silently clipping it.

## Two traps these scripts encode

- **Join on the manifest's `channel_map`, never the qptiff `Biomarker` string.**
  The scanner writes `E-cadherin`, `Ker8-18`, `Keratin 5`, `Pan-Cytokeratin`;
  matching raw strings silently drops four gate markers.
- **`evaluate_locked_pair` returns `target_redsea` at QC-retained length but the
  `full_*` masks and `call_idx` at full donor length.** Mixing them raises
  IndexError at best and misaligns silently at worst.
