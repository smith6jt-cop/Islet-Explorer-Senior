# REDSEA and its GPU path — adversarial review, 2026-08-27

**Status:** the CPU path is sound and faithful to upstream. **The GPU path cannot complete a
single donor** and must not be enabled until two defects are fixed. Written before spending
compute on a 20-donor rebuild, which is what it saved.

**Verdict: run production on CPU. `use_gpu = false`.**

---

## 1. The GPU path crashes on every donor (reproduced)

`redsea.py:510`, inside the `backend.on_gpu` branch:

```python
corrected_sum = d + comp_mode * e - alpha * (F @ e)
```

`d`, `e` and `F` are CuPy; `comp_mode` and `alpha` are **NumPy arrays**, because
`apply_no_spillover_mask` (`redsea.py:363-364`) promotes them to length-C vectors whenever any
channel is excluded. CuPy refuses implicit host↔device mixing for non-scalar operands. The CPU
branch is immune only because `redsea.py:526-528` wraps them in `np.asarray(...)`.

**This fires on every donor, not occasionally.** `config.ini` sets
`exclude_no_spillover = true`, and `marker_registry.json` has 9 markers that are not
`redsea_corrected` — `DAPI, FOXP3, Ki67, PCNA, SOX2, TCF_1, TOX, TP63, Collagen_IV` — so
`n_skipped > 0` always and the vector path is always taken.

Reproduced on the repo's own two-cell fixture:

```
CPU         scalar : OK  [[11.5, 0.8], [0.0, 0.2]]
CPU         vector : OK  [[11.5, 0.7], [0.0, 0.3]]
GPU         scalar : OK  [[11.5, 0.8], [0.0, 0.2]]   <- matches CPU exactly
GPU         vector : TypeError: Unsupported type <class 'numpy.ndarray'>
```

It fires at `redsea.py:902`, i.e. **after** rasterisation, the boundary band, all 59 channel
reads and bincounts, and the contact matrix — the entire expensive part of the donor, thrown
away. `run_redsea` uses `map_donors(..., on_error="raise")` (`redsea.py:1063`), so the cohort
aborts. Fix: cast both onto the device.

## 2. The largest donors exceed VRAM (reproduced)

`redsea.py:852-861` transfers whole images with **no tiling or chunking**: int32 mask + int32
nucleus + float32 channel = **12 bytes/pixel**, held live across the 59-channel loop. Against
18.29 GiB free on the RTX 4000 Ada (20 GiB card):

| donor | W × H | Gpx | 3 whole-image arrays | fits? |
|---|---|---|---|---|
| 6505 | 35520 × 56160 | 1.995 | **22.29 GiB** | **NO** |
| 6374 | 35520 × 51120 | 1.816 | **20.29 GiB** | **NO** |
| 6593 | 35520 × 44640 | 1.586 | 17.72 GiB | marginal |
| 6566 | 35520 × 43200 | 1.534 | 17.15 GiB | marginal |
| 6551 | 35520 × 42480 | 1.509 | 16.86 GiB | marginal |

That bound is **optimistic** — it excludes the band arrays, `band_idx` (int64), and the
`chg[band_idx_g]` gather temporary (computed twice, `redsea.py:864` and `:870`). Nothing is
ever freed: no `del`, no `free_all_blocks`, and all five index arrays stay live into
`compensate`. `get_backend` (`gpu.py:86-88`) advertises "size-gated GPU use" in its module
docstring and implements **no size gate at all** — only device presence.

## 3. There is no test coverage whatsoever

Eight `# pragma: no cover - GPU only` branches (5 in `gpu.py`, 3 in `redsea.py`) plus five
unmarked dual-path methods (`bincount`, `diags`, `clip_`, `asarray`, `csr_from_scipy`).
**No test in the repository constructs a `Backend`** — every `compensate` call in
`tests/current/test_redsea.py` omits `backend`, so `redsea.py:499` is never entered. `gpu.py`
has exactly one commit, the subtree import. This code has most likely never executed.

One test passing `backend=get_backend(True)` with a **vector** `alpha`, skipped when
`gpu_available()` is False, would have caught defect 1 in seconds.

## 4. Where the GPU would actually win, once fixed

Measured at 200 Mpx / 300k cells:

| | time |
|---|---|
| CPU bincount | 3.39 s |
| GPU bincount (resident) | 0.013 s |
| H2D transfer, 200 Mpx float32 | 0.075 s |

Extrapolated to donor 6505 (1.995 Gpx × 4 bincounts × 59 channels): **~2.2 h CPU vs ~1 min
GPU.** That is the prize, and it is worth fixing for.

The sparse matmul is **not** worth it: cold cuSPARSE start is ~3.5 s and `map_donors` uses a
spawn context, so every donor pays it fresh; and `compensate_compartments` (`redsea.py:398`)
takes no backend and recomputes `F @ edge` on CPU regardless, so the GPU result is only used
for a cross-check.

## 5. GPU/CPU numerical agreement — measured

| operation | result |
|---|---|
| `bincount` | **bit-identical**, max\|Δ\| = 0 |
| sparse `F @ edge` | **not** bit-identical — max\|Δ\| 1.9e-6, **max relative 2.0e-7** |

`bincount` is exact for a specific reason worth recording: every qptiff page is `uint16`, and
`read_channel` casts to float32, which represents 0..65535 exactly; CuPy's `bincount`
accumulates in float64, so every partial sum is exact and summation order is irrelevant.
**That guarantee dies the moment anything upstream makes channel data non-integral.**

The matmul divergence has a consequence: the cross-check at `redsea.py:925-940` compares
**CPU** `compensate_compartments` output against **GPU** `compensate` output at
`rtol=1e-9, atol=0.0`. Measured divergence is 2.0e-7 — 200× over that threshold — so the guard
would fail donors spuriously on GPU. Either give `compensate_compartments` the backend, or add
an `atol`.

## 6. Provenance gap

`use_gpu` and `gpu_device` are **not** in the REDSEA stage fingerprint (`pipeline.py:280-288`),
and `use_gpu` does not enter `run_id`. Since the matmul is not bit-identical, a GPU-built and a
CPU-built stage are content-addressed identically while carrying different bytes, and a
mixed-backend cohort would look uniform in the manifest. `config.py:537` also exposes
`PHENOCYCLER_USE_GPU`, so the backend can change with no config diff. Add `use_gpu` to the
fingerprint if the GPU path is ever enabled.

---

## Findings about the CPU path (which is what we are running)

**Fine, explicitly:** `norm_form = donor` is faithful to upstream — `F[i,j] = M[i,j]/B[j]`
where `B` is the contact-matrix column sum, matching redseapy's `cellBoundaryTotal`;
`downsample = 1.0`; `alpha = 1.0`; the clamp's position (sum-domain, pre-division — equivalent
to mean-domain, no sign error); and the `spawn` multiprocessing context, which correctly avoids
the fork+CUDA hazard.

**`edge_radius = 0` is correct, not a bug.** `skimage.morphology.diamond(0)` is a 1×1 all-true
footprint, so the dilation at `redsea.py:810` is the identity and `band == inner` exactly —
precisely the 1-px inner rim the docstring promises. (Two asides: the `binary_dilation` call is
therefore dead work on a ~2 Gpx boolean array, and it emits a skimage deprecation warning that
will break on 0.28.)

**Worth fixing or recording:**

- **`METHOD_HISTORY.md:201` cites the wrong line.** The clamp that reaches disk is
  `redsea.py:455` (in `compensate_compartments`), not `:529` — `:529`'s output feeds only a
  cross-check and a log line. The 0.5–19% figures stand.
- **`comp_mode = 1` inflates isolated cells.** With no neighbours, `F @ edge` is zero, so the
  result is `data + edge` — the rim counted twice with no compensating subtraction. Upstream
  deletes such cells (`redseapy/redsea.py:198-209`); we don't. The only isolated-cell test uses
  `comp_mode=0`, where passthrough is correct, so the production behaviour is **untested**.
  `n_isolated` is logged (`redsea.py:905-907`) but gated by nothing.
- **`gap_bridge = 1` bridges 2 px, not 1** — both sides grow by `grow+0.5`. Measured: gaps of
  1 px and 2 px are bridged, 3 px is not. Docstring off by 2×.
- **The clamp hides an unquantified over-subtraction.** `F` redistributes 100% of a cell's rim
  among its contacting neighbours even when most of that rim faces background, so
  over-subtraction is structurally expected. The signed pre-clip sum is never persisted, and
  nothing downstream distinguishes a clamped 0 from a true 0.
- **`clamp_frac` (`redsea.py:533`) over-counts** — computed over all GeoJSON rows including
  zero-area fragments, which register as "clamped" without anything having been clamped. It
  feeds a log line only; do not quote it as the clamp rate.
- **`apply_spillover_context_eligibility` perturbs surviving neighbours.** Excluded cells are
  zeroed in both masks (`redsea.py:674-677`), so their footprint becomes background — which
  `gap_bridge` may then annex, manufacturing contacts between cells that were physically
  separated, while their real rim no longer donates spillover and their neighbours are
  under-corrected.
- **A documented default contradicts config.** `geometry_qc.py:19-23` says an analysis-excluded
  polygon may still contribute spillover context, but `missing_geometry_policy = exclude` plus
  `geometry_qc.py:533-537` removes nucleus-less InstanSeg fragments from spillover context.
- **`save_intermediates` downcasts to float32** (`redsea.py:553-556`), so a recompensate run is
  not value-preserving against a full run. Dormant while the default is False.
- **`recompensate_from_intermediates` has no cross-check** and never passes the backend to
  `compensate_compartments` — if used as an OOM workaround, both guards are lost.

---

## Before the GPU path is enabled

1. Cast `alpha`/`comp_mode` onto the device at `redsea.py:510`.
2. Add the whole-image size gate `get_backend` already claims to have (~1.4 Gpx on this card),
   or tile the channel loop; free the index arrays before `compensate`.
3. Write one GPU test with a vector `alpha`, skipped when `gpu_available()` is False.
4. Give `compensate_compartments` the backend, or loosen `redsea.py:929`'s `atol=0.0`.
5. Add `use_gpu` to the REDSEA fingerprint (`pipeline.py:280-288`).
6. Paired CPU/GPU run on the smallest donor (6380, 0.382 Gpx) and diff the parquet before
   trusting the other 19.
