"""Single-tile feasibility check: does DeepCellTypes call THIS tissue sensibly?

One 2048-px tile from donor 6539, stock `predict()` with the package's own default
per-FOV normalization -- no hook, no halo, no histogram job. Then three threshold-free,
rank-based readouts against the run's own expression parquet, joined on `object_id`:

  1. measured CD3e / CD68 / CD20 vs the model's T-cell / Macrophage / B-cell probability
  2. measured E-cadherin / Pan-Cytokeratin vs its Epithelial probability
  3. where its endocrine calls fall relative to the QuPath islet geometry

Nothing here is a threshold, a gate or a cutoff: every readout is a rank statistic, and the
islet comparison uses QuPath GEOMETRY (`cell_region`), not any phenotype label.

If the model does not call this tissue sensibly, that is discoverable now rather than after
the whole harness exists.
"""
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import tifffile
from scipy.stats import spearmanr

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "core" / "phenocycler"))
sys.path.insert(0, str(HERE))

from phenocycler.redsea import qptiff_channels                    # noqa: E402
from phenocycler.artifacts import iter_qupath_geojson_features    # noqa: E402
from tiff_window import read_window                               # noqa: E402
from deepcell_types.utils import resolve_supported_marker         # noqa: E402
from deepcell_types import predict                                # noqa: E402

# Marker-dropout arm. DAPI is the panel's ONLY nuclear channel and the crosswalk withholds
# it by default (a panel-membership decision, recorded in deepcell_types_crosswalk.json).
# `--dapi` maps DAPI -> dsDNA and hands it to the model, so the cost of that exclusion is
# measured rather than assumed.
WITH_DAPI = "--dapi" in sys.argv
# Tile choice. The plan's per-FOV normalization argument assumes ~3,560 cells (0.88x the
# training FOV). The DENSEST tile on this slide holds 10,800 -- 2.7x it, off-distribution.
# `--median-tile` picks the tile whose cell count is closest to the slide median instead.
MEDIAN_TILE = "--median-tile" in sys.argv
# Both released checkpoints load through stock predict(); the SSL pretrain->finetune arm
# is the second one the plan fetches.
PTFT = "--ptft" in sys.argv

DONOR, TILE, MPP = "6539", 2048, 0.5
RUN = "30f6689e4d27d295f890"
MASK = "/data/codex_phenocycler/archive_2026-07-25/masks/6539.tif"
CKPT = pathlib.Path.home() / ".deepcell" / "models" / (
    "deepcell-types_2026-06-23_resmlp_ptft.pt" if PTFT else "deepcell-types_2026-06-15_resmlp.pt")
OUT = ROOT / "data" / "development" / "deepcell_types_review" / DONOR
OUT.mkdir(parents=True, exist_ok=True)
ARM = (("with-dapi" if WITH_DAPI else "no-dapi") + ("_median-tile" if MEDIAN_TILE else "")
       + ("_ptft" if PTFT else ""))

# Child classes count toward their parent -- the paper's own hierarchy, vendored in
# dct_figures/scoring.py as CELL_TYPE_HIERARCHY. Used here only to group probabilities.
TCELL = ["Tcell", "CD4T", "CD8T", "Treg", "NKT"]


def log(*a):
    print(*a, flush=True)


# --------------------------------------------------------------- inputs
man = json.loads((ROOT / "core/phenocycler/data/qupath_manifest.json").read_text())
im = next(i for i in man["images"] if i["donor_id"] == DONOR)
xw = json.loads((HERE / "deepcell_types_crosswalk.json").read_text())["entries"]
bc = json.loads((HERE / "broad_class_crosswalk.json").read_text())
EPI = [c for c, e in bc["entries"].items() if e["lineage"] == "Epithelial"]
ENDO = ["AlphaCell", "BetaCell", "Endocrine"]

tf, pages, ome, (H, W), px_um = qptiff_channels(im["qptiff"]["path"], 0)
log(f"ARM: {ARM}")
log(f"{DONOR}: {H}x{W}, {len(pages)} planes, {px_um:.6f} um/px (declaring mpp={MPP})")

# channel_name -> (plane index, registry marker key, DCT name); crosswalk nulls dropped.
sel = []
for c in im["channel_map"]:
    nm = c["channel_name"]
    if nm == "DAPI" and WITH_DAPI:
        dct = "dsDNA"
    elif nm in xw:
        dct = xw[nm]["target"]
        if dct is None:
            continue
    else:
        dct = resolve_supported_marker(nm)
        if dct is None:
            continue
    assert ome[c["channel_index"]] == nm, f"OME name mismatch at plane {c['channel_index']}"
    sel.append((c["channel_index"], c["marker"], dct))
log(f"{len(sel)} planes selected; {len(im['channel_map']) - len(sel)} dropped by the crosswalk")

# --------------------------------------------------------------- densest tile
mp = tifffile.TiffFile(MASK).pages[0]
counts = {}
for i in range(H // TILE):
    for j in range(W // TILE):
        sub = read_window(mp, i * TILE, (i + 1) * TILE, j * TILE, (j + 1) * TILE)
        counts[(i, j)] = np.unique(sub[sub > 0]).size
if MEDIAN_TILE:
    # Median over tiles this pipeline would actually process (>=1000 cells), then the
    # tile closest to it -- a typical field, not the extreme one.
    pop = {k: v for k, v in counts.items() if v >= 1000}
    med = np.median(list(pop.values()))
    best = min(pop, key=lambda k: abs(pop[k] - med))
    log(f"median-density tile (median over {len(pop)} populated tiles = {med:.0f})")
else:
    best = max(counts, key=counts.get)
bn = counts[best]
y0, x0 = best[0] * TILE, best[1] * TILE
log(f"tile: row {best[0]} col {best[1]} origin ({y0},{x0}) with {bn} cells")
mask = read_window(mp, y0, y0 + TILE, x0, x0 + TILE).astype(np.int32)

# label -> object_id, in the order rasterize_mask assigned them
t = time.time()
oids = [ft.get("id") for ft in iter_qupath_geojson_features(im["cell_geojson"]["path"])]
log(f"streamed {len(oids):,} object_ids in {time.time() - t:.0f}s")
assert mask.max() <= len(oids), "mask label exceeds geojson feature count"

# --------------------------------------------------------------- raw stack
raw = np.empty((len(sel), TILE, TILE), dtype=np.float32)
for k, (ch, _, _) in enumerate(sel):
    raw[k] = read_window(pages[ch], y0, y0 + TILE, x0, x0 + TILE)
names = [d for _, _, d in sel]
log(f"raw stack {raw.shape} {raw.nbytes / 1e6:.0f} MB")

# --------------------------------------------------------------- inference
import torch  # noqa: E402

dev = "cuda" if torch.cuda.is_available() else "cpu"
t = time.time()
res = predict(raw, mask, names, MPP, model_name=str(CKPT), device=dev,
              return_probabilities=True)
log(f"predict() on {len(res.cell_types):,} cells in {time.time() - t:.1f}s on {dev}")

# Probability columns must follow the model's own index order, not sorted order.
from deepcell_types.config import DCTConfig  # noqa: E402
ct_order = [k for k, _ in sorted(DCTConfig().ct2idx.items(), key=lambda kv: int(kv[1]))]

P = pd.DataFrame(res.probabilities, columns=ct_order)
P["object_id"] = [oids[i - 1] for i in res.cell_indices]
P["call"] = res.cell_types

# --------------------------------------------------------------- join
E = pd.read_parquet(ROOT / f"core/phenocycler/data/runs/{RUN}/expression/donor_id={DONOR}/data.parquet")
E["object_id"] = E["object_id"].astype(str)
P["object_id"] = P["object_id"].astype(str)
J = P.merge(E, on="object_id", how="inner", suffixes=("_p", "_e"))
log(f"joined {len(J):,} of {len(P):,} tile cells to the expression table\n")

J.to_parquet(OUT / f"feasibility_tile_{ARM}.parquet", index=False)

# --------------------------------------------------------------- readouts
def rho(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)          # NaN = unavailable, never zero
    if ok.sum() < 50:
        return np.nan, int(ok.sum())
    return spearmanr(x[ok], y[ok]).statistic, int(ok.sum())


log("=" * 72)
log("READOUT 1 - measured immune marker vs the model's probability for that cell type")
log(f"{'marker':<18}{'vs class(es)':<34}{'rho':>8}{'n':>8}")
for mk, grp, label in [("CD3e", TCELL, "Tcell+CD4T+CD8T+Treg+NKT"),
                       ("CD68", ["Macrophage"], "Macrophage"),
                       ("CD20", ["Bcell"], "Bcell"),
                       ("CD20", ["Bcell", "Plasma"], "Bcell+Plasma"),
                       ("CD68", ["Macrophage", "Monocyte"], "Macrophage+Monocyte"),
                       ("MPO", ["Neutrophil"], "Neutrophil")]:
    if mk not in J:
        continue
    r, n = rho(J[mk], J[grp].sum(axis=1))
    log(f"{mk:<18}{label:<34}{r:>8.3f}{n:>8,}")

log("")
log("READOUT 2 - measured epithelial marker vs epithelial probability")
log(f"{'marker':<18}{'vs class(es)':<34}{'rho':>8}{'n':>8}")
for mk in ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM"]:
    if mk not in J:
        continue
    for grp, label in [(["Epithelial"], "Epithelial (alone)"), (EPI, "Epithelial lineage (12)")]:
        r, n = rho(J[mk], J[grp].sum(axis=1))
        log(f"{mk:<18}{label:<34}{r:>8.3f}{n:>8,}")

log("")
log("READOUT 3 - endocrine calls vs QuPath islet GEOMETRY (cell_region), not any label")
if "cell_region" in J:
    log(f"cell_region values: {J['cell_region'].value_counts().to_dict()}")
    endo_p = J[ENDO].sum(axis=1)
    for reg, sub in J.groupby("cell_region"):
        e = endo_p[sub.index]
        log(f"  {str(reg):<22} n={len(sub):>6,}  median P(endocrine)={np.nanmedian(e):.4f}  "
            f"mean={np.nanmean(e):.4f}  called endocrine={(sub['call'].isin(ENDO)).mean():.1%}")
    # AUROC of P(endocrine) separating islet-core cells from the rest -- a real binary,
    # supplied by geometry rather than by any phenotype call.
    is_islet = J["cell_region"].astype(str).str.contains("slet", case=False, na=False)
    if is_islet.any() and (~is_islet).any():
        from scipy.stats import rankdata
        r_ = rankdata(endo_p.values)
        n1, n0 = int(is_islet.sum()), int((~is_islet).sum())
        auc = (r_[is_islet.values].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        log(f"  AUROC P(endocrine) islet vs non-islet: {auc:.3f}  (n_islet={n1:,}, n_other={n0:,})")

log("")
log("Top 15 calls in this tile:")
for k, v in J["call"].value_counts().head(15).items():
    log(f"  {k:<24}{v:>7,}  {v / len(J):>6.1%}")
log("")
log(f"wrote {OUT / f'feasibility_tile_{ARM}.parquet'}")
