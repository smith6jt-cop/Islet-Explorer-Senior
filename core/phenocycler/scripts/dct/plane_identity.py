"""Plane-identity gate: does each marker NAME point at the right image PLANE?

Reads one tile, measures per-cell mean intensity of every selected plane over the
mask, and Spearman-correlates it against that marker's own column in the run's
expression parquet, joined on object_id.  Every marker must correlate with itself
and with nothing else.  A wrong name->plane pairing mislabels every cell in the
pilot with every other gate green, and no other check in the plan can see it.
"""
import sys, json, time
sys.path.insert(0,'/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/core/phenocycler')
sys.path.insert(0,'/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/core/phenocycler/scripts/dct')
import numpy as np, pandas as pd, tifffile
from scipy.stats import spearmanr
from phenocycler.redsea import qptiff_channels
from phenocycler.artifacts import iter_qupath_geojson_features
from tiff_window import read_window

ROOT='/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior'
DONOR='6539'; TILE=2048
man=json.load(open(f'{ROOT}/core/phenocycler/data/qupath_manifest.json'))
im=[i for i in man['images'] if i['donor_id']==DONOR][0]
qp=im['qptiff']['path']; gj=im['cell_geojson']['path']
name2marker={c['channel_name']:c['marker'] for c in im['channel_map']}
xw=json.load(open(f'{ROOT}/core/phenocycler/scripts/dct/deepcell_types_crosswalk.json'))['entries']

tf,pages,ome,(H,W),px = qptiff_channels(qp,0)
print(f"{DONOR}: {H}x{W}, {len(pages)} planes, {px:.6f} um/px", flush=True)

# --- pick the densest tile from the archived mask -------------------------
mp = tifffile.TiffFile('/data/codex_phenocycler/archive_2026-07-25/masks/6539.tif').pages[0]
best=(None,-1)
for i in range(H//TILE):
    for j in range(W//TILE):
        sub=read_window(mp,i*TILE,(i+1)*TILE,j*TILE,(j+1)*TILE)
        n=np.unique(sub[sub>0]).size
        if n>best[1]: best=((i,j),n)
(ti,tj),ncell=best
y0,x0=ti*TILE,tj*TILE
print(f"densest 2048 tile: row {ti} col {tj}  origin ({y0},{x0})  {ncell} cells", flush=True)
mask=read_window(mp,y0,y0+TILE,x0,x0+TILE).astype(np.int64)

# --- label -> object_id, in the order rasterize_mask assigns them ---------
t=time.time(); oids=[]
for ft in iter_qupath_geojson_features(gj): oids.append(ft.get("id"))
print(f"streamed {len(oids):,} object_ids from geojson in {time.time()-t:.0f}s", flush=True)
labels=np.unique(mask[mask>0])
assert labels.max()<=len(oids), "mask label exceeds geojson feature count"
tile_oid=np.array([oids[l-1] for l in labels])

# --- per-cell mean intensity of every selected plane ----------------------
sel=[(nm,name2marker[nm]) for nm in ome if nm in name2marker
     and not (nm in xw and xw[nm]['target'] is None)]
print(f"{len(sel)} planes selected (crosswalk nulls dropped)", flush=True)
order=np.argsort(labels); lut=np.zeros(labels.max()+1,dtype=np.int64); lut[labels]=np.arange(labels.size)
flat=mask.ravel(); nzp=flat>0; idx=lut[flat[nzp]]
counts=np.bincount(idx,minlength=labels.size).astype(np.float64)
meas={}
for nm,mk in sel:
    ch=[c['channel_index'] for c in im['channel_map'] if c['channel_name']==nm][0]
    assert ome[ch]==nm, f"OME name mismatch at plane {ch}"
    w=read_window(pages[ch],y0,y0+TILE,x0,x0+TILE).astype(np.float64).ravel()[nzp]
    with np.errstate(invalid='ignore',divide='ignore'):
        meas[mk]=np.where(counts>0, np.bincount(idx,weights=w,minlength=labels.size)/np.where(counts>0,counts,1), np.nan)
M=pd.DataFrame(meas); M['object_id']=tile_oid

# --- join to the run's own expression table on object_id ------------------
E=pd.read_parquet(f'{ROOT}/core/phenocycler/data/runs/30f6689e4d27d295f890/expression/donor_id={DONOR}/data.parquet')
E['object_id']=E['object_id'].astype(str); M['object_id']=M['object_id'].astype(str)
J=M.merge(E,on='object_id',suffixes=('_tile','_run'),how='inner')
print(f"joined {len(J):,} of {len(M):,} tile cells to the expression table on object_id\n", flush=True)

cols=[m for m in meas if f'{m}_run' in J.columns or m in E.columns]
R=np.zeros((len(cols),len(cols))); N=np.zeros(len(cols),dtype=int)
for a,ma in enumerate(cols):
    x=J[f'{ma}_tile'] if f'{ma}_tile' in J else J[ma]
    for b,mb in enumerate(cols):
        y=J[f'{mb}_run'] if f'{mb}_run' in J else J[mb]
        xv=np.asarray(x,dtype=float); yv=np.asarray(y,dtype=float)
        ok=np.isfinite(xv)&np.isfinite(yv)          # expression table records missing as NaN, not 0
        R[a,b]=spearmanr(xv[ok],yv[ok]).statistic if ok.sum()>=50 else np.nan
        if a==b: N[a]=int(ok.sum())
from gates import assess_self_match

verdict = assess_self_match(cols, R, N, min_obs=50, min_self=0.30, margin=0.0)
print(f"{'marker':<18}{'self rho':>9}{'best other':>12}  {'which':<16}{'n':>7}")
for name, self_s, best, who, n, ok, why in sorted(verdict.rows, key=lambda r: (r[5], r[1])):
    note = "" if ok else f"   <-- FAIL: {why}"
    print(f"{name:<18}{self_s:>9.3f}{best:>12.3f}  {who:<16}{n:>7}{note}")
finite = [r[1] for r in verdict.rows if np.isfinite(r[1])]
print(f"\nself-correlation: median {np.median(finite):.3f}, min {min(finite):.3f}"
      if finite else "\nno finite self-correlations")
print(f"{sum(1 for r in verdict.rows if not r[5])} of {len(verdict.rows)} markers failed")
print("PLANE-IDENTITY GATE:", "PASS" if verdict.passed else "FAIL", "--", verdict.reason)
raise SystemExit(0 if verdict.passed else 1)
