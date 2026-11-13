"""Build the PhenoCycler-only, per-donor islet inflammation ordering from native outlines.

measure builds the shared density scale and independent donor graphs; attach copies the affected
cell AnnData files to a staging directory, then joins the validated islet information.
Run attach at each data host to avoid transferring complete sections. It neither infers
new cell identities nor assigns a cell to one of several overlapping peri-islet bands.
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import anndata as ad
from anndata.io import read_elem, write_elem
import h5py
import numpy as np
import pandas as pd
import shapely
from shapely.affinity import scale
from shapely.geometry import shape

PREVIOUS_VERSION = 'native-core-peri20-density-rank-v1'
VERSION = 'native-core-peri20-shared-density-v2'
BAND_UM = 20
QUAD_SEGS = 64
CELL_FIELDS = {'islet_core_area_um2':'core_area_um2', 'islet_peri_area_um2':'peri_area_um2',
               'core_immune_per_mm2':'core_immune_per_mm2','peri_immune_per_mm2':'peri_immune_per_mm2',
               'combined_immune_per_mm2':'combined_immune_per_mm2',
               'islet_pseudotime':'pseudotime', 'islet_pseudotime_rank':'pseudotime_rank',
               'islet_inflammation_rank_within_donor':'within_donor_density_rank',
               'islet_donor_graph_dpt':'donor_graph_dpt'}
PREV = {'islet_pseudotime':'islet_composition_dpt','islet_pseudotime_rank':'islet_composition_dpt_rank',
        'islet_inflammation':'islet_core_only_burden','core_immune_per_mm2':'islet_previous_core_immune_per_mm2'}

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def inflammation_order(density):
    """Average ranks of exact ties within the supplied reference population. Preserve missing."""
    s=pd.Series(density,dtype='float64')
    valid=s.dropna()
    if (~np.isfinite(valid) | valid.lt(0)).any():
        raise ValueError('Density must be finite and nonnegative when assessed')
    result=pd.Series(np.nan,index=s.index,dtype='float64')
    if valid.empty:
        return result, 'no assessed islets'
    if valid.nunique()==1:
        if valid.iloc[0]==0:
            result.loc[valid.index]=0
            return result,'all zero'
        return result,'constant positive; no relative ordering'
    r=valid.rank(method='average')
    result.loc[r.index]=(r-r.min())/(r.max()-r.min())
    return result,'ordered'

def regions(job):
    donor, fpath, img, qroot, troot = job
    fields=['islet_id','outline_wkb','outline_area_um2','geometry_reference','tissue_reference','n_core','n_peri','n_core_immune','n_peri_immune','geometry_valid','complete_20um_tissue_margin','core_pancreatic_tissue_fraction']
    frame=pd.read_parquet(fpath,columns=fields).set_index('islet_id')
    if not frame.index.is_unique: raise ValueError(f'{donor}: duplicate native islet IDs')
    def path_for(key,relocate):
        if frame[key].nunique()!=1: raise ValueError(f'{donor}: ambiguous geometry source')
        p=Path(frame[key].iloc[0])
        return Path(relocate)/p.name if relocate else p
    geo=path_for('geometry_reference',qroot); anno=path_for('tissue_reference',troot)
    if not np.isclose(img['pixel_size_um_x'],img['pixel_size_um_y'],rtol=0,atol=1e-9):
        raise ValueError(f'{donor}: geometry requires calibrated square pixels')
    px=img['pixel_size_um_x']
    features=json.loads(geo.read_bytes())['features']
    def microns(g):return scale(g,xfact=px,yfact=px,origin=(0,0))
    ids={str(f['id']) for f in features if f['properties'].get('size_class')=='Islet'}
    if ids!=set(frame.index): raise ValueError(f'{donor}: full native islet census mismatch')
    cores=shapely.from_wkb(frame.outline_wkb.to_numpy())
    valid=shapely.is_valid(cores) & ~shapely.is_empty(cores)
    if not np.array_equal(valid,frame.geometry_valid):
        raise ValueError(f'{donor}: geometry validity differs from cell membership census')
    if not np.allclose(shapely.area(cores),frame.outline_area_um2,rtol=2e-4,atol=1):
        raise ValueError(f'{donor}: native core area disagrees with membership source')
    clusters=[microns(shape(f['geometry'])) for f in features if f['properties'].get('size_class')=='EndocrineCluster']
    del features
    roi=json.loads(anno.read_bytes())['features']
    def annotate(name):
        return microns(shapely.union_all([shapely.make_valid(shape(f['geometry'])) for f in roi if (f['properties'].get('classification') or {}).get('name')==name]))
    tissue=annotate('Tissue').difference(annotate('pLN'))
    if not tissue.is_valid or tissue.area<=0: raise ValueError(f'{donor}: pancreatic tissue mask invalid')
    # Exactly the polygon blockers used in native 20 µm point membership; a peri band
    # can overlap another islet's band, but never an islet core or endocrine cluster.
    blocks=np.concatenate([cores[valid],np.asarray(clusters,dtype=object)])
    tree=shapely.STRtree(blocks)
    result=[]
    for core, ok, complete in zip(cores,valid,frame.complete_20um_tissue_margin,strict=True):
        if not ok:
            result.append((np.nan,np.nan));continue
        inner=core if complete else core.intersection(tissue)
        # Match the Euclidean distance predicate as a polygonal buffer. 64 segments/
        # quadrant bounds sagitta at a 20 µm radius to <0.002 µm.
        outer=core.buffer(BAND_UM,quad_segs=QUAD_SEGS)
        band=outer.difference(core)
        candidates=tree.query(band,predicate='intersects')
        # Exclude self (the ring shares an edge with its core), and union overlaps once.
        others=[shapely.make_valid(blocks[i]) for i in candidates if blocks[i] is not core]
        if others: band=band.difference(shapely.union_all(others))
        if not complete:band=band.intersection(tissue)
        result.append((inner.area,band.area))
    frame[['core_area_um2','peri_area_um2']]=np.asarray(result)
    if not np.allclose(frame.loc[valid,'core_area_um2']/shapely.area(cores[valid]),frame.loc[valid,'core_pancreatic_tissue_fraction'],atol=1e-6):
        raise ValueError(f'{donor}: clipped core differs from measured membership source')
    # Recount the actual native UUID membership independently of the region summary.
    mempath=Path(fpath).parent/'members.parquet'
    members=pd.read_parquet(mempath,columns=['islet','object_id','region','refined_class'])
    if members.duplicated(['islet','region','object_id']).any() or not members.islet.isin(frame.index).all():
        raise ValueError(f'{donor}: duplicate or unmapped membership')
    members['immune']=members.refined_class.eq('Immune')
    grouped=members.groupby(['islet','region'],observed=True).immune.agg(['size','sum'])
    for region in ('core','peri'):
        verified=grouped.xs(region,level=1).reindex(frame.index,fill_value=0)
        np.testing.assert_array_equal(verified['size'],frame['n_'+region])
        np.testing.assert_array_equal(verified['sum'],frame['n_'+region+'_immune'])
    return donor,frame[['core_area_um2','peri_area_um2','n_core','n_peri','n_core_immune','n_peri_immune']],{
        'native_geometry_sha256':sha(geo),'tissue_roi_sha256':sha(anno),'members_sha256':sha(mempath),
        'features_sha256':sha(fpath),'pixel_um':px,'boundary_clipped':int((~frame.complete_20um_tissue_margin).sum())}

def apply_order(o):
    for region in ('core','peri'):
        area=o[region+'_area_um2']
        supported=area.gt(0)&np.isfinite(area)&o['n_'+region].gt(0)
        o[region+'_density_supported']=supported
        o[region+'_immune_per_mm2']=(1e6*o['n_'+region+'_immune']/area).where(supported)
    o['density_eligible']=o.core_density_supported & o.peri_density_supported
    o['combined_area_um2']=o.core_area_um2+o.peri_area_um2
    o['n_combined_immune']=o.n_core_immune+o.n_peri_immune
    o['combined_immune_per_mm2']=(1e6*o.n_combined_immune/o.combined_area_um2).where(o.density_eligible)
    o['within_donor_density_rank']=np.nan
    o['burden']=pd.Series(None,index=o.index,dtype='object')
    o['density_order_status']=pd.Series(None,index=o.index,dtype='object')
    for donor,ids in o.groupby('donor',observed=True).groups.items():
        values=o.loc[ids,'combined_immune_per_mm2']
        rank,status=inflammation_order(values)
        o.loc[ids,'within_donor_density_rank']=rank
        o.loc[ids,'density_order_status']=status
        pos=values[values.gt(0)]
        o.loc[values[values.eq(0)].index,'burden']='Zero'
        o.loc[pos[pos.le(pos.median())].index,'burden']='Lower'
        o.loc[pos[pos.gt(pos.median())].index,'burden']='Upper'
    o['within_donor_burden']=pd.Categorical(o.burden,categories=['Zero','Lower','Upper'],ordered=True)
    o['density_order_status']=o.density_order_status.astype('category')
    density=o['combined_immune_per_mm2']
    maximum=density.max()
    o['pseudotime']=density / maximum if maximum > 0 else density.copy()
    # Use a single pooled reference for all donors. Tied measured zeros share zero.
    o['pseudotime_rank'],_=inflammation_order(density)
    positive=density[density.gt(0)]
    o['burden']=pd.Categorical(
        np.select([density.eq(0), density.gt(0) & density.le(positive.median()), density.gt(positive.median())],
                  ['Zero','Lower','Upper'], default=None),
        categories=['Zero','Lower','Upper'], ordered=True)
    return o

def measure(args):
    sentinel=args.output/'islet_inflammation_method.json'
    if sentinel.exists():
        receipt=json.loads(sentinel.read_text())
        if receipt.get('islets_sha256')!=sha(args.output/'pancreas_islets.h5ad'):
            raise ValueError('Existing islet stage is inconsistent; choose a new directory')
        raise FileExistsError('Islet staging already completed. Use a new output directory or run attach.')
    ad.settings.allow_write_nullable_strings=True
    args.output.mkdir(mode=0o700,parents=True,exist_ok=True)
    src=ad.read_h5ad(args.islets)
    if 'inflammation_order_version' in src.uns: raise ValueError('Use pre-update source; it must contain the reference composition DPT')
    o=src.obs.copy()
    if not o.index.is_unique or not np.array_equal(o.index.astype(str),o.islet_id.astype(str)): raise ValueError('Bad islet UUID index')
    images={str(int(im['donor_id'])):im for im in json.loads(args.manifest.read_text())['images']}
    jobs=[(d,args.features_root/f'donor_{d}/features.parquet',images[d],args.geometry_root,args.tissue_root) for d in o.donor.astype(str).unique()]
    parts=[]; source_hashes={}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for d,part,receipt in pool.map(regions,jobs):
            print(d,len(part),flush=True);part['donor_check']=d;parts.append(part);source_hashes[d]=receipt
    merged=pd.concat(parts,verify_integrity=True).reindex(o.index)
    if not merged.donor_check.astype(str).eq(o.donor.astype(str)).all(): raise ValueError('Geometry and AnnData donor UUIDs differ')
    for col in ('n_core','n_peri','n_core_immune','n_peri_immune'):np.testing.assert_array_equal(merged[col],o[col])
    o['composition_dpt']=o.pseudotime; o['composition_dpt_rank']=o.pseudotime_rank
    o['previous_core_only_burden']=o.burden
    o['previous_core_immune_per_mm2']=o.core_immune_per_mm2
    o['previous_core_density_rank']=o.pop('density_rank_in_donor')
    for col in ['core_area_um2','peri_area_um2']:o[col]=merged[col]
    src.obs=apply_order(o)
    src.uns['inflammation_order_version']=VERSION
    notes={'primary':'pseudotime is the combined measured immune density divided by the maximum assessed density in the complete cohort. Absolute density in immune cells/mm² is the primary disease comparison. pseudotime_rank uses average tied ranks rescaled once across all assessed islets. within_donor_density_rank retains the separate within-donor order.',
           'density':'(n_core_immune + n_peri_immune) / ((core_area_um2 + peri_area_um2)/1e6); no addition of compartment densities.',
           'regions':'Native islet core and exterior 20 µm band, clipped to annotated Tissue minus pLN. Band excludes every islet core and endocrine cluster, matching native region memberships.',
           'support':'Both regions require positive finite assessed area and at least one admitted cell. Zero immune cells is a measured zero; an empty region is unassessed. Zero inflammation maps to zero. A constant-positive donor receives its cohort score; its within-donor order is undefined. A cohort with a single positive density has a score of 1 and no relative order.',
           'burden':'Zero for observed zero combined density; Lower / Upper split positive combined densities at the pooled positive-islet median, ties at the median are Lower.',
           'old_dpt':'composition_dpt, composition_dpt_rank and X_islet_graph preserve the prior protein composition graph and its coordinates; unavailable graph positions do not remove density-ranked islets.',
           'band_buffer_quadrant_segments':QUAD_SEGS}
    notes['reference_max_immune_per_mm2']=float(src.obs.combined_immune_per_mm2.max())
    notes['positive_median_immune_per_mm2']=float(src.obs.loc[src.obs.combined_immune_per_mm2.gt(0),'combined_immune_per_mm2'].median())
    if 'trajectory_measurements' in src.uns:
        old=src.uns['trajectory_measurements'];old['previous_core_only_burden']=old.pop('burden','')
        old['burden']=notes['burden'];old['pseudotime']=notes['primary'];old['pseudotime_rank']=notes['primary']
    src.uns['inflammation_order_method']=notes
    out=args.output/'pancreas_islets.h5ad'
    tmp=out.with_suffix('.building.h5ad');src.write_h5ad(tmp,compression='gzip')
    check=ad.read_h5ad(tmp); np.testing.assert_allclose(check.obs.combined_immune_per_mm2,src.obs.combined_immune_per_mm2,equal_nan=True)
    np.testing.assert_allclose(check.obs.composition_dpt_rank,read_islet_column(args.islets,'pseudotime_rank'),equal_nan=True)
    np.testing.assert_allclose(check.obs.pseudotime_rank,src.obs.pseudotime_rank,equal_nan=True)
    tmp.replace(out)
    o=src.obs
    rows=[]
    for d,z in o.groupby('donor',observed=True):
        valid=z.dropna(subset=['pseudotime_rank'])
        if len(valid):
            assert valid.sort_values('combined_immune_per_mm2').pseudotime_rank.is_monotonic_increasing and valid.groupby('combined_immune_per_mm2').pseudotime_rank.nunique().max()<=1
        assert valid.loc[valid.combined_immune_per_mm2.eq(0),'pseudotime_rank'].eq(0).all()
        cmp=z[['pseudotime_rank','composition_dpt_rank']].dropna()
        rows.append(dict(donor=str(d),disease=str(z.disease.iloc[0]),islets=len(z),ranked=len(valid),
                         min_immune_per_mm2=valid.combined_immune_per_mm2.min(),median_immune_per_mm2=valid.combined_immune_per_mm2.median(),
                         max_immune_per_mm2=valid.combined_immune_per_mm2.max(),min_shared_score=z.pseudotime.min(),max_shared_score=z.pseudotime.max(),
                         zero=int(z.combined_immune_per_mm2.eq(0).sum()),
                         peri_positive_core_zero=int((z.peri_immune_per_mm2.gt(0)&z.core_immune_per_mm2.eq(0)).sum()),
                         graph_fitted=z.composition_dpt_rank.notna().sum(), shared=len(cmp),
                         difference_gt_0_1=int((cmp.pseudotime_rank-cmp.composition_dpt_rank).abs().gt(.1).sum()),
                         prior_dpt_rho=cmp.corr(method='spearman').iloc[0,1],
                         native_core_rank_rho=valid.pseudotime_rank.corr(valid.core_immune_per_mm2,method='spearman')))
    pd.DataFrame(rows).to_csv(args.output/'islet_inflammation_by_donor.csv',index=False)
    audit=['donor','disease','core_area_um2','peri_area_um2','combined_area_um2','n_core','n_peri','n_core_immune','n_peri_immune','n_combined_immune','core_immune_per_mm2','peri_immune_per_mm2','combined_immune_per_mm2','density_eligible','density_order_status','burden','pseudotime','pseudotime_rank','within_donor_density_rank','within_donor_burden','composition_dpt_rank']
    o[audit].to_parquet(args.output/'islet_inflammation_ranks.parquet')
    disease=o.groupby('disease',observed=True).agg(donors=('donor','nunique'), native_islets=('donor','size'),
        assessed_islets=('combined_immune_per_mm2','count'),zero_islets=('combined_immune_per_mm2',lambda x:x.eq(0).sum()),
        min_immune_per_mm2=('combined_immune_per_mm2','min'),median_immune_per_mm2=('combined_immune_per_mm2','median'),
        max_immune_per_mm2=('combined_immune_per_mm2','max'),min_shared_score=('pseudotime','min'),max_shared_score=('pseudotime','max'))
    disease.to_csv(args.output/'islet_inflammation_by_disease.csv')
    state={'version':VERSION,'method':notes,'islets_source_sha256':sha(args.islets),'islets_sha256':sha(out),'sources':source_hashes,
           'manifest_sha256':sha(args.manifest),'donors':len(rows),'native_islets':len(o),'ranked_islets':int(o.pseudotime_rank.notna().sum()),'composition_graph_fitted_islets':int(o.composition_dpt_rank.notna().sum()),'zero_immune_islets':int(o.combined_immune_per_mm2.eq(0).sum())}
    (args.output/'islet_inflammation_method.json').write_text(json.dumps(state,indent=2)+'\n')

def attach_donor_graph(directory):
    """Link independently fitted donor graphs to the native AnnData by donor and islet UUID."""
    directory=Path(directory)
    path=directory/'pancreas_islets.h5ad'
    data=ad.read_h5ad(path)
    method_file=directory/'donor_islet_trajectory_method.json'
    method=json.loads(method_file.read_text())
    for name,expected in method['output_sha256'].items():
        if sha(directory/name)!=expected:
            raise ValueError(f'Incorrect donor trajectory evidence: {name}')
    points=pd.read_parquet(directory/'donor_islet_trajectory_points.parquet').set_index('islet_id',verify_integrity=True)
    if len(points)!=data.n_obs or not data.obs_names.isin(points.index).all():
        raise ValueError('Donor graphs and native islets have different UUIDs')
    points=points.loc[data.obs_names]
    np.testing.assert_array_equal(points.donor.astype(str),data.obs.donor.astype(str))
    np.testing.assert_allclose(points.combined_immune_per_mm2,data.obs.combined_immune_per_mm2,equal_nan=True)
    np.testing.assert_allclose(points.composition_dpt,data.obs.composition_dpt,equal_nan=True)
    data.obs['donor_graph_dpt']=points.graph_dpt.to_numpy()
    data.obs['donor_graph_root']=points.graph_root.to_numpy()
    data.obs['donor_graph_cluster']=pd.Categorical(points.graph_cluster)
    data.obs['donor_graph_exclusion']=pd.Categorical(points.graph_exclusion)
    data.obsm['X_donor_islet_trajectory']=points[['graph_x','graph_y']].to_numpy()
    data.uns['donor_trajectory_method']={k:method[k] for k in ['version','graph','root','coordinates','interpretation']}
    data.obs['inflammation_score']=data.obs.pseudotime
    data.obs['cohort_density_rank']=data.obs.pseudotime_rank
    ad.settings.allow_write_nullable_strings=True
    with tempfile.NamedTemporaryFile(dir=directory,suffix='.h5ad',delete=False) as temp:
        staged=Path(temp.name)
    try:
        data.write_h5ad(staged,compression='gzip')
        np.testing.assert_allclose(read_islet_column(staged,'donor_graph_dpt'),points.graph_dpt,equal_nan=True)
        staged.replace(path)
    finally:
        staged.unlink(missing_ok=True)
    state_file=directory/'islet_inflammation_method.json'
    state=json.loads(state_file.read_text())
    state['islets_sha256']=sha(path)
    state['donor_graph_method_sha256']=sha(method_file)
    state_file.write_text(json.dumps(state,indent=2)+'\n')


def read_islet_column(path,col):
    with h5py.File(path,'r') as h:return np.asarray(read_elem(h['obs'][col]))

def attach_one(job):
    src,dst,islet_path=job
    dst.parent.mkdir(exist_ok=True,parents=True,mode=0o700)
    if src.resolve()==dst.resolve():raise ValueError('Choose a separate staging directory')
    # Existing targets must be complete and match the source islet table; never mutate a published file.
    if dst.exists():raise FileExistsError(dst)
    isl=ad.read_h5ad(islet_path).obs
    if not isl.index.is_unique:raise ValueError('Repeated islet IDs')
    handle=tempfile.NamedTemporaryFile(dir=dst.parent,prefix='.building-',suffix='.h5ad',delete=False);handle.close()
    tmp=Path(handle.name)
    try:
        shutil.copyfile(src,tmp)
        with h5py.File(tmp,'r+') as f:
            oc=f['obs'];us=f['uns']
            if read_elem(us['release_id'])!='io60-phenocycler-20261113-no-epcam':raise ValueError('Unexpected cell typing release')
            prior=read_elem(us['inflammation_order_version']) if 'inflammation_order_version' in us else None
            if prior not in (None,PREVIOUS_VERSION):raise ValueError(f'Unexpected source version: {prior}')
            cats=read_elem(oc['islet_id'])
            # Map categories once; never clip missing/ambiguous memberships into an islet.
            if not isinstance(cats,pd.Categorical):cats=pd.Categorical(cats)
            uniq=pd.Index(cats.categories)
            if not uniq.isin(isl.index).all():raise ValueError(f'{src.name}: unknown islet UUID')
            ix=isl.loc[uniq];codes=cats.codes;good=codes>=0
            lookup=read_elem(oc['donor'])
            donor=np.asarray(lookup).astype(str)[good]
            if not np.array_equal(donor,ix.donor.astype(str).to_numpy()[codes[good]]):raise ValueError('Cell linked to incorrect donor')
            def vector(field,dtype='float32'):
                x=np.full(len(codes),np.nan,dtype=dtype);x[good]=ix[field].to_numpy(dtype=dtype)[codes[good]];return x
            if prior is None:
                for old,new in PREV.items():oc.move(old,new)
            else:
                np.testing.assert_allclose(read_elem(oc['islet_pseudotime_rank']),vector('within_donor_density_rank'),equal_nan=True)
            np.testing.assert_allclose(read_elem(oc['islet_composition_dpt_rank']),vector('composition_dpt_rank'),equal_nan=True)
            np.testing.assert_allclose(read_elem(oc['islet_composition_dpt']),vector('composition_dpt'),equal_nan=True)
            for name,field in CELL_FIELDS.items():write_elem(oc,name,vector(field),dataset_kwargs={'compression':'lzf'})
            burden=np.full(len(codes),-1,dtype='int8');burden[good]=ix.burden.cat.codes.to_numpy()[codes[good]]
            write_elem(oc,'islet_inflammation',pd.Categorical.from_codes(burden,categories=ix.burden.cat.categories,ordered=True),dataset_kwargs={'compression':'lzf'})
            order=list(oc.attrs['column-order'])
            oc.attrs['column-order']=np.asarray(list(dict.fromkeys(order+list(PREV.values())+list(CELL_FIELDS)+['islet_inflammation'])),dtype=h5py.string_dtype())
            write_elem(us,'inflammation_order_version',VERSION)
            assert set(read_elem(oc['islet_inflammation']).categories)=={'Zero','Lower','Upper'}
            for name,field in CELL_FIELDS.items():np.testing.assert_allclose(read_elem(oc[name]),vector(field),equal_nan=True,rtol=0,atol=0)
            assert np.array_equal(read_elem(oc['islet_inflammation']).codes,burden)
            n=len(codes);linked=int(good.sum());ranked=int(np.isfinite(vector('pseudotime_rank')).sum())
        tmp.replace(dst)
        return {'file':str(dst.name),'cells':n,'linked':linked,'ranked':ranked,'source_sha256':sha(src),'sha256':sha(dst)}
    finally:
        tmp.unlink(missing_ok=True)

def attach(args):
    if args.source.resolve()==args.output.resolve():raise ValueError('Do not edit original cells in place')
    args.output.mkdir(exist_ok=True,parents=True,mode=0o700)
    files=[Path('pancreas_umap.h5ad'),Path('pancreas_random.h5ad'),*[Path('pancreas')/f.name for f in sorted((args.source/'pancreas').glob('*.h5ad'))]]
    if not files[2:]:raise ValueError('No full donor sections')
    with h5py.File(args.islets,'r') as h:
        if read_elem(h['uns']['inflammation_order_version'])!=VERSION:raise ValueError('Wrong islet rank version')
    meta=json.loads((args.source/'typing_release.json').read_text())
    sections=pd.read_csv(args.source/'pancreas_manifest.csv',dtype={'donor':str})
    if set(f.stem for f in files[2:])!=set(sections.donor):raise ValueError('Full donor section census mismatch')
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for file,r in zip(files,pool.map(attach_one,[(args.source/f,args.output/f,args.islets) for f in files]),strict=True):
            r['file']=str(file);rows.append(r);print(str(file),r['cells'],r['ranked'],flush=True)
    if sum(r['cells'] for r in rows[2:])!=meta['pancreas_cells']:raise ValueError('Attached cell census disagrees with typing release')
    pd.DataFrame(rows).to_csv(args.output/'islet_inflammation_cell_files.csv',index=False)
    details=json.loads((args.islets.parent/'islet_inflammation_method.json').read_text())
    if details['islets_sha256']!=sha(args.islets):raise ValueError('Islet provenance mismatch')
    meta['islet_inflammation_order']={k:details[k] for k in ['version','native_islets','ranked_islets','composition_graph_fitted_islets','zero_immune_islets','islets_sha256']}
    meta['islet_inflammation_order']['cell_file_manifest']='islet_inflammation_cell_files.csv'
    meta['islet_inflammation_order']['method']='islet_inflammation_method.json'
    meta['islet_inflammation_order'].update(donor_trajectories='donor_islet_trajectory_method.json',
        trajectory_method_sha256=details['donor_graph_method_sha256'], absolute_density='combined_immune_per_mm2',
        shared_score='pseudotime', donor_graph_dpt='donor_graph_dpt')
    meta['fitted_islets']=details['composition_graph_fitted_islets']
    (args.output/'typing_release.json').write_text(json.dumps(meta,indent=2)+'\n')
    stamp(args.output,args.islets)

def stamp(output,islets):
    # Post-computation provenance can be regenerated without rerunning all geometry.
    output,islets=Path(output),Path(islets)
    method_file=islets.parent/'islet_inflammation_method.json'
    record=json.loads(method_file.read_text())
    hashes={name:sha(islets.parent/name) for name in ('islet_inflammation_by_donor.csv','islet_inflammation_by_disease.csv','islet_inflammation_ranks.parquet')}
    code_hash=sha(Path(__file__).resolve())
    record['producer_sha256']=code_hash
    record['audit_sha256']=hashes
    method_file.write_text(json.dumps(record,indent=2)+'\n')
    receipt_file=output/'typing_release.json'
    receipt=json.loads(receipt_file.read_text())
    receipt['islet_inflammation_order'].update({'producer_sha256':code_hash,'method_sha256':sha(method_file),
        'cell_file_manifest_sha256':sha(output/'islet_inflammation_cell_files.csv')})
    receipt_file.write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='cmd',required=True)
    a=subs.add_parser('measure')
    for field in ('islets','features-root','manifest','output','geometry-root','tissue-root'):
        a.add_argument('--'+field,type=Path,required=field not in ('geometry-root','tissue-root'))
    a.add_argument('--workers',type=int,default=3)
    b=subs.add_parser('attach')
    for field in ('islets','source','output'):b.add_argument('--'+field,type=Path,required=True)
    b.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    if args.cmd=='measure':
        measure(args)
        from donor_islet_trajectories import build
        build(args.output/'pancreas_islets.h5ad',args.output)
        attach_donor_graph(args.output)
    else:
        attach(args)

if __name__=='__main__':main()
