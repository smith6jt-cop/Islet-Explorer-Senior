#!/usr/bin/env python3
"""Retrospective fixed spatial family: endocrine domains and conditional checkpoint alignment.

Independent donor endpoints only. Selected current aligned values and exact current
UUID eligibility; IF contributes no endocrine cell identities. The conditional nulls
are statistical expression-domain tests, not evidence of functional immune mechanisms.
"""
from pathlib import Path
import json, hashlib, sys, argparse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree
from scipy.stats import rankdata
import shapely
from shapely.geometry import shape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path('/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3')
OUT=REPO/'results/v3/islet_immune_reassessment_20261006'
TAB=OUT/'tables';FIG=OUT/'figures'
ASS=Path('/data/io60-calibration/pancreas_analysis_assignments_single_exclusions_20261005/manifest.json')
sys.path.insert(0,str(REPO/'core/phenocycler/scripts/typing'))
from build_islet_geometry import islet_export, _manifest
SEED=20261006;NNULL=199;FAMILY='spatial_v2_20261006'
MARKERS=['HLA_A','PD_L1','PD_1','Granzyme_B','CD107a','CD31','B3TUBB','Collagen_IV','E_cadherin','EpCAM']
ACOLS=['cell_id','cell_class','refined_class','immune_class','analysis_eligible','endocrine_object_class','endocrine_object_id','hormone_measurement_state','X_centroid','Y_centroid','disease_status']

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for z in iter(lambda:f.read(8*1024*1024),b''):h.update(z)
    return h.hexdigest()

def rankz(v):
    v=np.asarray(v,float);out=np.full(len(v),np.nan);ok=np.isfinite(v)
    if ok.sum()>1:
        q=(rankdata(v[ok],method='average')-.5)/ok.sum();sd=q.std()
        if sd>1e-10:out[ok]=(q-q.mean())/sd
        else:out[ok]=0
    return out

def geometry(native):
    px,st=_manifest();p=islet_export(st[native]);pixel=float(px[native]);d={}
    for f in json.loads(p.read_text())['features']:
        if f.get('properties',{}).get('size_class')!='Islet' or not f.get('geometry'):continue
        g=shape(f['geometry']);g=g if g.is_valid else g.buffer(0)
        if not g.is_empty:d[str(f['id'])]=shapely.affinity.scale(g,xfact=pixel,yfact=pixel,origin=(0,0))
    return p,d

def center(v,g):
    v=np.asarray(v,float)
    if v.ndim==1:v=v[:,None];one=True
    else:one=False
    n=np.bincount(g);out=v.copy()
    for j in range(v.shape[1]):out[:,j]-=(np.bincount(g,weights=v[:,j])/n)[g]
    return out[:,0] if one else out

def design(v,g=None):
    x=np.asarray(v,float)
    if g is not None:x=center(x,g)
    sd=x.std(axis=0);keep=sd>1e-8;x=x[:,keep]/sd[keep]
    return x,keep

def solver(x,w):
    sw=np.sqrt(w);pinv=np.linalg.pinv(x*sw[:,None],rcond=1e-9)
    return lambda y: pinv@(y*sw)

def neighbor_features(xy,tmask,imask,p,endo,radius,k=12):
    e=xy[endo];tree=cKDTree(xy)
    dist,ix=tree.query(e,k=min(k+1,len(xy)),workers=2)
    ix=ix[:,1:];dist=dist[:,1:]
    local=dist<=radius
    density=tree.query_ball_point(e,radius,return_length=True,workers=2)-1
    immune_count=cKDTree(xy[imask]).query_ball_point(e,radius,return_length=True,workers=2) if imask.any() else np.zeros(len(e))
    tis=np.where(tmask)[0];tt=cKDTree(xy[tmask]);neighbors=tt.query_ball_point(e,radius,workers=2)
    pd=rankz(p['PD_1'][tmask]);gb=rankz(p['Granzyme_B'][tmask]);cd=rankz(p['CD107a'][tmask]);cyto=(gb+cd)/2
    tc=np.asarray([len(v) for v in neighbors]);tp=np.full(len(e),np.nan);cy=np.full(len(e),np.nan)
    for i,ids in enumerate(neighbors):
        ids=np.asarray(ids,int)
        if not len(ids):tp[i]=0;cy[i]=0;continue
        d=np.linalg.norm(xy[tis[ids]]-e[i],axis=1);w=np.exp(-.5*(d/(radius/2))**2)
        good=np.isfinite(pd[ids]);goodc=np.isfinite(cyto[ids])
        if good.any():tp[i]=np.average(pd[ids][good],weights=w[good])
        if goodc.any():cy[i]=np.average(cyto[ids][goodc],weights=w[goodc])
    st=[]
    for m in ['CD31','B3TUBB','Collagen_IV']:
        v=rankz(p[m]);valid=local & np.isfinite(v[ix]);n=valid.sum(axis=1)
        st.append(np.divide(np.where(valid,v[ix],0).sum(axis=1),n,out=np.full(len(e),np.nan),where=n>0))
    return dict(tcount=tc,pd1=tp,cyto=cy,immune_count=immune_count,density=density,struct=np.asarray(st).T)

def h3_graph(f,geos,degree,maxdist):
    groups=[];edges_i=[];edges_j=[];edge_group=[];cell_group=np.full(len(f),-1,int);retained=[]
    for islet,z in f.groupby('islet'):
        if len(z)<40:continue
        for subtype,s in z.groupby('subtype'):
            if len(s)<10:continue
            ids=s.index.to_numpy(int);xy=s[['x','y']].to_numpy(float)
            dd,nn=cKDTree(xy).query(xy,k=min(degree+1,len(xy)))
            src=np.repeat(np.arange(len(xy)),nn.shape[1]-1);dst=nn[:,1:].ravel();d=dd[:,1:].ravel()
            ok=(src<dst)&(d<=maxdist)
            # Keep undirected unique edges; include reciprocal selections once.
            und=np.sort(np.column_stack([np.repeat(np.arange(len(xy)),nn.shape[1]-1),nn[:,1:].ravel()]),axis=1)
            und=np.unique(und[d<=maxdist],axis=0)
            if not len(und):continue
            lines=shapely.linestrings(xy[und]);inside=shapely.covers(geos[islet],lines);und=und[inside]
            if len(und)<10:continue
            gid=len(groups);groups.append(dict(islet=islet,subtype=subtype,ids=ids,n=len(ids),edges=len(und)))
            cell_group[ids]=gid;edges_i.extend(ids[und[:,0]]);edges_j.extend(ids[und[:,1]]);edge_group.extend([gid]*len(und));retained.extend(ids)
    return groups,np.asarray(edges_i,int),np.asarray(edges_j,int),np.asarray(edge_group,int),cell_group,np.unique(retained)

def h3(f,geos,degree,maxdist,rng):
    groups,ei,ej,eg,cg,retained=h3_graph(f,geos,degree,maxdist)
    invalid=dict(donor_effect=np.nan,status='unavailable',reason='fewer than5 exchangeable islets',n_cells=len(retained),n_islets=len(set(g['islet'] for g in groups)),n_groups=len(groups),n_edges=len(ei))
    if invalid['n_islets']<5:return invalid,[],None
    # Remove graph-ineligible cells from the nuisance fit; exact old->new indices.
    sub=f.loc[retained].copy().reset_index(drop=True);remap=np.full(len(f),-1,int);remap[retained]=np.arange(len(retained))
    ei=remap[ei];ej=remap[ej];cg=cg[retained]
    y=center(sub.hla.to_numpy(),cg)
    x,_=design(sub[['pos','pos2','nuc','area','tcount','immune','density','pd1','cyto','vascular','neural','collagen','ecad','epcam']].to_numpy(),cg)
    islet_index=pd.factorize(sub.islet,sort=True)[0];w=1/np.bincount(islet_index)[islet_index]
    fit=solver(x,w);yh=x@fit(y);res=y-yh
    members=[np.where(cg==i)[0] for i in range(len(groups))]
    gn=np.asarray([len(v) for v in members]);ge=np.bincount(eg,minlength=len(groups));g_islet=pd.factorize(np.asarray([g['islet'] for g in groups]),sort=True)[0]
    def metric(rr):
        ss=np.bincount(cg,weights=rr*rr,minlength=len(groups));prod=np.bincount(eg,weights=rr[ei]*rr[ej],minlength=len(groups))
        return np.divide(gn*prod,ge*ss,out=np.full(len(groups),np.nan),where=(ss>1e-8)&(ge>0))
    obs=metric(res);null=np.full((NNULL,len(groups)),np.nan)
    for b in range(NNULL):
        rp=res.copy()
        for ids in members:rp[ids]=rng.permutation(res[ids])
        yy=yh+rp;rr=yy-x@fit(yy);null[b]=metric(rr)
    avg=np.nanmean(null,axis=0);sd=np.nanstd(null,axis=0);exchangeable=np.isfinite(obs)&np.isfinite(avg)&(sd>1e-8)
    available_islets=np.unique(g_islet[exchangeable]);rows=[]
    for ii in available_islets:
        keep=(g_islet==ii)&exchangeable;ob=float(np.mean(obs[keep]));nu=null[:,keep].mean(axis=1);eff=ob-float(nu.mean())
        rows.append(dict(islet=groups[np.where(keep)[0][0]]['islet'],observed=ob,null_mean=float(nu.mean()),null_sd=float(nu.std()),effect=eff,n_cells=int(gn[keep].sum()),n_edges=int(ge[keep].sum()),n_subtypes=int(keep.sum())))
    effect=float(np.mean([r['effect'] for r in rows])) if len(rows)>=5 else np.nan
    out=dict(donor_effect=effect,status='measured' if np.isfinite(effect) else 'unavailable',reason='' if np.isfinite(effect) else 'fewer than5 variable null-supported islets',n_cells=len(sub),n_islets=len(rows),n_groups=int(exchangeable.sum()),n_edges=len(ei),nuisance_rank=int(np.linalg.matrix_rank(x)),null_replicates=NNULL,residual_variance=float(res.var()))
    if len(rows):
        ob=np.mean([r['observed'] for r in rows]);nu=np.asarray([np.mean([null[b,(g_islet==i)&exchangeable].mean() for i in available_islets]) for b in range(NNULL)])
        out['conditional_null_p']=(1+np.sum(np.abs(nu-nu.mean())>=abs(ob-nu.mean())-1e-12))/(NNULL+1)
    display=sub.copy();display['hla_residual']=res
    return out,rows,display

def h4(f,rng):
    keep=(f.tcount>0)&np.isfinite(f.pdl1)&np.isfinite(f.hla)&np.isfinite(f.pd1)&np.isfinite(f.cyto)
    z=f.loc[keep].copy().reset_index(drop=True)
    # Each block must have measured within-islet variation; not an independent cell replicate.
    sizes=z.groupby('islet').size();z=z.loc[z.islet.isin(sizes[sizes>=5].index)].reset_index(drop=True)
    block=pd.factorize(z.islet,sort=True)[0];nb=len(np.unique(block))
    out=dict(donor_effect=np.nan,status='unavailable',reason='requires8 independent islet blocks and80 cells',n_cells=len(z),n_islets=nb)
    if len(z)<80 or nb<8:return out,[]
    if z.pdl1.nunique()<2:out['reason']='measured PD-L1 response is constant';return out,[]
    h=rankz(z.hla);p=rankz(z.pd1);c=rankz(z.cyto);y=rankz(z.pdl1)
    nuisance=z[['pos','pos2','nuc','area','tcount','immune','density','vascular','neural','collagen','ecad','epcam']].to_numpy(float)
    # Subtype dummies; islet fixed effects removed through weighted within transformation.
    typ=pd.get_dummies(z.subtype,drop_first=True,dtype=float).to_numpy()
    reduced=np.column_stack([h,p,c,h*p,h*c,p*c,nuisance,typ]);triple=h*p*c
    xr,k=design(reduced,block);yy=center(y,block);target=center(triple,block);sd=target.std()
    if not np.isfinite(sd) or sd<1e-8:out['reason']='interaction does not vary';return out,[]
    target/=sd;w=1/np.bincount(block)[block];fr=solver(xr,w);t_res=target-xr@fr(target)
    if np.average(t_res*t_res,weights=w)<1e-8:out['reason']='interaction is aliased with pairwise/nuisance terms';return out,[]
    xf=np.column_stack([xr,target]);ff=solver(xf,w)
    if np.linalg.matrix_rank(xf)<xf.shape[1]:out['reason']='rank-deficient full model';return out,[]
    coef=ff(yy);pred=xr@fr(yy);res=yy-pred;null=[]
    for b in range(NNULL):
        signs=rng.choice([-1,1],nb);wild=center(pred+res*signs[block],block);null.append(float(ff(wild)[-1]))
    null=np.asarray(null);obs=float(coef[-1]);baseline=float(null.mean())
    out.update(donor_effect=obs-baseline,observed_coefficient=obs,null_mean=baseline,null_sd=float(null.std()),conditional_null_p=(1+np.sum(np.abs(null-baseline)>=abs(obs-baseline)-1e-12))/(NNULL+1),status='measured',reason='',nuisance_rank=int(np.linalg.matrix_rank(xr)),full_rank=int(np.linalg.matrix_rank(xf)),null_replicates=NNULL,interaction_partial_variance=float(np.average(t_res*t_res,weights=w)))
    original_sd=center(reduced,block).std(axis=0);rawcoef=np.zeros(reduced.shape[1]);rawcoef[k]=coef[:-1]/original_sd[k];rawtriple=coef[-1]/sd
    for hh,cc,tag in [(-1,-1,'lowHLA_lowCyto'),(1,-1,'highHLA_lowCyto'),(-1,1,'lowHLA_highCyto'),(1,1,'highHLA_highCyto')]:
        out['PD1_simple_slope_'+tag]=float(rawcoef[1]+rawcoef[3]*hh+rawcoef[5]*cc+rawtriple*hh*cc)
        out['PD1_pairwise_slope_'+tag]=float(rawcoef[1]+rawcoef[3]*hh+rawcoef[5]*cc)
    out['PD1_slope_high_high_minus_low_low']=out['PD1_simple_slope_highHLA_highCyto']-out['PD1_simple_slope_lowHLA_lowCyto']
    return out,[dict(islet=i,n_cells=int((z.islet==i).sum())) for i in sorted(z.islet.unique())]

def work(rec):
    donor=str(rec['donor_id']);native=rec['donor_id_native'];rng=np.random.default_rng(SEED+int(donor))
    assert sha(rec['path'])==rec['sha256'],'Current mask hash mismatch'
    a=pd.read_parquet(rec['path'],columns=ACOLS)
    assert a.cell_id.is_unique
    source=Path(rec['native_source_path']);schema=set(pq.ParquetFile(source).schema.names)
    sc=['object_id','X_centroid','Y_centroid','nucleus_area_um2','fundamentals_cell_area_um2']+['redsea_aligned_'+m for m in MARKERS if 'redsea_aligned_'+m in schema]
    s=pd.read_parquet(source,columns=sc);assert s.object_id.is_unique
    ix=pd.Index(s.object_id.astype(str)).get_indexer(a.cell_id.astype(str));assert (ix>=0).all(),'Exact UUID join required'
    s=s.iloc[ix].reset_index(drop=True);a=a.reset_index(drop=True)
    assert np.nanmax(abs(s[['X_centroid','Y_centroid']].to_numpy(float)-a[['X_centroid','Y_centroid']].to_numpy(float)))<1e-3
    eligible=a.analysis_eligible.fillna(False).to_numpy(bool)
    assert eligible.sum()==rec['typing_analysis_cells']
    a=a.loc[eligible].reset_index(drop=True);s=s.loc[eligible].reset_index(drop=True)
    xy=a[['X_centroid','Y_centroid']].to_numpy(float);assert np.isfinite(xy).all()
    p={m:s['redsea_aligned_'+m].to_numpy(float) if 'redsea_aligned_'+m in s else np.full(len(s),np.nan) for m in MARKERS}
    nuc=s.nucleus_area_um2.to_numpy(float);area=s.fundamentals_cell_area_um2.to_numpy(float)
    gpath,geos=geometry(native)
    disease=str(a.disease_status.dropna().mode().iloc[0]).replace('AAb+','AAb')
    base=dict(donor_id=donor,disease_status=disease,family=FAMILY,source_eligible_cells=len(a),native_T_cells=int(a.immune_class.eq('Imm-T-cell').sum()))
    provenance=dict(donor_id=donor,assignment_path=rec['path'],assignment_sha256=rec['sha256'],protein_path=str(source),protein_sha256=sha(source),geometry_path=str(gpath),geometry_sha256=sha(gpath),measurements=['redsea_aligned_'+m for m in MARKERS],null_replicates=NNULL)
    endpoints=[];sens=[];isrows=[];feature_file=TAB/f'spatial_H3_H4_endocrine_features_donor_{donor}.parquet'
    if not base['native_T_cells']:
        for h in ['H3','H4']:endpoints.append(dict(**base,hypothesis=h,donor_effect=np.nan,status='unavailable',reason='native T lineage scope unavailable; no zero imputation',n_cells=0,n_islets=0))
        return endpoints,sens,isrows,provenance
    mean_diameter=float(np.nanmean(2*np.sqrt(area[np.isfinite(area)&(area>0)]/np.pi)))
    displays=[]
    for identity in ['native','refined']:
        label=a.cell_class if identity=='native' else a.refined_class
        endo=label.isin(['Beta','Alpha','Delta']).to_numpy() & a.endocrine_object_class.eq('Islet').fillna(False).to_numpy(bool) & a.hormone_measurement_state.eq('usable_for_typing').fillna(False).to_numpy(bool)
        tmask=a.immune_class.eq('Imm-T-cell').fillna(False).to_numpy(bool)
        if identity=='refined':tmask = tmask & a.refined_class.eq('Immune').fillna(False).to_numpy(bool)
        imask=label.eq('Immune').fillna(False).to_numpy(bool)
        ei=np.where(endo)[0]
        if len(ei)<80 or not tmask.any():
            for h in ['H3','H4']:
                row=dict(**base,hypothesis=h,sensitivity=identity,donor_effect=np.nan,status='unavailable',reason='insufficient typed endocrine/T support',n_cells=len(ei),n_islets=0)
                if identity=='native':endpoints.append(row)
                else:sens.append(row)
            continue
        ids=a.iloc[ei].endocrine_object_id.astype(str).to_numpy();assert set(ids).issubset(geos)
        subtype=label.iloc[ei].astype(str).to_numpy();hla=np.full(len(ei),np.nan);pdl=np.full(len(ei),np.nan)
        for typ in ['Beta','Alpha','Delta']:
            target=np.where(subtype==typ)[0]
            hla[target]=rankz(p['HLA_A'][ei[target]]);pdl[target]=rankz(p['PD_L1'][ei[target]])
        depth=np.empty(len(ei));size=np.empty(len(ei))
        for uid in np.unique(ids):
            selected=np.where(ids==uid)[0];points=shapely.points(xy[ei[selected]]);assert shapely.covers(geos[uid],points).all(),'Assignment must agree with true outline'
            depth[selected]=shapely.distance(points,shapely.boundary(geos[uid]));size[selected]=geos[uid].area
        for radius in ([50,25,75] if identity=='native' else [50]):
            nf=neighbor_features(xy,tmask,imask,p,endo,radius)
            f=pd.DataFrame(dict(islet=ids,subtype=subtype,x=xy[ei,0],y=xy[ei,1],hla=hla,pdl1=pdl,pos=depth/np.sqrt(size/np.pi),nuc=np.log1p(nuc[ei]),area=np.log1p(area[ei]),tcount=np.log1p(nf['tcount']),immune=np.log1p(nf['immune_count']),density=np.log1p(nf['density']),pd1=nf['pd1'],cyto=nf['cyto'],vascular=nf['struct'][:,0],neural=nf['struct'][:,1],collagen=nf['struct'][:,2],ecad=rankz(p['E_cadherin'][ei]),epcam=rankz(p['EpCAM'][ei])))
            f['pos2']=f.pos**2;f['cell_id']=a.iloc[ei].cell_id.to_numpy();f['T_neighbor_count']=nf['tcount'];f['signed_dist_islet_um']=-depth;f['islet_area_um2']=size
            variables=['hla','pos','pos2','nuc','area','tcount','immune','density','pd1','cyto','vascular','neural','collagen','ecad','epcam']
            ff=f.loc[np.isfinite(f[variables].to_numpy(float)).all(axis=1)].reset_index(drop=True)
            if identity=='native' and radius==50:
                f['donor_id']=donor;f['disease_status']=disease;f.to_parquet(feature_file,index=False)
            for degree in ([4,2,8] if radius==50 and identity=='native' else [4]):
                hr,ir,display=h3(ff,geos,degree,2*mean_diameter,rng)
                row=dict(**base,**hr,hypothesis='H3',sensitivity=f'{identity}_radius{radius}_k{degree}',neighborhood_um=radius,graph_k=degree,edge_cap_um=2*mean_diameter)
                if identity=='native' and radius==50 and degree==4:
                    endpoints.append(row)
                    if display is not None:display['donor_id']=donor;display['disease_status']=disease;display.to_parquet(TAB/f'spatial_H3_residual_cells_donor_{donor}.parquet',index=False)
                    isrows.extend([dict(**base,**r,hypothesis='H3') for r in ir])
                else:sens.append(row)
            hr,ir=h4(ff,rng);row=dict(**base,**hr,hypothesis='H4',sensitivity=f'{identity}_radius{radius}',neighborhood_um=radius)
            if identity=='native' and radius==50:endpoints.append(row);isrows.extend([dict(**base,**r,hypothesis='H4') for r in ir])
            else:sens.append(row)
    print(f"H3/H4 donor {donor}: "+'; '.join(f"{z['hypothesis']} {z['status']} nIslet={z['n_islets']} effect={z['donor_effect']}" for z in endpoints),flush=True)
    return endpoints,sens,isrows,provenance

def plots(e):
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    colors={'ND':'#3a6e94','AAb':'#bd8a29','T1D':'#a74151'}
    for a,h,title in zip(ax,['H3','H4'],['HLA-A domain excess beyond contact','Conditional checkpoint alignment']):
        z=e.loc[e.hypothesis.eq(h)]
        for x,stage in enumerate(['ND','AAb','T1D']):
            s=z.loc[z.disease_status.eq(stage)&z.donor_effect.notna()].sort_values('donor_id');j=np.linspace(-.18,.18,len(s))
            for dx,(_,r) in zip(j,s.iterrows()):a.scatter(x+dx,r.donor_effect,color=colors[stage],s=24);a.annotate(r.donor_id,(x+dx,r.donor_effect),xytext=(2,2),textcoords='offset points',fontsize=6)
            a.text(x,.01,f"{len(s)}/{int(z.disease_status.eq(stage).sum())} supported",transform=a.get_xaxis_transform(),ha='center',fontsize=8)
        a.axhline(0,lw=.7,color='gray',ls='--');a.set_xticks(range(3),['ND','AAb','T1D']);a.set_title(title,fontsize=11);a.spines[['top','right']].set_visible(False)
        a.set_ylabel('Residual homotypic graph excess' if h=='H3' else 'Standardized triple interaction, null centered')
    fig.suptitle('Independent donor endpoints | selected current processed protein',fontsize=12)
    for ext in ['png','svg']:fig.savefig(FIG/f'spatial_protein_H3_H4_summary.{ext}',dpi=300)
    plt.close(fig)

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--workers',type=int,default=2);pa.add_argument('--donor',action='append');args=pa.parse_args()
    plan=OUT/'provenance/spatial_hypothesis_plan_v2.json';assert plan.exists(),'Fixed protocol required before endpoint calculation'
    records=json.loads(ASS.read_text())['donors']
    if args.donor:records=[r for r in records if str(r['donor_id']) in args.donor]
    all_e=[];all_s=[];all_i=[];prov=[]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for e,s,i,p in ex.map(work,records):
            all_e.extend(e);all_s.extend(s);all_i.extend(i);prov.append(p)
            pd.DataFrame(all_e).to_csv(TAB/'spatial_H3_H4_donor_endpoints.partial.csv',index=False)
    e=pd.DataFrame(all_e);s=pd.DataFrame(all_s);i=pd.DataFrame(all_i)
    e.to_csv(TAB/'spatial_H3_H4_donor_endpoints.csv',index=False);s.to_csv(TAB/'spatial_H3_H4_sensitivities.csv',index=False);i.to_csv(TAB/'spatial_H3_H4_islet_null_QC.csv',index=False)
    plots(e)
    (TAB/'spatial_H3_H4_manifest.json').write_text(json.dumps(dict(family=FAMILY,script_sha256=sha(__file__),plan_sha256=sha(plan),donors=prov,null='H3 residual exchangeability approximation; H4 block-wild pairwise baseline',primary_tests='Owned by six-family donor combiner; conditional cell/null p values not disease inference'),indent=2))

if __name__=='__main__':main()
