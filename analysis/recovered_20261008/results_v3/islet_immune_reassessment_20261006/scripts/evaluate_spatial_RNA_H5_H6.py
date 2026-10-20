#!/usr/bin/env python3
"""Frozen RNA spatial H5/H6: higher-order organizer motifs and linked node architecture.
Released native identity authority; raw markers are diagnostic only. No IF hormone identity.
"""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,json,hashlib,importlib.util,os,sys
import numpy as np,pandas as pd,h5py
import shapely
from scipy.spatial import cKDTree
from scipy.stats import rankdata
ROOT=Path('/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3')
R=ROOT/'results/v3/islet_immune_reassessment_20261006';T=R/'tables';F=R/'figures'
MAN=ROOT/'results/v3/complete_cell_typing_paper_style_20261005/source_tables/xenium_immune_subpanels_20261005/manifest.json'
SEED=20261006
LEGACY=['112','6347','6533','6579','6580']
BGENES=['MS4A1','CD79A','CD79B','CD19','CD22'];TGENES=['CD3E','CD3G','TRAC']
MEM=['TCF7','LEF1'];CYTO=['TOX','GZMB','PRF1','NKG7'];ORG=['CCL19','CCL21']
CONFIGS=[('primary',10,50.,'released'),('depth50',50,50.,'released'),('depth100',100,50.,'released'),
 ('depth200',200,50.,'released'),('radius25',10,25.,'released'),('radius75',10,75.,'released')]
sys.path.insert(0,str(ROOT/'core/phenocycler'))
sys.path.insert(0,str(ROOT/'core/phenocycler/scripts/restain'))
import pheno_annotations as PA
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def key(x):return str(int(str(x)))
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v)for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,(np.integer,np.floating,np.bool_)):x=x.item()
 if isinstance(x,float)and not np.isfinite(x):return None
 return x
def jsonwrite(p,x):Path(p).write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n')
def reader():
 p=ROOT/'analysis/xenium/scripts/regenerate_current_rna_umaps.py';s=importlib.util.spec_from_file_location('RNAreader',p)
 m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def program(frame,genes,depth_col='raw_library_counts'):
 if not all('RNA_'+g in frame for g in genes):return None
 depth=frame[depth_col].to_numpy(float)
 return np.mean(np.column_stack([np.log1p(1e4*frame['RNA_'+g].to_numpy(float)/depth)for g in genes]),axis=1)
def zscore(x):
 x=np.asarray(x,float);sd=x.std()
 return (x-x.mean())/sd if sd>1e-12 else np.zeros(len(x))
def numerical_columns(frame):
 out=[];names=[]
 for name in frame:
  x=frame[name].to_numpy(float)
  if not np.isfinite(x).all():continue
  if x.std()>1e-12:out.append(zscore(x));names.append(name)
 return np.column_stack(out)if out else np.empty((len(frame),0)),names
def weighted_center(x,blocks,w):
 out=np.asarray(x,float).copy()
 if out.ndim==1:out=out[:,None];one=True
 else:one=False
 for b in np.unique(blocks):
  k=blocks==b;out[k]-=np.average(out[k],axis=0,weights=w[k])
 return out[:,0]if one else out
def block_wild(X0,I,y,blocks,w,n_null=199):
 X=weighted_center(X0,blocks,w);I=weighted_center(I,blocks,w);y=weighted_center(y,blocks,w)
 sw=np.sqrt(w);A=X*sw[:,None];rank=int(np.linalg.matrix_rank(A))
 fittedX=np.linalg.lstsq(A,I*sw,rcond=None)[0];rI=I-X@fittedX
 fittedY=np.linalg.lstsq(A,y*sw,rcond=None)[0];rY=y-X@fittedY
 den=float(np.sum(w*rI*rI))
 if den<1e-10:return None
 beta=float(np.sum(w*rI*rY)/den)
 ub=np.unique(blocks)
 cc=np.array([np.sum(w[blocks==b]*rI[blocks==b]*rY[blocks==b])/den for b in ub])
 full=np.array([np.sum(w[blocks==b]*rI[blocks==b]*(rY[blocks==b]-beta*rI[blocks==b]))/den for b in ub])
 if len(ub)<=12:
  numbers=np.arange(2**len(ub),dtype=np.uint32);signs=2*((numbers[:,None]>>np.arange(len(ub)))&1).astype(float)-1
  exact=True
 else:signs=np.random.default_rng(SEED).choice([-1.,1.],size=(n_null,len(ub)));exact=False
 null=signs@cc;ci=beta+signs@full
 p=float(np.mean(null>=beta-1e-12))if exact else float((1+(null>=beta).sum())/(len(null)+1))
 lo,hi=np.quantile(ci,[.025,.975])
 return {'donor_effect':beta,'within_donor_p_positive':p,'CI95_low':float(lo),'CI95_high':float(hi),
  'null_mean':float(null.mean()),'null_sd':float(null.std(ddof=1)),'n_null':len(null),'exact_block_sign_enumeration':exact,
  'n_islet_blocks':len(ub),'n_baseline_rank_after_islet_centering':rank,
  'interaction_residual_sum_squares':den,'pairwise_baseline_R2_within_islet':float(1-np.sum(w*rY*rY)/np.sum(w*y*y))if np.sum(w*y*y)>0 else None,
  'null_policy':'Restricted pairwise-only baseline + whole-islet Rademacher wild residual signs; within-islet overlapping-neighborhood covariance and all covariates/positions remain fixed',
  'null_limit':'Wild residual test assumes independent islet blocks and approximately symmetric conditional block residuals; donor comparison remains independent donor inference'}

def impute_area_QC(frame):
 # Missing area is a model-QC covariate, never a cell-identity failure.
 frame=frame.copy()
 for field in ['nucleus_area_um2','cell_area_um2']:
  x=pd.to_numeric(frame[field],errors='coerce')if field in frame else pd.Series(np.nan,index=frame.index)
  bad=~np.isfinite(x)|x.le(0)
  frame[field+'_missing_QC']=bad.to_numpy(float)
  measured=x[~bad];fill=float(measured.median())if len(measured)else 0.
  frame[field]=x.where(~bad,fill)
 return frame

def h5_one(d,geometry,state,membership,report,footprint,config):
 name,depth,radius,route=config
 state=impute_area_QC(state)
 out={'hypothesis':'H5','donor_id':d,'disease_status':state.disease_status.iloc[0],'configuration':name,'depth_threshold':depth,'neighbor_radius_um':radius,
  'identity_route':route,'H5':np.nan,'status':'untestable','reason':'','n_supported_B_context':0,'n_supported_T_context':0}
 # Organizers remain inside actual masks. Physical cell neighbors may cross an islet boundary; no expanded mask exists.
 valid_geom=geometry[geometry.geometry_wkb.map(lambda x:footprint.covers(shapely.from_wkb(x)))].copy()
 allowed=set(valid_geom.islet_id);inside=state.native_islet_id.isin(allowed)
 k=state.raw_library_counts.ge(depth)&state.nucleus_count.gt(0)
 organ=k&inside&state.native_lineage.isin(['Stromal_Fibroblast_Stellate','Vascular'])
 # Identity is the current pancreas producer's released call, without a new detection gate.
 bs=k&state.native_lineage.eq('Immune_B_Plasma');ts=k&state.native_lineage.eq('Immune_T')
 out['identity_authority']='panc_20261003_bundle/release native lineage; raw marker support is diagnostic only'
 out['n_missing_nuclear_area_QC']=int(state.loc[k,'nucleus_area_um2_missing_QC'].sum())
 out['n_missing_cell_area_QC']=int(state.loc[k,'cell_area_um2_missing_QC'].sum())
 org=state[organ].copy();b=state[bs].copy();tt=state[ts].copy()
 out.update(n_organizers=len(org),n_supported_B_all=len(b),n_supported_T_all=len(tt),n_whole_masks_inside_RNA_capture=len(allowed))
 if len(org)<100:out['reason']='fewer_than100_inside_Islet_organizer_cells';return out,[],None
 score=program(org,ORG);bv=program(b,['LTB'])
 if score is None or bv is None:out['reason']='required_gene_not_measured';return out,[],None
 positive_org=(org[['RNA_CCL19','RNA_CCL21']].sum(axis=1)>0).to_numpy()
 if positive_org.sum()<20:out['reason']='fewer_than20_CCL19_or_CCL21_positive_organizers';return out,[],None
 vendor=Path(report['source_acquisition_path'])/'cells.parquet'
 full=pd.read_parquet(vendor,columns=['cell_id','x_centroid','y_centroid','nucleus_count','transcript_counts'])
 full.cell_id=full.cell_id.astype(str)
 if full.cell_id.duplicated().any():raise ValueError('Vendor duplicate barcode')
 calls=pd.read_parquet(report['native_calls_path'])
 calls=calls[calls.section_key.astype(str).eq(str(report['section_id']))].copy();calls.cell_id=calls.cell_id.astype(str)
 if calls.cell_id.duplicated().any()or not set(calls.cell_id)<=set(full.cell_id):raise ValueError('Native/vendor QC keys incomplete')
 admissible=calls.nucleus_count.gt(0)&calls.depth.ge(depth)&calls.n_genes_detected.ge(5)
 full=full[full.cell_id.isin(set(calls.loc[admissible,'cell_id']))].copy()
 fullxy=full[['x_centroid','y_centroid']].to_numpy(float);oxy=org[['x_centroid_um','y_centroid_um']].to_numpy(float)
 bt=cKDTree(b[['x_centroid_um','y_centroid_um']].to_numpy(float))if len(b)else None
 tt_tree=cKDTree(tt[['x_centroid_um','y_centroid_um']].to_numpy(float))if len(tt)else None
 all_tree=cKDTree(fullxy)
 bn=bt.query_ball_point(oxy,radius,workers=1)if bt is not None else [[]for _ in org.index]
 tn=tt_tree.query_ball_point(oxy,radius,workers=1)if tt_tree is not None else [[]for _ in org.index]
 deg=all_tree.query_ball_point(oxy,radius,return_length=True,workers=1)-1
 wide=all_tree.query_ball_point(oxy,100.,return_length=True,workers=1)-1
 area=np.array([footprint.intersection(shapely.Point(x,y).buffer(100.,quad_segs=8)).area for x,y in oxy])
 density=np.divide(wide,area,out=np.zeros(len(wide),float),where=area>0)
 bi=np.unique(np.concatenate([np.asarray(x,dtype=int)for x in bn]))if any(len(x)for x in bn)else np.array([],int)
 ti=np.unique(np.concatenate([np.asarray(x,dtype=int)for x in tn]))if any(len(x)for x in tn)else np.array([],int)
 out.update(n_supported_B_context=len(bi),n_supported_T_context=len(ti),n_B_LTB_positive_context=int((b.iloc[bi].RNA_LTB>0).sum())if len(bi)else 0,
  n_CCL_positive_organizers=int(positive_org.sum()),n_islets_organizer=int(org.native_islet_id.nunique()))
 if len(bi)<20 or len(ti)<20:out['reason']='fewer_than20_unique_supported_B_or_T_physical_neighbors';return out,[],None
 if out['n_B_LTB_positive_context']<10:out['reason']='fewer_than10_LTB_positive_supported_B_neighbors';return out,[],None
 bcount=np.array([len(x)for x in bn]);tcount=np.array([len(x)for x in tn])
 exposure=np.array([bv[x].sum()for x in bn])/np.maximum(deg,1)
 y=tcount/np.maximum(deg,1);oz=zscore(score);bz=zscore(exposure)
 kfit=(deg>=5)&np.isfinite(density)&(area>0)
 if kfit.sum()<100:out['reason']='fewer_than100_organizers_with_local_all_cell_support';return out,[],None
 frame=pd.DataFrame({'CCL_state':oz,'B_LTB_exposure':bz,'B_neighbor_count':np.log1p(bcount),'local_all_QC_degree':np.log1p(np.maximum(deg,0)),
  'local_all_QC_density100um':np.log1p(density*1e6),'raw_library_depth':np.log1p(org.raw_library_counts.to_numpy()),
  'nucleus_area':np.log1p(org.nucleus_area_um2.to_numpy()),'cell_area':np.log1p(org.cell_area_um2.to_numpy()),
  'nucleus_area_missing_QC':org.nucleus_area_um2_missing_QC.to_numpy(),'cell_area_missing_QC':org.cell_area_um2_missing_QC.to_numpy(),
  'continuous_signed_islet_distance':org.signed_distance_to_available_IF_islet_um.to_numpy(),
  'organizer_vascular':org.native_lineage.eq('Vascular').to_numpy(float)})
 controlcols=[]
 for field in ['genomic_control_counts','control_probe_counts','control_codeword_counts']:
  if field in org:
   col=field+'_per1000_gene_counts';frame[col]=org[field].to_numpy(float)/org.raw_library_counts.to_numpy(float)*1000;controlcols.append(field)
 out['control_fields_measured']=','.join(controlcols)
 if 'genomic_control_counts'not in controlcols:out['reason']='gDNA_control_unavailable';return out,[],None
 frame=frame.loc[kfit];X,names=numerical_columns(frame);blocks=org.native_islet_id.to_numpy()[kfit]
 if len(np.unique(blocks))<8:out['reason']='fewer_than8_independent_islet_blocks';return out,[],None
 if np.unique(y[kfit]).size<3 or np.count_nonzero(tcount[kfit])<20:out['reason']='insufficient_independent_T_adjacency_variation';return out,[],None
 if len(frame)<np.linalg.matrix_rank(X)+len(np.unique(blocks))+20:
  out['reason']='fewer_than20_residual_organizer_degrees_of_freedom_after_nuisance_and_islet_effects';return out,[],None
 result=block_wild(X,(oz*bz)[kfit],y[kfit],blocks,np.maximum(deg[kfit],1))
 if result is None:out['reason']='interaction_not_identifiable_beyond_pairwise_baseline';return out,[],None
 out.update(result,H5=result['donor_effect'],status='estimated',reason='',n_organizers_model=int(kfit.sum()),baseline_columns=','.join(names))
 nodes=org.loc[kfit,['donor_id','xenium_serial','cell_id','native_islet_id','x_centroid_um','y_centroid_um']].copy()
 nodes['CCL_state']=score[kfit];nodes['B_LTB_exposure']=exposure[kfit];nodes['n_B_neighbors']=bcount[kfit];nodes['n_T_neighbors']=tcount[kfit]
 nodes['all_QC_neighbor_degree']=deg[kfit];nodes['T_neighbor_fraction']=y[kfit];nodes['interaction_predictor']=(oz*bz)[kfit]
 return out,[],nodes

def h5():
 gm=pd.read_parquet(T/'postIF_islet_objects.parquet');gm.donor_id=gm.donor_id.map(key)
 sm=json.loads(MAN.read_text())['sections'];rows=[];nodeparts=[]
 meta=pd.read_csv(T/'rna_donor_metadata_corrected.csv',dtype={'donor_id':str})
 for d in meta.donor_id:
  path=T/('rna_postIF_immune_stromal_states_donor_'+d+'.parquet')
  if not path.is_file():
   rows.extend({'donor_id':d,'disease_status':meta.set_index('donor_id').loc[d,'disease_status'],'configuration':name,'depth_threshold':dep,
      'neighbor_radius_um':rad,'identity_route':route,'H5':np.nan,'status':'untestable','reason':'no_released_exact_serial_IF_Islet_mask'}
      for name,dep,rad,route in CONFIGS);continue
  s=pd.read_parquet(path);m=pd.read_parquet(T/('rna_postIF_islet_membership_donor_'+d+'.parquet'))
  g=gm[gm.donor_id.eq(d)].copy();reports=[r for r in sm if r['tissue']=='pancreas'and key(r['donor_id'])==d]
  if len(reports)!=1:raise ValueError('Ambiguous source section')
  report=reports[0]
  if not s.xenium_serial.astype(str).str.zfill(7).eq(str(report['source_serial']).zfill(7)).all():raise ValueError('Wrong serial states')
  # Reconstruct the documented capture footprint from complete vendor nucleated RNA cells; no mask expansion or clipping.
  v=pd.read_parquet(Path(report['source_acquisition_path'])/'cells.parquet',columns=['x_centroid','y_centroid','nucleus_count','transcript_counts'])
  v=v[v.nucleus_count.gt(0)&v.transcript_counts.gt(0)]
  footprint,prov=PA.xenium_tissue_outline(v[['x_centroid','y_centroid']].to_numpy(float))
  for conf in CONFIGS:
   out,_,nodes=h5_one(d,g,s,m,report,footprint,conf);out['source_state_sha256']=sha(path);out['capture_policy']='Whole own-Islet polygon inside documented RNA-cell tissue footprint; no clipping';rows.append(out)
   if nodes is not None and conf[0]=='primary':nodeparts.append(nodes)
  print('H5 ORGANIZER',d,[(r['configuration'],r['status'],r['reason'])for r in rows if r['donor_id']==d],flush=True)
  del s,m,v
 allrows=pd.DataFrame(rows);allrows.to_csv(T/'spatial_H5_organizer_donor_endpoints_all_configurations.csv',index=False)
 allrows['donor_effect']=allrows.get('donor_effect',allrows.H5);allrows['n_islets']=allrows.get('n_islet_blocks',allrows.get('n_islets_organizer',np.nan))
 primary=allrows[allrows.configuration.eq('primary')].copy();primary.to_csv(T/'spatial_H5_organizer_donor_endpoints.csv',index=False)
 allrows.to_csv(T/'spatial_H5_organizer_null_QC.csv',index=False)
 if nodeparts:pd.concat(nodeparts,ignore_index=True).to_parquet(T/'spatial_H5_organizer_nodes.parquet',index=False)
 return allrows

def residualize_state(score,covariates,labels):
 X,names=numerical_columns(covariates)
 lab=pd.get_dummies(pd.Series(labels).astype(str),drop_first=True,dtype=float).to_numpy()
 X=np.column_stack([np.ones(len(score)),X,lab])
 resid=score-X@np.linalg.lstsq(X,score,rcond=None)[0]
 return zscore(resid),int(np.linalg.matrix_rank(X)),names
def bins(x,n=5):
 return np.minimum((rankdata(x,method='average')-1)*n/max(len(x),1),n-1).astype(int)
def pln_config(d,section,state,config):
 name,depth,radius,route=config
 row={'donor_id':d,'section_id':section,'configuration':name,'depth_threshold':depth,'neighbor_radius_um':radius,'identity_route':route,
  'architecture_excess':np.nan,'status':'untestable','reason':''}
 tt=state.native_lineage.isin(['T_CD4','T_CD8','T_Treg','T_CD4_CD8_both','T_subset_unresolved'])
 k=tt&state.raw_library_counts.ge(depth)&state.nucleus_count.gt(0)&state.n_genes_detected.ge(5)
 z=impute_area_QC(state[k]);n=len(z);row['n_QC_T_cells']=n
 row['identity_authority']='pln_20261003 final native released T lineage; T_NK_unresolved remains unresolved; no marker inclusion gate'
 row['n_missing_nuclear_area_QC']=int(z.nucleus_area_um2_missing_QC.sum());row['n_missing_cell_area_QC']=int(z.cell_area_um2_missing_QC.sum())
 row['fraction_cell_argmax_agrees_with_cluster']=float(z.native_cell_argmax_agrees_with_cluster.mean())if n and 'native_cell_argmax_agrees_with_cluster'in z else np.nan
 if n<200:row['reason']='fewer_than200_QC_supported_T_cells';return row,None
 memory=program(z,MEM);cyto=program(z,CYTO)
 if memory is None or cyto is None:row['reason']='memory_or_cytotoxic_gene_unmeasured';return row,None
 nm=(z[['RNA_'+g for g in MEM]].sum(axis=1)>0).sum();nc=(z[['RNA_'+g for g in CYTO]].sum(axis=1)>0).sum()
 row.update(n_memory_program_positive=int(nm),n_TOX_cytotoxic_program_positive=int(nc),
  mean_memory_state=float(memory.mean()),mean_TOX_cytotoxic_state=float(cyto.mean()),mean_library_depth=float(z.raw_library_counts.mean()))
 if min(nm,nc)<20:row['reason']='fewer_than20_raw_positive_T_cells_per_program';return row,None
 if 'genomic_control_counts'not in z:row['reason']='gDNA_control_unavailable';return row,None
 cov=pd.DataFrame({'log_depth':np.log1p(z.raw_library_counts.to_numpy()),'log_nucleus_area':np.log1p(z.nucleus_area_um2.to_numpy()),
  'log_cell_area':np.log1p(z.cell_area_um2.to_numpy()),'nucleus_area_missing_QC':z.nucleus_area_um2_missing_QC.to_numpy(),
  'cell_area_missing_QC':z.cell_area_um2_missing_QC.to_numpy(),'gDNA_per1000_library':z.genomic_control_counts.to_numpy()/z.raw_library_counts.to_numpy()*1000})
 for f in ['control_probe_counts','control_codeword_counts']:
  if f in z:cov[f+'_per1000_library']=z[f].to_numpy()/z.raw_library_counts.to_numpy()*1000
 labels=z.native_lineage.astype(str).to_numpy()
 contrast,rnk,controls=residualize_state(memory-cyto,cov,labels)
 if contrast.std()<1e-10:row['reason']='state_contrast_no_residual_variation';return row,None
 xy=z[['x_centroid_um','y_centroid_um']].to_numpy(float)
 distances,indices=cKDTree(xy).query(xy,k=9,distance_upper_bound=radius,workers=1)
 ii=np.repeat(np.arange(n),8);jj=indices[:,1:9].ravel();valid=(jj<n)&(distances[:,1:9].ravel()<=radius)&(jj!=ii)
 ii,jj=ii[valid],jj[valid]
 row.update(n_directed_T_neighbor_edges=len(ii),n_nodes_with_T_neighbors=int(np.unique(ii).size),state_nuisance_rank=rnk,state_control_fields=','.join(controls),
  n_native_lineage_groups=int(pd.Series(labels).nunique()))
 if len(ii)<400 or np.unique(ii).size<200:row['reason']='insufficient_geometric_T_neighbor_support';return row,None
 labels_code=pd.factorize(labels,sort=True)[0];strata=labels_code*10000+bins(z.raw_library_counts.to_numpy(),5)*1000+bins(z.nucleus_area_um2.to_numpy(),4)*100+z.nucleus_area_um2_missing_QC.to_numpy(int)*10+(z.genomic_control_counts.to_numpy()>0).astype(int)
 _,inverse,counts=np.unique(strata,return_inverse=True,return_counts=True);fraction=float(np.mean(counts[inverse]>=2))
 row['fraction_cells_exchangeable_in_QC_strata']=fraction
 if fraction<.70:row['reason']='fewer_than70percent_T_cells_exchangeable_within_QC_strata';return row,None
 obs=float(np.mean(contrast[ii]*contrast[jj]))
 order=np.argsort(strata,kind='stable');rng=np.random.default_rng(SEED);nperm=199
 null=np.empty(nperm);permuted=np.empty(n)
 for r in range(nperm):
  shuffled=np.lexsort((rng.random(n),strata));permuted[order]=contrast[shuffled]
  null[r]=np.mean(permuted[ii]*permuted[jj])
 row.update(status='estimated',reason='',observed_T_state_assortativity=obs,null_mean_assortativity=float(null.mean()),
  null_sd_assortativity=float(null.std(ddof=1)),null_mean_MC_SE=float(null.std(ddof=1)/np.sqrt(nperm)),
  architecture_excess=float(obs-null.mean()),within_section_p_positive=float((1+(null>=obs).sum())/(nperm+1)),n_spatial_label_null=nperm,
  conditional_null='TCF7/LEF1 minus TOX/GZMB/PRF1/NKG7 residual score labels shuffled within section/native-T-lineage/depth-quintile/nucleus-area-quartile/gDNA-positive strata; coordinates and measured score distribution fixed',
  organization_interpretation='Program-contrast spatial assortativity beyond measured QC/native-lineage strata; not clone, exhaustion, antigen specificity or causal migration')
 example=z[['cell_id','x_centroid_um','y_centroid_um','native_lineage','raw_library_counts','nucleus_area_um2']].copy()
 example['memory_state']=memory;example['TOX_cytotoxic_state']=cyto;example['residual_program_contrast']=contrast
 return row,example if name=='primary'else None

def pln_worker(report):
 rd=reader();hp=Path(report['h5ad_path']);cp=Path(report['native_calls_path']);d=key(report['donor_id']);section=str(report['section_id'])
 for path,pin in [(hp,report['h5ad_sha256']),(cp,report['native_calls_sha256'])]:
  if sha(path)!=pin:raise ValueError('Pinned pLN input changed')
 X,ids,genes,raw=rd.read_matrix(hp)
 if not raw.get('raw_counts_integer_verified'):raise ValueError('pLN raw counts not verified integer')
 c=pd.read_parquet(cp);c=c[c.section_key.astype(str).eq(section)].copy();c.cell_id=c.cell_id.astype(str)
 if c.cell_id.duplicated().any()or set(c.cell_id)!=set(ids):raise ValueError('pLN barcode bijection failed')
 c=c.set_index('cell_id').reindex(ids);depth=np.asarray(X.sum(axis=1)).ravel().astype(float)
 if not np.array_equal(depth,c.depth.to_numpy())or not c.nucleus_measured.astype(bool).all():raise ValueError('pLN raw depth/nucleus mismatch')
 gi={str(g):i for i,g in enumerate(genes)};required=MEM+CYTO
 absent=sorted(set(required)-set(gi))
 if absent:return {'rows':[{'donor_id':d,'section_id':section,'configuration':name,'depth_threshold':dep,'neighbor_radius_um':rad,'identity_route':route,
    'disease_status':report['disease_status'],'status':'untestable','reason':'actual_panel_missing:'+','.join(absent),'architecture_excess':np.nan}for name,dep,rad,route in CONFIGS],
    'source_hashes':{str(hp):sha(hp),str(cp):sha(cp)},'example':None}
 values={g:np.asarray(X[:,gi[g]].toarray()).ravel()for g in sorted(set(required+BGENES+TGENES)&set(gi))}
 bn=sum((values[g]>0).astype(np.int8)for g in BGENES if g in values);tn=sum((values[g]>0).astype(np.int8)for g in TGENES if g in values)
 features=np.asarray(X.getnnz(axis=1)).ravel()
 choose=c.nucleus_count.gt(0).to_numpy()&(depth>=10)&(features>=5)&c.lineage.isin(['T_CD4','T_CD8','T_Treg','T_CD4_CD8_both','T_subset_unresolved']).to_numpy()
 pos=np.flatnonzero(choose)
 state=pd.DataFrame({'cell_id':ids[pos],'raw_library_counts':depth[pos],'nucleus_count':c.nucleus_count.to_numpy()[pos],
  'native_lineage':c.lineage.astype(str).to_numpy()[pos],'n_genes_detected':features[pos],
  'n_B_identity_genes_raw_positive':bn[pos],'n_T_identity_genes_raw_positive':tn[pos]})
 if 'cell_agrees_with_cluster_call'in c:state['native_cell_argmax_agrees_with_cluster']=c.cell_agrees_with_cluster_call.fillna(False).to_numpy(bool)[pos]
 with h5py.File(hp,'r')as h:
  if 'spatial'not in h['obsm']:raise ValueError('pLN spatial coordinates absent')
  xy=np.asarray(h['obsm/spatial'][()])
  if 'spatial_coordinate_unit'not in h['uns']:raise ValueError('pLN spatial units not certified')
  unit=h['uns/spatial_coordinate_unit'][()];unit=unit.decode()if isinstance(unit,bytes)else str(unit)
  if 'microm'not in unit and unit!='um':raise ValueError('pLN coordinate unit invalid')
  if xy.shape!=(len(ids),2)or not np.isfinite(xy).all():raise ValueError('pLN coordinate malformed')
  state['x_centroid_um']=xy[pos,0];state['y_centroid_um']=xy[pos,1]
  for f in ['cell_area','nucleus_area','genomic_control_counts']:
   if f in h['obs']:state[f+'_um2'if f.endswith('area')else f]=rd.read_h5_array(h['obs'][f])[pos]
 vendor=Path(report['source_acquisition_path'])/'cells.parquet'
 vhash=None
 if vendor.is_file():
  v=pd.read_parquet(vendor);v.cell_id=v.cell_id.astype(str)
  if v.cell_id.duplicated().any()or not set(ids)<=set(v.cell_id):raise ValueError('pLN vendor control keys incomplete')
  v=v.set_index('cell_id').reindex(ids)
  if not np.array_equal(v.transcript_counts.to_numpy(),depth):raise ValueError('pLN vendor/library counts mismatch')
  for f in ['genomic_control_counts','control_probe_counts','control_codeword_counts']:
   if f in v:state[f]=v[f].to_numpy()[pos]
  vhash=sha(vendor)
 for g,val in values.items():state['RNA_'+g]=val[pos]
 rows=[];example=None
 for conf in CONFIGS:
  row,ex=pln_config(d,section,state,conf);row.update(disease_status=report['disease_status'],source_serial=str(report['source_serial']).zfill(7),
   matrix_native_key_bijection=True,raw_library_native_counts_exact=True,source_h5ad=str(hp),source_h5ad_sha256=report['h5ad_sha256'],
   source_vendor_controls=str(vendor)if vendor.is_file()else None,source_vendor_control_sha256=vhash)
  rows.append(row)
  if ex is not None:example=ex
 return {'rows':rows,'source_hashes':{str(hp):sha(hp),str(cp):sha(cp)},'example':example}

def pln(workers=2):
 cache=T/'spatial_H6_pln_section_cache';cache.mkdir(exist_ok=True)
 reports=[r for r in json.loads(MAN.read_text())['sections']if r['tissue']=='pLN'];results=[];jobs=[]
 for report in reports:
  name=str(report['donor_id'])+'__'+str(report['section_id']);p=cache/(name+'.json')
  if p.is_file():
   old=json.loads(p.read_text())
   if old.get('script_sha256')==sha(__file__)and old.get('source_manifest_sha256')==sha(MAN):results.append(old);continue
  jobs.append((report,p))
 with ProcessPoolExecutor(max_workers=workers)as pool:
  fs={pool.submit(pln_worker,r):(r,p)for r,p in jobs}
  for future in as_completed(fs):
   report,p=fs[future];res=future.result();example=res.pop('example')
   if example is not None:
    example['donor_id']=key(report['donor_id']);example['section_id']=str(report['section_id']);example['source_serial']=str(report['source_serial']).zfill(7)
    example.to_parquet(T/('spatial_H6_pln_primary_state_section_'+str(report['section_id'])+'.parquet'),index=False)
   res.update(script_sha256=sha(__file__),source_manifest_sha256=sha(MAN));jsonwrite(p,res);results.append(res)
   print('H6 PLN',report['donor_id'],report['section_id'],[(r['configuration'],r['status'])for r in res['rows']],flush=True)
 q=pd.DataFrame([row for result in results for row in result['rows']]);q.to_csv(T/'spatial_H6_pln_section_QC.csv',index=False)
 out=[]
 for (d,name),z in q.groupby(['donor_id','configuration'],sort=True):
  usable=z[z.status.eq('estimated')]
  base={'donor_id':d,'disease_status':z.disease_status.iloc[0],'configuration':name,'depth_threshold':int(z.depth_threshold.iloc[0]),
   'neighbor_radius_um':float(z.neighbor_radius_um.iloc[0]),'identity_route':z.identity_route.iloc[0],
   'n_sections_available':len(z),'n_sections_estimated':len(usable),'section_aggregation':'Unweighted mean of independent section-specific organization metrics; coordinates are never merged',
   'architecture_excess':np.nan,'status':'untestable','reason':'No section with sufficient supported T cells and exchangeable measured state programs'}
  if len(usable):
   weights=usable.n_QC_T_cells.to_numpy(float)
   base.update(status='estimated',reason='',architecture_excess=float(usable.architecture_excess.mean()),
    observed_T_state_assortativity=float(usable.observed_T_state_assortativity.mean()),
    null_mean_assortativity=float(usable.null_mean_assortativity.mean()),null_sd_assortativity=float(usable.null_sd_assortativity.mean()),
    mean_memory_state=float(np.average(usable.mean_memory_state,weights=weights)),
    mean_TOX_cytotoxic_state=float(np.average(usable.mean_TOX_cytotoxic_state,weights=weights)),
    mean_library_depth=float(np.average(usable.mean_library_depth,weights=weights)),
    n_QC_T_cells=int(usable.n_QC_T_cells.sum()),
    organization_scope='Depth/area/control/native-lineage-conditioned T program organization, not donor mean expression or donor disease comparison')
  out.append(base)
 donor=pd.DataFrame(out);donor.to_csv(T/'spatial_H6_pln_architecture.csv',index=False)
 return donor

def rank_association(z,age_sex=False):
 y=rankdata(z.H3.to_numpy(float));x=rankdata(z.architecture_excess.to_numpy(float))
 columns=[np.ones(len(z))]
 disease=pd.get_dummies(z.disease_status,drop_first=True,dtype=float).to_numpy()
 if disease.shape[1]:columns.extend(disease.T)
 for f in ['mean_memory_state','mean_TOX_cytotoxic_state']:
  columns.append(rankdata(z[f].to_numpy(float)))
 if age_sex:
  columns.append(rankdata(z.age_years.to_numpy(float)));columns.append(z.sex_analysis.eq('F').to_numpy(float))
 X=np.column_stack(columns);rank=int(np.linalg.matrix_rank(X))
 if len(z)<rank+6:return None,{'n_pair':len(z),'conditioning_design_rank':rank,'status':'untestable','reason':'paired_donors_less_than_conditioning_design_rank_plus6'}
 rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
 if min(rx.std(),ry.std())<1e-10:return None,{'n_pair':len(z),'conditioning_design_rank':rank,'status':'untestable','reason':'residual_endpoint_has_no_variation'}
 r=float(np.corrcoef(rx,ry)[0,1])
 return (r,rx,ry,X),{'n_pair':len(z),'conditioning_design_rank':rank,'status':'estimated','reason':'','effect':r,'H6':r}
def linked(h3source,column='donor_effect',age_sex=False):
 arch=pd.read_csv(T/'spatial_H6_pln_architecture.csv',dtype={'donor_id':str})
 if h3source is None:
  pd.DataFrame([{'hypothesis':'H6','status':'awaiting_supported_H3_source','reason':'pLN architecture evaluated independently; H3 linkage requires exact paired donor endpoint'}]).to_csv(T/'spatial_H6_linked_architecture_tests.csv',index=False);return pd.DataFrame()
 p=Path(h3source);h=pd.read_csv(p,dtype={'donor_id':str});h.donor_id=h.donor_id.map(key)
 if column not in h:raise ValueError('H3 column absent:'+column)
 if 'configuration'in h:h=h[h.configuration.eq('primary')]
 if 'variant'in h and h.donor_id.duplicated().any():h=h[h.variant.eq('primary')]
 if h.donor_id.duplicated().any():raise ValueError('H3 source not unique per donor; explicitly select primary source')
 h=h[['donor_id',column]].rename(columns={column:'H3'})
 meta=pd.read_csv(T/'rna_donor_metadata_corrected.csv',dtype={'donor_id':str})
 rows=[];primary_pairs=None
 for config in [name for name,_,_,_ in CONFIGS]:
  for omit5 in [False,True]:
   for extra in ([False,True]if config=='primary'and not omit5 else[False]):
    z=arch[arch.configuration.eq(config)&arch.status.eq('estimated')&arch.disease_status.isin(['ND','T1D'])].merge(h,on='donor_id',how='inner',validate='one_to_one')
    z=z.merge(meta[['donor_id','age_years','sex_analysis']],on='donor_id',validate='one_to_one')
    if omit5:z=z[~z.donor_id.isin(LEGACY)]
    z=z.replace([np.inf,-np.inf],np.nan).dropna(subset=['H3','architecture_excess','mean_memory_state','mean_TOX_cytotoxic_state']+(['age_years','sex_analysis']if extra else[])).copy()
    base={'hypothesis':'H6','configuration':config,'exclude_five_legacy_panRNA_donors':omit5,'age_sex_extra_sensitivity':extra,
     'source_H3':str(p),'source_H3_sha256':sha(p),'H3_column':column,'source_pLN_architecture_sha256':sha(T/'spatial_H6_pln_architecture.csv'),
     'predicted_direction':'positive','p_value':np.nan,'CI95_low':np.nan,'CI95_high':np.nan,
     'test_policy':'Two-sided disease-stratified Freedman-Lane permutation of ranked H3 residuals; paired independent ND/T1D donors; condition disease and both mean immune states',
     'interpretation_limit':'Cross-sectional conditional association; no clonality, antigen specificity, trafficking or causal direction'}
    result,q=rank_association(z,extra);base.update(q)
    if result is not None:
     r,rx,ry,X=result;group=z.disease_status.to_numpy();rng=np.random.default_rng(SEED);null=np.empty(99999)
     slots=[np.flatnonzero(group==g)for g in np.unique(group)]
     for b in range(len(null)):
      perm=ry.copy()
      for idx in slots:perm[idx]=ry[rng.permutation(idx)]
      perm-=X@np.linalg.lstsq(X,perm,rcond=None)[0]
      null[b]=np.corrcoef(rx,perm)[0,1]
     pv=float((1+(np.abs(null)>=abs(r)-1e-12).sum())/(len(null)+1))
     boots=[];rng=np.random.default_rng(SEED)
     for b in range(10000):
      ii=np.concatenate([rng.choice(idx,len(idx),replace=True)for idx in slots])
      rr,qq=rank_association(z.iloc[ii].copy(),extra)
      if rr is not None:boots.append(rr[0])
     lo,hi=np.quantile(boots,[.025,.975])if boots else(np.nan,np.nan)
     base.update(p_value=pv,CI95_low=float(lo),CI95_high=float(hi),n_donor_permutations=len(null),
      n_valid_stratified_donor_bootstrap=len(boots),null_association_mean=float(null.mean()),null_association_sd=float(null.std()),
      n_ND=int((group=='ND').sum()),n_T1D=int((group=='T1D').sum()))
     if config=='primary'and not omit5 and not extra:
      z['H6_residual_rank_pLN']=rx;z['H6_residual_rank_H3']=ry;primary_pairs=z.copy()
    rows.append(base)
 tests=pd.DataFrame(rows);tests.to_csv(T/'spatial_H6_linked_architecture_tests.csv',index=False)
 if primary_pairs is not None:primary_pairs.to_csv(T/'spatial_H6_linked_primary_donors.csv',index=False)
 return tests

def h5_sensitivity(h5rows):
 from scipy.stats import permutation_test
 rows=[]
 for conf in h5rows.configuration.unique():
  for omit5 in [False,True]:
   z=h5rows[h5rows.configuration.eq(conf)&h5rows.status.eq('estimated')]
   if omit5:z=z[~z.donor_id.isin(LEGACY)]
   nd=z.loc[z.disease_status.eq('ND'),'H5'].dropna().to_numpy();td=z.loc[z.disease_status.eq('T1D'),'H5'].dropna().to_numpy()
   base={'hypothesis':'H5','configuration':conf,'exclude_five_legacy_panRNA_donors':omit5,'n_ND':len(nd),'n_T1D':len(td),
    'effect':np.nan,'p_value':np.nan,'CI95_low':np.nan,'CI95_high':np.nan,'status':'untestable','reason':'fewer_than4_ND_or4_T1D_donors',
    'predicted_direction':'positive','inferential_unit':'independent donor, equal weight','root_family':'Fixed six; root applies BH6 including untestable p=1'}
   if min(len(nd),len(td))>=4:
    effect=float(td.mean()-nd.mean());p=permutation_test((td,nd),lambda a,b,axis:a.mean(axis)-b.mean(axis),vectorized=True,
      permutation_type='independent',alternative='two-sided',n_resamples=99999,batch=5000,rng=np.random.default_rng(SEED))
    rng=np.random.default_rng(SEED);boot=rng.choice(td,(10000,len(td)),replace=True).mean(1)-rng.choice(nd,(10000,len(nd)),replace=True).mean(1)
    lo,hi=np.quantile(boot,[.025,.975]);base.update(status='estimated',reason='',effect=effect,p_value=float(p.pvalue),CI95_low=float(lo),CI95_high=float(hi))
   rows.append(base)
 return pd.DataFrame(rows)

def figures():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 h=pd.read_csv(T/'spatial_H5_organizer_donor_endpoints.csv',dtype={'donor_id':str})if(T/'spatial_H5_organizer_donor_endpoints.csv').is_file()else pd.DataFrame()
 fig,ax=plt.subplots(1,2,figsize=(11,4.7))
 if len(h):
  text=h.reason.value_counts().head(5)
  ax[0].barh(range(len(text)),text.to_numpy(),color='#477b8b');ax[0].set_yticks(range(len(text)),[str(x).replace('_',' ')for x in text.index],fontsize=7);ax[0].set_xlabel('Donors');ax[0].set_title('Independent motif support / identifiability')
  for j,disease in enumerate(['ND','AAb+','T1D']):
   z=h[h.disease_status.eq(disease)&h.status.eq('estimated')]
   ax[1].scatter(np.full(len(z),j)+np.linspace(-.05,.05,len(z)),z.H5,color={'ND':'#267687','AAb+':'#be8a30','T1D':'#af4d55'}[disease],s=35)
   for x,r in zip(np.full(len(z),j)+np.linspace(-.05,.05,len(z)),z.itertuples()):ax[1].annotate(r.donor_id,(x,r.H5),fontsize=7)
  ax[1].axhline(0,c='gray',ls=':');ax[1].set_xticks([0,1,2],['ND','AAb+','T1D']);ax[1].set_ylabel('Conditional organizer × B-LTB interaction');ax[1].set_title('Donor coefficients; unavailable stays unavailable')
 else:ax[0].text(.5,.5,'H5 foundation unavailable',ha='center');ax[1].axis('off')
 fig.suptitle('H5: third-order islet organizer motif beyond pairwise structure',fontsize=12)
 fig.tight_layout(rect=(0,.08,1,.94));fig.text(.01,.02,'Own true IF-Islet organizer masks; physical cell neighbors may cross boundaries. Released native B/T identities; marker QC only. Islet-block wild null. No expansion rings or cell-transfer identity.',fontsize=8)
 for ext in ['png','svg']:fig.savefig(F/('spatial_H5_organizer_motifs.'+ext),dpi=200)
 plt.close(fig)
 arch=pd.read_csv(T/'spatial_H6_pln_architecture.csv',dtype={'donor_id':str});a=arch[arch.configuration.eq('primary')]
 fig,ax=plt.subplots(1,2,figsize=(11,4.7));z=a[a.status.eq('estimated')]
 for disease,col in [('ND','#267687'),('T1D','#af4d55'),('AAb+','#be8a30')]:
  zz=z[z.disease_status.eq(disease)]
  ax[0].scatter(zz.null_mean_assortativity,zz.observed_T_state_assortativity,c=col,s=35,label=disease+' n='+str(len(zz)))
  for r in zz.itertuples():ax[0].annotate(r.donor_id,(r.null_mean_assortativity,r.observed_T_state_assortativity),fontsize=6,xytext=(2,3),textcoords='offset points')
 ax[0].set_xlabel('QC / native-lineage conditional-null assortativity');ax[0].set_ylabel('Observed T-program assortativity');ax[0].set_title('pLN architecture beyond state abundance');ax[0].legend(fontsize=7)
 pairs=T/'spatial_H6_linked_primary_donors.csv'
 if pairs.is_file():
  p=pd.read_csv(pairs,dtype={'donor_id':str})
  for disease,col in [('ND','#267687'),('T1D','#af4d55')]:
   zz=p[p.disease_status.eq(disease)];ax[1].scatter(zz.H6_residual_rank_pLN,zz.H6_residual_rank_H3,c=col,s=35)
   for r in zz.itertuples():ax[1].annotate(r.donor_id,(r.H6_residual_rank_pLN,r.H6_residual_rank_H3),fontsize=6,xytext=(2,3),textcoords='offset points')
  ax[1].set_xlabel('pLN residual architecture rank');ax[1].set_ylabel('Pancreatic H3 residual domain-excess rank');ax[1].set_title('Paired ND/T1D donors; disease and mean states adjusted')
 else:ax[1].text(.5,.5,'H3 source pending or paired support insufficient',ha='center',va='center',wrap=True);ax[1].set_axis_off()
 fig.suptitle('H6: linked pLN T-state and pancreatic islet domain organization',fontsize=12)
 fig.tight_layout(rect=(0,.08,1,.94));fig.text(.01,.02,'TCF7/LEF1 versus TOX/GZMB/PRF1/NKG7 residual-score architecture; native partition QC10counts/5genes, positive measured nucleus. Sections remain independent; donor is inferential unit.',fontsize=8)
 for ext in ['png','svg']:fig.savefig(F/('spatial_H6_linked_architecture.'+ext),dpi=200)
 plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);ap.add_argument('--pln-only',action='store_true');ap.add_argument('--link-only',action='store_true')
 ap.add_argument('--H3-source');ap.add_argument('--H3-column',default='donor_effect');args=ap.parse_args()
 T.mkdir(exist_ok=True);F.mkdir(exist_ok=True)
 plan={'id':'spatial_H5_H6_RNA_v3_released_identity','script_sha256':sha(__file__),'created_utc':datetime.now(timezone.utc).isoformat(),
  'user_correction':'User rejected ad hoc two-marker identity gates. Prior gated scarcity and H6 metrics withdrawn; native released identity with original partition QC replaces them. Higher-order null and physical estimands unchanged.',
  'source_shared_plan':str(R/'provenance/spatial_hypothesis_plan_v2.json'),'source_shared_plan_sha256':sha(R/'provenance/spatial_hypothesis_plan_v2.json'),
  'primary_depth':10,'sensitivities':CONFIGS,'primary_evidence':'Current released native identities only; partition QC >=10 raw counts,>=5 detected genes,measured nucleus_count>0. Raw markers and producer decidedness audit only, no new identity filter',
  'H5_response':'Released native T physical-neighbor fraction among all exact-serial nucleated QC RNA centroids within50um of organizer; organizer remains inside actual own Islet mask; no islet expansion ring',
  'H5_predictor':'Standardized organizer CCL19/CCL21 state x standardized physical-neighbor released B-cell LTB exposure',
  'H5_baseline':'Pairwise main states, B-neighbor count, local all-QC degree50um and tissue-footprint-corrected density100um, depth, nuclear/cell area, genomic/negative control burden, organizer vascular/fibro identity, continuous signed-islet distance; weighted islet fixed effects',
  'H5_estimand':'Weighted within-islet partial interaction coefficient; restricted pairwise-only wild residual null by whole islet block',
  'H5_minimum':'100 model organizers,20CCL-positive organizers,20unique B and20unique T physical neighbors,10LTB-positive B neighbors,8islet blocks,20residual organizer degrees of freedom beyond effective nuisance/islet rank and identifiable interaction',
  'H6_program':'TCF7/LEF1 versusTOX/GZMB/PRF1/NKG7 measured contrast, not exhaustion or clonality',
  'H6_graph':'Directed8-nearest-released-T neighbors capped50um; radius25/75 sensitivity',
  'H6_null':'199 score-label shuffles within section/native-T-lineage/depth quintile/nuclear-area quartile/gDNA-positive strata after QC/nuisance residualization; fixed positions and score distribution',
  'H6_minimum':'200 released T/section,20 raw-positive cells per program,400directed edges/200connected nodes,70% score-exchangeable cells; pairedND/T1D donors>=rank+6',
  'H6_link':'Partial rank association of paired donor architecture_excess with H3 donor_effect, conditioned disease and both mean states; two-sided donor Freedman-Lane within-disease residual permutation99999 and stratified donor bootstrap10000',
  'primary_family':'Root combines fixed six viaBH with untestablep1; historical mean-expression H5/H6 remain descriptive only',
  'limitations':['Native RNA stromal/vascular labels remain uncertain; no Pheno-stroma alias','Raw marker support is diagnostic and never an inclusion criterion','IF masks are exact-serial; additional locally available acquisitions are being called and audited',
    'Pairwise-baseline wild residual symmetry/islet independence and conditional residual exchangeability are approximations','Cross-tissue association is observational and cannot establish migration, antigen specificity or causal direction']}
 jsonwrite(T/'spatial_H5_H6_analysis_plan.json',plan)
 if not args.link_only:
  pln(args.workers)
  if not args.pln_only:
   if not(T/'qc_postIF_islet_manifest.json').is_file():raise ValueError('PostIF full state export incomplete; wait completion')
   h5rows=h5()
  elif(T/'spatial_H5_organizer_donor_endpoints_all_configurations.csv').is_file():h5rows=pd.read_csv(T/'spatial_H5_organizer_donor_endpoints_all_configurations.csv',dtype={'donor_id':str})
  else:h5rows=None
 else:h5rows=pd.read_csv(T/'spatial_H5_organizer_donor_endpoints_all_configurations.csv',dtype={'donor_id':str})if(T/'spatial_H5_organizer_donor_endpoints_all_configurations.csv').is_file()else None
 tests=linked(args.H3_source,args.H3_column)
 sensitivity=[tests]
 if h5rows is not None:sensitivity.append(h5_sensitivity(h5rows))
 pd.concat([x for x in sensitivity if len(x)],ignore_index=True).to_csv(T/'spatial_H5_H6_sensitivities.csv',index=False)if any(len(x)for x in sensitivity)else None
 figures()
 jsonwrite(T/'spatial_H5_H6_execution_manifest.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),
  'source_RNA_manifest_sha256':sha(MAN),'clinical_authority_source':str(T/'rna_donor_metadata_corrected.csv'),'primary_depth':10,'inferential_unit':'independent donor',
  'primary_family_owner':'root fixed six BH','source_native_files_modified':False,'H3_source':args.H3_source})
 print('SPATIAL RNA H5/H6 COMPLETE',flush=True)
if __name__=='__main__':main()
