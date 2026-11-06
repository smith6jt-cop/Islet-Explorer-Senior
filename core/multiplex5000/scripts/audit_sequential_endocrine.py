"""Prespecified all-category correspondence; geometry leads, overlap never gates."""
from pathlib import Path
import json,sys,concurrent.futures
import numpy as np,pandas as pd,shapely
from shapely.affinity import affine_transform
from scipy.spatial.distance import cdist
from multiplex5000.integration.serial_context import geometry_descriptor,assign_geometry,compare_neighbor_vectors
from multiplex5000.integration.endocrine_objects import verify_source,digest
O=None;S=[];PLAN={}

def load(s):
 p=Path(s['output']);q=json.loads((p/'verified.json').read_text());verify_source(p/'objects.parquet',q['output_sha256']);assert (p/'image_qc.json').exists()
 for key in ['n_multiply_owned_cells','n_cut_cells','n_construction_count_disagreements','n_single_constructor_geometry_disagreements']:
  if q[key]:raise ValueError(f'{s["role"]} {s["donor_id"]} fails {key}: {q[key]}')
 d=pd.read_parquet(p/'objects.parquet');d['geometry']=shapely.from_wkb(d.geometry_wkb)
 iq=pd.read_parquet(p/'image_object_qc.parquet');d=d.merge(iq[['structure_id','within_image','DAPI_positive_fraction']],on='structure_id',validate='one_to_one')
 return d

def anchors(objects,radius):
 a=objects[objects.object_class.eq('Islet')];axy=a[['x_um','y_um']].to_numpy();ids=a.structure_id.to_numpy();out={}
 for r in objects.itertuples():
  v=axy-np.array([r.x_um,r.y_um]);dist=np.linalg.norm(v,axis=1);i=np.flatnonzero((dist<=radius)&(ids!=r.structure_id));i=i[np.lexsort((ids[i],dist[i]))][:6];out[r.structure_id]=v[i]
 return out

def one(donor):
 ps=next(s for s in S if s['role']=='protein' and s['donor_id']==donor);rs=next(s for s in S if s['role']=='RNA' and s['donor_id']==donor)
 p=O/'results'/donor;p.mkdir(parents=True,exist_ok=True)
 if (p/'decision.json').exists():return json.loads((p/'decision.json').read_text())
 a=load(rs);excluded=a[~a.proposal_evidence_available].copy();a=a[a.proposal_evidence_available].copy();b=load(ps)
 verify_source(ps['transform_path'],ps['transform_sha256']);tf=json.loads(Path(ps['transform_path']).read_text());assert not tf.get('has_displacement')
 m=np.array(tf['affine']);params=[m[0,0],m[0,1],m[1,0],m[1,1],m[0,2],m[1,2]]
 a['geometry']=[affine_transform(g,params) for g in a.geometry];a['geometry_wkb']=[g.wkb for g in a.geometry]
 desc=pd.DataFrame([geometry_descriptor(g) for g in a.geometry]);
 for k in desc:a[k]=desc[k].to_numpy()
 a=a.sort_values('structure_id').reset_index(drop=True);b=b.sort_values('structure_id').reset_index(drop=True)
 a.drop(columns='geometry').to_parquet(p/'RNA_registered_objects.parquet',index=False)
 assignments={w:assign_geometry(a,b,position_weight=w,max_dist_um=200,area_ratio=2.5) for w in [.25,.5,.75]}
 pair=assignments[.25];pair['donor_id']=donor
 for w in [.5,.75]:
  other=dict(zip(assignments[w].moving_id,assignments[w].fixed_id));pair[f'same_target_weight_{w}']=[other.get(i)==j for i,j in zip(pair.moving_id,pair.fixed_id)]
 ai=a.set_index('structure_id');bi=b.set_index('structure_id')
 for prefix,tab,key in [('RNA',ai,'moving_id'),('protein',bi,'fixed_id')]:
  for col in ['object_class','area_um2','x_um','y_um','DAPI_positive_fraction','within_image']:
   pair[f'{prefix}_{col}']=pair[key].map(tab[col])
 pair['category_transition']=pair.RNA_object_class.ne(pair.protein_object_class)
 pair['same_cell_identity_claim']=False
 stability=pd.read_parquet(Path(rs['output'])/'seed_group_stability.parquet');members=pd.read_parquet(Path(rs['output'])/'object_members.parquet')
 st=members.merge(stability,on='cell_id',how='inner',validate='many_to_one').groupby('structure_id').stable_both.mean()
 pair['RNA_seed_fraction_stable_at_40_and_60um']=pair.moving_id.map(st)
 pair['RNA_proposal_sensitivity_flag']=pair.RNA_seed_fraction_stable_at_40_and_60um.lt(1)
 pair['intersection_area_um2']=[ai.loc[i,'geometry'].intersection(bi.loc[j,'geometry']).area for i,j in zip(pair.moving_id,pair.fixed_id)]
 pair['zero_projected_overlap']=pair.intersection_area_um2.lt(1e-8)
 pair['within_reviewed_images']=pair.RNA_within_image & pair.protein_within_image
 distance=cdist(b[['x_um','y_um']],b[['x_um','y_um']]);bindex={i:k for k,i in enumerate(b.structure_id)}
 evidence=[];controlrows=[];radiusstats=[];blocks=[];classstats=[]
 mx=float(np.median(pair.RNA_x_um));my=float(np.median(pair.RNA_y_um))
 for radius in [750,1500,2250]:
  av=anchors(a,radius);bv=anchors(b,radius); pools=[];actuals=[];rrows=[];poolclasses=[]
  for r in pair.itertuples():
   actual=compare_neighbor_vectors(av[r.moving_id],bv[r.fixed_id]);j=bindex[r.fixed_id]
   ratio=b.area_um2.to_numpy()/b.area_um2.iloc[j]
   ok=(distance[j]>=300)&(distance[j]<=2000)&(ratio>=1/2.5)&(ratio<=2.5)&b.object_class.eq(r.protein_object_class).to_numpy()
   candidates=np.flatnonzero(ok);candidates=candidates[np.lexsort((b.structure_id.to_numpy()[candidates],distance[j,candidates]))]
   cs=[];targetids=[]
   for k in candidates:
    target=b.structure_id.iloc[k]
    if abs(len(bv[target])-len(bv[r.fixed_id]))>2:continue
    control=compare_neighbor_vectors(av[r.moving_id],bv[target])
    if control['status']=='measured':cs.append(control['score']);targetids.append(target)
    if len(cs)>=20:break
   block=int(r.RNA_x_um>=mx)+2*int(r.RNA_y_um>=my)
   rr=dict(donor_id=donor,moving_id=r.moving_id,fixed_id=r.fixed_id,radius_um=radius,spatial_block=block,
       RNA_object_class=r.RNA_object_class,protein_object_class=r.protein_object_class,neighbor_score=actual['score'],neighbor_status=actual['status'],
       n_source_neighbors=actual['n_source_neighbors'],n_target_neighbors=actual['n_target_neighbors'],n_agreeing=actual['n_agreeing'],n_controls=len(cs),
       control_median=float(np.median(cs)) if cs else np.nan,control_p95=float(np.percentile(cs,95)) if cs else np.nan)
   rr['context_delta']=rr['neighbor_score']-rr['control_median'];rrows.append(rr)
   controlrows.extend(dict(donor_id=donor,moving_id=r.moving_id,fixed_id=r.fixed_id,radius_um=radius,control_target_id=t,score=v) for t,v in zip(targetids,cs))
   if actual['status']=='measured' and cs:
    actuals.append(actual['score']);pools.append(cs);poolclasses.append(r.RNA_object_class)
  ef=pd.DataFrame(rrows);evidence.extend(rrows)
  rng=np.random.default_rng(20261106);null=[float(np.median([rng.choice(v) for v in pools])) for _ in range(199)] if pools else []
  actualmedian=float(np.median(actuals)) if actuals else np.nan;p95=float(np.percentile(null,95)) if null else np.nan
  radiusstats.append(dict(donor_id=donor,radius_um=radius,n_contextual_queries=len(actuals),actual_median=actualmedian,control_median_p95=p95,exceeds_control=bool(len(actuals)>=10 and actualmedian>p95)))
  for cls in ['Islet','EndocrineCluster','EndocrineSingle']:
   ix=[i for i,c in enumerate(poolclasses) if c==cls];cp=[pools[i] for i in ix];ca=[actuals[i] for i in ix]
   crng=np.random.default_rng(20261106);cn=[float(np.median([crng.choice(v) for v in cp])) for _ in range(199)] if cp else []
   cm=float(np.median(ca)) if ca else np.nan;c95=float(np.percentile(cn,95)) if cn else np.nan
   cb=[]
   for _,cg in ef[ef.RNA_object_class.eq(cls)].dropna(subset=['context_delta']).groupby('spatial_block'):
    cb.append(bool(len(cg)>=3 and cg.context_delta.median()>0))
   classstats.append(dict(donor_id=donor,RNA_object_class=cls,radius_um=radius,n_contextual_queries=len(ca),actual_median=cm,control_median_p95=c95,
     exceeds_control=bool(len(ca)>=10 and cm>c95),n_positive_blocks=sum(cb)))
  for block,g in ef.groupby('spatial_block'):
   valid=g.dropna(subset=['context_delta']);blocks.append(dict(donor_id=donor,radius_um=radius,spatial_block=int(block),n_contextual_queries=len(valid),median_delta=float(valid.context_delta.median()),positive=bool(len(valid)>=3 and valid.context_delta.median()>0)))
 ev=pd.DataFrame(evidence);ev.to_parquet(p/'anatomy_evidence.parquet',index=False);pd.DataFrame(controlrows).to_parquet(p/'conditional_controls.parquet',index=False)
 primary=ev[ev.radius_um.eq(1500)].drop(columns=['donor_id','fixed_id','RNA_object_class','protein_object_class']);pair=pair.merge(primary,on='moving_id',validate='one_to_one')
 pair['geometry_weight_stable']=pair['same_target_weight_0.5']&pair['same_target_weight_0.75']
 pair['local_context_supported']=pair.neighbor_status.eq('measured')&pair.n_controls.gt(0)&pair.neighbor_score.gt(pair.control_p95)&pair.n_agreeing.ge(3)
 pair['ambiguous_geometry']=~pair.mutual_composite_minimum|pair.alternative_cost_margin.le(.05)|~pair.geometry_weight_stable
 pair['image_sampling_review_flag']=pair.RNA_DAPI_positive_fraction.eq(0)|pair.protein_DAPI_positive_fraction.eq(0)|pair.RNA_DAPI_positive_fraction.isna()|pair.protein_DAPI_positive_fraction.isna()
 pair['regional_candidate_for_review']=pair.local_context_supported&pair.within_reviewed_images&~pair.ambiguous_geometry&~pair.image_sampling_review_flag
 pair.to_parquet(p/'pairs.parquet',index=False);pair.drop(columns=[]).to_csv(p/'pairs.csv',index=False)
 for w,t in assignments.items():t.to_parquet(p/f'assignment_position_weight_{w}.parquet',index=False)
 qc=json.loads(Path(ps['transform_path']).with_name('qc.json').read_text());ncc=qc['footprint_ncc_z']
 npositive=sum(q['positive'] for q in blocks if q['radius_um']==1500);nr=sum(r['exceeds_control'] for r in radiusstats)
 support=ncc>=5 and nr>=2 and npositive>=3
 accepted_classes=[cls for cls in ['Islet','EndocrineCluster','EndocrineSingle'] if ncc>=5 and sum(c['exceeds_control'] for c in classstats if c['RNA_object_class']==cls)>=2 and next(c['n_positive_blocks'] for c in classstats if c['RNA_object_class']==cls and c['radius_um']==1500)>=3]
 decision=dict(donor_id=donor,historical_verdict=ps['historical_verdict'],n_RNA_objects=len(a),n_protein_objects=len(b),n_excluded_without_RNA_seed=len(excluded),n_pairs=len(pair),
  n_RNA_unmatched=len(a)-len(pair),n_protein_unmatched=len(b)-len(pair),n_zero_overlap=int(pair.zero_projected_overlap.sum()),n_category_transitions=int(pair.category_transition.sum()),
  n_geometry_weight_stable=int(pair.geometry_weight_stable.sum()),n_local_context_supported=int(pair.local_context_supported.sum()),n_regionally_supported_unambiguous_candidates=int(pair.regional_candidate_for_review.sum()),
  historical_NCC_z=ncc,n_supporting_radii=nr,n_positive_blocks=npositive,regional_rubric_met=bool(support),
  release_scope='exploratory_regional_candidates_with_review_flags' if support and accepted_classes else 'diagnostic_only_not_accepted_for_correspondence_release',
  pooled_regional_rubric_met=bool(support),release_accepted_for_regional_review=bool(support and accepted_classes),accepted_RNA_categories_for_regional_review=accepted_classes,
  production_accepted=False,individual_pairs_validated=False,radius_statistics=radiusstats,spatial_blocks=blocks,class_statistics=classstats,
  transform_sha256=ps['transform_sha256'],pair_sha256=digest(p/'pairs.parquet'))
 (p/'decision.json').write_text(json.dumps(decision,indent=2)+'\n');print('REASSESSED',donor,len(pair),'regional',support,'radii',nr,'blocks',npositive,flush=True)
 return decision

def main(argv=None):
 global O,S,PLAN
 import argparse
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',required=True,type=Path)
 args=ap.parse_args(argv);O=args.run;S=json.loads((O/'run_specs.json').read_text());PLAN=json.loads((O/'analysis_plan.json').read_text())
 assert json.loads((O/'visual_review.json').read_text())['complete'],'Actual-image visual review must precede correspondence'
 with concurrent.futures.ProcessPoolExecutor(max_workers=5) as ex:r=list(ex.map(one,[s['donor_id'] for s in S if s['role']=='protein']))
 (O/'release_decisions.json').write_text(json.dumps(r,indent=2)+'\n')
 pd.DataFrame([{k:v for k,v in x.items() if k not in ['radius_statistics','spatial_blocks','class_statistics']} for x in r]).to_csv(O/'donor_reassessment.csv',index=False)
 pair=pd.concat([pd.read_parquet(O/'results'/x['donor_id']/'pairs.parquet') for x in r]);pair.to_parquet(O/'all_pairs.parquet',index=False)
 pair.groupby(['donor_id','RNA_object_class','protein_object_class']).size().rename('n_candidates').reset_index().to_csv(O/'class_transitions.csv',index=False)

if __name__=='__main__':main()
