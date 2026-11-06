"""Post-image-review ambiguity diagnostic; never changes primary assignment.

For each assigned pair, examine other same-target-category objects within 200 um
of the transformed RNA center, including those outside the primary area gate.
Flag a closer alternative with neighborhood score within .05 of the assigned
score, or any alternative ahead by >.05. Report zero/.05/.10 sensitivity. These
thresholds are descriptive review triggers, not calibrated probabilities.
"""
from pathlib import Path
import json,concurrent.futures
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from audit_sequential_endocrine import anchors,compare_neighbor_vectors
O=None;S=[]
def one(d):
 donor=d['donor_id'];p=O/'results'/donor;a=pd.read_parquet(p/'RNA_registered_objects.parquet');ps=next(s for s in S if s['role']=='protein' and s['donor_id']==donor);b=pd.read_parquet(Path(ps['output'])/'objects.parquet')
 pairs=pd.read_parquet(p/'pairs.parquet');av=anchors(a,1500);bv=anchors(b,1500);dist=cdist(pairs[['RNA_x_um','RNA_y_um']],b[['x_um','y_um']]);rows=[];flags=[]
 for i,r in enumerate(pairs.itertuples()):
  candidates=b[(dist[i]<=200)&b.object_class.eq(r.protein_object_class)&b.structure_id.ne(r.fixed_id)];near=[]
  for q in candidates.itertuples():
   score=compare_neighbor_vectors(av[r.moving_id],bv[q.structure_id])['score']
   if not np.isfinite(score) or not np.isfinite(r.neighbor_score):continue
   distance=float(np.hypot(q.x_um-r.RNA_x_um,q.y_um-r.RNA_y_um));ratio=max(q.area_um2/r.RNA_area_um2,r.RNA_area_um2/q.area_um2)
   rr=dict(donor_id=donor,moving_id=r.moving_id,fixed_id=r.fixed_id,alternative_id=q.structure_id,alternative_score=score,assigned_score=r.neighbor_score,alternative_distance_um=distance,assigned_distance_um=r.distance_um,alternative_area_ratio=ratio,outside_primary_area_gate=bool(ratio>2.5));near.append(rr);rows.append(rr)
  f=dict(moving_id=r.moving_id,n_nearby_same_category_alternatives=len(near))
  for tolerance in [0,.05,.1]:
   f[f'anatomy_alternative_flag_tolerance_{tolerance}']=any((q['alternative_distance_um']<r.distance_um and q['alternative_score']>=r.neighbor_score-tolerance) or q['alternative_score']>r.neighbor_score+tolerance for q in near)
  flags.append(f)
 flag=pd.DataFrame(flags);newcols=[c for c in flag if c!='moving_id'];pairs=pairs.drop(columns=[c for c in newcols if c in pairs]).merge(flag,on='moving_id',validate='one_to_one')
 pairs['anatomy_alternative_review_flag']=pairs['anatomy_alternative_flag_tolerance_0.05']
 pairs['regional_candidate_before_alternative_review']=pairs.regional_candidate_for_review
 pairs['regional_candidate_for_review']=pairs.regional_candidate_before_alternative_review & ~pairs.anatomy_alternative_review_flag
 pairs.to_parquet(p/'pairs.parquet',index=False);pairs.to_csv(p/'pairs.csv',index=False);pd.DataFrame(rows).to_parquet(p/'local_anatomical_alternatives.parquet',index=False)
 import hashlib
 d['pair_sha256']=hashlib.sha256((p/'pairs.parquet').read_bytes()).hexdigest();d['n_review_candidates_before_secondary_alternatives']=int(pairs.regional_candidate_before_alternative_review.sum());d['n_anatomy_alternative_flags']=int(pairs.anatomy_alternative_review_flag.sum());d['n_regionally_supported_unambiguous_candidates']=int(pairs.regional_candidate_for_review.sum());d['secondary_review_policy']=__doc__
 (p/'decision.json').write_text(json.dumps(d,indent=2)+'\n');print('ALTERNATIVE_REVIEW',donor,d['n_review_candidates_before_secondary_alternatives'],'->',d['n_regionally_supported_unambiguous_candidates'],flush=True);return d
def main(argv=None):
 global O,S
 import argparse
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',required=True,type=Path)
 args=ap.parse_args(argv);O=args.run;S=json.loads((O/'run_specs.json').read_text())
 ds=json.loads((O/'release_decisions.json').read_text())
 if any('secondary_review_policy' in d for d in ds):raise RuntimeError('Secondary review already applied; preserve the receipt')
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as ex:r=list(ex.map(one,ds))
 (O/'release_decisions.json').write_text(json.dumps(r,indent=2)+'\n')
 pd.concat([pd.read_parquet(O/'results'/d['donor_id']/'pairs.parquet') for d in r]).to_parquet(O/'all_pairs.parquet',index=False)
 (O/'secondary_alternative_review.json').write_text(json.dumps(dict(policy=__doc__,timing='Added after actual-image review identified competing local structures; not part of the prespecified primary assignment/rubric',results=[{k:d[k] for k in ['donor_id','n_anatomy_alternative_flags','n_review_candidates_before_secondary_alternatives','n_regionally_supported_unambiguous_candidates']} for d in r]),indent=2)+'\n')

if __name__=='__main__':main()
