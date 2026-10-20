#!/usr/bin/env python3
"""One fixed six-test family; independent donor disease effects and conditional H6."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path('/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/islet_immune_reassessment_20261006')
T=R/'tables';F=R/'figures';SEED=20261006;FAMILY='spatial_v2_20261006'
COL={'ND':'#3a6e94','AAb':'#bd8a29','T1D':'#a74151'}
def read(name):return pd.read_csv(T/name,dtype={'donor_id':str})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def contrast(z):
    a=z.loc[z.disease_status.eq('T1D'),'donor_effect'].dropna().to_numpy(float)
    b=z.loc[z.disease_status.eq('ND'),'donor_effect'].dropna().to_numpy(float)
    out=dict(effect=np.nan,ci_low=np.nan,ci_high=np.nan,p=1.,n_test=len(a),n_ref=len(b),test_status='unavailable',reason='requires at least4 supported T1D and4ND donors',method='independent donor two-sided label permutation')
    if len(a)<4 or len(b)<4:return out
    effect=float(a.mean()-b.mean());rng=np.random.default_rng(SEED)
    perm=stats.permutation_test((a,b),lambda x,y,axis:np.mean(x,axis=axis)-np.mean(y,axis=axis),vectorized=True,batch=2000,n_resamples=99999,alternative='two-sided',rng=rng)
    boot=rng.choice(a,(10000,len(a)),replace=True).mean(axis=1)-rng.choice(b,(10000,len(b)),replace=True).mean(axis=1)
    lo,hi=np.quantile(boot,[.025,.975]);out.update(effect=effect,ci_low=float(lo),ci_high=float(hi),p=float(perm.pvalue),test_status='estimated',reason='')
    return out
def bh6(p):
    p=np.asarray(p,float);order=np.argsort(p);q=np.empty(6);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*6/np.arange(1,7))[::-1])[::-1]);return q
def primary_donors():
    h12=read('spatial_H1_H2_donor_endpoints.csv');h34=read('spatial_H3_H4_donor_endpoints.csv');h5=read('spatial_H5_organizer_donor_endpoints.csv');h6=read('spatial_H6_pln_architecture.csv')
    rows=[]
    for _,r in h12.iterrows():
        for h in ['H1','H2']:
            n=int(r[h+'_n_supported_islets']);value=r[h+'_endpoint'] if n>=5 else np.nan
            row=dict(hypothesis=h,donor_id=r.donor_id,disease_status=r.disease_status,donor_effect=value,n_islets=n,status='measured' if np.isfinite(value) else 'unavailable',reason=r[h+'_availability_reason'] if n>=5 else 'fewer than5 supported islets',endpoint_descriptive=r.get(h+'_endpoint_descriptive',r[h+'_endpoint']))
            for c in (['vascular_excess','neural_excess'] if h=='H1' else ['fragmentation_z','interface_z']):row[c]=r[h+'_'+c] if n>=5 else np.nan
            rows.append(row)
    for _,r in h34.iterrows():rows.append(r.to_dict())
    for _,r in h5.iterrows():
        rows.append(dict(hypothesis='H5',donor_id=r.donor_id,disease_status=r.disease_status,donor_effect=r.get('donor_effect',r.get('H5',np.nan)),n_islets=r.get('n_islet_blocks',r.get('n_islets_organizer',np.nan)),status=r.status,reason=r.get('reason','')))
    for _,r in h6.loc[h6.configuration.eq('primary')].iterrows():
        rows.append(dict(hypothesis='H6',donor_id=r.donor_id,disease_status=r.disease_status,donor_effect=r.architecture_excess,status=r.status,reason=r.get('reason',''),endpoint_scope='pLN spatial architecture only; H6 inference is paired pancreatic H3 association',mean_memory_state=r.mean_memory_state,mean_TOX_cytotoxic_state=r.mean_TOX_cytotoxic_state))
    z=pd.DataFrame(rows);z['disease_status']=z.disease_status.replace({'AAb+':'AAb'});z['family']=FAMILY
    assert not z.duplicated(['hypothesis','donor_id']).any(),'Duplicate independent donor endpoints'
    z.to_csv(T/'spatial_hypothesis_donor_endpoints.csv',index=False)
    return z
def sensitivity_donors(z):
    out=[]
    h12=read('spatial_H1_H2_sensitivities.csv')
    for c,s in h12.groupby('configuration'):
        if c=='native_primary':continue
        for h in ['H1','H2']:
            q=s.copy();q['donor_effect']=q[h+'_endpoint'].where(q[h+'_n_supported_islets'].ge(5));q.disease_status=q.disease_status.replace({'AAb+':'AAb'})
            r=contrast(q);r.update(hypothesis=h,sensitivity=c);out.append(r)
    h34=read('spatial_H3_H4_sensitivities.csv')
    for (h,c),s in h34.groupby(['hypothesis','sensitivity']):r=contrast(s);r.update(hypothesis=h,sensitivity=c);out.append(r)
    for h,s in z.loc[z.hypothesis.ne('H6')].groupby('hypothesis'):
        for excluded,label in [(['6579'],'exclude_pancreatitis6579'),(['6551'],'exclude_Peters_imaged6551')]:
            r=contrast(s.loc[~s.donor_id.isin(excluded)]);r.update(hypothesis=h,sensitivity=label);out.append(r)
        if h in ['H1','H2','H3','H4']:
            # These four legacy same-slide donors lack fully comparable native T
            # scope; remove together to assess the fully typed acquisition stratum.
            r=contrast(s.loc[~s.donor_id.isin(['112','6347','6533','6580'])]);r.update(hypothesis=h,sensitivity='fully_typed_acquisition_stratum');out.append(r)
    src=T/'spatial_H5_H6_sensitivities.csv'
    if src.exists():
        q=read(src.name).rename(columns={'p_value':'p','CI95_low':'ci_low','CI95_high':'ci_high','n_pair':'n_test','n_ND':'n_ref','status':'test_status'})
        for _,r in q.iterrows():
            d=r.to_dict();d['sensitivity']='RNA_'+str(r.get('configuration',''))+'_legacyout'+str(r.get('exclude_five_legacy_panRNA_donors',''))+'_agesex'+str(r.get('age_sex_extra_sensitivity',''));out.append(d)
    s=pd.DataFrame(out);s['family']=FAMILY;s['multiplicity_scope']='diagnostic sensitivity; not a new six-test discovery family';s.to_csv(T/'spatial_hypothesis_sensitivities.csv',index=False)
    return s
def main():
    plan=json.loads((R/'provenance/spatial_hypothesis_plan_v2.json').read_text());z=primary_donors();rows=[]
    for h in plan['hypotheses'][:5]:
        d=z.loc[z.hypothesis.eq(h['id'])];r=contrast(d);r.update(hypothesis=h['id'],name=h['name'],endpoint=h['estimand'],contrast='T1D minus ND',predicted_direction='negative' if h['id']=='H4' else 'positive',family=FAMILY)
        if h['id'] in ['H1','H2']:
            for c in (['vascular_excess','neural_excess'] if h['id']=='H1' else ['fragmentation_z','interface_z']):
                q=d[['disease_status',c]].rename(columns={c:'donor_effect'});cr=contrast(q);r[c+'_disease_effect']=cr['effect']
        rows.append(r)
    h6=read('spatial_H6_linked_architecture_tests.csv');q=h6.loc[h6.configuration.eq('primary') & ~h6.exclude_five_legacy_panRNA_donors.astype(bool) & ~h6.age_sex_extra_sensitivity.astype(bool)]
    assert len(q)==1,'Exactly one fixed H6 primary required'
    a=q.iloc[0];supported=a.status=='estimated' and np.isfinite(a.get('p_value',np.nan))
    rows.append(dict(hypothesis='H6',name=plan['hypotheses'][5]['name'],endpoint=plan['hypotheses'][5]['estimand'],effect=a.get('effect',np.nan),ci_low=a.get('CI95_low',np.nan),ci_high=a.get('CI95_high',np.nan),p=float(a.p_value) if supported else 1.,n_test=a.n_pair,n_ref=0,n_pairs=a.n_pair,nuisance_rank=a.conditioning_design_rank,test_status='estimated' if supported else 'unavailable',reason=a.get('reason',''),contrast='paired conditional donor association',method='within-disease Freedman-Lane residual permutation',predicted_direction='positive',family=FAMILY))
    tests=pd.DataFrame(rows);assert len(tests)==6 and set(tests.hypothesis)==set('H'+str(i) for i in range(1,7))
    tests['q_six']=bh6(tests.p.fillna(1));interpret=[]
    for _,r in tests.iterrows():
        if r.test_status=='unavailable':s='Untestable with current independent-donor support; '+str(r.reason)
        elif r.q_six>=.05:s='Inconclusive for the fixed prediction; nonsignificance is not refutation'
        else:
            predicted=(r.effect<0) if r.hypothesis=='H4' else (r.effect>0)
            s='Supports the specified conditional spatial prediction' if predicted else 'Contradicts the specified directional prediction'
            if predicted and r.hypothesis=='H2' and not (r.fragmentation_z_disease_effect>0 and r.interface_z_disease_effect>0):s='Joint score differs, but component directions do not jointly support beta fragmentation plus immune interface'
            if predicted and r.hypothesis=='H1':s+='; relative vascular-versus-neural alignment only'
            if predicted and r.hypothesis=='H4':s+='; higher-order term alone does not establish net checkpoint slope weakening'
        interpret.append(s)
    tests['interpretation']=interpret;tests.to_csv(T/'spatial_hypothesis_tests.csv',index=False)
    sensitivities=sensitivity_donors(z)
    fig,axes=plt.subplots(2,3,figsize=(15,8.3),layout='constrained')
    titles=['Boundary configuration','Beta topology','HLA-A domains','Checkpoint interaction','Organizer motifs','pLN–islet architecture']
    for ax,h,title in zip(axes.flat,['H1','H2','H3','H4','H5','H6'],titles):
        r=tests.loc[tests.hypothesis.eq(h)].iloc[0];d=z.loc[z.hypothesis.eq(h)]
        if h!='H6':
            for x,stage in enumerate(['ND','AAb','T1D']):
                s=d.loc[d.disease_status.eq(stage)&d.donor_effect.notna()].sort_values('donor_id');j=np.linspace(-.18,.18,len(s))
                for dx,(_,v) in zip(j,s.iterrows()):ax.scatter(x+dx,v.donor_effect,color=COL[stage],s=24);ax.annotate(v.donor_id,(x+dx,v.donor_effect),xytext=(2,2),textcoords='offset points',fontsize=6)
            ax.set_xticks(range(3),['ND','AAb','T1D']);ax.axhline(0,color='gray',lw=.8,ls='--');ax.set_ylabel('Conditional arrangement endpoint')
        else:
            p=T/'spatial_H6_linked_primary_donors.csv'
            if p.exists():
                q=read(p.name)
                for stage,s in q.groupby('disease_status'):
                    ax.scatter(s.H6_residual_rank_pLN,s.H6_residual_rank_H3,color=COL.get(stage,'gray'),s=25,label=stage)
                    for _,v in s.iterrows():ax.annotate(v.donor_id,(v.H6_residual_rank_pLN,v.H6_residual_rank_H3),xytext=(2,2),textcoords='offset points',fontsize=6)
                ax.legend(frameon=False,fontsize=8)
            ax.set_xlabel('pLN architecture, conditional rank residual');ax.set_ylabel('Pancreatic HLA-A domain rank residual')
        ax.set_title(h+' | '+title,fontsize=11);ax.spines[['top','right']].set_visible(False)
        if r.test_status=='unavailable':label='UNTESTABLE\n'+str(r.reason).replace('supported ','')
        else:label=f"effect={r.effect:.3g} | p={r.p:.3g} | BH6 q={r.q_six:.3g}\n"+(f"paired n={int(r.n_test)}" if h=='H6' else f"T1D n={int(r.n_test)}, ND n={int(r.n_ref)}")
        ax.text(.5,1.02,label,transform=ax.transAxes,fontsize=8,ha='center',va='bottom',color='#333333',wrap=True)
    fig.suptitle('Six replacement spatial hypotheses | fixed donor-level family\nAll cells inform local context; only independent donors inform disease significance',fontsize=13)
    for ext in ['png','svg']:fig.savefig(F/f'spatial_hypotheses_v2.{ext}',dpi=300)
    plt.close(fig)
    (T/'spatial_hypothesis_family_manifest.json').write_text(json.dumps(dict(family=FAMILY,script_sha256=sha(__file__),plan_sha256=sha(R/'provenance/spatial_hypothesis_plan_v2.json'),n_primary_tests=6,unavailable_p=1,statistics='Independent donor permutation and bootstrap; BH6. H6 conditional residual permutation by disease.',historical_family='six_hypothesis_tests.csv superseded as active hypotheses; descriptive archive only',source_hashes={p.name:sha(p) for p in [T/'spatial_H1_H2_donor_endpoints.csv',T/'spatial_H3_H4_donor_endpoints.csv',T/'spatial_H5_organizer_donor_endpoints.csv',T/'spatial_H6_linked_architecture_tests.csv']}),indent=2))
    print(tests[['hypothesis','effect','p','q_six','n_test','n_ref','test_status','interpretation']].to_string(index=False))

if __name__=='__main__':main()
