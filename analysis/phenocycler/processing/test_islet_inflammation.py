"""Focused invariants for tissue-supported and tied inflammation ranks."""
import numpy as np
import pandas as pd
from islet_inflammation import inflammation_order, apply_order

def test_zeros_and_ties():
    r,status=inflammation_order([0,0,1,1,9,np.nan])
    assert status=='ordered'
    np.testing.assert_allclose(r,[0,0,4/7,4/7,1,np.nan],equal_nan=True)
    for v,w,tag in [([0,0],[0,0],'all zero'),([0],[0],'all zero'),
                    ([7],[np.nan],'constant positive'),([5,5],[np.nan,np.nan],'constant positive'),
                    ([np.nan],[np.nan],'no assessed'),([0,5],[0,1],'ordered')]:
        result,status=inflammation_order(v)
        assert tag in status
        np.testing.assert_allclose(result,w,equal_nan=True)

def test_region_area_and_donor_isolation():
    z=pd.DataFrame(dict(donor=['x']*4+['y']*2,core_area_um2=[1e6]*6,peri_area_um2=[3e6]*6,
                        n_core=[4,4,4,0,4,4],n_peri=[8]*6,n_core_immune=[0,2,1,0,1,1],
                        n_peri_immune=[0,0,1,0,3,3]))
    o=apply_order(z)
    np.testing.assert_allclose(o.combined_immune_per_mm2,[0,.5,.5,np.nan,1,1],equal_nan=True)
    np.testing.assert_allclose(o.within_donor_density_rank,[0,1,1,np.nan,np.nan,np.nan],equal_nan=True)
    np.testing.assert_allclose(o.pseudotime,[0,.5,.5,np.nan,1,1],equal_nan=True)
    np.testing.assert_allclose(o.pseudotime_rank,[0,3/7,3/7,np.nan,1,1],equal_nan=True)
    assert list(o.burden.astype('string').dropna())==['Zero','Lower','Lower','Upper','Upper']

def test_native_band_excludes_core_neighbors_clusters_and_nontissue(tmp_path):
    import json, math
    import shapely
    from shapely.geometry import box, mapping
    from islet_inflammation import regions
    a,b=box(0,0,10,10),box(17,0,21,4)
    cluster=box(10,7,14,10)
    def candidate(n,g,k):return dict(id=n,geometry=mapping(g),properties={'size_class':k})
    def tissue(n,g):return dict(geometry=mapping(g),properties={'classification':{'name':n}})
    f=tmp_path/'islets.json';f.write_text(json.dumps({'features':[candidate('a',a,'Islet'),candidate('b',b,'Islet'),candidate('c',cluster,'EndocrineCluster')]}))
    # A right-sided tissue cut and an internal pLN exclusion cut into the A ring.
    t=tmp_path/'anno.json';t.write_text(json.dumps({'features':[tissue('Tissue',box(-100,-100,25,100)),tissue('pLN',box(-4,0,-3,4))]}))
    table=pd.DataFrame(dict(islet_id=['a','b'],outline_wkb=[a.wkb,b.wkb],outline_area_um2=[a.area,b.area],geometry_reference=str(f),tissue_reference=str(t),n_core=1,n_peri=1,n_core_immune=[1,0],n_peri_immune=[0,1],geometry_valid=True,complete_20um_tissue_margin=False,core_pancreatic_tissue_fraction=1.))
    features=tmp_path/'features.parquet';table.to_parquet(features,index=False)
    pd.DataFrame(dict(islet=['a','a','b','b'],object_id=['1','2','3','4'],region=['core','peri','core','peri'],refined_class=['Immune','Beta','Beta','Immune'])).to_parquet(tmp_path/'members.parquet',index=False)
    _,out,_=regions(('sample',features,{'pixel_size_um_x':1,'pixel_size_um_y':1},None,None))
    # Independent higher-resolution geometry; includes the union of blockers only once.
    expected=a.buffer(20,quad_segs=1024).intersection(box(-100,-100,25,100)).difference(shapely.union_all([a,b,cluster,box(-4,0,-3,4)]))
    assert abs(out.loc['a','peri_area_um2']-expected.area)<.2
    assert out.loc['a','core_area_um2']==100
