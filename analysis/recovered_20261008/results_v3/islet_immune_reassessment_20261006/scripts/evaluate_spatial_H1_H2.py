#!/usr/bin/env python3
"""Conditional arrangements of current protein cells around and within true islets."""
from pathlib import Path
import json, hashlib, sys, argparse
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
from shapely.geometry import shape
from shapely.affinity import scale
from shapely.strtree import STRtree
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO=Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
OUT=REPO/"results/v3/islet_immune_reassessment_20261006"
TAB=OUT/"tables"; FIG=OUT/"figures"; PROV=OUT/"provenance"
ASS=Path("/data/io60-calibration/pancreas_analysis_assignments_single_exclusions_20261005/manifest.json")
COMP=Path("/data/io60-calibration/composition_production_all_exclusions_20261005/endocrine_objects.parquet")
sys.path.insert(0,str(REPO/"core/phenocycler/scripts/typing"))
from build_islet_geometry import islet_export, _manifest
from subpanel_io import GEOJSON
from context import OBJECT_FOLDER
SEED=20261006; NPERM=199
COLORS={"ND":"#3a6e94","AAb":"#bd8a29","T1D":"#a74151"}
ENDO={"Beta","Alpha","Delta"}
CONFIGS=[("native_primary","cell_class",1.,6,True),
         ("refined_identity","refined_class",1.,6,True),
         ("half_cell_scale","cell_class",.5,6,True),
         ("double_cell_scale","cell_class",2.,6,True),
         ("graph_k3","cell_class",1.,3,True),
         ("graph_k10","cell_class",1.,10,True),
         ("unstratified_label_null","cell_class",1.,6,False)]
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()
def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,np.generic):return clean(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    if isinstance(v,Path):return str(v)
    return v
def savejson(path,value):path.write_text(json.dumps(clean(value),indent=2,allow_nan=False))
def savefig(fig,name):
    fig.savefig(FIG/(name+".png"),dpi=300,bbox_inches="tight")
    fig.savefig(FIG/(name+".svg"),bbox_inches="tight");plt.close(fig)
def label_seed(*parts):
    return int(hashlib.sha256("|".join(map(str,parts)).encode()).hexdigest()[:8],16)
def load_geometry(native,allowed):
    pix,stems=_manifest();pixel=float(pix[native]);stem=stems[native]
    src=islet_export(stem); out={}; audit={"geometry_path":str(src),"geometry_sha256":sha(src),"pixel_size_um":pixel,"source_features":0,"invalid_or_empty":0,"non_Islet":0,"not_composition_eligible":0}
    for f in json.loads(src.read_text())["features"]:
        audit["source_features"]+=1; uid=str(f.get("id","")); pr=f.get("properties",{})
        if pr.get("size_class")!="Islet":audit["non_Islet"]+=1;continue
        if uid not in allowed:audit["not_composition_eligible"]+=1;continue
        if not f.get("geometry"):audit["invalid_or_empty"]+=1;continue
        g=scale(shape(f["geometry"]),xfact=pixel,yfact=pixel,origin=(0,0))
        if g.is_empty or not g.is_valid:audit["invalid_or_empty"]+=1;continue
        if uid in out:raise ValueError("Duplicate islet UUID")
        shapely.prepare(g);out[uid]=g
    audit["valid_true_Islets"]=len(out)
    return out,stem,pixel,audit
def structure_geometry(stem,name,pixel):
    files=sorted(p for p in (Path(GEOJSON)/OBJECT_FOLDER[name]).glob("*.geojson") if stem in str(p))
    audit={"structure_name":name,"matches":len(files),"folder":str(Path(GEOJSON)/OBJECT_FOLDER[name]),"source_features":0,"invalid_or_empty":0,"unsupported":0}
    if len(files)!=1:return None,audit
    src=files[0];audit.update(path=str(src),sha256=sha(src));out=[]
    for f in json.loads(src.read_text())["features"]:
        audit["source_features"]+=1
        if not f.get("geometry"):audit["invalid_or_empty"]+=1;continue
        g=scale(shape(f["geometry"]),xfact=pixel,yfact=pixel,origin=(0,0))
        if g.is_empty or not g.is_valid:audit["invalid_or_empty"]+=1;continue
        if g.geom_type not in {"Polygon","MultiPolygon"}:audit["unsupported"]+=1;continue
        out.append(g)
    audit["valid_objects"]=len(out)
    return out,audit
def nearest_dist(tree,points):
    if tree is None or not len(points):return np.full(len(points),np.nan)
    ii,dd=tree.query_nearest(points,all_matches=False,return_distance=True)
    result=np.full(len(points),np.nan);result[ii[0]]=dd;return result
def core_groups(a,geos):
    """Authoritative current true-Islet parent UUID; validate each centroid geometry."""
    parents=a.endocrine_object_id.astype("string").fillna("").to_numpy(str)
    mask=a.endocrine_object_class.eq("Islet").fillna(False).to_numpy(bool)
    out={};audit={"n_parent_core":int(mask.sum()),"parent_absent_current_geometry":0,"parent_centroid_not_covered":0}
    for uid,ix in pd.Series(np.where(mask)[0]).groupby(parents[mask],sort=False):
        idx=ix.to_numpy(int)
        if uid not in geos:audit["parent_absent_current_geometry"]+=len(idx);continue
        p=shapely.points(a.iloc[idx][["X_centroid","Y_centroid"]].to_numpy(float))
        good=shapely.covers(geos[uid],p)
        audit["parent_centroid_not_covered"]+=int((~good).sum())
        out[uid]=idx[good]
    return out,audit
def ring_segments(g,exterior_only=False):
    starts=[];ends=[]
    for p in ([g] if g.geom_type=="Polygon" else list(g.geoms) if g.geom_type=="MultiPolygon" else []):
        for ring in ([p.exterior] if exterior_only else [p.exterior]+list(p.interiors)):
            coords=np.asarray(ring.coords,float);starts.append(coords[:-1]);ends.append(coords[1:])
    if not starts:return np.empty((0,2)),np.empty((0,2))
    st=np.concatenate(starts);en=np.concatenate(ends);valid=np.linalg.norm(en-st,axis=1)>0
    return st[valid],en[valid]
def exact_exterior_projection(g,points,max_distance):
    st,en=ring_segments(g,True);vec=en-st;length=np.linalg.norm(vec,axis=1)
    cumulative=np.r_[0,np.cumsum(length)[:-1]]
    tree=STRtree(shapely.linestrings(np.stack([st,en],axis=1)))
    ij,dd=tree.query_nearest(points,all_matches=True,return_distance=True,max_distance=max_distance)
    # GEOS line_locate_point resolves equal nearest segments by first ring segment.
    order=np.lexsort((ij[1],ij[0]));ii,jj=ij[:,order];dd=dd[order]
    _,first=np.unique(ii,return_index=True);ii=ii[first];jj=jj[first];dd=dd[first]
    xy=np.column_stack([shapely.get_x(points[ii]),shapely.get_y(points[ii])])
    fraction=np.clip(np.sum((xy-st[jj])*vec[jj],axis=1)/(length[jj]**2),0,1)
    arc=cumulative[jj]+fraction*length[jj]
    return ii,np.asarray(dd,float),np.asarray(arc,float)
def boundary_cells(a,geos,cellscale):
    """Exact nearest boundary-segment search; preserves polygon/parent/tie semantics."""
    ids=list(geos);gg=np.asarray(list(geos.values()),dtype=object);tree=STRtree(gg)
    xy=a[["X_centroid","Y_centroid"]].to_numpy(float);pts=shapely.points(xy)
    radius=-2*cellscale*np.log(1e-8)
    # Outside all polygons, nearest polygon distance equals nearest ring-segment distance.
    starts=[];ends=[];owners=[]
    for j,g in enumerate(gg):
        st,en=ring_segments(g)
        starts.append(st);ends.append(en);owners.append(np.full(len(st),j,int))
    st=np.concatenate(starts);en=np.concatenate(ends);owner=np.concatenate(owners)
    segmenttree=STRtree(shapely.linestrings(np.stack([st,en],axis=1)))
    ij,dd=segmenttree.query_nearest(pts,all_matches=False,return_distance=True,max_distance=radius)
    near=np.full(len(a),-1,int);dist=np.full(len(a),np.nan)
    near[ij[0]]=owner[ij[1]];dist[ij[0]]=dd
    # Exact covered-by test handles deep core and holes without nearest-polygon traversal.
    inside=tree.query(pts,predicate="covered_by")
    ambiguity=0
    if inside.shape[1]:
        order=np.lexsort((inside[1],inside[0]));ii,jj=inside[:,order]
        unique,first,counts=np.unique(ii,return_index=True,return_counts=True)
        near[unique]=jj[first];dist[unique]=0
        multi=unique[counts>1];ambiguity=len(multi)
        # Match the original STRtree's arbitrary tie policy only on rare overlapping cores.
        if len(multi):
            kk,zz=tree.query_nearest(pts[multi],all_matches=False,return_distance=True)
            near[multi[kk[0]]]=kk[1];dist[multi[kk[0]]]=zz
    mp={x:i for i,x in enumerate(ids)}
    parent=a.endocrine_object_id.astype("string").fillna("").map(mp).fillna(-1).to_numpy(int)
    parentcore=a.endocrine_object_class.eq("Islet").fillna(False).to_numpy(bool)&(parent>=0)
    near[parentcore]=parent[parentcore]
    candidates=np.isfinite(dist)&(dist<=radius)&(near>=0)
    # Outside equal-distance ties across different polygons: preserve original choice exactly.
    # All nearest segment matches are only queried at the candidate points without parent.
    unresolved=np.where(candidates&~parentcore)[0]
    outer_ties=0
    if len(unresolved):
        tt,td=segmenttree.query_nearest(pts[unresolved],all_matches=True,return_distance=True,max_distance=radius)
        sr=np.lexsort((owner[tt[1]],tt[0]));ti=tt[0,sr];to=owner[tt[1,sr]]
        change=np.r_[True,(ti[1:]!=ti[:-1])|(to[1:]!=to[:-1])]
        uniq_owner_count=np.bincount(ti[change],minlength=len(unresolved))
        tied=unresolved[uniq_owner_count>1];tied=tied[dist[tied]>0];outer_ties=len(tied)
        if len(tied):
            kk,zz=tree.query_nearest(pts[tied],all_matches=False,return_distance=True,max_distance=radius)
            near[tied[kk[0]]]=kk[1];dist[tied[kk[0]]]=zz
    groups={}
    for i,vals in pd.Series(np.where(candidates)[0]).groupby(near[candidates],sort=False):
        uid=ids[int(i)];g=geos[uid]
        if g.geom_type!="Polygon":continue
        idx=vals.to_numpy(int)
        selected,radial,arc=exact_exterior_projection(g,pts[idx],radius)
        keep=np.exp(-radial/(2*cellscale))>=1e-8
        groups[uid]=(idx[selected[keep]],radial[keep],arc[keep])
    return groups,{"n_numeric_truncation_cells":int((~candidates).sum()),"maximum_kernel_scale_um":2*cellscale,"maximum_per_cell_discarded_weight":1e-8,"radial_kernel":"exp(-absolute exterior boundary distance / donor mean cell diameter)","nearest_algorithm":"exact STRtree of original boundary segments; exact covered-by and original polygon-tie fallback; exact segment arclength, no simplification","n_boundary_segments":len(st),"n_overlap_core_ties":ambiguity,"n_cross_polygon_outer_ties":outer_ties}

def circular_scores(residual,profile):
    # Enumerate every nonzero shift; shift null tests orientation, not abundance.
    corr=np.fft.ifft(np.conj(np.fft.fft(residual))*np.fft.fft(profile)).real
    obs=float(np.dot(residual,profile));null=corr[1:]
    assert np.isclose(corr[0],obs,rtol=1e-7,atol=1e-9)
    return obs,null
def h1_islet(uid,g,group,labels,ell,vascular,neural,config,base):
    row=dict(base,hypothesis="H1",configuration=config,kernel_scale_um=ell,endpoint=np.nan,reason="",n_boundary_bins=0,n_core_cells=base["n_core_cells"])
    if g.geom_type!="Polygon":row["reason"]="multipart boundary has no single circular arclength";return row,None
    if group is None:row["reason"]="no eligible boundary support";return row,None
    if vascular is None or neural is None:row["reason"]="missing valid vascular or neural process source";return row,None
    idx,radial,arc=group;w=np.exp(-radial/ell);im=labels[idx]=="Immune"
    sw=float(w[im].sum());neff=float(sw**2/np.sum(w[im]**2)) if sw else 0.
    row.update(n_boundary_candidates=len(idx),n_immune_candidates=int(im.sum()),immune_kernel_mass=sw,immune_effective_count=neff,allcell_kernel_mass=float(w.sum()))
    if int(im.sum())<3 or neff<3 or sw<1:row["reason"]="requires >=3 actual/effective immune cells and >=1 immune-cell kernel mass";return row,None
    perimeter=float(g.exterior.length);n=max(64,int(np.ceil(perimeter/ell)));n=min(n,4096)
    row.update(n_boundary_bins=n,boundary_perimeter_um=perimeter,arc_step_um=perimeter/n,arc_cap_reached=n==4096)
    bins=np.minimum((arc/perimeter*n).astype(int),n-1)
    allhist=np.bincount(bins,weights=w,minlength=n);ihist=np.bincount(bins[im],weights=w[im],minlength=n)
    sigma=max(ell*n/perimeter,.5)
    allprofile=gaussian_filter1d(allhist,sigma,mode="wrap");immuneprofile=gaussian_filter1d(ihist,sigma,mode="wrap")
    r=immuneprofile/immuneprofile.sum()-allprofile/allprofile.sum()
    s=(np.arange(n)+.5)*perimeter/n;p=shapely.line_interpolate_point(g.exterior,s)
    v=np.exp(-nearest_dist(vascular,p)/ell);ner=np.exp(-nearest_dist(neural,p)/ell)
    row.update(vascular_profile_mean=float(v.mean()),neural_profile_mean=float(ner.mean()),vascular_profile_sd=float(v.std()),neural_profile_sd=float(ner.std()))
    if not np.isfinite(v).all() or not np.isfinite(ner).all():row["reason"]="nonfinite structural distance";return row,None
    if np.std(v-ner)<1e-8:row["reason"]="no varying vascular-minus-neural profile";return row,None
    obs,null=circular_scores(r,v-ner);ov,nv=circular_scores(r,v);on,nn=circular_scores(r,ner)
    sd=float(null.std(ddof=1))
    if sd<1e-10:row["reason"]="circular arrangement null has no variance";return row,None
    row.update(endpoint=obs-float(null.mean()),observed_score=obs,null_mean=float(null.mean()),null_sd=sd,
               null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)),
               n_null=len(null),conditional_rank_p=float((1+(null>=obs-1e-12).sum())/(len(null)+1)),
               vascular_excess=ov-float(nv.mean()),neural_excess=on-float(nn.mean()),
               null_standardized_excess=(obs-float(null.mean()))/sd,reason="available")
    cache={"arc_fraction":s/perimeter,"residual":r,"immune_profile":immuneprofile/immuneprofile.sum(),"background_profile":allprofile/allprofile.sum(),"vascular_profile":v,"neural_profile":ner,"null":null,"obs":obs,"xy":np.column_stack([shapely.get_x(p),shapely.get_y(p)]),"candidate_indices":idx,"radial":radial}
    return row,cache
def graph_edges(xy,g,k,ell):
    n=len(xy)
    if n<2:return np.array([],int),np.array([],int),np.array([],float),0
    d,j=cKDTree(xy).query(xy,k=min(k+1,n))
    if j.ndim==1:j=j[:,None]
    i=np.repeat(np.arange(n),j.shape[1]-1);jj=j[:,1:].reshape(-1)
    pairs=np.unique(np.sort(np.column_stack([i,jj]),axis=1),axis=0)
    pairs=pairs[pairs[:,0]!=pairs[:,1]]
    lines=shapely.linestrings(xy[pairs]);within=shapely.covers(g,lines)
    rejected=int((~within).sum());pairs=pairs[within]
    w=np.exp(-np.linalg.norm(xy[pairs[:,0]]-xy[pairs[:,1]],axis=1)/ell)
    return pairs[:,0],pairs[:,1],w,rejected
def fragment(labels,edge):
    i,j,w,_=edge;b=labels=="Beta"
    denom=float((w*(b[i].astype(float)+b[j].astype(float))).sum())
    if denom<=0:return np.nan
    return 1-float(2*w[b[i]&b[j]].sum())/denom
def interface(labels,edge):
    i,j,w,_=edge;b=labels=="Beta";immune=labels=="Immune"
    denom=float((w*(b[i].astype(float)+b[j].astype(float))).sum())
    if denom<=0:return np.nan
    return float(w[(b[i]&immune[j])|(b[j]&immune[i])].sum())/denom
def h2_islet(uid,g,idx,a,labels,ell,k,strat,config,base):
    row=dict(base,hypothesis="H2",configuration=config,kernel_scale_um=ell,graph_k=k,radial_stratified=strat,endpoint=np.nan,reason="")
    if len(idx)<2:row["reason"]="no admitted core cell graph";return row,None
    lab=labels[idx];xy=a.iloc[idx][["X_centroid","Y_centroid"]].to_numpy(float);end=np.isin(lab,list(ENDO))
    ei=np.where(end)[0];elab=lab[end];nB=int((elab=="Beta").sum());nnon=int((elab!="Beta").sum());nI=int((lab=="Immune").sum())
    row.update(n_beta=nB,n_alpha=int((elab=="Alpha").sum()),n_delta=int((elab=="Delta").sum()),n_immune=nI,n_endocrine=len(ei),beta_fraction=nB/len(ei) if len(ei) else np.nan)
    if len(ei)<10 or nB<3 or nnon<3:row["reason"]="requires >=10 endocrine, >=3 beta and >=3 other endocrine cells";return row,None
    if nI<1:row["reason"]="no native/refined immune cell in this own islet core";return row,None
    eg=graph_edges(xy[end],g,k,ell);ag=graph_edges(xy,g,k,ell)
    row.update(n_endocrine_edges=len(eg[0]),n_core_edges=len(ag[0]),n_endocrine_edges_rejected_outside=eg[3],n_core_edges_rejected_outside=ag[3])
    depth=np.asarray(shapely.distance(shapely.points(xy[end]),g.boundary),float)
    # Equal-count rank strata have no physical bands. Stable tie order uses source cell ID.
    ns=min(5,max(1,len(ei)//10)) if strat else 1
    order=np.lexsort((a.iloc[idx[ei]].cell_id.to_numpy(str),depth))
    strata=np.empty(len(ei),int)
    for z,ss in enumerate(np.array_split(order,ns)):strata[ss]=z
    blocks=[np.where(strata==z)[0] for z in range(ns)]
    movable=sum(len(v) for v in blocks if len(np.unique(elab[v]))>1)
    row.update(n_radial_rank_strata=ns,n_exchangeable_endocrine_positions=movable)
    if movable<2:row["reason"]="radial conditional labels have no exchangeable positions";return row,None
    observed=np.array([fragment(elab,eg),interface(lab,ag)])
    if not np.isfinite(observed).all():row["reason"]="no beta graph denominator";return row,None
    rng=np.random.default_rng(label_seed(SEED,base["donor_id"],uid,config));null=np.empty((NPERM,2))
    for z in range(NPERM):
        pp=elab.copy()
        for block in blocks:pp[block]=rng.permutation(pp[block])
        new=lab.copy();new[ei]=pp;null[z]=[fragment(pp,eg),interface(new,ag)]
    mean=null.mean(axis=0);sd=null.std(axis=0,ddof=1)
    row.update(observed_fragmentation=float(observed[0]),observed_beta_immune_interface=float(observed[1]),null_fragmentation_mean=float(mean[0]),null_interface_mean=float(mean[1]),null_fragmentation_sd=float(sd[0]),null_interface_sd=float(sd[1]),n_null=NPERM)
    if not np.isfinite(sd).all() or (sd<1e-10).any():row["reason"]="one or both topology null components have no variance";return row,None
    z=(observed-mean)/sd;zn=(null-mean)/sd;joint=zn.mean(axis=1);obs=float(z.mean())
    row.update(endpoint=obs,fragmentation_z=float(z[0]),interface_z=float(z[1]),fragmentation_excess=float(observed[0]-mean[0]),interface_excess=float(observed[1]-mean[1]),null_component_correlation=float(np.corrcoef(null.T)[0,1]),null_q025=float(np.quantile(joint,.025)),null_q975=float(np.quantile(joint,.975)),conditional_rank_p=float((1+(joint>=obs-1e-12).sum())/(NPERM+1)),both_components_positive=bool((z>0).all()),reason="available")
    cache={"xy":xy,"labels":lab,"endocrine_indices":ei,"endocrine_edges":eg,"core_edges":ag,"strata":strata,"null_components":null,"null_joint":joint,"obs":obs,"fragmentation_z":z[0],"interface_z":z[1]}
    return row,cache
def load_cells(rec):
    columns=["cell_id","donor_id","disease_status","cell_class","refined_class","immune_class","analysis_eligible","qc_analysis_eligible","endocrine_object_class","endocrine_object_id","X_centroid","Y_centroid"]
    assert sha(rec["path"])==rec["sha256"],"assignment changed"
    a=pd.read_parquet(rec["path"],columns=columns);assert not a.cell_id.duplicated().any()
    elig=a.analysis_eligible.fillna(False).to_numpy(bool)&a.qc_analysis_eligible.fillna(False).to_numpy(bool)
    assert int(elig.sum())==rec["typing_analysis_cells"]
    names=set(pq.ParquetFile(rec["native_source_path"]).schema.names)
    area="fundamentals_cell_area_um2" if "fundamentals_cell_area_um2" in names else "cell_area"
    p=pd.read_parquet(rec["native_source_path"],columns=["object_id","nucleus_area_um2",area])
    assert not p.object_id.duplicated().any()
    ix=pd.Index(p.object_id.astype(str)).get_indexer(a.cell_id.astype(str))
    assert (ix>=0).all(),"current cell UUID mismatch"
    nuc=pd.to_numeric(p.iloc[ix].nucleus_area_um2,errors="coerce").to_numpy(float)
    ca=pd.to_numeric(p.iloc[ix][area],errors="coerce").to_numpy(float)
    finitecoord=np.isfinite(a[["X_centroid","Y_centroid"]].to_numpy(float)).all(axis=1)
    nucok=np.isfinite(nuc)&(nuc>0);valid=elig&nucok&finitecoord
    audit={"source_cells":len(a),"typing_analysis_cells":int(elig.sum()),"invalid_nucleus_in_admitted":int((elig&~nucok).sum()),"invalid_coordinates_in_admitted":int((elig&~finitecoord).sum()),"n_spatial_admitted":int(valid.sum()),"area_source_field":area,"native_source_path":rec["native_source_path"],"native_source_sha256":sha(rec["native_source_path"])}
    a=a.loc[valid].reset_index(drop=True);ca=ca[valid]
    good=np.isfinite(ca)&(ca>0)
    if not good.any():raise ValueError("No positive cell area for continuous donor scale")
    ell=float(np.mean(2*np.sqrt(ca[good]/np.pi)))
    audit.update(mean_cell_equivalent_diameter_um=ell,n_cell_area_support=int(good.sum()),n_area_missing=int((~good).sum()))
    return a,ell,audit
def summarize(rows,rec,disease,config,base):
    out=dict(donor_id=str(rec["donor_id"]),disease_status=disease,configuration=config,**base)
    for h in ["H1","H2"]:
        rr=pd.DataFrame([x for x in rows if x["hypothesis"]==h]);supported=rr.loc[rr.reason.eq("available")] if len(rr) else rr
        out[h+"_n_islets"]=len(rr);out[h+"_n_supported_islets"]=len(supported)
        out[h+"_endpoint"]=float(supported.endpoint.mean()) if len(supported) else np.nan
        out[h+"_availability_reason"]="available" if len(supported) else "|".join(sorted(rr.reason.unique())) if len(rr) else "no valid true-Islet geometry"
        out[h+"_conditional_rank_fraction_below_05"]=float((supported.conditional_rank_p<=.05).mean()) if len(supported) else np.nan
        for field in (["vascular_excess","neural_excess","null_standardized_excess"] if h=="H1" else ["fragmentation_z","interface_z","fragmentation_excess","interface_excess","beta_fraction"]):
            out[h+"_"+field]=float(supported[field].mean()) if len(supported) else np.nan
    return out
def donor(rec,comp):
    print('H1/H2 '+str(rec['donor_id'])+': loading current cells',flush=True)
    a,ell,aa=load_cells(rec);disease=str(a.disease_status.mode().iloc[0]).replace("AAb+","AAb")
    cp=comp.loc[comp.donor_id.astype(str).eq(str(rec["donor_id_native"])) & comp.original_object_class.eq("Islet") & comp.composition_analysis_eligible.fillna(False)]
    geos,stem,pixel,ga=load_geometry(str(rec["donor_id_native"]),set(cp.endocrine_object_id.astype(str)))
    if not geos:raise ValueError("No valid current eligible true-Islet geometry")
    print('H1/H2 '+str(rec['donor_id'])+': exact boundary segment projection',flush=True)
    groups,cg=core_groups(a,geos);bg,ba=boundary_cells(a,geos,ell)
    print('H1/H2 '+str(rec['donor_id'])+': loading structural geometry',flush=True)
    vascular,va=structure_geometry(stem,"Blood vessel",pixel);neural,na=structure_geometry(stem,"Nerve",pixel)
    vt=STRtree(vascular) if vascular else None;nt=STRtree(neural) if neural else None
    rows=[];summary=[];diag=[];h1memo={};labels={c:a[c].astype("string").fillna("").to_numpy(str) for c in ["cell_class","refined_class"]}
    for config,labelcol,factor,k,strat in CONFIGS:
        print('H1/H2 '+str(rec['donor_id'])+': '+config,flush=True)
        cr=[];lc=labels[labelcol]
        for uid,g in geos.items():
            core=groups.get(uid,np.array([],int));lab=lc[core]
            base={"donor_id":str(rec["donor_id"]),"disease_status":disease,"islet_uuid":uid,"islet_area_um2":float(g.area),"geometry_type":g.geom_type,"n_core_cells":len(core),"n_core_endocrine":int(np.isin(lab,list(ENDO)).sum()),"n_core_immune":int((lab=="Immune").sum())}
            hkey=(uid,labelcol,factor)
            if hkey not in h1memo:h1memo[hkey]=h1_islet(uid,g,bg.get(uid),lc,ell*factor,vt,nt,config,base)
            h1,d1=h1memo[hkey];h1=dict(h1,configuration=config)
            h2,d2=h2_islet(uid,g,core,a,lc,ell*factor,k,strat,config,base)
            cr.extend([h1,h2])
            if config=="native_primary":
                if d1 is not None:diag.append(("H1",uid,g,h1,d1))
                if d2 is not None:diag.append(("H2",uid,g,h2,d2))
        rows.extend(cr);summary.append(summarize(cr,rec,disease,config,{"n_spatial_admitted":len(a),"mean_cell_equivalent_diameter_um":ell,"n_true_Islets":len(geos)}))
    # Example chosen by median area among supported islets, never by effect.
    examples={}
    for h in ["H1","H2"]:
        candidates=[x for x in diag if x[0]==h]
        if candidates:
            median=float(np.median([x[3]["islet_area_um2"] for x in candidates]))
            example=min(candidates,key=lambda x:(abs(x[3]["islet_area_um2"]-median),x[1]))
            examples[h]=example
    provenance={"donor_id":rec["donor_id"],"disease_status":disease,"assignment_path":rec["path"],"assignment_sha256":rec["sha256"],"cell_QC":aa,"islet_geometry":ga,"core_membership_QC":cg,"boundary_sampling_QC":ba,"vascular":va,"neural":na,"image_stem":stem,"label_source":"current selected REDSEA/UniFORM production native calls plus all-exclusion overlay; refined separate"}
    # Save cell-scale examples to avoid rereading all 25 donors for plotting.
    maps={}
    for h,ex in examples.items():
        _,uid,g,row,cache=ex;buf=g.buffer(4*ell);ix=np.unique(np.r_[groups.get(uid,[]),cache.get("candidate_indices",[])]).astype(int)
        coords=a.iloc[ix][["X_centroid","Y_centroid"]].to_numpy(float)
        inside=shapely.covers(buf,shapely.points(coords));ix=ix[inside]
        maps[h]={"donor_id":rec["donor_id"],"disease_status":disease,"islet_uuid":uid,"islet_wkb":g.wkb,"row":row,"cache":cache,"cell_xy":a.iloc[ix][["X_centroid","Y_centroid"]].to_numpy(float),"cell_labels":labels["cell_class"][ix],"cell_ids":a.iloc[ix].cell_id.to_numpy(str),
                 "vessel_wkb":[vascular[i].wkb for i in STRtree(vascular).query(buf,predicate="intersects")] if vascular else [],
                 "nerve_wkb":[neural[i].wkb for i in STRtree(neural).query(buf,predicate="intersects")] if neural else []}
    # numpy diagnostic arrays written losslessly; figures can be rerendered independently.
    import pickle
    with open(PROV/("spatial_H1_H2_example_donor_"+str(rec["donor_id"])+".pkl"),"wb") as f:pickle.dump(maps,f,protocol=5)
    pd.DataFrame(rows).to_csv(TAB/("spatial_H1_H2_islets_donor_"+str(rec["donor_id"])+".csv"),index=False)
    savejson(PROV/("spatial_H1_H2_sources_donor_"+str(rec["donor_id"])+".json"),provenance)
    print("Completed H1/H2 "+str(rec["donor_id"])+" "+disease+": "+str(len(geos))+" islets; primary supports "+str(summary[0]["H1_n_supported_islets"])+"/"+str(summary[0]["H2_n_supported_islets"]),flush=True)
    return rows,summary,provenance
def outline(ax,g,color,lw=.8):
    if g.geom_type=="Polygon":
        x,y=g.exterior.xy;ax.plot(x,y,color=color,lw=lw)
        for ring in g.interiors:x,y=ring.xy;ax.plot(x,y,color=color,lw=lw)
    elif g.geom_type=="MultiPolygon":
        for part in g.geoms:outline(ax,part,color,lw)
def donor_panel(ax,summary,key,title,ylabel):
    for j,stage in enumerate(["ND","AAb","T1D"]):
        s=summary.loc[summary.disease_status.eq(stage)].sort_values("donor_id");x=j+np.linspace(-.2,.2,len(s))
        for xx,(_,r) in zip(x,s.iterrows()):
            value=r[key]
            if np.isfinite(value):ax.scatter(xx,value,s=35,color=COLORS[stage]);ax.annotate(r.donor_id,(xx,value),xytext=(2,2),textcoords="offset points",fontsize=6)
    ax.set_xticks([0,1,2],["ND","AAb","T1D"]);ax.axhline(0,color="grey",ls="--",lw=.8);ax.set_title(title,fontsize=11);ax.set_ylabel(ylabel,fontsize=9);ax.spines[["top","right"]].set_visible(False)
def plot(summary):
    import pickle
    primary=summary.loc[summary.configuration.eq("native_primary")].copy()
    for h,name in [("H1","spatial_H1_neurovascular_interfaces"),("H2","spatial_H2_beta_topology")]:
        chosen=[]
        for stage in ["ND","T1D"]:
            s=primary.loc[primary.disease_status.eq(stage)&primary[h+"_n_supported_islets"].gt(0)]
            if len(s):
                m=float(s[h+"_n_supported_islets"].median())
                r=s.assign(delta=(s[h+"_n_supported_islets"]-m).abs()).sort_values(["delta","donor_id"]).iloc[0]
                p=PROV/("spatial_H1_H2_example_donor_"+str(r.donor_id)+".pkl")
                with open(p,"rb") as f:chosen.append(pickle.load(f)[h])
        fig,axes=plt.subplots(2,3,figsize=(16,9));fig.subplots_adjust(hspace=.52,wspace=.38)
        if h=="H1":
            donor_panel(axes[0,0],primary,"H1_endpoint","Vascular minus neural boundary alignment","Mean immune-exposure arrangement excess")
            donor_panel(axes[0,1],primary,"H1_vascular_excess","Vascular component","Mean immune-exposure arrangement excess")
            donor_panel(axes[0,2],primary,"H1_neural_excess","Neural component (direction diagnostic)","Mean immune-exposure arrangement excess")
        else:
            donor_panel(axes[0,0],primary,"H2_endpoint","Joint beta topology arrangement excess","Mean conditional standardized excess")
            donor_panel(axes[0,1],primary,"H2_fragmentation_z","Beta-neighbor fragmentation","Mean conditional standardized excess")
            donor_panel(axes[0,2],primary,"H2_interface_z","Beta–immune spatial adjacency","Mean conditional standardized excess")
        for n,ex in enumerate(chosen):
            ax=axes[1,n];g=shapely.from_wkb(ex["islet_wkb"])
            if h=="H1":
                for w in ex["vessel_wkb"]:outline(ax,shapely.from_wkb(w),"#4b8fbc",.5)
                for w in ex["nerve_wkb"]:outline(ax,shapely.from_wkb(w),"#c29a43",.6)
                for label,color in [("Immune","#a74151"),("Beta","#8f8f8f"),("Alpha","#bbbbbb"),("Delta","#d8d8d8")]:
                    q=ex["cell_xy"][ex["cell_labels"]==label];ax.scatter(q[:,0],q[:,1],s=10 if label=="Immune" else 3,color=color,label=label,zorder=3)
            else:
                cc=ex["cache"];xy=cc["xy"];ii,jj,ww,_=cc["core_edges"]
                for i,j,w in zip(ii,jj,ww):
                    if (cc["labels"][i]=="Beta" and cc["labels"][j]=="Immune") or (cc["labels"][j]=="Beta" and cc["labels"][i]=="Immune"):
                        ax.plot(xy[[i,j],0],xy[[i,j],1],color="#b1555e",lw=.7,alpha=.65)
                for label,color in [("Beta","#3f7e92"),("Alpha","#bd8a29"),("Delta","#77558b"),("Immune","#a74151")]:
                    q=xy[cc["labels"]==label];ax.scatter(q[:,0],q[:,1],s=11,color=color,label=label)
            outline(ax,g,"#222222",1.1);ax.set_aspect("equal");ax.invert_yaxis();ax.set_xlabel("Protein section x (um)",fontsize=9);ax.set_ylabel("Protein section y (um)",fontsize=9)
            ax.set_title(str(ex["donor_id"])+" "+ex["disease_status"]+" | median-area supported islet\n"+ex["islet_uuid"],fontsize=8)
            ax.legend(fontsize=6,frameon=False,ncol=2,loc="best")
        ax=axes[1,2]
        if chosen:
            ex=chosen[-1];cc=ex["cache"]
            if h=="H1":
                ax.plot(cc["arc_fraction"],cc["vascular_profile"],label="Vascular proximity",color="#4b8fbc")
                ax.plot(cc["arc_fraction"],cc["neural_profile"],label="Neural proximity",color="#c29a43")
                ax2=ax.twinx();ax2.plot(cc["arc_fraction"],cc["residual"],color="#a74151",label="Immune minus all-cell exposure")
                ax2.set_ylabel("Immune exposure residual",fontsize=8)
                ax.set_xlabel("Fraction of own islet outer arclength");ax.set_ylabel("exp(−distance / mean cell diameter)",fontsize=8)
                ax.set_title("Continuous boundary evidence; same example",fontsize=10);ax.legend(fontsize=7,loc="upper right")
            else:
                ax.hist(cc["null_joint"],bins=20,color="#b7c2c8");ax.axvline(cc["obs"],color="#a74151",lw=1.5)
                ax.set_title("Own-islet radial label null; same example",fontsize=10);ax.set_xlabel("Joint standardized topology score");ax.set_ylabel("Conditional null arrangements")
        else:ax.text(.5,.5,"No supported diagnostic example",ha="center",va="center",transform=ax.transAxes)
        fig.suptitle("H1: neurovascular boundary alignment" if h=="H1" else "H2: beta topology beyond cell depletion",fontsize=14)
        fig.text(.5,.012,"Every dot is one donor; current selected protein identities, nucleated eligible cells and true-Islet geometry. Examples chosen by support count then median area.\n"+
                 ("Continuous cell-scale kernels; exact circular structural-profile shifts preserve immune abundance, radial positions and all-cell exposure. Missing support is unavailable." if h=="H1" else "Protein-centroid neighbor graphs are spatial adjacency estimates. Radial label permutations preserve beta/alpha/delta counts and fixed immune positions; depleted or invariant islets are unavailable."),ha="center",fontsize=8)
        savefig(fig,name)
def main():
    parser=argparse.ArgumentParser();parser.add_argument("--plot-only",action="store_true");parser.add_argument("--donors",nargs="*");parser.add_argument("--resume",action="store_true");parser.add_argument("--workers",type=int,default=1);args=parser.parse_args()
    if args.plot_only:
        plot(pd.read_csv(TAB/"spatial_H1_H2_sensitivities.csv",dtype={"donor_id":str}));return
    for p in [TAB,FIG,PROV]:p.mkdir(parents=True,exist_ok=True)
    plan={"status":"frozen before computation","hypotheses":{
          "H1":{"estimand":"Equal-islet mean immune-versus-all-cell boundary exposure alignment with vascular minus neural process proximity; T1D minus ND positive",
                "kernel":"exp(-absolute outer-boundary distance / donor mean equivalent cell diameter), periodic arclength Gaussian smoothing at the same scale",
                "conditioning":"Current islet shape, admitted cell positions and counts, immune positions and radial distances, all-cell boundary exposure, structural profile values",
                "null":"All nonzero circular shifts of each own-islet vascular-minus-neural continuous profile; exact orientation null",
                "limitations":"Circular shifts model exchangeable orientation along arclength, not vascular causality; no single circular boundary for multipart islets; classifier process geometry differs from cell-body lineage"},
          "H2":{"estimand":"Equal-islet mean of conditional standardized beta-neighbor fragmentation and beta–immune weighted spatial adjacency excess; T1D minus ND positive, biological joint claim requires both components positive",
                "graph":"Own-islet symmetrized k-nearest centroid graph k=6; edges cannot leave own outline or cross holes; continuous length weight exp(-length / donor mean cell diameter)",
                "null":"199 seeded beta/alpha/delta label permutations within equal-count boundary-depth rank strata (>=10 endocrine per stratum), preserving subtype counts per stratum, all cell positions, immune identity and islet outline",
                "conditioning":"Beta abundance and native alpha/delta composition, within-islet radial architecture, fixed immune positions, outline geometry",
                "limitations":"Conditional on surviving mixed endocrine islets with core immune support; centroid adjacency is not segmented membrane contact; component null variance required"}},
          "analysis_population":"All 25 ordinary pancreatic protein donors with current all-exclusion assignments; only composition-eligible original Islet UUID outlines; finite coordinates and positive finite nucleus area required for every calculation",
          "identity":"Current production native calls derived from selected REDSEA/UniFORM values; refined calls are a sensitivity; no raw intensity retyping",
          "scale":"Mean equivalent eligible-cell diameter from fundamentals_cell_area_um2 in micrometres, no fixed new physical bands",
          "sparse_QC":{"H1":"at least 3 actual/effective immune cells and 1 immune-cell equivalent kernel mass in the continuous boundary kernel","H2":"at least 10 endocrine, 3 beta, 3 other endocrine, 1 own-core immune and variable component nulls"},
          "sensitivity_configurations":CONFIGS,"negative_control":"H1 neural process comparator and local all-cell exposure; nulls restricted to own islet, no pancreatic-to-node or between-islet cell links",
          "statistical_unit":"Donor; supported islets averaged equally within donor. Root common six-hypothesis donor testing and BH6; conditional per-islet ranks are QC, not independent disease tests",
          "acquisition_confounding":"Donor disease contrast requires root common-support acquisition sensitivities; within-islet arrangement null does not eliminate disease–acquisition aliasing",
          "missingness":"Unavailable endpoints saved with reasons, never zero-imputed; beta-depleted, immune-empty and constant-profile exclusions define conditional estimands",
          "seed":SEED,"source_assignment_manifest":str(ASS),"source_assignment_manifest_sha256":sha(ASS),"composition_source":str(COMP),"composition_sha256":sha(COMP),
          "version":"H1H2_spatial_v1_20261006"}
    plan["execution"]={"requested_parallel_workers":args.workers,"process_start_method":"spawn","manifest_order_aggregation":True,"code_sha256":sha(Path(__file__))}
    savejson(PROV/"spatial_H1_H2_plan.json",plan)
    records=json.loads(ASS.read_text())["donors"];comp=pd.read_parquet(COMP)
    if args.donors:records=[r for r in records if str(r["donor_id"]) in set(args.donors)]
    allrows=[];summaries=[];provenance=[];outcomes={};pending=[]
    for rec in records:
        cache=TAB/("spatial_H1_H2_islets_donor_"+str(rec["donor_id"])+".csv")
        ps=PROV/("spatial_H1_H2_sources_donor_"+str(rec["donor_id"])+".json")
        if args.resume and cache.exists() and ps.exists():
            p=json.loads(ps.read_text())
            assert p["assignment_sha256"]==rec["sha256"] and sha(rec["path"])==rec["sha256"]
            # Every measured/geometry source must still match; pure exact-search speedups preserve seeds and estimands.
            assert sha(p["cell_QC"]["native_source_path"])==p["cell_QC"]["native_source_sha256"]
            assert sha(p["islet_geometry"]["geometry_path"])==p["islet_geometry"]["geometry_sha256"]
            for structure in ["vascular","neural"]:
                if p[structure].get("path"):assert sha(p[structure]["path"])==p[structure]["sha256"]
            rows=pd.read_csv(cache,dtype={"donor_id":str,"islet_uuid":str}).to_dict("records")
            sums=[summarize([x for x in rows if x["configuration"]==config],rec,p["disease_status"],config,{"n_spatial_admitted":p["cell_QC"]["n_spatial_admitted"],"mean_cell_equivalent_diameter_um":p["cell_QC"]["mean_cell_equivalent_diameter_um"],"n_true_Islets":p["islet_geometry"]["valid_true_Islets"]}) for config,*_ in CONFIGS]
            outcomes[str(rec["donor_id"])]=(rows,sums,p)
            print("H1/H2 "+str(rec["donor_id"])+": validated cached sources and restored",flush=True)
        else:pending.append(rec)
    if args.workers>1 and pending:
        from concurrent.futures import ProcessPoolExecutor,as_completed
        import multiprocessing
        # Spawn avoids forking previously initialized Arrow/NumPy thread pools.
        with ProcessPoolExecutor(max_workers=min(args.workers,len(pending)),mp_context=multiprocessing.get_context("spawn")) as executor:
            futures={executor.submit(donor,rec,comp):str(rec["donor_id"]) for rec in pending}
            for future in as_completed(futures):
                d=futures[future];outcomes[d]=future.result()
    else:
        for rec in pending:outcomes[str(rec["donor_id"])]=donor(rec,comp)
    for rec in records:
        rows,sums,p=outcomes[str(rec["donor_id"])]
        allrows.extend(rows);summaries.extend(sums);provenance.append(p)
    s=pd.DataFrame(summaries);q=pd.DataFrame(allrows)
    s.to_csv(TAB/"spatial_H1_H2_sensitivities.csv",index=False)
    primary=s.loc[s.configuration.eq("native_primary")].copy();primary.to_csv(TAB/"spatial_H1_H2_donor_endpoints.csv",index=False)
    q.to_csv(TAB/"spatial_H1_H2_null_QC.csv",index=False)
    savejson(PROV/"spatial_H1_H2_source_manifest.json",{"plan":plan,"donors":provenance})
    plot(s)
    print("Saved H1/H2 donors="+str(len(primary))+" supported H1="+str(primary.H1_endpoint.notna().sum())+" H2="+str(primary.H2_endpoint.notna().sum()),flush=True)
if __name__=="__main__":main()
