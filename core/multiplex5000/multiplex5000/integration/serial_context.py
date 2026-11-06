"""Geometry-based candidate correspondence and focal-excluded neighbor evidence.

No projected intersection, biological label or expression enters these functions.
Neighbor scores are geometric diagnostics, not calibrated probabilities. Neither
assignment nor a diagnostic score accepts a production registration.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from shapely.geometry.polygon import orient


def geometry_descriptor(geometry):
    """Exact area moments plus size, solidity and circularity of a native object."""
    g=geometry
    if g is None or g.geom_type not in ('Polygon','MultiPolygon') or not g.is_valid or g.is_empty or g.area<=0:
        raise ValueError('Valid positive-area object required')
    center=np.array([g.centroid.x,g.centroid.y])
    moment=np.zeros((2,2),float)
    for part in g.geoms if g.geom_type=='MultiPolygon' else [g]:
        part=orient(part,sign=1.)
        for ring in [part.exterior,*part.interiors]:
            xy=np.asarray(ring.coords,dtype=float)-center
            x,y=xy[:-1].T;xx,yy=xy[1:].T
            cross=x*yy-xx*y
            moment[0,0]+=np.sum((x*x+x*xx+xx*xx)*cross)/12
            moment[1,1]+=np.sum((y*y+y*yy+yy*yy)*cross)/12
            moment[0,1]+=np.sum((2*x*y+x*yy+xx*y+2*xx*yy)*cross)/24
    moment[1,0]=moment[0,1]
    values=np.maximum(np.linalg.eigvalsh(moment/g.area),0)
    eccentricity=float(np.sqrt(max(0.,1-values[0]/values[-1]))) if values[-1]>0 else 0.
    return dict(x_um=float(center[0]),y_um=float(center[1]),area_um2=float(g.area),
        eccentricity=eccentricity,solidity=float(np.clip(g.area/g.convex_hull.area,0,1)),
        circularity=float(np.clip(4*np.pi*g.area/(g.length*g.length),0,1)))


def neighbor_vectors(xy,index,*,k=6,radius_um=1500.):
    """Return other objects relative to the focal center, with stable index ties."""
    xy=np.asarray(xy,float)
    if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():
        raise ValueError('Finite (n,2) positions required')
    if k<1 or not np.isfinite(radius_um) or radius_um<=0:raise ValueError('Positive finite neighborhood limits required')
    relative=xy-xy[index]
    distance=np.linalg.norm(relative,axis=1)
    ids=np.flatnonzero((np.arange(len(xy))!=index)&(distance<=radius_um))
    ids=ids[np.lexsort((ids,distance[ids]))][:k]
    return ids,relative[ids]


def compare_neighbor_vectors(source_vectors,target_vectors,*,tolerance_um=200.,minimum_neighbors=3):
    """Compare a focal-excluded arrangement using one-to-one partial matching.

    Each relative-vector agreement contributes max(0,1-error/tolerance).
    Missing or additional neighbors are allowed and reported. With too few
    neighbors the diagnostic is unavailable, not zero or an automatic failure.
    """
    a=np.asarray(source_vectors,float).reshape(-1,2)
    b=np.asarray(target_vectors,float).reshape(-1,2)
    if (not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(tolerance_um)
            or tolerance_um<=0 or minimum_neighbors<1):
        raise ValueError('Finite vectors and positive tolerance required')
    base=dict(n_source_neighbors=len(a),n_target_neighbors=len(b))
    if min(len(a),len(b))<minimum_neighbors:
        return dict(**base,status='insufficient_neighbors',score=np.nan,n_agreeing=0,
                    median_vector_error_um=np.nan,median_angle_error_deg=np.nan,links=[])
    distance=cdist(a,b)
    i,j=linear_sum_assignment(np.minimum(distance/tolerance_um,1.))
    good=distance[i,j]<tolerance_um
    ii,jj=i[good],j[good]
    agreement=np.maximum(0.,1-distance[i,j]/tolerance_um)
    angles=[]
    for u,v in zip(a[ii],b[jj]):
        norm=np.linalg.norm(u)*np.linalg.norm(v)
        if norm>0:angles.append(float(np.degrees(np.arccos(np.clip(np.dot(u,v)/norm,-1,1)))))
    return dict(**base,status='measured',score=float(2*agreement.sum()/(len(a)+len(b))),
        n_agreeing=int(good.sum()),median_vector_error_um=float(np.median(distance[ii,jj])) if good.any() else np.nan,
        median_angle_error_deg=float(np.median(angles)) if angles else np.nan,
        links=[dict(source_neighbor_index=int(u),target_neighbor_index=int(v),error_um=float(distance[u,v])) for u,v in zip(ii,jj)])


def geometry_costs(source,target,*,position_weight=.5,max_dist_um=200.,area_ratio=2.5):
    """Fixed weights for spatial distance and intrinsic morphology; no overlap."""
    if not 0<=position_weight<=1 or max_dist_um<=0 or area_ratio<=1:
        raise ValueError('Invalid position weight or candidate limits')
    cols=['x_um','y_um','area_um2','eccentricity','solidity','circularity']
    for table in [source,target]:
        if not np.isfinite(table[cols].to_numpy(float)).all() or (table.area_um2<=0).any():
            raise ValueError('Complete finite geometry descriptors required')
        if not table[['eccentricity','solidity','circularity']].apply(lambda s:s.between(0,1)).to_numpy().all():
            raise ValueError('Intrinsic shape descriptors must be fractions in [0,1]')
    distance=cdist(source[['x_um','y_um']],target[['x_um','y_um']])
    logarea=np.abs(np.log(source.area_um2.to_numpy()[:,None]/target.area_um2.to_numpy()[None,:]))
    delta=lambda name:np.abs(source[name].to_numpy()[:,None]-target[name].to_numpy()[None,:])
    intrinsic=(.70*logarea/np.log(area_ratio)+.20*delta('eccentricity')+
               .05*delta('solidity')+.05*delta('circularity'))
    feasible=(distance<=max_dist_um)&(logarea<=np.log(area_ratio))
    cost=position_weight*distance/max_dist_um+(1-position_weight)*intrinsic
    return cost,feasible,distance,intrinsic


def assign_geometry(source,target,*,position_weight=.5,max_dist_um=200.,area_ratio=2.5):
    """One-to-one candidate assignments with an unmatched option.

    Stable object IDs sort ties. An assignment is not an accepted correspondence.
    Costs are dimensionless and not calibrated probabilities. Unmatched sources
    cost 1; targets may remain unused. Feasible candidate costs are at most 1.
    """
    cols=['moving_id','fixed_id','cost','distance_um','intrinsic_cost','alternative_cost_margin','mutual_composite_minimum']
    if source.structure_id.duplicated().any() or target.structure_id.duplicated().any():
        raise ValueError('Unique native structure IDs required')
    source=source.sort_values('structure_id').reset_index(drop=True)
    target=target.sort_values('structure_id').reset_index(drop=True)
    if source.empty or target.empty:return pd.DataFrame(columns=cols)
    cost,feasible,distance,intrinsic=geometry_costs(source,target,position_weight=position_weight,
                                                max_dist_um=max_dist_um,area_ratio=area_ratio)
    masked=np.where(feasible,cost,1e6)
    rows,columns=linear_sum_assignment(np.concatenate([masked,np.ones((len(source),len(source)))],axis=1))
    out=[]
    for i,j in zip(rows,columns):
        if j>=len(target) or not feasible[i,j] or cost[i,j]>=1:continue
        row=masked[i].copy();row[j]=np.inf
        column=masked[:,j].copy();column[i]=np.inf
        alt=min(row.min(),column.min())
        out.append(dict(moving_id=source.structure_id.iloc[i],fixed_id=target.structure_id.iloc[j],
            cost=float(cost[i,j]),distance_um=float(distance[i,j]),intrinsic_cost=float(intrinsic[i,j]),
            alternative_cost_margin=float(alt-cost[i,j]) if alt<1e6 else np.nan,
            mutual_composite_minimum=bool(j==masked[i].argmin() and i==masked[:,j].argmin())))
    return pd.DataFrame(out,columns=cols)
