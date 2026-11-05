"""Area support for serial-section features; never cell correspondence.

The CRO N-glycan report, sections 9.4-9.5 (pp. 19-20), motivates fractional
coverage, area precision/recall/Dice/IoU and purity-yield reporting. Here both
sections are independent observations: neither mask is biological ground truth.
All geometry must already share a declared calibrated micrometre frame.
No thresholds in this module confer production acceptance.
"""
from __future__ import annotations

from itertools import combinations
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import shapely
from shapely.affinity import translate


def _area_object(geometry, name: str):
    if geometry is None:
        return None
    if (geometry.geom_type not in {"Polygon", "MultiPolygon"}
            or geometry.is_empty or not geometry.is_valid
            or not np.isfinite(geometry.area) or geometry.area <= 0):
        raise ValueError(f"{name}: valid positive-area object required; no silent repair")
    return geometry


def overlap_metrics(source, target) -> dict:
    """Exact common area and directional support, including interior holes.

    source_coverage is precision if source is treated as the prediction;
    target_coverage is recall. Those names are algebraic, not truth claims.
    Missing geometry remains unavailable; disjoint observed objects give zeros.
    """
    source=_area_object(source, "source");target=_area_object(target, "target")
    fields=("source_area_um2", "target_area_um2", "intersection_um2", "union_um2",
            "source_only_um2", "target_only_um2", "source_coverage", "target_coverage",
            "dice", "iou", "reciprocal_coverage")
    if source is None or target is None:
        return {"status":"unavailable_geometry", **{k:np.nan for k in fields}}
    a,b=float(source.area),float(target.area)
    common=float(source.intersection(target).area)
    common=min(common,a,b)  # bound floating-point intersection roundoff
    union=a+b-common
    return dict(status="measured",source_area_um2=a,target_area_um2=b,
                intersection_um2=common,union_um2=union,source_only_um2=a-common,
                target_only_um2=b-common,source_coverage=common/a,target_coverage=common/b,
                dice=2*common/(a+b),iou=common/union,
                reciprocal_coverage=min(common/a,common/b))


def fractional_links(source: Mapping[str, object], target: Mapping[str, object], *,
                     source_frame: str, target_frame: str,
                     target_available: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Every positive-area link and unmatched/ambiguous source area.

    Source-area fractions are NOT renormalized over observed targets. Area in
    overlapping targets is explicitly ambiguous and removed from each edge's
    exclusive fraction. This conserves source area without choosing a target by
    order or silently duplicating a measurement. A split source can have several
    non-overlapping links; none is promoted to a unique biological identity.
    With missing target geometry, covered_fraction describes only known targets;
    exclusivity, ambiguity and unassigned support remain unavailable.
    """
    if not source_frame or source_frame != target_frame or not source_frame.endswith(":um"):
        raise ValueError("both object sets must share an explicit calibrated :um frame")
    if any(not isinstance(k,str) or not k for k in [*source,*target]):
        raise ValueError("non-empty stable object identifiers are required")
    for k,g in source.items(): _area_object(g,k)
    for k,g in target.items(): _area_object(g,k)
    target_missing=not target_available or any(g is None for g in target.values())
    keys=[k for k in sorted(target) if target[k] is not None] if target_available else []
    gs=[target[k] for k in keys]
    tree=shapely.STRtree(gs) if gs else None
    links=[];summaries=[]
    for sid in sorted(source):
        sg=source[sid]
        if sg is None or not target_available:
            summaries.append(dict(source_id=sid,status="unavailable_geometry",n_overlapping_targets=0,
                source_area_um2=np.nan if sg is None else sg.area,covered_fraction=np.nan,
                exclusive_fraction=np.nan,ambiguous_fraction=np.nan,unassigned_fraction=np.nan))
            continue
        intersections={}
        for i in sorted(tree.query(sg)) if tree is not None else []:
            hit=sg.intersection(gs[i])
            if hit.area>0:intersections[keys[i]]=hit
        union=shapely.union_all(list(intersections.values()))
        contested=shapely.union_all([a.intersection(b) for a,b in combinations(intersections.values(),2)])
        covered=float(union.area/sg.area);ambiguous=float(contested.area/sg.area)
        status="partial_target_geometry" if target_missing else "measured"
        summaries.append(dict(source_id=sid,status=status,n_overlapping_targets=len(intersections),
             source_area_um2=sg.area,covered_fraction=covered,
             exclusive_fraction=np.nan if target_missing else max(0.,covered-ambiguous),
             ambiguous_fraction=np.nan if target_missing else ambiguous,
             unassigned_fraction=np.nan if target_missing else max(0.,1-covered)))
        for tid,hit in intersections.items():
            links.append(dict(source_id=sid,target_id=tid,**overlap_metrics(sg,target[tid]),
                target_geometry_status=status,
                exclusive_source_fraction=np.nan if target_missing else hit.difference(contested).area/sg.area,
                source_ambiguous_fraction=np.nan if target_missing else ambiguous,relationship="serial_section"))
    link_columns=["source_id","target_id","status","source_area_um2","target_area_um2",
                  "intersection_um2","union_um2","source_only_um2","target_only_um2",
                  "source_coverage","target_coverage","dice","iou","reciprocal_coverage",
                  "target_geometry_status","exclusive_source_fraction","source_ambiguous_fraction","relationship"]
    summary_columns=["source_id","status","n_overlapping_targets","source_area_um2","covered_fraction",
                     "exclusive_fraction","ambiguous_fraction","unassigned_fraction"]
    return pd.DataFrame(links,columns=link_columns),pd.DataFrame(summaries,columns=summary_columns)


def coverage_fidelity(coverage, area_um2, thresholds: Sequence[float]=(0.3,0.5,0.7)) -> pd.DataFrame:
    """CRO area precision/recall/Dice/IoU for a discretized mask.

    A representation check, NOT independent registration validation. Only
    observed finite coverages enter denominators; their count/area is reported.
    Footprints must be non-overlapping for these aggregate area metrics.
    """
    c=np.asarray(coverage,dtype=float)
    a=np.broadcast_to(np.asarray(area_um2,dtype=float),c.shape)
    if c.ndim!=1 or np.any(np.isinf(c)) or np.any(np.isfinite(c)&((c<0)|(c>1))):
        raise ValueError("coverage must be a vector of fractions or NaN")
    if np.any(~np.isfinite(a)) or np.any(a<=0):raise ValueError("finite positive footprint areas required")
    good=np.isfinite(c);c=c[good];a=a[good]
    out=[]
    def div(x,y):return x/y if y>0 else np.nan
    for threshold in thresholds:
        if not np.isfinite(threshold) or not 0<threshold<=1:raise ValueError("threshold must be in (0,1]")
        sel=c>=threshold
        tp=float(np.sum(c[sel]*a[sel]));fp=float(np.sum((1-c[sel])*a[sel]));fn=float(np.sum(c[~sel]*a[~sel]))
        out.append(dict(threshold=threshold,n_available=int(good.sum()),n_unavailable=int((~good).sum()),
            n_selected=int(sel.sum()),available_area_um2=float(a.sum()),selected_area_um2=float(a[sel].sum()),
            true_positive_area_um2=tp if len(c) else np.nan,false_positive_area_um2=fp if len(c) else np.nan,
            false_negative_area_um2=fn if len(c) else np.nan,precision=div(tp,tp+fp),recall=div(tp,tp+fn),
            dice=div(2*tp,2*tp+fp+fn),iou=div(tp,tp+fp+fn),
            status="measured_representation_only" if len(c) else "unavailable"))
    return pd.DataFrame(out)


def footprint_coverages(footprints, objects: Sequence[object] | None) -> np.ndarray:
    """Exact fractional cover of each sampling footprint by an object union.

    Supports affine-transformed MALDI squares and virtual resolution audits.
    Objects unavailable (None or a missing member) yields NaN; an observed empty
    object set yields zero. Intersections are unioned to avoid double counting.
    """
    footprints=list(footprints)
    for g in footprints:
        if _area_object(g,"footprint") is None:raise ValueError("footprint cannot be missing")
    if objects is None or any(g is None for g in objects):return np.full(len(footprints),np.nan)
    for g in objects:_area_object(g,"reference")
    out=np.zeros(len(footprints),float)
    if not len(objects):return out
    tree=shapely.STRtree(objects)
    query=tree.query(footprints,predicate="intersects")
    if query.shape[1]:
        for i in np.unique(query[0]):
            g=footprints[i]
            hits=[g.intersection(objects[j]) for j in query[1,query[0]==i]]
            out[i]=min(1.,max(0.,shapely.union_all(hits).area/g.area))
    return out


def paired_feature_summary(source_values, target_values, weights) -> dict:
    """Descriptive regional comparison with common finite feature support.

    Weights are geometric support, not calibrated probabilities; measured zeros
    remain values. No source measurement is copied onto a target cell or object.
    Effective weight count describes concentration, not independent sample size.
    """
    a=np.asarray(source_values,float);b=np.asarray(target_values,float);w=np.asarray(weights,float)
    if a.ndim!=1 or a.shape!=b.shape or a.shape!=w.shape:raise ValueError("equal-length vectors required")
    if np.any(np.isinf(w)) or np.any(np.isfinite(w)&((w<0)|(w>1))):raise ValueError("support weights must be fractions or NaN")
    good=np.isfinite(a)&np.isfinite(b)&np.isfinite(w)&(w>0)
    a=a[good];b=b[good];w=w[good];total=float(w.sum())
    return dict(n_input=int(len(good)),n_supported=int(good.sum()),weight_sum=total,
        source_weighted_mean=float(np.sum(w*a)/total) if total else np.nan,
        target_weighted_mean=float(np.sum(w*b)/total) if total else np.nan,
        source_minus_target_weighted_mean=float(np.sum(w*(a-b))/total) if total else np.nan,
        effective_weight_count=float(total**2/np.sum(w*w)) if total else np.nan)


def shift_sensitivity(source, target, radii_um: Sequence[float]=(25.,50.,100.)) -> list[dict]:
    """Fixed shifts in eight directions; sensitivity bounds, not confidence limits.

    Pair identity stays fixed. No correction is selected by the resulting overlap.
    Use independently estimated displacement models for calibrated uncertainty.
    """
    baseline=overlap_metrics(source,target);out=[]
    for radius in radii_um:
        if not np.isfinite(radius) or radius<=0:raise ValueError("positive finite shift radius required")
        values=[baseline]
        if source is not None and target is not None:
            values += [overlap_metrics(translate(source,xoff=radius*np.cos(theta),yoff=radius*np.sin(theta)),target)
                       for theta in np.arange(8)*np.pi/4]
        finite=baseline['status']=='measured'
        out.append(dict(radius_um=radius,
            reciprocal_coverage_min=min(x['reciprocal_coverage'] for x in values) if finite else np.nan,
            reciprocal_coverage_max=max(x['reciprocal_coverage'] for x in values) if finite else np.nan,
            dice_min=min(x['dice'] for x in values) if finite else np.nan,
            dice_max=max(x['dice'] for x in values) if finite else np.nan,
            interpretation="fixed_pair_shift_sensitivity_not_probability"))
    return out
