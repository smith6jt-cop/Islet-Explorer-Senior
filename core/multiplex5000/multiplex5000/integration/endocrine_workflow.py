"""Build Islets, EndocrineClusters and EndocrineSingles from native cell areas.

RNA evidence or protein classifier proposals feed the pinned live PubApp
whole-cell constructor. Fresh scratch projects preserve source images/qpdata.
Construction is candidate generation; image and correspondence review follow.
"""
from pathlib import Path
import argparse, concurrent.futures, copy, json, os, shutil, subprocess, sys, time
import ijson
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape,mapping,box
from shapely.affinity import scale

from .endocrine_objects import *
from .xenium_islets import read_section
from .serial_context import geometry_descriptor
import hashlib

def write_json(p,a): p.write_text(json.dumps(a,separators=(',',':'))+'\n')

def project(s,path,source_dir):
    p=path/'qupath';p.mkdir(exist_ok=True)
    lib=p/'scripts/lib';lib.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source_dir/'common.groovy',lib/'common.groovy')
    entry=copy.deepcopy(s['project_entry']);entry['entryID']=1
    # This fresh project has no source qpdata. The server supplies physical
    # calibration only; the constructor reads no image pixels.
    write_json(p/'project.qpproj',dict(version='0.7.0',createTimestamp=int(time.time()*1000),
        modifyTimestamp=int(time.time()*1000),lastID=1,images=[entry]))
    return p/'project.qpproj'

def protein(s,path):
    verify_source(s['cells_source']['path'],s['cells_source']['sha256'])
    verify_source(s['proposals_source']['path'],s['proposals_source']['sha256'])
    features=json.loads(Path(s['proposals_source']['path']).read_text())['features']
    features=[f for f in features if (f['properties'].get('classification') or {}).get('name')=='Islet']
    for f in features:f['properties']={'objectType':'annotation','classification':{'name':'Islet'}}
    props=np.array([shapely.make_valid(shape(f['geometry'])) for f in features],dtype=object)
    tree=shapely.STRtree(props)
    write_json(path/'proposals.geojson',dict(type='FeatureCollection',features=features))
    selected=[];geoms=[];native=[];total=0;nuclear=0
    with open(s['cells_source']['path'],'rb') as src:
        for f in ijson.items(src,'features.item',use_float=True):
            total+=1
            if not f.get('nucleusGeometry'):continue
            nuclear+=1
            g=shapely.make_valid(shape(f['geometry']))
            if not len(tree.query(g,predicate='intersects')):continue
            # Exclude no-nucleus objects completely, including construction context.
            selected.append(f);native.append(str(f['id']))
            geoms.append(scale(g,xfact=s['pixel_size_um'],yfact=s['pixel_size_um'],origin=(0,0)))
            if total%250000==0:print('protein',s['donor_id'],'read',total,'retained',len(selected),flush=True)
    write_json(path/'cells.geojson',dict(type='FeatureCollection',features=selected))
    pd.DataFrame(dict(cell_id=native,qupath_uuid=native,geometry_wkb=[g.wkb for g in geoms])).to_parquet(path/'native_cells.parquet',index=False)
    return dict(n_source_cells=total,n_nucleated_source_cells=nuclear,n_touching_nucleated_cells=len(selected),n_proposals=len(features))

def rna(s,path):
    cells=read_section(s)
    admitted=cells.nucleus_count.gt(0)&cells.depth.gt(0)
    cells.cell_id=cells.cell_id.astype(str)
    bp=s['vendor_boundaries'];verify_source(bp['path'],bp['sha256'])
    vertices=pd.read_parquet(bp['path'])
    available=set(vertices.cell_id.unique())
    missing=set(cells.loc[admitted,'cell_id'])-available
    if missing:raise ValueError(f'{len(missing)} admitted cells missing native segmentation')
    vertices=vertices[vertices.cell_id.isin(cells.loc[admitted,'cell_id'])]
    geos=vertex_polygons(vertices)
    proposals,seeds=rna_proposals(cells,geos,acquisition_id=s['acquisition_id'])
    write_json(path/'proposals.geojson',dict(type='FeatureCollection',features=proposal_features(proposals)))
    cellarr=np.array(geos.geometry.to_list(),dtype=object)
    _,ii=shapely.STRtree(cellarr).query(np.array(proposals.geometry.to_list(),dtype=object),predicate='intersects')
    chosen=geos.iloc[np.unique(ii)].copy()
    npath=Path(s['nucleus_boundaries']['path'])
    verify_source(npath,s['nucleus_boundaries']['sha256'])
    nv=pd.read_parquet(npath)
    nv=nv[nv.cell_id.isin(chosen.cell_id)]
    ng=vertex_polygons(nv).rename(columns={'geometry':'nucleus_geometry'})
    chosen=chosen.merge(ng,on='cell_id',how='left',validate='one_to_one')
    if chosen.nucleus_geometry.isna().any():raise ValueError('Measured nucleus polygon missing')
    chosen['qupath_uuid']=[native_uuid(s['acquisition_id'],v) for v in chosen.cell_id]
    selected=[dict(type='Feature',id=r.qupath_uuid,geometry=mapping(r.geometry),nucleusGeometry=mapping(r.nucleus_geometry),
                   properties=dict(objectType='cell')) for r in chosen.itertuples()]
    write_json(path/'cells.geojson',dict(type='FeatureCollection',features=selected))
    chosen.drop(columns=['geometry','nucleus_geometry']).assign(geometry_wkb=[g.wkb for g in chosen.geometry]).to_parquet(path/'native_cells.parquet',index=False)
    seeds.drop(columns='geometry').to_parquet(path/'RNA_seed_membership.parquet',index=False)
    proposals.drop(columns='geometry').assign(geometry_wkb=[g.wkb for g in proposals.geometry]).to_parquet(path/'RNA_proposals.parquet',index=False)
    sensitivity=[]
    for radius in [40,60]:
        labels=dbscan_labels(seeds[['x_um','y_um']].to_numpy(float),radius,1) if len(seeds) else []
        sensitivity.extend(dict(cell_id=cid,radius_um=radius,proposal_group=int(k)) for cid,k in zip(seeds.cell_id,labels))
    pd.DataFrame(sensitivity).to_parquet(path/'RNA_radius_membership_sensitivity.parquet',index=False)
    return dict(n_source_cells=len(cells),n_admitted_cells=int(admitted.sum()),n_RNA_seeds=len(seeds),
        n_touching_nucleated_cells=len(chosen),n_proposals=len(proposals),nucleus_boundary_source=str(npath),nucleus_boundary_sha256=digest(npath),
        **cells.attrs['input_accounting'])

def prepare(s,source_dir):
    start=time.time();path=Path(s['output']);path.mkdir(parents=True,exist_ok=True)
    if (path/'prepared.json').exists():return json.loads((path/'prepared.json').read_text())
    print('START',s['role'],s['donor_id'],flush=True)
    qc=protein(s,path) if s['role']=='protein' else rna(s,path)
    pp=project(s,path,source_dir)
    cfg=dict(pixel_size_um=s['pixel_size_um'],cells_geojson=str(path/'cells.geojson'),proposals_geojson=str(path/'proposals.geojson'),
        constructor_script=str(source_dir/'1b_islet_outlines_to_cells.groovy'),
        output_geojson=str(path/'constructed.geojson'),output_precise=str(path/'constructed_precise.json'),output_receipt=str(path/'constructor_receipt.json'))
    write_json(path/'constructor_config.json',cfg)
    result=dict(**s,**qc,project=str(pp),prepared_seconds=time.time()-start,
        prepared_hashes={x:digest(path/x) for x in ['cells.geojson','proposals.geojson','native_cells.parquet','constructor_config.json']})
    write_json(path/'prepared.json',result)
    print('PREPARED',s['role'],s['donor_id'],qc,'seconds',round(time.time()-start),flush=True)
    return result

def construct_prepared(s,source_dir,qupath):
    p=Path(s['output']);start=time.time()
    if (p/'verified.json').exists():return json.loads((p/'verified.json').read_text())
    prepared=json.loads((p/'prepared.json').read_text())
    project(s,p,source_dir) # fresh project JSON, never touches original qpdata
    cfg=json.loads((p/'constructor_config.json').read_text())
    for name,h in prepared['prepared_hashes'].items():verify_source(p/name,h)
    pins={str(source_dir/'1b_islet_outlines_to_cells.groovy'):'533331b7ec1d39a6fd35062daad46152163e78829c322b0f1266b8b523aa483a'}
    for x in [p/'cells.geojson',p/'proposals.geojson',p/'qupath/scripts/lib/common.groovy',
              Path(__file__).with_name('construct_whole_cells.groovy')]:
        pins[str(x)]=digest(x)
    print('CONSTRUCT',s['role'],s['donor_id'],flush=True)
    if (p/'execution.json').exists():
        execution=json.loads((p/'execution.json').read_text())
        assert execution['source_hashes']==pins
        for key,h in execution['outputs'].items():verify_source(cfg[key],h)
    else:
        execution=execute_constructor(config_path=p/'constructor_config.json',project_path=p/'qupath/project.qpproj',
            qupath=qupath,source_hashes=pins,log_path=p/'constructor.log')
        write_json(p/'execution.json',execution)
    precise=json.loads((p/'constructed_precise.json').read_text())
    features=[dict(type='Feature',id=r['id'],geometry=mapping(shapely.from_wkb(bytes.fromhex(r['geometry_wkb_hex']))),
        properties=dict(classification=dict(name=r['object_class']),measurements={'Whole cells':r['n_whole_cells']})) for r in precise]
    native=pd.read_parquet(p/'native_cells.parquet');native['geometry']=shapely.from_wkb(native.geometry_wkb)
    objects,links,qc=validate_constructed(features,native,pixel_size_um=s['pixel_size_um'])
    imported=pd.DataFrame(json.loads((p/'imported_cells_precise.json').read_text()))
    from shapely.affinity import scale
    imported['geometry']=[scale(shapely.from_wkb(bytes.fromhex(w)),xfact=s['pixel_size_um'],yfact=s['pixel_size_um'],origin=(0,0)) for w in imported.geometry_wkb_hex]
    imported=imported.merge(native[['qupath_uuid','cell_id']],on='qupath_uuid',validate='one_to_one')
    exact_objects,_,exact_qc=validate_constructed(features,imported,pixel_size_um=s['pixel_size_um'])
    qc['n_single_constructor_geometry_disagreements']=int((exact_objects.loc[exact_objects.object_class.eq('EndocrineSingle'),'single_relative_area_difference']>1e-10).sum())
    objects['single_constructor_relative_area_difference']=exact_objects.single_relative_area_difference
    qc['imported_cells_sha256']=digest(p/'imported_cells_precise.json')
    objects=objects.rename(columns={'structure_id':'constructor_object_id'})
    seedids=set(pd.read_parquet(p/'RNA_seed_membership.parquet').cell_id) if s['role']=='RNA' else set()
    if s['role']=='RNA':
        sc=links.loc[links.cell_id.isin(seedids)].groupby('object_id').size()
        objects['n_RNA_seed_cells']=objects.constructor_object_id.map(sc).fillna(0).astype(int)
        objects['proposal_evidence_available']=objects.n_RNA_seed_cells.gt(0)
    else:objects['proposal_evidence_available']=True
    members={k:sorted(z.cell_id.astype(str)) for k,z in links.groupby('object_id')}
    objects['structure_id']=[s['role']+':'+s['acquisition_id']+':'+hashlib.sha256(json.dumps(members.get(i,[])).encode()).hexdigest()[:24]
                            if members.get(i) else s['role']+':unavailable:'+i for i in objects.constructor_object_id]
    if objects.structure_id.duplicated().any():raise ValueError('Construction member sets collide')
    imap=dict(zip(objects.constructor_object_id,objects.structure_id));links['structure_id']=links.object_id.map(imap)
    for field in ['donor_id','source_serial','acquisition_id','role']:objects[field]=s[field];links[field]=s[field]
    objects['coordinate_frame']=('xenium:' if s['role']=='RNA' else '')+s['acquisition_id']+':um'
    objects['boundary_source']='pubapp_whole_cell_constructor'
    objects['proposal_source']='RNA_DBSCAN1_cell_hull' if s['role']=='RNA' else 'native_image_classifier'
    objects['individual_correspondence_validated']=False
    desc=pd.DataFrame([geometry_descriptor(shapely.from_wkb(w)) for w in objects.geometry_wkb])
    for col in desc:objects[col]=desc[col].to_numpy()
    objects.to_parquet(p/'objects.parquet',index=False);links.to_parquet(p/'object_members.parquet',index=False)
    raw_proposals=json.loads((p/'proposals.geojson').read_text())['features']
    prop=pd.DataFrame(dict(proposal_id=[str(f['id']) for f in raw_proposals],geometry=[scale(area_only(shapely.make_valid(shape(f['geometry']))),xfact=s['pixel_size_um'],yfact=s['pixel_size_um'],origin=(0,0)) for f in raw_proposals]))
    lineage=proposal_lineage(objects.assign(geometry=shapely.from_wkb(objects.geometry_wkb)),prop)
    lineage.to_parquet(p/'object_proposal_lineage.parquet',index=False)
    qc['n_objects_without_positive_area_proposal_lineage']=len(set(objects.structure_id)-set(lineage.structure_id))

    qc.update(donor_id=s['donor_id'],source_serial=s['source_serial'],role=s['role'],seconds=time.time()-start,
        class_counts=objects.object_class.value_counts().to_dict(),
        n_without_proposal_evidence=int((~objects.proposal_evidence_available).sum()),
        n_RNA_seeds_unassigned=len(seedids-set(links.cell_id)) if s['role']=='RNA' else None,
        source_hashes=pins,output_sha256=digest(p/'objects.parquet'),membership_sha256=digest(p/'object_members.parquet'))
    write_json(p/'verified.json',qc)
    print('VERIFIED',s['role'],s['donor_id'],qc['class_counts'],'QC', {k:v for k,v in qc.items() if k.startswith('n_')},flush=True)
    return qc


def run_manifest(manifest_path,output,*,qupath):
    manifest_path=Path(manifest_path);output=Path(output)
    m=json.loads(manifest_path.read_text())
    if m.get('schema')!='sequential_whole_endocrine_v1':raise ValueError('Unsupported whole-endocrine manifest schema')
    sections=m['sections'];keys=[(s['role'],s['acquisition_id']) for s in sections]
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate role/acquisition keys')
    # This entry point is explicitly bound to the scientifically verified live
    # PubApp constructor; a different upstream version needs a reviewed update.
    if m['constructor']['sha256']!='533331b7ec1d39a6fd35062daad46152163e78829c322b0f1266b8b523aa483a':
        raise ValueError('Unreviewed PubApp constructor generation')
    for key in ['constructor','common']:
        verify_source(m[key]['path'],m[key]['sha256'])
    for section in sections:
        if section['role'] not in ['RNA','protein']:raise ValueError('Unknown modality role')
        if section['role']=='RNA':
            if section['pixel_size_um']!=1:raise ValueError('RNA vendor geometry uses native micrometres')
            for key in ['vendor_boundaries','nucleus_boundaries']:verify_source(section[key]['path'],section[key]['sha256'])
    output.mkdir(parents=True,exist_ok=False)
    write_json(output/'input_manifest.json',m)
    source=output/'source_snapshot';source.mkdir()
    for key,name in [('constructor','1b_islet_outlines_to_cells.groovy'),('common','common.groovy')]:
        shutil.copy2(m[key]['path'],source/name)
    results=[]
    for section in sections:
        section={**section,'output':str(output/section['role']/section['source_serial'])}
        prepare(section,source);results.append(construct_prepared(section,source,qupath))
    write_json(output/'run_receipt.json',dict(status='constructed_candidates',schema=m['schema'],
        input_manifest_sha256=digest(manifest_path),module_sha256=digest(__file__),
        live_workflow_binding=m.get('live_workflow_binding'),results=results,
        individual_pairs_validated=False,production_accepted=False,
        next_required='Actual-image geometry QA followed by anatomy/geometry correspondence review'))
    return results


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--qupath',required=True,type=Path)
    a=ap.parse_args(argv);run_manifest(a.manifest,a.output,qupath=a.qupath)


if __name__=='__main__':main()
