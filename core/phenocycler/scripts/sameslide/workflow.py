#!/usr/bin/env python3
"""Manifest-driven four-slide development execution with validated, immutable attempts."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import traceback
from datetime import datetime, timezone
from uuid import uuid4

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
DEFAULT_OUT = REPO / 'data/development/sameslide_four'
sys.path.insert(0, str(HERE.parent / 'xenium_celltype'))
import gdna
from manifest import load_manifest, get_section, verify_sex_metadata
import evidence


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def identity(path):
    p = Path(path); before = p.stat(); digest = sha(p); after = p.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'input changed while hashing: {p}')
    return {'path':str(p.resolve()),'bytes':after.st_size,'mtime_ns':after.st_mtime_ns,'sha256':digest}


def write_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.'+uuid4().hex+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str, allow_nan=False)+'\n')
    tmp.replace(path)


def section_root(root, section):
    return Path(root)/f"donor_id={section['donor_id']}"/f"roi={section['roi']}"/f"serial={section['xenium_serial']}"


def valid_receipt(path, fingerprint=None):
    try:
        r = json.loads(Path(path).read_text())
        if r['status'] != 'validated_complete' or (fingerprint and r['fingerprint'] != fingerprint):
            return False
        for item in r['outputs']:
            p = Path(item['path'])
            if not p.is_file() or p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
                return False
        return bool(r['outputs']) and bool(r['validation'])
    except (OSError, ValueError, KeyError):
        return False


def stage_inputs(section, sex_source):
    b = Path(section['bundle'])
    files = [b/n for n in ('gene_panel.json','metrics_summary.csv','experiment.xenium',
                           'cells.parquet','cell_feature_matrix.h5')]
    files += [HERE/'workflow.py', HERE/'evidence.py', HERE/'manifest.py', Path(gdna.__file__)]
    return {'section':{k:section[k] for k in ('section_id','donor_id','roi','tissue','sex','design','bundle','xenium_serial')},
            'sex_source':sex_source, 'files':[identity(p) for p in files],
            'parameters':{'alpha':evidence.ALPHA,'min_genes':20,'min_adjusted_counts':50,
                          'density_quantile':.98,'prior_weight':'genomic_probe_count'},
            'software':{n:importlib.metadata.version(n) for n in ('numpy','scipy','pandas','anndata','h5py','pyarrow')}}


def run_rna(section, dest):
    import numpy as np
    import pandas as pd
    import h5py
    import anndata as ad
    from scipy import sparse
    b = Path(section['bundle'])
    comp = gdna.panel_composition(b)
    cells = pd.read_parquet(b/'cells.parquet')
    cells['cell_id'] = cells.cell_id.astype(str)
    if cells.cell_id.duplicated().any():
        raise ValueError('duplicate vendor cell IDs')
    # Vendor verification uses every vendor cell, before any biological exclusion.
    raw_rate = cells.genomic_control_counts.to_numpy(float)/comp['genomic_control']
    vendor = gdna.verify_against_vendor(float(raw_rate.mean()), comp['gene'], gdna.vendor_metrics(b))
    # Cells without nuclei contribute only to the vendor reconciliation above.
    bio = (np.isfinite(cells.nucleus_area) & (cells.nucleus_area > 0) &
           np.isfinite(cells.cell_area) & (cells.cell_area > 0)).to_numpy()
    rates = np.zeros(len(cells))
    rates[bio], _, area_table, prior_weight = gdna.background_rates(
        cells.loc[bio,'genomic_control_counts'].to_numpy(), cells.loc[bio,'cell_area'].to_numpy(),
        comp['genomic_control'])
    with h5py.File(b/'cell_feature_matrix.h5','r') as h:
        m = h['matrix']; names = m['features']['name'][:].astype(str)
        types = m['features']['feature_type'][:].astype(str)
        ids = m['barcodes'][:].astype(str)
        full = sparse.csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])[::-1])
    if len(set(ids)) != len(ids) or set(ids) != set(cells.cell_id):
        raise ValueError('matrix and full vendor cell population differ')
    positions = pd.Index(cells.cell_id).get_indexer(ids)
    cells = cells.iloc[positions].reset_index(drop=True); rates = rates[positions]; bio = bio[positions]
    keep_genes = types == 'Gene Expression'; genes = names[keep_genes]
    if len(set(genes)) != len(genes) or len(genes) != comp['gene']:
        raise ValueError('matrix gene features and declared panel disagree')
    counts = full[:,keep_genes].tocsr()
    if not np.issubdtype(counts.dtype, np.integer):
        counts = evidence.validate_counts(counts).astype(np.int64)
    qc, qc_meta = evidence.qc_table(counts,cells,rates,comp['gene'])
    qc.loc[~bio,'gdna_lambda'] = np.nan
    qc.loc[~bio,['expected_background','adjusted_qc_total']] = np.nan
    qc['gdna_raw'] = raw_rate[positions]
    qc.to_parquet(dest/'cell_qc.parquet',index=False)
    qc[~qc.qc_eligible].to_parquet(dest/'exclusions.parquet',index=False)
    eligible = qc.qc_eligible.to_numpy(); x = counts[eligible]; lam = rates[eligible]
    if not len(lam):
        raise ValueError('no RNA QC-eligible cells')
    support = evidence.detected(x,lam)
    summary = evidence.group_summary(x,support,lam,genes)
    summary.to_parquet(dest/'gene_summary.parquet',index=False)
    controls = evidence.control_report(x,support,genes,qc[eligible].reset_index(drop=True),section['sex'],section['tissue'])
    controls.to_csv(dest/'sex_and_curated_controls.csv',index=False)
    # Technical probes never enter identity fitting. Genomic controls are evaluated
    # descriptively because they were used to estimate lambda, not held out.
    tech_rows = []
    for feature_type in sorted(set(types)-{'Gene Expression'}):
        ix = np.flatnonzero(types == feature_type)
        tx = full[eligible][:,ix].tocsr()
        ts = evidence.detected(tx,lam)
        table = evidence.group_summary(tx,ts,lam,names[ix])
        table['feature_type'] = feature_type
        table['evaluation_role'] = ('fitted_control_not_held_out' if 'Genomic' in feature_type else 'held_out_technical')
        tech_rows.append(table)
    pd.concat(tech_rows,ignore_index=True).to_parquet(dest/'technical_controls.parquet',index=False)
    del full
    obs = qc.loc[eligible].set_index('cell_id').copy(); obs.index.name = None
    for col,value in [('donor_id',section['donor_id']),('roi',section['roi']),('sample',section['xenium_serial']),('sex',section['sex'])]:
        obs[col] = pd.Categorical([value]*len(obs))
    depth = np.asarray(x.sum(axis=1)).ravel(); target = float(np.median(depth[depth>0]))
    lognorm = x.astype(np.float32).multiply((target/depth)[:,None]).tocsr(); np.log1p(lognorm.data,out=lognorm.data)
    var = pd.DataFrame(index=pd.Index(genes,dtype=object))
    var['sex_control'] = [g in evidence.SEX_GENES for g in genes]
    var['curated_control'] = [g in evidence.CLEAN_NEG for g in genes]
    var['identity_eligible'] = ~(var.sex_control | var.curated_control)
    var['null_reference_eligible'] = var.identity_eligible
    var['below_gdna_background'] = summary.below_section_background.to_numpy()
    a = ad.AnnData(X=lognorm,obs=obs,var=var)
    a.layers['counts'] = x; a.layers['detected_gdna'] = support
    a.obsm['spatial'] = obs[['x_centroid','y_centroid']].to_numpy()
    a.uns['prep'] = {'donor_id':section['donor_id'],'roi':section['roi'],'serial':section['xenium_serial'],
                     'bundle':str(b),'qc':qc_meta,'nominal_detection_alpha':evidence.ALPHA,
                     'detection_method':'Poisson.sf(k-1, cell_gdna_lambda)',
                     'normalization':'median_depth_log1p; raw integer counts preserved',
                     'uncertainty':'Poisson intervals conditional on estimated background; no spatial-independence guarantee'}
    a.write_h5ad(dest/'rna.h5ad',compression='lzf')
    # Read actual stored arrays, not just an HDF5 header.
    del a
    check = ad.read_h5ad(dest/'rna.h5ad')
    if not np.array_equal(check.obs_names,ids[eligible]) or not np.array_equal(check.var_names,genes):
        raise ValueError('H5AD axes changed on round trip')
    if not np.issubdtype(check.layers['counts'].dtype,np.integer) or (check.layers['counts'] != x).nnz:
        raise ValueError('raw counts changed on round trip')
    if (check.layers['detected_gdna'] != support).nnz:
        raise ValueError('detection evidence changed on round trip')
    failures = controls[(controls.status == 'fail') & (controls.stratum == 'section')]
    report = {'section_id':section['section_id'],'panel':comp,'vendor_check':vendor,
        'prior_weight':prior_weight,'area_prior':area_table.to_dict('records'),'qc':qc_meta,
        'nominal_alpha':evidence.ALPHA,'sex':section['sex'],
        'female_Y_negative_test':'not_applicable' if section['sex']=='M' else ('failed' if len(failures) else 'passed'),
        'female_Y_failures':failures.probe.tolist(),
        'claims_gate':'failed_controls' if len(failures) else 'controls_evaluated',
        'raw_counts_preserved':True,'eligible_cell_ids_sha256':hashlib.sha256('\n'.join(ids[eligible]).encode()).hexdigest()}
    write_json(dest/'normalization.json',report)
    return {'raw_counts_roundtrip':True,'axes_roundtrip':True,'detection_roundtrip':True,
            'n_cells':int(eligible.sum()),'n_genes':len(genes),'vendor_agreement':True}


def execute(manifest_path, section_id, root=DEFAULT_OUT):
    m = load_manifest(manifest_path); section = get_section(m,section_id)
    sex = verify_sex_metadata(m)
    inputs = stage_inputs(section,sex)
    fingerprint = hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage = section_root(root,section)/'rna'; stage.mkdir(parents=True,exist_ok=True)
    with (stage/'stage.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        current = stage/'current.json'
        if valid_receipt(current,fingerprint):
            print(f'{section_id}: validated RNA stage reused',flush=True); return
        dest = stage/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]); dest.mkdir(parents=True)
        rec = {'section_id':section_id,'fingerprint':fingerprint,'started_utc':now(),
               'inputs':inputs,'job_id':os.environ.get('SLURM_JOB_ID'),'account':os.environ.get('SLURM_JOB_ACCOUNT'),
               'status':'running','attempt':str(dest)}
        write_json(dest/'receipt.json',rec)
        try:
            validation = run_rna(section,dest)
            outputs = [identity(p) for p in sorted(dest.iterdir()) if p.name != 'receipt.json']
            rec.update(status='validated_complete',finished_utc=now(),outputs=outputs,validation=validation)
            write_json(dest/'receipt.json',rec); write_json(current,rec)
        except Exception as exc:
            rec.update(status='failed',finished_utc=now(),error=str(exc),traceback=traceback.format_exc())
            write_json(dest/'receipt.json',rec)
            write_json(stage/'last_failure.json',rec)
            raise
        print(f'{section_id}: RNA validated {validation}',flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a = p.parse_args(argv); execute(a.manifest,a.section,a.out)


if __name__ == '__main__':
    main()
