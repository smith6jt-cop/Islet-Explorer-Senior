#!/usr/bin/env python3
"""Isolated paired selected-Cell control versus policy-own-compartment candidate.
No baseline call, registry, expression, shift, scale, or production path is overwritten.
All altered compartments get newly estimated UniFORM shifts and feature scales.
This is a review arm; it does not establish cell-type accuracy or select a new release.
"""
import os,sys
os.environ.setdefault("OMP_NUM_THREADS","6")
os.environ.setdefault("OPENBLAS_NUM_THREADS","6")
os.environ.setdefault("NUMBA_NUM_THREADS","6")
from pathlib import Path
import argparse,json,hashlib,subprocess,time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr
ROOT=Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
R=ROOT/"results/v3/islet_immune_reassessment_20261006"
DATA=ROOT/"core/phenocycler/data"
TYPING=ROOT/"core/phenocycler/scripts/typing"
CELLTYPE=ROOT/"core/phenocycler/scripts/celltype"
POLICY=ROOT/"results/v3/release_20261003/calibration_review/marker_redsea_ruling.csv"
REGISTRY=ROOT/"core/phenocycler/phenocycler/marker_registry.json"
GEOMETRY=ROOT/"data/development/typing"
BLOCK=Path("/data/io60-calibration/alignment_block_current")
ASSIGN=Path("/data/io60-calibration/pancreas_analysis_assignments_single_exclusions_20261005")
RELEASE=Path("/data/io60-calibration/pancreas_release_production")
ND=["0112","115","6374","6442","6461","6476","6482","6493","6539"]
PILOT=["0112","115","6442","6493"]
def canonical(d):return str(int(str(d)))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def writej(p,x):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 Path(p).write_text(json.dumps(x,indent=2,default=str,allow_nan=False)+"\n")
def read(p,columns=None):return pq.ParquetFile(p).read(columns=columns).to_pandas()
def pathone(root,d):
 files=sorted((Path(root)/f"donor_id={d}").glob("*.parquet"))
 if len(files)!=1:raise ValueError(f"{root} donor {d}: {len(files)} files, expected one")
 return files[0]
def source_audit():
 policy=pd.read_csv(POLICY).set_index("marker").measurement_compartment.to_dict()
 actual=[]
 for a in sorted((DATA/"audit/donors/expression_availability").glob("donor_id=*/data.parquet")):
  q=read(a);q["source_audit_path"]=str(a);q["source_audit_sha256"]=sha(a)
  q["ruling_compartment"]=q.marker.map(policy).fillna("not specified")
  q["policy_source_mismatch"]=(q.ruling_compartment.isin(["Nucleus","Membrane","Cytoplasm"]) & q.compartment.ne(q.ruling_compartment))
  actual.append(q)
 a=pd.concat(actual,ignore_index=True)
 a.to_csv(R/"tables/protein_actual_compartment_source_QC.csv",index=False)
 sample=read(R/"tables/protein_ND_allcell_nonhormone_sample_coordinates.parquet",["object_id","donor_id"])
 shifts=json.loads((ROOT/"data/development/celltype/shifts.json").read_text())["shifts"]
 proofs=[]
 markers=["CD3e","CD68","INS","Iba1","Pan_Cytokeratin"]
 for d in ND:
  ids=set(sample.loc[sample.donor_id.astype(str).eq(canonical(d)),"object_id"].iloc[:1000].astype(str))
  expr=read(pathone(DATA/"expression",d),["object_id"]+markers);expr=expr[expr.object_id.astype(str).isin(ids)]
  rr=read(pathone(DATA/"redsea",d),["object_id"]+[c+"__"+m for c in ["Cell","Membrane","Cytoplasm"] for m in markers])
  rr=rr[rr.object_id.astype(str).isin(ids)];q=expr.merge(rr,on="object_id",validate="one_to_one")
  if len(q)!=len(ids):raise AssertionError(f"{d}: sampled UUID universe changed")
  uf=read(ROOT/f"data/development/celltype/uniform/donor_id={d}/data.parquet",["object_id"]+markers)
  uf=uf[uf.object_id.astype(str).isin(ids)].rename(columns={m:"uniform_"+m for m in markers})
  rel=read(RELEASE/f"donor_id={d}/data.parquet",["object_id"]+["redsea_aligned_"+m for m in markers])
  rel=rel[rel.object_id.astype(str).isin(ids)]
  q=q.merge(uf,on="object_id",validate="one_to_one").merge(rel,on="object_id",validate="one_to_one")
  for m in markers:
   for c in ["Cell","Membrane","Cytoplasm"]:
    v=q[m].to_numpy(float);w=q[c+"__"+m].to_numpy(float);ok=np.isfinite(v)&np.isfinite(w)
    proofs.append(dict(donor_id=canonical(d),marker=m,comparison="current expression vs REDSEA "+c,n_UUIDs=len(q),n_joint_finite=int(ok.sum()),
     comparison_missing=int((~np.isfinite(w)).sum()),fraction_equal=float(np.isclose(v[ok],w[ok],rtol=1e-7,atol=1e-8).mean()) if ok.any() else np.nan,
     median_abs_difference=float(np.median(np.abs(v[ok]-w[ok]))) if ok.any() else np.nan))
   v=q["uniform_"+m].to_numpy(float);w=q[m].to_numpy(float)*float(shifts[m]["scale"][d]);ok=np.isfinite(v)&np.isfinite(w)
   proofs.append(dict(donor_id=canonical(d),marker=m,comparison="current UniFORM vs expression times recorded donor multiplier",n_UUIDs=len(q),n_joint_finite=int(ok.sum()),
    fraction_equal=float(np.isclose(v[ok],w[ok],rtol=1e-7,atol=1e-8).mean()),median_abs_difference=float(np.median(np.abs(v[ok]-w[ok])))))
   w=q["redsea_aligned_"+m].to_numpy(float);ok=np.isfinite(v)&np.isfinite(w)
   proofs.append(dict(donor_id=canonical(d),marker=m,comparison="selected release aligned vs current UniFORM",n_UUIDs=len(q),n_joint_finite=int(ok.sum()),
    fraction_equal=float(np.isclose(v[ok],w[ok],rtol=1e-7,atol=1e-8).mean()),median_abs_difference=float(np.median(np.abs(v[ok]-w[ok])))))
 p=pd.DataFrame(proofs);p.to_csv(R/"tables/protein_actual_compartment_numeric_parity_QC.csv",index=False)
 record=json.loads(Path("/data/io60-calibration/arms/redsea_shifted_frozen/run.json").read_text())
 checked=[]
 for d in ND:
  f=ROOT/f"data/development/celltype/uniform/donor_id={d}/data.parquet";h=sha(f)
  checked.append(dict(donor_id=canonical(d),path=str(f),actual_sha256=h,recorded_sha256=record["aligned_files"][d],matches=h==record["aligned_files"][d]))
 pd.DataFrame(checked).to_csv(R/"tables/protein_selected_aligned_hash_QC.csv",index=False)
 writej(R/"provenance/protein_actual_compartment_source_QC.json",dict(status="complete",actual_registry_compartments=pd.Series([m["compartment"] for m in json.loads(REGISTRY.read_text())["markers"]]).value_counts().to_dict(),
  marker_registry_path=str(REGISTRY),marker_registry_sha256=sha(REGISTRY),policy_path=str(POLICY),policy_sha256=sha(POLICY),
  conclusion="Selected expression reads REDSEA Cell for nonnuclear channels, despite own-compartment ruling; original nuclear measurements remain uncompensated Nucleus.",
  sampled_UUIDs_per_ND_donor=1000,sample_source="immutable admitted donor-balanced ND audit UUIDs",baseline_unchanged=True))
 print("actual-source audit saved",flush=True)
def prepare(out,donors,methods=None):
 from copy import deepcopy
 sys.path.insert(0,str(ROOT/"core/phenocycler"))
 from phenocycler.marker_registry import registry_from_dict
 from phenocycler.expression import build_selected_expression,write_donor_partition,CONTEXT_COLUMNS
 payload=json.loads(REGISTRY.read_text());policy=pd.read_csv(POLICY).set_index("marker").measurement_compartment.to_dict()
 own=deepcopy(payload)
 for m in own["markers"]:
  c=policy.get(m["name"])
  if c in ["Membrane","Cytoplasm","Nucleus"]:m["compartment"]=c
 regs={"Cell_control":registry_from_dict(payload),"own_compartment":registry_from_dict(own)}
 if methods is not None:
  if not methods or set(methods)-set(regs):raise ValueError("Unknown or empty requested measurement methods")
  regs={name:regs[name] for name in methods}
 out.mkdir(parents=True,exist_ok=True)
 writej(out/"registry_baseline_copy.json",payload);writej(out/"registry_own_compartment_candidate.json",own)
 led=[];measurement_rows=[]
 for d in donors:
  canonical_d=canonical(d)
  paths={s:pathone(DATA/s,d) for s in ["cells","redsea","geometry_qc","expression"]}
  ap=ASSIGN/f"donor_id={canonical_d}/data.parquet"
  assign=read(ap,["object_id","analysis_eligible","analysis_exclusion_reason"])
  bad=set(assign.loc[~assign.analysis_eligible.astype(bool),"object_id"].astype(str))
  # Preserve nuclear/source QC, add every mandatory current pancreas exclusion.
  qc=read(paths["geometry_qc"]);old=qc.analysis_eligible.astype(bool).to_numpy()
  mandatory=qc.object_id.astype(str).isin(bad).to_numpy()
  qc.loc[mandatory,"analysis_eligible"]=False;qc.loc[mandatory,"estimation_eligible"]=False
  availpath=DATA/f"audit/donors/expression_availability/donor_id={d}/data.parquet"
  availability=read(availpath);acquired=set(availability.loc[availability.available.astype(bool),"marker"])
  raw_cols=list(dict.fromkeys(["object_id"]+[c for c in CONTEXT_COLUMNS if c in pq.ParquetFile(paths["cells"]).schema.names]
                +[s.measurement_column for reg in regs.values() for s in reg.markers if s.measurement_source=="qupath_raw" and s.name in acquired]))
  cells=read(paths["cells"],raw_cols)
  corrected_cols=list(dict.fromkeys(["object_id"]+[s.measurement_column for reg in regs.values() for s in reg.markers if s.measurement_source=="redsea" and s.name in acquired]))
  corrected=read(paths["redsea"],corrected_cols)
  geo=read(GEOMETRY/f"revised_cells/donor_id={d}/data.parquet")
  geo.loc[geo.object_id.astype(str).isin(bad),"qc_analysis_eligible"]=False
  geo.loc[geo.object_id.astype(str).isin(bad),"qc_estimation_eligible"]=False
  geo.loc[geo.object_id.astype(str).isin(bad),"analysis_reason"]="immutable_current_analysis_exclusion"
  gp=out/f"geometry/revised_cells/donor_id={d}/data.parquet";gp.parent.mkdir(parents=True,exist_ok=True)
  if gp.exists():raise FileExistsError(gp)
  geo.to_parquet(gp,index=False)
  for arm,reg in regs.items():
   frame,av=build_selected_expression(cells,corrected,qc,donor_id=d,acquired_markers=acquired,registry=reg)
   # Do not preserve older QC-estimation-ineligible cells through a renamed eligible field.
   target=write_donor_partition(frame,out/arm/"expression",donor_id=d)
   avtarget=write_donor_partition(av,out/arm/"availability",donor_id=d)
   for m in av.itertuples():
    vals=pd.to_numeric(frame[m.marker],errors="coerce").to_numpy(float);eligible=frame.qc_analysis_eligible.to_numpy(bool)
    measurement_rows.append(dict(arm=arm,donor_id=canonical_d,marker=m.marker,source=m.source,compartment=m.compartment,
     acquired=bool(m.available),n_eligible=int(eligible.sum()),n_finite_eligible=int(np.sum(eligible&np.isfinite(vals))),
     n_unavailable_eligible=int(np.sum(eligible&~np.isfinite(vals))),n_zero_eligible=int(np.sum(eligible&np.isfinite(vals)&(vals==0)))))
   led.append(dict(donor_id=canonical_d,donor_id_native=d,arm=arm,rows=len(frame),n_baseline_analysis_eligible=int(old.sum()),
    n_current_mandatory_removed_from_source_QC=int((old&mandatory).sum()),n_new_eligible=int(frame.qc_analysis_eligible.sum()),
    expression=str(target),expression_sha256=sha(target),availability=str(avtarget),availability_sha256=sha(avtarget),
    source_files={k:dict(path=str(p),sha256=sha(p)) for k,p in paths.items()},assignment=dict(path=str(ap),sha256=sha(ap))))
   print("expression",arm,d,len(frame),flush=True)
  del corrected,cells,geo,qc
 # Exact current fragment/object geometry is reused only as geometry, never processed marker values.
 for name in ["crossed","qupath_fundamentals"]:
  target=out/"geometry"/name
  if not target.exists():target.symlink_to(GEOMETRY/name,target_is_directory=True)
 pd.DataFrame(measurement_rows).to_csv(out/"expression_measurement_QC.csv",index=False)
 writej(out/"source_manifest.json",dict(status="expression_complete",donors=donors,arms=list(regs),records=led,policy=str(POLICY),policy_sha256=sha(POLICY),
  baseline_registry_sha256=sha(REGISTRY),ruling_scope="33 Membrane,18 Cytoplasm,7 Nucleus; unmapped state channels explicitly retain baseline source",
  current_block_sidecars=str(BLOCK),current_block_manifest_sha256=sha(BLOCK/"manifest.json"),baseline_unchanged=True))
def uniform_worker(out,arm,donors,references_from=None):
 import panel,alignment_block as AB,run_step2 as S2,usability as U
 import matplotlib.pyplot as plt
 AB.OUT=BLOCK
 # Recompute the acquisition-group statistic on the new marker's own finite block population.
 # A missing cytoplasmic mean never changes another marker's admitted alignment population.
 rows=[];cohort=panel.cohort().set_index("donor_id");group=cohort.acq_group
 markers=panel.markers(include_variant=True)
 for d in donors:
  f=read(panel.expression_path(d));eligible=f.qc_analysis_eligible.to_numpy(bool)&AB.population(d,f.object_id,f.cell_region)
  for m in markers:
   if m not in set(panel.donor_markers(d)) or m not in f:continue
   v=pd.to_numeric(f.loc[eligible,m],errors="coerce").to_numpy(float);finite=v[np.isfinite(v)]
   rows.append(dict(donor_id=d,marker=m,n_eligible=int(eligible.sum()),n_finite=len(finite),p50=float(np.median(finite)) if len(finite) else np.nan))
 q=pd.DataFrame(rows);ratios=[]
 for m,g in q.groupby("marker"):
  ratios.append(dict(marker=m,**U.acquisition_group_ratio(g.set_index("donor_id").p50,group)))
 q=q.merge(pd.DataFrame(ratios),on="marker",how="left")
 q.to_parquet(out/arm/"marker_donor_usability.parquet",index=False)
 per=q.groupby("marker").agg(acq_ratio=("acq_ratio","first")).reset_index()
 chosen=[m for m in markers if any(m in set(panel.donor_markers(d)) for d in donors)]
 # A corrected source keeps reference donor identities, never old Cell multipliers
 # or floors. A diagnostic without a declared reference record chooses anew.
 references={}
 if references_from is not None:
  previous=json.loads(Path(references_from).read_text())
  if set(previous["donors"])!=set(donors):raise ValueError("Pinned reference source must cover the exact processing cohort")
  for m in chosen:
   if m not in previous["shifts"]:raise ValueError(f"Original reference missing for {m}")
   ref=previous["shifts"][m]["reference_donor"]
   vals=S2.marker_values(m,[ref],{}) if ref in donors else {}
   if ref not in vals:raise ValueError(f"Original reference donor {ref} unavailable for corrected {m}; no silent replacement")
   references[m]=ref
 shifts=S2.derive(chosen,donors,per,save_figs=True,figs=out/arm/"figures_step2",references=references)
 if references:
  for m,ref in references.items():
   if shifts[m]["reference_donor"]!=ref:raise ValueError(f"Pinned reference changed for {m}: {ref} -> {shifts[m]['reference_donor']}")
 def clean(v):
  if isinstance(v,(np.integer,)):return int(v)
  if isinstance(v,(np.floating,float)):return float(v) if np.isfinite(v) else "unavailable"
  if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
  if isinstance(v,list):return [clean(x) for x in v]
  return v
 writej(out/arm/"shifts.json",dict(uniform_commit=S2.UA.UNIFORM_COMMIT,num_bins=S2.UA.NUM_BINS,bar_mode_sd=S2.BAR_MODE_SD,
  donors=donors,markers=chosen,references_from=str(references_from) if references_from is not None else None,
  reference_source_sha256=sha(references_from) if references_from is not None else None,
  pinned_reference_donor_ids=references,expression_inputs={d:str(panel.expression_path(d)) for d in donors},
  shifts=clean(shifts),own_measurement_missingness="per-marker finite block population; no complete-case deletion across unrelated markers",
  usability_scope="Recomputed acquisition-group ratios; no inherited Cell acquisition ratio, shift, floor or span."))
 S2.apply_shift(shifts,chosen,donors,out=out/arm/"uniform")
 # A finite raw measurement without a measurable alignment histogram has no
 # estimated aligned unit. Preserve that evidence in expression, but never let
 # the raw value stand in a newly derived aligned-value column.
 unavailable_aligned=[]
 for d in donors:
  missing=[m for m in chosen if d in shifts[m].get("unavailable_donors",{})]
  if not missing:continue
  target=out/arm/f"uniform/donor_id={d}/data.parquet"
  frame=read(target)
  for m in missing:
   if m not in frame:continue
   unavailable_aligned.append(dict(donor_id=d,marker=m,n_raw_finite_inside_alignment_block=int(np.isfinite(frame[m].to_numpy(float)).sum()),
    reason=shifts[m]["unavailable_donors"][d],aligned_values_written_unavailable=True))
   frame[m]=np.nan
  frame.to_parquet(target,index=False)
 writej(out/arm/"unavailable_alignment_histograms.json",dict(status="complete",cases=unavailable_aligned,
  semantics="Raw expression preserved; an unavailable alignment estimate produces unavailable aligned values, never raw substitution"))
 import run_step3 as S3
 S3.ALIGNED=out/arm/"uniform";AB.OUT=BLOCK
 scale=S3.build_scale(donors,markers,n_per_donor=200000,seed=0,ratio=False)
 scale.sort_index().reset_index().to_parquet(out/arm/"scale.parquet",index=False)
 writej(out/arm/"scale_provenance.json",dict(producer="run_step3.build_scale",donors=donors,n_per_donor=200000,seed=0,nuclear_ratio=False,
  shifts_sha256=sha(out/arm/"shifts.json"),scale_sha256=sha(out/arm/"scale.parquet"),block_sidecars=str(BLOCK)))
 print("alignment and scale finished",arm,flush=True)
def subprocess_stage(out,arm,stage,cmd,env):
 p=out/arm/f"{stage}.log"
 with p.open("w") as log:
  r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,text=True)
 print(arm,stage,"returncode",r.returncode,flush=True)
 if r.returncode:raise RuntimeError(f"{arm} {stage} failed; see {p}")
 return dict(name=stage,command=cmd,log=str(p),returncode=r.returncode)
def typing(out,arm,donors,env):
 path=out/arm/"typing";path.mkdir(parents=True,exist_ok=True)
 da=["--donors",*donors];sources={f:sha(TYPING/f) for f in ["subpanel.py","build_cells.py","refine.py","immune_subtype.py"]}
 rec=dict(status="running",donors=donors,aligned_root=str(out/arm/"uniform"),scale=str(out/arm/"scale.parquet"),
  geometry_root=str(out/"geometry"),registry_criteria=str(GEOMETRY/"subpanel/subpanel_log.json"),
  registry_criteria_sha256=sha(GEOMETRY/"subpanel/subpanel_log.json"),source_scripts=sources,stages=[],baseline_unchanged=True,
  current_scripts_note="Current producer preserves Single exclusions; current script hashes are recorded independently of Oct3 baseline hashes.",
  representation_note="Refinement uses unshifted CD8 recovered from this method's matching shifts, converted CD8 feature units and a fresh representation; primary calls retain shifted CD8.")
 writej(path/"run.json",rec)
 stages=[
 ("subpanel",[sys.executable,str(TYPING/"subpanel.py"),*da,"--membership","geometry","--emit-detail","--refit","--work",str(path/"scratch/subpanel_context"),"--audit-out",str(path/"subpanel_audit.json"),"--out",str(path/"subpanel_rebuilt")]),
 ("build_cells_primary",[sys.executable,str(TYPING/"build_cells.py"),*da,"--calls",str(path/"subpanel_rebuilt"),"--out",str(path/"cells")]),
 ("refine",[sys.executable,str(Path(__file__).resolve()),"--worker","refine","--out",str(out),"--arm",arm,*da]),
 ("build_cells_with_refinement",[sys.executable,str(TYPING/"build_cells.py"),*da,"--calls",str(path/"subpanel_rebuilt"),"--refine",str(path/"refinement"),"--audit",str(path/"subpanel_audit.json"),"--out",str(path/"cells")]),
 ("immune_subtype",[sys.executable,str(TYPING/"immune_subtype.py"),*da])]
 for stage,cmd in stages:
  try:rec["stages"].append(subprocess_stage(out,arm,stage,cmd,env))
  except Exception:
   rec["status"]="failed";rec["failed_stage"]=stage;writej(path/"run.json",rec);raise
  writej(path/"run.json",rec)
 rec["status"]="complete";writej(path/"run.json",rec)
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument("--out",type=Path,default=R/"provenance/protein_own_compartment_candidate_pilot")
 ap.add_argument("--donors",nargs="+",default=PILOT)
 ap.add_argument("--audit-only",action="store_true")
 ap.add_argument("--prepare-only",action="store_true")
 ap.add_argument("--skip-audit",action="store_true")
 ap.add_argument("--skip-prepare",action="store_true")
 ap.add_argument("--worker",choices=["uniform","refine"])
 ap.add_argument("--arm",choices=["Cell_control","own_compartment"])
 ap.add_argument("--references-from",type=Path)
 ap.add_argument("--current-cohort",action="store_true",help="Allow only the exact eligible current QuPath manifest cohort for approved source repair")
 args=ap.parse_args();out=args.out.resolve();donors=list(dict.fromkeys("0112" if canonical(d)=="112" else canonical(d) for d in args.donors))
 if args.current_cohort:
  sys.path.insert(0,str(ROOT/"core/phenocycler"))
  from phenocycler.cohort import DONOR_EXCLUSIONS
  manifest=json.loads((DATA/"qupath_manifest.json").read_text())
  permitted={str(r["donor_id"]) for r in manifest["images"] if canonical(r["donor_id"]) not in {canonical(d) for d in DONOR_EXCLUSIONS}}
  if set(donors)!=permitted or len(donors)!=len(permitted):raise ValueError("Approved full-cohort source repair must use the exact current eligible QuPath manifest cohort")
 elif set(donors)-set(ND):raise ValueError("This candidate audit is restricted to the nine clinically ND donors")
 if ROOT==out or DATA in out.parents:raise ValueError("Candidate output must be isolated from stable production data")
 if args.worker:
  sys.path.insert(0,str(CELLTYPE));sys.path.insert(0,str(TYPING))
  if args.worker=="uniform":uniform_worker(out,args.arm,donors,args.references_from)
  else:
   import refine
   refine.SHIFTS=out/args.arm/"shifts.json"
   sys.argv=["refine.py","--out",str(out/args.arm/"typing/refinement"),"--donors",*donors,
    "--cd8-alignment","unshifted","--shifts",str(out/args.arm/"shifts.json")]
   refine.main()
  return
 if not args.skip_audit:source_audit()
 if args.audit_only:return
 if not args.skip_prepare:prepare(out,donors)
 if args.prepare_only:return
 for arm in ["Cell_control","own_compartment"]:
  env=os.environ.copy();env.update(IO60_EXPRESSION_ROOT=str(out/arm/"expression"),
   IO60_ALIGNED_ROOT=str(out/arm/"uniform"),IO60_SCALE_PATH=str(out/arm/"scale.parquet"),
   IO60_TYPING_ROOT=str(out/arm/"typing"),IO60_GEOMETRY_ROOT=str(out/"geometry"),
   IO60_SUBPANEL_LOG=str(GEOMETRY/"subpanel/subpanel_log.json"),
   IO60_REFINE_CD8_ALIGNMENT="unshifted",IO60_SHIFTS_PATH=str(out/arm/"shifts.json"))
  subprocess_stage(out,arm,"uniform",[sys.executable,str(Path(__file__).resolve()),"--worker","uniform","--out",str(out),"--arm",arm,"--donors",*donors],env)
  typing(out,arm,donors,env)
 writej(out/"comparison_status.json",dict(status="processing_complete",donors=donors,
  candidate_release_selected=False,baseline_unchanged=True,
  interpretation="Both arms independently calibrate the same donor cohort; this estimates a compartment repair plus reprocessing effect, without outcome-tuned UMAP or classifier relabeling."))
if __name__=="__main__":main()
