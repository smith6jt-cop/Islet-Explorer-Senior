"""A new, source-bound review workspace over completed sequential analyses.

Review notes never change scientific tables or activate a release. Images and
large tables remain in their original packages; this workspace records their
hashes and a portable export of the reviewer's progress.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from multiplex5000.artifacts import sha256_file as sha


SCHEMA = "multiplex5000-sequential-review/1"
DECISIONS = {"unreviewed", "evidence_checked", "needs_reassessment", "uncertain"}




def _safe_child(root, relative):
    path = (root / relative).resolve()
    if root.resolve() not in path.parents:
        raise ValueError(f"Artifact escapes its source package: {relative}")
    return path


def _records(frame):
    return json.loads(frame.to_json(orient="records"))


def build_workspace(recipe, output, *, verification, recipe_path):
    """Render source data plus editable notes into an already-created fresh output."""
    output = Path(output).resolve()
    base = Path(recipe["base_bundle"]["path"]).resolve()
    current = Path(recipe["current_islets"]["path"]).resolve()
    display = recipe["current_islets"]["display_directory"]
    artifacts = {}

    def link(path, *, expected=None):
        path = path.resolve()
        actual = sha(path)
        if expected is not None and actual != expected:
            raise ValueError(f"Changed review evidence: {path}")
        artifacts[str(path)] = {"sha256": actual, "bytes": path.stat().st_size}
        return quote(os.path.relpath(path, output), safe="/")

    manifest = {r["path"]: r for r in json.loads((current / "artifact_manifest.json").read_text())}

    def evidence(relative):
        if relative not in manifest:
            raise ValueError(f"Review evidence is absent from the pinned manifest: {relative}")
        return link(_safe_child(current, relative), expected=manifest[relative]["sha256"])

    current_file = current / "all_islet_relationships.parquet"
    base_file = base / "all_pairs.parquet"
    for file in [current_file, base_file, current / "all_donor_comparability.csv",
                 current / "artifact_manifest.json", base / "artifact_manifest.json"]:
        link(file)
    pairs = pd.read_parquet(current_file)
    cases = []
    for r in _records(pairs):
        pid = r["pair_id"]
        images = []
        if r["image_reviewed_in_current_scope"]:
            images = [
                {"label": "Current display (green contribution reduced)",
                 "url": evidence(f"{display}/figures/{pid}.jpg")},
                {"label": "Original image used for the recorded review",
                 "url": evidence(f"pair_figures/{pid}.jpg")},
            ]
        elif r.get("diagnostic_image_inspected"):
            images = [{"label": "Diagnostic image; numerical criteria remain unmet",
                       "url": evidence(f"{display}/diagnostic_figures/{pid}.jpg")}]
        detail = {k: r.get(k) for k in (
            "reason_code", "assessment_basis", "anatomy_observation", "geometry_observation",
            "competing_alternative_assessment", "numerical_screen_reasons", "reviewer", "review_type",
            "second_review_rationale", "distance_um", "RNA_area_um2", "protein_area_um2",
            "RNA_n_whole_cells", "protein_n_whole_cells", "RNA_n_endocrine_seed_cells",
            "historical_verdict", "historical_NCC_z", "zero_projected_overlap",
            "diagnostic_anatomy_observation", "diagnostic_scope")}
        cases.append({"id": "islet:" + pid, "label": pid, "generation": "current_islets",
                      "donor": str(r["donor_id"]), "moving_id": r["moving_id"], "fixed_id": r["fixed_id"],
                      "baseline": r["comparability_designation"], "image_available": bool(images),
                      "images": images, "detail": detail})
    # Keep the corrected multicell generation discoverable, without presenting
    # its 1,009 candidates as the current islet comparability decision set.
    for r in _records(pd.read_parquet(base_file)):
        key = json.dumps([str(r["donor_id"]), r["moving_id"], r["fixed_id"]], separators=(",", ":"))
        cid = "multicell:" + hashlib.sha256(key.encode()).hexdigest()[:24]
        cases.append({"id": cid, "label": str(r["donor_id"]) + " · " + r["RNA_object_class"] +
                      " → " + r["protein_object_class"], "generation": "corrected_multicell",
                      "donor": str(r["donor_id"]), "moving_id": r["moving_id"], "fixed_id": r["fixed_id"],
                      "baseline": "Exploratory review candidate" if r["eligible_for_exploratory_regional_review"]
                      else "Diagnostic candidate", "image_available": False, "images": [],
                      "detail": {k: r.get(k) for k in ("RNA_object_class", "protein_object_class",
                          "distance_um", "RNA_area_um2", "protein_area_um2", "zero_projected_overlap",
                          "ambiguous_geometry", "RNA_proposal_sensitivity_flag", "image_sampling_review_flag",
                          "anatomy_alternative_review_flag", "correspondence_policy")}})
    if len({r["id"] for r in cases}) != len(cases):
        raise ValueError("Duplicate review case identity")
    roster = _records(pd.read_csv(current / "all_donor_comparability.csv", dtype={"donor_id": str, "source_serial": str}))
    links = [{"label": "Current Islet-only comparison report", "url": evidence("REPORT.html")},
             {"label": "Corrected multicell report and subsequent revisions", "url": link(base / "REPORT.html")},
             {"label": "Current pair table", "url": evidence("all_islet_relationships.csv")},
             {"label": "Current comparable subset", "url": evidence("comparable_islet_relationships.csv")},
             {"label": "Native Single counts (within-section only)", "url": link(base / "native_single_counts.csv")}]
    related = []
    for item in recipe.get("related_reports", []):
        # These are dated navigation observations, not execution dependencies.
        p = Path(item["path"])
        actual = sha(p) if p.is_file() else None
        related.append({**item, "url": quote(os.path.relpath(p, output), safe="/"),
                        "current_sha256": actual, "changed_since_observation": actual != item.get("observed_sha256")})
    identity = {"recipe_sha256": sha(recipe_path), "cases": [{k: r[k] for k in
                ("id", "generation", "moving_id", "fixed_id", "baseline")} for r in cases], "artifacts": artifacts}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    data = {"schema": SCHEMA, "dataset_fingerprint": fingerprint, "created_utc": datetime.now(timezone.utc).isoformat(),
            "cases": cases, "roster": roster, "links": links, "related_reports": related,
            "source_artifacts": artifacts, "recipe_path": str(Path(recipe_path).resolve()),
            "recipe_sha256": sha(recipe_path), "counts": dict(Counter(r["generation"] for r in cases)),
            "scientific_designations_changed": False, "production_activated": False}
    (output / "review_data.json").write_text(json.dumps(data, indent=2) + "\n")
    embedded = json.dumps(data).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    (output / "REVIEW.html").write_text(HTML.replace("__DATA__", embedded))
    (output / "review_manifest.json").write_text(json.dumps({"schema": SCHEMA,
        "dataset_fingerprint": fingerprint, "files": {n: sha(output / n) for n in ("review_data.json", "REVIEW.html")},
        "recipe_sha256": sha(recipe_path), "source_artifacts": artifacts,
        "technical_verification_status": verification["status"], "scientific_designations_changed": False}, indent=2) + "\n")
    return {"review": str(output / "REVIEW.html"), "dataset_fingerprint": fingerprint,
            "counts": data["counts"], "source_artifacts_bound": len(artifacts), "new_reviewer_decisions": 0}


def validate_progress(data, exported):
    """Validate exported notes without treating them as scientific acceptance."""
    if not isinstance(exported, dict):
        raise ValueError("Review export must be a mapping")
    if exported.get("schema") != SCHEMA or exported.get("dataset_fingerprint") != data["dataset_fingerprint"]:
        raise ValueError("Review export belongs to different evidence; rebase requires explicit review")
    if exported.get("production_activated") or exported.get("scientific_designations_changed"):
        raise ValueError("Review notes cannot change scientific designations or activate a release")
    for field in ("history", "import_history"):
        if field in exported and (not isinstance(exported[field], list) or
                                  not all(isinstance(r, dict) for r in exported[field])):
            raise ValueError(f"Malformed {field}")
    decisions = exported.get("decisions")
    if not isinstance(decisions, dict):
        raise ValueError("Review decisions must be a mapping")
    ids = {r["id"] for r in data["cases"]}
    for key, record in decisions.items():
        if key not in ids or not isinstance(record, dict) or record.get("status") not in DECISIONS:
            raise ValueError(f"Invalid review decision: {key}")
        if not isinstance(record.get("note", ""), str) or not isinstance(record.get("reviewer", ""), str):
            raise ValueError(f"Malformed reviewer/note: {key}")
        if record["status"] != "unreviewed" and not record.get("reviewer", "").strip():
            raise ValueError(f"Reviewer required: {key}")
        if record["status"] != "unreviewed":
            stamp = record.get("updated_utc")
            try:
                if not isinstance(stamp, str) or datetime.fromisoformat(stamp.replace("Z", "+00:00")).tzinfo is None:
                    raise ValueError("Timestamp must contain a timezone")
            except (ValueError, TypeError) as error:
                raise ValueError(f"Valid review timestamp required: {key}") from error
    counts = Counter(r["status"] for r in decisions.values())
    counts["unreviewed"] += len(ids) - len(decisions)
    by_generation = {}
    for generation in sorted({r["generation"] for r in data["cases"]}):
        subset = [r for r in data["cases"] if r["generation"] == generation]
        by_generation[generation] = dict(Counter(decisions.get(r["id"], {}).get("status", "unreviewed") for r in subset))
    return {"status": "VALID_REVIEW_PROGRESS", "counts": dict(counts), "by_generation": by_generation,
            "scientific_designations_changed": False, "production_activated": False}


def review_status(directory, decisions):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "review_manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if sha(_safe_child(directory, name)) != expected:
            raise ValueError(f"Changed review workspace: {name}")
    data = json.loads((directory / "review_data.json").read_text())
    exported = json.loads(Path(decisions).read_text())
    result = validate_progress(data, exported)
    result.update({"review_dir": str(directory), "decisions": str(Path(decisions).resolve()),
                   "decisions_sha256": sha(decisions), "dataset_fingerprint": data["dataset_fingerprint"]})
    # Do not silently affirm progress against changed linked figures/tables.
    changed = [p for p, pin in manifest["source_artifacts"].items()
               if not Path(p).is_file() or sha(p) != pin["sha256"]]
    if changed:
        raise ValueError(f"Source evidence changed since workspace creation: {changed[:5]}")
    result["linked_evidence_hashes_verified"] = len(manifest["source_artifacts"])
    return result


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Sequential Xenium–PhenoCycler review</title>
<style>body{font:16px system-ui,sans-serif;color:#183337;background:#f3f5f1;margin:0}header,main{max-width:1450px;margin:auto;padding:24px}h1{font-size:30px;margin:4px 0 12px}h2{font-size:21px}p{line-height:1.5}a{color:#005c70}button,input,select,textarea{font:inherit;padding:9px;border:1px solid #a5b7b4;border-radius:5px;background:white}button{cursor:pointer}button.primary{background:#0f6162;color:white}button:disabled{opacity:.5}label{display:block;font-size:14px;margin:5px 0}nav,.filters,.actions{display:flex;gap:12px;flex-wrap:wrap;align-items:end}.banner{padding:14px;background:#fff1ce;border-left:4px solid #c68610}.layout{display:grid;grid-template-columns:320px 1fr;gap:22px}section,.panel{background:white;border:1px solid #d8e1db;border-radius:8px;padding:18px;margin-bottom:18px}#list{max-height:800px;overflow:auto}#list button{width:100%;text-align:left;margin:4px 0;font-size:14px}#list .selected{border:2px solid #0f6162;background:#e7f4ed}small{color:#4b6467}.ids{word-break:break-all;font:12px ui-monospace,monospace}.evidence{width:100%;display:block;border:1px solid #ddd}textarea{box-sizing:border-box;width:100%;min-height:110px}dt{font-weight:600;margin-top:9px}dd{margin:4px 0;white-space:pre-wrap;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}#message{white-space:pre-wrap}.warning{color:#8d320a}summary{cursor:pointer}#count,#progress{font-weight:600}details{margin:10px 0}@media(max-width:850px){.layout{grid-template-columns:1fr}#list{max-height:280px}header,main{padding:14px}}</style></head><body>
<header><small>MultiPlex5000 · sequential sections · review workspace</small><h1>Xenium–PhenoCycler evidence and progress</h1>
<p id="summary">The current Islet-only comparison and earlier corrected multicell generation remain separately searchable. Donors, sections, native object IDs and reviewed evidence stay attached to each case.</p>
<p class="banner">Recorded assessments include nonblinded AI image review. Your notes are a separate review record. They do not change comparability designations, establish cell-to-cell identity, accept a transform, or activate a release. “Not comparable” means insufficient or conflicting evidence. “Not assessable” means a failed or missing transform.</p>
<nav id="links"></nav><details><summary>Analysis history and downstream reports</summary><div id="related"></div></details>
<details><summary id="rosterLabel">Donor scope: paired and not assessable</summary><div style="overflow:auto"><table id="roster"></table></div></details>
</header><main><section><h2>Review progress</h2><div class="actions"><div><label for="reviewer">Your name or initials</label><input id="reviewer" autocomplete="off"></div><button id="export" class="primary">Export progress JSON</button><label>Resume saved progress<input id="import" type="file" accept="application/json,.json"></label></div><p id="progress"></p><p id="message" aria-live="polite">Notes save locally when you select Save note. Export progress before closing or moving this workspace.</p></section>
<section class="filters"><div><label for="generation">Generation</label><select id="generation"><option value="current_islets">Current Islet-only comparison</option><option value="corrected_multicell">Earlier corrected multicell candidates</option></select></div><div><label for="scope">Cases</label><select id="scope"><option value="images">With inspected image evidence</option><option value="all">All relationships</option><option value="comparable">Comparable / exploratory subset</option><option value="flagged">Needs reassessment / uncertain notes</option><option value="unreviewed">No new reviewer note</option></select></div><div><label for="donor">Donor</label><select id="donor"><option value="">All donors</option></select></div><div><label for="search">Native ID, pair ID or assessment</label><input id="search" type="search"></div><p id="count"></p></section>
<div class="layout"><aside class="panel"><div class="actions"><button id="prev">Previous page</button><button id="next">Next page</button></div><div id="list"></div></aside><section id="case"><h2 id="caseTitle"></h2><p id="baseline"></p><p class="ids" id="ids"></p><div id="images"></div><details open><summary>Recorded evidence and limitations</summary><dl id="detail"></dl></details><label for="decision">Your review status</label><select id="decision"><option value="unreviewed">Unreviewed</option><option value="evidence_checked">Evidence checked</option><option value="needs_reassessment">Needs scientific reassessment</option><option value="uncertain">Uncertain</option></select><label for="note">Observation, concern or next check</label><textarea id="note"></textarea><button class="primary" id="save">Save note</button><small> Saving a note does not rerun analysis.</small><details><summary>Saved note history</summary><pre id="history" style="white-space:pre-wrap"></pre></details></section></div></main>
<script id="data" type="application/json">__DATA__</script><script>
'use strict';const DATA=JSON.parse(document.getElementById('data').textContent);const $=id=>document.getElementById(id);const storageKey=DATA.schema+':'+DATA.dataset_fingerprint;const allowed=new Set(['unreviewed','evidence_checked','needs_reassessment','uncertain']);const ids=new Set(DATA.cases.map(r=>r.id));let state={schema:DATA.schema,dataset_fingerprint:DATA.dataset_fingerprint,decisions:{},history:[],import_history:[]};let selected=null,page=0,dirty=false;const pageSize=30;
function message(t,bad=false){$('message').textContent=t;$('message').className=bad?'warning':''}
function valid(x){if(!x||typeof x!=='object'||Array.isArray(x))throw Error('Review export must be an object');if(x.production_activated||x.scientific_designations_changed)throw Error('Notes cannot change scientific results or activate a release');if(x.schema!==DATA.schema||x.dataset_fingerprint!==DATA.dataset_fingerprint)throw Error('Different evidence fingerprint. This export cannot be applied to this workspace.');if(!x.decisions||Array.isArray(x.decisions)||typeof x.decisions!=='object')throw Error('Missing decision mapping');for(const [id,r] of Object.entries(x.decisions)){if(!ids.has(id)||!r||!allowed.has(r.status)||typeof(r.note??'')!=='string'||typeof(r.reviewer??'')!=='string'||(r.status!=='unreviewed'&&(!(r.reviewer||'').trim()||typeof r.updated_utc!=='string'||!/(Z|[+-]\d\d:\d\d)$/.test(r.updated_utc)||!Number.isFinite(Date.parse(r.updated_utc)))))throw Error('Invalid note: '+id)}if(x.history&&(!Array.isArray(x.history)||!x.history.every(r=>r&&typeof r==='object'&&!Array.isArray(r))))throw Error('Invalid note history');if(x.import_history&&(!Array.isArray(x.import_history)||!x.import_history.every(r=>r&&typeof r==='object'&&!Array.isArray(r))))throw Error('Invalid import history');return x}
try{const saved=localStorage.getItem(storageKey);if(saved)state=valid(JSON.parse(saved))}catch(e){message('Local progress could not be loaded: '+e.message,true)}
function persist(){try{localStorage.setItem(storageKey,JSON.stringify(state));return true}catch(e){message('Local saving unavailable. Your notes remain only in this open page. Export progress now. '+e.message,true);return false}}
function progress(){let rows=DATA.cases.filter(r=>r.generation===$('generation').value),done=rows.filter(r=>(state.decisions[r.id]?.status||'unreviewed')!=='unreviewed').length,flag=rows.filter(r=>['needs_reassessment','uncertain'].includes(state.decisions[r.id]?.status)).length;$('progress').textContent=done+' of '+rows.length+' cases have a new reviewer note; '+flag+' need follow-up in this generation.'}
const currentCases=DATA.cases.filter(r=>r.generation==='current_islets');$('summary').textContent='Current Islet-only comparison: '+currentCases.filter(r=>r.baseline==='Comparable as the same islet').length+' comparable relationships from '+currentCases.length+' candidates. The earlier corrected multicell generation remains separately searchable. Each case retains its donor, native IDs and evidence.';$('rosterLabel').textContent='All '+DATA.roster.length+' donors: paired and not-assessable scope';
for(const l of DATA.links){const a=document.createElement('a');a.href=l.url;a.textContent=l.label;$('links').append(a)}
for(const l of DATA.related_reports){const p=document.createElement('p'),a=document.createElement('a');a.href=l.url;a.textContent=l.label;p.append(a,document.createTextNode(' — '+l.status+'. '+l.scope+(l.changed_since_observation?' Changed since the dated observation; consult its current receipt.':'')));$('related').append(p)}
let header=$('roster').createTHead().insertRow();for(const s of ['Donor','RNA section key','Registration','Availability','Candidates','Comparable']){const c=document.createElement('th');c.textContent=s;header.append(c)}for(const r of DATA.roster){let row=$('roster').insertRow();for(const v of [r.donor_id,r.source_serial,r.historical_registration,r.assessment_status,r.candidate_pairs,r.comparable_as_same_islet])row.insertCell().textContent=v===null?'Not assessable':v}
for(const d of [...new Set(DATA.cases.map(r=>r.donor))].sort((a,b)=>a.localeCompare(b,undefined,{numeric:true}))){const o=new Option(d,d);$('donor').add(o)}
function filtered(){const query=$('search').value.toLowerCase(),scope=$('scope').value;return DATA.cases.filter(r=>r.generation===$('generation').value&&(!$('donor').value||r.donor===$('donor').value)&&(!query||JSON.stringify(r).toLowerCase().includes(query))&&(scope!=='images'||r.image_available)&&(scope!=='comparable'||['Comparable as the same islet','Exploratory review candidate'].includes(r.baseline))&&(scope!=='flagged'||['needs_reassessment','uncertain'].includes(state.decisions[r.id]?.status))&&(scope!=='unreviewed'||!state.decisions[r.id]||state.decisions[r.id].status==='unreviewed'))}
function select(r){if(dirty){message('Save your current note before selecting another case.',true);return}selected=r;$('case').hidden=false;$('caseTitle').textContent=r.label;$('baseline').textContent='Recorded assessment: '+r.baseline;$('ids').textContent='Donor '+r.donor+' | RNA '+r.moving_id+' | PhenoCycler '+r.fixed_id;$('images').replaceChildren();for(const im of r.images){const h=document.createElement('h3');h.textContent=im.label;const a=document.createElement('a');a.href=im.url;const img=document.createElement('img');img.src=im.url;img.alt=im.label+' for '+r.label;img.className='evidence';img.loading='lazy';a.append(img);$('images').append(h,a)}if(!r.images.length){const p=document.createElement('p');p.textContent='No pair-level inspected image is attached to this case in this generation. Absence of a displayed image is not a negative biological finding.';$('images').append(p)}$('detail').replaceChildren();for(const [k,v] of Object.entries(r.detail)){if(v===null||v==='')continue;const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=k.replaceAll('_',' ');dd.textContent=typeof v==='object'?JSON.stringify(v):String(v);$('detail').append(dt,dd)}const note=state.decisions[r.id]||{};$('decision').value=note.status||'unreviewed';$('note').value=note.note||'';if(note.reviewer&&!$('reviewer').value)$('reviewer').value=note.reviewer;$('history').textContent=JSON.stringify((state.history||[]).filter(h=>h.id===r.id),null,2);drawList()}
function drawList(){let rows=filtered();const max=Math.max(0,Math.ceil(rows.length/pageSize)-1);page=Math.min(page,max);$('count').textContent=rows.length+' matching cases';$('prev').disabled=page===0;$('next').disabled=page>=max;$('list').replaceChildren();for(const r of rows.slice(page*pageSize,(page+1)*pageSize)){const b=document.createElement('button');b.textContent=r.label+' · '+r.baseline+' · '+(state.decisions[r.id]?.status||'unreviewed');if(selected?.id===r.id)b.className='selected';b.onclick=()=>select(r);$('list').append(b)}progress()}
function filterChanged(){if(dirty){message('Your current note is unsaved. Save it before changing cases.',true)}page=0;drawList();if(!dirty){const r=filtered()[0];if(r)select(r);else $('case').hidden=true}}
for(const id of ['generation','scope','donor','search'])$(id).addEventListener(id==='search'?'input':'change',()=>{if(id==='generation'&&$('generation').value==='corrected_multicell'&&$('scope').value==='images')$('scope').value='all';filterChanged()});$('prev').onclick=()=>{page--;drawList()};$('next').onclick=()=>{page++;drawList()};for(const id of ['decision','note'])$(id).oninput=()=>{dirty=true};
$('save').onclick=()=>{if(!selected)return;const reviewer=$('reviewer').value.trim();if(!reviewer){message('Enter your name or initials before saving a note.',true);return}const previous=state.decisions[selected.id]||null;const next={status:$('decision').value,note:$('note').value,reviewer,updated_utc:new Date().toISOString()};state.history=state.history||[];state.history.push({id:selected.id,previous,next});state.decisions[selected.id]=next;dirty=false;if(persist())message('Note saved locally. Export progress to keep a separate copy.');select(selected)};
function exportPayload(){return {...state,exported_utc:new Date().toISOString(),scientific_designations_changed:false,production_activated:false}}
$('export').onclick=()=>{if(dirty){message('Save the current note before exporting.',true);return}const blob=new Blob([JSON.stringify(exportPayload(),null,2)+'\n'],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='sequential-review-progress-'+new Date().toISOString().replaceAll(':','-')+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);message('Progress exported. Keep the JSON file outside the immutable source packages.')};
function importProgress(x){valid(x);const conflicts=Object.keys(x.decisions).filter(id=>state.decisions[id]&&JSON.stringify(state.decisions[id])!==JSON.stringify(x.decisions[id]));if(conflicts.length)throw Error('Conflicting saved notes for '+conflicts.length+' cases. Current notes were preserved; use a separate workspace to compare exports.');const prior=JSON.parse(JSON.stringify(x));state.import_history=state.import_history||[];state.import_history.push({imported_utc:new Date().toISOString(),source_export:prior});state.decisions={...x.decisions,...state.decisions};state.history=[...new Map([...(x.history||[]),...(state.history||[])].map(h=>[JSON.stringify(h),h])).values()];const stored=persist();drawList();if(selected)select(selected);if(stored)message('Saved progress resumed. Imported history is retained.');}
$('import').onchange=async e=>{try{if(dirty)throw Error('Save the current note before importing');const f=e.target.files[0];if(f)importProgress(JSON.parse(await f.text()))}catch(err){message(err.message,true)}finally{e.target.value=''}};
window.addEventListener('beforeunload',e=>{if(dirty){e.preventDefault();e.returnValue=''}});filterChanged();
</script></body></html>'''
