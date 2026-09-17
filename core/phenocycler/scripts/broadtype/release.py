"""The A.4a release writer: one arm's calls -> a development-stamped identity release that
``phenocycler.integration.phenospace.load_release`` accepts.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/release.py \\
        --arm <producer> --donors 6436 6442 [--scope donor|frame] [--out data/development/broadtype/release]

Writes, under ``--out``::

    release.json                                   the contract phenospace validates
    label_space.json                               the 7 declared labels + Ambiguous + Unavailable
    calls/donor_id=<d>/roi=panc/data.parquet       one row per row of runs/<run>/expression/donor_id=<d>

Why this file exists.  The bake-off's arms write *screen records* (``records.CALL_COLUMNS``:
object_id, label, score, reason) and the gate writes *verdicts* (which declared labels were
verified on which donor).  Neither is something the integration layer may read: phenospace
refuses any producer vocabulary it has not been handed as a JSON declaration, and
``donor._evidence_mask`` refuses an ``evidence_positive`` with a single NA.  This module is the
one place that turns (calls, verdict, declarations) into that contract, so the integration and
app bridge consume exactly one artifact and never an arm's raw output.

What the release claims, and does not.  ``stamp = development`` and
``latch = NO_VERIFIED_TYPING (unchanged)``: the marker registry's latch is not lifted, no
production output is written, and nothing under ``data/runs/`` is touched.  A cell is
``evidence_positive`` only when the arm called a declared label AND the verdict verified that
label for that arm on that donor (held-out enrichment + spatial coherence + stability, see
``gate.judge_label``).  Every other cell -- Ambiguous, Unavailable, an unverified label -- is
``False``, never NA, and the per-cell ``reason`` says why:

    ''                                         called, label verified on this donor
    'label_not_verified'                       the arm called a declared label the verdict did not verify;
                                               lineage_native keeps the label, evidence_positive is False
    'qc_analysis_ineligible:<analysis_reason>' geometry QC excluded the cell from analysis; the cell is
                                               Unavailable regardless of what the arm said, because a
                                               measurement on a failed geometry is not a measurement
    'outside_screen_frame'                     frame-scope release only: the cell is in expression but
                                               not in the 16-core screen frame the arm was run on
    <the arm's own reason>                     Ambiguous / Unavailable as the arm recorded it

Precedence: QC ineligibility first (it is a property of the cell, not of the arm), then frame
coverage, then the arm's label, then verification.  A frame-scope release is allowed for
development so the downstream plumbing can be exercised before a whole-donor run exists; it
is not a whole-donor call and the release records ``scope`` so nothing can mistake it for one.

The state tier is respected: ``ki67_donor_rank`` is the within-donor rank of the Step-2
aligned Ki67 (a state marker, nucleus-measured; ``panel.NUCLEUS_MEASURED``) appended AFTER
identity, exactly as ``arm_base.state_columns`` does; it never enters a label.  ``cell_region``
and ``qupath_islet_num`` are QuPath audit columns carried for audit and refused as identity
inputs everywhere else in this lane.

Verdict shape.  ``verdict.py`` reads ``verdict/<arm>/donor_id=<d>.json`` (frame scope) or
``donor_id=<d>__<scope>.json`` (any other scope) as ``labels`` mapping label -> the
``gate.judge_label`` result ({"verified": bool, ...}, optionally classified as
{"status": "verified" | "unjudgeable" | "failed"}).  This module reads the same files the same
way, also accepting a ``labels`` list of such objects each carrying ``label``; a donor-scope
release prefers the scope-matched verdict and otherwise uses the frame verdict the A.4a
contract names, recording which.  ``verified`` must be a JSON boolean (or ``status`` one of the
three strings, and the two may not contradict); anything else, a missing file, or a verdict
naming a different arm or donor stops the writer.  Absence of evidence is never a pass.

No provenance string may carry a retired-lane token (``phenospace.RETIRED_INPUT_TOKENS``) or
``/assignments/``; the writer scans everything it writes before phenospace does.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from arm_base import within_donor_rank  # noqa: E402  (the one rank definition the lane uses)

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402
from phenocycler.integration.phenospace import RETIRED_INPUT_TOKENS  # noqa: E402

SCHEMA_VERSION = 1
MODALITY = "phenocycler"
ROI = "panc"                     # the run holds 20 pancreas sections; pLN sections are not in it
CELL_KEY = "object_id"
LATCH = "NO_VERIFIED_TYPING (unchanged)"
STAMP = "development"
STATE_RANK_MARKER = "Ki67"
COLUMNS = {
    "label": "lineage_native",
    "common_lineage": "common_lineage",
    "evidence_positive": "evidence_positive",
    "confidence": "confidence",
    "fine_type": "celltype_fine",
}
#: Column order of calls/donor_id=<d>/roi=panc/data.parquet.
CALL_COLUMNS = ("object_id", "lineage_native", "common_lineage", "evidence_positive", "confidence",
                "celltype_fine", "reason", "ki67_donor_rank", "cell_region", "qupath_islet_num")
STATUS_VERIFIED = "verified"     # verdict.py's label status vocabulary
REASON_QC = "qc_analysis_ineligible"
REASON_OUTSIDE_FRAME = "outside_screen_frame"
REASON_NOT_VERIFIED = "label_not_verified"
REASON_VOCABULARY = {
    "": "called; the label was verified for this arm on this donor",
    REASON_NOT_VERIFIED: "the arm called a declared label the verdict did not verify on this donor; "
                         "lineage_native keeps the label, evidence_positive is False",
    f"{REASON_QC}:<analysis_reason>": "geometry QC excluded the cell from analysis (geometry_qc analysis_reason); "
                                      "Unavailable regardless of the arm's output",
    REASON_OUTSIDE_FRAME: "frame-scope release: the cell is in expression but outside the screen frame",
    "<arm reason>": "Ambiguous / Unavailable exactly as the arm recorded it",
}


class ReleaseError(RuntimeError):
    """The inputs cannot be turned into a valid release; nothing partial is left behind."""


@dataclass
class Sources:
    """Everything the writer reads, so a test can point every root at a tmp tree.

    Defaults are the lane's pins.  ``repo_root`` is only used to write repo-relative paths
    into release.json; a path outside it is recorded absolute.
    """
    arms_root: pathlib.Path = field(default_factory=lambda: pins.OUT / "arms")
    verdict_root: pathlib.Path = field(default_factory=lambda: pins.OUT / "verdict")
    expression_root: pathlib.Path = field(default_factory=lambda: pins.EXPRESSION)
    expression_manifest: pathlib.Path = field(default_factory=lambda: pins.EXPRESSION_MANIFEST)
    geometry_qc_root: pathlib.Path = field(default_factory=lambda: pins.GEOMETRY_QC)
    uniform_root: pathlib.Path = field(default_factory=lambda: pins.UNIFORM)
    shifts_path: pathlib.Path = field(default_factory=lambda: pins.SHIFTS)
    declarations_path: pathlib.Path = field(default_factory=lambda: declarations.DECLARATIONS_JSON)
    run_id: str = pins.RUN_ID
    expression_content_id: str = pins.EXPRESSION_CONTENT_ID
    expression_output_sha256: str = pins.EXPRESSION_OUTPUT_SHA256
    repo_root: pathlib.Path = field(default_factory=lambda: pins.REPO)

    def calls_path(self, arm: str, donor: str, scope: str) -> pathlib.Path:
        # records.write_calls: config_id 'primary' -> calls__primary.parquet
        return self.arms_root / arm / f"donor_id={donor}" / scope / "calls__primary.parquet"

    def verdict_path(self, arm: str, donor: str, scope: str = "frame") -> pathlib.Path:
        """The verdict that speaks for these calls.

        ``verdict.py`` writes frame verdicts as ``donor_id=<d>.json`` and any other scope with a
        ``__<scope>`` suffix.  A scope-matched verdict is preferred when it exists; otherwise the
        frame verdict -- the screen's judgement, the file the A.4a contract names -- is used, and
        release.json records which one so the two can never be confused.
        """
        d = self.verdict_root / arm
        if scope != "frame":
            scoped = d / f"donor_id={donor}__{scope}.json"
            if scoped.exists():
                return scoped
        return d / f"donor_id={donor}.json"

    def verdict_summary_path(self) -> pathlib.Path:
        return self.verdict_root / "summary.json"

    def relpath(self, p: pathlib.Path) -> str:
        p = pathlib.Path(p).resolve()
        try:
            return str(p.relative_to(self.repo_root.resolve()))
        except ValueError:
            return str(p)


# --------------------------------------------------------------------------- small helpers
def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _strings(obj) -> list:
    """Every string anywhere in a JSON-like object (keys included)."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.append(str(k))
            out.extend(_strings(v))
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            out.extend(_strings(v))
    elif isinstance(obj, str):
        out.append(obj)
    return out


def scan_provenance(payload, *, what: str) -> None:
    """Refuse any string that names a retired identity input.

    phenospace checks the lane and the ``inputs[].path`` values; this checks EVERY string in
    what we write, so a retired token cannot hide in a verdict path or an arm name.
    """
    for s in _strings(payload):
        low = s.replace("\\", "/").lower()
        for tok in RETIRED_INPUT_TOKENS:
            if tok.lower() in low:
                raise ReleaseError(f"{what}: provenance string {s!r} names retired input token {tok!r}")
        if "/assignments/" in low:
            raise ReleaseError(f"{what}: provenance string {s!r} names runs/<id>/assignments, which is retired")


def release_id_for(arm: str, today: str | None = None) -> str:
    return f"broadtype-{arm}-dev-{today or time.strftime('%Y%m%d', time.gmtime())}"


def lane_for(arm: str) -> str:
    return f"broadtype-{arm}"


# --------------------------------------------------------------------------- the verdict
def verified_labels(payload: dict, *, arm: str, donor: str, declared: set) -> set:
    """The declared labels the verdict verified for ``arm`` on ``donor``.

    Accepts ``labels`` as a mapping label -> {"verified": bool, ...} (gate.judge_label keyed
    by label) or as a list of such objects carrying ``label``.  ``verified`` must be a JSON
    boolean.  A verdict that names another arm or donor, an undeclared label, or no labels
    at all is an error -- a silently empty set would turn every call into
    ``label_not_verified`` and look like a legitimate result.
    """
    if not isinstance(payload, dict):
        raise ReleaseError("verdict must be one JSON object")
    for key, want in (("arm", arm), ("donor_id", donor)):
        if key in payload and str(payload[key]) != str(want):
            raise ReleaseError(f"verdict names {key}={payload[key]!r}, expected {want!r}")
    labels = payload.get("labels")
    if isinstance(labels, dict):
        items = [({"label": k, **v} if isinstance(v, dict) else {"label": k, "verified": v})
                 for k, v in labels.items()]
    elif isinstance(labels, list):
        items = labels
    else:
        raise ReleaseError("verdict has no 'labels' mapping or list")
    out = set()
    for it in items:
        if not isinstance(it, dict) or "label" not in it:
            raise ReleaseError(f"verdict label entry is not an object with 'label': {it!r}")
        L = str(it["label"])
        if L not in declared:
            raise ReleaseError(f"verdict judges undeclared label {L!r}")
        v = it.get("verified", None)
        status = it.get("status", None)
        if isinstance(v, bool):
            ok = v
        elif status in (STATUS_VERIFIED, "unjudgeable", "failed"):
            ok = status == STATUS_VERIFIED          # verdict.py's classified form
        else:
            raise ReleaseError(f"verdict for {L!r} has non-boolean 'verified' {v!r} and no status")
        if isinstance(v, bool) and status is not None and (status == STATUS_VERIFIED) != v:
            raise ReleaseError(f"verdict for {L!r}: verified={v!r} contradicts status={status!r}")
        if ok:
            out.add(L)
    return out


# --------------------------------------------------------------------------- the label space
def label_space_payload(dec: dict, arm: str, *, declarations_sha256: str) -> dict:
    """The strict declaration phenospace loads; same shape as legacy_broad_space.json."""
    labels = []
    for L, spec in dec["labels"].items():
        labels.append({"label": L, "common_lineage": spec["common_lineage"], "reason": spec["reason"],
                       "markers": list(spec["markers"])})
    labels.append({"label": records.AMBIGUOUS, "common_resolution": "ambiguous",
                   "reason": "the arm did not resolve one declared label; the per-cell reason column says why"})
    labels.append({"label": records.UNAVAILABLE, "common_resolution": "unavailable",
                   "reason": "no identity call: geometry-QC ineligible, outside the screen frame, or the arm's "
                             "own recorded reason; never a negative"})
    return {
        "schema_version": SCHEMA_VERSION,
        "space_id": f"broadtype-{arm}-broad7-dev",
        "producer_lane": lane_for(arm),
        "declarations": {"method_version": dec["method_version"], "sha256": declarations_sha256,
                         "undeclared_labels": dict(dec.get("undeclared_labels", {}))},
        "labels": labels,
    }


# --------------------------------------------------------------------------- one donor
def _read_expression(src: Sources, donor: str) -> pd.DataFrame:
    p = src.expression_root / f"donor_id={donor}" / "data.parquet"
    if not p.exists():
        raise ReleaseError(f"donor {donor}: expression partition missing: {p}")
    cols = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    expr = pd.read_parquet(p, columns=cols)
    expr["object_id"] = expr["object_id"].astype(str)
    if expr["object_id"].duplicated().any():
        raise ReleaseError(f"donor {donor}: duplicate object_id in expression")
    return expr


def _read_qc_reason(src: Sources, donor: str, ids: pd.Series) -> pd.Series:
    """geometry_qc analysis_reason aligned to ``ids``; every cell must have one row."""
    p = src.geometry_qc_root / f"donor_id={donor}" / "data.parquet"
    if not p.exists():
        raise ReleaseError(f"donor {donor}: geometry_qc partition missing: {p}")
    gq = pd.read_parquet(p, columns=["object_id", "analysis_reason"])
    gq["object_id"] = gq["object_id"].astype(str)
    if gq["object_id"].duplicated().any():
        raise ReleaseError(f"donor {donor}: duplicate object_id in geometry_qc")
    lost = ~ids.isin(set(gq["object_id"]))
    if lost.any():
        raise ReleaseError(f"donor {donor}: {int(lost.sum())} expression cells have no geometry_qc row")
    aligned = gq.set_index("object_id")["analysis_reason"].reindex(ids.to_numpy())
    return aligned.fillna("").astype(str).reset_index(drop=True)


def _ki67_rank(src: Sources, donor: str, ids: pd.Series) -> np.ndarray:
    """Within-donor rank of the aligned Ki67 over the whole donor, aligned to ``ids``.

    Ranked over every cell in the uniform table (the donor), not over the released subset, so a
    frame-scope release carries the same number a whole-donor release would.  NaN where Ki67 is
    missing for the cell, and NaN for every cell when the donor's table has no Ki67 column.
    """
    p = src.uniform_root / f"donor_id={donor}" / "data.parquet"
    if not p.exists():
        raise ReleaseError(f"donor {donor}: uniform (Step-2 aligned) partition missing: {p}")
    import pyarrow.parquet as pq
    names = set(pq.ParquetFile(p).schema_arrow.names)
    if STATE_RANK_MARKER not in names:
        return np.full(len(ids), np.nan)
    tab = pd.read_parquet(p, columns=["object_id", STATE_RANK_MARKER])
    tab["object_id"] = tab["object_id"].astype(str)
    rank = within_donor_rank(tab[[STATE_RANK_MARKER]])[STATE_RANK_MARKER]
    rank.index = tab["object_id"].to_numpy()
    return rank.reindex(ids.to_numpy()).to_numpy(dtype=float)


def build_donor_calls(donor: str, arm: str, scope: str, src: Sources, dec: dict, verified: set) -> tuple:
    """One release row per expression row, in expression order.  Returns (frame, summary)."""
    declared = set(dec["labels"])
    expr = _read_expression(src, donor)
    ids = expr["object_id"]
    n = len(expr)

    cp = src.calls_path(arm, donor, scope)
    if not cp.exists():
        raise ReleaseError(f"donor {donor}: arm {arm!r} has no primary {scope}-scope calls at {cp}")
    calls = pd.read_parquet(cp)
    calls = records.validate_calls(calls[list(records.CALL_COLUMNS)].copy(), labels=declared)
    calls["object_id"] = calls["object_id"].astype(str)
    extra = ~calls["object_id"].isin(set(ids))
    if extra.any():
        raise ReleaseError(f"donor {donor}: {int(extra.sum())} arm calls are not in the run's expression")
    if scope == "donor" and len(calls) != n:
        raise ReleaseError(f"donor {donor}: whole-donor calls cover {len(calls):,} of {n:,} expression cells; "
                           "a donor-scope release may not leave a cell uncalled")

    out = expr[["object_id"]].merge(calls, on="object_id", how="left")   # keeps expression order
    covered = out["label"].notna().to_numpy()
    label = out["label"].astype(object).where(covered, records.UNAVAILABLE).astype(str)
    reason = out["reason"].fillna("").astype(str).where(covered, REASON_OUTSIDE_FRAME)
    score = pd.to_numeric(out["score"], errors="coerce").astype(float)

    # 1. geometry QC: a property of the cell, applied before anything the arm said
    elig = expr["qc_analysis_eligible"]
    inel = ~elig.fillna(False).astype(bool).to_numpy()
    if inel.any():
        qc_reason = _read_qc_reason(src, donor, ids)
        qr = qc_reason.to_numpy(dtype=object)
        qr = np.where(elig.isna().to_numpy(), "eligibility_unrecorded", qr)
        blank = inel & (qr == "")
        if blank.any():
            raise ReleaseError(f"donor {donor}: {int(blank.sum())} QC-ineligible cells carry no analysis_reason")
        label = label.where(~inel, records.UNAVAILABLE)
        reason = reason.where(~inel, pd.Series([f"{REASON_QC}:{r}" for r in qr], index=reason.index))
    # 2. verification: evidence only where the verdict verified the called label on THIS donor
    is_declared = label.isin(declared).to_numpy()
    evidence = is_declared & label.isin(verified).to_numpy()
    not_verified = is_declared & ~evidence
    reason = reason.where(~not_verified, REASON_NOT_VERIFIED)
    common = label.map(lambda L: dec["labels"][L]["common_lineage"] if L in dec["labels"] else "").astype(str)
    confidence = score.where(label.to_numpy() != records.UNAVAILABLE, np.nan)

    frame = pd.DataFrame({
        "object_id": ids.to_numpy(),
        "lineage_native": label.to_numpy(),
        "common_lineage": common.to_numpy(),
        "evidence_positive": evidence.astype(bool),
        "confidence": confidence.to_numpy(dtype=float),
        "celltype_fine": pd.array([pd.NA] * n, dtype="string"),
        "reason": reason.to_numpy(),
        "ki67_donor_rank": _ki67_rank(src, donor, ids),
        "cell_region": expr["cell_region"].fillna("").astype(str).to_numpy(),
        "qupath_islet_num": expr["islet_num"].fillna("").astype(str).to_numpy(),
    }, columns=list(CALL_COLUMNS))
    if frame["evidence_positive"].isna().any():
        raise ReleaseError(f"donor {donor}: evidence_positive has NA")   # unreachable by construction; kept as the contract
    summary = {
        "n_cells": int(n),
        "n_arm_calls": int(len(calls)),
        "n_called": int(is_declared.sum()),
        "n_evidence_positive": int(evidence.sum()),
        "n_label_not_verified": int(not_verified.sum()),
        "n_qc_ineligible": int(inel.sum()),
        "n_outside_frame": int((~covered).sum()),
        "labels": {k: int(v) for k, v in frame["lineage_native"].value_counts().sort_index().items()},
        "reasons": {k: int(v) for k, v in frame["reason"].value_counts().sort_index().head(40).items()},
        "verified_labels": sorted(verified),
        "ki67_rank_finite": int(np.isfinite(frame["ki67_donor_rank"].to_numpy()).sum()),
        "calls": src.relpath(cp),
        "calls_sha256": _sha256(cp),
    }
    return frame, summary


# --------------------------------------------------------------------------- the release
def write_release(arm: str, donors, *, scope: str = "donor", out: pathlib.Path | None = None,
                  src: Sources | None = None, today: str | None = None) -> dict:
    """Write release.json, label_space.json and calls/ under ``out``; return the release payload.

    Fails before writing anything if an input is missing or a provenance string is retired,
    and writes calls last so a half-written release never carries a valid release.json.
    """
    src = src or Sources()
    out = pathlib.Path(out) if out is not None else pins.OUT / "release"
    donors = [str(d) for d in donors]
    if scope not in ("donor", "frame"):
        raise ReleaseError(f"scope must be donor or frame, got {scope!r}")
    if not donors or len(set(donors)) != len(donors):
        raise ReleaseError("donors must be a non-empty, duplicate-free list")
    blocked = sorted(set(donors) & set(DONOR_EXCLUSIONS))
    if blocked:
        raise ReleaseError(f"release names centrally excluded donor(s) {blocked}; cohort.py owns that decision")
    for tok in RETIRED_INPUT_TOKENS:
        if tok.lower() in arm.lower():
            raise ReleaseError(f"arm name {arm!r} names retired input token {tok!r}")

    dec = json.loads(pathlib.Path(src.declarations_path).read_text())
    declared = set(dec["labels"])
    manifest = json.loads(pathlib.Path(src.expression_manifest).read_text())
    got_cid = str(manifest.get("content_id", ""))
    if got_cid != src.expression_content_id:
        raise ReleaseError(f"expression manifest content_id {got_cid!r} != pinned {src.expression_content_id!r}")
    row_counts = {str(d.get("donor_id")): int(d.get("row_count", -1))
                  for d in (manifest.get("output", {}) or {}).get("donors", []) if isinstance(d, dict)}
    for p, what in ((src.shifts_path, "shifts.json"), (src.verdict_summary_path(), "verdict summary")):
        if not pathlib.Path(p).exists():
            raise ReleaseError(f"{what} missing: {p}")

    verdicts, verified = {}, {}
    for d in donors:
        vp = src.verdict_path(arm, d, scope)
        if not vp.exists():
            raise ReleaseError(f"donor {d}: no verdict for arm {arm!r} at {vp}; evidence cannot be declared without it")
        vj = json.loads(vp.read_text())
        verified[d] = verified_labels(vj, arm=arm, donor=d, declared=declared)
        verdicts[d] = {"path": src.relpath(vp), "sha256": _sha256(vp),
                       "scope": str(vj.get("scope", "frame")) if isinstance(vj, dict) else "frame"}

    release_id = release_id_for(arm, today)
    lane = lane_for(arm)
    head = pins.git_head()
    dec_sha = _sha256(pathlib.Path(src.declarations_path))
    space = label_space_payload(dec, arm, declarations_sha256=dec_sha)

    frames, per_donor = {}, {}
    for d in donors:
        frames[d], per_donor[d] = build_donor_calls(d, arm, scope, src, dec, verified[d])
        want = row_counts.get(d)
        if want is not None and want >= 0 and want != per_donor[d]["n_cells"]:
            raise ReleaseError(f"donor {d}: expression has {per_donor[d]['n_cells']:,} rows but the manifest "
                               f"records {want:,}; the run moved under the pins")
        per_donor[d]["verdict"] = verdicts[d]
        undeclared = set(frames[d]["lineage_native"]) - declared - {records.AMBIGUOUS, records.UNAVAILABLE}
        if undeclared:
            raise ReleaseError(f"donor {d}: calls carry labels absent from label_space.json: {sorted(undeclared)}")

    payload = {
        "schema_version": SCHEMA_VERSION,
        "release_id": release_id,
        "modality": MODALITY,
        "producer": {"lane": lane, "commit": head["commit"], "dirty": bool(head["dirty"])},
        "label_space": "label_space.json",
        "calls_dir": "calls",
        "cell_key": CELL_KEY,
        "donors": donors,
        "inputs": [
            {"kind": "expression", "run_id": src.run_id, "content_id": src.expression_content_id,
             "output_content_sha256": src.expression_output_sha256, "path": src.relpath(src.expression_root)},
            {"kind": "uniform_aligned", "path": src.relpath(src.uniform_root), "sha256": _sha256(src.shifts_path),
             "shifts": src.relpath(src.shifts_path)},
            {"kind": "declarations", "path": src.relpath(src.declarations_path), "sha256": dec_sha,
             "method_version": dec["method_version"]},
        ],
        "columns": dict(COLUMNS),
        "latch": LATCH,
        "stamp": STAMP,
        "verdict_summary": src.relpath(src.verdict_summary_path()),
        "verdict_summary_sha256": _sha256(src.verdict_summary_path()),
        "scope": scope,
        "roi": ROI,
        "arm": {"name": arm, "scope": scope, "per_donor": per_donor},
        "reason_vocabulary": dict(REASON_VOCABULARY),
        "call_columns": list(CALL_COLUMNS),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if scope == "frame":
        payload["development_note"] = ("frame-scope release: cells outside the 16-core screen frame are "
                                       f"Unavailable({REASON_OUTSIDE_FRAME}); this is not a whole-donor call")
    scan_provenance(payload, what="release.json")
    scan_provenance(space, what="label_space.json")

    out.mkdir(parents=True, exist_ok=True)
    stale = out / "release.json"
    if stale.exists():
        stale.unlink()                       # no valid contract while calls are being rewritten
    (out / "label_space.json").write_text(json.dumps(space, indent=1) + "\n")
    for d, frame in frames.items():
        part = out / "calls" / f"donor_id={d}" / f"roi={ROI}"
        part.mkdir(parents=True, exist_ok=True)
        for old in part.glob("*.parquet"):
            old.unlink()
        frame.to_parquet(part / "data.parquet", index=False)
    (out / "release.json").write_text(json.dumps(payload, indent=1) + "\n")
    return payload


def verify_release(out: pathlib.Path, cfg=None):
    """Load what was written through phenospace itself, the way integration will."""
    from phenocycler.integration.phenospace import load_release
    if cfg is None:
        from phenocycler.config import load_config
        cfg = load_config(pins.PKG / "config.ini")
    cfg.phenocycler_identity_release = pathlib.Path(out).resolve()
    return load_release(cfg, MODALITY)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="A.4a development release writer for one broadtype arm")
    ap.add_argument("--arm", required=True, help="the producer arm (verdict.py's A.3 choice)")
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS))
    ap.add_argument("--scope", choices=("donor", "frame"), default="donor")
    ap.add_argument("--out", type=pathlib.Path, default=pins.OUT / "release")
    ap.add_argument("--arms-root", type=pathlib.Path, default=None,
                    help="override data/development/broadtype/arms (smoke tests against scratch)")
    ap.add_argument("--verdict-root", type=pathlib.Path, default=None,
                    help="override data/development/broadtype/verdict (smoke tests against scratch)")
    ap.add_argument("--today", default=None, help="YYYYMMDD for release_id; default UTC today")
    a = ap.parse_args(argv)
    pins.assert_pins()
    src = Sources()
    if a.arms_root is not None:
        src.arms_root = a.arms_root
    if a.verdict_root is not None:
        src.verdict_root = a.verdict_root
    t0 = time.time()
    payload = write_release(a.arm, a.donors, scope=a.scope, out=a.out, src=src, today=a.today)
    rel = verify_release(a.out)
    print(f"release {rel.release_id} lane={rel.producer_lane} donors={list(rel.donors)} "
          f"run={rel.run_id} content_id={rel.expression_content_id[:12]}… -> {a.out} ({time.time() - t0:.1f}s)")
    for d, s in payload["arm"]["per_donor"].items():
        print(f"  donor {d}: {s['n_cells']:,} cells, {s['n_called']:,} called, {s['n_evidence_positive']:,} evidence_positive, "
              f"{s['n_qc_ineligible']:,} qc_ineligible, {s['n_outside_frame']:,} outside_frame, "
              f"verified={s['verified_labels']}")
        print(f"     labels: {s['labels']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
