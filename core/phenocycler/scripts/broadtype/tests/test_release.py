"""The A.4a writer must produce what phenospace accepts, and refuse what it must refuse.

Everything here is synthetic (no cohort data): a tmp run tree in the shape
``tests/integration/conftest.py`` uses -- ``manifests/expression.json`` carrying ``content_id``,
a tiny ``expression`` partition -- plus tiny geometry_qc / uniform / shifts / arm-calls /
verdict files.  The load-bearing checks: ``load_release`` succeeds on the output; every
expression object_id is present exactly once, in expression order; ``evidence_positive`` has
no NA and is True exactly where a VERIFIED label was called on a QC-eligible cell; geometry-QC
ineligibility overrides the arm; a frame-scope release marks uncovered cells
``Unavailable(outside_screen_frame)``; ``attach_release_identity`` joins a matching geometry
frame; and a retired token anywhere in the provenance is refused before phenospace sees it.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import release  # noqa: E402
from records import AMBIGUOUS, UNAVAILABLE  # noqa: E402

from phenocycler.config import PipelineConfig  # noqa: E402
from phenocycler.integration.contract import coerce_cell_table, contract_profile, has_identity  # noqa: E402
from phenocycler.integration.phenospace import (  # noqa: E402
    PhenotypeSpaceError, attach_release_identity, load_phenotype_space, load_release,
)

ARM = "synthetic-arm"
DONOR = "6436"
RUN = "run-synthetic"
CONTENT_ID = "expression-content-synthetic"
IDS = [f"u{i:02d}" for i in range(12)]
DECLARATIONS = json.loads((HERE / "declarations.json").read_text())

#: What the synthetic arm says, per cell (label, score, reason).  Cell u04 is QC-ineligible
#: but the arm still called it; u09 has no Ki67 measurement.
ARM_CALLS = {
    "u00": ("Endocrine", 0.9, ""),
    "u01": ("Endocrine", 0.8, ""),
    "u02": ("Vascular", 0.7, ""),          # declared, NOT verified on this donor
    "u03": ("T_NK", 0.6, ""),              # declared, verified
    "u04": ("Endocrine", 0.95, ""),        # QC-ineligible -> Unavailable regardless
    "u05": (AMBIGUOUS, 0.4, "multi_lineage"),
    "u06": (UNAVAILABLE, np.nan, "arm_no_output"),
    "u07": ("Myeloid", 0.5, ""),           # declared, not verified
    "u08": ("Stromal", 0.55, ""),          # declared, verified
    "u09": ("Endocrine", 0.85, ""),
    "u10": ("B_Plasma", 0.65, ""),         # not verified
    "u11": ("Epithelial_nonendocrine", 0.75, ""),   # verified
}
VERIFIED = {"Endocrine", "T_NK", "Stromal", "Epithelial_nonendocrine"}
QC_INELIGIBLE = {"u04": "cell_area_below_min"}


def _write_json(path: pathlib.Path, payload) -> pathlib.Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return path


def build_sources(tmp_path: pathlib.Path, *, scope: str = "donor", calls_subset=None,
                  arm_labels: dict | None = None, verdict: dict | None = None,
                  arms_root_name: str = "arms") -> tuple:
    """A complete synthetic input tree.  Returns (Sources, PipelineConfig, data_dir)."""
    data = tmp_path / "data"
    run = data / "runs" / RUN
    _write_json(run / "manifests" / "expression.json", {
        "content_id": CONTENT_ID,
        "output": {"content_sha256": "synthetic-output-sha", "donors": [{"donor_id": DONOR, "row_count": len(IDS)}]},
    })
    rng = np.random.default_rng(0)
    ki67 = rng.uniform(0, 100, len(IDS))
    ki67[IDS.index("u09")] = np.nan
    expr = pd.DataFrame({
        "donor_id": DONOR, "object_id": IDS,
        "qc_analysis_eligible": [i not in QC_INELIGIBLE for i in IDS],
        "X_centroid": np.arange(len(IDS), dtype=float), "Y_centroid": np.arange(len(IDS), dtype=float) * 2,
        "cell_region": ["core" if k % 3 == 0 else "tissue" for k in range(len(IDS))],
        "islet_num": ["7" if k % 3 == 0 else "" for k in range(len(IDS))],
        "INS": rng.uniform(0, 10, len(IDS)), "CD31": rng.uniform(0, 10, len(IDS)), "Ki67": ki67,
    })
    p = run / "expression" / f"donor_id={DONOR}" / "data.parquet"
    p.parent.mkdir(parents=True)
    expr.to_parquet(p, index=False)
    # geometry_qc deliberately in a DIFFERENT row order: the join must be by object_id
    gq = pd.DataFrame({
        "object_id": IDS[::-1],
        "analysis_eligible": [i not in QC_INELIGIBLE for i in IDS[::-1]],
        "analysis_reason": [QC_INELIGIBLE.get(i, "") for i in IDS[::-1]],
    })
    p = run / "geometry_qc" / f"donor_id={DONOR}" / "data.parquet"
    p.parent.mkdir(parents=True)
    gq.to_parquet(p, index=False)
    dev = tmp_path / "development"
    uni = expr[["donor_id", "object_id", "qc_analysis_eligible", "INS", "CD31", "Ki67"]].copy()
    uni["Ki67"] = uni["Ki67"] + 3.0                    # a shift; ranks are shift-invariant
    p = dev / "celltype" / "uniform" / f"donor_id={DONOR}" / "data.parquet"
    p.parent.mkdir(parents=True)
    uni.to_parquet(p, index=False)
    shifts = _write_json(dev / "celltype" / "shifts.json", {"uniform_commit": "deadbeef", "shifts": {"Ki67": {}}})

    labels = arm_labels if arm_labels is not None else ARM_CALLS
    ids = list(labels) if calls_subset is None else [i for i in labels if i in calls_subset]
    calls = pd.DataFrame({
        "object_id": ids,
        "label": [labels[i][0] for i in ids],
        "score": [labels[i][1] for i in ids],
        "reason": [labels[i][2] for i in ids],
    })
    arms_root = dev / "broadtype" / arms_root_name
    cp = arms_root / ARM / f"donor_id={DONOR}" / scope / "calls__primary.parquet"
    cp.parent.mkdir(parents=True)
    calls.to_parquet(cp, index=False)

    verdict_root = dev / "broadtype" / "verdict"
    payload = verdict if verdict is not None else {
        "arm": ARM, "donor_id": DONOR,
        "labels": {L: {"verified": L in VERIFIED, "spatial": {"clear": L in VERIFIED}} for L in DECLARATIONS["labels"]},
    }
    _write_json(verdict_root / ARM / f"donor_id={DONOR}.json", payload)
    _write_json(verdict_root / "summary.json", {"arms": {ARM: {"survives": True}}})

    src = release.Sources(
        arms_root=arms_root, verdict_root=verdict_root,
        expression_root=run / "expression", expression_manifest=run / "manifests" / "expression.json",
        geometry_qc_root=run / "geometry_qc", uniform_root=dev / "celltype" / "uniform", shifts_path=shifts,
        declarations_path=HERE / "declarations.json",
        run_id=RUN, expression_content_id=CONTENT_ID, expression_output_sha256="synthetic-output-sha",
        repo_root=tmp_path,
    )
    cfg = PipelineConfig(data_dir=data)
    return src, cfg, data


def _geometry_frame(ids):
    return coerce_cell_table(
        pd.DataFrame({"cell_id": list(ids), "x_um": np.arange(len(ids), dtype=float),
                      "y_um": np.arange(len(ids), dtype=float) * 2}),
        modality="phenocycler", donor_id=DONOR, roi="panc", profile="geometry",
    )


def test_whole_donor_release_loads_and_declares_complete_evidence(tmp_path):
    src, cfg, _ = build_sources(tmp_path)
    out = tmp_path / "release"
    payload = release.write_release(ARM, [DONOR], scope="donor", out=out, src=src, today="20260905")

    rel = release.verify_release(out, cfg)
    assert rel.release_id == f"broadtype-{ARM}-dev-20260905"
    assert rel.producer_lane == f"broadtype-{ARM}" == rel.space.producer_lane
    assert rel.donors == (DONOR,)
    assert rel.run_id == RUN and rel.expression_content_id == CONTENT_ID
    assert rel.cell_key == "object_id"
    assert rel.label_column == "lineage_native" and rel.evidence_column == "evidence_positive"
    assert payload["latch"] == "NO_VERIFIED_TYPING (unchanged)" and payload["stamp"] == "development"
    assert payload["producer"]["commit"] and isinstance(payload["producer"]["dirty"], bool)
    kinds = [i["kind"] for i in payload["inputs"]]
    assert kinds == ["expression", "uniform_aligned", "declarations"]
    assert payload["inputs"][0]["output_content_sha256"] == "synthetic-output-sha"
    assert len(payload["inputs"][1]["sha256"]) == 64 and len(payload["inputs"][2]["sha256"]) == 64

    calls = pd.read_parquet(out / "calls" / f"donor_id={DONOR}" / "roi=panc" / "data.parquet")
    assert list(calls.columns) == list(release.CALL_COLUMNS)
    # every expression object_id exactly once, in expression order
    assert list(calls["object_id"]) == IDS
    # evidence: non-null boolean, True iff verified label on a QC-eligible cell
    assert calls["evidence_positive"].dtype == bool and not calls["evidence_positive"].isna().any()
    by = calls.set_index("object_id")
    for i, (lab, score, _reason) in ARM_CALLS.items():
        row = by.loc[i]
        if i in QC_INELIGIBLE:
            assert row["lineage_native"] == UNAVAILABLE
            assert row["reason"] == f"qc_analysis_ineligible:{QC_INELIGIBLE[i]}"
            assert not row["evidence_positive"] and np.isnan(row["confidence"])
        elif lab in VERIFIED:
            assert row["lineage_native"] == lab and row["evidence_positive"] and row["reason"] == ""
            assert row["confidence"] == pytest.approx(score)
        elif lab in DECLARATIONS["labels"]:
            assert row["lineage_native"] == lab and not row["evidence_positive"]
            assert row["reason"] == "label_not_verified" and row["confidence"] == pytest.approx(score)
        elif lab == AMBIGUOUS:
            assert row["lineage_native"] == AMBIGUOUS and row["reason"] == "multi_lineage"
            assert not row["evidence_positive"] and row["confidence"] == pytest.approx(score)
        else:
            assert row["lineage_native"] == UNAVAILABLE and row["reason"] == "arm_no_output"
            assert not row["evidence_positive"] and np.isnan(row["confidence"])
    # common_lineage comes from the declaration, '' for the non-mapped labels
    for L, spec in DECLARATIONS["labels"].items():
        sel = calls["lineage_native"] == L
        assert (calls.loc[sel, "common_lineage"] == spec["common_lineage"]).all()
    assert (calls.loc[calls["lineage_native"].isin([AMBIGUOUS, UNAVAILABLE]), "common_lineage"] == "").all()
    # celltype_fine NA; audit columns carried; Ki67 rank in (0,1), NaN where the measurement is missing
    assert calls["celltype_fine"].isna().all()
    assert by.loc["u00", "cell_region"] == "core" and by.loc["u00", "qupath_islet_num"] == "7"
    assert by.loc["u01", "qupath_islet_num"] == ""
    r = calls["ki67_donor_rank"]
    assert np.isnan(by.loc["u09", "ki67_donor_rank"]) and r.notna().sum() == len(IDS) - 1
    assert ((r.dropna() > 0) & (r.dropna() < 1)).all()
    assert payload["arm"]["per_donor"][DONOR]["n_evidence_positive"] == int(calls["evidence_positive"].sum()) == 6

    space = load_phenotype_space(out / "label_space.json")
    assert set(space.rules) == set(DECLARATIONS["labels"]) | {AMBIGUOUS, UNAVAILABLE}
    assert space.common_for("Epithelial_nonendocrine") == ("Exocrine", "mapped")
    assert space.common_for(AMBIGUOUS) == ("", "ambiguous")
    assert space.common_for(UNAVAILABLE) == ("", "unavailable")
    assert space.rules["Endocrine"].reason == DECLARATIONS["labels"]["Endocrine"]["reason"]


def test_attach_release_identity_joins_a_matching_geometry_frame(tmp_path):
    src, cfg, _ = build_sources(tmp_path)
    out = tmp_path / "release"
    release.write_release(ARM, [DONOR], scope="donor", out=out, src=src)
    rel = release.verify_release(out, cfg)
    attached = attach_release_identity(_geometry_frame(IDS[::-1]), rel, DONOR)   # geometry in another order
    assert contract_profile(attached) == "identity" and has_identity(attached)
    assert attached.attrs["identity_release_id"] == rel.release_id
    got = attached.set_index("cell_id")
    assert got.loc["u00", "lineage_common"] == "Endocrine" and got.loc["u11", "lineage_common"] == "Exocrine"
    assert got.loc["u05", "common_resolution"] == "ambiguous" and got.loc["u04", "common_resolution"] == "unavailable"
    assert not attached["evidence_positive"].isna().any()
    assert bool(got.loc["u03", "evidence_positive"]) and not bool(got.loc["u02", "evidence_positive"])
    # a geometry frame with one extra cell is a UUID mismatch, not a silent fill
    with pytest.raises(PhenotypeSpaceError, match="UUID mismatch"):
        attach_release_identity(_geometry_frame(IDS + ["u99"]), rel, DONOR)


def test_frame_scope_marks_uncovered_cells_outside_screen_frame(tmp_path):
    covered = {"u00", "u02", "u04", "u05", "u06"}
    src, cfg, _ = build_sources(tmp_path, scope="frame", calls_subset=covered)
    out = tmp_path / "release"
    payload = release.write_release(ARM, [DONOR], scope="frame", out=out, src=src)
    assert payload["scope"] == "frame" and "development_note" in payload
    calls = pd.read_parquet(out / "calls" / f"donor_id={DONOR}" / "roi=panc" / "data.parquet")
    assert list(calls["object_id"]) == IDS
    by = calls.set_index("object_id")
    outside = by.loc[[i for i in IDS if i not in covered]]
    assert (outside["lineage_native"] == UNAVAILABLE).all() and (outside["reason"] == "outside_screen_frame").all()
    assert not outside["evidence_positive"].any() and outside["confidence"].isna().all()
    assert by.loc["u00", "evidence_positive"] and by.loc["u04", "reason"].startswith("qc_analysis_ineligible:")
    assert not calls["evidence_positive"].isna().any()
    assert release.verify_release(out, cfg).donors == (DONOR,)
    # the same partial calls are NOT acceptable as a whole-donor release
    src2, _, _ = build_sources(tmp_path / "again", scope="donor", calls_subset=covered)
    with pytest.raises(release.ReleaseError, match="may not leave a cell uncalled"):
        release.write_release(ARM, [DONOR], scope="donor", out=tmp_path / "again" / "release", src=src2)


def test_retired_token_in_provenance_is_refused_before_writing(tmp_path):
    src, _, _ = build_sources(tmp_path, arms_root_name="restore_gated_redsea")
    out = tmp_path / "release"
    with pytest.raises(release.ReleaseError, match="retired"):
        release.write_release(ARM, [DONOR], scope="donor", out=out, src=src)
    assert not (out / "release.json").exists()
    src2, _, _ = build_sources(tmp_path / "b", arms_root_name="assignments")
    with pytest.raises(release.ReleaseError, match="assignments"):
        release.write_release(ARM, [DONOR], scope="donor", out=tmp_path / "b" / "release", src=src2)
    # and phenospace itself would have refused the same string, so the two checks agree
    payload = {"lane": "x", "inputs": [{"path": "data/restore_gated_redsea/calls.parquet"}]}
    with pytest.raises(release.ReleaseError):
        release.scan_provenance(payload, what="probe")


def test_missing_or_malformed_verdict_fails_closed(tmp_path):
    src, _, _ = build_sources(tmp_path)
    (src.verdict_root / ARM / f"donor_id={DONOR}.json").unlink()
    with pytest.raises(release.ReleaseError, match="no verdict"):
        release.write_release(ARM, [DONOR], out=tmp_path / "release", src=src)
    declared = set(DECLARATIONS["labels"])
    # list form is accepted; a non-boolean 'verified', another arm, or an undeclared label is not
    assert release.verified_labels({"labels": [{"label": "Endocrine", "verified": True},
                                               {"label": "Vascular", "verified": False}]},
                                   arm=ARM, donor=DONOR, declared=declared) == {"Endocrine"}
    with pytest.raises(release.ReleaseError, match="non-boolean"):
        release.verified_labels({"labels": {"Endocrine": {"verified": "yes"}}}, arm=ARM, donor=DONOR, declared=declared)
    with pytest.raises(release.ReleaseError, match="names arm"):
        release.verified_labels({"arm": "other", "labels": {}}, arm=ARM, donor=DONOR, declared=declared)
    with pytest.raises(release.ReleaseError, match="undeclared"):
        release.verified_labels({"labels": {"Neural": {"verified": True}}}, arm=ARM, donor=DONOR, declared=declared)
    with pytest.raises(release.ReleaseError, match="no 'labels'"):
        release.verified_labels({"verified": True}, arm=ARM, donor=DONOR, declared=declared)


def test_scope_suffixed_verdict_is_preferred_and_status_form_accepted(tmp_path):
    """verdict.py writes frame verdicts as donor_id=<d>.json and other scopes with a __<scope>
    suffix; a donor-scope release takes the donor-scope verdict when it exists and records
    which one it read, and falls back to the screen's frame verdict otherwise."""
    src, cfg, _ = build_sources(tmp_path)
    frame_v = src.verdict_root / ARM / f"donor_id={DONOR}.json"
    donor_v = src.verdict_root / ARM / f"donor_id={DONOR}__donor.json"
    assert src.verdict_path(ARM, DONOR, "donor") == frame_v          # no scoped file -> the frame verdict
    # the classified form verdict.py itself accepts: status strings instead of a boolean
    _write_json(donor_v, {"arm": ARM, "donor_id": DONOR, "scope": "donor",
                          "labels": {L: {"status": "verified" if L == "Myeloid" else "failed", "blocking_reasons": []}
                                     for L in DECLARATIONS["labels"]}})
    assert src.verdict_path(ARM, DONOR, "donor") == donor_v
    assert src.verdict_path(ARM, DONOR, "frame") == frame_v
    payload = release.write_release(ARM, [DONOR], scope="donor", out=tmp_path / "release", src=src)
    rec = payload["arm"]["per_donor"][DONOR]
    assert rec["verified_labels"] == ["Myeloid"] and rec["verdict"]["scope"] == "donor"
    assert rec["verdict"]["path"].endswith("donor_id=6436__donor.json")
    calls = pd.read_parquet(tmp_path / "release" / "calls" / f"donor_id={DONOR}" / "roi=panc" / "data.parquet")
    by = calls.set_index("object_id")
    assert bool(by.loc["u07", "evidence_positive"]) and not bool(by.loc["u00", "evidence_positive"])
    assert by.loc["u00", "reason"] == "label_not_verified"
    assert release.verify_release(tmp_path / "release", cfg).release_id.startswith(f"broadtype-{ARM}-dev-")
    declared = set(DECLARATIONS["labels"])
    with pytest.raises(release.ReleaseError, match="contradicts"):
        release.verified_labels({"labels": {"Endocrine": {"verified": True, "status": "failed"}}},
                                arm=ARM, donor=DONOR, declared=declared)


def test_undeclared_label_excluded_donor_and_moved_run_are_refused(tmp_path):
    bad = dict(ARM_CALLS)
    bad["u00"] = ("Neural", 0.9, "")
    src, _, _ = build_sources(tmp_path, arm_labels=bad)
    with pytest.raises(ValueError, match="undeclared labels"):
        release.write_release(ARM, [DONOR], out=tmp_path / "release", src=src)
    src2, _, _ = build_sources(tmp_path / "b")
    with pytest.raises(release.ReleaseError, match="excluded donor"):
        release.write_release(ARM, [DONOR, "6457"], out=tmp_path / "b" / "release", src=src2)
    src3, cfg3, _ = build_sources(tmp_path / "c")
    src3.expression_content_id = "something-else"
    with pytest.raises(release.ReleaseError, match="content_id"):
        release.write_release(ARM, [DONOR], out=tmp_path / "c" / "release", src=src3)
    # phenospace catches the same drift from its side when the manifest moves after writing
    src4, cfg4, data4 = build_sources(tmp_path / "d")
    out4 = tmp_path / "d" / "release"
    release.write_release(ARM, [DONOR], out=out4, src=src4)
    (data4 / "runs" / RUN / "manifests" / "expression.json").write_text(json.dumps({"content_id": "moved"}))
    with pytest.raises(PhenotypeSpaceError, match="content id"):
        load_release(PipelineConfig(data_dir=data4, phenocycler_identity_release=out4), "phenocycler")
