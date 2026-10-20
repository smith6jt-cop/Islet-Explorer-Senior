"""Arm E ``nimbus-rules`` without Nimbus: the label rule on a planted binary table, the plane
map's OME assertion, the FOV layout Nimbus expects, and (skip-if-missing) the worker's pins.

No cohort data.  The rule is the part of the arm that can be wrong in a way a run would not
reveal -- a negative that fails to block, an ablated marker that still votes, a missing cell
that becomes a label -- so each of those is planted and checked here.  Nimbus itself is only
import-checked, because a forward pass needs the checkpoint and a GPU.
"""
from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

BROADTYPE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BROADTYPE))
sys.path.insert(0, str(BROADTYPE / "arms"))

import records  # noqa: E402
import nimbus_rules as arm  # noqa: E402

LABELS = json.loads((BROADTYPE / "declarations.json").read_text())["labels"]
MARKERS = ["INS", "GCG", "SST", "E_cadherin", "Ker8_18", "CD31", "Vimentin", "SMA", "CD3e", "CD8", "CD20", "CD68"]


def _conf(rows: dict, markers=MARKERS) -> tuple:
    """rows: object_id -> {marker: confidence}; unmentioned markers are 0.05; None row = not scored."""
    ids = np.array(list(rows), dtype=object)
    C = pd.DataFrame(np.full((len(ids), len(markers)), 0.05), columns=markers)
    for i, (oid, vals) in enumerate(rows.items()):
        if vals is None:
            C.iloc[i] = np.nan
            continue
        for m, v in vals.items():
            C.loc[i, m] = v
    return ids, C


def test_rule_recovers_planted_labels_and_reasons():
    ids, C = _conf({
        "beta":      {"INS": 0.91},
        "alpha":     {"GCG": 0.77, "INS": 0.10},
        "exocrine":  {"E_cadherin": 0.80, "Ker8_18": 0.66},
        "endo_vasc": {"INS": 0.70, "CD31": 0.60},                  # two candidates -> multi_lineage
        "t_mac":     {"CD3e": 0.90, "CD68": 0.85},                 # each negates the other
        "nothing":   {},                                           # nothing positive
        "t":         {"CD3e": 0.55, "CD8": 0.51},
        "stromal_v": {"Vimentin": 0.62, "CD31": 0.58},             # CD31 negates Stromal; Vascular stands
        "gone":      None,                                         # not in any Nimbus row
    })
    calls, info = arm.calls_from_confidence(ids, C, MARKERS, LABELS, 0.5)
    records.validate_calls(calls, expected_ids=ids, labels=set(LABELS))
    got = calls.set_index("object_id")
    assert got.loc["beta", "label"] == "Endocrine" and got.loc["beta", "score"] == pytest.approx(0.91)
    assert got.loc["alpha", "label"] == "Endocrine" and got.loc["alpha", "score"] == pytest.approx(0.77)
    assert got.loc["exocrine", "label"] == "Epithelial_nonendocrine" and got.loc["exocrine", "score"] == pytest.approx(0.80)
    assert got.loc["t", "label"] == "T_NK"
    assert got.loc["stromal_v", "label"] == "Vascular"
    assert (got.loc["endo_vasc", ["label", "reason"]] == ["Ambiguous", "multi_lineage"]).all()
    assert (got.loc["t_mac", ["label", "reason"]] == ["Ambiguous", "all_candidates_negated"]).all()
    assert (got.loc["nothing", ["label", "reason"]] == ["Ambiguous", "no_positive_marker"]).all()
    assert (got.loc["gone", ["label", "reason"]] == ["Unavailable", "not_in_fov"]).all()
    assert np.isnan(got.loc["gone", "score"])
    assert info["n_not_in_fov"] == 1 and info["n_scored"] == 8
    assert info["positivity_rate"]["INS"] == pytest.approx(2 / 8)      # beta and endo_vasc; alpha has INS 0.10
    assert list(calls.columns) == list(records.CALL_COLUMNS)


def test_threshold_is_the_only_cut_and_is_inclusive():
    ids, C = _conf({"edge": {"CD31": 0.5}, "below": {"CD31": 0.4999}})
    calls, _ = arm.calls_from_confidence(ids, C, MARKERS, LABELS, 0.5)
    got = calls.set_index("object_id")
    assert got.loc["edge", "label"] == "Vascular"
    assert got.loc["below", "reason"] == "no_positive_marker"


def test_ablated_marker_neither_votes_nor_negates():
    # CD3e is in the table but NOT in the identity input: it may not vote for T_NK and may not
    # block Myeloid.  That is the transparency the held-out gate relies on.
    ids, C = _conf({"t_only_cd3": {"CD3e": 0.95}, "mac_with_cd3": {"CD68": 0.9, "CD3e": 0.9}})
    markers = [m for m in MARKERS if m != "CD3e"]
    calls, info = arm.calls_from_confidence(ids, C, markers, LABELS, 0.5)
    got = calls.set_index("object_id")
    assert got.loc["t_only_cd3", "reason"] == "no_positive_marker"
    assert got.loc["mac_with_cd3", "label"] == "Myeloid"
    assert "CD3e" not in info["positivity_rate"]


def test_a_label_with_no_marker_in_the_input_is_unreachable_not_defaulted():
    markers = ["INS", "CD31"]                       # only Endocrine and Vascular can be reached
    ids, C = _conf({"x": {"INS": 0.9}}, markers=markers)
    calls, info = arm.calls_from_confidence(ids, C, markers, LABELS, 0.5)
    assert calls.loc[0, "label"] == "Endocrine"
    assert set(info["unreachable_labels"]) == set(LABELS) - {"Endocrine", "Vascular"}


def test_partial_confidence_row_is_unavailable_not_negative():
    ids, C = _conf({"half": {"INS": 0.9}})
    C.loc[0, "CD31"] = np.nan
    calls, _ = arm.calls_from_confidence(ids, C, MARKERS, LABELS, 0.5)
    assert (calls.loc[0, ["label", "reason"]] == ["Unavailable", "marker_confidence_missing"]).all()


def test_plane_index_asserts_ome_name_agreement():
    cm = [{"channel_index": 0, "channel_name": "DAPI", "marker": "DAPI"},
          {"channel_index": 1, "channel_name": "E-cadherin", "marker": "E_cadherin"},
          {"channel_index": 2, "channel_name": "CD3e", "marker": "CD3e"}]
    ome = ["DAPI", "E-cadherin", "CD3e"]
    assert arm.plane_index(cm, ome, ["CD3e", "E_cadherin"]) == {"CD3e": 2, "E_cadherin": 1}
    with pytest.raises(RuntimeError, match="OME name mismatch"):
        arm.plane_index(cm, ["DAPI", "CD3e", "E-cadherin"], ["CD3e"])
    with pytest.raises(KeyError):
        arm.plane_index(cm, ome, ["CD68"])


def test_fov_layout_is_what_nimbus_reads(tmp_path):
    """Folder per FOV, ``<marker>.tiff`` per plane with the window's exact uint16 values, and a
    label TIFF whose int32 labels round-trip -- the contract ``MultiplexDataset`` reads."""
    import tifffile
    rng = np.random.default_rng(0)
    planes = {m: rng.integers(0, 4000, size=(96, 96), dtype=np.uint16) for m in ("CD3e", "INS", "SMA")}
    src = tmp_path / "planes.tif"
    with tifffile.TiffWriter(src) as tw:
        for m, a in planes.items():
            tw.write(a, description=f"<Biomarker>{m}</Biomarker>", contiguous=False, photometric="minisblack")
    tf = tifffile.TiffFile(src)
    pages = list(tf.pages)
    pidx = {"CD3e": 0, "INS": 1, "SMA": 2}
    window = {"y0": 10, "x0": 20, "size": 32}
    mask = np.zeros((32, 32), np.int32)
    mask[2:10, 2:10] = 1_938_318          # a GLOBAL label, far beyond the window's own count
    mask[20:30, 20:30] = 7
    fov = arm.write_fov(pages, pidx, window, mask, tmp_path / "fovs" / "w03", tmp_path / "masks" / "w03_whole_cell.tiff")
    tf.close()
    assert sorted(p.name for p in (tmp_path / "fovs" / "w03").iterdir()) == ["CD3e.tiff", "INS.tiff", "SMA.tiff"]
    for m, a in planes.items():
        got = tifffile.imread(tmp_path / "fovs" / "w03" / f"{m}.tiff")
        assert got.dtype == np.uint16 and np.array_equal(got, a[10:42, 20:52])
    back = tifffile.imread(fov["mask"])
    assert back.dtype == np.int32 and np.array_equal(back, mask)
    assert fov["name"] == "w03" and fov["n_labels"] == 2
    with pytest.raises(ValueError, match="clipped"):
        arm.write_fov(list(tifffile.TiffFile(src).pages), pidx, {"y0": 80, "x0": 0, "size": 32}, mask,
                      tmp_path / "fovs" / "w04", tmp_path / "masks" / "w04_whole_cell.tiff")


@pytest.mark.skipif(not arm.NIMBUS_PY.exists(), reason="external/.venv-nimbus is not installed here")
def test_worker_pins_match_the_environment():
    proc = subprocess.run([str(arm.NIMBUS_PY), str(arm.WORKER), "--check"], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    facts = json.loads(proc.stdout[proc.stdout.index("{"):])
    assert facts["nimbus_commit"] == arm.PIN["nimbus_commit"]
    assert facts["checkpoint_sha256"] == arm.PIN["checkpoint_sha256"]
    assert facts["torch"] == arm.PIN["torch"] and facts["numpy"] == arm.PIN["numpy"]


def test_params_and_pins_carry_no_frequency_terms():
    from criteria import FORBIDDEN_KEY_TOKENS
    def walk(d, path=""):
        for k, v in d.items():
            low = str(k).lower()
            assert not any(t in low for t in FORBIDDEN_KEY_TOKENS), f"{path}{k} looks like an expected-frequency term"
            if isinstance(v, dict):
                walk(v, f"{path}{k}.")
    walk(arm.PARAMS)
    assert arm.PARAMS["positivity_threshold"] == 0.5
    assert arm.PIN["nimbus_commit"].startswith("241281ef")
