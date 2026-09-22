"""Arm F2 ``kronos2-coral`` contract checks without cohort data, a GPU, or the model: the marker
mapping covers every declared identity marker, the novel-marker rows reuse existing KRONOS2 categories,
the channel plan is stable under ablation (so one CORAL store serves every configuration), the cache
key changes with the channel set and the model, and the arm module imports without torch/coral."""
import json
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "arms"))

import declarations  # noqa: E402
import kronos2_coral as arm  # noqa: E402

DEC = declarations.load()


def test_every_declared_identity_marker_has_a_kronos2_mapping():
    missing = [m for m in DEC["identity_markers"] if m not in arm.CANONICAL]
    assert not missing, missing
    assert arm.CANONICAL[arm.NUCLEAR] == "DAPI"


def test_no_state_marker_is_mapped():
    assert not set(arm.CANONICAL) & set(DEC["state_markers"])


@pytest.mark.data
def test_novel_rows_reuse_existing_kronos2_categories():
    """Reads the vendored KRONOS2 snapshot under ``external/``, which is not in version control."""
    meta = arm.KRONOS2_DIR / "marker_metadata.csv"
    if not meta.exists():
        pytest.skip("KRONOS2 snapshot not on disk")
    import pandas as pd
    k = pd.read_csv(meta)
    comps, fams = set(k["compartment"]), set(k["family"])
    names = {str(n).lower() for n in k["marker_name"]}
    for r in arm.NOVEL_ROWS:
        assert r["compartment"] in comps, r
        assert r["family"] in fams, r
        assert r["mean"] == "" and r["std"] == "", "stats are CORAL's to fill from the donor's windows"
        assert r["marker_name"].lower() not in names, f"{r['marker_name']} is in the vocabulary; it must not be re-registered"
    # every CANONICAL target that is not novel must be a vocabulary row (separator-insensitive)
    key = lambda s: "".join(ch for ch in str(s).lower() if ch.isalnum())  # noqa: E731
    vocab = {key(n) for n in k["marker_name"]}
    for m, c in arm.CANONICAL.items():
        if m not in arm.NOVEL:
            assert key(c) in vocab, (m, c)


def test_channel_plan_is_stable_under_ablation():
    full = sorted(DEC["identity_markers"])
    raw_a, can_a = arm.channel_plan("x", full)
    raw_b, can_b = arm.channel_plan("x", sorted([m for m in full if m != "CD3e"] + ["CD3e"]))
    assert raw_a == raw_b and can_a == can_b
    assert raw_a[0] == arm.NUCLEAR and len(raw_a) == len(full) + 1


def test_cache_key_tracks_channel_set_and_model():
    region = {"window": {"y0": 0, "x0": 0, "size": 2048}}
    assets = {"sha256": {"kronos2_vitb16_teacher.pth": "a" * 64}, "coral_commit": "b" * 40}
    k1 = arm.cache_key("6436", region, ["DAPI", "CD3e", "CD8"], assets)
    k2 = arm.cache_key("6436", region, ["DAPI", "CD3e"], assets)
    k3 = arm.cache_key("6436", region, ["CD8", "DAPI", "CD3e"], assets)
    k4 = arm.cache_key("6436", region, ["DAPI", "CD3e", "CD8"], {"sha256": {"kronos2_vitb16_teacher.pth": "c" * 64}, "coral_commit": "b" * 40})
    assert k1 != k2 and k1 == k3 and k1 != k4


def test_arm_module_imports_without_the_heavy_stack():
    assert "torch" not in sys.modules and "coral" not in sys.modules and "transformers" not in sys.modules


def test_default_params_declare_the_frozen_configuration():
    p = arm.DEFAULT_PARAMS
    assert p["patch_px"] == 64 and p["batch_size"] == 128 and p["scope_supported"] == ["frame"]
    assert "false_award_measurement" in p and "mask" in p
    json.dumps({k: v for k, v in p.items() if k != "run"})  # serialisable record
