"""Synthetic-only tests for the report lane: no run data, no release, no bundle is read."""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import identity_io as io  # noqa: E402
import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord, Manifest  # noqa: E402


def test_wrap_never_exceeds_width():
    long = "word " * 60
    for line in style.wrap(long, 40).split("\n"):
        assert len(line) <= 40
    assert style.wrap("a\nb", 10) == "a\nb"


def test_unavailable_patch_prefix_once():
    p = style.unavailable_patch("no release")
    assert p.get_label().startswith("Unavailable: ") and p.get_label().count("Unavailable") == 1


def test_axis_label_from_zero_omits_empty_unit():
    assert "()" not in style.axis_label_from_zero("area", "", False)
    assert "(um^2" in style.axis_label_from_zero("area", "um^2", False) or "um^2" in style.axis_label_from_zero("area", "um^2", False)


def test_status_display_and_order():
    assert style.status_display("AAB") == "Aab+"
    assert style.status_axis()[0] == "ND" and style.status_axis()[-1] == "T1D"


def test_stacked_bars_fractions_sum_to_one_and_ambiguous_hatched():
    common = {"Endocrine": "Endocrine", "T_NK": "T_NK"}
    colors = io.common_colors()
    table = pd.DataFrame({"Endocrine": [0.5, 0.2], "T_NK": [0.3, 0.3], io.AMBIGUOUS: [0.1, 0.4], io.UNAVAILABLE: [0.1, 0.1]}, index=["a", "b"])
    fig, ax = plt.subplots()
    handles = io.stacked_bars(ax, table, ["Endocrine", "T_NK", io.AMBIGUOUS, io.UNAVAILABLE], common, colors)
    tops = sorted(r.get_y() + r.get_height() for h in handles for r in h)
    assert abs(tops[-1] - 1.0) < 1e-9
    amb = [h for h in handles if h.get_label() == io.AMBIGUOUS][0]
    assert amb[0].get_hatch() == "///"
    plt.close(fig)


def test_label_color_sentinels():
    common = {"Epithelial_nonendocrine": "Exocrine"}
    colors = io.common_colors()
    assert io.label_color("Epithelial_nonendocrine", common, colors) == colors["Exocrine"]
    assert io.label_color(io.UNAVAILABLE, common, colors) == "#BBBBBB"
    assert io.label_color(io.AMBIGUOUS, common, colors) != io.label_color(io.UNAVAILABLE, common, colors)


def test_release_readers_return_none_without_release(tmp_path, monkeypatch):
    monkeypatch.setitem(io.RELEASE, "phenocycler", tmp_path / "nope")
    assert io.release_json("phenocycler") is None
    assert io.release_calls("phenocycler", "6436") is None
    assert io.common_of("phenocycler") == {}
    assert io.release_inputs("phenocycler") == []


def test_common_of_reads_label_space(tmp_path, monkeypatch):
    rel = tmp_path / "release"; rel.mkdir()
    (rel / "label_space.json").write_text(json.dumps({"labels": [
        {"label": "Endocrine", "common_lineage": "Endocrine"},
        {"label": "Epithelial_nonendocrine", "common_lineage": "Exocrine", "reason": "no acinar- or ductal-specific marker is declared on this panel"}]}))
    monkeypatch.setitem(io.RELEASE, "phenocycler", rel)
    c = io.common_of("phenocycler")
    assert c == {"Endocrine": "Endocrine", "Epithelial_nonendocrine": "Exocrine"}
    assert "acinar" not in json.dumps(list(c.values())).lower()


def test_manifest_merges_existing_entries(tmp_path):
    f1 = tmp_path / "a.txt"; f1.write_text("a")
    m = Manifest(tmp_path / "m.json")
    m.add_figure(FigureRecord(figure_id="P99", donor="6436", section="pheno", files=[str(f1)], inputs=[str(f1)], facts={"n": 1}))
    m.write()
    m2 = Manifest(tmp_path / "m.json")
    m2.add_figure(FigureRecord(figure_id="X99", donor="6436", section="xen", files=[str(f1)], inputs=[], facts={}))
    doc = json.loads(m2.write().read_text())
    assert set(doc["figures"]) == {"pheno/6436/P99", "xen/6436/X99"}
    assert doc["figures"]["pheno/6436/P99"]["facts"] == {"n": 1}


def test_not_available_panel_renders_and_names_reason(tmp_path, monkeypatch):
    import pheno_identity
    monkeypatch.setattr(pins, "OUT", tmp_path)
    prov = {"run_id": "run", "run_content_id": "c" * 40, "expression_content_id": "e" * 64, "expression_output_sha256": "o" * 64,
            "stage_content_ids": {}, "code": {"commit": "abc", "dirty": False}, "donors": {}, "identity_releases": {}, "written_utc": "now"}
    rec = pheno_identity._unavailable("P06", "6436", prov, "no PhenoCycler release yet")
    assert rec.facts["not_available"].startswith("no PhenoCycler release")
    assert any(f.endswith(".png") for f in rec.files) and any(f.endswith(".pdf") for f in rec.files)
    assert all(pathlib.Path(f).exists() for f in rec.files)


def test_pair_figures_declared_once_per_module():
    import multi, pheno_identity, xen_identity
    ids = []
    for mod in (pheno_identity, xen_identity, multi):
        ids += list(mod.PER_DONOR) + list(mod.PAIR)
        assert not set(mod.PER_DONOR) & set(mod.PAIR)
    assert len(ids) == len(set(ids)) == 18
    assert set(ids) >= {"P06", "P07", "P08", "P09", "P10", "P11", "X04", "X05", "X06", "X07", "X08", "X09", "M01", "M02", "M03", "M04", "M05", "M06"} - {"X10"}


def test_no_report_module_reads_the_explorer_bundle():
    for f in HERE.parent.glob("*.py"):
        assert "explorer_bundle" not in f.read_text() or f.name == "pins.py", f.name


def test_dev_line_is_present_on_every_title():
    fig, _ = style.figure(1, 1)
    style.title(fig, "T")
    assert any("NO_VERIFIED_TYPING" in t.get_text() for t in fig.texts)
    plt.close(fig)


def test_footer_names_release_ids_when_present():
    prov = {"run_id": "run", "run_content_id": "c" * 40, "expression_content_id": "e" * 64, "expression_output_sha256": "o" * 64,
            "stage_content_ids": {}, "code": {"commit": "abc", "dirty": True}, "donors": {}, "written_utc": "now",
            "identity_releases": {"phenocycler": "broadtype-tasp-dev-2026-09-06"}}
    t = style.footer_text(prov)
    assert "broadtype-tasp-dev-2026-09-06" in t and "(dirty)" in t
    prov["identity_releases"] = {}
    assert "none (QC figures need none)" in style.footer_text(prov)
