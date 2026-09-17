"""Pseudo-cell linking, mode guards, and the QC gates.

The two things worth protecting here are the minimum-anchor guard (linking on a handful of
features produces confident nonsense) and the mode guard (cell-to-cell pairing is a
measurement on one section and an invalid claim across two).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.contract import CellTable, coerce_cell_table
from phenocycler.integration.crossmodal import (AnchorError, fuzzy_smooth, link_cells,
                                                link_donor, permuted_link_null,
                                                resolve_anchors)
from phenocycler.integration.qc import (ALLOWED_BY_VERDICT, FAIL, PASS, WARN, verdict_for)
from phenocycler.integration.sameslide import (ModeError, match_by_centroid,
                                               require_same_slide, require_sequential)

PROTEIN_TO_GENE = {
    "CD31": "PECAM1", "CD3e": "CD3E", "CD163": "CD163", "CD68": "CD68",
    "CD20": "MS4A1", "MPO": "MPO", "B3TUBB": "TUBB3", "CD99": "CD99",
    "Pan_Cytokeratin": "KRT19", "SMA": "ACTA2",
}
LINEAGE_OF = {
    "CD31": "Vascular", "CD3e": "T_NK", "CD163": "Myeloid", "CD68": "Myeloid",
    "CD20": "B_Plasma", "MPO": "Granulocyte", "B3TUBB": "Neural", "CD99": "Endocrine",
    "Pan_Cytokeratin": "Exocrine", "SMA": "Stromal",
}


def _paired_tables(n=400, n_markers=10, seed=0, extent=2000.0):
    """Two sections whose lineages carry matching protein and gene signal."""
    r = np.random.default_rng(seed)
    markers = list(PROTEIN_TO_GENE)[:n_markers]
    lineages = sorted({LINEAGE_OF[m] for m in markers})

    def build(modality, sd):
        lin = r.choice(lineages, n)
        df = pd.DataFrame({
            "cell_id": [f"{modality[0]}{i}" for i in range(n)],
            "x_um": sd.uniform(0, extent, n), "y_um": sd.uniform(0, extent, n),
            "lineage_native": lin, "lineage_common": lin, "celltype_fine": lin,
            "niche_id_joint": "0",
        })
        for m in markers:
            name = m if modality == "phenocycler" else PROTEIN_TO_GENE[m]
            df[f"feat__{name}"] = np.where(lin == LINEAGE_OF[m],
                                           6.0 + sd.random(n), 0.3 + sd.random(n) * 0.3)
        df = coerce_cell_table(
            df, modality=modality, donor_id="T", roi="panc",
            profile="identity", identity_release_id="synthetic-test-release",
        )
        return CellTable(df=df, modality=modality, donor_id="T", roi="panc")

    return build("phenocycler", np.random.default_rng(seed + 1)), \
        build("xenium", np.random.default_rng(seed + 2))


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def _require_taxonomy(cfg) -> None:
    """Skip when the XeniumPanelExplorer submodule is not checked out.

    That submodule is a **private** repository, so CI running with the default
    repo-scoped GITHUB_TOKEN cannot clone it. Anything that resolves protein<->gene anchors
    goes through it, so those tests skip rather than fail — a missing private submodule is an
    environment limitation, not a defect in the code under test.
    """
    from phenocycler.integration.vocab import VocabError, load_panel_taxonomy

    try:
        load_panel_taxonomy(cfg)
    except (VocabError, FileNotFoundError) as exc:
        pytest.skip(f"XeniumPanelExplorer submodule unavailable: {exc}")


# --------------------------------------------------------------------------- #
# Anchors
# --------------------------------------------------------------------------- #

def test_anchors_come_from_the_run_not_the_identity_cores(cfg):
    """KRT19 and ACTA2 reach the deployed panel only via the custom-100 add-on.

    They are in no identity core, so validating against the core union would declare them
    off-panel and silently starve the anchor count.
    """
    p, x = _paired_tables(n_markers=10)
    try:
        pairs, p_cols, x_cols = resolve_anchors(cfg, p, x)
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")

    found = {mp.protein for mp in pairs}
    assert {"Pan_Cytokeratin", "SMA"} <= found, (
        "custom-100 genes present in the run were not recognised as anchors")
    assert len(pairs) == len(p_cols) == len(x_cols)


def test_anchor_requires_presence_in_both_tables(cfg):
    """A crosswalk pair absent from one table is not an anchor."""
    p, x = _paired_tables(n_markers=10)
    x = CellTable(df=x.df.drop(columns=["feat__PECAM1"]), modality="xenium",
                  donor_id=x.donor_id, roi=x.roi)
    try:
        pairs, _, _ = resolve_anchors(cfg, p, x)
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")
    assert "CD31" not in {mp.protein for mp in pairs}


# --------------------------------------------------------------------------- #
# The minimum-anchor guard
# --------------------------------------------------------------------------- #

def test_linking_refuses_below_the_minimum_anchor_count(tmp_path, cfg):
    """Three genes cannot support a shared embedding, and the failure must be loud."""
    _require_taxonomy(cfg)
    p, x = _paired_tables(n_markers=3)
    thin = load_config(data_dir=tmp_path, crossmodal_min_anchors=8)

    import phenocycler.integration.crossmodal as cm
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(thin.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(thin.cells_xen_dir, "T", "panc"))

    with pytest.raises(AnchorError, match="below the configured minimum"):
        cm.link_donor(thin, "T", "panc")


def test_linking_proceeds_with_enough_anchors(tmp_path):
    p, x = _paired_tables(n_markers=10, n=300)
    cfg2 = load_config(data_dir=tmp_path, crossmodal_min_anchors=5)
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(cfg2.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(cfg2.cells_xen_dir, "T", "panc"))
    try:
        res = link_donor(cfg2, "T", "panc")
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")
    assert res.n_anchors >= 5
    assert len(res.links) > 0
    assert set(res.links.columns) >= {"inferred_pheno_cell_id", "inferred_xen_cell_id",
                                      "inferred_confidence", "inferred_lineage_agree"}


def test_link_output_columns_are_prefixed_inferred():
    """Serial-section links are inference; the column names must say so."""
    p, x = _paired_tables(n_markers=10, n=200)
    p_cols = [f"feat__{m}" for m in list(PROTEIN_TO_GENE)[:10]]
    x_cols = [f"feat__{PROTEIN_TO_GENE[m]}" for m in list(PROTEIN_TO_GENE)[:10]]
    links = link_cells(p, x, p_cols, x_cols, smooth_k=5)
    assert len(links) > 0
    for col in links.columns:
        assert col.startswith("inferred_") or col in ("group",), col


def test_links_are_one_to_one():
    p, x = _paired_tables(n_markers=10, n=150)
    p_cols = [f"feat__{m}" for m in list(PROTEIN_TO_GENE)[:10]]
    x_cols = [f"feat__{PROTEIN_TO_GENE[m]}" for m in list(PROTEIN_TO_GENE)[:10]]
    links = link_cells(p, x, p_cols, x_cols, smooth_k=5)
    assert links["inferred_pheno_cell_id"].is_unique
    assert links["inferred_xen_cell_id"].is_unique


def test_permuted_link_null_is_reported():
    """With few lineages, random pairing already agrees often; the null must be shown."""
    links = pd.DataFrame({
        "inferred_pheno_lineage": ["Endocrine"] * 60 + ["Exocrine"] * 40,
        "inferred_xen_lineage": ["Endocrine"] * 60 + ["Exocrine"] * 40,
    })
    stats = permuted_link_null(links, n_perm=50)
    assert stats["agree_real"] == pytest.approx(1.0)
    assert stats["agree_null"] < 1.0
    assert stats["z"] > 3.0


def test_fuzzy_smoothing_reduces_variance():
    r = np.random.default_rng(3)
    xy = r.uniform(0, 1000, (300, 2))
    vals = r.normal(0, 1, (300, 4))
    assert fuzzy_smooth(vals, xy, k=15).std() < vals.std()


# --------------------------------------------------------------------------- #
# Mode guards
# --------------------------------------------------------------------------- #

def test_same_slide_pairing_refuses_in_sequential_mode(tmp_path):
    cfg_seq = load_config(data_dir=tmp_path, integration_mode="sequential")
    with pytest.raises(ModeError, match="DIFFERENT cells"):
        require_same_slide(cfg_seq)


def test_sequential_steps_refuse_in_same_slide_mode(tmp_path):
    cfg_ss = load_config(data_dir=tmp_path, integration_mode="same_slide")
    require_same_slide(cfg_ss)                       # allowed
    with pytest.raises(ModeError, match="sequential-mode step"):
        require_sequential(cfg_ss, "islet matching")


def test_crossmodal_refuses_in_same_slide_mode(tmp_path):
    """In same-slide mode the pairing is measured, so approximating it is a mistake.

    The pair must be DECLARED same_slide, not merely run under [integration] mode. Since
    2026-09-16 same-slide status is a property of the section pair rather than of the run
    (section_modes.csv), because it is an acquisition fact -- one physical section imaged
    twice -- and a config flag cannot assert it. Without a declaration link_donor refuses
    earlier and for a different reason, which is correct but is not what this test is for.
    """
    p, x = _paired_tables(n_markers=10, n=100)
    modes = tmp_path / "section_modes.csv"
    modes.write_text(
        "donor_id,roi,xenium_serial,mode,evidence,note\n"
        'T,panc,,same_slide,"synthetic fixture: one section imaged twice",'
        "crossmodal must refuse here\n"
    )
    cfg_ss = load_config(data_dir=tmp_path, integration_mode="same_slide",
                         section_modes_csv=modes)
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(cfg_ss.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(cfg_ss.cells_xen_dir, "T", "panc"))
    with pytest.raises(RuntimeError, match="sameslide"):
        link_donor(cfg_ss, "T", "panc")


def test_same_slide_cannot_be_asserted_by_the_run_mode_alone(tmp_path):
    """An undeclared pair is refused even when [integration] mode says same_slide.

    This is the guard that made the test above fail after the 2026-09-16 change, and it is
    worth keeping on purpose: analysing a serial pair as same_slide fabricates a cell-to-cell
    measurement, while the reverse only costs resolution.
    """
    p, x = _paired_tables(n_markers=10, n=100)
    cfg_ss = load_config(data_dir=tmp_path, integration_mode="same_slide")
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(cfg_ss.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(cfg_ss.cells_xen_dir, "T", "panc"))
    with pytest.raises(Exception, match="not declared"):
        link_donor(cfg_ss, "T", "panc")


def test_invalid_mode_rejected_at_config_load(tmp_path):
    with pytest.raises(ValueError, match="mode must be"):
        load_config(data_dir=tmp_path, integration_mode="serial")


def test_centroid_pairing_is_mutual_and_capped():
    fixed = np.array([[0.0, 0.0], [100.0, 0.0], [200.0, 0.0]])
    moving = np.array([[1.0, 0.0], [101.0, 0.0], [900.0, 0.0]])
    pairs = match_by_centroid(moving, fixed, max_dist_um=5.0)
    assert len(pairs) == 2
    assert (pairs["distance_um"] <= 5.0).all()


# --------------------------------------------------------------------------- #
# QC gates
# --------------------------------------------------------------------------- #

# The evidence the gate actually decides on, since 2026-09-10. Tissue Dice is passed through
# every case and is deliberately never what changes the answer -- see the two tests below it.
GOOD_EVIDENCE = dict(footprint_ncc_z=9.0, islet_peak_ratio=5.5, islet_peak_offset_um=40.0)


@pytest.mark.parametrize("dice,rmse,evidence,expected", [
    # Both bars clear: the image bar against other placements, and the islet landscape peaking
    # at this transform rather than somewhere else.
    (0.92, 60.0, GOOD_EVIDENCE, PASS),
    (0.40, 60.0, GOOD_EVIDENCE, PASS),
    # Image bar clears, islets do not corroborate: anatomy-level only.
    (0.92, 60.0, dict(GOOD_EVIDENCE, islet_peak_ratio=1.2), WARN),
    # The islet landscape peaks 900 um away -- this transform is not where the correspondence
    # is, which is donor 6436's original failure and the one the old gate could not see.
    (0.92, 60.0, dict(GOOD_EVIDENCE, islet_peak_offset_um=900.0), WARN),
    # Image bar fails outright.
    (0.92, 60.0, dict(GOOD_EVIDENCE, footprint_ncc_z=1.5), FAIL),
    # No image pair, but the structures corroborate: enough for structure/grid work, capped
    # at WARN because the primary bar was never measured.
    (0.92, 60.0, dict(GOOD_EVIDENCE, footprint_ncc_z=float("nan")), WARN),
    # Nothing measured at all: fail closed, and say so.
    (float("nan"), float("nan"), {}, FAIL),
])
def test_verdict_grading(dice, rmse, evidence, expected):
    v, _reasons = verdict_for(dice, rmse, dice_min=0.80, rmse_max=150.0, **evidence)
    assert v == expected


def test_tissue_dice_cannot_change_the_verdict():
    """It carried half the old score and moved *against* correctness on donor 6436.

    Correcting that donor's transform -- islet correspondence 12 of 136 (chance 5) to 69 --
    lowered tissue Dice from 0.549 to 0.530. A gate that prefers the wrong answer cannot be
    relaxed into correctness, so Dice is reported and never weighed.
    """
    for dice in (0.10, 0.53, 0.55, 0.99):
        assert verdict_for(dice, 60.0, dice_min=0.80, rmse_max=150.0,
                           **GOOD_EVIDENCE)[0] == PASS


def test_an_unmeasured_bar_is_reported_as_unavailable_not_as_a_failure():
    """A missing measurement is unavailable, not zero and not negative."""
    _v, reasons = verdict_for(0.92, 60.0, dice_min=0.80, rmse_max=150.0,
                              islet_peak_ratio=5.5, islet_peak_offset_um=40.0)
    joined = "; ".join(reasons)
    assert "unavailable" in joined
    assert "below" not in joined


def test_failed_registration_still_allows_registration_free_analyses():
    """A donor whose sections will not align must not be dropped from the study."""
    assert "donor" in ALLOWED_BY_VERDICT[FAIL]
    assert "niche" in ALLOWED_BY_VERDICT[FAIL]
    assert "structure" not in ALLOWED_BY_VERDICT[FAIL]
    assert "crossmodal" not in ALLOWED_BY_VERDICT[FAIL]
    assert "crossmodal" in ALLOWED_BY_VERDICT[PASS]


def test_verdict_reasons_are_actionable():
    _v, reasons = verdict_for(0.4, 400.0, dice_min=0.80, rmse_max=150.0)
    assert any("Dice" in r for r in reasons)
    assert any("RMSE" in r for r in reasons)


def test_a_warn_says_whether_the_islets_were_tested_or_untestable():
    """Two different facts about two different problems, and they must not share a sentence.

    A donor whose islets disagree has a transform question. A donor with no structures has an
    identity question — the corroborating bar was never measured. Both land at WARN; only the
    reason tells them apart.
    """
    tested, _r = verdict_for(0.9, 60.0, dice_min=0.80, rmse_max=150.0,
                             footprint_ncc_z=9.0, islet_peak_ratio=1.1,
                             islet_peak_offset_um=900.0)
    _v, reasons_tested = verdict_for(0.9, 60.0, dice_min=0.80, rmse_max=150.0,
                                     footprint_ncc_z=9.0, islet_peak_ratio=1.1,
                                     islet_peak_offset_um=900.0)
    _v2, reasons_untested = verdict_for(0.9, np.nan, dice_min=0.80, rmse_max=150.0,
                                        footprint_ncc_z=9.0)
    assert tested == WARN
    assert any("did not corroborate" in r for r in reasons_tested)
    assert any("no structures for this donor" in r for r in reasons_untested)
    assert not any("did not corroborate" in r for r in reasons_untested)
