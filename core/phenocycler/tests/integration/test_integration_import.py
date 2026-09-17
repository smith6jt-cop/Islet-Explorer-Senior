"""Xenium import: null handling, label fallback, and the same-slide docstring contract.

Regression guards for the three review findings on PR #5. Two of them fail *silently* — a
null label becomes the string "nan" and a fallback path quietly stops firing — so they need
tests rather than just a fix.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.import_xenium import _obs_str, import_sample

anndata = pytest.importorskip("anndata")


def _h5ad(obs: pd.DataFrame, tmp_path, n_genes: int = 6):
    """Write a minimal phenotyped-style h5ad; `import_sample` reads a path, not an object."""
    rng = np.random.default_rng(0)
    n = len(obs)
    X = rng.poisson(0.5, (n, n_genes)).astype("float32")
    a = anndata.AnnData(X=X, obs=obs,
                        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))
    a.layers["lognorm"] = np.log1p(X)
    a.obsm["spatial"] = np.column_stack([rng.uniform(0, 2000, n),
                                         rng.uniform(0, 2000, n)]).astype("float32")
    path = Path(tmp_path) / "0059865_phenotyped.h5ad"
    a.write_h5ad(path)
    return path


# --------------------------------------------------------------------------- #
# _obs_str null handling
# --------------------------------------------------------------------------- #

def test_obs_str_preserves_nulls_as_the_default():
    """`astype(str)` before `fillna` turns NaN into the literal string "nan"."""
    obs = pd.DataFrame({"celltype": pd.Series(["Endocrine", None, "Exocrine"], dtype=object)})
    out = _obs_str(obs, "celltype")
    assert list(out) == ["Endocrine", "", "Exocrine"]
    assert "nan" not in set(out)


def test_obs_str_handles_categorical_columns():
    """AnnData obs columns are routinely Categorical, where `fillna` raises outright."""
    obs = pd.DataFrame({"celltype": pd.Categorical(["Beta", None, "Alpha"])})
    out = _obs_str(obs, "celltype")
    assert list(out) == ["Beta", "", "Alpha"]


def test_obs_str_respects_a_custom_default():
    obs = pd.DataFrame({"x": pd.Series([None, "a"], dtype=object)})
    assert list(_obs_str(obs, "x", default="?")) == ["?", "a"]


def test_obs_str_missing_column_returns_default():
    obs = pd.DataFrame({"other": [1, 2]})
    assert list(_obs_str(obs, "celltype", default="-")) == ["-", "-"]


# --------------------------------------------------------------------------- #
# The second-order bug: a silently dead fallback
# --------------------------------------------------------------------------- #

def test_embedded_fine_labels_do_not_create_identity_without_a_release(tmp_path):
    """A label-shaped H5AD column is not release provenance."""
    # Categorical with unset values, which is how AnnData actually stores a null label
    # (h5py cannot serialise an object column containing None at all).
    n = 20
    obs = pd.DataFrame({
        "celltype_broad": pd.Categorical(["Endocrine"] * n),
        "celltype": pd.Categorical([None] * n, categories=["Beta", "Alpha"]),
        "celltype_fine": pd.Categorical(["Beta"] * 12 + ["Alpha"] * 8),
    }, index=[f"c{i}" for i in range(n)])

    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc")
    assert (df["celltype_fine"] == "").all()
    assert (df["lineage_native"] == "").all()


def test_null_labels_do_not_leak_the_string_nan(tmp_path):
    """A null broad label must become '' (and map to Unassigned), never the label "nan"."""
    n = 10
    obs = pd.DataFrame({
        "celltype_broad": pd.Categorical(["Endocrine"] * 6 + [None] * 4),
        "celltype": pd.Categorical(["Beta"] * 6 + [None] * 4),
    }, index=[f"c{i}" for i in range(n)])

    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc")
    assert "nan" not in set(df["lineage_native"])
    assert "nan" not in set(df["celltype_fine"])
    assert (df["lineage_native"] == "").all()
    assert (df["lineage_common"] == "").all()


def test_import_sample_stamps_the_join_key(tmp_path):
    """donor_id/roi reach the data here and nowhere else in the Xenium world."""
    obs = pd.DataFrame({"celltype_broad": ["Exocrine"] * 5, "celltype": ["Acinar"] * 5},
                       index=[f"c{i}" for i in range(5)])
    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc",
                       sample_key="0059865__panc")
    assert set(df["donor_id"]) == {"6539"}
    assert set(df["roi"]) == {"panc"}
    assert set(df["xenium_sample_key"]) == {"0059865__panc"}


def test_broad_lowsignal_is_not_reinterpreted_without_a_release(tmp_path):
    """Embedded fallback flags cannot open the identity contract either."""
    n = 10
    obs = pd.DataFrame({
        "celltype_broad": ["Exocrine"] * n,
        "celltype": ["Acinar"] * n,
        "broad_lowsignal": [True] * 3 + [False] * 7,
    }, index=[f"c{i}" for i in range(n)])
    df = import_sample(_h5ad(obs, tmp_path), donor_id="1", roi="panc")
    assert df["evidence_positive"].isna().all()


def test_arbitrary_release_id_string_is_not_an_import_escape_hatch(tmp_path):
    obs = pd.DataFrame({"celltype_broad": ["Exocrine"]}, index=["c0"])
    with pytest.raises(TypeError, match="identity_release_id"):
        import_sample(
            _h5ad(obs, tmp_path), donor_id="1", roi="panc",
            identity_release_id="not-a-validated-release",
        )


# --------------------------------------------------------------------------- #
# release coverage is a property of the SECTION, not the donor
# --------------------------------------------------------------------------- #

# The label space the declared Xenium release at data/development/xenium_typing/release
# already carries: Immune_T -> T_NK, Immune_B -> B_Plasma, Myeloid -> Myeloid.
_LABEL_SPACE = {
    "schema_version": 1,
    "space_id": "synthetic-space",
    "producer_lane": "synthetic-lane",
    "labels": [
        {"label": "Immune_T", "common_lineage": "T_NK"},
        {"label": "Immune_B", "common_lineage": "B_Plasma"},
        {"label": "Myeloid", "common_lineage": "Myeloid"},
    ],
}


def _write_release(tmp_path, *, donor="6442", roi_calls, release_id="synthetic",
                   dirname="release", donor_level_calls=None):
    """Write a declared xenium release.

    ``roi_calls`` maps roi -> list of (cell_id, label) pairs, or None for a partition
    directory that exists with no shard in it. ``donor_level_calls``, when given, writes one
    shard directly under ``donor_id=<donor>/`` instead -- the A.4 donor-level layout, which
    has no ROI partitions at all and so resolves for *every* roi asked of it.
    """
    import json

    release_dir = tmp_path / dirname
    calls = release_dir / "calls" / f"donor_id={donor}"
    calls.mkdir(parents=True)
    for roi, pairs in (roi_calls or {}).items():
        part = calls / f"roi={roi}"
        part.mkdir()
        if pairs is not None:
            _write_calls(part / "data_0.parquet", pairs)
    if donor_level_calls is not None:
        _write_calls(calls / "data_0.parquet", donor_level_calls)
    (release_dir / "label_space.json").write_text(json.dumps(_LABEL_SPACE))
    (release_dir / "release.json").write_text(json.dumps({
        "schema_version": 1,
        "release_id": release_id,
        "modality": "xenium",
        "producer": {"lane": "synthetic-lane", "commit": "abc123"},
        "label_space": "label_space.json",
        "calls_dir": "calls",
        "cell_key": "cell_id",
        "donors": [donor],
        "inputs": [{"kind": "expression", "run_id": "run1", "content_id": "c1",
                    "path": "expression"}],
        "columns": {"label": "phenotype", "common_lineage": "common_lineage"},
    }))
    return release_dir


def _write_calls(path, pairs):
    labels = [label for _cid, label in pairs]
    common = {rule["label"]: rule["common_lineage"] for rule in _LABEL_SPACE["labels"]}
    pd.DataFrame({
        "cell_id": [cid for cid, _label in pairs],
        "phenotype": labels,
        "common_lineage": [common[label] for label in labels],
    }).to_parquet(path, index=False)


def _panc_only_release(tmp_path, donor="6442"):
    """A declared xenium release that partitions ``donor`` by ``roi=panc`` only."""
    return _write_release(tmp_path, donor=donor, roi_calls={"panc": None},
                          release_id="synthetic-panc-only")


def _bundle(tmp_path, name, cell_ids):
    """A minimal Xenium output bundle; `_resolve_source` reads cells.parquet from it."""
    bundle = tmp_path / name
    bundle.mkdir(parents=True)
    n = len(cell_ids)
    pd.DataFrame({
        "cell_id": list(cell_ids),
        "x_centroid": [float(i) for i in range(n)],
        "y_centroid": [float(i) + 5.0 for i in range(n)],
        "cell_area": [10.0 + i for i in range(n)],
    }).to_parquet(bundle / "cells.parquet", index=False)
    return bundle


def _manifest(tmp_path, rows):
    """Write data/integration/manifest.csv and return the data dir."""
    data = tmp_path / "data"
    path = data / "integration" / "manifest.csv"
    path.parent.mkdir(parents=True)
    base = {"disease_status": "ND", "pair_status": "xenium_only", "pheno_image": "",
            "pheno_geojson": "", "xenium_serial": "0059947", "xenium_batch": "",
            "xenium_zarr": "", "xenium_h5ad": "", "section_gap_um": "", "block_id": ""}
    pd.DataFrame([{**base, **row} for row in rows]).to_csv(path, index=False)
    return data


def test_a_panc_only_release_does_not_abort_the_lymph_node_import(tmp_path):
    """The pLN sections are a different piece of tissue in a different cell-id namespace.

    ``attach_release_identity`` resolves calls by (donor, roi) and raises when that pair has
    none, so testing only ``release.donors`` sent donor 6442's lymph node into an attachment
    that raised and aborted the whole import -- after earlier donors had already been
    written. The section, not the donor, decides whether identity attaches.
    """
    from phenocycler.config import PipelineConfig
    from phenocycler.integration.contract import read_cell_table
    from phenocycler.integration.import_xenium import run_import_xenium

    data = tmp_path / "data"
    bundle = tmp_path / "output-XETG00298__0059947__pLN__20250401__152313"
    bundle.mkdir(parents=True)
    pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(4)],
        "x_centroid": [1.0, 2.0, 3.0, 4.0],
        "y_centroid": [5.0, 6.0, 7.0, 8.0],
        "cell_area": [10.0, 11.0, 12.0, 13.0],
    }).to_parquet(bundle / "cells.parquet", index=False)

    manifest = data / "integration" / "manifest.csv"
    manifest.parent.mkdir(parents=True)
    pd.DataFrame([{
        "donor_id": "6442", "roi": "pln_1", "disease_status": "ND",
        "pair_status": "xenium_only", "pheno_image": "", "pheno_geojson": "",
        "xenium_sample_key": "0059947__pln_1", "xenium_serial": "0059947",
        "xenium_batch": "", "xenium_bundle": str(bundle), "xenium_zarr": "",
        "xenium_h5ad": "", "section_gap_um": "", "block_id": "",
        "tissue": "pancreatic_lymph_node",
    }]).to_csv(manifest, index=False)

    cfg = PipelineConfig(data_dir=data,
                         xenium_identity_release=_panc_only_release(tmp_path))

    summary = run_import_xenium(cfg, roi="pln_1", only_paired=False)

    assert list(summary["identity"]) == ["geometry only"]
    assert int(summary["n_cells"].iloc[0]) == 4
    table = read_cell_table(cfg.cells_xen_dir, "6442", "pln_1", modality="xenium")
    assert table.tissue == "pancreatic_lymph_node"
    assert table.df["lineage_common"].eq("").all()
    assert table.df["evidence_positive"].isna().all()


# --------------------------------------------------------------------------- #
# the guard is STRUCTURAL: a release that covers roi=pln_1 must attach to roi=pln_1
# --------------------------------------------------------------------------- #

def _load_release(release_dir, data_dir):
    from phenocycler.config import PipelineConfig
    from phenocycler.integration.phenospace import load_release

    cfg = PipelineConfig(data_dir=data_dir, xenium_identity_release=release_dir)
    return load_release(cfg, "xenium")


def _two_section_cohort(tmp_path):
    """One serial, two separately segmented bundles, disjoint cell-id namespaces."""
    panc = _bundle(tmp_path, "output-XETG00298__0059947__Panc__20250401__152313",
                   [f"panc-{i}" for i in range(4)])
    pln = _bundle(tmp_path, "output-XETG00298__0059947__pLN__20250401__152313",
                  [f"pln-{i}" for i in range(4)])
    data = _manifest(tmp_path, [
        {"donor_id": "6442", "roi": "panc", "tissue": "pancreas",
         "xenium_sample_key": "0059947__panc", "xenium_bundle": str(panc)},
        {"donor_id": "6442", "roi": "pln_1", "tissue": "pancreatic_lymph_node",
         "xenium_sample_key": "0059947__pln_1", "xenium_bundle": str(pln)},
    ])
    return data, panc, pln


_PLN_CALLS = [("pln-0", "Immune_T"), ("pln-1", "Immune_B"),
              ("pln-2", "Myeloid"), ("pln-3", "Immune_T")]
_PANC_CALLS = [("panc-0", "Immune_T"), ("panc-1", "Immune_T"),
               ("panc-2", "Myeloid"), ("panc-3", "Immune_B")]


def test_a_release_covering_roi_pln_1_attaches_to_the_node_and_not_the_pancreas(tmp_path):
    """The immune vocabulary is licensed on pLN, so a pLN release must actually land.

    The guard that keeps pancreas calls off the lymph node is structural -- disjoint
    cell-id namespaces -- not a refusal to call immune cells in a node. A release that
    partitions roi=pln_1 has its own calls in the node's own namespace, and import must
    attach them, with no code change and no bleed into the pancreas section.
    """
    from phenocycler.config import PipelineConfig
    from phenocycler.integration.contract import read_cell_table
    from phenocycler.integration.import_xenium import run_import_xenium

    data, _panc, _pln = _two_section_cohort(tmp_path)
    release = _write_release(tmp_path, roi_calls={"pln_1": _PLN_CALLS},
                             release_id="synthetic-pln-only")
    cfg = PipelineConfig(data_dir=data, xenium_identity_release=release)

    summary = run_import_xenium(cfg, only_paired=False).set_index("roi")
    assert summary.loc["pln_1", "identity"] == "attached"
    assert summary.loc["panc", "identity"] == "geometry only"

    node = read_cell_table(cfg.cells_xen_dir, "6442", "pln_1", modality="xenium").df
    node = node.set_index("cell_id")
    assert list(node.loc[[c for c, _l in _PLN_CALLS], "lineage_native"]) == \
        [label for _c, label in _PLN_CALLS]
    assert list(node.loc[[c for c, _l in _PLN_CALLS], "lineage_common"]) == \
        ["T_NK", "B_Plasma", "Myeloid", "T_NK"]

    panc = read_cell_table(cfg.cells_xen_dir, "6442", "panc", modality="xenium").df
    assert panc["lineage_common"].eq("").all(), "a pLN release must not bleed into the pancreas"
    assert panc["lineage_native"].eq("").all()


def test_a_release_covering_roi_panc_attaches_to_the_pancreas_and_not_the_node(tmp_path):
    """The converse direction of the same partition test."""
    from phenocycler.config import PipelineConfig
    from phenocycler.integration.contract import read_cell_table
    from phenocycler.integration.import_xenium import run_import_xenium

    data, _panc, _pln = _two_section_cohort(tmp_path)
    release = _write_release(tmp_path, roi_calls={"panc": _PANC_CALLS},
                             release_id="synthetic-panc-calls")
    cfg = PipelineConfig(data_dir=data, xenium_identity_release=release)

    summary = run_import_xenium(cfg, only_paired=False).set_index("roi")
    assert summary.loc["panc", "identity"] == "attached"
    assert summary.loc["pln_1", "identity"] == "geometry only"

    panc = read_cell_table(cfg.cells_xen_dir, "6442", "panc", modality="xenium").df
    assert set(panc["lineage_common"]) == {"T_NK", "Myeloid", "B_Plasma"}
    node = read_cell_table(cfg.cells_xen_dir, "6442", "pln_1", modality="xenium").df
    assert node["lineage_common"].eq("").all()


# --------------------------------------------------------------------------- #
# forcing a pancreas release onto a node section: refused, and refused loudly
# --------------------------------------------------------------------------- #

def test_attaching_a_pancreas_release_to_a_node_section_is_refused_and_names_the_roi(tmp_path):
    """`attach_release_identity` stays strict; the import-side skip is not the protection."""
    from phenocycler.integration.import_xenium import import_sample
    from phenocycler.integration.phenospace import (
        PhenotypeSpaceError,
        attach_release_identity,
    )

    data, _panc, pln = _two_section_cohort(tmp_path)
    release = _load_release(
        _write_release(tmp_path, roi_calls={"panc": _PANC_CALLS},
                       release_id="synthetic-panc-calls"),
        data,
    )
    geometry = import_sample(pln, donor_id="6442", roi="pln_1",
                             tissue="pancreatic_lymph_node")

    with pytest.raises(PhenotypeSpaceError) as excinfo:
        attach_release_identity(geometry, release, "6442")
    message = str(excinfo.value)
    assert "pln_1" in message, "the refusal must name the roi it could not cover"
    assert "no calls" in message


def test_a_partial_cell_id_overlap_is_refused_rather_than_mostly_matched(tmp_path):
    """The measured failure mode if someone forces the attachment anyway.

    panc and pln_1 are separately segmented bundles, so their cell_id namespaces are
    disjoint and any apparent overlap is coincidence: measured 2026-09-17, 6442's 958,744
    pancreas calls hit 44 of its 177,799 pln_1 cells. A donor-level release -- one with no
    ROI partitions, so `calls_for` resolves it for every roi -- is the layout that can
    actually deliver pancreas calls to a node section, and the join must refuse the
    near-empty result instead of accepting a 3-of-100 overlap as mostly matched.
    """
    from phenocycler.integration.import_xenium import import_sample
    from phenocycler.integration.phenospace import (
        PhenotypeSpaceError,
        attach_release_identity,
    )

    node_ids = [f"pln-{i}" for i in range(100)]
    pln = _bundle(tmp_path, "output-XETG00298__0059947__pLN__20250401__152313", node_ids)
    data = _manifest(tmp_path, [
        {"donor_id": "6442", "roi": "pln_1", "tissue": "pancreatic_lymph_node",
         "xenium_sample_key": "0059947__pln_1", "xenium_bundle": str(pln)},
    ])
    # 200 pancreas calls, of which 3 happen to collide with a node cell_id.
    collisions = node_ids[:3]
    pancreas_calls = [(f"panc-{i}", "Immune_T") for i in range(200)] + \
        [(cid, "Immune_T") for cid in collisions]
    release = _load_release(
        _write_release(tmp_path, donor_level_calls=pancreas_calls, roi_calls={},
                       release_id="synthetic-donor-level"),
        data,
    )
    assert release.calls_for("6442", "pln_1"), (
        "a donor-level release resolves for every roi; that is what makes this reachable")

    geometry = import_sample(pln, donor_id="6442", roi="pln_1",
                             tissue="pancreatic_lymph_node")
    with pytest.raises(PhenotypeSpaceError) as excinfo:
        attach_release_identity(geometry, release, "6442")
    message = str(excinfo.value)
    assert "97 geometry cells lack calls" in message
    assert "200 calls lack geometry" in message
    assert "Refusing to change the frame" in message


# --------------------------------------------------------------------------- #
# same-slide docstring must match the implementation
# --------------------------------------------------------------------------- #

def test_sameslide_docstring_matches_what_pair_donor_does():
    """The docstring must describe the matchers that are actually reachable.

    It once claimed IoU was preferred while `pair_donor` only ever called centroids; it now
    claims both are available with `auto` selecting between them. Either way the risk is the
    same — the two drift apart silently, and a reader has no way to tell which criterion ran.
    """
    import inspect

    from phenocycler.integration import sameslide

    doc = inspect.getdoc(sameslide) or ""
    src = "".join(inspect.getsource(f) for f in (sameslide.pair_donor, sameslide._match,
                                                 sameslide.match_tables_by_iou))

    assert "match_by_centroid(" in src, "centroid matching must remain reachable"
    assert "match_by_iou(" in src, "IoU matching must be reachable from pair_donor"
    # And the docstring must not still describe IoU as unavailable.
    for stale in ("not yet wired in", "does not use it", "no same-slide data"):
        assert stale not in doc, f"module docstring still says {stale!r} about IoU"
    assert "auto" in doc and "iou" in doc.lower(), (
        "the docstring must say how the matcher is selected")


def test_pair_donor_records_which_matcher_ran():
    """`auto` falling back silently would make a polygon-loading failure look like a
    successful centroid run, and the two have different error characteristics."""
    import inspect

    from phenocycler.integration import sameslide

    src = inspect.getsource(sameslide.pair_donor)
    assert '"method": used' in src or '"method"' in src, (
        "stats must record which matcher produced the pairs")
