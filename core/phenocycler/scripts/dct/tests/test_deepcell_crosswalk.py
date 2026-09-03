"""The crosswalk must stay true to the *installed* DeepCellTypes package.

`deepcell_types.__version__` is '0.1.0' on every master commit, so the only thing that can
tell you the marker vocabulary or the alias table moved under you is a commit pin plus a
test that reads it back. These tests are the readback.

They need no checkpoint and no network: `list_supported_markers` / `resolve_supported_marker`
read the packaged `vocab.json`.
"""
import json
import os
import pathlib

import pytest

DCT = pathlib.Path(__file__).resolve().parents[1]
REPO = DCT.parents[3]
CROSSWALK = DCT / "deepcell_types_crosswalk.json"
MANIFEST = REPO / "core" / "phenocycler" / "data" / "qupath_manifest.json"

# Expected installed commits. Kept here AND in pyproject.toml's `dct` extra, which is the
# point: this test is what notices when the two drift.
PINNED = {
    "deepcell_types": "59d74a97711848dabc30fc01e242b6b1e4d3d3dc",
    "deepcell_auth": "523a4dd5d29941e076ffd378c6b032122589bdb2",
}


@pytest.fixture(scope="module")
def crosswalk():
    return json.loads(CROSSWALK.read_text())


@pytest.fixture(scope="module")
def dct_utils():
    """The package's own vocabulary helpers, forced onto the packaged vocab.json.

    DEEPCELL_TYPES_ZARR_PATH would silently swap the registry for an archive's, making
    every assertion below a statement about a different vocabulary.
    """
    utils = pytest.importorskip("deepcell_types.utils")
    old = os.environ.pop("DEEPCELL_TYPES_ZARR_PATH", None)
    try:
        yield utils
    finally:
        if old is not None:
            os.environ["DEEPCELL_TYPES_ZARR_PATH"] = old


def _dist_commit(dist_prefix):
    """The commit pip actually installed, from the dist-info direct_url.json."""
    import importlib.metadata as md

    for dist in md.distributions():
        if (dist.metadata["Name"] or "").replace("-", "_") != dist_prefix:
            continue
        raw = dist.read_text("direct_url.json")
        if raw is None:
            return None
        return json.loads(raw).get("vcs_info", {}).get("commit_id")
    return None


# --------------------------------------------------------------------------- pins


@pytest.mark.parametrize("pkg,commit", sorted(PINNED.items()))
def test_installed_commit_matches_the_pin(pkg, commit):
    pytest.importorskip(pkg)
    got = _dist_commit(pkg)
    if got is None:
        pytest.skip(f"{pkg} was not installed from git; no commit to check")
    assert got == commit, (
        f"{pkg} is installed at {got}, pinned at {commit}. Either the pin in "
        "pyproject.toml [dct] moved or someone reinstalled from master. Resolve before "
        "trusting any output — __version__ cannot tell these apart."
    )


def test_crosswalk_records_the_same_package_commit(crosswalk):
    assert crosswalk["package"]["name"] == "deepcell-types"
    assert crosswalk["package"]["commit"] == PINNED["deepcell_types"]


# --------------------------------------------------------------------- vocabulary


def test_vocabulary_is_the_released_one(dct_utils):
    """A checkpoint is refused unless the vocabulary matches it exactly, so guard the size."""
    assert len(dct_utils.list_supported_markers()) == 278
    assert len(dct_utils.list_supported_cell_types()) == 51
    assert dct_utils.list_model_versions() == ["2026-06-15", "2026-06-23-ptft"]


def test_every_rename_target_exists(crosswalk, dct_utils):
    supported = set(dct_utils.list_supported_markers())
    bad = {k: e["target"] for k, e in crosswalk["entries"].items()
           if e["target"] is not None and e["target"] not in supported}
    assert not bad, f"crosswalk targets absent from the model's registry: {bad}"


def test_rename_targets_are_canonical_not_aliases(crosswalk, dct_utils):
    """A target must be the canonical name, so the recorded value is what the model sees."""
    off = {k: (e["target"], dct_utils.resolve_supported_marker(e["target"]))
           for k, e in crosswalk["entries"].items()
           if e["target"] is not None
           and dct_utils.resolve_supported_marker(e["target"]) != e["target"]}
    assert not off, f"targets that are aliases rather than canonical names: {off}"


def test_no_entry_is_redundant(crosswalk, dct_utils):
    """The file's own invariant: it holds ONLY names the package cannot resolve itself.

    A redundant entry is how this file silently forks from upstream — the package learns an
    alias, we keep overriding it, and the two disagree with nothing to say so.
    """
    redundant = {k: dct_utils.resolve_supported_marker(k)
                 for k in crosswalk["entries"]
                 if dct_utils.resolve_supported_marker(k) is not None}
    assert not redundant, (
        f"the package already resolves these on its own; drop them: {redundant}"
    )


# ------------------------------------------------------------------------- schema


def test_schema_is_well_formed(crosswalk):
    assert crosswalk["schema_version"] == 1
    assert "channel_name" in crosswalk["key_space"]
    for name, e in crosswalk["entries"].items():
        assert set(e) == {"confidence", "reason", "target"}, f"{name} has unexpected keys"
        assert e["confidence"] in {"settled", "arguable"}, name
        assert e["reason"].strip(), f"{name} has no reason"
        assert e["target"] is None or isinstance(e["target"], str), name


def test_every_exclusion_says_why(crosswalk):
    """A null is a marker dropped from the panel. It never gets to be silent."""
    for name, e in crosswalk["entries"].items():
        if e["target"] is None:
            assert len(e["reason"].strip()) >= 15, f"{name} is excluded without a real reason"


def test_the_recorded_judgement_calls_are_still_flagged(crosswalk):
    """These four are substitutions, not synonyms; they must never quietly become 'settled'."""
    for name in ("INS", "Ker8-18", "TCF-1", "IDO1"):
        assert crosswalk["entries"][name]["confidence"] == "arguable", (
            f"{name} is a judgement call and has to stay marked as one"
        )


# ------------------------------------------------------- the panel, end to end


@pytest.mark.data
def test_the_whole_panel_resolves(crosswalk, dct_utils):
    """Every OME channel_name in the cohort is either translated, excluded, or resolved.

    A name that falls through all three is dropped by the model without comment, and the
    cell is typed on a smaller panel than the run record claims.
    """
    if not MANIFEST.exists():
        pytest.skip(f"{MANIFEST} not present (data is gitignored)")
    manifest = json.loads(MANIFEST.read_text())
    names = sorted({c["channel_name"] for im in manifest["images"] for c in im["channel_map"]})
    assert names, "manifest carried no channel names"

    entries = crosswalk["entries"]
    leftover = [n for n in names
                if n not in entries and dct_utils.resolve_supported_marker(n) is None]
    assert not leftover, f"names neither crosswalked nor resolvable: {leftover}"


# ------------------------------------------------- the DCT class -> broad class map

BROAD = DCT / "broad_class_crosswalk.json"


@pytest.fixture(scope="module")
def broad():
    return json.loads(BROAD.read_text())


def test_every_emittable_class_is_mapped(broad, dct_utils):
    """No class may fall through. A cell typed into a gap is a cell with no compartment."""
    vocab = set(dct_utils.list_supported_cell_types())
    mapped = set(broad["entries"])
    assert mapped == vocab, (
        f"unmapped classes: {sorted(vocab - mapped)}; stale entries: {sorted(mapped - vocab)}"
    )


def test_absent_entries_really_are_absent(broad, dct_utils):
    vocab = set(dct_utils.list_supported_cell_types())
    overlap = vocab & set(broad["absent_from_vocabulary"])
    assert not overlap, f"listed as absent but the model can emit them: {sorted(overlap)}"


def test_broad_classes_are_from_the_declared_set(broad):
    allowed = set(broad["broad_classes"])
    off = {k: e["broad_class"] for k, e in broad["entries"].items()
           if e["broad_class"] not in allowed}
    assert not off, f"broad_class outside the declared vocabulary: {off}"


def test_every_departure_from_the_lineage_default_is_argued(broad):
    """A broad_class that disagrees with its lineage's default must say why, in the file."""
    default = broad["lineage_to_broad"]
    for name, e in broad["entries"].items():
        if e["broad_class"] != default[e["lineage"]]:
            assert e.get("reason", "").strip(), (
                f"{name} overrides {e['lineage']}->{default[e['lineage']]} with "
                f"{e['broad_class']} and gives no reason"
            )


def test_erythrocytes_are_not_counted_as_immune(broad):
    """The study measures the immune environment. Upstream files red cells under Myeloid."""
    assert broad["entries"]["Erythrocyte"]["broad_class"] == "Other"


def test_the_langerhans_trap_stays_documented(broad):
    e = broad["entries"]["Langerhans"]
    assert e["broad_class"] == "Immune"
    assert "islet" in e["reason"].lower(), "the islet-of-Langerhans warning was dropped"


def test_expected_in_pancreas_is_flagged_report_only(broad):
    """If this ever becomes a filter, the policy text is the thing that says it must not."""
    assert "REPORT ONLY" in broad["policy"]["expected_in_pancreas"]
    assert all(isinstance(e["expected_in_pancreas"], bool) for e in broad["entries"].values())


def test_islet_endocrine_classes_are_expected_here(broad):
    for ct in ("AlphaCell", "BetaCell", "Endocrine"):
        e = broad["entries"][ct]
        assert e["expected_in_pancreas"] is True
        assert e["broad_class"] == "Epithelial"
