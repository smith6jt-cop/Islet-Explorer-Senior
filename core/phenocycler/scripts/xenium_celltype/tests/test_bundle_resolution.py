"""Which physical section a donor means — the bug that produced no error at all.

`xenium_paths.csv` has one row per (donor, roi, batch), and two donors carry more than one
row for `roi=panc`: 6533 and 112 were run in the 2024 nPOD cohort (serials 0041393 / 0041396)
and again in the 2026 Pheno-Xeno acquisition (0080870 / 0081366). Three separate call sites
selected with `.iloc[0]`, which is the 2024 row for both. A Pheno-Xeno analysis would have
run on the wrong physical section, completed normally, and reported plausible numbers.

These tests pin the two halves of the fix — the ambiguous case FAILS, and a named serial
resolves — plus the three things that had to be true around it for the fix to be worth
anything: the 2026 paths are absolute and must not be rewritten by a site root, the h5ad
cannot be paired with a section it has no way to name, and the R step gets R_LIBS_USER.
"""
import pathlib
import sys

import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
LANE = HERE.parent
sys.path.insert(0, str(LANE))
sys.path.insert(0, str(LANE / "hpg"))

import xenio                                          # noqa: E402

#: The real cohort table. These tests assert against the actual recorded rows, not a fixture:
#: the defect was a property of THIS file's contents, and a fixture that happened to hold one
#: row per donor would have passed against the broken code.
PATHS_CSV = LANE.parents[1] / "data" / "integration" / "xenium_paths.csv"

#: (donor, the 2024 nPOD serial `.iloc[0]` used to return, the 2026 Pheno-Xeno serial wanted)
PHENO_XENO = [("6533", "0041393", "0080870"), ("112", "0041396", "0081366")]


def test_the_two_pheno_xeno_donors_really_are_ambiguous():
    """If this fails the premise is gone and the rest of the file is testing nothing."""
    for donor, old, new in PHENO_XENO:
        serials = list(xenio.bundle_candidates(donor, PATHS_CSV, "panc").serial)
        assert sorted(serials) == sorted([old, new]), f"{donor}: {serials}"


@pytest.mark.parametrize("donor,old,new", PHENO_XENO)
def test_ambiguous_donor_fails_loud_instead_of_picking(donor, old, new):
    """The whole point. A silent wrong pick is the danger; an exception is the fix."""
    with pytest.raises(xenio.BundleAmbiguityError) as e:
        xenio.resolve_bundle(donor, PATHS_CSV, roi="panc")
    msg = str(e.value)
    # the error has to be actionable: both serials, and the flag that resolves it
    assert old in msg and new in msg, msg
    assert "--serial" in msg, msg


@pytest.mark.parametrize("donor,old,new", PHENO_XENO)
def test_naming_the_serial_selects_that_section(donor, old, new):
    """Both directions, because 'it returns the new one' would also pass if it ignored us."""
    assert xenio.resolve_bundle(donor, PATHS_CSV, roi="panc", serial=new).serial == new
    assert xenio.resolve_bundle(donor, PATHS_CSV, roi="panc", serial=old).serial == old

    ref = xenio.resolve_bundle(donor, PATHS_CSV, roi="panc", serial=new)
    assert ref.batch == "PhenoXeno"
    assert ref.n_candidates == 2
    # the 2026 bundle, not the 2024 one, is what the path points at
    assert new in str(ref.path) and old not in str(ref.path)


def test_a_whole_bundle_directory_name_is_accepted_as_a_serial():
    """Copy-pasting the bundle name is the obvious thing to try; it should not be an error."""
    ref = xenio.resolve_bundle(
        "6533", PATHS_CSV, roi="panc",
        serial="output-XETG00298__0080870__panc__20260910__211917")
    assert ref.serial == "0080870"


def test_an_unrecorded_serial_names_the_ones_that_exist():
    with pytest.raises(xenio.XeniumIOError) as e:
        xenio.resolve_bundle("6533", PATHS_CSV, roi="panc", serial="9999999")
    assert "0041393" in str(e.value) and "0080870" in str(e.value)


def test_unambiguous_donors_still_resolve_without_a_serial():
    """Failing loud must cost nothing for the donors that are not dual.

    The dual set is pinned rather than derived, so a donor QUIETLY becoming dual-bundle shows
    up here as a failure instead of silently widening where a serial is mandatory. It grew
    from two pairs to four on 2026-09-17, when the same-slide Pheno-Xeno pLN bundles
    (0080870__pLN / 0081366__pLN) were recorded in xenium_paths.csv beside the 2024 nPOD pLN
    runs. That is the same 2024-vs-2026 duality panc already had, now for pln_1.
    """
    xen = pd.read_csv(PATHS_CSV, dtype=str)
    xen = xen[xen.source_path.notna() & (xen.source_path != "NA")]
    dual = {(d, r) for (d, r), n in xen.groupby(["donor_id", "roi"]).size().items() if n > 1}
    assert dual == {("6533", "panc"), ("112", "panc"),
                    ("6533", "pln_1"), ("112", "pln_1")}, dual
    for (donor, roi), _ in xen.groupby(["donor_id", "roi"]):
        if (donor, roi) in dual:
            continue
        assert xenio.resolve_bundle(donor, PATHS_CSV, roi=roi).serial


def test_an_absolute_recorded_path_is_not_rewritten_by_a_root():
    """The 2026 bundles are recorded absolute; joining them to a site root invents a path."""
    abs_p = "/blue/maigan/smith6jt/Xenium_Analysis/data/output-XETG00298__0080870__panc__x"
    assert xenio.join_root(pathlib.Path("/data/xenium"), abs_p) == pathlib.Path(abs_p)
    # a relative path still gets the root, which is the case every nPOD donor relies on
    assert xenio.join_root(pathlib.Path("/data/xenium"), "HuPrime5k/x") == \
        pathlib.Path("/data/xenium/HuPrime5k/x")


def test_sites_does_not_strip_components_off_an_absolute_path():
    """`hpg/sites.py` rewrites recorded paths with root+strip; that must not touch an
    absolute one. It used to produce /orange/brusko/xenium/blue/maigan/... — nowhere."""
    import sites
    cfg = {"bundle_roots": [{"root": "/orange/brusko/xenium", "strip": 0},
                            {"root": "/orange/maigan/Xenium_T1D_PancPLN", "strip": 2}]}
    abs_p = "/blue/maigan/smith6jt/Xenium_Analysis/data/output-XETG00298__0080870__panc__x"
    assert sites._candidates(abs_p, cfg) == [pathlib.Path(abs_p)]
    rel = sites._candidates("HuPrime5kPanTissue/6KTRAF/MM036/output-XETG00298__0041393__P__x",
                            cfg)
    assert pathlib.Path("/orange/brusko/xenium/HuPrime5kPanTissue/6KTRAF/MM036/"
                        "output-XETG00298__0041393__P__x") in rel


@pytest.mark.parametrize("donor,old,new", PHENO_XENO)
def test_the_h5ad_is_refused_for_a_multi_section_donor(donor, old, new):
    """The hazard a serial-only fix would have left behind, and made worse.

    `<donor>_PancXenium.h5ad` is keyed on donor and roi with no serial, and BOTH these
    donors have one — they are 2024 cohort donors. Preferring it would take expression from
    the 2024 section and coordinates from the 2026 bundle and join them on `cell_id`. The
    ids overlap, so that join returns cells rather than raising: it would look ordinary.
    """
    reason = xenio.h5ad_refusal(donor, "panc", n_candidates=2)
    assert reason, "an h5ad with no serial must not be paired with one of two sections"
    assert "cell_id" in reason
    # and it must NOT fire for the ordinary single-section donors
    assert xenio.h5ad_refusal("6539", "panc", n_candidates=1) == ""


def test_serial_spec_grammar():
    assert xenio.parse_serial_spec(["0080870"]) == {"*": "0080870"}
    assert xenio.parse_serial_spec(["6533=0080870", "112=0081366"]) == {
        "6533": "0080870", "112": "0081366"}
    assert xenio.serial_for("6533", {"6533": "0080870"}) == "0080870"
    assert xenio.serial_for("6533", {"*": "0080870"}) == "0080870"
    assert xenio.serial_for("6539", {"6533": "0080870"}) is None
    # a bare serial silently applied across a cohort is exactly the class of bug being fixed
    with pytest.raises(xenio.XeniumIOError):
        xenio.parse_serial_spec(["0080870", "112=0081366"])
    with pytest.raises(xenio.XeniumIOError):
        xenio.parse_serial_spec(["6533=0080870", "6533=0041393"])


def test_build_one_resolves_the_serial_once_for_path_and_cells(monkeypatch):
    """The output path's serial and the cells written under it must come from ONE resolution.

    They were resolved independently — two `.iloc[0]` selections over the same CSV — so they
    agreed only by coincidence. This asserts the serial `load_donor` is asked for is the one
    the output path is named with.
    """
    import build_phenotyped_base as B

    seen = {}

    def fake_load_donor(donor, csv, *, roi, prefer, serial=None):
        seen["serial"] = serial
        raise RuntimeError("stop after resolution")

    monkeypatch.setattr(B.xenio, "load_donor", fake_load_donor)
    with pytest.raises(RuntimeError, match="stop after resolution"):
        B.build_one("6533", roi="panc", out_root=pathlib.Path("/tmp/nope"),
                    paths_csv=PATHS_CSV, serial="0080870")
    assert seen["serial"] == "0080870"


def test_build_one_refuses_an_ambiguous_donor():
    import build_phenotyped_base as B
    with pytest.raises(xenio.BundleAmbiguityError):
        B.build_one("6533", roi="panc", out_root=pathlib.Path("/tmp/nope"),
                    paths_csv=PATHS_CSV)


def test_preflight_reports_every_unresolvable_donor_before_building_any():
    """A cohort build is hours long; the refusal has to come in the first second."""
    import build_phenotyped_base as B
    with pytest.raises(B.BuildError) as e:
        B.preflight_serials(["6539", "6533", "112"], roi="panc",
                            paths_csv=PATHS_CSV, serials={})
    msg = str(e.value)
    assert "6533" in msg and "112" in msg
    assert "nothing was\nbuilt" in msg or "nothing was built" in msg.replace("\n", " ")
    # naming them makes it pass
    B.preflight_serials(["6539", "6533", "112"], roi="panc", paths_csv=PATHS_CSV,
                        serials={"6533": "0080870", "112": "0081366"})


@pytest.fixture
def fake_rscript(tmp_path):
    """A PATH carrying an `Rscript`, so these tests check the env, not this machine's R."""
    d = tmp_path / "bin"
    d.mkdir()
    exe = d / "Rscript"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    return str(d)


def test_rscript_is_given_R_LIBS_USER(fake_rscript):
    """KODAMA lives in R_LIBS_USER. `run_donor` spawned Rscript with the ambient environment,
    so an interactive run failed with `there is no package called 'KODAMA'` — after the PCA.

    The default is asserted against the path `hpg/env.sh` and `hpg/setup_envs.sh` install to,
    which is where `hpg/kodama.lock` records KODAMA 3.3 actually living.
    """
    import run_donor
    env = run_donor.r_environment({"PATH": fake_rscript})
    assert env["R_LIBS_USER"] == "/blue/maigan/smith6jt/R_libs"
    # an explicit one wins, so a caller that has sourced hpg/env.sh is not overridden
    env = run_donor.r_environment({"PATH": fake_rscript, "R_LIBS_USER": str(pathlib.Path.cwd())})
    assert env["R_LIBS_USER"] == str(pathlib.Path.cwd())


def test_r_environment_refuses_a_library_path_that_does_not_exist(fake_rscript):
    import run_donor
    with pytest.raises(RuntimeError, match="R_LIBS_USER"):
        run_donor.r_environment({"PATH": fake_rscript, "R_LIBS_USER": "/nonexistent/R_libs"})


def test_r_environment_refuses_when_Rscript_is_not_on_path():
    """R comes from the conda env; a missing interpreter should say so, not FileNotFoundError."""
    import run_donor
    with pytest.raises(RuntimeError, match="Rscript is not on PATH"):
        run_donor.r_environment({"PATH": "/nonexistent/bin",
                                 "R_LIBS_USER": "/blue/maigan/smith6jt/R_libs"})
