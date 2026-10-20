"""``stamp_if_values`` says what it did to every donor, including the ones it skipped.

The silence was the defect. The absent-arrays branch set NaN, ``continue``d and returned with
``markers`` empty and nothing printed, so a donor whose registration had PASSED and whose arrays
simply had not been produced yet was indistinguishable from one that was never registered --
which is how 6539 read as a registration failure for a week when all it needed was this step
re-running.

**Why these run in a subprocess.** This checkout's ``.venv`` has ``phenocycler`` installed
editable against a DIFFERENT worktree (``/home/smith6jt/wt-Q-ci/core/phenocycler``), so an
in-process ``import phenocycler`` in this suite does not reach the file under test. Each case
below runs with ``PYTHONPATH`` pointed at this tree and asserts, before anything else, that the
module it loaded is this one -- so the test cannot quietly pass against another checkout.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
LOCAL_PACKAGE = REPO / "core" / "multiplex5000"

PREAMBLE = textwrap.dedent(f'''
    import sys, types
    from pathlib import Path
    import numpy as np, pandas as pd
    import phenocycler.integration.postxen as P

    LOCAL = Path({str(LOCAL_PACKAGE)!r})
    assert Path(P.__file__).is_relative_to(LOCAL), (
        f"loaded {{P.__file__}}, which is not the file under test")

    class Cells:
        def __init__(self, df):
            self.df = df

    def fake_cfg(root):
        return types.SimpleNamespace(cells_xen_dir=Path(root))

    written = {{}}

    def install(monkey_cells=None, manifest_row=None):
        P.read_cell_table = lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError()) \\
            if monkey_cells is None else Cells(monkey_cells.copy())
        P.write_cell_table = lambda df, path: written.setdefault("path", str(path))
        P.partition_path = lambda base, donor, roi: Path(base) / f"donor_id={{donor}}"
        P.load_manifest = lambda cfg: pd.DataFrame()
        P.paired_rows = lambda man, roi=None, tissue=None: (
            pd.DataFrame() if manifest_row is None else pd.DataFrame([manifest_row]))
''').strip()


def _run(body: str, tmp_path: Path) -> subprocess.CompletedProcess:
    script = PREAMBLE + "\n" + textwrap.dedent(body)
    env = {**os.environ, "PYTHONPATH": str(LOCAL_PACKAGE),
           "PYTHONDONTWRITEBYTECODE": "1"}
    return subprocess.run([sys.executable, "-c", script], text=True, capture_output=True,
                          cwd=str(tmp_path), env=env)


def test_describe_stamp_names_every_outcome(tmp_path):
    """One line per donor, either way, and the skip line names the path it looked in."""
    result = _run('''
        stamped = P.describe_stamp({"donor_id": "6539", "roi": "panc",
                                    "markers": ["INS", "GCG"],
                                    "n_stamped_by_marker": {"INS": 319005, "GCG": 319005}})
        assert "6539/panc" in stamped and "INS 319,005" in stamped and "GCG 319,005" in stamped

        skipped = P.describe_stamp({
            "donor_id": "115", "roi": "panc", "markers": [],
            "skipped_reason": "no per-cell IF arrays for 115: 2 of 2 markers absent",
            "missing_arrays": ["/if/115_0041117_cell_INS.npy", "/if/115_0041117_cell_GCG.npy"]})
        assert "not stamped" in skipped
        assert "2 of 2 markers absent" in skipped
        assert "/if/115_0041117_cell_INS.npy" in skipped

        # The partial case a caller cannot see: markers is non-empty, one array is still absent.
        partial = P.describe_stamp({"donor_id": "6591", "roi": "panc", "markers": ["INS"],
                                    "n_stamped_by_marker": {"INS": 5},
                                    "missing_arrays": ["/if/6591_0041116_cell_GCG.npy"]})
        assert "stamped" in partial and "GCG absent" in partial
        assert "/if/6591_0041116_cell_GCG.npy" in partial

        errored = P.describe_stamp({"donor_id": "112", "roi": "panc", "markers": [],
                                    "error": "no Xenium geometry table",
                                    "looked_in": "/x/cells_xen"})
        assert "no Xenium geometry table" in errored and "/x/cells_xen" in errored
        print("OK")
    ''', tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_a_skipped_donor_prints_its_reason_and_the_path(tmp_path):
    result = _run(f'''
        root = Path({str(tmp_path)!r})
        cells = pd.DataFrame({{"cell_id": ["a", "b"]}})
        install(monkey_cells=cells, manifest_row=None)
        if_dir = root / "protein_rna"
        if_dir.mkdir(exist_ok=True)
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            out = P.stamp_if_values(fake_cfg(root), "6482", "panc", if_dir=if_dir,
                                    bundle_cells=root / "absent.parquet")
        line = buf.getvalue().strip()
        assert out["error"].startswith("bundle cells.parquet missing"), out
        assert line.startswith("[postxen] 6482/panc:"), line
        assert "absent.parquet" in line, line
        assert out["n_stamped"] == 0 and out["markers"] == []
        print("OK")
    ''', tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_announce_false_stays_silent_but_still_records(tmp_path):
    """A caller that prints its own line can turn this one off; the record is unchanged."""
    result = _run(f'''
        root = Path({str(tmp_path)!r})
        install(monkey_cells=pd.DataFrame({{"cell_id": ["a"]}}), manifest_row=None)
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            out = P.stamp_if_values(fake_cfg(root), "6482", "panc", if_dir=root / "if",
                                    bundle_cells=root / "absent.parquet", announce=False)
        assert buf.getvalue() == "", buf.getvalue()
        assert out["error"].startswith("bundle cells.parquet missing")
        print("OK")
    ''', tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_n_stamped_sums_every_marker_not_just_the_last(tmp_path):
    """``stamped`` was ASSIGNED inside the marker loop, so ``n_stamped`` reported GCG alone.

    The two markers are given DIFFERENT finite counts here, so a regression to "the last
    marker's count" would show as 3 rather than 5.
    """
    result = _run(f'''
        root = Path({str(tmp_path)!r})
        bundle = root / "bundle"; bundle.mkdir(exist_ok=True)
        order = np.array(["a", "b", "c"], dtype=object)
        pd.DataFrame({{"cell_id": order}}).to_parquet(bundle / "cells.parquet", index=False)
        if_dir = root / "if"; if_dir.mkdir(exist_ok=True)
        np.save(if_dir / "6539_0059865_cell_INS.npy", np.array([1.0, 2.0, np.nan]))
        np.save(if_dir / "6539_0059865_cell_GCG.npy", np.array([np.nan, 4.0, np.nan]))
        install(monkey_cells=pd.DataFrame({{"cell_id": ["c", "a", "b"]}}), manifest_row=None)
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            out = P.stamp_if_values(fake_cfg(root), "6539", "panc", if_dir=if_dir,
                                    bundle_cells=bundle / "cells.parquet")
        assert out["markers"] == ["INS", "GCG"], out
        assert out["n_stamped_by_marker"] == {{"INS": 2, "GCG": 1}}, out
        assert out["n_stamped"] == 3, out
        line = buf.getvalue().strip()
        assert "INS 2" in line and "GCG 1" in line, line
        print("OK")
    ''', tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_a_missing_array_is_nan_and_never_zero(tmp_path):
    """A marker whose array was never written is unavailable, not a dark cell."""
    result = _run(f'''
        root = Path({str(tmp_path)!r})
        bundle = root / "bundle"; bundle.mkdir(exist_ok=True)
        pd.DataFrame({{"cell_id": np.array(["a", "b"], dtype=object)}}).to_parquet(
            bundle / "cells.parquet", index=False)
        if_dir = root / "if"; if_dir.mkdir(exist_ok=True)
        np.save(if_dir / "6539_0059865_cell_INS.npy", np.array([1.0, 2.0]))
        holder = {{}}
        install(monkey_cells=pd.DataFrame({{"cell_id": ["a", "b"]}}), manifest_row=None)
        _write = P.write_cell_table
        def capture(df, path):
            holder["df"] = df.copy()
            return _write(df, path)
        P.write_cell_table = capture
        import io, contextlib
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            out = P.stamp_if_values(fake_cfg(root), "6539", "panc", if_dir=if_dir,
                                    bundle_cells=bundle / "cells.parquet")
        assert out["markers"] == ["INS"], out
        assert "missing_arrays" in out and out["missing_arrays"], out
        # The GLOB that found nothing is what is recorded -- not a path that exists, which is
        # the point: it tells the reader exactly where the array was expected.
        assert out["missing_arrays"][0].endswith("6539_*_cell_GCG.npy"), out["missing_arrays"]
        df = holder["df"]
        assert df["feat__if_GCG"].isna().all(), df["feat__if_GCG"].tolist()
        assert (df["feat__if_GCG"] == 0).sum() == 0
        assert "GCG absent" in buf.getvalue()
        print("OK")
    ''', tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
