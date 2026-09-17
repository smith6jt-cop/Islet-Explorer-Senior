"""``figure_manifest.json``: every input the report read, every file it wrote, each with a sha256.

    cd core/phenocycler && ../../.venv/bin/python scripts/report/manifest.py --verify

Two hash modes, both from ``phenocycler.artifacts`` so they are the pipeline's own:

* ``full`` -- ``sha256_file`` over every byte.  The default; used for every parquet, h5ad, h5,
  JSON and the GeoJSON exports (3-5 GB each; tens of seconds).
* ``sampled`` -- ``sha256_sampled_file``: three 1 MiB windows at deterministic offsets plus the
  file size, all folded into the digest.  Used for files above ``SAMPLED_ABOVE_BYTES`` (8 GiB) --
  in practice the two ~80 GB qptiffs -- because a full read of 160 GB per report run is not a
  cost worth paying for a figure lane.  **This is a weaker fingerprint** (a change confined to
  an unsampled region is not caught) and the manifest says so on the entry.  It is the same
  algorithm and the same parameters ``qupath_manifest.json`` records, so the two digests can be
  compared directly.

``verify`` re-hashes every listed input and output in its recorded mode and exits non-zero on
any drift, so a figure whose input moved cannot be mistaken for a current one.  Release ids are
recorded as ``null`` with the reason: the QC figures need no identity release.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from dataclasses import asdict, dataclass, field

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
from phenocycler.artifacts import FileFingerprint  # noqa: E402

SAMPLED_ABOVE_BYTES = 8 * 1024 ** 3
MANIFEST_NAME = "figure_manifest.json"


@dataclass
class FigureRecord:
    """What one figure is: its id, where it was written, what it read, and the numbers it printed."""
    figure_id: str
    donor: str                       # a donor id, or pins.PAIR_DIR for both-donor figures
    section: str                     # 'pheno' | 'xen' | '_pair'
    files: list                      # png, pdf, json sidecar
    inputs: list                     # paths read
    facts: dict = field(default_factory=dict)
    written_utc: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


def hash_mode_for(path: pathlib.Path, sampled_above: int = SAMPLED_ABOVE_BYTES) -> str:
    return "sampled" if path.stat().st_size > sampled_above else "full"


def fingerprint(path, *, sampled_above: int = SAMPLED_ABOVE_BYTES) -> dict:
    p = pathlib.Path(path).resolve()
    mode = hash_mode_for(p, sampled_above)
    fp = asdict(FileFingerprint.capture(p, hash_mode=mode))
    if mode == "sampled":
        fp["note"] = ("weaker fingerprint: 3 x 1 MiB windows + size (phenocycler.artifacts.sha256_sampled_file), "
                      "the same digest qupath_manifest.json records for this file; a change confined to an "
                      "unsampled region is not caught")
    return fp


class Manifest:
    """Accumulates inputs and figures across sections and donors, then writes one JSON."""

    def __init__(self, path: pathlib.Path | None = None, *, sampled_above: int = SAMPLED_ABOVE_BYTES,
                 provenance: dict | None = None):
        self.path = pathlib.Path(path) if path else pins.OUT / MANIFEST_NAME
        self.sampled_above = sampled_above
        self.doc = {"method_version": pins.METHOD_VERSION, "written_utc": "", "provenance": provenance or {},
                    "hash_modes": {"full": "sha256 over every byte",
                                   "sampled": f"files above {sampled_above} bytes; see note on each entry"},
                    "inputs": {}, "figures": {}}
        if self.path.exists():
            try:
                old = json.loads(self.path.read_text())
                self.doc["inputs"] = dict(old.get("inputs", {}))
                self.doc["figures"] = dict(old.get("figures", {}))
            except json.JSONDecodeError:
                pass

    def add_input(self, path) -> dict:
        """Fingerprint an input EVERY time it is read.

        The first version of this kept the first fingerprint it saw, which is wrong across a session
        that re-renders: a figure drawn after its input changed was recorded against the input's old
        hash, and ``--verify`` then reported drift on a manifest that was merely stale. Re-hashing on
        every read makes the manifest describe the inputs as they were at the LAST render, which is
        what ``--verify`` compares against. The cost is one hash per figure per input (sampled for
        files above ``sampled_above``).
        """
        p = pathlib.Path(path).resolve()
        self.doc["inputs"][str(p)] = fingerprint(p, sampled_above=self.sampled_above)
        return self.doc["inputs"][str(p)]

    def add_figure(self, rec: FigureRecord) -> None:
        for p in rec.inputs:
            self.add_input(p)
        entry = asdict(rec)
        entry["inputs"] = [str(pathlib.Path(p).resolve()) for p in rec.inputs]
        entry["files"] = {str(pathlib.Path(f).resolve()): fingerprint(f, sampled_above=self.sampled_above)
                          for f in rec.files}
        key = f"{rec.section}/{rec.donor}/{rec.figure_id}"
        self.doc["figures"][key] = entry

    def write(self) -> pathlib.Path:
        self.doc["written_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.doc, indent=1, default=str) + "\n")
        return self.path


def _recheck(entry: dict) -> tuple:
    p = pathlib.Path(entry["path"])
    if not p.exists():
        return False, "missing"
    fp = FileFingerprint.from_dict({k: entry[k] for k in ("path", "size_bytes", "mtime_ns", "sha256", "hash_mode",
                                                        "sample_size_bytes", "sample_count")})
    now = FileFingerprint.capture(p, hash_mode=fp.hash_mode, sample_size_bytes=fp.sample_size_bytes or 1 << 20,
                                  sample_count=fp.sample_count or 3)
    if now.sha256 != fp.sha256:
        return False, f"sha256 drift ({fp.hash_mode})"
    return True, "ok"


def verify(path: pathlib.Path | None = None, *, verbose: bool = True) -> tuple:
    """Re-hash every input and output.  Returns (ok, rows); ok is False on ANY drift or missing file."""
    path = pathlib.Path(path) if path else pins.OUT / MANIFEST_NAME
    if not path.exists():
        return False, [("manifest", str(path), "missing")]
    doc = json.loads(path.read_text())
    rows, ok = [], True
    for key, entry in doc.get("inputs", {}).items():
        good, why = _recheck(entry)
        ok &= good
        rows.append(("input", key, why))
    for fkey, fig in doc.get("figures", {}).items():
        for fpath, entry in fig.get("files", {}).items():
            good, why = _recheck(entry)
            ok &= good
            rows.append((f"figure {fkey}", fpath, why))
    if verbose:
        bad = [r for r in rows if r[2] != "ok"]
        print(f"figure manifest {path}: {len(rows)} files checked, {len(bad)} with drift/missing")
        for r in bad:
            print(f"  DRIFT {r[0]}: {r[1]} -- {r[2]}")
    return bool(ok), rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="figure manifest verify")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--manifest", default=None)
    a = ap.parse_args(argv)
    if a.verify:
        ok, _ = verify(pathlib.Path(a.manifest) if a.manifest else None)
        return 0 if ok else 1
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
