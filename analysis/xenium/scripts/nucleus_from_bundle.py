#!/usr/bin/env python3
"""Find a Xenium section's ``nucleus_count`` by SEARCHING the disk, not by trusting a path.

A cell with no nucleus is excluded from every calculation.  The pancreas h5ads under
``/data/xenium/h5ad_out`` carry no nucleus column at all, so the only handle on the rule is the
vendor bundle's ``cells.parquet``.  ``panc_xenium_cluster_name.nucleus_for`` looks that bundle up
through ``core/phenocycler/data/integration/xenium_paths.csv`` and gives up when the recorded
``source_path`` does not exist.  Recorded paths go stale: serial numbers stay put but tissue
tokens and directory spellings do not (``PanT`` for 115's pancreas, ``Region_1`` for 6591's,
``Panc_pLN`` for a section that carries both ROIs), and whole bundles get moved between
``/data/xenium/HuPrime5kPanTissue/``, ``/data/xenium_pln/`` and ``/data/xenium_bundles/``.

So this module never reads ``source_path`` as an instruction.  It indexes every
``cells.parquet`` it can find under the search roots, reads each bundle's own
``metrics_summary.csv`` for the cassette (donor) and region the vendor recorded in it, ranks the
candidates, and then decides by CONTENT: a bundle is this section's bundle only if its
``cell_id`` set is a superset of the h5ad's own ``cell_id`` set.  Nothing else is accepted.

That last rule is not ceremony.  Xenium ``cell_id`` strings are per-region, so two unrelated
sections share a few by coincidence: 6347's pancreas h5ad (1,393,487 cells) meets exactly 28 of
its ids in the 0041140 pLN bundle and 6580's (931,227) meets 31 in the 0041139 pLN bundle.  A
resolver that matched on "some ids are present" would hand back a nucleus rule built from 28
cells.  Containment is the only test that separates the right bundle from a coincidence.

When no bundle contains the section, the result says so and carries what was tried, plus any
section the vendor metrics describe but no bundle on this machine holds -- an absent measurement
is unavailable, never zero.

Use
---
    from nucleus_from_bundle import resolve_nucleus, nucleus_counts

    res = resolve_nucleus("6482", "panc")        # cell ids loaded from the h5ad
    res.status, res.bundle, res.n_zero           # 'resolved', <path>, 23,847
    res.mapping["aaaaaajp-1"]                    # 1

    counts, source = nucleus_counts("6482", "panc", cell_ids)   # drop-in for nucleus_for()

    $ python nucleus_from_bundle.py 6347 6579 6580 115 6482

Nothing here writes to the release, and nothing here decides what the release does with an
unresolved section.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Roots searched, in order.  The narrow ones first because a bundle found there is found in
# seconds; ``/data`` is the exhaustive sweep and is only walked when the narrow roots fail.
NARROW_ROOTS: tuple[Path, ...] = (
    Path("/data/xenium"),
    Path("/data/xenium_pln"),
    Path("/data/xenium_bundles"),
    Path("/data/QP_Xenium"),
)
WIDE_ROOTS: tuple[Path, ...] = (Path("/data"),)
H5AD_DIR = Path("/data/xenium/h5ad_out")

# Directories that cannot hold a vendor bundle and cost real time to walk.
PRUNE = {
    ".git", "__pycache__", ".cache", "node_modules", ".snapshots", "lost+found",
    "paraview_build", "alphafold", "alphafold_models", ".venv", "renv", "site-packages",
}
MAX_DEPTH = 8

# output-XETG00298__0041140__Panc_pLN__20241025__204655
BUNDLE_RE = re.compile(
    r"^output-(?P<instrument>[A-Za-z0-9]+)__(?P<serial>\d+)__(?P<region>.+?)__"
    r"(?P<date>\d{8})__(?P<time>\d+)$"
)

PANC_TOKENS = ("panc", "pan", "pant", "pancreas")
PLN_TOKENS = ("pln", "ln", "node")


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def _roi_family(roi: str) -> str:
    """``panc`` / ``pln`` / ``''`` -- the coarse ROI a token belongs to."""
    n = _norm(roi)
    if not n:
        return ""
    if "pln" in n or n.startswith("ln"):
        return "pln"
    if n.startswith("pan"):
        return "panc"
    return ""


def _region_families(region: str) -> set[str]:
    """Every ROI family a bundle's region token could serve.

    ``Panc_pLN`` is one section carrying both, and it is scored as both rather than guessed at.
    ``Region_1`` names no tissue, so it belongs to neither and survives on other evidence.
    """
    n = _norm(region)
    fams = set()
    if any(t in n for t in PLN_TOKENS):
        fams.add("pln")
    if any(n.startswith(t) or t in n for t in PANC_TOKENS):
        fams.add("panc")
    return fams


@dataclass
class Bundle:
    """One directory on disk that holds a ``cells.parquet``."""

    path: Path                      # the bundle directory
    cells: Path                     # its cells.parquet
    serial: str = ""
    region: str = ""
    run_date: str = ""
    cassette: str = ""              # vendor's cassette_name from metrics_summary.csv
    metrics_region: str = ""        # vendor's region_name
    num_cells_detected: int | None = None
    n_rows: int | None = None       # parquet row count (metadata only; no data read)
    has_nucleus_count: bool | None = None

    def as_dict(self) -> dict:
        d = {k: (str(v) if isinstance(v, Path) else v) for k, v in self.__dict__.items()}
        return d


@dataclass
class Section:
    """A section the vendor metrics describe, whether or not its bundle is on this machine."""

    metrics: Path
    serial: str = ""
    region: str = ""
    cassette: str = ""
    run_name: str = ""
    num_cells_detected: int | None = None
    bundle: Path | None = None      # the directory holding it, when one was found


@dataclass
class NucleusResolution:
    donor: str
    roi: str
    status: str                     # 'resolved' | 'unresolved'
    reason: str = ""
    bundle: str | None = None
    serial: str = ""
    region: str = ""
    n_query: int = 0
    n_bundle_rows: int = 0
    n_matched: int = 0
    n_missing: int = 0
    n_zero: int = 0
    counts: dict = field(default_factory=dict)          # nucleus_count value -> n cells
    mapping: dict = field(default_factory=dict)         # cell_id -> nucleus_count
    candidates: list = field(default_factory=list)      # what was tried, and why it was rejected
    orphan_sections: list = field(default_factory=list) # measured sections with no bundle here

    def report(self) -> str:
        head = f"{self.donor} {self.roi}: {self.status}"
        if self.status == "resolved":
            share = self.n_zero / self.n_query if self.n_query else float("nan")
            dist = ", ".join(f"{k}={v:,}" for k, v in sorted(self.counts.items()))
            return (
                f"{head}\n"
                f"  bundle          {self.bundle}\n"
                f"  serial/region   {self.serial} / {self.region}\n"
                f"  h5ad cells      {self.n_query:,}\n"
                f"  bundle rows     {self.n_bundle_rows:,}\n"
                f"  matched         {self.n_matched:,} ({self.n_missing:,} unmatched)\n"
                f"  nucleus_count   {dist}\n"
                f"  nucleus_count==0 {self.n_zero:,} cells -> excluded "
                f"({share:.4%} of the section)"
            )
        lines = [f"{head}\n  {self.reason}"]
        for c in self.candidates:
            lines.append(
                f"  tried {c['path']}  rows={c.get('n_rows')}  "
                f"ids_in_common={c.get('intersection')}  -> {c['verdict']}"
            )
        for s in self.orphan_sections:
            lines.append(
                f"  vendor metrics describe {s['serial']}__{s['region']} "
                f"(cassette {s['cassette']}, {s['num_cells_detected']:,} cells detected) but no "
                f"cells.parquet on this machine holds it: {s['metrics']}"
            )
        return "\n".join(lines)


def _walk_bundles(root: Path) -> list[Path]:
    """Every directory under ``root`` that holds a ``cells.parquet``.

    A directory that holds one is not descended into: a bundle's own subdirectories
    (``analysis``, ``cell_feature_matrix``, ``morphology_focus``) never hold another section.
    """
    found: list[Path] = []
    root = Path(root)
    if not root.is_dir():
        return found
    base_depth = len(root.parts)
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        here = Path(dirpath)
        if len(here.parts) - base_depth >= MAX_DEPTH:
            dirnames[:] = []
            continue
        dirnames[:] = [d for d in dirnames if d not in PRUNE and not d.startswith(".")]
        if "cells.parquet" in filenames:
            found.append(here)
            dirnames[:] = []
    return found


def _read_metrics(path: Path) -> dict:
    try:
        with path.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return {}
    return rows[0] if rows else {}


def _parquet_facts(cells: Path) -> tuple[int | None, bool | None]:
    """Row count and whether ``nucleus_count`` is a column -- parquet footer only."""
    try:
        import pyarrow.parquet as pq

        pf = pq.ParquetFile(cells)
        return pf.metadata.num_rows, "nucleus_count" in pf.schema_arrow.names
    except Exception:
        return None, None


def index_bundles(
    roots: tuple[Path, ...] = NARROW_ROOTS,
    cache: Path | None = None,
    refresh: bool = False,
) -> list[Bundle]:
    """Index every bundle under ``roots``.  Optionally cached, since the sweep costs minutes."""
    if cache is not None and cache.is_file() and not refresh:
        try:
            raw = json.loads(cache.read_text())
            return [
                Bundle(
                    path=Path(b["path"]), cells=Path(b["cells"]), serial=b["serial"],
                    region=b["region"], run_date=b["run_date"], cassette=b["cassette"],
                    metrics_region=b["metrics_region"],
                    num_cells_detected=b["num_cells_detected"], n_rows=b["n_rows"],
                    has_nucleus_count=b["has_nucleus_count"],
                )
                for b in raw
            ]
        except Exception:
            pass

    out: list[Bundle] = []
    seen: set[Path] = set()
    for root in roots:
        for d in _walk_bundles(Path(root)):
            real = d.resolve()
            if real in seen:
                continue
            seen.add(real)
            b = Bundle(path=d, cells=d / "cells.parquet")
            m = BUNDLE_RE.match(d.name)
            if m:
                b.serial = m.group("serial")
                b.region = m.group("region")
                b.run_date = m.group("date")
            mx = _read_metrics(d / "metrics_summary.csv")
            if mx:
                b.cassette = str(mx.get("cassette_name", "") or "")
                b.metrics_region = str(mx.get("region_name", "") or "")
                n = mx.get("num_cells_detected")
                b.num_cells_detected = int(float(n)) if n not in (None, "",) else None
            b.n_rows, b.has_nucleus_count = _parquet_facts(b.cells)
            out.append(b)
    if cache is not None:
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps([b.as_dict() for b in out], indent=1))
        except OSError:
            pass
    return out


def index_sections(roots: tuple[Path, ...] = NARROW_ROOTS) -> list[Section]:
    """Every section a ``metrics_summary.csv`` describes, bundle present or not.

    Vendor metrics survive a deleted bundle: ``/data/xenium_pln/_metrics/`` holds one per
    section, including sections whose ``cells.parquet`` is nowhere on this machine.  They are
    what lets an unresolved section be reported as a specific absence rather than a shrug.
    """
    out: list[Section] = []
    seen: set[Path] = set()
    for root in roots:
        root = Path(root)
        if not root.is_dir():
            continue
        base_depth = len(root.parts)
        for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
            here = Path(dirpath)
            if len(here.parts) - base_depth >= MAX_DEPTH:
                dirnames[:] = []
                continue
            dirnames[:] = [d for d in dirnames if d not in PRUNE and not d.startswith(".")]
            for name in filenames:
                if not name.endswith("metrics_summary.csv"):
                    continue
                p = here / name
                if p.resolve() in seen:
                    continue
                seen.add(p.resolve())
                mx = _read_metrics(p)
                if not mx:
                    continue
                serial, region = "", ""
                m = BUNDLE_RE.match(name.replace(".metrics_summary.csv", ""))
                if m:
                    serial, region = m.group("serial"), m.group("region")
                else:
                    m = BUNDLE_RE.match(here.name)
                    if m:
                        serial, region = m.group("serial"), m.group("region")
                n = mx.get("num_cells_detected")
                out.append(Section(
                    metrics=p, serial=serial,
                    region=region or str(mx.get("region_name", "") or ""),
                    cassette=str(mx.get("cassette_name", "") or ""),
                    run_name=str(mx.get("run_name", "") or ""),
                    num_cells_detected=int(float(n)) if n not in (None, "") else None,
                    bundle=here if (here / "cells.parquet").is_file() else None,
                ))
    return out


def _donor_evidence(b: Bundle, donor: str) -> int:
    """How strongly a bundle claims this donor, from its own name and its own metrics."""
    d = _norm(donor)
    score = 0
    if d and d in _norm(b.cassette):
        score += 3                         # the vendor's cassette_name, e.g. nPOD6347 / 6579
    if d and d in _norm(b.path.name):
        score += 2                         # the directory spelling, e.g. 6476_Pan
    return score


def _rank(bundles: list[Bundle], donor: str, roi: str, n_query: int) -> list[Bundle]:
    fam = _roi_family(roi)
    ranked = []
    for b in bundles:
        score = _donor_evidence(b, donor)
        fams = _region_families(b.region) | _region_families(b.metrics_region)
        if fam and fam in fams:
            score += 2
        elif fam and fams:
            score -= 1                     # the other tissue of the same donor
        if b.n_rows is not None and n_query:
            # A superset has at least as many rows, and the vendor's QC drops only a little,
            # so the right bundle sits just above the h5ad's own cell count.
            if b.n_rows >= n_query:
                score += 2
                slack = (b.n_rows - n_query) / n_query
                if slack < 0.10:
                    score += 2
            else:
                score -= 3
        ranked.append((score, -(b.n_rows or 0), b))
    ranked.sort(key=lambda t: (-t[0], t[1]))
    return [b for _, _, b in ranked]


def find_h5ads(donor: str, roi: str = "panc", h5ad_dir: Path = H5AD_DIR) -> list[Path]:
    """Every h5ad under ``h5ad_dir`` that is this donor's ``roi`` -- searched, not assumed.

    Donor 115 and 112 are written ``HDL115`` / ``HDL112``; a donor with several sections of one
    ROI carries a trailing letter (``6579a_pLN``, ``6579b_pLN``, ``6579c_pLN``).  Returned in
    name order so a caller that resolves them one at a time gets a stable order.
    """
    tissue = "pLN" if _roi_family(roi) == "pln" else "Panc"
    pats = [f"{donor}_{tissue}Xenium.h5ad", f"HDL{donor}_{tissue}Xenium.h5ad",
            f"{donor}?_{tissue}Xenium.h5ad", f"HDL{donor}?_{tissue}Xenium.h5ad"]
    hits: list[Path] = []
    for pat in pats:
        for h in sorted(Path(h5ad_dir).glob(pat)):
            if h not in hits:
                hits.append(h)
    return hits


def read_h5ad_cell_ids(path: Path) -> np.ndarray:
    """One h5ad's ``obs`` index, as strings."""
    import h5py

    with h5py.File(path, "r") as f:
        obs = f["obs"]
        key = obs.attrs.get("_index", "_index")
        key = key.decode() if isinstance(key, bytes) else key
        return obs[key][:].astype(str)


def h5ad_cell_ids(donor: str, roi: str = "panc", h5ad_dir: Path = H5AD_DIR) -> np.ndarray:
    """The section's own cell ids, when this (donor, roi) is a single section."""
    hits = find_h5ads(donor, roi, h5ad_dir)
    tissue = "pLN" if _roi_family(roi) == "pln" else "Panc"
    if not hits:
        raise FileNotFoundError(f"no {tissue} h5ad for donor {donor} under {h5ad_dir}")
    if len(hits) > 1:
        raise ValueError(
            f"donor {donor} has {len(hits)} {tissue} h5ads ({', '.join(h.name for h in hits)}); "
            "resolve one section at a time by passing cell_ids")
    return read_h5ad_cell_ids(hits[0])


def resolve_nucleus(
    donor: str,
    roi: str = "panc",
    cell_ids: np.ndarray | None = None,
    roots: tuple[Path, ...] | None = None,
    widen: bool = True,
    cache: Path | None = None,
) -> NucleusResolution:
    """``cell_id -> nucleus_count`` for one (donor, roi), found by searching.

    A candidate is accepted only when it contains EVERY id of the section.  Partial overlap is
    reported, never used: Xenium cell ids collide across sections by coincidence.
    """
    import pandas as pd

    if cell_ids is None:
        cell_ids = h5ad_cell_ids(donor, roi)
    cell_ids = np.asarray(cell_ids).astype(str)
    want = set(cell_ids.tolist())
    res = NucleusResolution(donor=str(donor), roi=roi, status="unresolved",
                            n_query=len(cell_ids))

    root_sets = [tuple(roots)] if roots else [NARROW_ROOTS] + ([WIDE_ROOTS] if widen else [])
    tried_paths: set[str] = set()
    for root_set in root_sets:
        bundles = index_bundles(root_set, cache=cache if root_set == NARROW_ROOTS else None)
        for b in _rank(bundles, str(donor), roi, len(cell_ids)):
            if str(b.cells) in tried_paths:
                continue
            if b.has_nucleus_count is False:
                tried_paths.add(str(b.cells))
                res.candidates.append({"path": str(b.cells), "n_rows": b.n_rows,
                                       "intersection": None,
                                       "verdict": "no nucleus_count column"})
                continue
            if b.n_rows is not None and b.n_rows < len(cell_ids):
                # Cannot contain the section: a superset has at least as many rows.  Recorded
                # rather than silently skipped, and not read.
                tried_paths.add(str(b.cells))
                res.candidates.append({"path": str(b.cells), "n_rows": b.n_rows,
                                       "intersection": None,
                                       "verdict": "too few rows to contain this section"})
                continue
            tried_paths.add(str(b.cells))
            try:
                frame = pd.read_parquet(b.cells, columns=["cell_id", "nucleus_count"])
            except Exception as exc:
                res.candidates.append({"path": str(b.cells), "n_rows": b.n_rows,
                                       "intersection": None,
                                       "verdict": f"unreadable: {exc}"})
                continue
            frame["cell_id"] = frame["cell_id"].astype(str)
            have = dict(zip(frame["cell_id"], frame["nucleus_count"].to_numpy()))
            missing = want - have.keys()
            if missing:
                res.candidates.append({
                    "path": str(b.cells), "n_rows": b.n_rows,
                    "intersection": len(want) - len(missing),
                    "verdict": f"not a superset ({len(missing):,} of {len(want):,} ids absent)"})
                continue
            counts = np.array([have[c] for c in cell_ids])
            res.status = "resolved"
            res.bundle = str(b.cells)
            res.serial, res.region = b.serial, b.region or b.metrics_region
            res.n_bundle_rows = int(b.n_rows or len(frame))
            res.n_matched = int(len(cell_ids))
            res.n_missing = 0
            res.counts = {int(k): int(v) for k, v in
                          sorted(Counter(counts.tolist()).items())}
            res.n_zero = int((counts == 0).sum())
            res.mapping = have
            res.candidates.append({"path": str(b.cells), "n_rows": b.n_rows,
                                   "intersection": len(want), "verdict": "superset -- accepted"})
            return res

    res.reason = (
        "UNMEASURED: no cells.parquet under "
        f"{', '.join(str(r) for p in root_sets for r in p)} contains every cell id of this "
        "section, so the nucleus rule is unavailable for it -- not passed, and not zero")
    for s in index_sections(NARROW_ROOTS):
        if _donor_evidence(Bundle(path=s.metrics, cells=s.metrics, cassette=s.cassette),
                           str(donor)) == 0:
            continue
        if s.bundle is not None:
            continue
        if _roi_family(roi) and _region_families(s.region) and \
                _roi_family(roi) not in _region_families(s.region):
            continue
        res.orphan_sections.append({
            "metrics": str(s.metrics), "serial": s.serial, "region": s.region,
            "cassette": s.cassette, "run_name": s.run_name,
            "num_cells_detected": s.num_cells_detected})
    return res


def nucleus_counts(
    donor: str, roi: str, cell_ids: np.ndarray, **kw
) -> tuple[np.ndarray | None, str]:
    """``(counts aligned to cell_ids, source)`` or ``(None, why-not)``.

    Same shape as ``panc_xenium_cluster_name.nucleus_for`` so a caller can swap the recorded-path
    lookup for a search without changing what it does with the answer.
    """
    res = resolve_nucleus(donor, roi, cell_ids=cell_ids, **kw)
    if res.status != "resolved":
        return None, res.reason
    counts = np.array([res.mapping[c] for c in np.asarray(cell_ids).astype(str)],
                      dtype=np.int32)
    return counts, f"bundle {res.bundle}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("donors", nargs="+", help="donor ids, e.g. 6347 6579 6580 115 6482")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--no-widen", action="store_true",
                    help="search only the Xenium roots, never the whole of /data")
    ap.add_argument("--cache", type=Path, default=None,
                    help="JSON index cache; the narrow sweep is reused from it")
    ap.add_argument("--json", type=Path, default=None, help="write the results here")
    args = ap.parse_args(argv)

    out = []
    for donor in args.donors:
        sections = find_h5ads(donor, args.roi)
        if not sections:
            # An h5ad this machine does not hold is an absence to name, not a crash: the CLI
            # keeps going so one missing donor does not hide the rest of the sweep.
            tissue = "pLN" if _roi_family(args.roi) == "pln" else "Panc"
            res = NucleusResolution(
                donor=str(donor), roi=args.roi, status="unresolved",
                reason=(f"UNMEASURED: no {tissue} h5ad for donor {donor} under {H5AD_DIR}, so "
                        "there is no section to resolve a nucleus rule for here"))
            print(res.report())
            print()
            out.append({k: v for k, v in res.__dict__.items() if k != "mapping"})
            continue
        for h in sections:
            # Several h5ads mean several sections of the same ROI (6579a/b/c_pLN); each is
            # resolved against its own cell ids rather than pooled.
            label = args.roi if len(sections) == 1 else f"{args.roi} [{h.stem}]"
            res = resolve_nucleus(donor, args.roi, cell_ids=read_h5ad_cell_ids(h),
                                  widen=not args.no_widen, cache=args.cache)
            res.roi = label
            print(res.report())
            print()
            d = {k: v for k, v in res.__dict__.items() if k != "mapping"}
            d["h5ad"] = str(h)
            out.append(d)
    if args.json:
        args.json.write_text(json.dumps(out, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
