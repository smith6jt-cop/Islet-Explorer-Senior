"""Resolve a donor's Xenium bundle and h5ad at whichever site this process is running on.

``xenio.py`` hard-codes ``/data/xenium``, which is true of this workstation and of nowhere
else. On HiPerGator the same cohort is split across two ``/orange`` trees with two different
layouts, and one of them -- ``/orange/brusko`` -- is the only copy of MM046/MM047, i.e. of
donors 6450, 6505, 6534 and 6623. A resolver that quietly skipped an unreadable donor would
drop two of the eleven Xenium runs and two of the cohort's three Aab+ donors, and the run
would still exit 0. So this module fails closed: a donor either resolves to a named file that
exists, or it raises with every path it tried.

``sites.json`` holds the layouts. ``XENIUM_SITE`` picks one; with nothing set, the site is
chosen by which one actually resolves, and the choice is recorded, never assumed.

Two things a site table cannot decide, both delegated to :mod:`xenio`:

* **which section.** ``xenium_paths.csv`` carries several bundles for 6533 and 112 at
  ``roi=panc`` (a 2024 nPOD run and a 2026 Pheno-Xeno run), so ``serial`` selects one and an
  unnamed ambiguity raises. This module used to take ``.iloc[0]`` — the 2024 run, silently.
* **an absolute recorded path.** The Pheno-Xeno bundles are recorded absolute because they do
  not live in the ``HuPrime5kPanTissue/<batch>/<run>/`` tree the ``strip``/``root`` rule
  rewrites. Stripping leading components off an absolute path produced
  ``/orange/brusko/xenium/blue/maigan/...``, which exists nowhere, so those two donors
  resolved at no site at all.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))                 # the lane root, for `xenio`

import xenio                                         # noqa: E402

SITES_JSON = Path(os.environ.get("XENIUM_SITES", HERE / "sites.json"))


class SiteResolutionError(RuntimeError):
    """A donor's data cannot be located — never downgraded to a skip."""


@dataclass(frozen=True)
class Resolved:
    donor_id: str
    roi: str
    site: str
    bundle: Path
    cells_parquet: Path
    h5ad: Path | None      # None means bundle-only, which is a fact about 4 donors
    expression_source: str  # "h5ad" | "bundle"
    serial: str = ""              # WHICH section; the bundle dir's XETG serial
    n_candidates: int = 1         # bundles recorded for this donor+roi, before selection
    h5ad_refused_reason: str = ""  # non-empty when an existing h5ad was deliberately not used


def _load_sites() -> dict:
    sites = json.loads(SITES_JSON.read_text())
    return {k: v for k, v in sites.items() if not k.startswith("_")}


def _candidates(rel: str, site_cfg: dict) -> list[Path]:
    # An absolute recorded path is already a location: it names no `HuPrime5kPanTissue/<batch>`
    # prefix to strip and no root to join, and doing either to it invents a path. It is
    # returned as itself and checked for existence like any other candidate.
    if Path(rel).is_absolute():
        return [Path(rel.rstrip("/"))]
    parts = Path(rel.strip("/")).parts
    out = []
    for entry in site_cfg["bundle_roots"]:
        strip = int(entry.get("strip", 0))
        if strip >= len(parts):
            continue
        out.append(Path(entry["root"]).joinpath(*parts[strip:]))
    return out


def _h5ad_for(donor_id: str, roi: str, h5ad_dir: Path) -> Path | None:
    tag = "Panc" if roi == "panc" else "pLN"
    if not h5ad_dir.is_dir():
        return None
    for p in sorted(h5ad_dir.glob(f"*_{tag}Xenium.h5ad")):
        if p.name.split("_")[0].lstrip("HDL").rstrip("abc") == donor_id:
            return p
    return None


def resolve(donor_id: str, xenium_paths_csv: Path, *, roi: str = "panc",
            site: str | None = None, serial: str | None = None) -> Resolved:
    """Locate one SECTION. Raises :class:`SiteResolutionError` rather than returning None.

    ``serial`` names which bundle when a donor+roi carries more than one; omitting it there
    raises rather than picking, because the bundles are different physical sections.
    """
    sites = _load_sites()
    site = site or os.environ.get("XENIUM_SITE") or None
    wanted = [site] if site else list(sites)
    if site and site not in sites:
        raise SiteResolutionError(f"XENIUM_SITE={site!r} is not in {SITES_JSON}")

    try:
        row = xenio.select_bundle_row(donor_id, xenium_paths_csv, roi=roi, serial=serial)
    except xenio.XeniumIOError as e:
        raise SiteResolutionError(str(e)) from e
    rel = str(row.source_path)
    n_candidates = int(row["n_candidates"])
    refused = xenio.h5ad_refusal(donor_id, roi, n_candidates)

    tried: list[str] = []
    for name in wanted:
        cfg = sites[name]
        for cand in _candidates(rel, cfg):
            cells = cand / "cells.parquet"
            tried.append(str(cells))
            if not cells.is_file():
                continue
            # The h5ad filename carries no serial, so for a multi-section donor it cannot be
            # paired with this bundle -- see xenio.h5ad_refusal. 6533 and 112 both HAVE one.
            h5 = None if refused else _h5ad_for(donor_id, roi, Path(cfg["h5ad_dir"]))
            matrix = cand / "cell_feature_matrix.h5"
            if h5 is None and not matrix.is_file():
                raise SiteResolutionError(
                    f"donor {donor_id} roi={roi}: found {cells} but neither a usable h5ad "
                    f"under {cfg['h5ad_dir']} nor {matrix}. Expression has no source; this "
                    f"donor cannot be run and must not be silently omitted."
                    + (f" (the h5ad there was refused: {refused})" if refused else ""))
            return Resolved(donor_id=donor_id, roi=roi, site=name, bundle=cand,
                            cells_parquet=cells, h5ad=h5,
                            expression_source="h5ad" if h5 is not None else "bundle",
                            serial=str(row.serial), n_candidates=n_candidates,
                            h5ad_refused_reason=refused if h5 is None else "")

    raise SiteResolutionError(
        f"donor {donor_id} roi={roi} serial={row.serial} resolved at no site. Recorded path "
        f"{rel!r}. Tried:\n  " + "\n  ".join(tried))
