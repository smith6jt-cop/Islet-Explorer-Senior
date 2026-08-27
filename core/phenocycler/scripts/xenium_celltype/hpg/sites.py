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
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
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


def _load_sites() -> dict:
    sites = json.loads(SITES_JSON.read_text())
    return {k: v for k, v in sites.items() if not k.startswith("_")}


def _candidates(rel: str, site_cfg: dict) -> list[Path]:
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
            site: str | None = None) -> Resolved:
    """Locate one donor. Raises :class:`SiteResolutionError` rather than returning None."""
    sites = _load_sites()
    site = site or os.environ.get("XENIUM_SITE") or None
    wanted = [site] if site else list(sites)
    if site and site not in sites:
        raise SiteResolutionError(f"XENIUM_SITE={site!r} is not in {SITES_JSON}")

    xen = pd.read_csv(xenium_paths_csv, dtype=str)
    hit = xen[(xen.donor_id == donor_id) & (xen.roi == roi)]
    hit = hit[hit.source_path.notna() & (hit.source_path != "NA")]
    if hit.empty:
        raise SiteResolutionError(
            f"donor {donor_id} roi={roi} has no bundle path in {xenium_paths_csv}")
    rel = str(hit.iloc[0].source_path)

    tried: list[str] = []
    for name in wanted:
        cfg = sites[name]
        for cand in _candidates(rel, cfg):
            cells = cand / "cells.parquet"
            tried.append(str(cells))
            if not cells.is_file():
                continue
            h5 = _h5ad_for(donor_id, roi, Path(cfg["h5ad_dir"]))
            matrix = cand / "cell_feature_matrix.h5"
            if h5 is None and not matrix.is_file():
                raise SiteResolutionError(
                    f"donor {donor_id} roi={roi}: found {cells} but neither an h5ad under "
                    f"{cfg['h5ad_dir']} nor {matrix}. Expression has no source; this donor "
                    f"cannot be run and must not be silently omitted.")
            return Resolved(donor_id=donor_id, roi=roi, site=name, bundle=cand,
                            cells_parquet=cells, h5ad=h5,
                            expression_source="h5ad" if h5 is not None else "bundle")

    raise SiteResolutionError(
        f"donor {donor_id} roi={roi} resolved at no site. Recorded relative path {rel!r}. "
        f"Tried:\n  " + "\n  ".join(tried))
