"""Step 4: tissue domains with UTAG, beside the per-cell calls.

UTAG (Kim et al., Nature Methods 2022; `ElementoLab/utag`) answers a different question from
Step 3 and that is why it is here.  Step 3 asks *what is this cell*; UTAG asks *what
neighbourhood is it in*, by message-passing marker expression over a spatial graph and
clustering the result.  What comes out are microanatomical domains -- islet, duct, vessel,
stroma, immune aggregate -- discovered from the data rather than drawn.

For this study that is the point.  The question is the immune environment **around islets**,
and islet proximity currently comes from one source: QuPath's `Islet_N` polygons and the
20 um expansion ring.  UTAG gives a second, independent account of the same tissue, built
from expression and position instead of from an annotator's outline -- so the two can be
compared rather than one assumed.  `cell_region` is deliberately not an input.

**Inputs are the aligned values, not the raw ones.**  Message passing averages a marker over
a cell's neighbours, so a donor-level scale difference would propagate into every domain.
Step 2's shift is what makes that safe.

**Distance.** ``max_dist`` is in the units of ``obsm['spatial']``, which here is micrometres
(the expression stage stores centroids in um).  UTAG's default of 20 was chosen for IMC at
1 um resolution; on this cohort the median cell area is 53 um^2, so a cell radius is ~4 um and
20 um reaches roughly two to three cell diameters.  Recorded rather than inherited silently.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402

#: The UTAG checkout this lane was built against.  Pinned the same way as UniFORM: the repo
#: is importable but not published, so the assertion is what stops a drifted checkout from
#: silently changing the domains.
UTAG_COMMIT = "056cb87ee90aa850450b806b075033bc8c50703a"
UTAG_PATH = panel.REPO / "external" / "utag"

ALIGNED = panel.OUT / "uniform"
OUT = panel.OUT / "domains"
MAX_DIST_UM = 20.0


def _import_utag():
    if not UTAG_PATH.exists():
        raise RuntimeError(
            f"UTAG checkout missing at {UTAG_PATH}.\n"
            f"  git clone https://github.com/smith6jt-cop/utag.git external/utag\n"
            f"  git -C external/utag checkout {UTAG_COMMIT}")
    head = subprocess.run(["git", "-C", str(UTAG_PATH), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    if head != UTAG_COMMIT:
        raise RuntimeError(f"UTAG checkout is at {head}, pinned to {UTAG_COMMIT}")
    sys.path.insert(0, str(UTAG_PATH))
    from utag import utag                                   # noqa: F401
    return utag


def donor_adata(donor: str, markers: list, max_cells: int | None = None, seed: int = 0):
    """One donor's aligned values as an AnnData with micrometre centroids in obsm['spatial'].

    ``cell_region`` and ``islet_num`` ride along in ``obs`` but are NOT variables: they are
    the independent account the domains get compared against, and feeding them in would make
    that comparison circular.
    """
    import anndata

    have = [m for m in markers if m in set(panel.donor_markers(donor))]
    cols = ["object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
            "cell_region", "islet_num"] + have
    df = pq.read_table(ALIGNED / f"donor_id={donor}" / "data.parquet",
                       columns=cols).to_pandas()
    df = df[df["qc_analysis_eligible"]]
    df = df[np.isfinite(df[have + ["X_centroid", "Y_centroid"]].to_numpy("float32")).all(1)]
    if max_cells and len(df) > max_cells:
        df = df.iloc[np.random.default_rng(seed).choice(len(df), max_cells, replace=False)]
    df = df.reset_index(drop=True)

    a = anndata.AnnData(np.arcsinh(df[have].to_numpy("float32")))    # variance-stabilised
    a.var_names = have
    a.obs = df[["object_id", "cell_region", "islet_num"]].copy()
    a.obs["donor_id"] = donor
    a.obsm["spatial"] = df[["X_centroid", "Y_centroid"]].to_numpy("float32")
    return a


def domain_islet_concordance(obs: pd.DataFrame, label: str) -> pd.DataFrame:
    """Does any UTAG domain correspond to the islets QuPath drew?

    The check Step 3 uses, pointed the other way: there the call was known and the geometry
    tested it; here the domain is unknown and the geometry names it.  Neither uses the other's
    inputs -- `cell_region` is polygon containment, the domains are expression and position.
    """
    base = (obs["cell_region"] == "core").mean()
    rows = []
    for d, g in obs.groupby(label, observed=True):
        frac = (g["cell_region"] == "core").mean()
        rows.append({"domain": str(d), "n": len(g), "core_frac": frac,
                     "peri_frac": (g["cell_region"] == "peri").mean(),
                     "base_rate": base,
                     "enrichment": frac / base if base > 0 else np.nan})
    return pd.DataFrame(rows).sort_values("enrichment", ascending=False)
