"""Step 0: carry QuPath's own signed islet distance across into this lane.

The run does not have it.  ``expression`` carries ``cell_region`` (``core`` / ``peri`` /
``tissue`` / ``other``) and ``islet_num``, both derived from the QuPath ``Parent`` string at
ingest, and nothing else about where a cell sits relative to an islet.  ``cell_region ==
'peri'`` is the **20 um expansion annotation**, so it answers "inside the ring or not" and
cannot answer "how far".

KINTSUGI's peri-islet analysis is built on a continuous signed distance -- its zone is
``0 <= d <= 50`` um, which the 20 um ring cannot express -- so that measurement has to come
from somewhere.  It already exists: QuPath computed it and wrote it into the cohort exports as
``Signed distance to annotation Islet um``, alongside ``Object ID``.  Ingest simply did not
select it into ``cells/``.

So this module does not *compute* a distance.  It carries across a measurement that was made
upstream, keyed by the same ``object_id`` as everything else, which is what ``CLAUDE.md``'s
"preserve the QuPath UUID" asks for and what makes the join to ``expression`` exact rather
than approximate.  ``data/Islets/*.geojson`` holds the islet polygons for all 20 donors and
could be used to recompute it, but a recomputation would be a *second* measurement with its
own boundary and rounding conventions; QuPath's is the one the annotations were drawn in and
the one KINTSUGI's published 6551 analysis used.

**The sign convention is checked, not assumed.**  ``validate`` compares the extracted value
against ``cell_region``, which was derived independently at ingest from a different column:
a ``core`` cell sits inside the islet annotation and must read <= 0, a ``peri`` cell sits in
the 20 um ring and must read within (0, 20].  The violation rate is reported rather than
suppressed, because a disagreement here means either the join or the convention is wrong and
both are silent failures otherwise.

Two things worth carrying:

* **The distances do not depend on the parenting.**  Donor 6539's QuPath parent assignment was
  repaired and re-exported on 2026-08-27 (``Cellmeasurements_6539parentFixed.csv``), taking it
  from 3 rings carrying cells to 569.  Measured over all 503,962 shared cells, **both**
  signed distances are bit-identical before and after: islet ``max |diff| = 0.0``, tissue
  non-null in both and equal.  The geometry was always right; only ``Parent`` was wrong.  So
  what this module extracted from the original export stands, and what the repair changes is
  ``cell_region`` -- which lives upstream in ``cells/``, not here.
* A *missing* measurement is null, never 0.  A distance of 0.0 means "on the boundary", and
  conflating that with "not measured" is exactly what ``CLAUDE.md`` forbids.
"""
from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))  # reuse the cohort/panel authority as-is
sys.path.insert(0, str(HERE.parents[1]))           # core/phenocycler, for `phenocycler`

import panel                                       # noqa: E402
from phenocycler.cells_parquet import (             # noqa: E402
    is_parquet_measurement,
    signed_distance_column,
)

#: The exact QuPath column names, in both spellings: QuPath 0.7 (the v3 re-export) writes
#: "annotation with Islet", the v2 batch exports "annotation Islet".  They carry a literal
#: MICRO SIGN; matching them by prefix or by ``LIKE`` would also match a differently-scoped
#: distance if one is ever added, so each is spelled out and the one present is asserted.
ISLET_COLS = ("Signed distance to annotation with Islet µm",
              "Signed distance to annotation Islet µm")
TISSUE_COLS = ("Signed distance to annotation with Tissue µm",
               "Signed distance to annotation Tissue µm")

#: Donors whose tissue distance is not a measurement.  **Empty, and 6539 is deliberately not
#: in it.**
#:
#: An earlier version of this file listed 6539 on the strength of
#: ``panel_availability.json``'s ``has_tissue_annotation: false``.  That field is inferred
#: from whether ``Annotation (Tissue)`` ever appears in the ``Parent`` column, which is a
#: statement about hierarchy resolution, not about whether the annotation exists.  Measured:
#: 6539's export carries a non-null ``Signed distance to annotation Tissue um`` for all
#: 503,962 of its cells, and those values are **identical** before and after the 2026-08-27
#: parent fix.  QuPath could not have computed them without the annotation.  So the tissue
#: distance was always available; what was missing was the parenting.
NO_TISSUE_ANNOTATION: dict[str, str] = {}

#: This lane writes beside the cosine lane, never into it.
LANE_OUT = panel.REPO / "data" / "development" / "kintsugi"
OUT = LANE_OUT / "distance"


def _sql_string(value: str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def csv_groups(cohort) -> dict[pathlib.Path, dict[str, str]]:
    """Group the manifest's images by the export they were measured in.

    The cohort spans two exports with different column counts (1,228 and 992), so they are
    read separately rather than unioned by name -- a union would make every read pay for the
    other file's columns.
    """
    groups: dict[pathlib.Path, dict[str, str]] = {}
    for image in cohort.images:
        # A v3 image carries its signed distances in the context Parquet, not the cells one.
        source = image.measurement_context or image.measurement_csv
        path = pathlib.Path(source.path)
        groups.setdefault(path, {})[image.image_id] = str(image.donor_id)
    return groups


def _donor_case(image_to_donor: dict[str, str]) -> str:
    branches = " ".join(
        f"WHEN {_sql_string(image)} THEN {_sql_string(donor)}"
        for image, donor in sorted(image_to_donor.items())
    )
    return (
        f'CASE "Image" {branches} ELSE '
        "error('unmapped QuPath Image: ' || coalesce(\"Image\", '<NULL>')) END"
    )


def read_source(path: pathlib.Path, *, sample: bool = False) -> str:
    """The DuckDB table function for one export: Parquet (v3) or CSV/TSV (v2)."""
    escaped = str(path).replace("'", "''")
    if is_parquet_measurement(path):
        return f"read_parquet('{escaped}')"
    extra = ", sample_size=1" if sample else ", ignore_errors=true"
    return f"read_csv('{escaped}', header=true, all_varchar=true{extra})"


def assert_columns_present(con, csv_path: pathlib.Path) -> tuple[str, str]:
    """Fail before the scan if the export does not carry the distance columns.

    Returns the (islet, tissue) column names actually present, in whichever spelling.
    ``Image`` is required only of a CSV export: a v3 context Parquet is per image and keyed
    by ``Object ID`` alone, and its image comes from the manifest.
    """
    names = {
        row[0]
        for row in con.execute(
            f"SELECT column_name FROM (DESCRIBE SELECT * FROM "
            f"{read_source(csv_path, sample=True)})"
        ).fetchall()
    }
    islet = signed_distance_column(names, "Islet")
    tissue = signed_distance_column(names, "Tissue")
    keys = ("Object ID",) if is_parquet_measurement(csv_path) else ("Object ID", "Image")
    missing = [
        *([] if islet else [" or ".join(repr(c) for c in ISLET_COLS)]),
        *([] if tissue else [" or ".join(repr(c) for c in TISSUE_COLS)]),
        *(c for c in keys if c not in names),
    ]
    if missing:
        raise RuntimeError(
            f"{csv_path} does not carry {missing}. This lane cannot substitute a "
            "recomputed distance for the one the annotations were drawn in."
        )
    return islet, tissue


def extract_sql(csv_path: pathlib.Path, image_to_donor: dict[str, str],
                out_dir: pathlib.Path, *, islet_col: str = ISLET_COLS[1],
                tissue_col: str = TISSUE_COLS[1]) -> str:
    """One column-selective pass over one export, partitioned by donor.

    ``all_varchar`` then ``TRY_CAST`` rather than letting the sniffer type 1,228 columns:
    the projection is pushed down, so only these four are ever converted, and a value that
    does not parse becomes NULL -- unavailable -- instead of failing the whole scan.

    ``islet_col`` / ``tissue_col`` are the spellings :func:`assert_columns_present` found.
    A v3 context Parquet is one image with no ``Image`` column, so its single manifest image
    names the donor directly.
    """
    out = str(out_dir).replace("'", "''")
    if is_parquet_measurement(csv_path):
        if len(image_to_donor) != 1:
            raise ValueError(
                f"{csv_path} is a per-image context Parquet but {len(image_to_donor)} "
                "manifest images point at it"
            )
        donor = _sql_string(next(iter(image_to_donor.values())))
        select_donor, where = f"{donor} AS donor_id", ""
    else:
        allowed = ", ".join(_sql_string(i) for i in sorted(image_to_donor))
        select_donor = f"{_donor_case(image_to_donor)} AS donor_id"
        where = f'WHERE "Image" IN ({allowed})'
    return f"""
COPY (
    SELECT {select_donor},
           "Object ID" AS object_id,
           TRY_CAST("{islet_col}" AS DOUBLE) AS signed_dist_islet_um,
           TRY_CAST("{tissue_col}" AS DOUBLE) AS signed_dist_tissue_um
    FROM {read_source(csv_path)}
    {where}
) TO '{out}' (FORMAT PARQUET, PARTITION_BY (donor_id),
              OVERWRITE_OR_IGNORE, FILENAME_PATTERN 'data_{{i}}', COMPRESSION zstd)
"""


#: Centroids sit anywhere inside a cell, and a cell straddling a boundary can have its
#: centroid on either side of it, so exact agreement is not the right bar.  What would matter
#: is a *systematic* disagreement -- a sign flip, a unit error, a mis-joined donor, or a ring
#: that is not the width it is named -- and those move these rates to tens of percent.
MAX_VIOLATION_RATE = 0.05

#: Slack on the measured ring width, for centroids of cells straddling the outer boundary.
RING_TOLERANCE_UM = 5.0


def assert_pinned_run_present() -> None:
    """Fail with the situation, not a bare FileNotFoundError, if the pinned run is gone.

    ``panel.RUN_ID`` names the run both this lane and the cosine lane are pinned to, and
    ``data/runs/`` is not something either lane writes.  If it has been removed, that is a
    fact about the workspace that should be stated once and clearly, rather than surfacing as
    a missing-parquet traceback from whichever donor happened to be read first.
    """
    if panel.EXPRESSION.exists():
        return
    runs = panel.PKG / "data" / "runs"
    present = sorted(d.name for d in runs.iterdir() if d.is_dir()) if runs.exists() else []
    latest = (runs / "LATEST").read_text().strip() if (runs / "LATEST").exists() else "<absent>"
    raise FileNotFoundError(
        f"the pinned run {panel.RUN_ID} is not on disk.\n"
        f"  expected: {panel.EXPRESSION}\n"
        f"  runs present: {present or '<none>'}\n"
        f"  data/runs/LATEST says: {latest}\n"
        "This lane and scripts/celltype are both pinned to that run, and its content hash is "
        "asserted by scripts/celltype/tests. Re-pinning is a deliberate change of provenance "
        "-- decide it explicitly rather than letting a different run be picked up silently."
    )

def validate(donor_id: str, run_expression: pathlib.Path | None = None) -> dict:
    """Check the extraction against ``cell_region``, which was derived independently.

    ``cell_region`` comes from the QuPath ``Parent`` string at ingest; the distance comes
    from a different column of the same export.  They agree only if the join is right and the
    sign convention is what this lane assumes, so this is a real check rather than a
    restatement.

    Two categories are counted separately because they mean different things:

    * **core_positive** -- a cell whose ``Parent`` is ``Islet_N`` reading a positive distance.
      A genuine sign or join failure.  Measured at 0.4-2.7% per donor, which is the rate at
      which a cell straddling the boundary has its centroid outside it.
    * **peri_nonpositive** -- a cell in ``Islet_N_exp20um`` reading <= 0.  **Expected, and not
      a defect.**  The QuPath column is the distance to the *nearest* islet annotation, while
      ``Parent`` names one specific islet, so where islets sit close together a cell in one
      islet's ring can be inside another.  Uniform at ~1-2% across every donor, which is what
      a geometric effect looks like; it is reported, not counted against the bar.
    * **peri_beyond_ring** -- a cell in the ring reading further out than that donor's own
      measured ring width.  Checked against :mod:`ring_width`, never against the nominal 20 um
      the annotation is *named* for: donor 6482's expansion is 25.1 um.
    """
    import pyarrow.parquet as pq

    import ring_width

    assert_pinned_run_present()

    expression = run_expression or panel.EXPRESSION
    exp = pq.read_table(expression / f"donor_id={donor_id}" / "data.parquet",
                        columns=["object_id", "cell_region"]).to_pandas()
    dist = pq.read_table(OUT / f"donor_id={donor_id}").to_pandas()

    merged = exp.merge(dist[["object_id", "signed_dist_islet_um"]], on="object_id", how="left")
    unmatched = int(merged["signed_dist_islet_um"].isna().sum())
    d = merged["signed_dist_islet_um"]

    ring = ring_width.measure(str(donor_id))
    outer = ring["ring_width_um_median"] + RING_TOLERANCE_UM

    core = merged["cell_region"] == "core"
    peri = merged["cell_region"] == "peri"
    core_pos = int((core & d.notna() & (d > 0)).sum())
    peri_nonpos = int((peri & d.notna() & (d <= 0)).sum())
    peri_beyond = int((peri & d.notna() & (d > outer)).sum())

    n_core, n_peri = int(core.sum()), int(peri.sum())
    return {
        "donor_id": str(donor_id),
        "n_expression": len(merged),
        "n_distance": len(dist),
        "unmatched_object_id": unmatched,
        "ring_width_um_median": ring["ring_width_um_median"],
        "n_rings_measured": ring["n_rings_measured"],
        "n_rings_in_geojson": ring["n_rings_in_geojson"],
        "ring_matches_nominal": ring["matches_nominal"],
        "n_core": n_core,
        "core_positive": core_pos,
        "core_positive_rate": (core_pos / n_core) if n_core else float("nan"),
        "n_peri": n_peri,
        "peri_nonpositive": peri_nonpos,
        "peri_nonpositive_rate": (peri_nonpos / n_peri) if n_peri else float("nan"),
        "peri_beyond_ring": peri_beyond,
        "peri_beyond_ring_rate": (peri_beyond / n_peri) if n_peri else float("nan"),
        "tissue_available": str(donor_id) not in NO_TISSUE_ANNOTATION,
        "tissue_unavailable_reason": NO_TISSUE_ANNOTATION.get(str(donor_id), ""),
    }
