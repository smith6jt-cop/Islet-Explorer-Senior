#!/usr/bin/env python
"""Turn a Xenium bundle into the h5ad the established typer expects.

`analysis/xenium/scripts/hierarchical_phenotype.py` is the reference Xenium typer -- two-tier
identity-core anchoring with scANVI propagation, canonical under that lane's do-not-regress
decisions #16/#17 -- and `phenocycler/integration/import_xenium.py` is already written against
its output. It has never run on this cohort, because what it consumes does not exist here: the
46 h5ads under `/data/xenium/h5ad_out/` carry raw counts and donor metadata and **nothing
else** -- no `obsm['spatial']`, no layers, no cell types. Four cohort donors have no h5ad at
all.

This builds that input. It is deliberately NOT a port of `notebooks/01_preprocessing_v2.ipynb`:
the typer needs exactly three things -- `layers['counts']`, `layers['lognorm']` and
`var['panel_for_embedding']` (the last stamped afterwards by `finalize_panel.py`) -- and it
fails loud on each. The notebook's remaining outputs (a flat scVI latent, a leiden, a UMAP)
are consumed by nothing downstream of the hierarchical rewire.

**Why `prefer="bundle"` for every donor, not just the four without an h5ad.** The dual source
is a fact about which files were transferred, not a methodological fork, and `xenio`'s own
cross-validation on donor 6539 settles that they agree: identical gene order, scores to 6e-7,
100% argmax agreement over 20,000 cells. What differs is the *population* -- the h5ad is an
upstream Seurat QC subset (317,610 of 6539's 319,005 cells). Reading the bundle everywhere
means all 20 donors start from the same cells and take the same QC recorded here, instead of
16 inheriting an undocumented upstream filter and 4 taking none. That asymmetry would sit
underneath every cross-donor comparison the study makes.

QC follows the Xenium lane's decision #10: absolute LOWER bounds (`min_genes`, `min_counts`)
because those are sanity floors for noise and segmentation artifacts, and a DENSITY-based
upper bound (`counts / cell_area` above a quantile) because in FFPE Xenium `cell_area` varies
5-10x and an absolute counts cap predominantly removes large *healthy* cells -- acinar,
ductal, large epithelial -- rather than doublets.

    cd core/phenocycler
    ../../.venv/bin/python scripts/xenium_celltype/build_phenotyped_base.py --donors 6539
    ../../.venv/bin/python scripts/xenium_celltype/build_phenotyped_base.py --workers 8

Writes `<out>/<serial>/<roi>/<serial>_phenotyped.h5ad`. The nesting is not cosmetic:
`integration/manifest.py` only falls back to the flat `<serial>/<serial>_phenotyped.h5ad` when
a serial carries exactly one ROI, and 25 of this cohort's 26 serials carry more than one
(every pancreas slide also carries pLN regions). The flat name would never be found.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))             # core/phenocycler, for `phenocycler`

import xenio                                          # noqa: E402

DEFAULT_OUT = Path("/data/xenium/processed")
DEFAULT_PATHS_CSV = HERE.parents[1] / "data" / "integration" / "xenium_paths.csv"

#: Decision #10. Lower bounds absolute (noise floors, not size-correlated); upper bound on
#: DENSITY, because Spearman(total_counts, cell_area) is ~0.78 in FFPE Xenium and an absolute
#: cap would drop large healthy cells rather than doublets.
MIN_GENES = 20
MIN_COUNTS = 50
DENSITY_QUANTILE = 0.98

#: The bundle directory name carries the slide serial, which is the sample id every downstream
#: consumer keys on -- `manifest.py`'s h5ad discovery, and the existing integration fixture.
#: The pattern now lives in `xenio` and is re-exported here: this module used to own a second
#: copy, and two regexes that are supposed to agree about what a serial is eventually do not.
_SERIAL = xenio.SERIAL_RE


class BuildError(RuntimeError):
    """A donor cannot be built -- never silently skipped."""


def serial_of(source_path: str) -> str:
    """The slide serial in a bundle path, raised as a `BuildError` for this lane's callers."""
    try:
        return xenio.serial_of(source_path)
    except xenio.XeniumIOError as e:
        raise BuildError(str(e)) from e


def bundle_extras(bundle: Path, cell_ids: np.ndarray) -> pd.DataFrame:
    """`cell_area` and friends for the cells we kept, in their order.

    `xenio.load_donor` reads only the three columns it needs. `cell_area` is read here rather
    than added there because it is a QC input, not a coordinate, and the density filter is
    this script's decision -- not something every consumer of `load_donor` should inherit.
    """
    cols = ["cell_id", "cell_area", "nucleus_area", "transcript_counts", "segmentation_method",
            # The vendor's own totals. `total_counts` INCLUDES control and genomic counts and
            # is the only record of them once the matrix is loaded gene-expression-only;
            # `sc.pp.calculate_qc_metrics` overwrites `obs['total_counts']` with a gene-only
            # row sum wherever it is called, so it is preserved here under its own name.
            "total_counts", "genomic_control_counts"]
    df = pd.read_parquet(bundle / "cells.parquet", columns=cols)
    df["cell_id"] = df.cell_id.astype(str)
    df = df.set_index("cell_id").reindex(pd.Index(cell_ids, name="cell_id"))
    if df.cell_area.isna().any():
        raise BuildError(f"{bundle}: cell_area missing for "
                         f"{int(df.cell_area.isna().sum())} of {len(df)} cells")
    return df


def qc_mask(counts: sp.csr_matrix, cell_area: np.ndarray, *,
            nucleus_area: np.ndarray | None = None,
            gdna_lambda: np.ndarray | None = None,
            n_panel_genes: int | None = None) -> tuple[np.ndarray, dict]:
    """Decision #10's filter, returning the mask AND what it removed.

    ``MIN_COUNTS`` is applied to **gDNA-corrected** counts when ``gdna_lambda`` is supplied.
    On these runs the genomic-DNA background is ~126 transcripts in a median 433-transcript
    cell, so a floor of 50 on raw counts admits cells that are majority background: measured,
    14.8% (0080870) and 18.6% (0081366) of the cells a raw floor of 50 keeps have more
    expected background than real signal. Correcting the quantity rather than raising the
    number costs 7.0% / 9.4% of cells; reaching the same place by raising the RAW floor to
    300 costs 31% / 34% and still admits them.

    The correction subtracts an EXPECTATION, not identified reads -- nothing says which counts
    are genomic. That is sound for a threshold and for aggregates, and must not be mistaken
    for a cleaned expression vector.
    """
    n_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    cell_area = np.asarray(cell_area, dtype=float)
    area_ok = np.isfinite(cell_area) & (cell_area > 0)
    nucleus_ok = (np.ones(len(cell_area), bool) if nucleus_area is None else
                  np.isfinite(nucleus_area) & (np.asarray(nucleus_area) > 0))
    density = np.divide(n_counts, cell_area, out=np.full(len(n_counts), np.nan),
                        where=area_ok)
    eligible = area_ok & nucleus_ok
    if not eligible.any():
        raise BuildError("no cells have valid area and a measured nucleus")
    cap = float(np.quantile(density[eligible], DENSITY_QUANTILE))

    info = {"min_genes": MIN_GENES, "min_counts": MIN_COUNTS,
            "density_quantile": DENSITY_QUANTILE, "density_cap_counts_per_um2": cap,
            "dropped_min_genes": int((n_genes < MIN_GENES).sum()),
            "dropped_density": int((density > cap).sum()),
            "dropped_invalid_area": int((~area_ok).sum()),
            "dropped_no_nucleus": int((~nucleus_ok).sum()),
            "nucleus_required": nucleus_area is not None}

    if gdna_lambda is None or n_panel_genes is None:
        counts_for_floor = n_counts
        info["min_counts_applied_to"] = "raw"
    else:
        expected_bg = np.asarray(gdna_lambda, dtype=float) * float(n_panel_genes)
        counts_for_floor = np.maximum(n_counts - expected_bg, 0.0)
        info.update({
            "min_counts_applied_to": "gdna_corrected",
            "expected_background_median": float(np.median(expected_bg)),
            "frac_cells_majority_background": float((expected_bg > 0.5 * n_counts).mean()),
            "frac_cells_background_exceeds_counts": float((expected_bg > n_counts).mean()),
        })
    info["dropped_min_counts"] = int((counts_for_floor < MIN_COUNTS).sum())

    keep = eligible & (n_genes >= MIN_GENES) & (counts_for_floor >= MIN_COUNTS) & (density <= cap)
    return keep, info


def _object_index(values) -> pd.Index:
    """An unnamed index of plain Python strings -- see the note in `build_one`."""
    return pd.Index(np.asarray(values, dtype=object), dtype=object)


def _object_categorical(values) -> pd.Categorical:
    """A categorical whose CATEGORIES are `object`, not pandas 3's inferred StringDtype."""
    arr = np.asarray(values, dtype=object)
    return pd.Categorical(arr, categories=pd.Index(np.unique(arr), dtype=object))


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_one(donor: str, *, roi: str, out_root: Path, paths_csv: Path,
              overwrite: bool = False, serial: str | None = None,
              bundle: Path | None = None) -> dict:
    """Build one SECTION. `serial` is required wherever a donor+roi carries more than one.

    The bundle is resolved ONCE, here, and the resolved serial is then handed to
    `xenio.load_donor` so both agree by construction. They used to resolve independently --
    two `.iloc[0]` selections over the same CSV -- which meant the output path's serial and
    the cells written under it were only ever equal by coincidence. For 6533 and 112 that
    coincidence also picked the wrong section: `.iloc[0]` is the 2024 nPOD run (0041393 /
    0041396), not the 2026 Pheno-Xeno bundle (0080870 / 0081366).
    """
    import anndata as ad

    ref = None if bundle is not None else xenio.resolve_bundle(donor, paths_csv, roi=roi, serial=serial)
    if bundle is not None:
        bundle = Path(bundle).resolve()
        found_serial = serial_of(str(bundle))
        if serial is not None and serial != found_serial:
            raise BuildError(f"serial {serial} disagrees with bundle {found_serial}")
        source_path, serial = str(bundle), found_serial
    else:
        source_path, serial = ref.source_path, ref.serial
    out = out_root / serial / roi / f"{serial}_phenotyped.h5ad"
    if out.exists() and not overwrite:
        return {"donor_id": donor, "serial": serial, "status": "exists", "path": str(out)}

    # the resolved serial, not the caller's possibly-absent one -- one section, named twice
    if bundle is not None:
        ids, genes, matrix = xenio._read_bundle_matrix(bundle)
        coordinates = pd.read_parquet(bundle / "cells.parquet").set_index("cell_id")
        coordinates.index = coordinates.index.astype(str)
        xy = coordinates.reindex(ids)[["x_centroid", "y_centroid"]].to_numpy()
        if not np.isfinite(xy).all():
            raise BuildError("bundle contains missing/non-finite cell coordinates")
        dd = xenio.DonorData(donor_id=donor, cell_id=ids, genes=genes, counts=matrix,
                            xy=xy, source="bundle", serial=serial, bundle=bundle,
                            n_dropped_no_spatial=0)
    else:
        dd = xenio.load_donor(donor, paths_csv, roi=roi, prefer="bundle", serial=serial)
        bundle = ref.path
    extras = bundle_extras(bundle, dd.cell_id)

    # Per-cell genomic-DNA background, from the bundle's own genomic control probes.
    # `summarise` re-derives the vendor's headline false-positive number and refuses if it
    # does not reproduce, so the correction below rests on a checked identity.
    import gdna as gdna_mod
    comp = gdna_mod.panel_composition(bundle)
    lam, lam_raw, area_tbl, shrink_a = gdna_mod.background_rates(
        extras["genomic_control_counts"].to_numpy(), extras.cell_area.to_numpy(),
        comp["genomic_control"])
    gdna_check = gdna_mod.verify_against_vendor(
        float(lam.mean()), comp["gene"], gdna_mod.vendor_metrics(bundle))

    keep, qc = qc_mask(dd.counts, extras.cell_area.to_numpy(),
                       nucleus_area=extras.nucleus_area.to_numpy(),
                       gdna_lambda=lam, n_panel_genes=comp["gene"])
    counts = dd.counts[keep].tocsr().astype(np.float32)
    lognorm = xenio.lognorm(dd.counts[keep].tocsr()).astype(np.float32)

    # Every string here is written as `object`, not pandas 3's default `StringDtype`.
    # anndata < 0.13 refuses to write a nullable StringArray unless
    # `settings.allow_write_nullable_strings` is toggled on, and this file has to be readable
    # by the typer's env (anndata 0.12.4) as well as by this one (0.13.2).  Categorical for
    # the low-cardinality columns, which is also how they should be stored.
    # Two traps here, both silent, both caught only because the round-trip was checked.
    #
    # 1. The indices are DELIBERATELY UNNAMED.  A named index is written as an ordinary
    #    column plus an `_index` attribute, and round-tripping that through a second writer
    #    (finalize_panel.py, under a different anndata) reset BOTH indices to positional
    #    strings -- var_names became '0','1','2'...  Every gene lookup then missed: the
    #    identity panel came back 0/5101 and donor_sex 'unknown', and nothing raised.  A
    #    typing run on that file would have been meaningless and looked fine.
    #
    # 2. Every string is forced to `object`, not pandas 3's inferred `StringDtype`.
    #    `pd.Index(np.array([...], dtype=object))` returns an ArrowStringArray under pandas
    #    3.0, and anndata < 0.13 REFUSES to write one unless
    #    `settings.allow_write_nullable_strings` is toggled.  Toggling it would make the file
    #    unreadable by the very env the typer runs in (scvi-env, anndata 0.12.4), so the
    #    dtype is pinned instead of the setting.
    #
    # Both are pinned by `tests/test_build_phenotyped_base.py`.
    n_keep = int(keep.sum())
    obs = pd.DataFrame(index=_object_index(dd.cell_id[keep]))
    for c in ("cell_area", "nucleus_area", "transcript_counts",
              "total_counts", "genomic_control_counts"):
        obs[c] = extras[c].to_numpy()[keep]
    obs["vendor_total_counts"] = extras["total_counts"].to_numpy()[keep]
    # The background every downstream detection decision has to know about. Carried on the
    # cell so no later stage has to re-open the bundle to find it.
    obs["gdna_lambda"] = lam[keep].astype(np.float32)
    obs["gdna_raw"] = lam_raw[keep].astype(np.float32)
    obs["gdna_expected_counts"] = (lam[keep] * comp["gene"]).astype(np.float32)
    obs["segmentation_method"] = _object_categorical(
        extras["segmentation_method"].to_numpy()[keep])
    obs["donor_id"] = _object_categorical([donor] * n_keep)
    obs["roi"] = _object_categorical([roi] * n_keep)
    obs["sample"] = _object_categorical([serial] * n_keep)

    adata = ad.AnnData(
        X=lognorm, obs=obs,
        var=pd.DataFrame(index=_object_index(dd.genes)))
    adata.layers["counts"] = counts
    adata.layers["lognorm"] = lognorm
    adata.obsm["spatial"] = dd.xy[keep].astype(np.float32)

    # Provenance. This is NOT a StageManifest and must not pretend to be one -- nothing here
    # writes under data/runs/. But it must be attributable, because it is the input every
    # Xenium cell type in this cohort descends from.
    adata.uns["prep"] = {
        "built_by": "core/phenocycler/scripts/xenium_celltype/build_phenotyped_base.py",
        "source": dd.source,
        "bundle": source_path,
        "bundle_resolved_path": str(bundle),
        "bundle_n_candidates": ref.n_candidates if ref is not None else 1,
        "bundle_selected_by": "explicit bundle" if ref is None else "resolved serial",
        "h5ad_refused_reason": dd.h5ad_refused_reason,
        "xenium_paths_csv_sha256": _sha256(paths_csv) if ref is not None else "",
        "donor_id": donor, "roi": roi, "serial": serial,
        "n_cells_in": int(dd.counts.shape[0]),
        "n_cells_out": int(keep.sum()),
        "n_dropped_no_spatial": int(dd.n_dropped_no_spatial),
        "gdna": {"n_genomic_probes": int(comp["genomic_control"]),
                 "n_panel_genes": int(comp["gene"]),
                 "shrink_a": float(shrink_a),
                 "lambda_mean": float(lam.mean()),
                 "lambda_raw_mean": float(lam_raw.mean()),
                 "vendor_check": gdna_check,
                 "area_table": area_tbl.to_dict("records")},
        "n_genes": int(len(dd.genes)),
        "qc": qc,
        "lognorm": "median-depth normalize_total then log1p (xenio.lognorm)",
    }

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    return {"donor_id": donor, "serial": serial, "status": "written", "path": str(out),
            "n_cells_in": int(dd.counts.shape[0]), "n_cells_out": int(keep.sum()), **qc}


def parse_serial_spec(spec) -> dict[str, str]:
    """`--serial` -> {donor: serial}, as `xenio` defines it, raised as this lane's BuildError."""
    try:
        return xenio.parse_serial_spec(spec)
    except xenio.XeniumIOError as e:
        raise BuildError(str(e)) from e


def preflight_serials(donors, *, roi: str, paths_csv: Path, serials: dict[str, str]) -> None:
    """Resolve every donor BEFORE building any of them, and report all failures at once.

    A cohort build takes hours. Discovering at donor 17 that 6533 needed a serial -- and
    discovering it as a traceback from inside a worker process -- wastes the other sixteen.
    Every donor is resolved here first, so an unnamed ambiguity is a refusal in the first
    second naming every donor it affects and every serial each could mean.
    """
    problems = []
    for d in donors:
        try:
            xenio.resolve_bundle(d, paths_csv, roi=roi, serial=xenio.serial_for(d, serials))
        except xenio.XeniumIOError as e:
            problems.append(f"  {d}: {e}")
    if problems:
        raise BuildError(
            "these donors do not resolve to exactly one Xenium bundle, so nothing was "
            "built:\n" + "\n".join(problems))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None,
                    help="default: the whole eligible cohort")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--bundle", type=Path, help="explicit bundle for exactly one donor")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--xenium-paths", type=Path, default=DEFAULT_PATHS_CSV)
    ap.add_argument("--workers", type=int, default=1,
                    help="gzip-write bound, so parallelism helps; each worker peaks ~3 GB")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--serial", nargs="*", default=None, metavar="[DONOR=]SERIAL",
                    help="Xenium slide serial naming WHICH section to build. Required for a "
                         "donor+roi carrying several bundles -- 6533 and 112 at roi=panc each "
                         "have a 2024 nPOD run and a 2026 Pheno-Xeno run. Either a bare "
                         "serial applied to every donor built (--serial 0080870), or "
                         "donor=serial pairs for a multi-donor run "
                         "(--serial 6533=0080870 112=0081366).")
    args = ap.parse_args(argv)

    if args.donors:
        donors = list(args.donors)
    else:
        import json as _json
        spec = _json.loads((HERE.parents[1] / "cohort-spec.json").read_text())
        donors = list(spec["expected_donors"])
    from phenocycler.cohort import ensure_eligible_donors
    donors = list(ensure_eligible_donors(donors, context="xenium phenotyped base"))
    serials = parse_serial_spec(args.serial)
    if args.bundle is not None:
        if len(donors) != 1:
            ap.error("--bundle requires exactly one --donors value")
    else:
        preflight_serials(donors, roi=args.roi, paths_csv=args.xenium_paths, serials=serials)

    print(f"building {len(donors)} donors  roi={args.roi}  -> {args.out}", flush=True)
    rows = []
    if args.workers <= 1:
        for d in donors:
            rows.append(build_one(d, roi=args.roi, out_root=args.out,
                                  paths_csv=args.xenium_paths, overwrite=args.overwrite,
                                  serial=xenio.serial_for(d, serials), bundle=args.bundle))
            print(f"  {rows[-1]}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            fut = {ex.submit(build_one, d, roi=args.roi, out_root=args.out,
                             paths_csv=args.xenium_paths, overwrite=args.overwrite,
                             serial=xenio.serial_for(d, serials), bundle=args.bundle): d
                   for d in donors}
            for f in as_completed(fut):
                rows.append(f.result())
                print(f"  {rows[-1]}", flush=True)

    summary = args.out / f"build_summary_{args.roi}.json"
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
