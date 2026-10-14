"""Loading Xenium expression + spatial for one donor, from either of the two sources.

There are two, and they are NOT interchangeable by donor:

* ``h5ad_out/<donor>_PancXenium.h5ad`` — 16 of the 20 run-cohort donors. Raw counts and
  donor metadata; **no** cell types, no ``obsm['spatial']``, no layers.
* the bundle's ``cell_feature_matrix.h5`` — the remaining four (6450, 6505, 6534, 6623,
  i.e. the whole of Xenium runs MM046 and MM047), which have no h5ad at all.

Dropping the bundle path would silently delete two of the eleven Xenium runs, and those two
carry two of the cohort's three Aab+ donors. Cross-validated on donor 6539, which has both:
identical gene order, scores agreeing to 6e-7, 100% argmax agreement over 20,000 cells. So
the choice of source is a fact about which files exist, not a methodological fork.

Spatial always comes from the bundle's ``cells.parquet`` (``x_centroid``/``y_centroid``, in
microns) — neither expression source carries coordinates. The h5ad is a strict subset of the
bundle's cells (317,610 of 319,005 for 6539; the difference is the Xenium QC drop), so the
join is an inner join on ``cell_id`` and the dropped count is reported, never assumed.

**A donor is not a section.** ``xenium_paths.csv`` has one row per (donor, roi, batch), and
two donors carry more than one row for ``roi=panc``: 6533 and 112 were each run in the 2024
nPOD cohort (serials 0041393 / 0041396) and again in the 2026 Pheno-Xeno acquisition
(0080870 / 0081366). Selecting with ``.iloc[0]`` returned the 2024 run for both, silently,
and a Pheno-Xeno analysis would have been performed on the wrong physical section with
nothing in the output able to show it.

So selection is explicit. ``select_bundle_row`` takes an optional ``serial`` and RAISES
:class:`BundleAmbiguityError` when a donor+roi carries several bundles and the caller named
none. There is no fail-safe default here and there must not be one: unlike
``integration/modes.py``, where falling back to ``sequential`` merely costs resolution, every
choice among bundles is a different set of cells, so a wrong guess is wrong data rather than
weak data. The caller names the serial or the caller gets an error.

**The h5ad has no serial.** ``<donor>_PancXenium.h5ad`` is keyed on donor and roi only, so
for a multi-section donor there is no way to tell which section it holds — and 6533 and 112
both HAVE one (they are 2024 cohort donors). Pairing that h5ad's expression with a 2026
bundle's coordinates joins two different physical sections on ``cell_id``; the ids are not
disjoint, so the join would not fail, it would just be meaningless. The h5ad is therefore
REFUSED for any donor+roi carrying more than one bundle, and ``h5ad_refused_reason`` records
why the bundle was read instead.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

XENIUM_ROOT = Path("/data/xenium")
H5AD_DIR = XENIUM_ROOT / "h5ad_out"

#: Feature rows in cell_feature_matrix.h5 that are real genes. The file also carries
#: negative-control probes/codewords and deprecated codewords — 5,101 genes out of 9,638
#: or 9,644 features depending on the run's control-set version.
GENE_FEATURE = "Gene Expression"

#: The bundle directory name carries the slide serial, which is the section identifier every
#: downstream consumer keys on. This is the one regex; ``build_phenotyped_base.serial_of``
#: delegates here so the two can never disagree about what a serial is.
SERIAL_RE = re.compile(r"output-XETG\d+__(\d+)__")


class XeniumIOError(RuntimeError):
    """A donor's expression or spatial data cannot be loaded — never silently skipped."""


class BundleAmbiguityError(XeniumIOError):
    """A donor+roi resolves to several bundles and the caller named no serial.

    Separate from :class:`XeniumIOError` because it is the one failure a caller can fix by
    passing an argument, and because a test that asserts "this fails loud" must be able to
    say which loud failure it means.
    """


@dataclass(frozen=True)
class BundleRef:
    """One resolved section: which serial, which batch, and where it is on this filesystem."""

    donor_id: str
    roi: str
    serial: str
    batch: str
    source_path: str          # exactly as recorded in xenium_paths.csv
    path: Path                # the recorded path joined to a root (or itself, if absolute)
    n_candidates: int         # how many bundles this donor+roi carries, before selection


@dataclass(frozen=True)
class DonorData:
    donor_id: str
    cell_id: np.ndarray
    genes: list[str]
    counts: sp.csr_matrix     # cells x genes, RAW
    xy: np.ndarray            # cells x 2, microns
    source: str               # "h5ad" | "bundle"
    n_dropped_no_spatial: int
    #: Which section these cells came from. Defaulted so existing positional construction
    #: (hpg/stage_prep.py builds one directly) keeps working.
    serial: str = ""
    bundle: Path | None = None
    h5ad_refused_reason: str = ""


def serial_of(source_path: str) -> str:
    """The slide serial in a bundle path, e.g. ``0080870``."""
    m = SERIAL_RE.search(str(source_path))
    if not m:
        raise XeniumIOError(f"no XETG serial in bundle path {source_path!r}")
    return m.group(1)


def _normalise_serial(serial: str) -> str:
    """Accept a bare serial or a whole bundle directory name, return the bare serial."""
    s = str(serial).strip().rstrip("/")
    m = SERIAL_RE.search(s)
    return m.group(1) if m else s


def bundle_candidates(donor_id: str, xenium_paths_csv: Path, roi: str = "panc") -> pd.DataFrame:
    """Every bundle row recorded for this donor+roi, with the serial parsed out."""
    xen = pd.read_csv(xenium_paths_csv, dtype=str)
    hit = xen[(xen.donor_id == donor_id) & (xen.roi == roi)]
    hit = hit[hit.source_path.notna() & (hit.source_path != "NA")].copy()
    if not hit.empty:
        hit["serial"] = [serial_of(p) for p in hit.source_path]
    return hit


def select_bundle_row(donor_id: str, xenium_paths_csv: Path, *, roi: str = "panc",
                      serial: str | None = None) -> pd.Series:
    """The ONE bundle row for this donor+roi+serial, or a loud failure.

    This is the whole point of the module's selection rule, so it is written once and the
    three call sites (``load_donor``, ``build_phenotyped_base.build_one``,
    ``hpg/sites.resolve``) all route through it. They used to carry three separate
    ``.iloc[0]`` selections, which is three chances to disagree about which section a donor
    means while every one of them reports success.
    """
    hit = bundle_candidates(donor_id, xenium_paths_csv, roi)
    if hit.empty:
        raise XeniumIOError(
            f"donor {donor_id} roi={roi} has no bundle path recorded in {xenium_paths_csv}")

    if serial:
        want = _normalise_serial(serial)
        chosen = hit[hit.serial == want]
        if chosen.empty:
            have = ", ".join(f"{r.serial} ({r.batch})" for _, r in hit.iterrows())
            raise XeniumIOError(
                f"donor {donor_id} roi={roi} has no bundle with serial {want!r}. "
                f"Recorded: {have}.")
        if len(chosen) > 1:
            raise BundleAmbiguityError(
                f"donor {donor_id} roi={roi} serial={want!r} matches {len(chosen)} rows in "
                f"{xenium_paths_csv}; the serial is supposed to identify one section.")
        row = chosen.iloc[0]
    elif len(hit) > 1:
        listing = "\n".join(
            f"    --serial {r.serial}   batch={r.batch}   {r.source_path}"
            for _, r in hit.iterrows())
        raise BundleAmbiguityError(
            f"donor {donor_id} roi={roi} carries {len(hit)} Xenium bundles and no serial was "
            f"named. Refusing to guess: each bundle is a DIFFERENT physical section, so a "
            f"wrong pick is wrong data, not merely weaker data. Name one:\n{listing}")
    else:
        row = hit.iloc[0]

    out = row.copy()
    out["n_candidates"] = int(len(hit))
    return out


def parse_serial_spec(spec) -> dict[str, str]:
    """A `--serial` command line -> ``{donor: serial}``; a bare serial keys ``"*"``.

    Two forms because two things are being asked for. One donor at a time wants
    ``--serial 0080870``; a cohort command that happens to contain an ambiguous donor wants
    ``--serial 6533=0080870`` and nothing said about the other twenty. Parsed here rather
    than per script so the grammar cannot drift between `run_donor`, `build_phenotyped_base`
    and `hpg/preflight`.
    """
    out: dict[str, str] = {}
    for item in (spec or []):
        donor, sep, serial = str(item).partition("=")
        key, val = (donor.strip(), serial.strip()) if sep else ("*", donor.strip())
        if not val:
            raise XeniumIOError(f"--serial {item!r} is empty; expected SERIAL or DONOR=SERIAL")
        if key in out and out[key] != val:
            raise XeniumIOError(f"--serial names {key!r} twice, as {out[key]!r} and {val!r}")
        out[key] = val
    if "*" in out and len(out) > 1:
        raise XeniumIOError(
            f"--serial mixes a bare serial ({out['*']!r}) with donor=serial pairs. Pick one "
            f"form: a bare serial applies to every donor, which is only meaningful for a "
            f"single-donor command.")
    return out


def serial_for(donor_id: str, serials: dict[str, str]) -> str | None:
    """This donor's serial from a parsed spec: its own entry, else the bare one, else None."""
    return serials.get(donor_id, serials.get("*"))


def join_root(root: Path, source_path: str) -> Path:
    """Join a recorded path to a site root, honouring an absolute recorded path.

    ``xenium_paths.csv`` records the nPOD cohort relative to a
    ``HuPrime5kPanTissue/<batch>/<run>/`` tree, but the two Pheno-Xeno bundles are recorded
    ABSOLUTE because they do not live in that tree at all. ``Path / "/abs"`` already returns
    the absolute path, so this used to work by accident; it is written down because
    ``hpg/sites.py`` strips leading components before joining, which silently mangles an
    absolute path into a root-relative one that does not exist.
    """
    p = str(source_path).rstrip("/")
    return Path(p) if Path(p).is_absolute() else Path(root) / p


def resolve_bundle(donor_id: str, xenium_paths_csv: Path, *, roi: str = "panc",
                   serial: str | None = None, root: Path | None = None) -> BundleRef:
    """Select one section and say where it is under ``root`` (default :data:`XENIUM_ROOT`)."""
    row = select_bundle_row(donor_id, xenium_paths_csv, roi=roi, serial=serial)
    source_path = str(row.source_path).rstrip("/")
    return BundleRef(
        donor_id=donor_id, roi=roi, serial=str(row.serial),
        batch=str(row.get("batch", "")), source_path=source_path,
        path=join_root(XENIUM_ROOT if root is None else root, source_path),
        n_candidates=int(row["n_candidates"]))


def _bundle_dir(donor_id: str, xenium_paths_csv: Path, roi: str = "panc",
                serial: str | None = None) -> Path:
    return resolve_bundle(donor_id, xenium_paths_csv, roi=roi, serial=serial).path


def _read_bundle_matrix(bundle: Path):
    import h5py
    p = bundle / "cell_feature_matrix.h5"
    if not p.exists():
        raise XeniumIOError(f"{p} is absent; transfer the bundle's cell_feature_matrix.h5")
    with h5py.File(p, "r") as f:
        g = f["matrix"]
        # stored CSC over (features, cells) -> the same buffers are CSR over (cells, features)
        X = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=tuple(g["shape"][:])[::-1])
        name = np.array([x.decode() for x in g["features"]["name"][:]])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        bar = np.array([x.decode() for x in g["barcodes"][:]])
    keep = ftype == GENE_FEATURE
    if not keep.any():
        raise XeniumIOError(f"{p} carries no '{GENE_FEATURE}' features")
    return bar, list(name[keep]), X[:, keep].tocsr()


def _h5ad_path(donor_id: str, roi: str = "panc", h5ad_dir: Path | None = None) -> Path | None:
    tag = "Panc" if roi == "panc" else "pLN"
    d = H5AD_DIR if h5ad_dir is None else Path(h5ad_dir)
    if not d.is_dir():
        return None
    for p in sorted(d.glob(f"*_{tag}Xenium.h5ad")):
        stem = p.name.split("_")[0]
        if stem.lstrip("HDL").rstrip("abc") == donor_id:
            return p
    return None


def h5ad_refusal(donor_id: str, roi: str, n_candidates: int) -> str:
    """Why the h5ad may not be paired with this bundle, or "" if it may.

    ``<donor>_PancXenium.h5ad`` carries no serial, so for a donor+roi with more than one
    bundle there is nothing in the filename, the obs, or the path that says WHICH section it
    holds. 6533 and 112 both have one — they are 2024 cohort donors — so a Pheno-Xeno run
    that preferred the h5ad would take expression from the 2024 section and coordinates from
    the 2026 bundle and join them on ``cell_id``. That join does not fail loudly: the ids
    overlap by construction, so it would return some cells and look ordinary.
    """
    if n_candidates > 1:
        return (f"donor {donor_id} roi={roi} carries {n_candidates} bundles and the h5ad "
                f"filename records no serial, so which section it holds is unknown; read "
                f"the named bundle instead of joining two sections on cell_id")
    return ""


def load_donor(donor_id: str, xenium_paths_csv: Path, *, roi: str = "panc",
               prefer: str = "h5ad", serial: str | None = None,
               bundle: Path | None = None) -> DonorData:
    """Counts + microns for one donor, from whichever source exists.

    ``prefer`` picks the source when both exist. It changes which cells are present (the
    h5ad is the upstream QC subset) but not the values — see the module docstring.

    ``serial`` names the section. It is REQUIRED for a donor+roi carrying more than one
    bundle (6533 and 112 at ``roi=panc``); without it the call raises
    :class:`BundleAmbiguityError` rather than returning one of them.
    """
    ref = resolve_bundle(donor_id, xenium_paths_csv, roi=roi, serial=serial)
    # ``bundle`` overrides where the files are, never WHICH section they are. The recorded
    # path in xenium_paths.csv is relative to a HuPrime5kPanTissue/<batch>/<run>/ tree that no
    # site lays out in full, and XENIUM_ROOT is this module's workstation default. hpg/sites.py
    # already owns per-site layout and fails closed; callers that consult it pass the bundle
    # it resolved so the two do not disagree. Section identity still comes from ``serial``.
    bundle_dir = Path(bundle) if bundle is not None else ref.path
    cells_p = bundle_dir / "cells.parquet"
    if not cells_p.exists():
        raise XeniumIOError(f"{cells_p} is absent; spatial coordinates live only there")
    cells = pd.read_parquet(cells_p, columns=["cell_id", "x_centroid", "y_centroid"])
    cells["cell_id"] = cells.cell_id.astype(str)

    refused = h5ad_refusal(donor_id, roi, ref.n_candidates)
    h5 = _h5ad_path(donor_id, roi) if (prefer == "h5ad" and not refused) else None
    if h5 is not None:
        import anndata as ad
        a = ad.read_h5ad(h5)
        ids = np.array(list(map(str, a.obs_names)))
        genes = list(map(str, a.var_names))
        counts = sp.csr_matrix(a.X)
        source = "h5ad"
    else:
        ids, genes, counts = _read_bundle_matrix(bundle_dir)
        source = "bundle"

    pos = pd.Series(np.arange(len(cells)), index=cells.cell_id.values)
    hit = pos.reindex(ids)
    ok = hit.notna().to_numpy()
    if not ok.any():
        raise XeniumIOError(f"donor {donor_id}: no cell_id from {source} matched cells.parquet")
    rows = hit[ok].to_numpy().astype(int)
    xy = cells.loc[:, ["x_centroid", "y_centroid"]].to_numpy()[rows]

    return DonorData(donor_id=donor_id, cell_id=ids[ok], genes=genes,
                     counts=counts[ok].tocsr(), xy=xy, source=source,
                     n_dropped_no_spatial=int((~ok).sum()),
                     serial=ref.serial, bundle=bundle_dir,
                     h5ad_refused_reason=refused if source == "bundle" else "")


def lognorm(counts: sp.csr_matrix, *, target: float | None = None) -> sp.csr_matrix:
    """Per-cell depth normalisation to `target` (default: the median depth) then log1p.

    For an argmax across lineages within a cell this is nearly a no-op — measured on 6539,
    raw and lognorm agree to within 0.2 percentage points on every lineage share, because a
    per-cell scalar cancels in a within-cell comparison. It matters here because the
    embedding compares cells to each other, where depth is a genuine confound.
    """
    depth = np.asarray(counts.sum(axis=1)).ravel()
    if target is None:
        target = float(np.median(depth[depth > 0])) if (depth > 0).any() else 1.0
    scale = np.divide(target, depth, out=np.zeros_like(depth, dtype=float), where=depth > 0)
    return sp.csr_matrix(counts.multiply(scale[:, None]).log1p())
