#!/usr/bin/env python3
"""Rebuild one donor's per-cell table from a revised QuPath export, REUSING REDSEA.

WHY THIS EXISTS INSTEAD OF A RE-INGEST. Written while the pipeline was still content-addressed:
`run_id` was a hash over the whole cohort manifest, so repointing three donors at their revised
exports re-keyed the run and forced all twenty donors through ingest, geometry, REDSEA and
expression again -- about 2.3 hours to recompute byte-identical results for seventeen donors that
did not change.

CORRECTED 2026-09-18: THAT FRAMEWORK WAS REMOVED ON 2026-09-16. There is no `run_id` and no
re-keying; stages write to stable roots and a stage whose output is present is skipped, so a
re-ingest would now cost only the donors actually rerun. The reason this module still exists is
the narrower one below -- the revised export did not re-segment anything, so REDSEA and the
UniFORM alignment can be REUSED per cell rather than recomputed -- and the twenty-donor cohort
output this lane reads still lives only under the historical `data/runs/<run_id>/` tree.

THE KEY FACT THAT MAKES REUSE CORRECT: the revised export did not re-segment anything. QuPath
scripts 0b and 1-4 re-ran against a NEW TISSUE ANNOTATION -- 0b removed the duplicate detections,
1_panc_objects re-derived the Islet annotations inside the new tissue boundary, and 3_signed_distance
recomputed the distances. Cell DETECTIONS keep their UUIDs from the original detection run (measured:
99.5-99.8% of revised object ids are already in the run, and the cell outline geojson is a SUPERSET
of every revised export). A cell that survives de-duplication is the same object, with the same
polygon, so its REDSEA-corrected and UniFORM-aligned marker values are unchanged. They are joined
here by `object_id` -- never centroid, never row order.

WHAT IS RECOMPUTED, because the new tissue annotation changed it:
  parent / islet_num     from the export's `Parent`
  cell_region            from the export's re-run signed distance, at the 20 um rule
  qc_analysis_eligible   from the geometry thresholds in config.ini
  qc_estimation_eligible the same, plus the no-cytoplasm rule, which gates ESTIMATION only

WHAT CANNOT BE RECOMPUTED HERE, and is therefore dropped from the eligibility rule: the three
raster checks (`raster_area_mismatch`, `raster_empty`, `nucleus_raster_empty`) need the whole-slide
label raster. On 6534 they excluded 51,517 cells in the old run, and `raster_area_mismatch` was
overwhelmingly the DUPLICATE artefact -- the loser of each pair lost every pixel to the survivor.
With the duplicates gone at source that reason should be largely empty, but "should be" is not
"measured", so the column records which rules were actually applied.

THE BACKLOG IS DISCOVERED, NOT NAMED. `data/RevisedQuPathExports` is the first place a donor's
cell matrix is looked for, and the maintainer drops new exports into it while this runs. A list of
donor ids in this file's argv would go stale the moment one lands, exactly as `cohort._rebuilt()`'s
literal tuple did. So the exports are read, each is resolved to a donor, and a donor is in the
backlog when its newest export is newer than its rebuilt partition -- or when it has none.

DONOR IDENTITY COMES FROM THE `Image` COLUMN, NEVER FROM THE FILENAME. The filenames do not agree
with the donor ids: `0112.tsv` and `HDL112_with_pLN.tsv` are the same image (0081366_HDL112_...),
`HDL115_with_pLN.tsv` is donor `115`, and `6533.tsv` and `6533_with_pLN.tsv` are one image too.
Parsing any of those spellings would be donor-name parsing, which this repo forbids for exactly
this reason. The `Image` value is matched against two tracked manifests -- `qupath_manifest.json`
and `data/donor_acquisition_runs.csv` -- and an export whose image neither names is REFUSED with
its image printed, never guessed at from its name.

`--list` prints that backlog with mtimes and does nothing else.

THE `*_with_pLN` EXPORTS CARRY A SECOND TISSUE, AND IT DOES NOT GO IN THE PANCREAS PARTITION.
Measured 2026-09-18 over the whole of each file (not a prefix):

    6533_with_pLN.tsv    Annotation (Tissue) 931,910  Islet_* 95,566  Annotation (pLN)  34,841  root 664
    6533.tsv             Annotation (Tissue) 931,910  Islet_* 95,566                            root 664
    HDL112_with_pLN.tsv  Annotation (Tissue) 641,639  Islet_* 117,166 Annotation (pLN) 617,396
    0112.tsv             Annotation (Tissue) 758,805  Islet_*       0
    HDL115_with_pLN.tsv  Annotation (Tissue) 656,835  Islet_* 129,699 Annotation (pLN) 261,761

`6533.tsv` is `6533_with_pLN.tsv` minus its pLN rows, to the cell. So dropping the pLN-parented
rows makes the newer export describe exactly the pancreas the older one did, and the count is
printed rather than absorbed. (`0112.tsv` is the odd one: 641,639 + 117,166 = 758,805, so its
islet cells are parented to `Tissue` -- that export has no islet annotations at all, and the
`_with_pLN` one supersedes it on substance and not only on mtime.)

The per-image annotation census is written to `revised_cells/section_parents.json`, which is the
evidence `scripts/senior/build_sections_layer.py` derives `section_relationship` from. It is a
record of which annotation classes the image carries, and it licenses no distance: pancreas and
lymphoid tissue in one image are two tissues that share a pixel grid.

    python core/phenocycler/scripts/typing/rebuild_donor_from_revised_export.py --list
    python core/phenocycler/scripts/typing/rebuild_donor_from_revised_export.py          # backlog
    python core/phenocycler/scripts/typing/rebuild_donor_from_revised_export.py 6505 6534
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S

sys.path.insert(0, str(S.REPO / "core/phenocycler"))
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

REVISED = S.REPO / "data/RevisedQuPathExports"
OUT = S.OUT / "revised_cells"

#: The two tracked image -> donor manifests. Neither is complete on its own: the QuPath manifest
#: names 0112 and 6533 and not 6579; the acquisition table names 6579 and not 0112 or 6533. They
#: are read together and any disagreement about one image is a refusal, not a preference.
MANIFEST = S.REPO / "core/phenocycler/data/qupath_manifest.json"
ACQUISITION = S.REPO / "data/donor_acquisition_runs.csv"

#: The per-image annotation census this lane writes for the sections layer.
SECTION_PARENTS = OUT / "section_parents.json"

#: QuPath appends the pyramid level to the project entry name in `Image`. Stripped so the bare
#: file name in either manifest matches an export written with or without it.
RESOLUTION_SUFFIX = re.compile(r" - resolution #\d+$")

#: `Parent` shapes, measured on every export in the directory.
ANNOTATION_PARENT = re.compile(r"^Annotation \((?P<name>.*)\)$")
ISLET_PARENT = re.compile(r"^Islet[_ ]?\d+$")
ROOT_PARENT = re.compile(r"^Root object")

#: An annotation class that delimits lymphoid tissue rather than pancreas. Matched on the class
#: name recorded in the export, case-insensitively -- `Annotation (pLN)` is the only spelling
#: measured here, and every class actually seen is written to section_parents.json, so a new one
#: shows up as an unclassified name instead of being silently counted as pancreas.
PLN_ANNOTATION = re.compile(r"^p\s*ln", re.IGNORECASE)

#: from core/phenocycler/config.ini [geometry_qc]. Read, not guessed -- see load_thresholds().
CONFIG = S.REPO / "core/phenocycler/config.ini"

#: Every islet expansion is 20 um. Derived from the signed distance rather than read from an
#: annotation, because 6482's annotated expansion measures 23.1 um and that is an error.
EXPANSION_UM = 20.0

IDENT = {
    "Object ID": "object_id",
    "Image": "image",
    "Parent": "parent_raw",
    "Classification": "qupath_classification",
    "Centroid X µm": "X_centroid",
    "Centroid Y µm": "Y_centroid",
}
SHAPE = {
    "Nucleus: Area µm^2": "nucleus_area_um2",
    "Cell: Area µm^2": "cell_area",
    "Cell: Solidity": "cell_solidity",
    "Nucleus/Cell area ratio": "nucleus_cell_area_ratio",
}
#: 3_signed_distance.groovy writes "with Islet"; an older run of it wrote the form without "with".
#: Prefer the new one and fall back, rather than assuming which export is in hand.
DIST = ["Signed distance to annotation with Islet µm", "Signed distance to annotation Islet µm"]


def load_thresholds() -> dict[str, float]:
    import configparser
    cp = configparser.ConfigParser()
    cp.read(CONFIG)
    g = cp["geometry_qc"]
    return {
        "min_cell_area": g.getfloat("min_cell_area"),
        "max_area_median_multiple": g.getfloat("max_area_median_multiple"),
        "min_cell_solidity": g.getfloat("min_cell_solidity"),
        "no_cytoplasm_nc_ratio": g.getfloat("no_cytoplasm_nc_ratio"),
    }


def read_export(path: pathlib.Path) -> pd.DataFrame:
    """Only the columns needed. The exports are 8-17 GB and 994 columns wide."""
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    index = {n: i for i, n in enumerate(header)}
    dist_col = next((c for c in DIST if c in index), None)
    if dist_col is None:
        raise KeyError(f"{path.name}: no signed-distance-to-Islet column; 3_signed_distance did not run")

    wanted = {**IDENT, **SHAPE, dist_col: "signed_dist_islet_um"}
    missing = [c for c in wanted if c not in index]
    if missing:
        raise KeyError(f"{path.name}: missing {missing}")

    picks = [index[c] for c in wanted]
    names = [wanted[c] for c in wanted]
    rows = []
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            rows.append([row[i] for i in picks])
    df = pd.DataFrame(rows, columns=names)
    for c in list(SHAPE.values()) + ["X_centroid", "Y_centroid", "signed_dist_islet_um"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def image_key(value: str) -> str:
    """The `Image` value reduced to the bare file name both manifests record."""
    return RESOLUTION_SUFFIX.sub("", str(value).strip())


def image_to_donor() -> dict[str, tuple[str, str]]:
    """`Image` file name -> (donor_id, which manifest said so).

    Both tracked manifests are read because neither covers the directory on its own, and a
    disagreement about one image raises rather than resolving to whichever was read last.
    """
    mapping: dict[str, tuple[str, str]] = {}

    def record(key: str, donor: str, source: str) -> None:
        key = image_key(key)
        if not key or not donor:
            return
        previous = mapping.get(key)
        if previous and previous[0] != donor:
            raise ValueError(
                f"{key}: {previous[1]} calls this donor {previous[0]} and {source} calls it "
                f"{donor}. An image belongs to one donor; refusing to pick.")
        mapping.setdefault(key, (donor, source))

    raw = json.loads(MANIFEST.read_text())
    for row in raw["images"]:
        record(row["image_id"], str(row["donor_id"]), MANIFEST.name)
    if ACQUISITION.exists():
        with ACQUISITION.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                record(row["image_file"], str(row["donor_id"]), ACQUISITION.name)
    return mapping


def read_first_image(path: pathlib.Path) -> str:
    """The `Image` value of the export's first data row, read without loading the file."""
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if "Image" not in header:
            raise KeyError(f"{path.name}: no Image column; donor cannot be resolved")
        column = header.index("Image")
        for row in reader:
            if len(row) > column and row[column].strip():
                return row[column].strip()
    raise ValueError(f"{path.name}: no data rows; donor cannot be resolved")


def partition_mtime(donor: str) -> float | None:
    target = OUT / f"donor_id={donor}/data.parquet"
    return target.stat().st_mtime if target.exists() else None


def survey() -> list[dict]:
    """One record per export in `RevisedQuPathExports`, with why it is or is not in the backlog.

    Nothing here is inferred from a file name.  The donor comes from the export's own `Image`
    value; the marker values it would have to join to come from `celltype/uniform`, and an
    export whose donor has no aligned partition is refused with that named, because joining on
    `object_id` needs rows on both sides and a missing measurement is unavailable, not zero.
    """
    mapping = image_to_donor()
    records = []
    for path in sorted(REVISED.glob("*.tsv")):
        record = {"tsv": path.name, "path": path, "bytes": path.stat().st_size,
                  "mtime": path.stat().st_mtime, "donor": None, "image": None,
                  "manifest": None, "state": None, "reason": ""}
        try:
            record["image"] = read_first_image(path)
        except (KeyError, ValueError, OSError) as exc:
            record["state"] = "refused"
            record["reason"] = f"unreadable: {exc}"
            records.append(record)
            continue
        resolved = mapping.get(image_key(record["image"]))
        if resolved is None:
            record["state"] = "refused"
            record["reason"] = (f"image {record['image']!r} is in neither {MANIFEST.name} nor "
                                f"{ACQUISITION.name}; record it there rather than parsing the "
                                "file name")
            records.append(record)
            continue
        record["donor"], record["manifest"] = resolved
        records.append(record)

    newest: dict[str, dict] = {}
    for record in records:
        if record["donor"] is None:
            continue
        held = newest.get(record["donor"])
        if held is None or record["mtime"] > held["mtime"]:
            newest[record["donor"]] = record

    for record in records:
        if record["state"] == "refused":
            continue
        donor = record["donor"]
        if donor in DONOR_EXCLUSIONS:
            record["state"] = "refused"
            record["reason"] = f"donor excluded: {DONOR_EXCLUSIONS[donor]}"
        elif newest[donor] is not record:
            record["state"] = "superseded"
            record["reason"] = f"a newer export for {donor}: {newest[donor]['tsv']}"
        elif not (S.ALIGNED / f"donor_id={donor}/data.parquet").exists():
            record["state"] = "refused"
            record["reason"] = (f"no aligned partition at {S.ALIGNED}/donor_id={donor}; this "
                                "donor has no REDSEA/UniFORM marker values to join to")
        else:
            built = partition_mtime(donor)
            if built is None:
                record["state"] = "backlog"
                record["reason"] = "no rebuilt partition"
            elif built < record["mtime"]:
                record["state"] = "backlog"
                record["reason"] = "the export is newer than the rebuilt partition"
            else:
                record["state"] = "fresh"
                record["reason"] = "the rebuilt partition is newer than the export"
    return records


def stamp(epoch: float | None) -> str:
    return "-" if epoch is None else time.strftime("%Y-%m-%d %H:%M", time.localtime(epoch))


def print_survey(records: list[dict]) -> None:
    print(f"{REVISED}  ({len(records)} exports)\n")
    print(f"  {'export':<24} {'GiB':>7}  {'written':<16} {'donor':>6}  "
          f"{'rebuilt':<16} {'state':<10} why")
    for record in sorted(records, key=lambda r: (r["state"] or "", r["tsv"])):
        donor = record["donor"] or "?"
        built = stamp(partition_mtime(donor)) if record["donor"] else "-"
        print(f"  {record['tsv']:<24} {record['bytes'] / 2**30:7.2f}  "
              f"{stamp(record['mtime']):<16} {donor:>6}  {built:<16} "
              f"{record['state']:<10} {record['reason']}")
    backlog = [r["donor"] for r in records if r["state"] == "backlog"]
    print(f"\n  backlog: {len(backlog)} donor(s): {' '.join(backlog) if backlog else 'none'}")


def classify_parents(parent: pd.Series) -> tuple[pd.Series, dict]:
    """(is_pln_annotation, census) from the export's `Parent` column.

    The census is what `section_relationship` is derived from downstream, so it records the
    annotation classes by name.  Nothing is inferred from an absence: a class that does not
    match the pLN spelling is counted under its own name, not folded into "pancreas".
    """
    text = parent.fillna("").astype(str).str.strip()
    annotation = text.str.extract(ANNOTATION_PARENT, expand=False)
    counts = annotation.dropna().value_counts().to_dict()
    census = {
        "annotation_parent_counts": {str(k): int(v) for k, v in sorted(counts.items())},
        "islet_parent_rows": int(text.str.match(ISLET_PARENT).sum()),
        "root_object_rows": int(text.str.match(ROOT_PARENT).sum()),
        "rows_in_export": int(len(text)),
    }
    census["pln_annotation_classes"] = sorted(
        name for name in census["annotation_parent_counts"] if PLN_ANNOTATION.match(name))
    census["other_parent_rows"] = int(
        len(text) - annotation.notna().sum() - census["islet_parent_rows"]
        - census["root_object_rows"])
    is_pln = annotation.notna() & annotation.fillna("").str.match(PLN_ANNOTATION)
    return is_pln.fillna(False), census


def record_section_parents(donor: str, record: dict, census: dict) -> None:
    """Merge this image's annotation census into revised_cells/section_parents.json."""
    payload = {"donors": {}}
    if SECTION_PARENTS.exists():
        try:
            payload = json.loads(SECTION_PARENTS.read_text())
        except json.JSONDecodeError:
            payload = {"donors": {}}
    payload.setdefault("donors", {})
    payload["evidence"] = (
        "QuPath Parent annotation classes per image, counted over the whole revised export. "
        "A filtering and sorting attribute only: two tissues in one image share a pixel grid "
        "and no anatomy, and this file licenses no cross-tissue distance.")
    payload["generated_at_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    payload["donors"][donor] = {
        "image": record["image"],
        "source_tsv": record["tsv"],
        "source_mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(record["mtime"])),
        **census,
    }
    SECTION_PARENTS.parent.mkdir(parents=True, exist_ok=True)
    SECTION_PARENTS.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def islet_number(parent: pd.Series) -> pd.Series:
    """`Islet_269` -> 269. Anything that is not an islet annotation is NA, not 0."""
    return pd.to_numeric(parent.str.extract(r"Islet[_ ]?(\d+)", expand=False), errors="coerce")


def region(dist: pd.Series) -> pd.Series:
    """core / expansion / tissue, from the signed distance alone.

    Negative or zero is inside the islet; (0, 20] is the peri-islet expansion; beyond that is
    ordinary tissue. Built this way so every expansion is exactly 20 um by construction.
    """
    out = pd.Series(pd.NA, index=dist.index, dtype="object")
    finite = np.isfinite(dist)
    out[finite & (dist <= 0)] = "core"
    out[finite & (dist > 0) & (dist <= EXPANSION_UM)] = "expansion"
    out[finite & (dist > EXPANSION_UM)] = "tissue"
    return out


def eligibility(df: pd.DataFrame, thr: dict[str, float]):
    """(analysis_eligible, analysis_reason, estimation_eligible, estimation_reason).

    THE PIPELINE KEEPS THREE SEPARATE ROLES -- analysis, estimation, spillover context -- and a
    rule that disqualifies a cell for one does NOT disqualify it for the others. Collapsing them
    into one flag is not a conservative simplification, it is a different answer: the first
    version of this function applied `no_cytoplasm` to ANALYSIS and excluded 100,059 cells for
    6534 that the run counts as analysis-eligible. Measured in the run's own geometry QC:

        analysis_reason            no_cytoplasm = 0
        estimation_reason          no_cytoplasm = 99,024
        spillover_context_reason   no_cytoplasm = 0

    A nucleus that fills its whole cell gives you nothing to ESTIMATE a background level from;
    it is still a cell, and still counts.

    A cell with no nucleus is excluded from every role -- we cannot be sure it is a cell at all.
    Order matters: the reason recorded is the FIRST rule the cell failed, matching how the
    geometry QC reports one reason per cell.
    """
    area = df["cell_area"]
    med = np.nanmedian(area[np.isfinite(area)])
    shared = [
        ("missing_nucleus_geometry", ~np.isfinite(df["nucleus_area_um2"])),
        ("cell_area_below_min", area < thr["min_cell_area"]),
        ("cell_area_above_max", area > thr["max_area_median_multiple"] * med),
        ("low_solidity", df["cell_solidity"] < thr["min_cell_solidity"]),
    ]

    def apply(checks):
        reason = pd.Series("", index=df.index, dtype="object")
        for name, failed in checks:
            reason[(reason == "") & failed.fillna(False)] = name
        return reason == "", reason

    a_ok, a_why = apply(shared)
    e_ok, e_why = apply(shared + [
        ("no_cytoplasm", df["nucleus_cell_area_ratio"] > thr["no_cytoplasm_nc_ratio"]),
    ])
    return a_ok, a_why, e_ok, e_why


def rebuild(record: dict, thr: dict[str, float]) -> None:
    started = time.time()
    donor, export = record["donor"], record["path"]

    df = read_export(export)
    seen = sorted(df["image"].dropna().astype(str).str.strip().unique())
    if len(seen) != 1:
        raise ValueError(f"{export.name}: {len(seen)} distinct Image values {seen}; this lane "
                         "rebuilds one image into one donor partition and will not split them")
    if image_key(seen[0]) != image_key(record["image"]):
        raise ValueError(f"{export.name}: resolved on {record['image']!r} but the file carries "
                         f"{seen[0]!r}")

    # THE SECOND TISSUE COMES OUT BEFORE ANYTHING IS DERIVED FROM IT. A cell parented to a pLN
    # annotation is in a lymph node: it has no islet to be at a distance from, and leaving it in
    # would put lymphoid cells in the pancreas denominator with nothing saying so. The census is
    # recorded either way, so an image with no pLN annotation is a measured zero.
    is_pln, census = classify_parents(df["parent_raw"])
    n_pln = int(is_pln.sum())
    if n_pln:
        df = df[~is_pln].reset_index(drop=True)
    census["rows_written_to_pancreas_partition"] = int(len(df))
    record_section_parents(donor, record, census)

    df.insert(0, "donor_id", donor)
    df["islet_num"] = islet_number(df["parent_raw"])
    df["cell_region"] = region(df["signed_dist_islet_um"])
    (df["qc_analysis_eligible"], df["analysis_reason"],
     df["qc_estimation_eligible"], df["estimation_reason"]) = eligibility(df, thr)

    # ---- REDSEA-corrected, UniFORM-aligned marker values, BY object_id -------------------
    aligned_path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    aligned = pq.read_table(aligned_path).to_pandas()
    markers = [c for c in aligned.columns
               if c not in {"donor_id", "object_id", "qc_analysis_eligible",
                            "X_centroid", "Y_centroid", "cell_region", "islet_num"}]
    merged = df.merge(aligned[["object_id"] + markers], on="object_id",
                      how="left", validate="one_to_one", indicator="_merge")

    # A revised cell with no aligned row has NO marker measurement. That is unavailable, not
    # zero -- the markers stay NaN and the count is reported rather than quietly dropped.
    #
    # COUNTED FROM THE MERGE INDICATOR, not from whether a marker column is NaN. The first
    # attempt probed markers[0] and reported 191,322 missing rows for 6534 when the true
    # number is 6,182: a marker is legitimately NaN for a cell that HAS a row, because not
    # every marker was acquired on every donor and REDSEA can leave one unset. "This column
    # is empty" and "this cell has no measurements" are different questions.
    unmeasured = int((merged["_merge"] == "left_only").sum())
    merged = merged.drop(columns="_merge")

    target = OUT / f"donor_id={donor}"
    target.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(target / "data.parquet", index=False)

    n = len(merged)
    elig = int(merged["qc_analysis_eligible"].sum())
    by_region = merged.loc[merged["qc_analysis_eligible"], "cell_region"].value_counts()
    pln = f" | pLN rows set aside {n_pln:>7,}" if n_pln else ""
    print(f"  {donor}: {n:>9,} cells | eligible {elig:>9,} ({100*elig/n:4.1f}%) | "
          f"core {by_region.get('core', 0):>7,} expansion {by_region.get('expansion', 0):>7,} "
          f"tissue {by_region.get('tissue', 0):>9,} | "
          f"no marker row {unmeasured:>6,} | islets {int(merged['islet_num'].nunique()):>5,}"
          f"{pln} | from {record['tsv']} | {time.time()-started:.0f}s", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("donors", nargs="*", default=[],
                    help="donor ids to rebuild; default is the discovered backlog")
    ap.add_argument("--list", action="store_true", dest="list_only",
                    help="print the backlog with mtimes and do nothing else")
    args = ap.parse_args()

    records = survey()
    if args.list_only:
        print_survey(records)
        return 0

    print_survey(records)
    by_donor = {r["donor"]: r for r in records if r["state"] in {"backlog", "fresh"}}
    if args.donors:
        chosen, missing = [], []
        for donor in args.donors:
            if donor in by_donor:
                chosen.append(by_donor[donor])
            else:
                missing.append(donor)
        if missing:
            refused = {r["donor"]: r["reason"] for r in records if r["state"] == "refused"}
            detail = "; ".join(f"{d}: {refused.get(d, 'no usable export in ' + str(REVISED))}"
                               for d in missing)
            raise SystemExit(f"[err] refusing to rebuild {missing}: {detail}")
    else:
        chosen = [r for r in records if r["state"] == "backlog"]
        if not chosen:
            print("\nnothing to rebuild: every resolvable export is already rebuilt")
            return 0

    thr = load_thresholds()
    print(f"\ngeometry thresholds from config.ini: {thr}")
    print(f"expansion fixed at {EXPANSION_UM:.0f} um, derived from the signed distance\n")
    OUT.mkdir(parents=True, exist_ok=True)
    for record in chosen:
        rebuild(record, thr)
    print(f"\n-> {OUT}")
    print(f"-> {SECTION_PARENTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
