#!/usr/bin/env python3
"""
Build the donor acquisition-batch covariate table from the qptiff headers.

Acquisition batch is aliased with disease status at the COARSE level in this cohort: no ND donor
was acquired on Fusion 2.4.1. The Fusion version therefore cannot be corrected on -- correcting
for it would remove the effect being measured -- it can only be disclosed, and it has to exist as
data to be disclosed at all. That is what this table is for.

Correction, where it is warranted, happens on the FINER variables instead: ``stain_group`` and
``acq_run`` in ``data/donor_acquisition_runs.csv``, which ``scripts/build_acquisition_runs.py``
recovers from the per-slide TIFF ``DateTime``. Four of its nine staining groups carry more than
one disease status, which is what makes correction possible at all. Where correction is not
warranted, stratify.

(Corrected 2026-09-18. This paragraph used to say the table had to exist "before Gate 5". The
Gate 2-5 apparatus is retired and on the never-revive list; the requirement is the maintainer's
standing rule on batch, not a gate. It also said the covariate "currently exists nowhere" --
both ``data/donor_acquisition_batch.csv`` and ``data/donor_acquisition_runs.csv`` exist now.)

Reads ONLY each qptiff's first-page description (no pixel data), so it is safe to run against /data.

An image the cohort specification declares takes its donor FROM the specification (the same rule
as ``build_acquisition_runs.py``): that is how HDL115, HDL112 and the same-slide sections
(`0080870_6533_pLNpanc_Scan3`, ...) get a row -- their names cannot be parsed, and until 2026-10-01
this generator stopped on them while the table on disk carried their rows by hand. An undeclared
file in --images is named `<donor_id>[panc]_Scan<n>.er.qptiff`, or is known non-cohort material and
skipped by name; any other filename stops the run rather than being parsed into a donor id.

    python scripts/build_acquisition_batch_table.py \
        --images /path/to/Images \
        --donor-metadata /path/to/donor_metadata_panc.xlsx \
        [--cohort-spec core/multiplex5000/cohort-spec.json]
    -> data/donor_acquisition_batch.csv
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
# The package lives in the `core/multiplex5000/` subtree (it was the `Phenocycler_Analysis` submodule
# until Phase 2 dissolved it), so the script runs without `pip install -e core/multiplex5000`.
# Index 0 pre-dates the subtree move and is kept: the checkout wins over any installed copy, so the
# exclusions enforced below are the ones readable in this tree rather than in some older wheel.
sys.path.insert(0, str(REPO / "core" / "multiplex5000"))
from phenocycler.cohort import DONOR_EXCLUSIONS                                # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
from build_acquisition_runs import COHORT_SPEC, declared_donors                # noqa: E402

# The 2026-07 cohort expansion: 15 donors acquired first, then these 7. A separate staining/imaging run
# with a physical panel change, i.e. a real acquisition batch. Cross-checked against panel_variant below.
BATCH2_DONORS = frozenset({"6450", "6505", "6523", "6534", "6566", "6591", "6623"})

# Export naming. A cohort scan is the donor id, an optional tissue suffix (batch 2 was named `panc`,
# with one `pacn` typo on 6534), then the scan number (6476 and 6579 are Scan2). The authoritative
# donor -> filename map is the array in senior/convert_viewer_ometiff.sh. Anchored on purpose: a
# partial match would invent a donor id out of an unrelated filename.
QPTIFF_NAME = r"(\d+)(?:panc|pacn)?_Scan\d+\.er\.qptiff"

# Known non-cohort material sharing the same image directory: HDL115pancLN is a lymph node from a
# subject outside the cohort. Recognised and skipped by name rather than parsed — the old loose
# pattern silently read it as donor "115".
NON_COHORT_NAME = r"[A-Za-z0-9]*LN_Scan\d+\.er\.qptiff"


# The same-slide PhenoCycler + Xenium sections: a slide number, the donor, `pLNpanc`. They were
# acquired on their own, on the 70-channel panel that carries BOTH CD38 and b-Catenin1, so they are
# a third acquisition batch -- recorded as `pLNpanc` in the table since 6533 and 112 joined it on
# 2026-09-16. The donor is never read out of this name; only the batch is.
SAME_SLIDE_NAME = r"\d+_(?:HDL)?\d+_pLNpanc_Scan\d+\.er\.qptiff"


def donor_for(base: str, declared: dict[str, str]) -> str | None:
    """The donor a qptiff belongs to. None = known non-cohort material, skipped. An unrecognised
    name stops the run: it is never parsed into a donor id."""
    if base in declared:
        return declared[base]
    if re.fullmatch(NON_COHORT_NAME, base):
        return None
    m = re.fullmatch(QPTIFF_NAME, base)                     # 6450panc_Scan1.er.qptiff -> 6450
    if not m:
        # Fail closed. The previous loose pattern took the first digit run after an optional
        # alpha prefix, so it mis-parsed rather than failed: Scan1.qptiff yielded donor "1" and
        # HDL115pancLN_Scan1.er.qptiff yielded "115" without anything declaring it. Anchoring the
        # whole name means an unrecognised file stops the run instead of inventing a donor.
        raise SystemExit(
            f"[err] qptiff filename is not declared in the cohort specification and does not match "
            f"the <donor_id>[panc]_Scan<n>.er.qptiff convention: {base!r}")
    return m.group(1)


def batch_for(donor: str, base: str) -> str:
    if re.fullmatch(SAME_SLIDE_NAME, base):
        return "pLNpanc"
    return "batch2" if donor in BATCH2_DONORS else "batch1"


#: What each batch's panel must be; a row that disagrees stops the run (the batch is hardcoded).
PANEL_OF_BATCH = {"batch1": {"CD38"}, "batch2": {"b_Catenin1", "neither"}, "pLNpanc": {"both"}}


def _first_page_description(path: str) -> str:
    import tifffile
    with tifffile.TiffFile(path) as t:
        return t.pages[0].description or ""


def _channel_names(path: str) -> list[str]:
    """Channel biomarker names, one per page, in acquisition order."""
    import tifffile
    names = []
    with tifffile.TiffFile(path) as t:
        for page in t.series[0].pages:
            d = page.description or ""
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", d)
            if m:
                names.append(m.group(1))
    return names


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--images", type=Path, required=True,
                    help="Directory of source .qptiff files (only first-page headers are read).")
    ap.add_argument("--donor-metadata", type=Path, required=True,
                    help="Excel file with donor_id and disease.status columns.")
    ap.add_argument("--cohort-spec", type=Path, default=COHORT_SPEC,
                    help="The cohort specification whose images take their donor from it.")
    ap.add_argument("--out", type=Path, default=REPO / "data" / "donor_acquisition_batch.csv")
    a = ap.parse_args(argv)
    declared = declared_donors(a.cohort_spec)

    # Validate the inputs up front, so a wrong path fails immediately and legibly rather than as a
    # KeyError from merging an empty frame (no .qptiff matches leaves the DataFrame with no columns)
    # or an unhandled FileNotFoundError from read_excel.
    if not a.images.is_dir():
        raise SystemExit(f"[err] --images is not a directory: {a.images}")
    if not a.donor_metadata.is_file():
        raise SystemExit(f"[err] --donor-metadata is not a file: {a.donor_metadata}")

    qptiffs = sorted(glob.glob(str(a.images / "*.qptiff")))
    if not qptiffs:
        raise SystemExit(f"[err] no .qptiff files found in {a.images}")

    rows = []
    for f in qptiffs:
        base = os.path.basename(f)
        donor = donor_for(base, declared)
        if donor is None:
            print(f"[skip] non-cohort sample: {base}", file=sys.stderr)
            continue
        desc = _first_page_description(f)
        soft = re.search(r"<AcquisitionSoftware>(.*?)</AcquisitionSoftware>", desc)
        chans = _channel_names(f)
        rows.append({
            "donor_id": donor,
            "image_file": base,
            "acquisition_software": soft.group(1) if soft else None,
            "n_channels": len(chans),
            "has_CD38": "CD38" in chans,
            "has_b_Catenin1": any(c.lower().startswith("b-catenin") or c.lower().startswith("b_catenin")
                                  for c in chans),
            "channels": ";".join(chans),
        })
    if not rows:
        raise SystemExit(f"[err] no cohort qptiff files in {a.images}; every match was non-cohort")

    t = pd.DataFrame(rows)

    # The spreadsheet writes HDL112 as 112 and the specification as 0112: one donor. The join is on
    # the id without leading zeroes; the table keeps the specification's spelling, which is the one
    # the partitions carry.
    meta = pd.read_excel(a.donor_metadata)[["donor_id", "disease.status"]]
    meta["_key"] = meta["donor_id"].astype(str).str.strip().str.lstrip("0")
    if meta["_key"].duplicated().any():
        raise SystemExit(f"[err] {a.donor_metadata}: a donor is listed twice")
    t["_key"] = t["donor_id"].str.lstrip("0")
    t = t.merge(meta[["_key", "disease.status"]], on="_key", how="left").drop(columns="_key")
    t["eligible"] = ~t["donor_id"].isin(DONOR_EXCLUSIONS)
    unknown = t.loc[t["eligible"] & t["disease.status"].isna(), "donor_id"].tolist()
    if unknown:                               # unavailable, never guessed -- and never silent
        print(f"[warn] no disease.status in {a.donor_metadata.name} for {unknown}", file=sys.stderr)
    # Panel variant: the batch-2 CD38 <-> b-Catenin1 swap, recorded rather than inferred downstream.
    t["panel_variant"] = t.apply(
        lambda r: "CD38" if r.has_CD38 and not r.has_b_Catenin1
        else ("b_Catenin1" if r.has_b_Catenin1 and not r.has_CD38
              else ("both" if r.has_CD38 and r.has_b_Catenin1 else "neither")), axis=1)
    # ACQUISITION BATCH — the covariate to analyse against. The two batches are separate staining and
    # imaging runs months apart with a physical panel change (batch 1 carries CD38, batch 2 carries
    # b-Catenin1 or neither), which is what makes them a batch rather than a label. `panel_variant`
    # and `acquisition_software` are each perfect proxies for it on this cohort — every batch-2 donor
    # is Fusion 2.4.1 and no batch-1 donor is — but the batch is the thing, not the software version.
    t["batch"] = [batch_for(d, f) for d, f in zip(t["donor_id"], t["image_file"])]
    off_panel = pd.Series([v not in PANEL_OF_BATCH[b] for v, b in zip(t["panel_variant"], t["batch"])], index=t.index)
    drift = t.loc[t.eligible & off_panel]
    if len(drift):                            # fail loudly rather than let the batch definition rot
        raise SystemExit(
            f"[err] panel_variant no longer agrees with the batch for {drift['donor_id'].tolist()} "
            f"({drift[['panel_variant', 'batch']].to_dict('records')}); "
            "the batch definition is hardcoded and must be re-derived if the cohort changes")

    t = t.sort_values("donor_id").reset_index(drop=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    t.drop(columns=["channels"]).to_csv(a.out, index=False)
    t.to_csv(a.out.with_name(a.out.stem + "_with_channels.csv"), index=False)

    e = t[t.eligible]
    ct = pd.crosstab(e["acquisition_software"], e["disease.status"])
    print(ct.to_string(), "\n")
    print(pd.crosstab(e["acquisition_software"], e["panel_variant"]).to_string(), "\n")

    # The specific confound: ND vs the rest, on Fusion 2.4.1 vs the rest.
    from scipy.stats import fisher_exact
    is241 = e["acquisition_software"].eq("Fusion 2.4.1")
    isnd = e["disease.status"].eq("ND")
    table = [[int((is241 & isnd).sum()), int((is241 & ~isnd).sum())],
             [int((~is241 & isnd).sum()), int((~is241 & ~isnd).sum())]]
    odds, p = fisher_exact(table)
    print(f"ND x Fusion-2.4.1 contingency (eligible donors only): {table}")
    print(f"Fisher exact p = {p:.4f}  -> batch is {'CONFOUNDED with' if p < 0.05 else 'independent of'} "
          "disease group; disclose it, do not correct it away")
    print(f"\n[done] {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
