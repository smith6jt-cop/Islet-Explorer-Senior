#!/usr/bin/env python3
"""Build the per-donor clinical table the app's Inflammation tab reads.

This script is named as the regeneration path in `data/app_data/DATA_PROVENANCE.md` and in the
Inflammation module's own on-screen error text, and it has never existed -- not in the tree and
not anywhere in git history.  So the table it was supposed to produce could not be refreshed and
the tab's nPOD panel has been permanently degraded.  This is that script.

Two sources, joined on the donor id:

  * `donor_metadata_panc.xlsx` -- the cohort's own sheet: status, age, sex, autoantibody count,
    T1D duration.  This is authoritative for cohort membership and for `disease.status`.
  * the nPOD "All-Download" export -- BMI, HbA1c, C-peptide and the four autoantibody assays,
    for every nPOD donor.  Deliberately gitignored: it carries free-text admission courses,
    clinical histories, medication lists and HLA types for 741 donors, of which this cohort
    needs a dozen numeric fields for twenty.  Only the derived CSV is committed.

Case-id mapping: donor `115` is filed upstream as `HDL115`; the script fails loudly rather than
dropping a donor it cannot resolve, because a silently missing donor is a silently missing row in
every downstream per-donor statistic.

`C_peptide_ng_mL` arrives censored on some donors (`<0.02`); the value and the censoring are kept
in separate columns rather than coercing the string to a number or the number to a string.

    python scripts/build_donor_metadata.py --out <domain data root>/donor_metadata.csv
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))
from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

COHORT_XLSX = pathlib.Path.home() / "IO60panc2nd" / "donor_metadata_panc.xlsx"
NPOD_GLOB = "All-Download_From_nPOD_*.xlsx"

#: nPOD column -> the name the app reads.
NPOD_FIELDS = {
    "age_years": "age_years_npod",
    "diabetes_hx_years": "t1d_duration_years",
    "BMI": "BMI",
    "HbA1c_percent": "HbA1c_percent",
    "GADA_Result": "AAb_GADA",
    "IA_2A_Result": "AAb_IA2A",
    "mIAA_Result": "AAb_mIAA",
    "ZnT8A_Result": "AAb_ZnT8A",
}

#: The app orders donor status ND < Aab+ < T1D; the cohort sheet spells the middle one "AAB".
STATUS = {"ND": "ND", "AAB": "Aab+", "T1D": "T1D"}


def _positive(value) -> bool | None:
    """An nPOD autoantibody result as a logical, or None when it was not assayed."""
    if pd.isna(value):
        return None
    return str(value).strip().lower() in {"positive", "pos", "true", "yes", "1"}


def _c_peptide(value) -> tuple[float | None, bool]:
    """(value, censored).  `<0.02` is a limit of detection, not the number 0.02."""
    if pd.isna(value):
        return None, False
    text = str(value).strip()
    if text.startswith("<"):
        return float(re.sub(r"^<\s*", "", text)), True
    try:
        return float(text), False
    except ValueError:
        return None, False


def _npod_index(npod: pd.DataFrame) -> dict[str, int]:
    """donor id -> row, tolerating the `HDL` prefix some cases carry upstream."""
    out: dict[str, int] = {}
    for position, case in enumerate(npod["case_id"].astype(str)):
        out[case] = position
        digits = re.sub(r"^\D+", "", case)
        if digits and digits not in out:
            out[digits] = position
    return out


def build(cohort_xlsx: pathlib.Path, npod_xlsx: pathlib.Path) -> pd.DataFrame:
    cohort = pd.read_excel(cohort_xlsx)
    cohort["donor_id"] = cohort["donor_id"].astype(str)
    cohort = cohort[~cohort["donor_id"].isin(DONOR_EXCLUSIONS)].copy()

    npod = pd.read_excel(npod_xlsx, sheet_name="data")
    index = _npod_index(npod)

    unresolved = [d for d in cohort["donor_id"] if d not in index]
    if unresolved:
        raise KeyError(
            f"{len(unresolved)} donor(s) not found in the nPOD export: {unresolved}. "
            "A dropped donor is a silently missing row in every per-donor statistic.")

    rows = []
    for _, donor in cohort.iterrows():
        source = npod.iloc[index[donor["donor_id"]]]
        peptide, censored = _c_peptide(source.get("C_peptide_ng_mL"))
        row = {
            "case_id": donor["donor_id"],
            "npod_case_id": str(source["case_id"]),
            "donor_status": STATUS.get(str(donor["disease.status"]), donor["disease.status"]),
            "age": donor["age"],
            "sex": donor["sex"],
            "n_aab": donor["aab_count"],
            "duration_years_cohort": donor.get("duration"),
            "C_peptide_ng_mL": peptide,
            "C_peptide_censored": censored,
        }
        for source_name, out_name in NPOD_FIELDS.items():
            value = source.get(source_name)
            row[out_name] = _positive(value) if out_name.startswith("AAb_") else value
        rows.append(row)

    frame = pd.DataFrame(rows)
    # The app bins age at 13 years; keep the bin beside the value so a plot cannot re-derive it
    # differently from the table.
    frame["age_group"] = pd.cut(frame["age"], [-1, 13, 200], labels=["<13", ">=13"])
    return frame.sort_values(["donor_status", "age"]).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=pathlib.Path, default=COHORT_XLSX)
    parser.add_argument("--npod", type=pathlib.Path, default=None)
    parser.add_argument("--out", type=pathlib.Path,
                        default=REPO / "data/development/typing/app/donor_metadata.csv")
    parser.add_argument("--all-imaged", action="store_true",
                        help="keep every donor in the cohort sheet rather than only the donors "
                             "the current run carries")
    args = parser.parse_args()

    npod = args.npod
    if npod is None:
        found = sorted(REPO.glob(NPOD_GLOB))
        if not found:
            raise FileNotFoundError(
                f"no nPOD export matching {NPOD_GLOB} in {REPO} -- it is gitignored on purpose; "
                "download a fresh one rather than committing it")
        npod = found[-1]

    frame = build(args.cohort, npod)
    if not args.all_imaged:
        # The cohort sheet covers every donor ever imaged; the app dataset covers the donors the
        # run actually carries.  Taking the run's own manifest as the authority means this table
        # cannot drift from the cells it describes.
        manifest = json.loads((REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
        in_run = {str(image["donor_id"]) for image in manifest["images"]}
        missing = sorted(in_run - set(frame["case_id"]))
        if missing:
            raise KeyError(f"run donors absent from the cohort sheet: {missing}")
        frame = frame[frame["case_id"].isin(in_run)].reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(f"{len(frame)} donors -> {args.out}")
    print(frame[["case_id", "donor_status", "age", "sex", "n_aab", "BMI",
                 "HbA1c_percent"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
