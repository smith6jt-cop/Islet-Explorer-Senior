"""Protein analysis admission: eligibility, exact exclusions, and explicit Single statuses."""
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ENDOCRINE_UNCERTAINTY = frozenset({"Endocrine (unresolved subtype)",
                                 "Endocrine (ambiguous subtype)"})
DEFAULT_EXCLUSIONS = Path("/data/io60-calibration/pancreas_analysis_assignments_single_exclusions_20261005")


def donor_key(value: object) -> str:
    text = str(value)
    return str(int(text)) if text.isdecimal() else text


def protein_analysis_mask(cells: pd.DataFrame, donor: str,
                          exclusions: pd.DataFrame) -> pd.Series:
    """Keep ordinary uncertain endocrine; never admit a keyed or explicitly excluded Single.

    Uncertainty in an ordinary endocrine subtype is not an exclusion. Eligibility remains an
    AND of the available source/current masks. Historical Single uncertainty is excluded only
    with its explicit object class/status or its exact entry in the exclusion registry.
    """
    admitted = cells.cell_class.notna() & ~cells.cell_class.astype(str).eq("Excluded")
    for column in ("qc_analysis_eligible", "analysis_eligible"):
        if column in cells:
            admitted &= cells[column].fillna(False).astype(bool)
    if "resolved_by" in cells:
        admitted &= ~cells.resolved_by.astype(str).eq("excluded_single_hormone")
    if "endocrine_object_class" in cells:
        single = cells.endocrine_object_class.astype(str).eq("EndocrineSingle")
        unsupported = cells.cell_class.astype(str).isin(ENDOCRINE_UNCERTAINTY)
        if "endocrine_subtype_status" in cells:
            unsupported |= cells.endocrine_subtype_status.astype(str).isin(
                {"no_usable_hormone", "tied_highest_hormones"})
        admitted &= ~(single & unsupported)
    if len(exclusions):
        if "donor_id" not in exclusions:
            raise ValueError("protein exclusions require donor_id and an exact cell key")
        key = next((k for k in ("object_id", "cell_id") if k in cells and k in exclusions), None)
        if key is None:
            raise ValueError("protein exclusions share no exact cell key with the source")
        selected = exclusions.loc[exclusions.donor_id.map(donor_key).eq(donor_key(donor)), key]
        admitted &= ~cells[key].astype(str).isin(selected.astype(str))
    return admitted


def load_protein_exclusions(path: Path, donor: str) -> pd.DataFrame:
    """Load the donor's explicit typing-excluded keys, without treating all Singles as excluded.

    An assignment's explicit mandatory_analysis_exclusion field is authoritative membership,
    including keys its current QC mask already excludes. Historical assignments without that
    field supply QC-eligible typing exclusions; their original QC failures are handled by the
    source mask. A keyed exclusion file is also accepted. Missing donor partitions stop the read.
    """
    path = Path(path)
    if path.is_dir():
        spellings = {str(donor), donor_key(donor)}
        if donor_key(donor) == "112":
            spellings.add("0112")
        matches = [path / f"donor_id={d}/data.parquet" for d in sorted(spellings)
                   if (path / f"donor_id={d}/data.parquet").is_file()]
        if len(matches) != 1:
            raise ValueError(f"{donor}: expected one mandatory exclusion partition under {path}")
        path = matches[0]
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix == ".parquet":
        parquet = pq.ParquetFile(path)
        columns = [c for c in ("donor_id", "object_id", "cell_id", "analysis_eligible",
                               "qc_analysis_eligible", "analysis_exclusion_reason",
                               "mandatory_analysis_exclusion")
                   if c in parquet.schema.names]
        frame = parquet.read(columns=columns).to_pandas()
    else:
        frame = pd.read_csv(path, dtype={"donor_id": str, "object_id": str, "cell_id": str})
    if "donor_id" not in frame:
        frame["donor_id"] = donor_key(donor)
    frame = frame.loc[frame.donor_id.map(donor_key).eq(donor_key(donor))].copy()
    if "mandatory_analysis_exclusion" in frame:
        frame = frame.loc[frame.mandatory_analysis_exclusion.fillna(False).astype(bool)]
    else:
        if "analysis_eligible" in frame:
            frame = frame.loc[~frame.analysis_eligible.fillna(False).astype(bool)]
        if "qc_analysis_eligible" in frame:
            frame = frame.loc[frame.qc_analysis_eligible.fillna(False).astype(bool)]
    key = next((k for k in ("object_id", "cell_id") if k in frame), None)
    if key is None or frame[key].isna().any() or frame.duplicated(["donor_id", key]).any():
        raise ValueError("mandatory protein exclusions require unique nonmissing donor/cell keys")
    return frame.reset_index(drop=True)
