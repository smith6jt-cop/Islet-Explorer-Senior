"""Count the immune-root marker evidence among pLN cells unresolved after broad review."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

IMMUNE_ROOT = ["CD68", "CD11c", "CD20", "CD3e", "CD11b"]
CYTOKERATIN = "Pan-Cytokeratin"
COMPLEMENTARY = ["CD79a", "CD163", "Iba1", "MPO", "CD66"]
ADDITIONAL = ["CD4", "CD8", "CD56", "CD206", "CD14", "CD209", "CD57"]


def _audit_donor(job):
    donor, review_file, marker_calls = job
    review = pd.read_parquet(review_file, filters=[("donor", "=", donor), ("candidate_state", "=", "no_match")],
                             columns=["object_id", "pln_lineage_review"])
    selected = pd.Index(review.loc[review.pln_lineage_review.isna(), "object_id"])
    if selected.hasnans or not selected.is_unique:
        raise ValueError(f"{donor}: residual object IDs are not complete and unique")
    markers = IMMUNE_ROOT + [CYTOKERATIN] + COMPLEMENTARY + ADDITIONAL
    native = pq.ParquetFile(marker_calls / f"{donor}__pln_1.parquet").read(
        columns=["object_id", *["marker_call:" + marker for marker in markers]]).to_pandas().set_index("object_id")
    if not native.index.is_unique or not selected.isin(native.index).all():
        raise ValueError(f"{donor}: residual cells disagree with the native marker-call index")
    calls = native.loc[selected]
    calls.columns = markers
    if not calls.isin([0, 1]).all().all():
        raise ValueError(f"{donor}: a marker is unassessed, so a complete-case census cannot be reported")
    positive = calls.eq(1)
    anchor = positive[IMMUNE_ROOT].any(axis=1)
    panck = positive[CYTOKERATIN]
    complementary = positive[COMPLEMENTARY].any(axis=1)
    if (anchor & ~panck).any():
        raise ValueError(f"{donor}: some unresolved cells satisfy the complete immune root")
    counts = {
        "immune_root_anchor_vetoed_panck": anchor & panck,
        "no_immune_root_anchor_panck_positive": ~anchor & panck,
        "no_immune_root_anchor_panck_negative": ~anchor & ~panck,
        "complementary_any": complementary,
        "complementary_without_root_anchor": complementary & ~anchor,
        "root_or_complementary": anchor | complementary,
        "root_or_complementary_or_additional": anchor | complementary | positive[ADDITIONAL].any(axis=1),
        "t_cell_pair": positive.CD3e & (positive.CD4 | positive.CD8),
        "b_cell_pair": positive.CD20 & positive.CD79a,
        "myeloid_pair": positive[["CD68", "Iba1", "CD163", "CD14"]].sum(axis=1).ge(2),
        "granulocyte_pair": positive.MPO & positive.CD66,
    }
    counts["any_paired_lineage_signature"] = (counts["t_cell_pair"] | counts["b_cell_pair"]
        | counts["myeloid_pair"] | counts["granulocyte_pair"])
    counts.update({f"positive_{marker}": positive[marker] for marker in markers})
    result = {"donor": donor, "residual_cells": len(selected),
              **{name: int(hit.sum()) for name, hit in counts.items()}}
    if sum(result[key] for key in ["immune_root_anchor_vetoed_panck", "no_immune_root_anchor_panck_positive",
                                   "no_immune_root_anchor_panck_negative"]) != len(selected):
        raise ValueError("Immune-root outcomes do not partition the residual cells")
    return result


def build_residual_marker_support(*, review_cells, marker_calls, donor_reference, output, workers=3):
    """Return one cell census per donor after verifying each single-marker marginal."""
    review_cells, marker_calls, donor_reference, output = map(Path, (review_cells, marker_calls, donor_reference, output))
    # Donor reference is the reviewed full-section marker output produced independently
    # from the original aligned intensities and scored using dedicated pLN thresholds.
    baseline = pd.read_parquet(donor_reference)
    baseline = baseline[baseline.review_status.ne("Cell and named cluster agree")]
    donors = sorted(baseline.donor.astype(str).unique())
    with ProcessPoolExecutor(max_workers=workers) as executor:
        rows = list(executor.map(_audit_donor, [(donor, review_cells, marker_calls) for donor in donors]))
    result = pd.DataFrame(rows).set_index("donor")
    for marker in IMMUNE_ROOT + [CYTOKERATIN] + COMPLEMENTARY + ADDITIONAL:
        acquired = marker.replace("-", "_")
        expected = baseline[baseline.marker.eq(acquired)].groupby("donor", observed=True)
        np.testing.assert_array_equal(result[f"positive_{marker}"], expected.n_positive.sum().reindex(result.index))
        np.testing.assert_array_equal(result.residual_cells, expected.n_measured.sum().reindex(result.index))
    disease = baseline[["donor", "disease"]].drop_duplicates().set_index("donor")
    result = result.join(disease, validate="one_to_one").reset_index()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".residual-markers-", dir=output.parent) as work:
        staged = Path(work)/output.name
        result.to_parquet(staged, index=False, compression="zstd")
        staged.replace(output)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["review-cells", "marker-calls", "donor-reference", "output"]:
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    table = build_residual_marker_support(**vars(parser.parse_args()))
    print(table[["residual_cells", "immune_root_anchor_vetoed_panck", "complementary_without_root_anchor",
                 "any_paired_lineage_signature"]].sum().to_string())
