#!/usr/bin/env python3
"""Prepare contact sheets and native-cell examples for manual protein diagnostic review."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from PIL import Image


def ref(path):
    path = path.resolve()
    with path.open("rb") as handle:
        return {"path": str(path), "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


def read(path, columns):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification", required=True, type=Path)
    args = parser.parse_args()
    root = args.verification.resolve()
    report_path = root / "remediation_checks.json"
    report = json.loads(report_path.read_text())
    assert report["status"] == "pass"
    receipt_path = Path(report["receipt"]["path"])
    assert ref(receipt_path)["sha256"] == report["receipt"]["sha256"]
    receipt = json.loads(receipt_path.read_text())
    stage = receipt_path.parent
    sources = {row["name"]: row for row in receipt["inputs"]}
    dest = root / "review_materials"
    dest.mkdir(exist_ok=False)
    (dest / "representatives").mkdir()
    (dest / "contact_sheets").mkdir()
    sections = []
    for row in report["rows"]:
        donor = row["donor_id"]
        mixed = read(stage / f"mixed/donor_id={donor}/data.parquet", [
            "object_id", "cell_uid", "acquisition_id", "segmentation_id", "tissue", "typing_state",
            "cell_class", "refined_class", "refine_state", "refine_note"])
        mixed["review_state"] = mixed.refine_state.fillna(mixed.typing_state)
        selected = mixed.groupby("review_state", sort=True, dropna=False).head(1).copy()
        assert selected.object_id.is_unique and len(selected) == mixed.review_state.nunique()
        aligned = read(Path(sources[f"aligned:{donor}"]["path"]), ["object_id", "X_centroid", "Y_centroid"])
        selected = selected.merge(aligned, on="object_id", validate="one_to_one")
        assert np.isfinite(selected[["X_centroid", "Y_centroid"]]).all().all()
        sidecar = read(stage / f"pancreas/refinement/donor_id={donor}/data.parquet", [
            "object_id", "concordance", "demoted", "refine_pass"])
        selected = selected.merge(sidecar, on="object_id", how="left", validate="one_to_one")
        target = dest / f"representatives/{donor}.json"
        target.write_text(selected.to_json(orient="records", indent=2) + "\n")
        sections.append({"donor_id": donor, "representatives": ref(target),
                         "review_states": selected.review_state.tolist(),
                         "call_plot": ref(root / f"release_audit/overlays/{donor}.png"),
                         "availability_plot": ref(root / f"refinement_states/overlays/{donor}.png"),
                         "coordinate_input": sources[f"aligned:{donor}"],
                         "plot_sample": ref(root / f"refinement_states/samples/{donor}.csv.gz")})
        print(donor, len(selected), "representative states", flush=True)
    sheets = []
    for start in range(0, len(sections), 3):
        group = sections[start:start + 3]
        canvas = Image.new("RGB", (1800, 1200), "white")
        for column, row in enumerate(group):
            for level, key in enumerate(("call_plot", "availability_plot")):
                with Image.open(row[key]["path"]) as plot:
                    assert plot.size == (600, 600)
                    canvas.paste(plot.convert("RGB"), (column * 600, level * 600))
        path = dest / f"contact_sheets/sheet_{start // 3 + 1:02}.png"
        canvas.save(path)
        sheets.append({**ref(path), "donor_ids": [row["donor_id"] for row in group]})
    manifest = {"schema": "protein-review-materials/1", "verification": ref(report_path),
                "producer": ref(Path(__file__)), "sections": sections, "contact_sheets": sheets,
                "sampling": "Plots contain the existing deterministic 12000-cell sample; examples select the first native cell in each review state from the full mixed table. Native coordinates are in the aligned table's declared micrometres.",
                "contact_sheet_layout": "Three sections per sheet; call maps above availability maps; no image resampling.",
                "visual_review": "pending; recorded separately after inspection",
                "scope": "Computational diagnostic review, not raw-marker image review or biological approval."}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
