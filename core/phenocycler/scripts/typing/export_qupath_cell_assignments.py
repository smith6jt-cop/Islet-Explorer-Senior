#!/usr/bin/env python3
"""Export QuPath-object-keyed native pancreas assignments, one CSV per donor/image.

The export keeps primary and refined identities separate and retains state/reason/evidence fields.
It is keyed to the exact source QuPath object UUID and image name; it does not rewrite a QuPath
project or collapse unresolved calls into a class.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


REPO = Path(__file__).resolve().parents[4]
EXPORT_FIELDS = {
    "object_id": "object_id",
    "image": "image",
    "cell_id": "cell_id",
    "donor_id": "donor_id",
    "primary_identity": "cell_class",
    "Classification": "cell_class",
    "refined_identity": "refined_class",
    "Refined classification": "refined_class",
    "immune_class": "immune_class",
    "immune_subtype": "immune_subtype",
    "primary_call_state": "primary_call_state",
    "primary_call_unavailability_reason": "primary_call_unavailability_reason",
    "refinement_state": "refinement_state",
    "immune_subtyping_state": "immune_subtyping_state",
    "top_identity_marker": "top_identity_marker",
    "top_expressed_marker": "top_expressed_marker",
    "argmax_scaled_no_supporting_only": "argmax_scaled_no_supporting_only",
    "immune_definition": "immune_definition",
    "resolved_by": "resolved_by",
    "refine_pass": "refine_pass",
    "refine_note": "refine_note",
    "concordance": "concordance",
    "concordance_before": "concordance_before",
    "immune_class_matches": "immune_class_matches",
    "immune_class_also": "immune_class_also",
    "immune_class_undecided": "immune_class_undecided",
    "immune_subtype_matches": "immune_subtype_matches",
    "immune_subtype_also": "immune_subtype_also",
    "immune_subtype_undecided": "immune_subtype_undecided",
    "qc_analysis_eligible": "qc_analysis_eligible",
    "analysis_reason": "analysis_reason",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def image_token(image: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", image).strip("._-")[:48] or "image"
    suffix = hashlib.sha256(image.encode("utf-8")).hexdigest()[:8]
    return f"{stem}_{suffix}"


def export_donor(path: Path, donor: str, out: Path) -> dict:
    source_names = set(pq.ParquetFile(path).schema.names)
    missing = sorted(set(EXPORT_FIELDS.values()) - source_names)
    if missing:
        raise ValueError(f"{donor}: candidate cell release missing QuPath export fields {missing}")
    source_columns = list(dict.fromkeys(EXPORT_FIELDS.values()))
    source = pd.read_parquet(path, columns=source_columns)
    frame = pd.DataFrame({name: source[column] for name, column in EXPORT_FIELDS.items()})
    frame["object_id"] = frame["object_id"].astype("string")
    frame["image"] = frame["image"].astype("string")
    frame["cell_id"] = frame["cell_id"].astype("string")
    if frame["object_id"].isna().any() or frame["object_id"].duplicated().any():
        raise ValueError(f"{donor}: object_id is missing or duplicated")
    if frame["image"].isna().any():
        raise ValueError(f"{donor}: a QuPath object has no source image label")
    if not (frame["object_id"] == frame["cell_id"]).all():
        raise ValueError(f"{donor}: release cell_id no longer matches QuPath object_id")

    donor_dir = out / f"donor_id={donor}"
    donor_dir.mkdir(parents=True, exist_ok=True)
    image_exports = []
    for image, subset in frame.groupby("image", sort=True, dropna=False):
        image = str(image)
        destination = donor_dir / f"cell_assignments_{image_token(image)}.csv"
        if destination.exists():
            raise FileExistsError(f"refusing to overwrite QuPath assignment export: {destination}")
        subset.to_csv(destination, index=False, na_rep="", lineterminator="\n")
        image_exports.append({"image": image, "rows": len(subset),
                              "path": str(destination.resolve()), "bytes": destination.stat().st_size,
                              "sha256": sha256(destination)})
    return {"donor_id": donor, "source_rows": len(frame),
            "images": image_exports, "source_path": str(path.resolve()),
            "source_sha256": sha256(path)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--donors", nargs="*", default=None,
                    help="optional validation subset; omitted means all release donors")
    args = ap.parse_args()
    source_manifest_path = args.release_root / "manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text())
    donors = [str(x["donor_id"]) for x in source_manifest.get("donors", [])]
    selected = donors if not args.donors else [d for d in donors if d in set(map(str, args.donors))]
    if args.donors and set(map(str, args.donors)) != set(selected):
        raise ValueError("requested QuPath export donors are not all present in the source release")
    if not selected:
        raise ValueError("no donors selected for QuPath assignment export")
    args.out.mkdir(parents=True, exist_ok=True)
    donor_outputs = []
    for donor in selected:
        path = args.release_root / f"donor_id={donor}" / "data.parquet"
        if not path.is_file():
            raise FileNotFoundError(path)
        donor_outputs.append(export_donor(path, donor, args.out))
        print(f"{donor}: exported {donor_outputs[-1]['source_rows']:,} QuPath-keyed assignments",
              flush=True)
    payload = {
        "title": "QuPath-keyed pancreas cell assignments",
        "status": "candidate_validation_build" if source_manifest.get("status") != "production" else "production",
        "source_release": str(args.release_root.resolve()),
        "source_release_manifest_sha256": sha256(source_manifest_path),
        "export_builder": str(Path(__file__).resolve()),
        "export_builder_sha256": sha256(Path(__file__).resolve()),
        "key": ["donor_id", "image", "object_id"],
        "qupath_join": "object_id is the exact source QuPath UUID; image is retained to prevent cross-image collisions",
        "identity_fields": "primary_identity and refined_identity are distinct; blank calls remain blank and state/reason fields explain availability",
        "confidence_fields": "no calibrated class probability is invented; native concordance, argmax-marker, match, undecided, and refinement evidence are exported as-is",
        "columns": list(EXPORT_FIELDS),
        "number_of_donors": len(donor_outputs),
        "rows": sum(x["source_rows"] for x in donor_outputs),
        "donors": donor_outputs,
    }
    (args.out / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {payload['rows']:,} rows across {len(donor_outputs)} donors -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
