"""Explicit section identities for the four-slide analysis; no donor-only fallback."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

DEFAULT_MANIFEST = Path(__file__).resolve().parents[2] / "config" / "sameslide_four.json"
REQUIRED_BUNDLE_FILES = (
    "gene_panel.json", "experiment.xenium", "metrics_summary.csv", "cells.parquet",
    "cell_feature_matrix.h5", "cell_boundaries.parquet", "nucleus_boundaries.parquet",
)


def canonical_donor(value):
    value = str(value).strip().upper()
    if value in {'112', '0112', 'HDL112', 'HDL0112'}:
        return '112'
    return str(int(value)) if value.isdecimal() else value


def verify_sex_metadata(manifest):
    import hashlib
    import pandas as pd
    path = Path(manifest['metadata_source'])
    df = pd.read_excel(path, dtype=str)
    df['donor_id'] = df.donor_id.map(canonical_donor)
    receipt = {'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'column': 'sex', 'donors': {}}
    for r in manifest['sections']:
        hits = df[df.donor_id == canonical_donor(r['donor_id'])]
        if len(hits) != 1 or hits.iloc[0]['sex'] != r['sex']:
            raise ValueError(f"{r['section_id']}: recorded sex disagrees with authoritative metadata")
        receipt['donors'][r['donor_id']] = r['sex']
    return receipt


def section_records(manifest: dict, design: str | None = None) -> list[dict]:
    return [r for r in manifest["sections"] if design is None or r["design"] == design]


def get_section(manifest: dict, section_id: str) -> dict:
    hits = [r for r in section_records(manifest) if r["section_id"] == section_id]
    if len(hits) != 1:
        raise ValueError(f"expected one explicitly named section {section_id!r}; found {len(hits)}")
    return hits[0]


def validate_manifest(manifest: dict, check_paths: bool = False) -> list[str]:
    """Return errors; pending protein inputs remain declared and never use legacy exports."""
    errors: list[str] = []
    if manifest.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    seen: set[str] = set()
    pairs: dict[tuple[str, str], set[str]] = {}
    for r in manifest.get("sections", []):
        required = ("section_id", "donor_id", "xenium_serial", "roi", "tissue", "sex",
                    "design", "bundle", "pixel_size_um", "acquisition_run")
        missing = [k for k in required if k not in r]
        if missing:
            errors.append(f"record missing {missing}")
            continue
        sid = r["section_id"]
        if r['donor_id'] != canonical_donor(r['donor_id']):
            errors.append(f'{sid}: donor alias must be canonicalized')
        if sid in seen:
            errors.append(f"duplicate section_id {sid}")
        seen.add(sid)
        if sid != f"{r['donor_id']}/{r['roi']}/{r['xenium_serial']}":
            errors.append(f"{sid}: section_id disagrees with explicit identity")
        if r["design"] not in ("same_slide", "xenium_only"):
            errors.append(f"{sid}: invalid design")
        if r["sex"] not in ("F", "M"):
            errors.append(f"{sid}: sex must be explicitly F or M")
        if (r["roi"], r["tissue"]) not in (("panc", "pancreas"), ("pln_1", "pancreatic_lymph_node")):
            errors.append(f"{sid}: inconsistent tissue and ROI")
        if f"__{r['xenium_serial']}__" not in Path(r["bundle"]).name:
            errors.append(f"{sid}: bundle serial mismatch")
        if not Path(r["bundle"]).is_absolute():
            errors.append(f"{sid}: bundle must be an absolute, explicitly resolved path")
        if float(r["pixel_size_um"]) <= 0:
            errors.append(f"{sid}: pixel_size_um must be positive")
        pairs.setdefault((r["donor_id"], r["roi"]), set()).add(r["design"])
        if r["design"] == "same_slide":
            ph = r.get("pheno") or {}
            if any(not ph.get(k) for k in ("image", "measurements", "geometry", "status")):
                errors.append(f"{sid}: missing explicit protein source")
            if ph.get("status") not in ("ready", "pending_v3_export"):
                errors.append(f"{sid}: unrecognised protein status")
            if ph.get("status") == "pending_v3_export" and ph.get("legacy_fallback"):
                errors.append(f"{sid}: pending v3 source must not fall back to legacy export")
        if check_paths:
            bundle = Path(r["bundle"])
            for name in REQUIRED_BUNDLE_FILES:
                p = bundle / name
                if not p.is_file() or p.stat().st_size == 0:
                    errors.append(f"{sid}: missing or empty {p}")
            ep = bundle / "experiment.xenium"
            if ep.is_file():
                experiment = json.loads(ep.read_text())
                if experiment["slide_id"] != r["xenium_serial"]:
                    errors.append(f"{sid}: experiment slide_id disagrees")
                if abs(experiment["pixel_size"] - r["pixel_size_um"]) > 1e-9:
                    errors.append(f"{sid}: experiment pixel size disagrees")
            if r["design"] == "same_slide":
                ph = r.get("pheno") or {}
                keys = ("image", "measurements", "geometry") if ph.get("status") == "ready" else ("image",)
                for k in keys:
                    if ph.get(k) and not Path(ph[k]).is_file():
                        errors.append(f"{sid}: missing protein {k}: {ph[k]}")
    if len(seen) != manifest.get("expected_sections", len(seen)):
        errors.append("section count differs from expected_sections")
    for key, designs in pairs.items():
        if designs != {"same_slide", "xenium_only"}:
            errors.append(f"{key}: requires both acquisition designs")
    return errors


def load_manifest(path: Path | str | None = None) -> dict:
    manifest = json.loads(Path(path or DEFAULT_MANIFEST).read_text())
    errors = validate_manifest(manifest)
    if errors:
        raise ValueError("invalid same-slide manifest:\n" + "\n".join(errors))
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--check-paths", action="store_true")
    args = ap.parse_args(argv)
    data = load_manifest(args.manifest)
    errors = validate_manifest(data, args.check_paths)
    print(json.dumps({"sections": len(data["sections"]), "errors": errors,
                      "pending_protein": [r["section_id"] for r in section_records(data, "same_slide")
                                          if r["pheno"]["status"] != "ready"]}, indent=2))
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())
