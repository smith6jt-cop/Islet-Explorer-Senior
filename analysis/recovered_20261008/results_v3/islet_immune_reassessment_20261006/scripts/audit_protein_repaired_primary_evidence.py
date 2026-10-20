#!/usr/bin/env python3
"""Report complete hormone evidence, fitted levels and routes without accuracy targets."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ENDO = {"Alpha", "Beta", "Delta", "Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)"}
HORMONES = ("INS", "GCG", "SST")

def read(path, columns=None):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args(); root = args.root.resolve()
    manifest = json.loads((root / "source_manifest.json").read_text())
    donor_rows = []; stratum_rows = []; count_rows = []; sources = []
    for method in manifest["arms"]:
        typing = root / method / "typing"
        audit_path = typing / "subpanel_audit.json"
        audits = {q["donor_id"]: q for q in json.loads(audit_path.read_text())}
        scale_path = root / method / "scale.parquet"
        scale = read(scale_path).set_index("marker")
        for donor in manifest["donors"]:
            path = typing / f"cells/donor_id={donor}/data.parquet"
            columns = ["object_id", "cell_class", "is_endocrine", "endocrine_route", "islet_stratum", "endocrine_object_class", "endocrine_subtype_status", "immune_rules_unavailable_count", "resolved_by", "excluded_reason"]
            calls = read(path, columns)
            aligned = root / method / f"uniform/donor_id={donor}/data.parquet"
            names = set(pq.ParquetFile(aligned).schema.names)
            values = read(aligned, ["object_id", *[m for m in HORMONES if m in names]])
            pair = calls.merge(values, on="object_id", how="left", validate="one_to_one")
            overlap = read(root / f"geometry/crossed/crossed_{donor}.parquet", ["object_id", "crossed_frac_hormone_fragment"])
            if overlap.object_id.duplicated().any():
                raise ValueError(f"Duplicate fragment-overlap UUIDs: {donor}")
            coverage = pair.object_id.map(overlap.set_index("object_id").crossed_frac_hormone_fragment)
            contaminated_or_unmeasured = coverage.isna() | coverage.gt(0)
            complete = np.ones(len(pair), dtype=bool)
            for marker in HORMONES:
                if marker not in pair or marker not in scale.index:
                    complete[:] = False; continue
                floor = float(scale.loc[marker, "floor"]); span = float(scale.loc[marker, "span"])
                # Reproduce declared scaling to account for evidence availability,
                # without interpreting output identities as validation.
                denominator = span if np.isfinite(span) and span > 0 else 1.0
                scaled = np.clip((pair[marker].to_numpy(np.float32) - floor) / denominator, 0, 1)
                scaled[contaminated_or_unmeasured.to_numpy()] = np.nan
                complete &= np.isfinite(scaled)
            a = audits[donor]
            if int((~complete).sum()) != a["n_incomplete_hormone"]:
                raise ValueError(f"Independent finite-three-hormone accounting differs from producer audit: {method} {donor}")
            broad = calls.cell_class.isin(ENDO)
            if not calls.is_endocrine.astype(bool).equals(broad):
                raise ValueError(f"Broad endocrine flag omits/includes unexpected subtype: {method} {donor}")
            single = calls.endocrine_object_class.fillna("").eq("EndocrineSingle")
            ordinary = broad & ~single
            donor_rows.append(dict(method=method, donor_id=donor, n_primary_calls=len(calls),
                n_incomplete_three_hormone=int((~complete).sum()), hormones_acquired=";".join(m for m in HORMONES if m in names),
                n_broad_endocrine=int(broad.sum()), n_ordinary_endocrine=int(ordinary.sum()),
                n_ordinary_unresolved=int((ordinary & calls.cell_class.eq("Endocrine (unresolved subtype)")).sum()),
                n_ordinary_ambiguous=int((ordinary & calls.cell_class.eq("Endocrine (ambiguous subtype)")).sum()),
                n_direct_single=int(single.sum()), n_unknown_immune_rule_evidence=int(calls.immune_rules_unavailable_count.gt(0).sum()),
                n_primary_excluded=int(calls.cell_class.eq("Excluded").sum()), scale_unavailable_markers=";".join(a["scale_unavailable_markers"]),
                n_by_cd99=a["n_by_cd99"], n_by_elongation_only=a["n_by_elongation_only"], n_by_hormone_fallback=a["n_by_hormone_fallback"]))
            for stratum in ["1", "2-3", "4-10", ">10"]:
                hit = calls.islet_stratum.astype(str).eq(stratum).to_numpy()
                level = a["islet_levels"][stratum]
                stratum_rows.append(dict(method=method, donor_id=donor, islet_stratum=stratum,
                    n_cells=int(hit.sum()), n_incomplete_three_hormone=int((hit & ~complete).sum()),
                    n_confirmed_reference=level["n_confirmed"], CD99_level=level["CD99"], nucleus_elongation_level=level["nucleus_elongation"],
                    CD99_level_available=level["CD99"] is not None, shape_level_available=level["nucleus_elongation"] is not None,
                    n_broad_endocrine=int((hit & broad.to_numpy()).sum()),
                    n_ordinary_unresolved=int((hit & ordinary.to_numpy() & calls.cell_class.eq("Endocrine (unresolved subtype)").to_numpy()).sum()),
                    n_ordinary_ambiguous=int((hit & ordinary.to_numpy() & calls.cell_class.eq("Endocrine (ambiguous subtype)").to_numpy()).sum())))
            for field in ["cell_class", "endocrine_route", "endocrine_subtype_status", "resolved_by", "excluded_reason"]:
                for value, count in calls[field].fillna("").astype(str).value_counts(dropna=False).items():
                    count_rows.append(dict(method=method, donor_id=donor, field=field, value=value, n_cells=int(count)))
            sources.append(dict(method=method, donor_id=donor, calls_path=str(path), calls_sha256=sha(path), aligned_path=str(aligned), aligned_sha256=sha(aligned)))
            print("verified repaired primary evidence accounting", method, donor, len(calls), flush=True)
    paths = [root / "repaired_primary_donor_evidence_QC.csv", root / "repaired_primary_stratum_evidence_QC.csv", root / "repaired_primary_routes_statuses_QC.csv"]
    for frame, path in zip([donor_rows, stratum_rows, count_rows], paths):
        pd.DataFrame(frame).to_csv(path, index=False)
    (root / "repaired_primary_evidence_QC.json").write_text(json.dumps(dict(status="passed", sources=sources,
        result_tables=[dict(path=str(p), sha256=sha(p)) for p in paths], accuracy_validated=False,
        interpretation="Evidence completeness, fitted-level availability and operational call/route accounting; no inferred fallback levels or accuracy from output frequencies"), indent=2) + "\n")

if __name__ == "__main__":
    main()
