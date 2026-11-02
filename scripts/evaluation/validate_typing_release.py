#!/usr/bin/env python3
"""Read-only full-cohort receipt/identity audit; writes only to a new evaluation root.

Run on a compute node for large releases. No fitting, gate acceptance, or label changes.
The declared output roster is never filtered to successful sections.
"""
import argparse
import hashlib
import json
import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=str, allow_nan=False) + "\n")


def read(path, wanted):
    names = pq.ParquetFile(path).schema_arrow.names
    return pq.ParquetFile(path).read(columns=[c for c in wanted if c in names]).to_pandas()


def receipt_audit(root, out):
    stage = root / "06_typing"
    receipt = json.loads((stage / "receipt.json").read_text())
    checks = []
    for role, items in (("input", receipt["inputs"]), ("output", receipt["outputs"])):
        for item in items:
            path = Path(item["path"]) if role == "input" else stage / item["path"]
            expected = item.get("sha256")
            actual = sha(path) if path.is_file() and expected else None
            checks.append({"role": role, "path": str(path), "expected": expected,
                           "actual": actual, "bytes": path.stat().st_size if path.is_file() else None,
                           "status": "pass" if actual and actual == expected else
                                     "unavailable" if expected is None else "fail"})
    pd.DataFrame(checks).to_csv(out / "receipt_checks.csv", index=False)
    write(out / "receipt.json", receipt)
    for stage_name in ("00_inventory", "03_protein_measurement", "04_rna_ingest", "05_normalization", "06_typing"):
        for filename in ("receipt.json", "gate/accepted.json", "gate/gate.json"):
            p = root / stage_name / filename
            if p.is_file():
                dest = out / "upstream" / stage_name / filename
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(p, dest)
    return receipt, checks


def spatial(out, name, frame, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Deterministic diagnostic sample only; accounting always uses all cells.
    frame = frame.iloc[np.linspace(0, len(frame) - 1, min(len(frame), 12000), dtype=int)]
    frame.to_csv(out / f"samples/{name}.csv.gz", index=False)
    fig, ax = plt.subplots(figsize=(5, 5), layout="constrained")
    for value, group in frame.groupby(label, dropna=False):
        ax.scatter(group.x_um, group.y_um, s=0.8, label=str(value), rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set(title=name, xlabel="x (um)", ylabel="y (um)")
    ax.legend(fontsize=5, loc="upper left", markerscale=3)
    fig.savefig(out / f"overlays/{name}.png", dpi=120)
    plt.close(fig)


def rna(root, out):
    acquisitions = pd.read_csv(root / "00_inventory/acquisitions.csv", keep_default_na=False)
    roster = acquisitions[acquisitions.platform == "xenium"]
    calls = pd.read_csv(root / "06_typing/rna_acquisitions.csv", keep_default_na=False)
    shutil.copyfile(root / "00_inventory/acquisitions.csv", out / "acquisitions.csv")
    shutil.copyfile(root / "00_inventory/relationships.csv", out / "relationships.csv")
    for name in ("rna_acquisitions.csv", "rna_controls.csv", "rna_classes.csv"):
        shutil.copyfile(root / "06_typing" / name, out / name)
    rows = []
    for entry in roster.to_dict("records"):
        found = calls[calls.acquisition_id == entry["acquisition_id"]]
        base = {"acquisition_id": entry["acquisition_id"], "donor_id": entry["donor_id"],
                "tissue": entry["tissue"], "identity_check": "blocked", "errors": []}
        if len(found) != 1:
            base["errors"].append("missing or duplicate typing acquisition")
            rows.append(base)
            continue
        row = found.iloc[0].to_dict()
        sample = row["sample_key"]
        base.update(row)
        try:
            cells = read(root / f"04_rna_ingest/cells/{sample}.parquet",
                         ["cell_uid", "cell_id", "typing_eligible", "nucleus_count", "x_um", "y_um"])
            p = root / f"06_typing/rna/{sample}/{sample}.calls.parquet"
            if not p.is_file():
                # A damaged/low-signal region is accounted without manufacturing calls.
                if row["status"] != "unavailable":
                    raise ValueError("candidate has no calls")
                labels = cells[["cell_uid"]].copy()
                labels["lineage"] = "Unavailable_region"
                base["identity_check"] = "unavailable_no_calls_expected"
            else:
                labels = read(p, ["cell_uid", "cell_id", "lineage"])
                if labels.cell_uid.isna().any() or labels.cell_uid.duplicated().any():
                    raise ValueError("null/duplicate typed cell_uid")
                if set(labels.cell_uid) != set(cells.cell_uid):
                    raise ValueError("typing cell_uid universe differs from ingest")
                base["identity_check"] = "pass"
            if cells.cell_uid.isna().any() or cells.cell_uid.duplicated().any():
                raise ValueError("null/duplicate ingested cell_uid")
            if len(cells) != int(row["cells"]):
                raise ValueError("segmented count differs from declared count")
            counts = ["called_cells", "broad_cells", "unresolved_cells", "uncalled_cells", "unavailable_cells"]
            if sum(int(row[c]) for c in counts) != len(cells):
                raise ValueError("state counts do not sum to segmented cells")
            joined = cells.merge(labels, on="cell_uid", validate="one_to_one", suffixes=("", "_call"))
            if "cell_id_call" in joined and not (joined.cell_id == joined.cell_id_call).all():
                raise ValueError("native ID changed")
            eligible = joined.typing_eligible.astype(bool)
            base["no_nucleus"] = int((joined.nucleus_count <= 0).sum())
            base["nucleus_admission_errors"] = int((eligible & (joined.nucleus_count <= 0)).sum())
            bad_call = ~joined.lineage.str.startswith("Unavailable_") & ~eligible
            base["ineligible_non_unavailable_calls"] = int(bad_call.sum())
            if bad_call.any() or base["nucleus_admission_errors"]:
                raise ValueError("nucleus-ineligible cell has a call/admission")
            if not np.isfinite(joined[["x_um", "y_um"]]).all().all():
                raise ValueError("nonfinite coordinates")
            report = root / f"06_typing/rna/{sample}/{sample}.report.json"
            if report.is_file():
                data = json.loads(report.read_text())
                # Preserve fitted diagnostics/coverage without copying large cell arrays.
                write(out / f"reports/{sample}.json", data)
            spatial(out, sample, joined[["cell_uid", "x_um", "y_um", "lineage"]], "lineage")
        except Exception as error:  # noqa: BLE001 -- retain each failed acquisition in the ledger
            base["identity_check"] = "fail"
            base["errors"].append(f"{type(error).__name__}: {error}")
        rows.append(base)
        print(entry["acquisition_id"], base["identity_check"], flush=True)
    write(out / "per_acquisition.json", rows)
    return rows


def spleen(root, out, receipt):
    rows = []
    for section in receipt["settings"]["typed"]:
        base = {"section_id": section, "panel": "spleen_pc67", "identity_check": "blocked", "errors": []}
        try:
            p = root / f"06_typing/data/cells/donor_id={section}/data.parquet"
            cells = read(p, ["cell_uid", "object_id", "cell_class", "refined_class", "status"])
            src = next(x for x in receipt["inputs"] if x["name"] == f"expression:{section}")
            expr = read(src["path"], ["object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid"])
            if cells.object_id.isna().any() or cells.object_id.duplicated().any() or cells.cell_uid.duplicated().any():
                raise ValueError("missing/duplicate native or qualified identifier")
            eligible = expr[expr.qc_analysis_eligible.astype(bool)]
            base.update(cells=len(cells), source_cells=len(expr), source_eligible=len(eligible),
                        source_ineligible=len(expr)-len(eligible))
            if set(eligible.object_id) != set(cells.object_id):
                raise ValueError("typed IDs differ from eligible source IDs")
            joined = cells.merge(eligible, on="object_id", validate="one_to_one").rename(
                columns={"X_centroid": "x_um", "Y_centroid": "y_um"})
            spatial(out, section, joined[["cell_uid", "x_um", "y_um", "cell_class"]], "cell_class")
            base["identity_check"] = "pass"
            base["class_counts"] = cells.refined_class.fillna("(none)").value_counts().to_dict()
        except Exception as error:  # noqa: BLE001 -- retain each failed acquisition in the ledger
            base["identity_check"] = "fail"
            base["errors"].append(f"{type(error).__name__}: {error}")
        rows.append(base)
        print(section, base["identity_check"], flush=True)
    for panel, sections in receipt["settings"]["not_typed"].items():
        rows += [{"section_id": s, "panel": panel, "identity_check": "blocked",
                  "reason": receipt["settings"]["not_typed_reason"]} for s in sections]
    for p in (root / "06_typing/tables").glob("*.csv"):
        dest = out / "tables" / p.name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(p, dest)
    write(out / "per_acquisition.json", rows)
    return rows


def mixed(root, out, receipt):
    """Reconcile the isolated legacy full25 mixed pancreas/embedded-pLN route."""
    rows = []
    inventory_source = next(x["path"] for x in receipt["inputs"] if x["name"] == "inventory_acquisitions")
    inventory = pd.read_csv(inventory_source, dtype=str, keep_default_na=False)
    for donor in receipt["settings"]["donors"]:
        native_source = next(x["path"] for x in receipt["inputs"] if x["name"] == f"aligned:{donor}")
        acquisition = inventory[(inventory.donor_id == donor) & (inventory.source == "qupath_manifest") & (inventory.tissue == "pancreas")]
        if len(acquisition) != 1:
            raise ValueError(f"{donor}: inventory acquisition is not unique")
        base = {"donor_id": donor, "acquisition_id": acquisition.iloc[0].acquisition_id,
                "segmentation_id": acquisition.iloc[0].segmentation_version,
                "identity_check": "blocked", "errors": []}
        try:
            cells = read(root / f"06_typing/mixed/donor_id={donor}/data.parquet",
                         ["object_id", "native_cell_id", "cell_uid", "acquisition_id", "segmentation_id", "qc_analysis_eligible",
                          "tissue", "cell_region", "typing_state", "cell_class", "refined_class", "candidate_state"])
            source = read(native_source, ["object_id", "qc_analysis_eligible", "cell_region", "X_centroid", "Y_centroid"])
            if cells.cell_uid.isna().any() or cells.cell_uid.duplicated().any() or cells.object_id.duplicated().any():
                raise ValueError("missing/duplicate typed identifier")
            if set(cells.object_id) != set(source.object_id) or not cells.native_cell_id.equals(cells.object_id):
                raise ValueError("typed native ID universe changed")
            joined = cells.merge(source, on="object_id", validate="one_to_one", suffixes=("", "_source"))
            if not joined.qc_analysis_eligible.equals(joined.qc_analysis_eligible_source):
                raise ValueError("QC eligibility changed")
            if not joined.cell_region.equals(joined.cell_region_source):
                raise ValueError("tissue membership changed")
            if not cells.acquisition_id.eq(base["acquisition_id"]).all() or not cells.segmentation_id.eq(base["segmentation_id"]).all():
                raise ValueError("acquisition/segmentation identity changed")
            eligible = joined.qc_analysis_eligible.astype(bool)
            node = joined.cell_region.eq("pLN")
            if not joined.loc[~eligible, "typing_state"].eq("unavailable_qc").all():
                raise ValueError("QC-ineligible cell received a typing state")
            if not joined.loc[eligible & ~node, "typing_state"].eq("pancreas_hierarchy").all():
                raise ValueError("eligible pancreas universe not covered")
            if joined.loc[node, "refined_class"].notna().any():
                raise ValueError("embedded pLN acquired a final pancreas class")
            joined = joined.rename(columns={"X_centroid": "x_um", "Y_centroid": "y_um"})
            if not np.isfinite(joined[["x_um", "y_um"]]).all().all():
                raise ValueError("nonfinite coordinates")
            joined["display_call"] = joined.refined_class.fillna(joined.typing_state)
            spatial(out, donor, joined[["cell_uid", "x_um", "y_um", "tissue", "typing_state", "display_call"]], "display_call")
            base.update(cells=len(cells), eligible=int(eligible.sum()), ineligible=int((~eligible).sum()),
                        embedded_pln=int(node.sum()), embedded_pln_eligible=int((node & eligible).sum()),
                        typing_states=cells.typing_state.value_counts().to_dict(), identity_check="pass")
        except Exception as error:  # noqa: BLE001 -- retain each failed acquisition in the ledger
            base["identity_check"] = "fail"
            base["errors"].append(f"{type(error).__name__}: {error}")
        rows.append(base)
        print(donor, base["identity_check"], flush=True)
    for name in ("typing_counts.json", "pancreas_scale.parquet", "pancreas/run.json", "pancreas/subpanel_audit.json"):
        src = root / "06_typing" / name
        if src.is_file():
            dest = out / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dest)
    write(out / "per_acquisition.json", rows)
    return rows


def pln(root, out, receipt):
    rows = []
    for row in pd.read_csv(root / "03_protein_measurement/sections.csv", dtype=str).to_dict("records"):
        entry = row["protein_entry"]
        cells = read(root / f"03_protein_measurement/cells/{entry}.parquet",
                     ["cell_uid", "object_id", "protein_eligible", "nucleus_area_um2", "x_um", "y_um"])
        calls = read(root / f"06_typing/protein/{entry}.parquet",
                     ["cell_uid", "object_id", "protein_eligible", "candidate_state"])
        errors = []
        if cells.cell_uid.duplicated().any() or calls.cell_uid.duplicated().any() or set(cells.cell_uid) != set(calls.cell_uid):
            errors.append("duplicate/mismatched cell_uid")
        joined = cells.merge(calls, on="cell_uid", validate="one_to_one", suffixes=("", "_call"))
        if not joined.object_id.equals(joined.object_id_call) or not joined.protein_eligible.equals(joined.protein_eligible_call):
            errors.append("native identity or eligibility changed")
        no_nucleus = joined.nucleus_area_um2.isna() | (joined.nucleus_area_um2 <= 0)
        if (joined.protein_eligible & no_nucleus).any():
            errors.append("nucleus-ineligible cell admitted")
        if not joined.loc[~joined.protein_eligible, "candidate_state"].eq("unavailable").all():
            errors.append("ineligible cell received candidate")
        if len(joined) != int(row["cells"]) or int(joined.protein_eligible.sum()) != int(row["eligible"]):
            errors.append("section counts do not reconcile")
        row.update(cells=len(joined), identity_check="fail" if errors else "pass", errors=errors)
        spatial(out, entry, joined[["cell_uid", "x_um", "y_um", "candidate_state"]], "candidate_state")
        rows.append(row)
    for name in ("protein_thresholds.csv", "protein_rule_counts.csv", "protein_candidate_status.csv"):
        shutil.copyfile(root / "06_typing" / name, out / name)
    write(out / "per_acquisition.json", rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--kind", required=True, choices=["rna", "spleen", "mixed", "pln"])
    args = parser.parse_args()
    if args.out.resolve().is_relative_to(args.root.resolve()):
        parser.error("evaluation output must be outside the scientific release")
    args.out.mkdir(parents=True, exist_ok=False)
    for name in ("samples", "overlays", "reports"):
        (args.out / name).mkdir()
    receipt, checks = receipt_audit(args.root, args.out)
    rows = (rna(args.root, args.out) if args.kind == "rna" else
            spleen(args.root, args.out, receipt) if args.kind == "spleen" else
            mixed(args.root, args.out, receipt) if args.kind == "mixed" else
            pln(args.root, args.out, receipt))
    summary = {"root": str(args.root), "producer": receipt["code"], "python": platform.python_version(),
               "script_sha256": sha(__file__), "input_output_files": len(checks),
               "hash_failures": [r for r in checks if r["status"] == "fail"],
               "unavailable_hashes": [r for r in checks if r["status"] == "unavailable"],
               "acquisitions": len(rows), "identity_failures": [r for r in rows if r["identity_check"] == "fail"],
               "scope": "computational reconciliation, no biological acceptance or production changes"}
    write(args.out / "summary.json", summary)


if __name__ == "__main__":
    main()
