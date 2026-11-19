"""Estimated Xenium islet outlines from explicitly declared RNA seed evidence.

Historical islet-only producer retained for replay. For all-size whole-cell
objects (including endocrine clusters and singles), use endocrine_workflow.
DBSCAN uses the existing structure caller. Convex hulls export its geometric
support as boundaries in native micrometres. These are expression-dependent
estimates, not IF segmentation, cell identities or accepted registrations.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from ..cohort import excluded_on
from .structures import dbscan_labels


@dataclass(frozen=True)
class BoundaryParams:
    eps_um: float = 50.0
    min_samples: int = 10

    def __post_init__(self):
        if not np.isfinite(self.eps_um) or self.eps_um <= 0:
            raise ValueError("eps_um must be finite and positive")
        if not isinstance(self.min_samples, int) or self.min_samples < 3:
            raise ValueError("min_samples must be an integer >=3 for area objects")


OBJECT_COLUMNS = ["structure_id", "acquisition_id", "coordinate_frame", "boundary_source",
                  "boundary_method", "boundary_status", "n_seed_cells", "area_um2",
                  "x_um", "y_um", "max_seed_span_um", "geometry_wkb"]


def cluster_boundaries(cells: pd.DataFrame, *, acquisition_id: str,
                       params: BoundaryParams = BoundaryParams()):
    """Return estimated objects, seed memberships and explicit QC accounting.

    ``seed_evidence`` is supplied by a declared RNA evidence producer. Nucleus
    and positive-depth admission are enforced here even if a caller's mask is
    permissive. Stable native IDs sort DBSCAN traversal and identify seed sets.
    No per-cell labels are assigned to non-seeds by spatial proximity.
    """
    from shapely.geometry import MultiPoint

    if not acquisition_id or not isinstance(acquisition_id, str):
        raise ValueError("An explicit acquisition ID is required")
    required = ["cell_id", "x_um", "y_um", "nucleus_count", "depth", "seed_evidence"]
    frame = cells[required].copy()
    if frame.cell_id.isna().any() or frame.cell_id.astype(str).eq("").any():
        raise ValueError("Missing native cell IDs")
    frame.cell_id = frame.cell_id.astype(str)
    if frame.cell_id.duplicated().any():
        raise ValueError("Native IDs must be unique within an acquisition")
    if not pd.api.types.is_bool_dtype(frame.seed_evidence) or frame.seed_evidence.isna().any():
        raise ValueError("seed_evidence must be an explicit nonmissing Boolean")
    xy = frame[["x_um", "y_um"]].to_numpy(float)
    if not np.isfinite(xy).all():
        raise ValueError("Finite native-micrometre coordinates required")
    admitted = (frame.nucleus_count.gt(0) & frame.depth.gt(0)).fillna(False)
    selected = frame.seed_evidence & admitted
    seeds = frame.loc[selected].sort_values("cell_id").reset_index(drop=True)
    coords = seeds[["x_um", "y_um"]].to_numpy(float)
    labels = (dbscan_labels(coords, params.eps_um, params.min_samples)
              if len(seeds) >= params.min_samples else np.full(len(seeds), -1, dtype=int))
    seeds["structure_id"] = pd.Series([None] * len(seeds), dtype="object")
    seeds["assignment_status"] = "dbscan_noise"
    objects = []
    for k in sorted(set(labels) - {-1}):
        ix = np.flatnonzero(labels == k)
        ids = seeds.cell_id.iloc[ix].tolist()
        digest = hashlib.sha256(json.dumps(ids, separators=(",", ":")).encode()).hexdigest()[:24]
        uid = f"rna_islet:{acquisition_id}:{digest}"
        g = MultiPoint(coords[ix]).convex_hull
        valid = g.geom_type == "Polygon" and g.is_valid and g.area > 0
        center = g.centroid
        span = np.ptp(coords[ix], axis=0)
        objects.append(dict(structure_id=uid, acquisition_id=acquisition_id,
            coordinate_frame=f"xenium:{acquisition_id}:um", boundary_source="rna_dbscan",
            boundary_method="seed_centroid_convex_hull",
            boundary_status="estimated" if valid else "unavailable_degenerate_cluster",
            n_seed_cells=len(ix), area_um2=float(g.area) if valid else np.nan,
            x_um=float(center.x), y_um=float(center.y),
            max_seed_span_um=float(span.max()), geometry_wkb=g.wkb if valid else None))
        seeds.loc[ix, "structure_id"] = uid
        seeds.loc[ix, "assignment_status"] = "cluster_member" if valid else "degenerate_cluster"
    seeds["acquisition_id"] = acquisition_id
    qc = dict(n_cells=len(frame), n_nucleus_positive_depth_cells=int(admitted.sum()),
              n_seed_evidence=int(frame.seed_evidence.sum()),
              n_seed_evidence_excluded_by_nucleus_or_depth=int((frame.seed_evidence & ~admitted).sum()),
              n_seeds=len(seeds), n_noise_seeds=int((labels < 0).sum()),
              n_clusters=len(objects), n_estimated_boundaries=sum(r["boundary_status"] == "estimated" for r in objects),
              largest_cluster_seed_fraction=max((r["n_seed_cells"] / len(seeds) for r in objects), default=0.),
              parameters=asdict(params), individual_islet_accuracy_validated=False)
    return pd.DataFrame(objects, columns=OBJECT_COLUMNS), seeds, qc


def _digest(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_section(spec: dict):
    """Read a pinned RNA-call/vendor-cell pair by exact native ID.

    The review adapter requires broad Endocrine, positive per-cell naming-gene
    support and the producer's lineage eligibility. No generic-label fallback.
    These calls are candidate evidence; the source release status is retained.
    """
    if excluded_on(str(spec["donor_id"]), "xenium"):
        raise ValueError("Section is excluded for the Xenium platform")
    if not re.fullmatch(r"[0-9]{7}", str(spec["source_serial"])):
        raise ValueError("An exact seven-digit acquisition serial is required")
    if spec["coordinate_units"] != "micrometres" or spec["roi"] != "panc":
        raise ValueError("Only explicitly declared pancreas coordinates in micrometres are supported")
    for role in ("cells", "calls"):
        path = Path(spec[f"{role}_path"])
        if _digest(path) != spec[f"{role}_sha256"]:
            raise ValueError(f"{role} SHA-256 mismatch: {path}")
    raw = pd.read_parquet(spec["cells_path"], columns=["cell_id", "x_centroid", "y_centroid", "nucleus_count"])
    wanted = ["cell_id", "source_serial", "donor_id", "roi", "lineage", "lineage_analysis_eligible",
              "n_naming_detected_for_cluster_call", "nucleus_count", "nucleus_measured", "depth"]
    calls = pd.read_parquet(spec["calls_path"], columns=wanted)
    for field, value in (("source_serial", spec["source_serial"]), ("donor_id", spec["donor_id"]), ("roi", spec["roi"])):
        if calls[field].isna().any() or not calls[field].astype(str).eq(str(value)).all():
            raise ValueError(f"Call-table {field} differs from declared section")
    for table in (raw, calls):
        if table.cell_id.isna().any() or table.cell_id.astype(str).duplicated().any():
            raise ValueError("Duplicate or missing native cell key")
        table.cell_id = table.cell_id.astype(str)
    raw_keys, call_keys = set(raw.cell_id), set(calls.cell_id)
    if not call_keys.issubset(raw_keys):
        raise ValueError("Native call keys are absent from the pinned vendor cells")
    if raw_keys != call_keys and spec.get("calls_scope") != "declared_subset_of_vendor_cells":
        raise ValueError("Different key sets require an explicitly declared call-table subset")
    merged = raw.merge(calls, on="cell_id", how="left", validate="one_to_one", indicator=True,
                       suffixes=("_vendor", "_calls"))
    present = merged._merge.eq("both")
    if not merged.loc[present, "nucleus_measured"].fillna(False).all():
        raise ValueError("Nucleus admission evidence unavailable")
    if not np.array_equal(merged.loc[present, "nucleus_count_vendor"].to_numpy(float),
                          merged.loc[present, "nucleus_count_calls"].to_numpy(float), equal_nan=True):
        raise ValueError("Nucleus counts disagree with the pinned vendor cells")
    result = merged.rename(columns={"x_centroid": "x_um", "y_centroid": "y_um", "nucleus_count_vendor": "nucleus_count"})
    result["seed_evidence"] = (result.lineage.eq("Endocrine") & result.lineage_analysis_eligible.fillna(False)
                               & result.n_naming_detected_for_cluster_call.gt(0)).fillna(False).astype(bool)
    result.attrs["input_accounting"] = dict(n_vendor_cells=len(raw), n_cells_with_calls=len(calls),
                                            n_vendor_cells_without_calls=int((~present).sum()))
    return result


def run_manifest(manifest_path: Path, output: Path):
    """Write new, source-labelled RNA boundary candidates; never overwrite a run."""
    import shapely
    from shapely.geometry import mapping

    manifest = json.loads(Path(manifest_path).read_text())
    if manifest.get("schema") != "xenium_islet_boundary_inputs_v1":
        raise ValueError("Unsupported input manifest schema")
    specs = manifest["sections"]
    ids = [r["acquisition_id"] for r in specs]
    if len(ids) != len(set(ids)):
        raise ValueError("Repeated acquisition ID in manifest")
    if output.exists():
        raise FileExistsError(f"Preserve the existing run; select a new output directory: {output}")
    output.mkdir(parents=True)
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    summary = []
    for spec in specs:
        cells = read_section(spec)
        for profile, parameters in manifest["profiles"].items():
            if not profile or Path(profile).name != profile or profile in (".", ".."):
                raise ValueError("Profile must be a simple directory name")
            objects, seeds, qc = cluster_boundaries(cells, acquisition_id=spec["acquisition_id"],
                                                   params=BoundaryParams(**parameters))
            qc.update(cells.attrs["input_accounting"])
            for key in ("donor_id", "source_serial", "roi", "source_release_id", "source_release_status"):
                objects[key] = spec[key]
                seeds[key] = spec[key]
            objects["profile"] = profile
            path = output / profile / spec["source_serial"]
            path.mkdir(parents=True, exist_ok=False)
            objects.to_parquet(path / "objects.parquet", index=False)
            seeds.to_parquet(path / "seed_membership.parquet", index=False)
            features = []
            for row in objects.to_dict("records"):
                wkb = row.pop("geometry_wkb")
                if wkb is not None:
                    features.append(dict(type="Feature", id=row["structure_id"], properties=row,
                                         geometry=mapping(shapely.from_wkb(wkb))))
            (path / "boundaries_um.geojson").write_text(json.dumps(dict(type="FeatureCollection",
                coordinate_frame=f"xenium:{spec['acquisition_id']}:um", coordinate_units="micrometres",
                boundary_source="rna_dbscan", individual_islet_accuracy_validated=False, features=features), allow_nan=False) + "\n")
            (path / "qc.json").write_text(json.dumps(qc, indent=2) + "\n")
            summary.append(dict(donor_id=spec["donor_id"], source_serial=spec["source_serial"],
                acquisition_id=spec["acquisition_id"], profile=profile,
                **{k: v for k, v in qc.items() if k != "parameters"}))
        print(spec["donor_id"], spec["source_serial"], "RNA boundary profiles complete", flush=True)
    pd.DataFrame(summary).to_csv(output / "section_summary.csv", index=False)
    files = [p for p in sorted(output.rglob("*")) if p.is_file()]
    record = dict(status="complete_candidate", manifest_sha256=_digest(Path(manifest_path)),
        module_sha256=_digest(Path(__file__)), sections=len(specs), profiles=list(manifest["profiles"]),
        registration_accepted=False, production_receipt_changed=False,
        artifacts=[dict(path=str(p.relative_to(output)), sha256=_digest(p)) for p in files])
    (output / "run_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return pd.DataFrame(summary)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    run_manifest(args.manifest, args.output)


if __name__ == "__main__":
    main()
