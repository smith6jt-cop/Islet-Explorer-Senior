#!/usr/bin/env python3
"""Fixed, exact-key measured-marker and centroid evidence for corrected protein calls.

This is a diagnostic of operational labels, not independent biological ground truth.
It neither changes labels nor fits thresholds, embeddings or representations.
"""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "6"
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

DEFAULT_ROOT = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008")
MARKERS = ("Caveolin", "CD34", "CD31", "Ker8_18", "E_cadherin", "Pan_Cytokeratin", "EpCAM",
           "INS", "GCG", "SST", "CD99", "Podoplanin", "B3TUBB", "CD56", "CD3e", "CD68",
           "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11c", "CD209", "CD4", "CD8", "FOXP3")
HORMONES = ("INS", "GCG", "SST")
ROI_DONORS = ("0112", "115", "6442", "6493")
POSITIVE = 0.10
PROGRAMS = {
    "vascular_Caveolin_AND_CD34": (("Caveolin",), ("CD34",)),
    "vascular_CD31_AND_Caveolin_or_CD34": (("CD31",), ("Caveolin", "CD34")),
    "epithelial_Ker8_18_AND_junction_or_cytokeratin": (("Ker8_18",), ("E_cadherin", "Pan_Cytokeratin")),
    "epithelial_EpCAM_AND_lineage_marker": (("EpCAM",), ("Ker8_18", "E_cadherin", "Pan_Cytokeratin")),
    "hormone_INS_AND_GCG": (("INS",), ("GCG",)),
    "hormone_INS_AND_SST": (("INS",), ("SST",)),
    "hormone_GCG_AND_SST": (("GCG",), ("SST",)),
    "T_cell_CD3e_AND_CD4_or_CD8": (("CD3e",), ("CD4", "CD8")),
    "B_cell_CD20_AND_CD79a": (("CD20",), ("CD79a",)),
    "myeloid_CD68_AND_Iba1": (("CD68",), ("Iba1",)),
}

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def read(path, columns):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def numeric_islets(frame):
    values = pd.to_numeric(frame.typing_islet_num, errors="coerce")
    counts = values.dropna().value_counts()
    # Deterministic annotation-number selection, with no marker intensity criterion.
    return sorted(int(v) for v, n in counts.items() if v > 0 and v == int(v) and n >= 20)[:2]

def program_state(frame, groups):
    """Three-state AND of OR groups; a finite positive witnesses an OR group.

    A missing alternative never invalidates a positive witness. A group is a known
    negative only if every alternative is observed and below the operational bar.
    """
    positive = np.ones(len(frame), dtype=bool)
    known_negative = np.zeros(len(frame), dtype=bool)
    for group in groups:
        values = frame[list(group)].to_numpy(float)
        finite = np.isfinite(values)
        witness = (finite & (values >= POSITIVE)).any(axis=1)
        negative = finite.all(axis=1) & ~witness
        positive &= witness
        known_negative |= negative
    known = positive | known_negative
    return positive, known

def summarize(frame, donor, label_field):
    marker_rows, program_rows = [], []
    for label, group in frame.groupby(label_field, sort=True):
        for marker in MARKERS:
            values = group[marker].to_numpy(float)
            finite = np.isfinite(values)
            observed = values[finite]
            marker_rows.append(dict(donor_id=donor, label_field=label_field, label=label, marker=marker,
                n_sample=len(group), n_observed=int(finite.sum()), n_unavailable=int((~finite).sum()),
                median=float(np.median(observed)) if len(observed) else None,
                p10=float(np.percentile(observed, 10)) if len(observed) else None,
                p90=float(np.percentile(observed, 90)) if len(observed) else None,
                n_operational_positive=int((observed >= POSITIVE).sum()),
                operational_positive_fraction=float((observed >= POSITIVE).mean()) if len(observed) else None))
        for name, groups in PROGRAMS.items():
            positive, known = program_state(group, groups)
            program_rows.append(dict(donor_id=donor, label_field=label_field, label=label, program=name,
                n_sample=len(group), n_decidable=int(known.sum()), n_unknown=int((~known).sum()),
                n_coexpressing=int(positive.sum()),
                coexpression_fraction=float(positive.sum() / known.sum()) if known.any() else None))
    return marker_rows, program_rows

def heatmap(summary, value, columns, out, title):
    # One equal donor/class sample per donor: each row is the median of donor-level
    # summaries, not a cell-frequency weighted claim about the complete cohort.
    matrix = summary.groupby(["label", columns])[value].median().unstack(columns)
    matrix = matrix.reindex(columns=list(MARKERS) if columns == "marker" else list(PROGRAMS))
    data = np.ma.masked_invalid(matrix.to_numpy(float))
    fig, ax = plt.subplots(figsize=(max(11, len(matrix.columns) * .46), max(5, len(matrix) * .43)))
    cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("#d5d5d5")
    plot = ax.imshow(data, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(matrix.index)), matrix.index, fontsize=9)
    ax.set_title(title, fontsize=11)
    fig.colorbar(plot, ax=ax, label=value.replace("_", " "), shrink=.8)
    fig.text(.01, .01, "Fixed UUID samples; gray = unavailable. Illustrative bar 0.10; fitted structural/hormone gates are not replayed.", fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(out.with_suffix(".png"), dpi=160)
    fig.savefig(out.with_suffix(".pdf")); plt.close(fig)

def roi_plot(frame, donor, roi, bounds, label_field, out):
    labels = sorted(frame[label_field].unique())
    cmap = plt.get_cmap("tab20")
    colors = {label: cmap(i % 20) for i, label in enumerate(labels)}
    fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharex=True, sharey=True)
    axes[0].scatter(frame.X_centroid, frame.Y_centroid, c=[colors[q] for q in frame[label_field]], s=5, linewidths=0)
    axes[0].set_title(label_field)
    panels = (("Caveolin / CD34", ("Caveolin", "CD34")), ("CD31", ("CD31",)),
              ("Ker8_18 / E_cadherin / Pan_Cytokeratin", ("Ker8_18", "E_cadherin", "Pan_Cytokeratin")),
              ("INS / GCG / SST", HORMONES))
    for ax, (name, markers) in zip(axes[1:], panels):
        values = frame[list(markers)].to_numpy(float)
        any_finite = np.isfinite(values).any(axis=1)
        # No zero substitution: unavailable points are displayed separately.
        intensity = np.full(len(frame), np.nan)
        intensity[any_finite] = np.nanmax(values[any_finite], axis=1)
        ax.scatter(frame.X_centroid[~any_finite], frame.Y_centroid[~any_finite], color="#c9c9c9", s=5, linewidths=0)
        plot = ax.scatter(frame.X_centroid[any_finite], frame.Y_centroid[any_finite], c=intensity[any_finite],
                          cmap="magma", vmin=0, vmax=1, s=5, linewidths=0)
        ax.set_title(name); fig.colorbar(plot, ax=ax, shrink=.62, label="declared scaled measurement")
    for ax in axes:
        ax.set_xlim(bounds[0], bounds[1]); ax.set_ylim(bounds[3], bounds[2]); ax.set_aspect("equal")
        ax.set_xlabel("exported X_centroid")
    axes[0].set_ylabel("exported Y_centroid")
    handles = [Line2D([0], [0], marker="o", ls="", color=colors[label], label=label, markersize=4) for label in labels]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .055), ncol=min(5, len(handles)), fontsize=8)
    fig.suptitle(f"Donor {donor}; fixed annotation islet {roi}; {len(frame):,} eligible cell centroids", fontsize=11)
    fig.text(.01, .015, "Centroid measurement maps; not IF image overlays. CD31 participates in refinement. Gray = unavailable.", fontsize=8)
    fig.tight_layout(rect=(0, .21, 1, .95))
    fig.savefig(out.with_suffix(".png"), dpi=160); fig.savefig(out.with_suffix(".pdf")); plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--assignments", type=Path, help="Fresh current-generation assignment root; exact UUID QC/admission intersection")
    ap.add_argument("--label-field", choices=("cell_class", "refined_class"), default="cell_class")
    ap.add_argument("--sample-per-class-donor", type=int, default=512)
    ap.add_argument("--roi-halfwidth", type=float, default=150, help="Halfwidth in exported centroid coordinate units")
    args = ap.parse_args(); root = args.root.resolve()
    if not 1 <= args.sample_per_class_donor <= 512 or args.roi_halfwidth <= 0:
        raise ValueError("Bounded diagnostic requires 1..512 UUIDs per donor/class and positive ROI halfwidth")
    typing = root / "own_compartment/typing"; run_path = typing / "run.json"
    run = json.loads(run_path.read_text())
    if run["status"] != "complete":
        raise ValueError("Corrected primary, fresh refinement and immune subtyping must be terminal")
    config_path = root / "generation_config.json"
    config = json.loads(config_path.read_text())
    if len(run["donors"]) != 25 or set(run["donors"]) != set(config["donors"]) or run["generation_config_sha256"] != sha(config_path):
        raise ValueError("Typing donor universe/configuration differs from frozen corrected all25 generation")
    assignment_binding = None
    assignment_records = {}
    if args.assignments:
        assignment_manifest = args.assignments / "manifest.json"
        assignment = json.loads(assignment_manifest.read_text())
        if assignment.get("state") != "complete" or assignment.get("policy") != "current_generation":
            raise ValueError("Program evidence requires a terminal fresh current_generation assignment view")
        assignment_records = {str(int(q["donor_id"])): q for q in assignment["donors"]}
        if len(assignment_records) != 25 or set(assignment_records) != {str(int(d)) for d in run["donors"]}:
            raise ValueError("Current assignment universe differs from exact corrected all25 typing donors")
        release_roots = {str(Path(q["native_source_path"]).parent.parent.resolve()) for q in assignment_records.values()}
        if len(release_roots) != 1: raise ValueError("Assignments mix native release generations")
        release_manifest = Path(next(iter(release_roots))) / "manifest.json"
        release = json.loads(release_manifest.read_text())
        if Path(release["calibration_run"]).resolve() != run_path.resolve() or release["calibration_run_sha256"] != sha(run_path):
            raise ValueError("Fresh assignment native release is not bound to this terminal corrected typing")
        assignment_binding = dict(manifest=str(assignment_manifest.resolve()), sha256=sha(assignment_manifest),
            release_manifest=str(release_manifest), release_manifest_sha256=sha(release_manifest))
    out = (args.out or root / ("program_ROI_evidence_" + args.label_field)).resolve()
    if out.exists(): raise FileExistsError("Preserve prior evidence; use a fresh output directory")
    out.mkdir(parents=True); shutil.copy2(__file__, out / "executed_source.py")
    donors = run["donors"]; scale_path = Path(run["scale"])
    scale = read(scale_path, ["marker", "floor", "span"]).set_index("marker")
    source_rows = []; sample_frames = []; marker_rows = []; program_rows = []; census = []; roi_rows = []
    for donor in donors:
        calls_path = typing / f"cells/donor_id={donor}/data.parquet"
        columns = ["object_id", "cell_class", "refined_class", "qc_analysis_eligible", "X_centroid", "Y_centroid", "typing_islet_num"]
        calls = read(calls_path, columns)
        if calls.object_id.isna().any() or calls.object_id.duplicated().any(): raise ValueError(f"Invalid typing UUIDs: {donor}")
        qc = calls.qc_analysis_eligible.fillna(False).astype(bool)
        eligible = qc & calls.cell_class.ne("Excluded")
        assignment_path = None
        if args.assignments:
            assignment_path = args.assignments / f"donor_id={str(int(donor))}/data.parquet"
            record = assignment_records[str(int(donor))]
            if assignment_path.resolve() != Path(record["path"]).resolve() or sha(assignment_path) != record["sha256"]:
                raise ValueError(f"Fresh current assignment partition checksum/path changed: {donor}")
            assignment = read(assignment_path, ["object_id", "analysis_eligible", "native_cell_class", "native_refined_class"])
            if assignment.object_id.duplicated().any(): raise ValueError(f"Duplicate assignment UUIDs: {donor}")
            keyed = assignment.set_index("object_id")
            admitted = calls.object_id.map(keyed.analysis_eligible)
            if admitted.isna().any(): raise ValueError(f"Current assignment misses a typed UUID: {donor}")
            for field in ("cell_class", "refined_class"):
                if not calls.object_id.map(keyed["native_" + field]).fillna("").eq(calls[field].fillna("")).all():
                    raise ValueError(f"Fresh assignment native labels differ from terminal corrected typing: {donor} {field}")
            eligible &= admitted.astype(bool)
        calls = calls.loc[eligible].copy()
        if calls[args.label_field].isna().any(): raise ValueError(f"Missing terminal label: {donor}")
        for label, n in calls[args.label_field].value_counts().items():
            census.append(dict(donor_id=donor, label_field=args.label_field, label=label, n_eligible=int(n)))
        sample = calls.sort_values("object_id").groupby(args.label_field, sort=True).head(args.sample_per_class_donor).copy()
        sample["sample_selected"] = True
        windows = []
        if donor in ROI_DONORS:
            if not np.isfinite(calls[["X_centroid", "Y_centroid"]].to_numpy(float)).all():
                raise ValueError(f"Fixed ROI donor has unavailable coordinates: {donor}")
            numeric = pd.to_numeric(calls.typing_islet_num, errors="coerce")
            for roi in numeric_islets(calls):
                center = calls.loc[numeric.eq(roi), ["X_centroid", "Y_centroid"]].median()
                half = args.roi_halfwidth
                bounds = [float(center.X_centroid-half), float(center.X_centroid+half), float(center.Y_centroid-half), float(center.Y_centroid+half)]
                window = calls.loc[calls.X_centroid.between(*bounds[:2]) & calls.Y_centroid.between(*bounds[2:])].copy()
                if len(window): windows.append((roi, bounds, window))
                roi_rows.append(dict(donor_id=donor, islet_num=roi, bounds=bounds, n_cells=len(window),
                    status="plotted" if len(window) else "unavailable_empty_fixed_window",
                    n_annotation_members=int(numeric.eq(roi).sum()), centroid_center=[float(center.X_centroid), float(center.Y_centroid)]))
        wanted = pd.concat([sample, *[q[2] for q in windows]], ignore_index=True).drop_duplicates("object_id")
        aligned_path = Path(run["aligned_root"]) / f"donor_id={donor}/data.parquet"
        names = set(pq.ParquetFile(aligned_path).schema.names)
        measured = read(aligned_path, ["object_id", *[m for m in MARKERS if m in names]])
        measured = measured.loc[measured.object_id.isin(wanted.object_id)].copy()
        if measured.object_id.duplicated().any(): raise ValueError(f"Duplicate aligned UUIDs: {donor}")
        wanted = wanted.merge(measured, on="object_id", how="left", validate="one_to_one", indicator=True)
        if not wanted._merge.eq("both").all(): raise ValueError(f"Unmatched measurement UUIDs: {donor}")
        wanted.drop(columns="_merge", inplace=True)
        fragment_path = Path(run["geometry_root"]) / f"crossed/crossed_{donor}.parquet"
        fragments = read(fragment_path, ["object_id", "crossed_frac_hormone_fragment"])
        if fragments.object_id.duplicated().any(): raise ValueError(f"Duplicate fragment UUIDs: {donor}")
        coverage = wanted.object_id.map(fragments.set_index("object_id").crossed_frac_hormone_fragment)
        wanted["hormone_fragment_coverage"] = coverage
        for marker in MARKERS:
            if marker not in wanted or marker not in scale.index:
                wanted[marker] = np.nan; continue
            floor = float(scale.loc[marker, "floor"]); span = float(scale.loc[marker, "span"])
            raw = wanted[marker].to_numpy(np.float32)
            scaled = np.full(len(raw), np.nan)
            if np.isfinite(floor) and np.isfinite(span) and span > 0:
                finite = np.isfinite(raw)
                scaled[finite] = np.clip((raw[finite]-floor)/span, 0, 1)
            wanted[marker] = scaled
        wanted.loc[coverage.isna() | coverage.gt(0), list(HORMONES)] = np.nan
        wanted["donor_id"] = donor
        sampled = wanted.loc[wanted.object_id.isin(sample.object_id)].copy()
        m, p = summarize(sampled, donor, args.label_field); marker_rows.extend(m); program_rows.extend(p)
        sample_frames.append(sampled)
        for roi, bounds, window in windows:
            selected = wanted.loc[wanted.object_id.isin(window.object_id)].copy()
            selected.to_parquet(out / f"ROI_{donor}_{roi}_exact_cells.parquet", index=False)
            roi_plot(selected, donor, roi, bounds, args.label_field, out / f"ROI_{donor}_{roi}")
        source_rows.append(dict(donor_id=donor, calls=str(calls_path), calls_sha256=sha(calls_path),
            aligned=str(aligned_path), aligned_sha256=sha(aligned_path), fragment_overlap=str(fragment_path), fragment_overlap_sha256=sha(fragment_path),
            assignment=str(assignment_path) if assignment_path else None, assignment_sha256=sha(assignment_path) if assignment_path else None,
            n_sample=len(sampled), n_eligible=len(calls),
            n_sample_fragment_coverage_unavailable=int(sampled.hormone_fragment_coverage.isna().sum()),
            n_sample_fragment_contaminated=int(sampled.hormone_fragment_coverage.gt(0).sum()),
            n_eligible_coordinates_unavailable=int((~np.isfinite(calls[["X_centroid", "Y_centroid"]].to_numpy(float)).all(axis=1)).sum())))
        print("fixed measured evidence", donor, len(sampled), len(windows), flush=True)
    samples = pd.concat(sample_frames, ignore_index=True)
    samples.to_parquet(out / "sample_exact_cells.parquet", index=False)
    marker = pd.DataFrame(marker_rows); program = pd.DataFrame(program_rows)
    marker.to_csv(out / "marker_evidence_by_donor_class.csv", index=False)
    program.to_csv(out / "coexpression_by_donor_class.csv", index=False)
    pd.DataFrame(census).to_csv(out / "eligible_label_census.csv", index=False)
    heatmap(marker, "operational_positive_fraction", "marker", out / "marker_evidence_heatmap", "Median donor-level operational-positive fraction by reported class")
    heatmap(program, "coexpression_fraction", "program", out / "coexpression_heatmap", "Median donor-level decidable coexpression fraction by reported class")
    files = [p for p in sorted(out.iterdir()) if p.is_file()]
    manifest = dict(status="complete", created_at_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        typing_run=str(run_path), typing_run_sha256=sha(run_path), source_sha256=sha(__file__), sources=source_rows,
        generation_config=str(config_path), generation_config_sha256=sha(config_path), assignment_binding=assignment_binding,
        scale=str(scale_path), scale_sha256=sha(scale_path), donors=donors, label_field=args.label_field,
        sample_selection="First lexicographically sorted UUIDs per donor/reported class, independent of marker intensity",
        sample_per_class_donor=args.sample_per_class_donor, n_sample_exact_UUIDs=len(samples),
        eligibility="Typed QC-eligible, non-Excluded primary cells" + (" intersected with fresh current-generation analysis_eligible by exact UUID" if args.assignments else ""),
        ROI_selection="First two numeric typing_islet_num IDs with >=20 eligible annotation members in fixed donors 0112/115/6442/6493",
        ROI_halfwidth_exported_coordinate_units=args.roi_halfwidth, ROI_windows=roi_rows,
        measurement_policy="Frozen own-compartment aligned floor/span; nonfinite remains unavailable; hormone fragment masking by exact UUID",
        CD8_units="Primary shifted own-compartment units; refinement itself used separate inverse-converted unshifted units",
        operational_positive_bar=POSITIVE, positive_bar_interpretation="Fixed illustrative/immune-rule diagnostic bar; does not reproduce fitted vascular, epithelial or hormone decision thresholds", programs=PROGRAMS,
        accuracy_validated=False, labels_modified=False, embedding_fitted=False,
        limitations=["Markers overlap the typing rules and provide operational evidence, not held-out biological ground truth",
            "CD31 is separate from the primary Caveolin/CD34 object gate but participates in refinement",
            "ROI maps display measured cell centroids rather than raw IF pixels or image overlays",
            "Bounded class-stratified samples do not estimate full-cohort cell frequencies",
            "Missing negative evidence remains unknown; no unobserved marker is set to zero"],
        outputs=[dict(path=str(p), sha256=sha(p)) for p in files])
    (out / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("FIXED PROGRAM AND CENTROID EVIDENCE COMPLETE; BIOLOGICAL VERDICT NOT INFERRED", flush=True)

if __name__ == "__main__":
    main()
