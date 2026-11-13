"""Sensitivity to measurement compartment in dedicated pLN, using the released 0.10 gate refit per section and arm."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
from dataclasses import replace
from functools import partial
from tempfile import NamedTemporaryFile
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from multiplex5000.tissues.evaluate import evaluate_rules_batch
from multiplex5000.protein_evidence import candidate_calls
from multiplex5000.tissues.lymph_node import load_rule_source

REFERENCE_CANDIDATES = {
    "CD3e": [
        "T_cell",
        "CD8_CTL",
        "Treg",
        "Tfh",
        "CD4_T_helper",
        "CD8_T_cytotoxic",
        "Effector_CD8",
        "Memory_CD4",
        "Memory_CD8",
    ],
    "CD20": ["B_cell", "Memory_B_cell", "Naive_B_cell", "GC_B_cell"],
    "CD68": ["Macrophage", "M1_macrophage", "M2_macrophage"],
    "Pan-Cytokeratin": [],
    "CD31": ["Endothelial", "Blood_vessel", "HEV", "Lymphatic_vessel"],
    "B3TUBB": ["Neural", "Nerve_fiber", "Schwann_cell"],
    "SMA": ["Smooth_muscle"],
}


def write_table(frame, path):
    with NamedTemporaryFile(
        dir=path.parent, prefix="." + path.name + ".", delete=False
    ) as tmp:
        pending = Path(tmp.name)
    try:
        frame.to_parquet(pending, index=False, compression="zstd")
        pending.replace(path)
    finally:
        pending.unlink(missing_ok=True)


def run_section(
    section, *, exports, stage, calls_dir, output, rule_source, compartments
):
    pa.set_cpu_count(1)
    donor, entry, source_file = section
    path = exports / Path(source_file).name
    root_rules = {r.phenotype: r.constraints for r in rule_source.rules}
    root_markers = sorted(compartments)
    features = ["Object ID", "Nucleus: Area µm^2"]
    features += ["Cell: " + m + ": Mean" for m in root_markers]
    features += [compartments[m] + ": " + m + ": Mean" for m in root_markers]
    raw = (
        pq.ParquetFile(path)
        .read(columns=features, use_threads=False)
        .to_pandas()
        .set_index("Object ID")
    )
    calls = (
        pq.ParquetFile(stage / "protein" / f"{entry}.parquet")
        .read(
            columns=["object_id", "protein_eligible"]
            + ["marker_call:" + m for m in root_markers],
            use_threads=False,
        )
        .to_pandas()
        .set_index("object_id")
    )
    status = pd.read_parquet(
        calls_dir / f"{donor}.parquet",
        columns=["analysis_eligible", "candidate_state", "candidate_phenotype"],
    )
    if (
        not status.index.is_unique
        or not raw.index.equals(calls.index)
        or not raw.index.equals(status.index)
    ):
        raise ValueError(
            f"{entry}: native object IDs differ between source measurements and calls"
        )
    eligible = status.analysis_eligible.to_numpy()
    unmatched = eligible & status.candidate_state.eq("no_match").to_numpy()
    if not np.array_equal(
        eligible, calls.protein_eligible.to_numpy()
    ) or not np.array_equal(eligible, raw["Nucleus: Area µm^2"].gt(0).to_numpy()):
        raise ValueError(f"{entry}: nucleus-admitted cells disagree")
    thresholds = (
        pd.read_csv(stage / "protein_thresholds.csv")
        .query("protein_entry == @entry")
        .set_index("marker")
    )
    states = {v: {} for v in ["released", "cell_raw", "own_raw"]}
    marker_rows = []
    for m in root_markers:
        t = thresholds.loc[m]
        arm = {}
        for name, column, scale in [
            ("released", "Cell", t.shift_scale),
            ("cell_raw", "Cell", 1.0),
            ("own_raw", compartments[m], 1.0),
        ]:
            a = candidate_calls(
                entry,
                m,
                raw[f"{column}: {m}: Mean"].to_numpy(),
                eligible,
                shift_scale=scale,
            )
            states[name][m] = a.states
            arm[name] = a
        if not np.array_equal(
            arm["released"].states, calls["marker_call:" + m].to_numpy()
        ) or not np.isclose(arm["released"].cutoff, t.cutoff):
            raise ValueError(
                f"{entry}: refitted {m} whole-cell calls differ from the release"
            )
        leaf_ref = status.candidate_state.eq("unique_candidate").to_numpy() & eligible
        for pop, mask in [
            ("all_admitted", eligible),
            ("source_unmatched", unmatched),
            (
                "source_marker_associated_candidates",
                leaf_ref
                & status.candidate_phenotype.isin(
                    REFERENCE_CANDIDATES.get(m, [])
                ).to_numpy(),
            ),
        ]:
            key = dict(
                donor=donor,
                marker=m,
                compartment=compartments[m],
                population=pop,
                cells=int(mask.sum()),
                reference_candidates=";".join(REFERENCE_CANDIDATES.get(m, []))
                if pop.startswith("source_marker")
                else "",
            )
            for name, ant in arm.items():
                key.update(
                    {
                        name + "_positive": int(np.sum((ant.states == 1) & mask)),
                        name + "_measured": int(np.sum((ant.states >= 0) & mask)),
                        name + "_cutoff": ant.cutoff,
                        name + "_floor": ant.floor,
                        name + "_top": ant.top,
                        name + "_shift": ant.shift_scale,
                        name + "_fit_background_n": ant.n_background_support,
                    }
                )
            both = mask & (arm["cell_raw"].states >= 0) & (arm["own_raw"].states >= 0)
            key.update(
                both_measured=int(both.sum()),
                gained=int(
                    np.sum(
                        both
                        & (arm["own_raw"].states == 1)
                        & (arm["cell_raw"].states == 0)
                    )
                ),
                lost=int(
                    np.sum(
                        both
                        & (arm["own_raw"].states == 0)
                        & (arm["cell_raw"].states == 1)
                    )
                ),
            )
            marker_rows.append(key)
    section_rows = []
    # Use the same three-state root evaluator as the source release. All required
    # positive and negative measurements, including all alternatives, must exist.
    matches = {}
    for arm_name, state in states.items():
        evaluated = evaluate_rules_batch(rule_source, state)
        matches[arm_name] = {
            name: (values == 1) & eligible for name, values in evaluated.states.items()
        }
    if any((mask & unmatched).any() for mask in matches["released"].values()):
        raise ValueError(f"{entry}: a source no-match cell satisfies a released root")
    native_output = output / "unmatched"
    qualifier = pd.DataFrame({"object_id": raw.index[unmatched], "donor": donor})
    for name, clauses in root_rules.items():
        term = "root_" + name.lower()
        for arm in ("cell_raw", "own_raw"):
            qualifier[term + "_" + arm] = matches[arm][name][unmatched]
        qualifier[term + "_own_all_markers_measured"] = np.logical_and.reduce(
            [states["own_raw"][m] >= 0 for m, op in clauses]
        )[unmatched]
    for arm in ("cell_raw", "own_raw"):
        qualifier[arm + "_roots"] = np.sum(
            [a[unmatched] for a in matches[arm].values()], axis=0
        ).astype("int8")
    write_table(qualifier, native_output / f"{donor}.parquet")
    for roots in matches.values():
        n_roots = np.sum(list(roots.values()), axis=0)
        roots.update(
            any_root=n_roots > 0, two_or_more_roots=n_roots > 1, no_root=n_roots == 0
        )
    for label in matches["released"]:
        z = {arm: roots[label] for arm, roots in matches.items()}
        for pop, mask in [("all_admitted", eligible), ("source_unmatched", unmatched)]:
            section_rows.append(
                dict(
                    donor=donor,
                    root=label,
                    population=pop,
                    cells=int(mask.sum()),
                    **{k: int(np.sum(v & mask)) for k, v in z.items()},
                    gained_vs_raw=int(np.sum(mask & ~z["cell_raw"] & z["own_raw"])),
                    lost_vs_raw=int(np.sum(mask & z["cell_raw"] & ~z["own_raw"])),
                )
            )
    return pd.DataFrame(section_rows), pd.DataFrame(marker_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exports",
        type=Path,
        required=True,
        help="Native pLN QuPath cell parquet directory",
    )
    parser.add_argument(
        "--release",
        type=Path,
        required=True,
        help="Dedicated release containing Stage 03 and 06",
    )
    parser.add_argument(
        "--calls",
        type=Path,
        required=True,
        help="November EpCAM-excluded candidate files, named by donor",
    )
    parser.add_argument(
        "--panel", type=Path, required=True, help="Dedicated pLN tissue rule manifest"
    )
    parser.add_argument(
        "--registry", type=Path, required=True, help="Marker compartment ruling CSV"
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Private per-object and donor summaries directory",
    )
    parser.add_argument(
        "--scanpy",
        type=Path,
        help="Install compact summary tables in this Scanpy bundle",
    )
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    source = load_rule_source(args.panel, dedicated=True)
    selected = tuple(
        source.by_phenotype[name] for name in source.candidate_lineage_roots
    )
    rule_source = replace(source, rules=selected)
    registry = pd.read_csv(args.registry).set_index("marker").measurement_compartment
    if not registry.index.is_unique:
        raise ValueError("Marker compartment registry contains duplicate markers")
    markers = sorted(set().union(*(rule.markers for rule in selected)))
    if "EpCAM" in markers:
        raise ValueError("Dedicated rule source retains EpCAM")
    compartments = {m: registry[source.panel_aliases.get(m, m)] for m in markers}
    roster = pd.read_csv(
        args.release / "03_protein_measurement/sections.csv", dtype={"donor_id": str}
    )
    if (
        roster.empty
        or not roster.donor_id.is_unique
        or not roster.protein_entry.is_unique
        or not roster.tissue.eq("pancreatic_lymph_node").all()
    ):
        raise ValueError("A dedicated pLN release with one entry per donor is required")
    sections = list(
        roster[["donor_id", "protein_entry", "source_path"]].itertuples(
            index=False, name=None
        )
    )
    (args.out / "unmatched").mkdir(exist_ok=True, parents=True)
    process = partial(
        run_section,
        exports=args.exports,
        stage=args.release / "06_typing",
        calls_dir=args.calls,
        output=args.out,
        rule_source=rule_source,
        compartments=compartments,
    )
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        files = list(pool.map(process, sections))
    for kind, parts in zip(("roots", "markers"), zip(*files)):
        result = pd.concat(parts, ignore_index=True)
        filename = f"pln_compartment_{kind}.parquet"
        write_table(result, args.out / filename)
        if args.scanpy is not None:
            args.scanpy.mkdir(exist_ok=True, parents=True)
            write_table(result, args.scanpy / filename)


if __name__ == "__main__":
    main()
