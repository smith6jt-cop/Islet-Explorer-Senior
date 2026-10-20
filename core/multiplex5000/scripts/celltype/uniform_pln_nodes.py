#!/usr/bin/env python3
"""UniFORM alignment of lymph-node tissue: the 19 dedicated pLN sections and the 5 different-block nodes.

Ruled 2026-10-18 ("Nodes aligned among pLN tissue"): a lymph node that shares a slide with a pancreas
but sits in a different block is not part of that pancreas's UniFORM sample (2026-10-02, same-block
basis); it is aligned with lymph-node tissue instead -- the 19 dedicated pLN sections of the IO60pLN
project -- on the markers every one of those panels acquired.

One measurement state for every sample: the raw QuPath means (``Cell: Mean``; ``Nucleus: Mean`` for
the nuclear markers, as the pancreas lane reads them). The dedicated sections have no REDSEA table,
so the nodes' REDSEA values would be a second measurement version in one alignment.

The shift, the refusal rule and the scale are the pancreas lane's (``run_step2.derive``,
``usability.acquisition_group_ratio``, ``run_step3.build_scale``), applied to samples instead of
donors. No reference sample is pinned: none exists for lymph-node tissue, so UniFORM chooses, and
the choice is recorded. These values are NOT in the pancreas alignment's units.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel  # noqa: E402
import score as SC  # noqa: E402
import uniform_adapter as UA  # noqa: E402
import usability as U  # noqa: E402
from run_step2 import BAR_MODE_SD, REFUSED  # noqa: E402

V3 = pathlib.Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
GEN = pathlib.Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008")
DEDICATED = V3 / "core/phenocycler/data/integration/cells_pln"
RAW_CELLS = V3 / "core/phenocycler/data/cells"
#: The stage table that measured each node cell's nucleus (the cell release covers pancreas rows only).
REDSEA_STAGE = V3 / "core/phenocycler/data/redsea"
AVAILABILITY = GEN / "own_compartment/availability"
SIDECARS = pathlib.Path("/data/io60-calibration/alignment_block_current")
#: The five different-block nodes (corrected_block_relationships.json, effective relationship).
NODES = {"0112": "2b123ea5-2aa3-41d0-a44b-2213b89c15e5", "115": "019f5ef5-fb89-4552-b159-e8b438259463",
         "6347": "ed8f8a13-54f3-4bf5-9f79-70cf1207ede1", "6533": "733e19b4-582c-4827-9a29-163d97a00783",
         "6580": "ee6522ed-8576-443a-aaf8-357463cb9ea4"}
#: Acquisition experiment of each sample, for the lane's group-ratio refusal rule.
SEVENTY_CHANNEL = {"0112", "6347", "6533", "6580"}
#: Dedicated-panel spelling -> cohort spelling.
RENAME = {"PDL1": "PD_L1"}
N_PER_SAMPLE = 200_000


def sha(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def measured(names: list[str], marker: str) -> str:
    """The column a marker is read from: the nucleus mean for the nuclear group, else the cell mean."""
    source = {v: k for k, v in RENAME.items()}.get(marker, marker)
    column = f"Nucleus__{source}" if marker in panel.NUCLEUS_MEASURED else source
    return column if column in names else ""


def dedicated_samples() -> dict[str, dict]:
    out = {}
    for path in sorted(DEDICATED.glob("donor_id=*/roi=pln_1/*.parquet")):
        donor = path.parts[-3].split("=", 1)[1]
        out[f"dedicated|{donor}|pln_1"] = {"donor_id": donor, "source": "dedicated_pLN_section",
                                           "annotation_uuid": "", "path": path, "group": "IO60pLN"}
    return out


def node_samples() -> dict[str, dict]:
    out = {}
    for donor, uuid in NODES.items():
        out[f"node|{donor}|{uuid}"] = {
            "donor_id": donor, "source": "different_block_node", "annotation_uuid": uuid,
            "path": RAW_CELLS / f"donor_id={donor}" / "data_0.parquet",
            "group": "pancreas_70_channel" if donor in SEVENTY_CHANNEL else "pancreas_standard"}
    return out


def read_dedicated(spec: dict, markers: list[str]) -> pd.DataFrame:
    f = pq.ParquetFile(spec["path"])
    names = f.schema_arrow.names
    columns = {m: measured(names, m) for m in markers}
    t = f.read(columns=["object_id", "X_centroid", "Y_centroid", "nucleus_area_um2", "qc_analysis_eligible",
                        *[c for c in columns.values() if c]]).to_pandas()
    eligible = t.qc_analysis_eligible.fillna(False).to_numpy(bool) & (pd.to_numeric(t.nucleus_area_um2) > 0).to_numpy()
    out = t[["object_id", "X_centroid", "Y_centroid"]].assign(eligible=eligible)
    for m, c in columns.items():
        out[m] = pd.to_numeric(t[c], errors="coerce").to_numpy(float) if c else np.nan
    return out


def read_node(spec: dict, markers: list[str]) -> pd.DataFrame:
    donor, uuid = spec["donor_id"], spec["annotation_uuid"]
    side = pd.read_parquet(SIDECARS / f"donor_id={donor}.parquet")
    keep = side.loc[side.block.eq("different-block") & side.pln_annotation_id.astype(str).eq(uuid), "object_id"]
    if keep.empty:
        raise ValueError(f"{donor}: no different-block cells for {uuid}")
    f = pq.ParquetFile(spec["path"])
    names = f.schema_arrow.names
    columns = {m: measured(names, m) for m in markers}
    t = f.read(columns=["object_id", "X_centroid", "Y_centroid", *[c for c in columns.values() if c]]).to_pandas()
    t = t[t.object_id.isin(set(keep))]
    qc = pq.read_table(GEN / f"own_compartment/expression/donor_id={donor}/data.parquet",
                       columns=["object_id", "qc_analysis_eligible"]).to_pandas()
    nucleus = pq.read_table(REDSEA_STAGE / f"donor_id={donor}" / "data_0.parquet",
                            columns=["object_id", "nucleus_area_px"]).to_pandas()
    t = t.merge(qc, on="object_id", how="left", validate="one_to_one").merge(
        nucleus, on="object_id", how="left", validate="one_to_one")
    if len(t) != len(keep):
        raise ValueError(f"{donor}: {len(keep)} sidecar cells but {len(t)} raw rows")
    if t.qc_analysis_eligible.isna().any() or t.nucleus_area_px.isna().all():
        raise ValueError(f"{donor}: node cells missing from the QC or nucleus tables")
    eligible = t.qc_analysis_eligible.fillna(False).to_numpy(bool) & (pd.to_numeric(t.nucleus_area_px) > 0).to_numpy()
    out = t[["object_id", "X_centroid", "Y_centroid"]].assign(eligible=eligible)
    for m, c in columns.items():
        out[m] = pd.to_numeric(t[c], errors="coerce").to_numpy(float) if c else np.nan
    return out


def acquired(spec: dict) -> set[str]:
    """Markers a sample's panel acquired: every column of a dedicated section; a node's donor record."""
    if spec["source"] == "dedicated_pLN_section":
        names = pq.ParquetFile(spec["path"]).schema_arrow.names
        return {RENAME.get(n, n) for n in names if not n.startswith("Nucleus__")}
    table = pd.read_parquet(AVAILABILITY / f"donor_id={spec['donor_id']}" / "data.parquet")
    return set(table.loc[table.available.astype(bool), "marker"])


def shared_markers(samples: dict) -> list[str]:
    registry = set(panel.markers(include_variant=True))
    return sorted(set.intersection(*(acquired(spec) & registry for spec in samples.values())))


def align(frames: dict, groups: dict, markers: list[str]) -> dict:
    reg, norm, _ = UA._import_uniform()
    result = {}
    for m in markers:
        vals, unavailable = {}, {}
        for s, f in frames.items():
            v = f.loc[f.eligible, m].to_numpy(float)
            if np.any(np.isfinite(v) & (v >= 1.0)):
                vals[s] = v
            else:
                unavailable[s] = "No finite eligible values at or above UniFORM's log floor (1.0)."
        have = sorted(vals)
        level = pd.Series({s: float(np.nanmedian(v[np.isfinite(v)])) for s, v in vals.items()})
        acq = U.acquisition_group_ratio(level, pd.Series(groups))["acq_ratio"]
        flagged = bool(np.isfinite(acq) and acq >= U.BARS["acq_ratio_flag"])
        if len(have) < 2:
            result[m] = {"reference_sample": have[0] if have else "unavailable", "scale": {s: 1.0 for s in have},
                         "provenance": "too_few_samples_to_align", "unavailable_samples": unavailable,
                         "acq_ratio": acq}
            continue
        counts, ranges = UA.marker_histograms(vals)
        before = UA.background_mode_spread(counts, ranges)
        hist = {m: {"counts": [counts[s] for s in have], "bin_edges": [None] * len(have)}}
        with contextlib.redirect_stdout(io.StringIO()):
            shift_map, refs, _ = reg.automatic_registration(hist, [m], [m], have, None, None, UA.NUM_BINS)
        bins = dict(zip(have, shift_map[m]))
        after = UA.background_mode_spread(counts, ranges, bins)
        log_shift = dict(zip(have, norm.calculate_shift_in_log_pixels({m: ranges}, [m], UA.NUM_BINS, shift_map)[m]))
        scale = {s: float(np.exp(-log_shift[s])) for s in have}
        refused = flagged and after["sd_log"] > BAR_MODE_SD
        if refused:
            scale = {s: 1.0 for s in have}
        over_bar = not refused and after["sd_log"] > BAR_MODE_SD
        result[m] = {"reference_sample": refs[0], "reference_source": "uniform_choice",
                     "shift_bins": {s: int(bins[s]) for s in have}, "shift_log": {s: float(log_shift[s]) for s in have},
                     "scale": scale, "unavailable_samples": unavailable,
                     "provenance": REFUSED if refused else ("automatic_registration_over_bar" if over_bar
                                                            else "automatic_registration"),
                     "over_bar": bool(over_bar), "acq_ratio": acq, "acq_ratio_flagged": flagged,
                     "mode_sd_before": before["sd_log"], "mode_sd_after": after["sd_log"],
                     "mode_fold_before": before["fold"], "mode_fold_after": after["fold"],
                     "global_min": ranges["global_min"], "global_max": ranges["global_max"]}
        tag = "REFUSED" if refused else ("ok" if not over_bar else "OVER BAR")
        print(f"  {m:<16} ref={refs[0]:<28} sd {before['sd_log']:.3f} -> {after['sd_log']:.3f}  {tag}", flush=True)
    return result


def build_scale(frames: dict, markers: list[str]) -> pd.DataFrame:
    """run_step3.build_scale's floor (aligned log1p mode, 512 bins) and top (p99.9), over all samples."""
    rng = np.random.default_rng(0)
    pool = {m: [] for m in markers}
    for s in sorted(frames):
        f = frames[s][frames[s].eligible]
        if len(f) > N_PER_SAMPLE:
            f = f.iloc[rng.choice(len(f), N_PER_SAMPLE, replace=False)]
        for m in markers:
            pool[m].append(f[m].to_numpy("float64"))
    floors = {}
    for m, chunks in pool.items():
        v = np.concatenate(chunks)
        v = v[np.isfinite(v) & (v > 0)]
        if len(v) < 1000:
            continue
        h, e = np.histogram(np.log1p(v), bins=512)
        floors[m] = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
    return SC.marker_scale({m: np.concatenate(c) for m, c in pool.items()}, floors).reset_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise SystemExit(f"{out} exists; choose a new output folder")
    samples = {**dedicated_samples(), **node_samples()}
    if len(samples) != 24:
        raise SystemExit(f"expected 19 dedicated sections and 5 nodes, found {len(samples)} samples")
    markers = shared_markers(samples)
    print(f"{len(samples)} samples, {len(markers)} shared markers", flush=True)
    frames = {}
    for s, spec in samples.items():
        frames[s] = (read_dedicated if spec["source"] == "dedicated_pLN_section" else read_node)(spec, markers)
        print(f"  read {s}: {len(frames[s]):,} cells, {int(frames[s].eligible.sum()):,} eligible", flush=True)
        if not frames[s].eligible.any():
            raise SystemExit(f"{s}: no eligible cell -- refusing to align without it")
    shifts = align(frames, {s: spec["group"] for s, spec in samples.items()}, markers)
    for s, f in frames.items():
        for m in markers:
            f[m] = f[m] * shifts[m]["scale"].get(s, np.nan)
    scale = build_scale(frames, markers)
    (out / "aligned").mkdir(parents=True)
    for s, f in frames.items():
        spec = samples[s]
        f.rename(columns={"eligible": "qc_analysis_eligible_with_nucleus"}).assign(
            sample_id=s, donor_id=spec["donor_id"], source=spec["source"], annotation_uuid=spec["annotation_uuid"]
        ).to_parquet(out / "aligned" / f"{s.replace('|', '__')}.parquet", index=False)
    (out / "shifts.json").write_text(json.dumps(shifts, indent=1, default=float) + "\n")
    scale.to_parquet(out / "scale.parquet", index=False)
    (out / "provenance.json").write_text(json.dumps({
        "created_UTC": datetime.now(timezone.utc).isoformat(), "ruling": "Nodes aligned among pLN tissue (2026-10-18)",
        "uniform_commit": UA.UNIFORM_COMMIT, "markers": markers, "measurement": "raw QuPath Cell: Mean; Nucleus: Mean for "
        + ", ".join(panel.NUCLEUS_MEASURED), "eligibility": "qc_analysis_eligible and nucleus area > 0",
        "reference": "uniform_choice (no pinned lymph-node reference exists)",
        "samples": {s: {**{k: str(v) for k, v in spec.items()}, "sha256": sha(spec["path"])} for s, spec in samples.items()},
        "node_inputs": {d: {"sidecar": str(SIDECARS / f"donor_id={d}.parquet"), "sidecar_sha256": sha(SIDECARS / f"donor_id={d}.parquet")}
                        for d in NODES},
        "code": {"script": str(pathlib.Path(__file__).resolve()), "sha256": sha(__file__)}}, indent=1) + "\n")
    print(f"-> {out}", flush=True)


if __name__ == "__main__":
    main()
