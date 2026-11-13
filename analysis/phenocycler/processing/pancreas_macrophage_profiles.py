"""Rebuild the primary- and final-parent macrophage sidecars from aligned native protein values."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MARKERS = ("CD68", "CD163", "M2Gal3")
VERSION = "pancreas-macrophage-gal3-two-classes-v2"
PROFILES = ("Macrophage_Gal3+", "Macrophage_Gal3-")
ANCHOR_COLUMNS = {"CD68": "macrophage_cd68_state", "CD163": "macrophage_cd163_state"}
RULE = "CD68 or CD163 at or above the 0.10 scaled gate establishes anchor positivity. A measured galectin-3 (M2Gal3/LGALS3) result assigns Gal3+ at or above the same gate, otherwise Gal3-."


def pool_previous(values):
    """Use the measured Gal3 sign of legacy labels while preserving unknown labels and missingness."""
    return values.astype("string").str.replace(
        r"^CD68[+-]CD163[+-]Macrophage_LGALS3([+-])(?:CD4[+-])?$",
        r"Macrophage_Gal3\1",
        regex=True,
    )


def marker_states(aligned, scale, threshold):
    values = aligned.reindex(columns=MARKERS).to_numpy(dtype=np.float32)
    parameters = scale.reindex(MARKERS)
    floor = parameters["floor"].to_numpy(float)
    span = parameters["span"].to_numpy(float)
    calibrated = np.isfinite(floor) & np.isfinite(span) & (span > 0)
    scaled = np.clip((values - floor) / np.where(calibrated, span, np.nan), 0.0, 1.0)
    measured = np.isfinite(values) & np.isfinite(scaled)
    positive = measured & (scaled >= threshold)
    called = positive[:, :2].any(axis=1) & measured[:, 2]
    result = pd.Series(
        pd.NA, index=aligned.index, dtype="string", name="macrophage_profile"
    )
    result.loc[called] = np.where(positive[called, 2], *PROFILES)
    states = pd.DataFrame(index=aligned.index)
    for position, column in enumerate(ANCHOR_COLUMNS.values()):
        states[column] = pd.Categorical.from_codes(
            np.where(measured[:, position], positive[:, position].astype(int), 2),
            categories=["below_gate", "positive", "unmeasured"],
        )
    return result, states


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def atomic(path, writer):
    path.parent.mkdir(exist_ok=True, parents=True)
    with NamedTemporaryFile(dir=path.parent, prefix=".profiles-", delete=False) as tmp:
        pending = Path(tmp.name)
    try:
        writer(pending)
        pending.replace(path)
    finally:
        pending.unlink(missing_ok=True)


def process(donor, *, aligned, stage, output, scale, gate):
    part = Path(f"donor_id={donor}") / "data.parquet"
    native = pq.ParquetFile(aligned / part)
    fields = [m for m in MARKERS if m in native.schema.names]
    raw = (
        native.read(columns=["object_id", *fields], use_threads=False)
        .to_pandas()
        .set_index("object_id")
    )
    if not raw.index.is_unique:
        raise ValueError(f"{donor}: duplicate aligned UUID")
    summaries = []
    for lane in ["immune_primary", "immune_final"]:
        target = stage / lane / part
        frame = pq.ParquetFile(target).read().to_pandas()
        classified = frame.immune_class.eq("Imm-macrophage").fillna(False)
        ids = pd.Index(frame.loc[classified, "object_id"])
        if not ids.is_unique or not ids.isin(raw.index).all():
            raise ValueError(
                f"{donor}: macrophage UUID outside the native measurements"
            )
        names, anchors = marker_states(raw.loc[ids], scale, gate[lane])
        recorded = (
            frame.loc[classified, "macrophage_profile"].astype("string").set_axis(ids)
        )
        prev = pool_previous(recorded)
        old_valid = prev.notna()
        if not prev[old_valid].eq(names[old_valid]).fillna(False).all():
            raise ValueError(
                f"{donor} {lane}: measured Gal3 sign disagrees with the previous call"
            )
        if not frame.loc[~classified, "macrophage_profile"].isna().all():
            raise ValueError(f"{donor}: profile outside the macrophage parent")
        frame["macrophage_profile"] = pd.Series(
            pd.NA, index=frame.index, dtype="string"
        )
        frame.loc[classified, "macrophage_profile"] = names.to_numpy()
        frame["macrophage_profile_status"] = "not_macrophage"
        frame.loc[classified, "macrophage_profile_status"] = np.where(
            names.notna(), "measured", "unavailable_marker_measurement"
        )
        for col in anchors:
            frame[col] = anchors[col].reindex(frame.object_id).array
        dest = output / lane / part
        atomic(dest, lambda p: frame.to_parquet(p, index=False, compression="zstd"))
        summaries.append(
            dict(
                donor=donor,
                lane=lane,
                parents=len(names),
                unavailable=int(names.isna().sum()),
                counts=names.value_counts().to_dict(),
                newly_measured=int((~old_valid & names.notna()).sum()),
                old_source_sha256=sha(target),
                output_sha256=sha(dest),
                source_bytes=target.stat().st_size,
            )
        )
    return summaries


def rebuild(stage, aligned, output, workers):
    if stage.resolve() == output.resolve():
        raise ValueError("Choose a separate staging directory")
    scale = pd.read_parquet(stage.parent / "pancreas_scale.parquet").set_index("marker")
    lanes = ["immune_primary", "immune_final"]
    covers = {
        lane: json.loads((stage / lane / "coverage.json").read_text()) for lane in lanes
    }
    donors = sorted(covers["immune_final"]["donors_built"])
    if set(donors) != set(covers["immune_primary"]["donors_built"]):
        raise ValueError("Parent-stage donor rosters differ")
    run = partial(
        process,
        aligned=aligned,
        stage=stage,
        output=output,
        scale=scale,
        gate={k: float(v["require_min"]) for k, v in covers.items()},
    )
    with ProcessPoolExecutor(workers) as pool:
        reports = [v for pair in pool.map(run, donors) for v in pair]
    for lane in lanes:
        coverage = covers[lane]
        counts = pd.read_csv(
            stage / lane / "per_donor_counts.csv", dtype={"donor_id": str}
        )
        rows = []
        for item in (r for r in reports if r["lane"] == lane):
            donor = item["donor"]
            by_donor = coverage["per_donor"][donor]
            by_donor["macrophage_profile"] = item["counts"]
            by_donor["macrophage_profile_unavailable"] = item["unavailable"]
            complete = item["parents"] > 0 and item["unavailable"] == 0
            for label in PROFILES:
                rows.append(
                    dict(
                        donor_id=donor,
                        level="macrophage_profile",
                        label=label,
                        count=item["counts"].get(label, 0) if complete else np.nan,
                        status="measured" if complete else "EMPTY_profile_unavailable",
                        reason=""
                        if complete
                        else "macrophage parent, a positive CD68/CD163 anchor, or measured M2Gal3 unavailable",
                        count_without_majority_term=np.nan,
                        supporting_only_note="",
                    )
                )
        result = pd.concat(
            [counts.loc[counts.level.ne("macrophage_profile")], pd.DataFrame(rows)],
            ignore_index=True,
        )
        atomic(
            output / lane / "per_donor_counts.csv",
            lambda p: result.to_csv(p, index=False),
        )
        atomic(
            output / lane / "coverage.json",
            lambda p: p.write_text(json.dumps(coverage, indent=2) + "\n"),
        )
    report = {
        "version": VERSION,
        "inputs": list(MARKERS),
        "rule": RULE,
        "profiles": list(PROFILES),
        "anchor_state_fields": ANCHOR_COLUMNS,
        "gate": {k: float(v["require_min"]) for k, v in covers.items()},
        "sections": reports,
    }
    atomic(
        output / "macrophage_profile_method.json",
        lambda p: p.write_text(json.dumps(report, indent=2) + "\n"),
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("stage", "aligned", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    rebuild(args.stage, args.aligned, args.output, args.workers)
