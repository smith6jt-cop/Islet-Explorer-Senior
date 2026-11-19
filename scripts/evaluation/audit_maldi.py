#!/usr/bin/env python3
"""Evaluate all declared MALDI sections without bypassing their inspection gate."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from multiplex5000 import maldi
from multiplex5000.parameter_observations import file_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    files = maldi.source_files(args.source)
    bindings = [{"path": str(f.path), "expected": f.sha256, "actual": file_hash(f.path)} for f in files]
    if any(x["expected"] != x["actual"] for x in bindings):
        raise ValueError("MALDI exports differ from SHA256SUMS")
    roster = pd.read_csv(args.inventory, keep_default_na=False, dtype={"donor_id": str})
    con = maldi.connect()
    features = maldi.features(con, files)
    pairing = maldi.pairing_problems(con, files)
    if any(v for k, v in pairing.items() if not k.endswith("_rows")):
        raise ValueError(f"MALDI pairing failed: {pairing}")
    frame = con.execute(maldi.spots_sql(files, features, roster.tissue_id.tolist())).fetchdf()
    names = [f for _, f in features]
    rows, detection = [], []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    (args.out / "overlays").mkdir()
    (args.out / "samples").mkdir()
    for section in roster.to_dict("records"):
        tid = section["tissue_id"]
        data = frame[frame.tissue_id == tid]
        values = data[names].to_numpy(float)
        summed = np.nansum(values, axis=1)
        on = maldi.on_tissue(summed)
        finite = np.isfinite(values)
        all_zero = finite.all(axis=1) & (values == 0).all(axis=1)
        all_missing = ~finite.any(axis=1)
        x, y = data.x_um.to_numpy(), data.y_um.to_numpy()
        # Report coordinate rounding residue; no new pass tolerance is introduced.
        offsets = np.column_stack((x - x.min(), y - y.min())) / maldi.SPOT_PITCH_UM
        residual = np.abs(offsets - np.rint(offsets)) * maldi.SPOT_PITCH_UM
        row = {"acquisition_id": section["acquisition_id"], "tissue_id": tid,
               "donor_id": section["donor_id"], "spots": len(data), "expected_spots": section["n_spots"],
               "unique_ids": int(data.spot_id.nunique()), "on_tissue": int(on.sum()),
               "all_zero": int(all_zero.sum()), "all_zero_on_tissue": int((all_zero & on).sum()),
               "all_missing": int(all_missing.sum()), "missing_values": int((~finite).sum()),
               "negative_values": int((values < 0).sum()), "grid_residue_max_um": float(residual.max()),
               "duplicate_xy": int(data.duplicated(["x_um", "y_um"]).sum()),
               "registration": "blocked_missing_inspected_tissue_and_neighbour",
               "object_summaries": "blocked_no_accepted_registration",
               "source_units": "exported intensity; calibration/unit declaration unavailable"}
        assert row["spots"] == row["expected_spots"] == row["unique_ids"]
        rows.append(row)
        for i, (role, name) in enumerate(features):
            denominator_all = int(finite[:, i].sum())
            denominator_on = int((finite[:, i] & on).sum())
            positives = int(((values[:, i] > 0) & finite[:, i]).sum())
            positives_on = int(((values[:, i] > 0) & finite[:, i] & on).sum())
            detection.append({"tissue_id": tid, "role": role, "feature": name,
                              "finite_all_spots": denominator_all, "finite_masked_spots": denominator_on,
                              "detected_all": positives, "detected_masked": positives_on,
                              "fraction_all_raster": positives / denominator_all if denominator_all else None,
                              "fraction_masked": positives_on / denominator_on if denominator_on else None})
        fig, ax = plt.subplots(figsize=(5, 4), layout="constrained")
        ax.scatter(x, y, c=np.log1p(summed), s=0.5, cmap="magma")
        ax.scatter(x[~on], y[~on], c="cyan", s=4, label="mask excluded")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set(title=f"{tid}: {len(data):,} spots, {int((~on).sum())} masked out",
               xlabel="instrument x (um)", ylabel="instrument y (um)")
        ax.legend(fontsize=7)
        fig.savefig(args.out / f"overlays/{tid}.png", dpi=120)
        plt.close(fig)
        chosen = np.unique(np.r_[np.linspace(0, len(data)-1, min(12, len(data)), dtype=int),
                                 np.flatnonzero(~on)[:12]])
        data.iloc[chosen].to_csv(args.out / f"samples/{tid}.csv", index=False)
    pd.DataFrame(rows).to_csv(args.out / "sections.csv", index=False)
    pd.DataFrame(detection).to_csv(args.out / "detection_denominators.csv", index=False)
    summary = {"bindings": bindings, "pairing": pairing, "sections": len(rows),
               "donors": int(roster.donor_id.nunique()), "spots": len(frame),
               "all_zero": sum(r["all_zero"] for r in rows),
               "all_zero_on_tissue": sum(r["all_zero_on_tissue"] for r in rows),
               "feature_counts": pd.Series([f[0] for f in features]).value_counts().to_dict(),
               "all_zero_probe": maldi.on_tissue(np.array([0., 0., 2., np.nan])).tolist(),
               "mask_denominator": "finite values at spots with any detected ion; not all raster spots",
               "registration_inspection": "all 32 declarations blank; no tissue inference from Tissue_id",
               "script_sha256": file_hash(Path(__file__))}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
