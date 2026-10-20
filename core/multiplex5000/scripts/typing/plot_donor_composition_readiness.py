#!/usr/bin/env python3
"""Plot donor-level object composition and feature availability for a composition build."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


GROUPS = ("endocrine", "vascular", "lymphatic", "neural", "stromal", "immune")
OBJECT_CLASSES = ("EndocrineSingle", "EndocrineCluster", "Islet")
CLASS_COLORS = {"EndocrineSingle": "#7aa6c2", "EndocrineCluster": "#e6ae61",
                "Islet": "#6f9d77"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composition", type=Path, required=True,
                        help="composition output directory")
    parser.add_argument("--out", type=Path, required=True,
                        help="figure basename; writes PNG and SVG")
    parser.add_argument("--trajectory-manifest", type=Path, default=None,
                        help="optional ND trajectory manifest for the feature-policy caption")
    args = parser.parse_args()
    manifest_path = args.composition / "manifest.json"
    objects_path = args.composition / "endocrine_objects.parquet"
    manifest = json.loads(manifest_path.read_text())
    if not objects_path.is_file():
        raise FileNotFoundError(f"cohort object table is missing: {objects_path}")
    objects = pd.read_parquet(objects_path)
    if objects.duplicated(["donor_id", "endocrine_object_id"]).any():
        raise ValueError("composition table repeats donor/object keys")
    donors = sorted(objects.donor_id.astype(str).unique())
    coverage = {str(r["donor_id"]): r for r in manifest.get("per_donor", [])}
    available = np.array([
        [bool(coverage.get(d, {}).get("available_groups", {}).get(g, False)) for g in GROUPS]
        for d in donors
    ], dtype=bool)
    counts = (objects.groupby([objects.donor_id.astype(str), "original_object_class"])
              .size().unstack(fill_value=0).reindex(index=donors, columns=OBJECT_CLASSES,
                                                    fill_value=0))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax_heat, ax_count) = plt.subplots(
        2, 1, figsize=(12.5, 8.8), sharex=True,
        gridspec_kw={"height_ratios": [1.1, 2.2], "hspace": 0.08})
    cmap = ListedColormap(["#e9e9e9", "#287a68"])
    ax_heat.imshow(available.astype(int).T, aspect="auto", interpolation="nearest", cmap=cmap,
                   vmin=0, vmax=1)
    ax_heat.set_yticks(np.arange(len(GROUPS)), labels=[g.title() for g in GROUPS])
    ax_heat.set_xticks(np.arange(len(donors)), labels=donors)
    ax_heat.tick_params(axis="x", labelbottom=False, length=0)
    ax_heat.tick_params(axis="y", length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_title("Measured feature groups by donor", loc="left", fontsize=11, pad=7)

    bottom = np.zeros(len(donors), dtype=float)
    for class_name in OBJECT_CLASSES:
        vals = counts[class_name].to_numpy(dtype=float)
        ax_count.bar(np.arange(len(donors)), vals, bottom=bottom,
                     color=CLASS_COLORS[class_name], width=0.82, label=class_name)
        bottom += vals
    ax_count.set_ylabel("Endocrine objects")
    ax_count.set_xticks(np.arange(len(donors)), labels=donors, rotation=65, ha="right")
    ax_count.tick_params(axis="y", length=0)
    for side in ("top", "right"):
        ax_count.spines[side].set_visible(False)
    ax_count.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_count.set_axisbelow(True)
    ax_count.legend(frameon=False, ncols=3, loc="upper left", bbox_to_anchor=(0, -0.24))

    status = str(manifest.get("analysis_status", "unspecified"))
    fig.suptitle(f"Endocrine-object composition readiness ({status})", x=0.08,
                 y=0.985, ha="left", fontsize=14, fontweight="bold")
    trajectory = (json.loads(args.trajectory_manifest.read_text())
                  if args.trajectory_manifest and args.trajectory_manifest.is_file() else {})
    unavailable = trajectory.get("unavailable_feature_groups", [])
    used = trajectory.get("feature_groups", [])
    if used:
        detail = (f"ND trajectory uses {len(used)} measured groups ({', '.join(used)}); "
                  f"unavailable and excluded from fitting: {', '.join(unavailable) or 'none'}.")
    else:
        detail = "Missing groups remain unavailable; no marker-intensity proxy is recoded as cell identity."
    fig.text(0.08, 0.945,
             "Grey = unavailable. " + detail,
             ha="left", va="top", fontsize=9, color="#555555")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.89, bottom=0.25)
    fig.savefig(args.out.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(args.out.with_suffix(".svg"), facecolor="white")
    plt.close(fig)
    print(f"wrote donor composition figure for {len(donors)} donors -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
