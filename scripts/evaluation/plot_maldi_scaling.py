#!/usr/bin/env python3
"""Show linear/log display scaling for every declared MALDI section; no normalization."""
import argparse
import html
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from multiplex5000 import maldi
from multiplex5000.parameter_observations import file_hash

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    sources = maldi.source_files(args.source)
    bindings = [{"path": str(f.path), "sha256": file_hash(f.path)} for f in sources]
    assert all(b["sha256"] == f.sha256 for b, f in zip(bindings, sources, strict=True))
    roster = pd.read_csv(args.inventory, keep_default_na=False, dtype={"donor_id": str})
    con = maldi.connect()
    features = maldi.features(con, sources)
    paired = maldi.pairing_problems(con, sources)
    assert not any(v for k, v in paired.items() if not k.endswith("_rows"))
    data = con.execute(maldi.spots_sql(sources, features, roster.tissue_id.tolist())).fetchdf()
    rows = []
    for section in roster.to_dict("records"):
        tid = section["tissue_id"]
        frame = data.loc[data.tissue_id.eq(tid)]
        expected = section.get("n_spots", section.get("expected_spots"))
        if expected is None:
            raise ValueError("inventory must provide an explicit expected spot count")
        assert len(frame) == expected == frame.spot_id.nunique()
        values = frame[[name for _, name in features]].to_numpy(float)
        assert np.isfinite(values).all() and (values >= 0).all()
        summed = values.sum(axis=1)
        on = maldi.on_tissue(summed)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout="constrained")
        for ax, color, title in zip(axes[:2], [summed, np.log1p(summed)],
                                     ["Linear display", "log1p display"], strict=True):
            points = ax.scatter(frame.x_um, frame.y_um, c=color, s=1.2, cmap="magma", vmin=0)
            bar = fig.colorbar(points, ax=ax, shrink=0.7)
            bar.set_label("Summed exported intensity" if title.startswith("Linear") else "ln(1 + summed exported intensity)")
            ax.set_title(title)
        axes[2].scatter(frame.x_um[on], frame.y_um[on], s=1.2, color="#777777", label="summed > 0")
        axes[2].scatter(frame.x_um[~on], frame.y_um[~on], s=3, color="#10aebd", label="all-zero, excluded")
        axes[2].set_title(f"Existing signal mask: {int((~on).sum()):,} excluded")
        axes[2].legend(fontsize=8)
        for ax in axes:
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set(xlabel="Instrument x (µm)", ylabel="Instrument y (µm)")
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.ticklabel_format(style="plain", useOffset=False)
            ax.tick_params(labelsize=9)
            if np.ptp(frame.x_um) / max(np.ptp(frame.y_um), 1.0) < 0.65:
                ax.xaxis.set_major_locator(MaxNLocator(2))
                ax.tick_params(axis="x", labelrotation=30)
        fig.suptitle(f"{tid} · {len(frame):,} spots · 50 µm sampling\n"
                     "Display transforms only; 48-feature sum, original coordinates, no registration", fontsize=12)
        path = args.out / f"{tid}.png"
        fig.savefig(path, dpi=120)
        plt.close(fig)
        rows.append({"acquisition_id": section["acquisition_id"], "tissue_id": tid, "spots": len(frame),
                     "mask_excluded": int((~on).sum()), "summed_intensity_min": float(summed.min()),
                     "summed_intensity_median": float(np.median(summed)), "summed_intensity_max": float(summed.max()),
                     "plot": path.name, "sha256": file_hash(path)})
    report = {"schema": "maldi-display-scaling/1", "scope": "display only; no scientific normalization or physical registration",
              "source_bindings": bindings, "inventory_sha256": file_hash(args.inventory),
              "script_sha256": file_hash(Path(__file__)), "intensity_units": "supplier-exported; calibration unavailable",
              "coordinate_units": "um", "sections": rows, "spots": len(data)}
    (args.out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    cards = "\n".join(f'<section><h2>{html.escape(r["tissue_id"])}</h2><a href="{r["plot"]}"><img src="{r["plot"]}" alt="Linear, log1p and mask views for {html.escape(r["tissue_id"])}"></a></section>' for r in rows)
    page = '<!doctype html><meta charset="utf-8"><title>MALDI scaling plots</title><style>body{font:16px system-ui;margin:2rem;max-width:1600px}img{width:100%;height:auto}section{margin:3rem 0}p{max-width:75ch}</style><h1>MALDI display scaling: all 32 sections</h1><p>Linear and log1p views use the same exported 48-feature intensity sum and original 50 µm raster coordinates. Colors are scaled within each section, with explicit color bars. No scientific normalization, tissue identity, registration or neighbour assignment is introduced. Cyan marks measured all-zero spots excluded by the existing signal mask; this mask is not an independently validated tissue outline.</p>' + cards
    (args.out / "index.html").write_text(page + "\n")
    print(json.dumps({"sections": len(rows), "spots": len(data), "index": str(args.out / "index.html")}))


if __name__ == "__main__":
    main()
