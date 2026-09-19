#!/usr/bin/env python3
"""Per-islet neighbourhood metrics -- the columns Spatial Cards A, B and C read.

The formulas are ``scripts/compute_neighborhood_metrics.py``'s, unchanged:

    peri_prop_<p>   = peri count of p / total peri cells
    enrich_z_<p>    = (observed - expected) / sqrt(expected),  expected = tissue_prop(p) x peri_total
    min_dist_<p>    = min distance from the islet's core centroid to a peri cell of p

Two choices this cohort forces, both stated rather than silent:

* **The Poisson baseline is the donor's own tissue, not the pooled cohort.**  The legacy script
  pools every donor into one ``tissue_props``, which was defensible for one cohort but here would
  make ``enrich_z`` read out donor composition: a donor whose whole pancreas is immune-rich would
  show enrichment at every islet, including islets that are immune-poor *for that donor*.  The
  question Card B asks -- is this islet's peri zone enriched over the tissue it sits in -- is the
  per-donor one.  The formula is untouched; only the population the expectation comes from is.

* **The named immune sub-metrics are not written.**  ``cd8_to_macro_ratio`` and
  ``tcell_density_peri`` count phenotypes (``CD8a Tcell``, ``CD4 Tcell``, ``Macrophage``) that the
  per-cell table this file reads does not name.  Writing them would put all-NaN columns behind app
  controls that look live.  ``immune_*`` is written, because ``Immune`` is a phenotype here.

  Corrected 2026-09-18: a second pass DOES name them -- ``immune_subtype.py``, over the ``Immune``
  class only, 20 donors -- but it writes to ``typing/immune/`` and **never modifies
  ``cell_class``**, and no app artifact reads that tree.  So the sub-metrics stay unwritten until
  something rules on how a subtype reaches an app artifact; the absence is a decision, not a gap.

A cell with no nucleus is already absent from the per-cell table's eligible rows; ``Excluded`` is
dropped here too so the denominator is the same one every other artifact uses.

    python core/phenocycler/scripts/typing/neighborhood.py
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"

#: The phenotypes the per-cell table calls immune.  One: ``immune_subtype.py``'s finer names live
#: in ``typing/immune/`` and are not joined into any app artifact.
IMMUNE_TYPES = ("Immune",)

# A CLASS WITH NO CELLS IN THIS COHORT, NAMED RATHER THAN DROPPED BY DEFAULT.
#
# ``Unresolved (no lineage evidence)`` is the ninth class.  ``refine.py`` would produce it, it is
# declared in ``run20_default.yml``, and no cell on disk carries it today -- so the vocabulary
# ``main()`` derives below, which is the classes actually OBSERVED in the cohort, has no entry
# for it and this file writes no ``peri_count_Unresolved`` / ``peri_prop_Unresolved`` /
# ``enrich_z_Unresolved`` / ``min_dist_Unresolved``.
#
# That absence is a DECISION, stated here, not pandas' default: the counts below now keep their
# null group (``dropna=False``), so nothing leaves a denominator silently, and the day the
# refinement is wired the class appears in the observed vocabulary and gets its four columns.
# Writing them NOW would add an all-null ``enrich_z_`` and ``min_dist_`` family for a class with
# no cells -- the same all-NaN-column-behind-a-live-looking-control this file already refuses
# for ``cd8_to_macro_ratio`` (see the module docstring).  Which of the two an artifact should
# carry is a question for whoever wires the pass, not a default to be picked here.

#: The app's word for the 20 um zone.
REGION_TO_APP = {"core": "core", "expansion": "peri", "tissue": "tissue"}


def _safe(name: str) -> str:
    """The column-name spelling every consumer of these metrics already assumes."""
    return name.replace(" ", "_").replace("+", "plus")


def _donor_cells(donor: str) -> pd.DataFrame:
    frame = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["phenotype", "cell_class", "cell_region_20um", "islet_key",
                 "X_centroid", "Y_centroid"],
    ).to_pandas()
    frame = frame[frame["cell_class"] != "Excluded"]
    frame["cell_region"] = frame["cell_region_20um"].map(REGION_TO_APP)
    return frame


def donor_metrics(donor: str, phenotypes: list[str]) -> pd.DataFrame:
    cells = _donor_cells(donor)

    # The baseline: this donor's whole tissue, every eligible cell in it.
    # ``dropna=False`` because this is a DENOMINATOR.  A cell whose phenotype did not map is
    # still a cell of this donor's tissue; dropping it would renormalise every other class
    # upward and inflate ``enrich_z`` for all of them.  Its null key is never fetched by the
    # ``.get(pheno, 0.0)`` below, so it raises no class's expectation -- it only keeps the total
    # honest.  A class with no cells in this donor gets ``0.0`` from that ``.get`` and then an
    # UNAVAILABLE z (``expected > 0`` is false), never a zero z: no expectation is not "expected
    # none".
    tissue_props = cells["phenotype"].value_counts(normalize=True, dropna=False).to_dict()

    islet_cells = cells[cells["islet_key"].notna()]
    peri = islet_cells[islet_cells["cell_region"] == "peri"]
    core = islet_cells[islet_cells["cell_region"] == "core"]
    peri_by_islet = dict(tuple(peri.groupby("islet_key")))
    core_by_islet = dict(tuple(core.groupby("islet_key")))

    records = []
    for islet in sorted(set(peri_by_islet) | set(core_by_islet)):
        p = peri_by_islet.get(islet, peri.iloc[:0])
        c = core_by_islet.get(islet, core.iloc[:0])
        # ``dropna=False`` again: ``p_total``/``c_total`` are ``len()``, so a dropped null
        # phenotype would leave the per-class counts no longer summing to the total they are
        # divided by.  Each class is then read with ``.get(pheno, 0)`` -- a class this islet has
        # none of is a printed 0, which is a measurement, and is kept distinct from
        # ``peri_prop_`` being NaN where the islet has no peri cell to be a proportion of.
        p_counts = p["phenotype"].value_counts(dropna=False)
        c_counts = c["phenotype"].value_counts(dropna=False)
        p_total, c_total = len(p), len(c)

        rec: dict[str, object] = {"case_id": int(donor), "islet_key": islet,
                                  "total_cells_peri": p_total, "total_cells_core": c_total}

        for pheno in phenotypes:
            safe = _safe(pheno)
            n = int(p_counts.get(pheno, 0))
            rec[f"peri_count_{safe}"] = n
            rec[f"peri_prop_{safe}"] = n / p_total if p_total else np.nan

            expected = tissue_props.get(pheno, 0.0) * p_total
            rec[f"enrich_z_{safe}"] = ((n - expected) / np.sqrt(expected)
                                       if expected > 0 else np.nan)

        peri_immune = sum(int(p_counts.get(t, 0)) for t in IMMUNE_TYPES)
        core_immune = sum(int(c_counts.get(t, 0)) for t in IMMUNE_TYPES)
        rec["immune_count_peri"] = peri_immune
        rec["immune_count_core"] = core_immune
        rec["immune_frac_peri"] = peri_immune / p_total if p_total else np.nan
        rec["immune_frac_core"] = core_immune / c_total if c_total else np.nan
        rec["immune_ratio"] = ((peri_immune / p_total) / (core_immune / c_total)
                               if core_immune > 0 and p_total and c_total else np.nan)

        # Distances are from the islet's core centroid, so an islet with no core cell has none.
        if p_total and c_total:
            cx = float(c["X_centroid"].mean())
            cy = float(c["Y_centroid"].mean())
            px = p["X_centroid"].to_numpy(float)
            py = p["Y_centroid"].to_numpy(float)
            labels = p["phenotype"].to_numpy()
            dist = np.sqrt((px - cx) ** 2 + (py - cy) ** 2)
            immune_mask = np.isin(labels, IMMUNE_TYPES)
            rec["min_dist_immune_mean"] = float(dist[immune_mask].min()) if immune_mask.any() else np.nan
            for pheno in phenotypes:
                mask = labels == pheno
                rec[f"min_dist_{_safe(pheno)}"] = float(dist[mask].min()) if mask.any() else np.nan
        else:
            rec["min_dist_immune_mean"] = np.nan
            for pheno in phenotypes:
                rec[f"min_dist_{_safe(pheno)}"] = np.nan

        records.append(rec)
    return pd.DataFrame(records)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    # One phenotype vocabulary across the cohort, so every islet has every column.
    phenotypes = sorted({p for donor in donors
                         for p in _donor_cells(donor)["phenotype"].dropna().unique()})
    print(f"phenotypes: {', '.join(phenotypes)}\n")

    frames = []
    for donor in donors:
        frame = donor_metrics(donor, phenotypes)
        peri_present = int((frame["total_cells_peri"] > 0).sum())
        print(f"  {donor}: {len(frame):>5,} islets, {peri_present:>5,} with peri cells, "
              f"median peri immune z "
              f"{np.nanmedian(frame['enrich_z_Immune']):+.2f}", flush=True)
        frames.append(frame)

    out = pd.concat(frames, ignore_index=True)
    target = APP / "neighborhood_metrics.csv"
    out.to_csv(target, index=False)

    families = {name: len([c for c in out.columns if c.startswith(name)])
                for name in ("peri_prop_", "peri_count_", "enrich_z_", "min_dist_")}
    print(f"\n{len(out):,} islets x {out.shape[1]} columns "
          + ", ".join(f"{n} {name}*" for name, n in families.items()))
    print(f"-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
