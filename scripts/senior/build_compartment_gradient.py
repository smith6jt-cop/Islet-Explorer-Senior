#!/usr/bin/env python3
"""The compartment gradient: islet core, peri-islet expansion, exocrine tissue, pLN.

Part 3 of the study question -- immune cells near and in islets against those in
the pancreatic lymph node and in the exocrine tissue.

What this builds
----------------
``data/development/integration/compartments/`` holds, per donor:

* ``compartment_composition.csv`` -- immune fraction and immune-subtype
  composition for every compartment, one row per donor x compartment x subtype.
* ``islet_size_strata.csv`` -- the same core-compartment numbers recomputed on
  islets of at least ``--min-islet-cells`` core cells, which is the stratum in
  which an islet-level statement is defensible at all.
* ``batch_eta_squared.csv`` -- eta-squared of every cross-donor quantity on
  ``stain_group`` and ``acq_run``, each against its own permutation null.
* ``manifest.json`` -- the reasons, the denominators and the withheld donors.

Three compartments are pancreas and one is not
----------------------------------------------
The pancreas compartments are **spatial**: core, the 20 um peri-islet expansion
and the remaining exocrine tissue are cut out of one section by one signed
distance.  The pLN compartment is **compositional only**.  A lymph node is a
different block from the pancreas beside it -- for most donors a different slide
entirely -- so it has no islet, and a signed distance to one is undefined there
rather than zero.  The pLN column joins to the pancreas columns on ``donor_id``
and ``donor_status`` and on nothing else.  No distance, adjacency,
neighbourhood or KDE crosses the tissue boundary, and nothing is co-plotted in
space across tissues.  ``refuse_cross_tissue`` makes that structural: asking
this module for a pLN distance raises.

  One exception is coming and is not inferred here.  The maintainer is
  annotating the pLN portions of pancreas slides as a distinct dataset, and some
  of those lymphoid structures are anatomically attached to the section, where
  an islet-to-lymphoid distance *is* meaningful.  That classification arrives as
  data.  ``LYMPHOID_ATTACHMENT_COLUMN`` names the column this reader will honour
  when it lands; until then every pancreas/pLN pair is treated as separate.

Compartments come from the signed distance, never from an annotation
--------------------------------------------------------------------
Every expansion is 20 um.  One donor's drawn expansion annotation measures
23.1 um against 18.8-19.3 for every other, and all of them are named
``_exp20um``, so the annotation cannot be trusted to mean what it is called.
The compartment is therefore cut from ``signed_dist_islet_um``:

    core       signed distance <= 0
    peri       0 < signed distance <= 20
    exocrine   signed distance > 20

QuPath's ``Parent`` is never read.  Where ``typing/islet_membership`` exists the
geometric outline membership is read back as an audit and the disagreement is
reported per donor -- it is an audit, not an override.

Two stratifications, both measured, both applied before anything pools
----------------------------------------------------------------------
1. The endocrine shape signal is real only above roughly ten cells per islet, so
   every core-compartment quantity is emitted a second time over islets with at
   least ``--min-islet-cells`` core cells.
2. 6539's islets are mostly one-cell classifier fragments.  It is measured, it
   is reported, and it is flagged ``fragmented_islets`` so that no cohort-pooled
   per-islet number silently carries it.

Neither stratification drops a donor.  A cohort-pooled per-islet number without
both is wrong; a donor withheld for being inconvenient is worse.

Labels
------
PhenoCycler pancreas and PhenoCycler pLN are named by different lanes, so their
labels are harmonised through
``core/phenocycler/data/integration/vocab_crosswalk.csv`` (generated; never
hand-edited) to the common lineages ``T_NK``, ``B_Plasma``, ``Myeloid`` and
``Granulocyte``.  A lineage the crosswalk does not classify is withheld by name
and never defaulted to non-immune.

Negative evidence is evidence, and "unavailable" is reserved for a measurement
that does not exist.  A donor whose pLN export carries no intensity measurement
has nothing to cluster and nothing to name; that is unavailable with a reason,
and it is never zero immune cells.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

CELLS = REPO / "data" / "development" / "typing" / "cells"
MEMBERSHIP = REPO / "data" / "development" / "typing" / "islet_membership"
PLN_PHENO = REPO / "data" / "development" / "typing" / "pln" / "pheno"
PLN_XENIUM = REPO / "data" / "development" / "xenium_typing" / "release_pln_immune_v2"
DONOR_META = REPO / "data" / "development" / "typing" / "app" / "donor_metadata.csv"
ACQ_RUNS = REPO / "data" / "donor_acquisition_runs.csv"
CROSSWALK = REPO / "core" / "phenocycler" / "data" / "integration" / "vocab_crosswalk.csv"
OUT = REPO / "data" / "development" / "integration" / "compartments"

#: The maintainer's pLN-attachment classification is not in the checkout yet.
#: When it arrives it arrives as a column on the pLN section table; this is the
#: name this reader will honour.  Nothing infers it.
LYMPHOID_ATTACHMENT_COLUMN = "attached_to_pancreas_section"

EXPANSION_UM = 20.0
STATUS_ORDER = ["ND", "Aab+", "T1D"]
PANCREAS_COMPARTMENTS = ["core", "peri", "exocrine"]
COMMON_LINEAGES = ["T_NK", "B_Plasma", "Myeloid", "Granulocyte"]

#: PhenoCycler pancreas fine immune calls -> common lineage.  NK is T_NK because
#: the crosswalk collapses T and NK on the Xenium side under one common lineage;
#: the fine call is kept in the subtype column, not thrown away.
PANC_IMMUNE_TO_COMMON = {
    "Imm-T-cell": "T_NK",
    "Imm-NK-cell": "T_NK",
    "Imm-B-cell": "B_Plasma",
    "Imm-plasma-cell": "B_Plasma",
    "Imm-macrophage": "Myeloid",
    "Imm-dendritic-cell": "Myeloid",
    "Imm-neutrophil": "Granulocyte",
}

#: PhenoCycler pLN cluster-named immune classes -> common lineage.
PLN_CLASS_TO_COMMON = {
    "T_cell": "T_NK",
    "NK_ILC": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}

#: Xenium pLN release labels -> common lineage.
XENIUM_LABEL_TO_COMMON = {
    "Immune_T": "T_NK",
    "Immune_NK_ILC": "T_NK",
    "Immune_B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}


class CrossTissueDistance(RuntimeError):
    """Raised when something asks for a distance that spans a tissue boundary."""


def refuse_cross_tissue(tissue: str, quantity: str) -> None:
    """Structural refusal, not a convention.

    A lymph node has no islet.  A signed distance to one is undefined there, not
    zero, so the refusal is an exception rather than a NaN that a later mean
    could swallow.
    """
    if tissue != "pancreas":
        raise CrossTissueDistance(
            f"{quantity} is a spatial quantity and {tissue!r} is not the pancreas "
            "section it was measured on. The pLN is a different block: it has no "
            "islet, so a signed distance to one is undefined, not zero. Join on "
            "donor_id and donor_status instead."
        )


def read_parquet_retrying(path: Path, columns=None, attempts: int = 5):
    """A concurrent workflow writes the pLN tree while this reads it.

    A partial read is a wrong number, not a small one, so retry and then withhold
    -- never record what was on disk halfway through a write.
    """
    import time

    last = None
    for attempt in range(attempts):
        try:
            return pd.read_parquet(path, columns=columns)
        except Exception as exc:  # pragma: no cover - timing dependent
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"could not read {path} after {attempts} attempts: {last}")


def eligible_donors() -> list[str]:
    donors = sorted(p.name.split("=", 1)[1] for p in CELLS.glob("donor_id=*"))
    return [d for d in donors if d not in DONOR_EXCLUSIONS]


def load_donor_meta() -> pd.DataFrame:
    meta = pd.read_csv(DONOR_META, dtype={"case_id": str})
    meta = meta.rename(columns={"case_id": "donor_id"})
    acq = pd.read_csv(ACQ_RUNS, dtype={"donor_id": str})
    acq = acq[acq["excluded"].astype(str).str.lower() != "true"]
    # One donor can carry several scans; the compartment quantities are per
    # donor, so the batch covariate is the donor's scan when it is unambiguous
    # and unrecorded when it is not.  Unrecorded is not a default value.
    agg = acq.groupby("donor_id").agg(
        acq_run=("acq_run", lambda s: s.iloc[0] if s.nunique() == 1 else pd.NA),
        stain_group=("stain_group", lambda s: s.iloc[0] if s.nunique() == 1 else pd.NA),
    )
    return meta.merge(agg, on="donor_id", how="left")


def load_crosswalk_immune() -> tuple[set[str], dict[str, str]]:
    """Which common lineages the crosswalk records as immune in the pLN."""
    cw = pd.read_csv(CROSSWALK)
    lin = cw[cw["kind"] == "lineage"]
    immune, scope = set(), {}
    for _, row in lin.iterrows():
        common = str(row["common"])
        if str(row["is_immune"]).strip().lower() == "yes":
            immune.add(common)
            scope[common] = str(row.get("immune_tissue_scope", ""))
    return immune, scope


# ---------------------------------------------------------------------------
# pancreas


def compartment_from_signed_distance(dist: pd.Series) -> pd.Series:
    """core / peri / exocrine, cut from the signed distance alone.

    The drawn expansion annotation is not consulted: one donor's is 23.1 um
    while every other is 18.8-19.3 um, and all of them are named _exp20um.
    A cell with no signed distance has an unavailable compartment, never a
    default one.
    """
    out = pd.Series(pd.NA, index=dist.index, dtype="object")
    have = dist.notna()
    out[have & (dist <= 0)] = "core"
    out[have & (dist > 0) & (dist <= EXPANSION_UM)] = "peri"
    out[have & (dist > EXPANSION_UM)] = "exocrine"
    return out


def pancreas_donor(donor: str, min_islet_cells: int) -> dict:
    path = CELLS / f"donor_id={donor}" / "data.parquet"
    cols = [
        "object_id", "cell_class", "immune_definition", "excluded_reason",
        "qc_analysis_eligible", "signed_dist_islet_um", "islet_num",
        "cell_region_20um", "acq_run", "stain_group",
    ]
    df = read_parquet_retrying(path, columns=cols)
    refuse_cross_tissue("pancreas", "signed_dist_islet_um")

    n_rows = len(df)
    # A cell with no nucleus never reached this tree; what is here and still
    # ineligible is excluded by name, and the name is kept.
    keep = df["qc_analysis_eligible"].fillna(False)
    keep &= df["excluded_reason"].fillna("").str.strip() == ""
    dropped = df.loc[~keep, "excluded_reason"].fillna("ineligible").value_counts().to_dict()
    df = df.loc[keep].copy()

    df["compartment"] = compartment_from_signed_distance(df["signed_dist_islet_um"])
    no_distance = int(df["compartment"].isna().sum())

    df["is_immune"] = df["cell_class"] == "Immune"
    fine = df["immune_definition"].fillna("").str.strip()
    df["fine"] = fine.where(df["is_immune"], "")
    df["common"] = df["fine"].map(PANC_IMMUNE_TO_COMMON)
    unmapped = sorted(set(df.loc[df["is_immune"] & df["common"].isna(), "fine"]) - {""})

    rows = []
    for comp in PANCREAS_COMPARTMENTS:
        sub = df[df["compartment"] == comp]
        rows.extend(_composition_rows(donor, "pancreas", "phenocycler", comp, sub, "all_islets"))

    # Stratum 2: only islets with at least min_islet_cells core cells.  An
    # islet-level statement below that is a statement about classifier
    # fragments, so it is computed separately rather than pooled in.
    core = df[df["compartment"] == "core"]
    sizes = core["islet_num"].dropna().value_counts()
    big = set(sizes[sizes >= min_islet_cells].index)
    core_big = core[core["islet_num"].isin(big)]
    strat_rows = _composition_rows(
        donor, "pancreas", "phenocycler", "core", core_big,
        f"islets_ge_{min_islet_cells}_core_cells",
    )

    # The peri-islet expansion inherits the same stratification: a 20 um
    # expansion drawn around a one-cell fragment is not a peri-islet expansion.
    peri = df[df["compartment"] == "peri"]
    peri_big = peri[peri["islet_num"].isin(big)]
    strat_rows.extend(_composition_rows(
        donor, "pancreas", "phenocycler", "peri", peri_big,
        f"islets_ge_{min_islet_cells}_core_cells",
    ))

    n_islets = int(len(sizes))
    n_big = int(len(big))
    frag = float((sizes <= 1).sum() / len(sizes)) if len(sizes) else float("nan")
    median = float(sizes.median()) if n_islets else float("nan")
    # The flag is the measured condition, not a judgement: a donor whose median
    # islet is smaller than the threshold at which an islet-level statement is
    # defensible has islets that are mostly classifier fragments.
    flag = "fragmented_islets" if (n_islets and median < min_islet_cells) else ""

    audit = _membership_audit(donor, df)

    return {
        "rows": rows,
        "strat_rows": strat_rows,
        "islets": {
            "donor_id": donor,
            "n_islets": n_islets,
            f"n_islets_ge_{min_islet_cells}": n_big,
            "median_core_cells_per_islet": median,
            "fraction_islets_one_cell": frag,
            "flag": flag,
        },
        "diagnostics": {
            "donor_id": donor,
            "n_rows": n_rows,
            "n_eligible": int(len(df)),
            "dropped": dropped,
            "cells_without_signed_distance": no_distance,
            "immune_calls_not_in_crosswalk": unmapped,
            "acq_run": _one(df["acq_run"]),
            "stain_group": _one(df["stain_group"]),
            "membership_audit": audit,
        },
    }


def _one(series: pd.Series):
    vals = series.dropna().unique()
    return str(vals[0]) if len(vals) == 1 else None


def _membership_audit(donor: str, df: pd.DataFrame) -> dict:
    """Read the geometric outline membership back, where it exists.

    This is an audit of the distance cut, not an override of it, and it is
    reported rather than silently reconciled.  QuPath's Parent column is not
    read here or anywhere in this module.
    """
    path = MEMBERSHIP / f"donor_id={donor}" / "data.parquet"
    if not path.exists():
        return {"state": "unavailable",
                "reason": "no islet_membership partition for this donor"}
    m = read_parquet_retrying(path, columns=["object_id", "in_islet"])
    j = df[["object_id", "compartment"]].merge(m, on="object_id", how="left")
    covered = j["in_islet"].notna()
    core_cut = j["compartment"] == "core"
    return {
        "state": "measured",
        "n_covered": int(covered.sum()),
        "n_core_by_distance": int(core_cut.sum()),
        "n_core_by_outline": int((j["in_islet"] == True).sum()),  # noqa: E712
        "n_disagree": int((covered & (core_cut != (j["in_islet"] == True))).sum()),  # noqa: E712
    }


def _composition_rows(donor, tissue, modality, compartment, sub, stratum) -> list[dict]:
    """One denominator per compartment, repeated across subtype rows.

    ``n_cells`` is the eligible denominator for the compartment.  Do not sum it
    across subtype rows -- the subtype counts sum to the immune count, and the
    immune count plus the non-immune count is the denominator.  That internal
    coherence is checked by the test suite.
    """
    n = int(len(sub))
    n_immune = int(sub["is_immune"].sum()) if n else 0
    rows = [{
        "donor_id": donor, "tissue": tissue, "modality": modality,
        "compartment": compartment, "stratum": stratum,
        "lineage": "Immune_total", "n_cells": n, "n_lineage": n_immune,
        "n_immune": n_immune,
        "fraction_of_compartment": (n_immune / n) if n else float("nan"),
        "fraction_of_immune": 1.0 if n_immune else float("nan"),
        "state": "measured" if n else "unavailable",
        "reason": "" if n else f"no eligible cell in compartment {compartment}",
    }]
    for lineage in COMMON_LINEAGES:
        k = int((sub["common"] == lineage).sum()) if n else 0
        rows.append({
            "donor_id": donor, "tissue": tissue, "modality": modality,
            "compartment": compartment, "stratum": stratum,
            "lineage": lineage, "n_cells": n, "n_lineage": k, "n_immune": n_immune,
            "fraction_of_compartment": (k / n) if n else float("nan"),
            "fraction_of_immune": (k / n_immune) if n_immune else float("nan"),
            "state": "measured" if n else "unavailable",
            "reason": "" if n else f"no eligible cell in compartment {compartment}",
        })
    return rows


# ---------------------------------------------------------------------------
# pLN -- compositional only


def pln_phenocycler(donor: str) -> tuple[list[dict], dict]:
    base = PLN_PHENO / f"donor_id={donor}"
    status = pd.read_csv(PLN_PHENO / "per_donor_status.csv", dtype={"donor_id": str})
    row = status[status["donor_id"] == donor]
    if not len(row):
        return [], {"donor_id": donor, "state": "unavailable",
                    "reason": "donor has no pLN section in the pheno tree"}
    state = str(row.iloc[0]["state"])
    if state != "typed":
        # The measured reason travels with the row so the card can name it.
        # "the export carries no intensity measurement, so there is nothing to
        # cluster and nothing to name" is a statement about the measurement,
        # and it is never rendered as zero immune cells.
        reason = str(row.iloc[0].get("reason", "") or f"pLN state is {state}")
        rows = [{
            "donor_id": donor, "tissue": "pancreatic_lymph_node",
            "modality": "phenocycler", "compartment": "pln",
            "stratum": "not_applicable_no_islets", "lineage": lineage,
            "n_cells": pd.NA, "n_lineage": pd.NA, "n_immune": pd.NA,
            "fraction_of_compartment": float("nan"),
            "fraction_of_immune": float("nan"),
            "state": state if state in {"unavailable", "withheld", "failed"} else "unavailable",
            "reason": reason,
        } for lineage in ["Immune_total"] + COMMON_LINEAGES]
        return rows, {"donor_id": donor, "state": state, "reason": reason}

    frames = []
    for part in sorted(base.glob("roi=*/data.parquet")):
        frames.append(read_parquet_retrying(part, columns=["is_immune", "immune_class"]))
    if not frames:
        return [], {"donor_id": donor, "state": "unavailable",
                    "reason": "pLN status says typed but no roi partition is on disk"}
    df = pd.concat(frames, ignore_index=True)
    cls = df["immune_class"].astype(str).str.strip()
    df["common"] = cls.map(PLN_CLASS_TO_COMMON)
    df["is_immune"] = df["is_immune"].fillna(False)
    unmapped = sorted(set(cls[df["is_immune"] & df["common"].isna()]) - {"", "nan"})
    rows = _composition_rows(donor, "pancreatic_lymph_node", "phenocycler",
                             "pln", df, "all_islets")
    for r in rows:
        r["stratum"] = "not_applicable_no_islets"
    return rows, {"donor_id": donor, "state": "measured",
                  "n_cells": int(len(df)), "classes_not_in_crosswalk": unmapped}


def pln_xenium() -> tuple[list[dict], dict]:
    """The RNA view of the same compartment, as its own modality column.

    It is never merged with the protein pLN numbers and never tuned against
    them.  Agreement is informative and disagreement is informative.
    Cell ids do not cross: a QuPath UUID and a 10x cell_id share no namespace,
    so nothing here joins on a cell.
    """
    path = PLN_XENIUM / "per_donor_summary.csv"
    if not path.exists():
        return [], {"state": "unavailable", "reason": "no Xenium pLN release on disk"}
    s = read_csv_retrying(path)
    s["donor_id"] = s["donor_id"].astype(str)
    rows = []
    col_map = {
        "T_NK": ["s2__Immune_T", "s2__Immune_NK_ILC"],
        "B_Plasma": ["s2__Immune_B_Plasma"],
        "Myeloid": ["s2__Myeloid"],
        "Granulocyte": ["s2__Granulocyte"],
    }
    for _, r in s.iterrows():
        # The denominator excludes cells with no nucleus, which are not cells we
        # can be sure of, and excludes nothing else by value.
        n = int(r["cells"]) - int(r.get("no_nucleus_excluded", 0) or 0)
        n_immune = int(r["immune"])
        rows.append({
            "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
            "modality": "xenium", "compartment": "pln",
            "stratum": "not_applicable_no_islets", "lineage": "Immune_total",
            "n_cells": n, "n_lineage": n_immune, "n_immune": n_immune,
            "fraction_of_compartment": n_immune / n if n else float("nan"),
            "fraction_of_immune": 1.0 if n_immune else float("nan"),
            "state": "measured", "reason": "",
        })
        for lineage, cols in col_map.items():
            have = [c for c in cols if c in s.columns]
            if not have:
                rows.append({
                    "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
                    "modality": "xenium", "compartment": "pln",
                    "stratum": "not_applicable_no_islets", "lineage": lineage,
                    "n_cells": n, "n_lineage": pd.NA, "n_immune": n_immune,
                    "fraction_of_compartment": float("nan"),
                    "fraction_of_immune": float("nan"),
                    "state": "unavailable",
                    "reason": f"{lineage} has no column in the Xenium pLN release",
                })
                continue
            k = int(sum(int(r[c] or 0) for c in have))
            rows.append({
                "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
                "modality": "xenium", "compartment": "pln",
                "stratum": "not_applicable_no_islets", "lineage": lineage,
                "n_cells": n, "n_lineage": k, "n_immune": n_immune,
                "fraction_of_compartment": k / n if n else float("nan"),
                "fraction_of_immune": k / n_immune if n_immune else float("nan"),
                "state": "measured", "reason": "",
            })
    return rows, {"state": "measured", "n_donors": int(s["donor_id"].nunique())}


def read_csv_retrying(path: Path, attempts: int = 5) -> pd.DataFrame:
    import time
    last = None
    for attempt in range(attempts):
        try:
            return pd.read_csv(path)
        except Exception as exc:  # pragma: no cover
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"could not read {path}: {last}")


# ---------------------------------------------------------------------------
# batch


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """Between-group share of the total sum of squares. NaN when undefined."""
    ok = ~pd.isna(values) & ~pd.isna(groups)
    v, g = np.asarray(values, float)[ok], np.asarray(groups)[ok]
    if v.size < 3 or len(set(g.tolist())) < 2:
        return float("nan")
    grand = v.mean()
    sst = float(((v - grand) ** 2).sum())
    if sst <= 0:
        return float("nan")
    ssb = 0.0
    for level in set(g.tolist()):
        sel = v[g == level]
        ssb += sel.size * (sel.mean() - grand) ** 2
    return float(ssb / sst)


def eta_with_null(values, groups, n_perm=2000, seed=0) -> dict:
    """Eta-squared against a permutation null.

    With 14-20 donors over 9-11 acquisition levels, a large eta-squared is the
    *expected* value under no effect at all -- the null median says so.  So the
    observed statistic is reported next to that null and its permutation p, and
    neither is allowed to gate anything.  A confounded quantity stays in.
    """
    obs = eta_squared(values, groups)
    if not np.isfinite(obs):
        return {"eta_squared": float("nan"), "null_median": float("nan"),
                "p_permutation": float("nan"), "n_perm": 0}
    rng = np.random.default_rng(seed)
    ok = ~pd.isna(values) & ~pd.isna(groups)
    v, g = np.asarray(values, float)[ok], np.asarray(groups)[ok]
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = eta_squared(rng.permutation(v), g)
    return {
        "eta_squared": obs,
        "null_median": float(np.nanmedian(null)),
        "p_permutation": float((np.nansum(null >= obs) + 1) / (n_perm + 1)),
        "n_perm": n_perm,
        "n_levels": int(len(set(g.tolist()))),
        "n_donors_used": int(v.size),
    }


def batch_table(comp: pd.DataFrame, meta: pd.DataFrame, n_perm: int) -> pd.DataFrame:
    m = meta.set_index("donor_id")
    out = []
    keys = comp[["tissue", "modality", "compartment", "stratum", "lineage"]].drop_duplicates()
    for _, k in keys.iterrows():
        sel = comp
        for col in ["tissue", "modality", "compartment", "stratum", "lineage"]:
            sel = sel[sel[col] == k[col]]
        sel = sel[sel["state"] == "measured"]
        if len(sel) < 3:
            continue
        vals = sel["fraction_of_compartment"].to_numpy(float)
        donors = sel["donor_id"].tolist()
        for grouping in ["stain_group", "acq_run", "donor_status"]:
            g = np.array([m[grouping].get(d, pd.NA) if d in m.index else pd.NA
                          for d in donors], dtype=object)
            stat = eta_with_null(vals, g, n_perm=n_perm)
            out.append({
                "quantity": f"{k['compartment']}_{k['lineage']}_fraction",
                "tissue": k["tissue"], "modality": k["modality"],
                "compartment": k["compartment"], "stratum": k["stratum"],
                "lineage": k["lineage"], "grouping": grouping,
                "n_donors_unrecorded": int(sum(pd.isna(x) for x in g)),
                **stat,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--min-islet-cells", type=int, default=10,
                    help="core cells an islet needs before an islet-level number "
                         "is defensible; the measured threshold is about 10")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--domain-root", default=str(REPO / "data" / "development" / "typing" / "app"),
                    help="the app domain whose cells these are. It is recorded in the "
                         "manifest so the card can refuse an artifact built for the "
                         "other domain: the two datasets share five donors and no "
                         "vocabulary, and a donor id is not a licence to read.")
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = load_donor_meta()
    immune_lineages, immune_scope = load_crosswalk_immune()
    donors = args.donors or eligible_donors()

    rows, strat, islets, diags, pln_notes = [], [], [], [], []
    for donor in donors:
        res = pancreas_donor(donor, args.min_islet_cells)
        rows.extend(res["rows"])
        strat.extend(res["strat_rows"])
        islets.append(res["islets"])
        diags.append(res["diagnostics"])
        pln_rows, note = pln_phenocycler(donor)
        rows.extend(pln_rows)
        pln_notes.append(note)

    xen_rows, xen_note = pln_xenium()
    cohort = set(donors)
    xen_in = [r for r in xen_rows if r["donor_id"] in cohort]
    xen_out = sorted({r["donor_id"] for r in xen_rows} - cohort)
    rows.extend(xen_in)

    comp = pd.DataFrame(rows)
    comp = comp.merge(meta[["donor_id", "donor_status", "stain_group", "acq_run"]],
                      on="donor_id", how="left")
    comp["donor_status"] = pd.Categorical(comp["donor_status"], STATUS_ORDER, ordered=True)
    comp = comp.sort_values(["tissue", "modality", "compartment", "stratum",
                             "donor_status", "donor_id", "lineage"])
    comp.to_csv(out / "compartment_composition.csv", index=False)

    st = pd.DataFrame(strat)
    st = st.merge(meta[["donor_id", "donor_status"]], on="donor_id", how="left")
    st.to_csv(out / "islet_size_strata.csv", index=False)
    pd.DataFrame(islets).to_csv(out / "islet_size_summary.csv", index=False)

    eta = batch_table(pd.concat([comp, st.assign(tissue="pancreas",
                                                 modality="phenocycler")],
                                ignore_index=True),
                      meta, args.n_perm)
    eta.to_csv(out / "batch_eta_squared.csv", index=False)

    manifest = {
        "generated_by": "scripts/senior/build_compartment_gradient.py",
        "source_domain_root": str(Path(args.domain_root).resolve()),
        "compartments": {
            "core": "signed_dist_islet_um <= 0",
            "peri": f"0 < signed_dist_islet_um <= {EXPANSION_UM}",
            "exocrine": f"signed_dist_islet_um > {EXPANSION_UM}",
            "pln": "compositional only; no distance exists in a lymph node",
        },
        "expansion_from": "signed distance, never the drawn annotation (one donor's "
                          "is 23.1 um against 18.8-19.3 for every other, and all are "
                          "named _exp20um)",
        "parent_column": "never read; membership is geometric",
        "cross_tissue": "pancreas and pLN join on donor_id and donor_status ONLY. "
                        "No distance, adjacency, neighbourhood or KDE crosses the "
                        "boundary. refuse_cross_tissue() makes that structural.",
        "pending_maintainer_classification": {
            "column": LYMPHOID_ATTACHMENT_COLUMN,
            "note": "lymphoid structures anatomically attached to a pancreas section "
                    "will license an islet-to-lymphoid distance. It arrives as data "
                    "and is not inferred here; until it lands every pancreas/pLN pair "
                    "is separate.",
        },
        "stratifications": {
            "min_islet_cells": args.min_islet_cells,
            "why": "the endocrine shape signal is real only above about 10 cells per "
                   "islet; 6539's islets are mostly one-cell classifier fragments",
        },
        "donor_exclusions_from": "phenocycler.cohort.DONOR_EXCLUSIONS",
        "donors_excluded": sorted(DONOR_EXCLUSIONS),
        "donors_built": donors,
        "immune_lineages_from_crosswalk": sorted(immune_lineages),
        "immune_tissue_scope": immune_scope,
        "pln_phenocycler": pln_notes,
        "pln_xenium": xen_note,
        "pln_xenium_donors_outside_cohort": xen_out,
        "pancreas_diagnostics": diags,
        "batch": "eta-squared on stain_group and acq_run, each against a "
                 "permutation null; coarse batch is aliased with disease status so "
                 "it is not corrected on. Nothing is gated on these numbers.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
