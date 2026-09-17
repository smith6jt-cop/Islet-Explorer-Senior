"""Readers for the identity-dependent figures: the two development releases, the verdicts, and the
integration stage outputs.  Every reader returns None (or an empty frame) with a stated reason when
its input does not exist yet, so a figure can say "not available: <why>" instead of failing."""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402

BROADTYPE = pins.REPO / "data" / "development" / "broadtype"
XENIUM_TYPING = pins.REPO / "data" / "development" / "xenium_typing"
RELEASE = {"phenocycler": BROADTYPE / "release", "xenium": XENIUM_TYPING / "release"}
INTEGRATION = pins.PKG / "data" / "integration"
AMBIGUOUS, UNAVAILABLE = "Ambiguous", "Unavailable"


def release_json(modality: str) -> dict | None:
    p = RELEASE[modality] / "release.json"
    return json.loads(p.read_text()) if p.exists() else None


def label_space(modality: str) -> dict:
    p = RELEASE[modality] / "label_space.json"
    if not p.exists():
        return {}
    return {r["label"]: r for r in json.loads(p.read_text())["labels"]}


def release_calls(modality: str, donor: str, columns=None) -> pd.DataFrame | None:
    p = RELEASE[modality] / "calls" / f"donor_id={donor}" / "roi=panc" / "data.parquet"
    if not p.exists():
        return None
    return pd.read_parquet(p, columns=columns)


def release_inputs(modality: str) -> list:
    d = RELEASE[modality]
    return [p for p in (d / "release.json", d / "label_space.json") if p.exists()]


def producer(modality: str) -> str | None:
    rj = release_json(modality)
    if rj is None:
        return None
    if modality == "phenocycler":
        return rj.get("arm", {}).get("name") or rj.get("producer", {}).get("lane", "").replace("broadtype-", "")
    return rj.get("arm", {}).get("name")


def verdict(modality: str, arm: str, donor: str) -> dict | None:
    if modality == "phenocycler":
        p = BROADTYPE / "verdict" / arm / f"donor_id={donor}.json"
    else:
        p = XENIUM_TYPING / "gate" / arm / f"donor_id={donor}.json"
    return json.loads(p.read_text()) if p.exists() else None


def verified_labels(modality: str, arm: str, donor: str) -> set:
    v = verdict(modality, arm, donor)
    if v is None:
        return set()
    if modality == "phenocycler":
        return {L for L, r in v.get("labels", {}).items() if r.get("verified")}
    return set(v.get("verified_lineages", []))


def verdict_summary() -> dict | None:
    p = BROADTYPE / "verdict" / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


def cells_table(modality: str, donor: str, columns=None) -> pd.DataFrame | None:
    sub = "cells_pheno" if modality == "phenocycler" else "cells_xen"
    p = INTEGRATION / sub / f"donor_id={donor}" / "roi=panc" / "data_0.parquet"
    if not p.exists():
        return None
    return pd.read_parquet(p, columns=columns)


def structures(modality: str, donor: str) -> pd.DataFrame | None:
    p = INTEGRATION / "structures" / modality / f"donor_id={donor}" / "roi=panc" / "islets.parquet"
    return pd.read_parquet(p) if p.exists() else None


def registration_qc() -> pd.DataFrame | None:
    p = INTEGRATION / "qc" / "registration_qc.csv"
    return pd.read_csv(p) if p.exists() else None


def composition() -> pd.DataFrame | None:
    p = INTEGRATION / "qc" / "composition.csv"
    return pd.read_csv(p) if p.exists() else None


def paired_islets() -> pd.DataFrame | None:
    p = INTEGRATION / "paired" / "islets.parquet"
    return pd.read_parquet(p) if p.exists() else None


def niche_profiles() -> pd.DataFrame | None:
    p = INTEGRATION / "niches" / "niche_profiles.csv"
    return pd.read_csv(p) if p.exists() else None


def registration_transform(donor: str) -> dict | None:
    p = INTEGRATION / "registration" / f"donor_id={donor}" / "roi=panc" / "transform.json"
    return json.loads(p.read_text()) if p.exists() else None


def restain_registration(donor: str) -> dict | None:
    """The same-section (IF <-> Xenium) registration record for the multi-platform text box."""
    p = pins.REPO / "data" / "development" / "restain" / "refine_panc_0.5um.json"
    if not p.exists():
        return None
    rows = [r for r in json.loads(p.read_text()) if str(r.get("donor_id")) == str(donor)]
    return rows[0] if rows else None


def protein_rna(donor: str) -> dict | None:
    p = pins.REPO / "data" / "development" / "restain" / "protein_rna_panc.json"
    if not p.exists():
        return None
    rows = [r for r in json.loads(p.read_text()) if str(r.get("donor_id")) == str(donor)]
    return rows[0] if rows else None


def common_of(modality: str) -> dict:
    """producer label -> common lineage, from the release's label space."""
    return {L: r.get("common_lineage", "") for L, r in label_space(modality).items() if "common_lineage" in r}


def comparable_lineages() -> tuple:
    try:
        sys.path.insert(0, str(pins.PKG))
        from phenocycler.integration.donor import comparable_for
        return tuple(comparable_for("pancreas"))
    except Exception:  # noqa: BLE001
        return ("Endocrine", "Stromal", "Vascular", "T_NK", "B_Plasma", "Myeloid", "Neural")


def common_colors() -> dict:
    try:
        sys.path.insert(0, str(pins.PKG))
        from phenocycler.integration.figures import COMMON_COLORS
        return dict(COMMON_COLORS)
    except Exception:  # noqa: BLE001
        return {"Endocrine": "#CCBB44", "Exocrine": "#4477AA", "Stromal": "#EE6677", "Vascular": "#66CCEE", "T_NK": "#228833",
                "B_Plasma": "#66AA55", "Myeloid": "#117733", "Granulocyte": "#E69F00", "Neural": "#B5838D", "Unassigned": "#BBBBBB"}


def label_color(label: str, common: dict, colors: dict) -> str:
    if label == UNAVAILABLE:
        return "#BBBBBB"
    if label == AMBIGUOUS:
        return "#888888"
    return colors.get(common.get(label, label), "#999999")


def stacked_bars(ax, table: pd.DataFrame, order: list, common: dict, colors: dict, *, x_positions=None, width=0.6, hatch_ambiguous=True):
    """``table``: index = bar name, columns = labels (fractions).  Ambiguous hatched, Unavailable grey."""
    x = np.arange(len(table)) if x_positions is None else np.asarray(x_positions)
    bottom = np.zeros(len(table))
    handles = []
    for L in order:
        if L not in table.columns:
            continue
        v = table[L].to_numpy(float)
        h = ax.bar(x, v, width, bottom=bottom, color=label_color(L, common, colors),
                   hatch="///" if (L == AMBIGUOUS and hatch_ambiguous) else None, edgecolor="white", linewidth=0.3, label=L)
        bottom += np.nan_to_num(v)
        handles.append(h)
    ax.set_xticks(x); ax.set_xticklabels(list(table.index), fontsize=8)
    ax.set_ylim(0, 1.0)
    return handles


# ------------------------------------------------------------------ lane modules under their own names
# Every lane has a ``pins.py`` (and broadtype/xenium_typing a ``gate.py``); importing one of them while the
# report's own ``pins`` sits in sys.modules hands the lane the WRONG pins.  These loaders import a lane's
# modules in isolation and restore the report's namespace afterwards.
_LANES: dict = {}
_SHADOWED = ("pins", "gate", "gates", "usability", "nullgate", "nullgate_fast", "criteria", "declarations", "records", "arm_base")


def _lane(name: str, modules: tuple) -> dict:
    import importlib
    if name in _LANES:
        return _LANES[name]
    lane = pins.PKG / "scripts" / name
    saved = {k: sys.modules.pop(k) for k in _SHADOWED if k in sys.modules}
    sys.path.insert(0, str(lane))
    try:
        out = {m: importlib.import_module(m) for m in modules}
    finally:
        sys.path.remove(str(lane))
        for k in _SHADOWED:
            sys.modules.pop(k, None)
        sys.modules.update(saved)
    _LANES[name] = out
    return out


def xenium_lane() -> dict:
    """``{'pins': scripts/xenium_typing/pins, 'gate': scripts/xenium_typing/gate}``."""
    return _lane("xenium_typing", ("pins", "gate"))


def broadtype_lane() -> dict:
    """``{'pins': scripts/broadtype/pins, 'gate': scripts/broadtype/gate}``."""
    return _lane("broadtype", ("pins", "gate"))
