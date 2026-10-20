"""Bake-off evidence figures (B01-B04): what every arm did on both donors, beside its nulls.

Read-only over the screen verdicts (``data/development/broadtype/verdict/<arm>/donor_id=<d>[__donor].json``)
and the Xenium gate records (``data/development/xenium_typing/gate/<arm>/donor_id=<d>.json``).  Nothing here
re-judges; the figures draw the numbers the verdicts recorded.

    B01  held-out marker enrichment per (arm, label, marker): h vs its null p95, both donors
    B02  spatial coherence per (arm, label): Moran's I vs its null p95
    B03  primary label composition per arm (stacked; Ambiguous hatched, Unavailable grey) + verified count
    B04  Xenium gate: award contrast vs null p95, Moran's I vs null, replicate stability, per arm x lineage
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import identity_io as io  # noqa: E402
import pins  # noqa: E402
import style  # noqa: E402
from manifest import FigureRecord  # noqa: E402

SEC = "bakeoff"
VERDICT = io.BROADTYPE / "verdict"
XGATE = io.XENIUM_TYPING / "gate"
LABELS = ("Endocrine", "Epithelial_nonendocrine", "Vascular", "Stromal", "T_NK", "B_Plasma", "Myeloid")
XLIN = ("Endocrine", "Exocrine", "Immune_B", "Immune_T", "Myeloid", "Neural", "Stromal", "Vascular")
OK, FAIL, UNJ = "#228833", "#EE6677", "#BBBBBB"


def _rec(fid, files, inputs, facts):
    return FigureRecord(figure_id=fid, donor="_pair", section=SEC, files=[str(f) for f in files], inputs=[str(p) for p in inputs], facts=facts)


def _side(out, stem, facts):
    p = out / f"{stem}.json"; p.write_text(json.dumps(facts, indent=1, default=str) + "\n"); return p


def verdicts(scope: str = "frame") -> dict:
    out = {}
    for p in sorted(VERDICT.glob("*/donor_id=*.json")):
        suffix = "" if scope == "frame" else f"__{scope}"
        if not p.name.endswith(f"{suffix}.json") or (scope == "frame" and "__" in p.stem):
            continue
        v = json.loads(p.read_text())
        if v.get("scope", "frame") != scope:
            continue
        out[(v["arm"], str(v["donor_id"]))] = v
    return out


def arms_of(vd: dict) -> list:
    arms = sorted({a for a, _ in vd})
    def key(a):
        n = sum(vd[(a, d)]["n_verified"] for d in pins.DONORS if (a, d) in vd)
        return (-n, a)
    return sorted(arms, key=key)


# ------------------------------------------------------------------ B01 held-out enrichment
def b01_heldout(prov: dict, scope: str = "frame") -> FigureRecord:
    vd = verdicts(scope)
    arms = arms_of(vd)
    gated = [L for L in LABELS if L != "Endocrine"]
    fig, axes = style.figure(len(arms), len(pins.DONORS), width=7.5, height=1.9)
    facts = {}
    for i, a in enumerate(arms):
        for j, d in enumerate(pins.DONORS):
            ax = axes[i, j]
            v = vd.get((a, d))
            if v is None:
                ax.axis("off"); style.note(ax, f"{a} {d}: not judged"); continue
            xs, hs, p95s, cols, names, x = [], [], [], [], [], 0
            for L in gated:
                r = v["labels"].get(L, {})
                for m in r.get("heldout", {}).get("per_marker", []):
                    xs.append(x); hs.append(m["h"] if m.get("judgeable") else np.nan); p95s.append(m["null_p"] if m.get("judgeable") else np.nan)
                    cols.append(OK if m.get("clear") else (FAIL if m.get("judgeable") else UNJ)); names.append(f"{m['marker']}"); x += 1
                x += 1
            if not xs:
                ax.set_yticks([]); ax.set_xticks([])
                style.note(ax, f"{a} - donor {d}: no held-out marker was judgeable ({v['n_verified']}/7 verified)\n"
                               + "  ".join(f"{L}:{v['labels'][L]['status']}" for L in LABELS), fontsize=6.5, width=95)
                facts[f"{a}|{d}"] = {L: {"status": v["labels"][L]["status"], "n_calls": v["labels"][L].get("n_primary_calls")} for L in LABELS}
                continue
            ax.bar(xs, np.nan_to_num(hs), color=cols, width=0.8)
            ax.plot(xs, p95s, "k_", markersize=6, mew=1.2, label="null p95 (200 size- and region-matched draws)")
            ax.axhline(1.0, color="#999999", lw=0.6, ls=":")
            ax.set_xticks(xs); ax.set_xticklabels(names, rotation=90, fontsize=5)
            top = float(np.nanmax(hs)) if np.isfinite(hs).any() else 4.0
            ax.set_ylim(0, max(4.0, top * 1.1))
            ax.set_ylabel("h", fontsize=7)
            ax.set_title(f"{a} - donor {d} ({style.status_display(pins.donor_status(d))}): {v['n_verified']}/7 verified\n"
                         + style.wrap("  ".join(f"{L}:{v['labels'][L]['status']}" for L in LABELS), 95), fontsize=6)
            facts[f"{a}|{d}"] = {L: {"status": v["labels"][L]["status"], "n_calls": v["labels"][L].get("n_primary_calls")} for L in LABELS}
    axes[0, 0].legend(fontsize=6, loc="upper right")
    fig.text(0.5, 0.955, "green = clears its null; red = judged, below null; grey = unjudgeable (< 20 tail cells in the frame); "
             "h = P(label | held-out marker tail) / P(label | any judged cell), the held-out marker removed from the arm's INPUT", ha="center", fontsize=7)
    style.title(fig, f"B01 held-out marker enrichment per arm, label and marker ({scope} scope) -- every bar beside its own null")
    style.footer(fig, prov); style.finish(fig, top=0.94, bottom=0.03)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, f"B01_heldout_enrichment_{scope}"); plt.close(fig)
    files.append(_side(out, f"B01_heldout_enrichment_{scope}", facts))
    return _rec("B01", files, sorted(VERDICT.glob("*/donor_id=*.json")), facts)


# ------------------------------------------------------------------ B02 spatial coherence
def b02_spatial(prov: dict, scope: str = "frame") -> FigureRecord:
    vd = verdicts(scope); arms = arms_of(vd)
    fig, axes = style.figure(1, len(pins.DONORS), width=7, height=4.5)
    facts = {}
    for j, d in enumerate(pins.DONORS):
        ax = axes[0, j]
        x = 0; xt, xl = [], []
        for a in arms:
            v = vd.get((a, d))
            if v is None:
                continue
            for L in LABELS:
                s = v["labels"][L].get("spatial") or {}
                if s.get("judgeable"):
                    ax.bar(x, s["I"], color=OK if s["clear"] else FAIL, width=0.8)
                    ax.plot([x], [s["null_p"]], "k_", markersize=6, mew=1.2)
                else:
                    ax.bar(x, 0.002, color=UNJ, width=0.8)
                xt.append(x); xl.append(f"{a[:14]} | {L[:10]}"); x += 1
                facts[f"{a}|{d}|{L}"] = {k: s.get(k) for k in ("I", "null_p", "n_calls", "clear", "judgeable")}
            x += 1
        ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=90, fontsize=4.5); ax.set_ylim(bottom=0)
        ax.set_ylabel("Moran's I of the one-hot label (k = 12)"); ax.set_title(f"donor {d} ({style.status_display(pins.donor_status(d))}); black tick = within-region label-permutation null p95", fontsize=8)
    style.title(fig, f"B02 spatial coherence per arm and label ({scope} scope); green clears the null, red does not, grey unjudgeable (< 50 calls)")
    style.footer(fig, prov); style.finish(fig, bottom=0.3)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, f"B02_spatial_coherence_{scope}"); plt.close(fig)
    files.append(_side(out, f"B02_spatial_coherence_{scope}", facts))
    return _rec("B02", files, sorted(VERDICT.glob("*/donor_id=*.json")), facts)


# ------------------------------------------------------------------ B03 composition per arm
def b03_composition(prov: dict, scope: str = "frame") -> FigureRecord:
    vd = verdicts(scope); arms = arms_of(vd)
    common = {L: L for L in LABELS}; common["Epithelial_nonendocrine"] = "Exocrine"
    colors = io.common_colors(); order = list(LABELS) + [io.AMBIGUOUS, io.UNAVAILABLE]
    fig, axes = style.figure(1, len(pins.DONORS), width=7, height=5)
    facts = {}
    for j, d in enumerate(pins.DONORS):
        ax = axes[0, j]
        rows = {}
        for a in arms:
            v = vd.get((a, d))
            if v is None:
                continue
            c = v.get("primary_label_counts") or {}; n = sum(c.values()) or 1
            rows[f"{a}\n{v['n_verified']}/7 verified"] = {L: c.get(L, 0) / n for L in order}
            facts[f"{a}|{d}"] = {"counts": c, "n_verified": v["n_verified"]}
        T = pd.DataFrame(rows).T
        io.stacked_bars(ax, T, order, common, colors)
        ax.tick_params(axis="x", labelsize=6, rotation=90); ax.set_ylabel("fraction of frame cells (primary run)")
        ax.set_title(f"donor {d} ({style.status_display(pins.donor_status(d))}), {n:,} frame cells per arm", fontsize=8)
    fig.legend(handles=[matplotlib.patches.Patch(color=io.label_color(L, common, colors), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in order], loc="lower center", bbox_to_anchor=(0.5, 0.04), ncol=5, fontsize=7)
    style.title(fig, f"B03 what each arm called on the same frozen frame ({scope} scope): composition is descriptive, never a criterion")
    style.footer(fig, prov); style.finish(fig, bottom=0.32)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, f"B03_composition_by_arm_{scope}"); plt.close(fig)
    files.append(_side(out, f"B03_composition_by_arm_{scope}", facts))
    return _rec("B03", files, sorted(VERDICT.glob("*/donor_id=*.json")), facts)


# ------------------------------------------------------------------ B04 Xenium gate
def b04_xenium(prov: dict) -> FigureRecord:
    recs = {}
    for p in sorted(XGATE.glob("*/donor_id=*.json")):
        v = json.loads(p.read_text()); recs[(v["arm"], str(v["donor_id"]))] = v
    arms = sorted({a for a, _ in recs})
    if not arms:
        fig, axes = style.figure(1, 1, width=7, height=3); axes[0, 0].axis("off"); style.note(axes[0, 0], "no Xenium gate records yet")
        style.title(fig, "B04 Xenium gate: not available"); style.footer(fig, prov); style.finish(fig)
        out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, "B04_xenium_gate"); plt.close(fig)
        return _rec("B04", files, [], {"not_available": "no gate records"})
    fig, axes = style.figure(len(arms), 2 * len(pins.DONORS), width=4.2, height=2.6)
    facts = {}
    for i, a in enumerate(arms):
        for j, d in enumerate(pins.DONORS):
            v = recs.get((a, d))
            axA, axB = axes[i, 2 * j], axes[i, 2 * j + 1]
            if v is None:
                axA.axis("off"); axB.axis("off"); style.note(axA, f"{a} {d}: not judged"); continue
            x = np.arange(len(XLIN))
            real = [v["lineages"][L]["award"].get("real_own") for L in XLIN]
            p95 = [v["lineages"][L]["award"].get("null_p95_own") for L in XLIN]
            col = [OK if v["lineages"][L]["verified"] else (FAIL if v["lineages"][L]["award"].get("judgeable") else UNJ) for L in XLIN]
            axA.bar(x, [np.nan_to_num(r if r is not None else np.nan) for r in real], color=col)
            axA.plot(x, [np.nan if q is None else q for q in p95], "k_", markersize=7, mew=1.2)
            axA.axhline(0, color="black", lw=0.5); axA.set_xticks(x); axA.set_xticklabels(XLIN, rotation=90, fontsize=5.5)
            axA.set_title(f"{a} - {d}: award contrast vs null p95 ({v['n_verified']}/8 verified)", fontsize=6.5); axA.set_ylabel("own-core contrast", fontsize=6)
            I = [v["lineages"][L]["spatial"].get("I") for L in XLIN]; Ip = [v["lineages"][L]["spatial"].get("null_p") for L in XLIN]
            stab = [v["lineages"][L]["stability"].get("invariant") for L in XLIN]
            axB.bar(x, [np.nan_to_num(q if q is not None else np.nan) for q in I], color=col)
            axB.plot(x, [np.nan if q is None else q for q in Ip], "k_", markersize=7, mew=1.2)
            for k, s in enumerate(stab):
                axB.text(k, 0, "stable" if s else ("unstable" if s is False else "unjudged"), rotation=90, fontsize=4.5, ha="center", va="bottom", color="black")
            axB.set_xticks(x); axB.set_xticklabels(XLIN, rotation=90, fontsize=5.5); axB.set_ylim(bottom=0)
            axB.set_title("Moran's I vs within-section permutation null; replicate stability printed", fontsize=6.5)
            facts[f"{a}|{d}"] = {L: {"verified": v["lineages"][L]["verified"], "n": v["label_counts"].get(L, 0), "awarded": v["lineages"][L]["award"].get("awarded"), "stable": stab[k]} for k, L in enumerate(XLIN)}
    style.title(fig, "B04 Xenium arms under one gate: green verified, red judged and failed, grey unjudgeable; every bar beside its null")
    style.footer(fig, prov); style.finish(fig, top=0.93, bottom=0.04)
    out = pins.OUT / "_pair" / SEC; files = style.savefig(fig, out, "B04_xenium_gate"); plt.close(fig)
    files.append(_side(out, "B04_xenium_gate", facts))
    return _rec("B04", files, sorted(XGATE.glob("*/donor_id=*.json")), facts)


def run(prov: dict, only=None) -> list:
    recs = []
    for fid, fn in (("B01", b01_heldout), ("B02", b02_spatial), ("B03", b03_composition)):
        if only and fid not in only:
            continue
        for scope in ("frame", "donor"):
            if verdicts(scope):
                r = fn(prov, scope); recs.append(r); print(f"[report] {fid} {scope}: {r.files[0]}", flush=True)
    if not only or "B04" in only:
        r = b04_xenium(prov); recs.append(r); print(f"[report] B04: {r.files[0]}", flush=True)
    return recs


if __name__ == "__main__":
    pins.assert_pins()
    from manifest import Manifest
    man = Manifest()
    for r in run(pins.provenance(), sys.argv[1:] or None):
        man.add_figure(r)
    man.write()
