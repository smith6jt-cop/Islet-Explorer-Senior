"""The B.2 Xenium release: a declared identity release that ``phenospace.load_release`` accepts.

    cd core/phenocycler && ../../.venv/bin/python scripts/xenium_typing/release.py --arm <arm> --donors 6436 6442

Cell set = EVERY bundle cell (``cells.parquet``; 695,113 / 958,744 -- exactly ``cells_xen``), because
``phenospace.attach_release_identity`` requires set equality with the geometry import.  Cells the base
builder dropped are ``Unavailable(qc_dropped:<';'-joined set>)`` with the QC recomputed per cell by
``build_phenotyped_base.qc_mask`` (the three masks overlap, so this is a set, never a partition).
``evidence_positive`` is a NON-NULL boolean (``donor._evidence_mask`` refuses NA): True iff the
lineage was VERIFIED for this arm on this donor in the gate json AND the cell is not low-signal.
``feat__if_INS`` / ``feat__if_GCG`` carry the re-stain lane's per-cell IF protein values where they
exist -- continuous evidence on the same section; no positivity is implied (Q4).  No provenance string
may name a retired input (``phenospace.RETIRED_INPUT_TOKENS``).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "xenium_celltype"))
import gate  # noqa: E402
import pins  # noqa: E402
sys.path.insert(0, str(pins.PKG))
from phenocycler.integration.phenospace import RETIRED_INPUT_TOKENS  # noqa: E402
from phenocycler.integration.vocab import XENIUM_TO_COMMON  # noqa: E402

SCHEMA_VERSION = 1
ROI = pins.ROI


def _sha(path: pathlib.Path) -> str:
    return pins.sha256_file(path)


def qc_dropped_reasons(donor: str, bundle_ids: np.ndarray) -> pd.Series:
    """'qc_dropped:<set>' per bundle cell the base excluded, recomputed with the builder's own rule."""
    import xenio  # scripts/xenium_celltype
    import build_phenotyped_base as bpb
    dd = xenio.load_donor(donor, pins.XENIUM_PATHS, roi=ROI, prefer="bundle")
    cells = pd.read_parquet(pins.cells_parquet(donor), columns=["cell_id", "cell_area"])
    area = cells.set_index(cells["cell_id"].astype(str))["cell_area"].reindex(np.asarray(dd.cell_id).astype(str)).to_numpy(float)
    counts = dd.counts
    n_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    density = np.where(area > 0, n_counts / np.where(area > 0, area, 1.0), 0.0)
    cap = float(np.quantile(density[np.isfinite(density)], bpb.DENSITY_QUANTILE)) if hasattr(bpb, "DENSITY_QUANTILE") else float(np.quantile(density[np.isfinite(density)], 0.98))
    min_genes = getattr(bpb, "MIN_GENES", 20)
    min_counts = getattr(bpb, "MIN_COUNTS", 50)
    fails = []
    for g, c, d_ in zip(n_genes, n_counts, density):
        f = []
        if g < min_genes:
            f.append("min_genes")
        if c < min_counts:
            f.append("min_counts")
        if d_ > cap:
            f.append("density")
        fails.append(";".join(f))
    s = pd.Series(fails, index=np.asarray(dd.cell_id).astype(str))
    return s.reindex(bundle_ids).fillna("not_in_matrix")


def label_space(arm: str) -> dict:
    labels = [{"label": L, "common_lineage": XENIUM_TO_COMMON[L],
               "reason": "Xenium broad lineage from the identity cores (core/panels snapshot)"} for L in pins.LINEAGES]
    labels += [{"label": gate.AMBIGUOUS, "common_resolution": "ambiguous", "reason": "the arm did not resolve a unique lineage"},
               {"label": gate.UNAVAILABLE, "common_resolution": "unavailable", "reason": "identity measurement unavailable (QC-dropped, not in the analysis, or the arm emitted nothing)"}]
    return {"schema_version": SCHEMA_VERSION, "space_id": f"xenium-typing-{arm}-broad8-v1",
            "producer_lane": f"xenium-typing-{arm}", "labels": labels}


def build_donor(arm: str, donor: str) -> tuple:
    calls = gate.read_calls(arm, donor).set_index("cell_id")
    gpath = pins.gate_json(arm, donor)
    if not gpath.exists():
        raise FileNotFoundError(f"{arm} donor {donor}: no gate verdict at {gpath}; run gate.judge first")
    verdict = json.loads(gpath.read_text())
    verified = set(verdict.get("verified_lineages", []))
    bundle = pd.read_parquet(pins.cells_parquet(donor))
    ids = bundle["cell_id"].astype(str).to_numpy()
    out = pd.DataFrame({"cell_id": ids})
    lab = calls["lineage_native"].reindex(ids)
    present = lab.notna().to_numpy()
    out["lineage_native"] = lab.fillna(gate.UNAVAILABLE).astype(str).to_numpy()
    out["celltype_fine"] = calls["celltype_fine"].reindex(ids).astype(object).where(present, None).to_numpy()
    out["confidence"] = calls["confidence"].reindex(ids).astype(float).to_numpy()
    out["reason"] = calls["reason"].reindex(ids).fillna("").astype(str).to_numpy()
    out["cycling"] = pd.array(calls["cycling"].reindex(ids).to_numpy(), dtype="boolean")
    qc = qc_dropped_reasons(donor, ids)
    dropped = ~present
    out.loc[dropped, "reason"] = np.where(qc[dropped].to_numpy() != "", "qc_dropped:" + qc[dropped].to_numpy(),
                                          "not_in_base")
    out.loc[dropped, "confidence"] = np.nan
    lowsig = out["reason"].eq("lowsignal").to_numpy()
    is_lin = out["lineage_native"].isin(pins.LINEAGES).to_numpy()
    out["common_lineage"] = out["lineage_native"].map(XENIUM_TO_COMMON).fillna("").to_numpy()
    out["evidence_positive"] = (is_lin & out["lineage_native"].isin(verified).to_numpy() & ~lowsig)
    unv = is_lin & ~out["lineage_native"].isin(verified).to_numpy() & (out["reason"] == "")
    out.loc[unv, "reason"] = "label_not_verified"
    for c in ("control_probe_counts", "total_counts"):
        if c not in bundle.columns:
            bundle[c] = np.nan
    out["control_probe_frac"] = (bundle["control_probe_counts"] / bundle["total_counts"].replace(0, np.nan)).to_numpy()
    for m in ("INS", "GCG"):
        p = pins.if_values_path(donor, m)
        col = f"feat__if_{m}"
        if p.exists():
            v = np.load(p)
            out[col] = v if len(v) == len(out) else np.nan
            if len(v) != len(out):
                print(f"[release] {donor}: IF {m} array length {len(v)} != bundle {len(out)}; written NaN", flush=True)
        else:
            out[col] = np.nan
    stats = {"n_cells": int(len(out)), "n_typed": int(present.sum()), "n_qc_dropped": int(dropped.sum()),
             "n_evidence_positive": int(out["evidence_positive"].sum()), "verified_lineages": sorted(verified),
             "labels": out["lineage_native"].value_counts().to_dict(),
             "calls_sha256": _sha(pins.arm_dir(arm, donor) / "calls.parquet"), "gate_json": str(gpath), "gate_sha256": _sha(gpath)}
    return out, stats


def write(arm: str, donors: list, out_dir: pathlib.Path) -> pathlib.Path:
    out_dir = pathlib.Path(out_dir)
    calls_dir = out_dir / "calls"
    per_donor = {}
    for d in donors:
        df, stats = build_donor(arm, d)
        assert df["evidence_positive"].isna().sum() == 0
        p = calls_dir / f"donor_id={d}" / f"roi={ROI}"
        p.mkdir(parents=True, exist_ok=True)
        df.to_parquet(p / "data.parquet", index=False)
        per_donor[d] = stats
    space = label_space(arm)
    (out_dir / "label_space.json").write_text(json.dumps(space, indent=1) + "\n")
    inputs = [{"kind": "base", "path": str(pins.base_path(d)), "sha256": _sha(pins.base_path(d)), "donor_id": d} for d in donors]
    inputs.append({"kind": "panel_snapshot", **pins.panel_snapshot()})
    rel = {"schema_version": SCHEMA_VERSION, "release_id": f"xenium-typing-{arm}-dev-{time.strftime('%Y%m%d', time.gmtime())}",
           "modality": "xenium", "producer": {"lane": space["producer_lane"], **pins.git_head()},
           "label_space": "label_space.json", "calls_dir": "calls", "cell_key": "cell_id", "roi": ROI,
           "donors": [str(d) for d in donors], "inputs": inputs,
           "columns": {"label": "lineage_native", "common_lineage": "common_lineage", "evidence_positive": "evidence_positive",
                       "confidence": "confidence", "fine_type": "celltype_fine"},
           "stamp": "development", "latch": "no production typing is declared; this is a development release",
           "arm": {"name": arm, "per_donor": per_donor}, "written_utc": pins.utc_now()}
    text = json.dumps(rel, indent=1) + "\n"
    low = (text + json.dumps(space)).lower()
    for tok in RETIRED_INPUT_TOKENS:
        if tok.lower() in low:
            raise RuntimeError(f"release provenance names a retired input token {tok!r}")
    (out_dir / "release.json").write_text(text)
    return out_dir


def verify(out_dir: pathlib.Path) -> None:
    from phenocycler.config import load_config
    from phenocycler.integration.phenospace import load_release
    cfg = load_config(pins.CONFIG_INI)
    object.__setattr__(cfg, "xenium_identity_release", str(out_dir)) if hasattr(cfg, "__dataclass_fields__") else setattr(cfg, "xenium_identity_release", str(out_dir))
    r = load_release(cfg, "xenium")
    print(f"[release] loads: {r.release_id} donors={r.donors} labels={len(r.space.rules)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True)
    ap.add_argument("--donors", nargs="+", default=["6436", "6442"])
    ap.add_argument("--out", default=str(pins.OUT / "release"))
    a = ap.parse_args(argv)
    out = write(a.arm, a.donors, pathlib.Path(a.out))
    verify(out)
    print(f"[release] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
