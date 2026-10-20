#!/usr/bin/env python3
"""What the same-section linkage supports, per donor, measured rather than asserted.

Two lanes claim a cell-level protein/RNA correspondence in this cohort, and they are not the
same claim:

**The post-Xenium IF re-stain.** The Xenium section is stained for INS and GCG afterwards, so
the IF image and the transcriptome are the SAME PHYSICAL SECTION. ``scripts/restain/protein_rna.py``
measures the two hormones inside every Xenium cell's OWN polygon, so there is no segmentation
transfer at all on that axis -- the uncertainty is registration, not fan-out. What there IS a
fan-out on is the other direction: the IF image has its own nuclei, and asking "which IF nucleus
is this Xenium cell" is many-to-one in both directions. That is measured here.

**The same-slide lane (6533, 112).** Full PhenoCycler then Xenium on one section, with the
protein CALL carried onto Xenium cells by largest shared area
(``core/multiplex5000/scripts/sameslide/lift_calls.py``, read not run here). That transfer is
many-to-one by construction and the lifted call is a protein statement about overlapping area,
not a measurement of the Xenium cell.

Everything this script prints is computed from files under the repository. Nothing is carried
across from prose. Where an input is absent the donor is reported as unsupported, not skipped.

Outputs (under ``data/development/integration/same_section/``):
  donor_linkage.csv    one row per donor: registration state, what the linkage supports
  fanout.csv           the many-to-one fan-out, both directions, per donor
  agreement.csv        per-cell protein-vs-RNA agreement where it is defensible
  linkage.json         the same numbers, with the run's provenance

    .venv/bin/python scripts/senior/same_section_linkage.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
VALIS = ROOT / "data/development/restain/valis"
INTEG = ROOT / "core/phenocycler/data/integration"
CELLS_XEN = INTEG / "cells_xen"
PROTEIN_RNA = ROOT / "data/development/restain/protein_rna"
#: The held-out gate for the two 0.5 um donors was written to a shared file rather than into
#: their donor directories, so a donor-directory-only reader silently finds no gate for them.
REFINE_0_5UM = ROOT / "data/development/restain/refine_panc_0.5um.json"
OUT = ROOT / "data/development/integration/same_section"

#: Held-out median nucleus offset above which the refined transform does not support a
#: cell-level statement. A nucleus is ~7 um across, so an offset at that scale means the
#: nearest nucleus after warping is routinely a DIFFERENT nucleus. This is a geometric
#: statement about the transform, not a tuned threshold and not fitted to any result.
PASS_P50_UM = 2.0

#: Cap for "these two nuclei are the same nucleus". Same reasoning: half a nucleus diameter.
NN_CAP_UM = 3.5


def donor_exclusions() -> set[str]:
    sys.path.insert(0, str(ROOT / "core/multiplex5000"))
    from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402
    return {str(d) for d in DONOR_EXCLUSIONS}


def held_out_p50(refinement: dict) -> tuple[float | None, float | None]:
    """Median and p90 held-out nucleus offset, pooled over the two spatial folds."""
    folds = refinement.get("held_out", {}).get("folds", [])
    p50 = [f["held_out"]["p50"] for f in folds if "held_out" in f]
    p90 = [f["held_out"]["p90"] for f in folds if "held_out" in f]
    if not p50:
        return None, None
    return float(np.mean(p50)), float(np.mean(p90))


def shared_gates() -> dict:
    """Gates recorded outside the donor directories, keyed by donor."""
    j = read_json(REFINE_0_5UM)
    return {str(r["donor_id"]): r for r in (j or [])}


def read_json(p: Path):
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def nn(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For every row of ``a``, the nearest row of ``b``: (distance, index)."""
    from scipy.spatial import cKDTree
    d, i = cKDTree(b).query(a, k=1, workers=-1)
    return d, i


def inventory(excluded: set[str]) -> pd.DataFrame:
    shared = shared_gates()
    sections = pd.read_csv(INTEG / "restain_sections.csv", comment="#")
    channels = {str(r.donor_id): r for r in sections.itertuples()}
    modes = pd.read_csv(INTEG / "section_modes.csv", comment="#")
    same_slide = {str(r.donor_id): r for r in modes.itertuples() if r.mode == "same_slide"}

    rows = []
    donors = sorted({p.name.split("=")[1] for p in VALIS.glob("donor_id=*")}
                    | set(same_slide))
    for d in donors:
        dd = VALIS / f"donor_id={d}"
        refined = read_json(dd / "refined.json")
        refinement = read_json(dd / "refinement.json") or shared.get(d)
        gate_source = ("refinement.json" if (dd / "refinement.json").exists()
                       else ("refine_panc_0.5um.json" if d in shared else None))
        p50 = p90 = None
        if refinement:
            p50, p90 = held_out_p50(refinement)
        xen = CELLS_XEN / f"donor_id={d}/roi=panc"
        xen_files = sorted(xen.glob("*.parquet"))
        if_cols = []
        n_xen_cells = None
        if xen_files:
            import pyarrow.parquet as pq
            pf = pq.ParquetFile(xen_files[0])
            names = pf.schema_arrow.names
            if_cols = [n for n in names if n.startswith("feat__if_")]
            n_xen_cells = pf.metadata.num_rows
        rows.append(dict(
            donor_id=d,
            centrally_excluded=d in excluded,
            same_slide=d in same_slide,
            has_if_nuclei=(dd / "if_nuclei_um.npy").exists(),
            has_xenium_nuclei=(dd / "xenium_nuclei_um.npy").exists(),
            has_refined_transform=refined is not None,
            has_refined_gate=refinement is not None,
            gate_source=gate_source,
            transform_winner=(refined or {}).get("meta", {}).get("candidate"),
            rotation_deg=(refined or {}).get("rotation_deg"),
            heldout_p50_um=p50,
            heldout_p90_um=p90,
            detection_ceiling=(refinement or {}).get("detection_ceiling"),
            n_if_nuclei=(refinement or {}).get("n_if_nuclei"),
            n_xenium_nuclei=(refinement or {}).get("n_xenium_nuclei"),
            n_xenium_cells=n_xen_cells,
            channel_map_recorded=d in channels,
            if_values_on_disk=bool(list(PROTEIN_RNA.glob(f"{d}_*_cell_INS.npy"))),
            if_stamped_on_cells=sorted(if_cols),
        ))
    return pd.DataFrame(rows)


def verdict(r: pd.Series) -> tuple[str, str]:
    """What this donor's linkage supports, and why."""
    if r.centrally_excluded:
        return "nothing", "centrally excluded from the cohort"
    if r.same_slide and not r.has_refined_transform:
        return "nothing", ("same-slide acquisition, but no refined transform in this checkout: "
                           "the lifted call cannot be reproduced or checked here")
    if not r.has_refined_transform:
        return "nothing", "no refined transform: registration did not produce one"
    if r.heldout_p50_um is None or not np.isfinite(r.heldout_p50_um):
        return "per-region", ("a transform exists but no held-out gate was run, so its "
                              "cell-level accuracy is unmeasured here")
    if r.heldout_p50_um > PASS_P50_UM:
        return "per-region", (f"held-out median nucleus offset {r.heldout_p50_um:.3f} um exceeds "
                              f"{PASS_P50_UM} um; the nearest warped nucleus is routinely a "
                              f"different nucleus")
    if not r.if_stamped_on_cells:
        return "per-region", ("registration passes, but no protein value is stamped on this "
                              "donor's Xenium cells, so there is nothing to join per cell")
    if not r.channel_map_recorded:
        return "per-region", ("no recorded channel-to-antibody mapping for this section: the "
                              "IF values cannot be named")
    return "per-cell", (f"held-out median nucleus offset {r.heldout_p50_um:.3f} um, protein "
                        f"stamped per Xenium polygon, channel mapping recorded")


def fanout(donors: list[str]) -> pd.DataFrame:
    rows = []
    for d in donors:
        dd = VALIS / f"donor_id={d}"
        fi, xn = dd / "if_nuclei_um.npy", dd / "xenium_nuclei_um.npy"
        if not (fi.exists() and xn.exists() and (dd / "refined.json").exists()):
            rows.append(dict(donor_id=d, measurable=False,
                             reason="no nuclei, or no refined transform to place them with"))
            continue
        # NUCLEI against NUCLEI. warped_if_um.npy is the valis-stage warp of Xenium CELL
        # centroids -- both wrong here: a cell centroid is displaced from its own nucleus,
        # and the valis warp is superseded by the refined one (measured on 6442: mean 1.37 um,
        # max 119 um apart). Warp the Xenium nuclei with the refined transform instead.
        sys.path.insert(0, str(ROOT / "core/multiplex5000"))
        from phenocycler.integration.transform import Transform  # noqa: E402
        ifn = np.load(fi)
        xnuc = Transform.load(dd, stem="refined").apply(np.load(xn))
        xids = np.load(dd / "xenium_nuclei_cell_ids.npy", allow_pickle=True)
        if len(xnuc) != len(xids):
            rows.append(dict(donor_id=d, measurable=False,
                             reason=f"nuclei {len(xnuc)} != nucleus cell ids {len(xids)}"))
            continue
        xcell = xnuc

        # Xenium nucleus -> nearest IF nucleus, and the reverse. Both are capped at NN_CAP_UM:
        # beyond that the pair is not a correspondence, it is the absence of one.
        d_xi, i_xi = nn(xnuc, ifn)
        d_ix, i_ix = nn(ifn, xnuc)
        x_ok, i_ok = d_xi <= NN_CAP_UM, d_ix <= NN_CAP_UM

        # Fan-in: how many Xenium cells claim the SAME IF nucleus.
        claims = np.bincount(i_xi[x_ok], minlength=len(ifn))
        per_x = claims[i_xi[x_ok]]
        # Fan-out: how many IF nuclei land inside the catchment of one Xenium cell.
        rev = np.bincount(i_ix[i_ok], minlength=len(xcell))

        rows.append(dict(
            donor_id=d, measurable=True, reason="",
            n_if_nuclei=int(len(ifn)), n_xenium_nuclei=int(len(xcell)),
            seg_ratio=float(len(ifn) / len(xcell)),
            xen_nuc_with_if_partner=float(x_ok.mean()),
            if_nuc_with_xen_partner=float(i_ok.mean()),
            nn_p50_um=float(np.median(d_xi)), nn_p90_um=float(np.percentile(d_xi, 90)),
            # one-to-one is per_x == 1; anything above is a shared IF nucleus
            share_one_to_one=float((per_x == 1).mean()),
            fanin_p50=float(np.median(per_x)), fanin_p95=float(np.percentile(per_x, 95)),
            fanin_max=int(per_x.max()) if per_x.size else 0,
            fanout_p50=float(np.median(rev[rev > 0])) if (rev > 0).any() else 0.0,
            fanout_p95=float(np.percentile(rev[rev > 0], 95)) if (rev > 0).any() else 0.0,
            fanout_max=int(rev.max()),
        ))
    return pd.DataFrame(rows)


def auroc(score: np.ndarray, pos: np.ndarray) -> float:
    """Rank AUROC. Threshold-free on purpose: the IF values are features, and naming a
    cut-off on them would be asserting positivity that nobody has defined."""
    r = pd.Series(score).rank().to_numpy()
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def nucleus_match_class(d: str) -> pd.Series:
    """Per Xenium cell: is its nucleus matched one-to-one to an IF nucleus?

    The IF values themselves are measured in the Xenium cell's OWN polygon, so they carry no
    transfer error. This classification is not about them -- it is the segmentation ambiguity
    a LIFTED CALL would have to pass through, applied to the same cells, so the agreement can
    be reported with and without it and the difference read directly.
    """
    dd = VALIS / f"donor_id={d}"
    sys.path.insert(0, str(ROOT / "core/multiplex5000"))
    from phenocycler.integration.transform import Transform  # noqa: E402
    ifn = np.load(dd / "if_nuclei_um.npy")
    xnuc = Transform.load(dd, stem="refined").apply(np.load(dd / "xenium_nuclei_um.npy"))
    xids = np.load(dd / "xenium_nuclei_cell_ids.npy", allow_pickle=True).astype(str)
    dist, idx = nn(xnuc, ifn)
    ok = dist <= NN_CAP_UM
    claims = np.bincount(idx[ok], minlength=len(ifn))
    cls = np.full(len(xids), "unmatched", dtype=object)
    cls[ok] = np.where(claims[idx[ok]] == 1, "one_to_one", "shared_if_nucleus")
    # A Xenium cell can carry more than one nucleus (measured: 958,918 nuclei over 958,744
    # cells on 6442). Keep the nearest nucleus for the cell rather than dropping either.
    t = pd.DataFrame({"cell_id": xids, "dist": dist, "cls": cls})
    t = t.sort_values("dist").drop_duplicates("cell_id", keep="first")
    return t.set_index("cell_id")["cls"]


def agreement(donors: list[str]) -> pd.DataFrame:
    """Per-cell protein-vs-RNA agreement, restricted to cells that carry a nucleus.

    The RNA call is the declared pancreas Xenium release (docs/PRODUCERS.md, producer 2),
    resolved by ``scripts/senior/xenium_release.resolve`` and joined to the ingest on
    ``cell_id``; the IF values come from the ingest, which is the only tree that carries
    them. Until 2026-09-23 the call was read off the ingest's own ``lineage_native``, the
    superseded scANVI release that config.ini:153 attached to it. A cell is *called* by the
    release module's own reading of its label (``classify_label``): ``Uncalled_*`` and
    ``Unavailable*`` are not calls and leave the comparison.

    The ``endocrine_fine`` strata are no longer produced. They split the endocrine call by
    a fine type the old arm inferred from surrogates; the declared release carries no fine
    endocrine type, because the 5K panel has no INS, GCG or SST. With no release column to
    read, the stratum is absent rather than re-derived.
    """
    sys.path.insert(0, str(ROOT))
    from scripts.senior import xenium_release  # noqa: E402
    release = xenium_release.resolve("pancreas")
    rows = []
    if not release.available:
        # Named, not silent: the frame carries the resolver's reason for every donor asked.
        return pd.DataFrame([dict(donor_id=d, stratum="unavailable",
                                  reason=f"no_current_xenium_pancreas_release ({release.reason})")
                             for d in donors])
    label = release.require_column("label")
    for d in donors:
        f = sorted((CELLS_XEN / f"donor_id={d}/roi=panc").glob("*.parquet"))
        if not f:
            continue
        cols = ["cell_id", "feat__if_INS", "feat__if_GCG"]
        df = pd.read_parquet(f[0], columns=cols)
        df["cell_id"] = df["cell_id"].astype(str)
        calls = release.read_calls(columns=["cell_id", label], donors=[d]).to_pandas()
        if calls.empty:
            rows.append(dict(donor_id=d, stratum="unavailable",
                             reason=f"{release.release_id}: no partition for this donor"))
            continue
        calls["cell_id"] = calls["cell_id"].astype(str)
        calls = calls.drop_duplicates("cell_id").set_index("cell_id")[label]
        df["call"] = calls.reindex(df["cell_id"]).to_numpy()
        dd = VALIS / f"donor_id={d}"
        nuc_ids = set(np.load(dd / "xenium_nuclei_cell_ids.npy", allow_pickle=True).tolist())
        has_nuc = df["cell_id"].isin(nuc_ids)
        n_all = len(df)
        df = df[has_nuc]                                   # no nucleus, no calculation
        imaged = df["feat__if_INS"].notna() & df["feat__if_GCG"].notna()
        sub = df[imaged]
        called = sub["call"].map(xenium_release.classify_label) == xenium_release.CALL
        sub = sub[called]
        endo = (sub["call"] == "Endocrine").to_numpy()
        mcls = nucleus_match_class(d).reindex(sub["cell_id"]).fillna("unmatched")
        mcls = mcls.to_numpy()
        base = dict(donor_id=d, release_id=release.release_id, call_column=label,
                    n_cells=n_all, n_with_nucleus=int(has_nuc.sum()),
                    n_nucleus_free_excluded=int(n_all - has_nuc.sum()),
                    n_imaged=int(imaged.sum()),
                    n_outside_imaged_area=int((~imaged).sum()),
                    n_called=int(len(sub)), n_endocrine=int(endo.sum()))
        for name, score in (("if_INS", sub["feat__if_INS"].to_numpy()),
                            ("if_GCG", sub["feat__if_GCG"].to_numpy())):
            base[f"auroc_{name}_endocrine"] = auroc(score, endo)
            base[f"median_{name}_endocrine"] = float(np.median(score[endo])) if endo.any() else np.nan
            base[f"median_{name}_other"] = float(np.median(score[~endo])) if (~endo).any() else np.nan
        rows.append(dict(**base, stratum="all", n_stratum=int(len(sub))))
        # The same agreement, split by the segmentation ambiguity. A gap between these rows
        # is the cost a many-to-one transfer would impose; no gap means the disagreement is
        # not a segmentation artefact.
        for cname in ("one_to_one", "shared_if_nucleus", "unmatched"):
            m = mcls == cname
            if m.sum() < 500:
                continue
            r = dict(base)
            r.update(stratum=f"nucleus_match:{cname}", n_stratum=int(m.sum()),
                     n_endocrine=int(endo[m].sum()))
            for name, col in (("if_INS", "feat__if_INS"), ("if_GCG", "feat__if_GCG")):
                sc = sub[col].to_numpy()[m]
                r[f"auroc_{name}_endocrine"] = auroc(sc, endo[m])
                r[f"median_{name}_endocrine"] = (float(np.median(sc[endo[m]]))
                                                 if endo[m].any() else np.nan)
                r[f"median_{name}_other"] = (float(np.median(sc[~endo[m]]))
                                             if (~endo[m]).any() else np.nan)
            rows.append(r)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    excluded = donor_exclusions()
    inv = inventory(excluded)
    v = inv.apply(verdict, axis=1, result_type="expand")
    inv["supports"], inv["because"] = v[0], v[1]
    inv.to_csv(args.out / "donor_linkage.csv", index=False)
    print(inv[["donor_id", "has_refined_transform", "heldout_p50_um", "supports"]]
          .to_string(index=False))

    fo = fanout(sorted(inv.donor_id[~inv.centrally_excluded]))
    fo.to_csv(args.out / "fanout.csv", index=False)
    print("\n== fan-out ==")
    print(fo.to_string(index=False))

    per_cell = sorted(inv.donor_id[inv.supports == "per-cell"])
    ag = agreement(per_cell)
    ag.to_csv(args.out / "agreement.csv", index=False)
    print("\n== agreement ==")
    print(ag.to_string(index=False))

    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = None
    (args.out / "linkage.json").write_text(json.dumps({
        "commit": commit,
        "pass_p50_um": PASS_P50_UM, "nn_cap_um": NN_CAP_UM,
        "donor_exclusions": sorted(excluded),
        "per_cell_donors": per_cell,
        "donor_linkage": json.loads(inv.to_json(orient="records")),
        "fanout": json.loads(fo.to_json(orient="records")),
        "agreement": json.loads(ag.to_json(orient="records")),
    }, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
