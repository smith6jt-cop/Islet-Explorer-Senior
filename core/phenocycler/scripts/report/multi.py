"""Sequential-slide multi-platform figures M01-M06 (PhenoCycler <-> Xenium serial sections, plus the
same-section post-Xenium IF).  The two registrations are never conflated: the serial-section verdict
comes from integration/qc/registration_qc.csv; the same-section IF<->Xenium registration is printed in
its own labelled box.  Stages the verdict does not permit are stated, never faked.
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

sys.path.insert(0, str(pins.PKG))
from phenocycler.config import load_config  # noqa: E402

sys.path.insert(0, str(pins.PKG / "scripts" / "xenium_typing"))
import pins as xpins  # noqa: E402

SEC = "multi"
RESTAIN = pins.REPO / "data" / "development" / "restain"
PAIR_RADIUS_UM = 6.0
FOCAL_INSULIN_CAVEAT = {"6436": "6436 (T1D): insulin is focal -- most islets are insulin-deficient; INS associations are driven by few cells"}


def cfg():
    return load_config(pins.PKG / "config.ini")


def _rec(fid, donor, files, inputs, facts, section=SEC):
    return FigureRecord(figure_id=fid, donor=str(donor), section=section, files=[str(f) for f in files],
                        inputs=[str(p) for p in inputs], facts=facts)


def _side(out, stem, facts):
    p = out / f"{stem}.json"; p.write_text(json.dumps(facts, indent=1, default=str) + "\n"); return p


def _out(donor, section=SEC):
    return pins.OUT / ("_pair" if section == "_pair" else str(donor)) / SEC


def _unavailable(fid, donor, prov, why, section=SEC):
    fig, axes = style.figure(1, 1, width=7, height=3)
    axes[0, 0].axis("off"); style.note(axes[0, 0], f"{fid}: not available\n{why}")
    style.title(fig, f"{fid} -- donor {donor}: not available"); style.footer(fig, prov, None if section == "_pair" else donor); style.finish(fig)
    files = style.savefig(fig, _out(donor, section), f"{fid}_not_available"); plt.close(fig)
    return _rec(fid, donor, files, [], {"not_available": why}, section=section)


def verdict_row(donor: str) -> dict | None:
    q = io.registration_qc()
    if q is None:
        return None
    r = q[(q["donor_id"].astype(str) == str(donor)) & (q["roi"].astype(str) == "panc")]
    return r.iloc[0].to_dict() if len(r) else None


def verdict_text(donor: str) -> str:
    r = verdict_row(donor)
    if r is None:
        return "serial-section registration QC: not run (integration 'qc' stage has no row for this donor)"
    def f(v, fmt="{:.3f}"):
        try:
            return fmt.format(float(v)) if pd.notna(v) else "n/a"
        except Exception:  # noqa: BLE001
            return str(v)
    return (f"SERIAL-SECTION registration (PhenoCycler <-> Xenium): verdict {r.get('verdict')}\n"
            f"tissue Dice {f(r.get('tissue_dice'))} (bar 0.80)   islet RMSE {f(r.get('islet_rmse_um'), '{:.1f}')} um (bar 150)   islet Dice {f(r.get('islet_dice'))}\n"
            f"reasons: {r.get('reasons')}\nallows: {r.get('allows')}\nflags: {r.get('sanity_flags')}")


def same_section_text(donor: str) -> str:
    p = RESTAIN / "valis" / f"donor_id={donor}" / "refined.json"
    q = io.restain_registration(donor)
    lines = ["SAME-SECTION registration (post-Xenium IF <-> Xenium, restain lane) -- a different registration, a different section pair"]
    if p.exists():
        r = json.loads(p.read_text())
        t = r.get("translation_um") or [np.nan, np.nan]
        lines.append(f"affine: scale {r.get('scale'):.4f}, rotation {r.get('rotation_deg'):.3f} deg, translation ({t[0]:.1f}, {t[1]:.1f}) um, mirrored {r.get('mirrored')}, displacement field {r.get('has_displacement')}")
    if q:
        lines.append(f"refine winner: {q.get('winner')} (nuclei ICP; detection ceiling {q.get('detection_ceiling')})")
        ho = q.get("held_out") if isinstance(q.get("held_out"), dict) else {}
        for f in ho.get("folds", []):
            h = f.get("held_out", {})
            lines.append(f"  fitted on {f.get('fitted_on')} half (n={f.get('n_fit'):,}), held-out residual p50 {h.get('p50', float('nan')):.2f} um, p90 {h.get('p90', float('nan')):.2f} um")
        if ho:
            lines.append(f"  half-to-half disagreement: translation {ho.get('translation_disagreement_um', float('nan')):.2f} um, scale {ho.get('scale_disagreement_ppm', float('nan')):.0f} ppm")
    if len(lines) == 1:
        lines.append("no record on disk")
    return "\n".join(lines)


def _common_calls(modality: str, donor: str) -> pd.DataFrame | None:
    c = io.release_calls(modality, donor)
    if c is None:
        return None
    key = "object_id" if modality == "phenocycler" else "cell_id"
    common = io.common_of(modality)
    c = c.rename(columns={key: "cell_id"})
    c["common"] = c["lineage_native"].map(common).fillna(c["lineage_native"].where(c["lineage_native"].isin([io.AMBIGUOUS, io.UNAVAILABLE]), io.AMBIGUOUS))
    return c


def _xy(modality: str, donor: str) -> pd.DataFrame:
    if modality == "phenocycler":
        t = pd.read_parquet(pins.cells_parquet(donor), columns=["object_id", "X_centroid", "Y_centroid"]).rename(columns={"object_id": "cell_id", "X_centroid": "x_um", "Y_centroid": "y_um"})
    else:
        t = pd.read_parquet(xpins.cells_parquet(donor), columns=["cell_id", "x_centroid", "y_centroid"]).rename(columns={"x_centroid": "x_um", "y_centroid": "y_um"})
    t["cell_id"] = t["cell_id"].astype(str)
    return t


def _raster(ax, df, order, colors, px=20.0, title=""):
    gx = (df["x_um"].to_numpy() / px).astype(int); gy = (df["y_um"].to_numpy() / px).astype(int)
    gx -= gx.min(); gy -= gy.min()
    W, H = gx.max() + 1, gy.max() + 1
    code = pd.Categorical(df["common"], categories=order).codes
    counts = np.zeros((H, W, len(order)), dtype=np.int32); np.add.at(counts, (gy, gx, code), 1)
    filled = counts.sum(axis=2) > 0; r = np.full((H, W), -1); r[filled] = counts.argmax(axis=2)[filled]
    rgb = np.ones((H, W, 3))
    for j, L in enumerate(order):
        rgb[r == j] = matplotlib.colors.to_rgb(colors.get(L, "#BBBBBB") if L not in (io.AMBIGUOUS, io.UNAVAILABLE) else ("#888888" if L == io.AMBIGUOUS else "#BBBBBB"))
    ax.imshow(rgb, origin="upper", interpolation="nearest"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(style.wrap(title), fontsize=8)
    style.scale_bar(ax, 1.0 / px, 1000.0) if hasattr(style, "scale_bar") else None


# ------------------------------------------------------------------ M01 registration overlay + the two registrations
def m01_registration(donor: str, prov: dict) -> FigureRecord:
    from phenocycler.integration.figures import registration_overlay
    out = _out(donor); out.mkdir(parents=True, exist_ok=True)
    inner = out / "M01_registration_overlay_inner.png"
    try:
        p = registration_overlay(cfg(), donor, out_path=inner)
    except Exception as exc:  # noqa: BLE001
        p = None; err = f"{type(exc).__name__}: {exc}"
    else:
        err = None
    fig = plt.figure(figsize=(16, 5.5)); gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1])
    axes = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]])
    ax = axes[0, 0]
    if p is not None and pathlib.Path(p).exists():
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(style.wrap("M01-A  integration.figures.registration_overlay (before / after / islet correspondence)"), fontsize=8)
    else:
        ax.axis("off"); style.note(ax, "M01-A not available: registration_overlay returned nothing\n(no transform.json for this donor: 'register' stage not run or not permitted)" + (f"\n{err}" if err else ""))
    ax = axes[0, 1]; ax.axis("off")
    ax.text(0.02, 0.98, style.wrap(verdict_text(donor), 96), va="top", ha="left", fontsize=7.5, family="monospace", transform=ax.transAxes, bbox=dict(boxstyle="round", fc="#FFF4E5", ec="#E69F00"))
    ax.text(0.02, 0.50, style.wrap(same_section_text(donor), 96), va="top", ha="left", fontsize=7.5, family="monospace", transform=ax.transAxes, bbox=dict(boxstyle="round", fc="#E8F1FF", ec="#4477AA"))
    style.title(fig, f"M01 registration and verdict -- donor {donor}"); style.footer(fig, prov, donor); style.finish(fig)
    files = style.savefig(fig, out, "M01_registration"); plt.close(fig)
    facts = {"serial_section": verdict_row(donor), "same_section_record": same_section_text(donor), "overlay_png": str(p) if p else None}
    files.append(_side(out, "M01_registration", facts))
    inputs = [x for x in (io.INTEGRATION / "qc" / "registration_qc.csv", RESTAIN / "valis" / f"donor_id={donor}" / "refined.json", RESTAIN / "refine_panc_0.5um.json") if x.exists()]
    return _rec("M01", donor, files, inputs, facts)


# ------------------------------------------------------------------ M02 side-by-side maps, own frames
def m02_maps(donor: str, prov: dict) -> FigureRecord:
    ph, xe = _common_calls("phenocycler", donor), _common_calls("xenium", donor)
    if ph is None and xe is None:
        return _unavailable("M02", donor, prov, "neither release exists yet")
    colors = io.common_colors()
    order = [L for L in colors if L != "Unassigned"] + [io.AMBIGUOUS, io.UNAVAILABLE]
    r = verdict_row(donor); verdict = (r or {}).get("verdict", "not run")
    fig, axes = style.figure(1, 2, width=7, height=6.5)
    counts = {}
    for ax, (mod, calls) in zip(axes[0], (("phenocycler", ph), ("xenium", xe))):
        if calls is None:
            ax.axis("off"); style.note(ax, f"{mod}: no release yet"); continue
        df = _xy(mod, donor).merge(calls[["cell_id", "common"]], on="cell_id", how="left"); df["common"] = df["common"].fillna(io.UNAVAILABLE)
        _raster(ax, df, order, colors, title=f"M02  {mod} in its OWN frame ({len(df):,} cells, 20 um majority raster of common_lineage)")
        counts[mod] = df["common"].value_counts().to_dict()
    fig.legend(handles=[matplotlib.patches.Patch(color=(colors.get(L) if L in colors else ("#888888" if L == io.AMBIGUOUS else "#BBBBBB")), label=L, hatch="///" if L == io.AMBIGUOUS else None) for L in order if any(L in c for c in counts.values())], loc="lower center", ncol=6, fontsize=7)
    style.title(fig, f"M02 common-lineage maps, each modality in its own frame (serial-section verdict {verdict}: the transform is drawn only in M01) -- donor {donor}")
    style.footer(fig, prov, donor); style.finish(fig, bottom=0.12)
    out = _out(donor); files = style.savefig(fig, out, "M02_side_by_side_maps"); plt.close(fig)
    facts = {"counts": counts, "verdict": verdict}
    files.append(_side(out, "M02_side_by_side_maps", facts))
    return _rec("M02", donor, files, io.release_inputs("phenocycler") + io.release_inputs("xenium"), facts)


# ------------------------------------------------------------------ M03 composition concordance (pair)
def m03_composition(donors: list, prov: dict) -> FigureRecord:
    comp = io.comparable_lineages()
    rows = []
    for d in donors:
        ph, xe = _common_calls("phenocycler", d), _common_calls("xenium", d)
        if ph is None or xe is None:
            continue
        for mode, mask_fn in (("evidence_positive", lambda c: c["evidence_positive"].astype(bool)), ("all_typed", lambda c: ~c["common"].isin([io.AMBIGUOUS, io.UNAVAILABLE]))):
            a, b = ph[mask_fn(ph)], xe[mask_fn(xe)]
            a, b = a[a["common"].isin(comp)], b[b["common"].isin(comp)]
            for L in comp:
                rows.append({"donor": d, "mode": mode, "lineage": L, "pheno": (a["common"] == L).sum() / max(len(a), 1), "xen": (b["common"] == L).sum() / max(len(b), 1), "n_pheno": int(len(a)), "n_xen": int(len(b))})
    if not rows:
        return _unavailable("M03", "_pair", prov, "both releases are needed for at least one donor", section="_pair")
    T = pd.DataFrame(rows)
    from scipy.stats import spearmanr
    rng = np.random.default_rng(0)
    fig, axes = style.figure(1, 2, width=6, height=5.5)
    facts = {}
    colors = io.common_colors()
    for ax, mode in zip(axes[0], ("evidence_positive", "all_typed")):
        t = T[T["mode"] == mode]
        for d, mk in zip(donors, ("o", "s")):
            tt = t[t["donor"] == d]
            ax.scatter(tt["pheno"], tt["xen"], s=50, marker=mk, c=[colors.get(L, "#999") for L in tt["lineage"]], edgecolors="black", linewidths=0.5, label=f"{d} ({style.status_display(pins.donor_status(d))})")
        rho = spearmanr(t["pheno"], t["xen"]).statistic if len(t) > 2 else np.nan
        null = []
        for _ in range(1000):
            perm = np.concatenate([rng.permutation(t.loc[t["donor"] == d, "xen"].to_numpy()) for d in donors if (t["donor"] == d).any()])
            null.append(spearmanr(t["pheno"], perm).statistic)
        p95 = float(np.nanquantile(null, 0.95))
        lim = max(float(t[["pheno", "xen"]].max().max()) * 1.1, 0.05)
        ax.plot([0, lim], [0, lim], color="#999999", lw=0.8, ls=":"); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("PhenoCycler fraction"); ax.set_ylabel("Xenium fraction")
        ax.set_title(style.wrap(f"M03  {mode.replace('_', ' ')}: Spearman rho {rho:.2f} vs within-donor label-permutation null p95 {p95:.2f} (n={len(t)} points, descriptive)"), fontsize=7.5)
        ax.legend(fontsize=6, loc="lower right")
        facts[mode] = {"rho": float(rho), "null_p95": p95, "points": t.to_dict(orient="records")}
    fig.legend(handles=[matplotlib.patches.Patch(color=colors.get(L, "#999"), label=L) for L in comp], loc="lower center", ncol=7, fontsize=7)
    style.title(fig, "M03 composition concordance over donor.comparable_for('pancreas') lineages, both donors"); style.footer(fig, prov); style.finish(fig, bottom=0.14)
    out = _out("_pair", "_pair"); files = style.savefig(fig, out, "M03_composition_concordance"); plt.close(fig)
    files.append(_side(out, "M03_composition_concordance", facts))
    return _rec("M03", "_pair", files, io.release_inputs("phenocycler") + io.release_inputs("xenium"), facts, section="_pair")


# ------------------------------------------------------------------ M04 islet-level
def m04_islets(donor: str, prov: dict) -> FigureRecord:
    r = verdict_row(donor); verdict = (r or {}).get("verdict", "not run")
    paired = io.paired_islets()
    paired = paired[paired["donor_id"].astype(str) == str(donor)] if paired is not None and "donor_id" in paired.columns else paired
    fig, axes = style.figure(1, 2, width=6.5, height=5)
    ax = axes[0, 0]
    facts = {"verdict": verdict}
    if paired is not None and len(paired):
        from phenocycler.integration.figures import islet_concordance_figure
        out = _out(donor); out.mkdir(parents=True, exist_ok=True)
        p = islet_concordance_figure(cfg(), out_path=out / "M04_islet_concordance_inner.png")
        if p:
            ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(style.wrap("M04-A  islet_concordance_figure (matched islets)"), fontsize=8)
        facts["n_pairs"] = int(len(paired))
    else:
        ax.axis("off"); style.note(ax, f"M04-A  structure matching not permitted at verdict {verdict}\n(paired/islets.parquet absent; ALLOWED_BY_VERDICT[FAIL] = donor, niche)")
    ax = axes[0, 1]
    drew = False
    for mod, col in (("phenocycler", "#4477AA"), ("xenium", "#EE6677")):
        s = io.structures(mod, donor)
        if s is None or not len(s):
            continue
        size_col = next((c for c in ("n_cells", "n_members", "size") if c in s.columns), None)
        if size_col is None:
            continue
        v = s[size_col].to_numpy(float)
        style.hist_from_zero(ax, v, bins=40, color=col, alpha=0.5, label=f"{mod}: {len(s):,} islets, median {np.median(v):.0f} cells") if "color" in style.hist_from_zero.__code__.co_varnames else ax.hist(v, bins=40, color=col, alpha=0.5, label=f"{mod}: {len(s):,} islets, median {np.median(v):.0f} cells")
        facts[mod] = {"n_islets": int(len(s)), "median_size": float(np.median(v))}
        drew = True
    if drew:
        ax.set_xlabel("islet size (cells, identity-derived S3 structures)"); ax.set_ylabel("islets"); ax.legend(fontsize=7)
        ax.set_title(style.wrap("M04-B  registration-free islet count and size per modality"), fontsize=8)
    else:
        ax.axis("off"); style.note(ax, "M04-B not available: no identity-derived structures on disk yet ('structures' stage not run)")
    style.title(fig, f"M04 islet level -- donor {donor} (serial-section verdict {verdict})"); style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "M04_islets"); plt.close(fig)
    files.append(_side(out, "M04_islets", facts))
    return _rec("M04", donor, files, [x for x in (io.INTEGRATION / "qc" / "registration_qc.csv",) if x.exists()], facts)


# ------------------------------------------------------------------ M05 same-section IF vs RNA
def m05_same_section(donor: str, prov: dict) -> FigureRecord:
    rec = io.protein_rna(donor)
    xe = io.release_calls("xenium", donor)
    if rec is None:
        return _unavailable("M05", donor, prov, "restain lane protein_rna_panc.json has no row for this donor")
    fig, axes = style.figure(1, 3, width=5.2, height=5.5)
    facts = {"caveat": FOCAL_INSULIN_CAVEAT.get(str(donor), "")}
    # A. Spearman rho per gene, INS and GCG
    ax = axes[0, 0]
    mk = rec["compartments"]["cell"]["markers"]
    genes = list(mk["INS"]["genes"]); x = np.arange(len(genes))
    for k, (m, col) in enumerate((("INS", "#CCBB44"), ("GCG", "#EE6677"))):
        rho = [mk[m]["genes"][g]["spearman_rho"] for g in genes]
        ax.bar(x + (k - 0.5) * 0.4, rho, 0.4, color=col, label=f"IF {m} vs transcript count (n={mk[m]['n_cells_measured']:,} cells)")
    hint = [mk["INS"]["genes"][g].get("subtype_hint", "") for g in genes]
    ax.set_xticks(x); ax.set_xticklabels([f"{g}\n{h}" for g, h in zip(genes, hint)], rotation=90, fontsize=6)
    ax.axhline(0, color="black", lw=0.6); ax.set_ylabel("Spearman rho per cell"); ax.legend(fontsize=6, loc="lower left")
    ax.set_title(style.wrap(f"M05-A  label-free association, {len(genes)}-gene endocrine identity core (same section; no cell is called INS+/GCG+)"), fontsize=7.5)
    facts["rho"] = {m: {g: mk[m]["genes"][g]["spearman_rho"] for g in genes} for m in ("INS", "GCG")}
    # B. IF hormone per Xenium label beside a label-permutation null
    ax = axes[0, 1]
    if xe is None:
        ax.axis("off"); style.note(ax, "M05-B not available: no Xenium release yet")
    else:
        common, colors = io.common_of("xenium"), io.common_colors()
        cells = pd.read_parquet(xpins.cells_parquet(donor), columns=["cell_id"]); cells["cell_id"] = cells["cell_id"].astype(str)
        lab = xe.set_index("cell_id")["lineage_native"].reindex(cells["cell_id"]).fillna(io.UNAVAILABLE).to_numpy()
        labs = [L for L in list(common) + [io.AMBIGUOUS] if (lab == L).sum() >= 50]
        rng = np.random.default_rng(0)
        for k, (m, col) in enumerate((("INS", "#CCBB44"), ("GCG", "#EE6677"))):
            v = np.load(xpins.if_values_path(donor, m)).astype(float)
            ok = np.isfinite(v)
            med = [np.median(v[ok & (lab == L)]) for L in labs]
            perm_meds = []
            for _ in range(200):
                pl = lab[ok][rng.permutation(ok.sum())]
                perm_meds.append([np.median(v[ok][pl == L]) for L in labs])
            pm = np.asarray(perm_meds); lo, hi = np.quantile(pm, 0.025, axis=0), np.quantile(pm, 0.975, axis=0)
            xx = np.arange(len(labs)) + (k - 0.5) * 0.4
            ax.bar(xx, med, 0.4, color=col, label=f"IF {m} median per label")
            ax.errorbar(xx, (lo + hi) / 2, yerr=(hi - lo) / 2, fmt="none", ecolor="black", elinewidth=0.8, capsize=2, label="label-permutation null 95% band" if k == 0 else None)
            facts[f"if_{m}_median_by_label"] = dict(zip(labs, map(float, med)))
        ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=7); ax.set_ylabel("IF intensity in the Xenium cell polygon (from 0)"); ax.set_ylim(bottom=0)
        ax.legend(fontsize=6); ax.set_title(style.wrap("M05-B  same-section IF hormone per Xenium released label, beside its null (Q4: no positive call)"), fontsize=7.5)
    # C. Xenium-labelled IF nucleus overlay
    ax = axes[0, 2]
    vd = RESTAIN / "valis" / f"donor_id={donor}"
    need = [vd / "if_nuclei_um.npy", vd / "warped_if_um.npy", vd / "xenium_cell_ids.npy"]
    if xe is None or not all(p.exists() for p in need):
        ax.axis("off"); style.note(ax, "M05-C not available: needs the Xenium release and the restain lane's IF nuclei + warped Xenium centroids")
    else:
        from scipy.spatial import cKDTree
        if_xy = np.load(need[0]); wx = np.load(need[1]); cid = np.load(need[2], allow_pickle=True).astype(str)
        tree = cKDTree(wx); d, j = tree.query(if_xy, distance_upper_bound=PAIR_RADIUS_UM)
        paired = np.isfinite(d)
        lab_by_cell = xe.set_index("cell_id")["lineage_native"]
        nl = np.full(len(if_xy), "unpaired", dtype=object); nl[paired] = lab_by_cell.reindex(cid[j[paired]]).fillna(io.UNAVAILABLE).to_numpy()
        common, colors = io.common_of("xenium"), io.common_colors()
        endo = [L for L in common if common[L] == "Endocrine"]
        cx, cy = (np.median(if_xy[np.isin(nl, endo)], axis=0) if np.isin(nl, endo).any() else np.median(if_xy, axis=0))
        half = 400.0
        w = (np.abs(if_xy[:, 0] - cx) < half) & (np.abs(if_xy[:, 1] - cy) < half)
        cols = [("#DDDDDD" if L == "unpaired" else io.label_color(L, common, colors)) for L in nl[w]]
        ax.scatter(if_xy[w, 0], if_xy[w, 1], s=3, c=cols, linewidths=0); ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xlabel("IF x (um)"); ax.set_ylabel("IF y (um)")
        ax.set_title(style.wrap(f"M05-C  IF nuclei ({w.sum():,} in an 800 um window) coloured by the Xenium label of the nearest warped Xenium cell within {PAIR_RADIUS_UM:.0f} um; {paired.mean():.0%} of {len(if_xy):,} IF nuclei paired; no protein call is made"), fontsize=7)
        facts["overlay"] = {"n_if_nuclei": int(len(if_xy)), "frac_paired": float(paired.mean()), "pair_radius_um": PAIR_RADIUS_UM, "labels_in_window": pd.Series(nl[w]).value_counts().to_dict()}
    style.title(fig, f"M05 same-section IF <-> Xenium -- donor {donor}. {FOCAL_INSULIN_CAVEAT.get(str(donor), '')}")
    style.footer(fig, prov, donor); style.finish(fig, top=0.88)
    out = _out(donor); files = style.savefig(fig, out, "M05_same_section_if_rna"); plt.close(fig)
    files.append(_side(out, "M05_same_section_if_rna", facts))
    inputs = [RESTAIN / "protein_rna_panc.json"] + [p for p in need if p.exists()] + io.release_inputs("xenium")
    return _rec("M05", donor, files, inputs, facts)


# ------------------------------------------------------------------ M06 niche-view concordance
def m06_niches(donor: str, prov: dict) -> FigureRecord:
    prof = io.niche_profiles()
    if prof is None or not len(prof):
        return _unavailable("M06", donor, prov, "registration-free niche profiles not on disk yet ('grid' stage not run)")
    dcol = next((c for c in ("donor_id", "donor") if c in prof.columns), None)
    t = prof[prof[dcol].astype(str) == str(donor)] if dcol else prof
    if not len(t):
        return _unavailable("M06", donor, prov, "niche_profiles.csv has no row for this donor")
    mcol = next((c for c in ("modality",) if c in t.columns), None)
    lin_cols = [c for c in t.columns if c in io.common_colors()]
    ncol = next((c for c in ("niche", "niche_id") if c in t.columns), None)
    if mcol is None or ncol is None or not lin_cols:
        return _unavailable("M06", donor, prov, f"niche_profiles.csv columns not understood: {list(t.columns)}")
    P = t[t[mcol] == "phenocycler"].set_index(ncol)[lin_cols]; X = t[t[mcol] == "xenium"].set_index(ncol)[lin_cols]
    common_n = [n for n in P.index if n in X.index]
    fig, axes = style.figure(1, 2, width=6.5, height=5)
    ax = axes[0, 0]
    colors = io.common_colors()
    io.stacked_bars(ax, pd.concat([P.loc[common_n].add_prefix("").rename(index=lambda n: f"P:{n}"), X.loc[common_n].rename(index=lambda n: f"X:{n}")]), lin_cols, {}, colors, hatch_ambiguous=False)
    ax.set_title(style.wrap("M06-A  niche composition in the common vocabulary, both modalities (registration-free)"), fontsize=8); ax.tick_params(axis="x", rotation=90)
    ax = axes[0, 1]
    if common_n:
        from scipy.stats import spearmanr
        real = np.mean([spearmanr(P.loc[n], X.loc[n]).statistic for n in common_n])
        rng = np.random.default_rng(0); null = []
        for _ in range(1000):
            perm = rng.permutation(len(lin_cols)); null.append(np.mean([spearmanr(P.loc[n].to_numpy(), X.loc[n].to_numpy()[perm]).statistic for n in common_n]))
        ax.hist(null, bins=30, color="#BBBBBB", label="lineage-permutation null"); ax.axvline(real, color="#d62728", label=f"observed mean rho {real:.2f}")
        ax.axvline(np.quantile(null, 0.95), color="black", ls=":", label=f"null p95 {np.quantile(null, 0.95):.2f}"); ax.legend(fontsize=7); ax.set_xlabel("mean per-niche Spearman rho")
        ax.set_title(style.wrap("M06-B  per-niche profile agreement vs null"), fontsize=8)
        facts = {"mean_rho": float(real), "null_p95": float(np.quantile(null, 0.95)), "n_niches": len(common_n)}
    else:
        ax.axis("off"); style.note(ax, "no niche shared by both modalities"); facts = {"n_niches": 0}
    style.title(fig, f"M06 niche-view concordance -- donor {donor}"); style.footer(fig, prov, donor); style.finish(fig)
    out = _out(donor); files = style.savefig(fig, out, "M06_niche_concordance"); plt.close(fig)
    files.append(_side(out, "M06_niche_concordance", facts))
    return _rec("M06", donor, files, [io.INTEGRATION / "niches" / "niche_profiles.csv"], facts)


PER_DONOR = {"M01": m01_registration, "M02": m02_maps, "M04": m04_islets, "M05": m05_same_section, "M06": m06_niches}
PAIR = {"M03": m03_composition}
