"""Who is inside the islets, given that QuPath already told us where they are.

The islets are a **prior, not a finding**.  Every image in this cohort was run through a
QuPath pixel classifier on CD99, INS, GCG and SST, and its output is the ``Islet_N`` polygons
and their 20 um expansion ring -- already joined per cell as ``cell_region``.  Rediscovering
them is not analysis: UTAG recovers them at 25x enrichment, which is a good check that the
data is sound and no new information at all.

That prior is strong and it is measurable.  Median within-donor percentile inside ``Islet_N``
minus outside it, over six donors: **CD99 +0.49 to +0.50 in every donor** -- more consistent
than INS, which runs +0.04 to +0.51 because T1D donors have lost their beta cells -- with
GCG +0.49 to +0.50 and CD56 +0.35 to +0.45.  Immune markers go the other way (CD3e -0.07 to
-0.29, CD11b -0.09 to -0.33), which is the point: inside an islet the non-endocrine cells are
the minority, and they are what this study is about.

So the work is not to call 25 million cells.  It is:

* **inside ``Islet_N``** -- take endocrine as the default and move a cell off it only on
  positive evidence of something else.  Islet cores are 1.5-6.2% of cells per donor.
* **in the 20 um ring** -- no default.  This is the immune environment around islets, the
  thing the study exists to measure, so every reference competes.

Endocrine SUBTYPE is a separate question and is answered here directly from the hormones,
because inside the prior we already know the cell is endocrine and the only question is which.

    python scripts/celltype/run_islet.py                 # whole cohort
    python scripts/celltype/run_islet.py --donors 6450
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402
import rules as R                                # noqa: E402
import score as SC                               # noqa: E402
import run_step3 as S3                           # noqa: E402

OUT = panel.OUT / "islet_residents"

#: References that describe an endocrine cell.  Inside the prior these are not scored -- the
#: pixel classifier already made that call, and re-deriving it is what produced Endo-beta at
#: 0.00% while the same cells sat inside polygons drawn from INS.
ENDOCRINE_PREFIX = ("Epi-endocrine", "Endo-alpha", "Endo-beta", "Endo-delta",
                    "Endo-hormone-negative")


def intruder_references(rules: pd.DataFrame, columns) -> list:
    """References that can legitimately displace the endocrine prior inside an islet.

    The question inside the prior is immune / vascular / lymphatic / other -- not what kind
    of epithelium.  So the whole epithelial branch is out of contention there, ancestors and
    siblings alike, because an endocrine cell satisfies those rules by construction: it is
    E-cadherin+ and carries Ker8/18 and Pan-CK. Left in, the generic ``Epithelial`` parent
    took 42% of core cells off the prior and ``Epi-exocrine`` another 26%.

    That is a scope decision, not a claim that no duct ever crosses an islet outline.  A
    Ker8/18-high cell inside a polygon drawn by an INS/GCG/SST/CD99 pixel classifier is a
    question about the outline, and this lane does not answer it -- the ring, where there is
    no prior, scores every reference including the epithelial ones.
    """
    par = dict(zip(rules["phenotype"], rules["parent"]))

    def branch(ph):
        cur, seen = ph, []
        while cur and cur != "all":
            seen.append(cur)
            cur = par.get(cur)
            if cur in seen:
                break
        return seen

    endo_branch = set()
    for ph in rules["phenotype"]:
        if str(ph).startswith(ENDOCRINE_PREFIX):
            endo_branch |= set(branch(ph))
    implied = {c for c in rules["phenotype"] if endo_branch & set(branch(c))}
    return [c for c in columns if c not in implied]

#: Hormones, for subtyping inside the prior.  Read directly: within an islet the question is
#: not "is this endocrine" but "which hormone", and that is what these markers answer.
HORMONE = {"Endo-beta": "INS", "Endo-alpha": "GCG", "Endo-delta": "SST"}


def native_hormone_resolver(donor, aligned, scale, gate):
    """Recover the applied alignment by cell ID and use the fitted floor and span."""
    def resolve(contenders):
        from multiplex5000 import roots
        from multiplex5000.artifacts import QuPathCohortManifest
        from multiplex5000.tissues.interior_dominance import measure_conflicts

        manifest = pathlib.Path(os.environ.get("IO60_QUPATH_MANIFEST")
                                or (roots.source("qupath_manifest") if hasattr(roots, "source") else None)
                                or roots.INTEGRATION.parent / "qupath_manifest.json")
        if not manifest.is_file():
            return {}
        images = [i for i in QuPathCohortManifest.read_json(manifest).images
                  if panel._same_donor(i.donor_id, donor)]
        if len(images) != 1:
            return {}
        markers = sorted({m for pair in contenders.values() for m in pair})
        selected = aligned.set_index("object_id", verify_integrity=True)
        selected.index = selected.index.astype(str)
        selected = selected.loc[list(contenders), markers]
        raw = pq.ParquetFile(panel.expression_path(donor)).read(
            columns=["object_id", *markers]).to_pandas().set_index("object_id", verify_integrity=True)
        raw.index = raw.index.astype(str)
        raw = raw.reindex(selected.index)
        calibration = {}
        for marker in markers:
            x, y = raw[marker].to_numpy(float), selected[marker].to_numpy(float)
            keep = np.isfinite(x) & np.isfinite(y) & (x > 0)
            if keep.any():
                factor = np.median(y[keep] / x[keep])
                if not np.allclose(x[keep] * factor, y[keep], rtol=1e-5, atol=1e-7):
                    raise ValueError(f"{donor} {marker}: aligned expression differs from its source scaling")
                calibration[marker] = (factor, scale.loc[marker, "floor"], scale.loc[marker, "span"])
        results = measure_conflicts(images[0], contenders, calibration, gate)
        for oid, result in results.items():
            result.update(image_content_id=images[0].content_id, contenders=contenders[oid], gate=gate)
        return results
    return resolve


def subtype_endocrine(Xs: np.ndarray, markers: list, min_scaled: float,
                      object_ids=None, resolve=None) -> pd.DataFrame:
    """Resolve two positive co-leaders using native interior area and median brightness."""
    n = len(Xs)
    names = [name for name, marker in HORMONE.items() if marker in markers]
    out = pd.DataFrame({"call": ["Endo-unresolved"] * n, "score": np.nan, "runner_up": "",
                        "margin": np.nan, "hormone_status": "missing_hormone_measurement",
                        "hormone_interior_evidence": None})
    if not names or not n:
        return out
    values = Xs[:, [markers.index(HORMONE[name]) for name in names]]
    scores = np.where(np.isfinite(values), values, -np.inf)
    order = np.argsort(-scores, axis=1, kind="stable")
    rows = np.arange(n)
    first = order[:, 0]
    best = scores[rows, first]
    finite = np.isfinite(best)
    positive = finite & (best >= min_scaled)
    ties = (scores == best[:, None]) & positive[:, None]
    unique = positive & (ties.sum(axis=1) == 1)
    out.loc[finite, "score"] = best[finite]
    out.loc[finite & ~positive, ["call", "hormone_status"]] = ["Endo-hormone-negative", "below_gate"]
    out.loc[unique, "call"] = np.asarray(names)[first[unique]]
    out.loc[unique, "hormone_status"] = "hormone_argmax"
    if len(names) > 1:
        second = order[:, 1]
        runner = unique & np.isfinite(scores[rows, second])
        out.loc[runner, "runner_up"] = np.asarray(names)[second[runner]]
        out.loc[runner, "margin"] = best[runner] - scores[rows[runner], second[runner]]
    else:
        out.loc[unique, "margin"] = best[unique]
    tied = positive & ~unique
    out.loc[tied, "margin"] = 0.0
    out.loc[tied, "hormone_status"] = "multiple_positive_leaders"
    pairs = np.flatnonzero(ties.sum(axis=1) == 2)
    out.loc[pairs, "hormone_status"] = "missing_native_evidence"
    if len(pairs) and resolve is not None and object_ids is not None:
        ids = np.asarray(object_ids, dtype=str)
        if len(ids) != n or len(set(ids[pairs])) != len(pairs):
            raise ValueError("Hormone conflicts require unique matching source object IDs")
        contenders = {ids[i]: tuple(HORMONE[names[j]] for j in np.flatnonzero(ties[i])) for i in pairs}
        measured = resolve(contenders)
        labels = {v: k for k, v in HORMONE.items()}
        for i in pairs:
            evidence = measured.get(ids[i])
            if evidence is None:
                continue
            out.at[i, "hormone_interior_evidence"] = json.dumps(evidence, sort_keys=True, allow_nan=False)
            winner = evidence["winner"]
            out.at[i, "hormone_status"] = evidence["status"]
            if winner is not None:
                a, b = contenders[ids[i]]
                if winner not in (a, b):
                    raise ValueError("Native hormone winner was not a contender")
                out.at[i, "call"] = labels[winner]
                out.at[i, "runner_up"] = labels[b if winner == a else a]
                out.at[i, "hormone_status"] = "hormone_interior_dominance"
    return out


def run_donor(donor: str, Rm: pd.DataFrame, scale: pd.DataFrame, reqs: dict,
              require_min: float, ratio: bool, intruders: list) -> pd.DataFrame:
    missing = [m for m in Rm.index if m not in set(panel.donor_markers(donor))]
    usable = [p for p in Rm.columns if not (Rm.loc[missing, p] != 0).any()]
    intruder = [p for p in usable if p in intruders]

    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    df, have, _ = S3._read_donor(donor, list(Rm.index), keep, ratio=False)
    df = df[df["qc_analysis_eligible"] & df["cell_region"].isin(("core", "peri"))]
    if df.empty:
        return pd.DataFrame()
    df = df.reset_index(drop=True)
    Xs = SC.apply_scale(df[have].to_numpy("float32"), have, scale)
    # The DAPI ratio is a per-reference decision -- see run_step3.NUCLEAR_RATIO.
    nuc_all = (SC.nuclear_only_references(Rm, panel.NUCLEUS_MEASURED)
               if ratio else np.zeros(len(Rm.columns), dtype=bool))
    Xs_nuc = None
    if nuc_all.any():
        dfr, _ = SC.nuclear_ratio(df.copy(), have, "DAPI", panel.NUCLEUS_MEASURED)
        Xs_nuc = SC.apply_scale(dfr[have].to_numpy("float32"), have,
                                S3._nuclear_scale(scale, [donor], have))
    nuc_by_ref = dict(zip(Rm.columns, nuc_all))
    core = (df["cell_region"] == "core").to_numpy()

    # The two regions are different questions and get different reference sets.  Inside
    # Islet_N the prior settles the epithelial one, so only intruders compete.  The ring has no
    # prior -- it is ordinary tissue, mostly acinar -- so every reference does.  Scoring the
    # ring against intruders alone left 49% of it Ambiguous: the acinar and ductal cells had
    # nowhere to go.
    def _call(rows, refs):
        if not rows.any():
            return None
        Xr = Xs[rows]
        Xr_nuc = Xs_nuc[rows] if Xs_nuc is not None else None
        nuc_ref = np.array([nuc_by_ref.get(r, False) for r in refs], dtype=bool)
        S, allowed = SC.score_and_gate(Xr, Xr_nuc, Rm.loc[have, refs].to_numpy("float32"),
                                       have, reqs, refs, require_min, nuc_ref)
        o = SC.call_from_scores(S, refs, allowed=allowed)
        o.index = np.flatnonzero(rows)
        return o

    parts = [x for x in (_call(core, intruder), _call(~core, usable)) if x is not None]
    out = pd.concat(parts).sort_index().reset_index(drop=True)

    # Inside Islet_N, a cell that shows no intruder is endocrine BY THE PRIOR, and is
    # subtyped by hormone rather than left Ambiguous.
    fell_through = core & (out["call"] == "Ambiguous").to_numpy()
    rows = np.flatnonzero(fell_through)
    sub = subtype_endocrine(Xs[rows], have, require_min, df.loc[rows, "object_id"],
                            native_hormone_resolver(donor, df, scale, require_min))
    for column in ("call", "runner_up", "score", "margin"):
        out.loc[rows, column] = sub[column].to_numpy()
    out.loc[rows, "reason"] = "endocrine_by_islet_prior"
    for column in ("hormone_status", "hormone_interior_evidence"):
        out.loc[rows, column] = sub[column].to_numpy()

    out.insert(0, "object_id", df["object_id"].to_numpy())
    out.insert(0, "donor_id", donor)
    out["cell_region"] = df["cell_region"].to_numpy()
    out["islet_num"] = df["islet_num"].to_numpy()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--require-min", type=float, default=S3.REQUIRE_MIN)
    ap.add_argument("--no-nuclear-ratio", action="store_true")
    args = ap.parse_args()

    donors = args.donors or panel.cohort()["donor_id"].tolist()
    flat = R.flatten(R.load_rules())
    Rm, _ = R.reference_matrix(flat, panel.markers(include_variant=True), S3.ENCODING,
                               strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    reqs = R.requirements(R.load_rules())
    ratio = not args.no_nuclear_ratio
    intruders = intruder_references(R.load_rules(), Rm.columns)
    excluded = [c for c in Rm.columns
                if c not in intruders and not str(c).startswith(ENDOCRINE_PREFIX)]
    print(f"Islet residents — {len(intruders)} references may displace the endocrine prior")
    print(f"  excluded as implied by it: {excluded}")

    scale = S3.build_scale(donors, list(Rm.index), ratio=ratio)
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for d in donors:
        t = time.time()
        c = run_donor(d, Rm, scale, reqs, args.require_min, ratio, intruders)
        if c.empty:
            print(f"  {d}: no islet cells"); continue
        p = OUT / f"donor_id={d}"; p.mkdir(parents=True, exist_ok=True)
        c.to_parquet(p / "data.parquet", index=False)
        frames.append(c)
        print(f"  {d}: {len(c):,} islet-associated cells "
              f"({(c.cell_region=='core').sum():,} core, {(c.cell_region=='peri').sum():,} "
              f"peri)   {time.time()-t:.0f}s", flush=True)
    calls = pd.concat(frames, ignore_index=True)

    for region in ("core", "peri"):
        g = calls[calls.cell_region == region]
        print(f"\n--- inside Islet_N ({region}), {len(g):,} cells ---"
              if region == "core" else
              f"\n--- the 20 um ring ({region}), {len(g):,} cells — the study's question ---")
        v = g["call"].value_counts(normalize=True)
        for c, f in v.head(14).items():
            print(f"    {c:<38}{f:7.2%}")
    calls.to_parquet(panel.OUT / "islet_residents_all.parquet", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
