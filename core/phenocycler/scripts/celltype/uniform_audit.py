"""Step 8e: two audits that must accompany every UniFORM refit.

UniFORM (``run_step2.py``) is the only normalization of record for the v3 re-run.  It aligns
each marker's background mode per donor with one rigid log-shift, and ``run_step3.build_scale``
then puts the aligned values on one cohort floor/top.  Two assumptions sit under that, and
neither is checked by the background-mode SD bar ``run_step2`` reports:

1. **Batch is gone** -- after alignment one cohort-wide level means the same thing on every
   donor.  Checked by the RESTORE residual-batch test (Chang et al., RESTORE), section (a).
2. **The background is invariant** -- the population whose mode is aligned is biologically
   the same across donors.  Where disease changes it (HLA class I rising in T1D exocrine
   tissue, interferon-induced markers), aligning the mode subtracts the biology.  Checked by
   the marker-invariance audit, section (b).

Both write tables and a markdown summary.  **Neither changes UniFORM**, refuses a shift, or
feeds a threshold anywhere: they are evidence for a person deciding whether a marker's shift
may stand, and nothing reads their output automatically.

(a) RESTORE residual-batch test
-------------------------------
For a marker ``i`` and a biologically exclusive partner ``j1``, cells positive for ``j1`` are
in-tissue negative controls for ``i``.  Per donor, ``i``'s threshold is

    t = max(mean + k_sigma * sd, q-quantile) * margin        (defaults 3, 0.95, 1.10)

of ``i`` over that donor's ``j1``-positive cells, and the pooled threshold is the same statistic
over every donor's ``j1``-positive cells together.  If batch effects are gone the two agree, so
``log(t_donor / t_pooled)`` and ``frac_donor - frac_pooled`` (the positive fraction of ``i`` under
each) carry no batch structure.  They are summarised per stain batch; a spread of batch medians
of the log ratio above ``tol_log_ratio`` is reported as residual batch.

``j1``-positivity is **reference-based**, never a frequency: a cell is ``j1``-positive when ``j1``
is strictly above a level that is either declared (``--partner-levels``) or fitted per donor as
the p99 of ``j1`` over that donor's reference population (``--reference-column``, a boolean
column -- the same p99-of-ordinary-tissue convention as ``typing/subpanel_io.reference_level``).
With neither, the audit refuses to run.  The ``j1``-positive set is decided once, on the aligned
values, and keyed by ``(donor_id, object_id)``; a per-donor shift is monotone, so the set is the
same before and after alignment and only ``i``'s location moves between the two runs.

Co-positivity is scored on a **different** exclusive partner ``j2``: the fraction of
``i``-positive cells that are also ``j2``-positive.  On ``j1`` itself it would be capped near 5 %
by construction (the threshold is at least the 95th percentile of ``i`` in ``j1``-positive
cells), so ``j2 == j1`` -- and ``j2 == i`` -- raise :class:`CircularPartnerError`.  A marker with
only one declared partner has no co-positivity row; nothing falls back to ``j1``.

The exclusivity ratio ``r = sigma2 / sigma1`` is the ratio of singular values of the
uncentred two-column matrix of the pair's values per donor, on **cohort-scale axes**
(``score.apply_scale``: clipped to [0, 1] at the cohort floor and top).  Mutually exclusive
populations lie along the two axes and give ``r`` near 1 when balanced; co-expressed markers
lie along the diagonal and give ``r`` near 0.  It is computed from the 2x2 Gram matrix, which
has the same singular values, so it is O(cells).  ``r`` depends on the balance of the two
populations as well as their exclusivity, so it is read across donors for one pair, not as a
pass/fail number on its own.

(b) Marker-invariance audit
---------------------------
From the shifts JSON ``run_step2`` writes (``shifts[marker]['shift_log'][donor]``; the applied
scale is ``exp(-shift_log)``), per marker fit ``shift_log ~ batch + disease`` by OLS and test
the disease term with an F statistic whose null distribution comes from permuting disease
labels **within batch**.  A batch holding one disease level contributes no permutations; if no
batch holds two, disease is fully confounded with batch, the term is not estimable, and the
marker is reported ``confounded`` -- which is *not* a pass.  A batch with no ND donor is named
in the design table whatever the marker.  Markers are flagged on the Benjamini-Hochberg q of
the permutation p across the markers tested.  Markers refused by ``run_step2``
(``refused_reagent_change``, ``too_few_donors_to_align``) have no shift to test and are listed
as such.

Excluded donors are derived from ``phenocycler.cohort.DONOR_EXCLUSIONS`` and dropped from both
audits; the dropped keys are recorded in the summary.

Inputs are QC-filtered on ``qc_analysis_eligible`` (a cell with no nucleus is never
analysis-eligible), which must be present.  A missing measurement is NaN and leaves every
denominator it would have entered; a donor without a marker gets an ``unavailable`` row, not a
zero.

Running it after a v3 UniFORM refit
-----------------------------------
From ``core/phenocycler/``, after ``run_step2.py --out <aligned> --shifts-out <shifts.json>``::

    python scripts/celltype/uniform_audit.py \\
        --expression <raw expression root> --apply-shifts \\
        --shifts <shifts.json> \\
        --donor-metadata <domain donor_metadata.csv> \\
        --acquisition-runs ../../data/donor_acquisition_runs.csv \\
        --reference-column <boolean reference-population column> \\
        --out ../../data/development/celltype/uniform_audit_v3

``--apply-shifts`` reads the pre-alignment values and runs the RESTORE test twice, as given
(``pre_alignment``) and multiplied by each donor's UniFORM scale (``aligned``): a batch spread
present before and gone after is the shift working, one present after is residual batch.
Without it, ``--expression`` is taken as already aligned (``run_step2``'s ``--out`` tree) and
the test runs once.  ``--pairs`` takes a CSV (``marker_a,marker_b``) or
``phenocycler/fundamental_pair_diagnostics.json``; the default is :data:`DEFAULT_PAIRS`.
``--scale`` takes a cohort floor/top parquet (``run_step3.build_scale``'s shape); without it one
is derived from the aligned values the same way ``build_scale`` does.  Read ``SUMMARY.md``
first; every table beside it is tidy and keyed on ``(marker, partner, donor_id)`` or ``marker``.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score as SC                                  # noqa: E402

from phenocycler.cohort import DONOR_EXCLUSIONS, filter_eligible_donors   # noqa: E402

#: Declared exclusive pairs, used when ``--pairs`` is not given.  Each is a pair of markers no
#: single cell should carry together; ``j1`` and ``j2`` for a marker are drawn from its partners
#: here.  An input, not a finding: edit the table, not the code, to change it.
DEFAULT_PAIRS: tuple[tuple[str, str], ...] = (
    ("INS", "GCG"), ("INS", "SST"), ("GCG", "SST"),
    ("CD3e", "CD20"), ("CD3e", "CD79a"), ("CD68", "CD3e"),
    ("Ker8_18", "CD45"),
)

#: Markers whose background is expected a priori to move with disease.  Printed first in the
#: summary so a reader looks at them; NEVER used to flag -- the flag comes from the test alone.
DISEASE_WATCHLIST: tuple[str, ...] = ("HLA_A", "HLA-A", "HLA_ABC", "HLA-ABC", "HLA_DR",
                                      "HLA-DR", "IDO1", "PD_L1")

#: Same 1 % convention as ``typing/subpanel_io.REFERENCE_QUANTILE``.
REFERENCE_QUANTILE = 99.0
K_SIGMA = 3.0
CONTROL_QUANTILE = 0.95
MARGIN = 1.10
#: Fewest j1-positive control cells a donor must contribute before its threshold is fitted.
MIN_CONTROL_CELLS = 100
#: Fewest cells with both values finite before a donor's exclusivity ratio is computed.
MIN_PAIR_CELLS = 100
#: A spread of per-batch median log(t_donor / t_pooled) above log(1.25) is reported.
TOL_LOG_RATIO = math.log(1.25)
N_PERM = 5000
ALPHA = 0.05

#: The shift provenances ``run_step2`` writes that carry a fitted shift.
SHIFTED = ("automatic_registration", "automatic_registration_over_bar")


class CircularPartnerError(ValueError):
    """A co-positivity partner equal to the control partner or to the marker itself."""


# ------------------------------------------------------------------------------------ inputs
def normalize_donor(d: object) -> str:
    """``0112`` and ``112`` are one donor, as ``run_step2.main`` normalizes them."""
    s = str(d).strip()
    return str(int(s)) if s.isdecimal() else s


def _parquet_files(path: pathlib.Path) -> list[pathlib.Path]:
    path = pathlib.Path(path)
    if path.is_dir():
        files = sorted(path.rglob("*.parquet"))
        if not files:
            raise FileNotFoundError(f"{path}: no parquet files under it")
        return files
    if not path.exists():
        raise FileNotFoundError(path)
    return [path]


def _partition_donor(file: pathlib.Path) -> str | None:
    for part in reversed(file.parts):
        if part.startswith("donor_id="):
            return part.split("=", 1)[1]
    return None


def load_expression(paths: Sequence[pathlib.Path | str], markers: Iterable[str],
                    reference_column: str | None = None) -> pd.DataFrame:
    """Read the QC-eligible cells of every donor, keyed on ``(donor_id, object_id)``.

    ``paths`` are parquet files or hive-partitioned roots (``donor_id=<d>/*.parquet``).  A
    marker a file does not carry is NaN for that file's cells -- unavailable, never zero.
    ``qc_analysis_eligible`` is required: it is what excludes a cell without a nucleus.
    Centrally excluded donors are dropped here and recorded in ``frame.attrs['dropped']``.
    """
    markers = list(dict.fromkeys(markers))
    frames = []
    for file in (f for p in paths for f in _parquet_files(pathlib.Path(p))):
        have = set(pq.read_schema(file).names)
        if "qc_analysis_eligible" not in have:
            raise KeyError(f"{file}: no qc_analysis_eligible column -- refusing to audit cells "
                           "that were never QC-filtered")
        if "object_id" not in have:
            raise KeyError(f"{file}: no object_id column")
        if reference_column is not None and reference_column not in have:
            raise KeyError(f"{file}: reference column {reference_column!r} is absent")
        cols = ["object_id", "qc_analysis_eligible"]
        cols += [c for c in ("donor_id", reference_column) if c and c in have]
        cols += [m for m in markers if m in have]
        t = pq.read_table(file, columns=cols).to_pandas()
        if "donor_id" not in t.columns:
            donor = _partition_donor(file)
            if donor is None:
                raise KeyError(f"{file}: no donor_id column and no donor_id= partition")
            t["donor_id"] = donor
        t = t.loc[t["qc_analysis_eligible"].fillna(False).astype(bool)].drop(
            columns="qc_analysis_eligible")
        for m in markers:
            t[m] = t[m].astype("float64") if m in t.columns else np.nan
        frames.append(t)
    frame = pd.concat(frames, ignore_index=True)
    frame["donor_id"] = frame["donor_id"].map(normalize_donor)
    donors = sorted(frame["donor_id"].unique())
    keep = set(filter_eligible_donors(donors))
    frame = frame.loc[frame["donor_id"].isin(keep)].reset_index(drop=True)
    frame.attrs["dropped"] = sorted(set(donors) - keep)
    if frame.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("duplicate (donor_id, object_id) keys across the expression inputs")
    return frame


def load_shifts(path: pathlib.Path | str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``run_step2``'s shifts JSON -> (tidy shifts, markers without a shift and why).

    The tidy frame has ``marker, donor_id, shift_log, scale, provenance``.  A donor whose
    scale is the string ``unavailable`` is left out: it has no shift, not a shift of 0.
    """
    meta = json.loads(pathlib.Path(path).read_text())
    shifts = meta["shifts"] if "shifts" in meta else meta
    rows, skipped = [], []
    for marker, rec in shifts.items():
        prov = str(rec.get("provenance", ""))
        if prov not in SHIFTED:
            skipped.append({"marker": marker, "provenance": prov})
            continue
        scale = rec.get("scale", {})
        for donor, s in rec.get("shift_log", {}).items():
            if not isinstance(s, (int, float)) or not math.isfinite(s):
                continue
            sc = scale.get(donor, math.exp(-s))
            if not isinstance(sc, (int, float)):
                continue
            rows.append({"marker": marker, "donor_id": normalize_donor(donor),
                         "shift_log": float(s), "scale": float(sc), "provenance": prov})
    tidy = pd.DataFrame(rows, columns=["marker", "donor_id", "shift_log", "scale", "provenance"])
    keep = set(filter_eligible_donors(tidy["donor_id"].unique()))
    tidy = tidy.loc[tidy["donor_id"].isin(keep)].reset_index(drop=True)
    return tidy, pd.DataFrame(skipped, columns=["marker", "provenance"])


def load_donor_metadata(metadata: pathlib.Path | str | pd.DataFrame,
                        runs: pathlib.Path | str | pd.DataFrame | None = None,
                        disease_col: str = "donor_status",
                        batch_col: str = "stain_group") -> pd.DataFrame:
    """``donor_id, disease, batch`` per donor, joining metadata to the acquisition runs.

    ``runs`` has ``data/donor_acquisition_runs.csv``'s shape (``donor_id, acq_run,
    stain_group``); when omitted, ``batch_col`` must be on ``metadata``.  A donor on the
    metadata but on neither table's batch column is kept with ``batch`` missing and reported
    by the audit, not assigned a batch.
    """
    def read(x):
        return x.copy() if isinstance(x, pd.DataFrame) else pd.read_csv(x, dtype=str)
    meta = read(metadata)
    if "donor_id" not in meta.columns:
        for alt in ("case_id", "donor"):
            if alt in meta.columns:
                meta = meta.rename(columns={alt: "donor_id"})
                break
    meta["donor_id"] = meta["donor_id"].map(normalize_donor)
    if runs is not None:
        r = read(runs)
        r["donor_id"] = r["donor_id"].map(normalize_donor)
        meta = meta.drop(columns=[c for c in r.columns if c != "donor_id" and c in meta.columns])
        meta = meta.merge(r.drop_duplicates("donor_id"), on="donor_id", how="left")
    for col in (disease_col, batch_col):
        if col not in meta.columns:
            raise KeyError(f"donor metadata has no {col!r} column")
    out = pd.DataFrame({"donor_id": meta["donor_id"],
                        "disease": meta[disease_col].astype("string").str.strip(),
                        "batch": meta[batch_col].astype("string").str.strip()})
    out = out.loc[out["donor_id"].isin(filter_eligible_donors(out["donor_id"]))]
    return out.drop_duplicates("donor_id").reset_index(drop=True)


def load_pairs(path: pathlib.Path | str | None = None) -> pd.DataFrame:
    """Exclusive pairs as ``marker_a, marker_b``, from a CSV, the fundamental-pair JSON, or
    :data:`DEFAULT_PAIRS`.  Order within a pair carries no meaning."""
    if path is None:
        pairs = list(DEFAULT_PAIRS)
    else:
        path = pathlib.Path(path)
        if path.suffix == ".json":
            doc = json.loads(path.read_text())
            pairs = [(p["target"], p["reference"]) for p in doc["pairs"]]
        else:
            t = pd.read_csv(path)
            pairs = list(zip(t.iloc[:, 0].astype(str), t.iloc[:, 1].astype(str), strict=True))
    out = pd.DataFrame(pairs, columns=["marker_a", "marker_b"])
    if (out["marker_a"] == out["marker_b"]).any():
        raise ValueError("a marker cannot be exclusive with itself")
    key = out.apply(lambda r: tuple(sorted((r.marker_a, r.marker_b))), axis=1)
    return out.loc[~key.duplicated()].reset_index(drop=True)


def partner_triples(pairs: pd.DataFrame) -> pd.DataFrame:
    """Every ``(marker i, control partner j1, co-positivity partner j2)`` the table allows.

    Both orientations of every pair are a ``(i, j1)``; ``j2`` is any other partner of ``i``.
    A ``(i, j1)`` with no other partner has one row with ``j2`` missing.
    """
    partners: dict[str, list[str]] = {}
    for a, b in pairs[["marker_a", "marker_b"]].itertuples(index=False):
        partners.setdefault(a, []).append(b)
        partners.setdefault(b, []).append(a)
    rows = []
    for i, js in partners.items():
        for j1 in js:
            others = [j for j in js if j != j1]
            for j2 in others or [None]:
                check_partners(i, j1, j2)
                rows.append({"marker": i, "partner": j1, "copos_partner": j2})
    return pd.DataFrame(rows, columns=["marker", "partner", "copos_partner"])


def check_partners(marker: str, partner: str, copos_partner: str | None) -> None:
    """Refuse the circular designs: co-positivity on ``j1`` itself, or on the marker."""
    if partner == marker:
        raise CircularPartnerError(f"{marker}: control partner is the marker itself")
    if copos_partner is not None and copos_partner in (partner, marker):
        raise CircularPartnerError(
            f"{marker}: co-positivity partner {copos_partner!r} must differ from the control "
            f"partner {partner!r} and the marker -- on the same pair co-positivity is capped "
            "near 5 % by construction")


# ------------------------------------------------------------------------ shared primitives
def apply_uniform_scale(frame: pd.DataFrame, shifts: pd.DataFrame,
                        markers: Iterable[str]) -> pd.DataFrame:
    """Multiply each ``(donor, marker)`` by its UniFORM scale, as ``run_step2.apply_shift``.

    A pair without a recorded scale is left as given, as ``apply_shift`` leaves a marker a
    donor never acquired or whose shift was refused.
    """
    out = frame.copy()
    lookup = shifts.set_index(["marker", "donor_id"])["scale"]
    for m in markers:
        if m not in out.columns:
            continue
        s = out["donor_id"].map(lambda d, m=m: lookup.get((m, d), 1.0)).to_numpy(float)
        out[m] = out[m].to_numpy(float) * s
    return out


def restore_threshold(values: np.ndarray, k_sigma: float = K_SIGMA,
                      quantile: float = CONTROL_QUANTILE, margin: float = MARGIN,
                      min_cells: int = MIN_CONTROL_CELLS) -> float:
    """RESTORE's negative-control threshold: ``max(mean + k*sd, q) * margin``, or NaN when
    fewer than ``min_cells`` finite control values are on hand."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < min_cells:
        return float("nan")
    return float(max(v.mean() + k_sigma * v.std(ddof=1), np.quantile(v, quantile)) * margin)


def reference_levels(frame: pd.DataFrame, markers: Iterable[str],
                     reference_column: str,
                     quantile: float = REFERENCE_QUANTILE) -> pd.DataFrame:
    """Per donor, the p99 of each marker over that donor's reference population.

    A level at or below 0.0 is refused (NaN): a population sitting at the floor returns the
    floor, which would call every measured cell positive while looking fitted -- the floor
    guard of ``subpanel_io.endocrine_level``.
    """
    rows = []
    for donor, t in frame.groupby("donor_id", sort=True):
        ref = t.loc[t[reference_column].fillna(False).astype(bool)]
        for m in markers:
            v = ref[m].to_numpy(float)
            v = v[np.isfinite(v)]
            level = float(np.percentile(v, quantile)) if v.size else float("nan")
            if not (np.isfinite(level) and level > 0.0):
                level = float("nan")
            rows.append({"donor_id": donor, "marker": m, "level": level,
                         "n_reference": int(v.size)})
    return pd.DataFrame(rows, columns=["donor_id", "marker", "level", "n_reference"])


def positive_masks(frame: pd.DataFrame, levels: pd.DataFrame,
                   markers: Iterable[str]) -> pd.DataFrame:
    """Boolean ``<marker>`` columns, keyed like ``frame``: strictly above the donor's level.

    ``levels`` has ``marker, level`` and optionally ``donor_id`` (a declared cohort level has
    none).  A cell with the value missing, or a donor with no level, is not positive.
    """
    out = pd.DataFrame(index=frame.index)
    per_donor = "donor_id" in levels.columns and levels["donor_id"].notna().any()
    for m in markers:
        lv = levels.loc[levels["marker"] == m]
        if per_donor:
            lvl = frame["donor_id"].map(lv.set_index("donor_id")["level"]).to_numpy(float)
        else:
            lvl = np.full(len(frame), float(lv["level"].iloc[0]) if len(lv) else np.nan)
        v = frame[m].to_numpy(float)
        with np.errstate(invalid="ignore"):
            out[m] = np.isfinite(v) & np.isfinite(lvl) & (v > lvl)
    return out


def exclusivity_ratio(x: np.ndarray, y: np.ndarray, min_cells: int = MIN_PAIR_CELLS) -> float:
    """``sigma2 / sigma1`` of the uncentred ``[x, y]`` matrix, from its 2x2 Gram matrix."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < min_cells:
        return float("nan")
    x, y = x[ok], y[ok]
    gram = np.array([[x @ x, x @ y], [x @ y, y @ y]])
    lam = np.clip(np.linalg.eigvalsh(gram), 0.0, None)
    if lam[1] <= 0:
        return float("nan")
    return float(np.sqrt(lam[0] / lam[1]))


def derive_cohort_scale(frame: pd.DataFrame, markers: Iterable[str]) -> pd.DataFrame:
    """Cohort floor/top per marker on aligned values, as ``run_step3.build_scale`` does: the
    floor is the pooled log1p mode, the top ``score.SCALE_Q``, via ``score.marker_scale``."""
    pool, floors = {}, {}
    for m in markers:
        v = frame[m].to_numpy(float)
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        pool[m] = v
        pos = v[v > 0]
        if pos.size >= 1000:
            h, e = np.histogram(np.log1p(pos), bins=512)
            floors[m] = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
    return SC.marker_scale(pool, floors)


# ----------------------------------------------------------------- (a) RESTORE residual batch
def restore_test(frame: pd.DataFrame, triples: pd.DataFrame, partner_pos: pd.DataFrame,
                 values_label: str = "as_given", k_sigma: float = K_SIGMA,
                 quantile: float = CONTROL_QUANTILE, margin: float = MARGIN,
                 min_cells: int = MIN_CONTROL_CELLS) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-donor vs pooled RESTORE thresholds, and co-positivity on a different partner.

    ``partner_pos`` holds boolean positivity columns (from :func:`positive_masks`) aligned to
    ``frame``'s index, fixed once on the aligned values.  Returns ``(thresholds, copos)``:
    ``thresholds`` per ``(marker, partner, donor_id)``, ``copos`` per ``(marker, partner,
    copos_partner, donor_id)``.
    """
    th_rows, cp_rows = [], []
    donors = sorted(frame["donor_id"].unique())
    for (i, j1), grp in triples.groupby(["marker", "partner"], sort=True):
        for j2 in grp["copos_partner"]:
            check_partners(i, j1, j2 if isinstance(j2, str) else None)
        if i not in frame.columns or j1 not in partner_pos.columns:
            continue
        vi = frame[i].to_numpy(float)
        ctrl = partner_pos[j1].to_numpy(bool)
        t_pool = restore_threshold(vi[ctrl], k_sigma, quantile, margin, min_cells)
        for d in donors:
            at = (frame["donor_id"] == d).to_numpy()
            v = vi[at]
            measured = np.isfinite(v)
            n_ctrl = int((ctrl[at] & measured).sum())
            t_d = restore_threshold(v[ctrl[at]], k_sigma, quantile, margin, min_cells)
            status = ("unavailable_marker" if not measured.any() else
                      "too_few_controls" if not np.isfinite(t_d) else "ok")
            with np.errstate(invalid="ignore"):
                pos_d = measured & (v > t_d)
                pos_p = measured & (v > t_pool)
            n = int(measured.sum())
            frac_d = pos_d.sum() / n if n and np.isfinite(t_d) else np.nan
            frac_p = pos_p.sum() / n if n and np.isfinite(t_pool) else np.nan
            ratio = (math.log(t_d / t_pool) if np.isfinite(t_d) and np.isfinite(t_pool)
                     and t_d > 0 and t_pool > 0 else np.nan)
            th_rows.append({"values": values_label, "marker": i, "partner": j1,
                            "donor_id": d, "status": status, "n_cells": n,
                            "n_controls": n_ctrl, "threshold_donor": t_d,
                            "threshold_pooled": t_pool, "log_threshold_ratio": ratio,
                            "frac_donor": frac_d, "frac_pooled": frac_p,
                            "frac_delta": (frac_d - frac_p) if np.isfinite(frac_d)
                            and np.isfinite(frac_p) else np.nan})
            for j2 in grp["copos_partner"]:
                if not isinstance(j2, str) or j2 not in partner_pos.columns:
                    continue
                p2 = partner_pos[j2].to_numpy(bool)[at]
                cp_rows.append({
                    "values": values_label, "marker": i, "partner": j1, "copos_partner": j2,
                    "donor_id": d, "n_pos_donor": int(pos_d.sum()),
                    "n_pos_pooled": int(pos_p.sum()),
                    "copos_donor": (p2 & pos_d).sum() / pos_d.sum() if pos_d.sum() else np.nan,
                    "copos_pooled": (p2 & pos_p).sum() / pos_p.sum() if pos_p.sum() else np.nan,
                })
    return pd.DataFrame(th_rows), pd.DataFrame(cp_rows)


def restore_batch_summary(thresholds: pd.DataFrame, donors: pd.DataFrame,
                          tol_log_ratio: float = TOL_LOG_RATIO) -> tuple[pd.DataFrame,
                                                                         pd.DataFrame]:
    """Per stain batch, and per ``(values, marker, partner)``, the residual-batch readout.

    ``donors`` has ``donor_id, batch``; a donor without a batch is summarised under
    ``<missing>``, never folded into another batch.  ``residual_batch`` is True when the
    spread of per-batch median log ratios exceeds ``tol_log_ratio``.
    """
    t = thresholds.merge(donors[["donor_id", "batch"]], on="donor_id", how="left")
    t["batch"] = t["batch"].astype("string").fillna("<missing>")
    ok = t.loc[t["status"] == "ok"]
    keys = ["values", "marker", "partner"]
    per_batch = (ok.groupby(keys + ["batch"])
                 .agg(n_donors=("donor_id", "nunique"),
                      median_log_ratio=("log_threshold_ratio", "median"),
                      median_frac_delta=("frac_delta", "median"),
                      median_frac_donor=("frac_donor", "median"),
                      median_frac_pooled=("frac_pooled", "median"))
                 .reset_index())
    by_key = {k: g for k, g in per_batch.groupby(keys, sort=True)}
    rows = []
    for key, g in t.groupby(keys, sort=True):
        b = by_key.get(key, per_batch.iloc[0:0])
        gok = g.loc[g["status"] == "ok"]
        spread = (float(b["median_log_ratio"].max() - b["median_log_ratio"].min())
                  if len(b) >= 2 else np.nan)
        rows.append({**dict(zip(keys, key, strict=True)),
                     "n_donors_ok": int(len(gok)), "n_donors_unavailable": int(len(g) - len(gok)),
                     "n_batches": int(len(b)),
                     "batch_spread_log_ratio": spread,
                     "donor_sd_log_ratio": float(gok["log_threshold_ratio"].std(ddof=1))
                     if len(gok) >= 2 else np.nan,
                     "max_abs_frac_delta": float(gok["frac_delta"].abs().max())
                     if len(gok) else np.nan,
                     "residual_batch": bool(np.isfinite(spread) and spread > tol_log_ratio)})
    return per_batch, pd.DataFrame(rows)


def exclusivity_table(frame: pd.DataFrame, pairs: pd.DataFrame, scale: pd.DataFrame,
                      min_cells: int = MIN_PAIR_CELLS) -> pd.DataFrame:
    """``r = sigma2 / sigma1`` per donor per declared pair, on cohort-scale axes."""
    rows = []
    for a, b in pairs[["marker_a", "marker_b"]].itertuples(index=False):
        if a not in frame.columns or b not in frame.columns \
                or a not in scale.index or b not in scale.index:
            continue
        for d, t in frame.groupby("donor_id", sort=True):
            x = t[[a, b]].to_numpy(float)
            ok = np.isfinite(x).all(axis=1)
            xs = SC.apply_scale(x[ok], [a, b], scale).astype(float)
            r = exclusivity_ratio(xs[:, 0], xs[:, 1], min_cells)
            rows.append({"marker_a": a, "marker_b": b, "donor_id": d,
                         "n_cells": int(ok.sum()), "r": r})
    return pd.DataFrame(rows, columns=["marker_a", "marker_b", "donor_id", "n_cells", "r"])


# ----------------------------------------------------------------- (b) marker invariance
def _dummies(labels: Sequence[str], reference: str | None = None) -> np.ndarray:
    levels = sorted(set(labels))
    if reference in levels:
        levels.remove(reference)
    else:
        levels = levels[1:]
    return np.array([[1.0 if x == lv else 0.0 for lv in levels] for x in labels]).reshape(
        len(labels), len(levels))


def _rss(X: np.ndarray, y: np.ndarray) -> float:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    return float(r @ r)


def _f_stat(y, Xr, Xf, df_num, df_den) -> float:
    rss_f = _rss(Xf, y)
    rss_r = _rss(Xr, y)
    if rss_f <= 1e-300:
        return float("inf") if rss_r > 1e-300 else float("nan")
    return ((rss_r - rss_f) / df_num) / (rss_f / df_den)


def design_table(donors: pd.DataFrame, reference_disease: str = "ND") -> pd.DataFrame:
    """Batch x disease donor counts, with each batch's confounding named."""
    rows = []
    for batch, g in donors.groupby("batch", sort=True, dropna=False):
        counts = g["disease"].value_counts().to_dict()
        rows.append({"batch": batch, "n_donors": len(g),
                     **{f"n_{k}": int(v) for k, v in sorted(counts.items())},
                     "n_disease_levels": len(counts),
                     "has_reference_disease": reference_disease in counts,
                     "within_batch_contrast": len(counts) >= 2})
    return pd.DataFrame(rows).fillna(0)


def _n_within_batch_arrangements(batches: np.ndarray, disease: np.ndarray) -> int:
    total = 1
    for b in np.unique(batches):
        lab = disease[batches == b]
        n = math.factorial(len(lab))
        for c in np.unique(lab, return_counts=True)[1]:
            n //= math.factorial(int(c))
        total *= n
    return total


def invariance_test(shifts: pd.DataFrame, donors: pd.DataFrame, n_perm: int = N_PERM,
                    seed: int = 0, reference_disease: str = "ND",
                    alpha: float = ALPHA) -> pd.DataFrame:
    """Per marker, ``shift_log ~ batch + disease`` with a within-batch permutation p.

    ``status`` is ``tested``, ``confounded`` (disease not estimable beyond batch: no batch
    holds two disease levels), ``saturated`` (no residual degrees of freedom) or
    ``too_few_donors``.  ``effect_<level>`` is the disease coefficient against
    ``reference_disease`` in ``shift_log`` units, batch held fixed; since the applied scale is
    ``exp(-shift_log)``, a positive effect means UniFORM scales that disease level's values
    DOWN relative to the reference level's.  ``flag`` is ``q < alpha`` (BH over tested markers).
    """
    rng = np.random.default_rng(seed)
    d = donors.dropna(subset=["disease", "batch"])
    rows = []
    for marker, g in shifts.groupby("marker", sort=True):
        t = g.merge(d, on="donor_id", how="inner")
        n = len(t)
        base = {"marker": marker, "n_donors": n,
                "n_donors_without_metadata": int(len(g) - n),
                "n_batches": int(t["batch"].nunique()),
                "disease_levels": ",".join(sorted(t["disease"].unique()))}
        if n < 3 or t["disease"].nunique() < 2:
            rows.append({**base, "status": "too_few_donors", "F": np.nan, "p_perm": np.nan})
            continue
        y = t["shift_log"].to_numpy(float)
        batch = t["batch"].astype(str).to_numpy()
        dis = t["disease"].astype(str).to_numpy()
        Xr = np.column_stack([np.ones(n), _dummies(list(batch))])
        Dd = _dummies(list(dis), reference_disease)
        Xf = np.column_stack([Xr, Dd])
        rr, rf = np.linalg.matrix_rank(Xr), np.linalg.matrix_rank(Xf)
        df_num, df_den = rf - rr, n - rf
        n_arr = _n_within_batch_arrangements(batch, dis)
        contrast = [b for b in np.unique(batch) if len(np.unique(dis[batch == b])) >= 2]
        base.update({"df_disease": int(df_num), "df_resid": int(df_den),
                     "n_within_batch_arrangements": int(min(n_arr, 10**12)),
                     "batches_with_contrast": ",".join(contrast)})
        if df_num == 0 or n_arr <= 1:
            rows.append({**base, "status": "confounded", "F": np.nan, "p_perm": np.nan})
            continue
        if df_den <= 0:
            rows.append({**base, "status": "saturated", "F": np.nan, "p_perm": np.nan})
            continue
        beta, *_ = np.linalg.lstsq(Xf, y, rcond=None)
        levels = [lv for lv in sorted(set(dis)) if lv != reference_disease] \
            if reference_disease in set(dis) else sorted(set(dis))[1:]
        effects = {f"effect_{lv}": float(b) for lv, b in zip(levels, beta[Xr.shape[1]:], strict=True)}
        f_obs = _f_stat(y, Xr, Xf, df_num, df_den)
        groups = [np.flatnonzero(batch == b) for b in np.unique(batch)]
        hits = 0
        for _ in range(n_perm):
            perm = dis.copy()
            for idx in groups:
                perm[idx] = rng.permutation(dis[idx])
            Xp = np.column_stack([Xr, _dummies(list(perm), reference_disease)])
            if np.linalg.matrix_rank(Xp) != rf:
                continue
            if _f_stat(y, Xr, Xp, df_num, df_den) >= f_obs - 1e-12:
                hits += 1
        p = (hits + 1) / (n_perm + 1)
        rows.append({**base, **effects, "status": "tested", "F": f_obs, "p_perm": p,
                     "min_attainable_p": max(1.0 / n_arr, 1.0 / (n_perm + 1))})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_bh"] = np.nan
    tested = out["status"] == "tested"
    if tested.any():
        out.loc[tested, "q_bh"] = bh(out.loc[tested, "p_perm"].to_numpy(float))
    out["flag"] = np.where(out["status"] == "tested",
                           np.where(out["q_bh"] < alpha, "disease_associated", "not_flagged"),
                           out["status"])
    out["watchlist"] = out["marker"].isin(DISEASE_WATCHLIST)
    return out


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q values."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


# ---------------------------------------------------------------------------------- report
def _md_table(df: pd.DataFrame, cols: list[str], max_rows: int = 60) -> str:
    if df.empty:
        return "_none_\n"
    d = df[[c for c in cols if c in df.columns]].head(max_rows)
    fmt = d.copy()
    for c in fmt.columns:
        if pd.api.types.is_float_dtype(fmt[c]):
            fmt[c] = fmt[c].map(lambda v: "" if pd.isna(v) else f"{v:.4g}")
    lines = ["| " + " | ".join(map(str, fmt.columns)) + " |",
             "|" + "---|" * len(fmt.columns)]
    lines += ["| " + " | ".join(map(str, r)) + " |" for r in fmt.itertuples(index=False)]
    more = f"\n_{len(df) - max_rows} more rows in the table beside this file._\n" \
        if len(df) > max_rows else "\n"
    return "\n".join(lines) + more


def summary_markdown(restore_summary: pd.DataFrame, copos: pd.DataFrame,
                     exclusivity: pd.DataFrame, invariance: pd.DataFrame,
                     design: pd.DataFrame, skipped: pd.DataFrame, params: dict) -> str:
    """The human-readable summary.  States findings; decides nothing."""
    out = ["# UniFORM audit (Step 8e)", "",
           "Evidence for a person reviewing a UniFORM refit. Nothing here changes a shift, "
           "a scale or a threshold.", "", "## Parameters and inputs", ""]
    out += [f"- `{k}`: {v}" for k, v in params.items()] + [""]

    out += ["## (a) RESTORE residual-batch test", "",
            "`batch_spread_log_ratio` is the spread of per-batch median log(t_donor / "
            "t_pooled); `residual_batch` is True above "
            f"{params.get('tol_log_ratio', TOL_LOG_RATIO):.3f}.", ""]
    rb = restore_summary.sort_values(["values", "residual_batch", "batch_spread_log_ratio"],
                                     ascending=[True, False, False]) \
        if not restore_summary.empty else restore_summary
    out.append(_md_table(rb, ["values", "marker", "partner", "n_donors_ok",
                              "n_donors_unavailable", "n_batches", "batch_spread_log_ratio",
                              "donor_sd_log_ratio", "max_abs_frac_delta", "residual_batch"]))
    out += ["### Co-positivity on a different partner (median over donors)", ""]
    if not copos.empty:
        cp = (copos.groupby(["values", "marker", "partner", "copos_partner"])
              [["copos_donor", "copos_pooled"]].median().reset_index())
        out.append(_md_table(cp, list(cp.columns)))
    else:
        out.append("_no marker has a second declared partner_\n")
    out += ["### Exclusivity ratio r = sigma2/sigma1 (cohort-scale axes)", "",
            "Near 1: the pair lies on the two axes; near 0: along the diagonal. Read across "
            "donors for one pair; population balance also lowers r.", ""]
    if not exclusivity.empty:
        ex = (exclusivity.groupby(["marker_a", "marker_b"])["r"]
              .agg(["count", "min", "median", "max"]).reset_index())
        out.append(_md_table(ex, list(ex.columns)))
    else:
        out.append("_none_\n")

    out += ["## (b) Marker-invariance audit", "", "### Design: batch x disease", ""]
    out.append(_md_table(design, list(design.columns)))
    no_ref = design.loc[~design["has_reference_disease"].astype(bool), "batch"].tolist() \
        if not design.empty else []
    if no_ref:
        out += [f"**Batches with no {params.get('reference_disease', 'ND')} donor:** "
                + ", ".join(map(str, no_ref)) + ". Disease and batch are confounded there; "
                "those donors contribute nothing to the within-batch permutation.", ""]
    if not design.empty and not design["within_batch_contrast"].astype(bool).any():
        out += ["**No batch holds two disease levels: the disease term is not estimable "
                "beyond batch for any marker.**", ""]
    if not invariance.empty:
        inv = invariance.assign(_w=~invariance["watchlist"]).sort_values(
            ["_w", "flag", "q_bh"]).drop(columns="_w")
        cols = ["marker", "watchlist", "status", "flag", "n_donors", "df_disease",
                "F", "p_perm", "q_bh", "min_attainable_p"] + \
            [c for c in inv.columns if c.startswith("effect_")] + ["batches_with_contrast"]
        out.append(_md_table(inv, cols, max_rows=200))
        flagged = invariance.loc[invariance["flag"] == "disease_associated", "marker"].tolist()
        conf = invariance.loc[invariance["flag"] == "confounded", "marker"].tolist()
        out += [f"Flagged disease-associated: {', '.join(flagged) or 'none'}.",
                f"Not testable (confounded): {', '.join(conf) or 'none'} -- not a pass.", ""]
    out += ["Watchlist markers are listed first because their background is expected to move "
            "with disease; the list is never used to flag.", ""]
    if not skipped.empty:
        out += ["### Markers without a fitted shift", ""]
        out.append(_md_table(skipped, ["marker", "provenance"], max_rows=200))
    return "\n".join(out)


# ------------------------------------------------------------------------------------- CLI
def run(expression: Sequence[pathlib.Path | str], shifts_path: pathlib.Path | str,
        donor_metadata: pathlib.Path | str, out: pathlib.Path | str,
        acquisition_runs: pathlib.Path | str | None = None,
        pairs_path: pathlib.Path | str | None = None,
        scale_path: pathlib.Path | str | None = None,
        reference_column: str | None = None,
        partner_levels: pathlib.Path | str | None = None,
        apply_shifts: bool = False, disease_col: str = "donor_status",
        batch_col: str = "stain_group", reference_disease: str = "ND",
        k_sigma: float = K_SIGMA, quantile: float = CONTROL_QUANTILE, margin: float = MARGIN,
        min_controls: int = MIN_CONTROL_CELLS, tol_log_ratio: float = TOL_LOG_RATIO,
        n_perm: int = N_PERM, seed: int = 0, alpha: float = ALPHA) -> dict:
    """Both audits end to end; writes the tables and ``SUMMARY.md`` under ``out``."""
    if reference_column is None and partner_levels is None:
        raise ValueError("partner positivity needs a reference: pass --reference-column or "
                         "--partner-levels; the audit does not call positivity from a "
                         "frequency")
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pairs = load_pairs(pairs_path)
    triples = partner_triples(pairs)
    markers = sorted(set(pairs["marker_a"]) | set(pairs["marker_b"]))
    shifts, skipped = load_shifts(shifts_path)
    donors = load_donor_metadata(donor_metadata, acquisition_runs, disease_col, batch_col)

    given = load_expression(expression, markers, reference_column)
    aligned = apply_uniform_scale(given, shifts, markers) if apply_shifts else given

    if partner_levels is not None:
        lv = pd.read_csv(partner_levels, dtype={"donor_id": str})
        if "donor_id" in lv.columns:
            lv["donor_id"] = lv["donor_id"].map(normalize_donor)
    else:
        lv = reference_levels(aligned, markers, reference_column)
    pos = positive_masks(aligned, lv, markers)

    runs = [("pre_alignment", given), ("aligned", aligned)] if apply_shifts \
        else [("as_given", given)]
    th, cp = [], []
    for label, frame in runs:
        a, b = restore_test(frame, triples, pos, label, k_sigma, quantile, margin, min_controls)
        th.append(a)
        cp.append(b)
    thresholds, copos = pd.concat(th, ignore_index=True), pd.concat(cp, ignore_index=True)
    per_batch, rsum = restore_batch_summary(thresholds, donors, tol_log_ratio)

    scale = pd.read_parquet(scale_path) if scale_path else derive_cohort_scale(aligned, markers)
    if "marker" in scale.columns:
        scale = scale.set_index("marker")
    excl = exclusivity_table(aligned, pairs, scale)

    audited = donors.loc[donors["donor_id"].isin(shifts["donor_id"])]
    design = design_table(audited, reference_disease)
    inv = invariance_test(shifts, donors, n_perm, seed, reference_disease, alpha)

    thresholds.to_parquet(out / "restore_thresholds.parquet", index=False)
    copos.to_parquet(out / "restore_copositivity.parquet", index=False)
    lv.to_csv(out / "partner_levels.csv", index=False)
    per_batch.to_csv(out / "restore_per_batch.csv", index=False)
    rsum.to_csv(out / "restore_summary.csv", index=False)
    excl.to_parquet(out / "exclusivity.parquet", index=False)
    design.to_csv(out / "invariance_design.csv", index=False)
    inv.to_csv(out / "invariance.csv", index=False)
    skipped.to_csv(out / "invariance_skipped_markers.csv", index=False)
    params = {"expression": [str(p) for p in expression], "shifts": str(shifts_path),
              "apply_shifts": apply_shifts, "pairs": str(pairs_path or "DEFAULT_PAIRS"),
              "partner_positivity": (f"declared levels {partner_levels}" if partner_levels
                                     else f"per-donor p{REFERENCE_QUANTILE:g} of "
                                          f"`{reference_column}`"),
              "scale": str(scale_path or "derived as run_step3.build_scale"),
              "k_sigma": k_sigma, "quantile": quantile, "margin": margin,
              "min_controls": min_controls, "tol_log_ratio": round(tol_log_ratio, 4),
              "n_perm": n_perm, "seed": seed, "alpha": alpha,
              "reference_disease": reference_disease,
              "excluded_donors_dropped": given.attrs.get("dropped", []),
              "central_exclusions": sorted(DONOR_EXCLUSIONS)}
    (out / "SUMMARY.md").write_text(summary_markdown(rsum, copos, excl, inv, design,
                                                     skipped, params))
    return {"thresholds": thresholds, "copositivity": copos, "restore_summary": rsum,
            "restore_per_batch": per_batch, "exclusivity": excl, "design": design,
            "invariance": inv, "skipped": skipped, "partner_levels": lv}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expression", nargs="+", required=True, type=pathlib.Path,
                    help="Parquet file(s) or donor_id= partitioned root(s)")
    ap.add_argument("--apply-shifts", action="store_true",
                    help="--expression is pre-alignment: audit it as given AND after the "
                         "shifts JSON's scale")
    ap.add_argument("--shifts", required=True, type=pathlib.Path,
                    help="run_step2's shifts JSON")
    ap.add_argument("--donor-metadata", required=True, type=pathlib.Path)
    ap.add_argument("--acquisition-runs", type=pathlib.Path, default=None,
                    help="data/donor_acquisition_runs.csv (donor_id, acq_run, stain_group)")
    ap.add_argument("--disease-column", default="donor_status")
    ap.add_argument("--batch-column", default="stain_group",
                    help="stain_group or acq_run -- never the Fusion version")
    ap.add_argument("--reference-disease", default="ND")
    ap.add_argument("--pairs", type=pathlib.Path, default=None,
                    help="CSV marker_a,marker_b or fundamental_pair_diagnostics.json")
    ap.add_argument("--scale", type=pathlib.Path, default=None,
                    help="cohort floor/top parquet (marker, floor, top, span)")
    ap.add_argument("--reference-column", default=None,
                    help="boolean column marking each donor's reference population")
    ap.add_argument("--partner-levels", type=pathlib.Path, default=None,
                    help="CSV marker,level[,donor_id] on the aligned scale")
    ap.add_argument("--k-sigma", type=float, default=K_SIGMA)
    ap.add_argument("--quantile", type=float, default=CONTROL_QUANTILE)
    ap.add_argument("--margin", type=float, default=MARGIN)
    ap.add_argument("--min-controls", type=int, default=MIN_CONTROL_CELLS)
    ap.add_argument("--tol-log-ratio", type=float, default=TOL_LOG_RATIO)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args(argv)
    res = run(a.expression, a.shifts, a.donor_metadata, a.out, a.acquisition_runs, a.pairs,
              a.scale, a.reference_column, a.partner_levels, a.apply_shifts,
              a.disease_column, a.batch_column, a.reference_disease, a.k_sigma, a.quantile,
              a.margin, a.min_controls, a.tol_log_ratio, a.n_perm, a.seed, a.alpha)
    rb = res["restore_summary"]
    inv = res["invariance"]
    print(f"RESTORE: {int(rb['residual_batch'].sum()) if len(rb) else 0} of {len(rb)} "
          f"(values, marker, partner) rows show residual batch")
    if len(inv):
        print("invariance: " + ", ".join(f"{k} {v}" for k, v in
                                         inv["flag"].value_counts().items()))
    print(f"wrote {a.out / 'SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
