#!/usr/bin/env python3
"""Stage 6 build driver -- run the decision-record adapters over what exists on disk.

``decision_record.py`` is a LIBRARY: it has no ``main``, no argparse and no
``if __name__ == "__main__"``. It exposes four adapters, each of which turns ONE row of ONE
lane into a :class:`DecisionRecord`, plus ``records_to_frame`` / ``assert_no_absence_as_zero``.
Nothing in the repository imports it. This file is the missing driver: it finds the lanes'
outputs, feeds them to the adapters unchanged, and writes one combined table.

What each adapter consumes, read off the code rather than assumed:

* ``from_group_call`` / ``from_group_calls``  -- rows of ``nullgate.annotate_groups``
  (columns ``group, n_cells, lineage, reason, real, null_mean, null_p95, z, runners_up``).
  Produced by ``run_donor.py`` (Stage 4). NOT PRESENT AT BUILD TIME: a Stage-4 run is in
  progress (KODAMA, both same-slide serials) and
  ``data/development/xenium_celltype/donor_id=*/roi=*/serial=*/`` still holds only the
  ``input/`` prep files -- no ``group_calls_<arm>.csv``. This adapter is skipped, not
  removed: re-run this driver once Stage 4 writes those files and it will pick them up
  through the same discovery, which is why the skip is printed rather than silent.
* ``from_cell_call`` / ``from_cell_calls``    -- rows of ``celltype/score.call_from_scores``
  (``call, score, runner_up, margin, reason, n_candidates``). The protein arm. RUNS, over
  both the native call tables and their lifted form (see below).
* ``from_broad_call`` / ``from_broad_obs``    -- ``adata.obs`` of ``hierarchical_phenotype``
  tier 1 (``celltype_broad, broad_confidence, broad_n_identity_genes, broad_lowsignal,
  broad_marker_margin``). No ``*_phenotyped.h5ad`` exists for any of this cohort's four
  bundles, so this adapter does not run either.
* ``from_composition_row`` / ``from_composition`` -- rows of ``run_donor.composition``
  (``lineage, status, n_cells, frac, ci_lo, ci_hi``). Same Stage-4 dependency as adapter 1;
  does not run.

Two protein inputs, both adapter 2, kept apart by the ``source`` column:

* ``sameslide_protein``  -- ``data/development/celltype/sameslide_protein/donor_id=*/
  data.parquet``, the lane's OWN output, one row per PhenoCycler object, carrying the cosine,
  the margin, the runner-up and ``n_candidates``.
* ``lifted_protein``     -- ``core/phenocycler/data/integration/sameslide/donor_id=*/roi=*/
  serial=*/protein_calls.parquet``, the same calls carried onto Xenium cells by
  ``sameslide/lift_calls.py``. That transfer keeps ``call`` and ``reason`` and drops
  ``score``/``runner_up``/``margin``/``n_candidates`` -- it is a categorical transfer by
  largest shared area, not a rescoring -- so every statistic column is NaN for this source BY
  CONSTRUCTION, and NaN is what gets written. The two uncertainty columns the lift does add
  (``call_share``, ``call_purity``) are areas of agreement, not a probability and not a
  cosine, so they are NOT laundered into ``confidence`` or into a field named ``cosine``.

Detection modes are NOT one of the four adapters -- no adapter reads
``eval_detection.py``'s output -- so :func:`detection_records` builds them through the plain
``DecisionRecord`` constructor, which is what enforces the rule. The decision recorded is
``eval_detection``'s own: selectivity is defined there as positive-control rate / clean-
negative rate, and the mode with the highest selectivity is the one that measurement selects.

    unit  = detect:<bundle>:<mode>
    award = the highest-selectivity mode; statistic = selectivity_gain over ``raw>=1``
    loser = unavailable, its own gain kept in ``statistic_value``, its name in ``runners_up``
            (the ``cleared_but_outranked`` status is NOT usable here: it requires an
            UNMEASURED statistic, and a losing mode's gain is measured)
    none  = bundles 0041393 and 0041396 have ``positive`` NaN for every mode, because their
            ``.classes.json`` carries only ``__pooled__`` -- no protein reference class met
            ``MIN_CLASS_CELLS``, so ``positive_rate`` divided by an empty denominator. Every
            mode is ``unavailable``; ``statistic_value`` is NaN and ``statistic_support`` is
            the MEASURED 0 -- zero modes had a measurable selectivity -- in the same row.

``n`` for a detection record is the number of cells the rates were computed over, summed from
the ``.classes.json`` sidecar (``class_tables`` writes one entry per context class, over every
row of the matrix, so the sum is the matrix's row count).

    python3 scripts/xenium_celltype/build_decision_records.py            # from core/phenocycler
    python3 scripts/xenium_celltype/build_decision_records.py --out <dir>
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from decision_record import (  # noqa: E402
    AWARDED, CLEARED_BUT_OUTRANKED, NO_GROUP_CLEARED, UNAVAILABLE,
    DecisionRecord, Statistic, assert_no_absence_as_zero, from_cell_calls,
    from_composition, from_group_calls, records_to_frame,
)

PHENO = HERE.parents[1]                       # core/phenocycler
REPO = HERE.parents[3]                        # monorepo root
PROVENANCE = ("arm", "source", "adapter", "donor_id", "serial", "input_path")

#: (donor_id, serial) for the same-slide pancreas pairs on disk.
SAMESLIDE = (("112", "0081366"), ("6533", "0080870"))
BUNDLES = ("0041393", "0041396", "0080870", "0081366")

SCHEMA = pa.schema(
    [pa.field("unit", pa.string()), pa.field("n", pa.int64()),
     pa.field("label", pa.string()), pa.field("reason", pa.string()),
     pa.field("statistic", pa.string()), pa.field("statistic_value", pa.float64()),
     pa.field("statistic_margin", pa.float64()), pa.field("statistic_lo", pa.float64()),
     pa.field("statistic_hi", pa.float64()), pa.field("statistic_support", pa.int64()),
     pa.field("null_kind", pa.string()), pa.field("null_mean", pa.float64()),
     pa.field("null_p95", pa.float64()), pa.field("null_z", pa.float64()),
     pa.field("confidence", pa.float64()), pa.field("runners_up", pa.string()),
     pa.field("status", pa.string())]
    + [pa.field(c, pa.string()) for c in PROVENANCE])

NUMERIC = ("n", "statistic_value", "statistic_margin", "statistic_lo", "statistic_hi",
           "statistic_support", "null_mean", "null_p95", "null_z", "confidence")


def _typed(frame: pd.DataFrame, **prov: str) -> pa.Table:
    """Fixed dtypes so every chunk appends to one parquet schema. Nothing is filled.

    ``statistic_support`` goes through nullable ``Int64`` precisely so a lane that reports NO
    support count stays null while a measured 0 stays 0 -- the distinction the module is for.
    """
    out = frame.copy()
    for c in ("unit", "label", "reason", "statistic", "null_kind", "runners_up", "status"):
        out[c] = out[c].astype("object").where(out[c].notna(), None)
        out[c] = [None if v is None else str(v) for v in out[c]]
    out["n"] = pd.to_numeric(out["n"], errors="raise").astype("int64")
    out["statistic_support"] = pd.array(
        [None if v is None or (isinstance(v, float) and np.isnan(v)) else int(v)
         for v in out["statistic_support"]], dtype="Int64")
    for c in ("statistic_value", "statistic_margin", "statistic_lo", "statistic_hi",
              "null_mean", "null_p95", "null_z", "confidence"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float64")
    for k, v in prov.items():
        out[k] = str(v)
    return pa.Table.from_pandas(out[list(SCHEMA.names)], schema=SCHEMA, preserve_index=False)


# ---------------------------------------------------------------------------
# adapter 2 -- the protein arm
# ---------------------------------------------------------------------------
def protein_chunks(path: pathlib.Path, id_col: str, unit_prefix: str, batch: int = 200_000):
    """Yield ``(frame_of_records, n_rows)`` for one call table, a batch at a time.

    Chunked because the four protein tables are 2.9M rows between them and the adapter
    materialises one frozen dataclass per row; a single list of those is an inspection tool,
    not a table (``from_broad_obs``'s docstring makes the same point about a donor).
    """
    pf = pq.ParquetFile(path)
    for rb in pf.iter_batches(batch_size=batch):
        df = rb.to_pandas()
        units = [f"{unit_prefix}:{v}" for v in df[id_col].astype(str)]
        frame = records_to_frame(from_cell_calls(df, units=units))
        assert_no_absence_as_zero(frame)
        yield frame, len(df)


# ---------------------------------------------------------------------------
# adapters 1 and 4 -- Stage 4, when it lands
# ---------------------------------------------------------------------------
def stage4_dirs() -> list[tuple[str, str, pathlib.Path]]:
    """``(donor, serial, dir)`` for every Stage-4 donor directory that exists."""
    root = REPO / "data/development/xenium_celltype"
    out = []
    for d in sorted(root.glob("donor_id=*/roi=*/serial=*")):
        out.append((d.parents[1].name.split("=", 1)[1], d.name.split("=", 1)[1], d))
    return out


def _unit_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in ("arm", "k") if c in df.columns]


def group_call_frames(dest: pathlib.Path, donor: str, serial: str):
    """``run_donor`` writes ``group_calls.csv``; ``verify_run``'s fixtures use
    ``group_calls_<arm>.csv``. Both are read, whichever exists."""
    for csv in sorted(list(dest.glob("group_calls.csv")) + list(dest.glob("group_calls_*.csv"))):
        df = pd.read_csv(csv)
        cols = _unit_cols(df)

        def unit(row, _c=cols, _s=serial, _d=donor):
            tag = ":".join(f"{c}{row[c]}" for c in _c)
            return f"group:{_d}:{_s}:{tag}:{row['group']}" if tag else \
                   f"group:{_d}:{_s}:{row['group']}"

        frame = records_to_frame(from_group_calls(df, unit=unit))
        assert_no_absence_as_zero(frame)
        yield csv, frame, len(df)


def composition_frames(dest: pathlib.Path, donor: str, serial: str):
    for csv in sorted(list(dest.glob("composition.csv")) + list(dest.glob("composition_*.csv"))):
        df = pd.read_csv(csv)
        cols = _unit_cols(df)
        recs = []
        if cols:
            for keys, sub in df.groupby(cols, sort=False):
                keys = keys if isinstance(keys, tuple) else (keys,)
                tag = ":".join(f"{c}{v}" for c, v in zip(cols, keys))
                recs.extend(from_composition(sub, unit=f"comp:{donor}:{serial}:{tag}"))
        else:
            recs.extend(from_composition(df, unit=f"comp:{donor}:{serial}"))
        frame = records_to_frame(recs)
        assert_no_absence_as_zero(frame)
        yield csv, frame, len(df)


# ---------------------------------------------------------------------------
# detection modes -- no adapter exists; the plain constructor enforces the rule
# ---------------------------------------------------------------------------
def detection_records(csv: pathlib.Path, bundle: str) -> list[DecisionRecord]:
    df = pd.read_csv(csv)
    side = csv.with_suffix(".classes.json")
    n = 0
    if side.exists():
        n = int(sum(v["n_cells"] for v in json.loads(side.read_text()).values()))
    sel = pd.to_numeric(df["selectivity"], errors="coerce")
    gain = pd.to_numeric(df["selectivity_gain"], errors="coerce")
    support = int(sel.notna().sum())              # a measured count: how many modes had one

    if support == 0:
        return [DecisionRecord(
            unit=f"detect:{bundle}:{m}", n=n, label=None,
            reason="no_positive_control_measured",
            statistic=Statistic(name="selectivity_gain", support=0),
            null=None, confidence=float("nan"), status=UNAVAILABLE)
            for m in df["mode"]]

    order = sel.sort_values(ascending=False, na_position="last").index
    win = order[0]
    second = gain.loc[order[1]] if len(order) > 1 else float("nan")
    margin = gain.loc[win] - second
    losers = [str(df.loc[i, "mode"]) for i in order[1:] if pd.notna(sel.loc[i])]

    recs = []
    for i, row in df.iterrows():
        mode = str(row["mode"])
        if i == win:
            recs.append(DecisionRecord(
                unit=f"detect:{bundle}:{mode}", n=n, label=mode,
                reason="highest_selectivity_over_clean_negatives",
                statistic=Statistic(name="selectivity_gain", value=gain.loc[i],
                                    margin=margin, support=support),
                null=None, confidence=float("nan"),
                runners_up=tuple(losers), status=AWARDED))
        else:
            recs.append(DecisionRecord(
                unit=f"detect:{bundle}:{mode}", n=n, label=None,
                reason=(f"outranked_by:{df.loc[win, 'mode']}" if pd.notna(sel.loc[i])
                        else "selectivity_not_measured"),
                statistic=Statistic(name="selectivity_gain", value=gain.loc[i],
                                    support=support),
                null=None, confidence=float("nan"),
                runners_up=(mode,), status=UNAVAILABLE))
    return recs


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path,
                    default=PHENO / "data/development/xenium_celltype")
    ap.add_argument("--name", default="decision_records")
    args = ap.parse_args(argv)

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)
    table_path = args.out / f"{args.name}.parquet"
    summary_path = args.out / f"{args.name}_summary.csv"
    check_path = args.out / f"{args.name}_nan_check.json"

    writer = pq.ParquetWriter(table_path, SCHEMA)
    tally: list[dict] = []

    def take(frame: pd.DataFrame, n_in: int, **prov):
        writer.write_table(_typed(frame, **prov))
        vc = frame["status"].value_counts()
        tally.append({**prov, "rows_in": n_in, "rows_out": len(frame),
                      "awarded": int(vc.get(AWARDED, 0)),
                      "unavailable": int(vc.get(UNAVAILABLE, 0)),
                      "cleared_but_outranked": int(vc.get(CLEARED_BUT_OUTRANKED, 0)),
                      "no_group_cleared": int(vc.get(NO_GROUP_CLEARED, 0)),
                      "statistic_measured": int(
                          pd.to_numeric(frame["statistic_value"],
                                        errors="coerce").notna().sum())})

    # -- detection modes (Stage 3) ---------------------------------------
    for b in BUNDLES:
        csv = REPO / f"data/development/xenium_celltype/detection_modes_{b}.csv"
        if not csv.exists():
            print(f"[skip] {csv} absent", flush=True)
            continue
        frame = records_to_frame(detection_records(csv, b))
        assert_no_absence_as_zero(frame)
        take(frame, len(frame), arm="detection", source="eval_detection",
             adapter="DecisionRecord(no adapter defined)", donor_id="", serial=b,
             input_path=str(csv))
        print(f"[detect] {b}: {len(frame)} records ({time.time()-t0:.0f}s)", flush=True)

    # -- protein arm (adapter 2) -----------------------------------------
    for donor, serial in SAMESLIDE:
        src = REPO / f"data/development/celltype/sameslide_protein/donor_id={donor}/data.parquet"
        if src.exists():
            tot = 0
            for frame, n_in in protein_chunks(src, "object_id", f"pheno:{donor}"):
                take(frame, n_in, arm="protein", source="sameslide_protein",
                     adapter="from_cell_calls", donor_id=donor, serial="",
                     input_path=str(src))
                tot += len(frame)
            print(f"[protein] donor {donor} native: {tot} records "
                  f"({time.time()-t0:.0f}s)", flush=True)
        else:
            print(f"[skip] {src} absent", flush=True)

        lifted = (PHENO / f"data/integration/sameslide/donor_id={donor}/roi=panc/"
                          f"serial={serial}/protein_calls.parquet")
        if lifted.exists():
            tot = 0
            for frame, n_in in protein_chunks(lifted, "cell_id", f"xenium:{serial}"):
                take(frame, n_in, arm="protein", source="lifted_protein",
                     adapter="from_cell_calls", donor_id=donor, serial=serial,
                     input_path=str(lifted))
                tot += len(frame)
            print(f"[protein] donor {donor} lifted ({serial}): {tot} records "
                  f"({time.time()-t0:.0f}s)", flush=True)
        else:
            print(f"[skip] {lifted} absent", flush=True)

    # -- Stage 4 (adapters 1 and 4), if it has landed --------------------
    seen4 = False
    for donor, serial, dest in stage4_dirs():
        for csv, frame, n_in in group_call_frames(dest, donor, serial):
            seen4 = True
            take(frame, n_in, arm="group", source="nullgate", adapter="from_group_calls",
                 donor_id=donor, serial=serial, input_path=str(csv))
            print(f"[group] {donor}/{serial} {csv.name}: {len(frame)} records", flush=True)
        for csv, frame, n_in in composition_frames(dest, donor, serial):
            seen4 = True
            take(frame, n_in, arm="composition", source="run_donor",
                 adapter="from_composition", donor_id=donor, serial=serial,
                 input_path=str(csv))
            print(f"[comp] {donor}/{serial} {csv.name}: {len(frame)} records", flush=True)
    if not seen4:
        print("[skip] no Stage-4 group_calls/composition on disk; adapters 1 and 4 "
              "(nullgate, composition) did not run", flush=True)

    writer.close()

    summary = (pd.DataFrame(tally)
               .groupby(["arm", "source", "adapter", "donor_id", "serial", "input_path"],
                        as_index=False).sum(numeric_only=True))
    summary.to_csv(summary_path, index=False)
    print("\n" + summary.to_string(index=False), flush=True)

    # -- read back and re-check the rule ---------------------------------
    back = pd.read_parquet(table_path)
    assert_no_absence_as_zero(back)
    check = {"path": str(table_path), "rows": int(len(back)), "columns": {}}
    for c in NUMERIC:
        col = back[c]
        v = pd.to_numeric(col, errors="coerce")
        check["columns"][c] = {
            "nulls": int(col.isna().sum()),
            "min": None if v.notna().sum() == 0 else float(v.min()),
            "max": None if v.notna().sum() == 0 else float(v.max()),
            "zeros": int((v == 0).sum()),
        }
    check["status_counts"] = {k: int(v) for k, v in back["status"].value_counts().items()}
    check["label_nulls"] = int(back["label"].isna().sum())
    check_path.write_text(json.dumps(check, indent=2))
    print("\n" + json.dumps(check, indent=2), flush=True)
    print(f"\n[done] {table_path} ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
