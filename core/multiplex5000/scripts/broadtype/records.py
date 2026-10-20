"""The common record every arm writes, so the gate reads one shape.

Calls parquet (one per arm x donor x configuration x replicate x ablation):
    object_id        str    the QuPath UUID -- the only join key
    label            str    one of declarations.LABELS, or 'Ambiguous' / 'Unavailable'
    score            float  the arm's own confidence-like number (NaN allowed for Unavailable)
    reason           str    '' when called; the recorded reason otherwise (never empty for
                            Ambiguous/Unavailable)
The JSONL record beside it carries provenance and the identity-input list the tiering test
reads.  An arm that emits no call for a frozen cell writes 'Unavailable' with reason
'arm_no_output' -- it never drops the cell.
"""
from __future__ import annotations

import json
import pathlib
import time
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

import pins

_git_head = pins.git_head

AMBIGUOUS = "Ambiguous"
UNAVAILABLE = "Unavailable"
CALL_COLUMNS = ("object_id", "label", "score", "reason")


@dataclass
class ArmRecord:
    arm: str
    donor_id: str
    config_id: str                      # 'primary', 'ablate:<marker>', 'seed:<n>', 'subsample:<n>'
    identity_markers: list              # the ONLY markers that entered identity (type markers)
    scope: str                          # 'frame' or 'donor'
    n_cells: int
    n_called: int
    n_ambiguous: int
    n_unavailable: int
    seconds: float
    env: str                            # interpreter path or R version string
    pinned: dict = field(default_factory=dict)      # run_id, content ids, external pins used
    params: dict = field(default_factory=dict)
    notes: str = ""
    code: dict = field(default_factory=_git_head)
    written_utc: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


def out_dir(arm: str, donor: str, scope: str) -> pathlib.Path:
    return pins.OUT / "arms" / arm / f"donor_id={donor}" / scope


def validate_calls(df: pd.DataFrame, expected_ids: np.ndarray | None = None,
                   labels: set | None = None) -> pd.DataFrame:
    missing = [c for c in CALL_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"calls missing columns {missing}")
    if df["object_id"].duplicated().any():
        raise ValueError("duplicate object_id in calls")
    bad = df.loc[df["label"].isin([AMBIGUOUS, UNAVAILABLE]) & (df["reason"].fillna("") == "")]
    if len(bad):
        raise ValueError(f"{len(bad)} Ambiguous/Unavailable calls carry no reason")
    if labels is not None:
        unknown = set(df["label"]) - set(labels) - {AMBIGUOUS, UNAVAILABLE}
        if unknown:
            raise ValueError(f"undeclared labels in calls: {sorted(unknown)}")
    if expected_ids is not None:
        want = pd.Index(expected_ids.astype(str))
        got = pd.Index(df["object_id"].astype(str))
        if not want.equals(got):
            lost = want.difference(got)
            extra = got.difference(want)
            raise ValueError(f"calls do not cover the frozen cell set: {len(lost)} missing, {len(extra)} extra")
    return df


def write_calls(df: pd.DataFrame, record: ArmRecord, *, labels: set | None = None,
                expected_ids: np.ndarray | None = None) -> pathlib.Path:
    d = out_dir(record.arm, record.donor_id, record.scope)
    d.mkdir(parents=True, exist_ok=True)
    df = validate_calls(df[list(CALL_COLUMNS)].copy(), expected_ids, labels)
    stem = record.config_id.replace(":", "_")
    path = d / f"calls__{stem}.parquet"
    df.to_parquet(path, index=False)
    with (d / "records.jsonl").open("a") as fh:
        fh.write(json.dumps({**asdict(record), "calls": str(path)}, default=str) + "\n")
    return path


def fill_unavailable(df: pd.DataFrame, expected_ids: np.ndarray, reason: str = "arm_no_output") -> pd.DataFrame:
    """Every frozen cell gets a row; cells the arm did not emit are Unavailable(reason)."""
    base = pd.DataFrame({"object_id": pd.Series(expected_ids).astype(str)})
    df = df.copy()
    df["object_id"] = df["object_id"].astype(str)
    out = base.merge(df, on="object_id", how="left")
    miss = out["label"].isna()
    out.loc[miss, "label"] = UNAVAILABLE
    out.loc[miss, "reason"] = reason
    out["score"] = out["score"].astype(float)
    out["reason"] = out["reason"].fillna("")
    return out[list(CALL_COLUMNS)]
