"""Step 4: measure the two lanes against each other on the cells they both call.

This is a **measurement of two methods**, not a release gate for either.  Neither lane has
ground truth: the cosine lane scores each cell against the rule table behind a requirement
gate, this lane clusters first and annotates the cluster.  Where they agree, the answer rests
on two mechanically independent groupings; where they disagree, the disagreement is the
finding, and which one is right is not settled by this file.

**What it joins.**  ``(donor_id, object_id)`` -- the QuPath UUID, never a centroid and never a
donor-name parse.  ``CLAUDE.md`` requires that key end to end and it is the reason the join
can be exact.

**Which cells.**  The cosine lane's per-cell ``calls/`` currently covers one donor, but
``islet_residents_all.parquet`` carries all 20 donors' islet-associated cells (core + peri)
in the same schema.  That is the population this study is about, so it is the population the
comparison runs on -- stated rather than presented as cohort-wide agreement.

**Vocabulary.**  Both lanes answer against the same ``phenotype_workflow.csv``, so labels are
directly comparable with no crosswalk.  That is the whole reason the rule table was held
fixed across the two lanes.  Agreement is reported at two depths, because they are different
claims: exact label, and the top-level branch (``Immune`` / ``Epithelial`` / ``Endothelial``
/ ``Neural`` / ``Mesenchymal``), which is what a nested pair like ``Immune`` vs
``Mac-M2-like`` should be scored as agreeing on.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
sys.path.insert(0, str(HERE))
import annotate                                     # noqa: E402
import panel                                        # noqa: E402
import rules                                        # noqa: E402

COSINE_ISLET = panel.OUT / "islet_residents_all.parquet"
REFUSALS = (annotate.AMBIGUOUS, annotate.UNAVAILABLE)


def branch_of(table: pd.DataFrame) -> dict[str, str]:
    """phenotype -> its top-level branch, so nested labels are not scored as disagreement."""
    anc = annotate.ancestry(table)
    return {p: (anc[p][-1] if anc[p] else p) for p in table["phenotype"]}


def compare(kintsugi_calls: pd.DataFrame, cosine_calls: pd.DataFrame,
            table: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    """Join on (donor_id, object_id) and report agreement at label and branch depth."""
    table = rules.load_rules() if table is None else table
    branch = branch_of(table)

    merged = kintsugi_calls.merge(
        cosine_calls, on=["donor_id", "object_id"], how="inner",
        suffixes=("_kintsugi", "_cosine"))
    if merged.empty:
        raise RuntimeError("no cell is called by both lanes; check the run and the donors")

    both = merged[~merged["call_kintsugi"].isin(REFUSALS)
                  & ~merged["call_cosine"].isin(REFUSALS)].copy()
    both["branch_kintsugi"] = both["call_kintsugi"].map(branch)
    both["branch_cosine"] = both["call_cosine"].map(branch)

    summary = {
        "n_joined": int(len(merged)),
        "n_kintsugi_refused": int(merged["call_kintsugi"].isin(REFUSALS).sum()),
        "n_cosine_refused": int(merged["call_cosine"].isin(REFUSALS).sum()),
        "n_both_called": int(len(both)),
        "label_agreement": float((both["call_kintsugi"] == both["call_cosine"]).mean())
        if len(both) else float("nan"),
        "branch_agreement": float((both["branch_kintsugi"] == both["branch_cosine"]).mean())
        if len(both) else float("nan"),
        "both_refused": int((merged["call_kintsugi"].isin(REFUSALS)
                             & merged["call_cosine"].isin(REFUSALS)).sum()),
    }
    matrix = pd.crosstab(merged["call_kintsugi"], merged["call_cosine"])
    return matrix, summary


def refusal_overlap(kintsugi_calls: pd.DataFrame, cosine_calls: pd.DataFrame) -> dict:
    """Do the two lanes refuse the SAME cells?

    Both lanes leave a large ``Ambiguous`` fraction and both attribute it to the same cause --
    pancreas is ~85% acinar and this 59-plex carries no acinar enzyme marker.  That they agree
    on the *explanation* is not evidence they are refusing the *same cells*, and this is what
    tells them apart.  A high overlap supports one shared cause; a low one means each lane is
    failing somewhere the other is not, and the acinar story covers less than it appears to.
    """
    merged = kintsugi_calls.merge(cosine_calls, on=["donor_id", "object_id"], how="inner",
                                  suffixes=("_kintsugi", "_cosine"))
    k = merged["call_kintsugi"].isin(REFUSALS).to_numpy()
    c = merged["call_cosine"].isin(REFUSALS).to_numpy()
    inter = int((k & c).sum())
    union = int((k | c).sum())
    return {
        "n_joined": int(len(merged)),
        "kintsugi_refused": int(k.sum()),
        "cosine_refused": int(c.sum()),
        "refused_by_both": inter,
        "jaccard": (inter / union) if union else float("nan"),
        "share_of_kintsugi_refusals_also_refused_by_cosine":
            (inter / int(k.sum())) if k.sum() else float("nan"),
        "share_of_cosine_refusals_also_refused_by_kintsugi":
            (inter / int(c.sum())) if c.sum() else float("nan"),
    }
