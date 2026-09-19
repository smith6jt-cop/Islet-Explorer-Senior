#!/usr/bin/env python3
"""Write the QuPath import package so the calls can be reviewed on the images.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the Nimbus/Leiden chain
(``substrate.py`` -> ``positivity.py`` -> ``call.py`` -> ``validate.py`` -> ``export_qupath.py``,
with ``metacluster.py``) that the typing README used to describe.  It did not produce the calls on
disk.  The current typing is the subpanel classification in ``subpanel.py`` / ``build_cells.py``;
see ``README.md``, "Retired, and not to be revived as fallbacks".  What follows is kept as a record
of what was measured and why, not as a lane to fall back to.

QuPath is the only place these calls can actually be judged, so this is the deliverable.  The
format is the repo's existing contract, read by ``scripts/groovy/import_cell_assignments.groovy``:
one ``cell_assignments_<image>.csv`` per image, matched to detections by the QuPath UUID and the
exact image name, plus a manifest carrying each file's size and SHA-256 which the importer checks
before it changes a single PathClass.

Uncertainty is exported, not hidden.  ``ambiguous`` and ``unavailable`` cells never receive the
best candidate as their class -- they get their own classes in QuPath -- while the candidate is
still carried in ``best_broad_type`` so a reviewer can see what the workflow was weighing.

Two per-cell fields make the review sharper than a flat label set:

``assignment_status``
    ``anchor`` -- the cell is itself positive for one of its lineage's markers, so the call rests
    on that cell's own evidence.  ``inferred`` -- the cell was named by the profile of the cluster
    it belongs to.  Reviewing anchors and inferred cells separately says whether the clustering is
    carrying the call or the marker is.
``confidence``
    the cell's own strongest Nimbus confidence among its assigned lineage's markers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import lineages as L  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"
OUT = DEV / "qupath_export"
MANIFESTS = DEV / "manifests"

COLUMNS = ("object_id", "image", "broad_type", "specific_type", "assignment_status", "confidence",
           "best_broad_type", "best_specific_type", "broad_assignment_status",
           "specific_assignment_status", "ranked_broad_probabilities",
           "ranked_specific_probabilities", "ranked_probabilities")
#: A cell is called an ``anchor`` when it is itself positive for one of its lineage's markers.
#: Nimbus's own decision boundary; used only to split anchors from inferred cells for review, and
#: never to make, withhold or change a call.
ANCHOR_CONF = 0.5


def image_name(donor: str, override: dict[str, str] | None = None) -> str:
    """The name QuPath will compare against.

    The importer refuses a row whose ``image`` column differs from
    ``getProjectEntry().getImageName()``, so this must be the name in the PROJECT, which is not
    always the ingest manifest's ``image_id`` -- a project entry can be renamed, and QuPath's own
    default depends on how the image was added (with or without the ``- resolution #N`` series
    suffix).  Pass ``--image-name <donor>=<exact name>`` when they differ; list the real ones with

        getProject().getImageList().each { println "'" + it.getImageName() + "'" }
    """
    if override and str(donor) in override:
        return override[str(donor)]
    man = json.loads((REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    return [i for i in man["images"] if str(i.get("donor_id")) == str(donor)][0]["image_id"]


def build(donor: str, override: dict[str, str] | None = None,
          include_unavailable: bool = False) -> pd.DataFrame:
    calls = pd.read_parquet(DEV / f"calls/donor_id={donor}/data.parquet")
    pos = pd.read_parquet(DEV / f"positivity/donor_id={donor}/data.parquet")
    markers = [c for c in pos.columns if c in set(L.identity_markers())]
    df = calls.merge(pos[["object_id", *markers]], on="object_id", how="left")

    idx = {m: i for i, m in enumerate(markers)}
    X = df[markers].to_numpy(np.float64)
    # per lineage, this cell's own strongest confidence among that lineage's markers
    lin_names, lin_cols = [], []
    for lin, ms in L.LINEAGE_MARKERS.items():
        use = [idx[m] for m in ms if m in idx]
        if use:
            lin_names.append(lin); lin_cols.append(use)
    P = np.vstack([np.nanmax(X[:, c], axis=1) for c in lin_cols]).T

    lineage = df["lineage"].to_numpy().astype(object)
    specific = df["cell_type"].to_numpy().astype(object)
    pos_of = {n: i for i, n in enumerate(lin_names)}
    own = np.array([P[i, pos_of[l]] if l in pos_of else np.nan for i, l in enumerate(lineage)])

    status = np.where(lineage == L.UNAVAILABLE, "unavailable",
              np.where(lineage == L.AMBIGUOUS, "ambiguous",
              np.where(own >= ANCHOR_CONF, "anchor", "inferred")))
    named = ~np.isin(lineage, [L.UNAVAILABLE, L.AMBIGUOUS])
    best_i = np.nanargmax(np.where(np.isfinite(P), P, -np.inf), axis=1)
    best_broad = np.array(lin_names, dtype=object)[best_i]

    def ranked(row: np.ndarray) -> str:
        order = np.argsort(-np.where(np.isfinite(row), row, -np.inf))
        return json.dumps([[lin_names[j], round(float(row[j]), 4)] for j in order[:5]
                           if np.isfinite(row[j])])

    out = pd.DataFrame({
        "object_id": df["object_id"].astype(str),
        "image": image_name(donor, override),
        "broad_type": np.where(named, lineage, ""),
        "specific_type": np.where(named, specific, ""),
        "assignment_status": status,
        "confidence": np.round(np.where(np.isfinite(own), own, np.nan), 4),
        "best_broad_type": best_broad,
        "best_specific_type": np.where(named, specific, best_broad),
        "broad_assignment_status": status,
        "specific_assignment_status": status,
    })
    rk = [ranked(P[i]) for i in range(len(P))]
    out["ranked_broad_probabilities"] = rk
    out["ranked_specific_probabilities"] = rk
    out["ranked_probabilities"] = rk
    out = out[list(COLUMNS)]
    if not include_unavailable:
        # Objects with no nucleus are filtered out before any typing happens, so by default they
        # are left out of the export as well: writing them back as "Phenocycler - Unavailable"
        # CHANGES their classification in QuPath, which is the opposite of leaving them alone.
        # Pass --include-unavailable to get the explicit class instead.
        keep = out["assignment_status"] != "unavailable"
        dropped = int((~keep).sum())
        out = out[keep]
        print(f"[{donor}] omitted {dropped:,} objects with no measurement "
              f"(QuPath leaves their existing classification untouched)")
    return out


def write(donors: list[str], override: dict[str, str] | None = None,
          include_unavailable: bool = False) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    MANIFESTS.mkdir(parents=True, exist_ok=True)
    files = []
    for donor in donors:
        frame = build(donor, override, include_unavailable)
        path = OUT / f"cell_assignments_donor_{donor}.csv"
        frame.to_csv(path, index=False)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        files.append({"path": str(path.resolve()), "donor_id": donor,
                      # the importer refuses anything but a full hash, and says so -- a sampled
                      # fingerprint would let a truncated or partially-rewritten CSV through
                      "hash_mode": "full",
                      "size_bytes": path.stat().st_size, "sha256": digest,
                      "rows": int(len(frame)),
                      "status_counts": frame["assignment_status"].value_counts().to_dict(),
                      "broad_counts": frame.loc[frame.broad_type != "", "broad_type"]
                                           .value_counts().to_dict()})
        print(f"[{donor}] {path.name}  rows={len(frame):,}  sha256={digest[:16]}…")
        for k, v in files[-1]["status_counts"].items():
            print(f"      {k:<12} {v:>9,}  {v/len(frame)*100:>5.1f}%")
    man = MANIFESTS / "qupath_export.json"
    man.write_text(json.dumps({"schema": 1, "producer": "scripts/typing", "files": files}, indent=2))
    print(f"\nmanifest: {man}")
    print(f"import with:\n  QuPath script --project <project.qpproj> --image \"<image name>\" \\\n"
          f"    {REPO}/core/phenocycler/scripts/groovy/import_cell_assignments.groovy \\\n"
          f"    --args \"{OUT},specific\"")
    return man


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--include-unavailable", action="store_true",
                    help="also export objects with no measurement, as an explicit Unavailable "
                         "class; by default they are omitted so QuPath leaves them alone")
    ap.add_argument("--image-name", nargs="*", default=[], metavar="DONOR=NAME",
                    help="exact QuPath project entry name for a donor, when it differs from the "
                         "ingest manifest's image_id")
    a = ap.parse_args()
    override = {}
    for item in a.image_name:
        if "=" not in item:
            raise SystemExit(f"--image-name expects DONOR=NAME, got {item!r}")
        k, v = item.split("=", 1)
        override[k.strip()] = v.strip()
    write(a.donors, override, a.include_unavailable)


if __name__ == "__main__":
    main()
