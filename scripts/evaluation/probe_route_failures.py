#!/usr/bin/env python3
"""Small adversarial probes of existing routes; writes evidence, changes no science."""
import json
import runpy
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from multiplex5000.integration.paired_typing import link_candidate_cells
from multiplex5000.integration.sameslide import match_by_centroid
from multiplex5000.stages import maldi_serial_registration as maldi

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs/dynamic_portability_20261102/fault_probes"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fixture = runpy.run_path(str(BASE / "core/multiplex5000/tests/current/test_paired_typing.py"))["_fixture"]
    results = []
    for case in ("duplicate_identifier", "wrong_coordinate_frame", "absent_gene_evidence"):
        pair, protein, pcalls, rna, rcalls, transform = fixture()
        if case == "duplicate_identifier":
            rna.loc[1, "cell_uid"] = rna.loc[0, "cell_uid"]
        elif case == "wrong_coordinate_frame":
            pair = replace(pair, xenium_frame="slide_pixels")
        else:
            rcalls = rcalls.drop(columns="gene:CD3E")
        try:
            link_candidate_cells(pair, protein, pcalls, rna, rcalls, transform)
            result = {"case": case, "status": "unexpected_acceptance"}
        except ValueError as exc:
            result = {"case": case, "status": "rejected", "reason": str(exc)}
        results.append(result)
    # Deliberately identical centroids: one arbitrary mutual-nearest link survives.
    tie = match_by_centroid(np.array([[0., 0.], [0., 0.]]), np.array([[0., 0.]]))
    results.append({"case": "segmentation_split_equal_centroid_tie", "status": "specification_gap",
                    "links": tie.to_dict("records"), "unmatched": 1,
                    "reason": "No tie/ambiguity or split/merge annotation in the returned table"})
    sensitivity = []
    for offset in (4.999, 5.0, 5.001):
        links = match_by_centroid(np.array([[offset, 0.]]), np.array([[0., 0.]]))
        sensitivity.append({"offset_um": offset, "links": len(links)})
    results.append({"case": "existing_5um_boundary_perturbation", "status": "specification_gap",
                    "observations": sensitivity,
                    "reason": "Synthetic sensitivity demonstration only; these offsets are not proposed acceptance tolerances or empirical registration errors"})
    annotations = OUT / "empty_annotations.geojson"
    annotations.write_text('{"type":"FeatureCollection","features":[]}\n')
    task = SimpleNamespace(input=lambda name: annotations, project=None,
                           output=lambda name: OUT / name)
    row = {"tissue_id": "synthetic", "donor_id": "synthetic", "tissue": "pancreas",
           "neighbour": {"tissue": "pancreas", "platform": "phenocycler", "acquisition_id": "synthetic:protein"}}
    with patch.object(maldi, "read_qptiff_channel", return_value=(np.ones((5, 5)), 50)):
        try:
            maldi._fixed(task, row, 0.5)
            results.append({"case": "missing_tissue_outline", "status": "unexpected_acceptance"})
        except maldi.MaldiRegistrationError as exc:
            results.append({"case": "missing_tissue_outline", "status": "rejected", "reason": str(exc)})
    row["neighbour"]["platform"] = "xenium"
    with patch.object(maldi, "_sections", return_value=[row]), patch.object(maldi, "_figure"):
        result = maldi.run(task)
    results.append({"case": "unsupported_xenium_neighbour", "status": "unavailable", "result": result,
                    "scope": "Routing branch probe with mocked inventory, not accepted registration"})
    (OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    if any(r["status"] == "unexpected_acceptance" for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
