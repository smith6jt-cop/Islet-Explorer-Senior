import json

import pytest

from analysis.cohesive.common import read_json
from analysis.cohesive.protocol import AMENDMENT, COMPUTABILITY, amended


def archived_shape():
    hypotheses = [{"id": key, "minimum_support": "At least5supported islets per donor",
                   "eligibility": "Native Beta/Alpha/Delta within own true Islet; native T lineage"}
                  for key in COMPUTABILITY]
    return {"family": "spatial_v2_20261006", "hypotheses": hypotheses,
            "statistics": {"minimum_disease_support": "At least 4 T1D and 4 ND donors"},
            "sensitivity": {"pancreatitis": "Exclude6579", "protein_identity": "refined sensitivity",
                            "cell_neighborhood": "50um primary"}}


def test_no_support_gate_or_6579_drop_survives_the_amendment():
    frozen = amended(archived_shape(), read_json(AMENDMENT))
    active = {key: value for key, value in frozen.items() if key != "amendment"}
    for hypothesis in active["hypotheses"]:
        assert "minimum_support" not in hypothesis
        assert hypothesis["computability"] == COMPUTABILITY[hypothesis["id"]]
        assert "Native" not in hypothesis["eligibility"] and "native T" not in hypothesis["eligibility"]
    assert "minimum_disease_support" not in active["statistics"]
    assert set(active["sensitivity"]) == {"cell_neighborhood"}
    assert "Exclude6579" not in json.dumps({k: v for k, v in active.items()
                                            if k not in ("hypotheses", "statistics")})
    assert "refined_class" in active["protein_identity"]


def test_amendment_refuses_another_protocol():
    other = dict(archived_shape(), family="spatial_v1")
    with pytest.raises(ValueError, match="another protocol"):
        amended(other, read_json(AMENDMENT))
