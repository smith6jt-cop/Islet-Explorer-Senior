"""Section-local protein rule evidence shared by pLN and spleen.

The tissue supplies its recorded hierarchy. Stage 03 supplies measured marker
values and nucleus admission; Stage 05 supplies the marker shifts. This step
fits the maintainer-selected 0.10 review candidate separately for each marker
and physical section, then evaluates the hierarchy without forcing a final
cell class. A missing marker or shift remains unavailable.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from ..protein_evidence import ProteinEvidence, candidate_calls
from .evaluate import BatchRuleEvaluation, evaluate_rules_batch
from .lymph_node import RuleSource


class ProteinSectionError(ValueError):
    """The declared section measurements and scales cannot be compared."""


@dataclass(frozen=True)
class ProteinSectionEvidence:
    markers: dict[str, ProteinEvidence]
    rules: BatchRuleEvaluation
    missing_markers_by_rule: dict[str, tuple[str, ...]]


def evaluate_protein_section(source: RuleSource, section: str,
                             measurements: Mapping[str, np.ndarray],
                             nucleus_eligible: np.ndarray,
                             shifts: Mapping[str, float | None]) -> ProteinSectionEvidence:
    """Produce per-cell marker and hierarchy states in this section's units.

    ``measurements`` is keyed by actual panel marker names; the rule source
    resolves its recorded aliases. Scales for another section or an absent
    marker are rejected. An unscaled measured marker is reported unavailable.
    Rule results are matches and evidence gaps, not a final served cell class.
    """
    if not measurements:
        raise ProteinSectionError("a protein section needs measured marker columns")
    if set(shifts) - set(measurements):
        raise ProteinSectionError("UniFORM scales name markers absent from this section")
    eligible = np.asarray(nucleus_eligible)
    if eligible.ndim != 1 or eligible.dtype != np.bool_:
        raise ProteinSectionError("nucleus admission must be one Boolean per cell")
    markers = {
        marker: candidate_calls(section, marker, values, eligible,
                                shift_scale=shifts.get(marker))
        for marker, values in sorted(measurements.items())
    }
    evaluated = evaluate_rules_batch(source, {marker: item.states
                                              for marker, item in markers.items()})
    # Even a rule with no marker clauses cannot type a cell without a nucleus.
    states = {name: np.where(eligible, values, -1).astype(np.int8)
              for name, values in evaluated.states.items()}
    rules = BatchRuleEvaluation(states, {name: int(np.count_nonzero(values == -1))
                                         for name, values in states.items()})
    return ProteinSectionEvidence(markers, rules,
                                  source.missing_by_phenotype(measurements))
