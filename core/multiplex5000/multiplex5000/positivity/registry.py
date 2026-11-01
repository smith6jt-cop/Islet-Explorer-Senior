"""registry.py -- marker roles and typing constants for one panel and one tissue, read from data, never from code.

Ported from Islet-Explorer scripts/uniform/typing/registry.py at 578767e (see ``positivity.SOURCE_COMMIT``). Changed in
the move: a registry is loaded by an explicit panel and tissue (the CODEX run's ``cohort.csv`` lookup is not carried
over), its files are package data under ``positivity/registry/``, and the loaded registry records its directory so the
labels file is read beside it.

Owner, 2026-10-28: "The lineage and shared lists are not static and will be panel and tissue specific which needs to be
clear for the pipeline as a whole." Every typing module asks this module which markers it may use and how; no marker
list is written into a typing module. A run names its panel and tissue, and the files are

    registry/<panel>__<tissue>.csv              one row per channel: role, compartment, lineage, measurement, ...
    registry/<panel>__<tissue>__constants.csv   one row per constant: value, source, reason

Roles (canonical marker biology for the tissue, owner 2026-10-28):
  defining     a positive value places a cell in the marker's compartment (and, inside it, the lineage)
  shared       canonically on several compartments; adds evidence inside a compartment a defining marker found,
               never places a cell
  state        measured and flagged (e.g. Ki67), never part of an identity
  none         acquired but not used in typing (e.g. ColIV, owner 2026-10-28: "excluded from typing completely")
  never_used   not acquired on every image of the run (owner 2026-10-28: "Don't use markers that are not in the panel
               for every image")
Measurement: one REDSEA compartment per marker for every cell, whole cell, membrane or nucleus (owner, 2026-10-28:
"We can't have two ways of typing. We need to use whole cell and/or membrane and/or nucleus for all measurements.").
Positivity is read on the whole cell for every marker, Ki67 on the nucleus; `differentiation` names the measurement a
marker is compared on when groups compete for a cell (owner, 2026-10-30: "Use membrane for differentiating, but then
the whole cell still needs to be used for positivity").
Lineage: 'Lineage' for a marker that identifies the lineage by itself (CD3e: T), 'Lineage:subtype' for a subtype marker
(INS: Endocrine:Beta, CD8a: T:CD8), several entries for a marker on several lineages (CD4: T:CD4;Myeloid). Typing code
names lineages, never markers; the markers come from these entries.
Level method: `background_method` fixes a marker's level method where the owner decided it (2026-10-28: SMA
negative_p99, CD20 negative_p999); `validation` names the check of such a level ('reported' = shown, not a gate);
`requires` lists markers that must also be positive before the marker's lineage names a cell (CD20 needs CD45).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REGISTRY_DIR = HERE / 'registry'
TYPING_ROLES = ('defining', 'shared')
COMPARTMENT_COLUMN = {'whole_cell': 'Cell', 'membrane': 'Membrane', 'nucleus': 'Nucleus'}


class RegistryError(ValueError):
    """A registry does not fit the images or the tissue it is used for."""


def split(s: str) -> list[str]:
    return [x for x in (s or '').split(';') if x]


@dataclass
class Registry:
    name: str
    table: pd.DataFrame
    constants: dict = field(default_factory=dict)
    constant_sources: dict = field(default_factory=dict)
    directory: Path = REGISTRY_DIR

    @property
    def spec(self) -> tuple[str, str, str]:
        """(panel, tissue, directory): what load() needs to read this registry again (e.g. in a worker process)."""
        panel, tissue = self.name.split('__')
        return panel, tissue, str(self.directory)

    # -- markers ----------------------------------------------------------------------------------------------------
    def role(self, marker: str) -> str:
        return self.table.loc[marker, 'role']

    def markers(self, *roles: str) -> list[str]:
        roles = roles or TYPING_ROLES
        return [m for m in self.table.index if self.table.loc[m, 'role'] in roles]

    def typing_markers(self) -> list[str]:
        """Markers with a typing role (defining or shared)."""
        return self.markers(*TYPING_ROLES)

    def measured_markers(self) -> list[str]:
        """Markers whose positivity is computed: typing markers and states."""
        return self.markers('defining', 'shared', 'state')

    def compartments(self) -> list[str]:
        return list(dict.fromkeys(self.table.loc[self.markers('defining'), 'compartment']))

    def defining(self, compartment: str) -> list[str]:
        t = self.table
        return [m for m in self.markers('defining') if t.loc[m, 'compartment'] == compartment]

    def shared(self, compartment: str) -> list[str]:
        t = self.table
        return [m for m in self.markers('shared') if compartment in split(t.loc[m, 'shared_by'])]

    def lineages(self, marker: str) -> list[str]:
        return split(self.table.loc[marker, 'lineage'])

    def lineage_entries(self, marker: str) -> list[tuple[str, str]]:
        """(lineage, subtype) pairs of a marker, e.g. 'T:CD4;Myeloid' -> [('T', 'CD4'), ('Myeloid', '')]."""
        return [tuple((e.split(':', 1) + [''])[:2]) for e in self.lineages(marker)]

    def broad(self, lineage: str) -> list[str]:
        """Defining markers that identify `lineage` by themselves: their only lineage entry is the bare lineage."""
        return [m for m in self.markers('defining') if self.lineage_entries(m) == [(lineage, '')]]

    def subtypes(self, lineage: str) -> dict[str, str]:
        """{marker: subtype} for the defining markers with an entry 'lineage:subtype', in registry order."""
        return {m: s for m in self.markers('defining') for lin, s in self.lineage_entries(m) if lin == lineage and s}

    def requires(self, marker: str) -> list[str]:
        """Markers that must also be positive before this marker's lineage names a cell (e.g. CD20 needs CD45)."""
        return split(self.table.loc[marker, 'requires']) if 'requires' in self.table.columns else []

    def background_method(self, marker: str) -> str:
        """The level method an owner decision fixed for this marker ('' = chosen among its class's methods)."""
        return self.table.loc[marker, 'background_method'] if 'background_method' in self.table.columns else ''

    def validation(self, marker: str) -> str:
        return self.table.loc[marker, 'validation'] if 'validation' in self.table.columns else ''

    def positivity_class(self, marker: str) -> str:
        return self.table.loc[marker, 'positivity_class']

    def measurement(self, marker: str) -> str:
        return self.table.loc[marker, 'measurement']

    def candidates(self, marker: str) -> list[str]:
        return split(self.table.loc[marker, 'measurement_candidates']) or [self.measurement(marker)]

    def differentiation(self, marker: str) -> str:
        """The measurement that decides between groups competing for a cell, where it is not the positivity one;
        '' otherwise (owner 2026-10-30)."""
        return self.table.loc[marker, 'differentiation'] if 'differentiation' in self.table.columns else ''

    def reading(self, marker: str) -> str:
        """The measurement a marker is compared on: its differentiation measurement, else its positivity one."""
        return self.differentiation(marker) or self.measurement(marker)

    def column(self, marker: str, measurement: str | None = None) -> str:
        """The REDSEA column holding this marker's measurement, e.g. 'Cell__INS'."""
        return f"{COMPARTMENT_COLUMN[measurement or self.measurement(marker)]}__{marker}"

    # -- constants --------------------------------------------------------------------------------------------------
    def const(self, name: str) -> float:
        if name not in self.constants:
            raise KeyError(f'{self.name}: constant {name} is not in the constants file')
        return self.constants[name]

    # -- checks -----------------------------------------------------------------------------------------------------
    def validate(self, image_columns: dict[str, set[str]]) -> None:
        """Fail if a marker with a role is missing on any image, or a never_used marker is marked present everywhere."""
        problems = []
        for m in self.measured_markers():
            for img, cols in image_columns.items():
                if self.column(m) not in cols:
                    problems.append(f'{m} ({self.column(m)}) missing on {img}')
            if self.table.loc[m, 'on_every_image'] != 'yes':
                problems.append(f'{m} has a role but is not acquired on every image')
        if problems:
            raise RegistryError(f'{self.name}: ' + '; '.join(problems))


def load(panel: str, tissue: str, directory: Path | str = REGISTRY_DIR) -> Registry:
    directory = Path(directory)
    name = f'{panel}__{tissue}'
    table = pd.read_csv(directory / f'{name}.csv', keep_default_na=False, dtype=str).set_index('marker')
    consts, sources = {}, {}
    cpath = directory / f'{name}__constants.csv'
    if cpath.exists():
        for r in csv.DictReader(cpath.open()):
            consts[r['name']] = float(r['value']) if r['value'] not in ('', 'pending') else float('nan')
            sources[r['name']] = r['source']
    return Registry(name, table, consts, sources, directory)
