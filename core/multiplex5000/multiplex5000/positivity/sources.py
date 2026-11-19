"""Where the typing method reads a section: Protein Measurement output, or a CODEX run directory.

The CODEX scripts read a run directory (``typing/redsea``, ``typing/context``, ``cells``, ``typing/coverage``); in IES
the same values come from Protein Measurement. A source gives the method what it read from those files, in the same
shape, so the computation is unchanged:

- ``levels_frame``: the analysis-eligible cells of a section with the requested REDSEA columns, their context
  (``context``, ``endocrine_object_class``, ``crossed_frac_*``) and centroids ``x_um``/``y_um``, as ``levels.load``
  built it (REDSEA rows, eligible cells, left-merged context and centroids);
- ``identity_frame``: the same cells with the structure columns ``in_<class>``, as ``identity.load_values`` built it
  (context merged one-to-one); a Protein Measurement source adds ``start_object``, the tissue object a cell starts in
  (a spleen region), kept apart from ``endocrine_object_class``, which ``levels`` reads as 'in an endocrine object';
- ``empty_nulls``: the empty-slide null of each (section, marker), the floor of an edge or gate level
  (``nulls``: glass outside every section on the section's slide);
- ``coverage``: the area-coverage frame of a section (``coverage.compute``), checked against ``redsea_ids``, every
  cell of the section's REDSEA table;
- ``established_calls``: an earlier typing's calls, read only by a labels row ``named_by: established``.

Context is the tissue's (objects a cell is in). A Protein Measurement source reads it from ``context_dir/<section>
.parquet`` when one is given; a section without one has every cell in ordinary tissue and in no object, so
object-based options are not evaluable there rather than invented.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import json
import numpy as np
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

#: the context columns of a cell (CODEX cell_context.py): its context and endocrine object, the fraction of its area
#: inside each structure class's objects, and whether its centroid is in one
OVERLAP_COLUMNS = ('crossed_frac_vessel', 'crossed_frac_nerve', 'crossed_frac_lymphatic')
STRUCTURE_IN_COLUMNS = ('in_Blood vessel', 'in_Nerve', 'in_Lymphatic')
CONTEXT_DEFAULTS = {'context': 'tissue', 'endocrine_object_class': None, 'start_object': None,
                    **{c: 0.0 for c in OVERLAP_COLUMNS}, **{c: False for c in STRUCTURE_IN_COLUMNS}}


class SourceError(ValueError):
    """A section's inputs are missing or do not describe one set of cells."""


@dataclass(frozen=True)
class CodexRunSource:
    """A CODEX Islet-Explorer run directory, read exactly as the pinned scripts read it (the parity check)."""

    run: Path

    def levels_frame(self, tag: str, redsea_columns: list[str], context_columns: list[str]) -> pd.DataFrame:
        run = Path(self.run)
        rs = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet',
                             columns=['object_id', 'qc_analysis_eligible'] + redsea_columns)
        ctx = pd.read_parquet(run / 'typing' / 'context' / f'{tag}.parquet',
                              columns=['object_id', 'context', 'endocrine_object_class'] + context_columns)
        xy = pd.read_parquet(run / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
        xy.columns = ['object_id', 'x_um', 'y_um']
        return rs[rs['qc_analysis_eligible']].merge(ctx, on='object_id', how='left').merge(xy, on='object_id', how='left')

    def identity_frame(self, tag: str, redsea_columns: list[str], context_columns: list[str]) -> pd.DataFrame:
        run = Path(self.run)
        rs = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet',
                             columns=['object_id', 'qc_analysis_eligible'] + redsea_columns)
        ctx = pd.read_parquet(run / 'typing' / 'context' / f'{tag}.parquet',
                              columns=['object_id', 'x_um', 'y_um', 'context', 'endocrine_object_class']
                              + context_columns)
        return rs[rs['qc_analysis_eligible']].merge(ctx, on='object_id', how='left', validate='one_to_one')

    def empty_nulls(self, quantile: float) -> dict:
        t = pd.read_csv(Path(self.run) / 'typing' / 'empty_area_nulls.csv')
        col = f'null_p{quantile:g}'.replace('.', '')
        return {(r.image, r.marker): float(getattr(r, col)) for r in t.itertuples()}

    def coverage(self, tag: str) -> pd.DataFrame:
        return pd.read_parquet(Path(self.run) / 'typing' / 'coverage' / f'{tag}.parquet')

    @contextmanager
    def native_window(self, tag):
        import tifffile

        root = Path(self.run)
        entries = json.loads((root / "project_manifest.json").read_text())["entries"]
        entry, = (e for e in entries if e["image_name"] == tag)
        metadata = root / "qupath/data" / str(entry["entry_id"]) / "server.json"
        channels = [c["name"] for c in json.loads(metadata.read_text())["metadata"]["channels"]]
        data = json.loads((root / "cell_geojson" / f"{tag}.geojson").read_text())
        features = data["features"] if isinstance(data, dict) else data
        features = [{"id": f["id"], "geometry": f["geometry"]} for f in features]
        del data
        with tifffile.TiffFile(entry["image_path"]) as tf:
            level = tf.series[0].levels[0]
            if len(level.shape) != 3 or level.shape[0] != len(channels) or len(set(channels)) != len(channels):
                raise SourceError(f"{tag}: native image and channel map differ")
            mapped = None
            if len(level.pages) == 1 and len(channels) > 1:
                mapped = tifffile.memmap(entry["image_path"], series=0, level=0, mode="r")
            def read(name):
                i = channels.index(name)
                return np.asarray(mapped[i]) if mapped is not None else level.pages[i].asarray()
            yield features, level.shape[-2:], read, channels

    def redsea_ids(self, tag: str) -> pd.Series:
        return pd.read_parquet(Path(self.run) / 'typing' / 'redsea' / f'{tag}.parquet', columns=['object_id'])['object_id']

    def established_calls(self, tag: str, relative: str) -> pd.DataFrame:
        return pd.read_parquet(Path(self.run) / relative / 'calls' / f'{tag}.parquet',
                               columns=['object_id', 'cell_class', 'excluded_reason'])


@dataclass(frozen=True)
class ProteinMeasurementSource:
    """Protein Measurement's published tables (``<stage>/data/{redsea,expression}/donor_id=<section>``).

    ``qc_analysis_eligible`` (expression) is the admission the method reads: a cell without a nucleus, or outside the
    recorded analysis exclusions, is not analysis-eligible and never enters a level, a reference or a call. Centroids
    are expression's ``X_centroid``/``Y_centroid`` (µm). ``context_dir`` and ``coverage_dir`` are optional: what the
    tissue's objects and the area-coverage step supply.

    Empty-slide nulls come from ``nulls.cohort_nulls``: the QuPath cohort ``manifest`` names each section's qptiff,
    ``tissue_dir`` holds every section's tissue outline in that qptiff's level-0 frame. They are cached in (or, when it
    exists, read from) ``empty_nulls_csv``. A source with neither refuses: an edge or gate level without its floor is
    not the method."""

    stage_dir: Path
    context_dir: Path | None = None
    coverage_dir: Path | None = None
    empty_nulls_csv: Path | None = None
    manifest: Path | None = None
    tissue_dir: Path | None = None

    def _partition(self, table: str, tag: str) -> Path:
        files = sorted((Path(self.stage_dir) / 'data' / table / f'donor_id={tag}').glob('*.parquet'))
        if len(files) != 1:
            raise SourceError(f'{tag}: expected one {table} partition, found {len(files)}')
        return files[0]

    def _cells(self, tag: str, redsea_columns: list[str]) -> pd.DataFrame:
        names = pq.ParquetFile(self._partition('redsea', tag)).schema_arrow.names
        missing = sorted(set(redsea_columns) - set(names))
        if missing:
            raise SourceError(f'{tag}: REDSEA has no {missing}')
        rs = pq.read_table(self._partition('redsea', tag), columns=['object_id'] + redsea_columns).to_pandas()
        ex = pq.read_table(self._partition('expression', tag),
                           columns=['object_id', 'qc_analysis_eligible', 'X_centroid', 'Y_centroid']).to_pandas()
        ex = ex.rename(columns={'X_centroid': 'x_um', 'Y_centroid': 'y_um'})
        d = rs.merge(ex, on='object_id', how='left', validate='one_to_one')
        if d['qc_analysis_eligible'].isna().any():
            raise SourceError(f'{tag}: REDSEA cells absent from expression')
        return d[d['qc_analysis_eligible'].astype(bool)]

    def _context(self, tag: str, d: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
        path = Path(self.context_dir) / f'{tag}.parquet' if self.context_dir else None
        d = d.copy()
        if path is not None and path.is_file():
            present = [c for c in columns if c in pq.ParquetFile(path).schema_arrow.names]
            ctx = pd.read_parquet(path, columns=['object_id'] + present)
            d = d.merge(ctx, on='object_id', how='left', validate='one_to_one')
        for c in columns:
            if c not in d:
                d[c] = CONTEXT_DEFAULTS[c]
        return d

    def levels_frame(self, tag: str, redsea_columns: list[str], context_columns: list[str]) -> pd.DataFrame:
        d = self._cells(tag, redsea_columns)
        d = self._context(tag, d, ['context', 'endocrine_object_class'] + context_columns)
        return d[['object_id', 'qc_analysis_eligible'] + redsea_columns + ['context', 'endocrine_object_class']
                 + context_columns + ['x_um', 'y_um']]

    def identity_frame(self, tag: str, redsea_columns: list[str], context_columns: list[str]) -> pd.DataFrame:
        d = self._cells(tag, redsea_columns)
        d = self._context(tag, d, ['context', 'endocrine_object_class', 'start_object'] + context_columns)
        return d[['object_id', 'qc_analysis_eligible'] + redsea_columns + ['x_um', 'y_um', 'context',
                                                                            'endocrine_object_class', 'start_object']
                 + context_columns]

    def empty_nulls(self, quantile: float) -> dict:
        from . import nulls

        cached = Path(self.empty_nulls_csv) if self.empty_nulls_csv else None
        if cached is not None and cached.is_file():
            t = pd.read_csv(cached)
        elif self.manifest is not None and self.tissue_dir is not None:
            t = nulls.cohort_nulls(self.manifest, self.tissue_dir)
            if cached is not None:
                cached.parent.mkdir(parents=True, exist_ok=True)
                t.to_csv(cached, index=False)
        else:
            raise SourceError('no empty-slide nulls: give empty_nulls_csv, or the manifest and tissue_dir to measure '
                              'them')
        col = f'null_p{quantile:g}'.replace('.', '')
        return {(r.image, r.marker): float(getattr(r, col)) for r in t.itertuples()}

    def coverage(self, tag: str) -> pd.DataFrame:
        if self.coverage_dir is None:
            raise SourceError(f'{tag}: no coverage directory; run coverage.compute first')
        return pd.read_parquet(Path(self.coverage_dir) / f'{tag}.parquet')

    def native_window(self, tag):
        if self.manifest is None:
            return None
        from ..artifacts import QuPathCohortManifest
        from ..stages.spleen_typing import native_window

        image, = (i for i in QuPathCohortManifest.read_json(Path(self.manifest)).images if i.donor_id == tag)
        return native_window(image, image.cell_geojson.path, image.qptiff.path)

    def redsea_ids(self, tag: str) -> pd.Series:
        return pq.read_table(self._partition('redsea', tag), columns=['object_id']).to_pandas()['object_id']

    def established_calls(self, tag: str, relative: str) -> pd.DataFrame:
        raise SourceError(f'{tag}: no established typing exists in Protein Measurement output ({relative})')
