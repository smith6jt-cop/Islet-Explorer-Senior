"""
Export clusters / phenotypes to Xenium Explorer (3.0+).

Xenium Explorer supports two import paths for coloring cells by a categorical
label (cluster id, phenotype, scVI label, etc.):

1. **Cell Groups CSV** -- a plain CSV the user loads via
   ``Cell`` panel -> ``Add cell categorization``.  The first column is
   ``cell_id`` and each remaining column is a categorical grouping.  This is
   the simplest and most portable format.

2. **analysis.zarr.zip** -- the native Xenium analysis bundle produced by
   ``xeniumranger``.  When placed next to the ``experiment.xenium`` file it
   is picked up automatically and its groupings appear in the
   ``Cell`` -> ``Clusters`` dropdown alongside any built-in ones.  The format
   used here matches the one emitted by ``sopa`` / ``spatialdata-io`` so it
   round-trips with those tools.

The module is intentionally self-contained (only numpy/pandas/anndata/zarr) so
it can be imported from any of the analysis notebooks without pulling in
``spatialdata-io``.
"""

from __future__ import annotations

import json
import tempfile
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from anndata import AnnData
except ImportError:  # pragma: no cover - anndata is a hard dep elsewhere
    AnnData = "AnnData"  # type: ignore[assignment]


PathLike = str | Path

# Xenium Explorer expects up to 8-bit RGB hex colors (``#RRGGBB``).
_DEFAULT_PALETTE = "tab20"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_cell_ids(adata: AnnData, cell_id_key: str | None) -> pd.Series:
    """Return the series of cell ids matched to ``adata.obs``.

    Xenium Explorer matches rows by the ``cell_id`` column of the original
    ``cells.parquet`` / ``cells.csv.gz``.  Users frequently store this under a
    column in ``.obs`` (``cell_id``, ``cell_ID``) or leave it as the AnnData
    index.  The caller can override via ``cell_id_key``.
    """
    if cell_id_key is not None:
        if cell_id_key not in adata.obs.columns:
            raise KeyError(
                f"cell_id_key={cell_id_key!r} not found in adata.obs"
            )
        return adata.obs[cell_id_key].astype(str)

    for candidate in ("cell_id", "cell_ID", "CellID", "cell"):
        if candidate in adata.obs.columns:
            return adata.obs[candidate].astype(str)

    return pd.Series(adata.obs_names.astype(str), index=adata.obs_names)


def _as_categorical_series(values: pd.Series) -> pd.Categorical:
    """Coerce any dtype (int cluster ids, strings, pandas Categorical) to
    ``pandas.Categorical`` with deterministic category order."""
    if isinstance(values.dtype, pd.CategoricalDtype):
        return values.cat.remove_unused_categories().values  # type: ignore[return-value]
    if pd.api.types.is_numeric_dtype(values):
        ordered = sorted(values.dropna().unique())
        return pd.Categorical(values.astype(str), categories=[str(c) for c in ordered])
    return pd.Categorical(values.astype(str))


# 20-color categorical fallback (tab20, kept verbatim) used when matplotlib
# is not available. Xenium Explorer accepts 8-bit ``#RRGGBB`` hex strings.
_TAB20_FALLBACK = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
    "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
]


def generate_color_palette(
    categories: Sequence[str],
    palette: str = _DEFAULT_PALETTE,
) -> dict:
    """Return a ``{category: '#RRGGBB'}`` mapping.

    Uses a matplotlib colormap when available (any categorical palette name
    like ``tab10``, ``tab20``, ``Set1`` works) and otherwise falls back to a
    built-in 20-color tab20 sequence cycled to fit ``categories``.
    """
    def _fallback():
        n = len(_TAB20_FALLBACK)
        return {str(cat): _TAB20_FALLBACK[i % n] for i, cat in enumerate(categories)}

    try:
        import matplotlib
        from matplotlib.colors import to_hex
    except ImportError:
        return _fallback()

    # `matplotlib.cm.get_cmap` was deprecated in 3.7 and REMOVED in 3.9; the registry
    # `matplotlib.colormaps` has been available since 3.5.  The old call raised
    # AttributeError, which the ImportError guard above does not catch -- so on any modern
    # matplotlib this helper did not degrade to the fallback, it crashed.  Both environments
    # here are past the removal (repo .venv 3.11.1, scvi-env 3.10.7).
    try:
        cmap = matplotlib.colormaps[palette]
    except (KeyError, AttributeError):
        # An unknown palette name is a caller error worth surfacing, but not worth losing
        # the export over -- the fallback is a real 20-colour categorical sequence.
        return _fallback()
    n = getattr(cmap, "N", 256)
    colors = {}
    for i, cat in enumerate(categories):
        colors[str(cat)] = to_hex(cmap(i % n))
    return colors


# ---------------------------------------------------------------------------
# CSV export (Cell Groups)
# ---------------------------------------------------------------------------

def export_groups_to_csv(
    adata: AnnData,
    group_keys: str | Iterable[str],
    output_path: PathLike,
    cell_id_key: str | None = None,
    rename: Mapping[str, str] | None = None,
) -> Path:
    """Write a Xenium Explorer "Cell Groups" CSV.

    The resulting file has ``cell_id`` as the first column followed by one
    column per entry in ``group_keys``.  Xenium Explorer will create a new
    categorical layer for each column when the CSV is loaded via
    ``Cell`` -> ``Add cell categorization``.

    Parameters
    ----------
    adata
        Annotated data object whose ``.obs`` holds the groupings.
    group_keys
        One or more ``.obs`` column names (e.g. ``"leiden"``, ``"celltype"``,
        ``"phenotype"``).
    output_path
        Destination ``.csv`` file.  Parent directories are created.
    cell_id_key
        Column in ``.obs`` that stores the original Xenium ``cell_id``.  When
        ``None`` common defaults are tried and otherwise the AnnData index is
        used.
    rename
        Optional ``{obs_column: column_name_in_csv}`` mapping, handy when the
        ``.obs`` column is cryptic (``leiden_0.5`` -> ``Clusters``).

    Returns
    -------
    Path
        The written CSV path.
    """
    if isinstance(group_keys, str):
        group_keys = [group_keys]
    group_keys = list(group_keys)

    missing = [k for k in group_keys if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"group_keys not found in adata.obs: {missing}")

    rename = dict(rename or {})

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = {"cell_id": _resolve_cell_ids(adata, cell_id_key).values}
    for key in group_keys:
        col_name = rename.get(key, key)
        data[col_name] = adata.obs[key].astype(str).values

    pd.DataFrame(data).to_csv(output_path, index=False)
    return output_path


# ---------------------------------------------------------------------------
# analysis.zarr.zip export
# ---------------------------------------------------------------------------

def _create_array(group, name: str, data: np.ndarray) -> None:
    """Version-agnostic shim for writing a 1-D array into a zarr group.

    zarr 2 exposes ``create_dataset(name, data=...)`` while zarr 3 moved to
    ``create_array(name, shape=..., dtype=...)`` followed by assignment.
    """
    if hasattr(group, "create_array"):
        arr = group.create_array(name=name, shape=data.shape, dtype=data.dtype)
        arr[...] = data
    else:  # zarr 2
        group.create_dataset(name, data=data)


def _write_cell_groups_zarr(
    store_root,  # zarr.Group
    groupings: Mapping[str, pd.Categorical],
    colors: Mapping[str, Mapping[str, str]],
) -> None:
    """Populate a zarr ``cell_groups`` sub-group in the layout Xenium
    Explorer expects.

    Layout written (mirrors sopa / spatialdata-io):

    ``cell_groups`` (group)
        attrs:
            major_version: 1
            minor_version: 0
            grouping_names: [name0, name1, ...]
            group_names: [[cat00, cat01, ...], [cat10, cat11, ...], ...]
        ``cell_groups/0`` (group)
            attrs:
                name, version, color_palette
            ``indices`` (1-D int32) - concatenated cell indices per category
            ``indptr``  (1-D int32) - start offsets per category (len == n_cats)
        ``cell_groups/1`` ...
    """
    grouping_names: list[str] = []
    group_names: list[list[str]] = []

    for i, (name, cats) in enumerate(groupings.items()):
        codes = np.asarray(cats.codes)
        categories = [str(c) for c in cats.categories]

        # Build an argsort of cell indices grouped by category so that
        # ``indices[indptr[k]:indptr[k+1]]`` yields the cells of category k.
        order = np.argsort(codes, kind="stable")
        sorted_codes = codes[order]
        # Skip any -1 codes (NaN categories) from the output.
        mask = sorted_codes >= 0
        order = order[mask]
        sorted_codes = sorted_codes[mask]

        indptr = np.zeros(len(categories) + 1, dtype=np.int32)
        if order.size:
            counts = np.bincount(sorted_codes, minlength=len(categories))
            indptr[1:] = np.cumsum(counts)

        sub = store_root.create_group(str(i))
        _create_array(sub, "indices", order.astype(np.int32))
        _create_array(sub, "indptr", indptr)

        palette = dict(colors.get(name) or generate_color_palette(categories))
        sub.attrs.update({
            "name": str(name),
            "major_version": 1,
            "minor_version": 0,
            "color_palette": [palette.get(c, "#cccccc") for c in categories],
        })

        grouping_names.append(str(name))
        group_names.append(categories)

    store_root.attrs.update({
        "major_version": 1,
        "minor_version": 0,
        "grouping_names": grouping_names,
        "group_names": group_names,
    })


def export_groups_to_zarr(
    adata: AnnData,
    group_keys: str | Iterable[str],
    output_path: PathLike,
    cell_id_key: str | None = None,
    rename: Mapping[str, str] | None = None,
    colors: Mapping[str, Mapping[str, str]] | None = None,
    palette: str = _DEFAULT_PALETTE,
) -> Path:
    """Write an ``analysis.zarr.zip`` that Xenium Explorer picks up next to
    ``experiment.xenium``.

    The written archive is **position-indexed** -- it records, per category,
    the integer row index into ``adata``.  This means the order of cells in
    ``adata`` must match the order that Xenium Explorer sees (typically the
    ``cells.parquet`` order).  When cells have been filtered/reordered, use
    :func:`export_groups_to_csv` instead, which is matched by ``cell_id``.

    Parameters
    ----------
    adata
        Annotated data object.
    group_keys
        ``.obs`` columns to export.
    output_path
        Destination path.  If it does not end in ``.zarr.zip`` the suffix is
        added.  Convention: ``analysis.zarr.zip``.
    cell_id_key
        Reserved for future index-by-cell-id support; currently unused but
        validated so CSV/zarr callsites share a signature.
    rename
        Optional ``{obs_column: name_in_explorer}`` mapping.
    colors
        Optional ``{grouping_name: {category: '#RRGGBB'}}`` mapping.  Missing
        entries fall back to ``palette``.
    palette
        Matplotlib colormap name used for default colors.

    Returns
    -------
    Path
        The written ``.zarr.zip`` path.
    """
    try:
        import zarr
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "export_groups_to_zarr requires the 'zarr' package"
        ) from exc

    if isinstance(group_keys, str):
        group_keys = [group_keys]
    group_keys = list(group_keys)

    missing = [k for k in group_keys if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"group_keys not found in adata.obs: {missing}")

    # Touch cell_id_key purely to validate callers.
    if cell_id_key is not None:
        _resolve_cell_ids(adata, cell_id_key)

    rename = dict(rename or {})
    colors = dict(colors or {})

    output_path = Path(output_path)
    if not str(output_path).endswith(".zarr.zip"):
        output_path = output_path.with_suffix("").with_suffix(".zarr.zip")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    groupings: dict[str, pd.Categorical] = {}
    for key in group_keys:
        name = rename.get(key, key)
        groupings[name] = _as_categorical_series(adata.obs[key])

    resolved_colors = {
        name: colors.get(name) or generate_color_palette(
            [str(c) for c in cats.categories], palette=palette
        )
        for name, cats in groupings.items()
    }

    # Xenium Explorer opens a v2-format zarr inside a plain zip.  Writing
    # directly through zarr's ZipStore proved unreliable across zarr 2/3
    # (zarr 3 flushes attrs twice, producing duplicate zip entries).  Writing
    # to a temp directory and zipping it up afterwards is simple and
    # version-agnostic.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp) / "analysis.zarr"
        try:  # zarr 3
            root = zarr.open_group(store=str(tmp_dir), mode="w", zarr_format=2)
        except TypeError:  # zarr 2
            root = zarr.open_group(store=str(tmp_dir), mode="w")
        cell_groups = root.create_group("cell_groups")
        _write_cell_groups_zarr(cell_groups, groupings, resolved_colors)

        with zipfile.ZipFile(output_path, mode="w", compression=zipfile.ZIP_STORED) as zf:
            for path in sorted(tmp_dir.rglob("*")):
                if path.is_file():
                    zf.write(path, arcname=path.relative_to(tmp_dir).as_posix())

    return output_path


# ---------------------------------------------------------------------------
# Top-level convenience
# ---------------------------------------------------------------------------

def export_for_xenium_explorer(
    adata: AnnData,
    group_keys: str | Iterable[str],
    output_dir: PathLike,
    sample_name: str,
    cell_id_key: str | None = None,
    rename: Mapping[str, str] | None = None,
    colors: Mapping[str, Mapping[str, str]] | None = None,
    palette: str = _DEFAULT_PALETTE,
    write_csv: bool = True,
    write_zarr: bool = True,
) -> dict:
    """Export one or more groupings in both Xenium Explorer formats.

    This is the function most notebooks should call.  It writes

    - ``<output_dir>/<sample_name>_cell_groups.csv`` for
      ``Cell -> Add cell categorization`` import.
    - ``<output_dir>/<sample_name>/analysis.zarr.zip`` for drag-and-drop
      alongside an ``experiment.xenium`` file.
    - ``<output_dir>/<sample_name>_color_palette.json`` recording the colors
      used so downstream figures can stay consistent.

    Returns a dict of ``{artifact: Path}`` that were written.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(group_keys, str):
        group_keys = [group_keys]
    group_keys = list(group_keys)

    artifacts: dict[str, Path] = {}

    if write_csv:
        csv_path = output_dir / f"{sample_name}_cell_groups.csv"
        artifacts["csv"] = export_groups_to_csv(
            adata,
            group_keys,
            csv_path,
            cell_id_key=cell_id_key,
            rename=rename,
        )

    if write_zarr:
        zarr_dir = output_dir / sample_name
        zarr_dir.mkdir(parents=True, exist_ok=True)
        zarr_path = zarr_dir / "analysis.zarr.zip"
        artifacts["zarr"] = export_groups_to_zarr(
            adata,
            group_keys,
            zarr_path,
            cell_id_key=cell_id_key,
            rename=rename,
            colors=colors,
            palette=palette,
        )

    # Emit the palette regardless so the same colors can be reused in figures.
    rename_map = dict(rename or {})
    colors = dict(colors or {})
    palette_dump: dict[str, dict] = {}
    for key in group_keys:
        name = rename_map.get(key, key)
        cats = _as_categorical_series(adata.obs[key])
        palette_dump[name] = colors.get(name) or generate_color_palette(
            [str(c) for c in cats.categories], palette=palette
        )
    palette_path = output_dir / f"{sample_name}_color_palette.json"
    palette_path.write_text(json.dumps(palette_dump, indent=2))
    artifacts["palette"] = palette_path

    return artifacts


__all__ = [
    "export_for_xenium_explorer",
    "export_groups_to_csv",
    "export_groups_to_zarr",
    "generate_color_palette",
]
