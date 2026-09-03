"""
Xenium Analysis Utilities.

Imports are lazy to keep the lightweight submodules (e.g.
``xenium_explorer_export``) usable in environments that do not have the full
scientific stack (scanpy/scvi-tools/squidpy).  Importing
``utils.xenium_explorer_export`` never triggers a scanpy import, but
``from utils import load_xenium_data`` does, on first access.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

_LAZY = {
    # analysis_utils — requires scanpy/squidpy
    "load_xenium_data": "utils.analysis_utils",
    "add_spatial_coordinates": "utils.analysis_utils",
    "calculate_qc_metrics": "utils.analysis_utils",
    "filter_cells_and_genes": "utils.analysis_utils",
    "normalize_and_hvg": "utils.analysis_utils",
    "plot_highly_variable_genes": "utils.analysis_utils",
    "plot_scvi_latent_diagnostic": "utils.analysis_utils",
    "umap_dot_size": "utils.analysis_utils",
    "high_contrast_palette": "utils.analysis_utils",
    "lineage_subsets_for": "utils.analysis_utils",
    "render_celltype_subset_panels": "utils.analysis_utils",
    "render_major_lineage_panel": "utils.analysis_utils",
    "run_standard_workflow": "utils.analysis_utils",
    "export_to_csv": "utils.analysis_utils",
    "create_summary_report": "utils.analysis_utils",
    # xenium_explorer_export — pure python + zarr
    "export_for_xenium_explorer": "utils.xenium_explorer_export",
    "export_groups_to_csv": "utils.xenium_explorer_export",
    "export_groups_to_zarr": "utils.xenium_explorer_export",
    "generate_color_palette": "utils.xenium_explorer_export",
}


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        module = importlib.import_module(_LAZY[name])
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'utils' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(list(_LAZY) + list(globals())))


if TYPE_CHECKING:  # help static analyzers see the public names
    from .analysis_utils import (  # noqa: F401
        add_spatial_coordinates,
        calculate_qc_metrics,
        create_summary_report,
        export_to_csv,
        filter_cells_and_genes,
        high_contrast_palette,
        lineage_subsets_for,
        load_xenium_data,
        normalize_and_hvg,
        plot_highly_variable_genes,
        plot_scvi_latent_diagnostic,
        render_celltype_subset_panels,
        render_major_lineage_panel,
        run_standard_workflow,
        umap_dot_size,
    )
    from .xenium_explorer_export import (  # noqa: F401
        export_for_xenium_explorer,
        export_groups_to_csv,
        export_groups_to_zarr,
        generate_color_palette,
    )


__all__ = list(_LAZY)
