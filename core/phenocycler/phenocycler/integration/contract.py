"""
The modality-agnostic cell-table contract.

PhenoCycler and Xenium arrive in completely different shapes: hive-partitioned parquet
keyed by a QuPath UUID on one side, an AnnData/SpatialData object keyed by a Xenium cell id
on the other. Rather than teach every downstream module about both, S1a/S1b normalise each
modality into ONE flat parquet schema and everything after that reads only this.

That is what lets ``structures``, ``register``, ``match``, ``grid`` and ``qc`` be written
once instead of twice, and it is what makes them testable without any imaging data — a
synthetic frame with the right columns is a valid input.

Schema::

    cell_id          str      modality-native id. PhenoCycler: QuPath object_id (UUID).
                              Xenium: cell_id from cells.parquet / adata.obs_names.
    donor_id         str      the cross-modality join key. ALWAYS str (the donor metadata
                              workbook stores it numeric; PhenoCycler derives it by regex).
    tissue           str      'pancreas' | 'pancreatic_lymph_node'. Set when the dataset is
                              IMPORTED, not inferred later — see below.
    roi              str      the section within the tissue: 'panc' | 'pln_1' .. 'pln_3'.
                              One donor can have a pancreas and up to three lymph nodes.
    modality         str      'phenocycler' | 'xenium'
    x_um, y_um       float32  centroid in MICRONS in the modality's own frame.
    x_um_reg,        float32  centroid in microns in the FIXED frame, written by S5.
    y_um_reg                  NaN until registration runs.
    area_um2         float32  segmented cell area.
    lineage_native   str      the modality's own broad label, unmodified.
    lineage_common   str      the harmonised label from vocab.py.
    celltype_fine    str      finest available subtype.
    confidence       float32  assign_margin (PhenoCycler) | broad_confidence (Xenium).
    evidence_positive boolean False for cells assigned by a default/fallback rule rather
                              than by positive evidence; null when no identity exists.
                              See below — this one matters.
    region           str      'core' | 'peri' | 'tissue' | 'other' where known.
    structure_id     str      islet id once structures.py has run, else ''.
    feat__<name>     float32  per-cell features. PhenoCycler: '<marker>_norm'.
                              Xenium: identity-core lineage scores.

``tissue`` is a property of the dataset, fixed at import, not something inferred from a
desired cohort downstream. Carrying it as a column lets a stage filter or pool only sections
that actually exist. The current PhenoCycler run contains 20 pancreas sections; Xenium pLN
catalog rows remain unpaired until corresponding measured PhenoCycler sections exist.

It also decides what is *meaningful*. Islets exist in pancreas and not in lymph node, so
structure extraction and islet matching are pancreas-only; the registration-free levels
(donor, niche) and the grid comparison apply to both. Modules consult
``phenocycler.integration.tissues`` rather than assuming pancreas.

``evidence_positive`` is not decoration. PhenoCycler's ``Epithelial`` is a *default sink* —
``lineage.py`` assigns it to every cell that fails all the positive gates (flagged there as
``epi_default``) — while Xenium's ``Exocrine`` is a positive call from an identity-core
score, and Xenium separately flags ``broad_lowsignal`` cells whose lineage came from an
argmax fallback rather than the classifier. Comparing raw composition fractions across that
asymmetry inflates PhenoCycler epithelium against Xenium exocrine for reasons that have
nothing to do with biology. Carrying one bool through the contract lets ``donor.py`` report
``frac_all`` and ``frac_evidence_positive`` side by side and headline the latter.
"""

from __future__ import annotations

import glob
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

MODALITIES = ("phenocycler", "xenium")

FEATURE_PREFIX = "feat__"

# Parquet-level contract identity. Column names alone cannot distinguish a current
# geometry-only table from an old identity table whose labels happened to be blank. Every
# writer stamps both values and every reader verifies them before exposing the frame.
CONTRACT_SCHEMA_VERSION = 1
CONTRACT_PROFILE_GEOMETRY = "geometry"
CONTRACT_PROFILE_IDENTITY = "identity"
CONTRACT_PROFILES = frozenset({CONTRACT_PROFILE_GEOMETRY, CONTRACT_PROFILE_IDENTITY})
CONTRACT_METADATA_SCHEMA_KEY = b"phenocycler_contract_schema"
CONTRACT_METADATA_PROFILE_KEY = b"phenocycler_contract_profile"
CONTRACT_METADATA_IDENTITY_KEY = b"identity_release_id"

#: (column, numpy dtype) for every non-feature column, in canonical order.
CELL_TABLE_COLUMNS: dict[str, str] = {
    "cell_id": "object",
    "donor_id": "object",
    "tissue": "object",
    "roi": "object",
    "modality": "object",
    "x_um": "float32",
    "y_um": "float32",
    "x_um_reg": "float32",
    "y_um_reg": "float32",
    "area_um2": "float32",
    "lineage_native": "object",
    "lineage_common": "object",
    "celltype_fine": "object",
    "confidence": "float32",
    "evidence_positive": "boolean",
    "region": "object",
    "structure_id": "object",
}

#: Columns a producer must supply; the rest are filled with typed nulls.
#:
#: ``lineage_native`` is deliberately NOT here. A cell table that carries only geometry and
#: intensity is a legitimate product — it is what ``runs/<id>/expression`` can supply while no
#: target/reference pair is verified — and the stages that need only coordinates must be able
#: to consume it. Identity-dependent stages gate on :func:`has_identity` instead, so the
#: question asked is "does this table actually carry labels?" rather than "is this stage on a
#: hardcoded list?".
REQUIRED_COLUMNS = ("cell_id", "donor_id", "modality", "x_um", "y_um")
IDENTITY_REQUIRED_COLUMNS = ("lineage_native",)


def has_identity(df: pd.DataFrame) -> bool:
    """True only for an identity profile carrying a declared release id.

    Label-shaped columns are not provenance. Historical caches contain them too, and an
    identity release may legitimately contain only Unavailable rows. The release id is the
    authority; content is validated separately.
    """

    return (
        df.attrs.get("contract_profile") == CONTRACT_PROFILE_IDENTITY
        and bool(str(df.attrs.get("identity_release_id", "")).strip())
    )


def contract_profile(df: pd.DataFrame) -> str:
    """Return the declared profile; undeclared frames are least-privilege geometry."""

    profile = str(df.attrs.get("contract_profile", "")).strip()
    return profile if profile in CONTRACT_PROFILES else CONTRACT_PROFILE_GEOMETRY


class ContractError(ValueError):
    """A cell table violates the contract."""


@dataclass(frozen=True)
class CellTable:
    """A validated cell table plus the identity of what it describes.

    Kept as a thin wrapper over a DataFrame rather than a subclass: pandas subclassing is
    fragile across operations, and every consumer here wants the plain frame.
    """

    df: pd.DataFrame
    modality: str
    donor_id: str
    roi: str = ""
    tissue: str = ""

    @property
    def n_cells(self) -> int:
        return len(self.df)

    @property
    def features(self) -> list[str]:
        return feature_columns(self.df)

    def xy(self, registered: bool = False) -> np.ndarray:
        """(n, 2) float64 centroids, in microns.

        ``registered=True`` returns the fixed-frame coordinates and raises if S5 has not
        run — silently falling back to unregistered coordinates would produce a
        plausible-looking but meaningless overlay, which is exactly the failure mode of the
        prototype this replaces.
        """
        if registered:
            if "x_um_reg" not in self.df.columns or self.df["x_um_reg"].isna().all():
                raise ContractError(
                    f"{self.modality}/{self.donor_id}: registered coordinates requested but "
                    "x_um_reg is empty — run the 'register' stage, then 'grid' (which "
                    "applies the stored transform and writes the registered coordinates)"
                )
            cols = ("x_um_reg", "y_um_reg")
        else:
            cols = ("x_um", "y_um")
        return self.df.loc[:, list(cols)].to_numpy(dtype=np.float64)


def feature_columns(df: pd.DataFrame) -> list[str]:
    """The ``feat__*`` columns, sorted, so feature order never depends on write order."""
    return sorted(c for c in df.columns if c.startswith(FEATURE_PREFIX))


def feature_name(column: str) -> str:
    """``feat__CD31`` -> ``CD31``."""
    return column[len(FEATURE_PREFIX):] if column.startswith(FEATURE_PREFIX) else column


def as_feature(name: str) -> str:
    """``CD31`` -> ``feat__CD31``."""
    return name if name.startswith(FEATURE_PREFIX) else f"{FEATURE_PREFIX}{name}"


def tissue_for_roi(roi: str) -> str:
    """Map an ROI name to its tissue. Thin alias for :func:`tissues.for_roi`."""
    from .tissues import for_roi

    return for_roi(roi)


def normalize_tissue(tissue: str) -> str:
    """Canonicalise a tissue name. Thin alias for :func:`tissues.normalize`."""
    from .tissues import normalize

    return normalize(tissue)


def normalize_donor_id(value) -> str:
    """Coerce a donor id to the canonical string form.

    The donor metadata workbook stores ``6539`` as a number, ``xenium_paths.csv`` as text,
    and PhenoCycler derives it with a regex over the image name. Left alone, ``6539`` and
    ``6539.0`` and ``' 6539'`` are three different join keys and the merge silently produces
    zero paired donors. Everything that touches a donor id goes through here.
    """
    if value is None:
        return ""
    if isinstance(value, float):
        if np.isnan(value):
            return ""
        if float(value).is_integer():
            return str(int(value))
        return str(value)
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    text = str(value).strip()
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def empty_cell_table() -> pd.DataFrame:
    """A zero-row frame with the full contract schema and correct dtypes."""
    out = pd.DataFrame({c: pd.Series(dtype=d) for c, d in CELL_TABLE_COLUMNS.items()})
    out.attrs.update(contract_profile=CONTRACT_PROFILE_GEOMETRY, identity_release_id="")
    return out


def coerce_cell_table(
    df: pd.DataFrame,
    *,
    modality: str,
    donor_id: str,
    profile: str,
    identity_release_id: str = "",
    roi: str = "",
    tissue: str = "",
) -> pd.DataFrame:
    """Fill missing optional columns, cast dtypes, and order columns canonically.

    Producers supply what they have; this fills the rest with typed nulls so every consumer
    can assume the full schema exists. Feature columns are preserved and sorted after the
    fixed columns.
    """
    if modality not in MODALITIES:
        raise ContractError(f"unknown modality {modality!r}; expected one of {MODALITIES}")
    if profile not in CONTRACT_PROFILES:
        raise ContractError(f"unknown contract profile {profile!r}; expected {sorted(CONTRACT_PROFILES)}")
    release_id = str(identity_release_id).strip()
    if profile == CONTRACT_PROFILE_IDENTITY:
        absent = [c for c in IDENTITY_REQUIRED_COLUMNS if c not in df.columns]
        if absent:
            raise ContractError(f"identity table for {modality}/{donor_id} is missing {absent}")
        if not release_id:
            raise ContractError("identity profile requires a non-empty identity_release_id")
    elif release_id:
        raise ContractError("geometry profile cannot carry identity_release_id")

    out = df.copy()
    missing = [c for c in REQUIRED_COLUMNS if c not in out.columns and c not in ("modality",)]
    # donor_id/modality are supplied as arguments, so only complain about genuinely absent data
    missing = [c for c in missing if c not in ("donor_id",)]
    if missing:
        raise ContractError(f"cell table for {modality}/{donor_id} is missing {missing}")

    out["modality"] = modality
    out["donor_id"] = normalize_donor_id(donor_id)
    if "roi" not in out.columns or roi:
        out["roi"] = roi
    # Tissue is declared at import. When a producer does not pass one, derive it from the ROI
    # rather than leaving it blank — a blank tissue would silently drop the row out of every
    # tissue-filtered analysis, which is a far worse failure than a wrong-but-visible label.
    if tissue:
        # Lenient on input, strict on output: a producer may pass 'pLN' or 'lymph node', but
        # what lands in the column is always the canonical identifier every consumer filters
        # on. An unrecognised name is rejected here rather than becoming a third tissue.
        out["tissue"] = normalize_tissue(tissue)
    elif "tissue" not in out.columns or (out["tissue"].astype(str) == "").all():
        out["tissue"] = tissue_for_roi(roi or (out["roi"].iloc[0] if len(out) else ""))

    for col, dtype in CELL_TABLE_COLUMNS.items():
        if col not in out.columns:
            if dtype == "boolean":
                # Evidence is a property of an identity call. Geometry-only rows have no
                # call, hence unavailable evidence rather than an invented True/False.
                out[col] = pd.Series(pd.NA, index=out.index, dtype="boolean")
            elif dtype.startswith("float"):
                out[col] = np.nan
            else:
                out[col] = ""
        if dtype.startswith("float"):
            out[col] = pd.to_numeric(out[col], errors="coerce").astype(dtype)
        elif dtype == "boolean":
            out[col] = out[col].astype("boolean")
        else:
            out[col] = out[col].fillna("").astype(str)

    feats = feature_columns(out)
    for c in feats:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float32")

    # Anything a producer added beyond the contract is kept, not dropped. The contract is a
    # floor that every consumer can rely on, not a ceiling: `immune_subclass`,
    # `endocrine_subtype`, `xenium_sample_key` and `spatial_niche` are all modality-specific
    # provenance that would be expensive to recompute and is cheap to carry.
    extras = [c for c in out.columns if c not in CELL_TABLE_COLUMNS and c not in feats]
    result = out.loc[:, list(CELL_TABLE_COLUMNS) + feats + sorted(extras)].reset_index(drop=True)
    result.attrs.update(contract_profile=profile, identity_release_id=release_id)
    validate_cell_table(result)
    return result


def validate_cell_table(
    df: pd.DataFrame,
    *,
    modality: Optional[str] = None,
    profile: Optional[str] = None,
) -> None:
    """Raise :class:`ContractError` if ``df`` is not a valid cell table.

    Checks the things that actually go wrong in practice: a missing column, coordinates
    that are secretly pixels, a mixed-modality frame, duplicate ids.
    """
    missing = [c for c in CELL_TABLE_COLUMNS if c not in df.columns]
    if missing:
        raise ContractError(f"cell table missing contract columns: {missing}")

    declared_profile = profile or contract_profile(df)
    if declared_profile not in CONTRACT_PROFILES:
        raise ContractError(f"invalid contract profile {declared_profile!r}")
    release_id = str(df.attrs.get("identity_release_id", "")).strip()
    if declared_profile == CONTRACT_PROFILE_IDENTITY and not release_id:
        raise ContractError("identity profile has no identity_release_id provenance")
    if declared_profile == CONTRACT_PROFILE_GEOMETRY and release_id:
        raise ContractError("geometry profile carries an identity_release_id")

    if len(df) == 0:
        return

    mods = set(df["modality"].unique())
    if not mods <= set(MODALITIES):
        raise ContractError(f"cell table has unknown modality values: {mods - set(MODALITIES)}")
    if len(mods) > 1:
        raise ContractError(f"cell table mixes modalities {mods}; keep one per file")
    if modality is not None and mods != {modality}:
        raise ContractError(f"expected modality {modality!r}, got {mods}")

    if df["cell_id"].duplicated().any():
        n = int(df["cell_id"].duplicated().sum())
        raise ContractError(f"cell table has {n} duplicate cell_id values")

    if df["donor_id"].eq("").any():
        raise ContractError("cell table has empty donor_id values")

    if df["tissue"].eq("").any():
        raise ContractError(
            "cell table has rows with no tissue — tissue is declared at import and every row "
            "must carry one, or the row silently vanishes from tissue-filtered analyses")

    if declared_profile == CONTRACT_PROFILE_GEOMETRY:
        if df["lineage_native"].fillna("").astype(str).str.strip().ne("").any():
            raise ContractError(
                "geometry profile carries lineage_native labels without identity provenance"
            )
        if df["evidence_positive"].notna().any():
            raise ContractError(
                "geometry-only cell table carries evidence_positive values — evidence is "
                "unavailable until an identity call exists"
            )

    for col in ("x_um", "y_um"):
        if not np.isfinite(df[col].to_numpy()).any():
            raise ContractError(f"cell table column {col} has no finite values")

    # Both modalities are natively in microns — PhenoCycler from QuPath's 'Centroid X um',
    # Xenium from obsm['spatial']. A tissue section is at most ~4 cm, so a coordinate span
    # far beyond that means someone handed us pixels. Catching it here is cheap; catching it
    # after registration means debugging a silently wrong overlay.
    span = float(np.nanmax(df["x_um"].to_numpy()) - np.nanmin(df["x_um"].to_numpy()))
    if span > 100_000:
        raise ContractError(
            f"x_um spans {span:,.0f} um (>10 cm) — coordinates look like pixels, not microns"
        )


def _contract_metadata(path: Path) -> dict[bytes, bytes]:
    """Read and validate the custom contract metadata on one parquet file."""

    import pyarrow.parquet as pq

    metadata = pq.read_schema(path).metadata or {}
    got_version = metadata.get(CONTRACT_METADATA_SCHEMA_KEY)
    got_profile = metadata.get(CONTRACT_METADATA_PROFILE_KEY)
    got_identity = metadata.get(CONTRACT_METADATA_IDENTITY_KEY, b"")
    if got_version != str(CONTRACT_SCHEMA_VERSION).encode():
        raise ContractError(
            f"{path} has contract schema version {got_version!r}; expected "
            f"{CONTRACT_SCHEMA_VERSION}. Re-export it; legacy integration caches are not "
            "current inputs."
        )
    if got_profile not in {
        CONTRACT_PROFILE_GEOMETRY.encode(),
        CONTRACT_PROFILE_IDENTITY.encode(),
    }:
        raise ContractError(f"{path} has invalid or missing contract profile {got_profile!r}")
    if got_profile == CONTRACT_PROFILE_IDENTITY.encode() and not got_identity:
        raise ContractError(f"{path} identity profile has no identity_release_id metadata")
    if got_profile == CONTRACT_PROFILE_GEOMETRY.encode() and got_identity:
        raise ContractError(f"{path} geometry profile carries identity_release_id metadata")
    return metadata


def write_cell_table(df: pd.DataFrame, path: Path) -> Path:
    """Write one versioned partition, creating parents and validating before persistence."""
    validate_cell_table(df)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    import pyarrow as pa
    import pyarrow.parquet as pq

    table = pa.Table.from_pandas(df, preserve_index=False)
    metadata = dict(table.schema.metadata or {})
    metadata[CONTRACT_METADATA_SCHEMA_KEY] = str(CONTRACT_SCHEMA_VERSION).encode()
    metadata[CONTRACT_METADATA_PROFILE_KEY] = contract_profile(df).encode()
    metadata[CONTRACT_METADATA_IDENTITY_KEY] = str(
        df.attrs.get("identity_release_id", "")
    ).encode()
    pq.write_table(table.replace_schema_metadata(metadata), path)
    return path


def partition_path(root: Path, donor_id: str, roi: str = "") -> Path:
    """Hive path for a donor (+roi), matching the core pipeline's ``donor_id=<id>`` layout.

    The Xenium side adds ``roi=<roi>`` because one donor can have a pancreas and up to three
    pancreatic lymph node regions, and (worse) a single Xenium slide serial can carry both —
    so donor alone does not identify a section.
    """
    base = Path(root) / f"donor_id={normalize_donor_id(donor_id)}"
    if roi:
        base = base / f"roi={roi}"
    return base / "data_0.parquet"


def read_cell_table(
    root: Path,
    donor_id: str,
    roi: str = "",
    *,
    modality: Optional[str] = None,
) -> CellTable:
    """Read one donor (+roi) partition and return a validated :class:`CellTable`."""
    root = Path(root)
    direct = partition_path(root, donor_id, roi)
    if direct.exists():
        files = [direct]
    else:
        pattern = str(root / f"donor_id={normalize_donor_id(donor_id)}" / "**" / "*.parquet")
        files = sorted(glob.glob(pattern, recursive=True))
        if roi:
            files = [f for f in files if f"roi={roi}" in f]
    if not files:
        raise FileNotFoundError(
            f"no cell table for donor {donor_id!r}"
            + (f" roi {roi!r}" if roi else "")
            + f" under {root}"
        )

    file_metadata = [_contract_metadata(Path(f)) for f in files]
    recorded_profiles = {m[CONTRACT_METADATA_PROFILE_KEY].decode() for m in file_metadata}
    release_ids = {m.get(CONTRACT_METADATA_IDENTITY_KEY, b"").decode() for m in file_metadata}
    if len(recorded_profiles) != 1:
        raise ContractError(
            f"cell-table partition mixes contract profiles {sorted(recorded_profiles)} under "
            f"{root}; rebuild it as one profile"
        )
    if len(release_ids) != 1:
        raise ContractError(f"cell-table partition mixes identity releases {sorted(release_ids)}")
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    df.attrs.update(
        contract_profile=next(iter(recorded_profiles)),
        identity_release_id=next(iter(release_ids)),
    )
    validate_cell_table(df, modality=modality)
    actual_profile = contract_profile(df)
    recorded_profile = next(iter(recorded_profiles))
    if actual_profile != recorded_profile:
        raise ContractError(
            f"cell-table metadata says profile {recorded_profile!r}, but contents justify "
            f"{actual_profile!r}: {files[0]}"
        )
    mod = str(df["modality"].iloc[0])
    got_roi = roi or (str(df["roi"].iloc[0]) if "roi" in df.columns and len(df) else "")
    got_tissue = (str(df["tissue"].iloc[0])
                  if "tissue" in df.columns and len(df) else tissue_for_roi(got_roi))
    return CellTable(df=df, modality=mod, donor_id=normalize_donor_id(donor_id),
                     roi=got_roi, tissue=got_tissue)


def discover_partitions(root: Path) -> list[tuple[str, str]]:
    """List ``(donor_id, roi)`` pairs present under ``root``. ``roi`` is '' when unpartitioned."""
    root = Path(root)
    if not root.exists():
        return []
    out: set[tuple[str, str]] = set()
    for donor_dir in sorted(root.glob("donor_id=*")):
        donor = donor_dir.name.split("=", 1)[1]
        roi_dirs = sorted(donor_dir.glob("roi=*"))
        if roi_dirs:
            out.update((donor, d.name.split("=", 1)[1]) for d in roi_dirs)
        elif any(donor_dir.glob("*.parquet")):
            out.add((donor, ""))
    return sorted(out)


def concat_cell_tables(tables: Iterable[CellTable]) -> pd.DataFrame:
    """Concatenate cell tables, keeping the union of feature columns."""
    frames = [t.df for t in tables]
    if not frames:
        return empty_cell_table()
    profiles = {contract_profile(frame) for frame in frames}
    release_ids = {str(frame.attrs.get("identity_release_id", "")) for frame in frames}
    if len(profiles) != 1 or len(release_ids) != 1:
        raise ContractError(
            f"cannot concatenate profiles {sorted(profiles)} / releases {sorted(release_ids)}"
        )
    out = pd.concat(frames, ignore_index=True, sort=False)
    out.attrs.update(
        contract_profile=next(iter(profiles)),
        identity_release_id=next(iter(release_ids)),
    )
    for c in feature_columns(out):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float32")
    validate_cell_table(out)
    return out
