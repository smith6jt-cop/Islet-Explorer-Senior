#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan (code tier 1):
nothing in the tree imports it, and the inputs named below are absent from this checkout.
It is kept as the record of how that artifact was built, not as a description of anything
current.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method.

Compile individual donor spreadsheets into a single streamlined Excel for the app.

Usage:
  python scripts/compile_donors.py --input donors.zip --output consolidated_streamlined.xlsx
  python scripts/compile_donors.py --input /path/to/donors_dir --output out.xlsx
  python scripts/compile_donors.py --input donors.zip --mapping mapping.yaml --output out.xlsx

Notes:
- Accepts CSV & XLSX. If XLSX has multiple sheets, uses 'Islets' if present or the first sheet.
- Tries to infer donor_id, donor_status, islet_id, and columns for area & markers.
- You can override detection with an optional mapping YAML (see DATA_COMPILATION.md).
"""
import argparse
import io
import zipfile
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import yaml  # optional
except Exception:  # pragma: no cover
    yaml = None


def _normalize(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


ENDO_MARKERS = ["INS", "GCG", "SST"]
IMMUNE_MARKERS = ["CD3E", "CD68", "CD163", "LGALS3"]
PREFERRED_MEASURES = ["pct", "frac", "dens", "count"]
MEASURE_SYNONYMS: Dict[str, List[str]] = {
    "pct": ["pct", "percent", "percentage", "perc", "pctg"],
    "frac": ["frac", "fraction", "prop", "proportion"],
    "dens": ["dens", "density", "per_area", "cells_per_um2", "cells_per_mm2"],
    "count": ["count", "n", "num", "cells"]
}


def find_marker_columns_for(df: pd.DataFrame, region: str, markers: List[str]) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    region_n = _normalize(region)
    cols = list(df.columns)
    cols_norm = [(_normalize(c), c) for c in cols]
    marker_norms = {m: _normalize(m) for m in markers}
    syn2canon = {}
    for canon, syns in MEASURE_SYNONYMS.items():
        for s in syns:
            syn2canon[_normalize(s)] = canon
    for norm, original in cols_norm:
        if region_n not in norm:
            continue
        canon_measure = None
        for syn_norm, canon in syn2canon.items():
            if syn_norm and syn_norm in norm:
                canon_measure = canon
                break
        if canon_measure is None:
            continue
        for m, mn in marker_norms.items():
            if mn and mn in norm:
                out.setdefault(canon_measure, {})
                out[canon_measure].setdefault(m, original)
    return out


def coalesce_area(df: pd.DataFrame, region: str) -> Tuple[Optional[str], pd.Series]:
    cand_exact = [f"islet_area_{region}_um2", f"{region}__islet_area_{region}_um2"]
    for c in cand_exact:
        if c in df.columns:
            return c, pd.to_numeric(df[c], errors="coerce")
    # fallback heuristic: columns containing 'area' & region tokens
    region_n = _normalize(region)
    areas = [c for c in df.columns if ('area' in c.lower() and region_n in _normalize(c))]
    if areas:
        return areas[0], pd.to_numeric(df[areas[0]], errors="coerce")
    # last resort: any area
    any_area = [c for c in df.columns if 'area' in c.lower()]
    if any_area:
        return any_area[0], pd.to_numeric(df[any_area[0]], errors="coerce")
    return None, pd.Series(np.nan, index=df.index)


def infer_col(df: pd.DataFrame, keys: List[str]) -> Optional[str]:
    cols = list(df.columns)
    cols_norm = { _normalize(c): c for c in cols }
    for k in keys:
        nk = _normalize(k)
        for cn, orig in cols_norm.items():
            if nk in cn:
                return orig
    return None


def load_table(path: str) -> Optional[pd.DataFrame]:
    try:
        if path.lower().endswith(".csv"):
            return pd.read_csv(path)
        if path.lower().endswith(".xlsx"):
            xl = pd.ExcelFile(path)
            sheet = 'Islets' if 'Islets' in xl.sheet_names else xl.sheet_names[0]
            return pd.read_excel(path, sheet_name=sheet)
    except Exception:
        return None
    return None


def guess_status_from_name(name: str) -> Optional[str]:
    s = name.lower()
    if 'aab' in s or 'aab+' in s or 'aab_pos' in s:
        return 'Aab+'
    if 't1d' in s or 'type1' in s:
        return 'T1D'
    if 'nd' in s or 'non-diabetic' in s or 'control' in s:
        return 'ND'
    return None


def compile_dir(root: str, mapping: Optional[Dict]=None) -> pd.DataFrame:
    records = []
    for base, _, files in os.walk(root):
        for fn in files:
            if not fn.lower().endswith((".csv",".xlsx")):
                continue
            path = os.path.join(base, fn)
            df = load_table(path)
            if df is None or df.empty:
                continue

            donor_id_col = infer_col(df, ["donor_id","donor","subject","case"]) or 'donor_id'
            islet_id_col = infer_col(df, ["islet_id","islet"]) or 'islet_id'
            donor_status_col = infer_col(df, ["donor_status","status","group","cohort"]) or None

            donor_status = None
            if donor_status_col and donor_status_col in df.columns:
                donor_status = df[donor_status_col].astype(str)
            else:
                donor_status = pd.Series(guess_status_from_name(path), index=df.index, dtype=object)

            # Area & diameter
            area_core_col, area_core = coalesce_area(df, 'core')
            area_union_col, area_union = coalesce_area(df, 'union')

            # Markers
            m_core_endo = find_marker_columns_for(df, 'core', ENDO_MARKERS)
            m_core_imm = find_marker_columns_for(df, 'core', IMMUNE_MARKERS)
            m_band_imm = find_marker_columns_for(df, 'band', IMMUNE_MARKERS)
            m_union_imm = find_marker_columns_for(df, 'union', IMMUNE_MARKERS)

            out = pd.DataFrame()
            out['donor_id'] = df[donor_id_col] if donor_id_col in df.columns else os.path.basename(base)
            out['islet_id'] = df[islet_id_col] if islet_id_col in df.columns else np.arange(len(df))
            out['donor_status'] = donor_status.fillna('NA')

            if area_core_col:
                out['islet_area_core_um2'] = area_core
            if area_union_col is not None and area_union_col != area_core_col:
                out['islet_area_union_um2'] = area_union

            # Copy detected marker columns through as-is (to preserve naming the app understands)
            def copy_map(mm: Dict[str, Dict[str,str]]):
                for meas, mapping in mm.items():
                    for _, col in mapping.items():
                        if col in df.columns:
                            out[col] = pd.to_numeric(df[col], errors='coerce')

            copy_map(m_core_endo)
            copy_map(m_core_imm)
            copy_map(m_band_imm)
            copy_map(m_union_imm)

            records.append(out)

    if not records:
        return pd.DataFrame()
    full = pd.concat(records, axis=0, ignore_index=True)
    # Clean donor_status values to canonical ND/Aab+/T1D/Other
    def canon_status(x: str) -> str:
        s = str(x).strip().lower()
        if s in {"nd","non-diabetic","control","ctrl","non diabetic"}:
            return "ND"
        if "aab" in s or "+" in s:
            return "Aab+"
        if s in {"t1d","type1","type 1","type-1"}:
            return "T1D"
        return x
    full['donor_status'] = full['donor_status'].astype(str).map(canon_status)
    return full


def extract_if_zip(input_path: str) -> str:
    if os.path.isdir(input_path):
        return input_path
    if input_path.lower().endswith('.zip'):
        tdir = tempfile.mkdtemp(prefix='donor_zip_')
        shutil.unpack_archive(input_path, tdir)
        return tdir
    raise FileNotFoundError(f"Unsupported input: {input_path}")


def _parse_region_name(s: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (islet_key, region_label) from strings like 'Islet_64_core'."""
    if not isinstance(s, str):
        return None, None
    m = re.match(r"^(Islet_\d+)_(core|band|union)$", s)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def compile_from_results_zip(zip_path: str) -> pd.DataFrame:
    """Compile from a results.zip containing TSV files like islet_markers_XXXX.tsv, islet_targets_XXXX.tsv, islet_composition_XXXX.tsv, LGALS3_XXXX.tsv"""
    rows = []
    with zipfile.ZipFile(zip_path) as zf:
        names = [n for n in zf.namelist() if n.lower().endswith('.tsv')]
        # Group files by donor token suffix
        by_donor: Dict[str, List[str]] = {}
        for n in names:
            base = os.path.basename(n)
            m = re.match(r"(islet_markers|islet_targets|islet_composition|LGALS3)_([A-Za-z0-9]+)\.tsv$", base)
            if not m:
                continue
            donor = m.group(2)
            by_donor.setdefault(donor, []).append(n)

        for donor, files in by_donor.items():
            # Load tables
            df_targets = df_markers = df_comp = df_lg = None
            for f in files:
                b = os.path.basename(f)
                with zf.open(f) as fh:
                    if b.startswith('islet_targets_'):
                        df_targets = pd.read_csv(fh, sep='\t')
                    elif b.startswith('islet_markers_'):
                        df_markers = pd.read_csv(fh, sep='\t')
                    elif b.startswith('islet_composition_'):
                        df_comp = pd.read_csv(fh, sep='\t')
                    elif b.startswith('LGALS3_'):
                        df_lg = pd.read_csv(fh, sep='\t')

            if df_targets is None or df_comp is None:
                # Require targets (for area) and composition (for mapping islets)
                continue

            # Areas from targets: region_um2 per islet & region
            areas = {}
            for _, r in df_targets.iterrows():
                ikey, region = _parse_region_name(r.get('region'))
                if not ikey or not region:
                    continue
                # Grab region_um2 once per (islet, region)
                areas.setdefault(ikey, {})
                if f'{region}_um2' not in areas[ikey]:
                    try:
                        areas[ikey][f'{region}_um2'] = float(r.get('region_um2', np.nan))
                    except Exception:
                        pass

            # Build base rows from composition (has islet_id and name like Islet_n)
            base = pd.DataFrame()
            base['islet_id'] = df_comp.get('islet_id')
            base['islet_key'] = df_comp.get('name')  # 'Islet_n'
            base['donor_id'] = donor
            base['donor_status'] = 'Unknown'
            base['cells_total'] = pd.to_numeric(df_comp.get('cells_total'), errors='coerce')
            # endocrine fractions from *_any / cells_total
            for marker, col in [('INS','Ins_any'), ('GCG','Glu_any'), ('SST','Stt_any')]:
                frac = None
                if col in df_comp.columns and 'cells_total' in df_comp.columns:
                    with np.errstate(divide='ignore', invalid='ignore'):
                        frac = pd.to_numeric(df_comp[col], errors='coerce') / base['cells_total']
                base[f'core__frac_{marker}'] = frac

            # Attach areas
            def get_area(row, region_label):
                d = areas.get(row['islet_key'], {})
                return d.get(f'{region_label}_um2', np.nan)
            base['islet_area_core_um2'] = base.apply(lambda r: get_area(r, 'core'), axis=1)
            base['islet_area_band_um2'] = base.apply(lambda r: get_area(r, 'band'), axis=1)
            base['islet_area_union_um2'] = base.apply(lambda r: get_area(r, 'union'), axis=1)

            # Markers (immune/myeloid) from islet_markers + LGALS3
            markers_df = None
            if df_markers is not None:
                markers_df = df_markers.copy()
            if df_lg is not None:
                markers_df = pd.concat([markers_df, df_lg], axis=0, ignore_index=True) if markers_df is not None else df_lg.copy()

            wide = pd.DataFrame(index=base.index)
            if markers_df is not None and not markers_df.empty:
                # Normalize and pivot
                def extract_key_region(s):
                    ikey, region = _parse_region_name(s)
                    return ikey, region
                markers_df = markers_df.rename(columns={'region_type':'type'})
                markers_df['islet_key'] = markers_df['region'].apply(lambda s: extract_key_region(s)[0])
                markers_df['region_label'] = markers_df['region'].apply(lambda s: extract_key_region(s)[1])
                # keep only rows matching base islets
                markers_df = markers_df[markers_df['islet_key'].isin(base['islet_key'])]
                markers_df['marker'] = markers_df['marker'].astype(str).str.upper()
                for measure, col in [('frac','pos_frac'), ('count','pos_count')]:
                    if col not in markers_df.columns:
                        continue
                    t = markers_df.pivot_table(index='islet_key', columns=['region_label','marker'], values=col, aggfunc='first')
                    # Flatten columns
                    t.columns = [f"{r}__{measure}_{m}" for r,m in t.columns]
                    # align to base by islet_key
                    t = t.reindex(base['islet_key']).reset_index(drop=True)
                    wide = pd.concat([wide, t], axis=1)

            # Object densities & counts from targets (Nerve/Capillary/Lymphatic)
            targ_wide = pd.DataFrame(index=base.index)
            if df_targets is not None and not df_targets.empty:
                tdf = df_targets.copy()
                tdf['islet_key'] = tdf['region'].apply(lambda s: _parse_region_name(s)[0])
                tdf['region_label'] = tdf['region'].apply(lambda s: _parse_region_name(s)[1])
                tdf = tdf[tdf['islet_key'].isin(base['islet_key'])]
                tdf['class'] = tdf['class'].astype(str)
                # area density
                if 'area_density' in tdf.columns:
                    td = tdf.pivot_table(index='islet_key', columns=['region_label','class'], values='area_density', aggfunc='first')
                    td.columns = [f"{r}__dens_{c}" for r,c in td.columns]
                    td = td.reindex(base['islet_key']).reset_index(drop=True)
                    targ_wide = pd.concat([targ_wide, td], axis=1)
                # counts
                if 'count' in tdf.columns:
                    tc = tdf.pivot_table(index='islet_key', columns=['region_label','class'], values='count', aggfunc='first')
                    tc.columns = [f"{r}__count_{c}" for r,c in tc.columns]
                    tc = tc.reindex(base['islet_key']).reset_index(drop=True)
                    targ_wide = pd.concat([targ_wide, tc], axis=1)

            donor_df = pd.concat([base.reset_index(drop=True), wide.reset_index(drop=True), targ_wide.reset_index(drop=True)], axis=1)
            rows.append(donor_df)

    if not rows:
        return pd.DataFrame()
    full = pd.concat(rows, axis=0, ignore_index=True)
    return full


def _normalize_donor_id(v) -> str:
    try:
        # Handle floats like 6356.0 -> '6356'
        if isinstance(v, float):
            if np.isnan(v):
                return ''
            return str(int(v))
        if isinstance(v, (np.integer, int)):
            return str(int(v))
        s = str(v).strip()
        # Strip trailing .0 if present
        if re.match(r"^\d+\.0$", s):
            return s.split('.')[0]
        return s
    except Exception:
        return str(v)


def apply_donor_status_mapping(df: pd.DataFrame, map_xlsx: str) -> pd.DataFrame:
    """Merge donor_status from a mapping Excel (expects columns like 'Case ID' and 'Donor Status')."""
    try:
        m = pd.read_excel(map_xlsx)
    except Exception:
        return df
    # Try to detect columns
    col_id = None
    for c in m.columns:
        if _normalize(c) in {"caseid","donor","donorid","id","case"}:
            col_id = c
            break
    col_status = None
    for c in m.columns:
        if _normalize(c) in {"donorstatus","status","group","cohort"}:
            col_status = c
            break
    if col_id is None or col_status is None:
        return df
    mm = m[[col_id, col_status]].copy()
    mm.columns = ["donor_id","donor_status_map"]
    mm["donor_id"] = mm["donor_id"].map(_normalize_donor_id)
    # Normalize compiled donor_id
    out = df.copy()
    out["donor_id"] = out["donor_id"].map(_normalize_donor_id)
    out = out.merge(mm, on="donor_id", how="left")
    # Prefer mapped status when available
    out["donor_status"] = out["donor_status_map"].combine_first(out.get("donor_status"))
    out = out.drop(columns=["donor_status_map"], errors="ignore")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Path to donor ZIP or a directory containing donor files')
    ap.add_argument('--output', required=True, help='Path to write consolidated Excel')
    ap.add_argument('--donor-map-xlsx', default=None, help='Optional Excel mapping of donor_id to donor_status (e.g., CODEX_Pancreas_Donors.xlsx)')
    ap.add_argument('--mapping', default=None, help='Optional YAML mapping file to override detection')
    args = ap.parse_args()

    df = pd.DataFrame()
    mapping = None
    if args.input.lower().endswith('.zip'):
        # Prefer specialized results.zip flow if TSVs present
        try:
            df = compile_from_results_zip(args.input)
        except Exception:
            df = pd.DataFrame()
        if df.empty:
            # fallback: extract & treat as directory of CSV/XLSX
            root = extract_if_zip(args.input)
            df = compile_dir(root, mapping)
    else:
        root = extract_if_zip(args.input)
        if args.mapping:
            if yaml is None:
                raise RuntimeError('pyyaml is not installed; install it or omit --mapping')
            with open(args.mapping, 'r') as fh:
                mapping = yaml.safe_load(fh)
        df = compile_dir(root, mapping)

    # Apply donor status mapping if provided
    if args.donor_map_xlsx:
        df = apply_donor_status_mapping(df, args.donor_map_xlsx)

    # Apply built-in mapping override requested by user (ensures only ND/Aab+/T1D)
    FIXED_STATUS = {
        # ND
        '6516': 'ND', '6356': 'ND', '6479': 'ND', '6548': 'ND', '0112': 'ND', '112': 'ND',
        # Aab+
        '6505': 'Aab+', '6450': 'Aab+', '6521': 'Aab+', '6538': 'Aab+', '6549': 'Aab+',
        # T1D
        '6563': 'T1D', '6533': 'T1D', '6534': 'T1D', '6550': 'T1D', '6551': 'T1D',
    }
    if not df.empty:
        df['donor_id'] = df['donor_id'].map(_normalize_donor_id)
        df['donor_status'] = df.apply(
            lambda r: FIXED_STATUS.get(r['donor_id'], r.get('donor_status', 'Unknown')),
            axis=1
        )
    if df.empty:
        raise SystemExit('No donor tables found or compiled result is empty.')

    with pd.ExcelWriter(args.output, engine='openpyxl') as xw:
        df.to_excel(xw, sheet_name='Islets', index=False)
    print(f'Wrote {len(df)} rows to {args.output}')


if __name__ == '__main__':
    main()
