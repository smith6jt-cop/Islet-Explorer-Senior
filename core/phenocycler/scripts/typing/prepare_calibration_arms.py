#!/usr/bin/env python3
"""Materialize the shifted/unshifted CD8 arms for the REDSEA/raw comparison.

Every unshifted arm starts from the matching shifted alignment and reverses only that
source's per-donor CD8 multiplier. Thus marker source, object keys, non-CD8 values,
geometry and missingness stay identical within each shifted/unshifted pair.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def donor_dirs(root: Path) -> dict[str, Path]:
    found = {p.name.split("=", 1)[1]: p for p in root.glob("donor_id=*") if p.is_dir()}
    return found


def materialize_unshifted(source: Path, shifts_path: Path, dest: Path,
                          force: bool = False) -> list[dict]:
    shifts = json.loads(shifts_path.read_text())
    cd8 = shifts["shifts"]["CD8"]
    donors = list(map(str, shifts["donors"]))
    dirs = donor_dirs(source)
    missing = sorted(set(donors) - set(dirs))
    extra = sorted(set(dirs) - set(donors))
    if missing or extra:
        raise ValueError(f"{source}: donor partitions differ; missing={missing}, extra={extra}")
    dest.mkdir(parents=True, exist_ok=True)
    ledger = []
    for donor in donors:
        inp = dirs[donor] / "data.parquet"
        outdir = dest / f"donor_id={donor}"
        out = outdir / "data.parquet"
        if out.exists() and not force:
            raise FileExistsError(f"{out} already exists; use --force to rebuild")
        multiplier = float(cd8.get("scale", {}).get(donor, 1.0))
        if multiplier == 0 or not multiplier > 0:
            raise ValueError(f"{donor}: invalid CD8 shift multiplier {multiplier}")
        pf = pq.ParquetFile(inp)
        names = pf.schema_arrow.names
        if "CD8" not in names or "object_id" not in names:
            raise ValueError(f"{inp}: missing CD8 or object_id")
        outdir.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(out, pf.schema_arrow, compression="zstd", use_dictionary=True)
        try:
            for batch in pf.iter_batches(batch_size=250_000):
                table = pa.Table.from_batches([batch])
                index = names.index("CD8")
                restored = pc.divide(table.column(index), pa.scalar(multiplier))
                table = table.set_column(index, "CD8", restored)
                writer.write_table(table)
        finally:
            writer.close()
        rows = pq.ParquetFile(out).metadata.num_rows
        if rows != pf.metadata.num_rows:
            out.unlink(missing_ok=True)
            raise AssertionError(f"{donor}: row count changed {pf.metadata.num_rows} -> {rows}")
        ledger.append({"donor_id": donor, "rows": rows, "cd8_multiplier_reversed": multiplier,
                       "source": str(inp), "output": str(out)})
        print(f"{donor}: {rows:,} rows; reversed CD8 multiplier {multiplier:.12g}", flush=True)
    return ledger


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--redsea-shifted", type=Path, required=True)
    ap.add_argument("--redsea-shifts", type=Path, required=True)
    ap.add_argument("--raw-shifted", type=Path, required=True)
    ap.add_argument("--raw-shifts", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    arms = [
        ("redsea", args.redsea_shifted, args.redsea_shifts),
        ("raw", args.raw_shifted, args.raw_shifts),
    ]
    manifest = {"comparison": "registry_REDSEA_vs_all_raw; CD8 shifted_vs_unshifted",
                "unshifted_operation": "divide aligned CD8 by the matching UniFORM donor multiplier",
                "arms": {}}
    for source_name, shifted_root, shifts_path in arms:
        manifest["arms"][f"{source_name}_shifted"] = {
            "root": str(shifted_root), "shifts": str(shifts_path),
            "shifts_sha256": sha256(shifts_path),
            "uniform_commit": json.loads(shifts_path.read_text()).get("uniform_commit")}
        target = args.out / f"{source_name}_unshifted"
        ledger = materialize_unshifted(shifted_root, shifts_path, target, args.force)
        manifest["arms"][f"{source_name}_unshifted"] = {
            "root": str(target), "shifts": str(shifts_path),
            "shifts_sha256": sha256(shifts_path),
            "uniform_commit": json.loads(shifts_path.read_text()).get("uniform_commit"),
            "operation": "only CD8 restored to pre-UniFORM values; all other aligned fields unchanged",
            "donors": ledger}
    args.out.mkdir(parents=True, exist_ok=True)
    out = args.out / "alignment_arms.json"
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(f"\nAlignment arm manifest -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
