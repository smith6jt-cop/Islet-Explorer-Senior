"""Measure the two physical anchors of the intensity scale, per donor per marker.

Random-access level-0 tile reads (never a whole channel -- that stalls the disk).
For each requested marker it reports:

  ceiling  -- fraction of tissue pixels at or near the 16-bit detector ceiling,
              i.e. how much of the top of the distribution is censored;
  zero     -- the modal background pixel value, which is camera offset plus
              autofluorescence, NOT digitiser zero.

Both are properties of one donor's own image for one marker. Nothing is shared
across donors or markers.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import tifffile

MANIFEST = Path(__file__).resolve().parents[2] / "data/qupath_manifest.json"
CEILING = 65535


def sample_tiles(fh, page, n: int, seed: int):
    """Decode n randomly chosen non-empty tiles from one level-0 page."""
    offsets = page.dataoffsets
    counts = page.databytecounts
    usable = [i for i, c in enumerate(counts) if c > 0]
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(usable))
    tiles = []
    for pos in order:
        idx = usable[int(pos)]
        fh.seek(offsets[idx])
        raw = fh.read(counts[idx])
        try:
            out = page.decode(raw, idx)
        except Exception:  # noqa: BLE001
            continue
        arr = np.asarray(out[0]).squeeze()
        if arr.ndim != 2 or arr.std() < 1.0:
            continue  # empty glass reads flat at the offset
        tiles.append(arr)
        if len(tiles) >= n:
            break
    return tiles


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True)
    ap.add_argument("--markers", nargs="+", required=True)
    ap.add_argument("--tiles", type=int, default=100)
    args = ap.parse_args()

    manifest = json.loads(MANIFEST.read_text())
    image = next(i for i in manifest["images"] if str(i["donor_id"]) == args.donor)
    index = {c["marker"]: c["channel_index"] for c in image["channel_map"]}
    path = image["qptiff"]["path"]

    print(f"donor {args.donor}  {Path(path).name}")
    print(f"{'marker':<16}{'pixels':>12}{'n==0':>13}{'%==0':>8}"
          f"{'n>=65535':>10}{'n>=64000':>10}{'p99':>8}{'p999':>8}{'max':>8}")
    with tifffile.TiffFile(path) as tf:
        pages = tf.series[0].pages
        fh = tf.filehandle
        for marker in args.markers:
            if marker not in index:
                print(f"{marker:<16}  not in this donor's panel")
                continue
            page = pages[index[marker]]
            tiles = sample_tiles(fh, page, args.tiles, seed=7)
            if not tiles:
                print(f"{marker:<16}  no tissue tiles found")
                continue
            v = np.concatenate([t.ravel() for t in tiles]).astype(np.uint16)
            nz = int((v == 0).sum())
            print(f"{marker:<16}{v.size:>12,}{nz:>13,}{100*nz/v.size:>7.2f}%"
                  f"{int((v >= CEILING).sum()):>10,}{int((v >= 64000).sum()):>10,}"
                  f"{int(np.quantile(v,0.99)):>8}{int(np.quantile(v,0.999)):>8}"
                  f"{int(v.max()):>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
