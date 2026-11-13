"""Protein identity resolution: owned native interior, one vote, no display ceiling."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import tifffile

from multiplex5000.tissues.interior_dominance import measure_conflicts, score_interior
from multiplex5000.artifacts import FileFingerprint


def test_exact_fifteen_points_uses_all_pixels_and_ties_abstain():
    # A has exactly 15 pp more than B (9 votes vs 6 out of 20); co-positive ties,
    # negatives and the third contender all stay in the denominator.
    x = {"A": [3]*9+[1]*6+[1]*4+[0], "B": [1]*9+[2]*6+[1]*4+[0],
         "C": [0]*20}
    # Set the below-B signal higher (but below B) so A has greater median brightness.
    x["A"][9:15] = [1.5]*6
    got = score_interior(x, dict.fromkeys(x, .5))
    assert got["votes"] == {"A": 9, "B": 6, "C": 0}
    assert got["winner"] == "A" and got["area_lead"] == .15
    x["A"] += [0]; x["B"] += [0]; x["C"] += [0]
    assert score_interior(x, dict.fromkeys(x, .5))["winner"] is None


def test_area_must_agree_with_unique_brightness_winner_and_finite_data():
    x = {"A": [4]*6 + [0]*4, "B": [3.9]*6 + [5]*4, "C": [0]*10}
    # A wins 6/10 pixels but B has a higher median when A's winning values alternate.
    x["A"][:3] = [1]*3
    x["B"][:3] = [.5]*3
    got = score_interior(x, dict.fromkeys(x, .3))
    assert got["area_lead"] == .2 and got["status"] == "brightness_does_not_confirm_area"
    x["C"][-1] = np.nan
    assert score_interior(x, dict.fromkeys(x, .3))["winner"] is None


def test_native_ownership_hole_and_late_neighbour_are_removed_before_erosion(tmp_path, monkeypatch):
    from multiplex5000 import redsea
    side = 28
    tiff = tmp_path/"image.tif"
    with tifffile.TiffWriter(tiff) as writer:
        for val in [100, 160]:
            writer.write(np.full((side, side), val, dtype=np.uint16))
    def ring(a,b,c,d):
        return [[a,b], [c,b], [c,d], [a,d], [a,b]]
    geom = tmp_path/"cells.geojson"
    features = [dict(id="target", geometry=dict(type="Polygon",
                     coordinates=[ring(2,2,22,22), ring(7,7,10,10)])),
                dict(id="later", geometry=dict(type="Polygon", coordinates=[ring(16,0,27,27)]))]
    geom.write_text(json.dumps(dict(type="FeatureCollection",features=features)))
    image = SimpleNamespace(cell_geojson=FileFingerprint.capture(geom), qptiff=FileFingerprint.capture(tiff))
    def pages(path, level):
        f = tifffile.TiffFile(path)
        return f, f.pages, ["A","B"], (side, side), 1.
    monkeypatch.setattr(redsea, "qptiff_channels", pages)
    monkeypatch.setattr(redsea, "validated_channel_markers", lambda image,names:names)
    rows = measure_conflicts(image, {"target":("A","B"),"absent":("A","B")},
                             {"A":(1,0,100), "B":(1,0,100)}, .1)
    # Owned x = 2..15, y = 2..22, less a 4x4 hole. Interior x=4..13,
    # y=4..20 less the hole expanded 2 pixels on all sides (8x8).
    assert rows["target"]["interior_pixels"] == 10*17-8*8
    assert rows["target"]["winner"] == "B" and rows["target"]["medians"]["B"] == 1.6
    assert rows["absent"]["status"] == "object_missing_in_segmentation"
