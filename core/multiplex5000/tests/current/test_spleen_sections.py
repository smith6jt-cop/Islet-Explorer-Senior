"""Spleen section keys: ``hdl090spl`` names HANDEL donor HDL090's spleen section.

The spleen PhenoCycler cohort's donors are HANDEL ids (``HDL090``), not bare digits, and its
sections must resolve to the spleen ROI rather than the pancreas default.
"""

from __future__ import annotations

import pytest

from multiplex5000.integration.tissues import for_roi
from multiplex5000.sections import SectionError, parse


@pytest.mark.parametrize("key", ["hdl090spl", "HDL090spl", "HDL090SPL"])
def test_hdl_spleen_key(key):
    s = parse(key)
    assert s.key == "hdl090spl"
    assert s.donor_id == "HDL090"          # joins on the donor's own spelling
    assert s.roi == "spl"
    assert not s.on_xenium_slide


def test_spleen_roi_is_the_spleen_tissue():
    assert for_roi(parse("hdl011spl").roi) == "spleen"


def test_bare_hdl_key_is_refused_not_pancreas():
    with pytest.raises(SectionError, match="region token"):
        parse("hdl090")


@pytest.mark.parametrize("key", ["abc090spl", "hd090spl", "hdlspl"])
def test_only_the_hdl_prefix_is_accepted(key):
    with pytest.raises(SectionError):
        parse(key)


@pytest.mark.parametrize("key,donor,roi", [("6539", "6539", "panc"), ("6539pln2", "6539", "pln_2"),
                                           ("6539xpln", "6539", "pln_1")])
def test_numeric_keys_unchanged(key, donor, roi):
    s = parse(key)
    assert (s.donor_id, s.roi) == (donor, roi)
