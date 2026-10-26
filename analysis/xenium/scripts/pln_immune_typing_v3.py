#!/usr/bin/env python3
"""Compatibility command for the package-owned pLN Xenium typing worker."""

from multiplex5000.tissues import lymph_node_xenium_worker as _worker
from multiplex5000.tissues.lymph_node_xenium_worker import *  # noqa: F403

if __name__ == "__main__":
    _worker.main()
