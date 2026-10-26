#!/usr/bin/env python3
"""Compatibility command for the package-owned pancreas Xenium worker.

The scientific implementation lives in ``multiplex5000.tissues.pancreas_xenium_worker``.
The existing release driver still invokes this path while it is migrated.
"""

from multiplex5000.tissues import pancreas_xenium_worker as _worker
from multiplex5000.tissues.pancreas_xenium_worker import *  # noqa: F403

if __name__ == "__main__":
    _worker.main()
