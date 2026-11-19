#!/usr/bin/env python3
"""Connect a real stdio MCP client to the project-only recording server."""
import argparse
import asyncio
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from multiplex5000.parameter_observations import CONTEXT_KEYS, SCHEMA, file_hash


async def check(project, checkout):
    global_db = Path.home() / ".kintsugi/parameter_learning.db"
    before = file_hash(global_db) if global_db.is_file() else None
    params = StdioServerParameters(command=sys.executable, args=[
        "-m", "multiplex5000.parameter_observations", "--project", str(project),
        "--kintsugi", str(checkout), "serve"], env=dict(os.environ))
    value = {"schema": SCHEMA, "observation_id": "mcp-protocol-verification", "revision": 1,
             "supersedes": None, "status": "evaluated", "operation": "recording_protocol_check",
             "context": {**dict.fromkeys(CONTEXT_KEYS), "acquisition_ids": [], "section_ids": [],
                         "segmentation_ids": [], "registry_hashes": {}, "panel_hashes": {},
                         "input_receipts": [], "producer_version": "protocol-check/1"},
             "subject": {"kind": "geometry", "identifiers": []}, "settings": [],
             "quality": {"before": None, "after": None}, "outcome": "synthetic protocol verification",
             "evidence": []}
    async with stdio_client(params) as (read, write):  # noqa: SIM117 -- session uses opened pipes
        async with ClientSession(read, write) as session:
            await session.initialize()
            names = sorted(t.name for t in (await session.list_tools()).tools)
            assert names == ["export_parameter_observations", "get_parameter_observation_history",
                             "record_parameter_observation"]
            async def call(name, args):
                result = await session.call_tool(name, args)
                assert not result.isError, result
                return json.loads(result.content[0].text)
            first = await call("record_parameter_observation", {"observation": value})
            retry = await call("record_parameter_observation", {"observation": value})
            assert first["record_id"] == retry["record_id"] and retry["duplicate"]
            value.update(revision=2, supersedes=1, status="blocked", outcome="synthetic revised observation")
            await call("record_parameter_observation", {"observation": value})
            history = await session.call_tool("get_parameter_observation_history", {})
            assert not history.isError
            export = await call("export_parameter_observations", {})
            assert export["records"] == 2
            bad = await session.call_tool("record_parameter_observation", {"observation": {**value, "status": "approved"}})
            assert bad.isError
    after = file_hash(global_db) if global_db.is_file() else None
    assert before == after
    from multiplex5000 import parameter_observations
    result = {"status": "pass", "transport": "stdio", "tools": names, "records": export["records"],
              "checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "server_sha256": file_hash(Path(parameter_observations.__file__)),
              "idempotency": "pass", "revision_history": "pass", "approval_rejected": True,
              "global_db_before": before, "global_db_after": after,
              "export_sha256": export["sha256"]}
    (project / "protocol_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--kintsugi", type=Path, required=True)
    args = parser.parse_args()
    asyncio.run(check(args.project.resolve(), args.kintsugi.resolve()))
