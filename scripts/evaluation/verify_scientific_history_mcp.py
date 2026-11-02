#!/usr/bin/env python3
"""Read, idempotently retry and export a real scientific project through MCP."""
import argparse
import asyncio
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from multiplex5000.parameter_observations import file_hash


async def main(project, checkout):
    original = project / "parameter_observations.json"
    payload = json.loads(original.read_text())
    # Use a populated, nested scientific record, not the synthetic protocol fixture.
    item = next(r for r in reversed(payload["records"]) if r["observation"]["operation"] == "rna_typing")
    expected_hash = file_hash(original)
    global_db = Path.home() / ".kintsugi/parameter_learning.db"
    global_before = file_hash(global_db) if global_db.is_file() else None
    params = StdioServerParameters(command=sys.executable, args=[
        "-m", "multiplex5000.parameter_observations", "--project", str(project),
        "--kintsugi", str(checkout), "serve"], env=dict(os.environ))
    async with stdio_client(params) as (read, write):  # noqa: SIM117 -- session uses opened pipes
        async with ClientSession(read, write) as session:
            await session.initialize()
            retry = await session.call_tool("record_parameter_observation", {"observation": item["observation"]})
            assert not retry.isError
            result = json.loads(retry.content[0].text)
            assert result["duplicate"] and result["record_id"] == item["record_id"]
            history = await session.call_tool("get_parameter_observation_history", {"observation_id": item["observation"]["observation_id"]})
            assert not history.isError
            assert any(item["observation"]["observation_id"] in c.text for c in history.content if c.type == "text")
            exported = await session.call_tool("export_parameter_observations", {"relative_path": "mcp_verified_export.json"})
            assert not exported.isError
            result = json.loads(exported.content[0].text)
            assert result["records"] == len(payload["records"]) and result["sha256"] == expected_hash
    global_after = file_hash(global_db) if global_db.is_file() else None
    assert global_before == global_after
    report = {"status": "pass", "checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "project": str(project), "scientific_observation": item["observation"]["observation_id"],
              "retry_duplicate": True, "read": "pass", "export_sha256": expected_hash,
              "records": len(payload["records"]), "global_before": global_before, "global_after": global_after}
    (project / "scientific_mcp_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--kintsugi", type=Path, required=True)
    args = parser.parse_args()
    asyncio.run(main(args.project.resolve(), args.kintsugi.resolve()))
