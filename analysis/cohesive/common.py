"""Small, fail-closed file bindings and transactional analysis tasks."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Callable, Mapping


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(path: Path | str):
    return json.loads(Path(path).read_text())


def json_bytes(value) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def write_json(path: Path, value) -> None:
    with path.open("xb") as handle:
        handle.write(json_bytes(value))


def replace_json(path: Path, value) -> None:
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(json_bytes(value))
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


# Left for the desktop, editors and services that share this workstation.
CPU_RESERVE = 8
MEMORY_FRACTION = 0.6


def available_gib() -> float:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 2**20
    raise ValueError("MemAvailable is not reported")


def worker_count(tasks: int, gib_per_task: float) -> int:
    """Processes for independent tasks, bounded by free cores and by memory available now."""
    cpus = max(1, (os.cpu_count() or 1) - CPU_RESERVE)
    memory = int(available_gib() * MEMORY_FRACTION // gib_per_task)
    return max(1, min(tasks, cpus, memory))


def peak_rss_gib() -> float:
    """This process's peak resident memory; exact per task because each task gets its own process."""
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20


def parallel_map(function: Callable, items: list, gib_per_task: float, label: str,
                 cost: list[float] | None = None) -> list:
    """Run function over items in fresh processes, one task per process, results in input order.

    With an estimated cost per item, the costliest start first so the longest task is not the tail.
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    order = sorted(range(len(items)), key=lambda i: -cost[i]) if cost else list(range(len(items)))
    workers = worker_count(len(items), gib_per_task)
    print(f"{label}: {len(items)} tasks on {workers} processes "
          f"({gib_per_task:g} GiB/task budget, {available_gib():.0f} GiB available)", flush=True)
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"),
                             max_tasks_per_child=1) as pool:
        done = list(pool.map(function, [items[i] for i in order]))
    results = [None] * len(items)
    for position, result in zip(order, done):
        results[position] = result
    return results


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def binding(path: Path | str) -> dict:
    path = Path(path).resolve(strict=True)
    before = path.stat()
    require(path.is_file(), f"Expected a file, not a directory: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    after = path.stat()
    require((before.st_ino, before.st_size, before.st_mtime_ns) ==
            (after.st_ino, after.st_size, after.st_mtime_ns),
            f"Input changed while hashing: {path}")
    return {"path": str(path), "sha256": digest.hexdigest(), "bytes": after.st_size}


HASH_THREADS = 16


def threaded(function: Callable, items) -> list:
    """Hash many files at once; hashlib releases the GIL on large blocks. Results in input order."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(HASH_THREADS) as pool:
        return list(pool.map(function, items))


def verify_binding(record: Mapping) -> dict:
    actual = binding(record["path"])
    require(actual["sha256"] == record["sha256"], f"SHA256 mismatch: {record['path']}")
    if "bytes" in record:
        require(actual["bytes"] == record["bytes"], f"Size mismatch: {record['path']}")
    return actual


def unique_bindings(records) -> list[dict]:
    selected = {}
    for record in records:
        path = str(Path(record["path"]).resolve())
        if path in selected:
            require(selected[path]["sha256"] == record["sha256"],
                    f"Contradictory bindings for {path}")
        else:
            selected[path] = dict(record, path=path)
    return [selected[path] for path in sorted(selected)]


@dataclass
class Context:
    inputs: Path
    protocol: Path
    output: Path
    completed: dict[str, Path]

    def task(self, name: str) -> Path:
        require(name in self.completed, f"Required task has not completed: {name}")
        return self.completed[name]

    @property
    def handoff(self) -> dict:
        return read_json(self.inputs)

    @property
    def original_repo(self) -> Path:
        return Path(self.handoff["effective_annotation_reader"]["api"]).resolve().parents[2]


def validate_outputs(directory: Path, receipt: dict) -> None:
    require(receipt["execution_status"] == "complete", f"Incomplete task: {directory}")
    declared = receipt["outputs"]
    require(bool(declared), f"Task has no declared output files: {directory}")
    actual_names = {
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file() and path != directory / "completion.json"
    }
    require(actual_names == {record["relative_path"] for record in declared},
            f"Output set changed: {directory}")
    for record in declared:
        path = (directory / record["relative_path"]).resolve()
        require(path.is_relative_to(directory.resolve()), f"Output escapes task: {path}")
        verify_binding({"path": str(path), "sha256": record["sha256"], "bytes": record["bytes"]})


def run_task(output: Path, name: str, inputs: list[dict], settings: dict,
             build: Callable[[Path], dict]) -> tuple[Path, dict]:
    verified = threaded(verify_binding, unique_bindings(inputs))
    specification = {"task": name, "inputs": verified, "settings": settings}
    fingerprint = hashlib.sha256(json_bytes(specification)).hexdigest()
    parent = output / "tasks" / name
    destination = parent / fingerprint
    if destination.exists():
        receipt = read_json(destination / "completion.json")
        require(receipt["specification"] == specification, f"Receipt mismatch: {destination}")
        validate_outputs(destination, receipt)
        print(f"{name}: reusing verified outputs", flush=True)
        return destination, receipt

    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(dir=parent, prefix=".incomplete-"))
    started = utc_now()
    published = False
    try:
        summary = build(staging)
        require(isinstance(summary, dict), f"{name} must return an explicit result summary")
        # Workers import code from disk when they start; code edited during the build would publish
        # results under bindings they did not run. Small, so always rechecked.
        for record in verified:
            if record["path"].endswith(".py"):
                verify_binding(record)
        files = sorted(path for path in staging.rglob("*") if path.is_file())
        require(bool(files), f"{name} produced no output files")
        require(all(path.resolve().is_relative_to(staging.resolve()) for path in files),
                f"{name} output contains an external symlink")
        outputs = []
        for path in files:
            record = binding(path)
            outputs.append({"relative_path": str(path.relative_to(staging)),
                            "sha256": record["sha256"], "bytes": record["bytes"]})
        receipt = {"execution_status": "complete", "started_UTC": started,
                   "completed_UTC": utc_now(), "specification": specification,
                   "outputs": outputs, "result": summary}
        write_json(staging / "completion.json", receipt)
        require(not destination.exists(), f"Another writer published {destination}")
        staging.rename(destination)
        published = True
        print(f"{name}: published {destination}", flush=True)
        return destination, receipt
    finally:
        if not published:
            import sys
            error = sys.exception()
            write_json(staging / "failure.json", {
                "execution_status": "failed", "started_UTC": started, "failed_UTC": utc_now(),
                "task": name, "error_type": type(error).__name__, "error": str(error),
                "partial_outputs_are_not_reusable": True,
            })
