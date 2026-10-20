"""One file checksum: phenocycler/artifacts.py's sha256_file, plus the named copies that must stay.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_one_file_checksum.py -q

Measured 2026-09-24: seventeen helpers across the tree hashed a whole file with sha256.  Sixteen
returned exactly what ``sha256_file`` returns -- every byte, the hex digest, an error on a missing
file -- and one did not (``build_manifest.sha256`` answers "dir" and "skipped_large").  A copy is
replaced only where its file already imports ``phenocycler`` at module level, which proves it runs
where the helper exists: the repo ``.venv``.  The others are kept and named below with the reason,
because scvi-env, the VALIS venv and the Nimbus venv cannot import ``phenocycler``, and a file
that never imports it gives no evidence of where it runs.  A new copy lands here and fails.
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]

THE_ONE = ("core/phenocycler/phenocycler/artifacts.py", "sha256_file")

#: Byte-equivalent copies that stay, and why each cannot import the one.
KEPT = {
    ("analysis/xenium/scripts/panc_xenium_release.py", "sha256_of"):
        "the Xenium lane is a separate system and this file never imports phenocycler",
    ("analysis/xenium/scripts/pln_immune_typing_v3.py", "sha256"):
        "the Xenium lane; imports phenocycler only inside one function",
    ("core/phenocycler/scripts/broadtype/arms/cosine_broad8.py", "_sha"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/phenocycler/scripts/broadtype/arms/kronos2_coral.py", "_sha256"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/phenocycler/scripts/broadtype/arms/kronos2_worker.py", "_sha256"):
        "a worker for an external model environment; never imports phenocycler",
    ("core/phenocycler/scripts/broadtype/arms/kronos_clusters.py", "_sha256"):
        "imports phenocycler only inside one function",
    ("core/phenocycler/scripts/broadtype/arms/nimbus_worker.py", "_sha256"):
        "runs under external/.venv-nimbus, which cannot import phenocycler",
    ("core/phenocycler/scripts/broadtype/declarations.py", "_sha"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/phenocycler/scripts/report/pins.py", "sha256_file"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/phenocycler/scripts/restain/islet_annotations.py", "_sha256"):
        "never imports phenocycler, and the restain lane also runs under the VALIS venv",
    ("core/phenocycler/scripts/sameslide/reprocess_6533.py", "digest"):
        "the HiPerGator same-slide lane; never imports phenocycler",
    ("core/phenocycler/scripts/sameslide/workflow.py", "sha"):
        "the HiPerGator same-slide driver; never imports phenocycler",
    ("core/phenocycler/scripts/xenium_celltype/build_phenotyped_base.py", "_sha256"):
        "a scvi-env file; imports phenocycler only inside one function",
    ("core/phenocycler/scripts/xenium_typing/pins.py", "sha256_file"):
        "a scvi-env lane; imports phenocycler only inside one function",
}

#: Shaped like a file checksum, and not the same function.  (artifacts.sha256_sampled_file hashes
#: samples of a large file, not every byte, and the detector below does not take it for one.)
NOT_EQUIVALENT = {
    ("scripts/retire/build_manifest.py", "sha256"):
        "answers 'dir' for a directory and 'skipped_large' above SHA_LIMIT",
}

#: The copies replaced on 2026-09-24: each binds the one to its old local name, because other
#: modules and tests call it through that name (RA._sha256, AX._sha256, RX._sha256).
REPLACED = (
    "core/phenocycler/scripts/broadtype/release.py",
    "core/phenocycler/scripts/restain/run_align.py",
    "core/phenocycler/scripts/restain/run_islet_annotation_export.py",
    "core/phenocycler/scripts/restain/run_registered_export.py",
)


def _hashes_a_whole_file(node: ast.FunctionDef) -> bool:
    """A small function that builds a sha256, feeds it a file's bytes, and returns the digest."""
    source = ast.unparse(node)
    reads = "open(" in source or "read_bytes(" in source
    returns_digest = any(isinstance(n, ast.Return) and n.value is not None
                         and "hexdigest()" in ast.unparse(n.value) for n in ast.walk(node))
    small = sum(1 for _ in ast.walk(node) if isinstance(_, ast.stmt)) <= 12
    return "sha256" in source and reads and returns_digest and small


def _helpers() -> set[tuple[str, str]]:
    files = subprocess.run(["git", "ls-files", "*.py"], cwd=REPO, capture_output=True,
                           text=True, check=True).stdout.split()
    found = set()
    for rel in files:
        if "/tests/" in rel or rel.startswith("tests/"):
            continue
        try:
            tree = ast.parse((REPO / rel).read_text())
        except (SyntaxError, UnicodeDecodeError, FileNotFoundError):
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _hashes_a_whole_file(node):
                found.add((rel, node.name))
    return found


def test_every_file_checksum_is_the_one_or_named_with_its_reason():
    found = _helpers()
    assert THE_ONE in found
    unnamed = found - {THE_ONE} - set(KEPT) - set(NOT_EQUIVALENT)
    assert not unnamed, (
        f"a file checksum no one has looked at: {sorted(unnamed)}. Import "
        f"phenocycler.artifacts.sha256_file if the file runs in the repo .venv, or name it in "
        f"KEPT with the reason it cannot")
    vanished = (set(KEPT) | set(NOT_EQUIVALENT)) - found
    assert not vanished, f"named copies that no longer exist; take them off the list: {vanished}"


def test_the_replaced_copies_use_the_one_under_their_old_name():
    for rel in REPLACED:
        tree = ast.parse((REPO / rel).read_text())
        assert not any(isinstance(n, ast.FunctionDef) and n.name == "_sha256" for n in tree.body), rel
        assert any(isinstance(n, ast.ImportFrom) and n.module == "phenocycler.artifacts"
                   and any(a.name == "sha256_file" and a.asname == "_sha256" for a in n.names)
                   for n in tree.body), rel
