"""One file checksum: phenocycler/artifacts.py's sha256_file, plus the named copies that must stay.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_one_file_checksum.py -q

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

THE_ONE = ("core/multiplex5000/multiplex5000/artifacts.py", "sha256_file")

#: Byte-equivalent copies that stay, and why each cannot import the one.
KEPT = {
    ("analysis/xenium/scripts/panc_xenium_release.py", "sha256_of"):
        "the Xenium lane is a separate system and this file never imports phenocycler",
    ("core/multiplex5000/scripts/broadtype/arms/cosine_broad8.py", "_sha"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/celltype/alignment_block.py", "_sha256"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/broadtype/arms/kronos2_coral.py", "_sha256"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/broadtype/arms/kronos2_worker.py", "_sha256"):
        "a worker for an external model environment; never imports phenocycler",
    ("core/multiplex5000/scripts/broadtype/arms/kronos_clusters.py", "_sha256"):
        "imports phenocycler only inside one function",
    ("core/multiplex5000/scripts/broadtype/arms/nimbus_worker.py", "_sha256"):
        "runs under external/.venv-nimbus, which cannot import phenocycler",
    ("core/multiplex5000/scripts/broadtype/declarations.py", "_sha"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/report/pins.py", "sha256_file"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/restain/islet_annotations.py", "_sha256"):
        "never imports phenocycler, and the restain lane also runs under the VALIS venv",
    ("core/multiplex5000/scripts/sameslide/reprocess_6533.py", "digest"):
        "the HiPerGator same-slide lane; never imports phenocycler",
    ("core/multiplex5000/scripts/sameslide/workflow.py", "sha"):
        "the HiPerGator same-slide driver; never imports phenocycler",
    ("core/multiplex5000/scripts/xenium_celltype/build_phenotyped_base.py", "_sha256"):
        "a scvi-env file; imports phenocycler only inside one function",
    ("core/multiplex5000/scripts/xenium_typing/pins.py", "sha256_file"):
        "a scvi-env lane; imports phenocycler only inside one function",
    ("core/multiplex5000/scripts/celltype/uniform_pln_nodes.py", "sha"):
        "never imports phenocycler; its interpreter is not evidenced",
    ("core/multiplex5000/scripts/typing/run_calibration_arm.py", "sha256"):
        "never imports phenocycler; its interpreter is not evidenced",
    # Recovered byte for byte from the v3 checkout (analysis/cohesive/recovered_v3_producers_20261018.json);
    # replacing a copy would break that record.
    ("analysis/xenium/scripts/materialize_protein_analysis_assignments.py", "digest"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/analyze_protein_literature_endpoints.py", "file_sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/build_embedded_pln_typing.py", "sha256"):
        "recovered from the v3 checkout, which records it by hash",
    ("core/multiplex5000/scripts/typing/build_endocrine_object_composition.py", "file_sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/build_nd_referenced_trajectories.py", "sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/build_pancreas_cell_release.py", "digest"):
        "recovered from the v3 checkout; the release it writes pins this file by hash",
    ("core/multiplex5000/scripts/typing/export_qupath_cell_assignments.py", "sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/prepare_calibration_arms.py", "sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    ("core/multiplex5000/scripts/typing/recompute_typing_scale.py", "sha256"):
        "recovered from the v3 checkout; never imports phenocycler",
    # Recovered byte for byte on 2026-10-20 (analysis/recovered_20261008/RECOVERY.json): a record of
    # what ran on 2026-10-08, not run from the repo; editing a copy would break the sha256 its
    # receipts record.
    **{(f"analysis/recovered_20261008/results_v3/islet_immune_reassessment_20261006/{rel}", "sha"):
       "recovered byte for byte; a current receipt records these bytes"
       for rel in (
           "cell_typing_validation/peer_review_typing_20261008/protein/build_conservative_protein_candidate.py",
           "cell_typing_validation/peer_review_typing_20261008/rna/"
           "combined_PAN_pLN_positive_annotation_writer_preparation_r1/materialize_exact_positive_annotations.py",
           "scripts/audit_protein_candidate_missing_negatives.py",
           "scripts/audit_protein_repaired_primary_evidence.py",
           "scripts/evaluate_corrected_protein_programs.py",
           "scripts/reprocess_protein_own_compartment_candidate.py",
           "scripts/verify_protein_candidate_expression_sources.py",
           "scripts/verify_protein_candidate_primary_integrity.py",
       )},
    **{(f"analysis/recovered_20261008/results_v3/islet_immune_reassessment_20261006/scripts/{name}", "sha"):
       "recovered byte for byte; the cohesive lane's reconcile stage snapshots it by hash"
       for name in ("combine_spatial_hypotheses_v2.py", "evaluate_spatial_H1_H2.py",
                    "evaluate_spatial_H3_H4.py", "evaluate_spatial_RNA_H5_H6.py")},
    ("analysis/recovered_20261008/typing-validation/current_defensible_handoff_20261014_183339/"
     "identity_review/read_current_rna.py", "sha256"):
        "recovered byte for byte; the cohesive lane loads it from the 10-14 handoff folder",
    ("analysis/shared/typing_release_view.py", "sha256"):
        "copied from the v3 checkout; the cohesive lane imports it there and checks its hash",
    ("analysis/xenium/scripts/build_complete_typing_paper_figures.py", "sha"):
        "copied from the v3 checkout with the test that loads it; never imports phenocycler",
}

#: Shaped like a file checksum, and not the same function.  (artifacts.sha256_sampled_file hashes
#: samples of a large file, not every byte, and the detector below does not take it for one.)
NOT_EQUIVALENT = {
    ("scripts/retire/build_manifest.py", "sha256"):
        "answers 'dir' for a directory and 'skipped_large' above SHA_LIMIT",
    ("analysis/cohesive/common.py", "binding"):
        "answers path, digest and size, and refuses a file that changed while it was hashed",
    ("core/multiplex5000/scripts/typing/refine.py", "_log"):
        "builds a donor's refinement log; its one digest is of the script's own source",
}

#: The copies replaced on 2026-09-24: each binds the one to its old local name, because other
#: modules and tests call it through that name (RA._sha256, AX._sha256, RX._sha256).
REPLACED = (
    "core/multiplex5000/scripts/broadtype/release.py",
    "core/multiplex5000/scripts/restain/run_align.py",
    "core/multiplex5000/scripts/restain/run_islet_annotation_export.py",
    "core/multiplex5000/scripts/restain/run_registered_export.py",
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
