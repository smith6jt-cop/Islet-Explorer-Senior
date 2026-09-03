#!/usr/bin/env python3
"""Vendor Xenium subpanels from XeniumPanelExplorer at a PINNED commit — tier-driven.

Reproducibility (CLAUDE.md #6) requires a pinned SHA — never fetch a floating
``main``. Since 2026-06-08 the upstream repo carries a machine-readable
``data/panel_roles.csv`` (one row per (tissue, subpanel): tier, phenotyping_use,
lineage_or_group, n_genes, phenotyping_source, note) + ``data/panel_gene_roles.csv``
(per-gene curation for the mixed panels). This script is now DRIVEN by that
manifest instead of a hand-maintained file list:

  PHASE 1  fetch the two manifests → data/reference/subpanels/{panel_roles,
           panel_gene_roles}.csv  (the snapshot root, sibling to SOURCE.txt).
  PHASE 2  parse panel_roles.csv and fetch, per tissue, every panel whose
           phenotyping_use is identity (broad_mask / embedding_extra),
           module_score (state), or fine_tier — into the local layout
             data/reference/subpanels/{shared,pancreas,thymus}/<stem>.csv
           rows tagged ``exclude`` / ``residual`` / ``custom`` are NOT vendored.
  PHASE 3  (scANVI-redesign 2026-06-08, @076d1192) fetch the per-panel
           IDENTITY_CORE files — the clean, CL-ontology-curated + specificity-
           audited identity gene set per cell-type panel (gene,subtype_hint,
           evidence,note). These REPLACE the fat subpanels in the broad CALL
           (lessons 1+2: a fat-panel mean is not an identity score; z-scoring it
           is the catch-all's silent killer). Layout:
             data/reference/subpanels/identity_core/{shared,pancreas,thymus}/<stem>.csv
  PHASE 4  fetch the provenance side-files: identity_audit_<t>.csv,
           identity_epithelial_general_<t>.csv (snapshot root) + the pinned CL
           ontology (ontology/cl-basic.obo + README.md) so the upstream
           gene-by-gene classification is reproducible offline.

Only the ``gene`` column of each subpanel CSV is consumed downstream by
load_subpanel_genes / state_panels; the identity_core consumer additionally
reads ``subtype_hint`` (fine-tier subtyping). The per-tissue "fat" detection
columns are app-only and ignored. The manifest's ``phenotyping_source`` path
quirk (shared rows: ``subpanels_shared/<f>``; per-tissue rows:
``data/tissues/<t>/subpanels/<f>``) is normalized to both the upstream raw path
and the local vendored relative path.

Usage:
    python scripts/fetch_subpanels.py                 # roles + identity/state/fine_tier + identity_core + extras (both tissues)
    python scripts/fetch_subpanels.py --roles-only    # just refresh the two manifests
    python scripts/fetch_subpanels.py --tissue pancreas
    python scripts/fetch_subpanels.py --dry-run       # print URLs, write nothing
"""
import argparse
import csv
import io
import sys
import urllib.error
import urllib.request
import os
from pathlib import Path

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
DEST = REPO / "data/development/xenium_panels/subpanels"

# Pinned commit — XeniumPanelExplorer @ 2026-06-08, PR #7 merged to main
# (076d1192): per-panel IDENTITY_CORE files (gene-by-gene CL-ontology + canonical
# + specificity audit) + extended panel_gene_roles.csv + identity_audit_<t>.csv +
# the pinned CL ontology. Re-pinned from 0828c565 (which predated the cores).
SHA = "076d1192a25946ac3f9cf941c98202d1f224579a"
RAW = f"https://raw.githubusercontent.com/smith6jt-cop/XeniumPanelExplorer/{SHA}/"

# Machine-readable manifests (PHASE 1). Vendored at the snapshot root.
ROLE_FILES = {
    "panel_roles.csv":      "data/panel_roles.csv",
    "panel_gene_roles.csv": "data/panel_gene_roles.csv",
}

# phenotyping_use values that the consumer phenotyper vendors + reads.
VENDORED_USES = {"broad_mask", "embedding_extra", "module_score", "fine_tier"}
# panels that carry an IDENTITY_CORE upstream (the cell-type identity panels):
# broad_mask (per-lineage) + embedding_extra (16_epithelial_general).
IDENTITY_CORE_USES = {"broad_mask", "embedding_extra"}
TISSUES = ("pancreas", "thymus")

# PHASE 4 provenance side-files (per tissue) + ontology. Upstream repo path ->
# local vendored relative path (under DEST).
EXTRA_FILES_PER_TISSUE = {
    "data/identity_audit_{t}.csv":               "identity_core/identity_audit_{t}.csv",
    "data/identity_epithelial_general_{t}.csv":  "identity_core/identity_epithelial_general_{t}.csv",
}
ONTOLOGY_FILES = {
    "data/ontology/cl-basic.obo": "ontology/cl-basic.obo",
    "data/ontology/README.md":    "ontology/README.md",
}


#: Where the panel taxonomy actually lives in the monorepo.  `core/panels` IS
#: XeniumPanelExplorer, imported as a subtree in consolidation Phase 2 -- so what this script
#: used to fetch over HTTP from a pinned SHA is already tracked, in git, one directory away.
#: Verified 2026-08-28 against `raw.githubusercontent.com/.../076d1192/`: `panel_roles.csv`,
#: `panel_gene_roles.csv` (1,613,132 bytes) and four identity-core files are **byte-identical**.
#:
#: So this reads locally instead of over the network.  One authority (`core/panels`, tracked,
#: pinned by canonical decision #17), one derived cache (`DEST`, regenerable, gitignored,
#: never hand-edited).  The alternative -- a second copy of the same CSVs committed under
#: `analysis/xenium/` -- is the drift this script was written to prevent, and a floating
#: `SHA` constant cannot pin a file the repo already contains.
SOURCE = Path(os.environ.get("XENIUM_PANELS_SRC", REPO / "core" / "panels"))


class PanelSourceError(RuntimeError):
    """The panel taxonomy is not where it should be -- never silently skipped."""


def _read_local(repo_path):
    """Read a XeniumPanelExplorer-relative path from the `core/panels` subtree.

    Returns None -- loudly -- for a path the upstream layout has but this one does not, so a
    missing panel surfaces as a named gap in the summary rather than an empty snapshot.  A
    missing *source tree* is a different failure and raises: it means the subtree is gone,
    not that one file moved.
    """
    if not SOURCE.is_dir():
        raise PanelSourceError(
            f"panel taxonomy source {SOURCE} does not exist. It should be the `core/panels` "
            "subtree (= XeniumPanelExplorer). Set XENIUM_PANELS_SRC to override."
        )
    p = SOURCE / repo_path
    if not p.is_file():
        print(f"  !! not in {SOURCE.name}: {repo_path}")
        return None
    return p.read_text(encoding="utf-8")


#: Kept as the name every phase calls, so the mapping functions below are untouched.
_http_get = _read_local


def _write_source_stamp():
    """Record WHICH `core/panels` this snapshot came from, content-addressed.

    A path is not a pin.  `git rev-parse HEAD:core/panels/data` is the tree object id of the
    taxonomy itself, so it changes if and only if a byte of it changes -- unlike the repo
    HEAD, which moves on every unrelated commit.  `tests/test_repo_health.py` compares this
    stamp against the live tree id, so a snapshot generated before an edit to `core/panels`
    fails a test instead of quietly disagreeing with its own source.
    """
    import subprocess
    try:
        tree = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "HEAD:core/panels/data"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception as e:  # noqa: BLE001 -- an unstamped snapshot must say so, not guess
        tree = f"UNKNOWN ({e!r})"
    DEST.mkdir(parents=True, exist_ok=True)
    (DEST / "SOURCE.txt").write_text(
        "Generated by analysis/xenium/scripts/fetch_subpanels.py -- do not hand-edit.\n"
        "source_tree:  core/panels/data\n"
        f"tree_id:      {tree}\n"
        f"upstream_pin: XeniumPanelExplorer@{SHA}\n"
        "note:         core/panels IS XeniumPanelExplorer, imported as a subtree in\n"
        "              consolidation Phase 2.  Verified byte-identical to the upstream pin\n"
        "              on 2026-08-28 for panel_roles.csv, panel_gene_roles.csv and the\n"
        "              identity cores.  Regenerate with `python scripts/fetch_subpanels.py`.\n"
    )
    print(f"stamped SOURCE.txt -> core/panels/data tree {tree}")


def _local_rel(phenotyping_source, stem):
    """Map a panel_roles `phenotyping_source` to the vendored local relative path.

    shared rows:      'subpanels_shared/<f>'           -> 'shared/<stem>.csv'
    pancreas rows:    'data/tissues/pancreas/subpanels' -> 'pancreas/<stem>.csv'
    thymus rows:      'data/tissues/thymus/subpanels'   -> 'thymus/<stem>.csv'
    """
    s = phenotyping_source.replace("\\", "/")
    if "subpanels_shared/" in s:
        return f"shared/{stem}.csv"
    if "tissues/pancreas/" in s:
        return f"pancreas/{stem}.csv"
    if "tissues/thymus/" in s:
        return f"thymus/{stem}.csv"
    # Fallback: keep the basename under a generic dir (shouldn't happen).
    return f"misc/{stem}.csv"


def _upstream_path(phenotyping_source):
    """Normalize a panel_roles `phenotyping_source` to the upstream repo path."""
    s = phenotyping_source.replace("\\", "/")
    if s.startswith("subpanels_shared/"):
        return "data/" + s
    return s  # pancreas/thymus rows already carry the 'data/tissues/...' prefix


def _identity_core_upstream(phenotyping_source, stem):
    """Upstream repo path of a panel's identity_core file (PHASE 3).

    shared rows   -> data/identity_core/<stem>.csv
    pancreas rows -> data/tissues/pancreas/identity_core/<stem>.csv
    thymus rows   -> data/tissues/thymus/identity_core/<stem>.csv
    """
    s = str(phenotyping_source).replace("\\", "/")
    if "subpanels_shared/" in s or s.startswith("subpanels_shared/"):
        return f"data/identity_core/{stem}.csv"
    if "tissues/pancreas/" in s:
        return f"data/tissues/pancreas/identity_core/{stem}.csv"
    if "tissues/thymus/" in s:
        return f"data/tissues/thymus/identity_core/{stem}.csv"
    return None


def _identity_core_local(phenotyping_source, stem):
    """Local vendored relative path of a panel's identity_core file."""
    s = str(phenotyping_source).replace("\\", "/")
    if "subpanels_shared/" in s or s.startswith("subpanels_shared/"):
        return f"identity_core/shared/{stem}.csv"
    if "tissues/pancreas/" in s:
        return f"identity_core/pancreas/{stem}.csv"
    if "tissues/thymus/" in s:
        return f"identity_core/thymus/{stem}.csv"
    return None


def _write(local_rel, data, dry_run, require_gene_col=True):
    dest = DEST / local_rel
    if dry_run:
        print(f"  [dry] {local_rel}")
        return True
    lines = data.splitlines()
    if not lines:
        print(f"  !! {local_rel}: empty")
        return False
    if require_gene_col and lines[0].split(",")[0].strip().strip('"') != "gene":
        print(f"  !! {local_rel}: first column is not 'gene' ({lines[0][:60]!r})")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    n_rows = sum(1 for ln in lines[1:] if ln.split(",")[0].strip())
    dest.write_text(data)
    print(f"  wrote {local_rel}  ({n_rows} rows)")
    return True


def fetch_roles(dry_run):
    """PHASE 1 — fetch the two manifests; return panel_roles text (or None)."""
    print("PHASE 1 — manifests")
    roles_text = None
    ok = True
    for local_rel, repo_path in ROLE_FILES.items():
        data = _http_get(repo_path)
        if data is None:
            ok = False
            continue
        if local_rel == "panel_roles.csv":
            roles_text = data
        _write(local_rel, data, dry_run, require_gene_col=False)
    return roles_text if ok else None


def fetch_panels(roles_text, tissues, dry_run):
    """PHASE 2 — vendor every identity/state/fine_tier panel from panel_roles.csv."""
    print("\nPHASE 2 — panels (tier-driven from panel_roles.csv)")
    rdr = csv.DictReader(io.StringIO(roles_text))
    seen = {}  # local_rel -> upstream_path (dedupe shared files fetched once)
    for row in rdr:
        if row["tissue"] not in tissues:
            continue
        if row["phenotyping_use"] not in VENDORED_USES:
            continue
        stem = row["subpanel"]
        local_rel = _local_rel(row["phenotyping_source"], stem)
        seen.setdefault(local_rel, _upstream_path(row["phenotyping_source"]))
    ok = fail = 0
    for local_rel, repo_path in sorted(seen.items()):
        data = _http_get(repo_path)
        if data is not None and _write(local_rel, data, dry_run):
            ok += 1
        else:
            fail += 1
    return ok, fail


def fetch_identity_cores(roles_text, tissues, dry_run):
    """PHASE 3 — vendor the per-panel identity_core files for the cell-type
    (broad_mask / embedding_extra) panels. These are the clean identity gene sets
    the broad CALL scores by ABSOLUTE expression (scANVI-redesign 2026-06-08)."""
    print("\nPHASE 3 — identity_core (clean per-panel identity gene sets)")
    rdr = csv.DictReader(io.StringIO(roles_text))
    seen = {}  # local_rel -> upstream_path
    for row in rdr:
        if row["tissue"] not in tissues:
            continue
        if row["phenotyping_use"] not in IDENTITY_CORE_USES:
            continue
        stem = row["subpanel"]
        up = _identity_core_upstream(row["phenotyping_source"], stem)
        loc = _identity_core_local(row["phenotyping_source"], stem)
        if up and loc:
            seen.setdefault(loc, up)
    ok = fail = 0
    for local_rel, repo_path in sorted(seen.items()):
        data = _http_get(repo_path)
        # identity_core files are header `gene,subtype_hint,evidence,note`.
        if data is not None and _write(local_rel, data, dry_run):
            ok += 1
        else:
            fail += 1
    return ok, fail


def fetch_extras(tissues, dry_run):
    """PHASE 4 — provenance side-files (identity_audit / epithelial_general per
    tissue) + the pinned CL ontology. Non-fatal: warn on miss."""
    print("\nPHASE 4 — provenance (identity audits + CL ontology)")
    ok = fail = 0
    for up_tmpl, loc_tmpl in EXTRA_FILES_PER_TISSUE.items():
        for t in tissues:
            data = _http_get(up_tmpl.format(t=t))
            if data is not None and _write(loc_tmpl.format(t=t), data, dry_run,
                                           require_gene_col=False):
                ok += 1
            else:
                fail += 1
    for repo_path, local_rel in ONTOLOGY_FILES.items():
        data = _http_get(repo_path)
        if data is not None and _write(local_rel, data, dry_run,
                                       require_gene_col=False):
            ok += 1
        else:
            fail += 1
    return ok, fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roles-only", action="store_true",
                    help="fetch ONLY panel_roles.csv + panel_gene_roles.csv")
    ap.add_argument("--tissue", choices=TISSUES, default=None,
                    help="restrict PHASE 2 to one tissue (default: both)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print what would be fetched, write nothing")
    args = ap.parse_args()

    tissues = (args.tissue,) if args.tissue else TISSUES
    print(f"source {SOURCE} (upstream pin {SHA[:8]}, verified byte-identical)")
    print(f"dest {DEST}" + ("  [DRY RUN]" if args.dry_run else ""))

    roles_text = fetch_roles(args.dry_run)
    if args.dry_run and roles_text is None:
        # On a dry run we still need the manifest to enumerate PHASE 2.
        roles_text = _http_get("data/panel_roles.csv")
    if roles_text is None:
        print("\n!! could not obtain panel_roles.csv — aborting")
        sys.exit(1)
    if args.roles_only:
        print("\n--roles-only: manifests done")
        sys.exit(0)

    ok2, fail2 = fetch_panels(roles_text, tissues, args.dry_run)
    print(f"\nPHASE 2: {ok2} ok, {fail2} failed")
    ok3, fail3 = fetch_identity_cores(roles_text, tissues, args.dry_run)
    print(f"\nPHASE 3: {ok3} ok, {fail3} failed")
    ok4, fail4 = fetch_extras(tissues, args.dry_run)
    print(f"\nPHASE 4: {ok4} ok, {fail4} failed")
    # PHASE 3 (identity_core) is the broad-call source of truth — a miss there is
    # fatal. PHASE 4 (provenance) misses are warnings only.
    if not args.dry_run and not (fail2 or fail3):
        # Stamp only a COMPLETE snapshot. A stamp on a partial one would assert a
        # provenance the files do not have, which is worse than no stamp at all.
        _write_source_stamp()
    sys.exit(1 if (fail2 or fail3) else 0)


if __name__ == "__main__":
    main()
