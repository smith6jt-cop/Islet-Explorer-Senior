# One entry point for testing, checking and (gated) rebuilding this repository.
#
# pytest -r: the last -r wins, so a bare -rs would hide the failure names that
# core/phenocycler/pyproject.toml:162 asks for with -ra. -rfEs names failures, errors and skips.
#
#   make test            every suite, each run in its own process, continuing past failures
#   make test-fast       the default lane without the data-state guards (-m "not data")
#   make test-data       only the data-state guards, in every suite that has one (red = stale)
#   make check           is every layer/app artifact newer than its inputs? (reads mtimes only)
#   make plan            the rebuild the driver WOULD run, in order; runs nothing
#   make geojson-outline REFUSES unless GO=1; the on-demand 61 GB GeoJSON outline pass, which
#                        writes data/development/typing/geometry (DONORS="115 6374" to limit)
#   make rebuild         REFUSES unless GO=1; then build_layers.py --run (skips fresh nodes)
#   make env-check       the .venv resolves phenocycler to THIS checkout and no script names another
#   make lock            regenerate requirements.lock from the .venv (pip freeze, editables excluded)
#
# Interpreters are spelled out because `python` and `Rscript` on PATH are not the ones this
# repository is built against: /usr/bin/python3 lacks pyarrow and duckdb, /usr/bin/Rscript is
# 4.6.1 without the app's library. analysis/xenium runs in the scvi-env conda environment.
#
# The ten script lanes run one at a time on purpose: each does sys.path.insert(0, <its own dir>)
# and same-named modules across lanes shadow each other when collected together
# (.github/workflows/phenocycler-ci.yml, job `scripts`).

SHELL   := /bin/bash
REPO    := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PY      := $(REPO)/.venv/bin/python
CONDA   := /home/smith6jt/miniconda3/bin/conda
RSCRIPT := /opt/R/4.5.1/bin/Rscript
PYENV   := PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
PYTEST  := $(PYENV) $(PY) -m pytest -p no:cacheprovider
LANES   := common kintsugi report restain typing xenium_celltype broadtype celltype dct xenium_typing
GO      ?=
# Tests that shell out to R find it here even from a clean shell (Rscript is not in /usr/bin).
export PATH := /opt/R/4.5.1/bin:/usr/local/bin:$(PATH)
# R in the C locale rejects a non-ASCII regex the app's tests use; a clean shell has no LANG.
export LANG ?= C.UTF-8
export LC_ALL ?= C.UTF-8

.PHONY: help test _test-all test-core test-current test-integration test-research test-fast test-data test-lanes test-xenium test-follicle \
        test-r test-r-parse test-r-app test-r-viewer test-r-panels test-scripts check plan geojson-outline rebuild env-check lock

help:
	@sed -n '6,15p' $(MAKEFILE_LIST) | sed 's/^# \{0,3\}//'

# `make test` keeps going after a failing suite so every count is reported, then exits non-zero
# if any suite failed (that is what -k does for independent targets).
test:
	@$(MAKE) --no-print-directory -k _test-all

_test-all: test-core test-lanes test-xenium test-follicle test-r test-scripts

# Three targets, not one recipe: make -k continues past a failing target but not past a failing
# line inside a recipe, and the default lane is red on stale artifacts today.
test-core: test-current test-integration test-research

test-current:
	@echo "== core/phenocycler tests/current"
	cd core/phenocycler && $(PYTEST) -q -rfEs tests/current

test-integration:
	@echo "== core/phenocycler tests/integration"
	cd core/phenocycler && $(PYTEST) -q -rfEs tests/integration

test-research:
	@echo "== core/phenocycler tests/research"
	cd core/phenocycler && $(PYTEST) -q -rfEs tests/research

test-fast:
	cd core/phenocycler && $(PYTEST) -q -rfEs -m "not data" tests/current

# The guards live in tests/current, in five of the script lanes and in analysis/xenium, so
# this runs `-m data` in each place the default `-m "not data"` deselects them: every lane,
# one at a time for the reason given at test-lanes, and analysis/xenium under scvi-env with
# its own rootdir. -k semantics by hand, as in test-lanes: every suite runs, then one exit
# status. pytest exits 5 when a selection is empty, which is what `-m data` is in a lane
# that has no guard today; that is not a failure, so 5 is passed over and every other
# non-zero status is kept. The loop stays over all of $(LANES) so a guard added to any lane
# is run here without editing this file.
test-data:
	@status=0; \
	echo "== core/phenocycler tests/current (-m data)"; \
	(cd core/phenocycler && $(PYTEST) -q -rfEs -m data tests/current) || status=1; \
	for lane in $(LANES); do \
	  echo "== lane $$lane (-m data)"; \
	  (cd core/phenocycler && $(PYTEST) -q -rfEs -m data scripts/$$lane/tests); rc=$$?; \
	  if [ $$rc -ne 0 ] && [ $$rc -ne 5 ]; then status=1; fi; \
	done; \
	echo "== analysis/xenium (scvi-env, -m data)"; \
	(cd analysis/xenium && $(PYENV) $(CONDA) run --no-capture-output -n scvi-env python -m pytest -p no:cacheprovider -rfEs -m data tests) || status=1; \
	exit $$status

test-lanes:
	@status=0; for lane in $(LANES); do \
	  echo "== lane $$lane"; \
	  (cd core/phenocycler && $(PYTEST) -q -rfEs scripts/$$lane/tests) || status=1; \
	done; \
	echo "== lane xenium_celltype/hpg"; \
	(cd core/phenocycler && $(PYTEST) -q -rfEs scripts/xenium_celltype/hpg/tests) || status=1; \
	exit $$status

# analysis/xenium/pyproject.toml already sets addopts = "-q"; a second -q would hide the summary.
test-xenium:
	@echo "== analysis/xenium (scvi-env)"
	cd analysis/xenium && $(PYENV) $(CONDA) run --no-capture-output -n scvi-env python -m pytest -p no:cacheprovider -rfEs tests

test-follicle:
	@echo "== analysis/follicle"
	cd analysis/follicle && $(PYTEST) -q tests

# Four targets, not one recipe, for the reason given at test-core.
test-r: test-r-parse test-r-app test-r-viewer test-r-panels

test-r-parse:
	@echo "== apps/islet parse check"
	$(RSCRIPT) apps/islet/tests/parse_check.R

test-r-app:
	@echo "== apps/islet testthat"
	cd apps/islet && $(RSCRIPT) tests/run_tests.R

test-r-viewer:
	@echo "== scripts/test_viewer_helpers.R"
	$(RSCRIPT) scripts/test_viewer_helpers.R

# tests/testthat.R calls test_check(), which looks for ./testthat relative to the cwd -- the
# layout R CMD check runs it from. From core/panels itself it finds no tests and reports success.
test-r-panels:
	@echo "== core/panels testthat"
	cd core/panels/tests && $(RSCRIPT) testthat.R

# scripts/test_*.py are the two real pytest modules (the two check_*.py beside them are
# hand-rolled check scripts with synthetic inputs, run as scripts, and the R checks live under
# apps/islet/tests). A new test_*.py under scripts/ is picked up by the glob.
test-scripts:
	@echo "== scripts/test_*.py"
	$(PYTEST) -q -rfEs scripts/test_*.py
	@echo "== scripts/check_*.py"
	$(PYENV) $(PY) scripts/check_build_tissuumaps_regions.py
	$(PYENV) $(PY) scripts/check_export_viewer_cell_markers.py

check:
	$(PY) scripts/senior/build_layers.py --check

plan:
	$(PY) scripts/senior/build_layers.py --plan

geojson-outline:
	@if [ "$(GO)" != "1" ]; then \
	  echo "refusing: 'make geojson-outline' streams 61 GB of GeoJSON and writes data/development/typing/geometry." >&2; \
	  echo "Run it deliberately: 'GO=1 make geojson-outline' (DONORS=\"115 6374\" to limit)." >&2; \
	  exit 2; \
	fi
	$(PY) core/phenocycler/scripts/typing/build_nonendocrine_outline.py --include-marginal $(if $(DONORS),--donors $(DONORS),)

rebuild:
	@if [ "$(GO)" != "1" ]; then \
	  echo "refusing: 'make rebuild' writes data/development/typing/app and /data/islet_layers." >&2; \
	  echo "Run 'make plan', post it, get a go, then run 'GO=1 make rebuild'." >&2; \
	  exit 2; \
	fi
	$(PY) scripts/senior/build_layers.py --run

env-check:
	$(PY) scripts/env_check.py

lock:
	$(PY) -m pip freeze --exclude-editable > requirements.lock
	@wc -l requirements.lock
