# Claude Code cloud environment — setup, env vars, network

This file is a companion to `CLAUDE.md`. It documents the cloud-environment configuration needed for Claude Code to build and test the Xenium Shiny app. The cloud env never connects to HiPerGator or to the user's lab network — it produces source code and tests against synthetic fixtures, then commits to git. Deployment to HiPerGator is the user's responsibility, not the agent's.

## 1. Setup script

Runs once when the environment is provisioned. Pinned to Ubuntu 22.04 (jammy) or 24.04 (noble); detects the codename and selects the correct PPM URL automatically. Total runtime is roughly 8–15 minutes on first run; subsequent sessions reuse the cached library.

```bash
#!/usr/bin/env bash
# setup.sh — provision Claude Code cloud env for the Xenium Shiny app
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

# -------------------------------------------------------------------
# 1. system deps
# -------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential cmake gfortran pkg-config ca-certificates \
  curl wget gnupg lsb-release software-properties-common \
  git \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libhdf5-dev hdf5-tools \
  libfontconfig1-dev libfreetype6-dev \
  libpng-dev libtiff5-dev libjpeg-dev \
  libcairo2-dev libxt-dev \
  libharfbuzz-dev libfribidi-dev \
  libgit2-dev \
  libglpk-dev libgmp3-dev libmpfr-dev libnlopt-dev \
  zlib1g-dev liblzma-dev libbz2-dev \
  pandoc \
  python3-venv

# -------------------------------------------------------------------
# 2. R 4.4+ (CRAN apt repo)
# -------------------------------------------------------------------
if ! command -v R >/dev/null 2>&1 || \
   awk -v ver="$(R --version 2>/dev/null | head -1 | awk '{print $3}')" \
       'BEGIN{exit !(ver < "4.4.0")}'; then
  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  CODENAME=$(lsb_release -cs)
  echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${CODENAME}-cran40/" \
    | sudo tee /etc/apt/sources.list.d/cran-r.list
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends r-base r-base-dev
fi

# -------------------------------------------------------------------
# 3. Quarto (for the export tab's static report)
# -------------------------------------------------------------------
if ! command -v quarto >/dev/null 2>&1; then
  QV="1.6.40"
  wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QV}/quarto-${QV}-linux-amd64.deb" \
    -O /tmp/quarto.deb
  sudo dpkg -i /tmp/quarto.deb
  rm /tmp/quarto.deb
fi

# -------------------------------------------------------------------
# 4. Site-wide ~/.Rprofile — Posit Public Package Manager binaries.
#    PPM URL format is /cran/__linux__/<codename>/latest; HTTPUserAgent
#    must be set for binary delivery.
# -------------------------------------------------------------------
CODENAME=$(lsb_release -cs)
mkdir -p "$HOME"
cat > "$HOME/.Rprofile" <<EOF
local({
  options(
    repos = c(
      RSPM = "https://packagemanager.posit.co/cran/__linux__/${CODENAME}/latest",
      CRAN = "https://cloud.r-project.org",
      BioC = "https://bioconductor.org/packages/release/bioc"
    ),
    HTTPUserAgent = sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version[["platform"]],
            R.version[["arch"]], R.version[["os"]])
    ),
    Ncpus = max(1L, parallel::detectCores() - 1L),
    download.file.method = "libcurl",
    timeout = 600
  )
  Sys.setenv(
    RENV_CONFIG_PAK_ENABLED = "TRUE",
    RENV_CONFIG_INSTALL_VERBOSE = "FALSE",
    LIBARROW_BINARY = "true",
    NOT_CRAN = "true"
  )
})
EOF

# -------------------------------------------------------------------
# 5. Persistent R library + renv cache
# -------------------------------------------------------------------
mkdir -p "$HOME/.R/library" "$HOME/.cache/R/renv"
echo 'R_LIBS_USER="$HOME/.R/library"' >> "$HOME/.Renviron"
echo "RENV_PATHS_CACHE=$HOME/.cache/R/renv" >> "$HOME/.Renviron"

# -------------------------------------------------------------------
# 6. Pre-install heavy packages (CRAN + Bioconductor + GitHub).
#    Use pak — it resolves system deps in one pass and uses PPM binaries.
# -------------------------------------------------------------------
R -q -e 'install.packages(c("pak","renv"))'

R -q -e 'pak::pak(c(
  "shiny","bslib","shinyWidgets","DT","plotly","ggplot2","patchwork",
  "Seurat","SeuratObject","clustree","harmony",
  "viridisLite","RColorBrewer","qs2","arrow","data.table",
  "fs","glue","future","future.apply","promises","later",
  "shinycssloaders","shinyjs","waiter","ggrastr",
  "shinyFiles","mclust","leidenAlg",
  "rmarkdown","knitr","quarto",
  "testthat","shinytest2","devtools","pkgload",
  "bioc::rhdf5","bioc::zellkonverter",
  "immunogenomics/presto"
))'

# -------------------------------------------------------------------
# 7. Sanity check
# -------------------------------------------------------------------
R -q -e 'cat("R:",as.character(getRversion()),"\n");
         pkgs <- c("Seurat","clustree","harmony","arrow","rhdf5","presto","qs2");
         stopifnot(all(sapply(pkgs, requireNamespace, quietly=TRUE)));
         cat("All required packages load OK\n")'
```

Notes:

- **PPM gives binary packages** — installing Seurat from source would need 30+ minutes plus an unpredictable mix of system libraries. Binaries make the same install ~2 minutes. The `__linux__/<codename>` URL pattern combined with the `HTTPUserAgent` option is what triggers binary delivery; both are required.
- **`pak` over `install.packages`** — `pak` resolves CRAN, Bioconductor, and GitHub sources in a single dependency graph, parallelizes downloads, and (crucially) reports system-library gaps clearly. Falls back gracefully if a package has no PPM binary.
- **`leidenAlg`** is on CRAN and provides Leiden clustering for Seurat 5; without it the app falls back to Louvain (the spec says so).
- **`presto`** is on GitHub (`immunogenomics/presto`), not CRAN. Required for fast Wilcoxon marker tests on 5K-gene data.

## 2. Environment variables

Set these at the env level (so they persist across shells); `~/.Renviron` covers R sessions.

| Variable | Value | Purpose |
|---|---|---|
| `R_LIBS_USER` | `$HOME/.R/library` | Persistent user library across sessions |
| `RENV_PATHS_CACHE` | `$HOME/.cache/R/renv` | Persist renv cache across sessions |
| `RENV_CONFIG_PAK_ENABLED` | `TRUE` | renv uses pak for installs |
| `LIBARROW_BINARY` | `true` | `arrow` installs as binary on Linux |
| `MAKEFLAGS` | `-j$(nproc)` | Parallel compilation when source builds are unavoidable |
| `NOT_CRAN` | `true` | Suppresses CRAN-incompatible test paths in some packages |
| `DEBIAN_FRONTEND` | `noninteractive` | Silent apt during setup |
| `TZ` | `America/New_York` | Match user's location; deterministic dates in tests |

Optional but recommended for the agent's working comfort:

| Variable | Value | Purpose |
|---|---|---|
| `R_PROFILE_USER` | `$HOME/.Rprofile` | Explicit (the script writes to this path) |
| `R_DEFAULT_PACKAGES` | `datasets,utils,grDevices,graphics,stats,methods` | Standard set; explicit for reproducibility |
| `OMP_NUM_THREADS` | `4` | Cap OpenMP threads — Seurat/Matrix occasionally over-subscribe |
| `OPENBLAS_NUM_THREADS` | `4` | Same, for BLAS |

Do **not** set:

- `GITHUB_TOKEN` unless the user explicitly provides one. The packages installed from GitHub here (`presto`) are public; anonymous installs work and rate limits are sufficient for one setup pass per env.
- API keys for any third-party service. The app does not call out.
- HiPerGator credentials. Out of scope.

## 3. Network allowlist

These are the domains the setup script and subsequent build/test cycles need to reach. List comes from the actual install plan; no speculative entries.

**Required for the setup script:**
- `cloud.r-project.org`, `cran.r-project.org` — CRAN apt repo and source fallback
- `packagemanager.posit.co` — PPM binary packages (the bulk of installs)
- `bioconductor.org` — for `rhdf5`, `zellkonverter`
- `*.bioconductor.org` — Bioconductor repository subdomains
- `github.com`, `api.github.com`, `codeload.github.com` — `presto` install, plus any future GitHub-hosted packages
- `objects.githubusercontent.com`, `raw.githubusercontent.com` — release artifacts and source tarballs from GitHub
- `quarto.org` — Quarto release index (only if the script falls back from the pinned download URL)
- `keyserver.ubuntu.com`, `keys.openpgp.org` — apt key servers (sometimes used by `apt-key` helpers)
- `archive.ubuntu.com`, `security.ubuntu.com`, `ports.ubuntu.com` — base Ubuntu apt mirrors
- `ppa.launchpad.net` — only if a PPA is added (not in this script, but listed defensively)

**Required during development and test runs:**
- `github.com`, `*.github.com` — git push/pull, issue links
- `packagemanager.posit.co` — `renv::restore()` resolves new deps if the lockfile changes
- `cloud.r-project.org` — fallback when PPM lacks a package

**Not required and should be excluded:**
- `*.10xgenomics.com` — the panel CSVs are committed to the repo; runtime fetch is wrong.
- `*.spotify.com`, `*.amazonaws.com`, third-party MCPs — none of these are touched by this app.
- The user's HiPerGator endpoints, Tailscale, lab fileservers — out of scope.

If the cloud env uses a wildcard-supporting allowlist, the minimal practical set is:

```
cloud.r-project.org
cran.r-project.org
packagemanager.posit.co
*.bioconductor.org
bioconductor.org
github.com
*.github.com
*.githubusercontent.com
quarto.org
*.ubuntu.com
keyserver.ubuntu.com
keys.openpgp.org
```

## 4. Persistence and cache hygiene

The setup script writes to `$HOME/.R/library` and `$HOME/.cache/R/renv`. If the cloud env exposes a workspace volume that survives session restarts (most do), put both under that volume by changing the script's paths — e.g., `/workspace/.R/library` and `/workspace/.cache/R/renv`. That way subsequent Claude Code sessions reuse the compiled library and avoid the 10-minute reinstall.

Do not commit either directory to git. `.gitignore` should include:

```
/.R/
/.cache/
/cache/
/renv/library/
/renv/staging/
/renv/python/
.Rproj.user/
.Rhistory
.RData
*.qs2
```

`renv.lock` and `renv/activate.R` **are** committed — that is what makes the project reproducible across the cloud env and HiPerGator.

## 5. What the agent should do if setup fails

If `apt-get install` fails because a system library is missing from the allowlist or a mirror, surface the exact error before improvising. If `pak` reports an unresolved system dependency (it usually names the apt package), add the apt package to step 1 of the setup script and re-run only step 6, not the whole script. Do not silently substitute a different R package for one that fails to install — the version pins in `CLAUDE.md` were chosen deliberately.
