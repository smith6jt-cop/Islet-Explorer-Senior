# Deploying Islet Explorer to UF RC PubApps

Public, HTTPS hosting of the Shiny app on UF Research Computing's **PubApps** web service
(`https://islet-explorer.rc.ufl.edu`). This runbook is the executable companion to the planning
notes; RC's own docs are the authority for anything RC-side:
[overview](https://docs.rc.ufl.edu/services/web_hosting/) ·
[deployment](https://docs.rc.ufl.edu/services/web_hosting/deployment/) ·
[Shiny](https://docs.rc.ufl.edu/services/web_hosting/shiny/) ·
[systemd](https://docs.rc.ufl.edu/services/web_hosting/systemd/) ·
[price list](https://it.ufl.edu/rc/get-started/price-list/).

## Model in one paragraph

RC gives you a RHEL 8 VM (base **2 vCPU / 16 GB / 1 TB `/pubapps/$project`**, rootless **Podman** +
conda, **no Docker**, **no access to /blue or /orange**). You build a conda env, run the app on an
**RC-assigned port** via `deploy/run.R` (`shiny::runApp(host="0.0.0.0", port=…)`), and keep it alive
with **systemd `--user`** + `loginctl enable-linger`. RC's load balancer reverse-proxies the free
`.rc.ufl.edu` subdomain (+ SSL) to your port. **No shiny-server, no nginx to manage for the app** —
but the Viewer's OME-TIFFs need HTTP Range serving, which `runApp` can't do, so a small **nginx
Podman sidecar** serves them.

## Scope decisions baked into these files

- **AI Assistant: OFF** (`ISLET_ENABLE_AI=FALSE`, the app default) → no UF Navigator API key on the
  public host. The app is safe to leave fully open.
- **Viewer: ON**, images hosted on PubApps → requires purchased storage + the nginx sidecar.
- **Shared code: monorepo-local `apps/R/ieCore`** → deploy the full clone, not an isolated
  copy of `apps/islet`; `r-yaml` loads the validated islet DOMAIN at startup.
- **Which dataset: an explicit choice, made in `deploy/env.sh`** (`ISLET_DOMAIN_CONFIG`) and
  repeated in the unit. It is not a default and must not be left to one. Until 2026-09-18 neither
  file named it, so a clean checkout got the code fallback — `config/islet_default.yml`, the legacy
  CODEX domain — because the file that selects the run cohort on the maintainer's workstation,
  `apps/islet/.Rprofile`, is gitignored and does not travel with a clone. Phase 2's staging and
  this setting are one decision: stage the domain you are going to name.
- **Layer artifacts: an explicit path** (`ISLET_LAYERS_DIR`). The code default is the absolute
  `/data/islet_layers`, which exists on the maintainer's workstation and on no PubApps VM.

## What's in this directory

| File | Purpose |
|---|---|
| `environment.yml` | conda/mamba runtime env (source of truth for deps) |
| `run.R` | app launcher (systemd `ExecStart` runs this) |
| `env.sh.example` | → copy to `env.sh`, systemd `EnvironmentFile` (no secrets) |
| `islet-explorer.service` | systemd `--user` unit for the app |
| `islet-images.container` | Podman Quadlet unit — nginx sidecar for OME-TIFFs |
| `nginx-images.conf` | Range-capable static config for the sidecar |
| `install_runtime_R.R` | OPTIONAL — GitHub pkgs (vitessceR/rdeck), only if needed |
| `build_tissuumaps_image.sh` | OPTIONAL — builds the viewer image from our TissUUmaps fork |
| `islet-tissuumaps.container` | OPTIONAL — Quadlet unit for the TissUUmaps viewer backend |
| `tissuumaps.cfg` | OPTIONAL — read-only Flask config for that container |

Placeholders to replace everywhere: `<project>` (RC service username / instance),
`<APP_PORT>` and `<IMG_PORT>` (RC-assigned).

---

## Phase 0 — Allocation + support ticket

1. Get a **PA-Instance (PubApps)** allocation (free 3-month Trial to pilot, or **$300/yr**). Add
   **PA-STORAGE** for the OME-TIFFs — **+2 TB ($280/yr)** recommended (3 TB total, fits all 22 donors)
   or +1 TB ($140/yr) for the current 15 — and, for concurrency headroom, **+2 PA-CPU ($88/yr)**.
   Recommended "runs well" ≈ **$668/yr**.
2. Open an RC support ticket specifying: **public/open** access; suggested domain
   **`islet-explorer.rc.ufl.edu`**; **two same-origin proxy routes** — `/` → `<APP_PORT>` (with
   **WebSocket upgrade** — Shiny needs it) and `/local_images/` → `<IMG_PORT>`; the storage bump;
   developer GatorLink usernames. RC returns the project username, both ports, and DNS/proxy setup.

## Phase 1 — Provision the env (on the VM)

```bash
# install miniforge/mamba per RC's Shiny doc, then:
mamba env create -p /pubapps/<project>/prod/conda -f deploy/environment.yml
/pubapps/<project>/prod/conda/bin/python -c "import anndata; print(anndata.__version__)"   # bridge check
```

## Phase 2 — Clone + stage data

```bash
mkdir -p /pubapps/<project>/prod/logs
git clone <repo-url> /pubapps/<project>/prod/islet-explorer
cd /pubapps/<project>/prod/islet-explorer
R CMD INSTALL apps/R/ieCore   # package/install contract check; the app sources this same tree
```
Stage data (nothing is in git; PubApps can't see /blue, so copy from the source machine, e.g. via a
HiPerGator hop or Globus):

**Stage the data root of the domain you named in `deploy/env.sh`, and only that one.** The two
domains are separate datasets sharing five donors and no vocabulary; staging both invites serving
one under the other's labels.

- **Run cohort (`config/run20_default.yml`, the value in `env.sh.example`) — 4.6 GB, measured
  2026-09-18** → `/pubapps/<project>/prod/islet-explorer/data/development/typing/app/`.
  `../../data/development/typing/app/…` resolves from `apps/islet/`, which is exactly 2 levels
  below the repo root. It holds `parquet/` (the DuckDB views, 18,384 islets over 20 donors),
  `cells/`, `donors/`, `gson/`, `json/`, the trajectory `adata_ins_root.h5ad` and the workbooks.
  The `*.superseded-*.h5ad` files are history and need not travel.
- **Legacy CODEX cohort (`config/islet_default.yml`) — ~1.9 GB** →
  `/pubapps/<project>/prod/islet-explorer/data/`: copy `data/app_data/` as a plain directory of
  real files (it was symlink-based historically; since 2026-08-16 it is a self-contained copy, so
  there is no symlink layout to preserve). `../../data/app_data/…` resolves the same way.

  > **`data/app_data/cells/` and top-level `data/cells/` are different things.** The former is
  > 5,214 per-islet CSVs the drill-down reads; the latter is 13 GB of hive-partitioned
  > `donor_id=…` pipeline output. Same basename, unrelated content — never merge them. This is
  > why the app keeps the `app_data/` segment rather than reading `data/` flat.
- **Layer artifacts (430 MB, measured 2026-09-18)** → whatever you set `ISLET_LAYERS_DIR` to,
  e.g. `/pubapps/<project>/islet_layers/`. Nine directories: `availability`, `compartments`,
  `distance_profiles`, `grid_composition`, `islet_summary` (143 MB), `pln` (240 MB), `restain`
  (47 MB), `sections`, `trajectory`. These are stamped with the domain they were built from, and
  five of the app's nine cards compare the artifact's recorded `source_domain_root` against the
  active domain root and refuse, by its own named reason, rather than draw (measured 2026-09-18;
  the exact five and their reasons are listed in `deploy/env.sh.example`). Every artifact carrying
  that field on 2026-09-18 records the run-cohort root, so they pair with `run20_default.yml`.
- **OME-TIFFs (~1.2 TB)** → `/pubapps/<project>/islet_ome_tiff/`: Globus-transfer every
  `<case>.ome.tiff` **and** its `<case>.offsets.json` sidecar.

Whichever domain you staged, do **not** copy the ~113 GB of pipeline scratch (`redsea_scratch/`,
`cells_redsea/`, top-level `cells/`, `phenotype/`, `restore*`, `singlecell_protein.h5ad`, …) — the
app never reads it.

## Phase 3 — Configure + start services

```bash
cd /pubapps/<project>/prod/islet-explorer
cp deploy/env.sh.example deploy/env.sh && chmod 600 deploy/env.sh
# edit deploy/env.sh: set <project> and <APP_PORT>, and CONFIRM the two dataset settings --
#   ISLET_DOMAIN_CONFIG  must name the domain whose data root you staged in Phase 2
#   ISLET_LAYERS_DIR     must be where you staged the nine layer directories
# Neither errors when wrong; see the comments in env.sh.example for what each one costs.

mkdir -p ~/.config/systemd/user ~/.config/containers/systemd
cp deploy/islet-explorer.service ~/.config/systemd/user/
cp deploy/islet-images.container ~/.config/containers/systemd/
# edit both copies: set <project>, <APP_PORT>, <IMG_PORT>

loginctl enable-linger "$USER"
systemctl --user daemon-reload
systemctl --user enable --now islet-explorer islet-images
systemctl --user status islet-explorer islet-images
```

## Phase 4 — Verify

```bash
# app up:
curl -s http://localhost:<APP_PORT>/ | head
# image server honors Range (must be 206 + Accept-Ranges, NOT 200):
curl -I -H "Range: bytes=0-1" http://localhost:<IMG_PORT>/local_images/115.ome.tiff
tail -f /pubapps/<project>/prod/logs/app.log      # watch R startup

# THE DATASET IS THE ONE YOU MEANT. 00_00_iecore.R prints one line at startup:
#   [DOMAIN] Loaded <app.name> from <resolved config path>
# Read it from the log. A blank ISLET_DOMAIN_CONFIG does not warn -- it falls back.
grep '^\[DOMAIN\]' /pubapps/<project>/prod/logs/app.log

# The layer root is the one you staged, and is not empty:
systemctl --user show islet-explorer -p Environment | tr ' ' '\n' | grep ISLET_
ls "$ISLET_LAYERS_DIR"      # expect the nine directories from Phase 2, none empty
```
Then browse `https://islet-explorer.rc.ufl.edu` and exercise each tab:
**Plot** (DuckDB load; confirm no AI panel) · **Trajectory** (the reticulate→`anndata::read_h5ad`
check) · **Statistics** (`lmerTest`/`emmeans`) · **Spatial** (DuckDB `tissue`/`cell_distance`) ·
**Drill-down** (per-islet CSV + GeoJSON) · **Viewer** (Avivator streams a donor's OME-TIFF via Range).
Change any input and confirm a plot updates → verifies **WebSockets** pass the proxy (if the page is
inert, WS upgrade isn't enabled — go back to RC).

## Viewer backend — optional TissUUmaps swap

The Viewer defaults to **Avivator** (OME-TIFFs streamed by the nginx sidecar, HTTP Range). To
swap in **TissUUmaps** — region drawing, GeoJSON import/export, per-cell marker overlays, and
no Range requirement — do this instead of (not alongside) the `islet-images` sidecar.

The TissUUmaps we run is **our fork, [smith6jt-cop/TissUUmaps](https://github.com/smith6jt-cop/TissUUmaps)**,
built locally — not the PyPI release or the published `cavenel/tissuumaps` image.

1. **Build a project per donor** on a machine that has the slides + `pyvips`
   (`pip install "pyvips[binary]"`, or system `libvips` + `pip install pyvips`):
   ```bash
   for c in 115 6476 6539 …; do
     python scripts/senior/build_tissuumaps_project.py \
       --image /pubapps/<project>/islet_ome_tiff/$c.ome.tiff \
       --out-dir /pubapps/<project>/tissuumaps \
       --regions data/Islets
   done
   ```
   `--regions` points at the directory of QuPath exports and picks each donor's file by case
   id (`.gz` included); the scale into full-res pixels is derived from the two pixel sizes and
   then checked against the image bounds, so watch for a `WARNING` on the `regions:` line.
   Idempotent — a killed run resumes, and channels extracted by older code or different
   settings are detected and redone automatically (no `--force` needed). Check the printed
   `channels:` line for the real disk ratio before doing all 22, and verify one donor visually
   (`--keep`) first: the layers are **8-bit** proxies of the 16-bit data.
2. **Ask RC for a third proxy route**: `/tissuumaps/` → `<TM_PORT>` (same origin as the app).
3. **Build the image from the fork** — there is no registry to pull it from, so this must
   succeed before the unit will start:
   ```bash
   bash deploy/build_tissuumaps_image.sh          # ~10-20 min cold: compiles openslide 3.4.1
   # TISSUUMAPS_REF=<sha> pins a commit; the script also tags islet-tissuumaps:<short-sha>
   ```
4. **Start the service:**
   ```bash
   cp deploy/islet-tissuumaps.container ~/.config/containers/systemd/   # edit <project>/<TM_PORT>
   systemctl --user daemon-reload && systemctl --user enable --now islet-tissuumaps
   curl -s -o /dev/null -w '%{http_code}\n' http://localhost:<TM_PORT>/ping     # -> 200
   ```
5. **Flip the app** in `deploy/env.sh`: `ISLET_VIEWER_BACKEND=tissuumaps`,
   `TISSUUMAPS_URL=/tissuumaps`, `TISSUUMAPS_PROJECT_DIR=/pubapps/<project>/tissuumaps`, then
   `systemctl --user restart islet-explorer`.

`deploy/tissuumaps.cfg` sets `READ_ONLY = True` and the volume is mounted `:ro` — required,
because the Flask app otherwise exposes a **writable file browser** over its slide directory.
Users can still draw regions in the browser and download the GeoJSON. Rationale and the
verified findings: [`docs/tissuumaps_evaluation.md`](../docs/tissuumaps_evaluation.md).

**Updating the viewer** is manual by design: a locally-built image never auto-updates. To pick
up new fork commits, re-run `deploy/build_tissuumaps_image.sh` and
`systemctl --user restart islet-tissuumaps`.

## Operating notes

- **Update the app:** `git -C /pubapps/<project>/prod/islet-explorer pull && systemctl --user restart islet-explorer`.
- **Reboot:** RC does not auto-start pubapps services on VM reboot — after one, `systemctl --user start islet-explorer islet-images`.
- **Logs:** app → `/pubapps/<project>/prod/logs/app.log`; sidecar → `podman logs islet-images`.
- **Local dev keeps the AI tab:** run with `ISLET_ENABLE_AI=TRUE`.
- **Security:** the previously-committed-on-disk UF Navigator key (`apps/islet/.Renviron`) should be
  **rotated/revoked** — it isn't used on the public host (AI off), but the old key is compromised.
- **Data sensitivity:** the app is fully open; if any data is unpublished/sensitive, ask RC to gate the
  `/` route behind **GatorLink SSO** (a ticket toggle).

---

## Local workstation: two apps intended to be served side by side (2026-08-16)

> **Status 2026-09-02.** The publication build moved to `/data/panc_CODEX/Islet-Explorer`, so the
> `islet-explorer` symlink below and the nginx `local_images` alias both dangle, and the
> `islet-explorer-merged` entry was never created: shiny-server currently serves no Islet Explorer.
> Repointing is a root action and a maintainer decision; nothing in this repo does it.

This workstation is intended to serve **two** builds through shiny-server, deliberately, and they
must not be conflated.

| Path | Points at | Purpose | State |
|---|---|---|---|
| `/srv/shiny-server/islet-explorer` | `~/panc_CODEX/Islet-Explorer/app/shiny_app` | **Publication-linked build. Do not repoint or retire it.** | dangling after the source moved to `/data/panc_CODEX/Islet-Explorer` |
| `/srv/shiny-server/islet-explorer-merged` | `~/IO60panc2nd/Islet-Explorer-Senior/apps/islet` | The merged monorepo app (this repo) | **not created yet** |

**Checked 2026-08-20:** `/srv/shiny-server/` contained only `index.html`, `islet-explorer` and
`sample-apps`. The second entry had not been made. The 2026-09-02 status above supersedes the
earlier claim that the first entry was live.

The two have **separate data roots** and neither reads the other's: the publication build reads
`~/panc_CODEX/Islet-Explorer/data/`, the merged build reads
`~/IO60panc2nd/Islet-Explorer-Senior/data/app_data/`. They are copies, not symlinks, so work on
one cannot perturb the other.

Creating the second entry needs root (the directory is root-owned), so it is **not** done by any
script in this repo:

```bash
sudo ln -sfn /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/apps/islet \
             /srv/shiny-server/islet-explorer-merged
sudo systemctl restart shiny-server
```

Then set `IMAGE_BASE_URL` in `apps/islet/.Renviron` for the Viewer — that file is gitignored and
does not travel with a clone. Optional: `SHOW_INFLAMMATION_TAB=1` to surface the Inflammation tab,
`ISLET_ENABLE_AI=TRUE` for the AI panel (keep it off on anything public — it needs a UF Navigator
key).

Verified 2026-08-16, against the CODEX domain and the state of the app on that date: the merged
app boots against the staged data (DuckDB views registered, mirai pool up, 10,736 islet spatial
records, 59 channels, 15 viewer images), serves HTTP 200, and renders all five tabs it had then.
That tab count is no longer current — `app.R` now registers nine `nav_panel_hidden` panels (Plot,
Trajectory, Statistics, Spatial, Neighborhoods, CrossModal, LymphNode, Viewer, plus Inflammation
behind `SHOW_INFLAMMATION_TAB=1`) and `tests/tab_smoke.R` drives eight of them. Nothing above has
been re-verified against the run cohort; treat the 2026-08-16 line as history, not as a current
check. Run it directly with:

```bash
cd apps/islet && Rscript -e 'shiny::runApp(".", port = 7788)'
```
