# Serving both apps locally

`deploy/DEPLOY.md` covers UF RC PubApps. This covers the **local** shiny-server + nginx host, where
two apps run side by side:

| app | what it is | URL | mounted at |
|---|---|---|---|
| **Junior** | the separate `panc_CODEX` app | `/islet-explorer/` | `/srv/shiny-server/islet-explorer` → `/data/panc_CODEX/Islet-Explorer/app/shiny_app` |
| **Senior** | this repo, PhenoCycler run cohort | `/islet-explorer-senior/` | `/srv/shiny-server/islet-explorer-senior` → `apps/islet` |

## Junior is working and nothing here touches it

Checked 2026-09-14: `GET /islet-explorer/` → **200**, and a Range request on one of its OME-TIFFs
→ **206 Partial Content**. The maintainer repointed the nginx alias to `/data/panc_CODEX/...`
after the project moved off `~/panc_CODEX`, and that fix is live.

**Junior's two nginx location blocks in `nginx-islet-explorer.conf` are byte-identical to the
installed file**, and no step below removes or repoints Junior's symlink. An earlier draft of this
document told you to delete that symlink and remount Senior in its place, which would have removed
Junior and served Senior at Junior's URL. That step is gone. If you kept a copy of the old
procedure, discard it — the give-away is any step that deletes `/srv/shiny-server/islet-explorer`
instead of creating `/srv/shiny-server/islet-explorer-senior`.

## The two privileged steps

Everything else is already done and needs no root.

```bash
# 1. Mount Senior as a SECOND app. Additive: Junior's symlink is untouched.
sudo ln -s /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/apps/islet \
           /srv/shiny-server/islet-explorer-senior

# 2. Install the nginx config carrying BOTH apps, then reload.
sudo cp /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/deploy/local/nginx-islet-explorer.conf \
        /etc/nginx/conf.d/islet-explorer.conf
sudo nginx -t && sudo systemctl reload nginx
```

The nginx file **replaces** `/etc/nginx/conf.d/islet-explorer.conf`; it must not be installed
alongside it. Two `server` blocks on `listen 8080; server_name _;` are a duplicate-default
conflict — nginx keeps the first and silently ignores the second's locations.

`systemctl reload nginx` is a graceful reload and does not interrupt Junior. Do **not** restart
shiny-server: adding a directory under `site_dir` is picked up per request, and a restart would
drop Junior's live sessions for no reason.

## Verify — both apps, and Range on both image routes

```bash
curl -s -o /dev/null -w 'junior app     %{http_code}\n' http://127.0.0.1:8080/islet-explorer/
curl -s -o /dev/null -w 'senior app     %{http_code}\n' http://127.0.0.1:8080/islet-explorer-senior/
curl -sI -H 'Range: bytes=0-99' \
  http://127.0.0.1:8080/islet-explorer/local_images/Aab_6450.ome.tiff | head -1
curl -sI -H 'Range: bytes=0-99' \
  http://127.0.0.1:8080/islet-explorer-senior/local_images/115.ome.tiff | head -1
```

Expect `200`, `200`, and **`HTTP/1.1 206 Partial Content`** twice. A `200` on either image line
means Range is not being served and Avivator will render a blank canvas with no error.

## Already done — no root needed

**`apps/islet/.Rprofile`** (gitignored, mode 0644). It does two things, and both matter:

- **`.libPaths()` prepends `/home/smith6jt/R/x86_64-pc-linux-gnu-library/4.5`.** shiny-server runs
  as `shiny`, which has none of the app's packages. Measured under a non-maintainer `HOME`:
  without this, `shiny` itself is unavailable and the app cannot start; `duckdb` is unavailable,
  so the app would fall back to Excel and lose the tissue, cell-distance and neighbourhood views;
  `lmerTest` is unavailable, so the Statistics tab's primary test fails.
- **Sets `ISLET_DOMAIN_CONFIG` only when unset**, so an explicit value still wins. The only change
  is the unset default: it was `islet_default.yml`, it is now `run20_default.yml`. The testthat
  suite sets its domain per test and is unaffected — measured 2026-09-18, `tests/testthat` runs
  737 assertions over 69 `test_that` blocks with one known error, `test-card-registry.R:176`.

This is `.Rprofile` and not `.Renviron` on purpose — Junior uses the same pattern, and both
alternatives were tried and rejected:

1. `R_LIBS_USER=${R_LIBS_USER-/home/...}` in a `.Renviron` is a **guaranteed no-op**. R assigns
   `R_LIBS_USER` a default *before* reading `.Renviron`, so the `${VAR-default}` form always finds
   a value and keeps it. Under `shiny` that default is `/home/shiny/R/...`, which is empty. It
   looks correct when tested as the maintainer only because that user's default already is the
   intended library — a test that cannot fail.
2. An app-local `.Renviron` **replaces** `~/.Renviron` rather than adding to it, which would have
   removed `KEY` and `BASE` from every R session started in `apps/islet`.

Junior's own `.Renviron` is mode `0600` and therefore unreadable by `shiny` — it is inert. **Do not
"fix" it while you are in here**: that would change Junior's running configuration, which is the
one thing this exercise is meant not to do.

**`ISLET_LAYERS_DIR` is NOT set here, and does not need to be on this host.** The code default is
the absolute `/data/islet_layers` (`ieCore::ie_layers_dir()`), which is where the nine layer
directories actually are on this workstation — 430 MB, `drwxrwxr-x smith6jt`, world-readable
through a world-traversable `/data`, so uid `shiny` reaches them. A host where that path does not
exist must set the variable; `deploy/env.sh.example` and `deploy/islet-explorer.service` carry it
for exactly that reason.

**Neither variable travels with a clone.** `.Rprofile` is gitignored (`.gitignore:31`), so this
section describes THIS workstation and nothing else. On any other host the two settings come from
`deploy/env.sh`, which is also gitignored — its template is `deploy/env.sh.example`, and since
2026-09-18 the template names both.

## Known, not fixed here

- **The Viewer lists 16 images, including donors 6457 and 6579.** `core/multiplex5000/multiplex5000/
  cohort.py` centrally excludes 6457; 6579's exclusion was lifted on 2026-09-17, but it has no
  typing yet. `discover_islet_assets()` is a plain
  listing of `www/local_images` with no cohort filter, so the dropdown offers donors that are
  excluded from analysis. The nginx alias deliberately serves the same directory — making it
  stricter than the app would only turn those entries into 404s. The real fix is a domain-aware
  filter in `discover_islet_assets()`, reading the active dataset's donors rather than restating
  the exclusion list. Pre-existing; unchanged by this deployment.
- **DuckDB's spill directory** defaults to `.tmp` relative to the app directory, which `shiny`
  cannot write. Only reached under memory pressure, but it would surface as a query-time failure.

## Two things that take the site down without an obvious cause

- **`chmod 700 /home/smith6jt`.** `shiny` and nginx both reach the app and the images through
  `o+x` on the home directory, which is currently `drwxrwxrwx`. Remove it and both apps stop.
- **Treating `shiny::runApp` as a deployment test.** It answers `200` with the whole body for a
  Range request, so the Viewer appears to work in development and fails only once deployed.
