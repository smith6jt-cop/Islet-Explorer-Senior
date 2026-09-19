#!/usr/bin/env bash
set -euo pipefail

# Install a local static Avivator under shiny/www/avivator.
# Prefers building from a local source tree if present; otherwise mirrors the public site.
#
# NOTE 2026-09-18: there is no `shiny/` directory in this repository. The app lives at
# `apps/islet/` and its static root is `apps/islet/www/`, so WWW_DIR below points at a path that
# does not exist and the install lands outside the app. Recorded as a finding; the path is not
# changed here.

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WWW_DIR="$APP_DIR/shiny/www"
DEST_DIR="$WWW_DIR/avivator"

mkdir -p "$DEST_DIR"

if [[ -f "$DEST_DIR/index.html" || -f "$DEST_DIR/dist/index.html" || -f "$DEST_DIR/sites/avivator/dist/index.html" ]]; then
  echo "Avivator already present at $DEST_DIR"
  exit 0
fi

# 1) Build from source if a source checkout is present in DEST_DIR
if [[ -f "$DEST_DIR/package.json" && -d "$DEST_DIR/sites/avivator" ]]; then
  echo "Found Avivator source at $DEST_DIR; attempting build (requires Node + pnpm) ..."
  if ! command -v pnpm >/dev/null 2>&1; then
    echo "pnpm is not installed. Install via: corepack enable && corepack prepare pnpm@latest --activate" >&2
    echo "Falling back to mirroring public site..." >&2
  else
    (cd "$DEST_DIR" && pnpm install && pnpm --filter @viv/avivator build)
    if [[ -f "$DEST_DIR/sites/avivator/dist/index.html" ]]; then
      echo "Build complete. You can serve from $DEST_DIR/sites/avivator/dist"
      exit 0
    else
      echo "Build did not produce dist; falling back to mirror." >&2
    fi
  fi
fi

# 2) Mirror the public Avivator site
echo "Mirroring public Avivator to $DEST_DIR (requires wget) ..."
if ! command -v wget >/dev/null 2>&1; then
  echo "wget not found. Install with: sudo apt-get install -y wget" >&2
  exit 1
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
wget \
  --quiet \
  --recursive \
  --no-parent \
  --page-requisites \
  --adjust-extension \
  --convert-links \
  --no-host-directories \
  --directory-prefix "$TMPDIR" \
  https://viv.gehlenborglab.org/avivator/

# Move mirrored content into DEST_DIR
rsync -a "$TMPDIR/avivator/" "$DEST_DIR/"

if [[ -f "$DEST_DIR/index.html" ]]; then
  echo "Avivator installed at $DEST_DIR"
else
  echo "Failed to mirror Avivator. Please check network and try again." >&2
  exit 1
fi

