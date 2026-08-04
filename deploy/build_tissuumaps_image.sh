#!/usr/bin/env bash
# deploy/build_tissuumaps_image.sh — build the TissUUmaps viewer image from OUR fork.
#
# The Viewer's TissUUmaps backend runs smith6jt-cop/TissUUmaps rather than the
# published cavenel/tissuumaps image, so the server we serve is the code we
# control. The .tmap project path the app uses works on either (the fork's
# additions are confined to _fnfilter, the /slide route and
# flask.standalone.addLayer) -- running the fork means fixes we make there ship
# without waiting on upstream.
#
# The fork's container/Dockerfile builds on Alpine and COMPILES openslide 3.4.1
# from source, so expect ~10-20 minutes for a cold build. It is cached after that.
#
# Usage:
#   bash deploy/build_tissuumaps_image.sh                  # build from master
#   TISSUUMAPS_REF=74afd3a bash deploy/build_tissuumaps_image.sh   # pin a commit
#   TISSUUMAPS_SRC=~/src/TissUUmaps bash deploy/build_tissuumaps_image.sh  # reuse a clone
#
# Produces two tags: localhost/islet-tissuumaps:<short-sha> and :latest, which is
# exactly what deploy/islet-tissuumaps.container references. The localhost/ prefix
# is written out rather than left implicit: podman normalizes an unqualified local
# build to localhost/<name>, but docker does not, so tagging it explicitly keeps
# the image name identical under either engine. Pin the unit to the sha tag if you
# want the running service frozen to a known build.
set -euo pipefail

REPO_URL="${TISSUUMAPS_REPO:-https://github.com/smith6jt-cop/TissUUmaps.git}"
REF="${TISSUUMAPS_REF:-master}"
SRC="${TISSUUMAPS_SRC:-}"
IMAGE="${TISSUUMAPS_IMAGE:-islet-tissuumaps}"
# Fully-qualified local name — must match Image= in islet-tissuumaps.container.
IMAGE_REF="localhost/${IMAGE}"

# RC PubApps is rootless podman, no docker. Prefer podman, fall back to docker.
ENGINE="${CONTAINER_ENGINE:-}"
if [[ -z "$ENGINE" ]]; then
  if command -v podman >/dev/null 2>&1; then ENGINE=podman
  elif command -v docker >/dev/null 2>&1; then ENGINE=docker
  else echo "ERROR: neither podman nor docker found" >&2; exit 1
  fi
fi
echo "[build] engine: $ENGINE"

CLEANUP=""
if [[ -z "$SRC" ]]; then
  SRC="$(mktemp -d)/TissUUmaps"
  CLEANUP="$SRC"
  echo "[build] cloning $REPO_URL @ $REF"
  git clone --quiet "$REPO_URL" "$SRC"
  git -C "$SRC" checkout --quiet "$REF"
else
  echo "[build] using existing clone: $SRC"
fi

SHA="$(git -C "$SRC" rev-parse --short HEAD)"
echo "[build] fork commit: $SHA ($(git -C "$SRC" log -1 --format=%s | cut -c1-60))"

# Build context is the repo root -- the Dockerfile COPYs ./container/ and ./tissuumaps/.
"$ENGINE" build \
  -t "${IMAGE_REF}:${SHA}" \
  -t "${IMAGE_REF}:latest" \
  -f "$SRC/container/Dockerfile" \
  "$SRC"

echo "[build] built ${IMAGE_REF}:${SHA} and ${IMAGE_REF}:latest"
"$ENGINE" images "$IMAGE_REF" || true

if [[ -n "$CLEANUP" ]]; then
  rm -rf "$(dirname "$CLEANUP")"
fi

cat <<EOF

Next:
  1. systemctl --user restart islet-tissuumaps      # if already installed
  2. curl -s -o /dev/null -w '%{http_code}\\n' http://localhost:<TM_PORT>/ping   # -> 200

To freeze the service to this build, set Image=${IMAGE_REF}:${SHA}
in ~/.config/containers/systemd/islet-tissuumaps.container.
EOF
