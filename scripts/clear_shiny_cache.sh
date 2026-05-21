#!/bin/bash
# Script to clear shiny-server cache and restart services
# Usage: ./clear_shiny_cache.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=== Clearing Shiny Server Cache ==="

# 1. Clear any temporary R files and caches
echo "Clearing R temp files..."
rm -rf /tmp/Rtmp*
rm -rf /var/lib/shiny-server/bookmarks/*
rm -rf ~/.Rcache

# 2. Clear nginx cache if it exists
echo "Clearing nginx cache..."
if [ -d /var/cache/nginx ]; then
    sudo rm -rf /var/cache/nginx/*
fi

# 3. Update modification time on app.R to force reload
echo "Touching app.R to force reload..."
touch "$REPO_ROOT/app/shiny_app/app.R"

# 4. If running as root, restart shiny-server
if [ "$EUID" -eq 0 ]; then
    echo "Restarting shiny-server..."
    systemctl restart shiny-server
    systemctl restart nginx
else
    echo "To restart services, run as root:"
    echo "sudo systemctl restart shiny-server"
    echo "sudo systemctl restart nginx"
fi

# 5. Clear any APP_LOADED_FLAG files
echo "Clearing app load flags..."
find /srv/shiny-server -name "APP_LOADED_FLAG" -delete 2>/dev/null || true
find "$REPO_ROOT" -name "APP_LOADED_FLAG" -delete 2>/dev/null || true

echo "=== Cache clearing complete ==="
echo "App version: $(grep 'APP_VERSION' "$REPO_ROOT/app/shiny_app/app.R")"