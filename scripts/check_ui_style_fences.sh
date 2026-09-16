#!/usr/bin/env bash
# Two fences for the Modernist design in apps/islet. Run from the repo root.
#
#   ./scripts/check_ui_style_fences.sh
#
# Ported from Junior (f827106) with Senior's own paths, Senior's own baseline and
# Senior's own re-measured accent collisions. Junior's numbers are NOT carried over:
# the two apps do not share a palette, and it turned out the fence matters more here.
#
# ---------------------------------------------------------------------------
# Fence 1 -- no inline skin styles.
#
#   An inline `style=` attribute beats any stylesheet without `!important`, so every
#   one that sets colour / background / border-radius / shadow is a hole in the design
#   system: modernist-shiny.css cannot reach it and the element keeps its 2019 look.
#   Layout declarations (display, flex, gap, height, overflow) are fine and are
#   deliberately NOT matched -- layout is the module's business, skin is the sheet's.
#
#   RATCHET, NOT A GATE. The skin reconcile of 2026-09-14 took it from 216 lines
#   across 14 files to the CURRENT BASELINE:
#
#        35 matching lines across 7 files
#         13 R/mod_statistics_server.R    9 R/mod_trajectory_server.R
#          5 R/mod_plot_server.R          3 R/mod_statistics_ui.R
#          2 R/mod_inflammation_server.R  2 R/00_globals.R
#          1 app.R
#
#   Junior's post-redesign figure is 36 across 7, and the residue is the same kind
#   in both: legend swatches, statistical-state callouts, and colours inside
#   generated HTML strings.
#
#   WHAT THE RESIDUE IS, AND WHY IT IS NOT SKIN. The rule the reconcile applied is
#   "strip chrome, keep data and state", and it is worth stating because the count
#   will never reach zero:
#     * LEGEND SWATCHES -- `div(style = "...; background-color: #0066CC; ...")` and
#       the sprintf'd ones. These ARE the data encoding, drawn as HTML instead of
#       inside the figure, and they must track the palette rather than the sheet.
#     * STATISTICAL-STATE CALLOUTS in mod_statistics_server.R -- the amber, red and
#       green result grounds. They encode what the test found, not how the page
#       looks.
#     * `00_globals.R` -- two colours inside an HTML STRING, not a tag attribute.
#     * `app.R:637` -- the phenotype-rule swatch, `background:%s` where `%s` is a
#       phenotype colour. It is left in the count rather than pattern-matched away,
#       because a regex that forgives "background:%s" would forgive every sprintf'd
#       skin style too. Expect the count to stop at 1 in app.R, not 0.
#
#   BLIND SPOT, measured and not fixed. The pattern is single-line: it matches
#   `style = "..."` but not `style = paste0("...",\n "...")`, where the declarations
#   sit on continuation lines. Four blue-gradient banners and three inline copies of
#   section_heading hid there through the first pass of this port -- the fence read
#   clean on files that were visibly blue. A multi-line scan is the check that
#   actually found them; grep for the palette's own hexes, not for `style =`.
#
# ---------------------------------------------------------------------------
# Fence 2 -- the chrome accent must never appear inside plot code.
#
#   Measured for Senior (CIEDE2000, sRGB/D65) over the 120 distinct colours its data
#   palettes actually declare -- the two domain configs, the four phenotype palettes
#   and the marker palette in the helpers, and PAL_LEIDEN:
#
#     IE_CRIMSON #b71533   deltaE  1.2 from #BE0032  Beta cell, Max Distinction
#                          deltaE  2.4 from #B2182B  Aab+, "Diverging" donor palette
#                          deltaE  7.8 from #BC4749  Blood vessel, run20 phenotypes
#                          deltaE  8.9 from #AA4455  CD8a T cell, Colourblind Safe
#                          deltaE  9.2 from #D32F2F  Beta cell, High Contrast
#                          -- FIVE data colours inside deltaE 10.
#     IE_GOLD    #c8a02e   deltaE  9.3 from #CCBB44  PAL_LEIDEN slot 4, in the very
#                          file that defines the accent (R/theme_modernist.R:62).
#
#   deltaE < 10 is not reliably discriminable side by side, let alone across a legend
#   and a plot. A crimson rule beside a Beta-cell point reads as a Beta-cell key.
#   So the accent stays lines, borders and underlines; a filled swatch is a legend.
#
#   IE_ACCENT_700 #800c21 is the deep step the modernist rule prescribes for accent
#   text at body size, and it is NOT fenced -- but it is not clean either: deltaE 8.0
#   from #7F180D (Macrophage, Max Distinction). Do not use it as a fill either.
#
#   Scope is wider than Junior's: Senior keeps plot code in spatial_helpers.R,
#   viewer_helpers.R and rdeck_helpers.R as well, and CELL_DIST_MARKER_COLORS lives
#   in spatial_helpers.R. R/theme_modernist.R is excluded because defining the two
#   constants is its job.
#
# Exit 0 = both clean. Exit 1 = a fence tripped. Exit 2 = could not run.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 2

FAIL=0
APP="apps/islet"
[ -d "$APP/R" ] || { echo "no $APP/R -- run from the repo root" >&2; exit 2; }

echo "== Fence 1: inline skin styles (ratchet; baseline 35 lines / 7 files) =="
HITS=$(grep -nE "style *=.*(background|linear-gradient|border-radius|box-shadow|(^|[^-])color: *#)" \
        "$APP"/R/*.R "$APP"/app.R 2>/dev/null)
if [ -n "$HITS" ]; then
  echo "$HITS" | awk -F: '{print $1}' | sort | uniq -c | sort -rn | sed 's/^/  /'
  N=$(echo "$HITS" | wc -l)
  echo "  TOTAL: $N inline skin declarations across $(echo "$HITS" | awk -F: '{print $1}' | sort -u | wc -l) files"
  if [ "$N" -gt 35 ]; then
    echo "  REGRESSION: $((N - 35)) above the 35-line baseline."
  elif [ "$N" -lt 35 ]; then
    echo "  Below baseline by $((35 - N)) -- update the baseline in this script's header."
  fi
  FAIL=1
else
  echo "  clean"
fi

echo
echo "== Fence 2: chrome accent inside plot code =="
HITS2=$(grep -nEi 'b71533|c8a02e' \
         "$APP"/R/mod_*_server.R "$APP"/R/*_helpers.R 2>/dev/null)
if [ -n "$HITS2" ]; then
  echo "$HITS2" | sed 's/^/  /'
  echo "  Chrome accent leaked into plot code. Crimson is deltaE 1.2 from the Max"
  echo "  Distinction Beta-cell colour and 2.4 from the Diverging Aab+ colour; gold is"
  echo "  9.3 from PAL_LEIDEN slot 4. Use a data palette, or IE_RULE / IE_MUTED for ink."
  FAIL=1
else
  echo "  clean"
fi

exit $FAIL
