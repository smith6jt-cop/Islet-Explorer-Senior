/**
 * Probe each pyramid level to find which resolution level(s) contain corrupt tiles.
 * Run on the problem image in QuPath script editor.
 *
 * This reads a sampling of tiles at each level to isolate where the ZipException occurs.
 * If level 0 is clean, the rewrite script will fix the issue.
 * If level 0 is corrupt, you'll need a Python tifffile-based repair.
 */

import qupath.lib.images.servers.ImageServer
import qupath.lib.regions.RegionRequest

def server = getCurrentServer()
def metadata = server.getMetadata()
def downsamples = server.getPreferredDownsamples()

println "=== Image Info ==="
println "Size: ${metadata.getWidth()} x ${metadata.getHeight()}"
println "Channels: ${metadata.getSizeC()}"
println "Pyramid levels: ${downsamples.length} (downsamples: ${downsamples})"
println "Tile size: ${metadata.getPreferredTileWidth()} x ${metadata.getPreferredTileHeight()}"
println ""

int tileW = metadata.getPreferredTileWidth() ?: 512
int tileH = metadata.getPreferredTileHeight() ?: 512

for (int level = 0; level < downsamples.length; level++) {
    double ds = downsamples[level]
    int levelW = (int)(metadata.getWidth() / ds)
    int levelH = (int)(metadata.getHeight() / ds)
    int nTilesX = (int)Math.ceil(levelW / (double)tileW)
    int nTilesY = (int)Math.ceil(levelH / (double)tileH)
    int totalTiles = nTilesX * nTilesY

    println "--- Level ${level} (downsample=${ds}, ${levelW}x${levelH}, ${totalTiles} tiles) ---"

    int tested = 0
    int failed = 0
    List<String> failedTiles = []

    // Test all tiles at this level
    for (int ty = 0; ty < nTilesY; ty++) {
        for (int tx = 0; tx < nTilesX; tx++) {
            int x = (int)(tx * tileW * ds)
            int y = (int)(ty * tileH * ds)
            int w = (int)(tileW * ds)
            int h = (int)(tileH * ds)

            // Clamp to image bounds
            if (x + w > metadata.getWidth()) w = metadata.getWidth() - x
            if (y + h > metadata.getHeight()) h = metadata.getHeight() - y

            try {
                def request = RegionRequest.createInstance(server.getPath(), ds, x, y, w, h)
                server.readRegion(request)
                tested++
            } catch (Exception e) {
                failed++
                tested++
                failedTiles << "  tile(${tx},${ty}) at pixel(${x},${y}): ${e.getCause()?.getMessage() ?: e.getMessage()}"
            }
        }

        // Progress update per row
        if ((ty + 1) % 10 == 0 || ty == nTilesY - 1) {
            print "  ... tested row ${ty + 1}/${nTilesY} (${tested} tiles, ${failed} failures)"
        }
    }

    if (failed == 0) {
        println "  PASSED: All ${tested} tiles OK"
    } else {
        println "  FAILED: ${failed}/${tested} tiles corrupt:"
        failedTiles.each { println it }
    }
    println ""
}

println "=== Done ==="
println "If only non-zero levels failed: run rewrite_ome_tiff.groovy to fix."
println "If level 0 failed: the base resolution has corrupt tiles — use Python tifffile to repair."