/**
 * update_channel_names_from_data_cn.groovy
 *
 * Reads channel names from server.json files inside data_cn (project root).
 * Matches entries to the current project by metadata.name (primary) or
 * URI filename (fallback), then applies channel names and colors.
 *
 * Expected server.json structure:
 * {
 *   "uri": "file:///path/to/image.ome.tiff",   <- optional
 *   "metadata": {
 *     "name": "2007_CC1A",                      <- matched against project image names
 *     "channels": [ {"name": "CD107a", "color": -65536}, ... ]
 *   }
 * }
 *
 * QuPath version: 0.5 / 0.6 compatible
 */

import qupath.lib.images.servers.ImageChannel
import qupath.lib.images.servers.ImageServerMetadata
import com.google.gson.JsonParser

// ── Configuration ─────────────────────────────────────────────────────────────

final String SOURCE_FOLDER   = "data_cn"
final String SOURCE_JSONNAME = "server.json"
final boolean SAVE_IMAGE_DATA = true

// ── End configuration ──────────────────────────────────────────────────────────

def project    = getProject()
def projectDir = project.getPath().getParent().toFile()
def sourceDir  = new File(projectDir, SOURCE_FOLDER)

if (!sourceDir.isDirectory()) {
    println "ERROR: Source folder not found: ${sourceDir.absolutePath}"
    return
}

// ── Step 1: Index server.json files from data_cn ──────────────────────────────
// Each entry is indexed under ALL candidate keys it can be matched by:
//   1. metadata.name  (e.g. "2007_CC1A")
//   2. URI filename   (e.g. "20_07_sp_cc1a.ome.tiff")
// Keys are stored lower-case for case-insensitive matching.

def sourceMap = [:]   // key (lower-case) -> list of [name, color] maps

sourceDir.listFiles()?.sort()?.each { entryDir ->
    if (!entryDir.isDirectory()) return

    def jsonFile = new File(entryDir, SOURCE_JSONNAME)
    if (!jsonFile.exists()) return

    try {
        def root = JsonParser.parseString(jsonFile.text).getAsJsonObject()

        // Must have metadata.channels
        if (!root.has("metadata") || !root.get("metadata").isJsonObject()) {
            println "[WARN] No 'metadata' object in ${jsonFile.absolutePath} -- skipping"
            return
        }
        def metadata = root.getAsJsonObject("metadata")
        if (!metadata.has("channels") || !metadata.get("channels").isJsonArray()) {
            println "[WARN] No 'channels' array in ${jsonFile.absolutePath} -- skipping"
            return
        }

        def channelList = []
        metadata.getAsJsonArray("channels").each { el ->
            def ch = el.getAsJsonObject()
            channelList << [
                name : ch.has("name")  ? ch.get("name").getAsString() : null,
                color: ch.has("color") ? ch.get("color").getAsInt()   : null
            ]
        }

        // Collect all candidate keys for this entry
        def keys = []

        // Key 1: metadata.name
        if (metadata.has("name") && metadata.get("name").isJsonPrimitive()) {
            keys << metadata.get("name").getAsString().toLowerCase()
        }

        // Key 2: filename from URI
        if (root.has("uri") && root.get("uri").isJsonPrimitive()) {
            String uri = root.get("uri").getAsString()
            try {
                keys << new File(new URI(uri).getPath()).getName().toLowerCase()
            } catch (Exception ignored) {
                keys << new File(uri).getName().toLowerCase()
            }
        }

        if (keys.isEmpty()) {
            println "[WARN] Could not determine any match key from ${jsonFile.absolutePath} -- skipping"
            return
        }

        keys.each { k -> sourceMap[k] = channelList }
        println "[INDEX] keys=${keys}  ->  ${channelList.size()} channels  (${entryDir.name})"

    } catch (Exception e) {
        println "[WARN] Failed to parse ${jsonFile.absolutePath}: ${e.message}"
    }
}

println "\nIndexed ${sourceMap.size()} keys from '${SOURCE_FOLDER}'\n"

if (sourceMap.isEmpty()) {
    println "ERROR: No valid ${SOURCE_JSONNAME} files found under ${sourceDir.absolutePath}"
    return
}

// ── Step 2: Apply channel names to current project entries ────────────────────

int updated = 0
int skipped = 0
int failed  = 0

project.getImageList().each { entry ->
    def rawName   = entry.getImageName()
    def imageName = rawName.toLowerCase()

    // Exact match
    def channelList = sourceMap[imageName]

    // Extension-stripped fallback (e.g. "2007_cc1a.ome.tiff" -> "2007_cc1a")
    if (channelList == null) {
        def noExt = imageName.replaceAll("\\.[^.]+\$", "")
        channelList = sourceMap[noExt]
    }

    // Partial match: project name is a substring of a source key or vice versa
    if (channelList == null) {
        def noExt = imageName.replaceAll("\\.[^.]+\$", "")
        channelList = sourceMap.find { k, v ->
            k.replaceAll("\\.[^.]+\$", "").contains(noExt) ||
            noExt.contains(k.replaceAll("\\.[^.]+\$", ""))
        }?.value
    }

    if (channelList == null) {
        println "[SKIP] ${rawName} -- no matching server.json found in ${SOURCE_FOLDER}"
        skipped++
        return
    }

    try {
        def imageData = entry.readImageData()
        def server    = imageData.getServer()
        def meta      = server.getMetadata()
        int nChannels = server.nChannels()

        if (channelList.size() != nChannels) {
            println "[WARN] ${rawName} -- source has ${channelList.size()} channels, image has ${nChannels}. Applying available names."
        }

        def newChannels = (0..<nChannels).collect { int i ->
            String name  = i < channelList.size() && channelList[i].name  ? channelList[i].name  : "Channel ${i+1}"
            int    color = i < channelList.size() && channelList[i].color ? channelList[i].color : ImageChannel.getDefaultChannelColor(i)
            new ImageChannel(name, color)
        }

        def newMeta = new ImageServerMetadata.Builder(meta)
                .channels(newChannels)
                .build()
        server.setMetadata(newMeta)

        if (SAVE_IMAGE_DATA) {
            entry.saveImageData(imageData)
        }

        println "[OK]   ${rawName}"
        newChannels.eachWithIndex { ch, i -> println "       Ch ${i}: ${ch.getName()}" }
        println ""
        updated++

    } catch (Exception e) {
        println "[FAIL] ${rawName} -- ${e.message}"
        e.printStackTrace()
        failed++
    }
}

println "---------------------------------------------------------"
println "Done.  Updated: ${updated}  |  Skipped: ${skipped}  |  Failed: ${failed}"
println "---------------------------------------------------------"
if (updated > 0 && SAVE_IMAGE_DATA) {
    println "Close and re-open each image to see updated channel names in the Brightness/Contrast panel."
}