/**
 * Export every tissue-region annotation as a GeoJSON FeatureCollection.
 *
 * Supersedes export_regions_geojson.groovy, which filtered to Follicle + PALS
 * only. That filter is why Trabeculae polygons do not exist for the
 * Phenocycler donors — and Trabeculae is the one region the cell-derived layer
 * cannot supply, because no marker in the panel identifies it.
 *
 * Writes to analysis/geojson_regions/ rather than analysis/geojson/, so the
 * existing exports are superseded rather than overwritten.
 *
 * Coordinates are image pixels, matching the old exports. Convert with the
 * per-image pixel size — note 1901HBMP004 and HDL073 are 0.374 µm/px even
 * though their OME-XML declares 0.5078 (a stale header from
 * convert_to_ome_tiff.py); the CODEX sections are 0.3774.
 *
 *   QuPath script -p project.qpproj scripts/export_all_regions_geojson.groovy
 *
 * Read-only: no --save needed, and the hierarchy is not modified.
 */

import qupath.lib.io.GsonTools
import com.google.gson.stream.JsonWriter

// Every region class the analysis consumes. SmallVessel and LargeVessel are
// excluded: vessel objects are taken from AnnotationsFinal.csv, which carries
// their centroid and area — all that count density and area fraction need.
//
// Both spellings of Trabecula are listed, and so is Peripheral_White_Pulp: the
// CODEX sections use the singular and carry the peripheral white pulp that
// build_region_layer_qupath.py unions into the perifollicular band. Those
// sections are not in this project today, so nothing here can drop them — but
// the filter has to match CLASS_TO_REGION or a future export silently would.
def REGION_CLASSES = ['Follicle', 'PALS', 'RedPulp',
                      'Trabeculae', 'Trabecula',
                      'Peripheral_White_Pulp'] as Set

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def annotations = hierarchy.getAnnotationObjects()

def filtered = annotations.findAll { anno ->
    REGION_CLASSES.contains(anno.getPathClass()?.toString())
}

def imageName = getProjectEntry().getImageName()

if (filtered.isEmpty()) {
    println "${imageName}: no region annotations — skipping"
    return
}

def project = getProject()
def outputDir = new File(project.getPath().getParent().toFile(), "analysis/geojson_regions")
outputDir.mkdirs()

def baseName = imageName
    .replaceAll(/\.ome\.tiff?$/, '')
    .replaceAll(/\.tiff?$/, '')
def outputFile = new File(outputDir, "${baseName}.geojson")

def gson = GsonTools.getInstance(true)
def writer = new JsonWriter(new java.io.FileWriter(outputFile))
writer.setIndent("  ")
writer.beginObject()
writer.name("type").value("FeatureCollection")
writer.name("features")
writer.beginArray()
for (obj in filtered) {
    gson.toJson(obj, obj.getClass(), writer)
}
writer.endArray()
writer.endObject()
writer.close()

def counts = REGION_CLASSES.collectEntries { cls ->
    [cls, filtered.count { it.getPathClass()?.toString() == cls }]
}.findAll { it.value > 0 }

println "${imageName}: ${counts} → ${outputFile.getName()}"
