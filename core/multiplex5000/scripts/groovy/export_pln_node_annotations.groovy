/*
 * Export the annotations of one QuPath image as GeoJSON, straight from its .qpdata.
 *
 * The pLN node outlines live only in the binary .qpdata of the /data/IO60pLN QuPath project and
 * had never been exported. This reads one data file with PathIO.readHierarchy, which needs no
 * image server -- the qptiff is never opened -- and writes every annotation on that image to
 * GeoJSON in pixel coordinates.
 *
 * It reads and does not write to the project. Run it with QuPath's CLI, which does not save data
 * files unless asked:
 *
 *   for e in 1 2 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do
 *     JAVA_TOOL_OPTIONS="-Xmx48g" /home/smith6jt/QuPath/bin/QuPath script \
 *       --args /data/IO60pLN/data/$e/data.qpdata \
 *       --args /data/IO60pLN/results/node_geojson/entry$e.geojson \
 *       core/multiplex5000/scripts/groovy/export_pln_node_annotations.groovy
 *   done
 *
 * The stdout ANNOT lines are the geometry audit: ROI type, vertex count, pixel area, bounds and
 * child count per annotation. Keep them -- the child count is the only place QuPath's own
 * containment is recorded, and it is the cross-check for the geometric membership, never its
 * source.
 */
import qupath.lib.io.PathIO
import qupath.lib.io.PathIO.GeoJsonExportOptions

def inPath = java.nio.file.Paths.get(args[0])
def outPath = new File(args[1])
outPath.getParentFile().mkdirs()

long t0 = System.currentTimeMillis()
def hierarchy = PathIO.readHierarchy(inPath)
long t1 = System.currentTimeMillis()
def annots = new ArrayList<>(hierarchy.getAnnotationObjects())
println "READ_MS=" + (t1 - t0) + " N_OBJECTS=" + hierarchy.nObjects() + " N_ANNOTATIONS=" + annots.size()
for (a in annots) {
    def r = a.getROI()
    def roiDesc = (r == null) ? "null" : String.valueOf(r.getRoiType()) + "/" + r.getClass().getSimpleName()
    def bounds = (r == null) ? "" : (r.getBoundsX() + "," + r.getBoundsY() + "," + r.getBoundsWidth() + "," + r.getBoundsHeight())
    println "ANNOT\tid=" + a.getID() + "\tclass=" + a.getPathClass() + "\tname=" + a.getName() +
            "\troi=" + roiDesc + "\tnPoints=" + (r == null ? -1 : r.getNumPoints()) +
            "\tareaPx=" + (r == null ? -1 : r.getArea()) + "\tbounds=" + bounds +
            "\tnChildren=" + a.nChildObjects()
}
PathIO.exportObjectsAsGeoJSON(outPath, annots, GeoJsonExportOptions.FEATURE_COLLECTION, GeoJsonExportOptions.EXCLUDE_MEASUREMENTS)
println "WROTE=" + outPath + " BYTES=" + outPath.length() + " TOTAL_MS=" + (System.currentTimeMillis() - t0)
