// import_gold_standard.groovy  (QuPath 0.5+)
//
// Load the blinded gold-standard sample into the open image so it can be annotated in QuPath.
//
// Every sampled detection gets:
//   * PathClass "GS_review"     -> visible, selectable, filterable in the hierarchy
//   * measurement "gs_review_id" -> the blinded review id, shown in the measurement table
// Nothing else is set: no pipeline call is imported, so the annotation stays blind.
//
// WORKFLOW
//   1. Run this script on the donor's image.
//   2. Analyze > Show detection measurements, sort by gs_review_id.
//   3. For each cell: inspect the channels, then set its class from the fixed vocabulary below
//      (right-click > Set class, or select + press the class's keyboard shortcut).
//      Vocabulary: Immune, Epithelial, Endothelial, Neural, Mesenchymal, Other, Unsure.
//      "Unsure" is a legitimate answer and must not be avoided — a forced guess is worse than a
//      recorded uncertainty, and the scorer reports the unsure rate separately.
//   4. Run export_gold_standard.groovy to write your answers back out.
//
// Two annotators must work in SEPARATE QuPath projects (or separate copies of the .qpdata) and must
// not see each other's classes — inter-annotator kappa is meaningless otherwise.
//
// EDIT `folder` to point at <data_dir>/gold_standard_v1/qupath.

// ---- config ----------------------------------------------------------------
def folder = new File("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/gold_standard_v1/qupath")
// ----------------------------------------------------------------------------

// Same colours as config.COMPARTMENT_COLORS, so a reviewer moving between the gold-standard project
// and a pipeline overlay is not reading two different colour languages.
def COLORS = [
    "Immune"     : [34, 136, 51],    "Epithelial" : [68, 119, 170],
    "Endothelial": [102, 204, 238],  "Neural"     : [181, 131, 141],
    "Mesenchymal": [170, 51, 119],   "Other"      : [187, 187, 187],
    "Unsure"     : [240, 200, 0],    "GS_review"  : [255, 0, 255],
]
COLORS.each { name, c -> getPathClass(name).setColor(c[0], c[1], c[2]) }

if (!folder.isDirectory()) { println "ERROR: folder not found: ${folder}"; return }
def csvFiles = folder.listFiles().findAll {
    it.isFile() && it.getName().startsWith("gold_standard_") && it.getName().endsWith(".csv")
}
if (csvFiles.isEmpty()) { println "ERROR: no gold_standard_*.csv in ${folder}"; return }

def detById = [:]
getDetectionObjects().each { detById[it.getID().toString()] = it }
println "Open image has ${detById.size()} detections."

int matched = 0, missing = 0
csvFiles.each { file ->
    def reader = new BufferedReader(new FileReader(file))
    reader.readLine()                       // header: object_id,review_id,image
    String line
    while ((line = reader.readLine()) != null) {
        def f = line.split(',')
        def obj = detById[f[0]]
        if (obj == null) { missing++; continue }      // belongs to another image
        obj.setPathClass(getPathClass("GS_review"))
        // review_id is GS#####; store the numeric part so the measurement table sorts correctly.
        obj.getMeasurementList().put("gs_review_id", f[1].replaceAll("[^0-9]", "").toDouble())
        matched++
    }
    reader.close()
}
fireHierarchyUpdate()
println "Flagged ${matched} detections for gold-standard review (${missing} rows for other images)."
println "Next: annotate each one, then run export_gold_standard.groovy."
