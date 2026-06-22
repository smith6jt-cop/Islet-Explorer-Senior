// Import per-image phenotype CSVs (CytoMAP-style) into their matching QuPath images.
// Each CSV: col0 = Centroid X µm, col1 = Centroid Y µm, col2 = Image name, cols 3+ = numeric
// measurements added to the detection at that centroid.
// Adapted for QuPath 0.7: hard-coded import folder, fast centroid hash-map matching,
// MeasurementList.put(...) (putMeasurement was removed in 0.7).

import qupath.lib.gui.measure.ObservableMeasurementTableData
import java.math.MathContext

// ---- config ----------------------------------------------------------------
// Quick verification: import only these donors (one small image = fast). [] / null = all.
def donorsToImport = ["6539"]
def folder = new File("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/phenotype/qupath_import")
// ----------------------------------------------------------------------------

def qupathProject = getProject()
if (qupathProject == null) { println "ERROR: no project open."; return }
if (!folder.isDirectory()) { println "ERROR: import folder not found: ${folder}"; return }

def csvFiles = folder.listFiles().findAll {
    it.isFile() && it.getName().startsWith("pheno_broad_") && it.getName().endsWith(".csv") &&
    (donorsToImport == null || donorsToImport.isEmpty() ||
     donorsToImport.any { d -> it.getName() == "pheno_broad_${d}.csv" })
}
if (csvFiles.isEmpty()) { println "ERROR: no matching pheno_broad_*.csv in ${folder}"; return }
println "Importing ${csvFiles.size()} CSV file(s) from ${folder}"

csvFiles.each { file ->
    def csvReader = new BufferedReader(new FileReader(file))
    def measurementNames = csvReader.readLine().split(',')   // header
    int length = measurementNames.size()
    println "Adding results from ${file.getName()}"

    def csv = []
    Set imageList = []
    String line
    while ((line = csvReader.readLine()) != null) {
        def toAdd = line.split(',')
        imageList.add(toAdd[2])
        csv << toAdd
    }
    csvReader.close()

    imageList.each { image ->
        def entry = qupathProject.getImageList().find { it.getImageName() == image }
        if (entry == null) { println "  NO PROJECT ENTRY FOR IMAGE ${image}"; return }

        def imageData = entry.readImageData()
        def hierarchy = imageData.getHierarchy()
        def objects = hierarchy.getDetectionObjects()
        def ob = new ObservableMeasurementTableData()
        ob.setImageData(imageData, objects)

        // Build a centroid -> object lookup ONCE (rounded to the CSV's QuPath precision),
        // then match each CSV row in O(1). Avoids a spatial query per cell.
        def lookup = [:]
        objects.each { o ->
            def k = round(ob.getNumericValue(o, "Centroid X µm")) + "_" +
                    round(ob.getNumericValue(o, "Centroid Y µm"))
            lookup[k] = o
        }

        def rows = csv.findAll { it[2] == image }
        int matched = 0, miss = 0
        rows.each { r ->
            def object = lookup[(r[0] as double) + "_" + (r[1] as double)]
            if (object == null) { miss++; return }
            def ml = object.getMeasurementList()
            for (int i = 3; i < length; i++) {
                ml.put(measurementNames[i], r[i] as double)
            }
            matched++
        }
        entry.saveImageData(imageData)
        println "  ${image}: matched ${matched}/${rows.size()} cells (${miss} unmatched)"
    }
}

println "Done!"

def round(double number) {
    BigDecimal bd = new BigDecimal(number)
    def result
    if (number < 100) {
        bd = bd.round(new MathContext(4))
        result = bd.doubleValue()
    } else if (number < 1000) {
        result = number.round(2)
    } else {
        result = number.round(1)
    }
    return result
}
