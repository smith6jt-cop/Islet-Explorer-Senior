// import_redsea_lever_measurements.groovy   (QuPath 0.5+)
// -----------------------------------------------------------------------------
// Load per-cell REDSEA lever values (uncorrected / baseline-alpha1 / receptivity-gate corrected
// E_cadherin & CD31, the over-subtraction %, the gate weight, and the E_cad+ neighbor fraction)
// onto the detections of the CURRENTLY OPEN image, matched by UUID (PathObject.getID() == object_id).
//
// Produced by scratchpad/redsea_qupath_export.py:
//   figures/redsea_validation/qupath_redsea_lever_<donor>.csv
//   header: object_id,image,Ecad_uncorr,Ecad_base,Ecad_gate,Ecad_removed_base_pct,
//           Ecad_gate_w,epi_frac,CD31_uncorr,CD31_base,CD31_gate
//
// AFTER RUNNING, verify visually:
//   * Measure > Show measurement maps  -> compare "Ecad_uncorr" vs "Ecad_base" vs "Ecad_gate".
//     Baseline (alpha=1) dims the SOLID-EPITHELIUM interior; the gate keeps it. Isolated / dim-interior
//     "victim" cells should look the SAME (cleaned) under both.
//   * Colour by "Ecad_removed_base_pct" -> the over-subtraction map (bright = heavily subtracted;
//     should track epithelium interior, i.e. "epi_frac").
//   * Colour by "Ecad_gate_w" -> receptivity weight (0 = protected real expresser .. 1 = fully cleaned).
// -----------------------------------------------------------------------------

// EDIT this path if your repo/donor differ:
def csvFile = new File("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/figures/redsea_validation/qupath_redsea_lever_6374.csv")

if (!csvFile.isFile()) { println "ERROR: CSV not found: ${csvFile}"; return }

// index detections of the open image by UUID
def detById = [:]
getDetectionObjects().each { detById[it.getID().toString()] = it }
println "Open image has ${detById.size()} detections."

def reader = new BufferedReader(new FileReader(csvFile))
def header = reader.readLine().split(',')
def measCols = (2..<header.length).toList()   // skip object_id (0) and image (1)
println "Importing REDSEA measurements: ${measCols.collect{ header[it] }}"

int matched = 0, missing = 0
String line
while ((line = reader.readLine()) != null) {
    def f = line.split(',', -1)               // -1 keeps empty trailing/interior fields
    def obj = detById[f[0]]
    if (obj == null) { missing++; continue }  // row belongs to a different image
    def ml = obj.getMeasurementList()
    measCols.each { j ->
        if (j < f.length && f[j] != null && !f[j].isEmpty()) {
            try { ml.put(header[j] as String, f[j] as double) } catch (ex) { /* non-numeric -> skip */ }
        }
    }
    matched++
}
reader.close()
fireHierarchyUpdate()
println "Added REDSEA lever measurements to ${matched} detections (${missing} unmatched / other image)."
println "Now: Measure > Show measurement maps, and compare Ecad_uncorr / Ecad_base / Ecad_gate."
