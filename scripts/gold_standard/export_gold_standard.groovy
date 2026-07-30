// export_gold_standard.groovy  (QuPath 0.5+)
//
// Write one annotator's gold-standard calls back out, for scripts/gold_standard/score.py.
//
// Exports every detection carrying a gs_review_id measurement, with the PathClass the annotator set.
// A cell still classed "GS_review" is UNANNOTATED and is exported as an empty compartment, so the
// scorer can report completion honestly rather than silently treating it as a call.
//
// Pass the immutable bundle directory and annotator ID as script arguments:
//   --args "/absolute/path/to/data/gold_standard_v2/<run_id>,A"
// Run once per image. Two reviewers must use different IDs and separate projects.

String bundlePath = args != null && args.size() > 0 ? args[0]?.trim() : ""
String annotator = args != null && args.size() > 1 ? args[1]?.trim() : ""
if (bundlePath == null || bundlePath.isEmpty() ||
        annotator == null || !(annotator ==~ /[A-Za-z0-9._-]+/)) {
    println "ERROR: pass the audit bundle directory and a safe annotator ID."
    return
}
def bundle = new File(bundlePath)
def manifestFile = new File(bundle, "manifest.json")
if (!manifestFile.isFile()) {
    println "ERROR: audit manifest not found under ${bundle}."
    return
}
def outDir = new File(bundle, "completed")

def VOCAB = ["Immune", "Epithelial", "Endothelial", "Neural", "Mesenchymal", "Other", "Unsure"] as Set

outDir.mkdirs()
def imageName = getCurrentImageData().getServer().getMetadata().getName()
def out = new File(outDir, "annotations_${annotator}_${imageName.replaceAll('[^A-Za-z0-9_.-]', '_')}.csv")
if (out.exists()) {
    println "ERROR: refusing to overwrite existing annotation export: ${out}"
    return
}

int n = 0, unannotated = 0, offVocab = 0
def csvField = { value ->
    '"' + (value == null ? "" : value.toString()).replace('"', '""') + '"'
}
out.withPrintWriter { w ->
    w.println("annotator,image,object_id,gs_review_id,compartment")
    getDetectionObjects().each { obj ->
        def id = obj.getMeasurementList().get("gs_review_id")
        if (id == null || Double.isNaN(id)) return
        def cls = obj.getPathClass()?.getName() ?: ""
        if (cls == "GS_review" || cls == "") { unannotated++; cls = "" }
        else if (!VOCAB.contains(cls)) { offVocab++ }
        w.println([
            annotator,
            imageName,
            obj.getID(),
            "GS${String.format('%05d', (int) id)}",
            cls,
        ].collect(csvField).join(","))
        n++
    }
}
println "Wrote ${n} rows -> ${out}"
println "  unannotated (still GS_review): ${unannotated}"
if (offVocab > 0)
    println "  WARNING: ${offVocab} rows carry a class outside the vocabulary; the scorer will reject them."
