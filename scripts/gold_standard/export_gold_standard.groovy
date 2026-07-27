// export_gold_standard.groovy  (QuPath 0.5+)
//
// Write one annotator's gold-standard calls back out, for scripts/gold_standard/score.py.
//
// Exports every detection carrying a gs_review_id measurement, with the PathClass the annotator set.
// A cell still classed "GS_review" is UNANNOTATED and is exported as an empty compartment, so the
// scorer can report completion honestly rather than silently treating it as a call.
//
// EDIT `annotator` and `outDir` before running. Run once per image; the file is appended per image so
// a multi-image project produces one CSV per annotator.

// ---- config ----------------------------------------------------------------
def annotator = "A"                     // "A" or "B" — must differ between the two reviewers
def outDir = new File("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/gold_standard_v1/completed")
// ----------------------------------------------------------------------------

def VOCAB = ["Immune", "Epithelial", "Endothelial", "Neural", "Mesenchymal", "Other", "Unsure"] as Set

outDir.mkdirs()
def imageName = getCurrentImageData().getServer().getMetadata().getName()
def out = new File(outDir, "annotations_${annotator}_${imageName.replaceAll('[^A-Za-z0-9_.-]', '_')}.csv")

int n = 0, unannotated = 0, offVocab = 0
out.withPrintWriter { w ->
    w.println("annotator,image,object_id,gs_review_id,compartment")
    getDetectionObjects().each { obj ->
        def id = obj.getMeasurementList().get("gs_review_id")
        if (id == null || Double.isNaN(id)) return
        def cls = obj.getPathClass()?.getName() ?: ""
        if (cls == "GS_review" || cls == "") { unannotated++; cls = "" }
        else if (!VOCAB.contains(cls)) { offVocab++ }
        w.println("${annotator},${imageName},${obj.getID()},GS${String.format('%05d', (int) id)},${cls}")
        n++
    }
}
println "Wrote ${n} rows -> ${out}"
println "  unannotated (still GS_review): ${unannotated}"
if (offVocab > 0)
    println "  WARNING: ${offVocab} rows carry a class outside the vocabulary; the scorer will reject them."
