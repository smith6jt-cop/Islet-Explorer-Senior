// set_broad_lineage_class.groovy
// After importing the Step-1 phenotype CSVs with "import clusters and cell neighborhoods.groovy"
// (which adds a numeric `broad_lineage` measurement, 1-6), run this to set each detection's
// PathClass from that code so "view phenotypes.groovy" can colour / toggle lineages on the image.
//
// Run per image (or "Run for project"). Legend matches data/phenotype/qupath_import/broad_lineage_legend.csv.

// NOTE: superseded by import_broad_lineage.groovy (UUID match → direct class). Kept in sync for
// the legacy measurement-based path. Codes match export_phenotype_for_qupath.py LNAMES order.
def legend = [1:'Epithelial', 2:'Fibroblast', 3:'Immune', 4:'Endocrine', 5:'Endothelial', 6:'Muscle']

def readMeas = { obj, String key ->
    def ml = obj.getMeasurementList()
    try { return ml.get(key) } catch (Throwable t) { }
    try { return ml.getMeasurementValue(key) } catch (Throwable t) { return Double.NaN }
}

int n = 0, miss = 0
getDetectionObjects().each { d ->
    double v = readMeas(d, 'broad_lineage')
    if (Double.isNaN(v)) { miss++; return }
    def name = legend[(int) Math.round(v)]
    if (name != null) { d.setPathClass(getPathClass(name)); n++ }
}
fireHierarchyUpdate()
println "Set broad_lineage class on ${n} detections (${miss} had no broad_lineage measurement)."
