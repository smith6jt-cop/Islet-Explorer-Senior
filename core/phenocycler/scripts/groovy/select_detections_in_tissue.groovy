/**
 * Select only the detections inside the "Tissue" annotation.
 *
 * BY CENTROID, NOT BY PARENT -- and that distinction is the whole point of this script.
 *
 * QuPath's `Parent` column looks like the obvious way to do this (`it.getParent()` up the
 * hierarchy until something is classified Tissue). It is wrong, and it was measured wrong on this
 * project on 2026-09-15: for the islet annotations the Parent labelled 11,556 cells `peri` that
 * sit up to 74.6 um INSIDE a core outline -- 2.7 % of islet-core cells cohort-wide and 35.3 % on
 * 6539. Anything built on the Parent inherits that. The centroid test asks the geometry directly
 * and agrees with the signed distance to sub-pixel precision.
 *
 * `contains(x, y)` is the same semantics the Python side uses (`shapely` `within`), so a cell
 * selected here is a cell the analysis will also call inside the tissue.
 *
 * Run in QuPath's script editor with the image open. Run for a whole project with
 * Run > Run for project.
 */

String TISSUE_CLASS = "Tissue"

def tissue = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() == TISSUE_CLASS && it.getROI() != null
}
if (tissue.isEmpty()) {
    throw new IllegalStateException("No '${TISSUE_CLASS}' annotation in this image.")
}
def rois = tissue.collect { it.getROI() }

def all = getDetectionObjects()
def inside = all.findAll { d ->
    def r = d.getROI()
    r != null && rois.any { it.contains(r.getCentroidX(), r.getCentroidY()) }
}

resetSelection()
selectObjects(inside)

def name = getProjectEntry()?.getImageName() ?: getCurrentServer()?.getMetadata()?.getName()
println "[tissue] ${name}"
println "[tissue] ${tissue.size()} '${TISSUE_CLASS}' annotation(s)"
println "[tissue] ${inside.size()} of ${all.size()} detections inside " +
        "(${all.size() - inside.size()} outside, ${String.format('%.2f', 100.0 * (all.size() - inside.size()) / Math.max(all.size(), 1))}%)"
