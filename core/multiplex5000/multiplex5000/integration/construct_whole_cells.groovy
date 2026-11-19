/* Execute the checksum-pinned PubApp constructor in a fresh scratch project.
 * No source project is opened, no source cells or images are saved or modified.
 * The scientific constructor and common library are supplied by the run receipt.
 */
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import org.locationtech.jts.io.WKBWriter
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryFactory
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.LinearRing
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.util.GeometryFixer
import org.locationtech.jts.geom.util.PolygonExtracter
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ImportCustomizer

def cfg = GsonTools.getInstance().fromJson(new File(args[0]).text, Map)
setPixelSizeMicrons((double)cfg.pixel_size_um, (double)cfg.pixel_size_um)
clearAllObjects()
// Parse native coordinates as JTS geometry, then create native QuPath cells.
// PathIO's generic GeoJSON ROI adapter rounds to .01 image units and may fail
// precision reduction for valid vendor polygons. The production constructor
// consumes native cells directly, so do the same without that extra rounding.
def gf = new GeometryFactory()
def ring = { coords -> gf.createLinearRing(coords.collect { xy -> new Coordinate((double)xy[0], (double)xy[1]) } as Coordinate[]) }
def polygon = { coords -> gf.createPolygon(ring(coords[0]), coords.drop(1).collect { ring(it) } as LinearRing[]) }
def nativeROI = { raw ->
    def g
    if (raw.type == 'Polygon') g = polygon(raw.coordinates)
    else if (raw.type == 'MultiPolygon') g = gf.createMultiPolygon(raw.coordinates.collect { polygon(it) } as Polygon[])
    else throw new IllegalArgumentException('Polygonal native cell/proposal required: ' + raw.type)
    if (!g.isValid()) g = GeometryFixer.fix(g)
    def polygons = []
    PolygonExtracter.getPolygons(g, polygons)
    g = g.getFactory().buildGeometry(polygons)
    return GeometryTools.geometryToROI(g, qupath.lib.regions.ImagePlane.getDefaultPlane())
}
def rawCells = GsonTools.getInstance().fromJson(new File(cfg.cells_geojson).text, Map)
def cells = rawCells.features.collect { f ->
    if (f.nucleusGeometry == null) throw new IllegalArgumentException('No nucleus in construction input')
    def cell = PathObjects.createCellObject(nativeROI(f.geometry), nativeROI(f.nucleusGeometry), null, null)
    cell.setID(UUID.fromString(f.id))
    return cell
}
addObjects(cells)
def rawProposals = GsonTools.getInstance().fromJson(new File(cfg.proposals_geojson).text, Map)
def proposals = rawProposals.features.collect { f ->
    def o = PathObjects.createAnnotationObject(nativeROI(f.geometry), getPathClass('Islet'))
    o.setID(UUID.fromString(f.id))
    return o
}
addObjects(proposals)
binding.setVariable('PROJECT_BASE_DIR', getProject().getPath().getParent().toString())
def compiler = new CompilerConfiguration()
compiler.addCompilationCustomizers(new ImportCustomizer().addStaticStars('qupath.lib.gui.scripting.QPEx'))
new GroovyShell(this.class.classLoader, binding, compiler).evaluate(new File(cfg.constructor_script))
def result = binding.hasVariable('V3_RESULT') ? binding.getVariable('V3_RESULT') : null
if (result == null || !result.ok) throw new IllegalStateException('Constructor did not report success: ' + result)
def annotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() in ['Islet', 'EndocrineCluster', 'EndocrineSingle'] }
exportObjectsToGeoJson(annotations, cfg.output_geojson, 'FEATURE_COLLECTION')
def writer = new WKBWriter()
def precise = annotations.collect { o ->
    [id: o.getID().toString(), object_class: o.getPathClass().toString(),
     n_whole_cells: o.getMeasurementList().get('Whole cells'),
     geometry_wkb_hex: WKBWriter.toHex(writer.write(o.getROI().getGeometry()))]
}
def imported = cells.collect { c -> [qupath_uuid:c.getID().toString(), geometry_wkb_hex:WKBWriter.toHex(writer.write(c.getROI().getGeometry()))] }
new File(new File(cfg.output_precise).getParentFile(), 'imported_cells_precise.json').text = GsonTools.getInstance().toJson(imported) + '\n'
new File(cfg.output_precise).text = GsonTools.getInstance().toJson(precise) + '\n'
new File(cfg.output_receipt).text = GsonTools.getInstance(true).toJson(result) + '\n'
println('COMPLETE_WHOLE_CELL_OBJECTS ' + annotations.size())
