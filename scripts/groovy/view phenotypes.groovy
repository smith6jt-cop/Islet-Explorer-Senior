// view phenotypes.groovy — show/hide detection classes on the image (QuPath 0.7).
// Run AFTER set_broad_lineage_class.groovy has set each detection's PathClass.
// (QuPath 0.7 replaced OverlayOptions.hiddenClassesProperty() with setPathClassHidden(PathClass, boolean).)

import javafx.application.Platform
import javafx.geometry.Insets
import javafx.scene.Scene
import javafx.scene.control.CheckBox
import javafx.scene.control.Label
import javafx.scene.control.ScrollPane
import javafx.scene.layout.BorderPane
import javafx.scene.layout.GridPane
import javafx.stage.Stage
import javafx.beans.value.ChangeListener
import qupath.lib.gui.QuPathGUI
import static qupath.lib.gui.scripting.QPEx.getCurrentViewer

def viewer = getCurrentViewer()
if (viewer == null) { println "ERROR: no image open."; return }
def overlay = viewer.getOverlayOptions()

// distinct, non-null detection classes (flat classes — no parent/base-class splitting)
def classes = getDetectionObjects().collect { it.getPathClass() }.findAll { it != null }.unique()
classes = classes.sort { it.toString() }
if (classes.isEmpty()) { println "No classified detections. Run set_broad_lineage_class.groovy first."; return }
println "Classes: ${classes.collect { it.toString() }}"

// start with everything visible
classes.each { overlay.setPathClassHidden(it, false) }
viewer.repaint()

Platform.runLater {
    def grid = new GridPane()
    grid.setPadding(new Insets(12)); grid.setVgap(5); grid.setHgap(10)
    int row = 0
    grid.add(new Label("Show / hide phenotype classes"), 0, row++)

    def allBox = new CheckBox("All"); allBox.setSelected(true)
    grid.add(allBox, 0, row++)

    def boxes = []
    classes.each { pc ->
        def cb = new CheckBox(pc.toString())
        cb.setSelected(true)
        cb.selectedProperty().addListener({ o, ov, nv ->
            overlay.setPathClassHidden(pc, !nv)
            viewer.repaint()
        } as ChangeListener)
        boxes << cb
        grid.add(cb, 0, row++)
    }

    // "All" toggles every class box (each fires its own listener)
    allBox.selectedProperty().addListener({ o, ov, nv ->
        boxes.each { it.setSelected(nv) }
    } as ChangeListener)

    def stage = new Stage()
    stage.initOwner(QuPathGUI.getInstance().getStage())
    stage.setScene(new Scene(new BorderPane(new ScrollPane(grid))))
    stage.setTitle("Select classes to display")
    stage.setWidth(280); stage.setHeight(360); stage.setResizable(true)
    stage.show()
}
