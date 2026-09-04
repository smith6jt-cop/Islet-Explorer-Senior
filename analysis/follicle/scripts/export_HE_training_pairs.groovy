/**
 * Export (hematoxylin image, binary mask) training pairs for the H&E follicle U-Net.
 *
 * For each project entry this writes, into data/:
 *   - images/<name>.png  — a SINGLE-CHANNEL hematoxylin (H) image (TYPE_BYTE_GRAY,
 *                          high H = dark), deconvolved from the downsampled RGB H&E.
 *   - tissue/<name>.png  — a uint8 tissue mask (255 = tissue) computed from RGB as
 *                          mean(R,G,B) < 220, so red pulp is correctly tissue.
 *   - masks/<name>.png   — the uint8 Follicle mask (255 = Follicle, 0 = background).
 *   - meta/<name>.json   — pixel-size metadata + the H-deconvolution transform.
 *
 *   QuPath script -p project.qpproj scripts/export_HE_training_pairs.groovy
 *
 * Resolution is chosen by TARGET_UM_PER_PX so all donors land on a common scale
 * regardless of scanner. Default 2 µm/px ≈ 4× downsample from a 0.5 µm scan.
 *
 * PARITY REQUIREMENT: the RGB→H deconvolution below MUST stay numerically
 * identical to `rgb_to_h` in scripts/HE_unet_utils.py (the inference-side
 * counterpart). Both use the reference H&E OD unit vectors
 *   H = [0.5626, 0.7201, 0.4062], E = [0.2159, 0.8012, 0.5581] (normalized),
 *   OD = -log10(max(channel, 1) / 255),
 *   cH = w·OD where w is the first row of (MMᵀ)⁻¹M, M = [Ĥ; Ê],
 *   cH clipped to [0, 2.0]; intensity = floor(255·(1 − cH/2.0)).
 * Any change to these constants must be made in BOTH files in the same commit.
 */
import javax.imageio.ImageIO
import java.awt.Color
import java.awt.RenderingHints
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import qupath.lib.regions.RegionRequest

double TARGET_UM_PER_PX = 2.0

def imageData = getCurrentImageData()
if (imageData == null) {
    println "No image data — skipping"
    return
}
def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()

def cal = server.getPixelCalibration()
double pxSizeUm = cal.getAveragedPixelSizeMicrons()
if (Double.isNaN(pxSizeUm) || pxSizeUm <= 0) {
    println "Missing pixel calibration — skipping ${server.getMetadata().getName()}"
    return
}
double downsample = TARGET_UM_PER_PX / pxSizeUm
int w = (int) Math.ceil(server.getWidth() / downsample)
int h = (int) Math.ceil(server.getHeight() / downsample)

def follicles = hierarchy.getAnnotationObjects().findAll {
    it.getPathClass()?.getName() == 'Follicle'
}
if (follicles.isEmpty()) {
    println "No Follicle annotations in ${server.getMetadata().getName()} — skipping"
    return
}

def outDir = new File(buildFilePath(PROJECT_BASE_DIR, 'data'))
new File(outDir, 'images').mkdirs()
new File(outDir, 'masks').mkdirs()
new File(outDir, 'tissue').mkdirs()
new File(outDir, 'meta').mkdirs()

// Strip any QuPath suffix (" - resolution #N", " - normalized", etc.) first, then extensions.
def baseName = server.getMetadata().getName()
    .replaceAll(\ .*$\, '')
    .replaceAll(\.ome\.tiff?$\, '')
    .replaceAll(\.tiff?$\, '')

// --- Read RGB region ------------------------------------------------------
def request = RegionRequest.createInstance(server.getPath(), downsample,
    0, 0, server.getWidth(), server.getHeight())
def rgbImage = server.readBufferedImage(request)

// Use the actual rendered image dimensions for every output so the PNGs match
// exactly (QuPath's downsampling can round differently from ceil(w/downsample)).
w = rgbImage.getWidth()
h = rgbImage.getHeight()

// --- RGB → hematoxylin deconvolution (parity with rgb_to_h) ---------------
double[] Href = [0.5626, 0.7201, 0.4062]
double[] Eref = [0.2159, 0.8012, 0.5581]
double hnorm = Math.sqrt(Href[0] * Href[0] + Href[1] * Href[1] + Href[2] * Href[2])
double enorm = Math.sqrt(Eref[0] * Eref[0] + Eref[1] * Eref[1] + Eref[2] * Eref[2])
double[] Hn = [Href[0] / hnorm, Href[1] / hnorm, Href[2] / hnorm]
double[] En = [Eref[0] / enorm, Eref[1] / enorm, Eref[2] / enorm]
double hh = Hn[0] * Hn[0] + Hn[1] * Hn[1] + Hn[2] * Hn[2]
double he = Hn[0] * En[0] + Hn[1] * En[1] + Hn[2] * En[2]
double ee = En[0] * En[0] + En[1] * En[1] + En[2] * En[2]
double det = hh * ee - he * he
double a00 = ee / det
double a01 = -he / det
double w0 = a00 * Hn[0] + a01 * En[0]
double w1 = a00 * Hn[1] + a01 * En[1]
double w2 = a00 * Hn[2] + a01 * En[2]
double OD_WHITE = 255.0
double CH_MAX = 2.0

// OD lookup over the 8-bit range — matches the per-pixel LUT in rgb_to_h.
double[] odLut = new double[256]
for (int v = 0; v < 256; v++) {
    odLut[v] = -Math.log10(Math.max(v, 1) / OD_WHITE)
}

int[] pixels = rgbImage.getRGB(0, 0, w, h, null, 0, w)
def hImage = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def tissueImage = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def hRaster = hImage.getRaster()
def tRaster = tissueImage.getRaster()
for (int yy = 0; yy < h; yy++) {
    for (int xx = 0; xx < w; xx++) {
        int argb = pixels[yy * w + xx]
        int r = (argb >> 16) & 0xff
        int g = (argb >> 8) & 0xff
        int b = argb & 0xff
        double cH = w0 * odLut[r] + w1 * odLut[g] + w2 * odLut[b]
        if (cH < 0.0) cH = 0.0
        if (cH > CH_MAX) cH = CH_MAX
        double inten = 255.0 * (1.0 - cH / CH_MAX)
        if (inten < 0.0) inten = 0.0
        if (inten > 255.0) inten = 255.0
        hRaster.setSample(xx, yy, 0, (int) Math.floor(inten))
        double mean = (r + g + b) / 3.0
        tRaster.setSample(xx, yy, 0, mean < 220 ? 255 : 0)
    }
}
ImageIO.write(hImage, 'PNG', new File(outDir, "images/${baseName}.png"))
ImageIO.write(tissueImage, 'PNG', new File(outDir, "tissue/${baseName}.png"))

// --- Mask -----------------------------------------------------------------
def mask = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def g2d = mask.createGraphics()
g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF)
g2d.setColor(Color.BLACK)
g2d.fillRect(0, 0, w, h)
g2d.setColor(Color.WHITE)
g2d.setTransform(AffineTransform.getScaleInstance(1.0 / downsample, 1.0 / downsample))
for (ann in follicles) {
    g2d.fill(ann.getROI().getShape())
}
g2d.dispose()
ImageIO.write(mask, 'PNG', new File(outDir, "masks/${baseName}.png"))

// --- Sidecar --------------------------------------------------------------
def meta = new File(outDir, "meta/${baseName}.json")
meta.text = """\
{
  "image": "${baseName}",
  "source_name": "${server.getMetadata().getName()}",
  "image_type": "hematoxylin",
  "channels": 1,
  "tissue_mask": "tissue/${baseName}.png",
  "h_transform": {
    "h_vector": [0.5626, 0.7201, 0.4062],
    "e_vector": [0.2159, 0.8012, 0.5581],
    "od_white": ${OD_WHITE},
    "ch_max": ${CH_MAX},
    "polarity": "high_H_dark"
  },
  "native_px_size_um": ${pxSizeUm},
  "target_px_size_um": ${TARGET_UM_PER_PX},
  "downsample": ${downsample},
  "width_native_px": ${server.getWidth()},
  "height_native_px": ${server.getHeight()},
  "width_export_px": ${w},
  "height_export_px": ${h},
  "n_follicle_annotations": ${follicles.size()}
}
"""

println "${baseName}: ${follicles.size()} follicles, ${w}×${h} @ ${TARGET_UM_PER_PX} µm/px (hematoxylin)"
