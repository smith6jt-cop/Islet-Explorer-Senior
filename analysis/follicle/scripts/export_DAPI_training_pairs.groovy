/**
 * Export (DAPI image, Follicle mask, region label map) training tuples for the
 * DAPI/Phenocycler follicle U-Net.
 *
 * For each project entry this writes, into data_dapi/:
 *   - images/<name>.png   — SINGLE-CHANNEL uint8 DAPI, normalised and INVERTED
 *                           (high DAPI = dark), matching the polarity of the
 *                           H&E pipeline's hematoxylin images.
 *   - tissue/<name>.png   — uint8 tissue mask (255 = tissue), Otsu on raw DAPI.
 *   - masks/<name>.png    — uint8 Follicle mask (255 = Follicle), derived from
 *                           the region map so the two cannot disagree.
 *   - regions/<name>.png  — uint8 region label index (see REGION_CLASSES).
 *   - meta/<name>.json    — pixel-size metadata + the DAPI transform + per-class areas.
 *
 *   QuPath script -p project.qpproj scripts/export_DAPI_training_pairs.groovy
 *   QuPath script -p project.qpproj -a TARGET_UM_PER_PX=2.0 -a EXCLUDE=HDL018,HDL021 \
 *       scripts/export_DAPI_training_pairs.groovy
 *
 * Read-only with respect to the project — no --save needed.
 *
 * images/ + masks/ + tissue/ + meta/ is exactly the contract
 * HE_unet_utils.discover_samples() already consumes, so the H&E and DAPI
 * cohorts can be pooled into one training set. regions/ is the DAPI-only extra.
 *
 * PARITY REQUIREMENT: the histogram / Otsu / percentile / normalise block below
 * MUST stay numerically identical to `hist256`, `otsu_from_hist256`,
 * `percentile_from_hist256` and `dapi_to_uint8` in scripts/DAPI_unet_utils.py.
 * Every step is integer or double arithmetic over a 256-bin histogram precisely
 * so the two languages agree to the byte:
 *   bin      = value >>> 8
 *   otsu     = argmax_t  wB(t)*wF(t)*(meanB(t) - meanF(t))^2      (first max wins)
 *   tissue   = raw >= ((otsu + 1) << 8)   // otsu's t is the last BACKGROUND bin
 *   lo / hi  = low edge of the bin containing the 1.0 / 99.5 percentile OF TISSUE
 *   intensity= floor(255 * (1 - clip((v - lo)/(hi - lo), 0, 1)))
 * Any change must be made in BOTH files in the same commit;
 * scripts/check_dapi_parity.py re-verifies it against a real export.
 *
 * MEMORY: reads ONE channel, never the full 25-33 channel stack. The largest
 * slide (HDL052, 28800x58320 native) is ~7312x14805 at 2 µm/px; peak heap is
 * roughly 1 GB. Reading all channels instead would need ~7 GB, which is why the
 * single-channel extraction below matters.
 */
import javax.imageio.ImageIO
import java.awt.Color
import java.awt.RenderingHints
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import java.awt.image.DataBufferInt
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.regions.RegionRequest

// --- Args -----------------------------------------------------------------
double TARGET_UM_PER_PX = 2.0
String OUT_SUBDIR = 'data_dapi'
def EXCLUDE = [] as Set
double LO_PCT = 1.0
double HI_PCT = 99.5

if (binding.hasVariable('args') && args != null) {
    for (a in args) {
        def kv = (a as String).split('=', 2)
        if (kv.length != 2) continue
        switch (kv[0].trim().toUpperCase()) {
            case 'TARGET_UM_PER_PX': TARGET_UM_PER_PX = Double.parseDouble(kv[1]); break
            case 'OUT_DIR':          OUT_SUBDIR = kv[1].trim(); break
            case 'EXCLUDE':          EXCLUDE = (kv[1].split(',').collect { it.trim() }) as Set; break
            case 'LO_PCT':           LO_PCT = Double.parseDouble(kv[1]); break
            case 'HI_PCT':           HI_PCT = Double.parseDouble(kv[1]); break
        }
    }
}

// Paint order: later wins. Index 0 is background/unlabeled. Follicle is painted
// LAST so `regions == 5` and the binary `masks` PNG agree pixel-for-pixel.
// Keep this list in sync with REGION_CLASSES in scripts/DAPI_unet_utils.py.
def REGION_CLASSES = ['Background', 'RedPulp', 'Trabeculae', 'LargeVessel', 'PALS', 'Follicle']
int FOLLICLE_INDEX = REGION_CLASSES.indexOf('Follicle')

// --- Image ----------------------------------------------------------------
def imageData = getCurrentImageData()
if (imageData == null) {
    println "No image data — skipping"
    return
}
def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()
def sourceName = server.getMetadata().getName()

def cal = server.getPixelCalibration()
double pxSizeUm = cal.getAveragedPixelSizeMicrons()
if (Double.isNaN(pxSizeUm) || pxSizeUm <= 0) {
    println "Missing pixel calibration — skipping ${sourceName}"
    return
}
double downsample = TARGET_UM_PER_PX / pxSizeUm

// Strip QuPath's " - resolution #N" suffix first, then the extension. The extra
// .qptiff rule is why this differs from export_HE_training_pairs.groovy: the
// three .qptiff entries (HDL052/HDL070/HDL079) are not matched by /\.tiff?$/
// and would otherwise keep ".er.qptiff" in every output filename.
def baseName = sourceName
    .replaceAll(/ .*$/, '')
    .replaceAll(/\.ome\.tiff?$/, '')
    .replaceAll(/\.qptiff$/, '')
    .replaceAll(/\.tiff?$/, '')

if (EXCLUDE.contains(baseName) || EXCLUDE.any { baseName.startsWith(it) }) {
    println "${baseName}: excluded by EXCLUDE arg — skipping"
    return
}

// --- Annotations ----------------------------------------------------------
def byClass = [:]
for (name in REGION_CLASSES) {
    if (name == 'Background') continue
    byClass[name] = hierarchy.getAnnotationObjects().findAll {
        it.getPathClass()?.getName() == name
    }
}
if (byClass['Follicle'].isEmpty()) {
    println "No Follicle annotations in ${sourceName} — skipping"
    return
}

// --- Locate the DAPI channel ---------------------------------------------
// Channel names have been lost on some of these images before (hence
// RestoreChannelNames.groovy), so resolve by name and fall back to index 0 with
// a warning rather than trusting either on its own.
def channels = server.getMetadata().getChannels()
int dapiIndex = -1
for (int i = 0; i < channels.size(); i++) {
    if (channels[i].getName()?.toUpperCase()?.startsWith('DAPI')) { dapiIndex = i; break }
}
if (dapiIndex < 0) {
    println "  WARNING: ${baseName}: no channel named DAPI among ${channels.size()} " +
            "channels; falling back to channel 0"
    dapiIndex = 0
} else if (dapiIndex != 0) {
    println "  note: ${baseName}: DAPI is channel ${dapiIndex}, not 0"
}
def dapiChannelName = channels[dapiIndex].getName()

// Single-channel server: reads ONE band off disk instead of the whole stack.
def dapiServer = new TransformedServerBuilder(server)
    .extractChannels(dapiChannelName)
    .build()

def request = RegionRequest.createInstance(dapiServer.getPath(), downsample,
    0, 0, server.getWidth(), server.getHeight())
def dapiImage = dapiServer.readRegion(request)

// Trust the rendered raster, not ceil(width/downsample) — QuPath's downsampling
// can round differently, and every output has to agree pixel-for-pixel.
int w = dapiImage.getWidth()
int h = dapiImage.getHeight()
def dapiRaster = dapiImage.getRaster()

// --- Region label map -----------------------------------------------------
// Painted into TYPE_INT_RGB, NOT TYPE_BYTE_GRAY: Java2D converts sRGB colours
// through a linear-grey colour space when filling a grey image, so
// new Color(3,3,3) would NOT land on sample value 3. INT_RGB stores the sRGB
// bytes verbatim, so the low byte is exactly the class index.
def labelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB)
def lg = labelImage.createGraphics()
lg.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF)
lg.setColor(new Color(0, 0, 0))
lg.fillRect(0, 0, w, h)
lg.setTransform(AffineTransform.getScaleInstance(1.0 / downsample, 1.0 / downsample))
for (int ci = 1; ci < REGION_CLASSES.size(); ci++) {
    def anns = byClass[REGION_CLASSES[ci]] ?: []
    if (anns.isEmpty()) continue
    lg.setColor(new Color(ci, ci, ci))
    for (ann in anns) {
        lg.fill(ann.getROI().getShape())
    }
}
lg.dispose()
int[] labelData = ((DataBufferInt) labelImage.getRaster().getDataBuffer()).getData()

// --- Histogram / Otsu / percentiles (parity block) ------------------------
// Row-at-a-time bulk reads: one native getSamples call per row rather than a
// per-pixel virtual call, which matters at ~10^8 pixels, while keeping the
// working buffer at one row instead of a full 433 MB int[].
int[] row = new int[w]
long[] hist = new long[256]
for (int yy = 0; yy < h; yy++) {
    dapiRaster.getSamples(0, yy, w, 1, 0, row)
    for (int xx = 0; xx < w; xx++) {
        hist[(row[xx] & 0xffff) >>> 8]++
    }
}

// THREE-class Otsu, lower threshold. A spleen DAPI histogram has three modes —
// glass, red pulp, white pulp — and 2-class Otsu lands between red pulp and
// white pulp, producing a "tissue" mask of white pulp only. t1 is the
// glass/tissue split; t2 is kept for provenance.
// The threshold bin itself belongs to the BACKGROUND side, so tissue is
// `bin > otsuBin`; using `>=` would keep the whole glass bin inside the mask.
int[] otsu3 = otsu3FromHist(hist)
int otsuBin = otsu3[0]
int otsuUpperBin = otsu3[1]
long[] tissueHist = new long[256]
long tissueTotal = 0L
for (int b = 0; b < 256; b++) {
    tissueHist[b] = (b > otsuBin) ? hist[b] : 0L
    tissueTotal += tissueHist[b]
}
long histTotal = 0L
for (int b = 0; b < 256; b++) histTotal += hist[b]
if (tissueTotal <= 0L) {
    // Everything fell below the threshold — normalise over the whole image
    // rather than dividing by zero. Degenerate slide, but do not crash.
    for (int b = 0; b < 256; b++) tissueHist[b] = hist[b]
    tissueTotal = histTotal
}
int lo = percentileFromHist(tissueHist, LO_PCT)
int hi = percentileFromHist(tissueHist, HI_PCT)
if (hi <= lo) hi = lo + 256
double tissueFrac = histTotal > 0L ? ((double) tissueTotal / histTotal) : 0.0d

// --- Single pass: DAPI -> uint8 (inverted), tissue, regions, Follicle mask -
def dapiOut = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def tissueOut = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def regionOut = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def maskOut = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def dOut = dapiOut.getRaster()
def tOut = tissueOut.getRaster()
def rOut = regionOut.getRaster()
def mOut = maskOut.getRaster()
double span = (double) (hi - lo)
int tissueCut = (otsuBin + 1) << 8   // tissue is bin > otsuBin (see above)
long[] classCount = new long[REGION_CLASSES.size()]
for (int yy = 0; yy < h; yy++) {
    int rowBase = yy * w
    dapiRaster.getSamples(0, yy, w, 1, 0, row)
    for (int xx = 0; xx < w; xx++) {
        int v = row[xx] & 0xffff
        double x = (v - lo) / span
        if (x < 0.0d) x = 0.0d
        if (x > 1.0d) x = 1.0d
        dOut.setSample(xx, yy, 0, (int) Math.floor(255.0d * (1.0d - x)))
        tOut.setSample(xx, yy, 0, v >= tissueCut ? 255 : 0)
        int cls = labelData[rowBase + xx] & 0xff
        rOut.setSample(xx, yy, 0, cls)
        // Derived from the label map, so masks and regions cannot disagree.
        mOut.setSample(xx, yy, 0, cls == FOLLICLE_INDEX ? 255 : 0)
        classCount[cls]++
    }
}
labelImage.flush()
labelImage = null
labelData = null

// --- Write ----------------------------------------------------------------
def outDir = new File(buildFilePath(PROJECT_BASE_DIR, OUT_SUBDIR))
for (sub in ['images', 'masks', 'tissue', 'regions', 'meta']) {
    new File(outDir, sub).mkdirs()
}
ImageIO.write(dapiOut, 'PNG', new File(outDir, "images/${baseName}.png"))
ImageIO.write(tissueOut, 'PNG', new File(outDir, "tissue/${baseName}.png"))
ImageIO.write(regionOut, 'PNG', new File(outDir, "regions/${baseName}.png"))
ImageIO.write(maskOut, 'PNG', new File(outDir, "masks/${baseName}.png"))

// --- Sidecar --------------------------------------------------------------
double pxAreaUm2 = TARGET_UM_PER_PX * TARGET_UM_PER_PX
def areaEntries = []
for (int ci = 0; ci < REGION_CLASSES.size(); ci++) {
    areaEntries << "\"${REGION_CLASSES[ci]}\": ${(classCount[ci] * pxAreaUm2) as long}"
}
def classListJson = REGION_CLASSES.collect { "\"${it}\"" }.join(', ')

def meta = new File(outDir, "meta/${baseName}.json")
meta.text = """\
{
  "image": "${baseName}",
  "source_name": "${sourceName}",
  "image_type": "dapi",
  "channels": 1,
  "dapi_channel_index": ${dapiIndex},
  "dapi_channel_name": "${dapiChannelName}",
  "tissue_mask": "tissue/${baseName}.png",
  "region_map": "regions/${baseName}.png",
  "region_classes": [${classListJson}],
  "dapi_transform": {
    "otsu_bin": ${otsuBin},
    "otsu_upper_bin": ${otsuUpperBin},
    "lo": ${lo},
    "hi": ${hi},
    "lo_pct": ${LO_PCT},
    "hi_pct": ${HI_PCT},
    "tissue_frac": ${String.format('%.6f', tissueFrac)},
    "invert": true,
    "polarity": "high_DAPI_dark"
  },
  "native_px_size_um": ${pxSizeUm},
  "target_px_size_um": ${TARGET_UM_PER_PX},
  "downsample": ${downsample},
  "width_native_px": ${server.getWidth()},
  "height_native_px": ${server.getHeight()},
  "width_export_px": ${w},
  "height_export_px": ${h},
  "n_follicle_annotations": ${byClass['Follicle'].size()},
  "class_area_um2": {${areaEntries.join(', ')}}
}
"""

def foll = classCount[FOLLICLE_INDEX]
def pctFoll = w * h > 0 ? (100.0d * foll / ((long) w * h)) : 0.0d
println "${baseName}: ${byClass['Follicle'].size()} Follicle annotations, " +
        "${String.format('%.2f', pctFoll)}% Follicle pixels, ${w}x${h} @ " +
        "${TARGET_UM_PER_PX} µm/px (DAPI ch${dapiIndex}, lo=${lo} hi=${hi} " +
        "otsu=${otsuBin}, tissue=${String.format('%.1f', tissueFrac * 100)}%)"
if (tissueFrac < 0.05d) {
    println "  WARNING: ${baseName}: tissue mask covers only " +
            "${String.format('%.1f', tissueFrac * 100)}% of the image. The " +
            "glass/tissue split likely failed, which biases negative-patch " +
            "sampling and the lo/hi window. Inspect tissue/${baseName}.png."
}

// --- Parity helpers -------------------------------------------------------
// Mirrors otsu_from_hist256 in DAPI_unet_utils.py. Strict '>' so the FIRST
// maximum wins, matching numpy's argmax.
static int otsuFromHist(long[] hist) {
    double total = 0.0d
    double sumAll = 0.0d
    for (int i = 0; i < 256; i++) {
        total += hist[i]
        sumAll += ((double) hist[i]) * i
    }
    if (total <= 0.0d) return 0
    double wB = 0.0d
    double sumB = 0.0d
    double best = -1.0d
    int bestT = 0
    for (int t = 0; t < 256; t++) {
        wB += hist[t]
        sumB += ((double) hist[t]) * t
        double wF = total - wB
        if (wB <= 0.0d || wF <= 0.0d) continue
        double meanB = sumB / wB
        double meanF = (sumAll - sumB) / wF
        double between = wB * wF * (meanB - meanF) * (meanB - meanF)
        if (between > best) { best = between; bestT = t }
    }
    return bestT
}

// Mirrors otsu3_from_hist256 in DAPI_unet_utils.py. Three classes, because a
// spleen DAPI histogram has three modes (glass / red pulp / white pulp) and a
// 2-class split lands inside the tissue. Returns [t1, t2]; t1 is the
// glass/tissue threshold. Strict '>' so the FIRST maximum wins.
static int[] otsu3FromHist(long[] hist) {
    double total = 0.0d
    double[] W = new double[257]
    double[] S = new double[257]
    for (int i = 0; i < 256; i++) {
        total += hist[i]
        W[i + 1] = W[i] + hist[i]
        S[i + 1] = S[i] + ((double) hist[i]) * i
    }
    if (total <= 0.0d) return [0, 0] as int[]
    double best = -1.0d
    int bt1 = 0
    int bt2 = 0
    for (int t1 = 0; t1 < 255; t1++) {
        double w0 = W[t1 + 1]
        double s0 = S[t1 + 1]
        if (w0 <= 0.0d) continue
        for (int t2 = t1 + 1; t2 < 255; t2++) {
            double w1 = W[t2 + 1] - W[t1 + 1]
            if (w1 <= 0.0d) continue
            double w2 = total - W[t2 + 1]
            if (w2 <= 0.0d) continue
            double s1 = S[t2 + 1] - S[t1 + 1]
            double s2 = S[256] - S[t2 + 1]
            double v = s0 * s0 / w0 + s1 * s1 / w1 + s2 * s2 / w2
            if (v > best) { best = v; bt1 = t1; bt2 = t2 }
        }
    }
    if (best < 0.0d) {
        int t = otsuFromHist(hist)
        return [t, t] as int[]
    }
    return [bt1, bt2] as int[]
}

// Mirrors percentile_from_hist256: first bin whose cumulative count reaches the
// target, returned as that bin's LOW EDGE so the value is an exact integer in
// both languages.
static int percentileFromHist(long[] hist, double pct) {
    long total = 0L
    for (int i = 0; i < 256; i++) total += hist[i]
    if (total <= 0L) return 0
    double target = pct / 100.0d * total
    long cum = 0L
    for (int b = 0; b < 256; b++) {
        cum += hist[b]
        if (cum >= target) return b << 8
    }
    return 255 << 8
}
