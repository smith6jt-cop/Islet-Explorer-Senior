/**
 * Re-export the current image as a new OME-TIFF with freshly generated pyramid levels.
 * This fixes corruption in specific pyramid tiles by reading pixel data through
 * the image server (which decodes from intact full-res data) and writing clean compressed tiles.
 *
 * Run in QuPath: Extensions > Script Editor > paste & run
 * Output is written alongside the original with a "_fixed" suffix.
 */

import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.regions.RegionRequest

def imageData = getCurrentImageData()
def server = imageData.getServer()
def metadata = server.getMetadata()

// Build output path next to the original
def originalPath = server.getURIs()[0].toString()
// Handle both file:/ and file:/// URI formats
def originalFile = new File(new URI(originalPath))
def outputName = originalFile.name.replaceFirst('(?i)\\.(ome\\.tiff?|tiff?)$', '_fixed.ome.tif')
def outputFile = new File(originalFile.parent, outputName)

print("Source: ${originalFile.name}")
print("Output: ${outputFile.absolutePath}")
print("Image size: ${metadata.getWidth()} x ${metadata.getHeight()}")
print("Channels: ${metadata.getSizeC()}")
print("Pixel type: ${metadata.getPixelType()}")

// Configure the writer
def builder = new OMEPyramidWriter.Builder(server)
    .tileSize(512)
    .downsamples(1, 4, 8, 16, 32)   // skip downsample=2 (corrupt level)
    .parallelize()

// Use lossless compression
builder.losslessCompression()

print("Writing new OME-TIFF (this may take a while for large images)...")
long startTime = System.currentTimeMillis()

builder.build().writePyramid(outputFile.absolutePath)

long elapsed = (System.currentTimeMillis() - startTime) / 1000
print("Done in ${elapsed} seconds: ${outputFile.absolutePath}")
print("")
print("Verify the new file, then you can replace the original if needed.")