
createAnnotationsFromPixelClassifier("Tissue", 10000.0, 10000.0)
selectObjectsByClassification("Tissue");
qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath("/home/smith6jt/QuPath/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1")
    .device("gpu0")
    .inputChannels([ColorTransforms.createChannelExtractor("DAPI"), ColorTransforms.createChannelExtractor("CD45RO"), ColorTransforms.createChannelExtractor("CD38"), ColorTransforms.createChannelExtractor("CD20"), ColorTransforms.createChannelExtractor("CD4"), ColorTransforms.createChannelExtractor("CD44"), ColorTransforms.createChannelExtractor("CD163"), ColorTransforms.createChannelExtractor("HLA-DR"), ColorTransforms.createChannelExtractor("CD68"), ColorTransforms.createChannelExtractor("CD8"), ColorTransforms.createChannelExtractor("CD66"), ColorTransforms.createChannelExtractor("CD141"), ColorTransforms.createChannelExtractor("CD57"), ColorTransforms.createChannelExtractor("Vimentin"), ColorTransforms.createChannelExtractor("CD3e"), ColorTransforms.createChannelExtractor("HLA-A"), ColorTransforms.createChannelExtractor("CD45")])
    .outputChannels()
    .tileDims(512)
    .interTilePadding(32)
    .nThreads(4)
    .makeMeasurements(false)
    .randomColors(false)
    .outputType("Default")
    .build()
    .detectObjects()

//
//runObjectClassifier("SPLEEN_latest")
//
//createAnnotationsFromPixelClassifier("SpleenRegions", 100.0, 1000.0)
