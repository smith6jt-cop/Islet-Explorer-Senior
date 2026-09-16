# Coordinate frames as IDENTITIES, not scales.
#
# Four pixel scales and three coordinate systems are in play across this project:
# Xenium 0.2125 um/px, post-Xenium IF 0.25 (and 0.5 for donors 115 and 6591),
# PhenoCycler 0.5077663810243286, and the DOMAIN's rounded 0.5078. Two tables can
# agree on um/px and still describe different places.
#
# The failure this file exists to prevent is silent. A Xenium cell at x_um = 5000 and a
# PhenoCycler cell at x_um = 5000 are not near each other: the measured section-to-section
# rotations across the cohort span -168.9 to +88.6 degrees, and 12 of 20 donors exceed 20.
# Plotted together they produce a picture that looks like data. So a frame is an identity
# -- (kind, donor, roi) -- and mixing two is an error, not a judgement call.
#
# Sections mounted on one slide are the case that motivates keying on roi as well as donor:
# a pancreas section and a pLN section cut from DIFFERENT BLOCKS share a raster and share no
# anatomy. `phenocycler_um:112:panc` and `phenocycler_um:112:pln_1` are different frames even
# though every number in them came off one image at one scale.

#' Frame kinds this registry recognises
#'
#' `hex_qr`, `signed_distance_um` and `none` are SCALE-FREE: a hex bin label and a signed
#' distance are not positions in anybody's coordinate system, so they carry no donor or roi.
#' The other three are spatial and require both.
#' @export
IE_FRAME_KINDS <- c("xenium_um", "if_um", "phenocycler_um",
                    "hex_qr", "signed_distance_um", "none")

#' @export
IE_FRAME_SCALE_FREE <- c("hex_qr", "signed_distance_um", "none")

#' Build a frame identity
#'
#' @param kind one of [IE_FRAME_KINDS]
#' @param donor donor id; required for spatial kinds, refused for scale-free ones
#' @param roi roi token (`panc`, `pln_1`, ...); same rule as `donor`
#' @return a single string, e.g. `"xenium_um:6436:panc"`
#' @export
ie_frame <- function(kind, donor = NULL, roi = NULL) {
  if (length(kind) != 1L || is.na(kind) || !nzchar(kind)) {
    stop("ie_frame(): kind must be one non-empty string", call. = FALSE)
  }
  if (!kind %in% IE_FRAME_KINDS) {
    stop("ie_frame(): unknown kind '", kind, "'. Known: ",
         paste(IE_FRAME_KINDS, collapse = ", "), call. = FALSE)
  }
  scale_free <- kind %in% IE_FRAME_SCALE_FREE
  if (scale_free) {
    if (!is.null(donor) || !is.null(roi)) {
      stop("ie_frame(): '", kind, "' is scale-free and takes no donor or roi. ",
           "A hex bin label and a signed distance are not positions.", call. = FALSE)
    }
    return(kind)
  }
  bad <- function(x) is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))
  if (bad(donor) || bad(roi)) {
    stop("ie_frame(): '", kind, "' is a spatial frame and requires both donor and roi. ",
         "Two sections from different blocks can share a scale and share no anatomy, ",
         "so the scale alone does not identify a frame.", call. = FALSE)
  }
  paste(kind, as.character(donor), as.character(roi), sep = ":")
}

#' Take a frame identity apart
#' @export
ie_frame_parse <- function(frame) {
  if (length(frame) != 1L || is.na(frame) || !nzchar(frame)) {
    stop("ie_frame_parse(): frame must be one non-empty string", call. = FALSE)
  }
  parts <- strsplit(frame, ":", fixed = TRUE)[[1]]
  kind <- parts[[1]]
  if (!kind %in% IE_FRAME_KINDS) {
    stop("ie_frame_parse(): unknown kind '", kind, "' in frame '", frame, "'", call. = FALSE)
  }
  if (kind %in% IE_FRAME_SCALE_FREE) {
    if (length(parts) != 1L) {
      stop("ie_frame_parse(): '", kind, "' is scale-free but carries ", length(parts) - 1L,
           " extra field(s): '", frame, "'", call. = FALSE)
    }
    return(list(kind = kind, donor = NA_character_, roi = NA_character_))
  }
  if (length(parts) != 3L) {
    stop("ie_frame_parse(): spatial frame must be 'kind:donor:roi'; got '", frame, "'",
         call. = FALSE)
  }
  list(kind = kind, donor = parts[[2]], roi = parts[[3]])
}

#' Attach a frame identity to a data frame
#' @export
ie_attach_frame <- function(df, frame) {
  ie_frame_parse(frame)  # refuses an unconstructible value before it can be attached
  attr(df, "ie_frame") <- frame
  df
}

#' The frame a data frame carries
#'
#' Errors when the attribute is absent. A frameless set of coordinates is a defect, not a
#' default -- there is no sensible frame to assume, and assuming one is how two tissues end
#' up on one pair of axes.
#' @export
ie_frame_of <- function(df, context = NULL) {
  frame <- attr(df, "ie_frame", exact = TRUE)
  if (is.null(frame)) {
    stop("ie_frame_of(): no ie_frame attribute",
         if (!is.null(context)) paste0(" on ", context) else "",
         ". Coordinates without a declared frame are a defect; attach one with ",
         "ie_attach_frame() at the point the data is read.", call. = FALSE)
  }
  frame
}

#' Refuse to proceed unless two data frames are in the same frame
#'
#' Call this INSIDE the plot-building reactive, not at the module edge: a later refactor
#' routes around an edge check, and the thing being protected is the plot.
#' @export
ie_assert_same_frame <- function(a, b, context = "this plot") {
  fa <- ie_frame_of(a, paste0(context, " (first argument)"))
  fb <- ie_frame_of(b, paste0(context, " (second argument)"))
  if (!identical(fa, fb)) {
    pa <- ie_frame_parse(fa); pb <- ie_frame_parse(fb)
    why <- if (!identical(pa$kind, pb$kind)) {
      "different modalities"
    } else if (!identical(pa$donor, pb$donor)) {
      "different donors"
    } else {
      paste0("different tissue sections (", pa$roi, " vs ", pb$roi,
             ") -- sections mounted on one slide are cut from different blocks: ",
             "they share a raster and share no anatomy")
    }
    stop("ie_assert_same_frame(): refusing to combine '", fa, "' with '", fb,
         "' in ", context, " -- ", why, ".", call. = FALSE)
  }
  invisible(fa)
}

#' The microns-per-pixel a producer recorded for a section
#'
#' There is deliberately no um/px constant in this registry. The value is a property of an
#' image, it varies within one modality (post-Xenium IF is 0.25 for most donors and 0.5 for
#' 115 and 6591), and the DOMAIN's 0.5078 is a rounded display constant, not a calibration.
#' Read what the producer wrote or stop.
#' @export
ie_um_per_px <- function(section_row, field = "um_per_px") {
  if (is.null(section_row) || !length(section_row) || !field %in% names(section_row)) {
    stop("ie_um_per_px(): no '", field, "' recorded for this section. The scale is a ",
         "property of the image and must come from its producer -- never from DOMAIN, ",
         "whose pixel_size_um is a rounded display constant.", call. = FALSE)
  }
  v <- suppressWarnings(as.numeric(section_row[[field]][1]))
  if (!is.finite(v) || v <= 0) {
    stop("ie_um_per_px(): recorded '", field, "' is not a positive number: ",
         format(section_row[[field]][1]), call. = FALSE)
  }
  v
}
