#' Fetch 10x Genomics' pre-designed Xenium panel gene lists into
#' `data/reference_panels/`. Runs once at setup; the app loads whatever
#' is present at startup (see [load_reference_panels()] in
#' R/load_panels.R), so re-running the script just refreshes the files.
#'
#' Each panel is normalized to a tidy CSV keyed on `gene`. Extra columns
#' from the source (Ensembl_ID, Num_Probesets, Codewords, Annotation,
#' Tissues, Cell types, etc.) are preserved as-is for downstream lookup.
#'
#' Sources resolved 2026-05-22 from the public 10x panel resources page
#' (`https://10xgen.com/...` shortlinks redirect to `cdn.10xgenomics.com`
#' CSVs). If a panel URL 404s after a 10x revision, fetch the shortlink
#' chain manually and update PANELS below.

PANELS <- list(
  list(id = "hPrime5K",
       display_name  = "Human Prime 5K Pan-Tissue & Pathways",
       species = "Human", tissue = "Pan-Tissue", expected_genes = 5001,
       url = "https://10xgen.com/prime-5k-human"),
  list(id = "mPrime5K",
       display_name  = "Mouse Prime 5K Pan-Tissue & Pathways",
       species = "Mouse", tissue = "Pan-Tissue", expected_genes = 5006,
       url = "https://10xgen.com/prime-5k-mouse"),
  list(id = "hBrain",  display_name = "Human Brain",
       species = "Human", tissue = "Brain",  expected_genes = 266,
       url = "https://10xgen.com/v1-human-brain"),
  list(id = "hBreast", display_name = "Human Breast",
       species = "Human", tissue = "Breast", expected_genes = 280,
       url = "https://10xgen.com/v1-human-breast"),
  list(id = "hColon",  display_name = "Human Colon",
       species = "Human", tissue = "Colon",  expected_genes = 322,
       url = "https://10xgen.com/v1-human-colon"),
  list(id = "hIO",     display_name = "Human Immuno-Oncology",
       species = "Human", tissue = "Multi-Tissue", expected_genes = 380,
       url = "https://10xgen.com/v1-human-immunooncology"),
  list(id = "hLung",   display_name = "Human Lung",
       species = "Human", tissue = "Lung",   expected_genes = 289,
       url = "https://10xgen.com/v1-human-lung"),
  list(id = "hMulti",  display_name = "Human Multi-Tissue & Cancer",
       species = "Human", tissue = "Multi-Tissue", expected_genes = 377,
       url = "https://10xgen.com/v1-human-multi"),
  list(id = "hSkin",   display_name = "Human Skin",
       species = "Human", tissue = "Skin",   expected_genes = 260,
       url = "https://10xgen.com/v1-human-skin"),
  list(id = "mBrain",  display_name = "Mouse Brain",
       species = "Mouse", tissue = "Brain",  expected_genes = 247,
       url = "https://10xgen.com/v1-mouse-brain"),
  list(id = "mMulti",  display_name = "Mouse Tissue Atlas",
       species = "Mouse", tissue = "Multi-Tissue", expected_genes = 379,
       url = "https://10xgen.com/v1-mouse-multi")
)

OUT_DIR <- "data/reference_panels"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Normalize whatever 10x ships into a `gene`-keyed data.frame. Different
# panels use `Genes`, `Gene`, or `gene_name`; some misspell "Ensembl"
# as "Ensemble". We rename gene -> `gene` and keep the rest as-is.
normalize_panel <- function(df) {
  if (!ncol(df)) stop("empty CSV")
  gene_col <- intersect(c("Genes", "Gene", "gene_name", "gene"), names(df))
  if (!length(gene_col)) {
    # fall back to the first column if nothing matches
    gene_col <- names(df)[1]
  }
  names(df)[names(df) == gene_col[1]] <- "gene"
  df$gene <- trimws(as.character(df$gene))
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  df <- df[!duplicated(df$gene), , drop = FALSE]
  df
}

manifest_rows <- list()
for (p in PANELS) {
  message("Fetching ", p$id, " <- ", p$url)
  tmp <- tempfile(fileext = ".csv")
  status <- tryCatch(
    utils::download.file(p$url, destfile = tmp, quiet = TRUE, mode = "wb"),
    error = function(e) -1L
  )
  if (!identical(status, 0L)) {
    warning(p$id, ": download failed (status ", status, "); skipping")
    next
  }
  df <- tryCatch(
    data.table::fread(tmp, na.strings = c("", "NA"), data.table = FALSE,
                      showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) {
    warning(p$id, ": parse failed or empty; skipping")
    next
  }
  df <- normalize_panel(df)
  out_path <- file.path(OUT_DIR, paste0(p$id, ".csv"))
  data.table::fwrite(df, out_path)
  n <- nrow(df)
  if (!is.null(p$expected_genes) && abs(n - p$expected_genes) > 5) {
    warning(p$id, ": expected ~", p$expected_genes, " genes, got ", n)
  }
  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    panel_id     = p$id,
    display_name = p$display_name,
    species      = p$species,
    tissue       = p$tissue,
    n_genes      = n,
    file         = basename(out_path),
    source_url   = p$url,
    stringsAsFactors = FALSE
  )
}

if (length(manifest_rows)) {
  manifest <- do.call(rbind, manifest_rows)
  manifest <- manifest[order(manifest$species, manifest$panel_id), ]
  data.table::fwrite(manifest, file.path(OUT_DIR, "manifest.csv"))
  message("Wrote ", nrow(manifest), " panels to ", OUT_DIR)
  print(manifest[, c("panel_id", "display_name", "n_genes")])
} else {
  warning("no panels were fetched")
}
