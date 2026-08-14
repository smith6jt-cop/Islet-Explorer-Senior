#!/usr/bin/env Rscript
#
# One-shot migration: data/panel_audit/  ->  three-layer tissue-aware layout.
#
#   data/
#   ├── reference_5k/xenium5k_genes.csv          (constant: 10x annotations)
#   ├── subpanels_shared/*.csv                   (biology-only subpanel defs)
#   └── tissues/pancreas/                        (tissue-specific)
#       ├── manifest.yml
#       ├── subpanel_summary.csv
#       ├── subpanels/{4 pancreas-specific CSVs}
#       ├── audit/xenium5k_audit.csv             (detection / exclude / n_cells)
#       ├── audit/xenium5k_excluded.csv
#       ├── audit/hIO_*.csv
#       └── custom_panel.csv                     (was custom_T1D_GWAS_panel.csv)
#
# Plus a stub `data/tissues/thymus/` so the tissue picker has a second
# tissue at startup.
#
# Re-runnable: writes into a `data/` tree alongside the existing
# `data/panel_audit/`. Verifies counts at the end. Does NOT delete the
# legacy directory — prints the rm command for the operator to run.

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

root <- normalizePath(".", winslash = "/")
src  <- file.path(root, "data", "panel_audit")
out_ref5k     <- file.path(root, "data", "reference_5k")
out_shared    <- file.path(root, "data", "subpanels_shared")
out_pancreas  <- file.path(root, "data", "tissues", "pancreas")
out_thymus    <- file.path(root, "data", "tissues", "thymus")

if (!dir.exists(src)) {
  stop("Source directory not found: ", src,
       "\nRun this script from the repository root.")
}

# Pancreas-specific subpanels — everything else is shared biology.
# Stellate cells are a pancreas-context feature; if thymus needs a
# fibroblast subpanel later it will get its own file.
pancreas_specific_subpanels <- c(
  "01_pancreas_endocrine",
  "02_pancreas_exocrine",
  "13_fibroblast_stellate",
  "38_insulin_biology_T1D_extended"
)

# Reference run IDs that produced the detection_pct_<run> columns in
# the current audit. Pancreas tissue's two Xenium runs.
pancreas_reference_runs <- c("0041323", "0041326")

# Biology / annotation columns that live in the constant 5K reference.
ref5k_cols <- c("gene", "gene_id", "full_name", "location",
                "cell_type", "cellchat_pathway")

# Columns that are tissue-audit data: move to tissues/<t>/audit/.
audit_cols_keep <- c("gene", "exclude_recommended",
                     paste0("detection_pct_", pancreas_reference_runs),
                     paste0("n_cells_",       pancreas_reference_runs),
                     "log2_detection_ratio")   # renamed from log2_..._326_over_323

# Columns we strip from subpanel CSVs (and from custom_panel.csv) on the
# way through. These all get re-joined at load time from the audit table
# (detection / exclude) or the reference 5K (biology).
strip_from_subpanels <- c(
  ref5k_cols[-1],                                # gene_id, full_name, location, cell_type, cellchat_pathway
  "exclude_recommended",
  paste0("detection_pct_", pancreas_reference_runs),
  paste0("n_cells_",       pancreas_reference_runs),
  "log2_detection_ratio_326_over_323",
  "log2_detection_ratio"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  invisible(p)
}

read_csv <- function(p) {
  data.table::fread(p, na.strings = c("", "NA"), data.table = FALSE,
                    showProgress = FALSE)
}

write_csv <- function(df, p) {
  data.table::fwrite(df, p, na = "")
}

drop_cols <- function(df, cols) {
  keep <- setdiff(names(df), cols)
  df[, keep, drop = FALSE]
}

# Rename log2_detection_ratio_<rb>_over_<ra> -> log2_detection_ratio
# (the audit always carries a single ratio between the first two runs).
canonicalize_log2_ratio <- function(df) {
  hit <- grep("^log2_detection_ratio_.+_over_.+$", names(df), value = TRUE)
  if (length(hit)) names(df)[match(hit[1], names(df))] <- "log2_detection_ratio"
  df
}

# ---------------------------------------------------------------------------
# 1. Build reference_5k/xenium5k_genes.csv (tissue-agnostic)
# ---------------------------------------------------------------------------

ensure_dir(out_ref5k)
xen5k_src <- read_csv(file.path(src, "xenium5k_in_audit.csv"))
have_ref <- intersect(ref5k_cols, names(xen5k_src))
ref5k    <- unique(xen5k_src[, have_ref, drop = FALSE])
write_csv(ref5k, file.path(out_ref5k, "xenium5k_genes.csv"))
cat(sprintf("[reference_5k] wrote xenium5k_genes.csv  (%d genes, %d cols)\n",
            nrow(ref5k), ncol(ref5k)))

# ---------------------------------------------------------------------------
# 2. Build tissues/pancreas/audit/xenium5k_audit.csv (tissue-specific)
#
# Audit gene set = union of every gene that has detection_pct in any
# source CSV (xenium5k_in_audit + every subpanel_NN_*.csv + custom panel).
# That guards against subpanels that reference genes not in the 4,992
# 5K-shared pool.
# ---------------------------------------------------------------------------

ensure_dir(file.path(out_pancreas, "audit"))

# Pull (gene + audit cols) from a source frame, normalising the ratio col.
pull_audit <- function(df) {
  df <- canonicalize_log2_ratio(df)
  have <- intersect(audit_cols_keep, names(df))
  if (!"gene" %in% have) return(NULL)
  df[, have, drop = FALSE]
}

audit_frames <- list()
audit_frames[["xenium5k_in_audit"]] <- pull_audit(xen5k_src)

subpanel_files <- list.files(src, pattern = "^subpanel_[0-9].*\\.csv$",
                             full.names = TRUE)
for (f in subpanel_files) {
  audit_frames[[basename(f)]] <- pull_audit(read_csv(f))
}

cust_src <- read_csv(file.path(src, "custom_T1D_GWAS_panel.csv"))
audit_frames[["custom_T1D_GWAS_panel"]] <- pull_audit(cust_src)

# Stack and de-duplicate by gene; prefer the first non-NA value per column
# (xenium5k_in_audit is listed first so its values win on conflict).
audit_frames <- Filter(function(d) !is.null(d) && nrow(d) > 0L, audit_frames)
all_cols <- unique(unlist(lapply(audit_frames, names)))
audit_frames <- lapply(audit_frames, function(d) {
  for (col in setdiff(all_cols, names(d))) d[[col]] <- NA
  d[, all_cols, drop = FALSE]
})
audit_stack <- do.call(rbind, audit_frames)
audit_stack <- audit_stack[!is.na(audit_stack$gene) &
                            nzchar(audit_stack$gene), , drop = FALSE]

# Collapse per gene: first non-NA wins for each column.
audit_dt <- data.table::as.data.table(audit_stack)
collapse_first_nonNA <- function(x) {
  i <- which(!is.na(x))
  if (length(i)) x[i[1]] else x[1]
}
audit_master <- as.data.frame(
  audit_dt[, lapply(.SD, collapse_first_nonNA), by = "gene"]
)
# Order columns as audit_cols_keep, then any extras.
ord <- c(intersect(audit_cols_keep, names(audit_master)),
         setdiff(names(audit_master), audit_cols_keep))
audit_master <- audit_master[, ord, drop = FALSE]
write_csv(audit_master, file.path(out_pancreas, "audit", "xenium5k_audit.csv"))
cat(sprintf("[pancreas/audit]  wrote xenium5k_audit.csv  (%d genes)\n",
            nrow(audit_master)))

# ---------------------------------------------------------------------------
# 3. Split subpanels: shared/ vs tissues/pancreas/subpanels/
# ---------------------------------------------------------------------------

ensure_dir(out_shared)
ensure_dir(file.path(out_pancreas, "subpanels"))

n_shared <- 0L; n_tissue <- 0L
for (f in subpanel_files) {
  key <- sub("^subpanel_", "", sub("\\.csv$", "", basename(f)))
  df  <- read_csv(f)
  df  <- canonicalize_log2_ratio(df)
  df  <- drop_cols(df, strip_from_subpanels)

  dest <- if (key %in% pancreas_specific_subpanels) {
    n_tissue <- n_tissue + 1L
    file.path(out_pancreas, "subpanels", paste0(key, ".csv"))
  } else {
    n_shared <- n_shared + 1L
    file.path(out_shared, paste0(key, ".csv"))
  }
  write_csv(df, dest)
}
cat(sprintf("[subpanels]       %d shared, %d pancreas-specific\n",
            n_shared, n_tissue))

# ---------------------------------------------------------------------------
# 4. Custom panel
# ---------------------------------------------------------------------------

cust_out <- drop_cols(canonicalize_log2_ratio(cust_src), strip_from_subpanels)
write_csv(cust_out, file.path(out_pancreas, "custom_panel.csv"))
cat(sprintf("[pancreas]        wrote custom_panel.csv  (%d genes, %d cols)\n",
            nrow(cust_out), ncol(cust_out)))

# ---------------------------------------------------------------------------
# 5. Excluded list + hIO_* supplementary tables
# ---------------------------------------------------------------------------

excluded <- read_csv(file.path(src, "xenium5k_already_excluded.csv"))
excluded <- canonicalize_log2_ratio(excluded)
write_csv(excluded, file.path(out_pancreas, "audit", "xenium5k_excluded.csv"))
cat(sprintf("[pancreas/audit]  wrote xenium5k_excluded.csv  (%d genes)\n",
            nrow(excluded)))

hIO_files <- list.files(src, pattern = "^hIO_.*\\.csv$", full.names = TRUE)
for (f in hIO_files) {
  df <- read_csv(f)
  write_csv(df, file.path(out_pancreas, "audit", basename(f)))
}
cat(sprintf("[pancreas/audit]  wrote %d hIO_*.csv files\n", length(hIO_files)))

# ---------------------------------------------------------------------------
# 6. Per-tissue subpanel summary (carries n_genes/pct from pancreas data)
# ---------------------------------------------------------------------------

summary_src <- read_csv(file.path(src, "subpanel_summary_v2.csv"))
write_csv(summary_src, file.path(out_pancreas, "subpanel_summary.csv"))
cat(sprintf("[pancreas]        wrote subpanel_summary.csv  (%d rows)\n",
            nrow(summary_src)))

# ---------------------------------------------------------------------------
# 7. Manifests
# ---------------------------------------------------------------------------

pancreas_manifest <- paste0(
  'tissue_id: pancreas\n',
  'display_name: "Human Pancreas"\n',
  'reference_runs:\n',
  paste(sprintf('  - "%s"', pancreas_reference_runs), collapse = "\n"),
  "\n",
  'custom_panel:\n',
  '  file: custom_panel.csv\n',
  '  display_name: "T1D-GWAS custom 100"\n',
  'notes: |\n',
  '  Audit columns: detection_pct_<run>, n_cells_<run>,\n',
  '  log2_detection_ratio (between first two reference_runs),\n',
  '  exclude_recommended.\n'
)
writeLines(pancreas_manifest,
           con = file.path(out_pancreas, "manifest.yml"))
cat("[pancreas]        wrote manifest.yml\n")

# Thymus stub: empty subpanels/audit, no custom_panel block.
ensure_dir(file.path(out_thymus, "subpanels"))
ensure_dir(file.path(out_thymus, "audit"))
thymus_manifest <- paste0(
  'tissue_id: thymus\n',
  'display_name: "Human Thymus"\n',
  'reference_runs: []\n',
  '# custom_panel: (omit until a thymus add-on is configured)\n',
  'notes: |\n',
  '  Stub tissue. Drop tissue-specific subpanels under subpanels/,\n',
  '  a 5K audit table at audit/xenium5k_audit.csv, then add\n',
  '  reference_runs to enable the detection scatter plot.\n'
)
writeLines(thymus_manifest, con = file.path(out_thymus, "manifest.yml"))
# Empty subpanel summary so the Overview table renders without erroring.
write_csv(data.frame(subpanel = character(),
                     n_genes  = integer(),
                     description = character(),
                     pct_of_5K_kept_pool = numeric()),
          file.path(out_thymus, "subpanel_summary.csv"))
cat("[thymus]          wrote manifest.yml + empty subpanel_summary.csv\n")

# ---------------------------------------------------------------------------
# 8. Sanity-check + report
# ---------------------------------------------------------------------------

cat("\n--- verification ---\n")
cat(sprintf("reference_5k genes:           %d (expected ~4992)\n",
            nrow(ref5k)))
cat(sprintf("pancreas audit genes:         %d (expected >= 4992)\n",
            nrow(audit_master)))
cat(sprintf("pancreas subpanels:           %d (expected %d)\n",
            length(list.files(file.path(out_pancreas, "subpanels"),
                              pattern = "\\.csv$")),
            length(pancreas_specific_subpanels)))
cat(sprintf("shared subpanels:             %d\n",
            length(list.files(out_shared, pattern = "\\.csv$"))))
cat(sprintf("pancreas custom panel genes:  %d (expected 100)\n",
            nrow(cust_out)))

cat("\nMigration complete. Once you have verified the new layout, remove the\n")
cat("legacy directory:\n")
cat(sprintf("  rm -rf %s\n", file.path("data", "panel_audit")))
