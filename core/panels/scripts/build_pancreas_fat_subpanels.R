#!/usr/bin/env Rscript
#
# Build self-contained "fat" pancreas subpanel CSVs that mirror the thymus
# tissue layout, so the files are accurate and consistent when read directly
# (e.g. by the Xenium_Analysis phenotyping repo) — not only after the app's
# load-time joins.
#
# The source gene lists are "thin": the shared pool (data/subpanels_shared/) and
# the four pancreas-specific files carry only gene,category,rationale,
# single_gene_panel_note; biology + detection are otherwise merged at runtime by
# R/load_panels.R, so an external reader of those raw files sees none of them.
# This script reads the thin gene lists and writes self-contained "fat" CSVs into
# data/tissues/pancreas/subpanels/, mirroring the thymus per-category schemas:
#
#   cell-type (01-20):  gene,category,rationale,exclude_recommended,
#                       in_5K_reference,<5 biology>,detection_pct_<r1>,
#                       detection_pct_<r2>,det_max,n_cells_<r1>,n_cells_<r2>,
#                       log2_detection_ratio
#   pathway  (21-49):   ...same minus in_5K_reference, plus
#                       pathway_source,pathway_term,was_uncategorized
#   residual (99*):     ...as cell-type minus det_max, plus
#                       HALLMARK,REACTOME,KEGG,n_pathways
#
# Pancreas keeps its own curated pathway taxonomy (glycolysis, OXPHOS, ...);
# we align file *structure*, not panel membership. pathway_source is "curated"
# and pathway_term is the panel slug (the filename already encodes identity).
#
# Re-runnable / idempotent: gene lists are read from the shared pool + the four
# pancreas-specific files; everything else is recomputed from the reference 5K
# table, the pancreas audit table, and the reannotated residual. Run from the
# repository root:  Rscript scripts/build_pancreas_fat_subpanels.R
#
# The pancreas manifest already sets use_shared_subpanels: false so the app reads
# these tissue-local files (one canonical set); a new tissue adopting this layout
# must set the same flag.

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

root         <- normalizePath(".", winslash = "/")
shared_dir   <- file.path(root, "data", "subpanels_shared")
panc_dir     <- file.path(root, "data", "tissues", "pancreas", "subpanels")
ref5k_path   <- file.path(root, "data", "reference_5k", "xenium5k_genes.csv")
audit_path   <- file.path(root, "data", "tissues", "pancreas", "audit",
                          "xenium5k_audit.csv")
summary_path <- file.path(root, "data", "tissues", "pancreas",
                          "subpanel_summary.csv")
reannot_path <- file.path(shared_dir, "99_uncategorized_REANNOTATED.csv")

# Pancreas reference runs that produced the detection_pct_<run> columns.
runs     <- c("0041323", "0041326")
det_cols <- paste0("detection_pct_", runs)
nc_cols  <- paste0("n_cells_", runs)
bio_cols <- c("gene_id", "full_name", "location", "cell_type", "cellchat_pathway")

# Pancreas-specific subpanel stems (live under tissues/pancreas/subpanels/);
# everything else in the 01-49 range is read from the shared pool.
pancreas_specific <- c("01_pancreas_endocrine", "02_pancreas_exocrine",
                       "13_fibroblast_stellate", "38_insulin_biology_T1D_extended")

# Residual stems (pancreas-specific buckets, currently parked in shared/).
residual_stems <- c("99_uncategorized_REANNOTATED",
                    "99c_residual_after_reannotation",
                    "99d_truly_unannotated")

stopifnot(dir.exists(shared_dir), dir.exists(panc_dir),
          file.exists(ref5k_path), file.exists(audit_path))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

read_csv  <- function(p) data.table::fread(p, na.strings = c("", "NA"),
                                           data.table = FALSE, showProgress = FALSE)
write_csv <- function(df, p) data.table::fwrite(df, p, na = "")

# Match the upstream audit's capitalised booleans ("True"/"False").
tf <- function(x) ifelse(is.na(x), NA_character_, ifelse(x, "True", "False"))

col_or_na <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NA_character_

# Refuse to silently drop curator notes if any source actually carries them.
guard_note <- function(df, src) {
  v <- col_or_na(df, "single_gene_panel_note")
  if (any(!is.na(v) & nzchar(v))) {
    stop("single_gene_panel_note carries data in '", src,
         "'; the thymus schema has no such column. Reconcile before dropping.")
  }
}

# Count pathway terms across the three ";"-separated library columns.
count_terms <- function(x) vapply(x, function(s) {
  if (is.na(s) || !nzchar(s)) return(0L)
  sum(nzchar(trimws(strsplit(s, ";", fixed = TRUE)[[1]])))
}, integer(1))

# ---------------------------------------------------------------------------
# Reference tables
# ---------------------------------------------------------------------------

ref5k   <- read_csv(ref5k_path)     # gene + 5 biology cols (constant, 4992 genes)
audit   <- read_csv(audit_path)     # gene, exclude_recommended, detection/n_cells, log2
reannot <- read_csv(reannot_path)   # gene, HALLMARK, KEGG, REACTOME (2409 residual genes)
residual_set <- unique(reannot$gene)

# Per-gene biology + detection block, aligned to `genes`.
annotate <- function(genes) {
  bi <- ref5k[match(genes, ref5k$gene), bio_cols, drop = FALSE]
  au <- audit[match(genes, audit$gene), , drop = FALSE]
  d1 <- suppressWarnings(as.numeric(au[[det_cols[1]]]))
  d2 <- suppressWarnings(as.numeric(au[[det_cols[2]]]))
  dmax <- suppressWarnings(pmax(d1, d2, na.rm = TRUE))
  dmax[is.infinite(dmax)] <- NA_real_
  list(bi = bi, exclude = au[["exclude_recommended"]],
       in5k = genes %in% ref5k$gene,
       d1 = d1, d2 = d2, dmax = dmax,
       n1 = au[[nc_cols[1]]], n2 = au[[nc_cols[2]]],
       log2 = au[["log2_detection_ratio"]])
}

# Shared leading columns (gene, curation, exclude). `with_in5k` toggles the
# in_5K_reference column (present on cell-type/residual, absent on pathway).
lead_block <- function(df, a, with_in5k) {
  base <- data.frame(
    gene                = df$gene,
    category            = col_or_na(df, "category"),
    rationale           = col_or_na(df, "rationale"),
    exclude_recommended = a$exclude,
    check.names = FALSE, stringsAsFactors = FALSE)
  if (with_in5k) base[["in_5K_reference"]] <- tf(a$in5k)
  base <- cbind(base, a$bi)
  base[[det_cols[1]]] <- a$d1
  base[[det_cols[2]]] <- a$d2
  base
}

build_celltype <- function(df) {
  a <- annotate(df$gene); out <- lead_block(df, a, with_in5k = TRUE)
  out[["det_max"]]              <- a$dmax
  out[[nc_cols[1]]]             <- a$n1
  out[[nc_cols[2]]]             <- a$n2
  out[["log2_detection_ratio"]] <- a$log2
  out
}

build_pathway <- function(df, term) {
  a <- annotate(df$gene); out <- lead_block(df, a, with_in5k = FALSE)
  out[["det_max"]]              <- a$dmax
  out[[nc_cols[1]]]             <- a$n1
  out[[nc_cols[2]]]             <- a$n2
  out[["log2_detection_ratio"]] <- a$log2
  out[["pathway_source"]]       <- "curated"
  out[["pathway_term"]]         <- term
  out[["was_uncategorized"]]    <- tf(df$gene %in% residual_set)
  out
}

build_residual <- function(df) {
  a <- annotate(df$gene); out <- lead_block(df, a, with_in5k = TRUE)
  out[[nc_cols[1]]]             <- a$n1
  out[[nc_cols[2]]]             <- a$n2
  out[["log2_detection_ratio"]] <- a$log2
  H <- col_or_na(df, "HALLMARK"); R <- col_or_na(df, "REACTOME"); K <- col_or_na(df, "KEGG")
  out[["HALLMARK"]]   <- H
  out[["REACTOME"]]   <- R
  out[["KEGG"]]       <- K
  out[["n_pathways"]] <- count_terms(H) + count_terms(R) + count_terms(K)
  out
}

# ---------------------------------------------------------------------------
# Drive: numbered subpanels (01-49), then residual buckets
# ---------------------------------------------------------------------------

shared_files    <- list.files(shared_dir, pattern = "\\.csv$", full.names = TRUE)
shared_numbered <- shared_files[grepl("^[0-9]{2}_", basename(shared_files)) &
                                !grepl("^99",      basename(shared_files))]
numbered_paths  <- c(shared_numbered,
                     file.path(panc_dir, paste0(pancreas_specific, ".csv")))

# Read every source once, run the note guard, and classify it. Defer all writes
# until after the 5K-only invariant is validated below.
sources <- list()
for (p in numbered_paths) {
  stem <- sub("\\.csv$", "", basename(p))
  df   <- read_csv(p); guard_note(df, stem)
  num  <- as.integer(sub("_.*", "", stem))
  sources[[stem]] <- list(df = df,
                          kind = if (num <= 20L) "celltype" else "pathway",
                          term = sub("^[0-9]+_", "", stem))
}
for (stem in residual_stems) {
  df <- read_csv(file.path(shared_dir, paste0(stem, ".csv")))
  guard_note(df, stem)
  sources[[stem]] <- list(df = df, kind = "residual", term = NA_character_)
}

# Fail fast: 5K-only membership is a hard invariant for these subpanels (the PR
# and downstream consumers treat it as a guarantee). Refuse to write anything if
# a source carries a gene outside the 5K reference.
offenders <- do.call(rbind, lapply(names(sources), function(stem) {
  bad <- setdiff(sources[[stem]]$df$gene, ref5k$gene)
  if (length(bad)) {
    data.frame(subpanel = stem, n_non5k = length(bad),
               examples = paste(utils::head(bad, 5L), collapse = ";"),
               stringsAsFactors = FALSE)
  }
}))
if (!is.null(offenders) && nrow(offenders)) {
  cat("\n[FAIL] subpanel genes NOT in the 5K reference (5K-only invariant):\n")
  print(offenders, row.names = FALSE)
  stop(sum(offenders$n_non5k), " subpanel gene(s) are not in the 5K reference; ",
       "refusing to write fat subpanels. Fix the source gene lists first.")
}

# Write phase — build each subpanel to its per-category schema and emit the CSV.
report <- list()
for (stem in names(sources)) {
  s   <- sources[[stem]]
  out <- switch(s$kind,
                celltype = build_celltype(s$df),
                pathway  = build_pathway(s$df, s$term),
                residual = build_residual(s$df))
  write_csv(out, file.path(panc_dir, paste0(stem, ".csv")))
  report[[stem]] <- data.frame(subpanel = stem, n_genes = nrow(out),
                               schema = s$kind, stringsAsFactors = FALSE)
}
rep <- do.call(rbind, report); rownames(rep) <- NULL

# ---------------------------------------------------------------------------
# Verification report
# ---------------------------------------------------------------------------

cat(sprintf("\nWrote %d subpanel files to %s\n", nrow(rep),
            sub(paste0(root, "/"), "", panc_dir)))
cat(sprintf("  by schema: %s\n",
            paste(sprintf("%s=%d", names(table(rep$schema)), table(rep$schema)),
                  collapse = ", ")))
cat("[OK] every subpanel gene is in the 5K reference (validated before writing)\n")

# Consistency: gene counts should match subpanel_summary.csv (numbered panels).
if (file.exists(summary_path)) {
  sm <- read_csv(summary_path)
  m  <- merge(rep[rep$schema != "residual", c("subpanel", "n_genes")],
              sm[, c("subpanel", "n_genes")],
              by = "subpanel", suffixes = c("_file", "_summary"))
  bad <- m[m$n_genes_file != m$n_genes_summary, ]
  if (nrow(bad)) {
    cat("\n[WARN] gene-count mismatches vs subpanel_summary.csv:\n")
    print(bad, row.names = FALSE)
  } else {
    cat(sprintf("[OK] all %d numbered subpanel gene counts match subpanel_summary.csv\n",
                nrow(m)))
  }
}

cat("\nDone. The pancreas manifest already sets `use_shared_subpanels: false`",
    "so the app reads these files; verify it if adapting this for a new tissue.\n")
