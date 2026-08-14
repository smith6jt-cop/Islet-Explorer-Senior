#' Tissue-aware panel loading.
#'
#' Three layers on disk:
#'   1. `data/reference_5k/xenium5k_genes.csv` — constant 10x Xenium Prime
#'      5K Human Pan-Tissue gene list with biology annotations
#'      (gene, gene_id, full_name, location, cell_type, cellchat_pathway).
#'   2. `data/subpanels_shared/*.csv` — biology-only subpanel definitions
#'      (gene, category, rationale, single_gene_panel_note [, HALLMARK,
#'      KEGG, REACTOME]).
#'   3. `data/tissues/<tissue>/` — tissue-specific input:
#'        manifest.yml        — tissue_id, reference_runs, optional custom_panel
#'        subpanel_summary.csv
#'        subpanels/*.csv     — tissue-specific subpanels (biology-only)
#'        audit/xenium5k_audit.csv — per-gene detection_pct_<run>,
#'                                    n_cells_<run>, log2_detection_ratio,
#'                                    exclude_recommended
#'        audit/xenium5k_excluded.csv (optional)
#'        audit/hIO_*.csv     (optional)
#'        custom_panel.csv    (optional, declared in manifest)
#'
#' At load time the composer left-joins the audit table and the 5K
#' reference biology columns onto every subpanel and the custom panel.
#' Downstream code receives a `panels` object with the same field shape
#' the pre-refactor code expected; the joined columns make subpanel
#' data.frames carry detection_pct_<run> and biology annotations inline.

# ---------------------------------------------------------------------------
# Atomic readers
# ---------------------------------------------------------------------------

.read_csv <- function(p) {
  data.table::fread(p, na.strings = c("", "NA"), data.table = FALSE,
                    showProgress = FALSE)
}

#' Read the constant 5K reference gene table.
#' @return data.frame keyed by `gene`.
load_reference_5k <- function(path = app_paths$reference_5k) {
  f <- file.path(path, "xenium5k_genes.csv")
  if (!file.exists(f)) {
    stop("reference 5K not found: ", f,
         "\nExpected the constant `xenium5k_genes.csv` to live here.")
  }
  .read_csv(f)
}

#' Normalize a reference-panel CSV to a `gene`-keyed data.frame.
#'
#' 10x panel CSVs ship with various spellings (`Genes`, `Gene`,
#' `gene_name`); user drop-ins can use anything. Renames the first
#' matching column to `gene`, trims whitespace, drops blanks/dupes.
.normalize_reference_panel_df <- function(df) {
  if (!ncol(df)) return(df)
  gene_col <- intersect(c("gene", "Gene", "Genes", "gene_name"), names(df))
  if (!length(gene_col)) gene_col <- names(df)[1]
  if (gene_col[1] != "gene") names(df)[names(df) == gene_col[1]] <- "gene"
  df$gene <- trimws(as.character(df$gene))
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  df[!duplicated(df$gene), , drop = FALSE]
}

#' Read the 10x pre-designed Xenium panels and any user-dropped panels.
#'
#' `scripts/fetch_10x_panels.R` writes a `manifest.csv` enumerating the
#' panels it fetched with display name, species, tissue, gene count, and
#' source URL. This loader is lenient: it scans the directory for every
#' `*.csv` (other than the manifest), and any file not covered by the
#' manifest is registered with filename-stem metadata so a user can drop
#' in their own panel CSV and it appears in the dropdowns alongside the
#' canonical 10x panels.
#'
#' @return list(manifest = data.frame, panels = named list of data.frames).
#'   Both are empty if the directory is absent.
load_reference_panels <- function(path = app_paths$reference_panels) {
  empty_manifest <- data.frame(
    panel_id = character(), display_name = character(),
    species = character(), tissue = character(),
    n_genes = integer(), file = character(),
    source_url = character(), stringsAsFactors = FALSE
  )
  empty <- list(manifest = empty_manifest,
                panels = stats::setNames(list(), character()))
  if (!dir.exists(path)) return(empty)

  mf_path <- file.path(path, "manifest.csv")
  manifest <- if (file.exists(mf_path)) {
    m <- .read_csv(mf_path)
    if (!nrow(m) || !"panel_id" %in% names(m)) empty_manifest else m
  } else empty_manifest
  # Backfill any optional columns so downstream code can rely on the shape.
  for (col in c("display_name", "species", "tissue", "source_url")) {
    if (!col %in% names(manifest)) manifest[[col]] <- NA_character_
  }
  if (!"n_genes" %in% names(manifest)) manifest$n_genes <- NA_integer_
  if (!"file"    %in% names(manifest)) manifest$file    <- NA_character_

  csvs <- list.files(path, pattern = "\\.csv$", full.names = FALSE)
  csvs <- setdiff(csvs, "manifest.csv")
  known <- !is.na(manifest$file) & manifest$file %in% csvs
  manifest <- manifest[known, , drop = FALSE]

  # Build manifest rows for any CSV that exists but isn't yet listed.
  extras <- setdiff(csvs, manifest$file)
  if (length(extras)) {
    add <- data.frame(
      panel_id     = sub("\\.csv$", "", extras),
      display_name = sub("\\.csv$", "", extras),
      species      = NA_character_,
      tissue       = NA_character_,
      n_genes      = NA_integer_,
      file         = extras,
      source_url   = NA_character_,
      stringsAsFactors = FALSE
    )
    manifest <- rbind(manifest, add[, names(manifest), drop = FALSE])
  }

  if (!nrow(manifest)) return(empty)

  panels <- stats::setNames(
    lapply(seq_len(nrow(manifest)), function(i) {
      fp <- file.path(path, manifest$file[i])
      if (!file.exists(fp)) return(NULL)
      .normalize_reference_panel_df(.read_csv(fp))
    }),
    manifest$panel_id
  )
  keep <- !vapply(panels, is.null, logical(1))
  manifest <- manifest[keep, , drop = FALSE]
  panels   <- panels[keep]

  # Fill in missing n_genes from the actual file contents.
  needs_count <- is.na(manifest$n_genes) | manifest$n_genes <= 0L
  if (any(needs_count)) {
    manifest$n_genes[needs_count] <- vapply(
      panels[needs_count], nrow, integer(1)
    )
  }

  list(manifest = manifest, panels = panels)
}

#' Read tissue-agnostic subpanel definitions.
#' @return named list of data.frames keyed by file stem.
load_shared_subpanels <- function(path = app_paths$subpanels_shared) {
  if (!dir.exists(path)) {
    return(stats::setNames(list(), character()))
  }
  files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) return(stats::setNames(list(), character()))
  keys <- sub("\\.csv$", "", basename(files))
  stats::setNames(lapply(files, .read_csv), keys)
}

#' Index every tissue directory under `root` by its manifest's `tissue_id`.
#'
#' The manifest is authoritative for tissue identity — folder names are
#' arbitrary discovery anchors, not the source of truth. Returns a named
#' character vector mapping `tissue_id` -> absolute folder path. Folders
#' missing a manifest, missing a `tissue_id`, or with duplicate ids are
#' surfaced via warning and skipped.
tissues_index <- function(root = app_paths$tissues_root) {
  if (!dir.exists(root)) return(stats::setNames(character(), character()))
  candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  ids   <- character()
  paths <- character()
  for (d in candidates) {
    mf <- file.path(d, "manifest.yml")
    if (!file.exists(mf)) next
    m <- tryCatch(yaml::read_yaml(mf), error = function(e) {
      warning("manifest.yml at ", mf, " failed to parse: ",
              conditionMessage(e), "; skipping")
      NULL
    })
    if (is.null(m)) next
    tid <- m$tissue_id
    if (is.null(tid) || !nzchar(tid)) {
      warning("manifest.yml at ", mf, " has no tissue_id; skipping")
      next
    }
    if (tid %in% ids) {
      warning("duplicate tissue_id '", tid, "' at ", d,
              "; already registered from ", paths[match(tid, ids)],
              ". Keeping the first.")
      next
    }
    ids   <- c(ids,   tid)
    paths <- c(paths, normalizePath(d, mustWork = FALSE))
  }
  stats::setNames(paths, ids)
}

#' List tissue ids known under `root`.
available_tissues <- function(root = app_paths$tissues_root) {
  unname(names(tissues_index(root)))
}

#' Resolve a tissue id to its folder path (manifest-driven).
tissue_dir_for <- function(tissue_id, root = app_paths$tissues_root) {
  idx <- tissues_index(root)
  if (!(tissue_id %in% names(idx))) {
    stop("tissue '", tissue_id, "' not found under ", root,
         ". Known: ", paste(names(idx), collapse = ", "))
  }
  unname(idx[tissue_id])
}

#' Read a tissue payload (manifest + subpanels + audit + optional custom).
#' @return list(manifest, summary, subpanels, audit, excluded, hIO, custom).
load_tissue <- function(tissue_id, root = app_paths$tissues_root) {
  tdir <- tissue_dir_for(tissue_id, root)
  manifest_path <- file.path(tdir, "manifest.yml")
  manifest <- yaml::read_yaml(manifest_path)

  summary_path <- file.path(tdir, "subpanel_summary.csv")
  summary_df <- if (file.exists(summary_path)) .read_csv(summary_path) else
    data.frame(subpanel = character(), n_genes = integer(),
               description = character(),
               pct_of_5K_kept_pool = numeric(),
               stringsAsFactors = FALSE)

  subp_dir <- file.path(tdir, "subpanels")
  subpanels <- if (dir.exists(subp_dir)) {
    files <- list.files(subp_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(files)) {
      stats::setNames(lapply(files, .read_csv),
                      sub("\\.csv$", "", basename(files)))
    } else stats::setNames(list(), character())
  } else stats::setNames(list(), character())

  audit_path <- file.path(tdir, "audit", "xenium5k_audit.csv")
  audit_df <- if (file.exists(audit_path)) .read_csv(audit_path) else
    data.frame(gene = character(), stringsAsFactors = FALSE)

  excluded_path <- file.path(tdir, "audit", "xenium5k_excluded.csv")
  excluded_df <- if (file.exists(excluded_path)) .read_csv(excluded_path) else
    data.frame(gene = character(), stringsAsFactors = FALSE)

  hIO_dir <- file.path(tdir, "audit")
  hIO <- if (dir.exists(hIO_dir)) {
    files <- list.files(hIO_dir, pattern = "^hIO_.*\\.csv$",
                        full.names = TRUE)
    if (length(files)) {
      stats::setNames(lapply(files, .read_csv),
                      sub("\\.csv$", "", basename(files)))
    } else stats::setNames(list(), character())
  } else stats::setNames(list(), character())

  custom_df <- NULL
  cp <- manifest$custom_panel
  if (!is.null(cp) && !is.null(cp$file)) {
    cp_path <- file.path(tdir, cp$file)
    if (file.exists(cp_path)) {
      custom_df <- .read_csv(cp_path)
    } else {
      warning("custom_panel.file declared in manifest but not found: ", cp_path)
    }
  }

  list(
    id        = tissue_id,
    manifest  = manifest,
    summary   = summary_df,
    subpanels = subpanels,
    audit     = audit_df,
    excluded  = excluded_df,
    hIO       = hIO,
    custom    = custom_df,
    dir       = tdir
  )
}

# ---------------------------------------------------------------------------
# Join helpers — restore the legacy per-subpanel schema by left-joining
# audit columns + reference_5k biology columns onto each gene table.
# ---------------------------------------------------------------------------

#' Audit-table column names for the given reference run IDs.
#' Used by validators that need to know what columns to expect.
audit_detection_cols <- function(reference_runs) {
  if (!length(reference_runs)) return(character())
  c(paste0("detection_pct_", reference_runs),
    paste0("n_cells_",       reference_runs),
    "log2_detection_ratio")
}

# Left-join `df` on `gene` against `ref` (data.frame). Columns from `ref`
# that already exist on `df` are skipped — caller-supplied values win.
.left_join_gene <- function(df, ref) {
  if (is.null(ref) || !nrow(ref) || !"gene" %in% names(df) ||
      !"gene" %in% names(ref)) {
    return(df)
  }
  add_cols <- setdiff(names(ref), c("gene", names(df)))
  if (!length(add_cols)) return(df)
  ix <- match(df$gene, ref$gene)
  for (col in add_cols) df[[col]] <- ref[[col]][ix]
  df
}

# Apply audit + ref5k joins to a gene table.
.annotate_with_tissue <- function(df, tissue, ref5k) {
  df <- .left_join_gene(df, tissue$audit)
  df <- .left_join_gene(df, ref5k)
  df
}

# ---------------------------------------------------------------------------
# Composer — returns the legacy `panels` shape
# ---------------------------------------------------------------------------

#' Compose the in-memory `panels` object for a tissue.
#'
#' Returns a list with the legacy keys downstream modules consume:
#'   summary, subpanels (named list), custom (data.frame or NULL),
#'   xenium5k (joined audit-as-reference table), excluded, hIO, meta.
#' Plus tissue-aware keys: tissue (id + manifest), reference_5k, reference_runs.
load_panels <- function(tissue_id = NULL,
                        reference_5k_path = app_paths$reference_5k,
                        reference_panels_path = app_paths$reference_panels,
                        subpanels_shared_path = app_paths$subpanels_shared,
                        tissues_root = app_paths$tissues_root) {

  ref5k <- load_reference_5k(reference_5k_path)
  refp  <- load_reference_panels(reference_panels_path)

  tids <- available_tissues(tissues_root)
  if (is.null(tissue_id)) {
    if (!length(tids)) stop("no tissues found under ", tissues_root)
    tissue_id <- tids[1]
  } else if (!(tissue_id %in% tids)) {
    stop("tissue '", tissue_id, "' not found. Available: ",
         paste(tids, collapse = ", "))
  }
  tissue <- load_tissue(tissue_id, root = tissues_root)

  # Tissues may opt out of the shared subpanel pool by setting
  # `use_shared_subpanels: false` in their manifest. Default TRUE so
  # tissues defined before this flag (e.g. pancreas) keep their behaviour.
  use_shared <- isTRUE(tissue$manifest$use_shared_subpanels %||% TRUE)
  shared <- if (use_shared) load_shared_subpanels(subpanels_shared_path) else
    stats::setNames(list(), character())

  # Merge shared + tissue-specific subpanels (tissue wins on collision),
  # then annotate every subpanel with audit + biology columns.
  merged_subp <- shared
  for (k in names(tissue$subpanels)) merged_subp[[k]] <- tissue$subpanels[[k]]
  annotated_subp <- lapply(merged_subp, .annotate_with_tissue,
                           tissue = tissue, ref5k = ref5k)

  # Custom panel — optional. Annotate when present.
  custom_df <- tissue$custom
  if (!is.null(custom_df) && nrow(custom_df)) {
    custom_df <- .annotate_with_tissue(custom_df, tissue, ref5k)
  }

  # xenium5k = audit table joined with reference_5k biology (preserves
  # the schema legacy code expected from `xenium5k_in_audit.csv`).
  xenium5k_df <- .left_join_gene(tissue$audit, ref5k)
  # If reference_5k has genes not present in audit, still surface them so
  # validators see the full 5K reference universe.
  missing_ref <- setdiff(ref5k$gene, xenium5k_df$gene)
  if (length(missing_ref)) {
    add <- ref5k[ref5k$gene %in% missing_ref, , drop = FALSE]
    for (col in setdiff(names(xenium5k_df), names(add))) add[[col]] <- NA
    xenium5k_df <- rbind(xenium5k_df, add[, names(xenium5k_df), drop = FALSE])
  }

  list(
    summary        = tissue$summary,
    subpanels      = annotated_subp,
    custom         = custom_df,
    xenium5k       = xenium5k_df,
    excluded       = tissue$excluded,
    hIO            = tissue$hIO,
    reference_5k   = ref5k,
    reference_panels          = refp$panels,
    reference_panels_manifest = refp$manifest,
    reference_runs = as.character(tissue$manifest$reference_runs %||% character()),
    tissue         = list(
      id           = tissue_id,
      display_name = tissue$manifest$display_name %||% tissue_id,
      manifest     = tissue$manifest,
      dir          = tissue$dir
    ),
    meta           = list(
      tissue_id        = tissue_id,
      reference_5k_dir = normalizePath(reference_5k_path, mustWork = FALSE),
      shared_dir       = normalizePath(subpanels_shared_path, mustWork = FALSE),
      tissue_dir       = normalizePath(tissue$dir, mustWork = FALSE),
      mtime            = Sys.time(),
      subpanel_keys    = sort(names(annotated_subp))
    )
  )
}

# ---------------------------------------------------------------------------
# Subpanel-key resolution + custom-panel display labels
# ---------------------------------------------------------------------------

#' Resolve a `subpanel` value from `summary` to a `subpanels` list key.
#'
#' The summary's `subpanel` column uses the long name
#' (`99d_truly_unannotated_subset_of_99c`) whereas the files on disk use
#' the short stem (`99d_truly_unannotated`). Match by prefix.
resolve_subpanel_key <- function(summary_subpanel, panels) {
  keys <- panels$meta$subpanel_keys
  hit  <- keys[startsWith(summary_subpanel, keys) |
               startsWith(keys, summary_subpanel)]
  if (length(hit) == 1L) hit else NA_character_
}

#' Stable opaque key for the custom-panel slot (tissue-agnostic).
custom_panel_slot_key <- function() "custom_panel"

#' Display label for the custom-panel slot.
#'
#' Priority:
#'   1. Uploaded panel — basename of the upload (no extension).
#'   2. Tissue manifest `custom_panel.display_name`.
#'   3. Generic fallback `"custom_panel"`.
custom_panel_label <- function(status, panels = NULL) {
  if (!is.null(status) && !is.null(status$source_name)) {
    return(tools::file_path_sans_ext(status$source_name))
  }
  if (!is.null(panels) && !is.null(panels$tissue$manifest$custom_panel$display_name)) {
    return(panels$tissue$manifest$custom_panel$display_name)
  }
  "custom_panel"
}

#' Pass-through display-label mapper for any panel key.
#'
#' Rewrites:
#'   * the custom-panel slot   -> uploaded filename or manifest display name
#'   * a `ref10x:<panel_id>` key -> the 10x panel's display name (from
#'                                  `panels$reference_panels_manifest`)
#' Every other key is returned as-is.
panel_display_label <- function(key, status, panels = NULL) {
  if (identical(key, custom_panel_slot_key())) {
    return(custom_panel_label(status, panels))
  }
  if (is.character(key) && length(key) == 1L && startsWith(key, "ref10x:") &&
      !is.null(panels)) {
    pid <- sub("^ref10x:", "", key)
    mf  <- panels$reference_panels_manifest
    if (!is.null(mf) && nrow(mf)) {
      hit <- which(mf$panel_id == pid)
      if (length(hit)) return(mf$display_name[hit[1]])
    }
    return(pid)
  }
  key
}

# ---------------------------------------------------------------------------
# Custom-panel uploads (session override)
# ---------------------------------------------------------------------------

#' Canonical column schema for the custom-panel slot when no upload has
#' been provided. Composed dynamically from the tissue's reference runs.
custom_panel_canonical_cols <- function(reference_runs = character()) {
  c("gene", "category", "rationale", "single_gene_panel_note",
    "exclude_recommended",
    audit_detection_cols(reference_runs))
}

#' Parse a user-uploaded CSV and enrich it to the canonical custom-panel
#' schema by filling missing columns from the active tissue's custom
#' panel (if any) and the tissue audit table for overlapping genes.
#'
#' @param csv_path Path to the uploaded CSV.
#' @param panels   Output of `load_panels()`. Enrichment sources:
#'                 `panels$custom` (tissue's default custom panel — may be
#'                 NULL when the tissue has none) then `panels$xenium5k`
#'                 (the joined audit reference).
#' @param source_name Optional label (defaults to `basename(csv_path)`).
#' @return list(df, n_genes, n_enriched_from_custom, n_enriched_from_5k,
#'              n_unmatched, source_name)
upload_custom_panel <- function(csv_path, panels, source_name = NULL) {
  if (!file.exists(csv_path)) {
    stop("upload_custom_panel: file not found: ", csv_path)
  }
  if (is.null(source_name)) source_name <- basename(csv_path)

  df <- data.table::fread(csv_path, na.strings = c("", "NA"),
                          data.table = FALSE, showProgress = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) {
    stop("Uploaded CSV is empty.")
  }

  gene_col <- grep("^gene$", names(df), ignore.case = TRUE, value = TRUE)
  if (length(gene_col) == 0L) {
    stop("Uploaded CSV must have a `gene` column (case-insensitive). ",
         "Found columns: ", paste(names(df), collapse = ", "))
  }
  if (gene_col[1] != "gene") names(df)[names(df) == gene_col[1]] <- "gene"

  df$gene <- trimws(as.character(df$gene))
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  df <- df[!duplicated(df$gene), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("Uploaded CSV has no valid (non-empty) gene rows.")
  }

  canon <- custom_panel_canonical_cols(panels$reference_runs %||% character())
  for (col in canon) if (!(col %in% names(df))) df[[col]] <- NA
  extra <- setdiff(names(df), canon)
  df <- df[, c(canon, extra), drop = FALSE]

  fill_cols <- setdiff(canon, "gene")
  src_custom <- panels$custom        # may be NULL when tissue has none
  src_5k     <- panels$xenium5k

  n_from_custom <- 0L
  n_from_5k     <- 0L
  n_unmatched   <- 0L

  for (i in seq_len(nrow(df))) {
    g <- df$gene[i]
    matched <- FALSE

    if (!is.null(src_custom) && nrow(src_custom) > 0L) {
      hit <- which(src_custom$gene == g)
      if (length(hit) >= 1L) {
        filled_here <- FALSE
        for (col in fill_cols) {
          if (!(col %in% names(src_custom))) next
          if (is.na(df[[col]][i]) && !is.na(src_custom[[col]][hit[1]])) {
            df[[col]][i] <- src_custom[[col]][hit[1]]
            filled_here <- TRUE
          }
        }
        if (filled_here) n_from_custom <- n_from_custom + 1L
        matched <- TRUE
      }
    }

    if (!is.null(src_5k) && nrow(src_5k) > 0L) {
      hit <- which(src_5k$gene == g)
      if (length(hit) >= 1L) {
        filled_here <- FALSE
        for (col in fill_cols) {
          if (!(col %in% names(src_5k))) next
          if (is.na(df[[col]][i]) && !is.na(src_5k[[col]][hit[1]])) {
            df[[col]][i] <- src_5k[[col]][hit[1]]
            filled_here <- TRUE
          }
        }
        if (filled_here && !matched) n_from_5k <- n_from_5k + 1L
        matched <- TRUE
      }
    }

    if (!matched) n_unmatched <- n_unmatched + 1L
  }

  list(
    df                     = df,
    n_genes                = nrow(df),
    n_enriched_from_custom = n_from_custom,
    n_enriched_from_5k     = n_from_5k,
    n_unmatched            = n_unmatched,
    source_name            = source_name
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
