#' Compare a loaded dataset's gene set against the panel reference.
#'
#' Reference set = `panels$reference_5k$gene` ∪ `panels$custom$gene` (when
#' a custom panel is configured for the active tissue). The validation
#' never blocks analysis; it returns a structured report the UI can render.
#'
#' @param genes character vector of gene symbols present in the dataset.
#' @param panels output of [load_panels()].
#' @return named list:
#'   * `n_data`      — genes in the dataset
#'   * `n_reference` — genes in the reference union
#'   * `intersection`, `missing_from_data`, `extra_in_data` — gene vectors
#'   * `n_*` — corresponding counts
#'   * `pct_reference_covered` — `intersection / n_reference`
#'   * `partial`    — TRUE if not a complete cover
#'   * `has_custom` — TRUE if a custom panel was part of the reference union
panel_validate <- function(genes, panels) {
  if (!length(genes)) {
    stop("`genes` is empty; nothing to validate.")
  }
  data_g <- unique(as.character(genes))
  ref5k_g  <- panels$reference_5k$gene %||% character()
  custom_g <- if (!is.null(panels$custom)) panels$custom$gene else character()
  ref_g  <- unique(c(ref5k_g, custom_g))

  inter <- intersect(data_g, ref_g)
  miss  <- setdiff(ref_g, data_g)
  extra <- setdiff(data_g, ref_g)

  list(
    n_data                 = length(data_g),
    n_reference            = length(ref_g),
    intersection           = inter,
    missing_from_data      = miss,
    extra_in_data          = extra,
    n_intersection         = length(inter),
    n_missing_from_data    = length(miss),
    n_extra_in_data        = length(extra),
    pct_reference_covered  = if (length(ref_g)) length(inter) / length(ref_g) else NA_real_,
    partial                = length(miss) > 0L || length(extra) > 0L,
    has_custom             = length(custom_g) > 0L
  )
}

#' Render a one-paragraph human-readable summary of [panel_validate()].
#'
#' @param custom_label Display label for the custom-panel slot when the
#'   reference union includes one. Pass the uploaded-panel stem (or the
#'   tissue manifest's display_name) to keep the summary consistent with
#'   the rest of the UI.
panel_validate_summary <- function(rep, custom_label = NULL) {
  ref_phrase <- if (isTRUE(rep$has_custom)) {
    sprintf("xenium5k ∪ %s", custom_label %||% "custom panel")
  } else {
    "xenium5k"
  }
  msg <- sprintf(
    paste0(
      "Loaded dataset: %d genes. Reference (%s): %d genes. ",
      "Intersection: %d (%.1f%% of reference). Missing from data: %d. ",
      "Extra in data (not in reference): %d."),
    rep$n_data, ref_phrase, rep$n_reference, rep$n_intersection,
    100 * rep$pct_reference_covered,
    rep$n_missing_from_data, rep$n_extra_in_data
  )
  if (rep$partial) {
    msg <- paste0(msg,
                  " The match is partial — analysis is not blocked, but ",
                  "downstream subpanel selection will only score genes that ",
                  "are present in the loaded data.")
  }
  msg
}
