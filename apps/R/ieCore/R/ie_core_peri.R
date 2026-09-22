# The Plot tab's Core+Peri measure is one proportion made from two: the core's and the
# peri-islet expansion's, each over its own cell total. Two absences meet in that arithmetic and
# they are not the same thing:
#
#   a side with NO CELLS                             a measured absence. The islet has no
#   (total 0, or total AND proportion both null)     peri-islet expansion, or nothing was counted
#                                                    there, so that side contributes nothing and
#                                                    the other stands alone. This is what
#                                                    mod_plot_server.R has always done -- "an islet
#                                                    with no peri zone still contributes its core,
#                                                    and vice versa" -- and it stays.
#   a side WITH CELLS and a null proportion,         a missing measurement. The producer recorded
#   or a side with a MEASURED proportion and a       one half of the pair and not the other
#   null total                                       (build_galectin_layer.py writes NaN where a
#                                                    count is blank or its two object-classifier
#                                                    passes disagree), and the islet's combined
#                                                    value is UNAVAILABLE: NA, never 0 over a kept
#                                                    denominator and never the other side alone.
#
# Until 2026-09-21 the module folded the second case into the first --
# `ifelse(is.finite(prop) & is.finite(total), prop * total, 0)` -- so the unmeasured side became 0
# cells, the denominator kept that side's cells, and the Plot tab drew the islet at exactly 0.
# On the 2026-09-14 galectin layer that was 719 of 15,557 islet-regions. Until 2026-09-22 a null
# total beside a finite proportion was read the same way as no cells there, and the islet was
# drawn at its other side's value alone. The arithmetic lives here so it can be tested without
# the module.
#
# THE BOUNDARY. Each column is read once, in .ie_read_side(), and a value that cannot be what
# the column holds is refused there, with a reason on the row: a total is a count (finite and
# non-negative) or null, and a refused total is a null total; a proportion is in [0, 1] or null,
# and a refused proportion makes the islet unavailable outright. Nothing downstream --
# Counts, Density, Percentage in the module -- tests for these again; an NA `n` carries
# through them.

#' Combine a core and a peri-islet proportion into one Core+Peri cell count and total
#'
#' Vectorised over islets; a scalar argument recycles. Returns a data.frame with one row per
#' islet and three columns:
#' \describe{
#'   \item{n}{the combined positive-cell count, `core_prop * core_total + peri_prop * peri_total`,
#'     a side with no cells contributing 0; `NA` where the islet is unavailable (below)}
#'   \item{total}{the combined cell total, a null side total counting as 0}
#'   \item{reason}{`NA` where `n` is a number; otherwise why it is not, one clause per refused
#'     or missing column, sides joined with `; `}
#' }
#' The caller derives its metric from these -- `100 * n / total` where `total > 0`, `round(n)`,
#' `n / area` -- and an `NA` `n` carries through every one of them. `total` is left as measured
#' even when `n` is `NA`: the denominator is a fact about the islet, and the absence is in the
#' numerator.
#'
#' **The rule, per side**, after each column has been read through `as.numeric()`:
#' \itemize{
#'   \item a total that is finite and non-negative is a count; anything else -- blank,
#'     unparseable, `NA`, `NaN`, `Inf`, `-Inf`, negative -- is a null total. A negative or
#'     infinite total is refused with a reason when that refusal is what makes the islet
#'     unavailable.
#'   \item a proportion that is finite and in `[0, 1]` is a measurement; `NA`/`NaN` is a null
#'     proportion; a finite value outside `[0, 1]` is not a proportion and the islet's combined
#'     value is `NA` with a reason, whatever the total beside it.
#'   \item a count of 0, or a null total beside a null proportion, is **no cells there**: the side
#'     contributes 0 to `n` and 0 to `total` and the other side stands alone (decided
#'     2026-09-22, preserving what the module always did for an islet with no peri-islet
#'     expansion).
#'   \item a positive count beside a null proportion, or a null total beside a measured
#'     proportion, is a **missing measurement**: the islet's `n` is `NA` with a reason.
#' }
#'
#' @param core_prop,peri_prop proportions in [0, 1] on each side, or null
#' @param core_total,peri_total cell counts on each side, or null
#' @return a data.frame with columns `n`, `total`, `reason`
#' @export
ie_combine_core_peri <- function(core_prop, peri_prop, core_total, peri_total) {
  num <- function(x) suppressWarnings(as.numeric(x))
  len <- max(length(core_prop), length(peri_prop), length(core_total), length(peri_total))
  if (len == 0L) {
    return(data.frame(n = numeric(0), total = numeric(0), reason = character(0),
                      stringsAsFactors = FALSE))
  }
  core <- .ie_read_side(rep_len(num(core_prop), len), rep_len(num(core_total), len), "core")
  peri <- .ie_read_side(rep_len(num(peri_prop), len), rep_len(num(peri_total), len),
                        "peri-islet expansion")

  n <- core$n + peri$n
  n[core$unavailable | peri$unavailable] <- NA_real_
  total <- core$total + peri$total
  reason <- ifelse(nzchar(core$reason) & nzchar(peri$reason),
                   paste(core$reason, peri$reason, sep = "; "),
                   paste0(core$reason, peri$reason))
  reason[!nzchar(reason)] <- NA_character_
  data.frame(n = n, total = total, reason = reason, stringsAsFactors = FALSE)
}

# One side's columns, read once. `prop` and `total` are already numeric (NA where the source
# was blank or unparseable) and the same length; `side` names the side in the reason.
.ie_read_side <- function(prop, total, side) {
  fmt <- function(x) format(x, trim = TRUE)

  # The count column: finite and non-negative, or null. `is.finite()` is FALSE for NA, NaN,
  # Inf and -Inf, so none of them can leak through `>= 0` as TRUE.
  is_count <- is.finite(total) & total >= 0
  refused_total <- !is.na(total) & !is_count           # -5, Inf, -Inf: a number, not a count
  has_cells <- is_count & total > 0

  # The proportion column: in [0, 1], or null.
  is_number <- is.finite(prop)
  refused_prop <- is_number & (prop < 0 | prop > 1)     # 1.5, -0.2: a number, not a proportion
  is_prop <- is_number & !refused_prop

  n <- ifelse(has_cells & is_prop, prop * total, 0)
  t <- ifelse(is_count, total, 0)

  # Why the islet is unavailable, when it is. The three are exclusive: a refused proportion
  # is a number, a lost proportion is not; a lost proportion sits beside a count, a lost
  # total does not.
  lost_prop <- has_cells & !is_number
  lost_total <- !is_count & is_prop
  reason <- rep("", length(total))
  reason[refused_prop] <- sprintf("%s proportion %s is not a proportion",
                                  side, fmt(prop[refused_prop]))
  reason[lost_prop] <- sprintf("%s proportion not measured with %s %s cells present",
                               side, fmt(total[lost_prop]), side)
  lost_blank <- lost_total & !refused_total
  lost_refused <- lost_total & refused_total
  reason[lost_blank] <- sprintf("%s total not measured with %s proportion %s present",
                                side, side, fmt(prop[lost_blank]))
  reason[lost_refused] <- sprintf("%s total %s is not a count, with %s proportion %s present",
                                  side, fmt(total[lost_refused]), side, fmt(prop[lost_refused]))

  list(n = n, total = t, unavailable = refused_prop | lost_prop | lost_total, reason = reason)
}
