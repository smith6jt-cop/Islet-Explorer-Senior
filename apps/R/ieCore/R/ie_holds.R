# Holds -- a hand-written file whose exact content gates whether a blocked stage is allowed
# to run. `core/phenocycler`'s NO_VERIFIED_TYPING latch is the model this borrows from: a
# blocked stage stays blocked until something says otherwise, and "something says otherwise"
# must be a fact this code can check, not an assumption a caller gets to make on its own
# authority.
#
# Lifting a hold is two claims, not one, and both must hold together:
#   1. the hold file says, exactly, "HOLD LIFTED" -- not "hold lifted", not that text with
#      extra prose around it, because an exact match is the one thing a stray edit or a
#      copy-pasted paragraph cannot produce by accident;
#   2. the retired pin that recorded the hold's reason is gone. A hold file that says
#      lifted while its retired pin still sits on disk is two records disagreeing, and this
#      code trusts neither of them in that state.

#' Assert a hold is lifted: the hold file matches exactly, and the retired pin is gone
#'
#' @param hold_file path to the file whose entire content must be exactly `"HOLD LIFTED"`
#' @param retired_pin path to the file that must be ABSENT once the hold is lifted
#' @param context a short label naming what is asking, used only in the error message
#' @return `invisible(TRUE)` when both conditions hold; otherwise `stop()`s
#' @export
ie_assert_hold_lifted <- function(hold_file, retired_pin, context = "this stage") {
  if (length(hold_file) != 1L || is.na(hold_file) || !nzchar(hold_file)) {
    stop("ie_assert_hold_lifted(): hold_file must be one non-empty string", call. = FALSE)
  }
  if (length(retired_pin) != 1L || is.na(retired_pin) || !nzchar(retired_pin)) {
    stop("ie_assert_hold_lifted(): retired_pin must be one non-empty string", call. = FALSE)
  }
  if (!file.exists(hold_file)) {
    stop("ie_assert_hold_lifted(): no hold file at '", hold_file, "' for ", context,
         ". A hold with no recorded state is not lifted.", call. = FALSE)
  }
  content <- paste(readLines(hold_file, warn = FALSE), collapse = "\n")
  content <- trimws(content)
  if (!identical(content, "HOLD LIFTED")) {
    stop("ie_assert_hold_lifted(): '", hold_file, "' does not read exactly 'HOLD LIFTED' ",
         "(read: '", content, "'). ", context, " stays blocked.", call. = FALSE)
  }
  if (file.exists(retired_pin)) {
    stop("ie_assert_hold_lifted(): '", hold_file, "' says lifted, but the retired pin '",
         retired_pin, "' is still present for ", context, ". Remove the pin before the hold ",
         "can be trusted; two records disagreeing is not a lifted hold.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Whether a hold is lifted, without stopping
#'
#' A thin, non-throwing wrapper over [ie_assert_hold_lifted()] for callers that need a
#' boolean -- a card deciding whether to render, for instance -- rather than an error to
#' propagate.
#' @export
ie_hold_is_lifted <- function(hold_file, retired_pin) {
  tryCatch({
    ie_assert_hold_lifted(hold_file, retired_pin)
    TRUE
  }, error = function(e) FALSE)
}
