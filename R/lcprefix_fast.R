#' Compute length of longest common prefix of two strings
#'
#' Fast alternative to \code{Biostrings::lcprefix} that avoids
#' S4 method dispatch overhead.
#'
#' @param a A single character string.
#' @param b A single character string.
#'
#' @return An integer: the number of leading characters that match.
#'
#' @keywords internal
lcprefix_fast <- function(a, b) {
  min_len <- min(nchar(a), nchar(b))
  if (min_len == 0L) return(0L)
  a_raw <- utf8ToInt(substr(a, 1L, min_len))
  b_raw <- utf8ToInt(substr(b, 1L, min_len))
  first_diff <- which(a_raw != b_raw)[1L]
  if (is.na(first_diff)) min_len else as.integer(first_diff - 1L)
}
