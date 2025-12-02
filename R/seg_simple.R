#' Segment a single indel sequence using Rcpp interface
#'
#' @param ins_or_del Character indicating insertion ("i") or deletion ("d")
#' @param string A single indel sequence to segment (character scalar)
#' @param context A single flanking context sequence (character scalar)
#'
#' @return A list with the following elements:
#'   \item{unit}{Repeat unit sequence (character)}
#'   \item{unit_length}{Unit length (integer)}
#'   \item{internal_rep}{Internal repeat region (character)}
#'   \item{internal_reps}{Internal repeat count (integer)}
#'   \item{spacer}{Spacer sequence (character)}
#'   \item{spacer_length}{Spacer length (integer)}
#'   \item{prime3_rep}{3' flanking repeat region (character)}
#'   \item{prime3_reps}{3' flanking repeat count (integer)}
#'   \item{original_reps}{Original repeat count (integer)}
#'
#' @details
#' This function segments an indel sequence by finding the optimal repeat unit
#' that best explains the sequence structure. The algorithm tries all possible
#' repeat unit sizes and selects the best segmentation based on:
#' \itemize{
#'   \item Highest 3' flanking repeat count (most important)
#'   \item Highest internal repeat count
#'   \item Lowest spacer length (prefer clean repeats)
#'   \item Lowest unit length (prefer simpler units)
#' }
#'
#' @examples
#' \dontrun{
#' # Simple AT repeat
#' result <- seg_simple("d", "ATATAT", "ATATGG")
#' print(result)
#'
#' # CG repeat
#' result <- seg_simple("d", "CGCGCG", "CGCGAA")
#' print(result)
#' }
#'
#' @export
seg_simple <- function(ins_or_del, string, context) {
  # Input validation
  if (length(string) != 1 || length(context) != 1) {
    stop("string and context must each be a single character value")
  }

  if (length(ins_or_del) != 1 || !ins_or_del %in% c("i", "d")) {
    stop("ins_or_del must be either 'i' (insertion) or 'd' (deletion)")
  }

  if (
    !is.character(string) || !is.character(context) || !is.character(ins_or_del)
  ) {
    stop("All arguments must be character values")
  }

  # Call the C++ function via Rcpp
  result <- segment_simple_cpp(ins_or_del, string, context)

  return(result)
}
