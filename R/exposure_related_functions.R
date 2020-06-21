#' @title Read an exposure matrix from a file
#'
#' @param file CSV file containing an exposure matrix.
#'
#' @param check.names Passed to \code{\link[utils]{read.csv}}.
#' \strong{IMPORTANT}: If \code{TRUE} this will replace the double
#' colon in identifiers of the form <tumor_type>::<sample_id>
#' with two periods (i.e. <tumor_type>..<sample_id>.
#' If \code{check.names} is true, generate a warning
#' if double colons were present.
#'
#' @return Matrix of exposures.
#'
#' @importFrom utils read.csv
#' 
#' @export
ReadExposure <- function(file, check.names = TRUE) {
  if (check.names) {
    headers <- read.csv(file, nrow = 1, header = FALSE, stringsAsFactors = FALSE)
    double.colon <- grep("::", unlist(headers)[-1], fixed = TRUE)
    if (length(double.colon) > 0) {
      warning(":: in sample ID replaced by ..; suggest calling with check.names = FALSE")
    }
  }
  retval <- read.csv(file, row.names = 1, check.names = check.names)
  if (any(duplicated(colnames(retval)))) {
    stop("There is duplicated column name in the input file")
  }
  return(data.matrix(retval))
}