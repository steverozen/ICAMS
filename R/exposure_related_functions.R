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

#' @title Write an exposure matrix to a file
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param file File to which to write the exposure matrix (as a CSV file).
#'
#' @importFrom utils write.csv
#' 
#' @export
WriteExposure <- function(exposure, file) {
  old.digits <- getOption("digits")
  options(digits = 22)
  write.csv(exposure, file, row.names = TRUE)
  on.exit(options(digits = old.digits)) 
}

#' Sort columns of an exposure matrix from largest to smallest (or vice versa)
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'   
#' @param decreasing If \code{TRUE}, sort from largest to smallest.
#' 
#' @return The original \code{exposure} with columns sorted.
#'
#' @export
SortExposure <- function(exposure, decreasing = TRUE) {
  retval <- exposure[, order(colSums(exposure), decreasing = decreasing)]
  return(retval)
}