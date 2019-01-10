#' @include utility_functions.R
NULL



#' Collapse a SNS 192 catalog matrix to a 96 catalog matrix
#'
#' @param cat192 A SNS 192 catalog matrix whose row names indicate the 192
#'   mutation types while its columns show the occurrences of each mutation type of
#'   different samples.
#'
#' @return A SNS 96 catalog matrix whose row names indicate the 96
#    mutation types while its columns show the occurrences of
#    each mutation type of different samples.
#' @export
Collapse192to96 <- function(cat192) {
  dt192 <- data.table(cat192)
  dt192$rn <- PyrTri(rownames(cat192))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[ , -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[.catalog.row.order96, , drop = FALSE]
}

#' Collapse a SNS 1536 catalog matrix to a 96 catalog matrix
#'
#' @param cat1536 A SNS 1536 catalog matrix whose row names indicate the 1536
#'   mutation types while its columns show the occurrences of each mutation type
#'   of different samples.
#' @return A SNS 96 catalog matrix whose row names indicate the 96
#'   mutation types while its columns show the occurrences of
#'   each mutation type of different samples.
#' @export
Collapse1536to96 <- function(cat1536) {
  dt <- data.table(cat1536)
  rn <- rownames(cat1536)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # sring for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[.catalog.row.order96, , drop = FALSE]
}
