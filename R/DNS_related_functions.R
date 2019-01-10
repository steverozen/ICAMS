#' @include VCF_related_functions.R utility_functions.R
NULL

#' CanonicalizeDNS
#'
#' @param ref.vec TODO
#' @param alt.vec TODO
#'
#' @return TODO
#' @export
CanonicalizeDNS <- function(ref.vec, alt.vec) {
  # TODO document

  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  Canonicalize1DNS <- function(DNS) {
    if (DNS %in% .catalog.row.order.DNS.78) {
      return(DNS)
    } else {
      ref <- substr(DNS, 1, 2)
      alt <- substr(DNS, 3, 4)
      out <- paste0(revc(ref), revc(alt))
    }
    stopifnot(out %in% .catalog.row.order.DNS.78)
    return(out)
  }
  ret <- sapply(paste0(ref.vec, alt.vec), FUN = Canonicalize1DNS)
  return(ret)
}

#' CanonicalizeQUAD
#'
#' @param quad TODO
#'
#' @return TODO
#' @export
CanonicalizeQUAD <- function(quad) {
  # TODO document

  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  Canonicalize1QUAD <- function(quad) {
    if (quad %in% .catalog.row.order.QUAD.136) {
      return(quad)
    } else {
      out <- revc(quad)
      stopifnot(out %in% .catalog.row.order.QUAD.136)
      return(out)
    }
  }

  ret <- sapply(quad, FUN = Canonicalize1QUAD)
  return(ret)
}

#' Collapse a DNS 144 catalog matrix to a DNS 78 catalog matrix
#'
#' @param catDNS144 A DNS 144 catalog matrix whose row names indicate the 192
#'   mutation types while its columns show the occurrences of each mutation type of
#'   different samples.
#' @import data.table
#' @return A DNS 78 catalog matrix whose row names indicate the 96 mutation
#'   types while its columns show the occurrences of each mutation type of different
#'   samples.
#' @export
Collapse144to78 <- function(catDNS144) {
  dt144 <- data.table(catDNS144)
  ref <- substr(rownames(catDNS144), 1, 2)
  alt <- substr(rownames(catDNS144), 3, 4)
  dt144$rn <- CanonicalizeDNS(ref, alt)
  dt78 <- dt144[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat78 <- as.matrix(dt78[ , -1])
  rownames(mat78) <- dt78$rn
  mat78 <- mat78[.catalog.row.order.DNS.78, , drop = FALSE]
}
