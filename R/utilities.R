#' PyrTri
#'
#' @param mutstring a mutation string
#'
#' @return a mutation string
#' @export
#' @keywords internal
PyrTri <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' PyrPenta
#'
#' @param mutstring a mutation string
#'
#' @return a mutation string
#' @export
#' @keywords internal
PyrPenta <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in string.vec
#'
#' @param string.vec a vector of type character.
#'
#' @return A vector of type characters with the reverse complement of
#'   of every string in string.vec.
#' @import Biostrings
#' @keywords internal
#' @export
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' RevcSNS96
#'
#' @param mutstring a mutation string
#'
#' @return a mutation string
#' @export
RevcSNS96 <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}
