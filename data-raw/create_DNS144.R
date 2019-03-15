OLDCanonicalizeDNS <- function(ref.vec, alt.vec) {
  # TODO document

  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  Canonicalize1DNS <- function(DNS) {
    if (DNS %in% catalog.row.order$DNS78) {
      return(DNS)
    } else {
      ref <- substr(DNS, 1, 2)
      alt <- substr(DNS, 3, 4)
      if (!ref %in% canonical.ref) {
        ref <- revc(ref)
        alt <- revc(alt)
      }
      out <- paste0(ref, alt)
      if (!out %in% catalog.row.order$DNS78) {
        stopifnot(ref == revc(ref))
        out <- paste0(ref, revc(alt))
      }
      stopifnot(out %in% catalog.row.order$DNS78)
      return(out)
    }
  }
  ret <- sapply(paste0(ref.vec, alt.vec), FUN = Canonicalize1DNS)
  return(ret)
}

xx <- function() {
  xxx <- c("A", "C", "G", "T")
  retv1 <- character(0)
  retv2 <- character(0)
  for (i1 in xxx) {
    for (i2 in xxx ) {
      for (i3 in xxx) {
        if (i1 == i3) next
        for (i4 in xxx) {
          if (i2 == i4) next
          retv1 <- c(retv1, paste0(i1, i2))
          retv2 <- c(retv2, paste0(i3, i4))
        }
      }
    }
  }
  return(list(ref=retv1, alt=retv2))
}

to.test <- paste0(foo$ref, foo$alt)
result <- OLDCanonicalizeDNS(foo$ref, foo$alt)
write(deparse(to.test), "data-raw/to.test.CanonicalizeDNS.txt")
write(deparse(result), "data-raw/result.CanonicalizeDNS.txt")
