#' Read Catalog Functions
#'
#' Read a catalog in PCAWG7 format from path
#'
#' \code{ReadCat96} Read a 96 SNS catalog from path
#'
#' \code{ReadCat192} Read a 192 SNS catalog from path
#'
#' \code{ReadCat1536} Read a 1536 SNS catalog from path
#'
#' \code{ReadCatDNS78} Read a 78 DNS catalog from path
#'
#' \code{ReadCatDNS144} Read a 144 DNS catalog from path
#'
#' \code{ReadCatQUAD136} Read a 136 QUAD catalog from path
#'
#' \code{ReadCatID} Read a ID (insertion/deletion) catalog from path
#' Please take note that the deletions Repeat Size ranges from 0 to 5+
#' in the catalog, but for plotting and end user documentation
#' it ranges from 1 to 6+.
#'
#' See also \code{\link{WriteCatalog}}
#' @param path Path to a catalog on disk in the "PCAWG7" format.
#' @param strict If TRUE, do additional checks on the input, and stop if the
#'   checks fail.
#' @return A catalog in canonical in-memory format.
#' @name ReadCatalog
NULL

#' Write Catalog Functions
#'
#' Write a mutation catalog to a file on disk
#'
#' \code{WriteCat96} Write a SNS 96 mutation catalog to a file on disk
#'
#' \code{WriteCat192} Write a SNS 192 mutation catalog to a file on disk
#'
#' \code{WriteCat1536} Write a SNS 1536 mutation catalog to a file on disk
#'
#' \code{WriteCatDNS78} Write a DNS 78 mutation catalog to a file on disk
#'
#' \code{WriteCatDNS144} Write a DNS 144 mutation catalog to a file on disk
#'
#' \code{WriteCatQUAD136} Write a 136 QUAD catalog from path
#'
#' \code{WriteCatID} Write a ID (insertion/deletion) catalog to a file on disk
#' Please take note that the deletions Repeat Size ranges from 0 to 5+
#' in the catalog, but for plotting and end user documentation
#' it ranges from 1 to 6+.
#'
#' See also \code{\link{ReadCatalog}}
#' @param ct A matrix of mutation catalog.
#' @param path The path of the file to be written on disk.
#' @param strict If TRUE, do additional checks on the input,
#'   and stop if the checks fail.
#' @name WriteCatalog
NULL

#' @rdname ReadCatalog
#' @export
ReadCat96 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[2] == "Trinucleotide")
  }
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[, -(1 : 2)]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.SNS.96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[catalog.row.order.SNS.96, ]
  return(out)
}

#' @rdname ReadCatalog
#' @include utility_functions.R
#' @export
ReadCat192 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  # cos.copy <- cos # For debugging, testing
  stopifnot(nrow(cos) == 192)
  if (strict) {
    stopifnot(names(cos)[2] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[3] == "Trinucleotide")
    stopifnot(names(cos)[1] == "Strand")
  }
  ref.gt.var       <- unlist(cos[, 2])
  before.ref.after <- unlist(cos[, 3])

  ## Find the rows labeled with "T", indicating the
  ## SNS is on the transcribed (which is the *antisense*) strand.
  transcribed.strand.pos <- which(cos[, 1] == 'T')

  before.ref.after[transcribed.strand.pos] <-
    revc(before.ref.after[transcribed.strand.pos])

  var <- substring(ref.gt.var, 3, 3)
  var[transcribed.strand.pos] <- revc(var[transcribed.strand.pos])

  tmp <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(tmp == catalog.row.order.SNS.192)
  }
  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[catalog.row.order.SNS.192, ]
  return(out)
}

#' @rdname ReadCatalog
#' @export
ReadCat1536 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 1536)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation.type"))
    stopifnot(names(cos)[2] == "Pentanucleotide")
  }
  names(cos)[1:2] <- c("Mutation type", "Pentanucleotide")
  ref.gt.var       <- cos[["Mutation type"]]
  before.ref.after <- cos[["Pentanucleotide"]]
  var <- substring(ref.gt.var, 3, 3)
  out <- as.matrix(cos[ , -(1 : 2)])
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.SNS.1536)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[catalog.row.order.SNS.1536, ]
  return(out)
}

#' @rdname ReadCatalog
#' @export
ReadCatDNS78 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2)]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, catalog.row.order.DNS.78))
  if ( (length(diff1) > 0)
       &&
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
       &&
       (sort(setdiff(catalog.row.order.DNS.78, rn) ==
             c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
  ) {
    cat("using temporary hack for old DNS canonicalization\n")
    # CGAA -> CGTT
    rn[rn == "CGAA"] <- "CGTT"

    # CGAC -> CGGT
    rn[rn == "CGAC"] <- "CGGT"

    # CGGA -> CGTC
    rn[rn == "CGGA"] <- "CGTC"

    # TAAC -> TAGT
    rn[rn == "TAAC"] <- "TAGT"

    # TAAG -> TACT
    rn[rn == "TAAG"] <- "TACT"

    # TACC -> TAGG
    rn[rn == "TACC"] <- "TAGG"
  }
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.DNS.78)
  }
  out <- out[catalog.row.order.DNS.78, ]
  return(out)
}

#' @rdname ReadCatalog
#' @export
ReadCatDNS144 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 144)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2)]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.DNS.144)
  }
  out <- out[catalog.row.order.DNS.144, ]
  return(out)
}

#' @rdname ReadCatalog
#' @export
ReadCatQUAD136 <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 136)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Quad", "quad", "QUAD"))
  }
  names(cos)[1] <- "Quad"
  out <- cos[, -1]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.DNS.136)
  }
  out <- out[catalog.row.order.DNS.136, ]
  return(out)
}

#' @rdname ReadCatalog
#' @export
ReadCatID <- function(path, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 83)
  cn <- names(cos)
  ex.cn <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
  # Repeat_MH_size is the size of repeat OR microhomology (MH)
  if (strict) { for (i in 1 : 4) { stopifnot(cn[i] == ex.cn[i]) } }
  names(cos)[1 : 4] <- ex.cn
  rn <- apply(cos[, 1 : 4], MARGIN = 1, paste, collapse = ":")
  # View(data.frame(mini=rn, good=catalog.row.order.ID))
  out <- as.matrix(cos[ , -(1 : 4)])
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == catalog.row.order.ID)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[catalog.row.order.ID, ]
  return(out)
}

#' @title Write a mutation catalog to a file on disk
#'
#' @description Called by exported functions to do the
#' actual writing of the catalog to disk.
#'
#' @param ct A matrix of mutation catalog.
#' @param path The path of the file to be written on disk.
#' @param num.row The number of rows in the file to be written.
#' @param row.order The row order to be used for writing the file.
#' @param row.header The row header to be used for writing the file.
#' @param strict If TRUE, do additional checks on the input, and stop if the
#'   checks fail.
#' @return A catalog in canonical in-memory format.
#' @keywords internal
WriteCat <- function(ct, path, num.row, row.order, row.header, strict) {
  mut.categories <- rownames(ct)
  stopifnot(num.row == nrow(ct))
  if (strict) {
    stopifnot(mut.categories == row.order)
  }
  ct <- ct[row.order, ]
  DT <- as.data.table(ct)
  fwrite(cbind(row.header, DT), file = path)
}

#' @rdname WriteCatalog
#' @export
WriteCat96 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 96, catalog.row.order.SNS.96, catalog.row.headers.SNS.96, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCat192 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 192, catalog.row.order.SNS.192, catalog.row.headers.SNS.192, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCat1536 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 1536, catalog.row.order.SNS.1536,
           catalog.row.headers.SNS.1536, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCatDNS78 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 78, catalog.row.order.DNS.78,
           catalog.row.headers.DNS.78, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCatDNS144 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 144, catalog.row.order.DNS.144,
           catalog.row.headers.DNS.144, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCatQUAD136 <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 136, catalog.row.order.DNS.136,
           catalog.row.headers.DNS.136, strict)
}

#' @rdname WriteCatalog
#' @export
WriteCatID <- function(ct, path, strict = TRUE) {
  WriteCat(ct, path, 83, catalog.row.order.ID, catalog.row.headers.ID, strict)
}
