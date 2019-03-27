#' Read catalog.
#'
#' Read a catalog in standardized format from path.
#'
#' See also \code{\link{WriteCatalog}}
#'
#' @param path Path to a catalog on disk in the standardized format.
#'
#' @param ref.genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @param region One of "genome", "exome".
#'
#' @param catalog.type One of "counts", "density",
#' "counts.signature", "density.signature".
#'
#' @param strict If TRUE, then stop if additional checks on the input fail.
#'
#' @return A catalog in standard in-memory format with attributes added.
#' See \code{\link{CreateCatalogAttribute}} for more details.
#'
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   ranges from 0 to 5+, but for plotting and end user documentation it ranges
#'   from 1 to 6+.
#'
#' @export
ReadCatalog <- function(path, ref.genome, region, catalog.type, strict = TRUE) {
  if (CheckCatalogAttribute(ref.genome, region, catalog.type)) {
    class.of.catalog <- CheckClassOfCatalogFromPath(path)
    UseMethod(generic = "ReadCatalog", object = class.of.catalog)
  }
}

#' Write a catalog
#'
#' Write a catalog to a file.
#'
#' See also \code{\link{ReadCatalog}}.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}};
#' see also \code{\link{CreateCatalogAttribute}}.
#'
#' @param path The path to the file to be created .
#'
#' @param strict If TRUE, then fail if additional checks on the input fail.
#'
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   ranges from 0 to 5+, but for plotting and end user documentation it ranges
#'   from 1 to 6+.
#'
#' @export
WriteCatalog <- function(catalog, path, strict = TRUE) {
  UseMethod(generic = "WriteCatalog")
}

ReadCatalog.SNS96 <- function(path, ref.genome, region,
                              catalog.type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 96)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation Type",
                                   "Mutation.type", "Mutation.Type"))
    stopifnot(names(cos)[2] == "Trinucleotide")
  }
  ref.gt.var       <- unlist(cos[, 1])
  before.ref.after <- unlist(cos[, 2])
  var <- substring(ref.gt.var, 3, 3)
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SNS96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SNS96, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

ReadCatalog.SNS192 <- function(path, ref.genome, region,
                               catalog.type, strict = TRUE) {
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
    stopifnot(tmp == ICAMS::catalog.row.order$SNS192)
  }
  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[ICAMS::catalog.row.order$SNS192, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

ReadCatalog.SNS1536 <- function(path, ref.genome, region,
                                catalog.type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 1536)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Mutation type", "Mutation Type",
                                   "Mutation.type", "Mutation.Type"))
    stopifnot(names(cos)[2] == "Pentanucleotide")
  }
  names(cos)[1:2] <- c("Mutation type", "Pentanucleotide")
  ref.gt.var       <- cos[["Mutation type"]]
  before.ref.after <- cos[["Pentanucleotide"]]
  var <- substring(ref.gt.var, 3, 3)
  out <- as.matrix(cos[ , -(1 : 2)], drop = FALSE)
  rownames(out) <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SNS1536)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SNS1536, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

ReadCatalog.DNS78 <- function(path, ref.genome, region,
                              catalog.type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, ICAMS::catalog.row.order$DNS78))
  if ( (length(diff1) > 0)
       &&
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
       &&
       (sort(setdiff(ICAMS::catalog.row.order$DNS78, rn) ==
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
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DNS78)
  }
  out <- out[ICAMS::catalog.row.order$DNS78, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

ReadCatalog.DNS144 <- function(path, ref.genome, region,
                               catalog.type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 144)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DNS144)
  }
  out <- out[ICAMS::catalog.row.order$DNS144, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

ReadCatalog.DNS136 <- function(path, ref.genome, region,
                               catalog.type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 136)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Quad", "quad", "QUAD"))
  }
  names(cos)[1] <- "Quad"
  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DNS136)
  }
  out <- out[ICAMS::catalog.row.order$DNS136, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, type))
}

ReadCatalog.ID <- function(path, ref.genome, region, type, strict = TRUE) {
  cos <- data.table::fread(path)
  stopifnot(nrow(cos) == 83)
  cn <- names(cos)
  ex.cn <- c("Type", "Subtype", "Indel_size", "Repeat_MH_size")
  # Repeat_MH_size is the size of repeat OR microhomology (MH)
  if (strict) { for (i in 1 : 4) { stopifnot(cn[i] == ex.cn[i]) } }
  names(cos)[1 : 4] <- ex.cn
  rn <- apply(cos[, 1 : 4], MARGIN = 1, paste, collapse = ":")
  # View(data.frame(mini=rn, good=ICAMS::catalog.row.order$ID))
  out <- as.matrix(cos[ , -(1 : 4)], drop = FALSE)
  rownames(out) <- rn
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$ID)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$ID, , drop = FALSE]
  return(CreateCatalogAttribute(out, ref.genome, region, catalog.type))
}

#' @title Write a catalog to a file.
#'
#' @description This internal function is called by exported functions to do the
#' actual writing of the catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{CreateCatalogAttribute}} for more details.
#'
#' @param path The path of the file to be written.
#'
#' @param num.row The number of rows in the file to be written.
#'
#' @param row.order The row order to be used for writing the file.
#'
#' @param row.header The row header to be used for writing the file.
#'
#' @param strict If TRUE, then stop if additional checks on the input fail.
#'
#' @keywords internal
WriteCat <- function(catalog, path, num.row, row.order, row.header, strict) {
  mut.categories <- rownames(catalog)
  stopifnot(num.row == nrow(catalog))
  if (strict) {
    stopifnot(mut.categories == row.order)
  }
  catalog <- catalog[row.order, , drop = FALSE]
  DT <- as.data.table(catalog)
  fwrite(cbind(row.header, DT), file = path)
}

WriteCatalog.SNS96Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 96, ICAMS::catalog.row.order$SNS96,
           catalog.row.headers.SNS.96, strict)
}

WriteCatalog.SNS192Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 192, ICAMS::catalog.row.order$SNS192,
           catalog.row.headers.SNS.192, strict)
}

WriteCatalog.SNS1536Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 1536, ICAMS::catalog.row.order$SNS1536,
           catalog.row.headers.SNS.1536, strict)
}

WriteCatalog.DNS78Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 78, ICAMS::catalog.row.order$DNS78,
           catalog.row.headers.DNS.78, strict)
}

WriteCatalog.DNS144Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 144, ICAMS::catalog.row.order$DNS144,
           catalog.row.headers.DNS.144, strict)
}

WriteCatalog.DNS136Catalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 136, ICAMS::catalog.row.order$DNS136,
           catalog.row.headers.DNS.136, strict)
}

WriteCatalog.IndelCatalog <- function(catalog, path, strict = TRUE) {
  WriteCat(catalog, path, 83, ICAMS::catalog.row.order$ID,
           catalog.row.headers.ID, strict)
}
