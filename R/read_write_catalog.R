#' Read catalog.
#'
#' Read a catalog in standardized format from path.
#'
#' See also \code{\link{WriteCatalog}}
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @param strict If TRUE, do additional checks on the input, and stop if the
#'   checks fail.
#'
#' @return A catalog as an S3 object; see \code{\link{as.catalog}}.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file, ref.genome = "hg19", 
#'                         region = "genome",
#'                         catalog.type = "counts")
ReadCatalog <- function(file, ref.genome, region, catalog.type, strict = TRUE) {
  StopIfRegionIllegal(region)
  StopIfCatalogTypeIllegal(catalog.type)
  class.of.catalog <- InferClassOfCatalogForRead(file) #
  UseMethod(generic = "ReadCatalog", object = class.of.catalog)
}

#' Write a catalog
#'
#' Write a catalog to a file.
#'
#' See also \code{\link{ReadCatalog}}.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}};
#' see also \code{\link{as.catalog}}.
#'
#' @param file The path to the file to be created.
#'
#' @param strict If TRUE, do additional checks on the input, and stop if the
#'   checks fail.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file, ref.genome = "hg19", 
#'                         region = "genome",
#'                         catalog.type = "counts")
#' WriteCatalog(catSBS96, file = file.path(tempdir(), "catSBS96.csv"))
WriteCatalog <- function(catalog, file, strict = TRUE) {
  UseMethod(generic = "WriteCatalog")
}

#' @export
ReadCatalog.SBS96Catalog <- function(file, ref.genome, region,
                                     catalog.type, strict = TRUE) {
  cos <- data.table::fread(file)
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
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SBS96)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.SBS192Catalog <- function(file, ref.genome, region,
                                      catalog.type, strict = TRUE) {
  
  StopIfTranscribedRegionIllegal(region)
  
  cos <- data.table::fread(file)
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
  ## SBS is on the transcribed (which is the *antisense*) strand.
  transcribed.strand.pos <- which(cos[, 1] == 'T')

  before.ref.after[transcribed.strand.pos] <-
    revc(before.ref.after[transcribed.strand.pos])

  var <- substring(ref.gt.var, 3, 3)
  var[transcribed.strand.pos] <- revc(var[transcribed.strand.pos])

  tmp <- paste0(before.ref.after, var)
  if (strict) {
    stopifnot(tmp == ICAMS::catalog.row.order$SBS192)
  }
  out <- cos[, -(1 : 3), drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- tmp
  out <- out[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.SBS1536Catalog <- function(file, ref.genome, region,
                                       catalog.type, strict = TRUE) {
  cos <- data.table::fread(file)
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
    stopifnot(rownames(out) == ICAMS::catalog.row.order$SBS1536)
  }
  if (ncol(out) == 1) colnames(out) <- colnames(cos)[3]
  out <- out[ICAMS::catalog.row.order$SBS1536, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.DBS78Catalog <- function(file, ref.genome, region,
                                     catalog.type, strict = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 78)
  if (strict) {
    stopifnot(names(cos)[1 : 2] == c("Ref", "Var"))
  }
  names(cos)[1 : 2] <- c("Ref", "Var")
  out <- cos[, -(1 : 2), drop = FALSE]
  out <- as.matrix(out)
  rn <- paste0(cos$Ref, cos$Var)
  diff1 <- sort(setdiff(rn, ICAMS::catalog.row.order$DBS78))
  if ( (length(diff1) > 0)
       &&
       (diff1 == c("CGAA", "CGAC", "CGGA", "TAAC", "TAAG", "TACC"))
       &&
       (sort(setdiff(ICAMS::catalog.row.order$DBS78, rn) ==
             c("CGGT", "CGTC", "CGTT", "TACT", "TAGG", "TAGT")))
  ) {
    cat("using temporary hack for old DBS canonicalization\n")
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
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS78)
  }
  out <- out[ICAMS::catalog.row.order$DBS78, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.DBS144Catalog <- function(file, ref.genome, region,
                                      catalog.type, strict = TRUE) {
  
  StopIfTranscribedRegionIllegal(region)
  
  cos <- data.table::fread(file)
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
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS144)
  }
  out <- out[ICAMS::catalog.row.order$DBS144, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.DBS136Catalog <- function(file, ref.genome, region,
                                      catalog.type, strict = TRUE) {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 136)
  if (strict) {
    stopifnot(names(cos)[1] %in% c("Quad", "quad", "QUAD"))
  }
  names(cos)[1] <- "Quad"
  out <- cos[, -1, drop = FALSE]
  out <- as.matrix(out)
  rownames(out) <- cos$Quad
  if (strict) {
    stopifnot(rownames(out) == ICAMS::catalog.row.order$DBS136)
  }
  out <- out[ICAMS::catalog.row.order$DBS136, , drop = FALSE]
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @export
ReadCatalog.IndelCatalog <- function(file, ref.genome, region,
                                     catalog.type, strict = TRUE) {
  cos <- data.table::fread(file)
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
  return(as.catalog(out, ref.genome, region, catalog.type))
}

#' @title Write a catalog to a file.
#'
#' @description This internal function is called by exported functions to do the
#' actual writing of the catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param file The path of the file to be written.
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
WriteCat <- function(catalog, file, num.row, row.order, row.header, strict) {
  mut.categories <- rownames(catalog)
  stopifnot(num.row == nrow(catalog))
  if (strict) {
    stopifnot(mut.categories == row.order)
  }
  catalog <- catalog[row.order, , drop = FALSE]
  DT <- as.data.table(catalog)
  fwrite(cbind(row.header, DT), file = file)
}

#' @export
WriteCatalog.SBS96Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 96, ICAMS::catalog.row.order$SBS96,
           catalog.row.headers.SBS.96, strict)
}

#' @export
WriteCatalog.SBS192Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 192, ICAMS::catalog.row.order$SBS192,
           catalog.row.headers.SBS.192, strict)
}

#' @export
WriteCatalog.SBS1536Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 1536, ICAMS::catalog.row.order$SBS1536,
           catalog.row.headers.SBS.1536, strict)
}

#' @export
WriteCatalog.DBS78Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 78, ICAMS::catalog.row.order$DBS78,
           catalog.row.headers.DBS.78, strict)
}

#' @export
WriteCatalog.DBS144Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 144, ICAMS::catalog.row.order$DBS144,
           catalog.row.headers.DBS.144, strict)
}

#' @export
WriteCatalog.DBS136Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 136, ICAMS::catalog.row.order$DBS136,
           catalog.row.headers.DBS.136, strict)
}

#' @export
WriteCatalog.IndelCatalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 83, ICAMS::catalog.row.order$ID,
           catalog.row.headers.ID, strict)
}
