#' Read catalog
#'
#' Read a catalog in standardized format from path.
#'
#' See also \code{\link{WriteCatalog}}
#'
#' @param file Path to a catalog on disk in a standardized format.
#' The recognized formats are:
#' 
#' * ICAMS formatted SBS96, SBS192, SBS1536, DBS78, DBS136, DBS144, ID, ID166 
#'   (see \code{\link{CatalogRowOrder}}).
#'   
#' * SigProfiler-formatted SBS96, DBS78 and ID83 catalogs;
#'   see \url{https://github.com/AlexandrovLab/SigProfilerExtractor}.
#' 
#' * COSMIC-formatted SBS96, SBS192 (a.k.a. TSB192), 
#'   DBS78, ID83 catalogs; 
#'   see \url{https://cancer.sanger.ac.uk/signatures/}.
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
#' @param strict Ignored and deprecated.
#'   
#' @param stop.on.error If TRUE, call \code{stop} on error; otherwise
#'   return a 1-column matrix of NA's with the attribute "error"
#'   containing error information. The number of rows may not
#'   be the correct number for the expected catalog type.
#'
#' @return A catalog as an S3 object; see \code{\link{as.catalog}}.
#'
#' @note In ID (small insertions and deletions) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @inheritSection MutectVCFFilesToCatalog Comments
#' 
#' @export
#' 
#' @md
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file)
#' 

ReadCatalog <- function(file, 
                        ref.genome    = NULL, 
                        region        = "unknown", 
                        catalog.type  = "counts", 
                        strict        = NULL,
                        stop.on.error = TRUE) {
  tryCatch(
    return(ReadCatalogInternal(
      file = file,
      ref.genome = ref.genome,
      region = region,
      catalog.type = catalog.type)),
    error = function(err) {
      return(
        ReadCatalogErrReturn(
          err.info      = err,
          nrow          = 1, # We do not know what type of catalog
          stop.on.error = stop.on.error # FALSE for ICAMS-server
          ))})
}

#' Internal read catalog function to be wrapped in a tryCatch
#' 
#' @inheritParams ReadCatalog
#' 
#' @importFrom stats na.omit
#' 
#' @keywords internal
ReadCatalogInternal <- function(file, ref.genome = NULL, region = "unknown", 
                                catalog.type = "counts") {
  StopIfRegionIllegal(region)
  StopIfCatalogTypeIllegal(catalog.type)
  ## The external catalog file is imported 
  ## as a matrix object
  ## and a catalog object in this step.
  dt <- data.table::fread(file)
  
  # In some rare cases, there may be all NA in some columns in dt.
  # So we remove the columns which have all NA in dt
  dt <- dt[, which(unlist(lapply(dt, function(x)!all(is.na(x))))), with = FALSE]
  
  # In some rare cases, there may be NA in dt, then the number of rows will not
  # be accurate to infer catalog type. So we remove the rows which have NA in dt
  dt <- stats::na.omit(dt)
  catalog <- InferCatalogInfo(dt)
  attr(catalog, "ref.genome") <- ref.genome
  cat2 <- as.catalog(catalog, ref.genome = ref.genome,
                     region = region, catalog.type = catalog.type)
  
  return(cat2)

}

#' Get error message and either stop or create a null error output for read catalog
#' 
#' @param err.info The information passed to the \code{tryCatch} \code{error}
#'   function argument.
#'   
#' @param nrow The number of rows to put in the 1-column NA return matrix.
#' 
#' @param stop.on.error If \code{TRUE} then call \code{stop()}.
#' 
#' @param do.message  If \code{TRUE} then \code{message} the error information.
#' 
#' @keywords internal

ReadCatalogErrReturn <- 
  function(err.info, nrow, stop.on.error = TRUE, do.message = TRUE) {
  if (!is.null(err.info$message)) err.info <- err.info$message
  if (stop.on.error) stop(err.info)
  if (do.message) message(err.info)
  null.out <- matrix(NA, ncol = 1, nrow = nrow)
  attr(null.out, "error") <- err.info
  return(null.out)
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
#' @note In ID (small insertions and deletions) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file)
#' WriteCatalog(catSBS96, file = file.path(tempdir(), "catSBS96.csv"))
WriteCatalog <- function(catalog, file, strict = TRUE) {
  UseMethod(generic = "WriteCatalog")
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
WriteCat <- function(catalog, file, num.row, row.order, row.header, strict,
                     sep = ",") {
  mut.categories <- rownames(catalog)
  stopifnot(num.row == nrow(catalog))
  if (strict) {
    stopifnot(mut.categories == row.order)
  }
  catalog <- catalog[row.order, , drop = FALSE]
  DT <- as.data.table(catalog)
  fwrite(cbind(row.header, DT), file = file, sep = sep)
}

#' @export
WriteCatalog.SBS96Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 96, ICAMS::catalog.row.order$SBS96,
           catalog.row.headers$SBS96, strict)
}

#' @keywords internal
WriteSBS96CatalogAsTsv <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 96, ICAMS::catalog.row.order$SBS96,
           catalog.row.headers.SBS.96.v1, strict, sep = "\t")
}

#' @export
WriteCatalog.SBS192Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 192, ICAMS::catalog.row.order$SBS192,
           catalog.row.headers$SBS192, strict)
}

#' @export
WriteCatalog.SBS1536Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 1536, ICAMS::catalog.row.order$SBS1536,
           catalog.row.headers$SBS1536, strict)
}

#' @export
WriteCatalog.DBS78Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 78, ICAMS::catalog.row.order$DBS78,
           catalog.row.headers$DBS78, strict)
}

#' @export
WriteCatalog.DBS144Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 144, ICAMS::catalog.row.order$DBS144,
           catalog.row.headers$DBS144, strict)
}

#' @export
WriteCatalog.DBS136Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 136, ICAMS::catalog.row.order$DBS136,
           catalog.row.headers$DBS136, strict)
}

#' @export
WriteCatalog.IndelCatalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 83, ICAMS::catalog.row.order$ID,
           catalog.row.headers$ID, strict)
}

#' @export
WriteCatalog.ID166Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 166, ICAMS::catalog.row.order$ID166,
           catalog.row.headers$ID166, strict)
}
