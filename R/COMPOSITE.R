#' @export
WriteCatalog.COMPOSITECatalog <-
  function(catalog, file, strict = TRUE) {
    WriteCat(catalog, file = file, 
             1697,
             ICAMS::catalog.row.order[["COMPOSITE"]],
             # rownames(catalog), 
             catalog.row.headers$COMPOSITE,
             strict = strict)
  }

#' @keywords internal
SplitCatCOMPOSITE <- function(catalog) {

  # Split COMPOSITE catalog to 3 catalogs
  # (1 SBS1536, 1 DBS78, 1 Indel83)
  # use function in ReadWriteCatalogs.R
  catList <- list("SBS1536" = catalog[1:1536,],
                  "DBS78"   = catalog[1537:1614,],
                  "ID"      = catalog[1615:1697,])

  return(catList)
}

#' Plot the SBS96 part of a SignatureAnalyzer COMPOSITE signature or catalog
#'
#' @param catalog Catalog or signature matrix
#'
#' @param name Name of file to print to.
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @keywords internal

Plot96PartOfCompositeToPDF <- function(catalog, name, type = "density") {
  cat1536 <- catalog[1:1536, ]
  stop("Not supported, need to change class to SBS1536Catalog here")
  cat96 <- Collapse192CatalogTo96(cat1536)  
  all.0 <- which(colSums(cat96) == 0)
  if (length(all.0) > 0 ) {
    cat96[ , all.0] <- 1
    cn <- colnames(cat96)
    cn[all.0] <- paste(cn[all.0], "WARNING all 0")
    colnames(cat96) <- cn
  }
  PlotCatalogToPdf(catalog = cat96, file = name)
}

#' Plot the a SignatureAnalyzer COMPOSITE signature or catalog into separate pdfs
#'
#' @param catalog Catalog or signature matrix
#'
#' @param filename.header Contain path and the beginning part of the file name.
#' The name of the pdf files will be:
#' \emph{filename.header}\code{.SBS.96.pdf}
#' \emph{filename.header}\code{.SBS.1536.pdf}
#' \emph{filename.header}\code{.DBS.78.pdf}
#' \emph{filename.header}\code{.ID.83.pdf}
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @param id A vector containing the identifiers of the samples
#' or signatures in \code{catalog}.
#'
#' @keywords internal
TestPlotCatCOMPOSITE <-
  function(catalog, filename.header, type, id = colnames(catalog)) {
  
  ## Read in COMPOSITE catalogue
  # test.COMPOSITE.sigs <- ReadCatalog(catalog)
  
  # Check
  # stopifnot(nrow(test.COMPOSITE.sigs) == 1697)
  spectra <- matrix(stats::rnbinom(n = 1697 * 2, mu = 1000, size = 3), ncol = 2)
  colnames(spectra) <- c("s1", "s1")
  spectra <- as.catalog(spectra, infer.rownames = TRUE)

  # Subset COMPOSITE catalogue
  test.SBS1536 <- spectra[1:1536,]
  test.DBS78   <- spectra[1537:1614,]
  test.ID      <- spectra[1615:1697,]
  stop("Not supported; need to change classes of the row subsets")
  
  ## Plot using ICAMS embedded plotting function
  
  ICAMS::PlotCatalogToPdf(test.SBS1536,
                          file = paste0(filename.header,".SBS.1536.pdf"))
  ICAMS::PlotCatalogToPdf(test.DBS78,
                          file = paste0(filename.header,".DBS.78.pdf"))
  ICAMS::PlotCatalogToPdf(test.ID,
                          file = paste0(filename.header,".ID.83.pdf"))
  }
