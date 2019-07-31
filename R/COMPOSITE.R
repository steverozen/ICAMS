#' @export
ReadCatalog.COMPOSITECatalog <-
  function(file, ref.genome = NULL, region = "unknown", 
           catalog.type = "counts", strict = TRUE)
    
  {
  retval <- read.csv(file, header = T, row.names = 1)
  return(as.catalog(retval,
                    ref.genome = ref.genome,
                    region = region,
                    catalog.type = catalog.type))
}


#' @export
WriteCatalog.COMPOSITECatalog <-
  function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file = file, 
           1697,
           ICAMS::catalog.row.order[["COMPOSITE"]],
           # For the on-disk representation of COMPOSITE
           # catalogs we just use the rownames from the
           # in-memory representation.
           rownames(catalog), 
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
  cat96 <- Collapse192CatalogTo96(cat1536)  
  all.0 <- which(colSums(cat96) == 0)
  if (length(all.0) > 0 ) {
    cat96[ , all.0] <- 1
    cn <- colnames(cat96)
    cn[all.0] <- paste(cn[all.0], "WARNING all 0")
    colnames(cat96) <- cn
  }
  PlotCatalogToPdf(catalog = cat96/sum(cat96), filename = name, type = type)
}

#' Plot the a SignatureAnalyzer COMPOSITE signature or catalog into separate pdfs
#'
#' @param catalog Catalog or signature matrix
#'
#' @param filename.header Contain path and the beginning part of the file name.
#' The name of the pdf files will be:
#' \code{filename.header}.SNS.96.pdf
#' \code{filename.header}.SNS.1536.pdf
#' \code{filename.header}.DNS.78.pdf
#' \code{filename.header}.ID.83.pdf
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @param id A vector containing the identifiers of the samples
#' or signatures in \code{catalog}.
#'
#' @keywords internl
TestPlotCatCOMPOSITE <-
  function(catalog, filename.header, type, id = colnames(catalog)) {
  
  ## Read in COMPOSITE catalogue
  test.COMPOSITE.sigs <- ReadCatalog(catalog)
  
  ## Check
  stopifnot(nrow(test.COMPOSITE.sigs) == 1697)
  # TODO WUYang: check whether the base context is in correct order
  
  ## Subsetting COMPOSITE catalogue
  test.SNS1536.sigs <- test.COMPOSITE.sigs[1:1536,]
  test.DNS78.sigs   <- test.COMPOSITE.sigs[1537:1614,]
  test.ID83.sigs    <- test.COMPOSITE.sigs[1615:1697,]
  
  ## Plot using ICAMS embedded plotting function
  
  ICAMS::PlotCatalogToPdf(test.SNS1536.sigs,
                          filename = paste0(filename.header,".SNS.1536.pdf"),
                          type = type,
                          id = id)
  ICAMS::PlotCatalogToPdf(test.DNS78.sigs,
                          filename = paste0(filename.header,".DNS.78.pdf"),
                          type = type,
                          id = id)
  ICAMS::PlotCatalogToPdf(test.ID83.sigs,
                          filename = paste0(filename.header,".ID.83.pdf"),
                          type = type,
                          id = id)
  }
