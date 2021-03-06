% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_catalog_formats.R
\name{ConvertCatalogToSigProfilerFormat}
\alias{ConvertCatalogToSigProfilerFormat}
\title{Covert an ICAMS Catalog to SigProfiler format}
\usage{
ConvertCatalogToSigProfilerFormat(input.catalog, file, sep = "\\t")
}
\arguments{
\item{input.catalog}{Either a character string, in which case this is the
path to a file containing a catalog in \code{\link[ICAMS]{ICAMS}}
format, or an in-memory \code{\link[ICAMS]{ICAMS}} catalog.}

\item{file}{The path of the file to be written.}

\item{sep}{Separator to use in the output file.}
}
\description{
Covert an ICAMS Catalog to SigProfiler format
}
\note{
This function can only transform SBS96, DBS78 and ID ICAMS catalog
to SigProfiler format.
}
\keyword{internal}
