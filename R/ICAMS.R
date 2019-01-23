#' ICAMS: In-depth Characterization and Analysis of Mutational Signatures
#'
#' This package has functions to read in VCF files from Strelka and GATK,
#' create SNS, DNS, ID catalogs and do different types of plotting.
#'
#' @section Reading catalogs:
#' Functions for reading a catalog in PCAWG7 format from path:
#' \code{\link{ReadCatalog}}
#'
#' @section Writing catalogs:
#' Functions for writting a mutation catalog to a file on disk:
#' \code{\link{WriteCatalog}}
#'
#' @section Collapsing catalogs:
#' Functions for collapsing a mutation catalog to a canonical one:
#' \code{\link{CollapseCatalog}}
#'
#' @section Plotting catalogs:
#' Functions for plotting the mutation catalog of one sample:
#' \code{\link{PlotCatalog}}
#'
#' Functions for plotting mutation catalog of different samples to a PDF file:
#' \code{\link{CatalogToPdf}}
#'
#' @docType package
#' @name ICAMS
NULL
