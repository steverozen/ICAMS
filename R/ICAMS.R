#' ICAMS: In-depth Characterization and Analysis of Mutational Signatures
#'
#' This package has functions to read in VCF files from Strelka and Mutect (in
#' the Broad GATK package), create, read, and write single nucleotide
#' substitutions (SNS), double nucleotide substitutions (DNS), insertions and
#' deletions (ID) catalogs and do different types of plotting.
#'
#' This alpha version only works with VCFs for human GRCh37, but will work for
#' arbitrary \strong{human} catalogs (assuming no major change in "opportunities"
#' between GRCh37 and GRCh38).
#'
#' @section Reading and splitting VCF files: \enumerate{
#'
#' \item \code{\link{ReadAndSplitStrelkaSNSVCFs}} Read and split Strelka single
#' nucleotide substitution (SNS) VCFs (not Strelka indel VCFS).
#'
#' \item \code{\link{ReadStrelkaIDVCFs}} Read Strelka indel (ID) VCFs (not Strelka
#' SNS VCFS).
#'
#' \item \code{\link{ReadAndSplitMutectVCFs}} Read and split Mutect VCFs, which
#' contain indels and double nucleotide substitutions (DNSs) as well and SNSs.
#'}
#'
#' @section Creating catalogs from VCF files:
#' \enumerate{
#' \item \code{\link{StrelkaSNSVCFFilesToCatalog}}, which creates 3 SNS catalogs (96,
#' 192, 1536) and 3 DNS catalogs (78, 136, 144) from the Strelka SNS VCFs.
#'
#' \item \code{\link{StrelkaIDVCFFilesToCatalog}}, which creates ID (indels) catalog
#' from the Strelka ID VCFs.
#'
#' \item \code{\link{MutectVCFFilesToCatalog}}, which creates 3 SNS catalogs (96,
#' 192, 1536), 3 DNS catalogs (78, 136, 144) and ID (indels) catalog from the
#' Mutect VCFs.
#' }
#'
#' @section Reading catalogs:
#' Functions for reading files that contain mutational
#' spectrum catalogs in standardized format. These
#' also work for reading mutational signature profiles.
#' \code{\link{ReadCatalog}}
#'
#' @section Writing catalogs:
#' Functions for writing a mutational spectrum catalog to a file on disk.
#' These also work for writing mutational signature profiles.
#' \code{\link{WriteCatalog}}
#'
#' @section Transforming catalogs:
#' Functions for transforming count spectra from a particular organism region to
#' an inferred count spectra based on the target nucleotide abundance.
#' \code{\link{TransformSpectra}}
#'
#' @section Collapsing catalogs:
#' Functions for collapsing a mutation catalog.
#' \code{\link{CollapseCatalog}}
#'
#' @section Plotting catalogs:
#' Functions for plotting mutation spectrum catalogs
#' to a PDF file. These also work for plotting
#' mutational signature profiles.
#' \code{\link{PlotCatalogToPdf}}
#'
#' @section Exported data: \enumerate{
#'
#' \item \code{\link{CatalogRowOrder}} Canonical order of row names in a catalog.
#'
#'\item \code{\link{TranscriptRanges}} Transcript ranges and strand information
#'for a particular organism. }
#' @docType package
#' @name ICAMS
NULL
