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
#' @section Reading VCF files:
#' \code{\link{ReadListOfStrelkaSNSVCFs}}, which only reads Strelka single nucleotide
#' substitution (SNS) VCFs, not Strelka
#' indel VCFS.\cr
#' \code{\link{ReadListOfStrelkaIDVCFs}}, which only reads Strelka indels (ID) VCFs, not Strelka
#' SNS VCFS.\cr
#' \code{\link{ReadListOfMutectVCFs}}, which reads Mutect VCFs, which contain indels and double
#' nucleotide substitutions (DNSs) as well and SNSs.
#'
#' @section Splitting of in-memory VCFs:
#' \code{\link{SplitListOfStrelkaSNSVCFs}}, which splits Strelka SNS VCFs
#' into SNS and inferred DNS VCFs. \cr
#'  \code{\link{SplitListOfMutectVCFs}}, which separates Mutect VCFs into their
#'  SNS, DNS, and indel components.
#'
#' @section Creating catalogs from VCF files:
#' \code{\link{StrelkaSNSVCFFilesToCatalog}}, which creates 3 SNS catalogs (96,
#' 192, 1536) and 3 DNS catalogs (78, 136, 144) from the Strelka SNS VCFs.\cr
#' \code{\link{StrelkaIDVCFFilesToCatalog}}, which creates ID (indels) catalog
#' from the Strelka ID VCFs.\cr
#' \code{\link{MutectVCFFilesToCatalog}}, which creates 3 SNS catalogs (96,
#' 192, 1536), 3 DNS catalogs (78, 136, 144) and ID (indels) catalog from the
#' Mutect VCFs.
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
#' Functions for collapsing a mutation catalog:
#' \code{\link{CollapseCatalog}}
#'
#' @section Plotting catalogs:
#' Functions for plotting mutation catalog of various samples to a PDF file:
#' \code{\link{PlotCatalogToPdf}}
#'
#' @docType package
#' @name ICAMS
NULL
