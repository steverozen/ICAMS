#' Global data used in ICAMS package
#'
#' Documentations for the global data used in ICAMS package
#'
#' \code{catalog.row.order96} The canonical order of row names in a SNS 96
#' catalog.
#'
#' \code{catalog.row.order192} The canonical order of row names in a SNS 192
#' catalog.
#'
#' \code{catalog.row.order1536} The canonical order of row names in a SNS 1536
#' catalog.
#'
#' \code{catalog.row.order.DNS.78} The canonical order of row names in a DNS 78
#' catalog.
#'
#' \code{catalog.row.order.DNS.144} The canonical order of row names in a DNS
#' 144 catalog.
#'
#' \code{catalog.row.order.QUAD.136} The canonical order of row names in a QUAD
#' 136 catalog.
#'
#' \code{catalog.row.order.ID} The canonical order of row names in a ID
#' (insertion and deletion) catalog. (Please take note that the deletions Repeat
#' Size ranges from 0 to 5+ in the catalog, but for plotting and end user
#' documentation it ranges from 1 to 6+.)
#'
#' \code{ct.96.row.headers} A data frame which contains the row headers
#' information for writing a SNS 96 catalog to disk in PCAWG7 format.
#'
#' \code{ct.192.row.headers} A data frame which contains the row headers
#' information for writing a SNS 192 catalog to disk in PCAWG7 format.
#'
#' \code{ct.1536.row.headers} A data frame which contains the row headers
#' information for writing a SNS 1536 catalog to disk in PCAWG7 format.
#'
#' \code{ct.DNS78.row.headers} A data frame which contains the row headers
#' information fo writing a DNS 78 catalog to disk in PCAWG7 format.
#'
#' \code{ct.DNS144.row.headers} A data frame which contains the row headers
#' information for writing a DNS 144 catalog to disk in PCAWG7 format.
#'
#' \code{ct.QUAD136.row.headers} A data frame which contains the row headers
#' information for writing a QUAD 136 catalog to disk in PCAWG7 format.
#'
#' \code{ct.ID.row.headers} A data frame which contains the row headers
#' information for writing a ID (insertion and deletion) catalog to disk in
#' PCAWG7 format. (Please take note that the deletions Repeat Size ranges from 0
#' to 5+ in the catalog, but for plotting and end user documentation it ranges
#' from 1 to 6+.)
#'
#' \code{abundance.2bp} A matrix containing dinucleotide abundance information
#' for human GRCh37. Its row names indicate 10 different types of 2 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS78}} and
#' \code{\link{CatDNS78ToPdf}}.
#'
#' \code{abundance.3bp} A matrix containing trinucleotide abundance information
#' for human GRCh37. Its row names indicate 32 different types of 3 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCat96}} and
#' \code{\link{Cat96ToPdf}}.
#'
#' \code{abundance.4bp} A matrix containing tetranucleotide abundance
#' information for human GRCh37. Its row names indicate 136 different types of 4
#' base pairs combinations while its column contains the occurrences of each
#' type. It can be used in plotting functions \code{\link{PlotCatQUAD136}} and
#' \code{\link{CatQUAD136ToPdf}}.
#'
#' \code{abundance.5bp} A matrix containing pentanucleotide abundance
#' information for human GRCh37. Its row names indicate 512 different types of 5
#' base pairs combinations while its column contains the occurrences of each
#' type. It can be used in plotting functions \code{\link{PlotCat1536}} and
#' \code{\link{Cat1536ToPdf}}.
#'
#' \code{trans.ranges.GRCh37} A data.table which contains transcript range and
#' strand information for human GRCh37. It is derived from a raw \strong{GFF3}
#' format file, from which only the following four gene types are kept to
#' facilitate transcriptional strand bias analysis: protein_coding,
#' retained_intron, processed_transcript and nonsense_mediated_decay. It
#' contains chromosome name, start, end position, strand information and gene
#' name and is keyed by chrom, chromStart, and chromEnd. It can be used in
#' function \code{\link{StrelkaVCFFilesToCatalog}}.
#'
#' \code{old.trans.ranges.GRCh37} A data.table which contains transcript range
#' and strand information for human GRCh37, which is derived from a raw
#' \strong{BED} format file and is keyed by chrom, chromStart, and chromEnd.
#' This is mostly for testing purpose, may be removed in the future.
#'
#' \code{to.reorder.192.for.plotting} A reordering of row names in a SNS 192
#' catalog for plotting purpose. It is used in plotting functions
#' \code{\link{PlotCat192}} and \code{\link{PlotCat192Strand}}.
#'
#' \code{to.reorder.144.for.plotting} A reordering of row names in a DNS 144
#' catalog for plotting purpose. It is used in plotting function
#' \code{\link{PlotCatDNS144}}.
#'
#' \code{order.for.QUAD136.plotting} An order of tetranucleotides for plotting
#' QUAD 136 catalog. It is used in plotting functions
#' \code{\link{PlotCatQUAD136}} and \code{\link{CatQUAD136ToPdf}}.
#'
#' \code{empty.cats} A list of 6 empty catalogs (SNS 96, SNS 192, SNS 1536, DNS
#' 78, DNS 144, QUAD 136). This is mainly used in the internal functions in
#' ICAMS package.
#'
#' @name data
NULL

#' Canonical Order of a SNS 96 Catalog
#'
#' The canonical order of row names in a SNS 96 catalog.
#' @format A string of characters with length 96.
"catalog.row.order96"

