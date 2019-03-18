#' Standard order of row names in a catalog.
#'
#' This data is designed for those
#' who need to create their own catalogs from formats not
#' supported by this package. The rownames denote the mutation
#' types.  For example, for SNS96 catalogs, the rowname
#'  AGAT represents a mutation from AGA > ATA.
#'
#' @format A list of character vectors indicating the standard
#'   orders of row names in catalogs.
#'
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   is in the range from 0 to 5+, but for plotting
#'   and end user documentation it ranges
#'   from 1 to 6+.
#'
#' @name CatalogRowOrder
"catalog.row.order"

#' Transcript ranges data
#'
#' Transcript ranges and strand information for a particular organism
#'
#' \code{trans.ranges.GRCh37} A data.table which contains transcript range and
#' strand information for \strong{Human} GRCh37. It is derived from a raw \strong{GFF3}
#' format file, from which only the following four gene types are kept to
#' facilitate transcriptional strand bias analysis: protein_coding,
#' retained_intron, processed_transcript and nonsense_mediated_decay. It
#' contains chromosome name, start, end position, strand information and gene
#' name and is keyed by chrom, chromStart, and chromEnd. It can be used in
#' function \code{\link{StrelkaSNSVCFFilesToCatalog}}.
#'
#' \code{trans.ranges.GRCh38} A data.table which contains transcript range and
#' strand information for \strong{Human} GRCh38. It is derived from a raw \strong{GFF3}
#' format file, from which only the following four gene types are kept to
#' facilitate transcriptional strand bias analysis: protein_coding,
#' retained_intron, processed_transcript and nonsense_mediated_decay. It
#' contains chromosome name, start, end position, strand information and gene
#' name and is keyed by chrom, chromStart, and chromEnd. It can be used in
#' function \code{\link{StrelkaSNSVCFFilesToCatalog}}.
#'
#' @format A data.table which contains transcript range and strand information
#'   for a particular organism.
#'
#' @name TranscriptRanges
NULL

#' @rdname CatalogRowOrder
"catalog.row.order"

#' @rdname TranscriptRanges
"trans.ranges.GRCh37"

#' @rdname TranscriptRanges
"trans.ranges.GRCh38"
