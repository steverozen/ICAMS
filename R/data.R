#' Standard order of row names in a catalog.
#'
#' This data is designed for those
#' who need to create their own catalogs from formats not
#' supported by this package. The rownames denote the mutation
#' types.  For example, for SBS96 catalogs, the rowname
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
#' Transcript ranges and strand information for a particular reference genome.
#'
#' \code{trans.ranges.GRCh37} A data.table which contains transcript range and
#' strand information for \strong{Human} GRCh37. It is derived from a raw
#' \strong{GFF3} format file
#' (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gff3.gz),
#' only genes that are associated with a CCDS ID are kept for
#' transcriptional strand bias analysis. Needed for
#' \code{\link{StrelkaSBSVCFFilesToCatalog}},
#' \code{\link{MutectVCFFilesToCatalog}},  \code{\link{VCFsToSBSCatalogs}} and
#' \code{\link{VCFsToDBSCatalogs}}.
#'
#' \code{trans.ranges.GRCh38} A data.table which contains transcript range and
#' strand information for \strong{Human} GRCh38. It is derived from a raw
#' \strong{GFF3} format file
#' (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gff3.gz),
#' only genes that are associated with a CCDS ID are kept for
#' transcriptional strand bias analysis. Needed for
#' \code{\link{StrelkaSBSVCFFilesToCatalog}},
#' \code{\link{MutectVCFFilesToCatalog}},  \code{\link{VCFsToSBSCatalogs}} and
#' \code{\link{VCFsToDBSCatalogs}}.
#'
#' @format A data.table which contains transcript range and strand information
#'   for a particular reference genome. It contains chromosome name, start, end
#'   position, strand information and gene name and is keyed by chrom,
#'   chromStart, and chromEnd. It uses one-based coordinate system.
#'
#' @name TranscriptRanges
NULL

#' @rdname TranscriptRanges
"trans.ranges.GRCh37"

#' @rdname TranscriptRanges
"trans.ranges.GRCh38"

# Quiets concerns of R CMD check about no visible binding for global variable
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("POS2", "POS", "bothstrand", "strand", ".", "CHROM",
                           "ALT", "count", "rn", "occurrences", "type", "strand",
                           "bothstrand", "chrom", "exome.start", "exome.end",
                           "count", "REF", "seq.21bases", "N", "pyr.mut", "nrn",
                           "mutation", "LOW", "ID", "REF.x", "REF.y", "ALT.x",
                           "ALT.y", "ref2alt", "minus1bs", "minus2bs", "plus1bs",
                           "plus2bs", "POS.plus.one", "HIGH", "POS.y", "VAF.x",
                           "VAF.y", "delete.flag", "trans.ranges.GRCh37", "cols"))
}
