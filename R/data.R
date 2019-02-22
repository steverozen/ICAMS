#' Canonical order of row names in a catalog
#' @format A string of characters indicating the canonical order of row names in
#'   a catalog.
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   ranges from 0 to 5+, but for plotting and end user documentation it ranges
#'   from 1 to 6+.
#' @name CatalogRowOrder
NULL

#' Row headers information for writing a catalog to disk in PCAWG7 format
#' @format A data frame which contains the row headers information for writing
#'   a catalog to disk in PCAWG7 format.
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   ranges from 0 to 5+, but for plotting and end user documentation it ranges
#'   from 1 to 6+.
#' @name CatalogRowHeaders
NULL

#' Nucleotide abundance file
#'
#' Nucleotide abundance information for a particular organism
#'
#' \code{abundance.2bp.genome.GRCh37}, \code{abundance.2bp.exome.GRCh37} A
#' matrix containing dinucleotide abundance information for \strong{Human}
#' GRCh37. Its row names indicate 10 different types of 2 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS78}} and
#' \code{\link{PlotCatDNS78ToPdf}}.
#'
#' \code{abundance.2bp.genome.GRCh38}, \code{abundance.2bp.exome.GRCh38} A
#' matrix containing dinucleotide abundance information for \strong{Human}
#' GRCh38. Its row names indicate 10 different types of 2 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS78}} and
#' \code{\link{PlotCatDNS78ToPdf}}.
#'
#' \code{abundance.2bp.genome.GRCm38}, \code{abundance.2bp.exome.GRCm38} A
#' matrix containing dinucleotide abundance information for \strong{Mouse}
#' GRCm38. Its row names indicate 10 different types of 2 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS78}} and
#' \code{\link{PlotCatDNS78ToPdf}}.
#'
#' \code{abundance.3bp.genome.GRCh37}, \code{abundance.3bp.exome.GRCh37} A
#' matrix containing trinucleotide abundance information for \strong{Human}
#' GRCh37. Its row names indicate 32 different types of 3 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS96}} and
#' \code{\link{PlotCatSNS96ToPdf}}.
#'
#' \code{abundance.3bp.genome.GRCh38}, \code{abundance.3bp.exome.GRCh38} A
#' matrix containing trinucleotide abundance information for \strong{Human}
#' GRCh38. Its row names indicate 32 different types of 3 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS96}} and
#' \code{\link{PlotCatSNS96ToPdf}}.
#'
#' \code{abundance.3bp.genome.GRCm37}, \code{abundance.3bp.exome.GRCm37} A
#' matrix containing trinucleotide abundance information for \strong{Mouse}
#' GRCm37. Its row names indicate 32 different types of 3 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS96}} and
#' \code{\link{PlotCatSNS96ToPdf}}.
#'
#' \code{abundance.4bp.genome.GRCh37}, \code{abundance.4bp.exome.GRCh37} A
#' matrix containing tetranucleotide abundance information for \strong{Human}
#' GRCh37. Its row names indicate 136 different types of 4 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS136}} and
#' \code{\link{PlotCatDNS136ToPdf}}.
#'
#' \code{abundance.4bp.genome.GRCh38}, \code{abundance.4bp.exome.GRCh38} A
#' matrix containing tetranucleotide abundance information for \strong{Human}
#' GRCh38. Its row names indicate 136 different types of 4 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS136}} and
#' \code{\link{PlotCatDNS136ToPdf}}.
#'
#' \code{abundance.4bp.genome.GRCm37}, \code{abundance.4bp.exome.GRCm37} A
#' matrix containing tetranucleotide abundance information for \strong{Mouse}
#' GRCm37. Its row names indicate 136 different types of 4 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatDNS136}} and
#' \code{\link{PlotCatDNS136ToPdf}}.
#'
#' \code{abundance.5bp.genome.GRCh37}, \code{abundance.5bp.exome.GRCh37} A
#' matrix containing pentanucleotide abundance information for \strong{Human}
#' GRCh37. Its row names indicate 512 different types of 5 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS1536}} and
#' \code{\link{PlotCatSNS1536ToPdf}}.
#'
#' \code{abundance.5bp.genome.GRCh38}, \code{abundance.5bp.exome.GRCh38} A
#' matrix containing pentanucleotide abundance information for \strong{Human}
#' GRCh38. Its row names indicate 512 different types of 5 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS1536}} and
#' \code{\link{PlotCatSNS1536ToPdf}}.
#'
#' \code{abundance.5bp.genome.GRCm37}, \code{abundance.5bp.exome.GRCm37} A
#' matrix containing pentanucleotide abundance information for \strong{Mouse}
#' GRCm37. Its row names indicate 512 different types of 5 base pairs
#' combinations while its column contains the occurrences of each type. It can
#' be used in plotting functions \code{\link{PlotCatSNS1536}} and
#' \code{\link{PlotCatSNS1536ToPdf}}.
#' @format A single-column
#' matrix containing the counts of particular sequences in
#' a genome or part of a genome. This include
#' 2-mers, 3-mers, 4-mers, 5-mers, stranded or strand-agnostic,
#' and genome-wide, in-transcript, or in-exome, for different
#' reference genome versions and for different organisms. The names should
#' be self explanatory.
#' @note In the ID (insertion and deletion) catalog, deletion repeat size
#'   ranges from 0 to 5+, but for plotting and end user documentation it ranges
#'   from 1 to 6+.
#' @name Abundance
NULL

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
#' \code{old.trans.ranges.GRCh37} A data.table which contains transcript range
#' and strand information for \strong{Human} GRCh37, which is derived from a raw
#' \strong{BED} format file and is keyed by chrom, chromStart, and chromEnd.
#' This is mostly for testing purpose, may be removed in the future.
#' @format A data.table which contains transcript range and strand information
#'   for a particular organism.
#' @name TranscriptRanges
NULL

#' @rdname CatalogRowOrder
"catalog.row.order.SNS.96"

#' @rdname CatalogRowOrder
"catalog.row.order.SNS.192"

#' @rdname CatalogRowOrder
"catalog.row.order.SNS.1536"

#' @rdname CatalogRowOrder
"catalog.row.order.DNS.78"

#' @rdname CatalogRowOrder
"catalog.row.order.DNS.144"

#' @rdname CatalogRowOrder
"catalog.row.order.DNS.136"

#' @rdname CatalogRowOrder
"catalog.row.order.ID"

#' @rdname CatalogRowHeaders
"catalog.row.headers.SNS.96"

#' @rdname CatalogRowHeaders
"catalog.row.headers.SNS.192"

#' @rdname CatalogRowHeaders
"catalog.row.headers.SNS.1536"

#' @rdname CatalogRowHeaders
"catalog.row.headers.DNS.78"

#' @rdname CatalogRowHeaders
"catalog.row.headers.DNS.144"

#' @rdname CatalogRowHeaders
"catalog.row.headers.DNS.136"

#' @rdname CatalogRowHeaders
"catalog.row.headers.ID"

#' @rdname Abundance
"abundance.2bp.exome.GRCh37"

#' @rdname Abundance
"abundance.2bp.genome.GRCh37"

#' @rdname Abundance
"abundance.3bp.exome.GRCh37"

#' @rdname Abundance
"abundance.3bp.genome.GRCh37"

#' @rdname Abundance
"abundance.4bp.exome.GRCh37"

#' @rdname Abundance
"abundance.4bp.genome.GRCh37"

#' @rdname Abundance
"abundance.5bp.exome.GRCh37"

#' @rdname Abundance
"abundance.5bp.genome.GRCh37"

#' @rdname Abundance
"abundance.2bp.exome.GRCh38"

#' @rdname Abundance
"abundance.2bp.genome.GRCh38"

#' @rdname Abundance
"abundance.3bp.exome.GRCh38"

#' @rdname Abundance
"abundance.3bp.genome.GRCh38"

#' @rdname Abundance
"abundance.4bp.exome.GRCh38"

#' @rdname Abundance
"abundance.4bp.genome.GRCh38"

#' @rdname Abundance
"abundance.5bp.exome.GRCh38"

#' @rdname Abundance
"abundance.5bp.genome.GRCh38"

#' @rdname Abundance
"abundance.2bp.exome.GRCm38"

#' @rdname Abundance
"abundance.2bp.genome.GRCm38"

#' @rdname Abundance
"abundance.3bp.exome.GRCm38"

#' @rdname Abundance
"abundance.3bp.genome.GRCm38"

#' @rdname Abundance
"abundance.4bp.exome.GRCm38"

#' @rdname Abundance
"abundance.4bp.genome.GRCm38"

#' @rdname Abundance
"abundance.5bp.exome.GRCm38"

#' @rdname Abundance
"abundance.5bp.genome.GRCm38"

#' @rdname TranscriptRanges
"trans.ranges.GRCh37"

#' @rdname TranscriptRanges
"old.trans.ranges.GRCh37"
