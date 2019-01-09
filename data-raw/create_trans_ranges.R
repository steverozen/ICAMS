ReadTranscriptRanges <- function(path) {
  # Read transcript ranges and strands from a bed format file.
  # Use this one for the new, cut down gff3 file (2018 11 24)
  #
  # Args:
  #   path: Path to the file with the transcript information with 1-based
  #         start end positions of genomic ranges.
  #
  # Returns:
  #   A data.table keyed by chrom, chromStart, and chromEnd

  d <- read.table(path)
  colnames(d) <- c("chrom", "chromStart", "chromEnd", "strand", "name")
  bed1 <- data.table(d)
  setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

ReadBedTranscriptRanges <- function(path) {
  # Read transcript ranges and strands from a bed format file.
  # Mostly for testing.
  #
  # Args:
  #   path: Path to the file with the transcript information (in bed format)
  #
  # Returns:
  #   A data.table keyed by chrom, chromStart, and chromEnd

  names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
  bed <- read.table(path, col.names = names, as.is = TRUE)

  # Delete duplicate entries in the BED file
  bed <- dplyr::distinct(bed, chrom, chromStart, chromEnd, strand, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  bed$chromStart <- bed$chromStart + 1L

  bed1 <- data.table(bed)
  setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

.trans.ranges <<-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh37.tsv")

.old.trans.ranges <<-
  ReadBedTranscriptRanges("data-raw/CCDS-intranscript-GRCh37p13_CURRENT_20131024.bed")
