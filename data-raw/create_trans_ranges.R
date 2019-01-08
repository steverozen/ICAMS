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

.trans.ranges <<-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh37.tsv")

