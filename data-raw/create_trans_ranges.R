# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh38.tsv")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh37.tsv")

