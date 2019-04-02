# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh38.csv")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCh37.csv")

trans.ranges.GRCm38 <-
  ReadTranscriptRanges("data-raw/TranscriptRanges.GRCm38.csv")

