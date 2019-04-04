# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh38.csv")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh37.csv")

exome.ranges.GRCh38 <-
  ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh38.csv")

exome.ranges.GRCh37 <-
  ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh37.csv")

