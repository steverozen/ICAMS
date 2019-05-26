# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh37.GENECODEv30.csv")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh38.GENECODEv30.csv")

exome.ranges.GRCh37 <-
  ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh37.SureSelectv6.csv")

exome.ranges.stranded.GRCh37 <-
  CreateExomeStrandedRanges("data-raw/ranges/ExomeRanges.GRCh37.SureSelectv6.csv")

exome.ranges.GRCh38 <-
  ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh38.SureSelectv6.csv")

exome.ranges.stranded.GRCh38 <-
  CreateExomeStrandedRanges("data-raw/ranges/ExomeRanges.GRCh38.SureSelectv6.csv")



