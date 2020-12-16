# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh37.GENECODEv30.csv")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh38.GENECODEv30.csv")

trans.ranges.GRCm38 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCm38.GENECODEvM21.csv")

if (FALSE) {
  exome.ranges.GRCh37 <-
    ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh37.SureSelectv6.csv")

  exome.ranges.stranded.GRCh37 <-
    CreateExomeStrandedRanges(file = "data-raw/ranges/ExomeRanges.GRCh37.SureSelectv6.csv",
                              trans.ranges = trans.ranges.GRCh37)

  exome.ranges.GRCh38 <-
    ReadBedRanges("data-raw/ranges/ExomeRanges.GRCh38.SureSelectv6.csv")

  exome.ranges.stranded.GRCh38 <-
    CreateExomeStrandedRanges(file = "data-raw/ranges/ExomeRanges.GRCh38.SureSelectv6.csv",
                              trans.ranges = trans.ranges.GRCh38)
  
  exome.ranges.GRCm38 <-
    ReadBedRanges("data-raw/ranges/ExomeRanges.GRCm38.SureSelectv1.csv")
  
  exome.ranges.stranded.GRCm38 <-
    CreateExomeStrandedRanges(file = "data-raw/ranges/ExomeRanges.GRCm38.SureSelectv1.csv",
                              trans.ranges = trans.ranges.GRCm38)
  
}
