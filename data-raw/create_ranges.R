# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh37 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh37.GENECODEv19.csv")

trans.ranges.GRCh38 <-
  ReadTranscriptRanges("data-raw/ranges/TranscriptRanges.GRCh38.GENECODEv28.csv")

#exome.ranges.GRCh37 <-
  #ReadBedRanges("data-raw/ranges/SureSelect_Human_All_Exon_V6.hg19.bed")

#exome.ranges.GRCh38 <-
  #ReadBedRanges("data-raw/ranges/SureSelect_Human_All_Exon_V6.GRCh38.bed")



