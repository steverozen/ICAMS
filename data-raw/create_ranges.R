# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

trans.ranges.GRCh37 <-
  CreateTransRanges("data-raw/ranges/gencode.GRCh37.annotation.gtf.gz")

trans.ranges.GRCh38 <-
  CreateTransRanges("data-raw/ranges/gencode.GRCh38.p12.annotation.gff3.gz")

exome.ranges.GRCh37 <-
  ReadBedRanges("data-raw/ranges/SureSelect_Human_All_Exon_V6.hg19.bed")

exome.ranges.GRCh38 <-
  ReadBedRanges("data-raw/ranges/SureSelect_Human_All_Exon_V6.GRCh38.bed")



