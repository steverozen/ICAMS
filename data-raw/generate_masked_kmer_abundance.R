# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

#########################################################################################
# Code for generating masked exome k-mer counts
tmp <- proc.time()
hs37d5_masked_exome_unstranded_2bp <-
  GetExomeKmerCounts(k = 2, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_unstranded_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_unstranded_2bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_2bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_unstranded_3bp <-
  GetExomeKmerCounts(k = 3, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_unstranded_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_unstranded_3bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_3bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_unstranded_4bp <-
  GetExomeKmerCounts(k = 4, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_unstranded_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_unstranded_4bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_4bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_unstranded_5bp <-
  GetExomeKmerCounts(k = 5, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_unstranded_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_unstranded_5bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_exome_unstranded_2bp <-
  GetExomeKmerCounts(k = 2, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_unstranded_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_unstranded_2bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_2bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_unstranded_3bp <-
  GetExomeKmerCounts(k = 3, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_unstranded_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_unstranded_3bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_3bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_unstranded_4bp <-
  GetExomeKmerCounts(k = 4, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_unstranded_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_unstranded_4bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_4bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_unstranded_5bp <-
  GetExomeKmerCounts(k = 5, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_unstranded_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_unstranded_5bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_5bp.csv")

########################################################################################

########################################################################################
# Code for generating masked stranded exome k-mer counts
tmp <- proc.time()
hs37d5_masked_exome_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg19",
                        stranded.ranges = exome.ranges.stranded.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_stranded_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_stranded_2bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_2bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg19",
                        stranded.ranges = exome.ranges.stranded.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_stranded_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_stranded_3bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_3bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg19",
                        stranded.ranges = exome.ranges.stranded.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_stranded_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_stranded_4bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_4bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg19",
                        stranded.ranges = exome.ranges.stranded.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_stranded_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_stranded_5bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_exome_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg38",
                        stranded.ranges = exome.ranges.stranded.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_stranded_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_stranded_2bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_2bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg38",
                        stranded.ranges = exome.ranges.stranded.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_stranded_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_stranded_3bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_3bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg38",
                        stranded.ranges = exome.ranges.stranded.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_stranded_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_stranded_4bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_4bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg38",
                        stranded.ranges = exome.ranges.stranded.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_stranded_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_stranded_5bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_5bp.csv")

#########################################################################################
# Code for generating masked genome stranded k-mer counts
tmp <- proc.time()
hs37d5_masked_genome_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg19",
                        stranded.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_stranded_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_stranded_2bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_2bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg19",
                        stranded.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_stranded_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_stranded_3bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_3bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg19",
                        stranded.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_stranded_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_stranded_4bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_4bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg19",
                        stranded.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_stranded_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_stranded_5bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_genome_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg38",
                        stranded.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_stranded_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_stranded_2bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_2bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg38",
                        stranded.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_stranded_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_stranded_3bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_3bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg38",
                        stranded.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_stranded_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_stranded_4bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_4bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg38",
                        stranded.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_stranded_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_stranded_5bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_5bp.csv")

#########################################################################################
# Code for generating masked genome unstranded k-mer counts
tmp <- proc.time()
hs37d5_masked_genome_unstranded_2bp <-
  GetGenomeKmerCounts(k = 2, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_unstranded_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_unstranded_2bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_2bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_unstranded_3bp <-
  GetGenomeKmerCounts(k = 3, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_unstranded_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_unstranded_3bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_3bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_unstranded_4bp <-
  GetGenomeKmerCounts(k = 4, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_unstranded_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_unstranded_4bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_4bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_unstranded_5bp <-
  GetGenomeKmerCounts(k = 5, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_unstranded_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_unstranded_5bp,
          "data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_genome_unstranded_2bp <-
  GetGenomeKmerCounts(k = 2, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_unstranded_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_unstranded_2bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_2bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_unstranded_3bp <-
  GetGenomeKmerCounts(k = 3, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_unstranded_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_unstranded_3bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_3bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_unstranded_4bp <-
  GetGenomeKmerCounts(k = 4, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_unstranded_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_unstranded_4bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_4bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_unstranded_5bp <-
  GetGenomeKmerCounts(k = 5, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_unstranded_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_unstranded_5bp,
          "data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_5bp.csv")
