# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

#########################################################################################
# Code for generating exome k-mer counts
tmp <- proc.time()
hs37d5_masked_exome_2bp <-
  GetExomeKmerCounts(k = 2, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_2bp,
          "data-raw/new_masked_abundance/hs37d5_masked_exome_2bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_3bp <-
  GetExomeKmerCounts(k = 3, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_3bp,
          "data-raw/new_masked_abundance/hs37d5_masked_exome_3bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_4bp <-
  GetExomeKmerCounts(k = 4, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_4bp,
          "data-raw/new_masked_abundance/hs37d5_masked_exome_4bp.csv")

tmp <- proc.time()
hs37d5_masked_exome_5bp <-
  GetExomeKmerCounts(k = 5, ref.genome = "hg19",
                     exome.ranges = exome.ranges.GRCh37,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_exome_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_exome_5bp,
          "data-raw/new_masked_abundance/hs37d5_masked_exome_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_exome_2bp <-
  GetExomeKmerCounts(k = 2, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_2bp,
          "data-raw/new_masked_abundance/GRCh38_masked_exome_2bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_3bp <-
  GetExomeKmerCounts(k = 3, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_3bp,
          "data-raw/new_masked_abundance/GRCh38_masked_exome_3bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_4bp <-
  GetExomeKmerCounts(k = 4, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_4bp,
          "data-raw/new_masked_abundance/GRCh38_masked_exome_4bp.csv")

tmp <- proc.time()
GRCh38_masked_exome_5bp <-
  GetExomeKmerCounts(k = 5, ref.genome = "hg38",
                     exome.ranges = exome.ranges.GRCh38,
                     filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_exome_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_exome_5bp,
          "data-raw/new_masked_abundance/GRCh38_masked_exome_5bp.csv")

#########################################################################################
# Code for generating stranded k-mer counts
tmp <- proc.time()
hs37d5_masked_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg19",
                        trans.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_stranded_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_stranded_2bp,
          "data-raw/new_masked_abundance/hs37d5_masked_stranded_2bp.csv")

tmp <- proc.time()
hs37d5_masked_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg19",
                        trans.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_stranded_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_stranded_3bp,
          "data-raw/new_masked_abundance/hs37d5_masked_stranded_3bp.csv")

tmp <- proc.time()
hs37d5_masked_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg19",
                        trans.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_stranded_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_stranded_4bp,
          "data-raw/new_masked_abundance/hs37d5_masked_stranded_4bp.csv")

tmp <- proc.time()
hs37d5_masked_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg19",
                        trans.ranges = trans.ranges.GRCh37,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_stranded_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_stranded_5bp,
          "data-raw/new_masked_abundance/hs37d5_masked_stranded_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_stranded_2bp <-
  GetStrandedKmerCounts(k = 2, ref.genome = "hg38",
                        trans.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_stranded_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_stranded_2bp,
          "data-raw/new_masked_abundance/GRCh38_masked_stranded_2bp.csv")

tmp <- proc.time()
GRCh38_masked_stranded_3bp <-
  GetStrandedKmerCounts(k = 3, ref.genome = "hg38",
                        trans.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_stranded_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_stranded_3bp,
          "data-raw/new_masked_abundance/GRCh38_masked_stranded_3bp.csv")

tmp <- proc.time()
GRCh38_masked_stranded_4bp <-
  GetStrandedKmerCounts(k = 4, ref.genome = "hg38",
                        trans.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_stranded_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_stranded_4bp,
          "data-raw/new_masked_abundance/GRCh38_masked_stranded_4bp.csv")

tmp <- proc.time()
GRCh38_masked_stranded_5bp <-
  GetStrandedKmerCounts(k = 5, ref.genome = "hg38",
                        trans.ranges = trans.ranges.GRCh38,
                        filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_stranded_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_stranded_5bp,
          "data-raw/new_masked_abundance/GRCh38_masked_stranded_5bp.csv")

#########################################################################################
# Code for generating genome k-mer counts
tmp <- proc.time()
hs37d5_masked_genome_2bp <-
  GetGenomeKmerCounts(k = 2, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_2bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_2bp,
          "data-raw/new_masked_abundance/hs37d5_masked_genome_2bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_3bp <-
  GetGenomeKmerCounts(k = 3, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_3bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_3bp,
          "data-raw/new_masked_abundance/hs37d5_masked_genome_3bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_4bp <-
  GetGenomeKmerCounts(k = 4, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_4bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_4bp,
          "data-raw/new_masked_abundance/hs37d5_masked_genome_4bp.csv")

tmp <- proc.time()
hs37d5_masked_genome_5bp <-
  GetGenomeKmerCounts(k = 5, ref.genome = "hg19",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg19.txt")
t.hs37d5_masked_genome_5bp <- proc.time() - tmp
write.csv(hs37d5_masked_genome_5bp,
          "data-raw/new_masked_abundance/hs37d5_masked_genome_5bp.csv")

#########################################################################################
tmp <- proc.time()
GRCh38_masked_genome_2bp <-
  GetGenomeKmerCounts(k = 2, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_2bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_2bp,
          "data-raw/new_masked_abundance/GRCh38_masked_genome_2bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_3bp <-
  GetGenomeKmerCounts(k = 3, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_3bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_3bp,
          "data-raw/new_masked_abundance/GRCh38_masked_genome_3bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_4bp <-
  GetGenomeKmerCounts(k = 4, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_4bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_4bp,
          "data-raw/new_masked_abundance/GRCh38_masked_genome_4bp.csv")

tmp <- proc.time()
GRCh38_masked_genome_5bp <-
  GetGenomeKmerCounts(k = 5, ref.genome = "hg38",
                      filter.path = "data-raw/simple_repeat_files/simpleRepeat.hg38.txt")
t.GRCh38_masked_genome_5bp <- proc.time() - tmp
write.csv(GRCh38_masked_genome_5bp,
          "data-raw/new_masked_abundance/GRCh38_masked_genome_5bp.csv")
