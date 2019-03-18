# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

abundance.2bp.genome.GRCh37 <-
  CreateDinucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_4bp.txt")
abundance.3bp.genome.GRCh37 <-
  CreateTrinucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_3bp.txt")
abundance.4bp.genome.GRCh37 <-
  CreateTetranucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_4bp.txt")
abundance.5bp.genome.GRCh37 <-
  CreatePentanucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_5bp.txt")

abundance.2bp.exome.GRCh37 <-
  CreateDinucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_AgilentV6_4bp.txt")
abundance.3bp.exome.GRCh37 <-
  CreateTrinucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_AgilentV6_3bp.txt")
abundance.4bp.exome.GRCh37 <-
  CreateTetranucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_AgilentV6_4bp.txt")
abundance.5bp.exome.GRCh37 <-
  CreatePentanucAbundance("data-raw/nucleotide_abundance/hs37d5_masked_AgilentV6_5bp.txt")

abundance.2bp.genome.GRCh38 <-
  CreateDinucAbundance("data-raw/nucleotide_abundance/hg38_masked_4bp.txt")
abundance.3bp.genome.GRCh38 <-
  CreateTrinucAbundance("data-raw/nucleotide_abundance/hg38_masked_3bp.txt")
abundance.4bp.genome.GRCh38 <-
  CreateTetranucAbundance("data-raw/nucleotide_abundance/hg38_masked_4bp.txt")
abundance.5bp.genome.GRCh38 <-
  CreatePentanucAbundance("data-raw/nucleotide_abundance/hg38_masked_5bp.txt")

abundance.2bp.exome.GRCh38 <-
  CreateDinucAbundance("data-raw/nucleotide_abundance/hg38_masked_AgilentV6_4bp.txt")
abundance.3bp.exome.GRCh38 <-
  CreateTrinucAbundance("data-raw/nucleotide_abundance/hg38_masked_AgilentV6_3bp.txt")
abundance.4bp.exome.GRCh38 <-
  CreateTetranucAbundance("data-raw/nucleotide_abundance/hg38_masked_AgilentV6_4bp.txt")
abundance.5bp.exome.GRCh38 <-
  CreatePentanucAbundance("data-raw/nucleotide_abundance/hg38_masked_AgilentV6_5bp.txt")
