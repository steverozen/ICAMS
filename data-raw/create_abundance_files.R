abundance.2bp.genome.GRCh37 <- CreateDinucAbundance("data-raw/hs37d5_masked_4bp.txt")
abundance.3bp.genome.GRCh37 <- CreateTrinucAbundance("data-raw/hs37d5_masked_3bp.txt")
abundance.4bp.genome.GRCh37 <- CreateTetranucAbundance("data-raw/hs37d5_masked_4bp.txt")
abundance.5bp.genome.GRCh37 <- CreatePentanucAbundance("data-raw/hs37d5_masked_5bp.txt")

abundance.2bp.exome.GRCh37 <- CreateDinucAbundance("data-raw/hs37d5_masked_AgilentV6_4bp.txt")
abundance.3bp.exome.GRCh37 <- CreateTrinucAbundance("data-raw/hs37d5_masked_AgilentV6_3bp.txt")
abundance.4bp.exome.GRCh37 <- CreateTetranucAbundance("data-raw/hs37d5_masked_AgilentV6_4bp.txt")
abundance.5bp.exome.GRCh37 <- CreatePentanucAbundance("data-raw/hs37d5_masked_AgilentV6_5bp.txt")
