# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

######################################################################################

abundance.2bp.genome.GRCh37 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_genome_2bp.csv")
abundance.3bp.genome.GRCh37 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_genome_3bp.csv")
abundance.4bp.genome.GRCh37 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/hs37d5_masked_genome_4bp.csv")
abundance.5bp.genome.GRCh37 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/hs37d5_masked_genome_5bp.csv")

abundance.2bp.stranded.GRCh37 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_stranded_2bp.csv")
abundance.3bp.stranded.GRCh37 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_stranded_3bp.csv")
abundance.4bp.stranded.GRCh37 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/hs37d5_masked_stranded_4bp.csv")
abundance.5bp.stranded.GRCh37 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/hs37d5_masked_stranded_5bp.csv")

abundance.2bp.exome.GRCh37 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_exome_2bp.csv")
abundance.3bp.exome.GRCh37 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/hs37d5_masked_exome_3bp.csv")
abundance.4bp.exome.GRCh37 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/hs37d5_masked_exome_4bp.csv")
abundance.5bp.exome.GRCh37 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/hs37d5_masked_exome_5bp.csv")

abundance.3bp.transcript.unstranded.GRCh37 <-
  Collapse192AbundanceTo96(abundance.3bp.stranded.GRCh37)
abundance.2bp.transcript.unstranded.GRCh37 <-
  Collapse144AbundanceTo78(abundance.2bp.stranded.GRCh37)
#######################################################################################
abundance.2bp.genome.GRCh38 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_genome_2bp.csv")
abundance.3bp.genome.GRCh38 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_genome_3bp.csv")
abundance.4bp.genome.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38_masked_genome_4bp.csv")
abundance.5bp.genome.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38_masked_genome_5bp.csv")

abundance.2bp.stranded.GRCh38 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_stranded_2bp.csv")
abundance.3bp.stranded.GRCh38 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_stranded_3bp.csv")
abundance.4bp.stranded.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38_masked_stranded_4bp.csv")
abundance.5bp.stranded.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38_masked_stranded_5bp.csv")

abundance.2bp.exome.GRCh38 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_exome_2bp.csv")
abundance.3bp.exome.GRCh38 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh38_masked_exome_3bp.csv")
abundance.4bp.exome.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38_masked_exome_4bp.csv")
abundance.5bp.exome.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38_masked_exome_5bp.csv")

abundance.3bp.transcript.unstranded.GRCh38 <-
  Collapse192AbundanceTo96(abundance.3bp.stranded.GRCh38)
abundance.2bp.transcript.unstranded.GRCh38 <-
  Collapse144AbundanceTo78(abundance.2bp.stranded.GRCh38)
#########################################################################################

abundance.2bp.flat <- rep(1, length(abundance.2bp.genome.GRCh37))
names(abundance.2bp.flat) <- names(abundance.2bp.genome.GRCh37)
abundance.2bp.flat.stranded <- rep(1, length(abundance.2bp.stranded.GRCh37))
names(abundance.2bp.flat.stranded) <- names(abundance.2bp.stranded.GRCh37)

abundance.3bp.flat <- rep(1, length(abundance.3bp.genome.GRCh37))
names(abundance.3bp.flat) <- names(abundance.3bp.genome.GRCh37)
abundance.3bp.flat.stranded <- rep(1, length(abundance.3bp.stranded.GRCh37))
names(abundance.3bp.flat.stranded) <- names(abundance.3bp.stranded.GRCh37)

abundance.4bp.flat <- rep(1, length(abundance.4bp.genome.GRCh37))
names(abundance.4bp.flat) <- names(abundance.4bp.genome.GRCh37)

abundance.5bp.flat <- rep(1, length(abundance.5bp.genome.GRCh37))
names(abundance.5bp.flat) <- names(abundance.5bp.genome.GRCh37)
