# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

######################################################################################

abundance.2bp.genome.unstranded.GRCh37 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_2bp.csv")
abundance.3bp.genome.unstranded.GRCh37 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_3bp.csv")
abundance.4bp.genome.unstranded.GRCh37 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_4bp.csv")
abundance.5bp.genome.unstranded.GRCh37 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_unstranded_5bp.csv")

abundance.2bp.genome.stranded.GRCh37 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_2bp.csv")
abundance.3bp.genome.stranded.GRCh37 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_genome_stranded_3bp.csv")

abundance.2bp.exome.unstranded.GRCh37 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_2bp.csv")
abundance.3bp.exome.unstranded.GRCh37 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_3bp.csv")
abundance.4bp.exome.unstranded.GRCh37 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_4bp.csv")
abundance.5bp.exome.unstranded.GRCh37 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_unstranded_5bp.csv")

abundance.2bp.exome.stranded.GRCh37 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_2bp.csv")
abundance.3bp.exome.stranded.GRCh37 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/GRCh37/hs37d5_masked_exome_stranded_3bp.csv")

abundance.3bp.transcript.unstranded.GRCh37 <-
  Collapse192AbundanceTo96(abundance.3bp.genome.stranded.GRCh37)
abundance.2bp.transcript.unstranded.GRCh37 <-
  Collapse144AbundanceTo78(abundance.2bp.genome.stranded.GRCh37)
#######################################################################################
abundance.2bp.genome.unstranded.GRCh38 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_2bp.csv")
abundance.3bp.genome.unstranded.GRCh38 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_3bp.csv")
abundance.4bp.genome.unstranded.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_4bp.csv")
abundance.5bp.genome.unstranded.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_unstranded_5bp.csv")

abundance.2bp.genome.stranded.GRCh38 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_2bp.csv")
abundance.3bp.genome.stranded.GRCh38 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_3bp.csv")
abundance.4bp.genome.stranded.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_4bp.csv")
abundance.5bp.genome.stranded.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_genome_stranded_5bp.csv")

abundance.2bp.exome.unstranded.GRCh38 <-
  CreateDinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_2bp.csv")
abundance.3bp.exome.unstranded.GRCh38 <-
  CreateTrinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_3bp.csv")
abundance.4bp.exome.unstranded.GRCh38 <-
  CreateTetranucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_4bp.csv")
abundance.5bp.exome.unstranded.GRCh38 <-
  CreatePentanucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_unstranded_5bp.csv")

abundance.2bp.exome.stranded.GRCh38 <-
  CreateStrandedDinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_2bp.csv")
abundance.3bp.exome.stranded.GRCh38 <-
  CreateStrandedTrinucAbundance("data-raw/new_masked_abundance/GRCh38/GRCh38_masked_exome_stranded_3bp.csv")

abundance.3bp.transcript.unstranded.GRCh38 <-
  Collapse192AbundanceTo96(abundance.3bp.genome.stranded.GRCh38)
abundance.2bp.transcript.unstranded.GRCh38 <-
  Collapse144AbundanceTo78(abundance.2bp.genome.stranded.GRCh38)
#########################################################################################

abundance.2bp.flat.unstranded <- rep(1, length(abundance.2bp.genome.unstranded.GRCh37))
names(abundance.2bp.flat.unstranded) <- names(abundance.2bp.genome.unstranded.GRCh37)
abundance.2bp.flat.stranded <- rep(1, length(abundance.2bp.genome.stranded.GRCh37))
names(abundance.2bp.flat.stranded) <- names(abundance.2bp.genome.stranded.GRCh37)

abundance.3bp.flat.unstranded <- rep(1, length(abundance.3bp.genome.unstranded.GRCh37))
names(abundance.3bp.flat.unstranded) <- names(abundance.3bp.genome.unstranded.GRCh37)
abundance.3bp.flat.stranded <- rep(1, length(abundance.3bp.genome.stranded.GRCh37))
names(abundance.3bp.flat.stranded) <- names(abundance.3bp.genome.stranded.GRCh37)

abundance.4bp.flat.unstranded <- rep(1, length(abundance.4bp.genome.unstranded.GRCh37))
names(abundance.4bp.flat.unstranded) <- names(abundance.4bp.genome.unstranded.GRCh37)

abundance.5bp.flat.unstranded <- rep(1, length(abundance.5bp.genome.unstranded.GRCh37))
names(abundance.5bp.flat.unstranded) <- names(abundance.5bp.genome.unstranded.GRCh37)


tmp.hg19 <- list(
  exome = list(
    "78"   = ICAMS:::abundance.2bp.exome.unstranded.GRCh37,
    "144"  = ICAMS:::abundance.2bp.exome.stranded.GRCh37,
    "96"   = ICAMS:::abundance.3bp.exome.unstranded.GRCh37,
    "192"  = ICAMS:::abundance.3bp.exome.stranded.GRCh37,
    "136"  = ICAMS:::abundance.4bp.exome.unstranded.GRCh37,
    "1536" = ICAMS:::abundance.5bp.exome.unstranded.GRCh37
  ),
  transcript = list(
    "78"   = ICAMS:::abundance.2bp.transcript.unstranded.GRCh37,
    "144"  = ICAMS:::abundance.2bp.genome.stranded.GRCh37,
    "96"   = ICAMS:::abundance.3bp.transcript.unstranded.GRCh37,
    "192"  = ICAMS:::abundance.3bp.genome.stranded.GRCh37
    # Nothing here for 136 or 1536
  ),
  genome = list(
    "96"   = ICAMS:::abundance.3bp.genome.unstranded.GRCh37,
    "78"   = ICAMS:::abundance.2bp.genome.unstranded.GRCh37,
    "136"  = ICAMS:::abundance.4bp.genome.unstranded.GRCh37,
    "1536" = ICAMS:::abundance.5bp.genome.unstranded.GRCh37
  )
)

tmp.GRCh38 <- list(
  exome = list(
    "78"   = ICAMS:::abundance.2bp.exome.unstranded.GRCh38,
    "144"  = ICAMS:::abundance.2bp.exome.stranded.GRCh38,
    "96"   = ICAMS:::abundance.3bp.exome.unstranded.GRCh38,
    "192"  = ICAMS:::abundance.3bp.exome.stranded.GRCh38,
    "136"  = ICAMS:::abundance.4bp.exome.unstranded.GRCh38,
    "1536" = ICAMS:::abundance.5bp.exome.unstranded.GRCh38
    ),
  transcript = list(
    "78"   = ICAMS:::abundance.2bp.transcript.unstranded.GRCh38,
    "144"  = ICAMS:::abundance.2bp.genome.stranded.GRCh38,
    "96"   = ICAMS:::abundance.3bp.transcript.unstranded.GRCh38,
    "192"  = ICAMS:::abundance.3bp.genome.stranded.GRCh38
    # Nothing here for 136 or 1536
  ),
  genome = list(
    "96"   = ICAMS:::abundance.3bp.genome.unstranded.GRCh38,
    "78"   = ICAMS:::abundance.2bp.genome.unstranded.GRCh38,
    "136"  = ICAMS:::abundance.4bp.genome.unstranded.GRCh38,
    "1536" = ICAMS:::abundance.5bp.genome.unstranded.GRCh38
    )
  )

all.abundance <- list(
  BSgenome.Hsapiens.1000genomes.hs37d5 = tmp.hg19,
  BSgenome.Hsapiens.UCSC.hg38 =   tmp.GRCh38
)
