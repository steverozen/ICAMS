# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/create_abundance_files.R")
source("data-raw/create_catalogs.R")
source("data-raw/create_order_for_DBS136_plotting.R")
source("data-raw/create_ranges.R")

usethis::use_data(catalog.row.order,
                  trans.ranges.GRCh37,
                  trans.ranges.GRCh38,
                  overwrite = TRUE)

usethis::use_data(to.reorder.SBS.192.for.plotting,
                  to.reorder.DBS.144.for.plotting,
                  order.for.DBS.136.plotting,
                  empty.cats,
                  exome.ranges.GRCh37,
                  exome.ranges.GRCh38,
                  exome.ranges.stranded.GRCh37,
                  exome.ranges.stranded.GRCh38,
                  catalog.row.headers.SBS.96,
                  catalog.row.headers.SBS.192,
                  catalog.row.headers.SBS.1536,
                  catalog.row.headers.DBS.78,
                  catalog.row.headers.DBS.144,
                  catalog.row.headers.DBS.136,
                  catalog.row.headers.ID,
                  abundance.2bp.genome.GRCh37,
                  abundance.3bp.genome.GRCh37,
                  abundance.4bp.genome.GRCh37,
                  abundance.5bp.genome.GRCh37,
                  abundance.2bp.stranded.GRCh37,
                  abundance.3bp.stranded.GRCh37,
                  abundance.4bp.stranded.GRCh37,
                  abundance.5bp.stranded.GRCh37,
                  abundance.2bp.exome.GRCh37,
                  abundance.3bp.exome.GRCh37,
                  abundance.4bp.exome.GRCh37,
                  abundance.5bp.exome.GRCh37,
                  abundance.2bp.genome.GRCh38,
                  abundance.3bp.genome.GRCh38,
                  abundance.4bp.genome.GRCh38,
                  abundance.5bp.genome.GRCh38,
                  abundance.2bp.stranded.GRCh38,
                  abundance.3bp.stranded.GRCh38,
                  abundance.4bp.stranded.GRCh38,
                  abundance.5bp.stranded.GRCh38,
                  abundance.2bp.exome.GRCh38,
                  abundance.3bp.exome.GRCh38,
                  abundance.4bp.exome.GRCh38,
                  abundance.5bp.exome.GRCh38,
                  abundance.2bp.flat,
                  abundance.3bp.flat,
                  abundance.4bp.flat,
                  abundance.5bp.flat,
                  abundance.2bp.flat.stranded,
                  abundance.3bp.flat.stranded,
                  abundance.3bp.transcript.unstranded.GRCh37,
                  abundance.2bp.transcript.unstranded.GRCh37,
                  abundance.3bp.transcript.unstranded.GRCh38,
                  abundance.2bp.transcript.unstranded.GRCh38,
                  homopolymer.ms.regex.pattern,
                  internal = TRUE,
                  overwrite = TRUE)

