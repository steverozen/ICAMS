# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/code/load_abundance_from_files.R")
source("data-raw/code/create_catalogs.R")
source("data-raw/code/create_order_for_DBS136_plotting.R")
source("data-raw/code/create_ranges.R")
source("data-raw/code/create_gene_expression_data.R")
source("data-raw/code/create_ICAMS_SigPro_ID.R")
source("data-raw/code/create_catalogs_COSMIC.R")

usethis::use_data(catalog.row.order,
                  trans.ranges.GRCh37,
                  trans.ranges.GRCh38,
                  trans.ranges.GRCm38,
                  all.abundance,
                  gene.expression.data.HepG2,
                  gene.expression.data.MCF10A,
                  overwrite = TRUE)

usethis::use_data(to.reorder.SBS.192.for.plotting,
                  to.reorder.DBS.144.for.plotting,
                  order.for.DBS.136.plotting,
                  empty.cats,
                  catalog.row.headers,
                  catalog.row.headers.sp,
                  catalog.row.headers.cosmic,
                  catalog.row.headers.SBS.96.v1,
                  abundance.2bp.genome.unstranded.GRCh37,
                  abundance.3bp.genome.unstranded.GRCh37,
                  abundance.4bp.genome.unstranded.GRCh37,
                  abundance.5bp.genome.unstranded.GRCh37,
                  abundance.2bp.transcript.stranded.GRCh37,
                  abundance.3bp.transcript.stranded.GRCh37,
                  abundance.2bp.exome.unstranded.GRCh37,
                  abundance.3bp.exome.unstranded.GRCh37,
                  abundance.4bp.exome.unstranded.GRCh37,
                  abundance.5bp.exome.unstranded.GRCh37,
                  abundance.2bp.exome.stranded.GRCh37,
                  abundance.3bp.exome.stranded.GRCh37,
                  abundance.2bp.genome.unstranded.GRCh38,
                  abundance.3bp.genome.unstranded.GRCh38,
                  abundance.4bp.genome.unstranded.GRCh38,
                  abundance.5bp.genome.unstranded.GRCh38,
                  abundance.2bp.transcript.stranded.GRCh38,
                  abundance.3bp.transcript.stranded.GRCh38,
                  abundance.2bp.exome.unstranded.GRCh38,
                  abundance.3bp.exome.unstranded.GRCh38,
                  abundance.4bp.exome.unstranded.GRCh38,
                  abundance.5bp.exome.unstranded.GRCh38,
                  abundance.2bp.exome.stranded.GRCh38,
                  abundance.3bp.exome.stranded.GRCh38,
                  abundance.2bp.genome.unstranded.GRCm38,
                  abundance.3bp.genome.unstranded.GRCm38,
                  abundance.4bp.genome.unstranded.GRCm38,
                  abundance.5bp.genome.unstranded.GRCm38,
                  abundance.2bp.transcript.stranded.GRCm38,
                  abundance.3bp.transcript.stranded.GRCm38,
                  abundance.2bp.exome.unstranded.GRCm38,
                  abundance.3bp.exome.unstranded.GRCm38,
                  abundance.4bp.exome.unstranded.GRCm38,
                  abundance.5bp.exome.unstranded.GRCm38,
                  abundance.2bp.exome.stranded.GRCm38,
                  abundance.3bp.exome.stranded.GRCm38,
                  abundance.2bp.flat.unstranded,
                  abundance.3bp.flat.unstranded,
                  abundance.4bp.flat.unstranded,
                  abundance.5bp.flat.unstranded,
                  abundance.2bp.flat.stranded,
                  abundance.3bp.transcript.unstranded.GRCh37,
                  abundance.2bp.transcript.unstranded.GRCh37,
                  abundance.3bp.transcript.unstranded.GRCh38,
                  abundance.2bp.transcript.unstranded.GRCh38,
                  abundance.3bp.transcript.unstranded.GRCm38,
                  abundance.2bp.transcript.unstranded.GRCm38,
                  homopolymer.ms.regex.pattern,
                  SigPro.to.ICAMS.ID,
                  ICAMS.to.SigPro.ID,
                  flat.abundance,
                  internal = TRUE,
                  overwrite = TRUE)

