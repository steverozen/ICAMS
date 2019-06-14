
# Do you have read.and.prep.192.duke.nus.catalog?
test.cat <- 
  read.and.prep.192.duke.nus.catalog(
    file.path(devtools::package_file("data-raw"),
              "possible-future-code/CrVI_MEFs.txt"))

c96 <- as.catalog(test.cat[[1]]) 
# Works, ref.genome = NULL, region = "unknown", abundance = NULL,
# catalog.type = "counts" (the default for as.catalog now)

# c192 <- as.catalog(test.cat[[2]]) Does not work because the 192 catalog
# read from Duke-NUS format does not have the right rownames

PlotCatalogToPdf(c96, 
                 file = file.path(devtools::package_file("data-raw"), 
                                  "c96.counts.pdf"))

# Does not work but could be made to work; I am guessing not urgent
# c96.sig <- TransformCatalog(c96.sig, target.ref.genome = NULL,
#                            target.catalog.type = "counts.signature")

# Nanhai is working on creating catalogs from VCFs when here is no
# transcrpt range information

## New
## 
## Steve note -- there is no need for a reference genome in
## TransformCatalog unless you are changing genomes, and in fact
## What you reall need are abundances. The ref.genome argument is
## just a convenience that lets you look up abundances.

BiocManager::install("BSgenome.Mmusculus.UCSC.mm8")
library(BSgenome.Mmusculus.UCSC.mm8)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

strelka.files <- 
  dir(devtools::package_file("data-raw/possible-future-code"),
      pattern = ".vcf", full.names = TRUE)


StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
  files = strelka.files,
  ref.genome = BSgenome.Mmusculus.UCSC.mm10,
  trans.ranges = NULL,
  region = "genome",
  output.file = "test.pdf")


