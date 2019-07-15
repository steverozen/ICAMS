## Resubmission
This is a resubmission. In this version I have:
* Moved two Bioconductor packages(BSgenome.Hsapiens.1000genomes.hs37d5, 
  BSgenome.Hsapiens.UCSC.hg38) from Imports to Suggests so that they can
  be used conditionally in "ICAMS" as requested by CRAN.

## Test environments
* Local Windows 10 install: R 3.6.0
* Travis-CI Ubuntu 14.04.5. LTS: R 3.6.0 and R devel
* local OS X Mojave install: R 3.6.0

## R CMD check results
There were no ERRORs, WARGNINGs, or NOTEs

## Downstream dependencies
None
