## Resubmission
This is a resubmission. In this version I have:
* In DESCRIPTION, replaced https://doi.org/10.1101/gr.230219.117
  with <doi:10.1101/gr.230219.117>.
* In DESCRIPTION (and README) explained the acronym ICAMS 
  (In-depth Characterization and Analysis of Mutational Signatures).
* Added small executable examples too all Rd-files.
* Ensured that all examples/vignettes/tests that write files write only
  to tempdir().
* Removed all instances of <<- from the code.

## Test environments
* Local Windows 10 install: R 3.6.0
* Travis-CI Ubuntu 14.04.5. LTS: R 3.6.0 and R devel
* local OS X Mojave install: R 3.6.0

## R CMD check results
There were no ERRORs, WARGNINGs, or NOTEs

## Downstream dependencies
None
