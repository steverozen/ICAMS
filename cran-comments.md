## Resubmission
This is a resubmission. In this version I have:
* Replaced all instances of cat() and print() from the code
  with message(), warning(), or stop() as appropriate.
  These are all either part of error handling or part of 
  tracing capability that is turned off by default.
* Changed 2 print() calls and 1 fwrite() call that were part of error handling to
  write.csv() calls writing to a file in tempdir() (grep tempfile */*.R);
  one of these was in CheckSeqContextInVCF as noted in the manual
  CRAN check.
* Replaced T with TRUE and F with FALSE everywhere.

## Test environments
* Local Windows 10 install: R 3.6.0
* Travis-CI Ubuntu 14.04.5. LTS: R 3.6.0 and R devel
* local OS X Mojave install: R 3.6.0

## R CMD check results
There were no ERRORs, WARGNINGs, or NOTEs

## Downstream dependencies
None
