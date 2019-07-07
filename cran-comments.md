## Resubmission
This is a resubmission. In this version I have:
* Ensured that every  
  opar <- par(...)
  is followed IMMEDIATEDLY by
  on.exit(par(opar))
  (not farther down in the code)

## Test environments
* Local Windows 10 install: R 3.6.0
* Travis-CI Ubuntu 14.04.5. LTS: R 3.6.0 and R devel
* local OS X Mojave install: R 3.6.0

## R CMD check results
There were no ERRORs, WARGNINGs, or NOTEs

## Downstream dependencies
None
