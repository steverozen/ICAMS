## Resubmission
This is a resubmission. In this version I have:
* Written package names, software names and API names in 
  single quotes (e.g. 'ICAMS') in the DESCRIPTION.
* In functions that change the user's par() settings,
  there is now an immediate call of on.exit() to reset
  to the original settings.

## Test environments
* Local Windows 10 install: R 3.6.0
* Travis-CI Ubuntu 14.04.5. LTS: R 3.6.0 and R devel
* local OS X Mojave install: R 3.6.0

## R CMD check results
There were no ERRORs, WARGNINGs, or NOTEs

## Downstream dependencies
None
