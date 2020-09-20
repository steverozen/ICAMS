## Resubmission  
This is a resubmission. In this version I have:    

* Made changes to some test functions in testthat to use packages in Suggests  
  conditionally as requested by CRAN.

## Test environments
* Local Windows 10 install: R 4.0.2
* Local OS X Mojave install: R 3.6.2 and R 4.0.2
* Travis-CI Ubuntu 16.04.6 LTS: R oldrelease (3.6.3), release (4.0.0) and devel
* Winbuilder: R oldrelease (3.6.3), release (4.0.2) and devel

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTES: 

* checking for future file timestamps ... NOTE   
  unable to verify current time
  
  The website CRAN used to check the time is currently not available.
  
* checking CRAN incoming feasibility ... NOTE    
  Maintainer: 'Steve Rozen <steverozen@gmail.com>'     
  Days since last update: 2     
  
  I made changes to some test functions in testthat to use packages in Suggests conditionally as    requested by CRAN.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
