test.cat <- 
  read.and.prep.192.duke.nus.catalog(
    file.path(devtools::package_file("data-raw"),
              "possible-future-code/CrVI_MEFs.txt"))

c96 <- as.catalog(test.cat[[1]])
