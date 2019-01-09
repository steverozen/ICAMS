ReadAbundance3Bp <- function(path) {
  # Read data from a nucleotide abundance file with 3 base pairs.
  #
  # Args:
  #   path: Path to the file with the nucleotide abundance information
  #         with 3 base pairs.
  #
  # Returns:
  #   A matrix whose row names indicate 32 different types of 3 base pairs
  #   combinations while its column contains the occurrences of each type.

  dt <- fread(path)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

ReadAbundance4Bp <- function(path) {
  # Read data from a nucleotide abundance file with 4 base pairs.
  #
  # Args:
  #   path: Path to the file with the nucleotide abundance information
  #         with 4 base pairs.
  #
  # Returns:
  #   A matrix whose row names indicate 10 different types of 2 base pairs
  #   combinations while its column contains the occurrences of each type.

  dt <- fread(path)
  colnames(dt) <- c("4bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 3) %in% canonical.ref,
           substr(dt[[1]], 2, 3),
           revc(substr(dt[[1]], 2, 3)))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

ReadAbundance5Bp <- function(path) {
  # Read data from a nucleotide abundance file with 5 base pairs.
  #
  # Args:
  #   path: Path to the file with the nucleotide abundance information
  #         with 5 base pairs.
  #
  # Returns:
  #   A matrix whose row names indicate 512 different types of 5 base pairs
  #   combinations while its column contains the occurrences of each type.

  dt <- fread(path)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

.abundance.2bp <<- ReadAbundance4Bp("data-raw/hs37d5_masked_4bp.txt")
.abundance.3bp <<- ReadAbundance3Bp("data-raw/hs37d5_masked_3bp.txt")
.abundance.5bp <<- ReadAbundance5Bp("data-raw/hs37d5_masked_5bp.txt")
