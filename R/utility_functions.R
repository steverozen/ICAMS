#' Standardize the Chromosome name annotations for a data frame
#'
#' @param df A data frame whose first column contains the Chromosome name.
#' @keywords internal
#' @return A data frame whose Chromosome names are only in the form of 1:22,
#'   "X" and "Y".
#' @export
StandardChromName <- function(df) {
  # Is there any row in df whose Chromosome names start with "GL"?
  if (sum(grepl("^GL", df[[1]])) > 0) {
    df <- df[-grep("^GL", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names are "Hs37D5"?
  if (sum(grepl("^Hs", df[[1]])) > 0) {
    df <- df[-grep("^Hs", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names contain "M"?
  if (sum(grepl("M", df[[1]])) > 0) {
    df <- df[-grep("M", df[[1]]), ]
  }

  # Remove the "chr" character in the Chromosome's name
  df[, 1] <- sub(pattern = "chr", replacement = "", df[[1]])

  return(df)
}

#' Create a Transcript Range file from the raw GFF3 File
#'
#' @param path The name/path of the raw GFF3 File, or a complete URL.
#' @keywords internal
#' @return A data frame which contains chromosome name, start, end position,
#'   strand information and gene name. Only the following four gene types are
#'   kept to facilitate transcriptional strand bias analysis: protein_coding,
#'   retained_intron, processed_transcript and nonsense_mediated_decay.
#' @export
CreateTransRange <- function(path) {
  df <- read.csv(path, header = FALSE, fill = TRUE, nrows = 20)
  # Count the number of comment lines
  n <- sum(grepl("#", df[, 1]))

  # Read in the raw GFF3 File while skipping the comment lines
  dt <- fread(path, header = FALSE, sep = "\t", fill = TRUE, skip = n)

  dt1 <- dt[dt$V3 == "gene", ]

  # Select out the four gene types for transcriptional strand bias analysis
  idx <-
    grepl("protein_coding", dt1$V9) |
    grepl("retained_intron", dt1$V9) |
    grepl("processed_transcript", dt1$V9) |
    grepl("nonsense_mediated_decay", dt1$V9)
  dt2 <- dt1[idx, ]

  # Split the 9th column of dt2 according to separator ";" and get a list
  list <- stringr::str_split(dt2$V9, ";")

  # Extract the character string which contains gene name information
  names <- sapply(list, stringr::str_subset, "gene_name")

  # Remove the "gene_name" characters
  names <- sub(pattern = "gene_name.", replacement = "", names)

  # Remove the quotation marks
  names <- gsub(pattern = '\"', replacement = "", names)

  # Remove the whitespace
  dt2$V9 <- gsub(pattern = "\\s", replacement = "", names)

  return(StandardChromName(dt2[, c(1, 4, 5, 7, 9)]))
}


#' PyrTri
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @export
#' @keywords internal
PyrTri <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' PyrPenta
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @export
#' @keywords internal
PyrPenta <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in string.vec
#'
#' @param string.vec a vector of type character.
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @return A vector of type characters with the reverse complement of
#'   of every string in string.vec.
#' @keywords internal
#' @export
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' RevcSNS96
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @export
RevcSNS96 <- function(mutstring) {
  # TODO (steve) document

  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}

#' RevcDNS144
#'
#' @param mutstring TODO
#'
#' @return TODO
#' @keywords internal
#' @export
RevcDNS144 <- function(mutstring) {
  # TODO (Nanhai) document

  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 2))
  target  <- revc(substr(mutstring, 3, 4))
  return(paste0(context, target))
}

#' Read transcript ranges and strands from a bed format file.
#' Use this one for the new, cut down gff3 file (2018 11 24)
#'
#' @param path Path to the file with the transcript information with 1-based
#'   start end positions of genomic ranges.
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @export
ReadTranscriptRanges <- function(path) {
  d <- utils::read.table(path)
  colnames(d) <- c("chrom", "chromStart", "chromEnd", "strand", "name")
  bed1 <- data.table(d)
  data.table::setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

#' Read transcript ranges and strands from a bed format file.
#' Mostly for testing.
#'
#' @param path Path to the file with the transcript information (in bed format).
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @export
ReadBedTranscriptRanges <- function(path) {
  names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
  bed <- utils::read.table(path, col.names = names, as.is = TRUE)

  # Delete duplicate entries in the BED file
  bed <- dplyr::distinct(bed, chrom, chromStart, chromEnd, strand, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  bed$chromStart <- bed$chromStart + 1L

  bed1 <- data.table(bed)
  data.table::setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

#' Read data from a nucleotide abundance file with 3 base pairs
#'
#' @param path Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A matrix whose row names indicate 32 different types of 3 base pairs
#'   combinations while its column contains the occurrences of each type.
#' @export
ReadAbundance3Bp <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}

#' Read data from a nucleotide abundance file with 4 base pairs
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#' @import data.table
#' @return A matrix whose row names indicate 10 different types of 2 base pairs
#'   combinations while its column contains the occurrences of each type.
#' @export
ReadAbundance4Bp <- function(path) {
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

#' Read data from a nucleotide abundance file with 5 base pairs
#'
#' @param path Path to the file with the nucleotide abundance information
#'   with 5 base pairs.
#' @import data.table
#' @return A matrix whose row names indicate 512 different types of 5 base
#'   pairs combinations while its column contains the occurrences of each type.
#' @export
ReadAbundance5Bp <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  mat <- as.matrix(dt1[, 2])
  rownames(mat) <- dt1[[1]]
  return(mat)
}
