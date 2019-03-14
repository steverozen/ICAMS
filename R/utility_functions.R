#' Collapse catalog functions
#'
#' Collapse a catalog matrix
#'
#' \code{Collapse192To96} Collapse a SNS 192 catalog matrix to a SNS 96 catalog matrix.
#'
#' \code{Collapse1536To96} Collapse a SNS 1536 catalog matrix to a SNS 96 catalog matrix.
#'
#' \code{Collapse144To78} Collapse a DNS 144 catalog matrix to a DNS 78 catalog matrix.
#' @param catalog A catalog matrix to be collapsed whose row names indicate the
#'   mutation types while its columns show the occurrences of each mutation
#'   type of different samples.
#' @return A canonical catalog matrix whose row names indicate the mutation
#'   types while its columns show the occurrences of each mutation type of
#'   different samples.
#' @name CollapseCatalog
NULL

#' @rdname CollapseCatalog
#' @export
Collapse192To96 <- function(catalog) {
  dt192 <- data.table(catalog)
  dt192$rn <- PyrTri(rownames(catalog))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[catalog.row.order$SNS96, , drop = FALSE]
}

#' @rdname CollapseCatalog
#' @export
Collapse1536To96 <- function(catalog) {
  dt <- data.table(catalog)
  rn <- rownames(catalog)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # sring for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[catalog.row.order$SNS96, , drop = FALSE]
}

#' @rdname CollapseCatalog
#' @export
Collapse144To78 <- function(catalog) {
  dt144 <- data.table(catalog)
  ref <- substr(rownames(catalog), 1, 2)
  alt <- substr(rownames(catalog), 3, 4)
  dt144$rn <- CanonicalizeDNS(ref, alt)
  dt78 <- dt144[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat78 <- as.matrix(dt78[ , -1])
  rownames(mat78) <- dt78$rn
  mat78 <- mat78[catalog.row.order$DNS78, , drop = FALSE]
}

#' Handle abundance (opportunity) specifications uniformly.
#'
#' @param abundance Either an abundance variable or string specifying an abundance.
#'
#' @param which.n The n for the n-mers, one of 2, 3, 4, 5 for 2-mers, 3-mers, etc.
#' @keywords internal
NormalizeAbundanceArg <- function(abundance, which.n) {
  if (class(abundance) %in% c("integer", "numeric")) {
    stopifnot(nchar(names(abundance)[1]) == which.n)
    return (abundance)
  }
  if (!which.n %in% 2:5) {
    stop("Argument which.n must be in 2:5, got ", which.n)
  }

  if (is.null(abundance)) {
    if (which.n == 2) {
      abundance = rep(1L, 10)
      names(abundance) <- names(ICAMS:::abundance.2bp.exome.GRCh37)
      return(abundance)
    } else if (which.n == 3) {
      abundance = rep(1L, 32)
      names(abundance) <- names(ICAMS:::abundance.3bp.exome.GRCh37)
      return(abundance)
    } else if (which.n == 4) {
      abundance = rep(1L, 136)
      names(abundance) <- names(ICAMS:::abundance.4bp.exome.GRCh37)
      return(abundance)
    } else if (which.n == 5) {
      abundance = rep(1L, 512)
      names(abundance) <- names(ICAMS:::abundance.5bp.exome.GRCh37)
      return(abundance)
    }
  }

  if (!abundance %in% c("GRCh37.genome", "GRCh37.exome",
                        "GRCh38.genome", "GRCh38.exome")) {
    stop ('abundance must be either an abundance matrix created by yourself
          or one of ("GRCh37.genome", "GRCh37.exome", "GRCh38.genome",
          "GRCh38.exome"), got ', abundance)
  }

  if (abundance == "GRCh37.genome") {
    if (which.n == 2) return(abundance.2bp.genome.GRCh37)
    if (which.n == 3) return(abundance.3bp.genome.GRCh37)
    if (which.n == 4) return(abundance.4bp.genome.GRCh37)
    if (which.n == 5) return(abundance.5bp.genome.GRCh37)
  }
  if (abundance == "GRCh37.exome"){
    if (which.n == 2) return(abundance.2bp.exome.GRCh37)
    if (which.n == 3) return(abundance.3bp.exome.GRCh37)
    if (which.n == 4) return(abundance.4bp.exome.GRCh37)
    if (which.n == 5) return(abundance.5bp.exome.GRCh37)
  }
  if (abundance == "GRCh38.genome") {
    if (which.n == 2) return(abundance.2bp.genome.GRCh38)
    if (which.n == 3) return(abundance.3bp.genome.GRCh38)
    if (which.n == 4) return(abundance.4bp.genome.GRCh38)
    if (which.n == 5) return(abundance.5bp.genome.GRCh38)
  }
  if (abundance == "GRCh38.exome"){
    if (which.n == 2) return(abundance.2bp.exome.GRCh38)
    if (which.n == 3) return(abundance.3bp.exome.GRCh38)
    if (which.n == 4) return(abundance.4bp.exome.GRCh38)
    if (which.n == 5) return(abundance.5bp.exome.GRCh38)
  }
  stop("Programming error: we should never get here")
}

#' Transform between count and density catalogs
#'  and signatures and between different
#' source sequence abundances.
#'
#' @details Only certain transformations are legal.
#' \enumerate{
#'
#' \item The type \code{"density"} must always be associated with a \code{NULL}
#' abundance.
#'
#' \item The type \code{"signature"} can be associated with a \code{NULL}
#' abundance.
#'
#' \item The type \code{"counts"} must \strong{not} be associated with
#'  the \code{NULL} abundance.
#'
#' \item Otherwise, the following are legal:
#' \enumerate{
#' \item \code{counts -> counts}
#' \item \code{counts -> density}
#' \item \code{counts -> signature}
#' \item \code{density -> counts} (in which case the semantics are to
#' infer the genome-wide or exome wide counts based on the
#' densities.)
#' \item \code{density -> signature}
#' \item \code{signature -> signature}
#' }
#'
#' }
#'
#'
#' @param catalog A catalog as described in \code{\link{ICAMS}}.
#'
#' @param source.abundance Either \code{NULL} or a numeric vector with one element
#' for each source sequence for the mutation types in \code{catalog}
#' or a string specifying such a vector, one of \code{"GRCh37.genome"},
#' \code{"GRCh37.exome"}, \code{"GRCh38.genome"}, or \code{"GRCh38.exome"}.
#' This is the
#' abundance upon which the counts, densities, or proportions
#' in \code{catalog} are based. For example, for SNS in
#' trinucleotide context, e.g. ACT > AGT or TAC > TTC, the source
#' sequences are ACT and TAC.
#'
#' @param target.abundance Same possibilities as \code{source.abundance}.
#'
#' @param which.n The length of the source sequences, one of 2:5.
#'
#' @param source.type A character specifying type of the input catalog,
#' one of \code{"counts"}, \code{"signature"} or \code{"density"}.
#'
#' @param target.type A character specifying type of the output catalog,
#' with the same possible values as \code{source.type}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}
#'
#' @export
TransformCatalog <-
  function(catalog,
           source.abundance = NULL,
           target.abundance = NULL, which.n,
           source.type, target.type = source.type) {

  stopifnot(source.type %in% c("counts", "signature", "density"))
  stopifnot(target.type %in% c("counts", "signature", "density"))
  if (source.type == "density" && !is.null(source.abundance)) {
    stop('The type "density" must always be associated with a NULL abundance.')
  }
  if (target.type == "density" && !is.null(target.abundance)) {
    stop('The type "density" must always be associated with a NULL abundance.')
  }

  if (target.type != source.type && source.type == "signature") {
    stop("Only a \"counts\" or \"density\" type catalog ",
         "can be transformed to a different type.")
  }

  if (target.type == source.type && source.type == "density") {
    stop("We cannot transform a \"density\" type catalog ",
         "to \"density\" type catalog as it would require ",
         "uniform source and target abundances.")
  }

  # Some error checking
  if (!nrow(catalog) %in% c(96, 192, 1536, 78, 136, 144)) {
    stop("This function can only transform catalogs from the type of ",
         "SNS96, SNS192, SNS1536, DNS78, DNS136, DNS144")
  }
  if (nrow(catalog) == 96 && which.n != 3) {
    stop("Argument which.n must be 3 for an SNS 96 catalog, got ", which.n)
  }
  if (nrow(catalog) == 192 && which.n != 3) {
    stop("Argument which.n must be 3 for an SNS 192 catalog, got ", which.n)
  }
  if (nrow(catalog) == 1536 && which.n != 5) {
    stop("Argument which.n must be 5 for an SNS 1536 catalog, got ", which.n)
  }
  if (nrow(catalog) == 78 && which.n != 2) {
    stop("Argument which.n must be 2 for a DNS 78 catalog, got ", which.n)
  }
  if (nrow(catalog) == 136 && which.n != 4) {
    stop("Argument which.n must be 4 for a DNS 136 catalog, got ", which.n)
  }
  if (nrow(catalog) == 144 && which.n != 2) {
    stop("Argument which.n must be 2 for a DNS 144 catalog, got ", which.n)
  }

  if (nrow(catalog) == 192) {
    if (source.type != "counts" || target.type %in% c("density", "counts")) {
      stop('For SNS 192 catalog, only transformation from "counts" to "signature" ',
           'is implemented at the current stage.')
    } else {
      return(apply(catalog, MARGIN = 2, function (x) x / sum(x)))
    }
  }

  source.abundance <- NormalizeAbundanceArg(source.abundance, which.n)
  target.abundance <- NormalizeAbundanceArg(target.abundance, which.n)
  stopifnot(all(names(source.abundance) == names(target.abundance)))

  factor <- target.abundance / source.abundance
  names(factor) <- names(target.abundance)
  out.catalog <- catalog

  # CAUTION: this function depends on how mutations are encoded in
  # the row names!
  transform.n.mer <- function(source.n.mer) {
    # For 96 and 192 SNS, source.n.mer is e.g. "ACT" (for the encoding of ACT >
    # AGT as "ACTG"); for SNS1536 the n-mer for AACAG > AATAG is AACAG, in the
    # encoding AACAGT. For DNS78 and DNS144 TGGA represents TG >GA, and the
    # source n-mer is TG. For DNS136, TTGA represents TTGA > TNNA, and the
    # source n-mer is TTGA.

    # First, get the rows with the given source.n.mer
    rows <- grep(paste("^", source.n.mer, sep=''), rownames(out.catalog))
    # Then update those rows using the factor for that source.n.mer

    out.catalog[rows, ] <<- out.catalog[rows, ] * factor[source.n.mer]
    # For density, factor = 1/source.abundance
  }

  # Extract the source sequences from catalog
  n.mers <- substr(rownames(catalog), 1, which.n)

  lapply(unique(n.mers), transform.n.mer)

  if (target.type == "signature") {
    out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
    return(out2)
  } else {
    return(out.catalog)
  }

}

#' Standardize the Chromosome name annotations for a data frame
#'
#' @param df A data frame whose first column contains the Chromosome name
#'
#' @return A data frame whose Chromosome names are only in the form of 1:22, "X"
#'   and "Y".
#' @keywords internal
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
#'
#' @return A data frame which contains chromosome name, start, end position,
#'   strand information and gene name. Only the following four gene types are
#'   kept to facilitate transcriptional strand bias analysis: protein_coding,
#'   retained_intron, processed_transcript and nonsense_mediated_decay.
#' @keywords internal
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
#' @keywords internal
PyrTri <- function(mutstring) {
  # TODO (Steve) document
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
#' @keywords internal
PyrPenta <- function(mutstring) {
  # TODO (Steve) document
  stopifnot(nchar(mutstring) == rep(6, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 3, 3) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 5)),
                  revc(substr(mutstring, 6, 6))),
           mutstring)
  return(output)
}

#' Reverse complement every string in \code{string.vec}.
#'
#' @param string.vec a vector of type character.
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @return A vector of type characters with the reverse complement of every
#'   string in \code{string.vec}.
#' @export
revc <- function(string.vec) {
  return(
    as.character(reverseComplement(DNAStringSet(string.vec)))
  )
}

#' @title Reverse complement strings that represent stranded SNSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded SNSs in trinucleotide context,
#' for example "AATC" represents AAT > ACT mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 3 characters
#' concatenated with the reverse complement of the
#' last character, e.g. "AATC" returns "ATTG".
#'
#' @keywords internal
RevcSNS96 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 3))
  target  <- revc(substr(mutstring, 4, 4))
  return(paste0(context, target))
}

#' @title Reverse complement strings that represent stranded DNSs
#'
#' @param mutstring A vector of 4-character strings representing
#' stranded DNSs, for example "AATC" represents AA > TC mutations.
#'
#' @return Return the vector of
#' reverse complements of the first 2 characters
#' concatenated with the reverse complement of the second
#' 2 characters, e.g. "AATC" returns "TTGA".
#'
#' @keywords internal
RevcDNS144 <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  context <- revc(substr(mutstring, 1, 2))
  target  <- revc(substr(mutstring, 3, 4))
  return(paste0(context, target))
}

#' Read transcript ranges and strands from a gff3 format file.
#' Use this one for the new, cut down gff3 file (2018 11 24)
#'
#' @param path Path to the file with the transcript information with 1-based
#'   start end positions of genomic ranges.
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @keywords internal
ReadTranscriptRanges <- function(path) {
  d <- utils::read.table(path)
  colnames(d) <- c("chrom", "chromStart", "chromEnd", "strand", "name")
  bed1 <- data.table(d)
  data.table::setkeyv(bed1, c("chrom", "chromStart", "chromEnd"))
  return(bed1)
}

#' Read transcript ranges and strands from a bed format file.
#'
#' This function is mostly for testing purpose, may be removed in the future.
#'
#' @param path Path to the file with the transcript information (in bed format).
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#' @export
#' @keywords internal
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

#' Create trinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 32 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#' @keywords internal
CreateTrinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 2) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create dinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#' @import data.table
#' @return A numeric vector whose names indicate 10 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#' @keywords internal
CreateDinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("4bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse(substr(dt[[1]], 2, 3) %in% canonical.ref,
           substr(dt[[1]], 2, 3),
           revc(substr(dt[[1]], 2, 3)))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create tetranucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#' @import data.table
#' @return A numeric vector whose names indicate 136 different types of 4 base pairs
#'   combinations while its values indicate the occurrences of each type.
#' @keywords internal
CreateTetranucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("4bp", "occurrences")
  dt$type <- CanonicalizeQUAD(dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create pentanucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information
#'   with 5 base pairs.
#' @import data.table
#' @return A numeric vector whose names indicate 512 different types of 5 base
#'   pairs combinations while its values indicate the occurrences of each type.
#' @keywords internal
CreatePentanucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("5bp", "occurrences")
  dt$type <-
    ifelse(substr(dt[[1]], 3, 3) %in% c("A", "G"), revc(dt[[1]]), dt[[1]])
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Take strings representing a genome and return the \link[BSgenome]{BSgenome} object.
#'
#' @param genome Either a variable containing a \link[BSgenome]{BSgenome} object
#'   or a character string acting as a genome identifier.
#'
#' @return If \code{genome} is \link[BSgenome]{BSgenome} object, return it.
#' Otherwise return the \link[BSgenome]{BSgenome} object identified by the
#' string \code{genome}.
#'
#' @keywords internal
NormalizeGenomeArg <- function(genome) {
  if (class(genome) == "character") {
    if (genome %in% c("GRCh38", "hg38")) {
      genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (genome %in% c("GRCh37", "hg19")) {
      genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else {
      stop("Unrecoginzed genome identifier:\n", genome,
           "\nNeed one of GRCh38, hg38, GRCh37, hg19")
    }
  }
  return(genome)
}
