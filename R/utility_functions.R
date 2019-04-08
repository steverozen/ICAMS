#' "Collapse" a catalog.
#'
#' "Collapse" a catalog. Do not use this function for
#' signature catalogs.
#'
#' \code{Collapse192To96} Collapse an SNS 192 catalog
#' to an SNS 96 catalog.
#'
#' \code{Collapse1536To96} Collapse an SNS 1536 catalog
#'  to an SNS 96 catalog.
#'
#' \code{Collapse144To78} Collapse a DNS 144 catalog
#' to a DNS 78 catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
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
  mat96 <- mat96[ICAMS::catalog.row.order$SNS96, , drop = FALSE]
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
  mat96 <- mat96[ICAMS::catalog.row.order$SNS96, , drop = FALSE]
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
  mat78 <- mat78[ICAMS::catalog.row.order$DNS78, , drop = FALSE]
}

#' Transform between count and density catalogs
#' and signatures.
#'
#' @details Only the following transformations are legal:
#'
#' \enumerate{
#' \item \code{counts -> counts}
#' \item \code{counts -> density}
#' \item \code{counts -> (counts.signature, density.signature)}
#' \item \code{density -> counts} (in which case the semantics are to
#' infer the genome-wide or exome-wide counts based on the
#' densities.)
#' \item \code{density -> (counts.signature, density.signature)}
#' \item \preformatted{(counts.signature, density.signature) ->
#'  (counts.signature, density.signature)}
#' (\code{density.signature -> density.signature} is a null operation.)
#' \item \code{density -> density} (A null operation.)
#' }
#'
#'
#' @param catalog An SNS or DNS catalog as described in \code{\link{ICAMS}};
#'  must \strong{not} be an ID (indel) catalog.
#'
#' @param target.ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param target.region One of "genome", "exome".
#'
#' @param target.catalog.type A character string acting as a catalog type
#'   identifier, one of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#'
#' @export
TransformCatalog <-
  function(catalog, target.ref.genome, target.region, target.catalog.type) {
    # Some error checking
    stopifnot(target.catalog.type %in% c("counts", "density",
                                         "counts.signature", "density.signature"))

    if (attributes(catalog)$type %in% c("counts.signature", "density.signature") &&
        !target.catalog.type %in% c("counts.signature", "density.signature")) {
      stop("Only a \"counts\" or \"density\" type catalog ",
           "can be transformed to a different type.")
    }

    if (attributes(catalog)$type == "density.signature" &&
        target.catalog.type == "density.signature") {
      return(catalog)
    }

    if (attributes(catalog)$type == "density" && target.catalog.type == "density") {
      return(catalog)
    }

    if (!nrow(catalog) %in% c(96, 192, 1536, 78, 136, 144)) {
      stop("This function can only transform catalogs from the type of ",
           "SNS96, SNS192, SNS1536, DNS78, DNS136, DNS144")
    }

    if (nrow(catalog) == 192) {
      if (attributes(catalog)$type != "counts" ||
          target.catalog.type %in% c("density", "counts", "density.signature")) {
        stop('For SNS 192 catalog, only transformation from "counts" to "counts.signature" ',
             'is implemented at the current stage.\n')
      } else {
        cat <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
        return(as.catalog(cat, target.ref.genome, target.region, target.catalog.type))
      }
    }

    if (nrow(catalog) == 144) {
      if (attributes(catalog)$type != "counts" ||
          target.catalog.type %in% c("density", "counts", "density.signature")) {
        stop('For DNS 144 catalog, only transformation from "counts" to "counts.signature" ',
             'is implemented at the current stage.\n')
      } else {
        cat <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
        return(as.catalog(cat, target.ref.genome, target.region, target.catalog.type))
      }
    }

    source.abundance <- attributes(catalog)$abundance
    cat <- CreateCatalogAbundance(catalog, target.ref.genome,
                                  target.region, target.catalog.type)
    target.abundance <- attributes(cat)$abundance
    stopifnot(names(source.abundance) == names(target.abundance))

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
    }

    lapply(names(source.abundance), transform.n.mer)

    if (target.catalog.type %in% c("counts.signature", "density.signature")) {
      out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
      return(as.catalog(out2, target.ref.genome,
                        target.region, target.catalog.type))
    } else {
      return(as.catalog(out.catalog, target.ref.genome,
                        target.region, target.catalog.type))
    }
  }

#' Standardize the Chromosome name annotations for a data frame.
#'
#' @param df A data frame whose first column contains the Chromosome name
#'
#' @return A data frame whose Chromosome names are only in the form of 1:22, "X"
#'   and "Y".
#'
#' @keywords internal
StandardChromName <- function(df) {
  # Is there any row in df whose Chromosome names have "GL"?
  if (sum(grepl("GL", df[[1]])) > 0) {
    df <- df[-grep("GL", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names have "KI"?
  if (sum(grepl("KI", df[[1]])) > 0) {
    df <- df[-grep("KI", df[[1]]), ]
  }

  # Is there any row in df whose Chromosome names have "random"?
  if (sum(grepl("random", df[[1]])) > 0) {
    df <- df[-grep("random", df[[1]]), ]
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
#' @return A data table which contains chromosome name, start, end position,
#'   strand information and gene name. It is keyed by chrom, chromStart, and
#'   chromEnd. Only the following four gene types are kept to facilitate
#'   transcriptional strand bias analysis: protein_coding, retained_intron,
#'   processed_transcript and nonsense_mediated_decay.
#'
#' @keywords internal
CreateTransRange <- function(path) {
  df <- read.csv(path, header = FALSE, fill = TRUE, nrows = 20)
  # Count the number of comment lines
  n <- sum(grepl("#", df[, 1]))

  # Read in the raw GFF3 File while skipping the comment lines
  dt <- data.table::fread(path, header = FALSE, sep = "\t", fill = TRUE, skip = n)

  dt1 <- dt[dt$V3 == "gene", ]

  # Select out the four gene types for transcriptional strand bias analysis
  idx <-
    grepl("protein_coding", dt1$V9) |
    grepl("retained_intron", dt1$V9) |
    grepl("processed_transcript", dt1$V9) |
    grepl("nonsense_mediated_decay", dt1$V9)
  dt2 <- dt1[idx, ]

  # Split the 9th column of dt2 according to separator ";" and get a list
  list <- stringi::stri_split_fixed(dt2$V9, ";")

  # Extract the character string which contains gene name information
  names <- sapply(list, stringi::stri_subset_fixed, "gene_name")

  # Remove the "gene_name" characters
  names <- sub(pattern = "gene_name.", replacement = "", names)

  # Remove the quotation marks
  names <- gsub(pattern = '\"', replacement = "", names)

  # Remove the whitespace
  dt2$V9 <- gsub(pattern = "\\s", replacement = "", names)

  # Select the necessary columns and standardize the chromosome names
  dt3 <- StandardChromName(dt2[, c(1, 4, 5, 7, 9)])

  colnames(dt3) <- c("chrom", "chromStart", "chromEnd", "strand", "name")
  chrOrder <-c((1:22), "X", "Y")
  dt3$chrom <- factor(dt3$chrom, chrOrder, ordered = TRUE)
  return(data.table::setkeyv(dt3, c("chrom", "chromStart", "chromEnd")))
}

#' @keywords internal
PyrTri <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1, 3)),
                  revc(substr(mutstring, 4, 4))),
           mutstring)
  return(output)
}

#' @keywords internal
PyrPenta <- function(mutstring) {
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
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#'
#' @return A vector of type characters with the reverse complement of every
#'   string in \code{string.vec}.
#'
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
#'
#' @keywords internal
ReadTranscriptRanges <- function(path) {
  df <- utils::read.csv(path)
  dt <- data.table(df)
  chrOrder <-c((1:22), "X", "Y")
  dt$chrom <- factor(dt$chrom, chrOrder, ordered = TRUE)
  data.table::setkeyv(dt, c("chrom", "chromStart", "chromEnd"))
  return(dt)
}

#' Read chromosome and position information from a bed format file.
#'
#' @param path Path to the file in bed format.
#'
#' @return A data.table keyed by chrom, chromStart, and chromEnd.
#'
#' @keywords internal
ReadBedRanges <- function(path) {
  df <- data.table::fread(path)
  df1 <- StandardChromName(df[, 1:3])
  colnames(df1) <- c("chrom", "chromStart", "chromEnd")

  # Delete duplicate entries in the BED file
  df2 <- dplyr::distinct(df1, chrom, chromStart, chromEnd, .keep_all = TRUE)

  # Bed file are 0 based start and 1 based end (an oversimplification).
  # We need to add 1L and not 1, otherwise the column turns to a double
  # we get a warning from data.table.
  df2$chromStart <- df2$chromStart + 1L

  dt <- data.table(df2)
  chrOrder <- c((1:22), "X", "Y")
  dt$chrom <- factor(dt$chrom, chrOrder, ordered = TRUE)
  data.table::setkeyv(dt, c("chrom", "chromStart", "chromEnd"))
  return(dt)
}

#' Create trinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 32 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
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

#' Create stranded trinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 3
#'   base pairs.
#'
#' @return A numeric vector whose names indicate 64 different types of 3 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedTrinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("3bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create dinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 10 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateDinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("2bp", "occurrences")
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  dt$type <-
    ifelse((dt[[1]]) %in% canonical.ref, dt[[1]], revc(dt[[1]]))
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create stranded dinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 2
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 16 different types of 2 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
#' @keywords internal
CreateStrandedDinucAbundance <- function(path) {
  dt <- fread(path)
  colnames(dt) <- c("2bp", "occurrences")
  dt$type <- dt[[1]]
  dt1 <- dt[, .(counts = sum(occurrences)), by = type]
  abundance <- dt1$counts
  names(abundance) <- dt1$type
  return(abundance)
}

#' Create tetranucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 4
#'   base pairs.
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 136 different types of 4 base pairs
#'   combinations while its values indicate the occurrences of each type.
#'
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
#'
#' @import data.table
#'
#' @return A numeric vector whose names indicate 512 different types of 5 base
#'   pairs combinations while its values indicate the occurrences of each type.
#'
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
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @return If \code{genome} is \link[BSgenome]{BSgenome} object, return it.
#' Otherwise return the \link[BSgenome]{BSgenome} object identified by the
#' string \code{genome}.
#'
#' @keywords internal
NormalizeGenomeArg <- function(ref.genome) {
  if (class(ref.genome) == "character") {
    if (ref.genome %in% c("GRCh38", "hg38")) {
      ref.genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (ref.genome %in% c("GRCh37", "hg19")) {
      ref.genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else {
      stop("Unrecoginzed genome identifier:\n", ref.genome,
           "\nNeed one of GRCh38, hg38, GRCh37, hg19")
    }
  }
  return(ref.genome)
}

#' Check attributes of catalog specified by user
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @return TRUE
#'
#' @keywords internal
CheckCatalogAttribute <- function(ref.genome, region, catalog.type) {
  if (!ref.genome %in% c("GRCh37", "hg19", "GRCh38", "hg38",
                         "BSgenome.Hsapiens.UCSC.hg38",
                         "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    stop("Unrecoginzed reference genome identifier: ", ref.genome,
         "\nNeed one of GRCh38, hg38, GRCh37, hg19",
         "BSgenome.Hsapiens.UCSC.hg38, ",
         "BSgenome.Hsapiens.1000genomes.hs37d5")
  }
  if (!region %in% c("genome", "exome", "transcript")) {
    stop("Unrecoginzed region identifier: ", region,
         "\nNeed one of genome, exome, transcript")
  }
  if (!catalog.type %in% c("counts", "density",
                           "counts.signature", "density.signature")) {
    stop("Unrecoginzed catalog type identifier: ", catalog.type,
         "\nNeed one of counts, density, counts.signature, density.signature")
  }
  return(TRUE)
}

#' Check the class of catalog from path
#'
#' @param path Path to a catalog on disk in the standardized format.
#'
#' @return An object with the corresponding class type of catalog.
#'
#' @keywords internal
CheckClassOfCatalogFromPath <- function(path) {
  cos <- data.table::fread(path)
  if (nrow(cos) == 96) {
    structure("ClassofCatalog", class = "SNS96")
  } else if (nrow(cos) == 192) {
    structure("ClassofCatalog", class = "SNS192")
  } else if (nrow(cos) == 1536) {
    structure("ClassofCatalog", class = "SNS1536")
  } else if (nrow(cos) == 78) {
    structure("ClassofCatalog", class = "DNS78")
  } else if (nrow(cos) == 144) {
    structure("ClassofCatalog", class = "DNS144")
  } else if (nrow(cos) == 136) {
    structure("ClassofCatalog", class = "DNS136")
  } else if (nrow(cos) == 83) {
    structure("ClassofCatalog", class = "ID")
  } else {
    stop("The catalog seems not to be a standard catalog supported by ICAMS",
         "number of rows is ", nrow(cos))
  }
}

#' Create the class attribute of a catalog
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return The original catalog with class attribute added.
#'
#' @keywords internal
CreateCatalogClass <- function(catalog) {
  if(!nrow(catalog) %in% c(96, 192, 1536, 78, 144, 136, 83)) {
    stop('This is not a catalog supported by ICAMS. The input catalog must
         be one type of "SNS96", "SNS192", "SNS1536", "DNS78", "DNS144",
         "DNS136", "ID(indel)"',
         'The number of rows of the input catalog is ', nrow(catalog))
  }
  if(nrow(catalog) == 96) {
    class(catalog) <- append(class(catalog), "SNS96Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 192) {
    class(catalog) <- append(class(catalog), "SNS192Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 1536) {
    class(catalog) <- append(class(catalog), "SNS1536Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 78) {
    class(catalog) <- append(class(catalog), "DNS78Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 144) {
    class(catalog) <- append(class(catalog), "DNS144Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 136) {
    class(catalog) <- append(class(catalog), "DNS136Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  if(nrow(catalog) == 83) {
    class(catalog) <- append(class(catalog), "IndelCatalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  return(catalog)
}

#' Create the abundance attribute of a catalog
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @return The original catalog with abundance attribute added.
#'
#' @keywords internal
CreateCatalogAbundance <- function(catalog, ref.genome, region, catalog.type) {
  if(!nrow(catalog) %in% c(96, 192, 1536, 78, 144, 136, 83)) {
    stop('This is not a catalog supported by ICAMS. The input catalog must
         be one type of "SNS96", "SNS192", "SNS1536", "DNS78", "DNS144",
         "DNS136", "ID(indel)"',
         'The number of rows of the input catalog is ', nrow(catalog))
  }

  if(nrow(catalog) == 96) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.3bp.flat
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.3bp.genome.GRCh37
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.3bp.exome.GRCh37
      }
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.3bp.genome.GRCh38
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.3bp.exome.GRCh38
      }
    }
  }

  if(nrow(catalog) == 192) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.3bp.flat.stranded
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
        attr(catalog, "abundance") <- abundance.3bp.stranded.GRCh37
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
        attr(catalog, "abundance") <- abundance.3bp.stranded.GRCh38
    }
  }

  if(nrow(catalog) == 1536) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.5bp.flat
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.5bp.genome.GRCh37
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.5bp.exome.GRCh37
      }
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.5bp.genome.GRCh38
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.5bp.exome.GRCh38
      }
    }
  }

  if(nrow(catalog) == 78) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.2bp.flat
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.2bp.genome.GRCh37
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.2bp.exome.GRCh37
      }
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.2bp.genome.GRCh38
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.2bp.exome.GRCh38
      }
    }
  }

  if(nrow(catalog) == 144) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.2bp.flat.stranded
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
        attr(catalog, "abundance") <- abundance.2bp.stranded.GRCh37
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
        attr(catalog, "abundance") <- abundance.2bp.stranded.GRCh38
    }
  }

  if(nrow(catalog) == 136) {
    if (catalog.type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.4bp.flat
      return(catalog)
    } else if (ref.genome %in%
               c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.4bp.genome.GRCh37
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.4bp.exome.GRCh37
      }
    } else if (ref.genome %in%
               c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
      if (region == "genome") {
        attr(catalog, "abundance") <- abundance.4bp.genome.GRCh38
      } else if (region == "exome") {
        attr(catalog, "abundance") <- abundance.4bp.exome.GRCh38
      }
    }
  }

  if(nrow(catalog) == 83) {
    attr(catalog, "abundance") <- NULL
  }

  return(catalog)
}

#' Preserve attributes of the input catalog
#'
#' @param pre.catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param new.catalog A new catalog which needs to inherit the necessary attributes from
#' \code{pre.catalog}.
#'
#' @return The new catalog that has inherited the necessary attributes from
#' \code{pre.catalog}
#'
#' @export
PreserveCatalogAttribute <- function(pre.catalog, new.catalog) {
  attr(new.catalog, "ref.genome") <- attributes(pre.catalog)$ref.genome
  attr(new.catalog, "region") <- attributes(pre.catalog)$region
  attr(new.catalog, "type") <- attributes(pre.catalog)$type
  attr(new.catalog, "abundance") <- attributes(pre.catalog)$abundance
  attr(new.catalog, "class") <- attributes(pre.catalog)$class
  return(new.catalog)
}

#' Create attributes of a catalog
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @return The original catalog with the following attributes added: ref.genome,
#'   region, type, abundance, class.
#'
#' @export
as.catalog <- function(catalog, ref.genome, region, catalog.type) {
  ref.genome <- NormalizeGenomeArg(ref.genome)@pkgname

  if (CheckCatalogAttribute(ref.genome, region, catalog.type)) {
    attr(catalog, "ref.genome") <- ref.genome
    attr(catalog, "region") <- region
    attr(catalog, "type") <- catalog.type
    catalog <- CreateCatalogAbundance(catalog, ref.genome, region, catalog.type)
    catalog <- CreateCatalogClass(catalog)
  }
}

#' Generate all possible k-mers of length k.
#'
#' @param k Length of kmers (k>=2)
#'
#' @return Character vector containing all possible k-mers.
#'
#' @keywords internal
GenerateKmer <- function(k) {
  base <- c("A", "C", "G", "T")
  list.of.base <- list()
  for (i in 1:k){
    list.of.base[[i]] <- base
  }
  permutation <- expand.grid(list.of.base, stringsAsFactors = FALSE)
  all.kmer.list <- character(nrow(permutation))
  for (i in 1:k) {
    all.kmer.list <- stringi::stri_c(all.kmer.list, permutation[[i]])
  }
  all.kmer.list <- stringi::stri_sort(all.kmer.list)
}

#' Generate an empty matrix of k-mer abundance
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return An empty matrix of k-mer abundance
#'
#' @keywords internal
GenerateEmptyKmerCounts <- function(k) {
  all.kmer.list <- GenerateKmer(k)
  kmer.counts <- matrix(0, nrow = length(all.kmer.list))
  rownames(kmer.counts) <- all.kmer.list
  colnames(kmer.counts) <- "Freq"
  return(kmer.counts)
}

#' Generate k-mer abundance from given nucleotide sequences
#'
#' @param sequences A vector of nucleotide sequences
#'
#' @param k Length of k-mers (k>=2)
#'
#' @return Matrix of the counts of each k-mer inside \code{sequences}
#'
#' @keywords  internal
GetSequenceKmerCounts <- function(sequences, k) {
  kmer.counts <- GenerateEmptyKmerCounts(k)

  for(start_idx in 1:k){
    temp.seqs <- substring(sequences, start_idx, nchar(sequences))
    temp.kmers <-
      stringi::stri_extract_all_regex(
        temp.seqs, pattern = paste(rep(".", each = k), collapse = ""))
    temp.kmers <- unlist(temp.kmers)
    temp.kmer.counts <- data.frame(table(temp.kmers))
    row.names(temp.kmer.counts) <- temp.kmer.counts[, 1]
    if (any(grepl("N", temp.kmer.counts[, 1]))) {
      temp.kmer.counts <- temp.kmer.counts[-grep("N", temp.kmer.counts[, 1]), ]
    }

    kmer.counts[row.names(temp.kmer.counts), ] <-
      kmer.counts[row.names(temp.kmer.counts), ]  +
      temp.kmer.counts$Freq
  }
  return(kmer.counts)
}

#' Generate k-mer abundance from a given genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simplerepeat masking is accepted now.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @note It will take about 30 minutes to generate k-mer abundance for
#' human reference genome (k <= 5).
#'
#' @return Matrix of the counts of each k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetGenomeKmerCounts <- function(k, ref.genome, filter.path) {
  kmer.counts <- GenerateEmptyKmerCounts(k)
  genome <- NormalizeGenomeArg(ref.genome)

  # Remove decoyed chromosomes and mitochondrial DNA
  chr.list <- seqnames(genome)[which(nchar(GenomeInfoDb::seqnames(genome)) <= 5)]
  if (any(grepl("M", chr.list))) {
    chr.list <- chr.list[-grep("M", chr.list)]
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
  }

  print("Start counting by chromosomes")

  for (idx in 1:length(chr.list)) {
    print(chr.list[idx])

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr.list[idx]), ]
      filter.bed <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      genome.bed <-
        GRanges(chr.list[idx],
                IRanges(1, as.numeric(GenomeInfoDb::seqlengths(genome)[idx])))
      filtered.genome.bed <- GenomicRanges::setdiff(genome.bed, filter.bed)
      genome.seq <- BSgenome::getSeq(genome, filtered.genome.bed,
                                     as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      genome.seq <- gsub(homopolymer.ms.regex.pattern, "N",genome.seq)

    } else {
      genome.seq <- BSgenome::getSeq(genome, chr.list[idx], as.character = TRUE)
    }

    kmer.counts <- kmer.counts + GetSequenceKmerCounts(genome.seq, k)
  }
  return(kmer.counts)
}

#' Generate stranded k-mer abundance from a given genome and gene annotation file
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simplerepeat masking is accepted now.
#'
#' @param trans.range A GFF3 trans.range.file
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @note It will take about 15 minutes to generate stranded k-mer abundance for
#'   human reference genome (k <= 5).
#'
#' @return Matrix of the counts of each stranded k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetStrandedKmerCounts <- function(k, ref.genome, trans.ranges, filter.path) {
  stranded.ranges <- StandardChromName(trans.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)

  # Check whether chromosome names in stranded.ranges are the same as in ref.genome
  if (!all(stranded.ranges$chrom %in% seqnames(genome))){
    stranded.ranges$chrom <- paste0("chr", stranded.ranges$chrom)
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!all(filter.df$V2 %in% seqnames(genome))){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }
  }

  print("Start counting by chromosomes")

  for (chr in unique(stranded.ranges$chrom)) {
    print(chr)
    temp.stranded.ranges <- stranded.ranges[stranded.ranges$chrom == chr, ]
    trans.range.bed <-
      with(temp.stranded.ranges,
           GRanges(chrom, IRanges(chromStart, chromEnd), strand = strand))

    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]
      filter.bed <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))
      filtered.trans.range.bed <- GenomicRanges::setdiff(trans.range.bed, filter.bed)
      stranded.seq <- BSgenome::getSeq(genome, filtered.trans.range.bed,
                                       as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      stranded.seq <- gsub(homopolymer.ms.regex.pattern, "N",stranded.seq)

    }else{
      stranded.seq <- BSgenome::getSeq(genome, filtered.trans.range.bed,
                                       as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(stranded.seq, k)
  }
  return(kmer.counts)
}

#' Generate exome k-mer abundance from a given reference genome
#'
#' @param k Length of k-mers (k>=2)
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param exome.range A keyed data table which has exome ranges information.
#' It has three columns: chrom, chromStart and chromEnd.
#'
#' @param filter.path If given, homopolymers will be masked from
#'   genome(sequence). Only simplerepeat masking is accepted now.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @note It will take about 1 minute to generate exome k-mer abundance for human
#'   reference genome (k <= 5).
#'
#' @return Matrix of the counts of exome k-mer across the \code{ref.genome}
#'
#' @keywords internal
GetExomeKmerCounts <- function(k, ref.genome, exome.ranges, filter.path) {
  exome.ranges <- StandardChromName(exome.ranges)
  genome <- NormalizeGenomeArg(ref.genome)
  kmer.counts <- GenerateEmptyKmerCounts(k)

  # Check whether chromosome names in exome.ranges are the same as in ref.genome
  if (!all(exome.ranges$chrom %in% seqnames(genome))){
    exome.ranges$chrom <- paste0("chr", exome.ranges$chrom)
  }

  if (!missing(filter.path)) {
    filter.df <- fread(filter.path, header = F, stringsAsFactors = F)
    filter.df <- filter.df[filter.df$V6 <= 6]
    filter.df <- StandardChromName(filter.df[, 2:ncol(filter.df)])
    # Check whether chromosome names in filter.df are the same as in ref.genome
    if (!all(filter.df$V2 %in% seqnames(genome))){
      filter.df$V2 <- paste0("chr", filter.df$V2)
    }

  }
  print("Start counting by chromosomes")

  for (chr in unique(exome.ranges$chrom)) {
    print(chr)
    temp.exome.ranges <- exome.ranges[exome.ranges$chrom == chr, ]
    exome.range.bed <-
      with(temp.exome.ranges, GRanges(chrom, IRanges(chromStart, chromEnd)))
    if (!missing(filter.path)) {
      chr.filter.df <- filter.df[which(filter.df$V2 == chr), ]
      filter.bed <- with(chr.filter.df, GRanges(V2, IRanges(V3 + 1, V4)))

      filtered.exome.range.bed <- GenomicRanges::setdiff(exome.range.bed, filter.bed)
      exome.seq <- BSgenome::getSeq(genome, filtered.exome.range.bed,
                                    as.character = TRUE)
      #Filter shorter homopolymer and microsatellites by regex
      exome.seq <- gsub(homopolymer.ms.regex.pattern, "N",exome.seq)

    }else{
      exome.seq <- BSgenome::getSeq(genome, exome.range.bed,
                                    as.character = TRUE)
    }
    kmer.counts <- kmer.counts + GetSequenceKmerCounts(exome.seq, k)
  }
  return(kmer.counts)
}


