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

#' TODO
#'
#' @param catalog TODO
#' @param target.ref.genome TODO
#' @param target.region TODO
#' @param target.type TODO
#'
#' @return TODO
#' @export
TransformCatalog <- function(catalog, target.ref.genome, target.region, target.type) {
  # Some error checking
  stopifnot(target.type %in% c("counts", "density",
                               "counts.signature", "density.signature"))

  if (attributes(catalog)$type %in% c("counts.signature", "density.signature") &&
      !target.type %in% c("counts.signature", "density.signature")) {
    stop("Only a \"counts\" or \"density\" type catalog ",
         "can be transformed to a different type.")
  }

  if (attributes(catalog)$type == "density" && target.type == "density") {
    stop("We cannot transform a \"density\" type catalog ",
         "to \"density\" type catalog as it would require ",
         "uniform source and target abundances.")
  }

  if (!nrow(catalog) %in% c(96, 192, 1536, 78, 136, 144)) {
    stop("This function can only transform catalogs from the type of ",
         "SNS96, SNS192, SNS1536, DNS78, DNS136, DNS144")
  }

  if (nrow(catalog) == 192) {
    if (attributes(catalog)$type != "counts" ||
        target.type %in% c("density", "counts", "density.signature")) {
      stop('For SNS 192 catalog, only transformation from "counts" to "counts.signature" ',
           'is implemented at the current stage.\n')
    } else {
      cat <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
      return(CreateCatalogAttribute(cat, target.ref.genome, target.region, target.type))
    }
  }

  if (nrow(catalog) == 144) {
    if (attributes(catalog)$type != "counts" ||
        target.type %in% c("density", "counts", "density.signature")) {
      stop('For DNS 144 catalog, only transformation from "counts" to "counts.signature" ',
           'is implemented at the current stage.\n')
    } else {
      cat <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
      return(CreateCatalogAttribute(cat, target.ref.genome, target.region, target.type))
    }
  }

  source.abundance <- attributes(catalog)$abundance
  cat <- CreateCatalogAbundance(catalog, target.ref.genome, target.region, target.type)
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

  if (target.type %in% c("counts.signature", "density.signature")) {
    out2 <- apply(out.catalog, MARGIN = 2, function (x) x / sum(x))
    return(CreateCatalogAttribute(out2, target.ref.genome,
                                  target.region, target.type))
  } else {
    return(CreateCatalogAttribute(out.catalog, target.ref.genome,
                                  target.region, target.type))
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
#'
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
#'
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

#' Create dinucleotide abundance
#'
#' @param path Path to the file with the nucleotide abundance information with 4
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
#' @param ref.genome A genome argument as described in \code{\link{ICAMS}}.
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
#' @param ref.genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param type A character string acting as a catalog type identifier, one of
#' "counts", "density", "signature".
#'
#' @return TRUE
#'
#' @keywords internal
CheckCatalogAttribute <- function(ref.genome, region, type) {
  if (!ref.genome %in% c("GRCh37", "hg19", "GRCh38", "hg38",
                         "BSgenome.Hsapiens.UCSC.hg38",
                         "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    stop("Unrecoginzed reference genome identifier: ", ref.genome,
         "\nNeed one of GRCh38, hg38, GRCh37, hg19",
         "BSgenome.Hsapiens.UCSC.hg38, ",
         "BSgenome.Hsapiens.1000genomes.hs37d5")
  }
  if (!region %in% c("genome", "exome", "transcription")) {
    stop("Unrecoginzed region identifier: ", region,
         "\nNeed one of genome, exome, transcription")
  }
  if (!type %in% c("counts", "density", "counts.signature", "density.signature")) {
    stop("Unrecoginzed catalog type identifier: ", type,
         "\nNeed one of counts, density, counts.signature, density.signature")
  }
  return(TRUE)
}

#' Check the class of catalog
#'
#' @param catalog An S3 object with class "catalog".
#'
#' @return An object with the corresponding class type of catalog.
#'
#' @keywords internal
CheckClassOfCatalog <- function(catalog) {
  if (nrow(catalog) == 96) {
    structure("ClassofCatalog", class = "SNS96")
  } else if (nrow(catalog) == 192) {
    structure("ClassofCatalog", class = "SNS192")
  } else if (nrow(catalog) == 1536) {
    structure("ClassofCatalog", class = "SNS1536")
  } else if (nrow(catalog) == 78) {
    structure("ClassofCatalog", class = "DNS78")
  } else if (nrow(catalog) == 144) {
    structure("ClassofCatalog", class = "DNS144")
  } else if (nrow(catalog) == 136) {
    structure("ClassofCatalog", class = "DNS136")
  } else if (nrow(catalog) == 83) {
    structure("ClassofCatalog", class = "ID")
  } else {
    stop("The catalog seems not to be a standard catalog supported by ICAMS",
         "number of rows is ", nrow(catalog))
  }
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
    class(catalog) <- append(class(catalog), "ID(indel)Catalog", after = 0)
    class(catalog) <- unique(attributes(catalog)$class)
  }
  return(catalog)
}

#' Create the abundance attribute of a catalog
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @param ref.genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param type A character string acting as a catalog type identifier, one of
#' "counts", "density", "counts.signature" and "density.signature".
#'
#' @return The original catalog with abundance attribute added.
#'
#' @keywords internal
CreateCatalogAbundance <- function(catalog, ref.genome, region, type) {
  if(!nrow(catalog) %in% c(96, 192, 1536, 78, 144, 136, 83)) {
    stop('This is not a catalog supported by ICAMS. The input catalog must
         be one type of "SNS96", "SNS192", "SNS1536", "DNS78", "DNS144",
         "DNS136", "ID(indel)"',
         'The number of rows of the input catalog is ', nrow(catalog))
  }

  if(nrow(catalog) == 96) {
    if (type %in% c("density", "density.signature")) {
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
    attr(catalog, "abundance") <- NULL
  }

  if(nrow(catalog) == 1536) {
    if (type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.3bp.flat
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
    if (type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.3bp.flat
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
    attr(catalog, "abundance") <- NULL
  }

  if(nrow(catalog) == 136) {
    if (type %in% c("density", "density.signature")) {
      attr(catalog, "abundance") <- abundance.3bp.flat
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
#' See \code{\link{CreateCatalogAttribute}} for more details.
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
#' @param ref.genome A genome argument as described in \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#'
#' @param type A character string acting as a catalog type identifier, one of
#' "counts", "density", "counts.signature" and "density.signature".
#'
#' @return The original catalog with the following attributes added: ref.genome,
#'   region, type, abundance, class.
#'
#' @export
CreateCatalogAttribute <- function(catalog, ref.genome, region, type) {
    ref.genome <- NormalizeGenomeArg(ref.genome)@pkgname

  if (CheckCatalogAttribute(ref.genome, region, type)) {
    attr(catalog, "ref.genome") <- ref.genome
    attr(catalog, "region") <- region
    attr(catalog, "type") <- type
    catalog <- CreateCatalogAbundance(catalog, ref.genome, region, type)
    catalog <- CreateCatalogClass(catalog)
  }
}
