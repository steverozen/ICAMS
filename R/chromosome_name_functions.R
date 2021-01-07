#' Standardize the chromosome name annotations for a data frame.
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

#' Standardize the chromosome name annotations for a data frame.
#'
#' @param df A data frame whose first column contains the Chromosome name
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @return A \strong{list} with the elements
#' * \code{df} a data frame with variants that had "legal" chromosome
#'   names (see below for illegal chromosome names).
#'   Leading "chr" strings are removed.
#' * \code{discarded.variants}: \strong{Non-NULL only if} there
#'   variants with illegal chromosome names; these are
#'   names that contain the strings "GL", "Hs", "KI", "M", "random".
#' @md
#'
#' @keywords internal
StandardChromNameNew <- function(df, name.of.VCF = NULL) {
  # Create an empty data frame for discarded variants
  discarded.variants <- df[0, ]

  # Is there any row in df whose Chromosome names have "GL"?
  if (sum(grepl("GL", df[[1]])) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", sum(grepl("GL", df[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'GL' and ",
            "were removed. ",
            "See discarded.variants in the return value for more details")
    df1 <- df[-grep("GL", df[[1]]), ]
    df1.to.remove <- df[grep("GL", df[[1]]), ]
    df1.to.remove$discarded.reason <- 'Chromosome name contains "GL"'
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df1.to.remove)
  } else {
    df1 <- df
  }

  # Is there any row in df whose Chromosome names have "KI"?
  if (sum(grepl("KI", df1[[1]])) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", sum(grepl("KI", df1[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'KI' and ",
            "were removed. ",
            "See discarded.variants in the return value for more details")
    df2 <- df1[-grep("KI", df1[[1]]), ]
    df2.to.remove <- df1[grep("KI", df1[[1]]), ]
    df2.to.remove$discarded.reason <- 'Chromosome name contains "KI"'
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df2.to.remove)
  } else {
    df2 <- df1
  }

  # Is there any row in df whose Chromosome names have "random"?
  if (sum(grepl("random", df2[[1]])) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", sum(grepl("random", df2[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'random' and ",
            "were removed. ",
            "See discarded.variants in the return value for more details")
    df3 <- df2[-grep("random", df2[[1]]), ]
    df3.to.remove <- df2[grep("random", df2[[1]]), ]
    df3.to.remove$discarded.reason <- 'Chromosome name contains "random"'
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df3.to.remove)
  } else {
    df3 <- df2
  }

  # Is there any row in df whose Chromosome names are "Hs37D5"?
  if (sum(grepl("^Hs", df3[[1]])) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", sum(grepl("^Hs", df3[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'Hs' and ",
            "were removed. ",
            "See discarded.variants in the return value for more details")
    df4 <- df3[-grep("^Hs", df3[[1]]), ]
    df4.to.remove <- df3[grep("^Hs", df3[[1]]), ]
    df4.to.remove$discarded.reason <- 'Chromosome name contains "Hs"'
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df4.to.remove)
  } else {
    df4 <- df3
  }

  # Is there any row in df whose Chromosome names contain "M"?
  if (sum(grepl("M", df4[[1]])) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", sum(grepl("M", df4[[1]])), " row out of ",
            nrow(df), " had chromosome names that contain 'M' and ",
            "were removed. ",
            "See discarded.variants in the return value for more details")
    df5 <- df4[-grep("M", df4[[1]]), ]
    df5.to.remove <- df4[grep("M", df4[[1]]), ]
    df5.to.remove$discarded.reason <- 'Chromosome name contains "M"'
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df5.to.remove)
  } else {
    df5 <- df4
  }

  if (nrow(discarded.variants) == 0) {
    return(list(df = df5))
  } else {
    return(list(df = df5, discarded.variants = discarded.variants))
  }
}

#' Check and, if possible, correct the chromosome names in a VCF \code{data.frame}.
#'
#' @param vcf.df A VCF as a \code{data.frame}. Check the names in column
#' \code{CHROM}.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param ref.genome The reference genome with the chromosome names to check
#' \code{vcf.df$CHROM} against; must be a Bioconductor
#' \code{BSgenome}, e.g.
#' \code{BSgenome.Hsapiens.UCSC.hg38}.
#'
#' @return If the \code{vcf.df$CHROM} values are correct or
#' can be corrected, then a vector of chromosome names
#' that can be used as a replacement for \code{vcf.df$CHROM}.
#' If the names in \code{vcf.df$CHROM} cannot be made to
#' be consistent with the chromosome names in \code{ref.genome},
#' then \code{stop}.
#'
#' @keywords internal
CheckAndFixChrNames <- function(vcf.df, ref.genome, name.of.VCF = NULL) {
  names.to.check <- unique(vcf.df$CHROM)
  # Check whether the naming of chromosomes in vcf.df is consistent
  if(!sum(grepl("^chr", names.to.check)) %in% c(0, length(names.to.check))) {
    stop("\nNaming of chromosomes in VCF ", dQuote(name.of.VCF),
         " is not consistent: ",
         paste(names.to.check, collapse = " "))
  }

  ref.genome.names <- seqnames(ref.genome)

  not.matched <- setdiff(names.to.check, ref.genome.names)

  # The names match -- we leave well-enough alone
  if (length(not.matched) == 0) return(vcf.df$CHROM)

  vcf.has.chr.prefix <- any(grepl(pattern = "^chr", names.to.check))
  ref.has.chr.prefix <- any(grepl(pattern = "^chr", ref.genome.names))

  new.chr.names <- vcf.df$CHROM
  if (ref.has.chr.prefix && !vcf.has.chr.prefix) {
    names.to.check <- paste0("chr", names.to.check)
    new.chr.names <- paste0("chr", new.chr.names)
    not.matched1 <- setdiff(names.to.check, ref.genome.names)
    if (length(not.matched1) == 0) return(new.chr.names)
  }

  if (!ref.has.chr.prefix && vcf.has.chr.prefix) {
    names.to.check <- gsub("chr", "", names.to.check)
    new.chr.names <- gsub("chr", "", new.chr.names)
    not.matched2 <- setdiff(names.to.check, ref.genome.names)
    if (length(not.matched2) == 0) return(new.chr.names)
  }

  organism <- BSgenome::organism(ref.genome)

  CheckForPossibleMatchedChrName <- function(chr1, chr2) {
    if (chr1 %in% names.to.check) {
      # If chr2 is already in names.to.check, then stop
      if (chr2 %in% names.to.check) {
        stopmessage <- function() {
          stop("\n", chr1, " and ", chr2, " both are chromosome names in VCF ",
               dQuote(name.of.VCF),
               ", which should not be the case for ", organism, ". Please check ",
               "your data or specify the correct ref.genome argument")
        }
        if (vcf.has.chr.prefix) {
          if (grepl(pattern = "^chr", chr1)) {
            stopmessage()
          } else {
            chr1 <- paste0("chr", chr1)
            chr2 <- paste0("chr", chr2)
            stopmessage()
          }
        } else {
          if (!grepl(pattern = "^chr", chr1)) {
            stopmessage()
          } else {
            chr1 <- gsub("chr", "", chr1)
            chr2 <- gsub("chr", "", chr2)
            stopmessage()
          }
        }
      }

      new.chr.names[new.chr.names == chr1] <<- chr2
      names.to.check <- setdiff(names.to.check, chr1)
      names.to.check <<- unique(c(names.to.check, chr2))
    }
  }

  if (organism == "Homo sapiens") {

    # Maybe the problem is that X and Y are encoded as chr23 and chr24
    CheckForPossibleMatchedChrName("chr23", "chrX")
    CheckForPossibleMatchedChrName("chr24", "chrY")

    # Maybe the problem is that X and Y are encoded as 23 and 24
    CheckForPossibleMatchedChrName("23", "X")
    CheckForPossibleMatchedChrName("24", "Y")
  }

  if (organism == "Mus musculus") {

    # Maybe the problem is that X and Y are encoded as chr20 and chr21
    CheckForPossibleMatchedChrName("chr20", "chrX")
    CheckForPossibleMatchedChrName("chr21", "chrY")

    # Maybe the problem is that X and Y are encoded as 20 and 21
    CheckForPossibleMatchedChrName("20", "X")
    CheckForPossibleMatchedChrName("21", "Y")
  }

  not.matched3 <- setdiff(names.to.check, ref.genome.names)
  if (length(not.matched3) == 0) return(new.chr.names)

  stop("\nChromosome names in VCF ", dQuote(name.of.VCF),
       " not in ref.genome for ", organism, ": ",
       # We report the _original_ list of not matched names
       paste(not.matched, collapse = " "))
}

#' Check and, if possible, correct the chromosome names in a trans.ranges \code{data.table}
#'
#' @param trans.ranges A \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param vcf.df A VCF as a \code{data.frame}. Check the names in column
#' \code{CHROM}.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param ref.genome The reference genome with the chromosome names to check
#' \code{vcf.df$CHROM} against; must be a Bioconductor
#' \code{BSgenome}, e.g.
#' \code{BSgenome.Hsapiens.UCSC.hg38}.
#'
#' @return If the \code{vcf.df$CHROM} values are correct or can be corrected,
#'   then a vector of chromosome names that can be used as a replacement for
#'   \code{trans.ranges$chrom}. If the names in \code{vcf.df$CHROM} cannot be
#'   made to be consistent with the chromosome names in
#'   \code{trans.ranges$chrom}, then \code{stop}.
#'
#' @keywords internal
CheckAndFixChrNamesForTransRanges <-
  function(trans.ranges, vcf.df, ref.genome, name.of.VCF = NULL) {
    names.to.check <- as.character(unique(trans.ranges$chrom))

    vcf.chr.names <- as.character(unique(vcf.df$CHROM))
    # Check whether the naming of chromosomes in vcf.df is consistent
    if(!sum(grepl("^chr", vcf.chr.names)) %in% c(0, length(vcf.chr.names))) {
      stop("\nNaming of chromosomes in VCF ", dQuote(name.of.VCF),
           " is not consistent: ",
           paste(vcf.chr.names, collapse = " "))
    }

    not.matched <- setdiff(vcf.chr.names, names.to.check)

    # The names match -- we leave well-enough alone
    if (length(not.matched) == 0) return(trans.ranges$chrom)

    vcf.has.chr.prefix <- any(grepl(pattern = "^chr", vcf.chr.names))
    trans.has.chr.prefix <- any(grepl(pattern = "^chr", names.to.check))

    new.chr.names <- as.character(trans.ranges$chrom)
    if (trans.has.chr.prefix && !vcf.has.chr.prefix) {
      names.to.check <- gsub("chr", "", names.to.check)
      new.chr.names <- gsub("chr", "", names.to.check)

      names.to.check <- paste0("chr", names.to.check)
      new.chr.names <- paste0("chr", new.chr.names)
      not.matched1 <- setdiff(vcf.chr.names, names.to.check)
      if (length(not.matched1) == 0) return(new.chr.names)
    }

    if (!trans.has.chr.prefix && vcf.has.chr.prefix) {
      names.to.check <- paste0("chr", names.to.check)
      new.chr.names <- paste0("chr", new.chr.names)
      not.matched2 <- setdiff(vcf.chr.names, names.to.check)
      if (length(not.matched2) == 0) return(new.chr.names)
    }

    organism <- BSgenome::organism(ref.genome)

    CheckForPossibleMatchedChrName <- function(chr1, chr2) {
      if (chr2 %in% vcf.chr.names) {
        # If chr1 is already in vcf.chr.names, then stop
        if (chr1 %in% vcf.chr.names) {
          stopmessage <- function() {
            stop("\n", chr2, " and ", chr1, " both are chromosome names in VCF ",
                 dQuote(name.of.VCF),
                 ", which should not be the case for ", organism, ". Please check ",
                 "your data or specify the correct ref.genome argument")
          }
          if (vcf.has.chr.prefix) {
            if (grepl(pattern = "^chr", chr2)) {
              stopmessage()
            } else {
              chr1 <- paste0("chr", chr1)
              chr2 <- paste0("chr", chr2)
              stopmessage()
            }
          } else {
            if (!grepl(pattern = "^chr", chr1)) {
              stopmessage()
            } else {
              chr1 <- gsub("chr", "", chr1)
              chr2 <- gsub("chr", "", chr2)
              stopmessage()
            }
          }
        }

        # Update the trans.ranges chromosome names to be consistent with
        # those in vcf.df
        new.chr.names[new.chr.names == chr1] <<- chr2
        names.to.check <- setdiff(names.to.check, chr1)
        names.to.check <<- unique(c(names.to.check, chr2))
      }
    }

    if (organism == "Homo sapiens") {

      # Maybe the problem is that X and Y are encoded as chr23 and chr24
      CheckForPossibleMatchedChrName("chrX", "chr23")
      CheckForPossibleMatchedChrName("chrY", "chr24")

      # Maybe the problem is that X and Y are encoded as 23 and 24
      CheckForPossibleMatchedChrName("X", "23")
      CheckForPossibleMatchedChrName("Y", "24")
    }

    if (organism == "Mus musculus") {

      # Maybe the problem is that X and Y are encoded as chr20 and chr21
      CheckForPossibleMatchedChrName("chrX", "chr20")
      CheckForPossibleMatchedChrName("chrY", "chr21")

      # Maybe the problem is that X and Y are encoded as 20 and 21
      CheckForPossibleMatchedChrName("X", "20")
      CheckForPossibleMatchedChrName("Y", "21")
    }

    not.matched3 <- setdiff(vcf.chr.names, names.to.check)
    if (length(not.matched3) == 0) return(new.chr.names)

    stop("\nChromosome names in VCF ", dQuote(name.of.VCF),
         " not in trans ranges for ", organism, ": ",
         # We report the _original_ list of not matched names
         paste(not.matched, collapse = " "))
  }
