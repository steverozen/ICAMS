#' @keywords internal
CheckAndRemoveDiscardedVariants <- function(
  vcf,
  name.of.VCF = NULL,
  chr.names.to.process = NULL
) {
  if (nrow(vcf) == 0) {
    return(list(df = vcf))
  }

  # Create an empty data frame for discarded variants
  discarded.variants <- vcf[0, ]

  # Remove rows with same REF and ALT
  idx <- which(vcf$REF == vcf$ALT)
  if (length(idx) > 0) {
    df.to.remove <- vcf[idx, ]
    df.to.remove$discarded.reason <- "Variant with same REF and ALT"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df.to.remove)
    vcf <- vcf[-idx, ]
  }

  # Remove rows with pound sign
  retval <- RemoveRowsWithPoundSignNew(df = vcf, name.of.VCF = name.of.VCF)
  df1 <- retval$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval$discarded.variants)

  # Remove rows with duplicated CHROM and POS
  retval1 <-
    RemoveRowsWithDuplicatedCHROMAndPOSNew(df = df1, name.of.VCF = name.of.VCF)
  df2 <- retval1$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval1$discarded.variants)

  if (is.null(chr.names.to.process)) {
    # Remove rows with unstandardized chromosome names
    retval2 <- StandardChromNameNew(df = df2, name.of.VCF = name.of.VCF)
    df3 <- retval2$df
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
  } else {
    # Only keep variants that are specified by chr.names.to.process
    retval2 <-
      SelectVariantsByChromName(
        df = df2,
        chr.names.to.process = chr.names.to.process,
        name.of.VCF = name.of.VCF
      )
    df3 <- retval2$df
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
  }

  # VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", df3$ALT, fixed = TRUE)
  if (length(multiple.alt) > 0) {
    warning(
      "VCF ",
      ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
      " has variants with multiple alternative alleles and were ",
      "discarded. See discarded.variants in the return value for more ",
      "details."
    )
    df4 <- df3[-multiple.alt, ]
    df4.to.remove <- df3[multiple.alt, ]
    df4.to.remove$discarded.reason <- "Variant with multiple alternative alleles"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df4.to.remove)
  } else {
    df4 <- df3
  }

  # Remove variants involving three or more nucleotides
  # (e.g. ACT > TGA or AACT > GGTA)
  other.df <- which(nchar(df4$REF) > 2 & nchar(df4$ALT) == nchar(df4$REF))

  if (length(other.df) > 0) {
    warning(
      "VCF ",
      ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
      " has variants involving three or more nucleotides and were ",
      "discarded. See discarded.variants in the return value for more ",
      "details."
    )
    df5 <- df4[-other.df, ]
    df5.to.remove <- df4[other.df, ]
    df5.to.remove$discarded.reason <- "Variant involves three or more nucleotides"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df5.to.remove)
  } else {
    df5 <- df4
  }

  # Remove complex indels
  complex.indels <- which(
    (nchar(df5$REF) > 0) & # exclude indels represented without base preceding the indel
      (nchar(df5$ALT) > 0) &
      (nchar(df5$REF) != nchar(df5$ALT)) &
      (substr(df5$REF, 1, 1) != substr(df5$ALT, 1, 1))
  )
  if (length(complex.indels) > 0) {
    warning(
      "VCF ",
      ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
      " has complex indels and were discarded. See discarded.variants ",
      "in the return value for more details."
    )
    df6 <- df5[-complex.indels, ]
    df6.to.remove <- df5[complex.indels, ]
    df6.to.remove$discarded.reason <- "Complex indel"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df6.to.remove)
  } else {
    df6 <- df5
  }

  # Remove wrong DBS variants that have same base in the same position in REF and ALT
  # (e.g. TA > TT or GT > CT)
  wrong.DBS.type1 <- dplyr::filter(
    df6,
    nchar(REF) == 2,
    nchar(ALT) == 2,
    substr(REF, 1, 1) == substr(ALT, 1, 1)
  )
  wrong.DBS.type2 <- dplyr::filter(
    df6,
    nchar(REF) == 2,
    nchar(ALT) == 2,
    substr(REF, 2, 2) == substr(ALT, 2, 2)
  )
  wrong.DBS <- dplyr::bind_rows(wrong.DBS.type1, wrong.DBS.type2)

  if (nrow(wrong.DBS) > 0) {
    warning(
      "VCF ",
      ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
      " has wrong DBS variants and were discarded. See discarded.variants ",
      "in the return value for more details."
    )
    wrong.DBS.pos <- wrong.DBS$POS
    wrong.DBS$discarded.reason <- "Wrong DBS variant"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, wrong.DBS)
    df7 <- dplyr::filter(df6, !POS %in% wrong.DBS.pos)
  } else {
    df7 <- df6
  }

  # Remove variants which have ambiguous REF bases (not A, C, G, T)
  ambiguous.refs <- which(!substr(df7$REF, 1, 1) %in% c("A", "C", "G", "T"))
  if (length(ambiguous.refs) > 0) {
    warning(
      "VCF ",
      ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
      " has ambiguous REF bases and were discarded. See discarded.variants ",
      "in the return value for more details."
    )
    df8 <- df7[-ambiguous.refs, ]
    df8.to.remove <- df7[ambiguous.refs, ]
    df8.to.remove$discarded.reason <- "Variant has ambiguous REF base"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df8.to.remove)
  } else {
    df8 <- df7
  }

  if (nrow(discarded.variants) == 0) {
    return(list(df = df8))
  } else {
    return(list(df = df8, discarded.variants = discarded.variants))
  }
}
