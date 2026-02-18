#' Check SBS mutation class in VCF with the corresponding SBS mutation matrix
#'
#' @param vcf An annotated SBS VCF with columns of SBS mutation
#'   classes added by \code{AddSBSClass}.
#'
#' @param mat The SBS mutation count matrix.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @keywords internal
CheckSBSClassInVCF <- function(vcf, mat, sample.id) {
  if (nrow(mat) %in% c(96, 1536)) {
    # One SBS mutation can be represented by more than 1 row in vcf
    # after annotation by AddTranscript if the mutation position falls in multiple
    # transcripts. When creating the SBS96 and SBS1536 mutation matrix,
    # we only need to count these mutations once.
    df <- dplyr::distinct(vcf, CHROM, ALT, POS, .keep_all = TRUE)

    if (nrow(df) != colSums(mat)) {
      stop(
        "In sample ",
        sample.id,
        ", the number of SBS",
        nrow(mat),
        " variants in the annotated VCF is not the same as the total ",
        "counts in mutation matrix."
      )
    }
  } else {
    # Only keep those mutations that fall within transcribed region
    # when generating SBS192 mutation matrix.
    df1 <- vcf[!is.na(trans.strand), ]

    # Discard variants that fall on transcripts on both strand.
    df2 <- df1[bothstrand == FALSE, ]

    # One SBS mutation can be represented by more than 1 row in df2 if the mutation
    # position falls into the range of multiple transcripts on the same strand. We
    # only need to count these mutations once.
    df3 <- dplyr::distinct(df2, CHROM, ALT, POS, .keep_all = TRUE)
    if (nrow(df3) != colSums(mat)) {
      stop(
        "In sample ",
        sample.id,
        ", the number of SBS",
        nrow(mat),
        " variants in the annotated VCF is not the same as the total ",
        "counts in mutation matrix."
      )
    }
  }
}

#' Add and check SBS class in an annotated VCF with the corresponding SBS
#' mutation matrix
#'
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateSBSVCF}}. It must
#'   *not* contain indels and must *not* contain DBS (double base
#'   substitutions), or triplet base substitutions etc., even if encoded as
#'   neighboring SBS.
#'
#' @param mat96 The SBS96 mutation count matrix.
#'
#' @param mat1536 The SBS1536 mutation count matrix.
#'
#' @param mat192 The SBS192 mutation count matrix.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @return The original \code{vcf} with three additional columns
#'   \code{SBS96.class}, \code{SBS192.class} and \code{SBS1536.class} added.
#'
#' @keywords internal
AddAndCheckSBSClassInVCF <-
  function(vcf, mat96, mat1536, mat192 = NULL, sample.id) {
    vcf1 <- AddSBSClass(vcf)
    CheckSBSClassInVCF(vcf1, mat96, sample.id)
    CheckSBSClassInVCF(vcf1, mat1536, sample.id)
    if (!is.null(mat192)) {
      CheckSBSClassInVCF(vcf1, mat192, sample.id)
    }
    return(vcf1)
  }

#' Check and return the SBS mutation matrix
#'
#' @inheritParams AddAndCheckSBSClassInVCF
#'
#' @param discarded.variants A \code{data.frame} which contains rows of SBS
#'   variants whose pentanucleotide context contains "N".
#'
#' @param return.annotated.vcf Whether to return the annotated VCF with
#'   additional columns showing the mutation class for each variant. Default is
#'   FALSE.
#'
#' @inheritSection CreateOneColSBSMatrix Value
#'
#' @keywords internal
CheckAndReturnSBSMatrix <-
  function(
    vcf,
    discarded.variants,
    mat96,
    mat1536,
    mat192 = NULL,
    return.annotated.vcf = FALSE,
    sample.id = "counts"
  ) {
    if (nrow(discarded.variants) == 0) {
      if (is.null(mat192)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catSBS96 = mat96, catSBS1536 = mat1536))
        } else {
          vcf.SBS.class <-
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(
            catSBS96 = mat96,
            catSBS1536 = mat1536,
            annotated.vcf = vcf.SBS.class
          ))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(
            catSBS96 = mat96,
            catSBS192 = mat192,
            catSBS1536 = mat1536
          ))
        } else {
          vcf.SBS.class <-
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(
            catSBS96 = mat96,
            catSBS192 = mat192,
            catSBS1536 = mat1536,
            annotated.vcf = vcf.SBS.class
          ))
        }
      }
    } else {
      if (is.null(mat192)) {
        if (return.annotated.vcf == FALSE) {
          return(list(
            catSBS96 = mat96,
            catSBS1536 = mat1536,
            discarded.variants = discarded.variants
          ))
        } else {
          vcf.SBS.class <-
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(
            catSBS96 = mat96,
            catSBS1536 = mat1536,
            annotated.vcf = vcf.SBS.class,
            discarded.variants = discarded.variants
          ))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(
            catSBS96 = mat96,
            catSBS192 = mat192,
            catSBS1536 = mat1536,
            discarded.variants = discarded.variants
          ))
        } else {
          vcf.SBS.class <-
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(
            catSBS96 = mat96,
            catSBS192 = mat192,
            catSBS1536 = mat1536,
            annotated.vcf = vcf.SBS.class,
            discarded.variants = discarded.variants
          ))
        }
      }
    }
  }

#' Create the matrix an SBS catalog for *one* sample from an in-memory VCF.
#'
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateSBSVCF}}. It must
#'   *not* contain indels and must *not* contain DBS (double base
#'   substitutions), or triplet base substitutions etc., even if encoded as
#'   neighboring SBS.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @param return.annotated.vcf Whether to return the annotated VCF with
#'   additional columns showing the mutation class for each variant. Default is
#'   FALSE.
#'
#' @import data.table
#'
#' @importFrom dplyr %>% group_by summarize
#'
#' @section Value: A list of three 1-column matrices with the names
#'   \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}. If transcript
#'   information is not available in \code{vcf}, \code{catSBS192} is not
#'   generated. Do not rely on the order of elements in the list. If
#'   \code{return.annotated.vcf} = TRUE, another element \code{annotated.vcf}
#'   will appear in the list. If there are SBS variants whose pentanucleotide
#'   context contains "N", they will be excluded in the analysis and an
#'   additional element \code{discarded.variants} will appear in the return
#'   list.
#'
#' @note catSBS192 only contains mutations in transcribed regions.
#'
#' @keywords internal
CreateOneColSBSMatrix <- function(
  vcf,
  sample.id = "count",
  return.annotated.vcf = FALSE
) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude DBSs
  # before calling the function. This function does not detect DBSs.

  CheckForEmptySBSVCF <- function(
    vcf,
    return.annotated.vcf,
    discarded.variants = NULL
  ) {
    if (0 == nrow(vcf)) {
      # Create 1-column matrix with all values being 0 and the correct row and
      # column labels.
      catSBS96 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$SBS96),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$SBS96, sample.id)
        )
      catSBS192 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$SBS192),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$SBS192, sample.id)
        )
      catSBS1536 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$SBS1536),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$SBS1536, sample.id)
        )

      if (return.annotated.vcf == FALSE) {
        list.to.return <-
          list(
            catSBS96 = catSBS96,
            catSBS192 = catSBS192,
            catSBS1536 = catSBS1536,
            discarded.variants = discarded.variants
          )
      } else {
        list.to.return <-
          list(
            catSBS96 = catSBS96,
            catSBS192 = catSBS192,
            catSBS1536 = catSBS1536,
            annotated.vcf = vcf,
            discarded.variants = discarded.variants
          )
      }
      # If discarded.variants is NULL, then remove this element
      list.to.return <- Filter(Negate(is.null), list.to.return)
      return(list.to.return)
    } else {
      return(FALSE)
    }
  }

  ret1 <- CheckForEmptySBSVCF(
    vcf = vcf,
    return.annotated.vcf = return.annotated.vcf
  )
  if (!is.logical(ret1)) {
    return(ret1)
  }

  stopifnot(nchar(vcf$ALT) == 1)
  stopifnot(nchar(vcf$REF) == 1)
  stopifnot(vcf$ALT != vcf$REF)

  discarded.variants <- vcf[0]
  mismatches <- which(vcf$REF != substr(vcf$seq.21bases, 11, 11))
  if (length(mismatches) != 0) {
    discarded.variants <- rbind(discarded.variants, vcf[mismatches, ])
    discarded.variants$discarded.reason <-
      paste0(
        'SBS variant whose reference base in ref.genome does not match the',
        ' reference base in the VCF file.'
      )
    message(
      "In sample ",
      sample.id,
      " ",
      length(mismatches),
      " row out of ",
      nrow(vcf),
      " had reference base in ref.genome that does not match the ",
      "reference base in the VCF file.\n",
      "Please check the ref.genome argument.\n",
      "See discarded.variants in the return value for more details"
    )
    vcf <- vcf[-mismatches, ]
  }

  # Delete the rows of SBS if the pentanucleotide context contains "N"
  idx <- grep("N", substr(vcf$seq.21bases, 9, 13))
  if (!length(idx) == 0) {
    discarded.variants <- rbind(discarded.variants, vcf[idx, ])
    discarded.variants$discarded.reason <-
      'SBS variant whose pentanucleotide context contains "N"'
    vcf <- vcf[-idx, ]
    warning(
      'Variants in the SBS vcf ',
      sample.id,
      ' whose pentanucleotide context contains "N" ',
      'have been deleted so as not to conflict with downstream processing. ',
      'See discarded.variants in the return value for more details.'
    )
  }

  ret2 <- CheckForEmptySBSVCF(
    vcf = vcf,
    return.annotated.vcf = return.annotated.vcf,
    discarded.variants = discarded.variants
  )
  if (!is.logical(ret2)) {
    return(ret2)
  }

  # Keep a copy of the original vcf
  vcf0 <- vcf

  # Create 2 new columns that show the 3072 and 1536 mutation type
  context <- substr(vcf$seq.21bases, 9, 13)
  vcf$mutation <- paste0(context, vcf$ALT)

  # PyrPenta maps to strand-agnostic category
  # e.g. ATGCT>T "ATGCTT" maps to AGCAT>A, "AGCATA"
  vcf$pyr.mut <- PyrPenta(vcf$mutation)

  # One SBS mutation can be represented by more than 1 row in vcf
  # after annotation by AddTranscript if the mutation position falls
  # in multiple transcripts. When creating the 1536 and 96 catalog,
  # we only need to count these mutations once.
  #vcf1 <- vcf[, .(REF = REF[1], pyr.mut = pyr.mut[1]),
  #            by = .(CHROM, ALT, POS)]
  vcf1 <- vcf %>%
    dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], pyr.mut = pyr.mut[1])

  # Create part of the 1536 catalog matrix but missing mutation
  # types have NA in the count column.
  tab1536 <- table(vcf1[, "pyr.mut"])
  stopifnot(setequal(
    setdiff(names(tab1536), ICAMS::catalog.row.order$SBS1536),
    c()
  ))
  dt1536 <- data.table(tab1536)

  colnames(dt1536) <- c("rn", "count")
  d <- data.table(rn = ICAMS::catalog.row.order$SBS1536)
  stopifnot(length(ICAMS::catalog.row.order$SBS1536) == 1536)
  x <- merge(d, dt1536, by = "rn", all.x = TRUE)
  x[is.na(count), count := 0]
  stopifnot(sum(x$count) == nrow(vcf1))
  mat1536 <- matrix(x$count)
  rownames(mat1536) <- x$rn
  mat1536 <- mat1536[ICAMS::catalog.row.order$SBS1536, , drop = FALSE]
  colnames(mat1536) <- sample.id

  # Create the 96 catalog matrix
  x[, nrn := paste0(substr(rn, 2, 4), substr(rn, 6, 6))]
  dt96 <- x[, sum(count), by = nrn]
  stopifnot(nrow(dt96) == 96)
  mat96 <- matrix(dt96$V1)
  rownames(mat96) <- dt96$nrn
  mat96 <- mat96[ICAMS::catalog.row.order$SBS96, , drop = FALSE]
  colnames(mat96) <- sample.id

  if (is.null(vcf$trans.strand)) {
    retval <-
      CheckAndReturnSBSMatrix(
        vcf = vcf0,
        discarded.variants = discarded.variants,
        mat96 = mat96,
        mat1536 = mat1536,
        mat192 = NULL,
        return.annotated.vcf = return.annotated.vcf,
        sample.id = sample.id
      )
    return(retval)
  }

  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 192 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One SBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 192 catalog, we only need to count these mutations once.
  # vcf3 <- vcf2[, .(REF = REF[1], mutation = mutation[1],
  #                 trans.strand = trans.strand[1]),
  #             by = .(CHROM, ALT, POS)]
  vcf3 <- vcf2 %>%
    dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(
      REF = REF[1],
      mutation = mutation[1],
      trans.strand = trans.strand[1]
    )

  # If vcf3 has empty rows, we will return 1-column SBS192 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    mat192 <-
      matrix(
        0,
        nrow = length(ICAMS::catalog.row.order$SBS192),
        ncol = 1,
        dimnames = list(ICAMS::catalog.row.order$SBS192, sample.id)
      )
    retval <-
      CheckAndReturnSBSMatrix(
        vcf = vcf0,
        discarded.variants = discarded.variants,
        mat96 = mat96,
        mat1536 = mat1536,
        mat192 = mat192,
        return.annotated.vcf = return.annotated.vcf,
        sample.id = sample.id
      )
    return(retval)
  }

  # Create the 192 catalog matrix
  tab192 <- table(
    paste0(substr(vcf3$mutation, 2, 4), substr(vcf3$mutation, 6, 6)),
    vcf3$trans.strand,
    useNA = "ifany"
  )
  stopifnot(sum(tab192) == nrow(vcf3))
  dt192 <- as.data.table(tab192)
  colnames(dt192) <- c("rn", "trans.strand", "count")
  dt192 <- dt192[!is.na(trans.strand)]
  dt192[trans.strand == "-", rn := RevcSBS96(rn)]
  dt192 <- dt192[, .(count = sum(count)), by = rn]
  x192 <- data.table(rn = ICAMS::catalog.row.order$SBS192)
  x <- merge(x192, dt192, by = "rn", all.x = TRUE)
  x[is.na(count), count := 0]
  mat192 <- matrix(x[, count])
  rownames(mat192) <- unlist(x[, 1])
  mat192 <- mat192[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  colnames(mat192) <- sample.id

  CheckAndReturnSBSMatrix(
    vcf = vcf0,
    discarded.variants = discarded.variants,
    mat96 = mat96,
    mat1536 = mat1536,
    mat192 = mat192,
    return.annotated.vcf = return.annotated.vcf,
    sample.id = sample.id
  )
}

#' Add sequence context and transcript information to an in-memory DBS VCF
#'
#' @param DBS.vcf An in-memory DBS VCF as a \code{data.frame}.
#'
#' @inheritParams AnnotateSBSVCF
#'
#' @return An in-memory DBS VCF as a \code{data.table}. This has been annotated
#'   with the sequence context (column name \code{seq.21bases}) and with
#'   transcript information in the form of a gene symbol (e.g. \code{"TP53"})
#'   and transcript strand. This information is in the columns
#'   \code{trans.start.pos}, \code{trans.end.pos} , \code{trans.strand},
#'   \code{trans.Ensembl.gene.ID} and \code{trans.gene.symbol} in the output.
#'   These columns are not added if \code{is.null(trans.ranges)}.
#'
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitVCFs(file, variant.caller = "strelka")
#' DBS.vcf <- list.of.vcfs$DBS[[1]]
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateDBSVCF <- function(
  DBS.vcf,
  ref.genome,
  trans.ranges = NULL,
  name.of.VCF = NULL
) {
  DBS.vcf <- AddSeqContext(
    df = DBS.vcf,
    ref.genome = ref.genome,
    name.of.VCF = name.of.VCF
  )
  CheckSeqContextInVCF(DBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    DBS.vcf <- AddTranscript(
      df = DBS.vcf,
      trans.ranges = trans.ranges,
      ref.genome = ref.genome,
      name.of.VCF = name.of.VCF
    )
  }
  return(as.data.table(DBS.vcf))
}

#' Add DBS mutation class to an annotated DBS VCF
#'
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateDBSVCF}}. It must
#'   *not* contain indels and must *not* contain SBS (single base
#'   substitutions), or triplet base substitutions etc.
#'
#' @return The original \code{vcf} with three additional columns
#'   \code{DBS78.class}, \code{DBS136.class} and \code{DBS144.class} added.
#'
#' @keywords internal
AddDBSClass <- function(vcf) {
  vcf$DBS78.class <- CanonicalizeDBS(vcf$REF, vcf$ALT)
  vcf$DBS136.class <- CanonicalizeQUAD(substr(vcf$seq.21bases, 10, 13))
  vcf$DBS144.class <- NA
  idx <- which(!is.na(vcf$trans.strand) & (vcf$bothstrand == FALSE))
  vcf$DBS144.class[idx] <- paste0(vcf$REF[idx], vcf$ALT[idx])
  idx1 <- which(vcf$trans.strand == "-" & (vcf$bothstrand == FALSE))
  vcf$DBS144.class[idx1] <- RevcDBS144(vcf$DBS144.class[idx1])
  return(vcf)
}

#' Check DBS mutation class in VCF with the corresponding DBS mutation matrix
#'
#' @param vcf An annotated DBS VCF with columns of DBS mutation
#'   classes added by \code{AddDBSClass}.
#'
#' @param mat The DBS mutation count matrix.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @keywords internal
CheckDBSClassInVCF <- function(vcf, mat, sample.id) {
  if (nrow(mat) %in% c(78, 136)) {
    # One DBS mutation can be represented by more than 1 row in vcf
    # after annotation by AddTranscript if the mutation position falls in multiple
    # transcripts. When creating the DBS78 and DBS136 mutation matrix,
    # we only need to count these mutations once.
    df <- dplyr::distinct(vcf, CHROM, ALT, POS, .keep_all = TRUE)

    if (nrow(df) != colSums(mat)) {
      stop(
        "In sample ",
        sample.id,
        ", the number of DBS",
        nrow(mat),
        " variants in the annotated VCF is not the same as the total ",
        "counts in mutation matrix."
      )
    }
  } else {
    # Only keep those mutations that fall within transcribed region
    # when generating DBS144 mutation matrix.
    df1 <- vcf[!is.na(trans.strand), ]

    # Discard variants that fall on transcripts on both strand.
    df2 <- df1[bothstrand == FALSE, ]

    # One DBS mutation can be represented by more than 1 row in df2 if the mutation
    # position falls into the range of multiple transcripts on the same strand. We
    # only need to count these mutations once.
    df3 <- dplyr::distinct(df2, CHROM, ALT, POS, .keep_all = TRUE)
    if (nrow(df3) != colSums(mat)) {
      stop(
        "In sample ",
        sample.id,
        ", the number of DBS",
        nrow(mat),
        " variants in the annotated VCF is not the same as the total ",
        "counts in mutation matrix."
      )
    }
  }
}

#' Add and check DBS class in an annotated VCF with the corresponding DBS
#' mutation matrix
#'
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateDBSVCF}}. It must
#'   *not* contain indels and must *not* contain SBS (single base
#'   substitutions), or triplet base substitutions etc.
#'
#' @param mat78 The DBS78 mutation count matrix.
#'
#' @param mat136 The DBS136 mutation count matrix.
#'
#' @param mat144 The DBS144 mutation count matrix.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @return The original \code{vcf} with three additional columns
#'   \code{DBS78.class}, \code{DBS136.class} and \code{DBS144.class} added.
#'
#' @keywords internal
AddAndCheckDBSClassInVCF <-
  function(vcf, mat78, mat136, mat144 = NULL, sample.id) {
    vcf1 <- AddDBSClass(vcf)
    CheckDBSClassInVCF(vcf1, mat78, sample.id)
    CheckDBSClassInVCF(vcf1, mat136, sample.id)
    if (!is.null(mat144)) {
      CheckDBSClassInVCF(vcf1, mat144, sample.id)
    }
    return(vcf1)
  }

#' Check and return the DBS mutation matrix
#'
#' @inheritParams AddAndCheckDBSClassInVCF
#'
#' @param discarded.variants A \code{data.frame} which contains rows of DBS
#'   variants whose tetranucleotide context contains "N".
#'
#' @param return.annotated.vcf Whether to return the annotated VCF with
#'   additional columns showing the mutation class for each variant. Default is
#'   FALSE.
#'
#' @inheritSection CreateOneColDBSMatrix Value
#'
#' @keywords internal
CheckAndReturnDBSMatrix <-
  function(
    vcf,
    discarded.variants,
    mat78,
    mat136,
    mat144 = NULL,
    return.annotated.vcf = FALSE,
    sample.id = "counts"
  ) {
    if (nrow(discarded.variants) == 0) {
      if (is.null(mat144)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS136 = mat136))
        } else {
          vcf.DBS.class <-
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(
            catDBS78 = mat78,
            catDBS136 = mat136,
            annotated.vcf = vcf.DBS.class
          ))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS144 = mat144, catDBS136 = mat136))
        } else {
          vcf.DBS.class <-
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(
            catDBS78 = mat78,
            catDBS144 = mat144,
            catDBS136 = mat136,
            annotated.vcf = vcf.DBS.class
          ))
        }
      }
    } else {
      if (is.null(mat144)) {
        if (return.annotated.vcf == FALSE) {
          return(list(
            catDBS78 = mat78,
            catDBS136 = mat136,
            discarded.variants = discarded.variants
          ))
        } else {
          vcf.DBS.class <-
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(
            catDBS78 = mat78,
            catDBS136 = mat136,
            annotated.vcf = vcf.DBS.class,
            discarded.variants = discarded.variants
          ))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(
            catDBS78 = mat78,
            catDBS144 = mat144,
            catDBS136 = mat136,
            discarded.variants = discarded.variants
          ))
        } else {
          vcf.DBS.class <-
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(
            catDBS78 = mat78,
            catDBS144 = mat144,
            catDBS136 = mat136,
            annotated.vcf = vcf.DBS.class,
            discarded.variants = discarded.variants
          ))
        }
      }
    }
  }

#' Create the matrix a DBS catalog for *one* sample from an in-memory VCF.
#'
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateDBSVCF}}. It must
#'   *not* contain indels and must *not* contain SBS (single base
#'   substitutions), or triplet base substitutions etc.
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @import data.table
#'
#' @importFrom dplyr %>% group_by summarize
#'
#' @section Value: A list of three 1-column matrices with the names \code{catDBS78},
#'   \code{catDBS136}, and \code{catDBS144}. If trans.ranges is NULL,
#'   \code{catDBS144} is not generated. Do not rely on the order of elements in
#'   the list. If \code{return.annotated.vcf} = TRUE, another element
#'   \code{annotated.vcf} will appear in the list. If there are DBS variants
#'   whose tetranucleotide context contains "N", they will be excluded in the
#'   analysis and an additional element \code{discarded.variants} will appear in
#'   the return list.
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @keywords internal
CreateOneColDBSMatrix <- function(
  vcf,
  sample.id = "count",
  return.annotated.vcf = FALSE
) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude SBSs
  # before calling the function. This function does not detect SBSs.

  CheckForEmptyDBSVCF <- function(
    vcf,
    return.annotated.vcf,
    discarded.variants = NULL
  ) {
    if (0 == nrow(vcf)) {
      # Create 1-column matrix with all values being 0 and the correct row labels.
      catDBS78 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$DBS78),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$DBS78, sample.id)
        )
      catDBS136 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$DBS136),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$DBS136, sample.id)
        )
      catDBS144 <-
        matrix(
          0,
          nrow = length(ICAMS::catalog.row.order$DBS144),
          ncol = 1,
          dimnames = list(ICAMS::catalog.row.order$DBS144, sample.id)
        )
      if (return.annotated.vcf == FALSE) {
        list.to.return <-
          list(
            catDBS78 = catDBS78,
            catDBS136 = catDBS136,
            catDBS144 = catDBS144,
            discarded.variants = discarded.variants
          )
      } else {
        list.to.return <-
          list(
            catDBS78 = catDBS78,
            catDBS136 = catDBS136,
            catDBS144 = catDBS144,
            annotated.vcf = vcf,
            discarded.variants = discarded.variants
          )
      }
      # Remove element discarded variants if it is NULL
      list.to.return <- Filter(Negate(is.null), list.to.return)
      return(list.to.return)
    } else {
      return(FALSE)
    }
  }

  ret1 <- CheckForEmptyDBSVCF(
    vcf = vcf,
    return.annotated.vcf = return.annotated.vcf
  )
  if (!is.logical(ret1)) {
    return(ret1)
  }

  stopifnot(nchar(vcf$ALT) == 2)
  stopifnot(nchar(vcf$REF) == 2)

  discarded.variants <- vcf[0]
  # Delete the rows of DBS if the tetranucleotide context contains "N"
  idx <- grep("N", substr(vcf$seq.21bases, 10, 13))
  if (!length(idx) == 0) {
    discarded.variants <- rbind(discarded.variants, vcf[idx, ])
    discarded.variants$discarded.reason <-
      'DBS variant whose tetranucleotide context contains "N"'
    vcf <- vcf[-idx, ]
    warning(
      'Variants in the DBS vcf ',
      sample.id,
      ' whose tetranucleotide context contains "N" ',
      'have been deleted so as not to conflict with downstream processing. ',
      'See discarded.variants in the return value for more details.'
    )
  }

  ret2 <- CheckForEmptyDBSVCF(
    vcf = vcf,
    return.annotated.vcf = return.annotated.vcf,
    discarded.variants = discarded.variants
  )
  if (!is.logical(ret2)) {
    return(ret2)
  }

  # One DBS mutation can be represented by more than 1 row in vcf after annotated by
  # AnnotateDBSVCF function if the mutation position falls into the range of
  # multiple transcripts. When creating the 78 and 136 catalog, we only need to
  # count these mutations once.
  # vcf1 <- vcf[, .(REF = REF[1], seq.21bases = seq.21bases[1]),
  #            by = .(CHROM, ALT, POS)]
  vcf1 <- vcf %>%
    dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], seq.21bases = seq.21bases[1])

  # Create the 78 DBS catalog matrix
  canon.DBS.78 <- CanonicalizeDBS(vcf1$REF, vcf1$ALT)
  tab.DBS.78 <- table(canon.DBS.78)
  row.order.78 <- data.table(rn = ICAMS::catalog.row.order$DBS78)
  DBS.dt.78 <- as.data.table(tab.DBS.78)

  # DBS.dt.78 has two columns, names canon.DBS.78 (from the table() function)
  # and N (the count)
  DBS.dt.78.2 <-
    merge(
      row.order.78,
      DBS.dt.78,
      by.x = "rn",
      by.y = "canon.DBS.78",
      all = TRUE
    )
  DBS.dt.78.2[is.na(N), N := 0]
  stopifnot(DBS.dt.78.2$rn == ICAMS::catalog.row.order$DBS78)
  DBS.mat.78 <- as.matrix(DBS.dt.78.2[, 2])
  rownames(DBS.mat.78) <- DBS.dt.78.2$rn
  colnames(DBS.mat.78) <- sample.id

  # Create the 136 DBS catalog matrix
  canon.DBS.136 <- CanonicalizeQUAD(substr(vcf1$seq.21bases, 10, 13))
  tab.DBS.136 <- table(canon.DBS.136)
  row.order.136 <- data.table(rn = ICAMS::catalog.row.order$DBS136)
  DBS.dt.136 <- as.data.table(tab.DBS.136)

  # DBS.dt.136 has two columns, names canon.DBS.136 (from the table() function)
  # and N (the count)
  DBS.dt.136.2 <-
    merge(
      row.order.136,
      DBS.dt.136,
      by.x = "rn",
      by.y = "canon.DBS.136",
      all = TRUE
    )
  DBS.dt.136.2[is.na(N), N := 0]
  stopifnot(DBS.dt.136.2$rn == ICAMS::catalog.row.order$DBS136)
  DBS.mat.136 <- as.matrix(DBS.dt.136.2[, 2])
  rownames(DBS.mat.136) <- DBS.dt.136.2$rn
  colnames(DBS.mat.136) <- sample.id

  if (is.null(vcf$trans.strand)) {
    retval <-
      CheckAndReturnDBSMatrix(
        vcf = vcf,
        discarded.variants = discarded.variants,
        mat78 = DBS.mat.78,
        mat136 = DBS.mat.136,
        mat144 = NULL,
        return.annotated.vcf = return.annotated.vcf,
        sample.id = sample.id
      )
    return(retval)
  }

  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 144 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One DBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 144 catalog, we only need to count these mutations once.
  # vcf3 <- vcf2[, .(REF = REF[1], trans.strand = trans.strand[1]),
  #              by = .(CHROM, ALT, POS)]
  vcf3 <- vcf2 %>%
    dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], trans.strand = trans.strand[1])

  # If vcf3 has empty rows, we will return 1-column DBS144 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    DBS.mat.144 <-
      matrix(
        0,
        nrow = length(ICAMS::catalog.row.order$DBS144),
        ncol = 1,
        dimnames = list(ICAMS::catalog.row.order$DBS144, sample.id)
      )
    retval <-
      CheckAndReturnDBSMatrix(
        vcf = vcf,
        discarded.variants = discarded.variants,
        mat78 = DBS.mat.78,
        mat136 = DBS.mat.136,
        mat144 = DBS.mat.144,
        return.annotated.vcf = return.annotated.vcf,
        sample.id = sample.id
      )
    return(retval)
  }

  # Create the 144 DBS catalog matrix
  # There are 144 stranded DBSs: 4 X 4 sources and 3 X 3 alternates;
  # 4 x 4 x 3 x 3 = 144.
  tab.DBS.144 <-
    table(paste0(vcf3$REF, vcf3$ALT), vcf3$trans.strand, useNA = "ifany")
  stopifnot(sum(tab.DBS.144) == nrow(vcf3))
  DBS.dt.144 <- as.data.table(tab.DBS.144)
  colnames(DBS.dt.144) <- c("rn", "trans.strand", "count")
  DBS.dt.144 <- DBS.dt.144[!is.na(trans.strand)]
  DBS.dt.144[trans.strand == "-", rn := RevcDBS144(rn)]
  DBS.dt.144 <- DBS.dt.144[, .(count = sum(count)), by = rn]
  row.order.144 <- data.table(rn = ICAMS::catalog.row.order$DBS144)

  # DBS.dt.144 has two columns, names rn and count
  DBS.dt.144.2 <- merge(row.order.144, DBS.dt.144, by = "rn", all.x = TRUE)
  DBS.dt.144.2[is.na(count), count := 0]
  stopifnot(DBS.dt.144.2$rn == ICAMS::catalog.row.order$DBS144)
  DBS.mat.144 <- as.matrix(DBS.dt.144.2[, 2])
  rownames(DBS.mat.144) <- DBS.dt.144.2$rn
  colnames(DBS.mat.144) <- sample.id

  CheckAndReturnDBSMatrix(
    vcf = vcf,
    discarded.variants = discarded.variants,
    mat78 = DBS.mat.78,
    mat136 = DBS.mat.136,
    mat144 = DBS.mat.144,
    return.annotated.vcf = return.annotated.vcf,
    sample.id = sample.id
  )
}

#' \strong{\[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka") instead\]}
#' Create SBS and DBS catalogs from Strelka SBS VCF files and plot them to PDF
#'
#' \strong{\[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka") instead\]}
#' Create 3 SBS catalogs (96, 192, 1536) and 3 DBS catalogs (78, 136, 144) from
#' the Strelka SBS VCFs specified by \code{files} and plot them to PDF. The
#' function will find and merge adjacent SBS pairs into DBS if their VAFs are
#' very similar. The default threshold value for VAF is 0.02.
#'
#' This function calls \code{\link{StrelkaSBSVCFFilesToCatalog}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @section Value:
#' A list containing the following objects:
#'
#' * \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}: Matrix of
#' 3 SBS catalogs (one each for 96, 192, and 1536).
#'
#' * \code{catDBS78}, \code{catDBS136}, \code{catDBS144}: Matrix of
#' 3 DBS catalogs (one each for 78, 136, and 144).
#'
#' * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#' that were excluded from the analysis. See the added extra column
#' \code{discarded.reason} for more details.
#'
#' * \code{annotated.vcfs}:
#' \strong{Non-NULL only if} \code{return.annotated.vcfs} = TRUE.
#' A list of elements:
#'     + \code{SBS}: SBS VCF annotated by \code{\link{AnnotateSBSVCF}} with
#'     three new columns \code{SBS96.class}, \code{SBS192.class} and
#'     \code{SBS1536.class} showing the mutation class for each SBS variant.
#'     + \code{DBS}: DBS VCF annotated by \code{\link{AnnotateDBSVCF}} with
#'     three new columns \code{DBS78.class}, \code{DBS136.class} and
#'     \code{DBS144.class} showing the mutation class for each DBS variant.
#'
#' If \code{trans.ranges} is not provided by user and cannot be inferred by
#' ICAMS, SBS 192 and DBS 144 catalog will not be generated. Each catalog has
#' attributes added. See \code{\link{as.catalog}} for more details.
#' @md
#'
#' @section Note: SBS 192 and DBS 144 catalogs include only mutations in
#'   transcribed regions.
#'
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#'
#' @export
#'
#' @examples
#' \dontrun{
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <-
#'     StrelkaSBSVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19",
#'                                             trans.ranges = trans.ranges.GRCh37,
#'                                             region = "genome",
#'                                             output.file =
#'                                             file.path(tempdir(), "StrelkaSBS"))}
#'}
StrelkaSBSVCFFilesToCatalogAndPlotToPdf <-
  function(
    files,
    ref.genome,
    trans.ranges = NULL,
    region = "unknown",
    names.of.VCFs = NULL,
    output.file = "",
    return.annotated.vcfs = FALSE,
    suppress.discarded.variants.warnings = TRUE
  ) {
    lifecycle::deprecate_soft(
      when = "3.0.0",
      what = "StrelkaSBSVCFFilesToCatalogAndPlotToPdf()",
      details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead'
    )
    catalogs0 <-
      StrelkaSBSVCFFilesToCatalog(
        files,
        ref.genome,
        trans.ranges,
        region,
        names.of.VCFs,
        return.annotated.vcfs,
        suppress.discarded.variants.warnings
      )

    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    if (output.file != "") {
      output.file <- paste0(output.file, ".")
    }

    for (name in names(catalogs)) {
      PlotCatalogToPdf(
        catalogs[[name]],
        file = paste0(output.file, name, ".pdf")
      )
      if (name == "catSBS192") {
        PlotCatalogToPdf(
          catalogs[[name]],
          file = paste0(output.file, "SBS12.pdf"),
          plot.SBS12 = TRUE
        )
      }
    }

    return(catalogs)
  }

#' Create ID (small insertions and deletions) catalog from Strelka ID VCF files and plot them to PDF
#'
#' \strong{Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka") instead}
#' Create ID (small insertions and deletions) catalog from the Strelka ID VCFs
#' specified by \code{files} and plot them to PDF
#'
#' This function calls \code{\link{StrelkaIDVCFFilesToCatalog}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @param output.file Optional. The base name of the PDF file to be produced;
#'   the file is ending in \code{catID.pdf}.
#'
#' @inheritSection StrelkaIDVCFFilesToCatalog Value
#'
#' @inheritSection VCFsToCatalogsAndPlotToPdf ID classification
#'
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @examples
#' \dontrun{
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catID <-
#'     StrelkaIDVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19",
#'                                            region = "genome",
#'                                            output.file =
#'                                            file.path(tempdir(), "StrelkaID"))}
#'}
StrelkaIDVCFFilesToCatalogAndPlotToPdf <-
  function(
    files,
    ref.genome,
    region = "unknown",
    names.of.VCFs = NULL,
    output.file = "",
    flag.mismatches = 0,
    return.annotated.vcfs = FALSE,
    suppress.discarded.variants.warnings = TRUE
  ) {
    lifecycle::deprecate_soft(
      when = "3.0.0",
      what = "StrelkaIDVCFFilesToCatalogAndPlotToPdf()",
      details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead'
    )

    list <-
      StrelkaIDVCFFilesToCatalog(
        files,
        ref.genome,
        region,
        names.of.VCFs,
        flag.mismatches,
        return.annotated.vcfs,
        suppress.discarded.variants.warnings
      )

    if (output.file != "") {
      output.file <- paste0(output.file, ".")
    }

    PlotCatalogToPdf(list$catalog, file = paste0(output.file, "catID", ".pdf"))

    return(list)
  }

#' \strong{\[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect") instead\]}
#' Create SBS, DBS and Indel catalogs from Mutect VCF files
#' and plot them to PDF
#'
#' \strong{\[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect") instead\]}
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the Mutect VCFs specified by \code{files} and plot them to
#' PDF
#'
#' This function calls \code{\link{MutectVCFFilesToCatalog}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link[BSgenome]{BSgenome}} object
#'   \enumerate{
#'     \item \code{BSgenome.Hsapiens.1000genomes.hs37d5}
#'     \item \code{BSgenome.Hsapiens.UCSC.hg38}
#'     \item \code{BSgenome.Mmusculus.UCSC.mm10}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param names.of.VCFs Optional. Character vector of names of the VCF files.
#'   The order of names in \code{names.of.VCFs} should match the order of VCF
#'   file paths in \code{files}. If \code{NULL}(default), this function will
#'   remove all of the path up to and including the last path separator (if any)
#'   in \code{files} and file paths without extensions (and the leading dot)
#'   will be used as the names of the VCF files.
#'
#' @param tumor.col.names Optional. Vector of column names or column indices in
#'   VCFs which contain the tumor sample information. The order of elements in
#'   \code{tumor.col.names} should match the order of VCFs specified in
#'   \code{files}. If \code{tumor.col.names} is equal to \code{NA}(default),
#'   this function will use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'
#' @param output.file Optional. The base name of the PDF files to be produced;
#'   multiple files will be generated, each ending in \eqn{x}\code{.pdf}, where
#'   \eqn{x} indicates the type of catalog plotted in the file.
#'
#' @param flag.mismatches Deprecated. If there are ID variants whose \code{REF}
#'   do not match the extracted sequence from \code{ref.genome}, the function
#'   will automatically discard these variants and an element
#'   \code{discarded.variants} will appear in the return value. See
#'   \code{\link{AnnotateIDVCF}} for more details.
#'
#' @param return.annotated.vcfs Logical. Whether to return the annotated VCFs
#'   with additional columns showing mutation class for each variant. Default is
#'   FALSE.
#'
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @section Value:
#' A list containing the following objects:
#'
#' * \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}: Matrix of
#' 3 SBS catalogs (one each for 96, 192, and 1536).
#'
#' * \code{catDBS78}, \code{catDBS136}, \code{catDBS144}: Matrix of
#' 3 DBS catalogs (one each for 78, 136, and 144).
#'
#' * \code{catID}: Matrix of ID (small insertions and deletions) catalog.
#'
#' * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#' that were excluded from the analysis. See the added extra column
#' \code{discarded.reason} for more details.
#'
#' * \code{annotated.vcfs}:
#' \strong{Non-NULL only if} \code{return.annotated.vcfs} = TRUE.
#' A list of elements:
#'     + \code{SBS}: SBS VCF annotated by \code{\link{AnnotateSBSVCF}} with
#'     three new columns \code{SBS96.class}, \code{SBS192.class} and
#'     \code{SBS1536.class} showing the mutation class for each SBS variant.
#'     + \code{DBS}: DBS VCF annotated by \code{\link{AnnotateDBSVCF}} with
#'     three new columns \code{DBS78.class}, \code{DBS136.class} and
#'     \code{DBS144.class} showing the mutation class for each DBS variant.
#'     + \code{ID}: ID VCF annotated by \code{\link{AnnotateIDVCF}} with one
#'     new column \code{ID.class} showing the mutation class for each
#'     ID variant.
#'
#' If \code{trans.ranges} is not provided by user and cannot be inferred by
#' ICAMS, SBS 192 and DBS 144 catalog will not be generated. Each catalog has
#' attributes added. See \code{\link{as.catalog}} for more details.
#' @md
#'
#' @inheritSection VCFsToCatalogsAndPlotToPdf ID classification
#'
#' @section Note:
#'  SBS 192 and DBS 144 catalogs include only mutations in transcribed regions.
#'  In ID (small insertions and deletions) catalogs, deletion repeat sizes range
#'  from 0 to 5+, but for plotting and end-user documentation deletion repeat
#'  sizes range from 1 to 6+.
#'
#' @section Comments:
#' To add or change attributes of the catalog, you can use function
#' \code{\link[base]{attr}}. \cr For example, \code{attr(catalog, "abundance")
#' <- custom.abundance}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <-
#'     MutectVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19",
#'                                         trans.ranges = trans.ranges.GRCh37,
#'                                         region = "genome",
#'                                         output.file =
#'                                         file.path(tempdir(), "Mutect"))}
#'}
MutectVCFFilesToCatalogAndPlotToPdf <-
  function(
    files,
    ref.genome,
    trans.ranges = NULL,
    region = "unknown",
    names.of.VCFs = NULL,
    tumor.col.names = NA,
    output.file = "",
    flag.mismatches = 0,
    return.annotated.vcfs = FALSE,
    suppress.discarded.variants.warnings = TRUE
  ) {
    lifecycle::deprecate_soft(
      when = "3.0.0",
      what = "MutectVCFFilesToCatalogAndPlotToPdf()",
      details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect")` instead'
    )

    catalogs0 <-
      MutectVCFFilesToCatalog(
        files,
        ref.genome,
        trans.ranges,
        region,
        names.of.VCFs,
        tumor.col.names,
        flag.mismatches,
        return.annotated.vcfs,
        suppress.discarded.variants.warnings
      )

    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    if (output.file != "") {
      output.file <- paste0(output.file, ".")
    }

    for (name in names(catalogs)) {
      PlotCatalogToPdf(
        catalogs[[name]],
        file = paste0(output.file, name, ".pdf")
      )
      if (name == "catSBS192") {
        PlotCatalogToPdf(
          catalogs[[name]],
          file = paste0(output.file, "SBS12.pdf"),
          plot.SBS12 = TRUE
        )
      }
    }

    return(catalogs0)
  }

#' Create SBS, DBS and Indel catalogs from VCFs and plot them to PDF
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the VCFs specified by \code{files} and plot them to
#' PDF
#'
#' This function calls \code{\link{VCFsToCatalogs}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @param output.dir The directory where the PDF files will be saved.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param variant.caller Name of the variant caller that produces the VCF, can
#'   be either \code{"strelka"}, \code{"mutect"}, \code{"freebayes"} or
#'   \code{"unknown"}. This information is needed to calculate the VAFs (variant
#'   allele frequencies). If variant caller is \code{"unknown"}(default) and
#'   \code{get.vaf.function} is NULL, then VAF and read depth will be NAs. If
#'   variant caller is \code{"mutect"}, do \strong{not} merge SBSs into DBS.
#'
#' @param num.of.cores The number of cores to use. Not available on Windows
#'   unless \code{num.of.cores = 1}.
#'
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link[BSgenome]{BSgenome}} object
#'   \enumerate{
#'     \item \code{BSgenome.Hsapiens.1000genomes.hs37d5}
#'     \item \code{BSgenome.Hsapiens.UCSC.hg38}
#'     \item \code{BSgenome.Mmusculus.UCSC.mm10}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param names.of.VCFs Optional. Character vector of names of the VCF files.
#'   The order of names in \code{names.of.VCFs} should match the order of VCF
#'   file paths in \code{files}. If \code{NULL}(default), this function will
#'   remove all of the path up to and including the last path separator (if any)
#'   in \code{files} and file paths without extensions (and the leading dot)
#'   will be used as the names of the VCF files.
#'
#' @param tumor.col.names Optional. Only applicable to \strong{Mutect} VCFs.
#'   Vector of column names or column indices in \strong{Mutect} VCFs which
#'   contain the tumor sample information. The order of elements in
#'   \code{tumor.col.names} should match the order of \strong{Mutect} VCFs
#'   specified in \code{files}. If \code{tumor.col.names} is equal to
#'   \code{NA}(default), this function will use the 10th column in all the
#'   \strong{Mutect} VCFs to calculate VAFs. See \code{\link{GetMutectVAF}} for
#'   more details.
#'
#' @param filter.status The character string in column \code{FILTER} of the VCF
#'   that indicates that a variant has passed all the variant caller's filters.
#'   Variants (lines in the VCF) for which the value in column \code{FILTER}
#'   does not equal \code{filter.status} are silently excluded from the output.
#'   The internal function \code{DefaultFilterStatus} tries to infer
#'   \code{filter.status} based on \code{variant.caller}. If
#'   \code{variant.caller} is "unknown", user must specify \code{filter.status}
#'   explicitly. If \code{filter.status = NULL}, all variants are retained. If
#'   there is no \code{FILTER} column in the VCF, all variants are retained with
#'   a warning.
#'
#' @param get.vaf.function Optional. Only applicable when \code{variant.caller} is
#' \strong{"unknown"}. Function to calculate VAF(variant allele frequency) and read
#'   depth information from original VCF. See \code{\link{GetMutectVAF}} as an example.
#'   If \code{NULL}(default) and \code{variant.caller} is "unknown", then VAF
#'   and read depth will be NAs.
#'
#' @param ... Optional arguments to \code{get.vaf.function}.
#'
#' @param max.vaf.diff \strong{Not} applicable if \code{variant.caller =
#'   "mutect"}. The maximum difference of VAF, default value is 0.02. If the
#'   absolute difference of VAFs for adjacent SBSs is bigger than
#'   \code{max.vaf.diff}, then these adjacent SBSs are likely to be "merely"
#'   asynchronous single base mutations, opposed to a simultaneous doublet
#'   mutation or variants involving more than two consecutive bases. Use
#'   negative value (e.g. -1) to suppress merging adjacent SBSs to DBS.
#'
#' @param base.filename Optional. The base name of the PDF files to be produced;
#'   multiple files will be generated, each ending in \eqn{x}\code{.pdf}, where
#'   \eqn{x} indicates the type of catalog plotted in the file.
#'
#' @param return.annotated.vcfs Logical. Whether to return the annotated VCFs
#'   with additional columns showing mutation class for each variant. Default is
#'   FALSE.
#'
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @param chr.names.to.process A character vector specifying the chromosome
#'   names in VCF whose variants will be kept and processed, other chromosome
#'   variants will be discarded. If NULL(default), all variants will be kept
#'   except those on chromosomes with names that contain strings "GL", "KI",
#'   "random", "Hs", "M", "JH", "fix", "alt".
#'
#' @section Value:
#' A list containing the following objects:
#'
#' * \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}: Matrix of
#' 3 SBS catalogs (one each for 96, 192, and 1536).
#'
#' * \code{catDBS78}, \code{catDBS136}, \code{catDBS144}: Matrix of
#' 3 DBS catalogs (one each for 78, 136, and 144).
#'
#' * \code{catID}: Matrix of ID (small insertions and deletions) catalog.
#'
#' * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#' that were excluded from the analysis. See the added extra column
#' \code{discarded.reason} for more details.
#'
#' * \code{annotated.vcfs}:
#' \strong{Non-NULL only if} \code{return.annotated.vcfs} = TRUE.
#' A list of elements:
#'     + \code{SBS}: SBS VCF annotated by \code{\link{AnnotateSBSVCF}} with
#'     three new columns \code{SBS96.class}, \code{SBS192.class} and
#'     \code{SBS1536.class} showing the mutation class for each SBS variant.
#'     + \code{DBS}: DBS VCF annotated by \code{\link{AnnotateDBSVCF}} with
#'     three new columns \code{DBS78.class}, \code{DBS136.class} and
#'     \code{DBS144.class} showing the mutation class for each DBS variant.
#'     + \code{ID}: ID VCF annotated by \code{\link{AnnotateIDVCF}} with one
#'     new column \code{ID.class} showing the mutation class for each
#'     ID variant.
#'
#' If \code{trans.ranges} is not provided by user and cannot be inferred by
#' ICAMS, SBS 192 and DBS 144 catalog will not be generated. Each catalog has
#' attributes added. See \code{\link{as.catalog}} for more details.
#' @md
#'
#' @section ID classification:
#' See \url{https://github.com/steverozen/ICAMS/blob/v3.0.9-branch/data-raw/PCAWG7_indel_classification_2021_09_03.xlsx}
#' for additional information on ID (small insertions and deletions) mutation
#' classification.
#'
#' See the documentation for \code{\link{Canonicalize1Del}} which first handles
#' deletions in homopolymers, then handles deletions in simple repeats with
#' longer repeat units, (e.g. \code{CACACACA}, see
#' \code{\link{FindMaxRepeatDel}}), and if the deletion is not in a simple
#' repeat, looks for microhomology (see \code{\link{FindDelMH}}).
#'
#' @section Note:
#'  SBS 192 and DBS 144 catalogs include only mutations in transcribed regions.
#'  In ID (small insertions and deletions) catalogs, deletion repeat sizes range
#'  from 0 to 5+, but for plotting and end-user documentation deletion repeat
#'  sizes range from 1 to 6+.
#'
#' @section Comments:
#' To add or change attributes of the catalog, you can use function
#' \code{\link[base]{attr}}. \cr For example, \code{attr(catalog, "abundance")
#' <- custom.abundance}.
#'
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <-
#'     VCFsToCatalogsAndPlotToPdf(file, ref.genome = "hg19",
#'                                output.dir = tempdir(),
#'                                variant.caller = "mutect",
#'                                region = "genome",
#'                                base.filename = "Mutect")}
VCFsToCatalogsAndPlotToPdf <-
  function(
    files,
    output.dir,
    ref.genome,
    variant.caller = "unknown",
    num.of.cores = 1,
    trans.ranges = NULL,
    region = "unknown",
    names.of.VCFs = NULL,
    tumor.col.names = NA,
    filter.status = DefaultFilterStatus(variant.caller),
    get.vaf.function = NULL,
    ...,
    max.vaf.diff = 0.02,
    base.filename = "",
    return.annotated.vcfs = FALSE,
    suppress.discarded.variants.warnings = TRUE,
    chr.names.to.process = NULL
  ) {
    num.of.cores <- AdjustNumberOfCores(num.of.cores)

    catalogs0 <-
      VCFsToCatalogs(
        files = files,
        ref.genome = ref.genome,
        variant.caller = variant.caller,
        num.of.cores = num.of.cores,
        trans.ranges = trans.ranges,
        region = region,
        names.of.VCFs = names.of.VCFs,
        tumor.col.names = tumor.col.names,
        filter.status = filter.status,
        get.vaf.function = get.vaf.function,
        ... = ...,
        max.vaf.diff = max.vaf.diff,
        return.annotated.vcfs = return.annotated.vcfs,
        suppress.discarded.variants.warnings = suppress.discarded.variants.warnings,
        chr.names.to.process = chr.names.to.process
      )

    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    if (base.filename != "") {
      base.filename <- paste0(base.filename, ".")
    }

    for (name in names(catalogs)) {
      non.empty.samples <- RetrieveNonEmptySamples(catalogs[[name]])
      # Only plot samples which have mutations for a specific mutation class
      if (!is.null(non.empty.samples)) {
        PlotCatalogToPdf(
          non.empty.samples,
          file = file.path(output.dir, paste0(base.filename, name, ".pdf"))
        )
        if (name == "catSBS192") {
          PlotCatalogToPdf(
            non.empty.samples,
            file = file.path(output.dir, paste0(base.filename, "SBS12.pdf")),
            plot.SBS12 = TRUE
          )
        }
      }
    }

    return(catalogs0)
  }

#' @keywords internal
RetrieveNonEmptySamples <- function(catalog) {
  tmp <- colSums(catalog)
  indices <- which(tmp > 0)
  if (length(indices) > 0) {
    return(catalog[, indices, drop = FALSE])
  } else {
    return(NULL)
  }
}

#' @keywords internal
CanonicalizeDBS <- function(ref.vec, alt.vec) {
  DBS <- paste0(ref.vec, alt.vec)
  idx <- which(!(DBS %in% ICAMS::catalog.row.order$DBS78))
  if (length(idx) == 0) {
    return(DBS)
  } else {
    out <- paste0(fastrc::fast_rc(ref.vec[idx]), fastrc::fast_rc(alt.vec[idx]))
    stopifnot(all(out %in% ICAMS::catalog.row.order$DBS78))
    DBS[idx] <- out
    return(DBS)
  }
}

#' @keywords internal
CanonicalizeQUAD <- function(quad) {
  idx <- which(!(quad %in% ICAMS::catalog.row.order$DBS136))
  if (length(idx) == 0) {
    return(quad)
  } else {
    out <- fastrc::fast_rc(quad[idx])
    stopifnot(all(out %in% ICAMS::catalog.row.order$DBS136))
    quad[idx] <- out
    return(quad)
  }
}

#' @keywords internal
CheckNamesOfVCFs <- function(files, names.of.VCFs) {
  stopifnot(inherits(names.of.VCFs, "character"))
  if (length(files) != length(names.of.VCFs)) {
    stop(
      "\nThe number of names in names.of.VCFs does not match ",
      "the number of VCF files"
    )
  }
}

#' @keywords internal
InferTransRanges <- function(ref.genome, trans.ranges = NULL) {
  if (!is.null(trans.ranges)) {
    return(trans.ranges)
  } else {
    if (IsGRCh37(ref.genome)) {
      return(ICAMS::trans.ranges.GRCh37)
    } else if (IsGRCh38(ref.genome)) {
      return(ICAMS::trans.ranges.GRCh38)
    } else if (IsGRCm38(ref.genome)) {
      return(ICAMS::trans.ranges.GRCm38)
    } else {
      return(trans.ranges)
    }
  }
}
