#' @include VCF_related_functions.R utilities.R
NULL

#' Create single nucleotide mutation catalog for *one* sample from
#' a Variant Call Format (VCF) file.
#'
#' @param vcf An in-memory VCF file annotated by the AddSequence and
#'   AddTranscript functions. It must *not* contain indels and must *not*
#'   contain DNS (double nucleotide substituions), or triplet base substituions
#'   etc., even if encoded as neighboring SNS.
#' @param sample.id Usually the sample id, but defaults to "count".
#' @include utilities.R
#' @import data.table
#' @return A list of three matrices containing the SNS mutation catalog:
#'   96, 192, 1536 catalog respectively.
#' @export
CreateOneColSNSCatalog <- function(vcf, sample.id = "count") {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude DNSs
  # before calling the function. This function does not detect DNSs.

  if (0 == nrow(vcf)) {
    return(list(cat96 = empty.cats$cat96, cat192 = empty.cats$cat192,
                cat1536 = empty.cats$cat1536))
  }

  stopifnot(nchar(vcf$ALT) == 1)
  stopifnot(nchar(vcf$REF) == 1)

  # Create 2 new columns that show the 3072 and 1536 mutation type
  context <- substr(vcf$seq.21context, 9, 13)
  vcf$mutation <- paste0(context, vcf$ALT)

  # PyrPenta maps to strand-agnostic category
  # eg ATGCT>T "ATGCTT" maps to AGCAT>A, "AGCATA"
  vcf$pyr.mut <- PyrPenta(vcf$mutation)

  # Create part of the 1536 catalog matrix but missing mutation
  # types have NA in the count column.
  tab1536 <- table(vcf[, "pyr.mut"])
  dt1536  <- data.table(tab1536)

  # TODO(steve): document better, but this deals with the case
  # in which not every mutation class was represented in the
  # VCF, in which case we will fill in with 0.
  colnames(dt1536) <- c("rn", "count")
  d <- data.table(rn = .catalog.row.order1536)
  stopifnot(length(.catalog.row.order1536) == 1536)
  x <- merge(d, dt1536, by = "rn", all.x = TRUE)
  x[is.na(count), count := 0]
  stopifnot(sum(x$count) == nrow(vcf))
  mat1536 <- matrix(x$count)
  rownames(mat1536) <- x$rn
  mat1536 <- mat1536[.catalog.row.order1536, , drop = FALSE]
  colnames(mat1536) <- sample.id

  # Create the 96 catalog matrix
  x[, nrn := paste0(substr(rn, 2, 4), substr(rn, 6, 6))]
  dt96 <- x[, sum(count), by = nrn]
  stopifnot(nrow(dt96) == 96)
  mat96 <- matrix(dt96$V1)
  rownames(mat96) <- dt96$nrn
  mat96 <- mat96[.catalog.row.order96, , drop = FALSE]
  colnames(mat96) <- sample.id

  # Create the 192 catalog matrix
  tab192  <- table(paste0(substr(vcf$mutation, 2, 4),
                          substr(vcf$mutation, 6, 6)),
                   vcf[, "strand"],
                   useNA = "ifany")
  stopifnot(sum(tab192) == nrow(vcf))
  dt192 <- as.data.table(tab192)
  colnames(dt192) <- c("rn", "strand", "count")
  dt192 <- dt192[!is.na(strand)]
  dt192[strand == "-", rn := RevcSNS96(rn)]
  dt192 <- dt192[ , .(count = sum(count)), by = rn]
  x192 <- data.table(rn = .catalog.row.order192)
  x <- merge(x192, dt192, by = "rn", all.x = TRUE)
  x[is.na(count), count := 0]
  mat192 <- matrix(x[, count])
  rownames(mat192) <- unlist(x[, 1])
  mat192 <- mat192[.catalog.row.order192, , drop = FALSE]
  colnames(mat192) <- sample.id

  return(list(cat96 = mat96, cat192 = mat192, cat1536 = mat1536))
}

#' Create a list of 3 catalogs (one each for 96, 192, 1536)
#' out of the contents of the VCFs in list.of.vcfs
#'
#' @param list.of.vcfs List vector of in-memory VCFs. The list names will be the
#'   sample ids in the output catalog.
#' @param genome Name of a particular reference genome (without quotations
#'   marks).
#' @param trans.ranges A data frame containing transcript ranges.
#'
#' @return A list of 3 catalogs, one each for 96, 192, 1536:
#'   cat96
#'   cat192
#'   cat1536
#' @export
SNSVCFsToCatalogs <- function(list.of.vcfs, genome, trans.ranges) {
  ncol <- length(list.of.vcfs)

  cat96 <- empty.cats$cat96
  cat192 <- empty.cats$cat192
  cat1536 <- empty.cats$cat1536

  for (i in 1 : ncol) {
    three.vcfs.df <- SplitSNSVCF(list.of.vcfs[[i]])
    SNS <- three.vcfs.df$SNS.vcf

    SNS <- AddSequence(SNS, seq = genome)
    CheckSeqContextInVCF(SNS, "seq.21context")
    SNS <- AddTranscript(SNS, trans.ranges)
    SNS.cat <- CreateOneColSNSCatalog(SNS)
    rm(SNS)
    cat96 <- cbind(cat96, SNS.cat$cat96)
    cat192 <- cbind(cat192, SNS.cat$cat192)
    cat1536 <- cbind(cat1536, SNS.cat$cat1536)
  }

  colnames(cat96) <- names(list.of.vcfs)
  colnames(cat192) <- names(list.of.vcfs)
  colnames(cat1536) <- names(list.of.vcfs)

  return(list(cat96 = cat96, cat192 = cat192, cat1536 = cat1536))
}

#' Collapse a SNS 192 catalog matrix to a 96 catalog matrix
#'
#' @param cat192 A SNS 192 catalog matrix whose row names indicate the 192
#'   mutation types while its columns show the occurrences of each mutation type of
#'   different samples.
#'
#' @return A SNS 96 catalog matrix whose row names indicate the 96
#    mutation types while its columns show the occurrences of
#    each mutation type of different samples.
#' @export
Collapse192to96 <- function(cat192) {
  dt192 <- data.table(cat192)
  dt192$rn <- PyrTri(rownames(cat192))
  dt96 <- dt192[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[ , -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[.catalog.row.order96, , drop = FALSE]
}

#' Collapse a SNS 1536 catalog matrix to a 96 catalog matrix
#'
#' @param cat1536 A SNS 1536 catalog matrix whose row names indicate the 1536
#'   mutation types while its columns show the occurrences of each mutation type
#'   of different samples.
#' @return A SNS 96 catalog matrix whose row names indicate the 96
#'   mutation types while its columns show the occurrences of
#'   each mutation type of different samples.
#' @export
Collapse1536to96 <- function(cat1536) {
  dt <- data.table(cat1536)
  rn <- rownames(cat1536)

  # The next gsub replaces the string representing a
  # single-base mutation in pentanucleotide with the corresponding
  # sring for that mutation in a trinucleotide context.
  dt$rn <- gsub(".(...).(.)", "\\1\\2", rn, perl = TRUE)
  dt96 <- dt[, lapply(.SD, sum), by = rn, .SDcols = ]
  mat96 <- as.matrix(dt96[, -1])
  rownames(mat96) <- dt96$rn
  mat96 <- mat96[.catalog.row.order96, , drop = FALSE]
}
