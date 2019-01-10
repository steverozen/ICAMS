#' @include VCF_related_functions.R utility_functions.R
NULL

#' Create single nucleotide mutation catalog for *one* sample from
#' a Variant Call Format (VCF) file.
#'
#' @param vcf An in-memory VCF file annotated by the AddSequence and
#'   AddTranscript functions. It must *not* contain indels and must *not*
#'   contain DNS (double nucleotide substituions), or triplet base substituions
#'   etc., even if encoded as neighboring SNS.
#' @param sample.id Usually the sample id, but defaults to "count".
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

#' Create double nucleotide catalog for *one* sample from
#' a Variant Call Format (VCF) file
#'
#' @param vcf An in-memory VCF file annotated by the AddSequence and
#'   AddTranscript functions. It must *not* contain indels and must
#'   *not* contain SNS (single nucleotide substituions), or triplet base
#'   substituions etc.
#' @param sample.id Usually the sample id, but defaults to "count".
#' @import data.table
#' @return A list of three matrices containing the DNS catalog:
#'   catDNS78, catDNS144, catQUAD136 respectively.
#' @export
CreateOneColDNSCatalog <- function(vcf, sample.id = "count") {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude SNSs
  # before calling the function. This function does not detect SNSs.

  if (0 == nrow(vcf)) {
    return(list(catDNS78 = empty.cats$catDNS78,
                catDNS144 = empty.cats$catDNS144,
                catQUAD136 = empty.cats$catQUAD136))
  }

  stopifnot(nchar(vcf$ALT) == 2)
  stopifnot(nchar(vcf$REF) == 2)

  # Create the 78 DNS catalog matrix
  canon.DNS.78 <- CanonicalizeDNS(vcf$REF, vcf$ALT)
  tab.DNS.78 <- table(canon.DNS.78)
  row.order.78 <- data.table(rn = .catalog.row.order.DNS.78)
  DNS.dt.78 <- as.data.table(tab.DNS.78)

  # DNS.dt.78 has two columns, names canon.DNS.78 (from the table() function)
  # and N (the count)
  DNS.dt.78.2 <-
    merge(row.order.78, DNS.dt.78,
          by.x = "rn", by.y = "canon.DNS.78", all = TRUE)
  DNS.dt.78.2[is.na(N), N := 0]
  stopifnot(DNS.dt.78.2$rn == .catalog.row.order.DNS.78)
  DNS.mat.78 <- as.matrix(DNS.dt.78.2[, 2])
  rownames(DNS.mat.78) <- DNS.dt.78.2$rn
  colnames(DNS.mat.78)<- sample.id

  # Create the 136 QUAD catalog matrix
  canon.QUAD.136 <- CanonicalizeQUAD(substr(vcf$seq.21context, 10, 13))
  tab.QUAD.136 <- table(canon.QUAD.136)
  row.order.136 <- data.table(rn = .catalog.row.order.QUAD.136)
  QUAD.dt.136 <- as.data.table(tab.QUAD.136)

  # QUAD.dt.136 has two columns, names canon.QUAD.136 (from the table() function)
  # and N (the count)
  QUAD.dt.136.2 <-
    merge(row.order.136, QUAD.dt.136,
          by.x = "rn", by.y = "canon.QUAD.136", all = TRUE)
  QUAD.dt.136.2[is.na(N), N := 0]
  stopifnot(QUAD.dt.136.2$rn == .catalog.row.order.QUAD.136)
  QUAD.mat.136 <- as.matrix(QUAD.dt.136.2[, 2])
  rownames(QUAD.mat.136) <- QUAD.dt.136.2$rn
  colnames(QUAD.mat.136)<- sample.id

  # Create the 144 DNS catalog matrix
  # There are 144 stranded DNSs: 4 X 4 sources and 3 X 3 alternates;
  # 4 x 4 x 3 x 3 = 144.
  tab.DNS.144  <-
    table(paste0(vcf$REF, vcf$ALT), vcf[, "strand"], useNA = "ifany")
  stopifnot(sum(tab.DNS.144) == nrow(vcf))
  DNS.dt.144 <- as.data.table(tab.DNS.144)
  colnames(DNS.dt.144) <- c("rn", "strand", "count")
  DNS.dt.144 <- DNS.dt.144[!is.na(strand)]
  DNS.dt.144[strand == "-", rn := RevcDNS144(rn)]
  DNS.dt.144 <- DNS.dt.144[, .(count = sum(count)), by = rn]
  row.order.144 <- data.table(rn = .catalog.row.order.DNS.144)

  # DNS.dt.144 has two columns, names rn and count
  DNS.dt.144.2 <- merge(row.order.144, DNS.dt.144, by = "rn", all.x = TRUE)
  DNS.dt.144.2[is.na(count), count := 0]
  stopifnot(DNS.dt.144.2$rn == .catalog.row.order.DNS.144)
  DNS.mat.144 <- as.matrix(DNS.dt.144.2[, 2])
  rownames(DNS.mat.144) <- DNS.dt.144.2$rn
  colnames(DNS.mat.144)<- sample.id

  return(list(catDNS78 = DNS.mat.78, catDNS144 = DNS.mat.144,
              catQUAD136 = QUAD.mat.136))
}

#' Create a list of 3 catalogs (one each for DNS78, DNS144 and QUAD136)
#' out of the contents of the VCFs in list.of.vcfs
#'
#' @param list.of.vcfs List vector of in-memory VCFs. The list names will be
#' the sample ids in the output catalog.
#' @param genome Name of a particular reference genome
#' (without quotations marks).
#' @param trans.ranges A data frame containing transcript ranges.
#'
#' @return A list of 3 catalogs, one each for DNS78, DNS144, QUAD136:
#'   catDNS78
#'   catDNS144
#'   catQUAD136
#' @export
DNSVCFsToCatalogs <- function(list.of.vcfs, genome, trans.ranges) {
  ncol <- length(list.of.vcfs)

  catDNS78 <- empty.cats$catDNS78
  catDNS144 <- empty.cats$catDNS144
  catQUAD136 <- empty.cats$catQUAD136

  for (i in 1 : ncol) {
    three.vcfs.df <- SplitSNSVCF(list.of.vcfs[[i]])
    DNS <- three.vcfs.df$DNS.vcf

    DNS <- AddSequence(DNS, seq = genome)
    DNS <- AddTranscript(DNS, trans.ranges)
    CheckSeqContextInVCF(DNS, "seq.21context")
    DNS.cat <- CreateOneColDNSCatalog(DNS)
    rm(DNS)
    catDNS78 <- cbind(catDNS78, DNS.cat$catDNS78)
    catDNS144 <- cbind(catDNS144, DNS.cat$catDNS144)
    catQUAD136 <- cbind(catQUAD136, DNS.cat$catQUAD136)
  }

  colnames(catDNS78) <- names(list.of.vcfs)
  colnames(catDNS144) <- names(list.of.vcfs)
  colnames(catQUAD136) <- names(list.of.vcfs)

  return(list(catDNS78  = catDNS78, catDNS144  = catDNS144,
              catQUAD136  = catQUAD136))
}

#' Create 3 SNS catalogs (96, 192, 1536) and 3 DNS catalogs (78, 136, 144)
#' in the VCFs specified by vector.of.file.paths
#'
#' @param vector.of.file.paths A vector containing the paths of the VCF files.
#' @param genome  Name of a particular reference genome
#'   (without quotations marks).
#' @param trans.ranges A data.table which contains transcript range and
#'   strand information.
#'
#' @return  A list of 3 SNS catalogs (one each for 96, 192, and 1536)
#'   and 3 DNS catalogs (one each for 78, 136, and 144)

#' @export
VCFFiles2Catalog <- function(vector.of.file.paths, genome, trans.ranges) {
  vcfs <- ReadListOfVCFs(vector.of.file.paths)
  return(c(SNSVCFsToCatalogs(vcfs, genome, trans.ranges),
           DNSVCFsToCatalogs(vcfs, genome, trans.ranges)))
}
