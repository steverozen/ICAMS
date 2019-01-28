#' @include utility_functions.R
NULL

#' Read in the data lines of a Variant Call Format (VCF) file created by
#'     Strelka version 1
#' @importFrom utils read.csv
#' @param path The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#' @keywords internal
ReadStrelkaVCF <- function(path) {
  df <- read.csv(path, header = FALSE, sep = "\t", quote = "",
                 col.names = paste0("c", 1 : 100), as.is = TRUE)

  # Delete the columns which are totally empty
  df <- df[!sapply(df, function(x) all(is.na(x)))]

  # Delete meta-information lines which start with "##"
  idx <- grep("^##", df[, 1])
  df1 <- df[-idx, ]

  # Extract the names of columns in the VCF file
  names <- c("CHROM", as.character(df1[1, ])[-1])
  df1 <- df1[-1, ]
  colnames(df1) <- names

  df1$POS <- as.integer(df1$POS)
  df1$VAF <- GetStrelkaVAF(df1)
  return(StandardChromName(df1))
}

#' Extract the VAFs (variant allele frequencies) from a VCF created by
#'     Strelka version 1
#'
#' @param strelka.vcf said VCF as a data.frame
#'
#' @return A vector of VAFs, one for each row of strelka.vcf
#' @export
GetStrelkaVAF <-function(strelka.vcf) {
  stopifnot(class(strelka.vcf) == "data.frame")
  if (!("TUMOR" %in% names(strelka.vcf)) ||
      !("FORMAT" %in% names(strelka.vcf))) {
    stop("strelka.vcf does not appear to a Strelka VCF, column names are",
         paste(colnames(strelka.vcf), collapse=" "))
  }
  TUMOR <- strelka.vcf[ , "TUMOR"]
  control <- unique(strelka.vcf[ , "FORMAT"])
  alt     <- strelka.vcf[ , "ALT"]
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  values <- strsplit(TUMOR, split=":", fixed=TRUE)
  vaf <- numeric(nrow(strelka.vcf))
  each.base.col <- c("AU", "CU", "GU", "TU")
  for (i in 1:length(vaf)) {
    row.i <- values[[i]]
    names(row.i) <- colnames
    all.read.counts <- row.i[each.base.col]
    x <- strsplit(all.read.counts, split=",", fixed=TRUE)
    tier1.counts <- lapply(X = x, FUN = function(x) x[1]) # Tier 1 calls
    tier1.counts <- as.numeric(unlist(tier1.counts))
    names(tier1.counts) <- each.base.col
    total.read.count <- sum(tier1.counts)
    alt.count <- tier1.counts[paste0(alt[i], "U")]
    vaf[i] <- alt.count/total.read.count
  }
  return(vaf)
}

#' Read in the data lines of a Variant Call Format (VCF) file created by
#'     MuTect
#' @importFrom utils read.csv
#' @param path The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#' @keywords internal
ReadMutectVCF <- function(path) {
  df <- read.csv(path, header = FALSE, sep = "\t", quote = "",
                 col.names = paste0("c", 1 : 100), as.is = TRUE)

  # Delete the columns which are totally empty
  df <- df[!sapply(df, function(x) all(is.na(x)))]

  # Delete meta-information lines which start with "##"
  idx <- grep("^##", df[, 1])
  df1 <- df[-idx, ]

  # Extract the names of columns in the VCF file
  names <- c("CHROM", as.character(df1[1, ])[-1])
  df1 <- df1[-1, ]
  colnames(df1) <- names

  df1$POS <- as.integer(df1$POS)
  df1$VAF <- GetMutectVAF(df1)
  return(StandardChromName(df1))
}

#' Extract the VAFs (variant allele frequencies) from a VCF created by
#'     MuTect
#'
#' @param mutect.vcf said VCF as a data.frame
#'
#' @return A vector of VAFs, one for each row of mutect.vcf
#' @export
GetMutectVAF <-function(mutect.vcf) {
  stopifnot(class(mutect.vcf) == "data.frame")

  # Select out the column which has the information for F1R2 and F2R1
  info <- mutect.vcf[[10]]

  # Get the raw data by splitting the character string
  raw <- strsplit(info, ":")

  # Define a function to extract the values for F1R2 and F2R1
  Extract <- function(x) as.integer(unlist(strsplit(x[5:6], ",")))

  values <- lapply(raw, FUN = Extract)

  # Define a function to calculate VAF according to F1R2 and F2R1 values
  CalculateVAF <- function(x) sum(x[2], x[4]) / sum(x)

  vaf <- sapply(values, FUN = CalculateVAF)

  return(vaf)
}

#' Add sequence context to a data frame with mutation records
#'
#' @param df An input data frame storing mutation records of a VCF file.
#' @param seq A particular reference genome.
#' @importFrom methods as
#' @importFrom BSgenome getSeq
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @return A data frame with a new column added to the input data frame,
#'     which contains sequence context information.
#' @keywords internal
#' @export
AddSequence <- function(df, seq = BSgenome.Hsapiens.1000genomes.hs37d5) {
  if (0 == nrow(df)) return(df)

  # Create a GRanges object with range width equals to 21
  Ranges <-
    as(data.frame(chrom = df$CHROM, start = df$POS - 10, end = df$POS + 10),
       "GRanges")

  # Extract sequence context from the reference genome
  df <- dplyr::mutate(df,
                      seq.21context = getSeq(seq,
                                             Ranges,
                                             as.character = TRUE))
  return(df)
}

#' Add transcript information to a data frame with mutation records
#'
#' @param df A data frame storing mutation records of a VCF file.
#' @param trans.ranges A data.table with the genomic ranges and
#'     strands of transcripts.
#' @import data.table
#' @return A data frame with new columns added to the input data frame,
#'     which contain the mutated gene's name, range and strand information.
#' @export
#' @keywords internal
AddTranscript <- function(df, trans.ranges) {
  if (nrow(df) == 0) {
    return(df)
  }

  # Find range overlaps between the df and trans.ranges
  df1 <- data.table(df)
  df1[, POS2 := POS]
  dt <- foverlaps(df1, trans.ranges,
                  by.x = c("CHROM", "POS", "POS2"),
                  type = "within", mult = "all")

  # Get the lower and upper bounds of the gene location range
  dt1 <- dt[, .(Start = min(chromStart), End = max(chromEnd),
                Name = name[1], strand = strand[1]),
            by = .(CHROM, ALT, POS)] # Note that is important to have
  # ALT in the by list because in a few cases
  # there are multiple ALT alleles at one POS.

  # Swap gene location according to strand information
  dt2 <- dt1[strand == "-", c("End", "Start") := .(Start, End)]

  return(cbind(df, dt2))
}

#' MakeVCFDNSdf TODO(steve) add average VAF
#'
#' Take DNS ranges and the original VCF and generate a VCF with
#' dinucleotide REF and ALT alleles. The output VCF has minimal columns:
#' just CHROM, POS, ID, REF, ALT.
#'
#' @param DNS.range.df Data frame with columns CHROM, LOW, HIGH
#' @param SNS.vcf.dt TODO
#' @import data.table
#' @return TODO
#' @export
MakeVCFDNSdf <- function(DNS.range.df, SNS.vcf.dt) {
  tmpvcf <- SNS.vcf.dt[ , c("CHROM", "POS", "REF", "ALT")]
  DNS.range.dt <- as.data.table(DNS.range.df)
  tmp1 <- merge(DNS.range.dt, tmpvcf,
                by.x = c("CHROM", "LOW"),
                by.y = c("CHROM", "POS"))
  tmp2 <- merge(tmp1, tmpvcf,
                by.x = c("CHROM", "HIGH"),
                by.y = c("CHROM", "POS"))
  tmp2[, POS := LOW]
  tmp2[, ID := "From merged SNSs"]
  tmp2[, REF := paste0(REF.x, REF.y)]
  tmp2[, ALT := paste0(ALT.x, ALT.y)]
  return(as.data.frame(tmp2[, c("CHROM", "POS", "ID", "REF", "ALT")]))
}

#' Split an in-memory VCF into SNS, DNS, and variants involving
#' > 2 consecutive bases
#'
#' SNSs are single nucleotide substitutions,
#' eg C>T, A<G,....  DNSs are double nucleotide substitutions,
#' eg CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT > TCTA, ...
#' @param vcf.df An in-memory data frame containing a VCF file contents.
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02.
#' @import data.table
#' @importFrom GenomicRanges reduce
#' @return A list of 3 in-memory objects with the elements:
#    SNS.vcf:   Data frame of pure SNS mutations -- no DNS or 3+BS mutations
#    DNS.vcf:   Data frame of pure DNS mutations -- no SNS or 3+BS mutations
#    ThreePlus: Data table with the key CHROM, LOW.POS, HIGH.POS and additional
#    information (reference sequence, alternative sequence, context, etc.)
#    Additional information not fully implemented at this point because of
#    limited immediate biological interest.
#' @export
SplitSNSVCF <- function(vcf.df, max.vaf.diff = 0.02) {
  stopifnot(class(vcf.df) == "data.frame")

  # Record the total number of input variants for later sanity checking.
  num.in <- nrow(vcf.df)

  # First we look for pairs of rows where the POS of one of the
  # rows is at POS + 1 of the other row. For example
  #
  # 1   200   foo    A    G
  # 1   201   foo    C    G
  #
  # Reprsents AC > GG
  #
  # But there could also be situations like this
  #
  # X   300  foo   C  T
  # X   301  foo   C  T
  # X   302  foo   A  G
  #
  # which represents CCA > TTG
  #
  # But first we just find pairs, so the
  # X chromosome example will appear as 2
  # pairs CC > and CA > TG (more below).

  vcf.dt <- data.table(vcf.df)
  vcf.dt[, POS.plus.one := POS + 1]
  dt2 <- merge(vcf.dt, vcf.dt,
               by.x = c("CHROM", "POS"),
               by.y = c("CHROM", "POS.plus.one"))

  # After this merge, each row contains one *pair*.
  # In each row, POS.y == POS - 1, and the neighboring SNS
  # are at postions POS and POS.y.
  dt2[, HIGH := POS]
  dt2[, LOW := POS.y]

  # Keep only SNS pairs that have very similar VAFs (variant allele frequencies).
  # If VAFs are not similar, the adjacent SNSs were likely to be "merely"
  # asynchronous single base mutations, and a simultaneous doublet mutation.
  non.SNS <- dt2[abs(VAF.x - VAF.y) <= max.vaf.diff]
  # If VAF.x or VAF.y is NA the row will not go into non.SNS.
  rm(dt2)

  if (nrow(non.SNS) == 0) {
    # Thre are no non.SNS mutations in the input.
    # Everything in vcf.df is an SNS. We are finished.
    empty <- vcf.df[-(1 : nrow(vcf.df)), ]
    return(list(SNS.vcf = vcf.df, DNS.vcf = empty,
                ThreePlus =
                  data.table(CHROM = character(),
                             LOW.POS = numeric(),
                             HIGH.POS = numeric())))
  }

  # Remove non SNS rows from the output VCF for the SNSs
  pairs.to.remove <-
    data.frame(non.SNS[, .(CHROM, POS = HIGH)])
  pairs.to.remove <-
    rbind(pairs.to.remove,
          data.frame(non.SNS[, .(CHROM, POS = LOW)]))
  dt.rm <- data.table(pairs.to.remove)
  dt.rm$delete.flag = TRUE
  out.SNS.dt <- merge(vcf.dt, dt.rm, by = c("CHROM", "POS"), all.x = TRUE)
  out.SNS.dt2 <- out.SNS.dt[is.na(delete.flag)]
  out.SNS.df <- as.data.frame(out.SNS.dt2[, delete.flag := NULL])
  num.SNS.out <- nrow(out.SNS.df)

  # Now separate doublets (DNS) from triplet and above base substitutions.
  # For ease of testing, keep only the genomic range information.
  non.SNS <- non.SNS[, c("CHROM", "LOW", "HIGH")]
  ranges <-
    as(data.frame(chrom = non.SNS$CHROM, start = non.SNS$LOW, end = non.SNS$HIGH),
       "GRanges")
  rranges <- reduce(ranges) # Merge overlapping ranges
  DNS.plus <- as.data.frame(rranges)
  if ((sum(DNS.plus$width) + num.SNS.out) != num.in) {
    if ((sum(DNS.plus$width) + num.SNS.out) > num.in) {
      cat("too many SNS\n")
      stop()
    } else {
      cat("possible site with multiple variant alleles involved in a DNS\n")
    }
  }
  DNSx <- DNS.plus[DNS.plus$width == 2, c("seqnames", "start", "end"), ]
  colnames(DNSx) <- c("CHROM", "LOW", "HIGH")
  DNS.vcf.df <- MakeVCFDNSdf(DNSx, vcf.dt)
  num.DNS.out <- nrow(DNS.vcf.df)

  other.ranges <- DNS.plus[DNS.plus$width > 2, ]
  num.other <- sum(other.ranges$width)

  if ((num.SNS.out + 2 * num.DNS.out + num.other) != num.in) {
    cat("Counts are off:", num.SNS.out, 2*num.DNS.out, num.other, "vs", num.in, "\n")
  }

  return(list(SNS.vcf = out.SNS.df, DNS.vcf = DNS.vcf.df,
              ThreePlus = other.ranges))
}

#' Check that the sequence context information is consistent with the value of
#' the column REF.
#'
#' @param vcf In-memory VCF as a data.frame; must be an SNS or DNS VCF.
#' @param column.to.use The column name as a string of the column in the VCF
#'   with the context information
#' @return Throws error with location information if the value of REF is
#'   inconsistent with the value of seq.21context. Assumes the first base of the
#'   reference allele is at position (size(<context string>)-1)/2, and generates
#'   error if this is not an integer. Indices are 1-based.
#' @keywords internal
CheckSeqContextInVCF <- function(vcf, column.to.use) {
  if (0 == nrow(vcf)) return()

  # Die if this is an indel VCF
  stopifnot(nchar(vcf$REF) == nchar(vcf$ALT))
  stopifnot(!any(vcf$REF == '-'))
  stopifnot(!any(vcf$ALT == '-'))
  cut.pos <- 1 + (nchar(vcf[, column.to.use]) - 1) / 2
  stopifnot(cut.pos == round(cut.pos))
  cut.from.ref <- substr(vcf[, column.to.use], cut.pos,
                         (cut.pos + nchar(vcf$REF)) - 1)
  error.rows <- which(vcf$REF != cut.from.ref)
  if (any(error.rows > 0)) {
    data.table::fwrite(x = as.data.table(vcf[error.rows, ]),
                       file = "error.rows.csv")
    cat("Seqence context of reference allele is inconsistent,
        see file error.rows.csv")
    stop()
  }
}

#' Read a list of Strelka VCF files from path
#'
#' @param vector.of.file.paths A vector containing the paths of the VCF files.
#'
#' @return A list of vcfs from vector.of.file.paths.
#' @export
ReadListOfStrelkaVCFs <- function(vector.of.file.paths) {
  vcfs <- lapply(vector.of.file.paths, FUN = ReadStrelkaVCF)
  names(vcfs) <- vector.of.file.paths
  return(vcfs)
}

#' Read a list of Mutect VCF files from path
#'
#' @param vector.of.file.paths A vector containing the paths of the VCF files.
#'
#' @return A list of vcfs from vector.of.file.paths.
#' @export
ReadListOfMutectVCFs <- function(vector.of.file.paths) {
  vcfs <- lapply(vector.of.file.paths, FUN = ReadMutectVCF)
  names(vcfs) <- vector.of.file.paths
  return(vcfs)
}

#' Create single nucleotide mutation catalog for *one* sample from
#' a Variant Call Format (VCF) file.
#'
#' @param vcf An in-memory VCF file annotated by the AddSequence and
#'   AddTranscript functions. It must *not* contain indels and must *not*
#'   contain DNS (double nucleotide substituions), or triplet base substitutions
#'   etc., even if encoded as neighboring SNS.
#' @param sample.id Usually the sample id, but defaults to "count".
#' @import data.table
#' @return A list of three matrices containing the SNS mutation catalog:
#'   96, 192, 1536 catalog respectively.
#' @keywords internal
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

#' Create SNS catalogs from VCFs
#'
#' Create a list of 3 catalogs (one each for 96, 192, 1536)
#' out of the contents of the VCFs in list.of.vcfs
#'
#' @param list.of.vcfs List of in-memory VCFs. The list names will be the
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
VCFsToSNSCatalogs <- function(list.of.vcfs, genome, trans.ranges) {
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
#'   substitutions etc.
#' @param sample.id Usually the sample id, but defaults to "count".
#' @import data.table
#' @return A list of three matrices containing the DNS catalog:
#'   catDNS78, catDNS144, catQUAD136 respectively.
#' @keywords internal
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

#' Create DNS catalogs from VCFs
#'
#' Create a list of 3 catalogs (one each for DNS78, DNS144 and QUAD136)
#' out of the contents of the VCFs in list.of.vcfs
#'
#' @param list.of.vcfs List of in-memory VCFs. The list names will be
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
VCFsToDNSCatalogs <- function(list.of.vcfs, genome, trans.ranges) {
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

#' Create SNS and DNS catalogs from VCF files
#'
#' Create 3 SNS catalogs (96, 192, 1536) and 3 DNS catalogs (78, 136, 144)
#' from the VCFs specified by vector.of.file.paths
#'
#' This function calls \code{\link{VCFsToSNSCatalogs}} and \code{\link{VCFsToDNSCatalogs}}
#' @param vector.of.file.paths A vector containing the paths of the VCF files.
#' @param genome  Name of a particular reference genome
#'   (without quotations marks).
#' @param trans.ranges A data.table which contains transcript range and
#'   strand information.
#' @return  A list of 3 SNS catalogs (one each for 96, 192, and 1536)
#'   and 3 DNS catalogs (one each for 78, 136, and 144)
#' @export
VCFFilesToCatalog <- function(vector.of.file.paths, genome, trans.ranges) {
  vcfs <- ReadListOfStrelkaVCFs(vector.of.file.paths)

  # TODO(nanhai): split.vcfs <- SplitListOfVCFs(vcfs)
  #
  # SplitListOfVCFs <- function(list.of.vcfs) {
  #   v1 <- lapply(list.of.vcfs, SplitSNSVCF.........
  #   SNS <- lapply(v1, funtion(x) x$SNS....)
  #   DNS <- la....
  #   return(list(SNS = SNS, DNS = DNS))
  # }
  #
  # return(c(VCFsToSNSCatalogs(split.vcfs$SNS, genome, trans.ranges),
  #     VCFsToDNSCatalogs(split.vcfs$DNS, genome, trans.ranges)))
  return(c(VCFsToSNSCatalogs(vcfs, genome, trans.ranges),
           VCFsToDNSCatalogs(vcfs, genome, trans.ranges)))
}

#' CanonicalizeDNS
#'
#' @param ref.vec TODO
#' @param alt.vec TODO
#'
#' @return TODO
#' @keywords internal
CanonicalizeDNS <- function(ref.vec, alt.vec) {
  # TODO document

  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  Canonicalize1DNS <- function(DNS) {
    if (DNS %in% .catalog.row.order.DNS.78) {
      return(DNS)
    } else {
      ref <- substr(DNS, 1, 2)
      alt <- substr(DNS, 3, 4)
      out <- paste0(revc(ref), revc(alt))
    }
    stopifnot(out %in% .catalog.row.order.DNS.78)
    return(out)
  }
  ret <- sapply(paste0(ref.vec, alt.vec), FUN = Canonicalize1DNS)
  return(ret)
}

#' CanonicalizeQUAD
#'
#' @param quad TODO
#'
#' @return TODO
#' @keywords internal
CanonicalizeQUAD <- function(quad) {
  # TODO document

  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  Canonicalize1QUAD <- function(quad) {
    if (quad %in% .catalog.row.order.QUAD.136) {
      return(quad)
    } else {
      out <- revc(quad)
      stopifnot(out %in% .catalog.row.order.QUAD.136)
      return(out)
    }
  }

  ret <- sapply(quad, FUN = Canonicalize1QUAD)
  return(ret)
}
