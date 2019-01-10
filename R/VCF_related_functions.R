#' Read in the data lines of a Variant Call Format (VCF) file
#' @importFrom utils read.csv
#' @param path The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#' @export
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
  return(df1)
}

#' Extract the VAFs (variant allele frequencies) from a VAF created by
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

#' Add sequence context to a data frame with mutation records
#'
#' @param df An input data frame storing mutation records of a VCF file.
#' @param seq A particular reference genome.
#'
#' @return A data frame with a new column added to the input data frame,
#'     which contains sequence context information.
#' @importFrom methods as
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @export
AddSequence <- function(df, seq = BSgenome.Hsapiens.1000genomes.hs37d5) {

  if (0 == nrow(df)) return(df)

  # Create a GRanges object with range width equals to 21
  Ranges <-
    as(data.frame(chrom = df$CHROM, start = df$POS - 10, end = df$POS + 10),
       "GRanges")

  # Extract sequence context from the reference genome
  df <- dplyr::mutate(df,
                      seq.21context = BSgenome::getSeq(seq,
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

#' Take DNS ranges and the original VCF and generate a VCF with
#' dinucleotide REF and ALT alleles. The output VCF has minimal columns:
#' just CHROM, POS, ID, REF, ALT.
#'
#' @param DNS.range.df Data frame with columns CHROM, LOW, HIGH
#' @param SNS.vcf.dt TODO
#' @import data.table
#' @return TODO
#' @export
# TODO(steve) add average VAF
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
  rranges <- GenomicRanges::reduce(ranges) # Merge overlapping ranges
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
#'
#' @return Throws error with location information if the value of REF is
#'   inconsistent with the value of seq.21context. Assumes the first base of the
#'   reference allele is at position (size(<context string>)-1)/2, and generates
#'   error if this is not an integer. Indices are 1-based.
#' @export
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

#' Read a list of VCF files from path
#'
#' @param vector.of.file.paths A vector containing the paths of the VCF files.
#'
#' @return A list of vcfs from vector.of.file.paths.
#' @export
ReadListOfVCFs <- function(vector.of.file.paths) {
  vcfs <- lapply(vector.of.file.paths, FUN = ReadStrelkaVCF)
  names(vcfs) <- vector.of.file.paths
  return(vcfs)
}

