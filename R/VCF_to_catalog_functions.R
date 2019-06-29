#' Extract the VAFs (variant allele frequencies) from a VCF file.
#'
#' @param vcf Said VCF as a data.frame.
#'
#' @return A vector of VAFs, one for each row of \code{vcf}.
#'
#' @name GetVAF
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' MakeDataFrameFromStrelkaSBSVCF <- 
#'   getFromNamespace("MakeDataFrameFromStrelkaSBSVCF", "ICAMS")
#' df <- MakeDataFrameFromStrelkaSBSVCF(file)
#' vaf <- GetStrelkaVAF(df)
NULL

#' Read in the data lines of an SBS VCF created by Strelka version 1
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#'
#' @keywords internal
MakeDataFrameFromStrelkaSBSVCF <- function(file) {
  df <- read.csv(file, header = FALSE, sep = "\t", quote = "",
                 col.names = paste0("c", 1 : 100), as.is = TRUE)
  
  # Delete the columns which are totally empty
  df <- df[!sapply(df, function(x) all(is.na(x)))]
  
  # Delete meta-information lines which start with "##"
  if (any(grepl("^##", df[, 1]))) {
    idx <- grep("^##", df[, 1])
    df1 <- df[-idx, ]
  } else {
    df1 <- df
  }
  
  # Extract the names of columns in the VCF file
  names <- c("CHROM", as.character(df1[1, ])[-1])
  df1 <- df1[-1, ]
  colnames(df1) <- names
  
  stopifnot(df1$REF != df1$ALT)
  df1$POS <- as.integer(df1$POS)
  
  # Is there any column in df1 with name "strand"?
  # If there is, change its name to "strand_old" so that it will
  # conflict with code in other parts of ICAMS package.
  if ("strand" %in% colnames(df1)) {
    colnames(df1)[which(colnames(df1) == "strand")] <- "strand_old"
    warning('There is column in VCF which has name "strand", ',
            'it has been renamed to "strand_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  
  # Is there any column in df1 with name "VAF"?
  # If there is, change its name to "VAF_old" so that it will
  # conflict with code in other parts of ICAMS package.
  if ("VAF" %in% colnames(df1)) {
    colnames(df1)[which(colnames(df1) == "VAF")] <- "VAF_old"
    warning('There is column in VCF which has name "VAF", ',
            'it has been renamed to "VAF_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  
  return(df1)
}

#' Read in the data lines of an SBS VCF created by Strelka version 1
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file with VAFs added.
#'
#' @keywords internal
ReadStrelkaSBSVCF <- function(file) {
  df <- MakeDataFrameFromStrelkaSBSVCF(file)
  df$VAF <- GetStrelkaVAF(df)
  return(StandardChromName(df))
}

#' Read in the data lines of an ID VCF created by Strelka version 1
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @keywords internal
ReadStrelkaIDVCF <- function(file) {
  df <- read.csv(file, header = FALSE, sep = "\t", quote = "",
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
  
  # Check whether the input VCF is a Strelka ID VCF
  if (!("TUMOR" %in% names(df1)) ||
      !("FORMAT" %in% names(df1))) {
    stop("\nvcf does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(df1), collapse=" "))
  }
  control <- unique(df1[ , "FORMAT"])
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (all(each.base.col %in% colnames)) {
    stop("\nvcf does not appear to be a Strelka ID VCF,", 
         "the value of column FORMAT is \n", 
         control)
  }

  df1$POS <- as.integer(df1$POS)
  return(StandardChromName(df1))
}

#' @rdname GetVAF
#'
#' @export
GetStrelkaVAF <-function(vcf) {
  stopifnot(class(vcf) == "data.frame")
  if (!("TUMOR" %in% names(vcf)) ||
      !("FORMAT" %in% names(vcf))) {
    stop("\nvcf does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(vcf), collapse=" "))
  }
  TUMOR <- vcf[ , "TUMOR"]
  control <- unique(vcf[ , "FORMAT"])
  alt     <- vcf[ , "ALT"]
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  values <- strsplit(TUMOR, split=":", fixed=TRUE)
  vaf <- numeric(nrow(vcf))
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (!all(each.base.col %in% colnames)) {
    stop("\nvcf does not appear to be a Strelka SBS VCF,", 
         "the value of column FORMAT is \n", 
         control)
  }
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

#' Read in the data lines of a Variant Call Format (VCF) file created by MuTect
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#'
#' @keywords internal
MakeDataFrameFromMutectVCF <- function(file) {
  df <- read.csv(file, header = FALSE, sep = "\t", quote = "",
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
  
  # Is there any column in df1 with name "strand"?
  # If there is, change its name to "strand_old" so that it will
  # conflict with code in other parts of ICAMS package.
  if ("strand" %in% colnames(df1)) {
    colnames(df1)[which(colnames(df1) == "strand")] <- "strand_old"
    warning('There is column in VCF which has name "strand", ',
            'it has been renamed to "strand_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  
  # Is there any column in df1 with name "VAF"?
  # If there is, change its name to "VAF_old" so that it will
  # conflict with code in other parts of ICAMS package.
  if ("VAF" %in% colnames(df1)) {
    colnames(df1)[which(colnames(df1) == "VAF")] <- "VAF_old"
    warning('There is column in VCF which has name "VAF", ',
            'it has been renamed to "VAF_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  
  return(StandardChromName(df1))
}

#' Read in the data lines of a Variant Call Format (VCF) file created by
#'     MuTect
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file with VAFs added.
#'
#' @keywords internal
ReadMutectVCF <- function(file) {
  df <- MakeDataFrameFromMutectVCF(file)
  df$VAF <- GetMutectVAF(df)
  return(StandardChromName(df))
}

#' @rdname GetVAF
#'
#' @export
GetMutectVAF <-function(vcf) {
  stopifnot(class(vcf) == "data.frame")
  if (!any(grepl("/1", unlist(vcf[1, ])))) {
    stop("vcf does not appear to be a Mutect VCF, please check the data")
  }

  # Select out the column which has the information for F1R2 and F2R1
  info <- vcf[[10]]

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

#' @title Split a mutect2 VCF into SBS, DBS, and ID VCFs, plus a list of other mutations
#'
#' @param vcf.df An in-memory data.frame representing a Mutect VCF, including
#'  VAFs, which are added by \code{\link{ReadMutectVCF}}.
#'
#' @return A list with the SBS, DBS, and ID portions of the VCF file, plus two
#' data.frames of other mutations
#'
#' @keywords internal
SplitOneMutectVCF <- function(vcf.df) {
  # Mutect VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)

  multiple.alt.df <- vcf.df[multiple.alt, ]
  df <- vcf.df[-multiple.alt, ]
  rm(multiple.alt, vcf.df)

  SBS.df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]

  DBS.df <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]

  other.df <- df[nchar(df$REF) > 2 & nchar(df$ALT) == nchar(df$REF), ]

  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]

  return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
              other=other.df, multiple.alt = multiple.alt.df))

}

#' Split each Mutect VCF into SBS, DBS, and ID VCFs (plus two
#' VCF-like data frame with left-over rows).
#'
#' @param list.of.vcfs List of VCFs as in-memory data.frames.
#'
#' @return A list with 3 in-memory VCFs and two left-over
#' VCF-like data frames with rows that were not incorporated
#' into the first 3 VCFs, as follows:
#'
#' \enumerate{
#'
#'  \item \code{SBS} VCF with only single base substitutions.
#'
#'  \item \code{DBS} VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  \item \code{ID} VCF with only small insertions and deletions.
#'
#'  \item \code{other.subs} VCF like data.frame with
#'  rows for coordinate substitutions involving
#'  3 or more nucleotides, e.g. ACT > TGA or AACT > GGTA.
#'
#'  \item \code{multiple.alternative.alleles} VCF like data.frame with
#'  rows for variants with multiple alternative alleles, for example
#'  ACT mutated to both AGT and ACT at the same position.
#'
#' }
#'
#' @keywords internal
SplitListOfMutectVCFs <- function(list.of.vcfs) {
  v1 <- lapply(list.of.vcfs, SplitOneMutectVCF)
  SBS <- lapply(v1, function(x) x$SBS)
  DBS <- lapply(v1, function(x) x$DBS)
  ID  <- lapply(v1, function(x) x$ID)
  other.subs <- lapply(v1, function(x) x$other.df)
  multiple.alternative.alleles <-
    lapply(v1, function(x) x$multiple.alt)

  return(list(SBS = SBS, DBS = DBS, ID = ID,
              other.subs = other.subs,
              multiple.alternative.alleles
              = multiple.alternative.alleles
  ))
}

#' Add sequence context to a data frame with mutation records
#'
#' @param df An input data frame storing mutation records of a VCF file.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param seq.context.width The number of preceding and following bases to be
#'   extracted around the mutated position from \code{ref.genome}. The default value is 10.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @importFrom BSgenome getSeq seqnames
#'
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#'
#' @import BSgenome.Hsapiens.UCSC.hg38
#' 
#' @import BSgenome.Mmusculus.UCSC.mm10
#'
#' @importFrom stats start end
#'
#' @return A copy of the input data.frame with a new column added
#'     that contains sequence context information.
#'
#' @keywords internal
AddSeqContext <- function(df, ref.genome, seq.context.width = 10) {
  if (0 == nrow(df)) return(df)
  ref.genome <- NormalizeGenomeArg(ref.genome)

  # Check if the format of sequence names in df and genome are the same.
  # Internally ICAMS uses human chromosomes labeled as "1", "2", ... "X"...
  # However, BSgenome.Hsapiens.UCSC.hg38 has chromosomes labeled
  # "chr1", "chr2", ....
  vcf.chr.names <- unique(df$CHROM)
  if (!all(vcf.chr.names %in% seqnames(ref.genome))) {
    tmp.chr <- paste0("chr", vcf.chr.names)
    if (!all(tmp.chr %in% seqnames(ref.genome))) {
      stop("Cannot match chromosome names:\n",
           sort(vcf.chr.names), "\nversus\n", sort(seqnames(ref.genome)))
    }

    chr.names <- paste0("chr", df$CHROM)
  } else {
    chr.names <- df$CHROM
  }
  # Create a GRanges object with the needed width.
  Ranges <-
    GRanges(chr.names,
            IRanges(start = df$POS - seq.context.width, # 10,
                    end = df$POS + seq.context.width) # 10
    )

  # Extract sequence context from the reference genome
  df$extracted.seq <- getSeq(ref.genome, Ranges, as.character = TRUE)

  names(df)[names(df) == "extracted.seq"] <-
    paste0("seq.", 2 * seq.context.width + 1, "bases")
  return(df)
}

#' Add transcript information to a data frame with mutation records
#'
#' @param df A data frame storing mutation records of a VCF file.
#'
#' @param trans.ranges A \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @import data.table
#'
#' @return A data frame with new columns added to the input data frame,
#'     which contain the mutated gene's name, range and strand information.
#'
#' @keywords internal
AddTranscript <- function(df, trans.ranges = NULL) {
  if (nrow(df) == 0) {
    return(df)
  }
  
  if (is.null(trans.ranges)) {
    return(data.table(df))
  }

  # Find range overlaps between the df and trans.ranges
  df1 <- data.table(df)
  df1[, POS2 := POS]
  dt <- foverlaps(df1, trans.ranges,
                  by.x = c("CHROM", "POS", "POS2"),
                  type = "within", mult = "all")

  # Find out mutations that fall on transcripts on both strands
  dt1 <- dt[, bothstrand := "+" %in% strand && "-" %in% strand,
            by = .(CHROM, ALT, POS)] # Note that is important to have
  # ALT in the by list because in a few cases
  # there are multiple ALT alleles at one POS.

  # Count the number of transcript ranges where a particular mutation
  # falls into
  dt2 <- dt1[, count := .N, by = .(CHROM, ALT, POS)]

  # Swap gene location according to strand information
  dt3 <- dt2[strand == "-", c("end", "start") := .(start, end)]

  # Reorder the columns of dt3
  df.colnames <- colnames(df)
  trans.ranges.colnames <- colnames(trans.ranges)[-1]
  setcolorder(dt3, neworder = c(df.colnames, trans.ranges.colnames))

  # Delete redundant column in dt3
  dt4 <- dt3[, POS2 := NULL]

  return(dt4)
}

#' MakeVCFDBSdf Take DBS ranges and the original VCF and generate a VCF with
#' dinucleotide REF and ALT alleles.
#'
#' @param DBS.range.df Data frame with columns CHROM, LOW, HIGH
#'
#' @param SBS.vcf.dt A data table containing the VCF from which
#' \code{DBS.range.df} was computed.
#'
#' @import data.table
#'
#' @return A minimal VCF with only the columns \code{CHROM},
#' \code{POS}, \code{ID}, \code{REF}, \code{ALT}.
#'
#' @keywords internal
MakeVCFDBSdf <- function(DBS.range.df, SBS.vcf.dt) {
  tmpvcf <- SBS.vcf.dt[ , c("CHROM", "POS", "REF", "ALT")]
  DBS.range.dt <- as.data.table(DBS.range.df)
  tmp1 <- merge(DBS.range.dt, tmpvcf,
                by.x = c("CHROM", "LOW"),
                by.y = c("CHROM", "POS"))
  tmp2 <- merge(tmp1, tmpvcf,
                by.x = c("CHROM", "HIGH"),
                by.y = c("CHROM", "POS"))
  tmp2[, POS := LOW]
  tmp2[, ID := "From merged SBSs"]
  tmp2[, REF := paste0(REF.x, REF.y)]
  tmp2[, ALT := paste0(ALT.x, ALT.y)]
  return(as.data.frame(tmp2[, c("CHROM", "POS", "ID", "REF", "ALT")]))
}

#' Split an in-memory Strelka VCF into SBS, DBS, and variants involving
#' > 2 consecutive bases
#'
#' SBSs are single base substitutions,
#' e.g. C>T, A<G,....  DBSs are double base substitutions,
#' e.g. CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT > TCTA, ...
#'
#' @param vcf.df An in-memory data frame containing a Strelka VCF file contents.
#'
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02.
#'
#' @import data.table
#'
#' @importFrom stats start end
#'
#' @importFrom GenomicRanges GRanges reduce
#'
#' @importFrom IRanges IRanges
#'
#' @return A list of 3 in-memory objects with the elements:
#' \enumerate{
#'    \item \code{SBS.vcf}:   Data frame of pure SBS mutations -- no DBS or 3+BS mutations
#'    \item \code{DBS.vcf}:   Data frame of pure DBS mutations -- no SBS or 3+BS mutations
#'    \item{ThreePlus}: Data table with the key CHROM, LOW.POS, HIGH.POS and additional
#'    information (reference sequence, alternative sequence, context, etc.)
#'    Additional information not fully implemented at this point because of
#'    limited immediate biological interest.
#'    }
#'
#' @keywords internal
SplitStrelkaSBSVCF <- function(vcf.df, max.vaf.diff = 0.02) {
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
  # In each row, POS.y == POS - 1, and the neighboring SBS
  # are at postions POS and POS.y.
  dt2[, HIGH := POS]
  dt2[, LOW := POS.y]

  # Keep only SBS pairs that have very similar VAFs (variant allele frequencies).
  # If VAFs are not similar, the adjacent SBSs were likely to be "merely"
  # asynchronous single base mutations, and a simultaneous doublet mutation.
  non.SBS <- dt2[abs(VAF.x - VAF.y) <= max.vaf.diff]
  # If VAF.x or VAF.y is NA the row will not go into non.SBS.
  rm(dt2)

  if (nrow(non.SBS) == 0) {
    # Thre are no non.SBS mutations in the input.
    # Everything in vcf.df is an SBS. We are finished.
    empty <- vcf.df[-(1 : nrow(vcf.df)), ]
    return(list(SBS.vcf = vcf.df, DBS.vcf = empty,
                ThreePlus =
                  data.table(CHROM = character(),
                             LOW.POS = numeric(),
                             HIGH.POS = numeric())))
  }

  # Remove non SBS rows from the output VCF for the SBSs
  pairs.to.remove <-
    data.frame(non.SBS[, .(CHROM, POS = HIGH)])
  pairs.to.remove <-
    rbind(pairs.to.remove,
          data.frame(non.SBS[, .(CHROM, POS = LOW)]))
  dt.rm <- data.table(pairs.to.remove)
  dt.rm$delete.flag = TRUE
  out.SBS.dt <- merge(vcf.dt, dt.rm, by = c("CHROM", "POS"), all.x = TRUE)
  out.SBS.dt2 <- out.SBS.dt[is.na(delete.flag)]
  out.SBS.df <- as.data.frame(out.SBS.dt2[, delete.flag := NULL])
  num.SBS.out <- nrow(out.SBS.df)

  # Now separate doublets (DBSs) from triplet and above base substitutions.
  # For ease of testing, keep only the genomic range information.
  non.SBS <- non.SBS[, c("CHROM", "LOW", "HIGH")]
  ranges <-
    GRanges(non.SBS$CHROM, IRanges(start = non.SBS$LOW, end = non.SBS$HIGH))
  rranges <- reduce(ranges) # Merge overlapping ranges
  DBS.plus <- as.data.frame(rranges)
  if ((sum(DBS.plus$width) + num.SBS.out) != num.in) {
    if ((sum(DBS.plus$width) + num.SBS.out) > num.in) {
      stop("Possible programming error or input problem: too many SBS")
    } else {
      warning("Possible site with multiple variant alleles involved in a DBS\n")
    }
  }
  DBSx <- DBS.plus[DBS.plus$width == 2, c("seqnames", "start", "end"), ]
  colnames(DBSx) <- c("CHROM", "LOW", "HIGH")
  DBS.vcf.df <- MakeVCFDBSdf(DBSx, vcf.dt)
  num.DBS.out <- nrow(DBS.vcf.df)

  other.ranges <- DBS.plus[DBS.plus$width > 2, ]
  num.other <- sum(other.ranges$width)

  if ((num.SBS.out + 2 * num.DBS.out + num.other) != num.in) {
    warning("Counts are off:", num.SBS.out, 2*num.DBS.out, num.other, "vs", num.in, "\n")
  }

  return(list(SBS.vcf = out.SBS.df, DBS.vcf = DBS.vcf.df,
              ThreePlus = other.ranges))
}

#' Split a list of in-memory Strelka SBS VCF into SBS, DBS, and variants involving
#' > 2 consecutive bases
#'
#' SBSs are single base substitutions,
#' e.g. C>T, A<G,....  DBSs are double base substitutions,
#' e.g. CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT > TCTA, ...
#'
#' @param list.of.vcfs A list of in-memory data frames containing Strelka SBS VCF file contents.
#'
#' @return A list of 3 in-memory objects with the elements:
#'    SBS.vcfs:  List of Data frames of pure SBS mutations -- no DBS or 3+BS mutations
#'    DBS.vcfs:  List of Data frames of pure DBS mutations -- no SBS or 3+BS mutations
#'    ThreePlus: List of Data tables with the key CHROM, LOW.POS, HIGH.POS and additional
#'    information (reference sequence, alternative sequence, context, etc.)
#'    Additional information not fully implemented at this point because of
#'    limited immediate biological interest.
#'
#' @keywords internal
SplitListOfStrelkaSBSVCFs <- function(list.of.vcfs) {
  split.vcfs<- lapply(list.of.vcfs, FUN = SplitStrelkaSBSVCF)
  SBS.vcfs <- lapply(split.vcfs, function(x) x$SBS.vcf)
  DBS.vcfs <- lapply(split.vcfs, function(x) x$DBS.vcf)
  ThreePlus <- lapply(split.vcfs, function(x) x$ThreePlus)
  return(list(SBS.vcfs = SBS.vcfs, DBS.vcfs = DBS.vcfs, ThreePlus = ThreePlus))
}

#' Check that the sequence context information is consistent with the value of
#' the column REF.
#'
#' @param vcf In-memory VCF as a data.frame; must be an SBS or DBS VCF.
#'
#' @param column.to.use The column name as a string of the column in the VCF
#'   with the context information.
#'
#' @return Throws error with location information if the value of REF is
#'   inconsistent with the value of seq.21bases. Assumes the first base of the
#'   reference allele is at position (size(<context string>)-1)/2, and generates
#'   error if this is not an integer. Indices are 1-based.
#'
#' @keywords internal
CheckSeqContextInVCF <- function(vcf, column.to.use) {
  if (0 == nrow(vcf)) return()

  # Die if this is an indel VCF
  stopifnot(nchar(vcf$REF) == nchar(vcf$ALT))
  stopifnot(!any(vcf$REF == '-'))
  stopifnot(!any(vcf$ALT == '-'))
  cut.pos <- 1 + (nchar(vcf$column.to.use) - 1) / 2
  stopifnot(cut.pos == round(cut.pos))
  cut.from.ref <- substr(vcf$column.to.use, cut.pos,
                         (cut.pos + nchar(vcf$REF)) - 1)
  error.rows <- which(vcf$REF != cut.from.ref)
  if (any(error.rows > 0)) {
    temp <- tempfile(fileext = ".csv")
    data.table::fwrite(x = as.data.table(vcf[error.rows, ]),
                       file = temp)
    stop("Seqence context of reference allele is inconsistent,",
         "see file ", temp)
  }
}

#' Read Strelka SBS (single base substitutions) VCF files.
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @return A list of vcfs from \code{files}.
#'
#' @keywords internal
ReadStrelkaSBSVCFs <- function(files) {
  vcfs <- lapply(files, FUN = ReadStrelkaSBSVCF)
  names(vcfs) <- tools::file_path_sans_ext(basename(files))
  return(vcfs)
}

#' Read and split Strelka SBS VCF files.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @return A list of 3 in-memory objects as follows:
#' \enumerate{
#'    \item \code{SBS.vcfs} List of data.frames of pure SBS mutations -- no DBS or 3+BS mutations.
#'
#'    \item \code{DBS.vcfs} List of data.frames of pure DBS mutations -- no SBS or 3+BS mutations.
#'
#'    \item \code{ThreePlus} List of data.tables with the key CHROM, LOW.POS, HIGH.POS. containing
#'    rows that that in the input that did not represent SBSs or DBSs.
#'
#'    }
#'
#' @seealso \code{\link{StrelkaSBSVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
ReadAndSplitStrelkaSBSVCFs <- function(files) {
  vcfs <- ReadStrelkaSBSVCFs(files)
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  return(split.vcfs)
}

#' Read Strelka ID (insertion and deletion) VCF files.
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @return A list of vcfs from \code{files}.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @seealso \code{\link{StrelkaIDVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadStrelkaIDVCFs(file)
ReadStrelkaIDVCFs <- function(files) {
  vcfs <- lapply(files, FUN = ReadStrelkaIDVCF)
  names(vcfs) <- tools::file_path_sans_ext(basename(files))
  return(vcfs)
}

#' Read Mutect VCF files.
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @return A list of vcfs from \code{files}.
#'
#' @keywords internal
ReadMutectVCFs <- function(files) {
  vcfs <- lapply(files, FUN = ReadMutectVCF)
  names(vcfs) <- tools::file_path_sans_ext(basename(files))
  return(vcfs)
}

#' Read and split Mutect VCF files.
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @return A list with 3 in-memory VCFs and two left-over
#' VCF-like data frames with rows that were not incorporated
#' into the first 3 VCFs, as follows:
#'
#' \enumerate{
#'
#'  \item \code{SBS} VCF with only single base substitutions.
#'
#'  \item \code{DBS} VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  \item \code{ID} VCF with only small insertions and deletions.
#'
#'  \item \code{other.subs} VCF like data.frame with
#'  rows for coordinate substitutions involving
#'  3 or more nucleotides, e.g. ACT > TGA or AACT > GGTA.
#'
#'  \item \code{multiple.alternative.alleles} VCF like data.frame with
#'  rows for variants with multiple alternative alleles, for example
#'  ACT mutated to both AGT and ACT at the same position.
#'
#' }
#'
#' @seealso \code{\link{MutectVCFFilesToCatalog}}
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitMutectVCFs(file)
ReadAndSplitMutectVCFs <- function(files) {
  vcfs <- ReadMutectVCFs(files)
  split.vcfs <- SplitListOfMutectVCFs(vcfs)
  return(split.vcfs)
}

#' Create the matrix an SBS catalog for *one* sample from an in-memory VCF.
#'
#' @param vcf An in-memory VCF file annotated by the AddSeqContext and
#'   AddTranscript functions. It must *not* contain indels and must *not*
#'   contain DBS (double base substitutions), or triplet base substitutions
#'   etc., even if encoded as neighboring SBS.
#'   
#' @param trans.ranges A \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.  
#'
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @import data.table
#'
#' @return A list of three 1-column matrices with the names
#' \code{catSBS96}, \code{catSBS192}, \code{catSBS1536}.
#'  If trans.ranges is NULL, \code{catSBS192} is not generated.
#'  Do not rely on the order of elements in the list.
#'
#' @note catSBS192 only contains mutations in transcribed regions.
#'
#' @keywords internal
CreateOneColSBSMatrix <- function(vcf, trans.ranges = NULL, 
                                   sample.id = "count") {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude DBSs
  # before calling the function. This function does not detect DBSs.

  if (0 == nrow(vcf)) {
    # Create 1-column matrix with all values being 0 and the correct row labels.
    catSBS96 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS96), ncol = 1)
    rownames(catSBS96) <- ICAMS::catalog.row.order$SBS96
    catSBS192 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS192), ncol = 1)
    rownames(catSBS192) <- ICAMS::catalog.row.order$SBS192
    catSBS1536 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS1536), ncol = 1)
    rownames(catSBS1536) <- ICAMS::catalog.row.order$SBS1536

    return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                catSBS1536 = catSBS1536))
  }

  stopifnot(nchar(vcf$ALT) == 1)
  stopifnot(nchar(vcf$REF) == 1)
  stopifnot(vcf$ALT != vcf$REF)
  # Create 2 new columns that show the 3072 and 1536 mutation type
  context <- substr(vcf$seq.21bases, 9, 13)
  vcf$mutation <- paste0(context, vcf$ALT)

  # PyrPenta maps to strand-agnostic category
  # e.g. ATGCT>T "ATGCTT" maps to AGCAT>A, "AGCATA"
  vcf$pyr.mut <- PyrPenta(vcf$mutation)

  # One SBS mutation can be represented by more than 1 row in vcf after annotated by
  # AddTranscript function if the mutation position falls into the range of
  # multiple transcripts. When creating the 1536 and 96 catalog, we only need to
  # count these mutations once.
  vcf1 <- vcf[, .(REF = REF[1], pyr.mut = pyr.mut[1]),
              by = .(CHROM, ALT, POS)]

  # Create part of the 1536 catalog matrix but missing mutation
  # types have NA in the count column.
  tab1536 <- table(vcf1[, "pyr.mut"])
  stopifnot(setequal(
    setdiff(names(tab1536), ICAMS::catalog.row.order$SBS1536),
    c()))
  dt1536  <- data.table(tab1536)

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
  
  if (is.null(trans.ranges)) {
    return(list(catSBS96 = mat96, catSBS1536 = mat1536))
  }
  
  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 192 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One SBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 192 catalog, we only need to count these mutations once.
  vcf3 <- vcf2[, .(REF = REF[1], mutation = mutation[1], strand = strand[1]),
              by = .(CHROM, ALT, POS)]

  # Create the 192 catalog matrix
  tab192  <- table(paste0(substr(vcf3$mutation, 2, 4),
                          substr(vcf3$mutation, 6, 6)),
                   vcf3$strand,
                   useNA = "ifany")
  stopifnot(sum(tab192) == nrow(vcf3))
  dt192 <- as.data.table(tab192)
  colnames(dt192) <- c("rn", "strand", "count")
  dt192 <- dt192[!is.na(strand)]
  dt192[strand == "-", rn := RevcSBS96(rn)]
  dt192 <- dt192[ , .(count = sum(count)), by = rn]
  x192 <- data.table(rn = ICAMS::catalog.row.order$SBS192)
  x <- merge(x192, dt192, by = "rn", all.x = TRUE)
  x[is.na(count), count := 0]
  mat192 <- matrix(x[, count])
  rownames(mat192) <- unlist(x[, 1])
  mat192 <- mat192[ICAMS::catalog.row.order$SBS192, , drop = FALSE]
  colnames(mat192) <- sample.id

  return(list(catSBS96 = mat96, catSBS192 = mat192, catSBS1536 = mat1536))
}

#' Create SBS catalogs from SBS VCFs
#'
#' Create a list of 3 catalogs (one each for 96, 192, 1536)
#' out of the contents in list.of.SBS.vcfs. The SBS VCFs must not contain
#' DBSs, indels, or other types of mutations.
#'
#' @param list.of.SBS.vcfs List of in-memory data frames of pure SBS mutations
#'   -- no DBS or 3+BS mutations. The list names will be the sample ids in the
#'   output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return A list of 3 SBS catalogs, one each for 96, 192, 1536: catSBS96
#'   catSBS192 catSBS1536. If trans.ranges = NULL, SBS 192 catalog will not be
#'   generated. Each catalog has attributes added. See \code{\link{as.catalog}}
#'   for more details.
#'
#' @note SBS 192 catalogs only contain mutations in transcribed regions.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.SBS.vcfs <- ReadAndSplitMutectVCFs(file)$SBS
#' catalogs.SBS <- VCFsToSBSCatalogs(list.of.SBS.vcfs, ref.genome = "hg19",
#'                                   trans.ranges = trans.ranges.GRCh37,
#'                                   region = "genome")
VCFsToSBSCatalogs <- function(list.of.SBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown") {
  ncol <- length(list.of.SBS.vcfs)

  catSBS96 <- empty.cats$catSBS96
  catSBS192 <- empty.cats$catSBS192
  catSBS1536 <- empty.cats$catSBS1536

  for (i in 1:ncol) {
    SBS <- list.of.SBS.vcfs[[i]]

    SBS <- AddSeqContext(SBS, ref.genome = ref.genome)

    # Delete the rows of SBS if the extracted sequence contains "N"
    idx <- grep("N", substr(SBS$seq.21bases, 9, 13))
    if (!length(idx) == 0) {
      SBS <- SBS[-idx, ]
      message(
        'Rows in the SBS vcf where surrounding sequence contains "N" ',
        'have been deleted so as not to conflict with downstrea processing')
    }

    CheckSeqContextInVCF(SBS, "seq.21bases")
    SBS <- AddTranscript(SBS, trans.ranges)
    SBS.cat <- CreateOneColSBSMatrix(SBS, trans.ranges)
    rm(SBS)
    catSBS96 <- cbind(catSBS96, SBS.cat$catSBS96)
    if (!is.null(trans.ranges)) {
      catSBS192 <- cbind(catSBS192, SBS.cat$catSBS192)
    }
    catSBS1536 <- cbind(catSBS1536, SBS.cat$catSBS1536)
  }

  colnames(catSBS96) <- names(list.of.SBS.vcfs)
  colnames(catSBS1536) <- names(list.of.SBS.vcfs)

  catSBS96 <-
    as.catalog(catSBS96, ref.genome = ref.genome,
               region = region, catalog.type = "counts")

  catSBS1536 <-
    as.catalog(catSBS1536, ref.genome = ref.genome,
               region = region, catalog.type = "counts",
               abundance = NULL)
  if (is.null(trans.ranges)) {
    return(list(catSBS96 = catSBS96, catSBS1536 = catSBS1536))
  }
  
  colnames(catSBS192) <- names(list.of.SBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catSBS192 <-
    as.catalog(catSBS192, ref.genome = ref.genome,
               region = in.transcript.region, 
               catalog.type = "counts",
               abundance = NULL)
  return(list(catSBS96 = catSBS96, catSBS192 = catSBS192, 
              catSBS1536 = catSBS1536))
}

#' Create double base catalog for *one* sample from
#' a Variant Call Format (VCF) file
#'
#' @param vcf An in-memory VCF file annotated by the AddSeqContext and
#'   AddTranscript functions. It must *not* contain indels and must
#'   *not* contain SBS (single base substitutions), or triplet base
#'   substitutions etc.
#' 
#' @param trans.ranges A \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.  
#'   
#' @param sample.id Usually the sample id, but defaults to "count".
#'
#' @import data.table
#'
#' @return A list of three 1-column matrices with the names
#' \code{catDBS78}, \code{catDBS144}, and \code{catDBS136}.
#'  If trans.ranges is NULL, \code{catDBS144} is
#'   not generated. Do not rely on the order of elements
#'   in the list.
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @keywords internal
CreateOneColDBSMatrix <- function(vcf, trans.ranges = NULL,
                                   sample.id = "count") {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude SBSs
  # before calling the function. This function does not detect SBSs.

  if (0 == nrow(vcf)) {
    # Create 1-column matrix with all values being 0 and the correct row labels.
    catDBS78 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS78), ncol = 1)
    rownames(catDBS78) <- ICAMS::catalog.row.order$DBS78
    catDBS144 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS144), ncol = 1)
    rownames(catDBS144) <- ICAMS::catalog.row.order$DBS144
    catDBS136 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS136), ncol = 1)
    rownames(catDBS136) <- ICAMS::catalog.row.order$DBS136

    return(list(catDBS78 = catDBS78,
                catDBS136 = catDBS136,
                catDBS144 = catDBS144))
  }

  stopifnot(nchar(vcf$ALT) == 2)
  stopifnot(nchar(vcf$REF) == 2)

  # One DBS mutation can be represented by more than 1 row in vcf after annotated by
  # AddTranscript function if the mutation position falls into the range of
  # multiple transcripts. When creating the 78 and 136 catalog, we only need to
  # count these mutations once.
  vcf1 <- vcf[, .(REF = REF[1], seq.21bases = seq.21bases[1]),
              by = .(CHROM, ALT, POS)]

  # Create the 78 DBS catalog matrix
  canon.DBS.78 <- CanonicalizeDBS(vcf1$REF, vcf1$ALT)
  tab.DBS.78 <- table(canon.DBS.78)
  row.order.78 <- data.table(rn = ICAMS::catalog.row.order$DBS78)
  DBS.dt.78 <- as.data.table(tab.DBS.78)

  # DBS.dt.78 has two columns, names canon.DBS.78 (from the table() function)
  # and N (the count)
  DBS.dt.78.2 <-
    merge(row.order.78, DBS.dt.78,
          by.x = "rn", by.y = "canon.DBS.78", all = TRUE)
  DBS.dt.78.2[is.na(N), N := 0]
  stopifnot(DBS.dt.78.2$rn == ICAMS::catalog.row.order$DBS78)
  DBS.mat.78 <- as.matrix(DBS.dt.78.2[, 2])
  rownames(DBS.mat.78) <- DBS.dt.78.2$rn
  colnames(DBS.mat.78)<- sample.id

  # Create the 136 DBS catalog matrix
  canon.DBS.136 <- CanonicalizeQUAD(substr(vcf1$seq.21bases, 10, 13))
  tab.DBS.136 <- table(canon.DBS.136)
  row.order.136 <- data.table(rn = ICAMS::catalog.row.order$DBS136)
  DBS.dt.136 <- as.data.table(tab.DBS.136)

  # DBS.dt.136 has two columns, names canon.DBS.136 (from the table() function)
  # and N (the count)
  DBS.dt.136.2 <-
    merge(row.order.136, DBS.dt.136,
          by.x = "rn", by.y = "canon.DBS.136", all = TRUE)
  DBS.dt.136.2[is.na(N), N := 0]
  stopifnot(DBS.dt.136.2$rn == ICAMS::catalog.row.order$DBS136)
  DBS.mat.136 <- as.matrix(DBS.dt.136.2[, 2])
  rownames(DBS.mat.136) <- DBS.dt.136.2$rn
  colnames(DBS.mat.136)<- sample.id
  
  if (is.null(trans.ranges)) {
    return(list(catDBS78 = DBS.mat.78, catDBS136 = DBS.mat.136))
  }
  
  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 144 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One DBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 144 catalog, we only need to count these mutations once.
  vcf3 <- vcf2[, .(REF = REF[1], strand = strand[1]),
               by = .(CHROM, ALT, POS)]

  # Create the 144 DBS catalog matrix
  # There are 144 stranded DBSs: 4 X 4 sources and 3 X 3 alternates;
  # 4 x 4 x 3 x 3 = 144.
  tab.DBS.144  <-
    table(paste0(vcf3$REF, vcf3$ALT), vcf3$strand, useNA = "ifany")
  stopifnot(sum(tab.DBS.144) == nrow(vcf3))
  DBS.dt.144 <- as.data.table(tab.DBS.144)
  colnames(DBS.dt.144) <- c("rn", "strand", "count")
  DBS.dt.144 <- DBS.dt.144[!is.na(strand)]
  DBS.dt.144[strand == "-", rn := RevcDBS144(rn)]
  DBS.dt.144 <- DBS.dt.144[, .(count = sum(count)), by = rn]
  row.order.144 <- data.table(rn = ICAMS::catalog.row.order$DBS144)

  # DBS.dt.144 has two columns, names rn and count
  DBS.dt.144.2 <- merge(row.order.144, DBS.dt.144, by = "rn", all.x = TRUE)
  DBS.dt.144.2[is.na(count), count := 0]
  stopifnot(DBS.dt.144.2$rn == ICAMS::catalog.row.order$DBS144)
  DBS.mat.144 <- as.matrix(DBS.dt.144.2[, 2])
  rownames(DBS.mat.144) <- DBS.dt.144.2$rn
  colnames(DBS.mat.144)<- sample.id

  return(list(catDBS78 = DBS.mat.78, catDBS144 = DBS.mat.144,
              catDBS136 = DBS.mat.136))
}

#' Create DBS catalogs from VCFs
#'
#' Create a list of 3 catalogs (one each for DBS78, DBS144 and DBS136)
#' out of the contents in list.of.DBS.vcfs. The VCFs must not contain
#' any type of mutation other then DBSs.
#'
#' @param list.of.DBS.vcfs List of in-memory data frames of pure DBS mutations
#'   -- no SBS or 3+BS mutations. The list names will be the sample ids in the
#'   output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return A list of 3 DBS catalogs, one each for 78, 144, 136: catDBS78
#'   catDBS144 catDBS136. If trans.ranges = NULL, DBS 144 catalog will not be
#'   generated. Each catalog has attributes added. See \code{\link{as.catalog}}
#'   for more details.
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.DBS.vcfs <- ReadAndSplitMutectVCFs(file)$DBS
#' catalogs.DBS <- VCFsToDBSCatalogs(list.of.DBS.vcfs, ref.genome = "hg19",
#'                                   trans.ranges = trans.ranges.GRCh37,
#'                                   region = "genome")
VCFsToDBSCatalogs <- function(list.of.DBS.vcfs, ref.genome, 
                              trans.ranges = NULL, region = "unknown") {
  ncol <- length(list.of.DBS.vcfs)

  catDBS78 <- empty.cats$catDBS78
  catDBS136 <- empty.cats$catDBS136
  catDBS144 <- empty.cats$catDBS144


  for (i in 1 : ncol) {
    DBS <- list.of.DBS.vcfs[[i]]

    DBS <- AddSeqContext(DBS, ref.genome = ref.genome)

    # Delete the rows of DBS if the extracted sequence contains "N"
    idx <- grep("N", substr(DBS$seq.21bases, 10, 13))
    if (!length(idx) == 0) {
      DBS <- DBS[-idx, ]
      message(
        'Rows in the DBS vcf where surrounding sequence contains "N" ',
        'have been deleted so as not to conflict with downstrea processing')
    }

    DBS <- AddTranscript(DBS, trans.ranges)
    CheckSeqContextInVCF(DBS, "seq.21bases")
    DBS.cat <- CreateOneColDBSMatrix(DBS, trans.ranges)
    rm(DBS)
    catDBS78 <- cbind(catDBS78, DBS.cat$catDBS78)
    catDBS136 <- cbind(catDBS136, DBS.cat$catDBS136)
    if (!is.null(trans.ranges)) {
      catDBS144 <- cbind(catDBS144, DBS.cat$catDBS144)
    }
  }

  colnames(catDBS78) <- names(list.of.DBS.vcfs)
  colnames(catDBS136) <- names(list.of.DBS.vcfs)

  catDBS78 <-
    as.catalog(catDBS78, ref.genome = ref.genome,
               region = region, catalog.type = "counts")
  
  catDBS136 <-
    as.catalog(catDBS136, 
               ref.genome = ref.genome,
               region = region,
               catalog.type = "counts",
               abundance = NULL)

  if (is.null(trans.ranges)) {
    return(list(catDBS78 = catDBS78, catDBS136 = catDBS136))
  }
  colnames(catDBS144) <- names(list.of.DBS.vcfs)
  in.transcript.region <- ifelse(region == "genome", "transcript", region)
  catDBS144 <-
    as.catalog(catDBS144, 
               ref.genome = ref.genome,
               region = in.transcript.region, 
               catalog.type = "counts",
               abundance = NULL)
  return(list(catDBS78 = catDBS78, catDBS136 = catDBS136, 
              catDBS144 = catDBS144))
}

#' Create SBS and DBS catalogs from Strelka SBS VCF files.
#'
#' Create 3 SBS catalogs (96, 192, 1536) and 3 DBS catalogs (78, 136, 144)
#' from the Strelka SBS VCFs specified by \code{files}
#'
#' This function calls \code{\link{VCFsToSBSCatalogs}} and
#' \code{\link{VCFsToDBSCatalogs}}.
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536) and 3 DBS
#'   catalogs (one each for 78, 136, and 144). If trans.ranges = NULL, SBS 192
#'   and DBS 144 catalog will not be generated. Each catalog has attributes
#'   added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catalogs <- StrelkaSBSVCFFilesToCatalog(file, ref.genome = "hg19",
#'                                         trans.ranges = trans.ranges.GRCh37,
#'                                         region = "genome")
StrelkaSBSVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown") {
  vcfs <- ReadStrelkaSBSVCFs(files)
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  return(c(VCFsToSBSCatalogs(split.vcfs$SBS.vcfs, ref.genome, trans.ranges, region),
           VCFsToDBSCatalogs(split.vcfs$DBS.vcfs, ref.genome, trans.ranges, region)))
  }

#' Create SBS and DBS catalogs from Strelka SBS VCF files and plot them to PDF
#'
#' Create 3 SBS catalogs (96, 192, 1536) and 3 DBS catalogs (78, 136, 144) from the
#' Strelka SBS VCFs specified by \code{files} and plot them to PDF
#'
#' This function calls \code{\link{StrelkaSBSVCFFilesToCatalog}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param output.file The name of the PDF file to be produced.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and their graphs plotted to PDF
#'   with specified file name. If trans.ranges = NULL, SBS 192 and DBS 144
#'   catalog will not be generated and plotted. Each catalog has attributes
#'   added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed regions.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catalogs <- 
#'   StrelkaSBSVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19",
#'                                           trans.ranges = trans.ranges.GRCh37,
#'                                           region = "genome",
#'                                           output.file = file.path(tempdir(), 
#'                                                                   "StrelkaSBS.pdf")) 
StrelkaSBSVCFFilesToCatalogAndPlotToPdf <-
  function(files, ref.genome, trans.ranges = NULL, 
           region = "unknown", output.file) {
    catalogs <-
      StrelkaSBSVCFFilesToCatalog(files, ref.genome,
                                  trans.ranges, region)

    PlotCatalogToPdf(catalogs$catSBS96, 
                     file = sub(".pdf", ".SBS96Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    if (!is.null(trans.ranges)) {
      PlotCatalogToPdf(catalogs$catSBS192, 
                       file = sub(".pdf", ".SBS192Catalog.pdf", 
                                  output.file, ignore.case = TRUE))
      PlotCatalogToPdf(catalogs$catSBS192,
                       file = sub(".pdf", ".SBS12Catalog.pdf", 
                                  output.file, ignore.case = TRUE),
                       plot.SBS12 = TRUE)
      PlotCatalogToPdf(catalogs$catDBS144, 
                       file = sub(".pdf", ".DBS144Catalog.pdf", 
                                  output.file, ignore.case = TRUE))
    }
    PlotCatalogToPdf(catalogs$catSBS1536,
                     file = sub(".pdf", ".SBS1536Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    PlotCatalogToPdf(catalogs$catDBS78,
                     file = sub(".pdf", ".DBS78Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    PlotCatalogToPdf(catalogs$catDBS136,
                     file = sub(".pdf", ".DBS136Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    
    return(catalogs)
  }

#' Create ID (indel) catalog from Strelka ID VCF files
#'
#' Create ID (indel) catalog from the Strelka ID VCFs specified by \code{files}
#'
#' This function calls \code{\link{VCFsToIDCatalogs}}
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return An ID (indel) catalog with attributes added. See
#'   \code{\link{as.catalog}} for more details.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catID <- StrelkaIDVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                           region = "genome")
StrelkaIDVCFFilesToCatalog <- function(files, ref.genome, region = "unknown") {
  vcfs <- ReadStrelkaIDVCFs(files)
  return(VCFsToIDCatalogs(vcfs, ref.genome, region))
}

#' Create ID (indel) catalog from Strelka ID VCF files and plot them to PDF
#'
#' Create ID (indel) catalog from the Strelka ID VCFs specified by \code{files}
#' and plot them to PDF
#'
#' This function calls \code{\link{StrelkaIDVCFFilesToCatalog}} and
#' \code{\link{PlotCatalogToPdf}}
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param output.file The name of the PDF file to be produced.
#'
#' @return An ID (indel) catalog and its graph plotted to PDF with specified
#'   file name. The ID (indel) catalog has attributes added. See
#'   \code{\link{as.catalog}} for more details.
#'
#' @note In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catID <- 
#'   StrelkaIDVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                          region = "genome",
#'                                          output.file = file.path(tempdir(), 
#'                                                                  "StrelkaID.pdf"))
StrelkaIDVCFFilesToCatalogAndPlotToPdf <-
  function(files, ref.genome, region = "unknown", output.file) {
    catalog <-
      StrelkaIDVCFFilesToCatalog(files, ref.genome, region)
    PlotCatalogToPdf(catalog, file = sub(".pdf", ".IndelCatalog.pdf", 
                                         output.file, ignore.case = TRUE))
    return(catalog)
  }

#' Create SBS, DBS and Indel catalogs from Mutect VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the Mutect VCFs specified by \code{files}
#'
#' This function calls \code{\link{VCFsToSBSCatalogs}},
#' \code{\link{VCFsToDBSCatalogs}} and \code{\link{VCFsToIDCatalogs}}
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and ID catalog. If trans.ranges =
#'   NULL, SBS 192 and DBS 144 catalog will not be generated. Each catalog has
#'   attributes added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catalogs <- MutectVCFFilesToCatalog(file, ref.genome = "hg19", 
#'                                     trans.ranges = trans.ranges.GRCh37,
#'                                     region = "genome")
MutectVCFFilesToCatalog <-
  function(files, ref.genome, trans.ranges = NULL, region = "unknown") {
  vcfs <- ReadMutectVCFs(files)
  if (IsGRCm38(ref.genome)) {
    warning("Mouse MuTect VCFs files not tested; use at your own risk")
  }
  split.vcfs <- SplitListOfMutectVCFs(vcfs)
  return(c(VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, trans.ranges, region),
           VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, trans.ranges, region),
           list(catID = VCFsToIDCatalogs(split.vcfs$ID, ref.genome, region))))
  }

#' Create SBS, DBS and Indel catalogs from Mutect VCF files and plot them to PDF
#'
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
#' @param trans.ranges a \code{\link[data.table]{data.table}} which contains
#'   transcript range and strand information. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param output.file The name of the PDF file to be produced.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144), Indel catalog and their graphs
#'   plotted to PDF with specified file name. If trans.ranges = NULL, SBS 192
#'   and DBS 144 catalog will not be generated and plotted. Each catalog has
#'   attributes added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' catalogs <- 
#'   MutectVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                       trans.ranges = trans.ranges.GRCh37,
#'                                       region = "genome",
#'                                       output.file = file.path(tempdir(), "Mutect.pdf"))
MutectVCFFilesToCatalogAndPlotToPdf <-
  function(files, ref.genome, trans.ranges = NULL, 
           region = "unknown", output.file) {
    catalogs <-
      MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, region)

    PlotCatalogToPdf(catalogs$catSBS96, 
                     file = sub(".pdf", ".SBS96Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    if (!is.null(trans.ranges)) {
      PlotCatalogToPdf(catalogs$catSBS192, 
                       file = sub(".pdf", ".SBS192Catalog.pdf", 
                                  output.file, ignore.case = TRUE))
      PlotCatalogToPdf(catalogs$catSBS192,
                       file = sub(".pdf", ".SBS12Catalog.pdf", 
                                  output.file, ignore.case = TRUE),
                       plot.SBS12 = TRUE)
      PlotCatalogToPdf(catalogs$catDBS144, 
                       file = sub(".pdf", ".DBS144Catalog.pdf", 
                                  output.file, ignore.case = TRUE))
    }
    PlotCatalogToPdf(catalogs$catSBS1536,
                     file = sub(".pdf", ".SBS1536Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    PlotCatalogToPdf(catalogs$catDBS78,
                     file = sub(".pdf", ".DBS78Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    PlotCatalogToPdf(catalogs$catDBS136,
                     file = sub(".pdf", ".DBS136Catalog.pdf", 
                                output.file, ignore.case = TRUE))
    PlotCatalogToPdf(catalogs$catID,
                     file = sub(".pdf", ".IndelCatalog.pdf", 
                                output.file, ignore.case = TRUE))

    return(catalogs)
  }

#' @keywords internal
CanonicalizeDBS <- function(ref.vec, alt.vec) {
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  Canonicalize1DBS <- function(DBS) {
    if (DBS %in% ICAMS::catalog.row.order$DBS78) {
      return(DBS)
    } else {
      ref <- substr(DBS, 1, 2)
      alt <- substr(DBS, 3, 4)
      out <- paste0(revc(ref), revc(alt))
    }
    stopifnot(out %in% ICAMS::catalog.row.order$DBS78)
    return(out)
  }
  ret <- sapply(paste0(ref.vec, alt.vec), FUN = Canonicalize1DBS)
  return(ret)
}

#' @keywords internal
CanonicalizeQUAD <- function(quad) {
  canonical.ref <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  Canonicalize1QUAD <- function(quad) {
    if (quad %in% ICAMS::catalog.row.order$DBS136) {
      return(quad)
    } else {
      out <- revc(quad)
      stopifnot(out %in% ICAMS::catalog.row.order$DBS136)
      return(out)
    }
  }

  ret <- sapply(quad, FUN = Canonicalize1QUAD)
  return(ret)
}
