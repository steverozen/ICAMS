#' Extract the VAFs (variant allele frequencies) from a VCF file.
#'
#' @param vcf Said VCF as a data.frame.
#' 
#' @param name.of.VCF Name of the VCF file. 
#'   
#' @param tumor.col.name Optional. Name of the column in VCF which contains the
#'   tumor sample information. It \strong{must} have quotation marks. If
#'   \code{tumor.col.name} is equal to \code{NA}(default), this function will
#'   use the 10th column to calculate VAFs.
#'
#' @return A vector of VAFs, one for each row of \code{vcf}.
#'
#' @name GetVAF
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' MakeDataFrameFromStrelkaSBSVCF <- 
#'   getFromNamespace("MakeDataFrameFromStrelkaSBSVCF", "ICAMS")
#' df <- MakeDataFrameFromStrelkaSBSVCF(file)
#' vaf <- GetStrelkaVAF(df)
NULL

#' @keywords internal
RemoveRowsWithPoundSign <- function(df, file) {
  pound.chrom.idx <- which(df$CHROM == "#CHROM")
  if (length(pound.chrom.idx) > 0) {
    warning("Removing ", length(pound.chrom.idx), 
            " rows with #CHROM from file ", file)
    df1 <- df[-pound.chrom.idx, ]
    return(df1)
  } else {
    return(df)
  }
}

#' @keywords internal
RemoveRowsWithDuplicatedCHROMAndPOS <- function(df, file) {
  dups <- which(duplicated(df[, c("CHROM", "POS")]))
  if (length(dups) > 0) {
    dups2 <- which(duplicated(df[ , c("CHROM", "POS")], fromLast = TRUE))
    warning("In ", file, " ", 2 * length(dups), " rows out of ",
            nrow(df), " had duplicate CHROM and POS and were removed: ",
            dups, dups2)
    df1 <- df[-c(dups, dups2), ]
    return(df1)
  } else {
    return(df)
  }
}

#' Is there any column in \code{df} with name "strand"?
#' If there is, change its name to "strand_old" so that it will
#' conflict with code in other parts of ICAMS package.
#' 
#' @keywords internal
RenameColumnsWithNameStrand <- function(df) {
  if ("strand" %in% colnames(df)) {
    colnames(df)[which(colnames(df) == "strand")] <- "strand_old"
    warning('There is column in VCF which has name "strand", ',
            'it has been renamed to "strand_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  return(df)
}

#' Is there any column in df1 with name "VAF"?
#' If there is, change its name to "VAF_old" so that it will
#' conflict with code in other parts of ICAMS package.
#' 
#' @keywords internal
RenameColumnsWithNameVAF <- function(df) {
  if ("VAF" %in% colnames(df)) {
    colnames(df)[which(colnames(df) == "VAF")] <- "VAF_old"
    warning('There is column in VCF which has name "VAF", ',
            'it has been renamed to "VAF_old" so as ',
            'not to conflict with code in other parts of ICAMS package.')
  }
  return(df)
}

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
                 col.names = paste0("c", 1:100), as.is = TRUE)
  
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
  
  df1 <- RenameColumnsWithNameStrand(df1)
  df1 <- RenameColumnsWithNameVAF(df1)
  
  df1 <- RemoveRowsWithPoundSign(df1, file)
  df1 <- RemoveRowsWithDuplicatedCHROMAndPOS(df1, file)
  
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
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @keywords internal
ReadStrelkaIDVCF <- function(file) {
  df <- read.csv(file, header = FALSE, sep = "\t", quote = "",
                 col.names = paste0("c", 1:100), as.is = TRUE)

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
    stop("\nVCF does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(df1), collapse=" "))
  }
  control <- unique(df1[ , "FORMAT"])
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (all(each.base.col %in% colnames)) {
    stop("\nVCF does not appear to be a Strelka ID VCF, ", 
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
  stopifnot("data.frame" %in% class(vcf))
  if (!("TUMOR" %in% names(vcf)) ||
      !("FORMAT" %in% names(vcf))) {
    stop("\nVCF does not appear to be a Strelka VCF, column names are \n",
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
    stop("\nVCF does not appear to be a Strelka SBS VCF, ", 
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

#' Read in the data lines of a Variant Call Format (VCF) file created by Mutect
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
                 col.names = paste0("c", 1:100), as.is = TRUE)
  
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
  
  df1 <- RenameColumnsWithNameStrand(df1)
  df1 <- RenameColumnsWithNameVAF(df1)
  
  df1 <- RemoveRowsWithPoundSign(df1, file)
  df1 <- RemoveRowsWithDuplicatedCHROMAndPOS(df1, file)
  
  return(StandardChromName(df1))
}

#' Read in the data lines of a Variant Call Format (VCF) file created by
#'     Mutect
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#' 
#' @param name.of.vcf Name of the VCF file. If \code{NULL}(default), this
#'   function will remove all of the path up to and including the last path
#'   separator (if any) in \code{file} and file path without extensions (and the
#'   leading dot) will be used as the name of the VCF file.
#'   
#' @param tumor.col.name Name of the column in VCF which contains the tumor
#'   sample information. It \strong{must} have quotation marks. If
#'   \code{tumor.col.name} is equal to \code{NA}(default), this function will
#'   use the 10th column to calculate VAFs. See \code{\link{GetMutectVAF}} for
#'   more details.
#'   
#' @return A data frame storing mutation records of a VCF file with VAFs added.
#'
#' @keywords internal
ReadMutectVCF <- 
  function(file, name.of.VCF = NULL, tumor.col.name = NA) {
  df <- MakeDataFrameFromMutectVCF(file)
  if (is.null(name.of.VCF)) {
    vcf.name <- tools::file_path_sans_ext(basename(file))
  } else {
    vcf.name <- name.of.VCF
  }
  
  df$VAF <- GetMutectVAF(df, vcf.name, tumor.col.name)
  return(StandardChromName(df))
}

#' @rdname GetVAF
#'
#' @export
GetMutectVAF <- function(vcf, name.of.VCF = NULL, tumor.col.name = NA) {
  stopifnot("data.frame" %in% class(vcf))
  if (!any(grepl("/1", unlist(vcf[1, ])))) {
    stop("\nVCF does not appear to be a Mutect VCF, please check the data")
  }
  
  if (!is.na(tumor.col.name)) {
    if (!tumor.col.name %in% colnames(vcf)) {
      stop("\n", dQuote(tumor.col.name), 
           " is not one of the column names in vcf ",
           ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)))
    }
  }

  type1 <- c("F1R2", "F2R1")
  type2 <- c("REF_F1R2", "ALT_F1R2", "REF_F2R1", "ALT_F2R1")
  
  ExtractInfo <- function(idx, type, vector1, vector2) {
    pos <- match(type, unlist(strsplit(vector1[idx], ":")))
    values <- unlist(strsplit(vector2[idx], ":"))[pos]
  }
  
  CalculateVAF <- function(idx, list) {
    values <- list[[idx]]
    x <- as.integer(unlist(strsplit(values, ",")))
    vaf <- sum(x[2], x[4]) / sum(x)
  }
  
  GetVAFs <- function(type, vector1, vector2) {
    info <- lapply(1:length(vector1), FUN = ExtractInfo, type = type,
                   vector1 = vector1, vector2 = vector2)
    vafs <- sapply(1:length(info), FUN = CalculateVAF, list = info)
  }
  
  CheckAndReturnVAFs <- function(vafs) {
    idx.zero.vaf <- which(vafs == 0)
    if(length(idx.zero.vaf) == 0) {
      return(vafs)
    } else {
      warning("\nThere are rows which have zero VAF value in vcf ",
              ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)), "\n",
              "Please check and specify the correct column name for tumor sample ",
              "using argument 'tumor.col.name'")
      return(vafs)
    }
  }
  
  GetAndReturnVAFs <- function(type, vector1, vector2) {
    vafs <- GetVAFs(type, vector1, vector2)
    CheckAndReturnVAFs(vafs)
  }
  
  if (all(type1 %in% unlist(strsplit(vcf$FORMAT[1], ":")))) {
    if(is.na(tumor.col.name)) {
      GetAndReturnVAFs(type1, vcf$FORMAT, vcf[[10]])
    } else {
      GetAndReturnVAFs(type1, vcf$FORMAT, vcf[[tumor.col.name]])
    }
  } else if (all(type2 %in% unlist(strsplit(vcf$FORMAT[1], ":")))) {
    if(is.na(tumor.col.name)) {
      GetAndReturnVAFs(type2, vcf$FORMAT, vcf[[10]])
    } else {
      GetAndReturnVAFs(type2, vcf$FORMAT, vcf[[tumor.col.name]])
    }
  }
}

#' @title Split a mutect2 VCF into SBS, DBS, and ID VCFs, plus a list of other mutations
#'
#' @param vcf.df An in-memory data.frame representing a Mutect VCF, including
#'  VAFs, which are added by \code{\link{ReadMutectVCF}}.
#'
#' @return A list of in-memory objects with the elements:
#' \enumerate{
#'    \item \code{SBS.vcf}:   Data frame of pure SBS mutations -- no DBS or 3+BS mutations
#'    \item \code{DBS.vcf}:   Data frame of pure DBS mutations -- no SBS or 3+BS mutations
#'    \item{ThreePlus}: Data table with the key CHROM, LOW.POS, HIGH.POS and additional
#'    information (reference sequence, alternative sequence, context, etc.)
#'    Additional information not fully implemented at this point because of
#'    limited immediate biological interest.
#'    \item{multiple.alt}: Rows that were removed before processing because they had
#'    more than one alternate allele.
#'    }
#'
#'
#' @keywords internal
SplitOneMutectVCF <- function(vcf.df) {
  # Mutect VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)
  multiple.alt.df <- vcf.df[multiple.alt, ]
  
  if (length(multiple.alt) != 0) {
    df <- vcf.df[-multiple.alt, ]
  } else {
    df <- vcf.df
  }

  SBS.df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]

  DBS.df <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]

  other.df <- df[nchar(df$REF) > 2 & nchar(df$ALT) == nchar(df$REF), ]

  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]
  complex.indels.to.remove <- 
    which(substr(ID.df$REF, 1, 1) != substr(ID.df$ALT, 1, 1))
  complex.indels <- ID.df[complex.indels.to.remove, ]
  if (length(complex.indels.to.remove > 0)) {
    ID.df <- ID.df[-complex.indels.to.remove, ]
  }
  
  other.df2 <- rbind(other.df, complex.indels)

  return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
              other = other.df2, multiple.alt = multiple.alt.df))

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
  other.subs <- lapply(v1, function(x) x$other)
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
  chr.names <- CheckAndFixChrNames(vcf.df = df, ref.genome = ref.genome)
  
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
  
  # Rename some of the columns in dt3
  setnames(dt3, 
           old = c("start", "end", "strand", "Ensembl.gene.ID", "gene.symbol"), 
           new = c("trans.start.pos", "trans.end.pos", "trans.strand", 
                   "trans.Ensembl.gene.ID", "trans.gene.symbol"))

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
#' @return A list of in-memory objects with the elements:
#' \enumerate{
#'    \item \code{SBS.vcf}:   Data frame of pure SBS mutations -- no DBS or 3+BS mutations
#'    \item \code{DBS.vcf}:   Data frame of pure DBS mutations -- no SBS or 3+BS mutations
#'    \item{ThreePlus}: Data table with the key CHROM, LOW.POS, HIGH.POS and additional
#'    information (reference sequence, alternative sequence, context, etc.)
#'    Additional information not fully implemented at this point because of
#'    limited immediate biological interest.
#'    \item{multiple.alt}: Rows that were removed before processing because they had
#'    more than one alternate allele.
#'    }
#'
#' @keywords internal
SplitStrelkaSBSVCF <- function(vcf.df, max.vaf.diff = 0.02) {
  stopifnot("data.frame" %in% class(vcf.df))
  
  # Strelka SBS VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)
  multiple.alt.df <- vcf.df[multiple.alt, ]
  
  if (length(multiple.alt) != 0) {
    vcf.df <- vcf.df[-multiple.alt, ]
  }

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
    empty <- vcf.df[-(1:nrow(vcf.df)), ]
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
  out.SBS.df <- 
    as.data.frame(out.SBS.dt2[, c("POS.plus.one", "delete.flag") := NULL])
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
              ThreePlus = other.ranges, multiple.alt = multiple.alt.df))
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
#' @return A list of in-memory objects with the elements: \enumerate{
#'    \item \code{SBS.vcfs}:  List of Data frames of pure SBS mutations -- no DBS or 3+BS mutations
#'    \item \code{DBS.vcfs}:  List of Data frames of pure DBS mutations -- no SBS or 3+BS mutations
#'    \item \code{ThreePlus}: List of Data tables with the key CHROM, LOW.POS, HIGH.POS and additional
#'    information (reference sequence, alternative sequence, context, etc.)
#'    Additional information not fully implemented at this point because of
#'    limited immediate biological interest.
#'    \item \code{multiple.alt} Rows with multiple alternate alleles (removed from
#'    \code{SBS.vcfs} etc.)\
#'    }
#'
#' @keywords internal
SplitListOfStrelkaSBSVCFs <- function(list.of.vcfs) {
  split.vcfs <- lapply(list.of.vcfs, FUN = SplitStrelkaSBSVCF)
  SBS.vcfs   <- lapply(split.vcfs, function(x) x$SBS.vcf)
  DBS.vcfs   <- lapply(split.vcfs, function(x) x$DBS.vcf)
  ThreePlus  <- lapply(split.vcfs, function(x) x$ThreePlus)
  mult.alt   <- lapply(split.vcfs, function(x) x$multiple.alt)
  return(list(SBS.vcfs = SBS.vcfs,
              DBS.vcfs = DBS.vcfs,
              ThreePlus = ThreePlus,
              multiple.alt = mult.alt))
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
#' @importFrom utils write.csv
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
    write.csv(vcf[error.rows, ], file = temp)
    stop("Seqence context of reference allele is inconsistent,",
         "see file ", temp)
  }
}

#' Read Strelka SBS (single base substitutions) VCF files.
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @param names.of.VCFs Character vector of names of the VCF files. The order
#'   of names in \code{names.of.VCFs} should match the order of VCF file paths
#'   in \code{files}. If \code{NULL}(default), this function will remove all of
#'   the path up to and including the last path separator (if any) and file
#'   paths without extensions (and the leading dot) will be used as the names of
#'   the VCF files.
#'   
#' @return A list of vcfs from \code{files}.
#'
#' @keywords internal
ReadStrelkaSBSVCFs <- function(files, names.of.VCFs = NULL) {
  vcfs <- lapply(files, FUN = ReadStrelkaSBSVCF)
  if (is.null(names.of.VCFs)) {
    names(vcfs) <- tools::file_path_sans_ext(basename(files))
  } else {
    CheckNamesOfVCFs(files, names.of.VCFs)
    names(vcfs) <- names.of.VCFs
  }
  
  return(vcfs)
}

#' Read Mutect VCF files.
#'
#' @param files Character vector of file paths to the VCF files.
#'
#' @param names.of.VCFs Character vector of names of the VCF files. The order of
#'   names in \code{names.of.VCFs} should match the order of VCF file paths in
#'   \code{files}. If \code{NULL}(default), this function will remove all of the
#'   path up to and including the last path separator (if any) in \code{files}
#'   and file paths without extensions (and the leading dot) will be used as the
#'   names of the VCF files.
#'   
#' @param tumor.col.names Character vector of column names in VCFs which contain
#'   the tumor sample information. The order of names in \code{tumor.col.names}
#'   should match the order of VCFs specified in \code{files}. If
#'   \code{tumor.col.names} is equal to \code{NA}(default), this function will
#'   use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'   
#' @return A list of vcfs from \code{files}.
#'
#' @keywords internal
ReadMutectVCFs <- 
  function(files, names.of.VCFs = NULL, tumor.col.names = NA) {
  if (is.null(names.of.VCFs)) {
    vcfs.names <- tools::file_path_sans_ext(basename(files))
  } else {
    CheckNamesOfVCFs(files, names.of.VCFs)
    vcfs.names <- names.of.VCFs
  }
  num.of.files <- length(files)
  if (all(is.na(tumor.col.names))) {
    tumor.col.names <- rep(NA, num.of.files)
  }
  
  GetMutectVCFs <- function(idx, files, vector1, vector2) {
    ReadMutectVCF(file = files[idx], name.of.VCF = vector1[idx],
                  tumor.col.name = vector2[idx])
  }
  
  vcfs <- lapply(1:num.of.files, FUN = GetMutectVCFs, 
                 files, vcfs.names, tumor.col.names)
  names(vcfs) <- vcfs.names
  return(vcfs)
}

#' Add sequence context and transcript information to an in-memory SBS VCF.
#' 
#' @param SBS.vcf An in-memory SBS VCF as a \code{data.frame}.
#' 
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link{BSgenome}} object 
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'
#' @return An in-memory SBS VCF as a \code{data.table}. This has been annotated
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
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateSBSVCF <- function(SBS.vcf, ref.genome, trans.ranges = NULL) {
  SBS.vcf <- AddSeqContext(SBS.vcf, ref.genome = ref.genome)
  
  # Delete the rows of SBS if the extracted sequence contains "N"
  idx <- grep("N", substr(SBS.vcf$seq.21bases, 9, 13))
  if (!length(idx) == 0) {
    SBS.vcf <- SBS.vcf[-idx, ]
    message(
      'Rows in the SBS vcf where surrounding sequence contains "N" ',
      'have been deleted so as not to conflict with downstream processing')
  }
  
  CheckSeqContextInVCF(SBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    SBS.vcf <- AddTranscript(SBS.vcf, trans.ranges)
  }
  return(as.data.table(SBS.vcf))
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
CreateOneColSBSMatrix <- function(vcf, sample.id = "count") {
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
  mismatches <- which(vcf$REF != substr(vcf$seq.21bases, 11, 11))
  if (length(mismatches) != 0) {
    stop("\nSample ", sample.id, 
         ":\nThe reference base in ref.genome does not match the ", 
         "reference base in ", length(mismatches),
         "rows in the VCF file.\n",
         "Please check the ref.genome argument.")
  }
  
  # Create 2 new columns that show the 3072 and 1536 mutation type
  context <- substr(vcf$seq.21bases, 9, 13)
  vcf$mutation <- paste0(context, vcf$ALT)

  # PyrPenta maps to strand-agnostic category
  # e.g. ATGCT>T "ATGCTT" maps to AGCAT>A, "AGCATA"
  vcf$pyr.mut <- PyrPenta(vcf$mutation)

  # One SBS mutation can be represented by more than 1 row in vcf 
  # after annotation by AddTranscript if the mutation position falls 
  # is in multiple transcripts. When creating the 1536 and 96 catalog,
  # we only need to count these mutations once.
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
  
  if (is.null(vcf$trans.strand)) {
     return(list(catSBS96 = mat96, catSBS1536 = mat1536))
  }
  
  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 192 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One SBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 192 catalog, we only need to count these mutations once.
  vcf3 <- vcf2[, .(REF = REF[1], mutation = mutation[1], 
                   trans.strand = trans.strand[1]),
               by = .(CHROM, ALT, POS)]

  # Create the 192 catalog matrix
  tab192  <- table(paste0(substr(vcf3$mutation, 2, 4),
                          substr(vcf3$mutation, 6, 6)),
                   vcf3$trans.strand,
                   useNA = "ifany")
  stopifnot(sum(tab192) == nrow(vcf3))
  dt192 <- as.data.table(tab192)
  colnames(dt192) <- c("rn", "trans.strand", "count")
  dt192 <- dt192[!is.na(trans.strand)]
  dt192[trans.strand == "-", rn := RevcSBS96(rn)]
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

#' Add sequence context and transcript information to an in-memory DBS VCF.
#' 
#' @param DBS.vcf An in-memory DBS VCF as a \code{data.frame}.
#' 
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
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
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' DBS.vcf <- list.of.vcfs$DBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateDBSVCF <- function(DBS.vcf, ref.genome, trans.ranges = NULL) {
  DBS.vcf <- AddSeqContext(DBS.vcf, ref.genome = ref.genome)
  
  # Delete the rows of DBS if the extracted sequence contains "N"
  idx <- grep("N", substr(DBS.vcf$seq.21bases, 10, 13))
  if (!length(idx) == 0) {
    DBS.vcf <- DBS.vcf[-idx, ]
    message(
      'Rows in the DBS vcf where surrounding sequence contains "N" ',
      'have been deleted so as not to conflict with downstream processing')
  }
  
  CheckSeqContextInVCF(DBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    DBS.vcf <- AddTranscript(DBS.vcf, trans.ranges)
  }
  return(as.data.table(DBS.vcf))
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
#' @return A list of three 1-column matrices with the names \code{catDBS78},
#'   \code{catDBS144}, and \code{catDBS136}. If trans.ranges is NULL,
#'   \code{catDBS144} is not generated. Do not rely on the order of elements in
#'   the list.
#'
#' @note DBS 144 catalog only contains mutations in transcribed regions.
#'
#' @keywords internal
CreateOneColDBSMatrix <- function(vcf, sample.id = "count") {
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
  # AnnotateDBSVCF function if the mutation position falls into the range of
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
  
  if (is.null(vcf$trans.strand)) {
    return(list(catDBS78 = DBS.mat.78, catDBS136 = DBS.mat.136))
  }
  
  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 144 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One DBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 144 catalog, we only need to count these mutations once.
  vcf3 <- vcf2[, .(REF = REF[1], trans.strand = trans.strand[1]),
               by = .(CHROM, ALT, POS)]

  # Create the 144 DBS catalog matrix
  # There are 144 stranded DBSs: 4 X 4 sources and 3 X 3 alternates;
  # 4 x 4 x 3 x 3 = 144.
  tab.DBS.144  <-
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
  colnames(DBS.mat.144)<- sample.id

  return(list(catDBS78 = DBS.mat.78, catDBS144 = DBS.mat.144,
              catDBS136 = DBS.mat.136))
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
#' @inheritParams MutectVCFFilesToCatalogAndPlotToPdf
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and their graphs plotted to PDF
#'   with specified file name. If trans.ranges = NULL, SBS 192 and DBS 144
#'   catalog will not be generated and plotted. Each catalog has attributes
#'   added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed regions.
#' 
#' @inheritSection MutectVCFFilesToCatalogAndPlotToPdf Comments
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     StrelkaSBSVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19",
#'                                             trans.ranges = trans.ranges.GRCh37,
#'                                             region = "genome",
#'                                             output.file = 
#'                                             file.path(tempdir(), "StrelkaSBS"))}
StrelkaSBSVCFFilesToCatalogAndPlotToPdf <- function(files, 
                                                    ref.genome, 
                                                    trans.ranges = NULL, 
                                                    region = "unknown", 
                                                    names.of.VCFs = NULL, 
                                                    output.file = "") {
    
    catalogs <-
      StrelkaSBSVCFFilesToCatalog(files, ref.genome,
                                  trans.ranges, region, names.of.VCFs)
    
    if (output.file != "") output.file <- paste0(output.file, ".")
    
    for (name in names(catalogs)) {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, name, ".pdf"))
      if (name == "catSBS192") {
        PlotCatalogToPdf(catalogs[[name]],
                         file = paste0(output.file, "SBS12.pdf"),
                         plot.SBS12 = TRUE)
      }
    }
    
    return(catalogs)
  }

#' Create ID (small insertion and deletion) catalog from Strelka ID VCF files
#' and plot them to PDF
#'
#' Create ID (small insertion and deletion) catalog from the Strelka ID VCFs
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
#' @return A list whose first element is an ID (small insertion and deletion)
#'   catalog with its graph plotted to PDF with specified file name. The ID
#'   catalog has class attribute "IndelCatalog" added. See
#'   \code{\link{as.catalog}} for more details. The second element of the
#'   returned list is a list of further annotated VCFs.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catID <- 
#'     StrelkaIDVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                            region = "genome",
#'                                            output.file = 
#'                                            file.path(tempdir(), "StrelkaID"))}
#'                                                                    
StrelkaIDVCFFilesToCatalogAndPlotToPdf <- function(files, 
                                                   ref.genome, 
                                                   region = "unknown", 
                                                   names.of.VCFs = NULL, 
                                                   output.file = "") {
    
    list <-
      StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs)
    
    if (output.file != "") output.file <- paste0(output.file, ".")
    
    PlotCatalogToPdf(list$catalog, file = paste0(output.file, "catID", ".pdf"))
    
    return(list)
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
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link{BSgenome}} object 
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
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
#' @param tumor.col.names Optional. Character vector of column names in VCFs
#'   which contain the tumor sample information. The order of names in
#'   \code{tumor.col.names} should match the order of VCFs specified in
#'   \code{files}. If \code{tumor.col.names} is equal to \code{NA}(default),
#'   this function will use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'
#' @param output.file Optional. The base name of the PDF files to be produced;
#'   multiple files will be generated, each ending in \eqn{x}\code{.pdf}, where
#'   \eqn{x} indicates the type of catalog plotted in the file.
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144), Indel catalog and their graphs
#'   plotted to PDF with specified file name. If trans.ranges = NULL, SBS 192
#'   and DBS 144 catalog will not be generated and plotted. Each catalog has
#'   attributes added. See \code{\link{as.catalog}} for more details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
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
#'                       "Mutect.GRCh37.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     MutectVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                         trans.ranges = trans.ranges.GRCh37,
#'                                         region = "genome",
#'                                         output.file = 
#'                                         file.path(tempdir(), "Mutect"))}
MutectVCFFilesToCatalogAndPlotToPdf <- function(files, 
                                                ref.genome, 
                                                trans.ranges = NULL, 
                                                region = "unknown", 
                                                names.of.VCFs = NULL, 
                                                tumor.col.names = NA,
                                                output.file = "") {
    
    catalogs <-
      MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                              region, names.of.VCFs, tumor.col.names)
    
    if (output.file != "") output.file <- paste0(output.file, ".")
    
    for (name in names(catalogs)) {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, name, ".pdf"))
      if (name == "catSBS192") {
        PlotCatalogToPdf(catalogs[[name]],
                         file = paste0(output.file, "SBS12.pdf"),
                         plot.SBS12 = TRUE)
    }
    }
    
    return(catalogs)
}

#' Create a zip file which contains catalogs and plot PDFs from Mutect VCF files
#'
#' Create 3 SBS catalogs (96, 192, 1536), 3 DBS catalogs (78, 136, 144) and
#' Indel catalog from the Mutect VCFs specified by \code{dir}, save the catalogs
#' as CSV files, plot them to PDF and generate a zip archive of all the output files.
#'
#' This function calls \code{\link{MutectVCFFilesToCatalog}},
#' \code{\link{PlotCatalogToPdf}}, \code{\link{WriteCatalog}} and
#' \code{\link[zip]{zipr}}.
#'
#' @param dir Pathname of the directory which contains the Mutect VCF files.
#'   Each Mutect VCF \strong{must} have a file extension ".vcf" (case
#'   insensitive) and share the \strong{same} \code{ref.genome} and
#'   \code{region}.
#'   
#' @param zipfile Pathname of the zip file to be created.    
#'   
#' @param ref.genome  A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}. 
#'   
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link{BSgenome}} object 
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
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
#'   The order of names in \code{names.of.VCFs} should match the order of VCFs
#'   listed in \code{dir}. If \code{NULL}(default), this function will remove
#'   all of the path up to and including the last path separator (if any) in
#'   \code{dir} and file paths without extensions (and the leading dot) will be
#'   used as the names of the VCF files.
#'   
#' @param tumor.col.names Optional. Character vector of column names in VCFs which contain
#'   the tumor sample information. The order of names in \code{tumor.col.names}
#'   should match the order of VCFs listed in \code{dir}. If
#'   \code{tumor.col.names} is equal to \code{NA}(default), this function will
#'   use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'   
#' @param base.filename Optional. The base name of the CSV and PDF files to be
#'   produced; multiple files will be generated, each ending in
#'   \eqn{x}\code{.csv} or \eqn{x}\code{.pdf}, where \eqn{x} indicates the type
#'   of catalog.
#'
#' @importFrom utils glob2rx 
#' 
#' @importFrom zip zipr 
#'
#' @return  A list of 3 SBS catalogs (one each for 96, 192, and 1536), 3 DBS
#'   catalogs (one each for 78, 136, and 144) and Indel catalog. If trans.ranges
#'   = NULL, SBS 192 and DBS 144 catalog will not be generated and plotted. Each
#'   catalog has attributes added. See \code{\link{as.catalog}} for more
#'   details.
#'
#' @note SBS 192 and DBS 144 catalogs include only mutations in transcribed
#'   regions. In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+.
#' 
#' @section Comments:
#' To add or change attributes of the catalog, you can use function \code{\link[base]{attr}}. \cr
#' For example, \code{attr(catalog, "abundance") <- custom.abundance}.
#' 
#' @export
#' 
#' @examples 
#' dir <- c(system.file("extdata/Mutect-vcf",
#'                      package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catalogs <- 
#'     MutectVCFFilesToZipFile(dir, 
#'                             zipfile = paste0(tempdir(), "/test.zip"),
#'                             ref.genome = "hg19", 
#'                             trans.ranges = trans.ranges.GRCh37,
#'                             region = "genome",
#'                             base.filename = "Mutect")
#'   unlink(paste0(tempdir(), "/test.zip"))}
MutectVCFFilesToZipFile <- function(dir,
                                    zipfile, 
                                    ref.genome, 
                                    trans.ranges = NULL, 
                                    region = "unknown", 
                                    names.of.VCFs = NULL, 
                                    tumor.col.names = NA,
                                    base.filename = ""){
  files <- list.files(path = dir, pattern = "\\.vcf$", 
                      full.names = TRUE, ignore.case = TRUE)
  catalogs <-
    MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                            region, names.of.VCFs, tumor.col.names)
  
  if (base.filename != "") {
    output.file <- paste0(tempdir(), "\\", base.filename, ".")
  } else {
    output.file <- paste0(tempdir(), "\\", base.filename)
  }
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192") {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
  invisible(catalogs)
}

#' @keywords internal
CanonicalizeDBS <- function(ref.vec, alt.vec) {

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

#' @keywords internal
CheckNamesOfVCFs <- function(files, names.of.VCFs) {
  stopifnot(inherits(names.of.VCFs, "character"))
  if (length(files) != length(names.of.VCFs)) {
    stop("\nThe number of names in names.of.VCFs does not match ",
         "the number of VCF files")
  }
}

#' @keywords internal
InferTransRanges <- function(ref.genome, trans.ranges) {
  if (!is.null(trans.ranges)) {
    return(trans.ranges)
  } else {
    if (IsGRCh37(ref.genome)) {
      return(trans.ranges.GRCh37)
    } else if (IsGRCh38(ref.genome)) {
      return(trans.ranges.GRCh38)
    } else if (IsGRCm38(ref.genome)) {
      return(trans.ranges.GRCm38)
    } else {
      return(trans.ranges)
    }
  }
}