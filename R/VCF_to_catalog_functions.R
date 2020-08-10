#' Extract the VAFs (variant allele frequencies) and read depth information from
#' a VCF file
#' 
#' @param vcf Said VCF as a data.frame.
#' 
#' @param name.of.VCF Name of the VCF file. 
#'   
#' @param tumor.col.name Optional. Only applicable to \strong{Mutect} VCF. Name
#'   of the column in \strong{Mutect} VCF which contains the tumor sample
#'   information. It \strong{must} have quotation marks. If
#'   \code{tumor.col.name} is equal to \code{NA}(default), this function will
#'   use the 10th column to calculate VAFs.
#'
#' @return The original \code{vcf} with two additional columns added which
#'   contain the VAF(variant allele frequency) and read depth information.
#'
#' @name GetVAF
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' MakeDataFrameFromVCF <- getFromNamespace("MakeDataFrameFromVCF", "ICAMS")
#' df <- MakeDataFrameFromVCF(file)$df
#' df1 <- GetStrelkaVAF(df)
NULL

#' @keywords internal
RemoveRowsWithPoundSign <- function(df, file) {
  pound.chrom.idx <- which(df$CHROM == "#CHROM")
  if (length(pound.chrom.idx) > 0) {
    warning("In ", file, " ", length(pound.chrom.idx), " row out of ",
            nrow(df), " had value #CHROM in column 'CHROM' and were removed. ",
            "See discarded.variants in the return value for more details")
    df1 <- df[-pound.chrom.idx, ]
    return(list(df = df1, discarded.variants = df[pound.chrom.idx, ]))
  } else {
    return(list(df = df))
  }
}

#' @keywords internal
RemoveRowsWithDuplicatedCHROMAndPOS <- function(df, file) {
  dups <- which(duplicated(df[, c("CHROM", "POS")]))
  if (length(dups) > 0) {
    dups2 <- which(duplicated(df[ , c("CHROM", "POS")], fromLast = TRUE))
    warning("In ", file, " ", 2 * length(dups), " row out of ",
            nrow(df), " had duplicate CHROM and POS and were removed. ",
            "See the discarded variants in the return value for more details")
    df1 <- df[-c(dups, dups2), ]
    return(list(df = df1, discarded.variants = df[c(dups, dups2), ]))
  } else {
    return(list(df = df))
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
#' @param name.of.VCF Name of the VCF file. If \code{NULL}(default), this
#'   function will remove all of the path up to and including the last path
#'   separator (if any) in \code{file} and file path without extensions (and the
#'   leading dot) will be used as the name of the VCF file.
#'   
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @return A \strong{list} whose first element \code{df} is a data frame storing
#'   data lines of the VCF file with two additional columns added which contain
#'   the VAF(variant allele frequency) and read depth information. A second element
#'   \code{discarded.variants} \strong{only} appears if there are variants that 
#'   are excluded from the analysis.
#'
#' @keywords internal
ReadStrelkaSBSVCF <- function(file, name.of.VCF = NULL, 
                              suppress.discarded.variants.warnings = TRUE) {
  
  retval <- MakeDataFrameFromVCF(file, suppress.discarded.variants.warnings)
  
  if (is.null(name.of.VCF)) {
    vcf.name <- tools::file_path_sans_ext(basename(file))
  } else {
    vcf.name <- name.of.VCF
  }
  
  # Create an empty data frame for discarded variants
  discarded.variants <- retval$df[0, ]
  
  if (!is.null(retval$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval$discarded.variants)
  }
  
  if (suppress.discarded.variants.warnings == TRUE) {
    retval2 <- suppressWarnings(StandardChromName(retval$df, file)) 
  } else {
    retval2 <- StandardChromName(retval$df, file)
  }
  if (!is.null(retval2$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
  }
  
  df1 <- GetStrelkaVAF(retval2$df, vcf.name)
  
  if (nrow(discarded.variants) == 0) {
    return(list(df = df1))
  } else {
    return(list(df = df1, discarded.variants = discarded.variants))
  }
}

#' Read in the data lines of a Variant Call Format (VCF) file
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#' 
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @return A list of elements:
#' * \code{df}: A data frame storing data lines from the original VCF file.
#' * \code{discarded.variants}: \strong{Only appearing when} there are variants
#' that are excluded in the analysis.
#' @md
#'
#' @keywords internal
MakeDataFrameFromVCF <- 
  function(file, suppress.discarded.variants.warnings = TRUE) {
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
  
  # Create an empty data frame for discarded variants
  discarded.variants <- df1[0, ]
  
  if (suppress.discarded.variants.warnings == TRUE) {
    retval <- suppressWarnings(RemoveRowsWithPoundSign(df1, file))
  } else {
    retval <- RemoveRowsWithPoundSign(df1, file)
  }
  
  if (!is.null(retval$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval$discarded.variants)
  }
  
  if (suppress.discarded.variants.warnings == TRUE) {
    retval2 <- 
      suppressWarnings(RemoveRowsWithDuplicatedCHROMAndPOS(retval$df, file))
  } else {
    retval2 <- RemoveRowsWithDuplicatedCHROMAndPOS(retval$df, file)
  }
  
  if (!is.null(retval2$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
  }
  
  if (nrow(discarded.variants) == 0) {
    return(list(df = retval2$df))
  } else {
    return(list(df = retval2$df, discarded.variants = discarded.variants))
  }
}

#' Read in the data lines of an ID VCF created by Strelka version 1
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#' 
#' @param name.of.VCF Name of the VCF file. If \code{NULL}(default), this
#'   function will remove all of the path up to and including the last path
#'   separator (if any) in \code{file} and file path without extensions (and the
#'   leading dot) will be used as the name of the VCF file.
#'   
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @return A \strong{list} whose first element \code{df} is a data frame storing
#'   data lines of the VCF file. A second element \code{discarded.variants}
#'   \strong{only} appears if there are variants that are excluded from the
#'   analysis.
#'   
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @keywords internal
ReadStrelkaIDVCF <- function(file, name.of.VCF = NULL,
                             suppress.discarded.variants.warnings = TRUE) {
  retval <- MakeDataFrameFromVCF(file, suppress.discarded.variants.warnings)
    
  # Get the name of VCF
  if (is.null(name.of.VCF)) {
    vcf.name <- tools::file_path_sans_ext(basename(file))
  } else {
    vcf.name <- name.of.VCF
  }
  
  df1 <- retval$df
  # Check whether the input VCF is a Strelka ID VCF
  if (!("TUMOR" %in% names(df1)) ||
      !("FORMAT" %in% names(df1))) {
    stop("\nVCF ", dQuote(vcf.name),
         " does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(df1), collapse=" "))
  }
  control <- unique(df1[ , "FORMAT"])
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (all(each.base.col %in% colnames)) {
    stop("\nVCF ", dQuote(vcf.name),
         " does not appear to be a Strelka ID VCF, ", 
         "the value of column FORMAT is \n", 
         control)
  }
  
  # Create an empty data frame for discarded variants
  discarded.variants <- df1[0, ]
  if (!is.null(retval$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval$discarded.variants)
  }
  
  if (suppress.discarded.variants.warnings == TRUE) {
    retval2 <- suppressWarnings(StandardChromName(df1, file)) 
  } else {
    retval2 <- StandardChromName(df1, file)
  }
  if (!is.null(retval2$discarded.variants)) {
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
  }

  if (nrow(discarded.variants) == 0) {
    return(list(df = retval2$df))
  } else {
    return(list(df = retval2$df, discarded.variants = discarded.variants))
  }
}

#' @rdname GetVAF
#'
#' @export
GetStrelkaVAF <-function(vcf, name.of.VCF = NULL) {
  stopifnot("data.frame" %in% class(vcf))
  if (!("TUMOR" %in% names(vcf)) ||
      !("FORMAT" %in% names(vcf))) {
    stop("\nVCF ", 
         ifelse(is.null(name.of.VCF), "", paste0(dQuote(name.of.VCF), " ")), 
         "does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(vcf), collapse=" "))
  }
  
  TUMOR <- vcf[ , "TUMOR"]
  control <- unique(vcf[ , "FORMAT"])
  alt     <- vcf[ , "ALT"]
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  values <- strsplit(TUMOR, split=":", fixed=TRUE)
  vaf <- numeric(nrow(vcf))
  read.depth <- integer(nrow(vcf))
  
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (!all(each.base.col %in% colnames)) {
    stop("\nVCF ", 
         ifelse(is.null(name.of.VCF), "", paste0(dQuote(name.of.VCF), " ")), 
         "does not appear to be a Strelka SBS VCF, ", 
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
    read.depth[i] <- total.read.count
  }
  
  return(cbind(vcf, VAF = vaf, read.depth = read.depth))
}

#' Read in the data lines of a Variant Call Format (VCF) file created by
#'     Mutect
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#' 
#' @param name.of.VCF Name of the VCF file. If \code{NULL}(default), this
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
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#' 
#' @return A \strong{list} whose first element \code{df} is a data frame storing
#'   data lines of the VCF file with two additional columns added which contain
#'   the VAF(variant allele frequency) and read depth information. A second element
#'   \code{discarded.variants} \strong{only} appears if there are variants that 
#'   are excluded from the analysis.
#'   
#' @keywords internal
ReadMutectVCF <- 
  function(file, name.of.VCF = NULL, tumor.col.name = NA,
           suppress.discarded.variants.warnings = TRUE) {
    retval <- MakeDataFrameFromVCF(file, suppress.discarded.variants.warnings)
    if (is.null(name.of.VCF)) {
      vcf.name <- tools::file_path_sans_ext(basename(file))
    } else {
      vcf.name <- name.of.VCF
    }
    
    # Create an empty data frame for discarded variants
    discarded.variants <- retval$df[0, ]
    
    if (!is.null(retval$discarded.variants)) {
      discarded.variants <- 
        dplyr::bind_rows(discarded.variants, retval$discarded.variants)
    }
    
    if (suppress.discarded.variants.warnings == TRUE) {
      retval2 <- suppressWarnings(StandardChromName(retval$df, file)) 
    } else {
      retval2 <- StandardChromName(retval$df, file)
    }
    if (!is.null(retval2$discarded.variants)) {
      discarded.variants <- 
        dplyr::bind_rows(discarded.variants, retval2$discarded.variants)
    }
    
    df1 <- GetMutectVAF(retval2$df, vcf.name, tumor.col.name)
    
    if (nrow(discarded.variants) == 0) {
      return(list(df = df1))
    } else {
      return(list(df = df1, discarded.variants = discarded.variants))
    }
  }

#' @rdname GetVAF
#'
#' @export
GetMutectVAF <- function(vcf, name.of.VCF = NULL, tumor.col.name = NA) {
  stopifnot("data.frame" %in% class(vcf))
  
  # Specify the possible variable names in Mutect VCF that stores count of reads
  # information
  type1 <- c("F1R2", "F2R1")
  type2 <- c("REF_F1R2", "ALT_F1R2", "REF_F2R1", "ALT_F2R1")
  
  if (!all(sapply(type1, FUN = grepl, x = vcf$FORMAT[1])) &&
      !all(sapply(type2, FUN = grepl, x = vcf$FORMAT[1]))) {
    
    warning("\nVCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " does not appear to be a Mutect VCF, please check the data")
    
    vcf$VAF <- NA
    vcf$read.depth <- NA
    return(vcf)
    
    #stop("\nVCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
    #     " does not appear to be a Mutect VCF, please check the data")
  }

  
  #if (!any(grepl("/1", unlist(vcf[1, ]), fixed = TRUE)) && 
  #    !any(grepl("|1", unlist(vcf[1, ]), fixed = TRUE))) {
  #  stop("\nVCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
  #       " does not appear to be a Mutect VCF, please check the data")
  #}
  
  if (!is.na(tumor.col.name)) {
    if (!tumor.col.name %in% colnames(vcf)) {
      stop("\n", dQuote(tumor.col.name), 
           " is not one of the column names in vcf ",
           ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)))
    }
  }

  ExtractInfo <- function(idx, type, vector1, vector2) {
    pos <- match(type, unlist(strsplit(vector1[idx], ":")))
    values <- unlist(strsplit(vector2[idx], ":"))[pos]
  }
  
  CalculateVAF <- function(idx, list) {
    values <- list[[idx]]
    x <- as.integer(unlist(strsplit(values, ",")))
    return(data.frame(VAF = sum(x[2], x[4]) / sum(x), read.depth = sum(x)))
  }
  
  GetVAFs <- function(type, vector1, vector2) {
    info <- lapply(1:length(vector1), FUN = ExtractInfo, type = type,
                   vector1 = vector1, vector2 = vector2)
    vafs <- lapply(1:length(info), FUN = CalculateVAF, list = info)
    do.call("rbind", vafs)
  }
  
  CheckAndReturnVAFs <- function(vafs) {
    idx.zero.vaf <- which(vafs$VAF == 0)
    if(length(idx.zero.vaf) == 0) {
      return(cbind(vcf, vafs))
    } else {
      zero.vaf.row <- length(idx.zero.vaf)
      total.vaf.row <- nrow(vafs)
      warning("\nThere are ", zero.vaf.row, " out of total ", total.vaf.row, 
              " rows which have zero VAF value in vcf ",
              ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)), "\n",
              "Please check the data and if necessary, specify the correct ", 
              "column name for tumor sample using argument 'tumor.col.name'")
      return(cbind(vcf, vafs))
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

#' @rdname GetVAF
#'
#' @export
GetFreebayesVAF <- function(vcf, name.of.VCF = NULL) {
  # From header of freebayes vcf file:
  # SRF = "Number of reference observations on the forward strand"
  # SRR = "Number of reference observations on the reverse strand"
  # SAF = "Number of alternate observations on the forward strand"
  # SAR = "Number of alternate observations on the reverse strand"
  key.words <- c("SRF", "SRR", "SAF", "SAR")
  
  # Check whether the vcf is indeed a freebayes vcf
  if(!all(sapply(key.words, FUN = grepl, x = vcf$INFO[1]))) {
    stop("\nVCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
         " does not appear to be a freebayes VCF, please check the data")
  }
  
  info.list <- strsplit(vcf$INFO, split = ";")
  CalculateVAF <- function(vector, key.words) {
    # Get the index of those items in vector that contain key.words
    idx <- sapply(key.words, FUN = grep, x = vector)
    
    # Get the necessary information for calculating VAF
    info <- vector[idx]
    info1 <- unlist(strsplit(info, split = "="))
    info2 <- as.integer(info1[c(2, 4, 6, 8)])
    names(info2) <- key.words
    
    return(data.frame(VAF = sum(info2[3:4]) / sum(info2), 
                      read.depth = sum(info2)))
  }
  
  list1 <- lapply(info.list, FUN = CalculateVAF, key.words = key.words)
  vafs <- do.call("rbind", list1)
  return(cbind(vcf, vafs))  
}

#' Read in the data lines of a Variant Call Format (VCF) file
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#' 
#' @param variant.caller Name of the variant caller that produces the VCF, can be either
#' \code{strelka}, \code{mutect} or \code{freebayes}. This information is needed to 
#' calculate the VAFs(variant allel frequencies). If \code{NULL}(default), then VAF and
#' read depth information will not be added to the original VCF.
#' 
#' @param name.of.VCF Name of the VCF file. If \code{NULL}(default), this
#'   function will remove all of the path up to and including the last path
#'   separator (if any) in \code{file} and file path without extensions (and the
#'   leading dot) will be used as the name of the VCF file.
#'   
#' @param tumor.col.name Optional. Only applicable to \strong{Mutect} VCF. Name
#'   of the column in \strong{Mutect} VCF which contains the tumor sample
#'   information. It \strong{must} have quotation marks. If
#'   \code{tumor.col.name} is equal to \code{NA}(default), this function will
#'   use the 10th column to calculate VAFs. See \code{\link{GetMutectVAF}} for
#'   more details.
#'   
#' @return A data frame storing data lines of the VCF file with two additional
#'   columns added which contain the VAF(variant allele frequency) and read
#'   depth information.
#'   
#' @keywords internal
ReadVCF <- 
  function(file, variant.caller = NULL, name.of.VCF = NULL, tumor.col.name = NA) {
    df1 <- df <- MakeDataFrameFromVCF(file)$df
    df1$VAF <- NA
    df1$read.depth <- NA
    
    if (is.null(variant.caller)) {
      return(StandardChromName(df1))
    } 
    
    # Check whether the variant caller is supported by ICAMS
    if (!variant.caller %in% c("strelka", "mutect", "freebayes")) {
      stop(paste0("\nVariant caller", variant.caller, "is not supported by",
                  " ICAMS, please specify either ", dQuote("strelka"), ", ", 
                  dQuote("mutect"), " or ", dQuote("freebayes")))
    }
    
    # Get the name of VCF
    if (is.null(name.of.VCF)) {
      vcf.name <- tools::file_path_sans_ext(basename(file))
    } else {
      vcf.name <- name.of.VCF
    }
    
    if (variant.caller == "strelka") {
      # Check whether the input VCF is a Strelka VCF
      if (!("TUMOR" %in% names(df)) ||
          !("FORMAT" %in% names(df))) {
        stop("\nVCF ", dQuote(vcf.name),
             " does not appear to be a Strelka VCF, column names are \n",
             paste(colnames(df), collapse=" "))
      }
      
      # Check for any SBS in df and only calcuate VAF for those SBS variants
      SBS.idx <- which(nchar(df$REF) == 1 & nchar(df$ALT) == 1)
      if (length(SBS.idx) == 0) {
        return(StandardChromName(df1))
      } else {
        SBS.df <- df[SBS.idx, ]
        SBS.df1 <- GetStrelkaVAF(SBS.df, vcf.name)
        df1[SBS.idx, ]$VAF <- SBS.df1$VAF
        df1[SBS.idx, ]$read.depth <- SBS.df1$read.depth
        return(StandardChromName(df1))
      }
    }  
    
    if (variant.caller == "mutect") {
      df2 <- GetMutectVAF(df, vcf.name, tumor.col.name)
      return(StandardChromName(df2))
    }
    
    if (variant.caller == "freebayes") {
      # Check for any SBS in df and only calcuate VAF for those SBS variants
      SBS.idx <- which(nchar(df$REF) == 1 & nchar(df$ALT) == 1)
      if (length(SBS.idx) == 0) {
        return(StandardChromName(df1))
      } else {
        SBS.df <- df[SBS.idx, ]
        SBS.df1 <- GetFreebayesVAF(SBS.df, vcf.name)
        df1[SBS.idx, ]$VAF <- SBS.df1$VAF
        df1[SBS.idx, ]$read.depth <- SBS.df1$read.depth
        return(StandardChromName(df1))
      }
    }
  }

#' Read VCF files
#'
#' @param files Character vector of file paths to the VCF files.
#' 
#' @param variant.caller Name of the variant caller that produces \strong{all}
#'   the VCFs specified by \code{files}, can be either \code{strelka},
#'   \code{mutect} or \code{freebayes}. This information is needed to calculate
#'   the VAFs(variant allel frequencies). If \code{NULL}(default), then VAF and
#'   read depth information will not be added to the original VCFs.
#'
#' @param names.of.VCFs Character vector of names of the VCF files. The order
#'   of names in \code{names.of.VCFs} should match the order of VCF file paths
#'   in \code{files}. If \code{NULL}(default), this function will remove all of
#'   the path up to and including the last path separator (if any) and file
#'   paths without extensions (and the leading dot) will be used as the names of
#'   the VCF files.
#'
#' @param tumor.col.names Optional. Only applicable to \strong{Mutect} VCFs.
#'   Character vector of column names in \strong{Mutect} VCFs which contain the
#'   tumor sample information. The order of names in \code{tumor.col.names}
#'   should match the order of \strong{Mutect} VCFs specified in \code{files}.
#'   If \code{tumor.col.names} is equal to \code{NA}(default), this function
#'   will use the 10th column in all the \strong{Mutect} VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'   
#' @return A list of data frames storing data lines of the VCF files with two
#'   additional columns added which contain the VAF(variant allele frequency)
#'   and read depth information.
#'
#' @keywords internal
ReadVCFs <- function(files, variant.caller = NULL, names.of.VCFs = NULL, 
                     tumor.col.names = NA) {
  if (is.null(names.of.VCFs)) {
    vcfs.names <- tools::file_path_sans_ext(basename(files))
  } else {
    # Check whether the number of VCFs match the number of names
    # in names.of.VCFs
    CheckNamesOfVCFs(files, names.of.VCFs)
    vcfs.names <- names.of.VCFs
  }
  
  num.of.files <- length(files)
  if (all(is.na(tumor.col.names))) {
    tumor.col.names <- rep(NA, num.of.files)
  }
  
  ReadVCF1 <- function(idx, files, variant.caller, vector1, vector2) {
    ReadVCF(file = files[idx], variant.caller = variant.caller,
            name.of.VCF = vector1[idx], tumor.col.name = vector2[idx])
  }
  
  vcfs <- lapply(1:num.of.files, FUN = ReadVCF1, files = files, 
                 variant.caller = variant.caller, 
                 vector1 = vcfs.names, vector2 = tumor.col.names)
  names(vcfs) <- vcfs.names
  return(vcfs)
}

#' @keywords internal
CheckAndReturnSplitOneMutectVCF <- 
  function(SBS.df, DBS.df, ID.df, other.subs.df, multiple.alt.df) {
    if (nrow(other.subs.df) == 0) {
      if (nrow(multiple.alt.df) == 0) {
        return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df))
      } else {
        return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df, 
                    multiple.alt = multiple.alt.df))
      }
    } else {
      if (nrow(multiple.alt.df) == 0) {
        return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
                    other.subs = other.subs.df))
      } else {
        return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
                    other.subs = other.subs.df, 
                    multiple.alt = multiple.alt.df))
      }
    }
  }

#' @title Split a mutect2 VCF into SBS, DBS, and ID VCFs, plus a list of other mutations
#'
#' @param vcf.df An in-memory data.frame representing a Mutect VCF, including
#'  VAFs, which are added by \code{\link{ReadMutectVCF}}.
#'  
#' @param name.of.VCF Name of the VCF file.
#'
#' @return A list with 3 in-memory VCFs and two left-over
#' VCF-like data frames with rows that were not incorporated
#' into the first 3 VCFs, as follows:
#'
#'  * \code{SBS}: VCF with only single base substitutions.
#'
#'  * \code{DBS}: VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  * \code{ID}: VCF with only small insertions and deletions.
#'
#'  * \code{other.subs}: \strong{Only appearing when} there is VCF like
#'  data.frame with rows for coordinate substitutions involving 3 or more
#'  nucleotides (e.g. ACT > TGA or AACT > GGTA) and rows for complex indels.
#'
#'  * \code{multiple.alt}: \strong{Only appearing when} there is VCF like
#'  data.frame with rows for variants with multiple alternative alleles, for
#'  example ACA mutated to both AGA and ATA at the same position.
#'  @md
#'
#' @keywords internal
SplitOneMutectVCF <- function(vcf.df, name.of.VCF = NULL) {
  # Mutect VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)
  multiple.alt.df <- vcf.df[multiple.alt, ]
  
  if (length(multiple.alt) != 0) {
    df <- vcf.df[-multiple.alt, ]
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants with multiple alternative alleles and were ",
            "discarded. See element multiple.alt in the return value for more ",
            "details.")
  } else {
    df <- vcf.df
  }

  SBS.df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]

  DBS.df <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]

  other.df <- df[nchar(df$REF) > 2 & nchar(df$ALT) == nchar(df$REF), ]
  
  if (nrow(other.df) > 0) {
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants involving three or more nucleotides and were ",
            "discarded. See element other.subs in the return value for ", 
            "more details.")
  }
    
  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]
  complex.indels.to.remove <- 
    which(substr(ID.df$REF, 1, 1) != substr(ID.df$ALT, 1, 1))
  complex.indels <- ID.df[complex.indels.to.remove, ]
  if (length(complex.indels.to.remove) > 0) {
    ID.df <- ID.df[-complex.indels.to.remove, ]
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has complex indels and were discarded. See element other.subs ",
            "in the return value for more details.")
  }
  
  other.df2 <- rbind(other.df, complex.indels)
  
  CheckAndReturnSplitOneMutectVCF(SBS.df, DBS.df, ID.df, 
                                  other.df2, multiple.alt.df)
}

#' @keywords internal
CheckAndReturnSplitListOfMutectVCFs <-
  function(SBS.list, DBS.list, ID.list, other.subs.list, multiple.alt.list,
           not.analyzed.list) {
    # Remove NULL elements from the list
    other.subs.list2 <- Filter(Negate(is.null), other.subs.list)
    multiple.alt.list2 <- Filter(Negate(is.null), multiple.alt.list)
    not.analyzed.list2 <- Filter(Negate(is.null), not.analyzed.list)
    
    if (length(other.subs.list2) == 0) {
      if (length(multiple.alt.list2) == 0) {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list))
        } else {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      not.analyzed = not.analyzed.list2))
        }
      } else {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      multiple.alt = multiple.alt.list2))
        } else {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      multiple.alt = multiple.alt.list2,
                      not.analyzed = not.analyzed.list2))
        }
      }
    } else {
      if (length(multiple.alt.list2) == 0) {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      other.subs = other.subs.list2))
        } else {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      other.subs = other.subs.list2,
                      not.analyzed = not.analyzed.list2))
        }
      } else {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      other.subs = other.subs.list2,
                      multiple.alt = multiple.alt.list2))
        } else {
          return(list(SBS = SBS.list, DBS = DBS.list, ID = ID.list,
                      other.subs = other.subs.list2,
                      multiple.alt = multiple.alt.list2,
                      not.analyzed = not.analyzed.list2))
        }
      }
    }
  }

#' Split each Mutect VCF into SBS, DBS, and ID VCFs (plus two
#' VCF-like data frame with left-over rows).
#'
#' @param list.of.vcfs List of VCFs as in-memory data.frames.
#' 
#' @inheritParams ReadAndSplitMutectVCFs
#'
#' @inheritSection ReadAndSplitMutectVCFs Value
#'
#' @keywords internal
SplitListOfMutectVCFs <- 
  function(list.of.vcfs,
           suppress.discarded.variants.warnings = TRUE) {
    names.of.VCFs <- names(list.of.vcfs)
    list.of.vcfs.df <- lapply(list.of.vcfs, function(f1) f1$df)
    list.of.discarded.variants <- 
      lapply(list.of.vcfs, function(f2) f2$discarded.variants)
    
    GetSplitMutectVCFs <- function(idx, list.of.vcfs) {
      split.vcfs <- SplitOneMutectVCF(list.of.vcfs[[idx]], 
                                      name.of.VCF = names(list.of.vcfs)[idx])
      return(split.vcfs)
    }
    num.of.vcfs <- length(list.of.vcfs.df)
    if (suppress.discarded.variants.warnings == TRUE) {
      v1 <- suppressWarnings(lapply(1:num.of.vcfs, GetSplitMutectVCFs,
                                    list.of.vcfs = list.of.vcfs.df))
    } else {
      v1 <- lapply(1:num.of.vcfs, GetSplitMutectVCFs,
                   list.of.vcfs = list.of.vcfs.df)
    }
    names(v1) <- names.of.VCFs
    SBS <- lapply(v1, function(x) x$SBS)
    DBS <- lapply(v1, function(x) x$DBS)
    ID  <- lapply(v1, function(x) x$ID)
    other.subs <- lapply(v1, function(x) x$other.subs)
    multiple.alt <- lapply(v1, function(x) x$multiple.alt)
    
    CheckAndReturnSplitListOfMutectVCFs(SBS, DBS, ID, other.subs, 
                                        multiple.alt, 
                                        list.of.discarded.variants)
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
    GenomicRanges::GRanges(chr.names,
                           IRanges::IRanges(start = df$POS - seq.context.width, # 10,
                                            end = df$POS + seq.context.width) # 10
    )

  # Extract sequence context from the reference genome
  df$extracted.seq <- BSgenome::getSeq(ref.genome, Ranges, as.character = TRUE)

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
#' @return A minimal VCF with only the columns \code{CHROM}, \code{POS},
#'   \code{ID}, \code{REF}, \code{ALT}, \code{VAF}, \code{read.depth}.
#'   
#' @keywords internal
MakeVCFDBSdf <- function(DBS.range.df, SBS.vcf.dt) {
  tmpvcf <- SBS.vcf.dt[ , c("CHROM", "POS", "REF", "ALT", "TUMOR", "VAF")]
  DBS.range.dt <- as.data.table(DBS.range.df)
  tmp1 <- merge(DBS.range.dt, tmpvcf,
                by.x = c("CHROM", "LOW"),
                by.y = c("CHROM", "POS"))
  tmp2 <- merge(tmp1, tmpvcf,
                by.x = c("CHROM", "HIGH"),
                by.y = c("CHROM", "POS"))
  # Calculate the read depth for tier1
  tmp2[, DP.x := as.integer(sapply(strsplit(TUMOR.x, ":"), "[", 1))]
  tmp2[, DP.y := as.integer(sapply(strsplit(TUMOR.y, ":"), "[", 1))]
  tmp2[, read.depth := pmin(DP.x, DP.y)]
  
  tmp2[, VAF := rowMeans(cbind(VAF.x, VAF.y))]
  tmp2[, POS := LOW]
  tmp2[, ID := "From merged SBSs"]
  tmp2[, REF := paste0(REF.x, REF.y)]
  tmp2[, ALT := paste0(ALT.x, ALT.y)]
  return(as.data.frame(tmp2[, c("CHROM", "POS", "ID", "REF", "ALT", 
                                "VAF", "read.depth")]))
}

#' @keywords internal
CheckAndReturnSplitStrelkaSBSVCF <- 
  function(SBS.df, DBS.df, ThreePlus.df, multiple.alt.df) {
    if (nrow(ThreePlus.df) == 0) {
      if (nrow(multiple.alt.df) == 0) {
        return(list(SBS.vcf = SBS.df, DBS.vcf = DBS.df))
      } else {
        return(list(SBS.vcf = SBS.df, DBS.vcf = DBS.df, 
                    multiple.alt = multiple.alt.df))
      }
    } else {
      if (nrow(multiple.alt.df) == 0) {
        return(list(SBS.vcf = SBS.df, DBS.vcf = DBS.df,
                    ThreePlus = ThreePlus.df))
      } else {
        return(list(SBS.vcf = SBS.df, DBS.vcf = DBS.df,
                    ThreePlus = ThreePlus.df, 
                    multiple.alt = multiple.alt.df))
      }
    }
  }

#' Split an in-memory Strelka VCF into SBS, DBS, and variants involving
#' > 2 consecutive bases
#'
#' SBSs are single base substitutions,
#' e.g. C>T, A>G,....  DBSs are double base substitutions,
#' e.g. CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT>TCTA, ...
#'
#' @param vcf.df An in-memory data frame containing a Strelka VCF file contents.
#'
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02.
#' 
#' @param name.of.VCF Name of the VCF file.
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
#' 
#' \enumerate{
#'    \item \code{SBS.vcf}: Data frame of pure SBS mutations -- no DBS or 3+BS
#'    mutations.
#'
#'    \item \code{DBS.vcf}: Data frame of pure DBS mutations -- no SBS or 3+BS
#'    mutations.
#'
#'    \item \code{ThreePlus}: Data table with the key CHROM, LOW.POS, HIGH.POS
#'    and additional information (reference sequence, alternative sequence,
#'    context, etc.) Additional information not fully implemented at this point
#'    because of limited immediate biological interest.
#'
#'    \item \code{multiple.alt} Rows with multiple alternate alleles (removed
#'    from \code{SBS.vcf} etc.)
#'    
#'    }
#'
#' @keywords internal
SplitStrelkaSBSVCF <- 
  function(vcf.df, max.vaf.diff = 0.02, name.of.VCF = NULL) {
  stopifnot("data.frame" %in% class(vcf.df))
  
  # Strelka SBS VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)
  multiple.alt.df <- vcf.df[multiple.alt, ]
  
  if (length(multiple.alt) != 0) {
    vcf.df <- vcf.df[-multiple.alt, ]
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants with multiple alternative alleles and were ",
            "discarded. See element multiple.alt in the return value for more ",
            "details.")
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
  # pairs CC > TT and CA > TG (more below).

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
  # If VAFs are not similar, the adjacent SBSs are likely to be "merely"
  # asynchronous single base mutations, opposed to a simultaneous doublet mutation.
  non.SBS <- dt2[abs(VAF.x - VAF.y) <= max.vaf.diff]
  # TODO: if (any(is.na(VAF.x)) || any(is.na(VAF.y)))
  # If VAF.x or VAF.y is NA the row will not go into non.SBS.
  rm(dt2)

  if (nrow(non.SBS) == 0) {
    # There are no non.SBS mutations in the input.
    # Everything in vcf.df is an SBS. We are finished.
    empty <- vcf.df[-(1:nrow(vcf.df)), ]
    return(list(SBS.vcf = vcf.df, DBS.vcf = empty,
                ThreePlus =
                  data.table(CHROM = character(),
                             LOW.POS = numeric(),
                             HIGH.POS = numeric()),
                multiple.alt = empty))
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
    GenomicRanges::GRanges(non.SBS$CHROM, 
                           IRanges::IRanges(start = non.SBS$LOW, end = non.SBS$HIGH))
  rranges <- GenomicRanges::reduce(ranges) # Merge overlapping ranges
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
  
  CheckAndReturnSplitStrelkaSBSVCF(out.SBS.df, DBS.vcf.df, 
                                   other.ranges, multiple.alt.df)
}

#' @keywords internal
CheckAndReturnSplitListOfStrelkaSBSVCFs <-
  function(SBS.list, DBS.list, ThreePlus.list, multiple.alt.list,
           not.analyzed.list) {
    # Remove NULL elements from the list
    ThreePlus.list2 <- Filter(Negate(is.null), ThreePlus.list)
    multiple.alt.list2 <- Filter(Negate(is.null), multiple.alt.list)
    not.analyzed.list2 <- Filter(Negate(is.null), not.analyzed.list)
    
    if (length(ThreePlus.list2) == 0) {
      if (length(multiple.alt.list2) == 0) {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list))
        } else {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      not.analyzed = not.analyzed.list2))
        }
      } else {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      multiple.alt = multiple.alt.list2))
        } else {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      multiple.alt = multiple.alt.list2,
                      not.analyzed = not.analyzed.list2))
        }
      }
    } else {
      if (length(multiple.alt.list2) == 0) {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      ThreePlus = ThreePlus.list2))
        } else {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      ThreePlus = ThreePlus.list2,
                      not.analyzed = not.analyzed.list2))
        }
      } else {
        if (length(not.analyzed.list2) == 0) {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      ThreePlus = ThreePlus.list2,
                      multiple.alt = multiple.alt.list2))
        } else {
          return(list(SBS.vcfs = SBS.list, DBS.vcfs = DBS.list,
                      ThreePlus = ThreePlus.list2,
                      multiple.alt = multiple.alt.list2,
                      not.analyzed = not.analyzed.list2))
        }
      }
    }
  }

#' Split a list of in-memory Strelka SBS VCF into SBS, DBS, and variants involving
#' > 2 consecutive bases
#'
#' SBSs are single base substitutions,
#' e.g. C>T, A<G,....  DBSs are double base substitutions,
#' e.g. CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT>TCTA, ...
#'
#' @param list.of.vcfs A list of in-memory data frames containing Strelka SBS
#'   VCF file contents.
#'   
#' @inheritParams ReadAndSplitStrelkaSBSVCFs
#'   
#' @inheritSection ReadAndSplitStrelkaSBSVCFs Value
#'
#' @keywords internal
SplitListOfStrelkaSBSVCFs <- 
  function(list.of.vcfs, suppress.discarded.variants.warnings = TRUE) {
    names.of.VCFs <- names(list.of.vcfs)
    list.of.vcfs.df <- lapply(list.of.vcfs, function(f1) f1$df)
    list.of.discarded.variants <- 
      lapply(list.of.vcfs, function(f2) f2$discarded.variants)
    
    GetSplitStrelkaSBSVCFs <- function(idx, list.of.vcfs) {
      split.vcfs <- SplitStrelkaSBSVCF(list.of.vcfs[[idx]], 
                                       name.of.VCF = names(list.of.vcfs)[idx])
      return(split.vcfs)
    }
    num.of.vcfs <- length(list.of.vcfs.df)
    if (suppress.discarded.variants.warnings == TRUE) {
      split.vcfs <- 
        suppressWarnings(lapply(1:num.of.vcfs, GetSplitStrelkaSBSVCFs,
                                list.of.vcfs = list.of.vcfs.df))
    } else {
      split.vcfs <- lapply(1:num.of.vcfs, GetSplitStrelkaSBSVCFs,
                           list.of.vcfs = list.of.vcfs.df)
    }
    names(split.vcfs) <- names.of.VCFs
    SBS.vcfs   <- lapply(split.vcfs, function(x) x$SBS.vcf)
    DBS.vcfs   <- lapply(split.vcfs, function(x) x$DBS.vcf)
    ThreePlus  <- lapply(split.vcfs, function(x) x$ThreePlus)
    mult.alt   <- lapply(split.vcfs, function(x) x$multiple.alt)
    CheckAndReturnSplitListOfStrelkaSBSVCFs(SBS.vcfs, DBS.vcfs, ThreePlus,
                                            mult.alt, list.of.discarded.variants)
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
#' @inheritParams ReadMutectVCFs
#'    
#' @inheritSection ReadMutectVCFs Value
#'
#' @keywords internal
ReadStrelkaSBSVCFs <- function(files, names.of.VCFs = NULL,
                               suppress.discarded.variants.warnings = TRUE) {
  vcfs <- 
    lapply(files, FUN = ReadStrelkaSBSVCF, name.of.VCF = names.of.VCFs,
           suppress.discarded.variants.warnings = 
             suppress.discarded.variants.warnings)
  if (is.null(names.of.VCFs)) {
    names(vcfs) <- tools::file_path_sans_ext(basename(files))
  } else {
    # Check whether the number of VCFs match the number of names
    # in names.of.VCFs
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
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'   
#' @section Value: A list of \strong{lists}. Each list has a first element \code{df}
#'   which is a data frame that stores data lines of a VCF with additional
#'   columns \code{VAF} (variant allele frequency) and \code{read.depth} added.
#'   A second element \code{discarded.variants} \strong{only} appears if there
#'   are variants that are excluded from the analysis.
#'   
#' @keywords internal
ReadMutectVCFs <- 
  function(files, names.of.VCFs = NULL, tumor.col.names = NA,
           suppress.discarded.variants.warnings = TRUE) {
  if (is.null(names.of.VCFs)) {
    vcfs.names <- tools::file_path_sans_ext(basename(files))
  } else {
    # Check whether the number of VCFs match the number of names
    # in names.of.VCFs
    CheckNamesOfVCFs(files, names.of.VCFs)
    vcfs.names <- names.of.VCFs
  }
  num.of.files <- length(files)
  if (all(is.na(tumor.col.names))) {
    tumor.col.names <- rep(NA, num.of.files)
  }
  
  GetMutectVCFs <- function(idx, files, names.of.VCFs, tumor.col.names,
                            suppress.discarded.variants.warnings) {
    ReadMutectVCF(file = files[idx], name.of.VCF = names.of.VCFs[idx],
                  tumor.col.name = tumor.col.names[idx],
                  suppress.discarded.variants.warnings)
  }
  
  vcfs <- lapply(1:num.of.files, FUN = GetMutectVCFs, 
                 files, vcfs.names, tumor.col.names, 
                 suppress.discarded.variants.warnings)
  names(vcfs) <- vcfs.names
  return(vcfs)
}

#' Add sequence context and transcript information to an in-memory SBS VCF
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
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateSBSVCF <- function(SBS.vcf, ref.genome, trans.ranges = NULL) {
  SBS.vcf <- AddSeqContext(SBS.vcf, ref.genome = ref.genome)
  CheckSeqContextInVCF(SBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    SBS.vcf <- AddTranscript(SBS.vcf, trans.ranges)
  }
  return(as.data.table(SBS.vcf))
}


#' Add SBS mutation class to an annotated SBS VCF
#' 
#' @param vcf An in-memory VCF file annotated with sequence context and
#'   transcript information by function \code{\link{AnnotateSBSVCF}}. It must
#'   *not* contain indels and must *not* contain DBS (double base
#'   substitutions), or triplet base substitutions etc., even if encoded as
#'   neighboring SBS.
#'
#' @return The original \code{vcf} with three additional columns
#'   \code{SBS96.class}, \code{SBS192.class} and \code{SBS1536.class} added.
#'   
#' @keywords internal
AddSBSClass <- function(vcf) {
  col.names <- colnames(vcf)
  vcf$SBS1536.class <- paste0(substr(vcf$seq.21bases, 9, 13), vcf$ALT)
  vcf$SBS1536.class <- PyrPenta(vcf$SBS1536.class)
  vcf$SBS96.class <- paste0(substr(vcf$SBS1536.class, 2, 4),
                            substr(vcf$SBS1536.class, 6, 6))
  vcf$SBS192.class <- NA
  idx <- which(!is.na(vcf$trans.strand) & (vcf$bothstrand == FALSE))
  vcf$SBS192.class[idx] <- vcf$SBS96.class[idx]
  idx1 <- which(vcf$trans.strand == "-" & (vcf$bothstrand == FALSE))
  vcf$SBS192.class[idx1] <- RevcSBS96(vcf$SBS192.class[idx1])
  
  # Reorder the columns of vcf
  new.col.names <- c(col.names, "SBS96.class", "SBS192.class", "SBS1536.class")
  data.table::setDT(vcf)
  setcolorder(vcf, new.col.names)
  return(vcf)
}

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
      stop("In sample ", sample.id, ", the number of SBS", nrow(mat), 
           " variants in the annotated VCF is not the same as the total ", 
           "counts in mutation matrix.")
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
      stop("In sample ", sample.id, ", the number of SBS", nrow(mat), 
           " variants in the annotated VCF is not the same as the total ", 
           "counts in mutation matrix.")
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
  function(vcf, discarded.variants, mat96, mat1536, mat192 = NULL, 
           return.annotated.vcf = FALSE, sample.id = "counts") {
    
    if (nrow(discarded.variants) == 0) {
      if (is.null(mat192)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catSBS96 = mat96, catSBS1536 = mat1536))
        } else {
          vcf.SBS.class <- 
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, sample.id)
          return(list(catSBS96 = mat96, catSBS1536 = mat1536,
                      annotated.vcf = vcf.SBS.class))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(catSBS96 = mat96, catSBS192 = mat192, 
                      catSBS1536 = mat1536))
        } else {
          vcf.SBS.class <- 
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(catSBS96 = mat96, catSBS192 = mat192, catSBS1536 = mat1536,
                      annotated.vcf = vcf.SBS.class))
        }
      }
    } else {
      if (is.null(mat192)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catSBS96 = mat96, catSBS1536 = mat1536,
                      discarded.variants = discarded.variants))
        } else {
          vcf.SBS.class <- 
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, sample.id)
          return(list(catSBS96 = mat96, catSBS1536 = mat1536,
                      annotated.vcf = vcf.SBS.class,
                      discarded.variants = discarded.variants))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(catSBS96 = mat96, catSBS192 = mat192, 
                      catSBS1536 = mat1536, 
                      discarded.variants = discarded.variants))
        } else {
          vcf.SBS.class <- 
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
          return(list(catSBS96 = mat96, catSBS192 = mat192, catSBS1536 = mat1536,
                      annotated.vcf = vcf.SBS.class,
                      discarded.variants = discarded.variants))
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
CreateOneColSBSMatrix <- function(vcf, sample.id = "count",
                                  return.annotated.vcf = FALSE) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude DBSs
  # before calling the function. This function does not detect DBSs.

  if (0 == nrow(vcf)) {
    # Create 1-column matrix with all values being 0 and the correct row and
    # column labels.
    catSBS96 <- 
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS96), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$SBS96, sample.id))
    catSBS192 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS192), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$SBS192, sample.id))
    catSBS1536 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS1536), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$SBS1536, sample.id))
    
    if (return.annotated.vcf == FALSE) {
      return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                  catSBS1536 = catSBS1536))
    } else {
      return(list(catSBS96 = catSBS96, catSBS192 = catSBS192,
                  catSBS1536 = catSBS1536, annotated.vcf = vcf))
    }
  }

  stopifnot(nchar(vcf$ALT) == 1)
  stopifnot(nchar(vcf$REF) == 1)
  stopifnot(vcf$ALT != vcf$REF)
  mismatches <- which(vcf$REF != substr(vcf$seq.21bases, 11, 11))
  if (length(mismatches) != 0) {
    stop("\nSample ", sample.id, 
         ":\nThe reference base in ref.genome does not match the ", 
         "reference base in ", length(mismatches),
         " rows in the VCF file.\n",
         "Please check the ref.genome argument.")
  }
  
  discarded.variants <- vcf[0]
  # Delete the rows of SBS if the pentanucleotide context contains "N"
  idx <- grep("N", substr(vcf$seq.21bases, 9, 13))
    if (!length(idx) == 0) {
      discarded.variants <- rbind(discarded.variants, vcf[idx, ])
      vcf <- vcf[-idx, ]
      warning(
        'Variants in the SBS vcf whose pentanucleotide context contains "N" ',
        'have been deleted so as not to conflict with downstream processing. ',
        'See discarded.variants in the return value for more details.')
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
    CheckAndReturnSBSMatrix(vcf0, discarded.variants, mat96, 
                            mat1536, return.annotated.vcf, sample.id)
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
  
  # If vcf3 has empty rows, we will return 1-column SBS192 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    mat192 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS192), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$SBS192, sample.id))
    CheckAndReturnSBSMatrix(vcf0, discarded.variants, mat96, mat1536, mat192,
                            return.annotated.vcf, sample.id)
  }

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
  
  CheckAndReturnSBSMatrix(vcf0, discarded.variants, mat96, mat1536, mat192,
                          return.annotated.vcf, sample.id)
}

#' Add sequence context and transcript information to an in-memory DBS VCF
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
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' DBS.vcf <- list.of.vcfs$DBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateDBSVCF <- function(DBS.vcf, ref.genome, trans.ranges = NULL) {
  DBS.vcf <- AddSeqContext(DBS.vcf, ref.genome = ref.genome)
  CheckSeqContextInVCF(DBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    DBS.vcf <- AddTranscript(DBS.vcf, trans.ranges)
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
      stop("In sample ", sample.id, ", the number of DBS", nrow(mat), 
           " variants in the annotated VCF is not the same as the total ", 
           "counts in mutation matrix.")
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
      stop("In sample ", sample.id, ", the number of DBS", nrow(mat), 
           " variants in the annotated VCF is not the same as the total ", 
           "counts in mutation matrix.")
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
  function(vcf, discarded.variants, mat78, mat136, mat144 = NULL, 
           return.annotated.vcf = FALSE, sample.id = "counts") {
    
    if (nrow(discarded.variants) == 0) {
      if (is.null(mat144)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS136 = mat136))
        } else {
          vcf.DBS.class <- 
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, sample.id)
          return(list(catDBS78 = mat78, catDBS136 = mat136,
                      annotated.vcf = vcf.DBS.class))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS144 = mat144, 
                      catDBS136 = mat136))
        } else {
          vcf.DBS.class <- 
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(catDBS78 = mat78, catDBS144 = mat144, catDBS136 = mat136,
                      annotated.vcf = vcf.DBS.class))
        }
      }
    } else {
      if (is.null(mat144)) {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS136 = mat136,
                      discarded.variants = discarded.variants))
        } else {
          vcf.DBS.class <- 
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, sample.id)
          return(list(catDBS78 = mat78, catDBS136 = mat136,
                      annotated.vcf = vcf.DBS.class,
                      discarded.variants = discarded.variants))
        }
      } else {
        if (return.annotated.vcf == FALSE) {
          return(list(catDBS78 = mat78, catDBS144 = mat144, 
                      catDBS136 = mat136, 
                      discarded.variants = discarded.variants))
        } else {
          vcf.DBS.class <- 
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
          return(list(catDBS78 = mat78, catDBS144 = mat144, catDBS136 = mat136,
                      annotated.vcf = vcf.DBS.class,
                      discarded.variants = discarded.variants))
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
CreateOneColDBSMatrix <- function(vcf, sample.id = "count",
                                  return.annotated.vcf = FALSE) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude SBSs
  # before calling the function. This function does not detect SBSs.

  if (0 == nrow(vcf)) {
    # Create 1-column matrix with all values being 0 and the correct row labels.
    catDBS78 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS78), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$DBS78, sample.id))
    catDBS136 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS136), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$DBS136, sample.id))
    catDBS144 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS144), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$DBS144), sample.id)
    if (return.annotated.vcf == FALSE) {
      return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                  catDBS144 = catDBS144))
    } else {
      return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                  catDBS144 = catDBS144, annotated.vcf = vcf))
    }
  }

  stopifnot(nchar(vcf$ALT) == 2)
  stopifnot(nchar(vcf$REF) == 2)
  
  discarded.variants <- vcf[0]
  # Delete the rows of DBS if the tetranucleotide context contains "N"
  idx <- grep("N", substr(vcf$seq.21bases, 10, 13))
  if (!length(idx) == 0) {
    discarded.variants <- rbind(discarded.variants, vcf[idx, ])
    vcf <- vcf[-idx, ]
    warning(
      'Variants in the DBS vcf whose tetranucleotide context contains "N" ',
      'have been deleted so as not to conflict with downstream processing. ',
      'See discarded.variants in the return value for more details.')
  }

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
    CheckAndReturnDBSMatrix(vcf, discarded.variants, DBS.mat.78, DBS.mat.136, 
                            return.annotated.vcf, sample.id)
  }
  
  # There may be some mutations in vcf which fall on transcripts on both
  # strands. We do not consider those mutations when generating the 144 catalog.
  vcf2 <- vcf[bothstrand == FALSE, ]

  # One DBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 144 catalog, we only need to count these mutations once.
  vcf3 <- vcf2[, .(REF = REF[1], trans.strand = trans.strand[1]),
               by = .(CHROM, ALT, POS)]
  
  # If vcf3 has empty rows, we will return 1-column DBS144 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    DBS.mat.144 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS144), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$DBS144), sample.id)
    CheckAndReturnDBSMatrix(vcf, discarded.variants, DBS.mat.78, DBS.mat.136, 
                            DBS.mat.144, return.annotated.vcf, sample.id)
  }

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

  CheckAndReturnDBSMatrix(vcf, discarded.variants, DBS.mat.78, DBS.mat.136, 
                          DBS.mat.144, return.annotated.vcf, sample.id)
}

#' Create SBS and DBS catalogs from Strelka SBS VCF files and plot them to PDF
#'
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
#' * \code{discarded.variants}: 
#' \strong{Only appearing when} there are variants that were excluded in the
#' analysis. 
#' A list of elements:
#'     + \code{SBS}: SBS variants whose pentanucleotide context contains "N".
#'     + \code{DBS}: DBS variants whose tetranucleotide context contains "N".
#'     + \code{ThreePlus}: Variants involving three or more nucleotides (e.g. ACT >
#'       TGA or AACT > GGTA).
#'     + \code{multiple.alt}: Variants with multiple alternative alleles, for
#'       example ACA mutated to both AGA and ATA at the same position.
#'     + \code{not.analyzed}: Variants discarded immediately after reading in
#'     the VCFs:
#'         - Duplicated "CHROM" and "POS" values.
#'         - Chromosome names that contain "#".
#'         - Chromosome names that contain "GL".
#'         - Chromosome names that contain "KI".
#'         - Chromosome names that contain "random".
#'         - Chromosome names that contain "Hs".
#'         - Chromosome names that contain "M".
#'   
#' * \code{annotated.vcfs}: 
#' \strong{Only appearing when} \code{return.annotated.vcfs} = TRUE.
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
StrelkaSBSVCFFilesToCatalogAndPlotToPdf <- 
  function(files, 
           ref.genome, 
           trans.ranges = NULL, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           output.file = "",
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    
    catalogs0 <- 
      StrelkaSBSVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                                  region, names.of.VCFs, 
                                  return.annotated.vcfs, 
                                  suppress.discarded.variants.warnings)
    
    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
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
#' @inheritSection StrelkaIDVCFFilesToCatalog Value
#' 
#' @inheritSection VCFsToIDCatalogs ID classification
#' 
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   catID <- 
#'     StrelkaIDVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                            region = "genome",
#'                                            output.file = 
#'                                            file.path(tempdir(), "StrelkaID"))}
#'                                                                    
StrelkaIDVCFFilesToCatalogAndPlotToPdf <- 
  function(files, 
           ref.genome, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           output.file = "",
           flag.mismatches = 0,
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    
    list <-
      StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs,
                                 flag.mismatches, return.annotated.vcfs,
                                 suppress.discarded.variants.warnings)
    
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
#' * \code{catID}: Matrix of ID (small insertion and deletion) catalog.
#' 
#' * \code{discarded.variants}: 
#' \strong{Only appearing when} there are variants that were excluded in the
#' analysis. 
#' A list of elements:
#'     + \code{SBS}: SBS variants whose pentanucleotide context contains "N".
#'     + \code{DBS}: DBS variants whose tetranucleotide context contains "N".
#'     + \code{ID}: ID variants discarded that can belong to the following
#'       categories:
#'         - Variants which have empty REF or ALT allels.
#'         - Variants whose REF allels do not match the extracted sequence from
#'          \code{ref.genome}.
#'         - Variants which cannot be categorized according to the canonical
#'         representation. See catalog.row.order$ID for the canonical
#'         representation.
#'     + \code{other.subs}: Variants involving three or more nucleotides (e.g. ACT >
#'       TGA or AACT > GGTA) and complex indels.
#'     + \code{multiple.alt}: Variants with multiple alternative alleles, for
#'       example ACA mutated to both AGA and ATA at the same position.
#'     + \code{not.analyzed}: Variants discarded immediately after reading in
#'     the VCFs:
#'         - Duplicated "CHROM" and "POS" values.
#'         - Chromosome names that contain "#".
#'         - Chromosome names that contain "GL".
#'         - Chromosome names that contain "KI".
#'         - Chromosome names that contain "random".
#'         - Chromosome names that contain "Hs".
#'         - Chromosome names that contain "M".
#'   
#' * \code{annotated.vcfs}: 
#' \strong{Only appearing when} \code{return.annotated.vcfs} = TRUE.
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
#' See
#' \url{https://github.com/steverozen/ICAMS/raw/master/data-raw/PCAWG7_indel_classification_2017_12_08.xlsx}
#' for additional information on ID (small insertion and deletion) mutation
#' classification. 
#' 
#' See the documentation for \code{\link{Canonicalize1Del}} which first handles
#' deletions in homopolymers, then handles deletions in simple repeats with
#' longer repeat units, (e.g. \code{CACACACA}, see
#' \code{\link{FindMaxRepeatDel}}), and if the deletion is not in a simple
#' repeat, looks for microhomology (see \code{\link{FindDelMH}}).
#' 
#' See the code for unexported function \code{\link{CanonicalizeID}}
#' and the functions it calls for handling of insertions.
#'
#' @section Note:
#'  SBS 192 and DBS 144 catalogs include only mutations in transcribed regions.
#'  In ID (small insertion and deletion) catalogs, deletion repeat sizes range
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
#'     MutectVCFFilesToCatalogAndPlotToPdf(file, ref.genome = "hg19", 
#'                                         trans.ranges = trans.ranges.GRCh37,
#'                                         region = "genome",
#'                                         output.file = 
#'                                         file.path(tempdir(), "Mutect"))}
MutectVCFFilesToCatalogAndPlotToPdf <- 
  function(files, 
           ref.genome, 
           trans.ranges = NULL, 
           region = "unknown", 
           names.of.VCFs = NULL, 
           tumor.col.names = NA,
           output.file = "",
           flag.mismatches = 0,
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    
    catalogs0 <-
      MutectVCFFilesToCatalog(files, ref.genome, trans.ranges, 
                              region, names.of.VCFs, tumor.col.names,
                              flag.mismatches, return.annotated.vcfs,
                              suppress.discarded.variants.warnings)
    
    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
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
    
    return(catalogs0)
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
