#' Extract the VAFs (variant allele frequencies) and read depth information from a VCF file
#'
#' @param vcf Said VCF as a data.frame.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param tumor.col.name Optional. Only applicable to \strong{Mutect} VCF. Name
#'   or index of the column in \strong{Mutect} VCF which contains the tumor
#'   sample information. It \strong{must} have quotation marks if specifying the
#'   column name. If \code{tumor.col.name} is equal to \code{NA}(default), this
#'   function will use the 10th column to calculate VAFs.
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
#' df <- MakeDataFrameFromVCF(file)
#' df1 <- GetStrelkaVAF(df)
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
RemoveRowsWithPoundSignNew <- function(df, name.of.VCF = NULL) {
  pound.chrom.idx <- which(df$CHROM == "#CHROM")
  if (length(pound.chrom.idx) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", length(pound.chrom.idx), " row out of ",
            nrow(df), " had value #CHROM in column 'CHROM' and were removed. ",
            "See discarded.variants in the return value for more details")
    df1 <- df[-pound.chrom.idx, ]
	df1.to.remove <- df[pound.chrom.idx, ]
	df1.to.remove$discarded.reason <- 'Chromosome name is "#CHROM"'
    return(list(df = df1, discarded.variants = df1.to.remove))
  } else {
    return(list(df = df))
  }
}

#' @keywords internal
RemoveRowsWithDuplicatedCHROMAndPOS <- function(df, file) {
  dups <- which(duplicated(df[, c("CHROM", "POS")]))
  if (length(dups) > 0) {
    dups2 <- which(duplicated(df[ , c("CHROM", "POS")], fromLast = TRUE))
    warning("In ", file, " ", 2 * length(dups), " rows out of ",
            nrow(df), " had duplicate CHROM and POS and were removed: ",
            dups2, " ", dups)
    df1 <- df[-c(dups, dups2), ]
    return(df1)
  } else {
    return(df)
  }
}

#' @keywords internal
RemoveRowsWithDuplicatedCHROMAndPOSNew <- function(df, name.of.VCF = NULL) {
  discarded.variants <- df[0, ]
  # Find out variants which have the same "CHROM", "POS", "REF", "ALT"
  dups <- which(duplicated(df[, c("CHROM", "POS", "REF", "ALT")]))

  if (length(dups) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", 2 * length(dups), " row out of ",
            nrow(df), " had same CHROM, POS, REF and ALT and only one copy is kept. ",
            "See discarded.variants in the return value for more details")
    df.to.remove <- df[dups, ]
    df.to.remove$discarded.reason <- "Variant with same CHROM, POS, REF and ALT as another variant"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df.to.remove)
    df1 <- df[-dups, ]
  } else {
    df1 <- df
  }

  # Find out variants which have the same "CHROM", "POS", "REF" but different "ALT"
  dups2 <- which(duplicated(df1[, c("CHROM", "POS", "REF")]))

  if (length(dups2) > 0) {
    warning("In VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " ", 2 * length(dups2), " row out of ",
            nrow(df), " had same CHROM, POS, REF but different ALT and were removed. ",
            "See discarded.variants in the return value for more details")
    dups3 <- which(duplicated(df1[, c("CHROM", "POS", "REF")], fromLast = TRUE))
    df1.to.remove <- df1[c(dups2, dups3), ]
    df1.to.remove$discarded.reason <- "Variant with same CHROM, POS, REF but different ALT"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df1.to.remove)
    df2 <- df1[-c(dups2, dups3), ]
  } else {
    df2 <- df1
  }

  if (nrow(discarded.variants) > 0) {
    return(list(df = df2, discarded.variants = discarded.variants))
  } else {
    return(list(df = df2))
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
    warning('A column named "strand" in the VCF ',
            'was renamed to "strand_old" ',
            'avoid conflict with a newly added column named "strand".')
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
    warning('A column named "VAF" in the VCF ',
            'was renamed to "VAF_old" ',
            'avoid conflict with a newly added column named "VAF".')
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
#' @inheritSection ReadMutectVCF Value
#'
#' @keywords internal
ReadStrelkaSBSVCF <- function(file, name.of.VCF = NULL) {

  df <- MakeDataFrameFromVCF(file)

  if (is.null(name.of.VCF)) {
    vcf.name <- tools::file_path_sans_ext(basename(file))
  } else {
    vcf.name <- name.of.VCF
  }

  if (nrow(df) == 0) {
    return(df)
  } else {
    df1 <- GetStrelkaVAF(vcf = df, name.of.VCF = vcf.name)
    return(df1)
  }
}

#' Read a VCF file into a data frame with minimal processing.
#'
#' @details Header lines beginning "##" are removed, and column
#'   "#CHROM" is renamed to "CHROM". Other column names are
#'   unchanged. Columns "#CHROM", "POS", "REF", and "ALT" must
#'   be in the input.
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @export
#'
#' @return A data frame storing mutation records of a VCF file.
#'
#' @examples
#' file <- c(system.file("extdata/Strelka-SBS-vcf",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' df <- SimpleReadVCF(file)
SimpleReadVCF <- function(file) {
  return(MakeDataFrameFromVCF(file = file))
}

#' Read in the data lines of a Variant Call Format (VCF) file
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @return A data frame storing mutation records of a VCF file.
#'
#' @keywords internal
MakeDataFrameFromVCF <- function(file) { # name.of.VCF = NULL) {

  # Suppress the warning when the VCF is totally empty
  tryCatch({
    df1 <-
      suppressWarnings(data.table::fread(file, na.strings = "",
                                         skip = "#CHROM", fill = TRUE))

    if (nrow(df1) == 0) {
      return(df1)
    }

    required.col.names <- c("#CHROM", "POS", "REF", "ALT")
    col.names.exist <- required.col.names %in% colnames(df1)
    col.names.not.available <- required.col.names[!col.names.exist]

    if (!all(col.names.exist)) {
      stop("some columns required in VCF are not available ",
           paste(col.names.not.available, collapse = " "))
    }

    },
    error = function(err.info) {
      if (!is.null(err.info$message)) {
      stop(file, " does not appear to be a VCF file.\nDetails: ",
           err.info$message)
      }
    })

  # Extract the names of columns in the VCF file
  names <- c("CHROM", colnames(df1)[-1])
  colnames(df1) <- names

  df1$CHROM <- as.character(df1$CHROM)

  return(df1)
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
#' @return A data frame storing data lines of the VCF file.
#'
#' @inheritSection VCFsToIDCatalogs Note
#'
#' @keywords internal
ReadStrelkaIDVCF <- function(file, name.of.VCF = NULL) {
  df1 <- MakeDataFrameFromVCF(file)

  # Get the name of VCF
  if (is.null(name.of.VCF)) {
    vcf.name <- tools::file_path_sans_ext(basename(file))
  } else {
    vcf.name <- name.of.VCF
  }

  if (nrow(df1) == 0) {
    return(df1)
  }

  # Check whether the input VCF is a Strelka ID VCF
  if (!("TUMOR" %in% names(df1)) ||
      !("FORMAT" %in% names(df1))) {
    stop("\nVCF ", dQuote(vcf.name),
         " does not appear to be a Strelka VCF, column names are \n",
         paste(colnames(df1), collapse=" "))
  }
  control <- unique(df1[["FORMAT"]])
  stopifnot(length(control) == 1)
  colnames <- unlist(strsplit(control, split=":", fixed=TRUE))
  each.base.col <- c("AU", "CU", "GU", "TU")
  if (all(each.base.col %in% colnames)) {
    stop("\nVCF ", dQuote(vcf.name),
         " does not appear to be a Strelka ID VCF, ",
         "the value of column FORMAT is \n",
         control)
  }

  return(df1)
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

  vcf <- RenameColumnsWithNameVAF(vcf)

  TUMOR <- vcf[["TUMOR"]]
  control <- unique(vcf[["FORMAT"]])
  alt     <- vcf[["ALT"]]
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

#' Read in the data lines of a Variant Call Format (VCF) file created by Mutect
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
#' @param tumor.col.name Name or index of the column in VCF which contains the
#'   tumor sample information. It \strong{must} have quotation marks if
#'   specifying the column name. If \code{tumor.col.name} is equal to
#'   \code{NA}(default), this function will use the 10th column to calculate
#'   VAFs. See \code{\link{GetMutectVAF}} for more details.
#'
#' @section Value: A data frame storing data lines of a VCF file with two
#'   additional columns added which contain the VAF(variant allele frequency)
#'   and read depth information.
#'
#' @keywords internal
ReadMutectVCF <-
  function(file, name.of.VCF = NULL, tumor.col.name = NA) {
    df <- MakeDataFrameFromVCF(file)
    if (is.null(name.of.VCF)) {
      vcf.name <- tools::file_path_sans_ext(basename(file))
    } else {
      vcf.name <- name.of.VCF
    }

    if (nrow(df) == 0) {
      return(df)
    } else {
      df1 <- GetMutectVAF(vcf = df, name.of.VCF = vcf.name,
                          tumor.col.name = tumor.col.name)
      return(df1)
    }
  }

#' @rdname GetVAF
#'
#' @export
GetMutectVAF <- function(vcf, name.of.VCF = NULL, tumor.col.name = NA) {
  stopifnot("data.frame" %in% class(vcf))

  if (nrow(vcf) == 0) {
    #vcf$VAF <- NA
    #vcf$read.depth <- NA
    return(vcf)
  }

  vcf <- RenameColumnsWithNameVAF(vcf)

  # Specify the possible variable names in Mutect VCF that stores count of reads
  # information

  # From header of Mutect VCF:
  # F1R2 = "Count of reads in F1R2 pair orientation supporting each allele"
  # F2R1 = "Count of reads in F2R1 pair orientation supporting each allele"
  # ALT_F1R2 = "Count of reads in F1R2 pair orientation supporting the alternate allele"
  # ALT_F2R1 = "Count of reads in F2R1 pair orientation supporting the alternate allele"
  # REF_F1R2 = "Count of reads in F1R2 pair orientation supporting the reference allele"
  # REF_F2R1 = "Count of reads in F2R1 pair orientation supporting the reference allele"
  type1 <- c("F1R2", "F2R1")
  type2 <- c("REF_F1R2", "ALT_F1R2", "REF_F2R1", "ALT_F2R1")

  vcf.format <- unlist(stringi::stri_split_fixed(vcf$FORMAT[1], ":"))
  is.type1 <- all(type1 %in% vcf.format)
  is.type2 <- all(type2 %in% vcf.format)
  if (!is.type1 && !is.type2) {
    warning("\nVCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " does not appear to be a Mutect VCF, please check the data")
    vcf$VAF <- NA
    vcf$read.depth <- NA
    return(vcf)
  }

  if (!is.na(tumor.col.name)) {
    if (is.character(tumor.col.name)) {
      if (!tumor.col.name %in% colnames(vcf)) {
        stop("\n", dQuote(tumor.col.name),
             " is not one of the column names in vcf ",
             ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)))
      }
    } else if (is.numeric(tumor.col.name)) {
        if (!tumor.col.name %in% 1:ncol(vcf)) {
          stop("\n", tumor.col.name,
               " is not one of the column indices in vcf ",
               ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)))
        }
    } else {
      stop("\n",
           "tumor.col.name should either be the colname name or column index in vcf ",
           ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)))
    }
  }

  GetTumorColumn <- function(vcf, tumor.col.name) {
    if (is.na(tumor.col.name)) {
      return(vcf[[10]])
    } else {
      return(vcf[[tumor.col.name]])
    }
  }
  tumor.col <- GetTumorColumn(vcf, tumor.col.name)

  ExtracVAFAndReadDepth <- function(tumor.col, type) {
    format.info <- stringi::stri_split_fixed(vcf$FORMAT, ":")
    read.counts.idx <- lapply(format.info, function(x) match(type, x))
    tumor.info.list <- stringi::stri_split_fixed(tumor.col, ":")

    Extract <- function(idx, tumor.info.list, read.counts.idx) {
      x <- tumor.info.list[[idx]]
      idx <- read.counts.idx[[idx]]
      as.integer(unlist(strsplit(x[idx], ",")))
    }
    num <- nrow(vcf)
    read.counts.info <- lapply(1:num, FUN = Extract,
                               tumor.info.list = tumor.info.list,
                               read.counts.idx = read.counts.idx)
    vafs <- sapply(read.counts.info, function(x) {
      vaf <- sum(x[c(2, 4)]) / sum(x)
    })
    read.depth <- sapply(read.counts.info, function(x) sum(x))
    return(data.frame(VAF = vafs, read.depth = read.depth))
  }

  vafs <- NULL
  if (is.type1) {
    vafs <- ExtracVAFAndReadDepth(tumor.col, type1)
  } else if (is.type2) {
    vafs <- ExtracVAFAndReadDepth(tumor.col, type2)
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
  CheckAndReturnVAFs(vafs)
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

  vcf <- RenameColumnsWithNameVAF(vcf)

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


#' @param vcf An in-memory VCF data frame.
#'
#' @param mc.cores The number of cores to use. Not available on Windows
#'   unless \code{mc.cores = 1}.
#'
#' @importFrom parallel mclapply
#'
#' @export
#' 
#' @note 
#' \code{\link{GetPCAWGConsensusVAF}} is analogous to \code{\link{GetMutectVAF}},
#' calculating VAF and read depth from PCAWG7 consensus vcfs
#' 
#' @rdname GetVAF
#'
GetPCAWGConsensusVAF <- function(vcf, mc.cores = 1) {

  vcf <- RenameColumnsWithNameVAF(vcf)

  info <- vcf$INFO
  tmp <- stringi::stri_split_fixed(info, ";")
  alt.counts <- parallel::mclapply(tmp, FUN = function(x) {
    idx <- grep("t_alt_count", x, fixed = TRUE)
    if (length(idx) == 0) {
      return(as.integer(NA))
    } else {
      alt.info <- x[idx]
      alt.count <- gsub("t_alt_count=", "", alt.info)
      return(as.integer(alt.count))
    }
  }, mc.cores = mc.cores)
  alt.counts1 <- unlist(alt.counts)

  ref.counts <- parallel::mclapply(tmp, FUN = function(x) {
    idx <- grep("t_ref_count", x, fixed = TRUE)
    if (length(idx) == 0) {
      return(as.integer(NA))
    } else {
      ref.info <- x[idx]
      ref.count <- gsub("t_ref_count=", "", ref.info)
      return(as.integer(ref.count))
    }
  }, mc.cores = mc.cores)
  ref.counts1 <- unlist(ref.counts)

  read.depth <- alt.counts1 + ref.counts1
  vaf <- alt.counts1/read.depth
  vcf$VAF <- vaf
  vcf$read.depth <- read.depth
  return(vcf)
}

#' Read in the data lines of a Variant Call Format (VCF) file
#'
#' @importFrom utils read.csv
#'
#' @param file The name/path of the VCF file, or a complete URL.
#'
#' @param variant.caller Name of the variant caller that produces the VCF, can
#'   be either \code{"strelka"}, \code{"mutect"}, \code{"freebayes"} or
#'   \code{"unknown"}. This information is needed to calculate the VAFs (variant
#'   allele frequencies). If \code{"unknown"}(default) and
#'   \code{get.vaf.function} is NULL, then VAF and read depth will be NAs.
#'
#' @param name.of.VCF Name of the VCF file. If \code{NULL}(default), this
#'   function will remove all of the path up to and including the last path
#'   separator (if any) in \code{file} and file path without extensions (and the
#'   leading dot) will be used as the name of the VCF file.
#'
#' @param tumor.col.name Optional. Only applicable to \strong{Mutect} VCF. Name
#'   or index of the column in \strong{Mutect} VCF which contains the tumor
#'   sample information. It \strong{must} have quotation marks if specifying the
#'   column name. If \code{tumor.col.name} is equal to \code{NA}(default), this
#'   function will use the 10th column to calculate VAFs. See
#'   \code{\link{GetMutectVAF}} for more details.
#'
#' @param filter.status The character string in column \code{FILTER} of the VCF
#'   that indicates that a variant has passed all the variant caller's filters.
#'   Variants (lines in the VCF) for which the value in column \code{FILTER}
#'   does not equal \code{filter.status} are silently excluded from the output.
#'   If \code{NULL}, all variants are retained. In almost all cases, the default
#'   value of \code{"PASS"} is what the user would want. If \code{filter.status}
#'   is not \code{NULL} but the input file does not contain the column \code{FILTER}
#'   \code{filter.status} is ignored with a warning.
#'
#' @param get.vaf.function Optional. Only applicable when \code{variant.caller} is
#' \strong{"unknown"}. Function to calculate VAF(variant allele frequency) and read
#'   depth information from original VCF. See \code{\link{GetMutectVAF}} as an example.
#'   If \code{NULL}(default) and \code{variant.caller} is "unknown", then VAF
#'   and read depth will be NAs.
#'
#' @param ... Optional arguments to \code{get.vaf.function}.
#'
#' @param num.of.cores The number of cores to use. Not available on Windows
#'   unless \code{num.of.cores = 1}.
#'
#' @return A data frame storing data lines of the VCF file with two additional
#'   columns added which contain the VAF(variant allele frequency) and read
#'   depth information.
#'
#' @keywords internal
ReadVCF <-
  function(file, variant.caller = "unknown", name.of.VCF = NULL, tumor.col.name = NA,
           filter.status = "PASS", get.vaf.function = NULL, ...) {
    df0 <- MakeDataFrameFromVCF(file) # , name.of.VCF = name.of.VCF)

    if (nrow(df0) == 0) {
      return(df0)
    }

    # Remove rows that don't have the specified filter status
    if (is.null(filter.status)) {
      df <- df0
    } else {
      # Check whether df0 has column name "FILTER"
      if (!"FILTER" %in% colnames(df0)) {
        warning("\nThere is no column FILTER in the file ", file,
                "\nargument filter.status is ignored and all variants will be retained")
        df <- df0
      } else {
        df <- dplyr::filter(df0, FILTER == filter.status)
      }
    }
    rm(filter.status)

    if (nrow(df) == 0) {
      return(df)
    }

    df1 <- df
    df1$VAF <- as.numeric(NA)
    df1$read.depth <- as.numeric(NA)

    if (variant.caller == "unknown") {
      if (is.null(get.vaf.function)) {
        return(df1)
      } else {
        df2 <- get.vaf.function(df, ...)
        return(df2)
      }
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

      # Check for any SBS in df and only calculate VAF for those SBS variants
      SBS.idx0 <- which(nchar(df$REF) == 1 & nchar(df$ALT) == 1)
      SBS.multiple.alt <-
        which(nchar(df$REF) == 1 & grepl(",", df$ALT, fixed = TRUE))
      SBS.idx <- c(SBS.idx0, SBS.multiple.alt)
      if (length(SBS.idx) == 0) {
        return(df)
      } else {
        SBS.df <- df[SBS.idx, ]
        SBS.df1 <- GetStrelkaVAF(vcf = SBS.df, name.of.VCF = vcf.name)
        df1[SBS.idx, ]$VAF <- SBS.df1$VAF
        df1[SBS.idx, ]$read.depth <- SBS.df1$read.depth
        return(df1)
      }
    }

    if (variant.caller == "mutect") {
      df2 <- GetMutectVAF(vcf = df, name.of.VCF = vcf.name,
                          tumor.col.name = tumor.col.name)
      return(df2)
    }

    if (variant.caller == "freebayes") {
      # Check for any SBS in df and only calculate VAF for SBS variants
      SBS.idx0 <- which(nchar(df$REF) == 1 & nchar(df$ALT) == 1)
      SBS.multiple.alt <-
        which(nchar(df$REF) == 1 & grepl(",", df$ALT, fixed = TRUE))
      SBS.idx <- c(SBS.idx0, SBS.multiple.alt)
      if (length(SBS.idx) == 0) {
        return(df)
      } else {
        SBS.df <- df[SBS.idx, ]
        SBS.df1 <- GetFreebayesVAF(vcf = SBS.df, name.of.VCF = vcf.name)
        df1[SBS.idx, ]$VAF <- SBS.df1$VAF
        df1[SBS.idx, ]$read.depth <- SBS.df1$read.depth
        return(df1)
      }
    }
  }

#' Read VCF files
#'
#' @inheritParams ReadAndSplitVCFs
#'
#' @importFrom parallel mclapply
#'
#' @return A list of data frames storing data lines of the VCF files with two
#'   additional columns added which contain the VAF(variant allele frequency)
#'   and read depth information.
#'
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadVCFs(file, variant.caller = "mutect")
ReadVCFs <- function(files, variant.caller = "unknown", num.of.cores = 1,
                     names.of.VCFs = NULL,
                     tumor.col.names = NA, filter.status = "PASS",
                     get.vaf.function = NULL, ...) {
  num.of.cores <- AdjustNumberOfCores(num.of.cores)

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
    ReadVCF(file = files[idx],
            variant.caller = variant.caller,
            name.of.VCF = vector1[idx],
            tumor.col.name = vector2[idx],
            filter.status = filter.status,
            get.vaf.function = get.vaf.function,
            ... = ...)
  }

  vcfs <- parallel::mclapply(1:num.of.files, FUN = ReadVCF1, files = files,
                             variant.caller = variant.caller,
                             vector1 = vcfs.names, vector2 = tumor.col.names,
                             mc.cores = num.of.cores)
  names(vcfs) <- vcfs.names
  return(vcfs)
}

#' @keywords internal
CheckAndRemoveDiscardedVariants <- function(vcf, name.of.VCF = NULL) {
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

  # Remove rows with unstandardized chromosome names
  retval2 <- StandardChromNameNew(df = df2, name.of.VCF = name.of.VCF)
  df3 <- retval2$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval2$discarded.variants)

  # VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", df3$ALT, fixed = TRUE)
  if (length(multiple.alt) > 0) {
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants with multiple alternative alleles and were ",
            "discarded. See discarded.variants in the return value for more ",
            "details.")
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
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants involving three or more nucleotides and were ",
            "discarded. See discarded.variants in the return value for more ",
            "details.")
    df5 <- df4[-other.df, ]
    df5.to.remove <- df4[other.df, ]
    df5.to.remove$discarded.reason <- "Variant involves three or more nucleotides"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, df5.to.remove)
  } else {
    df5 <- df4
  }

  # Remove complex indels
  complex.indels <- which((nchar(df5$REF) != nchar(df5$ALT)) &
                            (substr(df5$REF, 1, 1) != substr(df5$ALT, 1, 1)))
  if (length(complex.indels) > 0) {
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has complex indels and were discarded. See discarded.variants ",
            "in the return value for more details.")
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
  wrong.DBS.type1 <- dplyr::filter(df6, nchar(REF) == 2, nchar(ALT) == 2,
                                   substr(REF, 1, 1) == substr(ALT, 1, 1))
  wrong.DBS.type2 <- dplyr::filter(df6, nchar(REF) == 2, nchar(ALT) == 2,
                                   substr(REF, 2, 2) == substr(ALT, 2, 2))
  wrong.DBS <- dplyr::bind_rows(wrong.DBS.type1, wrong.DBS.type2)

  if (nrow(wrong.DBS) > 0) {
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has wrong DBS variants and were discarded. See discarded.variants ",
            "in the return value for more details.")
    wrong.DBS.pos <- wrong.DBS$POS
    wrong.DBS$discarded.reason <- "Wrong DBS variant"
    discarded.variants <-
      dplyr::bind_rows(discarded.variants, wrong.DBS)
    df7 <- dplyr::filter(df6, !POS %in% wrong.DBS.pos)
  } else {
    df7 <- df6
  }

  if (nrow(discarded.variants) == 0) {
    return(list(df = df7))
  } else {
    return(list(df = df7, discarded.variants = discarded.variants))
  }
}

#' @title Split a mutect2 VCF into SBS, DBS, and ID VCFs, plus a list of other mutations
#'
#' @param vcf.df An in-memory data.frame representing a Mutect VCF, including
#'  VAFs, which are added by \code{\link{ReadMutectVCF}}.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @return A list with 3 in-memory VCFs and discarded variants that were not
#'   incorporated into the first 3 VCFs:
#'
#'  * \code{SBS}: VCF with only single base substitutions.
#'
#'  * \code{DBS}: VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  * \code{ID}: VCF with only small insertions and deletions.
#'
#'  * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#'  that were excluded from the analysis. See the added extra column
#'  \code{discarded.reason} for more details.
#'  @md
#'
#' @keywords internal
SplitOneMutectVCF <- function(vcf.df, name.of.VCF = NULL) {
  if (nrow(vcf.df) == 0) {
    return(list(SBS = vcf.df, DBS = vcf.df, ID = vcf.df))
  }

  # Create an empty data frame for discarded variants
  discarded.variants <- vcf.df[0, ]

  # Check and remove discarded variants
  retval <-
    CheckAndRemoveDiscardedVariants(vcf = vcf.df, name.of.VCF = name.of.VCF)
  df <- retval$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval$discarded.variants)

  SBS.df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]
  DBS.df <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]
  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]

  if (nrow(discarded.variants) == 0) {
    return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df))
  } else {
    return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
                discarded.variants = discarded.variants))
  }
}

#' Split each Mutect VCF into SBS, DBS, and ID VCFs (plus
#' VCF-like data frame with left-over rows)
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

    GetSplitMutectVCFs <- function(idx, list.of.vcfs) {
      split.vcfs <- SplitOneMutectVCF(list.of.vcfs[[idx]],
                                      name.of.VCF = names(list.of.vcfs)[idx])
      return(split.vcfs)
    }
    num.of.vcfs <- length(list.of.vcfs)
    if (suppress.discarded.variants.warnings == TRUE) {
      v1 <- suppressWarnings(lapply(1:num.of.vcfs, GetSplitMutectVCFs,
                                    list.of.vcfs = list.of.vcfs))
    } else {
      v1 <- lapply(1:num.of.vcfs, GetSplitMutectVCFs,
                   list.of.vcfs = list.of.vcfs)
    }
    names(v1) <- names.of.VCFs
    SBS <- lapply(v1, function(x) x$SBS)
    DBS <- lapply(v1, function(x) x$DBS)
    ID  <- lapply(v1, function(x) x$ID)
    discarded.variants <- lapply(v1, function(x) x$discarded.variants)

    # Remove NULL elements from discarded.variants
    discarded.variants1 <- Filter(Negate(is.null), discarded.variants)

    if (length(discarded.variants1) == 0) {
      return(list(SBS = SBS, DBS = DBS, ID = ID))
    } else {
      return(list(SBS = SBS, DBS = DBS, ID = ID,
                  discarded.variants = discarded.variants1))
    }
  }

#' Split an in-memory SBS VCF into pure SBSs, pure DBSs, and variants involving
#' > 2 consecutive bases
#'
#' SBSs are single base substitutions,
#' e.g. C>T, A>G,....  DBSs are double base substitutions,
#' e.g. CC>TT, AT>GG, ...  Variants involving > 2 consecutive
#' bases are rare, so this function just records them. These
#' would be variants such ATG>CCT, AGAT>TCTA, ...
#'
#' @param vcf.df An in-memory data frame containing an SBS VCF file contents.
#'
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02. If
#'   the absolute difference of VAFs for adjacent SBSs is bigger than
#'   \code{max.vaf.diff}, then these adjacent SBSs are likely to be "merely"
#'   asynchronous single base mutations, opposed to a simultaneous doublet
#'   mutation or variants involving more than two consecutive bases.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param always.merge.SBS If \code{TRUE} merge adjacent SBSs as DBSs
#'   regardless of VAFs and regardless of the value of \code{max.vaf.diff}.
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
#'    \item \code{discarded.variants}: \strong{Non-NULL only if} there are
#'    variants that were excluded from the analysis. See the added extra column
#'    \code{discarded.reason} for more details.
#'
#'    }
#'
#' @keywords internal
SplitSBSVCF <- function(vcf.df,
                        max.vaf.diff = 0.02,
                        name.of.VCF = NULL,
                        always.merge.SBS) {
  stopifnot("data.frame" %in% class(vcf.df))

  if (nrow(vcf.df) == 0) {
    return(list(SBS.vcf = vcf.df, DBS.vcf = vcf.df))
  }

  # Create an empty data frame for discarded variants
  discarded.variants <- vcf.df[0, ]

  # Check and remove discarded variants
  retval <-
    CheckAndRemoveDiscardedVariants(vcf = vcf.df, name.of.VCF = name.of.VCF)
  vcf.df <- retval$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval$discarded.variants)

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

  if (!always.merge.SBS) {
    # Keep only SBS pairs that have very similar VAFs (variant allele
    # frequencies). If VAFs are not similar, the adjacent SBSs are likely to be
    # "merely" asynchronous single base mutations, opposed to a simultaneous
    # doublet mutation.
    non.SBS <- dt2[abs(VAF.x - VAF.y) <= max.vaf.diff]
  } else {
    non.SBS <-dt2
  }
  rm(dt2)

  if (nrow(non.SBS) == 0) {
    # There are no non.SBS mutations in the input.
    # Everything in vcf.df is an SBS. We are finished.
    empty <- vcf.df[-(1:nrow(vcf.df)), ]
    if (nrow(discarded.variants) == 0) {
      return(list(SBS.vcf = vcf.df, DBS.vcf = empty))
    } else {
      return(list(SBS.vcf = vcf.df, DBS.vcf = empty,
                  discarded.variants = discarded.variants))
    }

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
  DBSx$CHROM <- as.character(DBSx$CHROM)
  DBS.vcf.df <- MakeVCFDBSdf(DBSx, vcf.dt)
  num.DBS.out <- nrow(DBS.vcf.df)

  other.ranges <- DBS.plus[DBS.plus$width > 2, ]
  if (nrow(other.ranges) > 0) {
    colnames(other.ranges)[1:3] <- c("CHROM", "LOW.POS", "HIGH.POS")
    other.ranges$discarded.reason <- "Variants that do not represent SBS or DBS"
    if (nrow(discarded.variants) == 0) {
      discarded.variants <- other.ranges
    } else {
      discarded.variants <- dplyr::bind_rows(discarded.variants, other.ranges)
    }
    warning("VCF ", ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
            " has variants involving three or more nucleotides and were ",
            "discarded. See discarded.variants in the return value for more ",
            "details.")
  }

  num.other <- sum(other.ranges$width)

  if ((num.SBS.out + 2 * num.DBS.out + num.other) != num.in) {
    warning("Counts are off:", num.SBS.out, 2*num.DBS.out, num.other, "vs", num.in, "\n")
  }

  if (nrow(discarded.variants) == 0) {
    return(list(SBS.vcf = out.SBS.df, DBS.vcf = DBS.vcf.df))
  } else {
    return(list(SBS.vcf = out.SBS.df, DBS.vcf = DBS.vcf.df,
                discarded.variants = discarded.variants))
  }
}

#' @title Split a VCF into SBS, DBS, and ID VCFs, plus a list of other mutations
#'
#' @param vcf.df An in-memory data.frame representing a VCF, including
#'  VAFs, which are added by \code{\link{ReadVCF}}.
#'
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02. If
#'   the absolute difference of VAFs for adjacent SBSs is bigger than
#'   \code{max.vaf.diff}, then these adjacent SBSs are likely to be "merely"
#'   asynchronous single base mutations, opposed to a simultaneous doublet
#'   mutation or variants involving more than two consecutive bases.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param always.merge.SBS If \code{TRUE} merge adjacent SBSs as DBSs
#'   regardless of VAFs and regardless of the value of \code{max.vaf.diff}.
#'
#' @return A list with 3 in-memory VCFs and discarded variants that were not
#'   incorporated into the first 3 VCFs:
#'
#'  * \code{SBS}: VCF with only single base substitutions.
#'
#'  * \code{DBS}: VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  * \code{ID}: VCF with only small insertions and deletions.
#'
#'  * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#'  that were excluded from the analysis. See the added extra column
#'  \code{discarded.reason} for more details.
#'  @md
#'
#' @keywords internal
SplitOneVCF <- function(vcf.df,
                        max.vaf.diff     = 0.02,
                        name.of.VCF      = NULL,
                        always.merge.SBS = FALSE) {
  if (nrow(vcf.df) == 0) {
    return(list(SBS = vcf.df, DBS = vcf.df, ID = vcf.df))
  }

  # Create an empty data frame for discarded variants
  discarded.variants <- vcf.df[0, ]

  # Check and remove discarded variants
  retval <-
    CheckAndRemoveDiscardedVariants(vcf = vcf.df, name.of.VCF = name.of.VCF)
  df <- retval$df
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, retval$discarded.variants)

  SBS.df0 <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]

  # Try to get DBS from adjacent SBSs according to similar VAFs
  split.dfs <- SplitSBSVCF(vcf.df           = SBS.df0,
                           max.vaf.diff     = max.vaf.diff,
                           name.of.VCF      = name.of.VCF,
                           always.merge.SBS = always.merge.SBS
                           )
  SBS.df <- split.dfs$SBS.vcf
  DBS.df0 <- split.dfs$DBS.vcf
  discarded.variants <-
    dplyr::bind_rows(discarded.variants, split.dfs$discarded.variants)

  DBS.df1 <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]
  DBS.df <- dplyr::bind_rows(DBS.df0, DBS.df1)
  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]

  if (nrow(discarded.variants) == 0) {
    return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df))
  } else {
    return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
                discarded.variants = discarded.variants))
  }
}

#' Split each VCF into SBS, DBS, and ID VCFs (plus
#' VCF-like data frame with left-over rows)
#'
#' @param list.of.vcfs List of VCFs as in-memory data frames. The VCFs should
#'   have \code{VAF} and \code{read.depth} information added. See
#'   \code{ReadVCFs} for more details.
#'
#' @param variant.caller Name of the variant caller that produces the VCF, can
#'   be either \code{"strelka"}, \code{"mutect"}, \code{"freebayes"} or
#'   \code{"unknown"}. If variant caller is \code{"mutect"}, do \strong{not} merge
#'   SBSs into DBS.
#'
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02. If
#'   the absolute difference of VAFs for adjacent SBSs is bigger than
#'   \code{max.vaf.diff}, then these adjacent SBSs are likely to be "merely"
#'   asynchronous single base mutations, opposed to a simultaneous doublet
#'   mutation or variants involving more than two consecutive bases.
#'
#' @param num.of.cores The number of cores to use. Not available on Windows
#'   unless \code{num.of.cores = 1}.
#'
#' @param always.merge.SBS If \code{TRUE} merge adjacent SBSs as DBSs
#'   regardless of VAFs and regardless of the value of \code{max.vaf.diff}
#'   and regardless of the value of \code{get.vaf.function}. It is an
#'   error to set this to \code{TRUE} when \code{variant.caller = "mutect"}.
#'
#' @inheritParams ReadAndSplitMutectVCFs
#'
#' @inheritSection ReadAndSplitMutectVCFs Value
#'
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Mutect-vcf",
#'                       "Mutect.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadVCFs(file, variant.caller = "mutect")
#' split.vcfs <- SplitListOfVCFs(list.of.vcfs, variant.caller = "mutect")
SplitListOfVCFs <-function(list.of.vcfs,
                           variant.caller,
                           max.vaf.diff = 0.02,
                           num.of.cores = 1,
                           suppress.discarded.variants.warnings = TRUE,
                           always.merge.SBS                     = FALSE) {
  names.of.VCFs <- names(list.of.vcfs)

  GetSplitVCFs <- function(idx, list.of.vcfs, variant.caller) {
    if (variant.caller == "mutect") {
      if (always.merge.SBS) {
        stop("always.merge.SBS must be FALSE when variant.caller = \"mutect\"")
      }
      split.vcfs <- SplitOneMutectVCF(vcf.df = list.of.vcfs[[idx]],
                                      name.of.VCF = names(list.of.vcfs)[idx])
    } else {
      split.vcfs <- SplitOneVCF(list.of.vcfs[[idx]],
                                max.vaf.diff     = max.vaf.diff,
                                name.of.VCF      = names(list.of.vcfs)[idx],
                                always.merge.SBS = always.merge.SBS)
    }

    return(split.vcfs)
  }
  num.of.vcfs <- length(list.of.vcfs)
  if (suppress.discarded.variants.warnings == TRUE) {
    v1 <- suppressWarnings(parallel::mclapply(1:num.of.vcfs, GetSplitVCFs,
                                              list.of.vcfs = list.of.vcfs,
                                              variant.caller = variant.caller,
                                              mc.cores = num.of.cores))
  } else {
    v1 <- parallel::mclapply(1:num.of.vcfs, GetSplitVCFs,
                             list.of.vcfs = list.of.vcfs,
                             variant.caller = variant.caller,
                             mc.cores = num.of.cores)
  }
  names(v1) <- names.of.VCFs
  SBS <- lapply(v1, function(x) x$SBS)
  DBS <- lapply(v1, function(x) x$DBS)
  ID  <- lapply(v1, function(x) x$ID)
  discarded.variants <- lapply(v1, function(x) x$discarded.variants)

  # Remove NULL elements from discarded.variants
  discarded.variants1 <- Filter(Negate(is.null), discarded.variants)

  if (length(discarded.variants1) == 0) {
    return(list(SBS = SBS, DBS = DBS, ID = ID))
  } else {
    return(list(SBS = SBS, DBS = DBS, ID = ID,
                discarded.variants = discarded.variants1))
  }
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
AddSeqContext <-
  function(df, ref.genome, seq.context.width = 10, name.of.VCF = NULL) {
  if (0 == nrow(df)) return(df)
  ref.genome <- NormalizeGenomeArg(ref.genome)

  # Check if the format of sequence names in df and genome are the same
  chr.names <- CheckAndFixChrNames(vcf.df = df,
                                   ref.genome = ref.genome,
                                   name.of.VCF = name.of.VCF)

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
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @import data.table
#'
#' @importFrom dplyr %>% group_by mutate
#'
#' @return A data frame with new columns added to the input data frame,
#'     which contain the mutated gene's name, range and strand information.
#'
#' @keywords internal
AddTranscript <-
  function(df, trans.ranges = NULL, ref.genome, name.of.VCF = NULL) {
    if (nrow(df) == 0) {
      return(df)
    }

    # Sometimes CHROM is numeric, but this breaks foverlaps, below.
    df$CHROM <- as.character(df$CHROM)

    if (is.null(trans.ranges)) {
      return(data.table(df))
    }

    df <- RenameColumnsWithNameStrand(df)

    ref.genome <- NormalizeGenomeArg(ref.genome = ref.genome)

    # Check whether the chromosome name format of trans.ranges matches with that
    # in df. If not, change chromosome name format in trans.ranges
    new.chr.names <-
      CheckAndFixChrNamesForTransRanges(trans.ranges = trans.ranges,
                                        vcf.df = df,
                                        ref.genome = ref.genome,
                                        name.of.VCF = name.of.VCF)
    trans.ranges$chrom <- new.chr.names

    # We need to set key for trans.ranges for using data.table::foverlaps
    if (!data.table::haskey(trans.ranges)) {
      data.table::setkeyv(trans.ranges, c("chrom", "start", "end"))
    }

    # Find range overlaps between the df and trans.ranges
    df1 <- data.table(df)
    df1[, POS2 := POS]
    dt <- data.table::foverlaps(df1, trans.ranges,
                                by.x = c("CHROM", "POS", "POS2"),
                                type = "within", mult = "all")

    # Find out mutations that fall on transcripts on both strands
    #dt1 <- dt[, bothstrand := "+" %in% strand && "-" %in% strand,
    #          by = .(CHROM, ALT, POS)] # Note that is important to have
    # ALT in the by list because in a few cases
    # there are multiple ALT alleles at one POS.

    dt1 <- dt %>% dplyr::group_by(CHROM, ALT, POS) %>%
      dplyr::mutate(bothstrand = "+" %in% strand && "-" %in% strand)
    data.table::setDT(dt1)

    # Count the number of transcript ranges where a particular mutation
    # falls into
    dt2 <- dt1 %>% dplyr::group_by(CHROM, ALT, POS) %>% dplyr::mutate(count = dplyr::n())
    #dt2 <- dt1[, count := .N, by = .(CHROM, ALT, POS)]
    data.table::setDT(dt2)

    # Swap gene location according to strand information
    dt3 <- dt2[strand == "-", c("end", "start") := .(start, end)]

    # Reorder the columns of dt3
    df.colnames <- colnames(df)
    trans.ranges.colnames <- colnames(trans.ranges)[-1]
    data.table::setcolorder(dt3, neworder = c(df.colnames, trans.ranges.colnames))

    # Rename some of the columns in dt3
    data.table::setnames(dt3,
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
  # tmpvcf <- SBS.vcf.dt[ , c("CHROM", "POS", "REF", "ALT", "VAF", "read.depth")]
  tmpvcf <- SBS.vcf.dt
  DBS.range.dt <- as.data.table(DBS.range.df)
  tmp1 <- merge(DBS.range.dt, tmpvcf,
                by.x = c("CHROM", "LOW"),
                by.y = c("CHROM", "POS"))
  tmp2 <- merge(tmp1, tmpvcf,
                by.x = c("CHROM", "HIGH"),
                by.y = c("CHROM", "POS"))
  # Calculate the read depth for DBS from merged SBS
  # tmp2[, DP.x := as.integer(sapply(strsplit(TUMOR.x, ":"), "[", 1))]
  # tmp2[, DP.y := as.integer(sapply(strsplit(TUMOR.y, ":"), "[", 1))]
  tmp2[, read.depth := pmin(read.depth.x, read.depth.y)]

  tmp2[, VAF := rowMeans(cbind(VAF.x, VAF.y))]
  tmp2[, POS := LOW]
  tmp2[, remark.for.DBS := "From merged SBSs"]
  tmp2[, REF := paste0(REF.x, REF.y)]
  tmp2[, ALT := paste0(ALT.x, ALT.y)]

  # Delete some of the columns
  tmp2[, c("read.depth.x", "read.depth.y", "VAF.x", "VAF.y", "LOW", "HIGH",
           "REF.x", "REF.y", "ALT.x", "ALT.y", "POS.plus.one.x",
           "POS.plus.one.y") := NULL]

  old.col.names <- setdiff(colnames(SBS.vcf.dt), "POS.plus.one")
  col.names.order1 <-
    c("CHROM", "POS", "REF", "ALT", "VAF", "read.depth", "remark.for.DBS")
  col.names.order2 <- setdiff(old.col.names, col.names.order1)

  # Retrieve the column information for DBS from original SBS VCF
  for (name in col.names.order2) {
    name1 <- paste0(name, c(".x", ".y"))

    # Get the unique column information for DBS
    GetUniqueInformation <- function(x) {
      y <- paste(unique(unlist(x)), collapse = ",")
      class(y) <- class(unique(unlist(x)))
      return(y)
    }
    tmp2[, (name) := apply(X = .SD, MARGIN = 1,
                           FUN = GetUniqueInformation), .SDcols = name1]

    # Delete the redundant columns
    tmp2[, (name1) := NULL]
  }

  # Return the DBS data table according to specific column orders
  col.names.order <- c(col.names.order1, col.names.order2)
  return(tmp2[, ..col.names.order])
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
#' @param max.vaf.diff The maximum difference of VAF, default value is 0.02. If
#'   the absolute difference of VAFs for adjacent SBSs is bigger than
#'   \code{max.vaf.diff}, then these adjacent SBSs are likely to be "merely"
#'   asynchronous single base mutations, opposed to a simultaneous doublet
#'   mutation or variants involving more than two consecutive bases.
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
#'    \item \code{discarded.variants}: \strong{Non-NULL only if} there are
#'    variants that were excluded from the analysis. See the added extra column
#'    \code{discarded.reason} for more details.
#'
#'    }
#'
#' @keywords internal
SplitStrelkaSBSVCF <- function(vcf.df, max.vaf.diff = 0.02, name.of.VCF = NULL, always.merge.SBS = FALSE) {
  stopifnot("data.frame" %in% class(vcf.df))

  retval <- SplitSBSVCF(vcf.df = vcf.df,
                        max.vaf.diff = max.vaf.diff,
                        name.of.VCF = name.of.VCF,
                        always.merge.SBS = always.merge.SBS)
  return(retval)
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

    GetSplitStrelkaSBSVCFs <- function(idx, list.of.vcfs) {
      split.vcfs <- SplitStrelkaSBSVCF(list.of.vcfs[[idx]],
                                       name.of.VCF = names(list.of.vcfs)[idx])
      return(split.vcfs)
    }
    num.of.vcfs <- length(list.of.vcfs)
    if (suppress.discarded.variants.warnings == TRUE) {
      split.vcfs <-
        suppressWarnings(lapply(1:num.of.vcfs, GetSplitStrelkaSBSVCFs,
                                list.of.vcfs = list.of.vcfs))
    } else {
      split.vcfs <- lapply(1:num.of.vcfs, GetSplitStrelkaSBSVCFs,
                           list.of.vcfs = list.of.vcfs)
    }
    names(split.vcfs) <- names.of.VCFs
    SBS.vcfs   <- lapply(split.vcfs, function(x) x$SBS.vcf)
    DBS.vcfs   <- lapply(split.vcfs, function(x) x$DBS.vcf)
    discarded.variants <- lapply(split.vcfs, function(x) x$discarded.variants)

    # Remove NULL elements from discarded.variants
    discarded.variants1 <- Filter(Negate(is.null), discarded.variants)

    if (length(discarded.variants1) == 0) {
      return(list(SBS.vcfs = SBS.vcfs, DBS.vcfs = DBS.vcfs))
    } else {
      return(list(SBS.vcfs = SBS.vcfs, DBS.vcfs = DBS.vcfs,
                  discarded.variants = discarded.variants1))
    }
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

  vcf <- data.table::as.data.table(vcf)
  # use .. notation to find column.to.use as a vector of column positions,
  # like it would work in data.frame
  cut.pos <- 1 + (nchar(unlist(vcf[, ..column.to.use])) - 1) / 2
  stopifnot(cut.pos == round(cut.pos))
  cut.from.ref <- substr(unlist(vcf[, ..column.to.use]), cut.pos,
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
ReadStrelkaSBSVCFs <- function(files, names.of.VCFs = NULL) {

  vcfs <-
    lapply(files, FUN = ReadStrelkaSBSVCF, name.of.VCF = names.of.VCFs)

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
#' @param tumor.col.names Vector of column names or column indices in VCFs which
#'   contain the tumor sample information. The order of elements in
#'   \code{tumor.col.names} should match the order of VCFs specified in
#'   \code{files}. If \code{tumor.col.names} is equal to \code{NA}(default),
#'   this function will use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link{GetMutectVAF}} for more details.
#'
#' @section Value: A list of data frames which store data lines of VCF files
#'   with two additional columns added which contain the VAF(variant allele
#'   frequency) and read depth information.
#'
#' @keywords internal
ReadMutectVCFs <-
  function(files, names.of.VCFs = NULL, tumor.col.names = NA) {
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

  GetMutectVCFs <- function(idx, files, names.of.VCFs, tumor.col.names) {
    ReadMutectVCF(file = files[idx], name.of.VCF = names.of.VCFs[idx],
                  tumor.col.name = tumor.col.names[idx])
  }

  vcfs <- lapply(1:num.of.files, FUN = GetMutectVCFs,
                 files = files, names.of.VCFs = vcfs.names,
                 tumor.col.names = tumor.col.names)

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
#' @param name.of.VCF Name of the VCF file.
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
AnnotateSBSVCF <- function(SBS.vcf, ref.genome,
                           trans.ranges = NULL, name.of.VCF = NULL) {
  SBS.vcf <- AddSeqContext(df = SBS.vcf, ref.genome = ref.genome,
                           name.of.VCF = name.of.VCF)
  CheckSeqContextInVCF(SBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    SBS.vcf <- AddTranscript(df = SBS.vcf, trans.ranges = trans.ranges,
                             ref.genome = ref.genome,
                             name.of.VCF = name.of.VCF)
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
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
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
            AddAndCheckSBSClassInVCF(vcf, mat96, mat1536, mat192, sample.id)
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
CreateOneColSBSMatrix <- function(vcf, sample.id = "count",
                                  return.annotated.vcf = FALSE) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude DBSs
  # before calling the function. This function does not detect DBSs.

  CheckForEmptySBSVCF <- function(vcf, return.annotated.vcf) {
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
    } else {
      return(FALSE)
    }
  }

  ret1 <- CheckForEmptySBSVCF(vcf = vcf,
                              return.annotated.vcf = return.annotated.vcf)
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
      paste0('SBS variant whose reference base in ref.genome does not match the',
             ' reference base in the VCF file.')
    message("In sample ", sample.id, " ", length(mismatches), " row out of ",
            nrow(vcf), " had reference base in ref.genome that does not match the ",
            "reference base in the VCF file.\n",
            "Please check the ref.genome argument.\n",
            "See discarded.variants in the return value for more details")
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
        'Variants in the SBS vcf ', sample.id,
        ' whose pentanucleotide context contains "N" ',
        'have been deleted so as not to conflict with downstream processing. ',
        'See discarded.variants in the return value for more details.')
    }

  ret2 <- CheckForEmptySBSVCF(vcf = vcf,
                              return.annotated.vcf = return.annotated.vcf)
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
  vcf1 <- vcf %>% dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], pyr.mut = pyr.mut[1])

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
    retval <-
      CheckAndReturnSBSMatrix(vcf = vcf0, discarded.variants = discarded.variants,
                              mat96 = mat96, mat1536 = mat1536, mat192 = NULL,
                              return.annotated.vcf = return.annotated.vcf,
                              sample.id = sample.id)
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
  vcf3 <- vcf2 %>% dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], mutation = mutation[1], trans.strand = trans.strand[1])

  # If vcf3 has empty rows, we will return 1-column SBS192 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    mat192 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$SBS192), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$SBS192, sample.id))
    retval <-
      CheckAndReturnSBSMatrix(vcf = vcf0, discarded.variants = discarded.variants,
                              mat96 = mat96, mat1536 = mat1536, mat192 = mat192,
                              return.annotated.vcf = return.annotated.vcf,
                              sample.id = sample.id)
    return(retval)
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

  CheckAndReturnSBSMatrix(vcf = vcf0, discarded.variants = discarded.variants,
                          mat96 = mat96, mat1536 = mat1536, mat192 = mat192,
                          return.annotated.vcf = return.annotated.vcf,
                          sample.id = sample.id)
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
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' DBS.vcf <- list.of.vcfs$DBS.vcfs[[1]]
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.DBS.vcf <- AnnotateDBSVCF(DBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)}
AnnotateDBSVCF <- function(DBS.vcf, ref.genome,
                           trans.ranges = NULL, name.of.VCF = NULL) {
  DBS.vcf <- AddSeqContext(df = DBS.vcf, ref.genome = ref.genome,
                           name.of.VCF = name.of.VCF)
  CheckSeqContextInVCF(DBS.vcf, "seq.21bases")
  trans.ranges <- InferTransRanges(ref.genome, trans.ranges)
  if (!is.null(trans.ranges)) {
    DBS.vcf <- AddTranscript(df = DBS.vcf, trans.ranges = trans.ranges,
                             ref.genome = ref.genome,
                             name.of.VCF = name.of.VCF)
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
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
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
            AddAndCheckDBSClassInVCF(vcf, mat78, mat136, mat144, sample.id)
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
CreateOneColDBSMatrix <- function(vcf, sample.id = "count",
                                  return.annotated.vcf = FALSE) {
  # Error checking:
  # This function cannot handle insertion, deletions, or complex indels,
  # Therefore we check for this problem; but we need to exclude SBSs
  # before calling the function. This function does not detect SBSs.

  CheckForEmptyDBSVCF <- function(vcf, return.annotated.vcf) {
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
               dimnames = list(ICAMS::catalog.row.order$DBS144, sample.id))
      if (return.annotated.vcf == FALSE) {
        return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                    catDBS144 = catDBS144))
      } else {
        return(list(catDBS78 = catDBS78, catDBS136 = catDBS136,
                    catDBS144 = catDBS144, annotated.vcf = vcf))
      }
    } else {
      return(FALSE)
    }
  }

  ret1 <- CheckForEmptyDBSVCF(vcf = vcf,
                              return.annotated.vcf = return.annotated.vcf)
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
      'Variants in the DBS vcf ', sample.id,
      ' whose tetranucleotide context contains "N" ',
      'have been deleted so as not to conflict with downstream processing. ',
      'See discarded.variants in the return value for more details.')
  }

  ret2 <- CheckForEmptyDBSVCF(vcf = vcf,
                              return.annotated.vcf = return.annotated.vcf)
  if (!is.logical(ret2)) {
    return(ret2)
  }

  # One DBS mutation can be represented by more than 1 row in vcf after annotated by
  # AnnotateDBSVCF function if the mutation position falls into the range of
  # multiple transcripts. When creating the 78 and 136 catalog, we only need to
  # count these mutations once.
  # vcf1 <- vcf[, .(REF = REF[1], seq.21bases = seq.21bases[1]),
  #            by = .(CHROM, ALT, POS)]
  vcf1 <- vcf %>% dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], seq.21bases = seq.21bases[1])

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
    retval <-
      CheckAndReturnDBSMatrix(vcf = vcf, discarded.variants = discarded.variants,
                              mat78 = DBS.mat.78, mat136 = DBS.mat.136,
                              mat144 = NULL,
                              return.annotated.vcf = return.annotated.vcf,
                              sample.id = sample.id)
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
  vcf3 <- vcf2 %>% dplyr::group_by(CHROM, ALT, POS) %>%
    dplyr::summarise(REF = REF[1], trans.strand = trans.strand[1])

  # If vcf3 has empty rows, we will return 1-column DBS144 matrix with all
  # values being 0 and the correct row labels
  if (nrow(vcf3) == 0) {
    DBS.mat.144 <-
      matrix(0, nrow = length(ICAMS::catalog.row.order$DBS144), ncol = 1,
             dimnames = list(ICAMS::catalog.row.order$DBS144, sample.id))
    retval <-
      CheckAndReturnDBSMatrix(vcf = vcf, discarded.variants = discarded.variants,
                              mat78 = DBS.mat.78, mat136 = DBS.mat.136,
                              mat144 = DBS.mat.144,
                              return.annotated.vcf = return.annotated.vcf,
                              sample.id = sample.id)
    return(retval)
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

  CheckAndReturnDBSMatrix(vcf = vcf, discarded.variants = discarded.variants,
                          mat78 = DBS.mat.78, mat136 = DBS.mat.136,
                          mat144 = DBS.mat.144,
                          return.annotated.vcf = return.annotated.vcf,
                          sample.id = sample.id)
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
    lifecycle::deprecate_soft(when = "3.0.0", 
                              what = "StrelkaSBSVCFFilesToCatalogAndPlotToPdf()",
                              details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead')
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

#' \strong{[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka") instead]}
#' Create ID (small insertions and deletions) catalog from Strelka ID VCF files
#' and plot them to PDF 
#'
#' \strong{[Deprecated, use VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka") instead]}
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
    lifecycle::deprecate_soft(when = "3.0.0", 
                              what = "StrelkaIDVCFFilesToCatalogAndPlotToPdf()",
                              details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead')

    list <-
      StrelkaIDVCFFilesToCatalog(files, ref.genome, region, names.of.VCFs,
                                 flag.mismatches, return.annotated.vcfs,
                                 suppress.discarded.variants.warnings)

    if (output.file != "") output.file <- paste0(output.file, ".")

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
#'   \code{\link{BSgenome}} object
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
#' @section ID classification:
#' See \url{https://github.com/steverozen/ICAMS/raw/master/data-raw/PCAWG7_indel_classification_2017_12_08.xlsx}
#' for additional information on ID (small insertions and deletions) mutation
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
    lifecycle::deprecate_soft(when = "3.0.0", 
                              what = "MutectVCFFilesToCatalogAndPlotToPdf()",
                              details = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect")` instead')
    
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
#'   \code{\link{BSgenome}} object
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
#'   If \code{NULL}, all variants are retained. In almost all cases, the default
#'   value of \code{"PASS"} is what the user would want.
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
#'   absolute difference of VAFs for adjacent SBSs is bigger than \code{max.vaf.diff},
#'   then these adjacent SBSs are likely to be "merely" asynchronous single base
#'   mutations, opposed to a simultaneous doublet mutation or variants involving
#'   more than two consecutive bases.
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
#' See \url{https://github.com/steverozen/ICAMS/raw/master/data-raw/PCAWG7_indel_classification_2017_12_08.xlsx}
#' for additional information on ID (small insertions and deletions) mutation
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
  function(files,
           output.dir,
           ref.genome,
           variant.caller = "unknown",
           num.of.cores = 1,
           trans.ranges = NULL,
           region = "unknown",
           names.of.VCFs = NULL,
           tumor.col.names = NA,
           filter.status = "PASS",
           get.vaf.function = NULL,
           ...,
           max.vaf.diff = 0.02,
           base.filename = "",
           return.annotated.vcfs = FALSE,
           suppress.discarded.variants.warnings = TRUE) {
    num.of.cores <- AdjustNumberOfCores(num.of.cores)

    catalogs0 <-
      VCFsToCatalogs(files = files,
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
                     suppress.discarded.variants.warnings =
                       suppress.discarded.variants.warnings)

    # Retrieve the catalog matrix from catalogs0
    catalogs <- catalogs0
    catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL
    if (base.filename != "") base.filename <- paste0(base.filename, ".")

    for (name in names(catalogs)) {
      non.empty.samples <- RetrieveNonEmptySamples(catalogs[[name]])
      # Only plot samples which have mutations for a specific mutation class
      if (!is.null(non.empty.samples)) {
        PlotCatalogToPdf(non.empty.samples,
                         file = file.path(output.dir,
                                          paste0(base.filename, name, ".pdf")))
        if (name == "catSBS192") {
          PlotCatalogToPdf(non.empty.samples,
                           file = file.path(output.dir,
                                            paste0(base.filename, "SBS12.pdf")),
                           plot.SBS12 = TRUE)
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
    out <- paste0(revc(ref.vec[idx]), revc(alt.vec[idx]))
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
    out <- revc(quad[idx])
    stopifnot(all(out %in% ICAMS::catalog.row.order$DBS136))
    quad[idx] <- out
    return(quad)
  }
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
