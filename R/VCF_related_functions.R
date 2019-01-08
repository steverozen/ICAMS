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

