### ICAMS related functions needed only for processing PCAWG7 data 


# source("ICAMS_dependencies.R")
# source("ICAMS_function.R")
# source("testDBS.R")


# Example

library(data.table)

d1 <- fread("data-raw/draft_code/sup-table-5.csv")

# assume microbiopsy_id is the ID we want




#' Read a PCAWG7 "simple" file as if it were a VCF.
#' 
#' @param path Path to the simple file.
#' 
#' @return A list of VCFs as in-memory data.frame's, one for
#' each sample in the input file.
#' 
#' @importFrom data.table fread
ReadSimpleAsVCF <- function(path) {
  
  simple.col.order <- 
    c("CancerType",
      "SampleID",
      "DataCollection",
      "GenomeVersion",
      "MutationType",
      "CHROM",
      "POS",
      "POSEnd",
      "REF",
      "ALT",
      "MoreInfo")
  
  vaf.col.order <- 
    c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "SampleID",
      "VAF",
      "DataCollection",
      "POSEnd",
      "CancerType",
      "MutationType",
      "GenomeVersion",
      "MoreInfo"
    )
  
  dt <- fread(
    path,
    col.names = simple.col.order
  )
  
  dt <- as.data.frame(dt)
  
  if (nrow(dt) == 0) {
    # File is empty
    return(list())
  }
  
  if (length(unique(dt$GenomeVersion)) != 1) {
    cat(paste(unique(dt$GenomeVersion), collapse="|"), "\n")
    stop()
  }
  
  GetVaf <- function(x) {
    if (length(x) == 0) return(NA)
    v <- grep(pattern = "vaf=", x = x, fixed = TRUE, value = TRUE)
    if (length(v) != 1) {
      if (!grepl("^dbsnp", x, perl = TRUE)) {
        cat("funny record, no vaf, returning NA\n", x, "\n")
      }
      return(NA)
    }
    return(sub("vaf=", "", x = v,fixed = TRUE))
  }
  
  vaf.plus <- strsplit(dt$MoreInfo, split = ";", fixed = TRUE)
  vaf <- sapply(X = vaf.plus, FUN = GetVaf)
  dt$VAF <- as.numeric(vaf)
  
  dt <- dt[ , vaf.col.order ] # Re-order the columns
  
  list.of.vcf <- split(dt, dt$SampleID)
  
  return(list.of.vcf)
}

#' Read a PCAWG7 "simple" file as if it were a VCF.
#' 
#' @param path Path to the simple file.
#' 
#' @param use.sample.id If TRUE translate from aliquot identifier
#' sample identifier.
#' 
#' @return Three lists, each with on element for each
#' sample in the VCF in \code{path}.
#' 
#' \enumerate{
#'
#'  \item \code{SNS.vcfs} A list of VCFs as a data.frame.
#'  
#'  \item \code{DNS.vcfs} A list of VCFs as data.frame.
#'  
#'  \item \code{ThreePlus} A list of VCF-like data frames
#'  that contain left over rows that represent coordinate mutations
#'  from > 2 adjacent nucleoides, e.g. ATC > CAG.
#' 
#' }
#' 
#' @importFrom ICAMS SplitListOfStrelkaSNSVCFs

ReadAndSplitPCAWG7SNS <- function(path, use.sample.id = TRUE) {
  
  SBS.vcfs <- ReadSimpleAsVCF("BTCA-SG.simple.gz")
  
  split.vcfs <- ICAMS:::SplitListOfStrelkaSNSVCFs(SBS.vcfs)
  
  #TODO(Steve): remove extra columns (?)
  
  if (use.sample.id) {
    names(split.vcfs$SNS.vcfs)  <- AliquotID2SampleID(names(SBS.vcfs))
    names(split.vcfs$DNS.vcfs)  <- AliquotID2SampleID(names(SBS.vcfs))
    names(split.vcfs$ThreePlus) <- AliquotID2SampleID(names(SBS.vcfs))
  } else {
    names(split.vcfs$SNS.vcfs)  <- names(SBS.vcfs)
    names(split.vcfs$DNS.vcfs)  <- names(SBS.vcfs)
    names(split.vcfs$ThreePlus) <- names(SBS.vcfs)
  }
  return(split.vcfs) 
}


#' Create SNS and DNS catalogs from Strelka SNS VCF files
#'
#' Create 3 SNS catalogs (96, 192, 1536) and 3 DNS catalogs (78, 136, 144)
#' from the Strelka SNS VCFs specified by vector.of.file.paths
#'
#' This function calls \code{\link{VCFsToSNSCatalogs}} and
#' \code{\link{VCFsToDNSCatalogs}}
#' 
#' @param file The name of a single "simple" file.
#' 
#' @param genome  Name of a particular reference genome
#'   (without quotations marks).
#'   
#' @param trans.ranges A data.table which contains
#'  transcript ranges for \code{genome}.
#'   
#' @return  A list of 3 SNS catalogs (one each for 96, 192, and 1536)
#'   and 3 DNS catalogs (one each for 78, 136, and 144)
#'   
#' @importFrom ICAMS VCFsToSNSCatalogs VCFsToDNSCatalogs
#' @export
PCAWG7SimpleFilesToCatalogs <-
  function(file, genome, trans.ranges) {
    split.vcfs <- ReadAndSplitPCAWG7SNS(file)
    return(
      c(ICAMS:::VCFsToSNSCatalogs(
        split.vcfs$SNS.vcfs, genome, trans.ranges),
        ICAMS:::VCFsToDNSCatalogs(
          split.vcfs$DNS.vcfs, genome, trans.ranges)))
  }

#' Test various functionality
#' 
TestSampleAliquot <- function() {
  # These are BOCA greylist
  x1 <- 
    c("f8467ec8-2d61-ba21-e040-11ac0c483584",
      "f856fa85-fdb8-c0b0-e040-11ac0d480b4e", 
      "f8696c79-b165-92a6-e040-11ac0c4804bf")
  x2 <- AliquotID2SampleID(x1)
  x3 <- SampleID2AliquotID(x2)
  stopifnot(x1 == x3)
  
  CaTypeAliquot2CaTypeSampleID <- function(ca.type.and.aliquot.id) {
    out.list <- strsplit(ca.type.and.aliquot.id, split="::", fixed = TRUE)
    out2 <- unlist(out.list)
    out3 <- matrix(out2, ncol = 2, byrow = T)
    out4 <- paste0(out3[ ,1], "::", AliquotID2SampleID(out3[ ,2]))
    return(out4)
  }

  x4 <- paste0(c("FOO", "BAR", "BAM"), "::", x1)
  x5 <- CaTypeAliquot2CaTypeSampleID(x4)
  stopifnot(x5 == c("FOO::SP116515", "BAR::SP116499", "BAM::SP116505"))
  cat("ok\n")
}
