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

#' ReadSimpleAsVCF
#'
#' @param path TODO
#' @param verbose TODO
#'
#' @return TODO
#' @export
ReadSimpleAsVCF <- function(path, verbose = 0) {
  dt <- fread(
    path,
    col.names = simple.col.order
  )

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
  dt[ , VAF := as.numeric(vaf)]

  dt <- dt[ , ..vaf.col.order ]

  list.of.vcf <- split(as.data.frame(dt), dt$SampleID)

  return(list.of.vcf)
}

#' ReadIDSimpleAsVCF
#'
#' @param path TODO
#' @param verbose TODO
#'
#' @return TODO
#' @export
ReadIDSimpleAsVCF <- function(path, verbose = 0) {
  dt <- fread(
    path,
    col.names = simple.col.order
  )

  if (nrow(dt) == 0) {
    # File is empty
    return(list())
  }

  if (length(unique(dt$GenomeVersion)) != 1) {
    cat(paste(unique(dt$GenomeVersion), collapse="|"), "\n")
    stop()
  }

  df <- as.data.frame(dt)

  stopifnot(df$REF != "-" | df$MutationType == 'INS')
  stopifnot(df$ALT != "-" | df$MutationTYpe == 'DEL')
  df[df$REF == "-", "REF"] <- ""
  df[df$ALT == "-", "ALT"] <- ""

  list.of.vcf <- split(df, df$SampleID)

  return(list.of.vcf)
}

#' CreateSampleId
#'
#' Get catalog sorting function from msigtools
#'
#' @param filename TODO
#' @param aliquot.id TODO
#'
#' @return TODO
#' @export
CreateSampleId <- function(filename, aliquot.id) {
  filename <- sub("\\..*", "", filename, perl = TRUE)
  return(paste0(filename, "::", aliquot.id))
}

#' OneSBSVCFAllCatalogs
#'
#' TODO(steve) merge this in where needed
#'
#' @param vcf TODO
#' @param genome TODO
#' @param sample.id TODO
#' @param verbose TODO
#'
#' @return TODO
#' @export
OneSBSVCFAllCatalogs <- function(vcf,
                                 genome,
                                 sample.id = "count",
                                 verbose = 0) {
  three.vcfs.df <- SplitSNSVCF(vcf)
  SNS <- three.vcfs.df$SNS.vcf
  DNS <- three.vcfs.df$DNS.vcf
  if (verbose > 0) {
    cat("length SNS and DNS vcfs", nrow(SNS), nrow(DNS), "\n")
  }

  SNS <- AddSequence(SNS, seq = genome)
  CheckSeqContextInVCF(SNS, "seq.21context")
  SNS <- AddTranscript(SNS, .trans.ranges)
  SNS.cat <- CreateOneColSNSCatalog(SNS, colname = sample.id)
  rm(SNS)

  DNS <- AddSequence(DNS, seq = genome)
  DNS <- AddTranscript(DNS, .trans.ranges)
  CheckSeqContextInVCF(DNS, "seq.21context")
  DNS.cat <- CreateOneColDNSCatalog(DNS, colname = sample.id)
  rm(DNS)

  ret <- SNS.cat
  SNS.cat[["catDNS78"]] <- DNS.cat
  return(ret)
}

#' DirectoryOfSBSSimple2ListOfVCF
#'
#' @param path.to.dir TODO
#' @param cts TODO
#' @param verbose TODO
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @return TODO
#' @export
DirectoryOfSBSSimple2ListOfVCF <-
  function(path.to.dir, cts = empty.cats, verbose = 3) {
    genome <- BSgenome.Hsapiens.1000genomes.hs37d5

    in.files <- list.files(path.to.dir, pattern = "*.simple*")
    out.list <- list()
    for (fl in in.files) {
      filepath <- paste0(path.to.dir, fl)
      cat("Processing file", filepath, "\n")
      out.list <- RCurl::merge.list(out.list, ReadSimpleAsVCF(filepath))
    }
    return(out.list)
  }

#' DirectoryOfSBSSimple2Catalog
#'
#' @param path.to.dir TODO
#' @param cts TODO
#' @param verbose TODO
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @return TODO
#' @export
DirectoryOfSBSSimple2Catalog <-
  function(path.to.dir, cts = empty.cats, verbose = 3) {
    genome <- BSgenome.Hsapiens.1000genomes.hs37d5

    in.files <- list.files(path.to.dir, pattern = "*.simple*")
    for (fl in in.files) {
      filepath <- paste0(path.to.dir, fl)
      cat("Processing file", filepath, "\n")
      vcf.list <- ReadSimpleAsVCF(filepath)
      for (sample.id in names(vcf.list)) {
        id2 <- CreateSampleId(fl, sample.id)
        if (verbose > 0 ) cat("Sample", id2, "\n")

        three.vcfs.df <- SplitSNSVCF(vcf.list[[sample.id]])
        SNS <- three.vcfs.df$SNS.vcf
        DNS <- three.vcfs.df$DNS.vcf
        cat("length SNS and DNS vcfs", nrow(SNS), nrow(DNS), "\n")

        SNS <- AddSequence(SNS, seq = genome)
        CheckSeqContextInVCF(SNS, "seq.21context")
        SNS <- AddTranscript(SNS, .trans.ranges)
        SNS.cat <- CreateOneColSNSCatalog(SNS, colname = id2)
        rm(SNS)

        DNS <- AddSequence(DNS, seq = genome)
        DNS <- AddTranscript(DNS, .trans.ranges)
        CheckSeqContextInVCF(DNS, "seq.21context")
        DNS.cat <- CreateOneColDNSCatalog(DNS, colname = id2)
        rm(DNS)

        cts$cat96 <- cbind(cts$cat96, SNS.cat$cat96)
        cts$cat192 <- cbind(cts$cat192, SNS.cat$cat192)
        cts$cat1536 <- cbind(cts$cat1536, SNS.cat$cat1536)
        cts$catDNS78 <- cbind(cts$catDNS78, DNS.cat)
      }
    }
    return(cts)
  }

#' InitAliquotSample
#'
#' In InitAliquotSample, note the hack to add one more mapping. One aliquot ID
#' is not in PCAWG as per ZHANG Junjun's email
#' https://mail.google.com/mail/u/0/#search/Junjun.Zhang%40oicr.on.ca++pcawg_wgs/FMfcgxmTmlJdCfGhxbDHpsJCwknCNTDw
#' [1] "CMDI-UK::f92cf0a2-0172-a1cc-e040-11ac0c486b1a" >
#' setdiff(StripCaType(colnames(ofcat96)), StripCaType(ycol)) [1]
#' "CMDI-UK::SP116883" was the associated sample id. See
#' https://mail.google.com/mail/u/0/#search/SP116883/FMfcgxmTmlJdCfGhxbDHpsJCwknCNTDw
#'
#' @return TODO
#' @export
InitAliquotSample <- function() {
  if (!exists(".ali2samp")) {
    t1 <-
      fread("pcawg_wgs_tumour_specimen-2017-03-02.tsv.gz",
            header = FALSE,
            col.names = c("AliquotID", "SampleID", "TumourInfo", "GreyOrWhite"))
    L <-
      list(t1,
           list("f92cf0a2-0172-a1cc-e040-11ac0c486b1a", "SP116883", "CMDI-UK", "Excluded"))
    #.ali2samp <<- copy(rbindlist(L))
    #.samp2ali <<- copy(.ali2samp)
    setkey(.ali2samp, AliquotID)
    setkey(.samp2ali, SampleID)
  }
}

#' AliquotID2SampleID
#'
#' @param aliquot.id TODO
#'
#' @return TODO
#' @export
AliquotID2SampleID <- function(aliquot.id) {
  InitAliquotSample()
  rows <- .ali2samp[aliquot.id]
  stopifnot(rows$AliquotID == aliquot.id)
  return(rows[ , SampleID])
}

#' SampleID2AliquotID
#'
#' @param sample.id TODO
#'
#' @return TODO
#' @export
SampleID2AliquotID <- function(sample.id) {
  InitAliquotSample()
  rows <- .samp2ali[sample.id]
  stopifnot(rows$SampleID == sample.id)
  return(rows[ , AliquotID])
}

#' StripCaType
#'
#' @param x TODO
#'
#' @return TODO
#' @export
StripCaType <- function(x) {
  return(sub("^.*::","", x, perl = TRUE))
}

#' SplitCaTypeAndID
#'
#' @param x TODO
#'
#' @return TODO
#' @export
SplitCaTypeAndID <- function(x) {
  x <- tstrsplit(x, "::", fixed = TRUE)
  names(x) <- c("CAType", "SampleID")
  return(as.data.frame(x, stringsAsFactors = FALSE))
}

#' CaTypeAliquot2CaTypeSampleID
#'
#' @param ca.type.and.aliquot.id TODO
#'
#' @return TODO
#' @export
CaTypeAliquot2CaTypeSampleID <- function(ca.type.and.aliquot.id) {
  out.list <- strsplit(ca.type.and.aliquot.id, split="::", fixed = TRUE)
  out2 <- unlist(out.list)
  out3 <- matrix(out2, ncol = 2, byrow = T)
  out4 <- paste0(out3[ ,1], "::", AliquotID2SampleID(out3[ ,2]))
  return(out4)
}

#' WriteMinimalVCFList
#'
#' Each VCF in list.of.VCFs needs to have columns
#' CHROM, POS, REF, ALT, VAF, STRAND ??????.
#' They will be extracted and re-ordered if necessary
#' to be consistent. The header line will begin with
#' #CHROM
#'
#' @param list.of.VCFs TODO
#' @param path TODO
#'
#' @return TODO
#' @export
WriteMinimalVCFList <- function(list.of.VCFs, path) {
}
