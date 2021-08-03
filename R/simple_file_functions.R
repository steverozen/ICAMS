#' Generate reconstructed VCFs from indel (small insertions and deletions) simple files
#'
#' @param files Character vector of file paths to the indel simple files.
#' 
#' @param output.dir The directory where the reconstructed VCFs will be saved.
#' 
#' @param num.parallel.files The (maximum) number of files to run in
#'   parallel. On Microsoft Windows machines it is silently changed to 1. Each
#'   file in turn can require multiple cores, as governed by
#'   \code{mc.cores.per.file}.
#'
#' @param mc.cores.per.file The maximum number of cores to use for each
#'   file. On Microsoft Windows machines it is silently changed to 1.
#'   
#' @keywords internal
#'
GenerateVCFsFromIndelSimpleFiles <- 
  function(files, output.dir, num.parallel.files = 1, mc.cores.per.file = 1) {
    time.used <- system.time(
      retval <- parallel::mclapply(files, FUN = function(x) {
        GenerateVCFsFromIndelSimpleFile1(file = x,
                                         output.dir = output.dir,
                                         max.mc.cores = mc.cores.per.file)
      }, mc.cores = AdjustNumberOfCores(num.parallel.files))
    )
    return(time.used)
  }

#' Generate reconstructed VCFs from indel (small insertions and deletions) simple file
#'
#' @param file The name/path of the simple indel file, or a complete URL.
#' 
#' @param output.dir The directory where the reconstructed VCFs will be saved.
#' 
#' @param max.mc.cores The maximum number of cores to use. On Microsoft Windows
#'   machines it is silently changed to 1.
#'
#' @keywords internal
GenerateVCFsFromIndelSimpleFile1 <- function(file, output.dir, max.mc.cores = 1) {
  
  dt <- data.table::fread(file)
  
  simple.col.names <- c("CancerType",
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
  
  colnames(dt) <- simple.col.names
  dt$ID <- dt$QUAL <- dt$INFO <- "."
  dt$FILTER <- "PASS"
  dt$FORMAT <- "GT:AD:AF:DP"
  
  vcf.col.names <- c("CHROM",
                     "POS",
                     "ID",
                     "REF",
                     "ALT",
                     "QUAL",
                     "FILTER",
                     "INFO",
                     "FORMAT",
                     "SampleID",
                     "DataCollection",
                     "POSEnd",
                     "CancerType",
                     "MutationType",
                     "GenomeVersion",
                     "MoreInfo")
  data.table::setcolorder(dt, neworder = vcf.col.names)
  
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  out <- parallel::mclapply(unique(dt$SampleID), FUN = function(x) {
    sample.id <- x
    genome.version <- unique(dt$GenomeVersion)
    if (genome.version %in% c("GRCh37", "hg19", "hs37d5")) {
      genome.version <- "GRCh37"
    } else if (genome.version %in% c("GRCh38", "hg38")) {
      genome.version <- "GRCh38"
    }
    
    ref.genome <- NormalizeGenomeArg(genome.version)
    dt1 <- dplyr::filter(dt, SampleID == sample.id)
    
    # Canonicalize the simple file INDEL to VCF format INDEL
    # Check if the format of sequence names in dt1 and genome are the same
    chr.names <- CheckAndFixChrNames(vcf.df = dt1, ref.genome = ref.genome)
    
    ref.genome.seqlengths <- ref.genome@seqinfo@seqlengths
    names(ref.genome.seqlengths) <- ref.genome@seqinfo@seqnames 
    
    # Change the format of small deletions
    del.idx <- which(dt1$ALT == "-") # Get all the indices of small dels
    if (length(del.idx) > 0) {
      pos.del <- dt1[del.idx, ]$POS - 1
      chrom.del <- chr.names[del.idx]
      
      # Check whether there are variants whose chrom.del and pos.del
      # falls outside the ref.genome
      ret <- sapply(1:length(chrom.del), FUN = function(x) {
        chrom <- chrom.del[x]
        pos <- pos.del[x]
        return(pos > ref.genome.seqlengths[chrom])
      })
      
      idx.to.remove <- which(ret)
      
      if (length(idx.to.remove) > 0) {
        idx.to.remove.in.dt <- del.idx[idx.to.remove]
        chrom.del <- chrom.del[-idx.to.remove]
        pos.del <- pos.del[-idx.to.remove]
        del.idx <- del.idx[-idx.to.remove]
      }
      
      ref.del <- BSgenome::getSeq(x = ref.genome,
                                  names = chrom.del, 
                                  start = pos.del, 
                                  width = 1,
                                  as.character = TRUE)
      dt1[del.idx, ]$POS <- pos.del
      dt1[del.idx, ]$REF <- paste0(ref.del, dt1[del.idx, ]$REF)
      dt1[del.idx, ]$ALT <- ref.del
      if (length(idx.to.remove) > 0) {
        message("Remove deletion variants whose CHROM and POS are beyond the boundaries of ref.genome")
        print(dt1[idx.to.remove.in.dt, ])
        dt1 <- dt1[-idx.to.remove.in.dt, ]
        chr.names <- chr.names[-idx.to.remove.in.dt]
      } else {
        dt1 <- dt1
      }
      
    }
    
    # Find the indices of small insertions
    ins.idx <- which(dt1$REF == "-")
    if (length(ins.idx) > 0) {
      pos.ins <- dt1[ins.idx, ]$POS - 1
      chrom.ins <- chr.names[ins.idx]
      
      # Check whether there are variants whose chrom.del and pos.del
      # falls outside the ref.genome
      ret <- sapply(1:length(chrom.ins), FUN = function(x) {
        chrom <- chrom.ins[x]
        pos <- pos.ins[x]
        return(pos > ref.genome.seqlengths[chrom])
      })
      
      idx.to.remove <- which(ret)
      
      if (length(idx.to.remove) > 0) {
        idx.to.remove.in.dt <- ins.idx[idx.to.remove]
        chrom.ins <- chrom.ins[-idx.to.remove]
        pos.ins <- pos.ins[-idx.to.remove]
        ins.idx <- ins.idx[-idx.to.remove]
      }
      
      ref.ins <- BSgenome::getSeq(x = ref.genome,
                                  names = chrom.ins, start = pos.ins, width = 1,
                                  as.character = TRUE)
      dt1[ins.idx, ]$POS <- pos.ins
      dt1[ins.idx, ]$REF <- ref.ins
      dt1[ins.idx, ]$ALT <- paste0(ref.ins, dt1[ins.idx, ]$ALT)
      
      if (length(idx.to.remove) > 0) {
        message("Remove insertion variants whose CHROM and POS are beyond the boundaries of ref.genome")
        print(dt1[idx.to.remove.in.dt, ])
        dt1 <- dt1[-idx.to.remove.in.dt, ]
      } else {
        dt1 <- dt1
      }
    }
    
    # Change the first column name from "CHROM" to "#CHROM"
    colnames(dt1)[1] <- "#CHROM"
    vcf.name <- paste0(unique(dt$CancerType), ".", sample.id, ".", 
                       genome.version, ".indel.vcf")
    vcf.output.file.path <- file.path(output.dir, vcf.name)
    data.table::fwrite(list("##VCF reconstructed from simple file"), 
                       file = vcf.output.file.path)
    
    data.table::fwrite(dt1, file = vcf.output.file.path, append = TRUE,
                       col.names = TRUE, sep = "\t")
  }, mc.cores = AdjustNumberOfCores(max.mc.cores)) # end of mclapply for sample id
}

#' @keywords internal
GeneratePCAWGAliquotID <- function (file){ 
  dt <- data.table::fread(file)
  aliquot.ids.info <- dt$tumor_wgs_aliquot_id
  SP.ids.info <- dt$tumor_wgs_icgc_specimen_id
  
  aliquot.ids.list <- stringi::stri_split_fixed(aliquot.ids.info, ",")
  SP.ids.list <- stringi::stri_split_fixed(SP.ids.info, ",")
  
  aliquot.ids <- unlist(aliquot.ids.list)
  SP.ids <- unlist(SP.ids.list)
  
  names(aliquot.ids) <- SP.ids
  return(aliquot.ids)
}

#' @keywords internal
GeneratePCAWGAliquotID2 <- function(file) {
  dt <- data.table::fread(file)
  aliquot.ids.info <- dt$aliquot_id
  SP.ids.info <- dt$icgc_specimen_id
  
  aliquot.ids.list <- stringi::stri_split_fixed(aliquot.ids.info, ",")
  SP.ids.list <- stringi::stri_split_fixed(SP.ids.info, ",")
  
  aliquot.ids <- unlist(aliquot.ids.list)
  SP.ids <- unlist(SP.ids.list)
  
  names(aliquot.ids) <- SP.ids
  return(aliquot.ids)
}