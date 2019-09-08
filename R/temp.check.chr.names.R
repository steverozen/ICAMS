#' Check and, if possible, correct the chromosome names in a VCF \code{data.frame}.
#' 
#' @param vcf.df A VCF as a \code{data.frame}. Check the names in column
#' \code{CHROM}.
#' 
#' @param ref.genome The reference genome with the chromosome names to check
#' \code{vcf.df$CHROM} against; must be a Bioconductor 
#' \code{\link[BSgenome]{BSgenome}}, e.g.
#' \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}.
#' 
#' @return If the \code{vcf.df$CHROM} values are correct or
#' can be corrected, then a vector of chromosome names
#' that can be used as a replacement for \code{vcf.df$CHROM}.
#' If the names in \code{vcf.df$CHROM} cannot be made to
#' be consistent with the chromosome names in in \code{ref.genome},
#' then \code{stop}.
#' 
#' @export

CheckAndFixChrNames <- function(vcf.df, ref.genome) {
  
  names.to.check <- unique(vcf.df$CHROM)

  ref.genome.names <- seqnames(ref.genome)
  
  not.matched <- setdiff(names.to.check, ref.genome.names)
  
  # The names match -- we leave well-enough alone
  if (length(not.matched) == 0) return(vcf.df$CHROM)
  
  chr.index <- grep(pattern = "^chr", names.to.check)
  vcf.has.chr.prefix <- length(chr.index) != 0
  
  chr.index <- grep(pattern = "^chr", ref.genome.names)
  ref.has.chr.prefix <- length(chr.index) != 0

  new.chr.names <- vcf.df$CHROM
  if (ref.has.chr.prefix && !vcf.has.chr.prefix) {
    names.to.check <- paste0("chr", names.to.check)
    new.chr.names <- paste0("chr", new.chr.names)
    not.matched1 <- setdiff(names.to.check, ref.genome.names)
    if (length(not.matched1) == 0) return(new.chr.names)
  }

  organism <- BSgenome::organism(ref.genome)
  
  if (organism == "Homo sapiens") {
    
    # Maybe the problem is that X and Y are encoded as 23 and 24
    if ("chr23" %in% names.to.check) {
      new.chr.names[new.chr.names == "chr23"] <- "chrX"
      names.to.check <- setdiff(names.to.check, "chr23")
      names.to.check <- unique(c(names.to.check, "chrX"))
    }
    if ("chr24" %in% names.to.check) {
      new.chr.names[new.chr.names == "chr24"] <- "chrY"
      names.to.check <- setdiff(names.to.check, "chr24")
      names.to.check <- unique(c(names.to.check, "chrY"))
    }
    
    # Maybe the problem is that X and Y are encoded as 23 and 24
    if ("23" %in% names.to.check) {
      new.chr.names[new.chr.names == "23"] <- "X"
      names.to.check <- setdiff(names.to.check, "23")
      names.to.check <- unique(c(names.to.check, "X"))
    }
    if ("24" %in% names.to.check) {
      new.chr.names[new.chr.names == "24"] <- "Y"
      names.to.check <- setdiff(names.to.check, "24")
      names.to.check <- unique(c(names.to.check, "Y"))
    }
  }
  
  not.matched3 <- setdiff(names.to.check, ref.genome.names)
  if (length(not.matched3) == 0) return(new.chr.names)
  
  stop("Chromosome names in input not in ref.genome for ",
       organism, ": ", 
       # We report the _original_ list of not matched names
       paste(not.matched, collapse = " "))
}
