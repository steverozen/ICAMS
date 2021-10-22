myAnnotateIDVCF<-function (ID.vcf, ref.genome, flag.mismatches = 0, name.of.VCF = NULL,
                           
                           suppress.discarded.variants.warnings = TRUE)
  
{
  
  if (nrow(ID.vcf) == 0) {
    
    return(list(annotated.vcf = ID.vcf))
    
  }
  
  discarded.variants <- ID.vcf[0, ]
  
  if (suppress.discarded.variants.warnings == TRUE) {
    
    retval <- suppressWarnings(ICAMS:::CheckAndRemoveDiscardedVariants(vcf = ID.vcf,
                                                                       
                                                                       name.of.VCF = name.of.VCF))
    
  }
  
  else {
    
    retval <- ICAMS:::CheckAndRemoveDiscardedVariants(vcf = ID.vcf,
                                                      
                                                      name.of.VCF = name.of.VCF)
    
  }
  
  df <- retval$df
  
  discarded.variants <- dplyr::bind_rows(discarded.variants,
                                         
                                         retval$discarded.variants)
  
  ref.genome <- ICAMS:::NormalizeGenomeArg(ref.genome)
  
  idx <- which(nchar(df$REF) == nchar(df$ALT))
  
  if (length(idx) > 0) {
    
    df1 <- df[-idx, ]
    
    df1.to.remove <- df[idx, ]
    
    df1.to.remove$discarded.reason <- "ID variant with same number of bases for REF and ALT alleles"
    
    discarded.variants <- dplyr::bind_rows(discarded.variants,
                                           
                                           df1.to.remove)
    
  }
  
  else {
    
    df1 <- df
    
  }
  
  idx1 <- which(df1$REF == "" | df1$ALT == "")
  
  if (length(idx1) > 0) {
    
    df2 <- df1[-idx1, ]
    
    df2.to.remove <- df1[idx1, ]
    
    df2.to.remove$discarded.reason <- "Variant has empty REF or ALT alleles"
    
    discarded.variants <- dplyr::bind_rows(discarded.variants,
                                           
                                           df2.to.remove)
    
  }
  
  else {
    
    df2 <- df1
    
  }
  
  complex.indels.to.remove <- which(substr(df2$REF, 1, 1) !=
                                      
                                      substr(df2$ALT, 1, 1))
  
  if (length(complex.indels.to.remove) > 0) {
    
    df3 <- df2[-complex.indels.to.remove, ]
    
    df3.to.remove <- df2[complex.indels.to.remove, ]
    
    df3.to.remove$discarded.reason <- "Complex indel"
    
    discarded.variants <- dplyr::bind_rows(discarded.variants,
                                           
                                           df3.to.remove)
    
  }
  
  else {
    
    df3 <- df2
    
  }
  
  stopifnot(substr(df3$REF, 1, 1) == substr(df3$ALT, 1, 1))
  
  var.width <- abs(nchar(df3$ALT) - nchar(df3$REF))
  
  is.del <- nchar(df3$ALT) <= nchar(df3$REF)
  
  var.width.in.genome <- ifelse(is.del, var.width, 0)
  
  df3$seq.context.width <- var.width * 6
  
  
  
  ## minimum context width = 20
  
  df3$seq.context.width[df3$seq.context.width < 20]<-20
  
  
  
  chr.names <- ICAMS:::CheckAndFixChrNames(vcf.df = df3, ref.genome = ref.genome,
                                           
                                           name.of.VCF = name.of.VCF)
  
  Ranges <- GenomicRanges::GRanges(chr.names, IRanges::IRanges(start = df3$POS - df3$seq.context.width,
                                                               
                                                               end = df3$POS + var.width.in.genome + df3$seq.context.width))
  
  df3$seq.context <- BSgenome::getSeq(ref.genome, Ranges, as.character = TRUE)
  
  seq.to.check <- substr(df3$seq.context, df3$seq.context.width +
                           
                           1, df3$seq.context.width + var.width.in.genome + 1)
  
  mismatches <- which(seq.to.check != df3$REF)
  
  if (length(mismatches) > 0) {
    
    df3$seq.to.check <- seq.to.check
    
    df4 <- df3[-mismatches, ]
    
    df4.to.remove <- df3[mismatches, ]
    
    df4.to.remove$discarded.reason <- "ID variant whose REF alleles do not match the extracted sequence from ref.genome"
    
    discarded.variants <- dplyr::bind_rows(discarded.variants,
                                           
                                           df4.to.remove)
    
  }
  
  else {
    
    df4 <- df3
    
  }
  
  if (nrow(discarded.variants) > 0) {
    
    if (suppress.discarded.variants.warnings == TRUE) {
      
      return(list(annotated.vcf = df4, discarded.variants = discarded.variants))
      
    }
    
    else {
      
      warning("\nSome ID variants were discarded, see element discarded.variants",
              
              " in the return value for more details")
      
      return(list(annotated.vcf = df4, discarded.variants = discarded.variants))
      
    }
    
  }
  
  else {
    
    return(list(annotated.vcf = df4))
    
  }
  
}