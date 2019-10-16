dist2TSS <- function(annotated.SBS.vcf) {
  
  df <- annotated.SBS.vcf
  dist2TSSbin <- function(df, pos = "POS", TSS = "trans.start.pos") {
    df <- data.frame(df)
    df$dist2TSS <- abs(df[, pos] - df[, TSS])
    df$distBins <- NULL
    df$distBins[df$dist2TSS <= 1e+05] <- 1
    df$distBins[df$dist2TSS <= 2e+05 & df$dist2TSS > 1e+05] <- 2
    df$distBins[df$dist2TSS <= 3e+05 & df$dist2TSS > 2e+05] <- 3
    df$distBins[df$dist2TSS <= 4e+05 & df$dist2TSS > 3e+05] <- 4
    df$distBins[df$dist2TSS <= 5e+05 & df$dist2TSS > 4e+05] <- 5
    df$distBins[df$dist2TSS <= 6e+05 & df$dist2TSS > 5e+05] <- 6
    df$distBins[df$dist2TSS <= 7e+05 & df$dist2TSS > 6e+05] <- 7
    df$distBins[df$dist2TSS <= 8e+05 & df$dist2TSS > 7e+05] <- 8
    df$distBins[df$dist2TSS <= 9e+05 & df$dist2TSS > 8e+05] <- 9
    df$distBins[df$dist2TSS <= 1e+06 & df$dist2TSS > 9e+05] <- 10
    return(df)
  }
  
  df <- dist2TSSbin(annotated.SBS.vcf)
  df <- subset(df, df$distBins != "NA")
  
  df$mutation <- paste0(df$REF, ">", df$ALT)
  i <- which(df$trans.strand == "-")
  df$mutation[i] <- paste0(ICAMS:::revc(df$REF[i]), ">", ICAMS:::revc(df$ALT[i]))
  
  output <- matrix(data = NA, nrow = 10, ncol = 12)
  rownames(output) <- paste0("group", 1:10)
  colnames(output)<-c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")
  
  i <- as.data.frame(table(df$mutation[df$distBins == 1]))
  ## ensure all mutation classes are represented
  rownames(i)<-i$Var1
  i<-i[colnames(output),]
  i$Var1<-colnames(output)
  i$Freq[is.na(i$Freq)]<-0
  
  colnames(output) <- names(i)
  output[1, ] <- i
  
}