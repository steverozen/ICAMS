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
  
  output[1, ] <- i$Freq
  for (class in colnames(output)) {
    for (row in 2:10) {
      output[row, class] <- sum(df$mutation[df$distBins == row] == class)
    }
  }
  
  output <- as.matrix(output)
  output <- cbind(output, Sense = rowSums(output[, c(4:6, 10:12)]))
  output <- cbind(output, antiSense = rowSums(output[, c(1:3, 7:9)]))
  
  return(output)
}

Plotdist2TSS <- function(output, main="") {
  output <- data.frame(output)
  output2 <- output[,c("antiSense","Sense")]
  rownames(output2) <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", 
                         "8-9", "9-10")
  colnames(output2) <- c("antiSense", "Sense")
  output2 <- t(output2)
  
  ## logistic regression. Should regress to AS/S ratio instead of trans.strand? Not working as 
  ## this is probably regressing the wrong parameters.
  
  ## df$newStrand <- NA
  ## df$newStrand[df$trans.strand == '+'] <- 'sense'
  ## df$newStrand[df$trans.strand == '-'] <- 'antisense'
  ## fit<-glm(as.factor(df$newStrand) ~ df$dist2TSS, family = 'binomial')
  ## pValue <- capture.output(cat('p =', round(coef(summary(fit))[2])))
  
  barplot(output2, beside = TRUE, ylim = c(0, max(output2)*1.2), 
                           col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.9, 
          ylab="# of mutations",main=main) 
  legend("topright", c("antiSense", "Sense"), bty="n", 
                                        fill=c("grey80", "grey32"),cex=0.8)
  # text (27.4, max(output2)*0.74, pValue, cex=0.8)
}

## Plotting individual mutation classes. Mutation class can be specified in this function.
## Use PlotStrandBiasTSSAll to plot all mutations classes.

PlotStrandBiasTSS <- function(matrix, class) {
  
  stopifnot(class %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  dt <- data.table(matrix)
  revc1 <- function(string) {
    return(paste0(ICAMS:::revc(substr(string, 1, 1)), '>', 
                  ICAMS:::revc(substr(string, 3, 3))))
  }
  output <- dt[, c(revc1(class), class), with = FALSE]
  colnames(output) <- c("antiSense", "Sense")
  barplot(t(output), beside=TRUE, ylim = c(0, max(output) + 10), 
          col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.4, ylab="# of mutations",
          names.arg=paste0('group', 1:10),main=paste0(class, "mutations"))
  legend("topright", c("antiSense", "Sense"), bty="n", 
         fill=c("grey80", "grey32"))
  return(output)
}

PlotStrandBiasTSSAll <- function(matrix) {
  
  dt <- data.table(matrix)
  revc1 <- function(string) {
    return(paste0(ICAMS:::revc(substr(string, 1, 1)), '>', 
                  ICAMS:::revc(substr(string, 3, 3))))
  }
  
  for (class in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) {
    output <- dt[, c(revc1(class), class), with = FALSE]
    barplot(t(output), beside=TRUE, ylim = c(0, max(output) *1.2), 
            col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.4, ylab="# of mutations",
            names.arg=paste0('group', 1:10))
    legend("topright", c("antiSense", "Sense"), bty="n", 
           fill=c("grey80", "grey32"))
    title(c(class, "mutations"))
  }
}

#' Plot distance to TSS to a PDF file. Unfinished, not working yet.
#'
#' @param file The name of the PDF file to be produced.

## dist2PDF <- function(file) {
##   grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
##   (barplot(output2, beside = TRUE, ylim = c(0, max(output) + 30), legend = TRUE))
##   grDevices::dev.off()
##   invisible(TRUE)
##   
## }
