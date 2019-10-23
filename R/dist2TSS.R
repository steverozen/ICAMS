dist2TSS <- function(annotated.SBS.vcf, plot.type) {
  
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
  df$mutation[i] <- paste0(revc(df$REF[i]), ">", revc(df$ALT[i]))
  
  revc1 <- function(string) {
    start <- revc(substr(string, 1, 1))
    end <- revc(substr(string, 3, 3))
    return(paste0(start, ">", end))
  }
  df <- df[df$mutation %in% c(plot.type, revc1(plot.type)), ]
  
  ## logistic regression. Should regress to AS/S ratio instead of trans.strand? Not working as 
  ## this is probably regressing the wrong parameters.

  type <- c('C>A', 'C>T', 'C>G', "T>A", "T>C", 'T>G')
  df$class <- 0
  df$class[df$mutation %in% type] <- 1
  
  fit<-glm(as.factor(df$class) ~ df$dist2TSS, family = 'binomial')
  pValue <- coef(summary(fit))[2, 4]
  
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
  
  return(list(matrix = output, pvalue = pValue))
}

Plotdist2TSS <- function(output, plot.type) {
  output1 <- data.frame(output$matrix)
  output2 <- output1[,c("antiSense","Sense")]
  rownames(output2) <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", 
                         "8-9", "9-10")
  colnames(output2) <- c("antiSense", "Sense")
  output2 <- t(output2)
  
  barplot(output2, beside = TRUE, ylim = c(0, max(output2)*1.2), 
                           col = c("grey80", "grey32"), cex.axis = 0.8, cex.names = 0.9, 
          ylab="# of mutations", main=plot.type) 
  legend("topright", c("antiSense", "Sense"), bty="n", 
                                        fill=c("grey80", "grey32"),cex=0.8)
  text (28.9, max(output2)*0.74, paste0('p value =', 
                                        formatC(output$pvalue, format = "e", digits = 2)), 
                                        cex=0.8, xpd = NA)
}

#' Plot transcriptional strand bias with respect to distance to transcription
#' start site
#' 
#' @param annotated.SBS.vcf An annotated SBS VCF as a data.table which contains
#'   sequence context and transcript information. Please refer to
#'   \code{\link{AnnotateSBSVCF}} for more details.
#'
#' @param plot.type A character string indicating one mutation type to be
#'   plotted. It should be one of "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#' 
#' @importFrom stats coef
#' 
#' @return \code{invisible(TRUE)}
#' 
#' @export
#'
#' @examples 
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)
#'   PlotTransBiasDist2TSS(annotated.SBS.vcf, plot.type = c("C>A", "C>T"))
#' }
PlotTransBiasDist2TSS <- function (annotated.SBS.vcf, plot.type){
  output <- dist2TSS(annotated.SBS.vcf, plot.type)
  Plotdist2TSS(output, plot.type)
  return(invisible(TRUE))
}

#' Plot transcriptional strand bias with respect to distance to transcription
#' start site to a PDF file
#'
#' @param annotated.SBS.vcf An annotated SBS VCF as a data.table which contains
#'   sequence context and transcript information. Please refer to
#'   \code{\link{AnnotateSBSVCF}} for more details.
#'
#' @param plot.type A vector of character indicating types to be plotted. It
#'   should be within "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'
#' @param file The name of output file
#'    
#' @importFrom stats glm
#' 
#' @return \code{invisible(TRUE)}
#' 
#' @export
#'
#' @examples
#' file <- c(system.file("extdata",
#'                       "Strelka.SBS.GRCh37.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)
#'   PlotTransBiasDist2TSSToPDF(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                              plot.type = c("C>A","C>G","C>T","T>A","T>C"),
#'                              file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasDist2TSSToPDF <- function(annotated.SBS.vcf, plot.type, file) {
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  opar <- par(mfrow = c(4, 2), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  on.exit(par(opar))
  num <- length(plot.type)
  type <- plot.type
  for (i in 1:num) {
    PlotTransBiasDist2TSS(annotated.SBS.vcf = annotated.SBS.vcf, 
                          plot.type = plot.type[i])
  }
  
  grDevices::dev.off()
  invisible(TRUE)
  
}
