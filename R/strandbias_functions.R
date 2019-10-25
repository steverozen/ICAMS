StrandBiasAsExpressionLevel <- 
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col, 
           TPM.col, num.of.bins) {
    
    TPM <- data.table(expression.level[, c(Ensembl.gene.ID.col, TPM.col)])
    names(TPM) <- c("trans.Ensembl.gene.ID", "TPM")
    
    # Delete genes expressed on both strands and NA Ensembl.gene.ID
    vcf <- annotated.SBS.vcf[-which(annotated.SBS.vcf$bothstrand == TRUE), ]
    vcf <- vcf[-which(is.na(vcf$trans.Ensembl.gene.ID) == TRUE), ]
    df <- merge(vcf, TPM, by = "trans.Ensembl.gene.ID", all.x = TRUE)
    df <- df[-which(is.na(df$TPM) == TRUE), ]
    
    # One SBS mutation can be represented by more than 1 row in df if the
    # mutation position falls into the range of multiple transcripts on the same
    # strand. We need to sum up the TPM values of multiple transcripts for one particular
    # mutation.
    df1 <- df[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                  seq.21bases = seq.21bases[1],
                  TPM = sum(TPM)), by = .(CHROM, POS)]
    
    df1$Exp_Level <- NA
    cutoffs <- 
      stats::quantile(df1$TPM, c(1:(num.of.bins - 1)/num.of.bins), na.rm = T)
    df1[, "Exp_Level"] <- 
      cut(df1$TPM, breaks = c(-Inf, cutoffs, Inf), labels = (1:num.of.bins))
    
    # Orient the ref.context and the var.context column
    df1$mutation <- paste0(substr(df1$seq.21bases, 10, 12), df1$ALT)
    df1[trans.strand == "-", `:=`(mutation, RevcSBS96(mutation))]
    df1$mutation <- 
      paste0((substr(df1$mutation, 2, 2)), ">", substr(df1$mutation, 4, 4))
    
    # Create a table for logistic regression
    dt <- df1[, c("mutation", "TPM")]
    type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    dt <- dt[mutation %in% type, `:=`(class, 1)]
    dt <- dt[!mutation %in% type, `:=`(class, 0)]
    logit.model <- stats::glm(class ~ TPM, family = binomial, data = dt)
    p.value <- summary(logit.model)$coefficients[2, 4]
    
    # Plot transcriptional strand bias as a function of gene expression
    result <- matrix(data = 0, nrow = num.of.bins, ncol = 12)
    rownames(result) <- c(1:num.of.bins)
    
    mutation.type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                       "G>T", "G>C", "G>A", "A>T", "A>G", "A>C")
    colnames(result) <- mutation.type
    for (x in 1:num.of.bins) {
      for (j in 1:12) {
        type <- mutation.type[j]
        result[x, type] <- nrow(df1[mutation == type & Exp_Level == x, ])
      }
    }
    
    return(list(plotmatrix = result, logit.df = dt, p.value = p.value))
  }


PlotTransBiasExp1 <- function(list, type, num.of.bins, ymax = NULL) {
  result <- list$plotmatrix
  foo1 <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (!type %in% foo1) {
    stop("please input: 'C>A','C>G','C>T','T>A','T>C','T>G'")
  }
  plotCombo <- cbind(Target = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), 
                     rev = c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C"), 
                     colLight = c("#4040FF",  "black", "#FF4040", 
                                  "#838383", "#40FF40", "#FF667F"), 
                     colDark = c("#00CCFF", "grey40", "#FF99CC", 
                                 "#C6C6C6", "#B0FFB0", "#FCB9C9"))
  i <- which(plotCombo[, "Target"] == type)
  tmp <- t(result[, c(plotCombo[i, 1], plotCombo[i, 2])])
  colnames(tmp) <- c(1:num.of.bins)
  if (is.null(ymax)) {
    ymax <- max(tmp)
  }
  bp <- barplot(tmp, beside = T, main = plotCombo[i, "Target"], 
                ylim = c(0, 1.5 * ymax), col = plotCombo[i, 3:4])
  legend("topright", legend = c("Untranscribed", "Transcribed"), 
         fill = plotCombo[i, 3:4], bty = "n", cex = 0.8)
  
  # drawing the triangle that represent expression value
  segments(min(bp), -1.25 * ymax * 0.35, max(bp), -1.25 * ymax * 0.35, xpd = NA)
  segments(min(bp), -1.25 * ymax * 0.35, max(bp), -1.25 * ymax * 0.25, xpd = NA)
  segments(max(bp), -1.25 * ymax * 0.35, max(bp), -1.25 * ymax * 0.25, xpd = NA)
  text(max(bp) * 1.06, -1.25 * ymax * 0.27, labels = "EXP", xpd = NA, cex = 0.7)
  text(max(bp) * 1.06, -1.25 * ymax * 0.32, labels = "level", xpd = NA, cex = 0.7)
  
  # Carrying out logistic regression
  dt <- data.table(list$logit.df)
  revc1 <- function(string) {
    start <- revc(substr(string, 1, 1))
    end <- revc(substr(string, 3, 3))
    return(paste0(start, ">", end))
  }
  dt <- dt[mutation %in% c(type, revc1(type)), ]
  dt <- dt[mutation == type, `:=`(class, 1)]
  dt <- dt[mutation != type, `:=`(class, 0)]
  logit.model <- glm(class ~ TPM, family = binomial, data = dt)
  p.value <- summary(logit.model)$coefficients[2, 4]
  text(0.585 * max(bp), 1.15 * ymax, labels = paste0("p = ", signif(p.value, 3)), 
       cex = 0.8, pos = 4, offset = TRUE)
}


#' Plot transcription strand bias with respect to gene expression level.
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param expression.level A \code{data.frame} which contains the transcription
#'   level of genes. \cr See \code{\link{gene.expression.level.example.GRCh37}}
#'   for more details.
#'   
#' @param Ensembl.gene.ID.col Name of column which has the Ensembl gene ID
#'   information in \code{experession.level}.
#' 
#' @param TPM.col Name of column which has the TPM (Transcripts Per Kilobase Million)
#'   information in \code{experession.level}.
#'   
#' @param num.of.bins The number of bins that will be plotted in the graph.
#' 
#' @param plot.type A character string indicating one mutation type to be
#'   plotted. It should be one of "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'   
#' @param ymax Limit for the y axis. If not specified, it defaults to NULL and
#'   the y axis limit equals 1.5 times of the maximum mutation counts in a
#'   specific mutation type.
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
#'   PlotTransBiasExp(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                    expression.level = gene.expression.level.example.GRCh37, 
#'                    Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
#'                    num.of.bins = 4, plot.type = "C>A")
#' }
PlotTransBiasExp <-
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col, TPM.col,
           num.of.bins, plot.type, ymax = NULL) { 
    
    list1 <- 
      StrandBiasAsExpressionLevel(annotated.SBS.vcf, expression.level, 
                                  Ensembl.gene.ID.col, TPM.col, num.of.bins)
    
    PlotTransBiasExp1(list = list1, type = plot.type, 
                      num.of.bins = num.of.bins, ymax = ymax)
    invisible(TRUE)
  }

#' Plot Transcription Strand Bias on Expression level to PDF.
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#' 
#' @param file The name of output file.
#'   
#' @param expression.level A \code{data.frame} which contains the transcription
#'   level of genes. \cr See \code{\link{gene.expression.level.example.GRCh37}}
#'   for more details.
#'   
#' @param Ensembl.gene.ID.col Name of column which has the Ensembl gene ID
#'   information in \code{experession.level}.
#' 
#' @param TPM.col Name of column which has the TPM (Transcripts Per Kilobase Million)
#'   information in \code{experession.level}.
#'   
#' @param num.of.bins The number of bins that will be plotted in the graph.
#' 
#' @param plot.type A vector of character indicating types to be plotted. It
#'   can be one or more types from "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'   The default is to print all the six mutation types.
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
#'   PlotTransBiasExpToPdf(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                         expression.level = gene.expression.level.example.GRCh37, 
#'                         Ensembl.gene.ID.col = "Ensembl.gene.ID", TPM.col = "TPM",
#'                         num.of.bins = 4, 
#'                         plot.type = c("C>A","C>G","C>T","T>A","T>C"), 
#'                         file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasExpToPdf <- 
  function(annotated.SBS.vcf, file, expression.level, Ensembl.gene.ID.col,
           TPM.col, num.of.bins, 
           plot.type = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) {
    list <- 
      StrandBiasAsExpressionLevel(annotated.SBS.vcf, expression.level, 
                                  Ensembl.gene.ID.col, TPM.col, 
                                  num.of.bins)
    
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
    
    opar <- par(mfrow = c(4, 3), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
    on.exit(par(opar))
    num <- length(plot.type)
    type <- plot.type
    for (i in 1:num) {
      PlotTransBiasExp1(list = list, type = plot.type[i], 
                        num.of.bins = num.of.bins, ymax = max(list$plotmatrix))
    }
    
    grDevices::dev.off()
    invisible(TRUE)
    
  }

##################################################################################
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

#' Plot transcription strand bias with respect to distance to transcription
#' start site.
#' 
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param plot.type A character string indicating one mutation type to be
#'   plotted. It should be one of "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#' 
#' @importFrom stats coef
#' 
#' @return \code{invisible(TRUE)}
#' 
#' @references Conaway, J. W. & Conaway, R. C. Transcription Elongation and
#'   Human Disease. \emph{Annu. Rev. Biochem}. 68, 301–319 (1999),
#'   https://doi.org/10.1146/annurev.biochem.68.1.301
#'   
#' @references Hu, J., Adar, S., Selby, C. P., Lieb, J. D. & Sancar, A.
#'   Genome-wide analysis of human global and transcription-coupled excision
#'   repair of UV damage at single-nucleotide resolution. \emph{Genes Dev}. 29,
#'   948–960 (2015), https://doi.org/10.1101/gad.261271.115
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
#'   PlotTransBiasDist2TSS(annotated.SBS.vcf, plot.type = "C>T")
#' }
PlotTransBiasDist2TSS <- function (annotated.SBS.vcf, plot.type){
  output <- dist2TSS(annotated.SBS.vcf, plot.type)
  Plotdist2TSS(output, plot.type)
  return(invisible(TRUE))
}

#' Plot transcription strand bias with respect to distance to transcription
#' start site to a PDF file.
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param plot.type A vector of character indicating types to be plotted. It
#'   should be within "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'
#' @param file The name of output file.
#'    
#' @importFrom stats glm
#' 
#' @return \code{invisible(TRUE)}
#' 
#' @references Conaway, J. W. & Conaway, R. C. Transcription Elongation and
#'   Human Disease. \emph{Annu. Rev. Biochem}. 68, 301–319 (1999),
#'   https://doi.org/10.1146/annurev.biochem.68.1.301
#'   
#' @references Hu, J., Adar, S., Selby, C. P., Lieb, J. D. & Sancar, A.
#'   Genome-wide analysis of human global and transcription-coupled excision
#'   repair of UV damage at single-nucleotide resolution. \emph{Genes Dev}. 29,
#'   948–960 (2015), https://doi.org/10.1101/gad.261271.115
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

