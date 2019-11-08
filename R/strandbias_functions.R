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
#'   PlotTransBiasGeneExp(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                        expression.level = gene.expression.level.example.GRCh37, 
#'                        Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                        TPM.col = "TPM", num.of.bins = 4, plot.type = "C>A")
#' }
PlotTransBiasGeneExp <-
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col, TPM.col,
           num.of.bins, plot.type, ymax = NULL) { 
    
    list1 <- 
      StrandBiasAsExpressionLevel(annotated.SBS.vcf, expression.level, 
                                  Ensembl.gene.ID.col, TPM.col, num.of.bins)
    
    PlotGeneExp(list = list1, type = plot.type, 
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
#'   PlotTransBiasGeneExpToPdf(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                             expression.level = gene.expression.level.example.GRCh37, 
#'                             Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                             TPM.col = "TPM", num.of.bins = 4, 
#'                             plot.type = c("C>A","C>G","C>T","T>A","T>C"), 
#'                             file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasGeneExpToPdf <- 
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
      PlotGeneExp(list = list, type = plot.type[i], 
                  num.of.bins = num.of.bins, ymax = max(list$plotmatrix))
    }
    
    grDevices::dev.off()
    invisible(TRUE)
    
  }

StrandBiasAsExpressionLevel <- 
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col, 
           TPM.col, num.of.bins) {
    
    TPM <- data.table(expression.level[, c(Ensembl.gene.ID.col, TPM.col)])
    names(TPM) <- c("trans.Ensembl.gene.ID", "TPM")
    
    # Delete genes expressed on both strands and rows with NA trans.strand
    vcf <- annotated.SBS.vcf[annotated.SBS.vcf$bothstrand == FALSE, ]
    vcf <- vcf[!is.na(trans.strand), ]
    
    # Construct gene expression levels according to num.of.bins
    dt <- data.table(TPM)
    dt1 <- dt[trans.Ensembl.gene.ID %in% trans.ranges.GRCh37$Ensembl.gene.ID, ]
    # dt1 <- dt[TPM != 0, ]
    dt1$exp.level <- NA
    cutoffs <- 
      stats::quantile(dt1$TPM, c(1:(num.of.bins - 1)/num.of.bins), na.rm = TRUE)
    dt1[, "exp.level"] <- 
      cut(dt1$TPM, breaks = c(-Inf, cutoffs, Inf), labels = (1:num.of.bins))
    
    df <- merge(vcf, dt1, by = "trans.Ensembl.gene.ID", all.x = TRUE)
    df <- df[!is.na(TPM), ]
    df$exp.level <- as.integer(df$exp.level)
    
    # One SBS mutation can be represented by more than 1 row in df if the
    # mutation position falls into the range of multiple transcripts on the same
    # strand. We choose the maximum gene expression level and TPM value of
    # multiple transcripts for one particular mutation.
    df1 <- df[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                  seq.21bases = seq.21bases[1], exp.level = max(exp.level), 
                  TPM = max(TPM)), by = .(CHROM, POS)]
    
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
        result[x, type] <- nrow(df1[mutation == type & exp.level == x, ])
      }
    }
    
    return(list(plotmatrix = result, logit.df = dt, pvalue.overall = p.value))
  }

PlotGeneExp <- function(list, type, num.of.bins, ymax = NULL) {
  result <- list$plotmatrix
  foo1 <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (!type %in% foo1) {
    stop("please input: 'C>A','C>G','C>T','T>A','T>C','T>G'")
  }
  plotCombo <- PlotStrandBiasColorMatrix()
  i <- which(plotCombo[, "Target"] == type)
  tmp <- t(result[, c(plotCombo[i, 2], plotCombo[i, 1])])
  colnames(tmp) <- c(1:num.of.bins)
  
  if (is.null(ymax)) {
    ymax <- max(tmp)
  }
  bp <- barplot(tmp, beside = TRUE, main = plotCombo[i, "Target"], 
                ylim = c(0, ifelse(ymax != 0, 1.5 * ymax, 5)), 
                col = plotCombo[i, 4:3], axisnames = FALSE,
                ylab = "counts")
  legend.list <- legend("topright", legend = c("Transcribed", "Untranscribed"), 
                        fill = plotCombo[i, 4:3], bty = "n", cex = 0.8)
  text(-min(bp) + 0.25, ifelse(ymax != 0, -ymax * 0.2, -1), 
       labels = "Expression", adj = 0.5, xpd = NA)
  
  # Draw a triangle that represents expression value
  if (ymax != 0) {
    y.pos <- c(-ymax * 0.2, -ymax * 0.2, -ymax * 0.1)
  } else {
    y.pos <- c(-1, -1, -0.5)
  }
  polygon(c(min(bp), max(bp), max(bp)), y.pos, xpd = NA)
  text(min(bp), ifelse(ymax != 0, -ymax * 0.3, -1.5), 
       labels = "Low", xpd = NA)
  text(max(bp), ifelse(ymax != 0, -ymax * 0.3, -1.5), 
       labels = "High", xpd = NA)
  
  # Carrying out logistic regression
  dt <- data.table(list$logit.df)
  dt <- dt[mutation %in% c(type, revcSBS6(type)), ]
  dt <- dt[mutation == type, `:=`(class, 1)]
  dt <- dt[mutation != type, `:=`(class, 0)]
  if (nrow(dt) != 0 ) {
    logit.model <- glm(class ~ TPM, family = binomial, data = dt)
    p.value <- summary(logit.model)$coefficients[2, 4]
    text(legend.list$rect$left * 1.005, 1.15 * ymax, 
         labels = paste0("p value = ", signif(p.value, 3)), 
         cex = 0.8, pos = 4)
  }
  
}

##################################################################################
#' Plot transcription strand bias with respect to distance to transcription
#' start site.
#' 
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param plot.type A character string indicating one mutation type to be
#'   plotted. It should be one of "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'   
#' @param ymax Limit for the y axis. If not specified, it defaults to NULL and
#'   the y axis limit equals 1.2 times of the maximum mutation counts in a
#'   specific mutation type.
#' 
#' @importFrom stats coef
#' 
#' @return \code{invisible(TRUE)}
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
#'   PlotTransBiasDistToTSS(annotated.SBS.vcf, plot.type = "C>T")
#' }
PlotTransBiasDistToTSS <- 
  function (annotated.SBS.vcf, plot.type, ymax = NULL){
  result <- DistToTSS(annotated.SBS.vcf)
  PlotDistToTSS(result, plot.type, ymax)
  invisible(TRUE)
}

#' Plot transcription strand bias with respect to distance to transcription
#' start site to a PDF file.
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param file The name of output file.
#' 
#' @param plot.type A vector of character indicating types to be plotted. It
#'   can be one or more types from "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'   The default is to print all the six mutation types.
#'    
#' @importFrom stats glm
#' 
#' @return \code{invisible(TRUE)}
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
#'   PlotTransBiasDistToTSSToPdf(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                              file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasDistToTSSToPdf <- 
  function(annotated.SBS.vcf, file, 
           plot.type = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) {
    
  list <- DistToTSS(annotated.SBS.vcf)  
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  opar <- par(mfrow = c(4, 2), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  on.exit(par(opar))
  num <- length(plot.type)
  type <- plot.type
  
  for (i in 1:num) {
    PlotDistToTSS(list = list, 
                 plot.type = plot.type[i], 
                 ymax = max(list$plotmatrix[, c(1:12)]))
  }
  
  grDevices::dev.off()
  invisible(TRUE)
  
}

DistToTSS <- function(annotated.SBS.vcf) {
  # Delete rows where a mutation falls on both strands or has NA trans.strand value
  vcf <- annotated.SBS.vcf[annotated.SBS.vcf$bothstrand == FALSE, ]
  vcf <- vcf[!is.na(trans.strand), ]
  
  distToTSSbin <- function(df, pos = "POS", TSS = "trans.start.pos") {
    df <- data.frame(df)
    df$distToTSS <- abs(df[, pos] - df[, TSS])
    df$distBins <- NULL
    df$distBins[df$distToTSS <= 1e+05] <- 1
    df$distBins[df$distToTSS <= 2e+05 & df$distToTSS > 1e+05] <- 2
    df$distBins[df$distToTSS <= 3e+05 & df$distToTSS > 2e+05] <- 3
    df$distBins[df$distToTSS <= 4e+05 & df$distToTSS > 3e+05] <- 4
    df$distBins[df$distToTSS <= 5e+05 & df$distToTSS > 4e+05] <- 5
    df$distBins[df$distToTSS <= 6e+05 & df$distToTSS > 5e+05] <- 6
    df$distBins[df$distToTSS <= 7e+05 & df$distToTSS > 6e+05] <- 7
    df$distBins[df$distToTSS <= 8e+05 & df$distToTSS > 7e+05] <- 8
    df$distBins[df$distToTSS <= 9e+05 & df$distToTSS > 8e+05] <- 9
    df$distBins[df$distToTSS <= 1e+06 & df$distToTSS > 9e+05] <- 10
    df$distBins[df$distToTSS > 1e+06 ] <- 11
    return(df)
  }
  
  vcf <- data.table(distToTSSbin(vcf))
  
  # One SBS mutation can be represented by more than 1 row in vcf if the
  # mutation position falls into the range of multiple transcripts on the same
  # strand. We choose the minimum distBins and distToTSS of
  # multiple transcripts for one particular mutation.
  df <- vcf[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                seq.21bases = seq.21bases[1], distToTSS = min(distToTSS), 
                distBins = min(distBins)), by = .(CHROM, POS)]
  
  df$mutation <- paste0(df$REF, ">", df$ALT)
  i <- which(df$trans.strand == "-")
  df$mutation[i] <- paste0(revc(df$REF[i]), ">", revc(df$ALT[i]))
  
  # Carry out logistic regression. 
  type <- c("C>A", "C>T", "C>G", "T>A", "T>C", "T>G")
  df$class <- 0
  df$class[df$mutation %in% type] <- 1
  
  fit <- glm(as.factor(df$class) ~ df$distToTSS, family = 'binomial')
  pvalue <- coef(summary(fit))[2, 4]
  
  output <- matrix(data = NA, nrow = 11, ncol = 12)
  rownames(output) <- paste0("group", 1:11)
  colnames(output) <-
    c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")
  
  i <- as.data.frame(table(df$mutation[df$distBins == 1]))
  # Ensure all mutation classes are represented
  rownames(i) <- i$Var1
  i <- i[colnames(output), ]
  i$Var1 <- colnames(output)
  i$Freq[is.na(i$Freq)] <- 0
  
  output[1, ] <- i$Freq
  for (class in colnames(output)) {
    for (row in 2:11) {
      output[row, class] <- sum(df$mutation[df$distBins == row] == class)
    }
  }
  
  output <- as.matrix(output)
  output <- cbind(output, Untranscribed = rowSums(output[, c(4:6, 10:12)]))
  output <- cbind(output, Transcribed = rowSums(output[, c(1:3, 7:9)]))
  
  return(list(plotmatrix = output, logit.df = df, 
              pvalue.overall = pvalue))
}

PlotDistToTSS <- function(list, plot.type, ymax = NULL) {
  output1 <- list$plotmatrix
  if (plot.type == "combined") {
    output2 <- output1[, c("Transcribed", "Untranscribed")]
  }
  output2 <- output1[, c(revcSBS6(plot.type), plot.type)]
  rownames(output2) <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", 
                         "8-9", "9-10", ">10")
  output2 <- t(output2)
  
  color.matrix <- PlotStrandBiasColorMatrix()
  colors <- color.matrix[color.matrix[, "Target"] %in% plot.type, c(4, 3)]
  
  if (is.null(ymax)) {
    ymax <- max(output2)
  }
  
  bp <- barplot(output2, beside = TRUE, 
                ylim = c(0, ifelse (ymax != 0, ymax * 1.2, 5)), 
                col = colors, cex.axis = 0.8, cex.names = 0.9,
                xlab = expression('Distance to TSS (10'^5*' bp)'),
                ylab = "counts", main = plot.type) 
  legend.list <- legend("topright", c("Transcribed", "Untranscribed"), bty = "n", 
                 fill = colors,  cex = 0.8)
        
  
  # Carry out logistic regression.
  df <- list$logit.df
  df1 <- df[df$mutation %in% c(plot.type, revcSBS6(plot.type)), ]
  if (nrow(df1) != 0) {
    fit <- glm(as.factor(df1$class) ~ df1$distToTSS, family = "binomial")
    pvalue <- coef(summary(fit))[2, 4]
    text(legend.list$rect$left * 1.005, ymax * 0.93, 
         paste0("p value = ", signif(pvalue, 3)), pos = 4, cex = 0.8, xpd = NA)
  }
  
}

revcSBS6 <- function(string) {
  start <- revc(substr(string, 1, 1))
  end <- revc(substr(string, 3, 3))
  return(paste0(start, ">", end))
}

PlotStrandBiasColorMatrix <- function() {
  matrix <- cbind(Target = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), 
                  rev = c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C"), 
                  colLight = c("#87CEFA",  "#A9A9A9", "#FFB6C1", 
                               "#DCDCDC", "#CCFFCC", "#FFB2BF"), 
                  colDark = c("#0000ff", "#000000", "#ff4040", 
                              "#838383", "#40ff40", "#ff667f")) 
  return(matrix)
}
