#' Plot transcription strand bias with respect to gene expression values.
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by \code{\link{AnnotateSBSVCF}}.
#'
#' @param expression.data A \code{\link{data.table}} which contains the
#'   expression values of genes. \cr See \code{\link{GeneExpressionData}} for more
#'   details.
#'   
#' @param Ensembl.gene.ID.col Name of column which has the Ensembl gene ID
#'   information in \code{expression.data}.
#' 
#' @param expression.value.col Name of column which has the gene expression
#'   values in \code{expression.data}.
#'   
#' @param num.of.bins The number of bins that will be plotted on the graph.
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
#'                        expression.data = gene.expression.data.HepG2, 
#'                        Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                        expression.value.col = "TPM", num.of.bins = 4, 
#'                        plot.type = "C>A")
#' }
PlotTransBiasGeneExp <-
  function(annotated.SBS.vcf, expression.data, Ensembl.gene.ID.col, 
           expression.value.col, num.of.bins, plot.type, ymax = NULL) { 
    list1 <- 
      StrandBiasGeneExp(annotated.SBS.vcf, expression.data, 
                        Ensembl.gene.ID.col, expression.value.col, 
                        num.of.bins)
    
    PlotGeneExp(list = list1, type = plot.type, 
                num.of.bins = num.of.bins, ymax = ymax)
    invisible(TRUE)
  }

#' Plot transcription strand bias with respect to gene expression values to a
#' PDF file.
#' 
#' @inheritParams PlotTransBiasGeneExp
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
#'                             expression.data = gene.expression.data.HepG2, 
#'                             Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                             expression.value.col = "TPM", num.of.bins = 4, 
#'                             plot.type = c("C>A","C>G","C>T","T>A","T>C"), 
#'                             file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasGeneExpToPdf <- 
  function(annotated.SBS.vcf, file, expression.data, Ensembl.gene.ID.col,
           expression.value.col, num.of.bins, 
           plot.type = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) {
    list <- StrandBiasGeneExp(annotated.SBS.vcf, expression.data, 
                              Ensembl.gene.ID.col, expression.value.col, 
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

StrandBiasGeneExp <- 
  function(annotated.SBS.vcf, expression.data, Ensembl.gene.ID.col, 
           expression.value.col, num.of.bins) {
    setDT(expression.data)
    data <- expression.data[, c(Ensembl.gene.ID.col, expression.value.col), 
                           with = FALSE]
    names(data) <- c("trans.Ensembl.gene.ID", "exp.value")
    
    # Delete genes expressed on both strands and rows with NA trans.strand
    vcf <- annotated.SBS.vcf[annotated.SBS.vcf$bothstrand == FALSE, ]
    vcf <- vcf[!is.na(trans.strand), ]
    
    # Construct gene expression levels according to num.of.bins
    dt <- data.table(data)
    dt1 <- dt[trans.Ensembl.gene.ID %in% trans.ranges.GRCh37$Ensembl.gene.ID, ]
    
    dt1$exp.level <- NA
    cutoffs <- 
      stats::quantile(dt1$exp.value, c(1:(num.of.bins - 1)/num.of.bins), na.rm = TRUE)
    dt1[, "exp.level"] <- 
      cut(dt1$exp.value, breaks = c(-Inf, cutoffs, Inf), labels = (1:num.of.bins))
    
    df <- merge(vcf, dt1, by = "trans.Ensembl.gene.ID", all.x = TRUE)
    df <- df[!is.na(exp.value), ]
    df$exp.level <- as.integer(df$exp.level)
    
    # One SBS mutation can be represented by more than 1 row in df if the
    # mutation position falls into the range of multiple transcripts on the same
    # strand. We choose the maximum gene expression level and exp.value of
    # multiple transcripts for one particular mutation.
    df1 <- df[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                  seq.21bases = seq.21bases[1], 
                  trans.Ensembl.gene.ID = trans.Ensembl.gene.ID[1],
                  trans.gene.symbol = trans.gene.symbol[1],
                  exp.level = max(exp.level), 
                  exp.value = max(exp.value)), by = .(CHROM, POS)]
    
    # Orient the ref.context and the var.context column
    df1$mutation <- paste0(substr(df1$seq.21bases, 10, 12), df1$ALT)
    df1[trans.strand == "-", `:=`(mutation, RevcSBS96(mutation))]
    df1$mutation <- 
      paste0((substr(df1$mutation, 2, 2)), ">", substr(df1$mutation, 4, 4))
    
    # Create a table for logistic regression
    # dt <- df1[, c("mutation", "exp.value")]
    dt <- df1
    type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    dt <- dt[mutation %in% type, `:=`(class, 1)]
    dt <- dt[!mutation %in% type, `:=`(class, 0)]
    # max.expression <- max(dt$exp.value)
    # dt1 <- dt[exp.value != max.expression, ]
    dt1 <- dt[exp.value != outlier(exp.value), ]
    logit.model <- stats::glm(class ~ exp.value, family = binomial, 
                              data = dt1)
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
  allowed.type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (!type %in% allowed.type) {
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
    #values <- sort(dt$exp.value, decreasing = TRUE)[1:2]
    #dt1 <- dt[!(exp.value %in% values), ]
    dt1 <- dt[exp.value != outlier(exp.value), ]
    logit.model <- glm(class ~ exp.value, family = binomial, data = dt1)
    p.value <- summary(logit.model)$coefficients[2, 4]
    text(legend.list$rect$left * 1.005, 1.15 * ymax, 
         labels = paste0("p value = ", signif(p.value, 2)), 
         cex = 0.8, pos = 4)
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