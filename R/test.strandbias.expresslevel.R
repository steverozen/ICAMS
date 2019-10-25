StrandBiasAsExpressionLevel <- 
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col, 
           TPM.col, num.of.bins) {
           
    TPM <- data.table(expression.level[, c(Ensembl.gene.ID.col, TPM.col)])
    names(TPM) <- c("trans.Ensembl.gene.ID", "TPM")
    
    # Delete genes expressed on both strands and NA Ensembl.gene.ID---some weird
    # gene names --needs to discuss further!!!
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
#' @param file The name of output file.
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
#' @param plot.type A vector of character indicating types to be plotted. It
#'   should be within "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
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
  function(annotated.SBS.vcf, expression.level, Ensembl.gene.ID.col,
           TPM.col, num.of.bins, plot.type, file) {
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
