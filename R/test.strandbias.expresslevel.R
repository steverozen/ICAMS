StrandBiasAsExpressionLevel <- 
  function(annotated.SBS.vcf, expression.level, ref.genome, n # num.bins 
           ) {
    
    ref.genome <- NormalizeGenomeArg(ref.genome)
    if (ref.genome@pkgname == "BSgenome.Hsapiens.1000genomes.hs37d5") {
      trans.ranges <- trans.ranges.GRCh37
    } else if (ref.genome@pkgname == "BSgenome.Hsapiens.UCSC.hg38)") {
      trans.ranges <- trans.ranges.GRCh38
    } else if (ref.genome@pkgname == "BSgenome.Mmusculus.UCSC.mm10") {
      trans.ranges <- trans.ranges.GRCm38
    }
    
    TPM <- data.table(expression.level)
    names(TPM) <- c("trans.gene.symbol", "TPM")
    
    # Exclude the lncRNA 
    TPM <- TPM[trans.gene.symbol %in% trans.ranges$gene.symbol, ]
    
    # Remove duplicate entries---why multiple genes show twice in our RNA-seq
    # data---needs to discuss further!!!
    TPM <- unique(TPM, by = "trans.gene.symbol")
    
    # Delete genes expressed on both strands and NA gene symbol---some wierd
    # gene names --needs to discuss further!!!
    vcf <- annotated.SBS.vcf[-which(annotated.SBS.vcf$bothstrand == TRUE), ]
    vcf <- vcf[-which(is.na(vcf$trans.gene.symbol) == TRUE), ]
    df <- merge(vcf, TPM, by = "trans.gene.symbol", all.x = TRUE)
    df <- df[-which(is.na(df$TPM) == TRUE), ]
    df$Exp_Level <- NA
    cutoffs <- stats::quantile(df$TPM, c(1:(n - 1)/n), na.rm = T)
    df[, "Exp_Level"] <- 
      cut(df$TPM, breaks = c(-Inf, cutoffs, Inf), labels = (1:n))
    
    # Orient the ref.context and the var.context column
    df$mutation <- paste0(substr(df$seq.21bases, 10, 12), df$ALT)
    df[trans.strand == "-", `:=`(mutation, RevcSBS96(mutation))]
    df$mutation <- 
      paste0((substr(df$mutation, 2, 2)), ">", substr(df$mutation, 4, 4))
    
    # Create a table table for logistic regression
    dt <- df[, c("mutation", "TPM")]
    type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    dt <- dt[mutation %in% type, `:=`(class, 1)]
    dt <- dt[!mutation %in% type, `:=`(class, 0)]
    logit.model <- stats::glm(class ~ TPM, family = binomial, data = dt)
    p.value <- summary(logit.model)$coefficients[2, 4]
    
    # Plot transcriptional strand bias as a function of gene expression
    result <- matrix(data = 0, nrow = n, ncol = 12)
    rownames(result) <- c(1:n)
    
    
    mutation.type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                       "G>T", "G>C", "G>A", "A>T", "A>G", "A>C")
    #i <- table(df$mutation[df$Exp_Level == 1])
    colnames(result) <- mutation.type
    for (x in 1:n) {
      for (j in 1:12) {
        type <- mutation.type[j]
        result[x, type] <- nrow(df[mutation == type & Exp_Level == x, ])
      }
    }
    
    return(list(plotmatrix = result, logit.df = dt, p.value = p.value))
  }


PlotTransBiasExp1 <- function(list, type, n, ymax = NULL) {
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
  i = which(plotCombo[, "Target"] == type)
  tmp <- t(result[, c(plotCombo[i, 1], plotCombo[i, 2])])
  colnames(tmp) <- c(1:n)
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
  text(0.75 * max(bp), 1.15 * ymax, labels = paste0("p = ", signif(p.value, 3)), 
       cex = 0.8, pos = 4, offset = TRUE)
}


#' Plot transcriptional strand bias with respect to gene expression level
#'
#' @param annotated.SBS.vcf An annotated SBS VCF as a data.table which contains
#'   sequence context and transcript information. Please refer to
#'   \code{\link{AnnotateSBSVCF}} for more details.
#'   
#' @param expression.level A two columns data.frame which contains the
#'   transcription level of genes. The first column should be the gene symbols.
#'   The second column should be the numeric expression level e.g.TPM
#'   (Transcripts Per Kilobase Million).
#'   
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'   
#' @param n The number of bins that will be plotted in the graph.
#' 
#' @param plot.type A character string indicating one mutation type to be
#'   plotted. It should be one of "C>A", "C>G", "C>T", "T>A", "T>C", "T>G".
#'   
#' @param ymax Limit for the y axis. If not specified, it defaults to NULL and
#'   the y axis limit equals to 1.5 times of the maximum mutation counts in a
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
#'                    expression.level = sample.gene.expression.value.GRCh37, 
#'                    ref.genome = "hg19", n = 4, plot.type = "C>A")
#' }
PlotTransBiasExp <-
  function(annotated.SBS.vcf, expression.level, ref.genome, 
                             n, plot.type, ymax = NULL) { # change n to num.bins
  list1 <- StrandBiasAsExpressionLevel(annotated.SBS.vcf, expression.level, 
                                       ref.genome, n)
  PlotTransBiasExp1(list = list1, type = plot.type, n = n, ymax = ymax)
  invisible(TRUE)
}



#' Plot Transcriptional Strand Bias on Expression level to PDF
#'
#' @param file The name of output file
#' 
#' @param annotated.SBS.vcf An annotated SBS VCF as a data.table which contains
#'   sequence context and transcript information. Please refer to
#'   \code{\link{AnnotateSBSVCF}} for more details.
#'   
#' @param expression.level A two columns data.frame which contains the
#'   transcription level of genes. The first column should be the gene symbols.
#'   The second column should be the numeric expression level e.g.TPM
#'   (Transcripts Per Kilobase Million).
#'   
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'   
#' @param n The number of bins that will be plotted in the graph.
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
#'                         expression.level = sample.gene.expression.value.GRCh37, 
#'                         ref.genome = "hg19", n = 4, 
#'                         plot.type = c("C>A","C>G","C>T","T>A","T>C"), 
#'                         file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasExpToPdf <- function(annotated.SBS.vcf, expression.level, 
                                  ref.genome, n, plot.type, file) {
  list <- StrandBiasAsExpressionLevel(annotated.SBS.vcf, expression.level, 
                                      ref.genome, n)
  
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  opar <- par(mfrow = c(4, 3), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  on.exit(par(opar))
  num <- length(plot.type)
  type <- plot.type
  for (i in 1:num) {
    PlotTransBiasExp1(list = list, type = plot.type[i], 
                      n = n, ymax = max(list$plotmatrix))
  }
  
  grDevices::dev.off()
  invisible(TRUE)
  
}
