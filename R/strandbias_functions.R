#' Plot transcription strand bias with respect to gene expression values
#'
#' @param annotated.SBS.vcf An SBS VCF annotated by
#'   \code{\link{AnnotateSBSVCF}}. It \strong{must} have transcript range
#'   information added.
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
#' @param damaged.base One of \code{NULL}, \code{"purine"} or
#'   \code{"pyrimidine"}. This function allocates approximately
#'   equal numbers of mutations from \code{damaged.base} into
#'   each of \code{num.of.bins} bin by expression level. E.g.
#'   if \code{damaged.base} is \code{"purine"}, then mutations from
#'   A and G will be allocated in approximately equal numbers to
#'   each expression-level bin. The rationale for the name \code{damaged.base}
#'   is that the direction of strand bias is a result of whether the damage
#'   occurs on a purine or pyrimidine.
#'   If \code{NULL}, the function attempts to infer the \code{damaged.base}
#'   based on mutation counts.
#'   
#' @param ymax Limit for the y axis. If not specified, it defaults to NULL and
#'   the y axis limit equals 1.5 times of the maximum mutation counts in a
#'   specific mutation type.
#'   
#' @importFrom stats glm
#'   
#' @return A list whose first element is a logic value indicating whether the
#'   plot is successful. The second element is a named numeric vector containing
#'   the p-values printed on the plot.
#'   
#' @section Note: 
#' The p-values are calculated by logistic regression using function
#' \code{\link[stats]{glm}}. The dependent variable is labeled "1" and "0" if
#' the mutation from \code{annotated.SBS.vcf} falls onto the untranscribed and
#' transcribed strand respectively. The independent variable is the binary
#' logarithm of the gene expression value from \code{expression.data} plus one,
#' i.e. \eqn{log_2 (x + 1)}{log2 (x + 1)} where \eqn{x} stands for gene
#' expression value.
#'   
#' @export
#'
#' @examples 
#' file <- c(system.file("extdata/Strelka-SBS-vcf/",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)
#'   PlotTransBiasGeneExp(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                        expression.data = gene.expression.data.HepG2, 
#'                        Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                        expression.value.col = "TPM", 
#'                        num.of.bins = 4, plot.type = "C>A")
#' }
PlotTransBiasGeneExp <-
  function(annotated.SBS.vcf, expression.data, Ensembl.gene.ID.col, 
           expression.value.col, num.of.bins, plot.type, damaged.base = NULL,
           ymax = NULL) { 
    list <- StrandBiasGeneExp(annotated.SBS.vcf, expression.data, 
                              Ensembl.gene.ID.col, expression.value.col, 
                              num.of.bins, damaged.base)
    
    PlotGeneExp(list = list, type = plot.type, 
                num.of.bins = num.of.bins, ymax = ymax)
    return(list(plot.success = TRUE, p.values = list$p.values))
  }

#' Plot transcription strand bias with respect to gene expression values to a
#' PDF file
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
#' @inherit PlotTransBiasGeneExp return
#' 
#' @inheritSection PlotTransBiasGeneExp Note
#' 
#' @export
#'
#' @examples
#' file <- c(system.file("extdata/Strelka-SBS-vcf/",
#'                       "Strelka.SBS.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs(file)
#' SBS.vcf <- list.of.vcfs$SBS.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.SBS.vcf <- AnnotateSBSVCF(SBS.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37)
#'   PlotTransBiasGeneExpToPdf(annotated.SBS.vcf = annotated.SBS.vcf, 
#'                             expression.data = gene.expression.data.HepG2, 
#'                             Ensembl.gene.ID.col = "Ensembl.gene.ID", 
#'                             expression.value.col = "TPM", 
#'                             num.of.bins = 4, 
#'                             plot.type = c("C>A","C>G","C>T","T>A","T>C"), 
#'                             file = file.path(tempdir(), "test.pdf"))
#' }
PlotTransBiasGeneExpToPdf <- 
  function(annotated.SBS.vcf, file, expression.data, Ensembl.gene.ID.col,
           expression.value.col, num.of.bins, 
           plot.type = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
           damaged.base = NULL) {
    list <- StrandBiasGeneExp(annotated.SBS.vcf, expression.data, 
                              Ensembl.gene.ID.col, expression.value.col, 
                              num.of.bins, damaged.base)
    
    # Setting the width and length for A4 size plotting
    grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
    
    opar <- par(mfrow = c(4, 3), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
    on.exit(par(opar))
    num <- length(plot.type)
    
    # Calculate ymax on the plot based on plot.type
    ymax <- max(list$plotmatrix[, c(plot.type, revcSBS6(plot.type))])
    
    for (i in 1:num) {
      PlotGeneExp(list = list, type = plot.type[i], 
                  num.of.bins = num.of.bins, ymax = ymax)
    }
    
    grDevices::dev.off()
    return(list(plot.success = TRUE, p.values = list$p.values))
    
  }

#' @keywords internal
CalculateExpressionLevel <- function(dt, num.of.bins, type, damaged.base) {
  dt.purine <- dt[mutation == revcSBS6(type), ]
  dt.pyrimidine <- dt[mutation == type, ]
  
  if (is.null(damaged.base)) {
    if (nrow(dt.purine) >= nrow(dt.pyrimidine)) {
      dt1 <- dt.purine
      dt2 <- dt.pyrimidine
    } else {
      dt1 <- dt.pyrimidine
      dt2 <- dt.purine
    }
  } else if (damaged.base == "purine") {
    dt1 <- dt.purine
    dt2 <- dt.pyrimidine
  } else if (damaged.base == "pyrimidine") {
    dt1 <- dt.pyrimidine
    dt2 <- dt.purine
  } else {
    stop('\nThe input for damaged.base must be one of NULL, "purine"',  
         '\n"pyrimidine". Please check the value. ')
  }
  
  setorder(dt1, exp.value)
  setorder(dt2, exp.value)
  if (num.of.bins == 1) {
    dt[, exp.level := num.of.bins]
    return(dt)
  } else {
    if (nrow(dt1) <= num.of.bins) {
      dt1$exp.level <- 1:nrow(dt1)
      max.exp.value <- max(dt1$exp.value)
      dt3 <- dt2[exp.value > max.exp.value, ]
      break.points1 <- c(0, dt1[1:nrow(dt1), exp.value])
      dt3[, exp.level := nrow(dt1) + 
            cut(1:nrow(dt3), breaks = num.of.bins - nrow(dt1), labels = FALSE)]
      index <- cumsum(table(dt3$exp.level))
      break.points2 <- c(dt3[index, exp.value])
      break.points <- c(break.points1, break.points2)
      
      GetExpLevel <- function(i, exp.value, break.points) {
        lower <- break.points[i]
        upper <- break.points[i + 1]
        total.match <- sum(exp.value >lower & exp.value <= upper)
        return(rep(i, total.match))
      }
      setorder(dt, exp.value)
      list <- lapply(1:num.of.bins, FUN = GetExpLevel, exp.value = dt$exp.value, 
                     break.points = break.points)
      exp.level <- do.call("c", list)
      dt$exp.level <- exp.level
      return(dt)
    } else {
      dt1[, exp.level := cut(1:nrow(dt1), breaks = num.of.bins, labels = FALSE)]
      idx <- cumsum(table(dt1$exp.level)) + 1
      break.points <- c(0, dt1[idx[1:(num.of.bins - 1)], exp.value], 
                        max(dt$exp.value) + 1)
      dup.idx <- which(duplicated(break.points))
      
      GetExpLevel1 <- function(i, exp.value, break.points) {
        lower <- break.points[i]
        upper <- break.points[i + 1]
        total.match <- sum(exp.value >= lower & exp.value < upper)
        return(rep(i, total.match))
      }
      
      setorder(dt, exp.value)
      list <- lapply(1:num.of.bins, FUN = GetExpLevel1, exp.value = dt$exp.value, 
                     break.points = break.points)
      exp.level <- do.call("c", list)
      
      if (length(dup.idx) != 0) {
        idx.max <- max(dup.idx)
        idx2 <- c(0, cumsum(table(dt1$exp.level) * 2))
        for (i in 1:(idx.max - 1)) {
          exp.level[(idx2[i] + 1):idx2[i + 1]] <- i
        }
      }
      
      dt$exp.level <- exp.level
      return(dt)
    }
  }
}

#' @keywords internal
CalculatePValues <- function(dt) {
  # Construct a matrix which can later store the p-values
  p.values <- rep(NA, 7)
  names(p.values) <- c("overall", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  setDT(dt)
  dt <- dt[mutation %in% type, class := 1]
  dt <- dt[!mutation %in% type, class := 0]
  
  logit.model <- stats::glm(class ~ log2(exp.value + 1), 
                            family = binomial, 
                            data = dt)
  p.values["overall"] <- summary(logit.model)$coefficients[2, 4]
  
  for (i in 1:6) {
    dt1 <- dt[mutation %in% c(type[i], revcSBS6(type[i])), ]
    if (nrow(dt1) != 0) {
      logit.model1 <- stats::glm(class ~ log2(exp.value + 1), 
                                 family = binomial, 
                                 data = dt1)
      # Calculate the number of unique expression values in dt1
      num.exp.value <- length(unique(dt1$exp.value))
      if (num.exp.value == 1) {
        # If there is only one unique expression value in dt1, then this
        # predictor variable will be dropped from the logistic regression model
        p.values[type[i]] <- NA
      } else {
        p.values[type[i]] <- summary(logit.model1)$coefficients[2, 4]
      }
    }
  }
  
  return(p.values)
}

#' @keywords internal
revcSBS6 <- function(string) {
  start <- revc(substr(string, 1, 1))
  end <- revc(substr(string, 3, 3))
  return(paste0(start, ">", end))
}

#' @keywords internal
PlotStrandBiasColorMatrix <- function() {
  matrix <- cbind(Target = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), 
                  rev = c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C"), 
                  colLight = c("#87CEFA",  "#A9A9A9", "#FFB6C1", 
                               "#DCDCDC", "#CCFFCC", "#FFB2BF"), 
                  colDark = c("#0000ff", "#000000", "#ff4040", 
                              "#838383", "#40ff40", "#ff667f")) 
  return(matrix)
}

#' @keywords internal
StrandBiasGeneExp <- 
  function(annotated.SBS.vcf, expression.data, Ensembl.gene.ID.col, 
           expression.value.col, num.of.bins, damaged.base) {
    dt1 <- expression.data[, c(Ensembl.gene.ID.col, expression.value.col),
                           with = FALSE]
    idx <- which(colnames(dt1) == Ensembl.gene.ID.col)
    colnames(dt1)[idx] <- "trans.Ensembl.gene.ID"
    idx1 <- which(colnames(dt1) == expression.value.col)
    colnames(dt1)[idx1] <- "exp.value"
    
    # Delete mutations which fall to transcripts on both strands and rows with
    # NA trans.strand
    vcf <- annotated.SBS.vcf[annotated.SBS.vcf$bothstrand == FALSE, ]
    vcf <- vcf[!is.na(trans.strand), ]
    
    df <- merge(vcf, dt1, by = "trans.Ensembl.gene.ID", all.x = TRUE)
    df <- df[!is.na(exp.value), ]
    
    # One SBS mutation can be represented by more than 1 row in df if the
    # mutation position falls into the range of multiple transcripts on the same
    # strand. We choose the maximum gene expression level and exp.value of
    # multiple transcripts for one particular mutation.
    df1 <- df[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                  seq.21bases = seq.21bases[1], 
                  trans.Ensembl.gene.ID = trans.Ensembl.gene.ID[1],
                  trans.gene.symbol = trans.gene.symbol[1],
                  trans.start.pos = trans.start.pos[1], 
                  trans.end.pos = trans.end.pos[1],
                  exp.value = max(exp.value)), by = .(CHROM, POS)]
    
    # Orient the ref.context and the var.context column
    df1$mutation <- paste0(substr(df1$seq.21bases, 10, 12), df1$ALT)
    df1[trans.strand == "-", `:=`(mutation, RevcSBS96(mutation))]
    df1$mutation <- 
      paste0((substr(df1$mutation, 2, 2)), ">", substr(df1$mutation, 4, 4))
    
    # Carry out logistic regression and get the p-values
    p.values <- CalculatePValues(df1)
    
    # Construct a matrix that later can be used to plot transcriptional strand
    # bias as a function of gene expression
    result <- matrix(data = 0, nrow = num.of.bins, ncol = 12)
    rownames(result) <- c(1:num.of.bins)
    
    mutation.type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                       "G>T", "G>C", "G>A", "A>T", "A>G", "A>C")
    colnames(result) <- mutation.type
    
    for (i in 1:6) {
      type <- mutation.type[i]
      df2 <- df1[mutation %in% c(type, revcSBS6(type)), ]
      if (nrow(df2) == 0) {
        result[, type] <- 0
        result[, revcSBS6(type)] <- 0
      } else {
        df2 <- CalculateExpressionLevel(df2, num.of.bins, type, damaged.base)
        for (j in 1:num.of.bins) {
          result[j, type] <- nrow(df2[mutation == type & exp.level == j, ])
          result[j, revcSBS6(type)] <- 
            nrow(df2[mutation == revcSBS6(type) & exp.level == j, ])
        }
      }
    }
    
    return(list(plotmatrix = result, vcf.df = df1, 
                p.values = p.values))
  }

#' @keywords internal
PlotGeneExp <- function(list, type, num.of.bins, ymax = NULL) {
  result <- list$plotmatrix
  allowed.type <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (!type %in% allowed.type) {
    stop("please input: 'C>A','C>G','C>T','T>A','T>C','T>G'")
  }
  plotCombo <- PlotStrandBiasColorMatrix()
  i <- which(plotCombo[, "Target"] == type)
  if (num.of.bins == 1) {
    tmp <- as.matrix(result[, c(plotCombo[i, 2], plotCombo[i, 1])])
  } else {
    tmp <- t(result[, c(plotCombo[i, 2], plotCombo[i, 1])])
  }
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
  text(mean(c(min(bp), max(bp))), ifelse(ymax != 0, -ymax * 0.3, -1.5), 
       labels = "Expression", xpd = NA)
    
  if (!is.na(list$p.values[type])) {
    text(legend.list$rect$left * 1.005, 1.15 * ymax, 
         labels = paste0("p-value = ", signif(list$p.values[type], 2)), 
         cex = 0.8, pos = 4)
  }
}