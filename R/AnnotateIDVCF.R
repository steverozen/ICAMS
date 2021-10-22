#' @title Add sequence context to an in-memory ID (insertion/deletion) VCF, and
#'   confirm that they match the given reference genome
#'
#' @param ID.vcf An in-memory ID (insertion/deletion) VCF as a
#'   \code{data.frame}. This function expects that there is a "context base" to
#'   the left, for example REF = ACG, ALT = A (deletion of CG) or REF = A, ALT =
#'   ACC (insertion of CC).
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param flag.mismatches Deprecated. If there are ID variants whose \code{REF}
#'   do not match the extracted sequence from \code{ref.genome}, the function
#'   will automatically discard these variants. See element
#'   \code{discarded.variants} in the return value for more details.
#'   
#' @param name.of.VCF Name of the VCF file.
#' 
#' @param suppress.discarded.variants.warnings Logical. Whether to suppress
#'   warning messages showing information about the discarded variants. Default
#'   is TRUE.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @importFrom BSgenome getSeq seqnames
#'
#' @importFrom stats start end
#' 
#' @importFrom utils write.csv
#' 
#' @importFrom dplyr bind_rows
#'
#' @return A list of elements:
#'   * \code{annotated.vcf}: The original VCF data
#'   frame with two new columns added to the input data frame:
#'       + \code{seq.context}: The sequence embedding the variant.
#'       + \code{seq.context.width}: The width of \code{seq.context} to the left.
#'   * \code{discarded.variants}: \strong{Non-NULL only if} there are variants
#'   that were excluded from the analysis. See the added extra column
#'   \code{discarded.reason} for more details.
#' @md
#' 
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMS"))
#' ID.vcf <- ReadAndSplitVCFs(file, variant.caller = "strelka")$ID[[1]]
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   list <- AnnotateIDVCF(ID.vcf, ref.genome = "hg19")
#'   annotated.ID.vcf <- list$annotated.vcf}
AnnotateIDVCF <- 
  function(ID.vcf, ref.genome, flag.mismatches = 0, name.of.VCF = NULL,
           suppress.discarded.variants.warnings = TRUE) {
    if (nrow(ID.vcf) == 0) {
      return(list(annotated.vcf = ID.vcf))
    }
    
    # Create an empty data frame for discarded variants
    discarded.variants <- ID.vcf[0, ]
    
    # Check and remove discarded variants
    if (suppress.discarded.variants.warnings == TRUE){
      retval <- 
        suppressWarnings(CheckAndRemoveDiscardedVariants(vcf = ID.vcf, 
                                                         name.of.VCF = name.of.VCF)) 
    } else {
      retval <- 
        CheckAndRemoveDiscardedVariants(vcf = ID.vcf, name.of.VCF = name.of.VCF)
    }
    df <- retval$df
    discarded.variants <- 
      dplyr::bind_rows(discarded.variants, retval$discarded.variants)
    
    ref.genome <- NormalizeGenomeArg(ref.genome)
    
    # Remove variants which have the same number of bases 
    # for REF and ALT alleles
    idx <- which(nchar(df$REF) == nchar(df$ALT))
    if (length(idx) > 0) {
      df1 <- df[-idx, ]
      df1.to.remove <- df[idx, ]
      df1.to.remove$discarded.reason <- 
        "ID variant with same number of bases for REF and ALT alleles"
      discarded.variants <- dplyr::bind_rows(discarded.variants, df1.to.remove)
    } else {
      df1 <- df
    }
    
    # stopifnot(nchar(df$REF) != nchar(df$ALT)) # This has to be an indel, maybe a complex indel
    
    # Remove variants which have empty REF or ALT alleles
    idx1 <- which(df1$REF == "" | df1$ALT == "")
    if (length(idx1) > 0) {
      df2 <- df1[-idx1, ]
      df2.to.remove <- df1[idx1, ]
      df2.to.remove$discarded.reason <- "Variant has empty REF or ALT alleles"
      discarded.variants <- dplyr::bind_rows(discarded.variants, df2.to.remove)
    } else {
      df2 <- df1
    }
    
    # We expect either eg ref = ACG, alt = A (deletion of CG) or
    # ref = A, alt = ACC (insertion of CC)
    complex.indels.to.remove <- 
      which(substr(df2$REF, 1, 1) != substr(df2$ALT, 1, 1))
    if (length(complex.indels.to.remove) > 0) {
      df3 <- df2[-complex.indels.to.remove, ]
      df3.to.remove <- df2[complex.indels.to.remove, ]
      df3.to.remove$discarded.reason <- "Complex indel"
      discarded.variants <- 
        dplyr::bind_rows(discarded.variants, df3.to.remove)
    } else {
      df3 <- df2
    }
    stopifnot(substr(df3$REF, 1, 1) == substr(df3$ALT, 1, 1))
    
    # First, figure out how much sequence context is needed.
    var.width <- abs(nchar(df3$ALT) - nchar(df3$REF))
    
    is.del <- nchar(df3$ALT) <= nchar(df3$REF)
    var.width.in.genome <- ifelse(is.del, var.width, 0)
    
    # Set the minimum seq.context.width to be 21, this is to facilitate
    # extended sequence context analysis
    df3$seq.context.width <- ifelse(var.width * 6 < 21, 21, var.width * 6)
    
    # 6 because we need to find out if the insertion or deletion is embedded
    # in up to 5 additional repeats of the inserted or deleted sequence.
    # Then add 1 to avoid possible future issues.
    
    # Extract sequence context from the reference genome
    
    # Check if the format of sequence names in df and genome are the same
    chr.names <- CheckAndFixChrNames(vcf.df = df3, ref.genome = ref.genome,
                                     name.of.VCF = name.of.VCF)
    
    # Create a GRanges object with the needed width.
    Ranges <-
      GRanges(chr.names,
              IRanges(start = df3$POS - df3$seq.context.width, # 10,
                      end = df3$POS + var.width.in.genome + 
                        df3$seq.context.width) # 10
      )
    
    df3$seq.context <- BSgenome::getSeq(ref.genome, Ranges, as.character = TRUE)
    
    seq.to.check <-
      substr(df3$seq.context, df3$seq.context.width + 1,
             df3$seq.context.width + var.width.in.genome + 1)
    
    mismatches <- which(seq.to.check != df3$REF)
    
    if (length(mismatches) > 0) {
      df3$seq.to.check <- seq.to.check
      df4 <- df3[-mismatches, ]
      df4.to.remove <- df3[mismatches, ]
      df4.to.remove$discarded.reason <- 
        "ID variant whose REF alleles do not match the extracted sequence from ref.genome"
      discarded.variants <- 
        dplyr::bind_rows(discarded.variants, df4.to.remove)
    } else {
      df4 <- df3
    }
    
    if (nrow(discarded.variants) > 0) {
      if (suppress.discarded.variants.warnings == TRUE) {
        return(list(annotated.vcf = df4, discarded.variants = discarded.variants))
      } else {
        warning("\nSome ID variants were discarded, see element discarded.variants", 
                " in the return value for more details")
        return(list(annotated.vcf = df4, discarded.variants = discarded.variants))
      }
    } else {
      return(list(annotated.vcf = df4))
    }
  }
