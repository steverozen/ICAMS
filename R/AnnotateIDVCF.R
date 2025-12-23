library(GenomicRanges)

#' @title Add sequence context and transcript information to an in-memory ID
#'   (insertion/deletion) VCF, and confirm that they match the given reference
#'   genome
#'
#' @param ID.vcf An in-memory ID (insertion/deletion) VCF as a
#'   \code{data.frame}. This function expects that there is a "context base" to
#'   the left, for example REF = ACG, ALT = A (deletion of CG) or REF = A, ALT =
#'   ACC (insertion of CC).
#'
#' @param ref.genome Can be a string or a BSgenome. If a string, it should
#' a well-known name for reference genome in BSgenome
#'
#' @param flag.mismatches Deprecated. If there are ID variants whose \code{REF}
#'   do not match the extracted sequence from \code{ref.genome}, the function
#'   will automatically discard these variants. See element
#'   \code{discarded.variants} in the return value for more details.
#'
#' @param name.of.VCF Name of the VCF file.
#'
#' @param suppress.discarded.variants.warnings If TRUE, do warn when variants
#' that cannot be processed are discarded.
#'
#' @param explain_indels If TRUE generate message on stdout showing how the
#' indel was categorized.
#'
#' @param context_width_multiplier Used to guess how much sequence on each side
#' of an indel is needed to categorize it.
#'
#' @param add_transcript_ranges If TRUE add transcript ranges to the output VCF.
#
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
#'   frame with new columns added to the input data frame, including:
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
  function(
    ID.vcf,
    ref.genome,
    flag.mismatches = 0,
    name.of.VCF = NULL,
    suppress.discarded.variants.warnings = TRUE,
    explain_indels = 1,
    context_width_multiplier = 20L,
    add_transcript_ranges = TRUE
  ) {
    return_list = justify_id_vcf(
      ID.vcf = ID.vcf,
      ref.genome = ref.genome,
      name.of.VCF = name.of.VCF,
      suppress.discarded.variants.warnings,
      explain_indels = explain_indels,
      context_width_multiplier = context_width_multiplier
    )

    justified_vcf = return_list$annotated.vcf
    discarded_variants = return_list$discarded.variants

    df5 = justified_vcf
    if (add_transcript_ranges) {
      trans.ranges <- InferTransRanges(ref.genome)
      if (!is.null(trans.ranges)) {
        df5 <- AddTranscript(
          # overwrite previous df5
          df = justified_vcf,
          trans.ranges = trans.ranges,
          ref.genome = ref.genome,
          name.of.VCF = name.of.VCF
        )
      }
    }

    results_as_list = categorize_indels_in_vcf(
      df5
    )
    indel_info_df = data.table::rbindlist(results_as_list, fill = TRUE)

    df6 = cbind(data.table::as.data.table(df5), indel_info_df)

    return(list(
      annotated.vcf = df6,
      discarded.variants = discarded_variants
    ))
  }
