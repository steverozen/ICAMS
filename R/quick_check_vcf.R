#' Filter, deduplicate, and check an annotated indel VCF
#'
#' Helper used by \code{\link{annot_vcf_to_83_catalog}},
#' \code{\link{annot_vcf_to_89_catalog}}, and
#' \code{\link{annot_vcf_to_476_catalog}} to perform shared
#' preprocessing: rename \code{#CHROM} to \code{CHROM}, optionally
#' filter to PASS variants, warn about positions with differing ALT
#' alleles, and deduplicate by position.
#'
#' @param annot_vcf A data frame with at least columns
#'   \code{CHROM} (or \code{#CHROM}), \code{POS}, and \code{ALT}.
#'   If \code{FILTER_PASS} is \code{TRUE}, a \code{FILTER} column
#'   is also required.
#'
#' @param FILTER_PASS If \code{TRUE}, retain only rows where the
#'   \code{FILTER} column equals \code{"PASS"}.
#'
#' @param do_message If \code{TRUE}, emit diagnostic messages showing
#'   row counts at each processing step.
#'
#' @return A data frame deduplicated by position, with a \code{pos_id}
#'   column added.
#'
#' @importFrom dplyr %>% filter group_by select
#'
#' @keywords internal
quick_check_vcf <- function(
  annot_vcf,
  FILTER_PASS = FALSE,
  do_message = FALSE
) {
  if (colnames(annot_vcf)[1] == "#CHROM") {
    colnames(annot_vcf)[1] <- "CHROM"
  }

  if (do_message) {
    message("initial annot_vcf rows = ", nrow(annot_vcf))
  }

  if (FILTER_PASS) {
    annot_vcf %>%
      filter(FILTER == "PASS") -> annot_vcf
  }

  if (do_message) {
    message("num PASS rows = ", nrow(annot_vcf))
  }

  annot_vcf %>%
    dplyr::mutate(pos_id = paste0(CHROM, "-", POS)) -> vcf_with_pos_id

  vcf_with_pos_id %>%
    group_by(pos_id) %>%
    filter(dplyr::n_distinct(ALT) > 1) %>%
    dplyr::select(pos_id, REF, ALT) -> multiple_alts

  if (nrow(multiple_alts) > 0) {
    warning(
      "Differences in 'ALT'; only 1 ALT value chosen arbitrarily at the following positions: ",
      paste(capture.output(print(multiple_alts)), collapse = '\n')
    )
  }

  vcf_with_pos_id %>% dplyr::distinct(pos_id, .keep_all = TRUE) -> cleaner_vcf
  if (do_message) {
    message("num PASS && unique rows = ", nrow(cleaner_vcf))
  }

  cleaner_vcf
}
