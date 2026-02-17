#' Convert an annotated indel VCF to a Koh 89-category catalog
#'
#' @description Take an annotated indel VCF data frame (with column
#'   \code{Koh_89} as produced by indel classification functions) and
#'   produce a single-column matrix of mutation counts in the
#'   89-category Koh classification scheme.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Optionally filters to PASS variants.
#'   \item Removes duplicate positions (warns if ALT alleles differ).
#'   \item Tallies counts per Koh 89 category and returns a matrix
#'     with one row per category (using
#'     \code{ICAMS::catalog.row.order$ID89}).
#' }
#'
#' @param annot_vcf A data frame with at least columns
#'   \code{CHROM}, \code{POS}, \code{ALT}, and \code{Koh_89}.
#'   If \code{FILTER_PASS} is \code{TRUE}, a \code{FILTER} column
#'   is also required.
#'
#' @param sample_id A character string used as the column name in the
#'   returned matrix.
#'
#' @param FILTER_PASS If \code{TRUE}, retain only rows where the
#'   \code{FILTER} column equals \code{"PASS"}.
#'
#' @param do_message If \code{TRUE}, emit diagnostic messages showing
#'   row counts at each processing step.
#'
#' @return A single-column matrix with 89 rows (one per Koh category)
#'   and integer mutation counts. Row names are the Koh 89 category
#'   strings; the column name is \code{sample_id}.
#'
#' @importFrom dplyr %>% filter group_by pull mutate if_else
#'
#' @export
annot_vcf_to_89_catalog <- function(
  annot_vcf,
  sample_id = "no_sample_id_provided",
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
    pull(pos_id) %>%
    unique() -> multiple_alts

  if (length(multiple_alts) > 0) {
    warning(
      "Differences in 'ALT'; only 1 ALT value chosen arbitrarily at the following positions: ",
      paste(capture.output(print(multiple_alts)), collapse = '\n')
    )
  }

  vcf_with_pos_id %>% dplyr::distinct(pos_id, .keep_all = TRUE) -> cleaner_vcf
  if (do_message) {
    message("num PASS && unique rows = ", nrow(cleaner_vcf))
  }

  cleaner_vcf %>%
    dplyr::count(Koh_89) -> compacted_vcf

  if (do_message) {
    message("num PASS && unique mutations = ", sum(compacted_vcf$n))
  }

  data.table::data.table(Koh_89 = ICAMS::catalog.row.order$ID89) %>%
    dplyr::left_join(compacted_vcf) %>%
    mutate(n = if_else(is.na(n), 0L, n)) -> almost
  almost <- as.data.frame(almost)
  rownames(almost) <- pull(almost, Koh_89)
  colnames(almost)[2] <- sample_id

  almost[, -1, drop = FALSE]
}
