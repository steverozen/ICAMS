#' Convert an annotated indel VCF to a Koh 476-category catalog
#'
#' @description Take an annotated indel VCF data frame (with columns
#'   \code{Koh_476} and \code{R} as produced by indel classification
#'   functions) and produce a single-column matrix of mutation counts
#'   in the 476-category Koh classification scheme.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Optionally filters to PASS variants.
#'   \item Removes duplicate positions (warns if ALT alleles differ).
#'   \item Collapses single-base indels with repeat count \eqn{\ge 9}
#'     into an \code{"R(9,)"} bin.
#'   \item Tallies counts per Koh 476 category and returns a matrix
#'     with one row per category (using
#'     \code{ICAMS::catalog.row.order$ID476}).
#' }
#'
#' @param annot_vcf A data frame with at least columns
#'   \code{CHROM}, \code{POS}, \code{ALT}, \code{Koh_476}, and
#'   \code{R} (repeat count). If \code{FILTER_PASS} is \code{TRUE},
#'   a \code{FILTER} column is also required.
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
#' @return A single-column matrix with 476 rows (one per Koh
#'   category) and integer mutation counts. Row names are the Koh 476
#'   category strings; the column name is \code{sample_id}.
#'
#' @importFrom dplyr %>% pull mutate if_else
#'
#' @export
annot_vcf_to_476_catalog <- function(
  annot_vcf,
  sample_id = "no_sample_id_provided",
  FILTER_PASS = FALSE,
  do_message = FALSE
) {
  cleaner_vcf <- quick_check_vcf(annot_vcf, FILTER_PASS, do_message)

  cleaner_vcf %>%
    dplyr::mutate(
      Koh_476 = if_else(
        R >= 9 &
          stringr::str_detect(
            Koh_476,
            "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)"
          ),
        stringr::str_replace(Koh_476, "R\\d+", "R(9,)"),
        Koh_476
      )
    ) -> update_type_strings

  if (do_message) {
    message(
      "num PASS && unique rows after updating type string = ",
      nrow(update_type_strings)
    )
  }

  update_type_strings %>%
    dplyr::count(Koh_476) -> compacted_vcf

  if (do_message) {
    message("num PASS && unique mutations = ", sum(compacted_vcf$n))
  }

  data.table::data.table(Koh_476 = ICAMS::catalog.row.order$ID476) %>%
    dplyr::left_join(compacted_vcf) %>%
    mutate(n = if_else(is.na(n), 0L, n)) -> almost
  almost <- as.data.frame(almost)
  rownames(almost) <- pull(almost, Koh_476)
  colnames(almost)[2] <- sample_id

  almost[, -1, drop = FALSE]
}
