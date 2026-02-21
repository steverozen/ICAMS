#' Convert an annotated indel VCF to a Koh 89-category catalog
#'
#' @description Take an annotated indel VCF data frame (with column
#'   \code{Koh_89} as produced by indel classification functions) and
#'   produce a single-column data frame of mutation counts in the
#'   89-category Koh classification scheme.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Optionally filters to PASS variants.
#'   \item Removes duplicate positions (warns if ALT alleles differ).
#'   \item Tallies counts per Koh 89 category and returns a data frame
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
#'   returned data frame.
#'
#' @param FILTER_PASS If \code{TRUE}, retain only rows where the
#'   \code{FILTER} column equals \code{"PASS"}.
#'
#' @param clip_le_9 Only keep variants with "R" <= 9, to approximate
#'   PCAWG indel calling.
#'
#' @param do_message If \code{TRUE}, emit diagnostic messages showing
#'   row counts at each processing step.
#'
#' @return A single-column data frame with 89 rows (one per Koh category)
#'   and integer mutation counts. Row names are the Koh 89 category
#'   strings; the column name is \code{sample_id}.
#'
#' @export
annot_vcf_to_89_catalog <- function(
  annot_vcf,
  sample_id = "no_sample_id_provided",
  FILTER_PASS = FALSE,
  do_message = FALSE,
  clip_le_9 = FALSE
) {
  zero_catalog <- function() {
    rn <- ICAMS::catalog.row.order$ID89
    m <- data.frame(x = rep(0L, length(rn)), row.names = rn)
    colnames(m) <- sample_id
    m
  }

  if (nrow(annot_vcf) == 0) return(zero_catalog())

  cleaner_vcf <- quick_check_vcf(annot_vcf, FILTER_PASS, do_message)

  if (nrow(cleaner_vcf) == 0) return(zero_catalog())

  # Replaced dplyr with data.table for performance:
  # if (clip_le_9) {
  #   cleaner_vcf <- dplyr::filter(cleaner_vcf, R <= 9)
  # }
  if (clip_le_9) {
    dt <- data.table::as.data.table(cleaner_vcf)
    cleaner_vcf <- as.data.frame(dt[R <= 9])
    if (do_message) {
      message("num rows after R <= 9 filter = ", nrow(cleaner_vcf))
    }
  }

  # Replaced dplyr with data.table for performance:
  # cleaner_vcf %>%
  #   dplyr::count(Koh_89) -> compacted_vcf
  dt <- data.table::as.data.table(cleaner_vcf)
  compacted_vcf <- dt[, .N, by = Koh_89]

  if (do_message) {
    message("num PASS && unique mutations = ", sum(compacted_vcf$N))
  }

  # Replaced dplyr with data.table for performance:
  # data.table::data.table(Koh_89 = ICAMS::catalog.row.order$ID89) %>%
  #   dplyr::left_join(compacted_vcf, by = "Koh_89") %>%
  #   mutate(n = if_else(is.na(n), 0L, n)) -> almost
  all_cats <- data.table::data.table(Koh_89 = ICAMS::catalog.row.order$ID89)
  almost <- compacted_vcf[all_cats, on = "Koh_89"]
  data.table::setnafill(almost, fill = 0L, cols = "N")

  almost <- as.data.frame(almost)
  rownames(almost) <- almost$Koh_89
  colnames(almost)[colnames(almost) == "N"] <- sample_id

  almost[, "Koh_89" != colnames(almost), drop = FALSE]
}
