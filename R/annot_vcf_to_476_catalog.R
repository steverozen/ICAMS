#' Convert an annotated indel VCF to a Koh 476-category catalog
#'
#' @description Take an annotated indel VCF data frame (with columns
#'   \code{Koh_476} and \code{R} as produced by indel classification
#'   functions) and produce a single-column data frame of mutation counts
#'   in the 476-category Koh classification scheme.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Optionally filters to PASS variants.
#'   \item Removes duplicate positions (warns if ALT alleles differ).
#'   \item Collapses single-base indels with repeat count \eqn{\ge 9}
#'     into an \code{"R(9,)"} bin.
#'   \item Tallies counts per Koh 476 category and returns a data frame
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
#' @return A single-column data frame with 476 rows (one per Koh
#'   category) and integer mutation counts. Row names are the Koh 476
#'   category strings; the column name is \code{sample_id}.
#'
#' @export
annot_vcf_to_476_catalog <- function(
  annot_vcf,
  sample_id = "no_sample_id_provided",
  FILTER_PASS = FALSE,
  do_message = FALSE,
  clip_le_9 = FALSE
) {
  zero_catalog <- function() {
    rn <- ICAMS::catalog.row.order$ID476
    m <- data.frame(x = rep(0L, length(rn)), row.names = rn)
    colnames(m) <- sample_id
    m
  }

  if (nrow(annot_vcf) == 0) return(zero_catalog())

  cleaner_vcf <- quick_check_vcf(annot_vcf, FILTER_PASS, do_message)

  if (nrow(cleaner_vcf) == 0) return(zero_catalog())

  # Replaced dplyr with data.table for performance:
  # cleaner_vcf %>%
  #   dplyr::mutate(
  #     Koh_476 = if_else(
  #       R >= 9 &
  #         stringr::str_detect(
  #           Koh_476,
  #           "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)"
  #         ),
  #       stringr::str_replace(Koh_476, "R\\d+", "R(9,)"),
  #       Koh_476
  #     )
  #   ) -> update_type_strings
  dt <- data.table::as.data.table(cleaner_vcf)
  mask <- dt$R >= 9 &
    stringr::str_detect(
      dt$Koh_476,
      "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)"
    )
  data.table::set(
    dt, which(mask), "Koh_476",
    stringr::str_replace(dt$Koh_476[mask], "R\\d+", "R(9,)")
  )
  update_type_strings <- dt

  if (do_message) {
    message(
      "num PASS && unique rows after updating type string = ",
      nrow(update_type_strings)
    )
  }

  if (nrow(update_type_strings) == 0) return(zero_catalog())

  # Replaced dplyr with data.table for performance:
  # if (clip_le_9) {
  #   update_type_strings <- dplyr::filter(update_type_strings, R <= 9)
  # }
  if (clip_le_9) {
    update_type_strings <- update_type_strings[R <= 9]
    if (do_message) {
      message(
        "num rows after R <= 9 filter = ",
        nrow(update_type_strings)
      )
    }
  }

  # Replaced dplyr with data.table for performance:
  # update_type_strings %>%
  #   dplyr::count(Koh_476) -> compacted_vcf
  compacted_vcf <- update_type_strings[, .N, by = Koh_476]

  if (do_message) {
    message("num PASS && unique mutations = ", sum(compacted_vcf$N))
  }

  # Replaced dplyr with data.table for performance:
  # data.table::data.table(Koh_476 = ICAMS::catalog.row.order$ID476) %>%
  #   dplyr::left_join(compacted_vcf, by = "Koh_476") %>%
  #   mutate(n = if_else(is.na(n), 0L, n)) -> almost
  all_cats <- data.table::data.table(Koh_476 = ICAMS::catalog.row.order$ID476)
  almost <- compacted_vcf[all_cats, on = "Koh_476"]
  data.table::setnafill(almost, fill = 0L, cols = "N")

  almost <- as.data.frame(almost)
  rownames(almost) <- almost$Koh_476
  colnames(almost)[colnames(almost) == "N"] <- sample_id

  almost[, "Koh_476" != colnames(almost), drop = FALSE]
}
