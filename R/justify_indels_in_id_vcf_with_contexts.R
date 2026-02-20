#' @title Justify indels in a VCF and update positions accordingly
#'
#' @description
#' For each indel in the input VCF, this function:
#' \enumerate{
#'   \item Justifies the indel position (moves it as far left as possible)
#'   \item Calculates how much the position was shifted
#'   \item Updates the POS column by decrementing it by the shift amount
#'   \item Updates the seq.context.width column accordingly
#'   \item Adds a new column 'pos_shift' showing how much each position was moved
#'   \item Categorizes the justified indel
#' }
#'
#' @section Important invariant:
#' The seq.context string is never modified. Both POS and seq.context.width are
#' decremented by the same amount (pos_shift), which maintains the relationship:
#' genomic position POS corresponds to position (seq.context.width + 1) in seq.context.
#'
#' For example, if seq.context was extracted from genomic positions
#' (POS - seq.context.width, POS + var.width + seq.context.width), then after
#' justification by shift S:
#' \itemize{
#'   \item seq.context remains unchanged
#'   \item POS becomes POS - S
#'   \item seq.context.width becomes seq.context.width - S
#'   \item The variant is now at position (seq.context.width - S) + 1 in seq.context
#' }
#'
#' @param vcf A data.frame representing a VCF, which must contain the following columns:
#'   \itemize{
#'     \item seq.context: Ample surrounding sequence on each side of the variants
#'     \item REF: The reference alleles (includes one unaltered base at the start)
#'     \item ALT: The alternative alleles
#'     \item POS: The genomic positions
#'     \item seq.context.width: The width of seq.context to the left
#'   }
#'
#' @param explain_indels 0 stay silent,
#' if 1, print explanation of there was a change in POS,
#' if 3, aways print
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item All original columns from the input vcf
#'     \item Updated POS values (decremented by pos_shift)
#'     \item Updated REF and ALT values if the positon of the indel has been moved
#'     \item Updated seq.context.width values (decremented by pos_shift)
#'     \item New column 'pos_shift': Amount the position was moved (usually 0)
#'   }
#'
#' @importFrom stringi stri_sub_replace
#'
#' @importFrom data.table rbindlist as.data.table
#'
#' @keywords internal
#'
justify_indels_in_id_vcf_with_contexts <- function(vcf, explain_indels = 1) {
  requireNamespace("stringi")

  if (nrow(vcf) == 0) {
    return(vcf)
  }

  # Verify required columns
  required_cols <- c("seq.context", "REF", "ALT", "POS", "seq.context.width")
  missing_cols <- setdiff(required_cols, colnames(vcf))
  if (length(missing_cols) > 0) {
    stop("VCF missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Initialize column for position shift
  vcf$pos_shift <- 0L

  # Process each row

  for (i in seq_len(nrow(vcf))) {
    # Extract information for this indel
    orig_ref <- vcf$REF[i]
    orig_alt <- vcf$ALT[i]
    seq_context <- vcf$seq.context[i]
    orig_pos <- vcf$seq.context.width[i] + 1 # Position in seq.context

    # Remove the common prefix base
    # REF and ALT are expected to have the same first base (the context base)
    if (!(nchar(orig_ref) == 1 || nchar(orig_alt) == 1)) {
      stop("We should not take the substr's below")
    }
    ref <- substr(orig_ref, 2, nchar(orig_ref))
    alt <- substr(orig_alt, 2, nchar(orig_alt))
    pos <- orig_pos + 1 # Position after removing common base

    # Determine if insertion or deletion and justify
    if (nchar(alt) < nchar(ref)) {
      # Deletion
      tmp_long <- seq_context
      end_pos <- pos + nchar(ref) - 1
      tmp_short <- stringi::stri_sub_replace(
        seq_context,
        from = pos,
        to = end_pos,
        replacement = ""
      )

      tmp <- justify_indel(tmp_long, tmp_short, pos, ref)
      new_ref = tmp$del_str_plus
      new_alt = stringi::stri_sub(new_ref, 1, 1)
      if (tmp$edge_warning) {
        warning(
          "There was not enough context to completly justify the deletion at ",
          vcf$CHOM[i],
          vcf$POS[i]
        )
      }

      # Calculate position shift (how much we moved left)
      pos_shift <- pos - tmp$leftmost_pos
    } else if (nchar(alt) > nchar(ref)) {
      # Insertion
      tmp_short <- seq_context
      tmp_long <- stringi::stri_sub_replace(
        tmp_short,
        pos,
        pos - 1,
        replacement = alt
      )

      tmp <- justify_indel(tmp_long, tmp_short, pos, alt)
      new_alt = tmp$del_str_plus
      new_ref = stringi::stri_sub(new_alt, 1, 1)
      if (tmp$edge_warning) {
        warning(
          "There was not enough context to completly justify the insertion at ",
          vcf$CHOM[i],
          vcf$POS[i]
        )
      }

      # Calculate position shift (how much we moved left)
      pos_shift <- pos - tmp$leftmost_pos

      # Categorize using the justified position
    } else {
      # Should not happen (same length) - already filtered out upstream
      result <- indel_all_na_return("Non-indel")
      pos_shift <- 0L
    }

    if (
      explain_indels == 2 ||
        (explain_indels == 1 && pos_shift > 0) ||
        tmp$edge_warning
    ) {
      explain_indel_justification(
        long_str = tmp_long,
        short_str = tmp_short,
        pos_before = pos,
        pos_after = tmp$leftmost_pos,
        is_ins = nchar(alt) > nchar(ref)
      )
      if (tmp$edge_warning) {
        message(
          "!!! We are not sure how far to the left we could ",
          " have moved the indel if we had more context ",
          "on the left"
        )
      }
    }

    # Store the position shift (use set() to avoid data.table column copy)
    data.table::set(vcf, i = as.integer(i), j = "pos_shift", value = as.integer(pos_shift))

    # Update POS: decrement by the shift amount
    # When we justify an indel and move it left by pos_shift bases,
    # the genomic position also moves left (decreases)
    data.table::set(vcf, i = as.integer(i), j = "POS", value = vcf$POS[i] - pos_shift)

    if (pos_shift == 0) {
      stopifnot(orig_alt == new_alt)
      stopifnot(orig_ref == new_ref)
    } else {
      data.table::set(vcf, i = as.integer(i), j = "ALT", value = new_alt)
      data.table::set(vcf, i = as.integer(i), j = "REF", value = new_ref)
    }
    # Update seq.context.width: decrement by the shift amount
    # The justified position is now closer to the left edge of seq.context
    data.table::set(vcf, i = as.integer(i), j = "seq.context.width", value = vcf$seq.context.width[i] - pos_shift)

    if (explain_indels > 1) {
      message(sprintf(
        "Row %d: shifted by %d bases (POS: %d -> %d)",
        i,
        pos_shift,
        vcf$POS[i] + pos_shift,
        vcf$POS[i]
      ))
    }
  }

  return(vcf)
}
