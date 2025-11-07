#' Print an indel with a dotted deletion marker.
#'
#' @param long_str Character string containing the pre-edit sequence.
#' @param short_str Character string containing the post-edit sequence.
#' @param pos One-based position of the edit within `long_str`.
#' @param ins_or_del Single character, `"i"` for insertion or `"d"` for deletion.
#'
#' @keywords internal
show_indel <- function(long_str, short_str, pos, ins_or_del) {
  stopifnot(
    is.character(long_str),
    length(long_str) == 1,
    is.character(short_str),
    length(short_str) == 1,
    is.numeric(pos),
    length(pos) == 1,
    is.character(ins_or_del),
    length(ins_or_del) == 1
  )

  ins_or_del <- match.arg(ins_or_del, c("i", "d"))
  pos_int <- as.integer(pos)
  if (pos_int != pos) {
    stop("`pos` must be an integer index.")
  }

  del_len <- nchar(long_str) - nchar(short_str)
  if (del_len < 0) {
    stop("`short_str` must be derived by deleting characters from `long_str`.")
  }

  prefix <- if (pos_int <= 1) "" else substr(short_str, 1, pos_int - 1)
  suffix <- if (pos_int <= nchar(short_str)) {
    substr(short_str, pos_int, nchar(short_str))
  } else {
    ""
  }

  expected_prefix <- if (pos_int <= 1) "" else substr(long_str, 1, pos_int - 1)
  expected_suffix_start <- pos_int + del_len
  expected_suffix <- if (expected_suffix_start > nchar(long_str)) {
    ""
  } else {
    substr(long_str, expected_suffix_start, nchar(long_str))
  }

  if (!identical(prefix, expected_prefix)) {
    stop(
      "`short_str` prefix does not align with `long_str` at position ",
      pos_int - 1,
      "."
    )
  }
  if (!identical(suffix, expected_suffix)) {
    stop("`short_str` suffix does not align with `long_str` after deletion.")
  }

  dotted_line <- paste0(prefix, strrep(".", del_len), suffix)

  if (ins_or_del == "d") {
    message(long_str)
    message(dotted_line)
  } else {
    message(dotted_line)
    message(long_str)
  }
}
