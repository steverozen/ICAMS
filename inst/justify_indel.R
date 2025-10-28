# justify_indel

# To test the code in this file try test_VCFsToCatalogs.R

# Take the view that there is a deletion in long_str
# at pos that creates short_str.
#
# pos is the 1-based position in long string at which to make a deletion to get short_str
# dually, pos is the 1-based position in short string before which one can
# make an insertion to get long_str.
#
# Move pos as far to left as possible so that a deletion
# at that position still results in an edit of long_str
# to short_str
#
# This can also be interpreted an inserion in short_str
# immediately in front of pos that generates long_str

justify_indel = function(long_str, short_str, pos, expected_delta = NULL) {
  library(stringi)

  del_len = nchar(long_str) - nchar(short_str)

  if (!is.null(expected_delta)) {
    if (expected_delta != substr(long_str, pos, pos + del_len - 1)) {
      message(
        "expected_delta = ",
        expected_delta,
        " != substr(long_str, pos, pos + del_len - 1)) = ",
        substr(long_str, pos, pos + del_len - 1)
      )
      browser()
      stop()
    }
  }

  stopifnot(
    stri_sub_replace(
      long_str,
      from = pos,
      to = pos + del_len - 1,
      replacement = ""
    ) ==
      short_str
  )
  for (test_pos in (pos - 1):0) {
    if (
      stri_sub_replace(
        long_str,
        test_pos,
        test_pos + del_len - 1,
        replacement = ""
      ) !=
        short_str
    ) {
      break
    }
  }
  leftmost_pos = test_pos + 1
  return(
    list(
      leftmost_pos = leftmost_pos,
      del_str = stri_sub(
        long_str,
        leftmost_pos,
        leftmost_pos + del_len - 1
      )
    )
  )
}

if (FALSE) {
  justify_indel("abc", "ac", 2)
  justify_indel("abbbc", "abbc", 4)
  justify_indel("xycagcaguv", "xycaguv", 4)
  justify_indel("xycagcaguv", "xycaguv", 5)
  justify_indel("xycagcaguv", "xycaguv", 6)
  justify_indel("xycagcaguv", "xycaguv", 7) # error
  justify_indel("xycagcagcagcaguv", "xycagcagcaguv", 9)
}
