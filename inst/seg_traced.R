#' Segment a single indel sequence with detailed tracing
#'
#' This is a wrapper around the standalone segment_traced C++ program that
#' processes a single indel-context pair. Use this for detailed analysis of
#' individual sequences, or use \code{segment_traced()} for vectorized processing.
#'
#' @param string A single indel sequence to segment (character scalar)
#' @param context A single flanking context sequence (character scalar)
#' @param trace Logical. If TRUE, print detailed trace output to console
#' @param binary_path Path to the segment_traced binary. Defaults to
#'   "./segment_traced" in the current working directory.
#'
#' @return A list with the following elements:
#'   \item{unit}{Repeat unit sequence (character)}
#'   \item{unit_length}{Unit length (integer)}
#'   \item{internal_rep}{Internal repeat region (character)}
#'   \item{internal_reps}{Internal repeat count (integer)}
#'   \item{spacer}{Spacer sequence (character)}
#'   \item{spacer_length}{Spacer length (integer)}
#'   \item{prime3_rep}{3' flanking repeat region (character)}
#'   \item{prime3_reps}{3' flanking repeat count (integer)}
#'
#' @details
#' This function requires the segment_traced binary to be compiled. To compile:
#' \code{g++ -std=c++11 -o segment_traced segment_traced.cpp}
#'
#' The trace output shows the step-by-step decision-making process including:
#' \itemize{
#'   \item Each candidate repeat unit tested
#'   \item LTRS (Longest Tandem Repeat Search) results
#'   \item Scoring for each segmentation
#'   \item Comparison logic between segmentations
#'   \item Final best segmentation selection
#' }
#'
#' @examples
#' \dontrun{
#' # Simple AT repeat with trace
#' result <- seg_traced("ATATAT", "ATATGG", trace = TRUE)
#' }
#'
#' @export
seg_traced <- function(
  ins_or_del,
  string,
  context,
  trace = FALSE,
  binary_path = devtools::package_file("segment_traced")
) {
  # Input validation
  if (length(string) != 1 || length(context) != 1) {
    stop("string and context must each be a single character value")
  }

  if (!file.exists(binary_path)) {
    stop(paste0(
      "segment_traced binary not found at: ",
      binary_path,
      "\n",
      "Compile with: g++ -std=c++11 -o segment_traced segment_traced.cpp"
    ))
  }

  # Build command
  args <- c(string, context)
  if (trace) {
    args <- c(args, "--trace")
  }

  # Call the binary and capture output
  output <- system2(binary_path, args = args, stdout = TRUE, stderr = TRUE)

  # Parse the output to extract numeric results
  # Look for the summary section which contains the key values
  unit_length_line <- grep("Unit length:", output, value = TRUE)
  internal_bases_line <- grep("Internal repeats:", output, value = TRUE)
  spacer_len_line <- grep("Spacer:", output, value = TRUE)
  prime3_bases_line <- grep("3' context repeats:", output, value = TRUE)

  # Extract numeric values
  unit_len <- as.integer(sub(
    ".*Unit length: ([0-9]+).*",
    "\\1",
    unit_length_line[1]
  ))
  internal_bases <- as.integer(sub(
    ".*Internal repeats: ([0-9]+) bases.*",
    "\\1",
    internal_bases_line[1]
  ))
  spacer_len <- as.integer(sub(
    ".*Spacer: ([0-9]+) bases.*",
    "\\1",
    spacer_len_line[1]
  ))
  prime3_bases <- as.integer(sub(
    ".*3' context repeats: ([0-9]+) bases.*",
    "\\1",
    prime3_bases_line[1]
  ))

  # Extract string components (same logic as code.cpp)
  # unit = string[0:unit_length]
  unit <- substr(string, 1, unit_len)

  # internal_rep = string[unit_length : unit_length + internal_bases]
  if (internal_bases > 0) {
    internal_rep <- substr(string, unit_len + 1, unit_len + internal_bases)
  } else {
    internal_rep <- ""
  }

  # spacer = string[unit_length + internal_bases : unit_length + internal_bases + spacer_len]
  if (spacer_len > 0) {
    spacer_start <- unit_len + internal_bases + 1
    spacer <- substr(string, spacer_start, spacer_start + spacer_len - 1)
  } else {
    spacer <- ""
  }

  # prime3_rep = context[0:prime3_bases]
  if (prime3_bases > 0) {
    prime3_rep <- substr(context, 1, prime3_bases)
  } else {
    prime3_rep <- ""
  }

  # Convert bases to repeat counts (divide by unit length)
  # This matches the logic in code.cpp lines 103-106
  if (unit_len != 0) {
    internal_reps <- internal_bases %/% unit_len
    prime3_reps <- prime3_bases %/% unit_len
  } else {
    internal_reps <- 0
    prime3_reps <- 0
  }

  original_reps = ifelse(
    ins_or_del == "i",
    prime3_reps,
    ifelse(
      spacer_len == 0,
      internal_reps + 1 + prime3_reps,
      prime3_reps
    )
  )

  # Return list matching the structure of the original segment() function
  return(list(
    unit = unit,
    unit_length = unit_len,
    internal_rep = internal_rep,
    internal_reps = internal_reps,
    spacer = spacer,
    spacer_length = spacer_len,
    prime3_rep = prime3_rep,
    prime3_reps = prime3_reps,
    original_reps = original_reps
  ))
}
