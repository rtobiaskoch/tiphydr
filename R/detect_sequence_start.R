#' Find the modal start position of a pattern across all sequences
#'
#' Searches for the most common (modal) start position of a pattern
#' (typically "ATG" start codon) across a DNAStringSet. Useful for
#' determining the consensus CDS start position before translation.
#'
#' @param xstring_set A `DNAStringSet` containing sequences to analyze.
#' @param pattern A character string representing the pattern to search for.
#'   Default is "ATG" (start codon).
#' @param verbose Logical; if TRUE (default), prints a summary message
#'   showing the detected position and counts of sequences with/without pattern.
#'
#' @return A list with two components:
#'   - `position`: integer, the modal (most common) start position
#'   - `table`: table object, frequency table of all positions found
#'
#' @details
#' For each sequence in the DNAStringSet, finds the first occurrence of the
#' pattern using search_pattern_in_xstring(). Sequences without the pattern
#' are excluded from the calculation. The modal position is the most frequently
#' occurring position.
#'
#' If no sequences contain the pattern, an error is raised.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   dna <- Biostrings::DNAStringSet(c(
#'     "Seq1" = "ATGAAATAG",
#'     "Seq2" = "ATGCCCTAA"
#'   ))
#'   result <- detect_sequence_start(dna, pattern = "ATG", verbose = TRUE)
#'   # position: 1
#' }
#'
#' @importFrom Biostrings DNAStringSet
detect_sequence_start <- function(xstring_set, pattern = "ATG", verbose = TRUE) {
  # --- Input validation ------------------------------------------------------
  if (!is(xstring_set, "DNAStringSet")) {
    stop("xstring_set must be a DNAStringSet")
  }

  if (!is.character(pattern) || length(pattern) != 1 || nchar(pattern) == 0) {
    stop("pattern must be a non-empty character string")
  }

  # --- Find first match position in each sequence ---------------------------
  # Loop through sequences and collect first match position for each
  first_positions <- purrr::map_int(
    seq_along(xstring_set),
    function(i) {
      positions <- search_pattern_in_xstring(xstring_set[[i]], pattern)
      # Return first position if found, NA otherwise
      if (length(positions) > 0) positions[1] else NA_integer_
    }
  )

  # --- Exclude sequences missing the pattern --------------------------------
  # Get only the positions that were found (not NA)
  found_positions <- first_positions[!is.na(first_positions)]

  # If no sequences contain the pattern, error
  if (length(found_positions) == 0L) {
    stop("No sequences contain pattern '", pattern, "'")
  }

  # --- Calculate mode of valid positions ------------------------------------
  # Create a frequency table
  position_table <- table(found_positions)

  # Find the modal position (the one with highest frequency)
  modal_position <- as.integer(names(position_table)[which.max(position_table)])

  # --- Prepare output -------------------------------------------------------
  result <- list(
    position = modal_position,
    table = position_table
  )

  # --- Optional verbose message -----------------------------------------------
  if (verbose) {
    n_found <- length(found_positions)
    n_missing <- sum(is.na(first_positions))
    message(
      "Start codon '", pattern, "' detected at position ",
      modal_position, " (found in ", n_found, " sequences",
      if (n_missing > 0) paste0("; missing in ", n_missing, " sequences") else "",
      ")"
    )
  }

  result
}
