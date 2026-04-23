#' Find pattern match positions in a single XString
#'
#' Generic wrapper for Biostrings::matchPattern() that extracts start positions
#' of all matches for a pattern in a single sequence (DNAString, AAString, etc).
#' Returns an integer vector of positions (1-indexed).
#'
#' @param xstring A Biostrings::XString object (DNAString, AAString, etc)
#' @param pattern A character string representing the pattern to search for.
#'   For AAString, pattern is automatically converted to uppercase.
#'
#' @return Integer vector of start positions of matches.
#'   Returns empty integer(0) if no matches found.
#'
#' @details
#' Uses Biostrings::matchPattern() to find all non-overlapping matches.
#' For amino acids, the pattern is converted to uppercase before matching
#' to ensure case-insensitive matching.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   dna <- Biostrings::DNAString("ATGATGATG")
#'   positions <- search_pattern_in_xstring(dna, "ATG")
#'   # Result: c(1, 4, 7)
#' }
#'
#' @importFrom Biostrings matchPattern start
search_pattern_in_xstring <- function(xstring, pattern) {
  # --- Input validation ------------------------------------------------------
  # Ensure pattern is character
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("pattern must be a character string of length 1")
  }

  # Check if this is an AAString (amino acids)
  is_aa <- is(xstring, "AAString")

  # Convert pattern to uppercase if we're searching an AAString
  if (is_aa) {
    pattern <- toupper(pattern)
  }

  # --- Find matches using Biostrings ----------------------------------------
  # matchPattern returns a Views object containing match ranges
  matches <- Biostrings::matchPattern(pattern, xstring)

  # Extract start positions from the matches
  # If no matches, this returns integer(0)
  if (length(matches) == 0L) {
    return(integer(0))
  }

  # Get start position for each match
  positions <- Biostrings::start(matches)

  positions
}
