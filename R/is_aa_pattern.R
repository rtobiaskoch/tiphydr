#' Detect if a pattern is amino acids or nucleotides
#'
#' @description
#' A pure helper function that determines whether a character pattern
#' represents amino acid codes or nucleotide codes. For ambiguous patterns
#' that contain characters valid in both alphabets (e.g., 'A'), prioritizes
#' amino acid classification.
#'
#' @param pattern A character string representing a sequence pattern
#'
#' @return Logical. TRUE if pattern contains only standard amino acid codes,
#'   FALSE if pattern contains nucleotides, invalid characters, or is empty.
#'
#' @details
#' Standard amino acid codes recognized: ACDEFGHIKLMNPQRSTVWY
#' The function is case-insensitive.
#'
#' For ambiguous characters (e.g., 'A' is both AA and nucleotide),
#' amino acid classification is prioritized, returning TRUE.
#'
#' @examples
#' is_aa_pattern("K")  # TRUE
#' is_aa_pattern("ATG")  # FALSE (nucleotide start codon)
#' is_aa_pattern("A")  # TRUE (prioritizes AA)
#' is_aa_pattern("XYZ")  # FALSE (invalid characters)
#'
#' @export
is_aa_pattern <- function(pattern) {
  # Validate input: must be character and non-empty
  if (!is.character(pattern) || length(pattern) == 0 || nchar(pattern) == 0) {
    return(FALSE)
  }

  # Standard amino acid codes (uppercase)
  aa_codes <- "ACDEFGHIKLMNPQRSTVWY"
  # Standard nucleotide codes (uppercase) - T, U, G, C, A
  nt_codes <- "ACGTU"

  # Convert to uppercase for comparison
  pattern_upper <- toupper(pattern)

  # Split pattern into individual characters
  pattern_chars <- unlist(strsplit(pattern_upper, ""))

  # Check if all characters are valid in either alphabet
  aa_chars <- unlist(strsplit(aa_codes, ""))
  nt_chars <- unlist(strsplit(nt_codes, ""))

  # All chars must be valid (in either AA or NT alphabet)
  if (!all(pattern_chars %in% union(aa_chars, nt_chars))) {
    return(FALSE)
  }

  # Check if there are any AA-only characters (not in NT alphabet)
  # AA-only chars are: DEFHIKLMNPQRSVWY (excludes A, C, G, T, U)
  aa_only_chars <- setdiff(aa_chars, nt_chars)
  has_aa_only <- any(pattern_chars %in% aa_only_chars)

  # Check if there are any NT-only characters (only U, T, or G when not in other contexts)
  # NT-only chars are: T, U (G is tricky - it's in both, but if we see T or U it's likely NT)
  nt_only_chars <- c("T", "U")
  has_nt_only <- any(pattern_chars %in% nt_only_chars)

  # If has AA-only chars, it's definitely AA
  if (has_aa_only) {
    return(TRUE)
  }

  # If has NT-only chars (T or U), it's definitely NT
  if (has_nt_only) {
    return(FALSE)
  }

  # If only ambiguous chars (A, C, G) with no AA-only or NT-only markers,
  # prioritize AA classification (return TRUE)
  return(TRUE)
}
