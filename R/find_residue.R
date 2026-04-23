#' Search for DNA or AA patterns in sequence data with auto-translation
#'
#' Main public function for unified pattern finding across DNAStringSet or
#' AAStringSet objects. Automatically handles translation when searching for
#' AA patterns in DNA sequences.
#'
#' @param xstring_set A `DNAStringSet` or `AAStringSet` containing sequences.
#' @param pattern Character vector of one or more patterns to search for.
#'   Can be nucleotide or amino acid codes.
#' @param start_position Optional integer >= 1. If provided and xstring_set
#'   is DNA, translation starts from this position. If NULL (default) and
#'   searching for AA pattern in DNA, auto-detects start position using
#'   detect_sequence_start().
#' @param verbose Logical; if TRUE (default), prints diagnostic messages
#'   when auto-translating or detecting start positions.
#'
#' @return A tibble with columns:
#'   - `sequence_name`: character, sequence identifier
#'   - `location`: integer, position of match within the (translated) sequence
#'   - `pattern`: character, the pattern that was matched
#'
#'   Rows are in long format (one row per match per pattern per sequence).
#'   Returns empty tibble with correct column structure if no matches found.
#'
#' @details
#' **Pattern Detection Logic:**
#' - If xstring_set is `AAStringSet`: searches for AA patterns as-is
#' - If xstring_set is `DNAStringSet` and all patterns are AA:
#'   - Auto-detects start position (or uses provided start_position)
#'   - Translates DNA to AA using translate_from_position()
#'   - Searches translated protein sequences
#' - If xstring_set is `DNAStringSet` and patterns are nucleotides:
#'   - Searches DNA directly
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   dna <- Biostrings::DNAStringSet(c(
#'     "Seq1" = "ATGAAATAG",
#'     "Seq2" = "ATGCCCTAA"
#'   ))
#'   # Find start codons
#'   find_residue(dna, pattern = "ATG")
#'
#'   # Auto-translate and find lysine (K)
#'   find_residue(dna, pattern = "K")
#' }
#'
#' @importFrom dplyr tibble bind_rows
#' @importFrom purrr map_df
#' @importFrom Biostrings DNAStringSet AAStringSet
#' @importFrom methods is
find_residue <- function(xstring_set,
                         pattern,
                         start_position = NULL,
                         verbose = TRUE) {

  # --- Input validation ------------------------------------------------------
  if (!is(xstring_set, "DNAStringSet") && !is(xstring_set, "AAStringSet")) {
    stop("xstring_set must be a DNAStringSet or AAStringSet")
  }

  if (!is.character(pattern) || length(pattern) == 0) {
    stop("pattern must be a non-empty character vector")
  }

  is_dna <- is(xstring_set, "DNAStringSet")
  is_protein <- is(xstring_set, "AAStringSet")

  # --- Determine search mode (DNA vs protein vs auto-translate) ----------------
  # Check if all patterns are AA patterns
  all_patterns_are_aa <- all(purrr::map_lgl(pattern, is_aa_pattern))

  if (is_protein) {
    # Protein input: search directly
    search_set <- xstring_set
    do_translate <- FALSE
  } else if (is_dna && all_patterns_are_aa) {
    # DNA input with AA patterns: translate
    do_translate <- TRUE

    # Determine start position
    if (is.null(start_position)) {
      # Auto-detect using modal start position
      start_info <- detect_sequence_start(
        xstring_set,
        pattern = "ATG",
        verbose = verbose
      )
      start_position <- start_info$position
    } else {
      if (!is.numeric(start_position) || start_position < 1) {
        stop("start_position must be a numeric >= 1")
      }
      if (verbose) {
        message("Translating from position ", start_position)
      }
    }

    # Translate all sequences from determined start position
    search_set <- translate_from_position(xstring_set, start_position = start_position)
  } else {
    # DNA input with nucleotide patterns: search DNA directly
    search_set <- xstring_set
    do_translate <- FALSE
  }

  # --- Search for patterns in (possibly translated) sequences ----------------
  seq_names <- names(search_set)

  # Build result rows for each sequence and each pattern
  result_rows <- purrr::map_df(
    seq_along(search_set),
    function(seq_idx) {
      seq <- search_set[[seq_idx]]
      seq_name <- seq_names[seq_idx]

      # For each pattern, find all matches
      pattern_rows <- purrr::map_df(
        pattern,
        function(pat) {
          positions <- search_pattern_in_xstring(seq, pat)

          if (length(positions) == 0L) {
            # No matches for this pattern in this sequence
            return(dplyr::tibble(
              sequence_name = character(),
              location = integer(),
              pattern = character()
            ))
          }

          # Create one row per match
          dplyr::tibble(
            sequence_name = rep(seq_name, length(positions)),
            location = positions,
            pattern = rep(pat, length(positions))
          )
        }
      )

      pattern_rows
    }
  )

  # --- Ensure correct column types and names even if empty -------------------
  if (nrow(result_rows) == 0L) {
    result_rows <- dplyr::tibble(
      sequence_name = character(),
      location = integer(),
      pattern = character()
    )
  }

  result_rows
}
