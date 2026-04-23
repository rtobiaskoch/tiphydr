#' Translate DNA to amino acids starting from a specified position
#'
#' Extracts a subsequence from a DNAStringSet starting at a given position,
#' then translates to amino acids. Fuzzy nucleotides (e.g., N, Y, R) in codons
#' are represented as 'X' (unknown amino acid) in the output.
#'
#' @param xstring_set A `DNAStringSet` object containing sequences to translate.
#' @param start_position An integer >= 1 specifying the nucleotide position
#'   (1-indexed) to start translation. Defaults to 1 (start of sequence).
#'
#' @return An `AAStringSet` with the same names as `xstring_set`, containing
#'   the translated amino acid sequences.
#'
#' @examples
#' \dontrun{
#'   dna <- Biostrings::DNAStringSet(
#'     c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA")
#'   )
#'   aa <- translate_from_position(dna, start_position = 1)
#'   # Result: AAStringSet with sequences "MK*" and "MP*"
#' }
#'
#' @export
translate_from_position <- function(xstring_set, start_position = 1) {
  # --- Input validation ------------------------------------------------------
  if (!is(xstring_set, "DNAStringSet")) {
    stop("xstring_set must be a DNAStringSet")
  }

  if (!is.integer(start_position) && !is.numeric(start_position)) {
    stop("start_position must be numeric")
  }

  if (start_position < 1) {
    stop("start_position must be >= 1")
  }

  # --- Extract and translate -------------------------------------------------
  # Extract subsequence from start_position to the end of each sequence
  subseqs <- Biostrings::subseq(xstring_set, start = start_position)

  # Translate to amino acids with fuzzy codon handling
  # if.fuzzy.codon = "X" converts any codon with ambiguous nucleotides to "X"
  translated <- Biostrings::translate(subseqs, if.fuzzy.codon = "X")

  translated
}
