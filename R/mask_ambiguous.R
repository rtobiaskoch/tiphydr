#' Mask ambiguous characters in a DNA sequence set
#'
#' Replaces every character that is not an unambiguous base (A/C/G/T by default)
#' with a single replacement character (N by default). This includes IUPAC
#' ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) and gap characters
#' (`-`, `.`).
#'
#' Use this instead of a bare `gsub()`: calling `gsub()` on a DNAStringSet
#' coerces it to a plain character vector and loses the Biostrings class. This
#' wrapper does the substitution on the character representation and coerces the
#' result back to a DNAStringSet, preserving sequence names.
#'
#' @param seqs DNAStringSet; sequences to clean.
#' @param replacement length-1 character used to replace ambiguous characters
#'   (default "N").
#' @param keep character; the set of characters to preserve, given as a single
#'   string (default "ACGT"). Matching is case-insensitive. To preserve
#'   alignment gaps, pass e.g. `keep = "ACGT-"`.
#'
#' @return DNAStringSet with ambiguous characters replaced; names preserved.
#' @export
mask_ambiguous <- function(seqs, replacement = "N", keep = "ACGT") {
  if (!methods::is(seqs, "DNAStringSet")) {
    stop("seqs must be a DNAStringSet")
  }
  if (length(replacement) != 1L || nchar(replacement) != 1L) {
    stop("replacement must be a single character")
  }

  # Build a negated character class from the keep-set (both cases), so anything
  # outside it is replaced. e.g. keep = "ACGT" -> pattern "[^ACGTacgt]"
  pattern <- paste0("[^", keep, tolower(keep), "]")

  # gsub() returns a plain character vector (names preserved); coerce back so
  # the caller keeps a DNAStringSet for downstream Biostrings operations.
  Biostrings::DNAStringSet(
    gsub(pattern, replacement, as.character(seqs))
  )
}
