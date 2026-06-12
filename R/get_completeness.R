#' Genome completeness from unambiguous base content
#'
#' Computes, per sequence, the fraction of positions that are unambiguous
#' nucleotides (A/C/G/T) out of the sequence's total length. This is a
#' reference-free proxy for how complete a consensus genome is: gaps (`-`, `.`)
#' and ambiguity codes (`N`, `R`, `Y`, ...) count as missing. It lets you derive
#' a coverage-like quality metric directly from the genome itself, which is
#' useful when per-sample read coverage (e.g. `samtools coverage`) is
#' unavailable - for instance for genomes already deposited in a database.
#'
#' Zero-length sequences return a completeness of 0 (rather than `NaN`).
#'
#' @param sequences A [Biostrings::DNAStringSet].
#'
#' @return A data.frame with one row per input sequence and two columns:
#'   `seq` (the sequence name) and `completeness` (numeric fraction in `[0, 1]`).
#'
#' @examples
#' \dontrun{
#' ss <- Biostrings::readDNAStringSet("genomes.fasta")
#' get_completeness(ss)
#' }
#'
#' @importFrom Biostrings alphabetFrequency width
#' @export
get_completeness <- function(sequences) {

  if (!methods::is(sequences, "DNAStringSet")) {
    stop("`sequences` must be a DNAStringSet, got: ", class(sequences)[1])
  }

  # baseOnly = TRUE collapses the alphabet to columns A, C, G, T, other - so
  # every gap/ambiguity character lands in "other" and is excluded from the
  # numerator below.
  counts <- Biostrings::alphabetFrequency(sequences, baseOnly = TRUE)

  atcg   <- rowSums(counts[, c("A", "C", "G", "T"), drop = FALSE])
  widths <- Biostrings::width(sequences)

  # Guard against zero-length sequences (some consensus pipelines emit empties).
  completeness <- ifelse(widths > 0, atcg / widths, 0)

  data.frame(
    seq          = names(sequences),
    completeness = completeness,
    stringsAsFactors = FALSE
  )
}
