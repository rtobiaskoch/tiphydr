#' Nucleotide composition at alignment columns
#'
#' Counts A/C/G/T, gaps and ambiguity codes at the requested alignment
#' columns, and gives the percent of each unambiguous base among the A/C/G/T
#' calls only (gaps and ambiguous codes are reported as counts, not in the
#' percent denominator). Positions are alignment columns, so they only mean
#' the same site across sequences if the input is aligned.
#'
#' @param alignment A [Biostrings::DNAMultipleAlignment], or an aligned
#'   [Biostrings::DNAStringSet] (all sequences the same width).
#' @param positions Integer vector of 1-based alignment columns. `NULL`
#'   (default) uses every column.
#'
#' @return A tibble with one row per position: `position`, counts `A`, `C`,
#'   `G`, `T`, `gap`, `ambiguous`, `n_acgt`, and `pct_A`, `pct_C`, `pct_G`,
#'   `pct_T`.
#'
#' @examples
#' aln <- Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "ACNT", s3 = "A-GA"))
#' nt_freq(aln, positions = c(2, 4))
#'
#' @importFrom Biostrings consensusMatrix width
#' @export
nt_freq <- function(alignment, positions = NULL) {

  # DNAMultipleAlignment -> DNAStringSet so one code path handles both inputs
  if (methods::is(alignment, "DNAMultipleAlignment")) {
    alignment <- methods::as(alignment, "DNAStringSet")
  }
  if (!methods::is(alignment, "DNAStringSet")) {
    stop("`alignment` must be a DNAMultipleAlignment or DNAStringSet, got: ",
         class(alignment)[1])
  }

  # Guard: columns only mean the same site across sequences if aligned
  aln_width <- unique(Biostrings::width(alignment))
  if (length(aln_width) != 1) {
    stop("`alignment` sequences must all be the same width (aligned).")
  }

  if (is.null(positions)) positions <- seq_len(aln_width)
  if (any(positions < 1 | positions > aln_width)) {
    stop("`positions` must be within 1 and ", aln_width, ".")
  }

  # rows = IUPAC letters + gap, cols = alignment columns, values = counts
  counts <- Biostrings::consensusMatrix(alignment)[, positions, drop = FALSE]

  # t() -> one row per position; as_tibble drops the stray row names
  tibble::as_tibble(t(counts)) |>
    dplyr::transmute(
      position  = positions,
      A, C, G, T,
      gap       = `-`,
      ambiguous = length(alignment) - A - C - G - T - gap,
      n_acgt    = A + C + G + T
    ) |>
    dplyr::mutate(dplyr::across(
      c(A, C, G, T),
      ~ round(100 * .x / n_acgt, 2),
      .names = "pct_{.col}"
    ))
}
