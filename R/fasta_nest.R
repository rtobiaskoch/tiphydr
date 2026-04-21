#' Nest a DNAStringSet as a list-column in a metadata dataframe
#'
#' Analogous to an sf geometry column: one sequence per metadata row.
#' The resulting `sequence` list-column is compatible with dplyr::filter(),
#' mutate(), and group_by(). Use fasta_unnest() to reverse.
#'
#' Matching is done by splitting each sequence name on `delim` and checking
#' whether the row's `id_col` value exactly equals any resulting part.
#' This avoids regex metacharacter issues common in FASTA sequence names.
#'
#' @param df dataframe with metadata
#' @param biostring DNAStringSet of sequences
#' @param id_col string; name of column in df to match against sequence name parts
#' @param delim string; delimiter used to split sequence names before matching
#'   (e.g. "|" splits "WNV|2021|CO|CSU123" into c("WNV","2021","CO","CSU123"))
#' @param verbose logical; if TRUE prints match summary (default TRUE)
#'
#' @return tibble with all original columns plus a `sequence` list-column of
#'   length-1 DNAStringSet objects (NULL where no match found)
#' @export
fasta_nest <- function(df, biostring, id_col, delim, verbose = TRUE) {

  # Validate that id_col exists in the dataframe before proceeding
  if (!id_col %in% names(df)) {
    stop("id_col '", id_col, "' not found in df columns: ",
         paste(names(df), collapse = ", "))
  }

  # names() pulls the sequence identifiers from the DNAStringSet (e.g. "WNV|2021|CO|NY10_001")
  seq_names   <- names(biostring)

  # Pre-split all sequence names once — reused for every row lookup.
  # stringr::fixed() treats delim as a literal string, not a regex pattern.
  # This is critical because "|" is a regex metacharacter that would otherwise
  # match everything.
  split_names <- stringr::str_split(seq_names, stringr::fixed(delim))

  # --- Inner matcher: returns a length-1 DNAStringSet or NULL ---
  match_sequence <- function(id_value) {
    # Check if id_value appears as an exact part in any split sequence name
    matched_idx <- which(
      purrr::map_lgl(split_names, ~ id_value %in% .x)
    )

    if (length(matched_idx) == 0L) return(NULL)

    # Warn but recover gracefully when an id matches more than one sequence
    if (length(matched_idx) > 1L) {
      warning("id '", id_value, "' matched ", length(matched_idx),
              " sequence names; using first: ", seq_names[matched_idx[1]])
    }

    # Subset the DNAStringSet to return exactly one sequence, preserving its name
    biostring[matched_idx[1]]
  }

  # df[[id_col]] pulls the matching column vector from the metadata tibble
  id_values <- df[[id_col]]

  # Apply match_sequence() to every row's id value; result is a list of
  # DNAStringSet (length 1) or NULL entries
  sequences <- purrr::map(id_values, match_sequence)

  # --- Summary counts for messaging and warning ---
  n_matched   <- sum(!purrr::map_lgl(sequences, is.null))
  n_unmatched <- nrow(df) - n_matched

  # Warn once, listing all unmatched ids, so users can diagnose naming issues
  if (n_unmatched > 0L) {
    unmatched_ids <- id_values[purrr::map_lgl(sequences, is.null)]
    warning(n_unmatched, " rows had no matching sequence: ",
            paste(unmatched_ids, collapse = ", "))
  }

  if (verbose) {
    message("Matched ", n_matched, " of ", nrow(df), " rows to sequences")
  }

  # dplyr::mutate adds the sequence list-column without modifying the input df
  dplyr::mutate(df, sequence = sequences)
}
