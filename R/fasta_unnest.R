#' Extract metadata and sequences from a nested sequence dataframe
#'
#' Inverse of fasta_nest(). Rows with NULL sequences (no FASTA match) are
#' dropped with a message. Returns metadata and sequences as separate elements.
#'
#' @param nested_df tibble with a `sequence` list-column as produced by
#'   fasta_nest()
#' @param verbose logical; if TRUE prints count of rows dropped (default TRUE)
#'
#' @return named list: metadata (tibble, sequence column removed) and
#'   biostring (DNAStringSet)
#' @export
fasta_unnest <- function(nested_df, verbose = TRUE) {

  if (!"sequence" %in% names(nested_df)) {
    stop("nested_df must have a 'sequence' list-column (produced by fasta_nest())")
  }

  # purrr::map_lgl checks each element of the list-column for NULL
  null_mask <- purrr::map_lgl(nested_df$sequence, is.null)
  n_null    <- sum(null_mask)

  if (n_null > 0L) {
    if (verbose) message("Dropping ", n_null, " rows with NULL sequences")
    nested_df <- nested_df[!null_mask, ]
  }

  # do.call(c, ...) concatenates the list of length-1 DNAStringSets into one
  biostring <- do.call(c, nested_df$sequence)

  # dplyr::select removes the sequence list-column, leaving clean metadata
  metadata  <- dplyr::select(nested_df, -sequence)

  list(metadata = metadata, biostring = biostring)
}
