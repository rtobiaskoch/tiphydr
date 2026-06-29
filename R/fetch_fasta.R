#' Fetch one or more sequences from NCBI as a FASTA file
#'
#' Uses `rentrez::entrez_fetch()` to retrieve sequences from NCBI and write
#' them to a single FASTA file. Accepts one or more accession numbers; all
#' results are concatenated into `filename`.
#'
#' @param accession character vector of NCBI accession numbers (e.g. "AF481864")
#' @param filename  path to the output FASTA file (created or overwritten)
#' @param db        NCBI database to query (default `"nucleotide"`)
#' @param verbose   logical; print progress messages (default `TRUE`)
#'
#' @return A `DNAStringSet` of the fetched sequences, invisibly. The primary
#'   side effect is writing `filename`.
#' @export
fetch_fasta <- function(accession, filename, db = "nucleotide", verbose = TRUE) {

  # --- input validation -------------------------------------------------------
  if (!is.character(accession) || length(accession) == 0) {
    stop("accession must be a non-empty character vector")
  }
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) {
    stop("filename must be a single non-empty string")
  }

  if (verbose) {
    message(
      "Fetching ", length(accession), " accession(s) from NCBI [", db, "]: ",
      paste(accession, collapse = ", ")
    )
  }

  # entrez_fetch accepts a character vector of accession strings as `id`,
  # fetching all in a single API call and returning concatenated FASTA text
  raw_fasta <- rentrez::entrez_fetch(
    db      = db,
    id      = accession,
    rettype = "fasta",
    retmode = "text"
  )

  if (!nzchar(trimws(raw_fasta))) {
    stop("NCBI returned no sequences for: ", paste(accession, collapse = ", "))
  }

  # Write the raw FASTA text to disk before parsing
  writeLines(raw_fasta, filename)

  seqs <- Biostrings::readDNAStringSet(filename)

  if (verbose) {
    message("Written ", length(seqs), " sequence(s) to ", filename)
  }

  invisible(seqs)
}
