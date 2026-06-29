# Internal: combine base query + named filter vector into a single NCBI query string
build_entrez_query <- function(query, filters = NULL) {
  if (is.null(filters)) return(query)
  filter_terms <- paste0(filters, "[", names(filters), "]")
  paste(c(query, filter_terms), collapse = " AND ")
}

#' Fetch all FASTA sequences matching an NCBI search query
#'
#' Searches NCBI using a query string plus optional field-tagged filters, then
#' downloads all matching sequences in batches using NCBI web history. All
#' results are written to a single FASTA file.
#'
#' NCBI query syntax: field tags are written as `value[Field]`, e.g.
#' `"West Nile virus[Organism]"`. Pass additional constraints via `filters`
#' rather than embedding them in `query` directly.
#'
#' Set the `ENTREZ_KEY` environment variable to your NCBI API key to raise
#' the rate limit from 3 to 10 requests per second.
#'
#' @param query      character; base NCBI search term (e.g. `"West Nile virus[Organism]"`)
#' @param filename   character; path to output FASTA file (created or overwritten)
#' @param filters    named character vector of additional NCBI field filters.
#'                   Names are field tags; values are the search terms.
#'                   e.g. `c(Title = "complete genome", PDAT = "2010:2024")`
#'                   Each entry is combined with `AND` as `value[Name]`.
#' @param db         character; NCBI database (default `"nucleotide"`)
#' @param batch_size integer; sequences per API call (default 500; max 10000)
#' @param retmax     numeric; maximum total sequences to fetch (default `Inf` = all)
#' @param verbose    logical; print progress messages (default `TRUE`)
#'
#' @return A `DNAStringSet` of the fetched sequences, invisibly. Primary side
#'   effect is writing `filename`.
#' @export
fetch_fasta_query <- function(
    query,
    filename,
    filters    = NULL,
    db         = "nucleotide",
    batch_size = 500L,
    retmax     = Inf,
    verbose    = TRUE
) {

  # --- input validation -------------------------------------------------------
  if (!is.character(query) || length(query) != 1L || !nzchar(query)) {
    stop("query must be a single non-empty string")
  }
  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) {
    stop("filename must be a single non-empty string")
  }
  if (!is.null(filters)) {
    if (!is.character(filters) || is.null(names(filters)) || any(!nzchar(names(filters)))) {
      stop("filters must be a named character vector (names are NCBI field tags, e.g. c(Organism = 'WNV'))")
    }
  }
  if (!is.numeric(retmax) || retmax <= 0) {
    stop("retmax must be a positive number or Inf")
  }

  # --- build full query -------------------------------------------------------
  full_query <- build_entrez_query(query, filters)

  if (verbose) message("Query: ", full_query)

  # --- search with web history ------------------------------------------------
  # use_history = TRUE posts results to NCBI server; avoids re-sending ID lists
  # retmax = 0 here: we only need the total count + web_history handle, not IDs
  search_result <- rentrez::entrez_search(
    db          = db,
    term        = full_query,
    use_history = TRUE,
    retmax      = 0L
  )

  total_available <- as.integer(search_result$count)
  if (total_available == 0L) {
    stop("No sequences found for query: ", full_query)
  }

  n_to_fetch <- min(retmax, total_available)
  if (verbose) {
    message("Found ", total_available, " sequence(s); fetching ", n_to_fetch)
    if (is.finite(retmax) && retmax < total_available) {
      message("  (capped at retmax = ", retmax, ")")
    }
  }

  # --- rate limiting ----------------------------------------------------------
  # NCBI: 3 req/sec without API key, 10 req/sec with ENTREZ_KEY env var
  # rentrez reads ENTREZ_KEY automatically for auth; sleep handles client-side pacing
  sleep_secs <- if (nzchar(Sys.getenv("ENTREZ_KEY"))) 0.11 else 0.34

  # --- batch fetch ------------------------------------------------------------
  # Overwrite any existing file before first append
  file.create(filename)

  starts <- seq(0L, n_to_fetch - 1L, by = batch_size)

  for (retstart in starts) {
    this_batch <- min(batch_size, n_to_fetch - retstart)

    if (verbose) {
      message(
        "  Fetching ", retstart + 1L, "–", retstart + this_batch,
        " of ", n_to_fetch, " ..."
      )
    }

    # retstart passed via ... to the underlying NCBI eutils call
    raw_fasta <- rentrez::entrez_fetch(
      db          = db,
      web_history = search_result$web_history,
      rettype     = "fasta",
      retmode     = "text",
      retmax      = this_batch,
      retstart    = retstart
    )

    cat(raw_fasta, file = filename, append = TRUE)

    # Sleep between batches only (not after the last one)
    if (retstart + this_batch < n_to_fetch) Sys.sleep(sleep_secs)
  }

  seqs <- Biostrings::readDNAStringSet(filename)

  if (verbose) {
    message("Written ", length(seqs), " sequence(s) to ", filename)
  }

  invisible(seqs)
}
