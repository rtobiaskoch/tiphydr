#' Read one or more FASTA files into DNAStringSet objects
#'
#' Wraps `Biostrings::readDNAStringSet()` with multi-file support. When a
#' single path is supplied the DNAStringSet is returned directly; when multiple
#' paths are supplied a named list is returned, with list names derived from
#' the file basenames (extension stripped).
#'
#' @param paths   character vector of file paths to FASTA files
#' @param verbose logical; if TRUE prints sequence count per file (default TRUE)
#'
#' @return Single path: a `DNAStringSet`. Multiple paths: a named list of
#'   `DNAStringSet` objects, list names derived from filenames without extension.
#' @export
fasta_read <- function(paths, verbose = TRUE) {

  # --- input validation ------------------------------------------------------
  if (!is.character(paths) || length(paths) == 0) {
    stop("paths must be a non-empty character vector")
  }

  # Guard against empty strings (e.g. system.file() returning "" when a file
  # is not found inside a package) — these give a misleading "Files not found:"
  # error with no filename shown
  blank_paths <- paths[!nzchar(paths)]
  if (length(blank_paths) > 0) {
    stop("paths contains ", length(blank_paths),
         " empty string(s). Did system.file() fail to find the file?")
  }

  # Check all files exist before attempting any reads
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }

  # --- inner reader ----------------------------------------------------------
  # Pure function: reads one FASTA, optionally messages, returns DNAStringSet
  read_one <- function(path) {
    seqs <- Biostrings::readDNAStringSet(path)  # from Biostrings package
    if (verbose) {
      message("Read ", length(seqs), " sequences from ", basename(path))
    }
    seqs
  }

  # --- dispatch on number of paths -------------------------------------------
  # Single file: return DNAStringSet directly (not wrapped in a list)
  if (length(paths) == 1L) {
    return(read_one(paths))
  }

  # Multiple files: read all with purrr::map, name by filename sans extension
  result        <- purrr::map(paths, read_one)          # purrr for functional style
  names(result) <- tools::file_path_sans_ext(           # tools:: strips ".fasta" etc.
    basename(paths)
  )

  # Warn when the same filename appears more than once — silently duplicated
  # list names are hard to debug downstream (consistent with fasta_combine pattern)
  dupe_names <- names(result)[duplicated(names(result))]
  if (length(dupe_names) > 0) {
    warning("Duplicate list names produced from paths: ",
            paste(unique(dupe_names), collapse = ", "))
  }

  result
}
