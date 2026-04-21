#' Combine a list of DNAStringSet objects into one
#'
#' @param fasta_list list of DNAStringSet objects (e.g. output of fasta_read()
#'   with multiple paths)
#' @param verbose logical; if TRUE prints combined sequence count (default TRUE)
#'
#' @return DNAStringSet
#' @export
fasta_combine <- function(fasta_list, verbose = TRUE) {

  # Fix 1: separate type and length checks for clearer error messages
  if (!is.list(fasta_list)) {
    stop("fasta_list must be a list of DNAStringSet objects, got: ", class(fasta_list))
  }
  if (length(fasta_list) == 0) {
    stop("fasta_list must be a non-empty list")
  }

  # Fix 2: validate that every element is a DNAStringSet
  not_dna <- !purrr::map_lgl(fasta_list, is, "DNAStringSet")
  if (any(not_dna)) {
    stop("All elements of fasta_list must be DNAStringSet objects. ",
         "Non-DNAStringSet elements at positions: ",
         paste(which(not_dna), collapse = ", "))
  }

  combined <- do.call(c, fasta_list)

  dupe_names <- names(combined)[duplicated(names(combined))]
  if (length(dupe_names) > 0) {
    warning("Duplicate sequence names found: ", paste(unique(dupe_names), collapse = ", "))
  }

  if (verbose) {
    message("Combined ", length(fasta_list), " DNAStringSets into ",
            length(combined), " total sequences")
  }

  combined
}
