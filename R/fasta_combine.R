#' Combine a list of DNAStringSet objects into one
#'
#' @param fasta_list list of DNAStringSet objects (e.g. output of fasta_read()
#'   with multiple paths)
#' @param verbose logical; if TRUE prints combined sequence count (default TRUE)
#'
#' @return DNAStringSet
#' @export
fasta_combine <- function(fasta_list, verbose = TRUE) {

  if (!is.list(fasta_list) || length(fasta_list) == 0) {
    stop("fasta_list must be a non-empty list of DNAStringSet objects")
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
