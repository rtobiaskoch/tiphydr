#' Swap FASTA sequence names using a metadata lookup
#'
#' Replaces each name in `fasta` with the corresponding entry in `replace_col`,
#' matched by position against `match_col`. Sequence names with no match in
#' `match_col` are left unchanged (with a warning).
#'
#' @param fasta      Biostrings::XStringSet (e.g. DNAStringSet) whose `names()`
#'                    are matched against `match_col`.
#' @param match_col  Character vector of ids matched against `names(fasta)`.
#'                    A match_col entry matches a fasta name if it occurs
#'                    anywhere within it (e.g. "accession" matches
#'                    "accession|date|deme"), so exact full-name values still
#'                    work as before.
#' @param replace_col Character vector of replacement names, same length and
#'                    order as `match_col`.
#'
#' @return `fasta` with `names()` swapped where a match was found.
#' @export
fasta_name_switch <- function(fasta, match_col, replace_col) {
  if (!methods::is(fasta, "XStringSet")) {
    stop("fasta must be a Biostrings::XStringSet, got: ", class(fasta)[1])
  }

  if (length(match_col) != length(replace_col)) {
    stop("match_col and replace_col must have the same length.")
  }

  # duplicate match_col ids would make the lookup ambiguous (which replace_col wins?)
  dupes <- match_col[duplicated(match_col)]
  if (length(dupes) > 0) {
    stop("duplicate match_col values: ", paste(unique(dupes), collapse = ", "))
  }

  fasta_names <- names(fasta)

  # for each fasta name, find the (single) match_col entry contained within it
  matches <- purrr::map_int(fasta_names, function(nm) {
    hits <- which(stringr::str_detect(nm, stringr::fixed(match_col)))
    if (length(hits) == 0) return(NA_integer_)
    if (length(hits) > 1) {
      stop("name '", nm, "' matches multiple match_col entries: ",
           paste(match_col[hits], collapse = ", "))
    }
    hits
  })

  unmatched <- fasta_names[is.na(matches)]
  if (length(unmatched) > 0) {
    warning(
      "no match_col entry for: ", paste(unmatched, collapse = ", "),
      " — names left unchanged"
    )
  }

  # replace matched names with corresponding replace_col entries; leave unmatched names as-is
  names(fasta)[!is.na(matches)] <- replace_col[matches[!is.na(matches)]]

  fasta
}
