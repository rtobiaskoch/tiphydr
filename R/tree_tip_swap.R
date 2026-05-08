#' Rename tip labels in a phylogenetic tree.
#'
#' Searches each value in `id` as a literal substring within `tree$tip.label`
#' and replaces any matched tip with the parallel value from `newnames`. This
#' allows short identifiers (e.g. accession IDs) to match longer tip labels
#' formatted as `VIRUS|YEAR|STATE|ACCESSION`. Tips with no match are left
#' unchanged.
#'
#' Matching uses `grepl(..., fixed = TRUE)` — literal substring, not regex.
#' This avoids misinterpretation of `|` and other metacharacters that commonly
#' appear in tip label formatting.
#'
#' @param tree     ape::phylo object.
#' @param id       Character vector of substrings to search for within
#'                 `tree$tip.label`.
#' @param newnames Character vector of replacement labels (same length as `id`).
#' @param match    Logical (default TRUE). When TRUE, errors if any pattern in
#'                 `id` does not match any tip label.
#'
#' @return ape::phylo with updated tip labels.
#' @export
tree_tip_swap <- function(tree, id, newnames, match = TRUE) {
  if (length(id) != length(newnames)) {
    stop("id and newnames must have the same length.")
  }

  # duplicate ids would produce ambiguous 1-to-many mappings
  dupes <- id[duplicated(id)]
  if (length(dupes) > 0) {
    stop("duplicate id values: ", paste(unique(dupes), collapse = ", "))
  }

  # for each id pattern, find which tip indices contain it as a substring
  # fixed = TRUE: literal match so | in tip labels is not treated as regex OR
  match_idx <- purrr::map(id, ~ which(grepl(.x, tree$tip.label, fixed = TRUE)))

  # a pattern matching multiple tips is ambiguous — refuse to guess
  multi_match <- id[purrr::map_lgl(match_idx, ~ length(.x) > 1)]
  if (length(multi_match) > 0) {
    stop(
      "id patterns match multiple tips (ambiguous): ",
      paste(multi_match, collapse = ", ")
    )
  }

  # patterns with zero hits
  no_match_mask <- purrr::map_lgl(match_idx, ~ length(.x) == 0)
  no_match <- id[no_match_mask]

  if (match && length(no_match) > 0) {
    stop(
      "id patterns not found in tree tips: ",
      paste(no_match, collapse = ", ")
    )
  }

  # build lookup from the resolved (full) tip label to the replacement name
  # only include patterns that matched exactly one tip
  single_mask <- !no_match_mask
  resolved_tips <- tree$tip.label[purrr::map_int(match_idx[single_mask], 1L)]
  lookup <- stats::setNames(newnames[single_mask], resolved_tips)

  # replace matched tips; leave unmatched tips unchanged
  tree$tip.label <- ifelse(
    tree$tip.label %in% names(lookup),
    lookup[tree$tip.label],
    tree$tip.label
  )

  tree
}
