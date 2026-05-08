#' Rename tip labels in a phylogenetic tree.
#'
#' Matches each value in `id` against `tree$tip.label` and replaces the
#' matched tip with the parallel value from `newnames`. Tips with no match
#' are left unchanged.
#'
#' @param tree     ape::phylo object.
#' @param id       Character vector of values matching `tree$tip.label`.
#' @param newnames Character vector of replacement labels (same length as `id`).
#' @param match    Logical (default TRUE). When TRUE, errors if any value in
#'                 `id` is absent from `tree$tip.label`.
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

  # when match = TRUE every id must resolve to an existing tip
  if (match) {
    missing_from_tree <- id[!id %in% tree$tip.label]
    if (length(missing_from_tree) > 0) {
      stop(
        "id values not found in tree tips: ",
        paste(missing_from_tree, collapse = ", ")
      )
    }
  }

  # named vector: old_label -> new_label
  lookup <- stats::setNames(newnames, id)

  # replace matched tips; keep originals where no match exists
  tree$tip.label <- ifelse(
    tree$tip.label %in% names(lookup),
    lookup[tree$tip.label],
    tree$tip.label
  )

  tree
}
