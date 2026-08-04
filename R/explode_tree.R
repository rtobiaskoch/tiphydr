#' Split a dated tree into introduction clades
#'
#' "Explodes" a time-scaled tree into its constituent introductions: edges where
#' the inferred ancestral deme changes. It reads a dated tree and a node x deme
#' probability table (e.g. TreeTime DTA output), detects introductions into every
#' deme using the continued-transmission definition (see
#' \code{\link{detect_introductions}}), and returns both a tidy summary and the
#' subtree for each multi-tip introduction clade.
#'
#' All node demes come from \code{node_probs}; every node's calendar date is read
#' directly off the dating tool's own native output (see
#' \code{\link{node_dates_from_timetree}} for why \code{annotated_tree_path} must
#' be that native, pre-round-trip file and must share \code{tree}'s exact
#' topology).
#'
#' @param tree An \code{ape::phylo} time tree.
#' @param node_probs A data frame with a \code{node} column (matching ape node
#'   numbering) and one column per deme holding state probabilities (e.g. TreeTime
#'   DTA output loaded with \code{read.delim()}). Optional; supply either this
#'   with \code{annotated_tree_path}, or both \code{clade_dwell} and
#'   \code{tip_membership}.
#' @param annotated_tree_path Path to the dating tool's native annotated NEXUS
#'   output, passed through to \code{\link{node_dates_from_timetree}}. Required
#'   when \code{node_probs} is supplied.
#' @param clade_dwell A data frame with clade dwell times (simmap source). Optional;
#'   must be supplied with \code{tip_membership} if using the simmap path.
#' @param tip_membership A data frame with tip membership information (simmap source).
#'   Optional; must be supplied with \code{clade_dwell} if using the simmap path.
#' @param confidence Minimum node probability for an established introduction
#'   (default 0.5). Passed to \code{\link{detect_introductions}}. When internal
#'   transitions exist but none clear it, a warning reports the max available
#'   probability and only singleton introductions are returned (empty
#'   \code{trees}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{introductions}}{A tibble, one row per kept tip (see
#'       \code{\link{detect_introductions}}), including \code{intro_node} (the
#'       ape node id of the introduction) so callers can join back to the
#'       source clade table.}
#'     \item{\code{trees}}{A named list of \code{ape::phylo} subtrees, one per
#'       multi-tip introduction clade, named by \code{intro_clade_id}. Singleton
#'       introductions (\code{clade_size == 1}) appear in \code{introductions}
#'       only.}
#'   }
#' @export
explode_tree <- function(
  tree,
  node_probs = NULL,
  annotated_tree_path = NULL,
  clade_dwell = NULL,
  tip_membership = NULL,
  confidence = 0.5
) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be an ape::phylo object")
  }

  has_node <- !is.null(node_probs)
  has_simmap <- !is.null(clade_dwell) || !is.null(tip_membership)

  if (has_node && has_simmap) {
    stop(
      "Supply exactly one source: node_probs (+ annotated_tree_path), OR ",
      "both clade_dwell and tip_membership -- not both at once."
    )
  }
  if (!has_node && !has_simmap) {
    stop(
      "Supply exactly one source: node_probs (+ annotated_tree_path), OR ",
      "both clade_dwell and tip_membership."
    )
  }
  if (has_simmap && (is.null(clade_dwell) || is.null(tip_membership))) {
    stop("clade_dwell and tip_membership must both be supplied together.")
  }

  if (has_node) {
    if (is.null(annotated_tree_path)) {
      stop("annotated_tree_path is required when using node_probs.")
    }
    if (!is.data.frame(node_probs) || !"node" %in% names(node_probs)) {
      stop("node_probs must be a data frame with a 'node' column")
    }
    tree_df <- build_tree_df(tree, node_probs, annotated_tree_path)
    intros <- detect_introductions(tree_df, confidence = confidence)
  } else {
    stop("simmap source not yet implemented") # replaced in Task 4-6
  }

  # Build one subtree per multi-tip clade: descend to the introduction node, then
  # keep only the tips that passed the continued-transmission filter.
  multi <- intros |>
    dplyr::filter(.data$clade_size > 1) |>
    dplyr::arrange(.data$intro_clade_id)

  if (nrow(multi) == 0L) {
    trees <- list()
  } else {
    trees <- multi |>
      dplyr::group_by(.data$intro_clade_id) |>
      dplyr::group_map(function(rows, key) {
        clade <- ape::extract.clade(tree, rows$intro_node[[1]])
        ape::keep.tip(clade, rows$tipname)
      }) |>
      stats::setNames(unique(multi$intro_clade_id))
  }

  list(
    introductions = intros,
    trees = trees
  )
}
