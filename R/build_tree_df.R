#' Assemble an annotated tree table for introduction detection
#'
#' Converts a tree to a tidy \code{tbl_tree} (one row per node, with \code{parent}
#' edges) and joins on the inferred deme, its probability, and a calendar date
#' for every node. This is the single object consumed by
#' \code{\link{detect_introductions}}.
#'
#' Node states (tips and internal nodes) come from \code{node_probs}; node dates
#' are read directly off the dating tool's own native output via
#' \code{\link{node_dates_from_timetree}} -- see that function's docstring for
#' why \code{annotated_tree_path} must be the dating tool's native, pre-round-trip
#' file and must share \code{tree}'s exact topology.
#'
#' @param tree An \code{ape::phylo} time tree.
#' @param node_probs Node x deme probability table (see
#'   \code{\link{most_likely_deme}}). Node ids must match ape node numbering.
#' @param annotated_tree_path Path to the dating tool's native annotated NEXUS
#'   output, passed through to \code{\link{node_dates_from_timetree}}.
#'
#' @return A \code{tbl_tree} (tibble) with columns \code{node}, \code{parent},
#'   \code{label}, \code{inferred_state}, \code{confidence_state},
#'   \code{inferred_date}, and \code{is_tip}. The \code{tbl_tree} class is kept so
#'   that \code{tidytree} traversal (\code{parent}, \code{ancestor},
#'   \code{offspring}) works downstream.
#' @export
build_tree_df <- function(tree, node_probs, annotated_tree_path) {

  ntip <- ape::Ntip(tree)

  # Topology table: node, parent, branch.length, label. Keeps class tbl_tree.
  base <- tidytree::as_tibble(tree)

  states <- most_likely_deme(node_probs)
  dates  <- node_dates_from_timetree(tree, annotated_tree_path)

  base |>
    dplyr::left_join(states, by = "node") |>
    dplyr::left_join(dates, by = "node") |>
    dplyr::mutate(is_tip = .data$node <= ntip)
}
