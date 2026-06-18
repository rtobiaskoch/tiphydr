#' Assemble an annotated tree table for introduction detection
#'
#' Converts a tree to a tidy \code{tbl_tree} (one row per node, with \code{parent}
#' edges) and joins on the inferred deme, its probability, and an inferred calendar
#' date for every node. This is the single object consumed by
#' \code{\link{detect_introductions}}.
#'
#' Node states (tips and internal nodes) come from \code{node_probs}; node dates
#' come from the time tree via \code{\link{node_dates_from_timetree}}, with tip
#' dates taken from \code{tip_dates}.
#'
#' @param tree An \code{ape::phylo} time tree.
#' @param node_probs Node x deme probability table (see
#'   \code{\link{most_likely_deme}}). Node ids must match ape node numbering.
#' @param tip_dates Tibble with \code{tip_label} and \code{date} (\code{Date}).
#'
#' @return A \code{tbl_tree} (tibble) with columns \code{node}, \code{parent},
#'   \code{label}, \code{inferred_state}, \code{confidence_state},
#'   \code{inferred_date}, and \code{is_tip}. The \code{tbl_tree} class is kept so
#'   that \code{tidytree} traversal (\code{parent}, \code{ancestor},
#'   \code{offspring}) works downstream.
#' @export
build_tree_df <- function(tree, node_probs, tip_dates) {

  ntip <- ape::Ntip(tree)

  # Topology table: node, parent, branch.length, label. Keeps class tbl_tree.
  base <- tidytree::as_tibble(tree)

  states <- most_likely_deme(node_probs)
  dates  <- node_dates_from_timetree(tree, tip_dates)

  base |>
    dplyr::left_join(states, by = "node") |>
    dplyr::left_join(dates, by = "node") |>
    dplyr::mutate(is_tip = .data$node <= ntip)
}
