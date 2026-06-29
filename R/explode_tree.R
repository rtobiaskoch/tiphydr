#' Split a dated tree into introduction clades
#'
#' "Explodes" a time-scaled tree into its constituent introductions: edges where
#' the inferred ancestral deme changes. It reads a dated tree and a node x deme
#' probability table (e.g. TreeTime DTA output), detects introductions into every
#' deme using the continued-transmission definition (see
#' \code{\link{detect_introductions}}), and returns both a tidy summary and the
#' subtree for each multi-tip introduction clade.
#'
#' Tip calendar dates are parsed from the tree's tip labels (a delimited field);
#' all node demes come from \code{node_probs_path}; internal-node dates are
#' inferred from the time tree (see \code{\link{node_dates_from_timetree}}).
#'
#' @param tree An \code{ape::phylo} time tree whose tip labels embed a date
#'   field (e.g. \code{strain|date|deme|...}).
#' @param node_probs A data frame with a \code{node} column (matching ape node
#'   numbering) and one column per deme holding state probabilities (e.g. TreeTime
#'   DTA output loaded with \code{read.delim()}).
#' @param confidence Minimum node probability for an established introduction
#'   (default 0.5). Passed to \code{\link{detect_introductions}}. When internal
#'   transitions exist but none clear it, a warning reports the max available
#'   probability and only singleton introductions are returned (empty
#'   \code{trees}).
#' @param delim Field delimiter in the tip labels (default \code{"|"}).
#' @param date_field Index of the date field within each tip label (default 2).
#' @param date_format Date format of the date field, passed to
#'   \code{\link[base]{as.Date}} (default \code{"\%Y-\%m-\%d"}). An unparseable
#'   date yields \code{NA} and a clear error rather than aborting mid-parse.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{introductions}}{A tibble, one row per kept tip (see
#'       \code{\link{detect_introductions}}), without the internal \code{intro_node}.}
#'     \item{\code{trees}}{A named list of \code{ape::phylo} subtrees, one per
#'       multi-tip introduction clade, named by \code{intro_clade_id}. Singleton
#'       introductions (\code{clade_size == 1}) appear in \code{introductions}
#'       only.}
#'   }
#' @export
explode_tree <- function(
  tree,
  node_probs,
  confidence = 0.5,
  delim = "|",
  date_field = 2L,
  date_format = "%Y-%m-%d"
) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be an ape::phylo object")
  }
  if (!is.data.frame(node_probs) || !"node" %in% names(node_probs)) {
    stop("node_probs must be a data frame with a 'node' column")
  }

  # Parse the date field out of each tip label. Split wide enough to isolate the
  # date column (the remainder collects in the final column).
  fields <- stringr::str_split_fixed(
    tree$tip.label,
    stringr::fixed(delim),
    n = date_field + 1L
  )
  tip_dates <- tibble::tibble(
    tip_label = tree$tip.label,
    # format = ... makes a bad date NA (caught below) instead of erroring here
    date = as.Date(fields[, date_field], format = date_format)
  )
  if (anyNA(tip_dates$date)) {
    stop(
      "could not parse a valid date from field ",
      date_field,
      " of these tip labels: ",
      paste(tip_dates$tip_label[is.na(tip_dates$date)], collapse = ", ")
    )
  }

  tree_df <- build_tree_df(tree, node_probs, tip_dates)
  intros <- detect_introductions(tree_df, confidence = confidence)

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
