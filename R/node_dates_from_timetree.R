#' Infer calendar dates for every node of a dated (time-scaled) tree
#'
#' A time tree's branch lengths are proportional to elapsed time, but the unit
#' (days vs decimal years) is often unclear and the root has no calendar anchor.
#' This fits a molecular-clock linear model — tip calendar date as a function of
#' root-to-node distance — and predicts a date for every internal node. Because
#' both the slope (the time unit) and the intercept (the root date) are estimated
#' from the tip dates, the result is correct regardless of the branch-length unit.
#'
#' Tip dates are returned as their observed values, not the fitted ones.
#'
#' The fit's \eqn{R^2} is a clocklikeness / temporal-signal check. For a genuine
#' time tree the tips fall almost exactly on the line (\eqn{R^2 \approx 1}); a low
#' value flags that the input is probably not time-scaled (e.g. a phylogram in
#' substitutions/site), violates a strict clock, or contains an erroneous tip
#' date. A warning is raised below \code{min_r2} but dates are still returned.
#'
#' @param tree An \code{ape::phylo} object with branch lengths (a time tree).
#' @param tip_dates A tibble with columns \code{tip_label} (matching
#'   \code{tree$tip.label}) and \code{date} (\code{Date}).
#' @param min_r2 Minimum acceptable \eqn{R^2} of the date~depth fit before a
#'   warning is emitted (default 0.95).
#'
#' @return A tibble with columns \code{node} (ape node id) and
#'   \code{inferred_date} (\code{Date}) for all tips and internal nodes.
#' @export
node_dates_from_timetree <- function(tree, tip_dates, min_r2 = 0.95) {

  if (is.null(tree$edge.length)) {
    stop("tree must have branch lengths (a time-scaled / dated tree)")
  }

  # Root-to-node distance for every node, indexed by ape node id.
  depth <- ape::node.depth.edgelength(tree)
  ntip  <- ape::Ntip(tree)

  # Align observed tip dates to tree tip order (ape tip ids are 1:ntip).
  tips <- tibble::tibble(
    node      = seq_len(ntip),
    tip_label = tree$tip.label,
    depth     = depth[seq_len(ntip)]
  ) |>
    dplyr::left_join(tip_dates, by = "tip_label")

  if (anyNA(tips$date)) {
    missing <- tips$tip_label[is.na(tips$date)]
    stop("no date supplied for tips: ", paste(missing, collapse = ", "))
  }

  # No spread in tip depth (e.g. contemporaneous sampling) means the slope is
  # unidentifiable — the regression cannot recover the time scale.
  if (stats::var(tips$depth) == 0) {
    stop("all tips are equidistant from the root (no temporal spread in depth); ",
         "cannot estimate node dates")
  }

  # Molecular-clock fit: numeric date ~ depth. Slope absorbs the time unit.
  fit <- stats::lm(as.numeric(date) ~ depth, data = tips)

  # Clocklikeness / temporal-signal check: tips should lie on the line.
  r2 <- summary(fit)$r.squared
  if (r2 < min_r2) {
    warning(
      sprintf(
        paste0("weak temporal signal: date~depth R2 = %.3f (< %.2f). The tree ",
               "may not be time-scaled or may violate a strict clock; inferred ",
               "internal-node dates could be unreliable."),
        r2, min_r2
      ),
      call. = FALSE
    )
  }

  # Predict a date for every node from its depth (unname: predict() carries names).
  all_nodes <- tibble::tibble(node = seq_along(depth), depth = depth)
  predicted <- unname(stats::predict(fit, newdata = all_nodes))

  all_nodes |>
    dplyr::mutate(
      inferred_date = as.Date(round(predicted), origin = "1970-01-01")
    ) |>
    # Replace fitted tip dates with the observed dates.
    dplyr::rows_update(
      dplyr::select(tips, node, inferred_date = date),
      by = "node"
    ) |>
    dplyr::select(node, inferred_date)
}
