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
    intros <- .introductions_from_simmap(clade_dwell, tip_membership, confidence)
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

#' Build a per-tip introductions tibble from simmap-aggregated clade tables
#'
#' Two-step confidence filter, mirroring detect_introductions()'s
#' `confidence_state >= confidence` gate but split across the two axes the
#' simmap data actually has: clade-level (posterior_support -- is this
#' introduction event real at all) and tip-level (membership_prob -- is this
#' specific tip's assignment to it reliable). See
#' docs/superpowers/specs/2026-08-04-explode-tree-simmap-source-design.md for
#' the full rationale.
#'
#' @param clade_dwell Tibble, one row per clade: intro_node, deme,
#'   posterior_support, inferred_intro_source,
#'   inferred_intro_source_probability, clade_size, inferred_intro_date,
#'   last_sample_date (the shape of dta_compute.R's _dta_clade_dwell.tsv).
#' @param tip_membership Tibble, one row per (tipname, intro_node): tipname,
#'   intro_node, deme, membership_prob (the shape of
#'   _dta_tip_membership.tsv) -- the FULL table, not pre-filtered to
#'   is_modal == TRUE rows, since the modal pick here is recomputed after
#'   establishment filtering, not read off the file's own is_modal column.
#' @param confidence Minimum value, shared by both filter steps: clade-level
#'   posterior_support and (post-establishment) tip-level membership_prob.
#' @return A tibble with the same columns as detect_introductions()'s output:
#'   tipname, deme, intro_clade_id, inferred_intro_date, last_sample_date,
#'   inferred_intro_source, inferred_intro_source_probability,
#'   intro_confidence_state, clade_size, intro_node.
.introductions_from_simmap <- function(clade_dwell, tip_membership, confidence) {
  established <- dplyr::filter(clade_dwell, .data$posterior_support >= confidence)

  if (nrow(clade_dwell) > 0L && nrow(established) == 0L) {
    warning(
      sprintf(
        paste0(
          "%d candidate clade(s) found in clade_dwell but none reach ",
          "confidence %.2f (max posterior_support available %.2f); ",
          "explode_tree() will return zero introductions. Lower ",
          "`confidence` to include them."
        ),
        nrow(clade_dwell),
        confidence,
        max(clade_dwell$posterior_support)
      ),
      call. = FALSE
    )
  }

  if (nrow(established) == 0L) {
    return(tibble::tibble(
      tipname = character(0),
      deme = character(0),
      intro_clade_id = integer(0),
      inferred_intro_date = as.Date(character(0)),
      last_sample_date = as.Date(character(0)),
      inferred_intro_source = character(0),
      inferred_intro_source_probability = numeric(0),
      intro_confidence_state = numeric(0),
      clade_size = integer(0),
      intro_node = integer(0)
    ))
  }

  # Restrict tip candidates to established clades BEFORE ranking -- a tip's
  # global argmax (raw is_modal) can point at a clade that didn't establish.
  candidates <- tip_membership |>
    dplyr::inner_join(
      dplyr::select(established, "intro_node", "deme"),
      by = c("intro_node", "deme")
    )

  # Recompute each tip's modal pick among survivors only.
  modal <- candidates |>
    dplyr::group_by(.data$tipname) |>
    dplyr::slice_max(.data$membership_prob, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$membership_prob >= confidence)

  modal |>
    dplyr::left_join(established, by = c("intro_node", "deme")) |>
    dplyr::rename(intro_confidence_state = "posterior_support") |>
    dplyr::select(
      "tipname",
      "deme",
      "intro_clade_id",
      "inferred_intro_date",
      "last_sample_date",
      "inferred_intro_source",
      "inferred_intro_source_probability",
      "intro_confidence_state",
      "clade_size",
      "intro_node"
    )
}
