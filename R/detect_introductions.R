#' Detect deme introductions from an annotated tree table
#'
#' An introduction is an edge where a node's inferred deme differs from its
#' parent's deme — the lineage entered a new deme. This detects introductions into
#' every deme, in two forms, mirroring the established hand-parsing workflow:
#'
#' \itemize{
#'   \item \strong{Established introductions}: an internal node (confidence
#'     \eqn{\ge} \code{confidence}) whose deme differs from its parent's. Its
#'     descendant tips are kept only under the \emph{continued-transmission}
#'     rule — every ancestor on the path from the tip back to the introduction
#'     node must also be the destination deme, i.e. the lineage never left.
#'   \item \strong{Singleton introductions}: a tip whose deme differs from its
#'     parent's deme (an introduction with no sampled onward transmission).
#' }
#'
#' Each independent transition into a deme is its own clade, so reintroductions
#' are counted separately. Clades are numbered chronologically by introduction
#' date (\code{intro_clade_id}).
#'
#' @param tree_df A \code{tbl_tree} from \code{\link{build_tree_df}}.
#' @param confidence Minimum node probability for an established (internal-node)
#'   introduction. Default 0.5. Singleton (tip) introductions are not filtered.
#'   If internal transitions exist but none clear \code{confidence}, a warning
#'   reports the maximum available probability so the threshold can be tuned.
#'
#' @return A tibble with one row per kept tip: \code{tipname}, \code{deme}
#'   (destination), \code{intro_clade_id}, \code{inferred_intro_date},
#'   \code{last_sample_date} (most recent tip date in the clade),
#'   \code{inferred_intro_source} (parent deme), \code{inferred_intro_source_probability},
#'   \code{clade_size}, and \code{intro_node} (ape node id of the introduction,
#'   used to extract subtrees).
#' @export
detect_introductions <- function(tree_df, confidence = 0.51) {
  # Per-node lookups keyed by node id (character keys for safe subsetting).
  state <- stats::setNames(tree_df$inferred_state, tree_df$node)
  conf <- stats::setNames(tree_df$confidence_state, tree_df$node)
  date <- stats::setNames(tree_df$inferred_date, tree_df$node)

  # Attach each node's parent attributes. Root's parent is itself -> no transition.
  # Drop the tbl_tree class first (plain tibble) so dplyr does not try to
  # reconstruct an invalid <phylo> from the mutated table.
  df <- tibble::as_tibble(tree_df) |>
    dplyr::mutate(
      # unname: subsetting a named lookup vector carries the key as a name
      parent_state = unname(state[as.character(.data$parent)]),
      parent_conf = unname(conf[as.character(.data$parent)])
    )

  # Edges where the deme changed = introductions.
  transitions <- dplyr::filter(df, .data$inferred_state != .data$parent_state)

  # ── Established introductions (internal nodes above confidence) ───────────────
  internal_transitions <- dplyr::filter(transitions, !.data$is_tip)
  established <- dplyr::filter(
    internal_transitions,
    .data$confidence_state >= confidence
  )

  # Explain a (possibly empty) trees list: internal transitions exist but the
  # confidence filter removed all of them, so only singletons can be returned.
  if (nrow(internal_transitions) > 0L && nrow(established) == 0L) {
    warning(
      sprintf(
        paste0(
          "%d internal-node transition(s) found but none reach confidence ",
          "%.2f (max available %.2f); no established (multi-tip) clades — ",
          "only singleton introductions returned. Lower `confidence` to ",
          "include them."
        ),
        nrow(internal_transitions),
        confidence,
        max(internal_transitions$confidence_state)
      ),
      call. = FALSE
    )
  }

  established_tips <- purrr::pmap_dfr(
    list(
      established$node,
      established$inferred_state,
      established$inferred_date,
      established$parent_state,
      established$parent_conf
    ),
    function(intro_node, dest, intro_date, source_state, source_prob) {
      # Descendant tips of the introduction node.
      tips <- tidytree::offspring(tree_df, intro_node, tiponly = TRUE)
      if (nrow(tips) == 0L) {
        return(NULL)
      }

      # Continued-transmission filter: keep a tip only if every ancestor at or
      # after the introduction date is the destination deme (never left).
      keep <- purrr::map_lgl(tips$node, function(tip_node) {
        # as_tibble drops tbl_tree so the filter does not trigger a phylo rebuild
        anc <- tibble::as_tibble(tidytree::ancestor(tree_df, tip_node))
        anc <- dplyr::filter(anc, .data$inferred_date >= intro_date)
        all(anc$inferred_state == dest)
      })

      tips <- tips[keep, , drop = FALSE]
      if (nrow(tips) == 0L) {
        return(NULL)
      }

      tibble::tibble(
        tipname = tips$label,
        deme = dest,
        inferred_intro_date = intro_date,
        inferred_intro_source = source_state,
        inferred_intro_source_probability = source_prob,
        sample_date = tips$inferred_date,
        intro_node = intro_node
      )
    }
  )

  # ── Singleton introductions (tips whose parent is a different deme) ────────────
  singletons <- transitions |>
    dplyr::filter(.data$is_tip) |>
    dplyr::transmute(
      tipname = .data$label,
      deme = .data$inferred_state,
      inferred_intro_date = .data$inferred_date,
      inferred_intro_source = .data$parent_state,
      inferred_intro_source_probability = .data$parent_conf,
      sample_date = .data$inferred_date,
      intro_node = .data$node
    )

  intros <- dplyr::bind_rows(established_tips, singletons)
  if (nrow(intros) == 0L) {
    return(intros)
  }

  # One clade id per introduction node, numbered chronologically by intro date
  # (ties broken by node id so each distinct clade gets a distinct id).
  clade_order <- intros |>
    dplyr::distinct(.data$intro_node, .data$inferred_intro_date) |>
    dplyr::arrange(.data$inferred_intro_date, .data$intro_node) |>
    dplyr::mutate(intro_clade_id = dplyr::row_number()) |>
    dplyr::select("intro_node", "intro_clade_id")

  intros |>
    dplyr::left_join(clade_order, by = "intro_node") |>
    dplyr::add_count(.data$intro_node, name = "clade_size") |>
    # last_sample_date: most recent tip date within each clade
    dplyr::mutate(
      last_sample_date = max(.data$sample_date),
      .by = "intro_node"
    ) |>
    dplyr::arrange(.data$intro_clade_id) |>
    dplyr::select(
      "tipname",
      "deme",
      "intro_clade_id",
      "inferred_intro_date",
      "last_sample_date",
      "inferred_intro_source",
      "inferred_intro_source_probability",
      "clade_size",
      "intro_node"
    )
}
