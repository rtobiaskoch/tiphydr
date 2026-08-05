#' Plot a phylogenetic tree coloured by a trait (discrete or continuous)
#'
#' Colours a tidy tree table by any trait column. By default it follows the
#' nextstrain branch-colouring convention: each branch takes its parent
#' (ancestral) node's trait, while tips and internal nodes are coloured by
#' their own value. Set \code{branch_state = "node"} to colour each branch by
#' its own endpoint instead. The
#' trait type is detected automatically — a numeric column gets a continuous
#' \pkg{viridis} gradient, anything else (character/factor) gets a discrete
#' RColorBrewer palette. Returns a \pkg{ggplot2} object that can be extended
#' with additional layers.
#'
#' If a per-node confidence column is present it is double-encoded as point size
#' and alpha (larger, more opaque = higher confidence). When absent, plain
#' points are drawn — so the function works on any \code{tbl_tree}, not just the
#' output of \code{\link{build_tree_df}}.
#'
#' Integer-coded categories (e.g. demes stored as 1/2/3) will render as a
#' gradient; pass them as a factor to get discrete colours.
#'
#' @param tree_df A \code{tbl_tree} with \code{node}, \code{parent}, and the
#'   \code{trait} column (e.g. from \code{\link{build_tree_df}}).
#' @param trait Character scalar — name of the column to colour by.
#' @param confidence Character scalar — name of an optional numeric confidence
#'   column (default \code{"confidence_state"}). Silently ignored if absent.
#' @param palette RColorBrewer palette name for discrete traits. Default
#'   \code{"Dark2"}. Ignored for continuous traits.
#' @param branch_state Which node's trait paints each branch. \code{"parent"}
#'   (default) is the nextstrain convention — colour flows down from the
#'   ancestor. \code{"node"} paints each branch with its own child endpoint's
#'   trait, so a branch always matches the point it leads to.
#'
#'   Neither is universally right, because a single colour cannot represent a
#'   branch on which the state changes partway. They differ in which half they
#'   misrepresent: \code{"parent"} shows the state the lineage arrived in and
#'   hides a transition that happens on the branch; \code{"node"} shows the
#'   state it ended in and hides where it came from. \code{"node"} is the more
#'   defensible choice for terminal branches specifically, since a tip's state
#'   is observed data rather than a reconstruction — painting it with an
#'   inferred ancestral state puts inference on top of an observation.
#'
#' @return A \code{ggplot} object. Additional \pkg{ggplot2} layers can be added
#'   with \code{+}.
#' @export
plot_tree_trait <- function(tree_df, trait,
                            confidence = "confidence_state",
                            palette    = "Dark2",
                            branch_state = c("parent", "node")) {
  branch_state <- match.arg(branch_state)
  # --- Validate: only require what is actually used --------------------------
  required_cols <- c("node", "parent")
  missing_cols  <- setdiff(required_cols, names(tree_df))
  if (length(missing_cols) > 0) {
    stop(
      "tree_df must be a tidy tree table with columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ", paste(missing_cols, collapse = ", "),
      ".\nBuild one with build_tree_df(tree, node_probs, tip_dates)."
    )
  }
  if (!trait %in% names(tree_df)) {
    stop(
      "Column '", trait, "' not found in tree_df. Available columns: ",
      paste(names(tree_df), collapse = ", ")
    )
  }

  # --- Detect trait type: numeric => continuous, else discrete --------------
  is_continuous <- is.numeric(tree_df[[trait]])

  # --- Detect optional confidence encoding ----------------------------------
  # Present + numeric => double-encode as size/alpha; otherwise plain points.
  has_confidence <- confidence %in% names(tree_df) &&
    is.numeric(tree_df[[confidence]])

  # --- Compute parent-trait for ancestral-state branch colouring ------------
  # Each branch is coloured by its parent node's trait (nextstrain convention:
  # colour flows from the ancestor). Join on ape node numbering via `parent`.
  parent_states <- dplyr::select(
    tree_df,
    node,
    parent_trait = !!rlang::sym(trait)
  )
  node_data <- dplyr::left_join(
    tree_df,
    parent_states,
    by = c("parent" = "node")
  )

  # --- Bin confidence into 4 quartile levels (only if used) -----------------
  if (has_confidence) {
    node_data <- dplyr::mutate(
      node_data,
      conf_level = cut(
        .data[[confidence]],
        breaks         = c(0, 0.25, 0.5, 0.75, 1.0),
        labels         = c("0-25%", "25-50%", "50-75%", "75-100%"),
        include.lowest = TRUE
      )
    )
  }

  # --- Point aesthetics: colour by own trait, size/alpha by confidence ------
  point_aes <- if (has_confidence) {
    ggplot2::aes(color = .data[[trait]], size = .data$conf_level,
                 alpha = .data$conf_level)
  } else {
    ggplot2::aes(color = .data[[trait]])
  }

  # --- Convert to phylo for ggtree (annotations reattached via %<+%) ---------
  phylo <- tidytree::as.phylo(tree_df)

  # Which column paints the branches -- see the branch_state @param.
  branch_col <- if (branch_state == "parent") "parent_trait" else trait

  p <- ggtree::ggtree(phylo)
  p <- ggtree::`%<+%`(p, node_data) +
    ggplot2::aes(color = .data[[branch_col]]) +
    ggtree::geom_tippoint(point_aes) +
    ggtree::geom_nodepoint(point_aes)

  # --- Colour scale: viridis for continuous, brewer for discrete ------------
  # Branches (parent_trait) and points (trait) share one colour scale.
  if (is_continuous) {
    p <- p + ggplot2::scale_color_viridis_c(name = trait, na.value = "grey80")
  } else {
    # brewer.pal minimum is 3; slice to the actual number of levels.
    # na.omit: root node has NA parent_trait — not a real level.
    trait_levels <- unique(stats::na.omit(node_data[[trait]]))
    n_colors <- max(3L, length(trait_levels))
    colors <- RColorBrewer::brewer.pal(n_colors, palette)[seq_along(trait_levels)]
    names(colors) <- trait_levels
    p <- p + ggplot2::scale_color_manual(
      values = colors, name = trait, na.translate = FALSE
    )
  }

  # --- Confidence size/alpha scales (only if used) --------------------------
  if (has_confidence) {
    p <- p +
      ggplot2::scale_size_manual(
        values = c("0-25%" = 1, "25-50%" = 2, "50-75%" = 3, "75-100%" = 4),
        name = "Confidence", na.translate = FALSE
      ) +
      ggplot2::scale_alpha_manual(
        values = c("0-25%" = 0.25, "25-50%" = 0.5, "50-75%" = 0.75, "75-100%" = 1.0),
        name = "Confidence", na.translate = FALSE
      )
  }

  p + ggplot2::theme_classic()
}
