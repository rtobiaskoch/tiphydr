#' Plot a phylogenetic tree coloured by a discrete trait
#'
#' Visualises the output of \code{\link{build_tree_df}} with branches coloured
#' by the parent (ancestral) node's trait — the nextstrain branch-colouring
#' convention — and tips and internal nodes coloured by their own trait value.
#' Returns a \pkg{ggplot2} object that can be extended with additional layers.
#'
#' @param tree_df A \code{tbl_tree} returned by \code{\link{build_tree_df}}.
#' @param trait Character scalar — name of a discrete (character or factor)
#'   column in \code{tree_df} to colour by.
#' @param palette RColorBrewer palette name. Default \code{"Dark2"}.
#'
#' @return A \code{ggplot} object. Additional \pkg{ggplot2} layers can be added
#'   with \code{+}.
#' @export
plot_tree_trait <- function(tree_df, trait, palette = "Dark2") {
  # --- Validate: was tree_df produced by build_tree_df()? --------------------
  # Checked against sentinel columns that build_tree_df() always produces.
  sentinel_cols <- c("node", "parent", "label", "branch.length", "is_tip", "confidence_state")
  missing_cols <- setdiff(sentinel_cols, names(tree_df))
  if (length(missing_cols) > 0) {
    stop(
      "tree_df does not look like output from build_tree_df(). ",
      "Missing columns: ",
      paste(missing_cols, collapse = ", "),
      ".\n",
      "To prepare your tree, run:\n",
      "  build_tree_df(tree, node_probs, tip_dates)"
    )
  }

  # --- Validate: trait column exists and is discrete -------------------------
  if (!trait %in% names(tree_df)) {
    stop(
      "Column '",
      trait,
      "' not found in tree_df. ",
      "Available columns: ",
      paste(names(tree_df), collapse = ", ")
    )
  }
  if (is.numeric(tree_df[[trait]])) {
    stop(
      "Column '",
      trait,
      "' is numeric. plot_tree_trait supports only discrete ",
      "(character or factor) traits. Convert to factor first if needed."
    )
  }

  # --- Compute parent-trait for ancestral-state branch colouring ------------
  # Each branch segment is coloured by its parent node's trait (nextstrain
  # convention: colour flows from the ancestor, revealing where a lineage came from).
  # The join uses ape node numbering: tree_df$parent holds the parent node id.
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

  # --- Bin confidence_state into 4 quartile levels ---------------------------
  # Used to double-encode confidence via point size and alpha — larger and more
  # opaque = higher probability ancestral state assignment.
  conf_breaks <- c(0, 0.25, 0.5, 0.75, 1.0)
  conf_labels <- c("0-25%", "25-50%", "50-75%", "75-100%")
  node_data <- dplyr::mutate(
    node_data,
    conf_level = cut(
      confidence_state,
      breaks        = conf_breaks,
      labels        = conf_labels,
      include.lowest = TRUE
    )
  )

  # --- Build discrete colour palette -----------------------------------------
  # brewer.pal minimum is 3; slice down to the actual number of trait levels.
  # na.omit: root node has NA for parent_trait — not a trait level.
  trait_levels <- unique(stats::na.omit(node_data[[trait]]))
  n_colors <- max(3L, length(trait_levels))
  colors <- RColorBrewer::brewer.pal(n_colors, palette)[seq_len(length(
    trait_levels
  ))]
  names(colors) <- trait_levels

  # --- Convert tbl_tree to phylo for ggtree ----------------------------------
  # as.phylo strips annotation columns; node_data is reattached via %<+% below.
  phylo <- tidytree::as.phylo(tree_df)

  # --- Build plot ------------------------------------------------------------
  # Branches:       parent node's trait (ancestral-state convention)
  # Tip points:     tip's own trait value
  # Internal nodes: node's own trait value
  # Root branch parent_trait = NA — dropped silently by na.translate = FALSE
  p_base <- ggtree::ggtree(phylo)
  ggtree::`%<+%`(p_base, node_data) +
    ggplot2::aes(color = parent_trait) +
    ggtree::geom_tippoint(
      ggplot2::aes(color = .data[[trait]], size = conf_level, alpha = conf_level)
    ) +
    ggtree::geom_nodepoint(
      ggplot2::aes(color = .data[[trait]], size = conf_level, alpha = conf_level)
    ) +
    ggplot2::scale_color_manual(
      values       = colors,
      name         = trait,
      na.translate = FALSE
    ) +
    ggplot2::scale_size_manual(
      values = c("0-25%" = 1, "25-50%" = 2, "50-75%" = 3, "75-100%" = 4),
      name   = "Confidence",
      na.translate = FALSE
    ) +
    ggplot2::scale_alpha_manual(
      values = c("0-25%" = 0.25, "25-50%" = 0.5, "50-75%" = 0.75, "75-100%" = 1.0),
      name   = "Confidence",
      na.translate = FALSE
    ) +
    ggplot2::theme_classic()
}
