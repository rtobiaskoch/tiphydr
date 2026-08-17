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
#' Points are drawn Nextstrain-style: shape 21 with an independent border and
#' fill. The border is a flat, fully-opaque light grey so points read as
#' distinct circles against the (alpha-dimmed) branches behind them; fill
#' carries the trait colour. If a per-node confidence column is present it is
#' double-encoded as point size and fill alpha (larger, more opaque = higher
#' confidence) — the alpha is baked directly into the fill colour rather than
#' mapped as its own aesthetic, so it dims the fill only, not the border. When
#' confidence is absent, fill is a flat alpha of 0.6 and only size is plain —
#' so the function works on any \code{tbl_tree}, not just the output of
#' \code{\link{build_tree_df}}.
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
#'   \code{"Dark2"}. Ignored for continuous traits, and ignored when
#'   \code{named_palette} is supplied.
#' @param named_palette Optional named character vector (trait level -> hex
#'   colour) for discrete traits, e.g. a project-stable palette instead of the
#'   auto-generated Dark2-by-appearance-order one. Ignored for continuous
#'   traits. Passing your own palette here -- rather than appending
#'   \code{+ scale_color_manual(values = ...)} after the call -- matters
#'   because branches AND point fills need to agree: fills are LITERAL
#'   per-row colours baked in during this call (see the point styling
#'   paragraph above), so a scale added afterward only repaints the branches,
#'   leaving points on the old palette.
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
                            branch_state = c("parent", "node"),
                            named_palette = NULL) {
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

  # --- Colour lookup for the trait, reused below to bake literal point fills -
  # brewer.pal minimum is 3; slice to the actual number of levels. na.omit:
  # root node has NA parent_trait -- not a real level.
  if (is_continuous) {
    trait_range <- range(node_data[[trait]], na.rm = TRUE)
    trait_color_of <- function(x) {
      scales::gradient_n_pal(viridisLite::viridis(256))(
        scales::rescale(x, from = trait_range)
      )
    }
  } else if (!is.null(named_palette)) {
    colors <- named_palette
    trait_color_of <- function(x) unname(colors[as.character(x)])
  } else {
    trait_levels <- unique(stats::na.omit(node_data[[trait]]))
    n_colors <- max(3L, length(trait_levels))
    colors <- RColorBrewer::brewer.pal(n_colors, palette)[seq_along(trait_levels)]
    names(colors) <- trait_levels
    trait_color_of <- function(x) unname(colors[as.character(x)])
  }

  # --- Point aesthetics: Nextstrain-style shape 21, border/fill independent -
  # Fill is a LITERAL per-row colour (not a mapped + scaled aesthetic), with
  # confidence (when present) baked directly into its alpha channel -- ggplot's
  # own `alpha` aesthetic would dim the border too, since it is not
  # fill-specific on shape 21. The confidence -> alpha values match the
  # quartile breakpoints the old scale_alpha_manual() used exactly. Border is
  # set to a flat, fully-opaque light grey as a literal geom param below (not
  # here), independent of the trait or confidence.
  if (has_confidence) {
    conf_alpha <- c(
      "0-25%" = 0.25, "25-50%" = 0.5, "50-75%" = 0.75, "75-100%" = 1.0
    )
    node_data$.fill <- scales::alpha(
      trait_color_of(node_data[[trait]]),
      unname(conf_alpha[as.character(node_data$conf_level)])
    )
    point_aes <- ggplot2::aes(fill = I(.data$.fill), size = .data$conf_level)
  } else {
    node_data$.fill <- scales::alpha(trait_color_of(node_data[[trait]]), 0.6)
    point_aes <- ggplot2::aes(fill = I(.data$.fill))
  }

  # --- Convert to phylo for ggtree (annotations reattached via %<+%) ---------
  phylo <- tidytree::as.phylo(tree_df)

  # Which column paints the branches -- see the branch_state @param.
  branch_col <- if (branch_state == "parent") "parent_trait" else trait

  # Branch alpha 0.6 matches the point fill's flat-case alpha, dimming
  # branches so the (fully-opaque-bordered) points stand out against them.
  p <- ggtree::ggtree(phylo, alpha = 0.6)
  p <- ggtree::`%<+%`(p, node_data) +
    ggplot2::aes(color = .data[[branch_col]]) +
    ggtree::geom_tippoint(point_aes, shape = 21, colour = "grey90", stroke = 0.3) +
    ggtree::geom_nodepoint(point_aes, shape = 21, colour = "grey90", stroke = 0.3)

  # --- Colour scale for branches (parent_trait): viridis / brewer -----------
  # Points no longer use this scale directly (their fill is literal, baked
  # from the same trait_color_of() lookup above), but branches still do, and
  # that is what provides the trait legend.
  if (is_continuous) {
    p <- p + ggplot2::scale_color_viridis_c(name = trait, na.value = "grey80")
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = colors, name = trait, na.translate = FALSE
    )
  }

  # --- Confidence size scale (only if used) ----------------------------------
  if (has_confidence) {
    p <- p +
      ggplot2::scale_size_manual(
        values = c("0-25%" = 1, "25-50%" = 2, "50-75%" = 3, "75-100%" = 4),
        name = "Confidence", na.translate = FALSE
      )
  }

  p + ggplot2::theme_classic()
}
